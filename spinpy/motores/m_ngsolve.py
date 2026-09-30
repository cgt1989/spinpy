"""
m_ngsolve.py — Motor NGSolve (nucleo C++ con interfaz Python; LGPL-2.1).

La malla de spinpy pasa a netgen sin transformarla (`AddPoints`,
`AddElements`): hexaedros de un voxel, o tetraedros por sus esquinas con
espacio P2, y las caras de base y techo como elementos de borde con nombre.
Las formas son simbolicas:

  lineal     sigma(u) : eps(v) dx = -p v_z ds(techo)
  no lineal  Variation(W(F) dx - p u_z ds(techo)): NGSolve deriva la energia
             y la linealiza (AssembleLinearization); Newton propio con
             predictor consistente del plato.

Resolvedores: directo = PARDISO de MKL (paquete `mkl` de pip; sin el,
`sparsecholesky` de NGSolve). Iterativo: con P2, CG + BDDC (el
precondicionador de subestructuracion de NGSolve); con hex8 (orden 1) BDDC
degenera en un directo, asi que se usa CG + pyamg con modos rigidos sobre la
matriz que monta NGSolve.
"""

from __future__ import annotations

import os

import numpy as np

from ._comun import (Cronometro, ErrorMotor, lame, modos_rigidos_gdl,
                     resolver_scipy, solucion)

NOMBRE = "NGSolve"
LICENCIA = "LGPL-2.1"


def version():
    import ngsolve
    return ngsolve.__version__


def _ng():
    import ngsolve as ng
    ng.SetNumThreads(int(os.environ.get("SPINPY_HILOS",
                                        os.cpu_count() or 1)))
    return ng


def _malla(p):
    ng = _ng()
    from netgen.meshing import FaceDescriptor
    from netgen.meshing import Mesh as NGMesh
    nodos, elems = p["nodos"], p["elems"]
    if elems.shape[1] == 8:
        esq, e3, orden, nv = np.arange(nodos.shape[0]), elems, 1, 4
    else:
        esq = np.unique(elems[:, :4])
        e3, orden, nv = np.searchsorted(esq, elems[:, :4]), 2, 3
    m = NGMesh(dim=3)
    m.AddPoints(np.ascontiguousarray(nodos[esq]))
    m.AddElements(dim=3, index=1,
                  data=np.ascontiguousarray(e3.astype(np.int32)), base=0)
    for k, (nom, caras) in enumerate((("base", p["caras_base"]),
                                      ("techo", p["caras_techo"])), start=1):
        m.Add(FaceDescriptor(surfnr=k, domin=1, domout=0, bc=k))
        m.SetBCName(k - 1, nom)
        c = np.searchsorted(esq, caras[:, :nv]).astype(np.int32)
        m.AddElements(dim=2, index=k, data=np.ascontiguousarray(c), base=0)
    return ng, ng.Mesh(m), esq, orden


def _dx(ng, orden):
    # hex8: Gauss 2x2x2, como la app y FEBio. P2: la regla por omision.
    if orden == 1:
        return ng.dx(intrules={ng.HEX: ng.IntegrationRule(ng.HEX, 3)})
    return ng.dx


def _vertice(fes, ng, v):
    return fes.GetDofNrs(ng.NodeId(ng.VERTEX, int(v)))


def _espacio(ng, p, mesh, orden, esq, plato):
    m = p["meta"]
    z = "base" + ("|techo" if plato else "")
    if m["apoyo"] == "empotrado":
        fes = ng.VectorH1(mesh, order=orden, dirichlet="base",
                          dirichletz="techo" if plato else "")
    else:
        fes = ng.VectorH1(mesh, order=orden, dirichletz=z)
    libres = fes.FreeDofs()
    if m["apoyo"] != "empotrado":
        for nodo, comps in ((m["ancla_xy"], (0, 1)), (m["ancla_y"], (1,))):
            d = _vertice(fes, ng, np.searchsorted(esq, nodo))
            for c in comps:
                libres.Clear(d[c])
    return fes, libres


def _techo_z(ng, fes, mesh):
    todos = np.nonzero(np.array(list(fes.GetDofs(mesh.Boundaries("techo"))),
                                bool))[0]
    d0 = _vertice(fes, ng, 0)
    paso = d0[1] - d0[0]
    if paso == fes.ndof // 3:          # componentes por bloques
        return todos[todos // paso == 2]
    return todos[(todos - d0[0]) % paso == 2]


def _a_spinpy(ng, p, gfu, mesh, esq):
    """u (N, 3) en el orden de spinpy: vertices por GDL; los intermedios
    por evaluacion del campo P2 (continuo) en el punto medio de la arista."""
    fes = gfu.space
    u = np.zeros((p["nodos"].shape[0], 3))
    vec = gfu.vec.FV().NumPy()
    d0 = _vertice(fes, ng, 0)
    paso = d0[1] - d0[0]
    base = np.arange(mesh.nv) + d0[0]
    for c in range(3):
        u[esq, c] = vec[base + c * paso]
    med = np.setdiff1d(np.arange(p["nodos"].shape[0]), esq)
    if med.size:
        X = p["nodos"][med]
        u[med] = gfu(mesh(X[:, 0], X[:, 1], X[:, 2]))
    return u


def _inversa(ng, mat, libres):
    try:
        return mat.Inverse(libres, inverse="pardiso"), "PARDISO (MKL)"
    except Exception:
        return (mat.Inverse(libres, inverse="sparsecholesky"),
                "sparsecholesky (NGSolve)")


def _energia(ng, material, F, lam, mu):
    Id = ng.Id(3)
    C = F.trans * F
    if material == "svk":
        Eg = 0.5 * (C - Id)
        return lam / 2 * ng.Trace(Eg) ** 2 + mu * ng.InnerProduct(Eg, Eg)
    J = ng.Det(F)
    return (mu / 2 * (ng.Trace(C) - 3) - mu * ng.log(J)
            + lam / 2 * ng.log(J) ** 2)


def resolver(p):
    m = p["meta"]
    t = Cronometro()
    t("malla")
    ng, mesh, esq, orden = _malla(p)
    plato = m["control"] == "plato"
    fes, libres = _espacio(ng, p, mesh, orden, esq, plato)
    techo_z = _techo_z(ng, fes, mesh)
    L = np.array(list(libres), bool)
    lam, mu = lame(m["E"], m["nu"])
    u, v = fes.TnT()
    carga = ng.Parameter(0.0)
    gfu = ng.GridFunction(fes)
    modo = m.get("solver", "auto")
    if modo == "auto":
        modo = "iterativo" if (m["analisis"] == "lineal"
                               and fes.ndof > 3e5) else "directo"
    info = {}
    if m["analisis"] == "lineal":
        t("montaje")

        def eps(w):
            return ng.Sym(ng.Grad(w))

        a = ng.BilinearForm(fes, symmetric=True)
        a += ng.InnerProduct(2 * mu * eps(u) + lam * ng.Trace(eps(u))
                             * ng.Id(3), eps(v)) * _dx(ng, orden)
        pre = ng.Preconditioner(a, "bddc") if (modo == "iterativo"
                                               and orden == 2) else None
        a.Assemble()
        f = ng.LinearForm(fes)
        if not plato:
            f += -m["cargas"][-1] * v[2] * ng.ds("techo")
        f.Assemble()
        t("solucion")
        if plato:
            gfu.vec.FV().NumPy()[techo_z] = -m["cargas"][-1] * m["H"]
        r = f.vec.CreateVector()
        r.data = f.vec - a.mat * gfu.vec
        du = gfu.vec.CreateVector()
        if modo == "directo":
            inv, info["solver"] = _inversa(ng, a.mat, libres)
            du.data = inv * r
        elif orden == 2:
            # BDDC proyecta sobre los GDL libres del espacio: `libres` ES ese
            # BitArray (FreeDofs devuelve una referencia), anclas incluidas.
            inv = ng.solvers.CGSolver(a.mat, pre.mat,
                                      tol=m.get("tol", 1e-10),
                                      maxiter=3000, printrates=False)
            du.data = inv * r
            info.update(solver="CG + BDDC (NGSolve)",
                        iteraciones=int(inv.iterations))
        else:
            import scipy.sparse as sp
            fi, co, va = a.mat.COO()
            A = sp.csr_matrix((va.NumPy(), (fi.NumPy(), co.NumPy())),
                              shape=(a.mat.height, a.mat.width))
            I = np.nonzero(L)[0]
            xyz, comp = _xyz_comp(ng, p, fes, mesh, esq)
            x, info = resolver_scipy(A[I][:, I], r.FV().NumPy()[I],
                                     B=modos_rigidos_gdl(xyz[I], comp[I]),
                                     tol=m.get("tol", 1e-10))
            du[:] = 0.0
            du.FV().NumPy()[I] = x
            info["solver"] = "CG + pyamg (matriz de NGSolve)"
        gfu.vec.data += du
        fint = gfu.vec.CreateVector()
        fint.data = a.mat * gfu.vec
        rr = (fint.FV().NumPy() - f.vec.FV().NumPy())[L]
        fn = np.linalg.norm(fint.FV().NumPy()[techo_z])
        info.setdefault("residuo_rel", float(np.linalg.norm(rr) / fn))
        F = [float(-fint.FV().NumPy()[techo_z].sum())]
        iters = None
    else:
        t("montaje")
        a = ng.BilinearForm(fes, symmetric=True)
        a += ng.Variation(_energia(ng, m["material"],
                                   ng.Id(3) + ng.Grad(u), lam, mu)
                          * _dx(ng, orden))
        if not plato:
            a += ng.Variation(carga * u[2] * ng.ds("techo"))
        res = gfu.vec.CreateVector()
        du = gfu.vec.CreateVector()
        dd = gfu.vec.CreateVector()
        Kdd = gfu.vec.CreateVector()
        fint = gfu.vec.CreateVector()
        F, iters = [], []
        t("solucion")
        for c in m["cargas"]:
            dd[:] = 0.0
            if plato:
                dd.FV().NumPy()[techo_z] = -c * m["H"] - \
                    gfu.vec.FV().NumPy()[techo_z]
            else:
                carga.Set(c)
            for k in range(m.get("max_iter", 30)):
                a.Apply(gfu.vec, res)
                r = res.FV().NumPy()
                a.AssembleLinearization(gfu.vec)
                if k == 0 and plato:
                    Kdd.data = a.mat * dd
                    res.data += Kdd
                    gfu.vec.data += dd
                else:
                    ref = max(np.linalg.norm(r[techo_z]),
                              c * m["A_osea_techo"] if not plato else 0.0,
                              1e-300)
                    if k > 0 and np.linalg.norm(r[L]) / ref < \
                            m.get("tol_newton", 1e-10):
                        break
                inv, nombre = _inversa(ng, a.mat, libres)
                du.data = inv * res
                gfu.vec.data -= du
            else:
                raise ErrorMotor(f"Newton no convergio en la carga {c:g}")
            iters.append(k)
            if plato:
                F.append(float(-r[techo_z].sum()))
            else:
                # sin el termino de carga: fuerzas internas en el techo
                carga.Set(0.0)
                a.Apply(gfu.vec, fint)
                carga.Set(c)
                F.append(float(-fint.FV().NumPy()[techo_z].sum()))
        info["solver"] = "Newton (propio) + " + nombre
    t("post")
    uu = _a_spinpy(ng, p, gfu, mesh, esq)
    tiempos = t.fin()
    return solucion(uu, F, {"motor": "ngsolve", "tiempos": tiempos,
                            "iteraciones_newton": iters,
                            "n_gdl": int(fes.ndof), **info})


def _xyz_comp(ng, p, fes, mesh, esq):
    """Posicion y componente de cada GDL (orden 1: todos en vertices)."""
    xyz = np.zeros((fes.ndof, 3))
    comp = np.zeros(fes.ndof, np.int64)
    d0 = _vertice(fes, ng, 0)
    paso = d0[1] - d0[0]
    base = np.arange(mesh.nv) + d0[0]
    for c in range(3):
        xyz[base + c * paso] = p["nodos"][esq]
        comp[base + c * paso] = c
    return xyz, comp
