"""
m_ngsolve.py: Motor NGSolve (nucleo C++ con interfaz Python; LGPL-2.1).

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


def _techo_z(ng, fes, p, esq):
    """GDL z de los VERTICES del techo.

    El espacio H1 de NGSolve es JERARQUICO: con orden 2 cada arista lleva
    una funcion burbuja que vale cero en los vertices, y solo las funciones
    de vertice forman particion de la unidad. De ahi dos reglas que no son las
    de un espacio lagrangiano nodal:
      * el plato impone el valor en los GDL de vertice; los de burbuja del
        techo quedan fijos a CERO (u_z constante en el plano del techo);
      * la reaccion es la suma de las fuerzas internas en los GDL de vertice
        (trabajo virtual con v = e_z en el techo, cuyo interpolante tiene
        burbujas nulas). Sumar tambien las burbujas daba 0,035 N en lugar de
        0,04 N en un bloque TET10 (prueba del bloque 29).
    Con orden 1 todos los GDL son de vertice y la regla coincide con la
    nodal.
    """
    vt = np.searchsorted(esq, np.intersect1d(p["techo_nodos"], esq))
    d0 = _vertice(fes, ng, 0)
    return vt + d0[0] + 2 * (d0[1] - d0[0])


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
    """Directo por omision: el Cholesky disperso de NGSolve.

    MEDIDO (comparativa_motores/, 4 hilos): `sparsecholesky` fue mas rapido
    y uso menos de la mitad de memoria que PARDISO de MKL en todos los casos:
    hex8 a 48^3 5,3 s / 0,7 GB frente a 7,6 s / 1,5 GB; a 80^3 (662 000 GDL)
    33 s / 3,4 GB frente a 52 s / 7,9 GB; TET10 a 48^3 (791 000 GDL) 45 s /
    3,2 GB frente a 69 s / 8,3 GB. Por eso MKL no hace falta; PARDISO queda
    como opcion con SPINPY_PARDISO=1.
    """
    if os.environ.get("SPINPY_PARDISO") and not os.environ.get(
            "SPINPY_SIN_PARDISO"):
        try:
            return mat.Inverse(libres, inverse="pardiso"), "PARDISO (MKL)"
        except Exception:
            pass
    return (mat.Inverse(libres, inverse="sparsecholesky"),
            "sparsecholesky (NGSolve)")


def _modo_auto(ng, m, orden, ndof):
    """Resolvedor por omision, con lo medido en comparativa_motores/:

    TET10 (P2): CG + BDDC. En el espinodoide a 32^3 (315 000 GDL) 18 s y
    1,3 GB; el Cholesky, 14 s y 1,0 GB; a 48^3 (791 000 GDL) 59 s y 3,9 GB
    frente a 45 s y 3,2 GB. Se deja BDDC porque su memoria crece mas despacio
    que la del factor al subir la resolucion.
    hex8 (orden 1): el Cholesky si cabe en la mitad de la memoria del equipo
    (a 365 000 GDL, 14 s frente a 70 s del CG + pyamg); si no, el iterativo.
    No lineal: directo (una factorizacion por iteracion de Newton).
    """
    if m["analisis"] != "lineal":
        return "directo"
    if orden == 2:
        return "iterativo"
    from .. import tiempos
    try:
        from ..fem import memoria_equipo_MB
        tot = memoria_equipo_MB()
    except Exception:
        tot = None
    mem = tiempos.memoria_fem("hex8", ndof)
    return "directo" if (tot is None or mem < 0.5 * tot) else "iterativo"


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
    techo_z = _techo_z(ng, fes, p, esq)
    L = np.array(list(libres), bool)
    lam, mu = lame(m["E"], m["nu"])
    u, v = fes.TnT()
    carga = ng.Parameter(0.0)
    gfu = ng.GridFunction(fes)
    modo = m.get("solver", "auto")
    if modo == "auto":
        modo = _modo_auto(ng, m, orden, fes.ndof)
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
        # `Compile()` reordena el arbol de expresiones de la energia (sin
        # compilador de C: no es `realcompile`). MEDIDO en el espinodoide a
        # 24^3, cuatro pasos SVK: 16,0 s sin compilar y 8,3 s compilado, con
        # las mismas fuerzas a todas las cifras.
        a += ng.Variation(_energia(ng, m["material"],
                                   ng.Id(3) + ng.Grad(u), lam, mu).Compile()
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
