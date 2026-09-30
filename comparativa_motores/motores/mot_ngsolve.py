"""
mot_ngsolve.py — Adaptador de NGSolve (nucleo en C++, interfaz Python).

La malla de spinpy se pasa a netgen tal cual (`AddPoints`, `AddElements`):
hexaedros de un voxel u hojas TET10 por sus esquinas, y las caras de base y
techo como elementos de borde con nombre. Las formas se escriben de modo
simbolico y NGSolve las integra y deriva:

  lineal     sigma(u) : eps(v) dx = -p v_z ds(techo);
  no lineal  Variation(W(F) dx), con W de St. Venant-Kirchhoff o el
             neo-Hookeano de FEBio; NGSolve linealiza la energia
             automaticamente (AssembleLinearization) y el Newton es de
             aqui (cinco lineas), con control por desplazamiento del plato.

Resolvedores: directo = PARDISO de MKL (`pip install mkl`; sin MKL,
`sparsecholesky` propio de NGSolve). Iterativo: con P2, CG + BDDC, el
precondicionador de dominio de NGSolve; con hex8 (orden 1) BDDC degenera en
un directo, asi que se usa CG + pyamg con modos rigidos sobre la matriz que
monta NGSolve.
"""

from __future__ import annotations

import numpy as np
import ngsolve as ng
from netgen.meshing import FaceDescriptor
from netgen.meshing import Mesh as NGMesh

from ._base import (Cronometro, modos_rigidos_gdl, principal,
                    resolver_scipy)

ng.SetNumThreads(int(__import__("os").environ.get("OMP_NUM_THREADS", "4")))


def _malla(p):
    nodos, elems = p["nodos"], p["elems"]
    if elems.shape[1] == 8:
        esq, e3, orden = np.arange(nodos.shape[0]), elems, 1
    else:
        esq = np.unique(elems[:, :4])
        e3, orden = np.searchsorted(esq, elems[:, :4]), 2
    cb, ct = p["caras_base"], p["caras_techo"]
    nv = 4 if elems.shape[1] == 8 else 3
    cb = np.searchsorted(esq, cb[:, :nv])
    ct = np.searchsorted(esq, ct[:, :nv])
    m = NGMesh(dim=3)
    m.AddPoints(np.ascontiguousarray(nodos[esq]))
    m.AddElements(dim=3, index=1,
                  data=np.ascontiguousarray(e3.astype(np.int32)), base=0)
    for k, (nom, caras) in enumerate((("base", cb), ("techo", ct)), start=1):
        m.Add(FaceDescriptor(surfnr=k, domin=1, domout=0, bc=k))
        m.SetBCName(k - 1, nom)
        m.AddElements(dim=2, index=k,
                      data=np.ascontiguousarray(caras.astype(np.int32)),
                      base=0)
    return ng.Mesh(m), esq, orden


def _medios(p, esq):
    n = p["nodos"].shape[0]
    return np.setdiff1d(np.arange(n), esq)


def _dx(orden):
    # hex8: Gauss 2x2x2, como la app y FEBio. tet10: la regla por omision.
    if orden == 1:
        return ng.dx(intrules={ng.HEX: ng.IntegrationRule(ng.HEX, 3)})
    return ng.dx


def _espacio(p, mesh, orden, extra_z=""):
    m = p["meta"]
    if m["apoyo"] == "empotrado":
        fes = ng.VectorH1(mesh, order=orden, dirichlet="base",
                          dirichletz=extra_z)
    else:
        fes = ng.VectorH1(mesh, order=orden,
                          dirichletz="base" + ("|" + extra_z if extra_z
                                               else ""))
    libres = fes.FreeDofs()
    if m["apoyo"] != "empotrado":
        esq = p["_esq"]
        for nodo, comps in ((m["ancla_xy"], (0, 1)), (m["ancla_y"], (1,))):
            v = int(np.searchsorted(esq, nodo))
            d = fes.GetDofNrs(ng.NodeId(ng.VERTEX, v))
            for c in comps:
                libres.Clear(d[c])
    return fes, libres


def _a_spinpy(p, gfu, mesh, esq):
    """u (N, 3) en el orden de spinpy: vertices por GDL, intermedios por
    evaluacion puntual (el campo P2 es continuo; en el punto medio de una
    arista vale lo que su GDL de arista representa)."""
    n = p["nodos"].shape[0]
    fes = gfu.space
    u = np.zeros((n, 3))
    vec = gfu.vec.FV().NumPy()
    nv = mesh.nv
    d = np.array([fes.GetDofNrs(ng.NodeId(ng.VERTEX, 0))])
    paso = d[0, 1] - d[0, 0]           # separacion entre componentes
    base = np.arange(nv) + d[0, 0]
    for c in range(3):
        u[esq, c] = vec[base + c * paso]
    med = _medios(p, esq)
    if med.size:
        X = p["nodos"][med]
        pts = mesh(X[:, 0], X[:, 1], X[:, 2])
        u[med] = gfu(pts)
    return u


def _inversa(a, fes, libres, m, orden, pre=None):
    modo = m.get("solver", "iterativo")
    if modo == "directo":
        try:
            return a.mat.Inverse(libres, inverse="pardiso"), \
                "PARDISO (MKL) via NGSolve"
        except Exception:
            return a.mat.Inverse(libres, inverse="sparsecholesky"), \
                "sparsecholesky (NGSolve)"
    return None, None


def resolver(p):
    m = p["meta"]
    t = Cronometro()
    t("malla")
    mesh, esq, orden = _malla(p)
    p["_esq"] = esq
    if m["analisis"] == "nl":
        return _no_lineal(p, t, mesh, esq, orden)
    if m["analisis"] != "lineal":
        raise NotImplementedError(m["analisis"])
    fes, libres = _espacio(p, mesh, orden)
    u, v = fes.TnT()
    E, nu = m["E"], m["nu"]
    mu, lam = E / 2 / (1 + nu), E * nu / (1 + nu) / (1 - 2 * nu)

    def eps(w):
        return ng.Sym(ng.Grad(w))

    def sig(w):
        return 2 * mu * eps(w) + lam * ng.Trace(eps(w)) * ng.Id(3)

    t("montaje")
    a = ng.BilinearForm(fes, symmetric=True)
    a += ng.InnerProduct(sig(u), eps(v)) * _dx(orden)
    iterativo = m.get("solver", "iterativo") == "iterativo"
    pre = None
    if iterativo and orden == 2:
        pre = ng.Preconditioner(a, "bddc")
    a.Assemble()
    f = ng.LinearForm(fes)
    f += -m["p"] * v[2] * ng.ds("techo")
    f.Assemble()
    gfu = ng.GridFunction(fes)
    t("solucion")
    info = {}
    if not iterativo:
        inv, nombre = _inversa(a, fes, libres, m, orden)
        gfu.vec.data = inv * f.vec
    elif orden == 2:
        inv = ng.solvers.CGSolver(a.mat, pre.mat, tol=m.get("tol", 1e-10),
                                  maxiter=2000, printrates=False)
        gfu.vec.data = inv * f.vec
        nombre = "CG + BDDC (NGSolve)"
        info["iteraciones"] = int(inv.iterations)
    else:
        import scipy.sparse as sp
        fila, col, val = a.mat.COO()
        A = sp.csr_matrix((val.NumPy(), (fila.NumPy(), col.NumPy())),
                          shape=(a.mat.height, a.mat.width))
        I = np.nonzero(np.array(list(libres), bool))[0]
        xyz, comp = _xyz_comp(p, fes, mesh, esq)
        x, info = resolver_scipy(A[I][:, I], f.vec.FV().NumPy()[I],
                                 B=modos_rigidos_gdl(xyz[I], comp[I]),
                                 tol=m.get("tol", 1e-10))
        gfu.vec.FV().NumPy()[I] = x
        nombre = "CG + pyamg (matriz de NGSolve)"
    if "residuo_rel" not in info:
        r = f.vec.CreateVector()
        r.data = f.vec - a.mat * gfu.vec
        rn = np.asarray(r.FV().NumPy())[np.array(list(libres), bool)]
        fn = np.asarray(f.vec.FV().NumPy())[np.array(list(libres), bool)]
        info["residuo_rel"] = float(np.linalg.norm(rn) / np.linalg.norm(fn))
    info["solver"] = nombre
    t("post")
    uu = _a_spinpy(p, gfu, mesh, esq)
    tiempos = t.fin()
    return uu, {"motor": "ngsolve", **info, "tiempos": tiempos,
                "n_gdl": int(fes.ndof)}, None


def _xyz_comp(p, fes, mesh, esq):
    """Posicion y componente de cada GDL (solo orden 1: todos en vertices)."""
    n = fes.ndof
    xyz = np.zeros((n, 3))
    comp = np.zeros(n, np.int64)
    d0 = fes.GetDofNrs(ng.NodeId(ng.VERTEX, 0))
    paso = d0[1] - d0[0]
    base = np.arange(mesh.nv) + d0[0]
    for c in range(3):
        xyz[base + c * paso] = p["nodos"][esq]
        comp[base + c * paso] = c
    return xyz, comp


def _energia(material, F, lam, mu):
    Id = ng.Id(3)
    C = F.trans * F
    if material == "svk":
        Eg = 0.5 * (C - Id)
        return lam / 2 * ng.Trace(Eg) ** 2 + mu * ng.InnerProduct(Eg, Eg)
    J = ng.Det(F)
    return mu / 2 * (ng.Trace(C) - 3) - mu * ng.log(J) \
        + lam / 2 * ng.log(J) ** 2


def _no_lineal(p, t, mesh, esq, orden):
    m = p["meta"]
    E, nu = m["E"], m["nu"]
    mu, lam = E / 2 / (1 + nu), E * nu / (1 + nu) / (1 - 2 * nu)
    fes, libres = _espacio(p, mesh, orden, extra_z="techo")
    u, _v = fes.TnT()
    t("montaje")
    a = ng.BilinearForm(fes, symmetric=True)
    a += ng.Variation(_energia(m["material"], ng.Id(3) + ng.Grad(u), lam, mu)
                      * _dx(orden))
    gfu = ng.GridFunction(fes)
    todos = fes.GetDofs(mesh.Boundaries("techo"))
    d0 = fes.GetDofNrs(ng.NodeId(ng.VERTEX, 0))
    paso = d0[1] - d0[0]
    ndof_c = fes.ndof // 3
    techo = np.nonzero(np.array(list(todos), bool))[0]
    techo_z = techo[(techo // ndof_c) == 2] if paso == ndof_c else \
        techo[(techo - d0[0]) % paso == 2]
    L = np.array(list(libres), bool)
    res = gfu.vec.CreateVector()
    du = gfu.vec.CreateVector()
    fuerzas, iters = [], []
    t("solucion")
    dd = gfu.vec.CreateVector()
    Kdd = gfu.vec.CreateVector()
    for eps in m["eps_plato"]:
        # Predictor consistente (ver mot_skfem): el incremento del plato
        # entra por el sistema linealizado.
        dd[:] = 0.0
        dd.FV().NumPy()[techo_z] = -eps * m["H"] - \
            gfu.vec.FV().NumPy()[techo_z]
        for k in range(m.get("max_iter", 30)):
            a.Apply(gfu.vec, res)
            r = res.FV().NumPy()
            a.AssembleLinearization(gfu.vec)
            if k == 0:
                Kdd.data = a.mat * dd
                res.data += Kdd
                gfu.vec.data += dd
            else:
                ref = max(np.linalg.norm(r[techo_z]), 1e-300)
                if np.linalg.norm(r[L]) / ref < m.get("tol_newton", 1e-10):
                    break
            inv = a.mat.Inverse(libres, inverse="pardiso")
            du.data = inv * res
            gfu.vec.data -= du
        else:
            raise RuntimeError(f"Newton no convergio a eps={eps}")
        iters.append(k)
        fuerzas.append(float(-r[techo_z].sum()))
    t("post")
    uu = _a_spinpy(p, gfu, mesh, esq)
    tiempos = t.fin()
    return uu, {"motor": "ngsolve", "solver": "Newton (propio) + PARDISO",
                "tiempos": tiempos, "iteraciones_newton": iters,
                "n_gdl": int(fes.ndof)}, {"F_reac": np.array(fuerzas)}


if __name__ == "__main__":
    principal(resolver)
