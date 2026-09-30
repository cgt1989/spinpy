"""
m_skfem.py — Motor scikit-fem (Python puro sobre numpy/scipy; BSD-3).

scikit-fem monta formas bilineales y lineales escritas en Python y devuelve
matrices de scipy; no trae resolvedor ni Newton. Aqui:
  lineal     `linear_elasticity` + traccion en el techo; SuperLU o CG + pyamg
             con modos rigidos (`resolver_scipy`).
  no lineal  formas de Lagrange total escritas a mano (St. Venant-Kirchhoff y
             neo-Hookeano de FEBio) y Newton-Raphson propio con predictor
             consistente del plato.
Elementos: hex8 = ElementHex1 (Gauss 2x2x2); tet10 = ElementTetP2 sobre las
esquinas (aristas rectas).
"""

from __future__ import annotations

import numpy as np

from ._comun import (Cronometro, ErrorMotor, fijos_basicos, lame,
                     modos_rigidos_gdl, resolver_scipy, solucion)

NOMBRE = "scikit-fem"
LICENCIA = "BSD-3-Clause"
# C3D8 (spinpy) -> orden de vertices de MeshHex de scikit-fem
PERM_HEX = [0, 4, 3, 1, 7, 5, 2, 6]


def version():
    import skfem
    return skfem.__version__


def _malla(p):
    from skfem import ElementHex1, ElementTetP2, MeshHex, MeshTet
    nodos, elems = p["nodos"], p["elems"]
    if elems.shape[1] == 8:
        return (MeshHex(nodos.T.copy(), elems[:, PERM_HEX].T.copy()),
                np.arange(nodos.shape[0]), ElementHex1(), 3)
    esq = np.unique(elems[:, :4])
    t = np.searchsorted(esq, elems[:, :4])
    return MeshTet(nodos[esq].T.copy(), t.T.copy()), esq, ElementTetP2(), 2


def _mapa_gdl(basis, mesh, esq, nodos):
    """gdl[k, c]: GDL de scikit-fem del nodo k de spinpy, componente c."""
    from ._comun import emparejar
    n = nodos.shape[0]
    gdl = -np.ones((n, 3), np.int64)
    gdl[esq, :] = basis.nodal_dofs.T
    if basis.edge_dofs.size:
        med = mesh.p[:, mesh.edges].mean(axis=1).T
        otros = np.setdiff1d(np.arange(n), esq)
        gdl[otros, :] = basis.edge_dofs[:, emparejar(med, nodos[otros])].T
    if (gdl < 0).any():
        raise ErrorMotor("GDL sin asignar")
    return gdl


def mm(A, B):
    """Producto matricial de campos tensoriales (3, 3, elementos, puntos)."""
    return np.einsum("ij...,jk...->ik...", A, B)


def _formas(material, lam, mu):
    from skfem import BilinearForm, LinearForm
    from skfem.helpers import ddot, det, eye, grad, inv, transpose

    def cin(w):
        G = grad(w["w"])
        F = G + eye(np.ones_like(G[0, 0]), 3)
        return F, mm(transpose(F), F)

    def tension(F, C):
        Id = eye(np.ones_like(C[0, 0]), 3)
        if material == "svk":
            Eg = 0.5 * (C - Id)
            return lam * (Eg[0, 0] + Eg[1, 1] + Eg[2, 2]) * Id + 2 * mu * Eg
        Ci = inv(C)
        return mu * (Id - Ci) + lam * np.log(det(F)) * Ci

    def dS(F, C, dE):
        Id = eye(np.ones_like(C[0, 0]), 3)
        if material == "svk":
            return lam * (dE[0, 0] + dE[1, 1] + dE[2, 2]) * Id + 2 * mu * dE
        Ci = inv(C)
        lnJ = np.log(det(F))
        return (2 * (mu - lam * lnJ) * mm(Ci, mm(dE, Ci))
                + lam * ddot(Ci, dE) * Ci)

    @BilinearForm
    def tangente(du, v, w):
        F, C = cin(w)
        Gu = grad(du)
        dE = 0.5 * (mm(transpose(F), Gu) + mm(transpose(Gu), F))
        return ddot(mm(Gu, tension(F, C)) + mm(F, dS(F, C, dE)), grad(v))

    @LinearForm
    def interna(v, w):
        F, C = cin(w)
        return ddot(mm(F, tension(F, C)), grad(v))

    return tangente, interna


def resolver(p):
    from skfem import (Basis, ElementVector, FacetBasis, LinearForm, asm)
    from skfem.models.elasticity import linear_elasticity
    m = p["meta"]
    t = Cronometro()
    t("malla")
    mesh, esq, el, intorder = _malla(p)
    basis = Basis(mesh, ElementVector(el), intorder=intorder)
    gdl = _mapa_gdl(basis, mesh, esq, p["nodos"])
    z1 = float(p["nodos"][:, 2].max())
    tol = 1e-9 * max(m["H"], 1.0)
    techo_f = mesh.facets_satisfying(lambda x: x[2] > z1 - tol,
                                     boundaries_only=True)
    D = np.unique([gdl[n, c] for n, c in fijos_basicos(p)])
    techo_z = gdl[p["techo_nodos"], 2]
    plato = m["control"] == "plato"
    if plato:
        D = np.unique(np.concatenate([D, techo_z]))
    I = np.setdiff1d(np.arange(basis.N), D)
    xyz = np.zeros((basis.N, 3))
    comp = np.zeros(basis.N, np.int64)
    for c in range(3):
        xyz[gdl[:, c]] = p["nodos"]
        comp[gdl[:, c]] = c
    B = modos_rigidos_gdl(xyz[I], comp[I])

    @LinearForm
    def presion_unitaria(v, w):
        return -1.0 * v[2]

    f1 = asm(presion_unitaria, FacetBasis(mesh, basis.elem, facets=techo_f,
                                          intorder=4))
    lam, mu = lame(m["E"], m["nu"])
    modo = m.get("solver", "auto")
    x = np.zeros(basis.N)
    info = {}
    if m["analisis"] == "lineal":
        t("montaje")
        K = asm(linear_elasticity(lam, mu), basis)
        t("solucion")
        c = m["cargas"][-1]
        if plato:
            x[techo_z] = -c * m["H"]
            b = -(K @ x)[I]
        else:
            b = c * f1[I]
        if modo == "auto":
            modo = "directo" if I.size <= 6000 else "iterativo"
        x[I], info = resolver_scipy(K[I][:, I].tocsr(), b, B=B, modo=modo,
                                    tol=m.get("tol", 1e-10))
        F = [float(-(K @ x)[techo_z].sum())]
        iters = None
    else:
        tangente, interna = _formas(m["material"], lam, mu)
        if modo == "auto":
            modo = "directo" if I.size <= 60000 else "iterativo"
        F, iters = [], []
        t("solucion")
        for c in m["cargas"]:
            dd = np.zeros(basis.N)
            fext = np.zeros(basis.N)
            if plato:
                dd[techo_z] = -c * m["H"] - x[techo_z]
            else:
                fext = c * f1
            for k in range(m.get("max_iter", 30)):
                w = basis.interpolate(x)
                fint = asm(interna, basis, w=w)
                r = fint - fext
                K = asm(tangente, basis, w=w)
                if k == 0 and plato:
                    rhs = -(r + K @ dd)[I]
                    x += dd
                else:
                    ref = max(np.linalg.norm(fint[techo_z]),
                              np.linalg.norm(fext), 1e-300)
                    if k > 0 and np.linalg.norm(r[I]) / ref < \
                            m.get("tol_newton", 1e-10):
                        break
                    rhs = -r[I]
                du, info = resolver_scipy(K[I][:, I].tocsr(), rhs, B=B,
                                          modo=modo, tol=1e-12)
                x[I] += du
            else:
                raise ErrorMotor(f"Newton no convergio en la carga {c:g}")
            iters.append(k)
            F.append(float(-fint[techo_z].sum()))
        info["solver"] = "Newton (propio) + " + info.get("solver", "")
    tiempos = t.fin()
    return solucion(x[gdl], F, {"motor": "skfem", "tiempos": tiempos,
                                "iteraciones_newton": iters,
                                "n_gdl": int(basis.N), **info})
