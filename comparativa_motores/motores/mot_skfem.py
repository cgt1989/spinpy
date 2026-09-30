"""
skfem.py — Adaptador de scikit-fem (Python puro sobre numpy/scipy).

scikit-fem monta formas bilineales y lineales escritas en Python y devuelve
matrices de scipy; no trae resolvedor ni Newton propios. Por eso aqui:
  * lineal: montaje con `linear_elasticity` y resolucion con `resolver_scipy`
    (SuperLU o CG + pyamg con modos rigidos, la receta de la app);
  * no lineal: formas de Lagrange total (St. Venant-Kirchhoff y neo-Hookeano
    de FEBio) escritas a mano, con Newton-Raphson y control por
    desplazamiento del plato. Las formulas estan en el INFORME, seccion 2.

Elementos: hex8 = `ElementHex1` (trilineal, Gauss 2x2x2); tet10 =
`ElementTetP2` sobre las esquinas (aristas rectas: los nodos intermedios son
los de spinpy).
"""

from __future__ import annotations

import numpy as np
from skfem import (Basis, BilinearForm, ElementHex1, ElementTetP2,
                   ElementVector, FacetBasis, LinearForm, MeshHex, MeshTet,
                   asm, condense)
from skfem.helpers import ddot, det, eye, grad, inv, transpose
from skfem.models.elasticity import lame_parameters, linear_elasticity

from ._base import Cronometro, modos_rigidos_gdl, principal, resolver_scipy

# C3D8 (spinpy) -> orden de vertices de MeshHex de scikit-fem
PERM_HEX = [0, 4, 3, 1, 7, 5, 2, 6]


def _malla(p):
    nodos, elems = p["nodos"], p["elems"]
    if elems.shape[1] == 8:
        return MeshHex(nodos.T.copy(), elems[:, PERM_HEX].T.copy()), \
            np.arange(nodos.shape[0]), ElementHex1(), 3
    esq = np.unique(elems[:, :4])
    t = np.searchsorted(esq, elems[:, :4])
    return MeshTet(nodos[esq].T.copy(), t.T.copy()), esq, ElementTetP2(), 2


def _mapa_gdl(basis, mesh, esq, nodos):
    """gdl[k, c]: GDL de scikit-fem del nodo k de spinpy, componente c."""
    from ._base import emparejar
    n = nodos.shape[0]
    gdl = -np.ones((n, 3), np.int64)
    gdl[esq, :] = basis.nodal_dofs.T
    if basis.edge_dofs.size:
        med = mesh.p[:, mesh.edges].mean(axis=1).T          # (Ne, 3)
        otros = np.setdiff1d(np.arange(n), esq)
        i = emparejar(med, nodos[otros])
        gdl[otros, :] = basis.edge_dofs[:, i].T
    if (gdl < 0).any():
        raise RuntimeError("GDL sin asignar")
    return gdl


def _preparar(p):
    m = p["meta"]
    t = Cronometro()
    t("malla")
    mesh, esq, el, intorder = _malla(p)
    basis = Basis(mesh, ElementVector(el), intorder=intorder)
    gdl = _mapa_gdl(basis, mesh, esq, p["nodos"])
    z = p["nodos"][:, 2]
    tol = 1e-9 * max(m["H"], 1.0)
    z1 = float(z.max())
    techo_f = mesh.facets_satisfying(lambda x: x[2] > z1 - tol,
                                     boundaries_only=True)
    fijos = []
    base = p["base_nodos"]
    if m["apoyo"] == "empotrado":
        fijos += [gdl[base, c] for c in range(3)]
    else:
        fijos += [gdl[base, 2], gdl[[m["ancla_xy"]], 0],
                  gdl[[m["ancla_xy"]], 1], gdl[[m["ancla_y"]], 1]]
    D = np.unique(np.concatenate(fijos))
    return mesh, basis, gdl, techo_f, D, t


def _xyz_comp(basis, gdl, nodos):
    n = basis.N
    xyz = np.zeros((n, 3))
    comp = np.zeros(n, np.int64)
    for c in range(3):
        xyz[gdl[:, c]] = nodos
        comp[gdl[:, c]] = c
    return xyz, comp


def resolver(p):
    m = p["meta"]
    if m["analisis"] == "nl":
        return _no_lineal(p)
    if m["analisis"] != "lineal":
        raise NotImplementedError(m["analisis"])
    mesh, basis, gdl, techo_f, D, t = _preparar(p)
    t("montaje")
    lam, mu = lame_parameters(m["E"], m["nu"])
    K = asm(linear_elasticity(lam, mu), basis)
    fb = FacetBasis(mesh, basis.elem, facets=techo_f,
                    intorder=4)
    carga = m["p"]

    @LinearForm
    def presion(v, w):
        return -carga * v[2]

    f = asm(presion, fb)
    t("solucion")
    A, b, x, I = condense(K, f, D=D)
    xyz, comp = _xyz_comp(basis, gdl, p["nodos"])
    sol, info = resolver_scipy(A, b, B=modos_rigidos_gdl(xyz[I], comp[I]),
                               modo=m.get("solver", "iterativo"),
                               tol=m.get("tol", 1e-10))
    x[I] = sol
    tiempos = t.fin()
    u = x[gdl]
    return u, {"motor": "skfem", **info, "tiempos": tiempos,
               "n_gdl": int(basis.N)}, None


# ---------------------------------------------------------------------------
# No lineal: Lagrange total, control por desplazamiento del plato
# ---------------------------------------------------------------------------

def mm(A, B):
    """Producto matricial de campos tensoriales (3, 3, elementos, puntos)."""
    return np.einsum("ij...,jk...->ik...", A, B)


def _formas(material, lam, mu):
    def cinematica(w):
        G = grad(w["w"])
        F = G + eye(np.ones_like(G[0, 0]), 3)
        C = mm(transpose(F), F)
        return F, C

    def tension(F, C):
        Id = eye(np.ones_like(C[0, 0]), 3)
        if material == "svk":
            Eg = 0.5 * (C - Id)
            tr = Eg[0, 0] + Eg[1, 1] + Eg[2, 2]
            return lam * tr * Id + 2 * mu * Eg
        Ci = inv(C)
        lnJ = np.log(det(F))
        return mu * (Id - Ci) + lam * lnJ * Ci

    def dS(F, C, dE):
        Id = eye(np.ones_like(C[0, 0]), 3)
        if material == "svk":
            tr = dE[0, 0] + dE[1, 1] + dE[2, 2]
            return lam * tr * Id + 2 * mu * dE
        Ci = inv(C)
        lnJ = np.log(det(F))
        return (2 * (mu - lam * lnJ) * mm(Ci, mm(dE, Ci))
                + lam * ddot(Ci, dE) * Ci)

    @BilinearForm
    def tangente(du, v, w):
        F, C = cinematica(w)
        S = tension(F, C)
        Gu = grad(du)
        dE = 0.5 * (mm(transpose(F), Gu) + mm(transpose(Gu), F))
        dP = mm(Gu, S) + mm(F, dS(F, C, dE))
        return ddot(dP, grad(v))

    @LinearForm
    def interna(v, w):
        F, C = cinematica(w)
        return ddot(mm(F, tension(F, C)), grad(v))

    return tangente, interna


def _no_lineal(p):
    m = p["meta"]
    mesh, basis, gdl, _, D, t = _preparar(p)
    lam, mu = lame_parameters(m["E"], m["nu"])
    tangente, interna = _formas(m["material"], lam, mu)
    techo_z = gdl[p["techo_nodos"], 2]
    D = np.unique(np.concatenate([D, techo_z]))
    I = np.setdiff1d(np.arange(basis.N), D)
    xyz, comp = _xyz_comp(basis, gdl, p["nodos"])
    B = modos_rigidos_gdl(xyz[I], comp[I])
    x = np.zeros(basis.N)
    fuerzas, iters, hist = [], [], []
    t("solucion")
    for eps in m["eps_plato"]:
        # Predictor CONSISTENTE: el incremento del plato entra por el sistema
        # linealizado (K_ff du_f = -r_f - K_fd du_d), no solo en los nodos del
        # techo. Imponerlo solo alli deja una primera iteracion muy
        # distorsionada, y con SVK (energia no convexa) Newton puede acabar en
        # otra rama de equilibrio: medido en el bloque al 20 %.
        dd = np.zeros(basis.N)
        dd[techo_z] = -eps * m["H"] - x[techo_z]
        for k in range(m.get("max_iter", 30)):
            w = basis.interpolate(x)
            r = asm(interna, basis, w=w)
            K = asm(tangente, basis, w=w)
            if k == 0:
                rhs = -(r + K @ dd)[I]
                x += dd
            else:
                nr = np.linalg.norm(r[I])
                ref = max(np.linalg.norm(r[techo_z]), 1e-300)
                hist.append(nr / ref)
                if nr / ref < m.get("tol_newton", 1e-10):
                    break
                rhs = -r[I]
            du, _ = resolver_scipy(K[I][:, I], rhs, B=B,
                                   modo=m.get("solver", "directo"),
                                   tol=1e-12)
            x[I] += du
        else:
            raise RuntimeError(f"Newton no convergio a eps={eps}")
        iters.append(k)
        fuerzas.append(float(-r[techo_z].sum()))
    tiempos = t.fin()
    return x[gdl], {"motor": "skfem", "solver": "Newton (propio) + "
                    + m.get("solver", "directo"), "tiempos": tiempos,
                    "iteraciones_newton": iters, "n_gdl": int(basis.N),
                    "historia_residuo": hist}, {"F_reac": np.array(fuerzas)}


if __name__ == "__main__":
    principal(resolver)
