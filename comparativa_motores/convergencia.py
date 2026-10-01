"""
convergencia.py: Orden de convergencia de cada motor frente a una solucion
exacta, con refinamiento h (malla) y p (grado del polinomio).

    cd comparativa_motores
    python convergencia.py h          # hex8, TET4 y TET10; n = 4, 8, 16, 32
    python convergencia.py p          # TET de grado 1 a 8 sobre una malla fija
    FENICSX_PYTHON=/ruta/al/python/de/fx python convergencia.py h p

El resto de la comparativa mide la diferencia ENTRE motores sobre una misma
malla; eso no dice si la solucion discreta se acerca a la del problema
continuo, ni a que ritmo. Aqui se resuelve un problema con solucion exacta
conocida y se mide el error de discretizacion al refinar. La teoria a priori
(Ciarlet; Brenner y Scott) predice, para una solucion suave y elementos de
grado k,

    |u - u_h|_H1 <= C h^k |u|_H(k+1),     ||u - u_h||_L2 <= C h^(k+1) |u|_H(k+1)

y, con la malla fija y el grado p creciente, un error que baja como
exp(-b p) (Babuska y Suri). En escala logaritmica las dos cosas son rectas:
log(e) frente a log(h) con pendiente k + 1 (L2) o k (H1), y log(e) frente a p.

PROBLEMA. Cubo unidad (mm), E = 20 GPa, nu = 0,30, sin fuerzas de volumen y
con el desplazamiento exacto impuesto en todo el contorno. La solucion exacta
es de Papkovich-Neuber (`exacta`): cumple las ecuaciones de Navier sin
fuerza de volumen, depende de nu y no es polinomica, asi que ningun espacio
de elementos finitos la contiene.

MISMO PROBLEMA DISCRETO. Con refinamiento h, todos los motores reciben la
misma malla (de `malla_cubo`) y los mismos valores nodales en el contorno
(interpolacion lagrangiana; NGSolve, que usa una base jerarquica, convierte
el valor del punto medio de cada arista en el coeficiente de su burbuja). El
error se calcula con el postproceso comun (`errores`), a partir del
desplazamiento nodal de cada motor, con una cuadratura de grado alto. Con
refinamiento p, NGSolve y FEniCSx integran el error con sus propias
herramientas y cada uno impone el contorno con su interpolacion.

Salida: `resultados/convergencia.jsonl`, una linea por (estudio, motor,
elemento, malla, grado). La figura la hace `figuras.py`.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
RES = AQUI / "resultados"
SALIDA = RES / "convergencia.jsonl"

E_MPA, NU = 20000.0, 0.30
#: Numeros de onda de la solucion exacta; c = sqrt(a^2 + b^2) = 2.
A_, B_, D_ = 1.2, 1.6, 1.5
C_ = float(np.hypot(A_, B_))

#: Mallas del refinamiento h (celdas por lado) y del refinamiento p.
N_H = (4, 8, 16, 32)
N_P, GRADOS_P = 4, tuple(range(1, 9))
TOL_ITER = 1e-12
#: Por encima, NGSolve (P2) y FEniCSx pasan de su directo a su iterativo.
GDL_DIRECTO = 4e5


# ---------------------------------------------------------------------------
# Solucion exacta
# ---------------------------------------------------------------------------

def exacta(x, y, z, M, E=E_MPA, nu=NU):
    """Desplazamiento exacto y su gradiente G[i][j] = du_i/dx_j.

    Papkovich-Neuber: 2 mu u = 4 (1 - nu) psi - grad(x . psi + phi), con
    psi = (0, 0, sin(a x) cos(b y) e^(c z)) y phi = e^(d y) cos(d z), ambas
    armonicas (c^2 = a^2 + b^2). Cumple mu lap u + (lam + mu) grad div u = 0
    (comprobado simbolicamente y por diferencias en `comprobar_exacta`).

    `M` es el modulo con sin, cos y exp: numpy, ufl o ngsolve. Asi la misma
    expresion sirve para el postproceso comun y para los motores que integran
    el error por su cuenta (refinamiento p).
    """
    mu2 = E / (1 + nu)                       # 2 mu
    k = 3 - 4 * nu
    a, b, c, d = A_, B_, C_, D_
    s, co = M.sin(a * x), M.cos(a * x)
    sb, cb = M.sin(b * y), M.cos(b * y)
    ez, ey = M.exp(c * z), M.exp(d * y)
    sd, cd = M.sin(d * z), M.cos(d * z)
    u = [-a * z * co * cb * ez / mu2,
         (b * z * s * sb * ez - d * ey * cd) / mu2,
         ((k - c * z) * s * cb * ez + d * ey * sd) / mu2]
    G = [[a * a * z * s * cb * ez / mu2,
          a * b * z * co * sb * ez / mu2,
          -a * co * cb * ez * (1 + c * z) / mu2],
         [a * b * z * co * sb * ez / mu2,
          (b * b * z * s * cb * ez - d * d * ey * cd) / mu2,
          (b * s * sb * ez * (1 + c * z) + d * d * ey * sd) / mu2],
         [a * (k - c * z) * co * cb * ez / mu2,
          (-b * (k - c * z) * s * sb * ez + d * d * ey * sd) / mu2,
          ((c * (k - c * z) - c) * s * cb * ez + d * d * ey * cd) / mu2]]
    return u, G


def u_exacta(X):
    """u exacta (N, 3) en los puntos X (N, 3)."""
    u, _ = exacta(X[:, 0], X[:, 1], X[:, 2], np)
    return np.stack(u, axis=1)


def comprobar_exacta(n=200, h=1e-4, semilla=0):
    """Residuo de Navier y del gradiente por diferencias centradas.

    Devuelve (max |div sigma| / max |sigma|/L, error relativo del gradiente);
    los dos deben ser del orden del truncamiento de la diferencia (~1e-7).
    """
    lam = E_MPA * NU / ((1 + NU) * (1 - 2 * NU))
    mu = E_MPA / (2 * (1 + NU))
    X = np.random.default_rng(semilla).uniform(0.1, 0.9, (n, 3))

    def grad(Y):
        _, G = exacta(Y[:, 0], Y[:, 1], Y[:, 2], np)
        return np.array(G).transpose(2, 0, 1)            # (n, 3, 3)

    def sigma(Y):
        G = grad(Y)
        eps = 0.5 * (G + G.transpose(0, 2, 1))
        tr = np.trace(eps, axis1=1, axis2=2)
        return lam * tr[:, None, None] * np.eye(3) + 2 * mu * eps

    div = np.zeros((n, 3))
    Gfd = np.zeros((n, 3, 3))
    for j in range(3):
        e = np.zeros(3)
        e[j] = h
        div += (sigma(X + e)[:, :, j] - sigma(X - e)[:, :, j]) / (2 * h)
        Gfd[:, :, j] = (u_exacta(X + e) - u_exacta(X - e)) / (2 * h)
    r_nav = float(np.abs(div).max() / np.abs(sigma(X)).max())
    r_G = float(np.abs(Gfd - grad(X)).max() / np.abs(grad(X)).max())
    return r_nav, r_G


# ---------------------------------------------------------------------------
# Mallas estructuradas del cubo unidad
# ---------------------------------------------------------------------------

#: Esquinas de un hexaedro en el orden C3D8 (el de la app).
ESQUINAS = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                     [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
#: Aristas de un TET10 en el orden C3D10.
ARISTAS_TET = ((0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3))


def malla_cubo(n, tipo):
    """Malla de [0, 1]^3 con n celdas por lado. Devuelve (nodos, elems).

    hex8   n^3 hexaedros C3D8.
    tet4   cada celda en seis tetraedros de Kuhn (los seis caminos de la
           esquina (0,0,0) a la (1,1,1)); la particion es conforme porque
           todas las celdas se cortan igual.
    tet10  la de tet4 con un nodo en el punto medio de cada arista (C3D10).
    """
    m = n + 1
    ii, jj, kk = np.meshgrid(np.arange(m), np.arange(m), np.arange(m),
                             indexing="ij")
    ijk = np.stack([ii.ravel(order="F"), jj.ravel(order="F"),
                    kk.ravel(order="F")], axis=1)
    nodos = ijk / float(n)

    def nod(c):
        return c[..., 0] + m * c[..., 1] + m * m * c[..., 2]

    ci, cj, ck = np.meshgrid(np.arange(n), np.arange(n), np.arange(n),
                             indexing="ij")
    celda = np.stack([ci.ravel(order="F"), cj.ravel(order="F"),
                      ck.ravel(order="F")], axis=1)
    if tipo == "hex8":
        return nodos, nod(celda[:, None, :] + ESQUINAS[None, :, :])
    tets = []
    for perm in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1),
                 (2, 1, 0)):
        v = [np.zeros(3, int)]
        for eje in perm:
            w = v[-1].copy()
            w[eje] = 1
            v.append(w)
        tets.append(nod(celda[:, None, :] + np.array(v)[None, :, :]))
    elems = np.concatenate(tets, axis=0)
    x = nodos[elems]
    vol = np.einsum("ij,ij->i", np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0]),
                    x[:, 3] - x[:, 0])
    neg = vol < 0
    elems[neg, 1], elems[neg, 2] = elems[neg, 2], elems[neg, 1].copy()
    if tipo == "tet4":
        return nodos, elems
    pares = np.concatenate([np.sort(elems[:, list(a)], axis=1)
                            for a in ARISTAS_TET], axis=0)
    unicas, inv = np.unique(pares, axis=0, return_inverse=True)
    medios = 0.5 * (nodos[unicas[:, 0]] + nodos[unicas[:, 1]])
    medio_e = inv.reshape(6, -1).T + nodos.shape[0]
    return np.vstack([nodos, medios]), np.hstack([elems, medio_e])


def en_contorno(nodos, tol=1e-12):
    return np.any((nodos < tol) | (nodos > 1 - tol), axis=1)


def caras_contorno(elems):
    """Caras exteriores (esquinas) con la normal hacia fuera."""
    if elems.shape[1] == 8:
        loc = ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5),
               (2, 3, 7, 6), (3, 0, 4, 7))
    else:
        loc = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))
    caras = np.concatenate([elems[:, list(f)] for f in loc], axis=0)
    _, idx, cuenta = np.unique(np.sort(caras, axis=1), axis=0,
                               return_index=True, return_counts=True)
    return caras[idx[cuenta == 1]]


def emparejar(coord_motor, coord_ref, tol=1e-9):
    """Indice i tal que coord_motor[i[k]] == coord_ref[k]."""
    from scipy.spatial import cKDTree
    d, i = cKDTree(coord_motor).query(coord_ref)
    if d.max() > tol:
        raise RuntimeError(f"nodos sin pareja (distancia {d.max():.2e})")
    return i


def lame(E=E_MPA, nu=NU):
    return E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))


# ---------------------------------------------------------------------------
# Postproceso comun: error en L2 y en la seminorma H1
# ---------------------------------------------------------------------------

def _gauss01(q):
    x, w = np.polynomial.legendre.leggauss(q)
    return 0.5 * (x + 1), 0.5 * w


def _cuadratura(tipo, q=6):
    """Puntos (nq, 3) y pesos del elemento de referencia.

    hex: Gauss-Legendre q^3 en [-1, 1]^3. tet: Gauss colapsado (Duffy) sobre
    el tetraedro unidad, exacto para polinomios de grado 2q - 3.
    """
    x, w = _gauss01(q)
    U, V, W = np.meshgrid(x, x, x, indexing="ij")
    P = np.meshgrid(w, w, w, indexing="ij")
    peso = (P[0] * P[1] * P[2]).ravel()
    if tipo == "hex8":
        return np.stack([2 * U.ravel() - 1, 2 * V.ravel() - 1,
                         2 * W.ravel() - 1], axis=1), 8 * peso
    u, v, w_ = U.ravel(), V.ravel(), W.ravel()
    xi = np.stack([u, v * (1 - u), w_ * (1 - u) * (1 - v)], axis=1)
    return xi, peso * (1 - u) ** 2 * (1 - v)


def _forma(tipo, xi):
    """N (nq, nn) y dN/dxi (nq, nn, 3) en los puntos `xi`."""
    if tipo == "hex8":
        s = 2 * ESQUINAS - 1
        f = 1 + xi[:, None, :] * s[None, :, :]            # (nq, 8, 3)
        N = f.prod(axis=2) / 8
        dN = np.stack([s[None, :, 0] * f[..., 1] * f[..., 2],
                       s[None, :, 1] * f[..., 0] * f[..., 2],
                       s[None, :, 2] * f[..., 0] * f[..., 1]], axis=2) / 8
        return N, dN
    L = np.stack([1 - xi.sum(axis=1), xi[:, 0], xi[:, 1], xi[:, 2]], axis=1)
    dL = np.array([[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    if tipo == "tet4":
        return L, np.broadcast_to(dL, (xi.shape[0], 4, 3)).copy()
    N = np.zeros((xi.shape[0], 10))
    dN = np.zeros((xi.shape[0], 10, 3))
    for i in range(4):
        N[:, i] = L[:, i] * (2 * L[:, i] - 1)
        dN[:, i] = (4 * L[:, i] - 1)[:, None] * dL[i]
    for k, (i, j) in enumerate(ARISTAS_TET):
        N[:, 4 + k] = 4 * L[:, i] * L[:, j]
        dN[:, 4 + k] = 4 * (L[:, i][:, None] * dL[j] + L[:, j][:, None]
                            * dL[i])
    return N, dN


def errores(nodos, elems, u, tipo, q=6, bloque=4000):
    """Error absoluto y relativo de u (nodal) frente a la solucion exacta.

    Interpola u con las funciones de forma del elemento (isoparametrico) e
    integra |u - u_ex|^2 y |grad u - grad u_ex|^2 con una cuadratura de
    grado alto. Devuelve e_L2, e_H1 (seminorma), y las normas de u_ex.
    """
    xi, w = _cuadratura(tipo, q)
    N, dN = _forma(tipo, xi)
    acc = np.zeros(4)
    for i0 in range(0, elems.shape[0], bloque):
        e = elems[i0:i0 + bloque]
        X = nodos[e]                                     # (ne, nn, 3)
        Ue = u[e]
        xq = np.einsum("qa,eai->eqi", N, X)
        J = np.einsum("eai,qaj->eqij", X, dN)            # dx_i/dxi_j
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)
        dNx = np.einsum("qaj,eqji->eqai", dN, Jinv)      # dN_a/dx_i
        uh = np.einsum("qa,eac->eqc", N, Ue)
        Gh = np.einsum("eac,eqai->eqci", Ue, dNx)
        ue, Ge = exacta(xq[..., 0], xq[..., 1], xq[..., 2], np)
        ue = np.stack(ue, axis=-1)
        Ge = np.stack([np.stack(f, axis=-1) for f in Ge], axis=-2)
        dv = w[None, :] * detJ
        acc += [np.sum(dv * ((uh - ue) ** 2).sum(-1)),
                np.sum(dv * ((Gh - Ge) ** 2).sum((-1, -2))),
                np.sum(dv * (ue ** 2).sum(-1)),
                np.sum(dv * (Ge ** 2).sum((-1, -2)))]
    eL2, eH1, nL2, nH1 = np.sqrt(acc)
    return {"e_L2": float(eL2), "e_H1": float(eH1),
            "e_L2_rel": float(eL2 / nL2), "e_H1_rel": float(eH1 / nH1),
            "norma_L2": float(nL2), "norma_H1": float(nH1)}


# ---------------------------------------------------------------------------
# Motores (refinamiento h). Cada uno devuelve (u nodal, meta).
# ---------------------------------------------------------------------------

def _libres(nodos):
    return np.nonzero(~en_contorno(nodos))[0]


def _solve_scipy(A, b, B, n_dir=20000):
    """SuperLU hasta n_dir incognitas; por encima, CG + pyamg (modos rigidos)
    con tolerancia 1e-12, muy por debajo del menor error de discretizacion."""
    from scipy.sparse.linalg import splu
    A = A.tocsr()
    if A.shape[0] <= n_dir:
        x = splu(A.tocsc()).solve(b)
        nombre = "SuperLU (scipy)"
    else:
        import pyamg
        ml = pyamg.smoothed_aggregation_solver(A, B=B, max_coarse=500)
        x = ml.solve(b, tol=TOL_ITER, maxiter=3000, accel="cg")
        nombre = "CG + pyamg (modos rigidos)"
    r = float(np.linalg.norm(b - A @ x) / np.linalg.norm(b))
    return x, nombre, r


def _modos(xyz, comp):
    sys.path.insert(0, str(AQUI.parent))
    from spinpy.motores._comun import modos_rigidos_gdl
    return modos_rigidos_gdl(xyz, comp)


def motor_app(nodos, elems, uD):
    """Matriz elemental de la app (`elastic.hex8_ke`, Gauss 2x2x2) y su
    receta de resolucion: LU por debajo de 6000 GDL, CG + AMG por encima con
    la tolerancia de la app (1e-8)."""
    sys.path.insert(0, str(AQUI.parent))
    import scipy.sparse as sp
    from scipy.sparse.linalg import splu

    from spinpy.elastic import hex8_ke
    from spinpy.resistencia import UMBRAL_DIRECTO, _modos_rigidos
    if elems.shape[1] != 8:
        raise NotImplementedError("la app solo tiene hex8")
    h = float(nodos[elems[0, 1], 0] - nodos[elems[0, 0], 0])
    ke = hex8_ke(h, h, h, E_MPA, NU)
    edof = (3 * elems[:, :, None] + np.arange(3)).reshape(-1, 24)
    nd = 3 * nodos.shape[0]
    K = sp.coo_matrix((np.tile(ke.ravel(), edof.shape[0]),
                       (np.repeat(edof, 24, axis=1).ravel(),
                        np.tile(edof, (1, 24)).ravel())),
                      shape=(nd, nd)).tocsr()
    lib = _libres(nodos)
    I = (3 * lib[:, None] + np.arange(3)).ravel()
    x = np.zeros(nd)
    fr = np.nonzero(en_contorno(nodos))[0]
    D = (3 * fr[:, None] + np.arange(3)).ravel()
    x[D] = uD[fr].ravel()
    Kii = K[I][:, I].tocsr()
    b = -(K[I][:, D] @ x[D])
    if I.size <= UMBRAL_DIRECTO:
        x[I] = splu(Kii.tocsc()).solve(b)
        nombre = "LU directo (app)"
    else:
        import pyamg
        ml = pyamg.smoothed_aggregation_solver(
            Kii, B=_modos_rigidos(nodos)[I], max_coarse=500)
        x[I] = ml.solve(b, tol=1e-8, maxiter=500, accel="cg")
        nombre = "AMG (con modos rigidos) + CG (app)"
    r = float(np.linalg.norm(b - Kii @ x[I]) / np.linalg.norm(b))
    return x.reshape(-1, 3), {"solver": nombre, "residuo_rel": r,
                              "n_gdl": int(nd)}


def motor_skfem(nodos, elems, uD):
    from skfem import (Basis, ElementHex1, ElementTetP1, ElementTetP2,
                       ElementVector, MeshHex, MeshTet, asm)
    from skfem.models.elasticity import linear_elasticity
    sys.path.insert(0, str(AQUI.parent))
    from spinpy.motores.m_skfem import PERM_HEX, _mapa_gdl
    if elems.shape[1] == 8:
        esq = np.arange(nodos.shape[0])
        mesh = MeshHex(nodos.T.copy(), elems[:, PERM_HEX].T.copy())
        el, io = ElementHex1(), 3
    else:
        esq = np.unique(elems[:, :4])
        t = np.searchsorted(esq, elems[:, :4])
        mesh = MeshTet(nodos[esq].T.copy(), t.T.copy())
        el, io = ((ElementTetP2(), 4) if elems.shape[1] == 10
                  else (ElementTetP1(), 2))
    basis = Basis(mesh, ElementVector(el), intorder=io)
    gdl = _mapa_gdl(basis, mesh, esq, nodos)
    lam, mu = lame()
    K = asm(linear_elasticity(lam, mu), basis)
    fr = np.nonzero(en_contorno(nodos))[0]
    D = gdl[fr].ravel()
    x = np.zeros(basis.N)
    x[D] = uD[fr].ravel()
    I = np.setdiff1d(np.arange(basis.N), D)
    xyz = np.zeros((basis.N, 3))
    comp = np.zeros(basis.N, np.int64)
    for c in range(3):
        xyz[gdl[:, c]] = nodos
        comp[gdl[:, c]] = c
    x[I], nombre, r = _solve_scipy(K[I][:, I], -(K[I][:, D] @ x[D]),
                                   _modos(xyz[I], comp[I]))
    return x[gdl], {"solver": nombre, "residuo_rel": r,
                    "n_gdl": int(basis.N)}


def motor_ngsolve(nodos, elems, uD):
    import ngsolve as ng
    from netgen.meshing import FaceDescriptor
    from netgen.meshing import Mesh as NGMesh
    ng.SetNumThreads(int(os.environ.get("SPINPY_HILOS", os.cpu_count() or 1)))
    hexa = elems.shape[1] == 8
    esq = np.arange(nodos.shape[0]) if hexa else np.unique(elems[:, :4])
    e3 = elems if hexa else np.searchsorted(esq, elems[:, :4])
    orden = 2 if elems.shape[1] == 10 else 1
    m = NGMesh(dim=3)
    m.AddPoints(np.ascontiguousarray(nodos[esq]))
    m.AddElements(dim=3, index=1, data=np.ascontiguousarray(
        e3.astype(np.int32)), base=0)
    m.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    m.SetBCName(0, "borde")
    m.AddElements(dim=2, index=1, data=np.ascontiguousarray(
        caras_contorno(e3).astype(np.int32)), base=0)
    mesh = ng.Mesh(m)
    fes = ng.VectorH1(mesh, order=orden, dirichlet="borde")
    u, v = fes.TnT()
    lam, mu = lame()

    def eps(w):
        return ng.Sym(ng.Grad(w))

    dx = (ng.dx(intrules={ng.HEX: ng.IntegrationRule(ng.HEX, 3)}) if hexa
          else ng.dx)
    a = ng.BilinearForm(fes, symmetric=True)
    a += ng.InnerProduct(2 * mu * eps(u) + lam * ng.Trace(eps(u)) * ng.Id(3),
                         eps(v)) * dx
    # Con TET10 a n = 32 (8e5 GDL) el Cholesky del cubo MACIZO agota los
    # 15 GB (medido: std::bad_alloc; el espinodoide, poroso, tiene mucho
    # menos relleno), y BDDC tambien, al factorizar su problema grueso. Ahi
    # se resuelve la matriz de NGSolve con CG + pyamg y los modos rigidos
    # (la ruta de su adaptador con hex8), a 1e-12.
    iterativo = orden == 2 and fes.ndof > GDL_DIRECTO
    a.Assemble()
    gfu = ng.GridFunction(fes)
    vec = gfu.vec.FV().NumPy()
    d0 = fes.GetDofNrs(ng.NodeId(ng.VERTEX, 0))
    paso = d0[1] - d0[0]
    base = np.arange(mesh.nv) + d0[0]
    fr_v = en_contorno(nodos[esq])
    for c in range(3):
        vec[base[fr_v] + c * paso] = uD[esq[fr_v], c]
    if orden == 2:
        # Base jerarquica: el GDL de arista es el coeficiente de una burbuja
        # que vale beta en el punto medio. Valor nodal del punto medio ->
        # coeficiente = (u_m - (u_a + u_b) / 2) / beta.
        med = np.setdiff1d(np.arange(nodos.shape[0]), esq)
        ed = np.array([[v_.nr for v_ in e.vertices] for e in mesh.edges])
        Xm = 0.5 * (nodos[esq[ed[:, 0]]] + nodos[esq[ed[:, 1]]])
        idx_m = med[emparejar(nodos[med], Xm)]
        dofs_e = np.array([fes.GetDofNrs(ng.NodeId(ng.EDGE, e))
                           for e in range(mesh.nedge)])
        prueba = ng.GridFunction(fes)
        prueba.vec.FV().NumPy()[dofs_e[0, 0]] = 1.0
        beta = float(prueba(mesh(*Xm[0]))[0])
        fr_e = en_contorno(Xm)
        for c in range(3):
            vm = uD[idx_m[fr_e], c]
            va = uD[esq[ed[fr_e, 0]], c]
            vb = uD[esq[ed[fr_e, 1]], c]
            vec[dofs_e[fr_e, c]] = (vm - 0.5 * (va + vb)) / beta
    r = gfu.vec.CreateVector()
    r.data = -1.0 * (a.mat * gfu.vec)
    if iterativo:
        import scipy.sparse as sp
        fi, co, va = a.mat.COO()
        A = sp.csr_matrix((va.NumPy(), (fi.NumPy(), co.NumPy())),
                          shape=(fes.ndof, fes.ndof))
        I = np.nonzero(np.array(list(fes.FreeDofs()), bool))[0]
        # Modos rigidos en la base jerarquica: los campos lineales los
        # representan las funciones de vertice; las burbujas valen cero.
        xyz = np.zeros((fes.ndof, 3))
        comp = np.full(fes.ndof, -1)
        for c in range(3):
            xyz[base + c * paso] = nodos[esq]
            comp[base + c * paso] = c
        B = _modos(xyz, comp)
        x, _, res = _solve_scipy(A[I][:, I], r.FV().NumPy()[I], B[I], n_dir=0)
        vec[I] += x
        nombre = "CG + pyamg (matriz de NGSolve, modos rigidos)"
    else:
        inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        nombre = "sparsecholesky (NGSolve)"
        du = gfu.vec.CreateVector()
        du.data = inv * r
        gfu.vec.data += du
    U = np.zeros((nodos.shape[0], 3))
    for c in range(3):
        U[esq, c] = vec[base + c * paso]
    med = np.setdiff1d(np.arange(nodos.shape[0]), esq)
    if med.size:
        X = nodos[med]
        U[med] = gfu(mesh(X[:, 0], X[:, 1], X[:, 2]))
    # Error integrado por NGSolve, para contrastar el postproceso comun.
    ue, Ge = exacta(ng.x, ng.y, ng.z, ng)
    ucf = ng.CF(tuple(ue))
    Gcf = ng.CF(tuple(g for fila in Ge for g in fila), dims=(3, 3))
    q = 2 * orden + 8
    eL2 = np.sqrt(ng.Integrate(ng.InnerProduct(gfu - ucf, gfu - ucf), mesh,
                               order=q))
    eH1 = np.sqrt(ng.Integrate(ng.InnerProduct(ng.Grad(gfu) - Gcf,
                                               ng.Grad(gfu) - Gcf), mesh,
                               order=q))
    return U, {"solver": nombre, "n_gdl": int(fes.ndof),
               "e_L2_propio": float(eL2), "e_H1_propio": float(eH1)}


def motor_sfepy(nodos, elems, uD):
    import logging
    logging.disable(logging.WARNING)
    from sfepy.base.base import output
    output.set_output(quiet=True)
    import sfepy.discrete as sd
    from sfepy.discrete import FieldVariable, Function
    from sfepy.discrete.conditions import Conditions, EssentialBC
    from sfepy.discrete.fem import FEDomain, Field, Mesh
    from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
    from sfepy.solvers.ls import ScipyDirect
    from sfepy.solvers.nls import Newton
    from sfepy.terms import Term
    sys.path.insert(0, str(AQUI.parent))
    from spinpy.motores.m_sfepy import _solver_lineal
    hexa = elems.shape[1] == 8
    esq = np.arange(nodos.shape[0]) if hexa else np.unique(elems[:, :4])
    conn = elems if hexa else np.searchsorted(esq, elems[:, :4])
    orden = 2 if elems.shape[1] == 10 else 1
    mesh = Mesh.from_data("cubo", nodos[esq], None,
                          [conn.astype(np.int32)],
                          [np.ones(conn.shape[0], np.int32)],
                          ["3_8" if hexa else "3_4"])
    dom = FEDomain("d", mesh)
    omega = dom.create_region("Omega", "all")
    gamma = dom.create_region("Gamma", "vertices of surface", "facet")
    field = Field.from_args("u", np.float64, 3, omega, approx_order=orden)
    u = FieldVariable("u", "unknown", field)
    v = FieldVariable("v", "test", field, primary_var_name="u")
    fun = Function("uD", lambda ts, coors, **kw: u_exacta(coors))
    ebc = EssentialBC("borde", gamma, {"u.all": fun})
    mat = sd.Material("m", D=stiffness_from_youngpoisson(3, E_MPA, NU))
    term = Term.new("dw_lin_elastic(m.D, v, u)",
                    sd.Integral("i", order=2 * orden), omega, m=mat, v=v,
                    u=u)
    pb = sd.Problem("p", equations=sd.Equations([sd.Equation("eq", term)]))
    pb.set_bcs(ebcs=Conditions([ebc]))
    n_gdl = int(field.n_nod * 3)
    if n_gdl <= 20000:
        ls, nombre = ScipyDirect({"method": "superlu"}), \
            "SuperLU (ls.scipy_direct)"
    else:
        ls, nombre = _solver_lineal("iterativo", field, u, TOL_ITER)
    status = {}
    pb.set_solver(Newton({"i_max": 1, "eps_a": 1e-300, "is_linear": True},
                         lin_solver=ls, status=status))
    st = pb.solve(save_results=False)
    U = np.asarray(st()).reshape(-1, 3)[emparejar(field.get_coor(), nodos)]
    return U, {"solver": nombre, "n_gdl": n_gdl}


def motor_fenicsx(nodos, elems, uD):
    """Corre en el entorno de conda-forge (proceso hijo, ver `_hijo`)."""
    import basix.ufl
    import ufl
    from dolfinx import fem, mesh as dmesh
    from dolfinx.fem.petsc import LinearProblem
    from mpi4py import MPI
    hexa = elems.shape[1] == 8
    esq = np.arange(nodos.shape[0]) if hexa else np.unique(elems[:, :4])
    if hexa:
        celda, cells = "hexahedron", elems[:, [0, 1, 3, 2, 4, 5, 7, 6]]
    else:
        celda, cells = "tetrahedron", np.searchsorted(esq, elems[:, :4])
    orden = 2 if elems.shape[1] == 10 else 1
    dom = ufl.Mesh(basix.ufl.element("Lagrange", celda, 1, shape=(3,)))
    msh = dmesh.create_mesh(MPI.COMM_SELF, cells.astype(np.int64), dom,
                            nodos[esq])
    V = fem.functionspace(msh, ("Lagrange", orden, (3,)))
    fd = msh.topology.dim - 1
    msh.topology.create_connectivity(fd, fd + 1)
    caras = dmesh.exterior_facet_indices(msh.topology)
    uD_f = fem.Function(V)
    uD_f.interpolate(lambda X: u_exacta(X.T).T)
    bc = fem.dirichletbc(uD_f, fem.locate_dofs_topological(V, fd, caras))
    lam, mu = lame()
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    def eps(w):
        return ufl.sym(ufl.grad(w))

    a = ufl.inner(2 * mu * eps(u) + lam * ufl.tr(eps(u)) * ufl.Identity(3),
                  eps(v)) * ufl.dx
    L = ufl.inner(fem.Constant(msh, np.zeros(3)), v) * ufl.dx
    n_gdl = V.dofmap.index_map.size_global * 3
    if n_gdl <= GDL_DIRECTO:
        nombre = "Cholesky MUMPS (PETSc)"
        uh = LinearProblem(a, L, bcs=[bc], petsc_options={
            "ksp_type": "preonly", "pc_type": "cholesky",
            "pc_factor_mat_solver_type": "mumps"},
            petsc_options_prefix="conv_").solve()
        if isinstance(uh, tuple):
            uh = uh[0]
    else:
        uh, nombre = _fenicsx_gamg(V, a, L, bc)
    X = V.tabulate_dof_coordinates()
    U = uh.x.array.reshape(-1, 3)[emparejar(X, nodos)]
    return U, {"solver": nombre, "n_gdl": int(n_gdl)}


def _fenicsx_gamg(V, a, L, bc):
    """CG + GAMG con los seis modos rigidos como espacio nulo cercano (la
    receta de `spinpy.motores.m_fenicsx`), residuo no precondicionado."""
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    from mpi4py import MPI
    from petsc4py import PETSc
    A = fp.assemble_matrix(fem.form(a), bcs=[bc])
    A.assemble()
    b = fp.assemble_vector(fem.form(L))
    fp.apply_lifting(b, [fem.form(a)], bcs=[[bc]])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    fp.set_bc(b, [bc])
    X = V.tabulate_dof_coordinates()[:V.dofmap.index_map.size_local]
    vecs = []
    for k in range(6):
        m = np.zeros((X.shape[0], 3))
        if k < 3:
            m[:, k] = 1.0
        else:
            i, j = ((0, 1), (1, 2), (2, 0))[k - 3]
            m[:, i], m[:, j] = -X[:, j], X[:, i]
        w = A.createVecLeft()
        w.setArray(m.ravel())
        for q in vecs:
            w.axpy(-w.dot(q), q)
        w.normalize()
        vecs.append(w)
    A.setNearNullSpace(PETSc.NullSpace().create(vectors=vecs))
    ksp = PETSc.KSP().create(MPI.COMM_SELF)
    ksp.setOperators(A)
    ksp.setType("cg")
    ksp.getPC().setType("gamg")
    ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
    ksp.setTolerances(rtol=TOL_ITER, max_it=5000)
    uh = fem.Function(V)
    ksp.solve(b, uh.x.petsc_vec)
    if ksp.getConvergedReason() <= 0:
        raise RuntimeError(f"GAMG no convergio ({ksp.getConvergedReason()})")
    uh.x.scatter_forward()
    return uh, "CG + GAMG (PETSc, modos rigidos)"


MOTORES_H = {"app": motor_app, "ngsolve": motor_ngsolve,
             "fenicsx": motor_fenicsx, "skfem": motor_skfem,
             "sfepy": motor_sfepy}


# ---------------------------------------------------------------------------
# Refinamiento p (malla fija, grado creciente)
# ---------------------------------------------------------------------------

def p_ngsolve(nodos, elems, grados):
    import ngsolve as ng
    from netgen.meshing import FaceDescriptor
    from netgen.meshing import Mesh as NGMesh
    ng.SetNumThreads(int(os.environ.get("SPINPY_HILOS", os.cpu_count() or 1)))
    m = NGMesh(dim=3)
    m.AddPoints(np.ascontiguousarray(nodos))
    m.AddElements(dim=3, index=1, data=np.ascontiguousarray(
        elems[:, :4].astype(np.int32)), base=0)
    m.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    m.SetBCName(0, "borde")
    m.AddElements(dim=2, index=1, data=np.ascontiguousarray(
        caras_contorno(elems[:, :4]).astype(np.int32)), base=0)
    mesh = ng.Mesh(m)
    ue, Ge = exacta(ng.x, ng.y, ng.z, ng)
    ucf = ng.CF(tuple(ue))
    Gcf = ng.CF(tuple(g for fila in Ge for g in fila), dims=(3, 3))
    lam, mu = lame()
    filas = []
    for p in grados:
        t0 = time.perf_counter()
        fes = ng.VectorH1(mesh, order=p, dirichlet="borde")
        u, v = fes.TnT()
        eps_u, eps_v = ng.Sym(ng.Grad(u)), ng.Sym(ng.Grad(v))
        a = ng.BilinearForm(fes, symmetric=True)
        a += ng.InnerProduct(2 * mu * eps_u + lam * ng.Trace(eps_u)
                             * ng.Id(3), eps_v) * ng.dx
        a.Assemble()
        gfu = ng.GridFunction(fes)
        gfu.Set(ucf, ng.BND)
        r = gfu.vec.CreateVector()
        r.data = -1.0 * (a.mat * gfu.vec)
        gfu.vec.data += a.mat.Inverse(fes.FreeDofs(),
                                      inverse="sparsecholesky") * r
        q = 2 * p + 8
        acc = [ng.Integrate(ng.InnerProduct(gfu - ucf, gfu - ucf), mesh,
                            order=q),
               ng.Integrate(ng.InnerProduct(ng.Grad(gfu) - Gcf,
                                            ng.Grad(gfu) - Gcf), mesh,
                            order=q),
               ng.Integrate(ng.InnerProduct(ucf, ucf), mesh, order=q),
               ng.Integrate(ng.InnerProduct(Gcf, Gcf), mesh, order=q)]
        eL2, eH1, nL2, nH1 = np.sqrt(acc)
        filas.append({"grado": p, "n_gdl": int(fes.ndof),
                      "solver": "sparsecholesky (NGSolve)",
                      "contorno": "Set(BND) de NGSolve",
                      "tiempo_s": time.perf_counter() - t0,
                      "e_L2": float(eL2), "e_H1": float(eH1),
                      "e_L2_rel": float(eL2 / nL2),
                      "e_H1_rel": float(eH1 / nH1)})
    return filas


def p_fenicsx(nodos, elems, grados):
    import basix
    import basix.ufl
    import ufl
    from dolfinx import fem, mesh as dmesh
    from dolfinx.fem.petsc import LinearProblem
    from mpi4py import MPI
    dom = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1,
                                     shape=(3,)))
    msh = dmesh.create_mesh(MPI.COMM_SELF, elems[:, :4].astype(np.int64),
                            dom, nodos)
    fd = msh.topology.dim - 1
    msh.topology.create_connectivity(fd, fd + 1)
    caras = dmesh.exterior_facet_indices(msh.topology)
    x = ufl.SpatialCoordinate(msh)
    ue, Ge = exacta(x[0], x[1], x[2], ufl)
    ucf, Gcf = ufl.as_vector(ue), ufl.as_matrix(Ge)
    lam, mu = lame()
    filas = []
    for p in grados:
        t0 = time.perf_counter()
        el = basix.ufl.element("Lagrange", "tetrahedron", p, shape=(3,),
                               lagrange_variant=basix.LagrangeVariant
                               .gll_warped)
        V = fem.functionspace(msh, el)
        uD = fem.Function(V)
        uD.interpolate(lambda X: u_exacta(X.T).T)
        bc = fem.dirichletbc(uD, fem.locate_dofs_topological(V, fd, caras))
        u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
        eu, ev = ufl.sym(ufl.grad(u)), ufl.sym(ufl.grad(v))
        a = ufl.inner(2 * mu * eu + lam * ufl.tr(eu) * ufl.Identity(3),
                      ev) * ufl.dx
        L = ufl.inner(fem.Constant(msh, np.zeros(3)), v) * ufl.dx
        uh = LinearProblem(a, L, bcs=[bc], petsc_options={
            "ksp_type": "preonly", "pc_type": "cholesky",
            "pc_factor_mat_solver_type": "mumps"},
            petsc_options_prefix=f"p{p}_").solve()
        if isinstance(uh, tuple):
            uh = uh[0]
        dxq = ufl.dx(metadata={"quadrature_degree": 2 * p + 8})
        acc = [fem.assemble_scalar(fem.form(f * dxq)) for f in (
            ufl.inner(uh - ucf, uh - ucf),
            ufl.inner(ufl.grad(uh) - Gcf, ufl.grad(uh) - Gcf),
            ufl.inner(ucf, ucf), ufl.inner(Gcf, Gcf))]
        eL2, eH1, nL2, nH1 = np.sqrt(np.real(acc))
        filas.append({"grado": p,
                      "n_gdl": int(V.dofmap.index_map.size_global * 3),
                      "solver": "Cholesky MUMPS (PETSc)",
                      "contorno": "interpolacion GLL de DOLFINx",
                      "tiempo_s": time.perf_counter() - t0,
                      "e_L2": float(eL2), "e_H1": float(eH1),
                      "e_L2_rel": float(eL2 / nL2),
                      "e_H1_rel": float(eH1 / nH1)})
    return filas


MOTORES_P = {"ngsolve": p_ngsolve, "fenicsx": p_fenicsx}


# ---------------------------------------------------------------------------
# Campana
# ---------------------------------------------------------------------------

def _version(motor):
    if motor == "app":
        sys.path.insert(0, str(AQUI.parent))
        from spinpy import __version__
        return __version__
    dist = {"ngsolve": "ngsolve", "fenicsx": "fenics-dolfinx",
            "skfem": "scikit-fem", "sfepy": "sfepy"}[motor]
    from importlib.metadata import version
    return version(dist)


def _hijo(entrada, salida):
    """Una corrida en ESTE interprete (el de FEniCSx); entrada y salida .npz.
    """
    d = dict(np.load(entrada, allow_pickle=True))
    tarea = str(d["tarea"])
    t0 = time.perf_counter()
    if tarea == "h":
        U, meta = MOTORES_H[str(d["motor"])](d["nodos"], d["elems"],
                                             d["uD"])
        meta["tiempo_s"] = time.perf_counter() - t0
        meta["version"] = _version(str(d["motor"]))
        np.savez(salida, U=U, meta=json.dumps(meta))
    else:
        filas = MOTORES_P[str(d["motor"])](d["nodos"], d["elems"],
                                           [int(g) for g in d["grados"]])
        np.savez(salida, meta=json.dumps({"filas": filas,
                                          "version": _version(
                                              str(d["motor"]))}))


def _en_fenicsx(py, **datos):
    tmp = tempfile.mkdtemp(prefix="conv_")
    ent, sal = Path(tmp) / "ent.npz", Path(tmp) / "sal.npz"
    np.savez(ent, **datos)
    p = subprocess.run([py, str(Path(__file__).resolve()), "--hijo",
                        str(ent), str(sal)], capture_output=True, text=True,
                       env=dict(os.environ, OMP_NUM_THREADS="4"))
    if p.returncode != 0 or not sal.exists():
        raise RuntimeError(p.stderr[-2000:])
    d = dict(np.load(sal, allow_pickle=True))
    meta = json.loads(str(d["meta"]))
    return (d.get("U"), meta)


def _guardar(fila):
    RES.mkdir(exist_ok=True)
    with open(SALIDA, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(fila) + "\n")
    print(json.dumps({k: fila[k] for k in fila if k in (
        "estudio", "motor", "elemento", "n", "grado", "n_gdl", "e_L2_rel",
        "e_H1_rel", "tiempo_s", "du_rel_ref")}), flush=True)


def estudio_h(motores, py_fx, ns=N_H, n_max_python_p2=16,
              tipos=("hex8", "tet4", "tet10")):
    """Refinamiento h. scikit-fem y SfePy montan en Python y su iterativo
    tarda de 15 a 17 minutos con TET10 a 8e5 GDL (Tabla 2): con TET10 se
    corren hasta n = `n_max_python_p2`."""
    for tipo in tipos:
        for n in ns:
            nodos, elems = malla_cubo(n, tipo)
            uD = u_exacta(nodos)
            ref = None
            # La referencia es el primer motor con resolvedor directo en
            # todos los tamanos: NGSolve (Cholesky) si se corre.
            for mot in sorted(motores, key=lambda m: m != "ngsolve"):
                if mot == "app" and tipo != "hex8":
                    continue
                if (mot in ("skfem", "sfepy") and tipo == "tet10"
                        and n > n_max_python_p2):
                    continue
                t0 = time.perf_counter()
                try:
                    if mot == "fenicsx":
                        U, meta = _en_fenicsx(py_fx, tarea="h", motor=mot,
                                              nodos=nodos, elems=elems,
                                              uD=uD)
                    else:
                        U, meta = MOTORES_H[mot](nodos, elems, uD)
                        meta["version"] = _version(mot)
                    meta.setdefault("tiempo_s", time.perf_counter() - t0)
                except Exception as e:                   # noqa: BLE001
                    _guardar({"estudio": "h", "motor": mot,
                              "elemento": tipo, "n": n, "h": 1.0 / n,
                              "ok": False, "motivo": repr(e)[:300]})
                    continue
                fila = {"estudio": "h", "motor": mot, "elemento": tipo,
                        "n": n, "h": 1.0 / n, "ok": True, **meta,
                        **errores(nodos, elems, U, tipo)}
                fila["du_contorno"] = float(
                    np.abs(U[en_contorno(nodos)]
                           - uD[en_contorno(nodos)]).max()
                    / np.abs(uD).max())
                if ref is None:
                    ref = U
                fila["du_rel_ref"] = float(np.abs(U - ref).max()
                                           / np.abs(ref).max())
                _guardar(fila)


def estudio_p(motores, py_fx, n=N_P, grados=GRADOS_P):
    nodos, elems = malla_cubo(n, "tet4")
    for mot in motores:
        if mot not in MOTORES_P:
            continue
        try:
            if mot == "fenicsx":
                _, meta = _en_fenicsx(py_fx, tarea="p", motor=mot,
                                      nodos=nodos, elems=elems,
                                      grados=np.array(grados))
                filas, ver = meta["filas"], meta["version"]
            else:
                filas, ver = MOTORES_P[mot](nodos, elems, grados), \
                    _version(mot)
        except Exception as e:                           # noqa: BLE001
            _guardar({"estudio": "p", "motor": mot, "n": n, "ok": False,
                      "motivo": repr(e)[:300]})
            continue
        for f in filas:
            _guardar({"estudio": "p", "motor": mot, "elemento": "tet",
                      "n": n, "h": 1.0 / n, "ok": True, "version": ver,
                      **f})


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "--hijo":
        _hijo(sys.argv[2], sys.argv[3])
        return
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("estudios", nargs="+", choices=("h", "p"))
    ap.add_argument("--motores", nargs="+",
                    default=["app", "ngsolve", "fenicsx", "skfem", "sfepy"])
    ap.add_argument("--python-fenicsx",
                    default=os.environ.get("FENICSX_PYTHON"))
    ap.add_argument("--elementos", nargs="+",
                    default=["hex8", "tet4", "tet10"])
    ap.add_argument("--n", nargs="+", type=int, default=list(N_H))
    ap.add_argument("--nuevo", action="store_true",
                    help="borra resultados/convergencia.jsonl antes")
    a = ap.parse_args()
    motores = [m for m in a.motores
               if m != "fenicsx" or a.python_fenicsx]
    r_nav, r_G = comprobar_exacta()
    print(f"solucion exacta: residuo de Navier {r_nav:.1e}, "
          f"gradiente {r_G:.1e} (diferencias centradas)")
    if a.nuevo and SALIDA.exists():
        SALIDA.unlink()
    if "h" in a.estudios:
        estudio_h(motores, a.python_fenicsx, ns=a.n, tipos=a.elementos)
    if "p" in a.estudios:
        estudio_p(motores, a.python_fenicsx)


if __name__ == "__main__":
    main()
