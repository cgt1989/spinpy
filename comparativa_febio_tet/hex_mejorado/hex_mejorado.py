"""hex_mejorado.py: Prototipo: ladrillos (hex8) con superficie ajustada y
modos incompatibles, como alternativa barata a TET10.

Misma topologia y mismos grados de libertad que `solido.malla_hex`. Tres
piezas independientes:

  suavizado   Taubin lambda/mu sobre el grafo de aristas de las caras libres
              (`malla_hex_suave(..., iters=k)`). Resultado: no mejora la
              rigidez con puntales de ~2 voxeles (ver INFORME.md).
  proyeccion  los nodos de la superficie libre se llevan por Newton sobre la
              isosuperficie f = T de un campo continuo
              (`malla_hex_suave(..., campo=(f, grad, T))`), con el desplazamiento
              acotado a `tope` voxeles por componente y los nodos de las caras
              del cubo retenidos en su plano.
  incompatible  hexaedro de Wilson-Taylor (9 modos internos condensados,
              `ke_lote_I`): corrige la rigidez excesiva a flexion del hex8 sin
              anadir grados de libertad. En voxeles regulares sigue siendo UNA
              matriz 24x24 para todos los elementos.

En ambos casos, los elementos con jacobiano pobre (< `jmin` del regular en
algun punto de Gauss) recuperan su posicion original a mitades.

`ensayo` es el de la app (1 MPa, apoyo deslizante, presion uniforme en el
techo, mm-N-MPa) con rigidez por elemento solo para los hexaedros deformados.
E_app usa el desplazamiento equivalente por trabajo (media ponderada por area
del techo), la definicion de FEBio en `febio.py`; von Mises se evalua en el
centro del elemento, como en `resistencia`.
"""
from __future__ import annotations

import time
import numpy as np
from scipy import sparse
import pyamg

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from spinpy.solido import malla_hex, ESQUINAS
from spinpy.resistencia import (_matriz_D, _modos_rigidos, _solo_portante,
                                capa_superficie, percentil_ponderado)
from spinpy.elastic import hex8_ke

XN = np.array([-1, 1, 1, -1, -1, 1, 1, -1], float)
YN = np.array([-1, -1, 1, 1, -1, -1, 1, 1], float)
ZN = np.array([-1, -1, -1, -1, 1, 1, 1, 1], float)
G = 1 / np.sqrt(3)
PG = np.array([[a, b, c] for a in (-G, G) for b in (-G, G) for c in (-G, G)])

# caras del hexaedro (orden C3D8), salientes
CARAS = np.array([[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4],
                  [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]])
# vecino por cara: (eje, signo)
VEC = [(2, -1), (2, 1), (1, -1), (0, 1), (1, 1), (0, -1)]


def dN(xi, eta, zet):
    return 0.125 * np.vstack([XN * (1 + eta * YN) * (1 + zet * ZN),
                              YN * (1 + xi * XN) * (1 + zet * ZN),
                              ZN * (1 + xi * XN) * (1 + eta * YN)])   # 3x8


DN_G = np.stack([dN(*p) for p in PG])            # (8, 3, 8)
DN_C = dN(0, 0, 0)


def jacobianos(X):
    """det J en los 8 puntos de Gauss y en el centro, normalizado por h^3/8."""
    J = np.einsum("gik,ekj->egij", DN_G, X)
    d = np.linalg.det(J)
    Jc = np.einsum("ik,ekj->eij", DN_C, X)
    return d, np.linalg.det(Jc)


def superficie_libre(BW):
    """Caras libres (quad, ids de nodo de la rejilla) y nodos de superficie."""
    nx, ny, nz = BW.shape
    nnx, nny = nx + 1, ny + 1
    P = np.pad(BW, 1, constant_values=True)   # fuera del cubo = material
    i, j, k = np.nonzero(BW)
    quads = []
    for f, (e, s) in enumerate(VEC):
        idx = [i + 1, j + 1, k + 1]
        idx[e] = idx[e] + s
        libre = ~P[tuple(idx)]
        ii, jj, kk = i[libre], j[libre], k[libre]
        q = np.stack([(ii + ESQUINAS[c, 0]) + nnx * (jj + ESQUINAS[c, 1])
                      + nnx * nny * (kk + ESQUINAS[c, 2]) for c in CARAS[f]],
                     axis=1)
        quads.append(q)
    return np.concatenate(quads)


def malla_hex_suave(BW, sp, iters=20, lam=0.5, mu=-0.53, tope=0.5,
                    jmin=0.2, verbose=False, campo=None, newton=10):
    """campo = (f, grad, T): proyecta los nodos de superficie sobre f = T
    (Newton a lo largo del gradiente) en lugar de suavizar con Taubin."""
    BW = np.asarray(BW, bool)
    sp = np.repeat(float(sp), 3) if np.ndim(sp) == 0 else np.asarray(sp, float)
    nodos, elems, inf = malla_hex(BW, sp)
    nx, ny, nz = BW.shape
    nnx, nny = nx + 1, ny + 1
    # ids globales de rejilla de los nodos compactos
    g = (np.round(nodos / sp).astype(np.int64))
    gid = g[:, 0] + nnx * g[:, 1] + nnx * nny * g[:, 2]
    orden = np.argsort(gid)
    q = superficie_libre(BW)
    q = orden[np.searchsorted(gid[orden], q)]           # a ids compactos
    # aristas del grafo de superficie
    ar = np.concatenate([q[:, [0, 1]], q[:, [1, 2]], q[:, [2, 3]], q[:, [3, 0]]])
    ar = np.unique(np.sort(ar, axis=1), axis=0)
    n = nodos.shape[0]
    A = sparse.coo_matrix((np.ones(2 * len(ar)),
                           (np.r_[ar[:, 0], ar[:, 1]], np.r_[ar[:, 1], ar[:, 0]])),
                          shape=(n, n)).tocsr()
    grado = np.asarray(A.sum(1)).ravel()
    sup = grado > 0
    W = sparse.diags(np.where(sup, 1 / np.maximum(grado, 1), 0)) @ A
    # restricciones de plano: coordenada fija si el nodo esta en una cara
    lo = np.zeros(3); hi = np.array([nx, ny, nz]) * sp
    fijo = np.zeros((n, 3), bool)
    for e in range(3):
        fijo[:, e] = (np.abs(nodos[:, e] - lo[e]) < 1e-9 * sp[e]) | \
                     (np.abs(nodos[:, e] - hi[e]) < 1e-9 * sp[e])
    X0 = nodos.copy()
    X = nodos.copy()
    movil = sup[:, None] & ~fijo
    if campo is not None:
        fz, gz, T = campo
        idx = np.nonzero(sup)[0]
        Y = X0[idx].copy()
        fj = fijo[idx]
        for it in range(newton):
            v = fz(Y) - T
            gr = gz(Y)
            gr[fj] = 0.0
            paso = (v / np.maximum((gr ** 2).sum(1), 1e-30))[:, None] * gr
            Y = Y - paso
            Y = X0[idx] + np.clip(Y - X0[idx], -tope * sp, tope * sp)
        X[idx] = Y
        iters = 0
    for it in range(iters):
        for f in (lam, mu):
            L = W @ X - X
            X = X + f * L * movil
            X = X0 + np.clip(X - X0, -tope * sp, tope * sp)
    # control de calidad: retroceso por elemento
    vol0 = np.prod(sp)
    for r in range(30):
        d, dc = jacobianos(X[elems])
        malo = (d.min(axis=1) < jmin * vol0 / 8)
        if not malo.any():
            break
        nm = np.unique(elems[malo])
        X[nm] = X0[nm] + 0.5 * (X[nm] - X0[nm])
    else:
        raise RuntimeError("no se pudo recuperar la calidad")
    d, dc = jacobianos(X[elems])
    vol = (d * 1.0).sum(axis=1)                        # pesos Gauss = 1
    deform = np.any(np.abs(X[elems] - X0[elems]) > 1e-12, axis=(1, 2))
    info = dict(inf, tipo="hex8s", iters=iters, tope=tope, retrocesos=r,
                n_sup=int(sup.sum()), n_deform=int(deform.sum()),
                V=float(vol.sum()), V_vox=float(elems.shape[0] * vol0),
                jmin=float(d.min() / (vol0 / 8)))
    info["dV_pct"] = 100 * (info["V"] / info["V_vox"] - 1)
    return X, elems, vol, deform, info


def ke_lote(X, D):
    """Rigidez 24x24 de cada hexaedro de X (e, 8, 3), Gauss 2x2x2."""
    ne = X.shape[0]
    Ke = np.zeros((ne, 24, 24))
    for gp in range(8):
        J = np.einsum("ik,ekj->eij", DN_G[gp], X)
        det = np.linalg.det(J)
        dNdx = np.linalg.solve(J, np.broadcast_to(DN_G[gp], (ne, 3, 8)))
        B = matriz_B(dNdx)
        Ke += np.einsum("eai,ab,ebj->eij", B, D, B) * det[:, None, None]
    return 0.5 * (Ke + Ke.transpose(0, 2, 1))


def matriz_B(dNdx):
    ne = dNdx.shape[0]
    B = np.zeros((ne, 6, 24))
    B[:, 0, 0::3] = dNdx[:, 0]; B[:, 1, 1::3] = dNdx[:, 1]; B[:, 2, 2::3] = dNdx[:, 2]
    B[:, 3, 1::3] = dNdx[:, 2]; B[:, 3, 2::3] = dNdx[:, 1]
    B[:, 4, 0::3] = dNdx[:, 2]; B[:, 4, 2::3] = dNdx[:, 0]
    B[:, 5, 0::3] = dNdx[:, 1]; B[:, 5, 1::3] = dNdx[:, 0]
    return B


def ensayo(X, elems, vol, deform, sp, forma, E=20e3, nu=0.3, sigma0=1.0,
           tol=1e-10, incompatible=False):
    """Ensayo de la app en mm-N-MPa: presion sigma0 en el techo, deslizante."""
    sp = np.asarray(sp, float)
    D = _matriz_D(E, nu)
    n = X.shape[0]; ndof = 3 * n; ne = elems.shape[0]
    edof = (3 * elems[:, :, None] + np.arange(3)).reshape(ne, 24)
    fke = ke_lote_I if incompatible else ke_lote
    X0 = np.array(ESQUINAS, float)[None] * sp
    ke0 = fke(X0, D)[0]
    vals = np.broadcast_to(ke0.ravel(), (ne, 576)).copy()
    idx = np.nonzero(deform)[0]
    for a in range(0, idx.size, 20000):
        s = idx[a:a + 20000]
        vals[s] = fke(X[elems[s]], D).reshape(-1, 576)
    K = sparse.coo_matrix((vals.ravel(), (np.repeat(edof, 24, 1).ravel(),
                                           np.tile(edof, (1, 24)).ravel())),
                          shape=(ndof, ndof)).tocsr()
    H = forma[2] * sp[2]
    A = forma[0] * sp[0] * forma[1] * sp[1]
    z = X[:, 2]
    tz = 1e-9 * H
    # techo: caras superiores (nodos 4..7) en el plano z = H; fuerza
    # consistente con 2x2 Gauss en la cara (cuadrilatero plano en z)
    et = np.nonzero(np.all(np.abs(z[elems[:, 4:]] - H) < tz, axis=1))[0]
    F = np.zeros(ndof)
    Q = X[elems[et][:, 4:8], :2]                       # (m, 4, 2)
    xs = np.array([-1, 1, 1, -1.]); ys = np.array([-1, -1, 1, 1.])
    for a_ in (-G, G):
        for b_ in (-G, G):
            Nq = 0.25 * (1 + a_ * xs) * (1 + b_ * ys)
            dxi = 0.25 * np.vstack([xs * (1 + b_ * ys), ys * (1 + a_ * xs)])
            J = np.einsum("ik,ekj->eij", dxi, Q)
            det = np.linalg.det(J)
            for c in range(4):
                np.add.at(F, 3 * elems[et, 4 + c] + 2, -sigma0 * Nq[c] * det)
    F_total = -F.sum()
    base = np.nonzero(np.abs(z) < tz)[0]
    cb = X[base]
    a = base[np.argmin(cb[:, 0] + cb[:, 1])]
    b = base[np.argmax(cb[:, 0] - cb[:, 1])]
    fijos = np.unique(np.r_[3 * base + 2, 3 * a, 3 * a + 1, 3 * b + 1])
    libres = np.setdiff1d(np.arange(ndof), fijos)
    Kff = K[libres][:, libres].tocsr()
    B = _modos_rigidos(X)[libres]
    t0 = time.perf_counter()
    ml = pyamg.smoothed_aggregation_solver(Kff, B=B, max_coarse=500)
    u = np.zeros(ndof)
    res = []
    u[libres] = ml.solve(F[libres], tol=tol, maxiter=1000, accel="cg",
                         residuals=res)
    t_sol = time.perf_counter() - t0
    rel = np.linalg.norm(F[libres] - Kff @ u[libres]) / np.linalg.norm(F[libres])
    # desplazamiento equivalente por trabajo = media ponderada por area
    uz = -(F @ u) / F_total
    E_app = abs((F_total / A) / (uz / H))
    # tension en el centro de cada elemento
    Jc = np.einsum("ik,ekj->eij", DN_C, X[elems])
    dNdx = np.linalg.solve(Jc, np.broadcast_to(DN_C, (ne, 3, 8)))
    eps = np.einsum("eai,ei->ea", matriz_B(dNdx), u[edof])
    s = eps @ D.T
    vm = np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2
                        + (s[:, 2] - s[:, 0]) ** 2)
                 + 3 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2))
    return dict(E_app=E_app, vm=vm / sigma0, res=rel, t_sol=t_sol,
                n_gdl=int(libres.size), it=len(res), F=F_total)


def correr(BW, sp, suave=True, incompatible=False, **kw):
    BW = _solo_portante(np.asarray(BW, bool))
    sp3 = np.repeat(float(sp), 3)
    t0 = time.perf_counter()
    if suave:
        X, el, vol, dfm, inf = malla_hex_suave(BW, sp3, **kw)
    else:
        X, el, inf = malla_hex(BW, sp3)
        vol = np.full(el.shape[0], np.prod(sp3)); dfm = np.zeros(el.shape[0], bool)
        inf = dict(inf, V=float(vol.sum()), dV_pct=0.0, jmin=1.0, n_deform=0)
    t_m = time.perf_counter() - t0
    r = ensayo(X, el, vol, dfm, sp3, BW.shape, incompatible=incompatible)
    capa = capa_superficie(BW)[np.nonzero(BW)]
    vs, ws = r["vm"][capa], vol[capa]
    return dict(E_app=r["E_app"], vm_max=float(vs.max()) if vs.size else np.nan,
                vm_p99=float(percentil_ponderado(vs, ws, 99, "hazen")) if vs.size else np.nan,
                dV_pct=inf["dV_pct"], jmin=inf["jmin"], n_def=inf["n_deform"],
                n_el=int(el.shape[0]), gdl=r["n_gdl"], res=r["res"],
                t_malla=t_m, t_sol=r["t_sol"], it=r["it"])


# ---------------------------------------------------------------------------
# Hexaedro con modos incompatibles (Wilson-Taylor, QM6 en 3D): 9 modos
# internos (1 - xi^2), (1 - eta^2), (1 - zeta^2) por componente, condensados.
# Derivadas de los modos con J0 del centro y factor detJ0/detJ (Taylor 1976)
# para pasar la prueba de la parcela en elementos distorsionados.
# ---------------------------------------------------------------------------

def _dM(xi, eta, zet):
    return np.diag([-2 * xi, -2 * eta, -2 * zet])        # d(1-s^2)/ds, 3x3


def matriz_Ba(dMdx):
    """B de los 9 modos internos: dMdx (e, 3 derivadas, 3 modos)."""
    ne = dMdx.shape[0]
    B = np.zeros((ne, 6, 9))
    B[:, 0, 0::3] = dMdx[:, 0]; B[:, 1, 1::3] = dMdx[:, 1]; B[:, 2, 2::3] = dMdx[:, 2]
    B[:, 3, 1::3] = dMdx[:, 2]; B[:, 3, 2::3] = dMdx[:, 1]
    B[:, 4, 0::3] = dMdx[:, 2]; B[:, 4, 2::3] = dMdx[:, 0]
    B[:, 5, 0::3] = dMdx[:, 1]; B[:, 5, 1::3] = dMdx[:, 0]
    return B


def ke_lote_I(X, D):
    ne = X.shape[0]
    J0 = np.einsum("ik,ekj->eij", DN_C, X)
    d0 = np.linalg.det(J0)
    Kuu = np.zeros((ne, 24, 24)); Kua = np.zeros((ne, 24, 9)); Kaa = np.zeros((ne, 9, 9))
    for gp in range(8):
        J = np.einsum("ik,ekj->eij", DN_G[gp], X)
        det = np.linalg.det(J)
        Bu = matriz_B(np.linalg.solve(J, np.broadcast_to(DN_G[gp], (ne, 3, 8))))
        dM = np.broadcast_to(_dM(*PG[gp]), (ne, 3, 3))
        Ba = matriz_Ba(np.linalg.solve(J0, dM)) * (d0 / det)[:, None, None]
        DBu = np.einsum("ab,ebj->eaj", D, Bu); DBa = np.einsum("ab,ebj->eaj", D, Ba)
        w = det[:, None, None]
        Kuu += np.einsum("eai,eaj->eij", Bu, DBu) * w
        Kua += np.einsum("eai,eaj->eij", Bu, DBa) * w
        Kaa += np.einsum("eai,eaj->eij", Ba, DBa) * w
    K = Kuu - Kua @ np.linalg.solve(Kaa, Kua.transpose(0, 2, 1))
    return 0.5 * (K + K.transpose(0, 2, 1))
