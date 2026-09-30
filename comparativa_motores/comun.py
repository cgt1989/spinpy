"""
comun.py: Casos de prueba y postproceso comun de la comparativa de motores.

Un CASO es un problema discreto completo (malla, material, apoyo, carga)
escrito en un .npz que leen los cinco adaptadores de `motores/`. El
POSTPROCESO se hace aqui, con las funciones de la app, a partir del
desplazamiento nodal que devuelve cada motor: asi E_app, von Mises, el p99 de
la capa superficial y el criterio de Pistoia se calculan con la MISMA
formula para los cinco, y lo que se compara es la solucion del sistema.

Unidades de los casos: mm, N, MPa (E en MPa, tension en MPa, u en mm). El
postproceso pasa a Pa solo para llamar a `resistencia.criterio_pistoia`, que
trabaja en las unidades de E_s y es indiferente a ellas salvo por el
cociente adimensional eps_eff.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))

from spinpy import febio                                       # noqa: E402
from spinpy.resistencia import (_solo_portante, capa_superficie,  # noqa: E402
                                criterio_pistoia)
from spinpy.solido import malla_hex                            # noqa: E402

E_MPA, NU = 20000.0, 0.30
SIGMA_REF_MPA = 1.0          # 1 MPa sobre la seccion bruta: el ensayo de la app


# ---------------------------------------------------------------------------
# Construccion de casos
# ---------------------------------------------------------------------------

def _anclas(nodos, base):
    """Las dos anclas del apoyo deslizante, con la regla de la app."""
    cb = nodos[base]
    a = int(base[int(np.argmin(cb[:, 0] + cb[:, 1]))])
    b = int(base[int(np.argmax(cb[:, 0] - cb[:, 1]))])
    return a, b


def _caras_z_hex(nodos, elems):
    z = nodos[:, 2]
    tol = 1e-9 * max(float(np.ptp(z)), 1.0)
    techo = elems[np.all(z[elems[:, 4:8]] >= z.max() - tol, axis=1), 4:8]
    base = elems[np.all(z[elems[:, 0:4]] <= z.min() + tol, axis=1), 0:4]
    return base, techo


def caso_hex(BW, spacing, nombre, apoyo="deslizante", **meta):
    """Ensayo de compresion en z sobre la malla de voxeles (hex8) de la app.

    La malla es la de `resistencia.ensayo_compresion`: solo el hueso portante
    (`_solo_portante`), un hexaedro por voxel, nodos en el orden de
    `np.unique` sobre la rejilla. `malla_hex` produce exactamente ese orden,
    de modo que el `u` de la app y el de los motores son comparables GDL a GDL.
    """
    BW = np.asarray(BW, bool)
    spacing = np.broadcast_to(np.asarray(spacing, float), (3,)).copy()
    B = _solo_portante(BW)
    nodos, elems, _ = malla_hex(B, spacing)
    caras_base, caras_techo = _caras_z_hex(nodos, elems)
    z = nodos[:, 2]
    tol = 1e-9 * max(float(np.ptp(z)), 1.0)
    base = np.nonzero(z <= z.min() + tol)[0]
    techo = np.nonzero(z >= z.max() - tol)[0]
    a, b = _anclas(nodos, base)
    A_bruta = float(BW.shape[0] * spacing[0] * BW.shape[1] * spacing[1])
    A_osea = float(caras_techo.shape[0] * spacing[0] * spacing[1])
    sup = capa_superficie(B)[np.nonzero(B)]
    m = {"nombre": nombre, "tipo": "hex8", "apoyo": apoyo, "E": E_MPA,
         "nu": NU, "sigma_ref": SIGMA_REF_MPA, "A_bruta": A_bruta,
         "A_osea_techo": A_osea, "p": SIGMA_REF_MPA * A_bruta / A_osea,
         "H": float(np.ptp(z)), "ancla_xy": a, "ancla_y": b,
         "forma": list(BW.shape), "spacing": spacing.tolist(),
         "BVTV": float(BW.mean()), "frac_portante": float(B.sum() / BW.sum()),
         "n_nodos": int(nodos.shape[0]), "n_elems": int(elems.shape[0]),
         "n_gdl": int(3 * nodos.shape[0]), "analisis": "lineal"}
    m.update(meta)
    return {"nodos": nodos, "elems": elems, "caras_techo": caras_techo,
            "caras_base": caras_base, "base_nodos": base,
            "techo_nodos": techo, "superficie": sup,
            "vol_elem": np.full(elems.shape[0], float(np.prod(spacing))),
            "BW": B, "meta": m}


def caso_tet(malla, nombre, apoyo="deslizante", **meta):
    """Ensayo de compresion sobre una malla TET10 de `febio.mallar`.

    Se reutiliza el post-proceso de la integracion con FEBio
    (`febio.preparar_tet10`): base y techo en sus planos exactos, caras tri6
    del techo con la normal hacia fuera, capa superficial sin las caras del
    cubo. Los motores reciben las ESQUINAS y los nodos intermedios; los que
    construyen su propio espacio P2 sobre las esquinas (todos menos el
    adaptador propio) generan los mismos nodos, porque las aristas son rectas
    y el nodo intermedio esta en su punto medio (`solido.verificar_tet10`).
    """
    nodos, elems = np.array(malla["nodos"], float), malla["elems"]
    # ARISTAS RECTAS. `febio.preparar_tet10` devuelve a su plano exacto los
    # nodos a menos de 1e-3 h de un plano del cubo; si una esquina se mueve y
    # la otra no, el nodo intermedio deja de estar en el punto medio (medido:
    # hasta 7e-5 del lado del cubo). Un motor que construye P2 sobre las
    # esquinas usaria entonces otra geometria que el que recibe los diez
    # nodos. Se recoloca cada intermedio en su punto medio (un nodo compartido
    # recibe el mismo valor desde todos sus elementos) y se declara cuanto se
    # movio.
    mov = 0.0
    for k, (i, j) in enumerate(((0, 1), (1, 2), (0, 2), (0, 3), (1, 3),
                                (2, 3))):
        med = 0.5 * (nodos[elems[:, i]] + nodos[elems[:, j]])
        mov = max(mov, float(np.abs(nodos[elems[:, 4 + k]] - med).max()))
        nodos[elems[:, 4 + k]] = med
    z = nodos[:, 2]
    z0, z1 = malla["z0"], malla["z1"]
    tol = 1e-6 * float(malla["spacing"][2])
    base = np.nonzero(np.abs(z - z0) <= tol)[0]
    techo = np.nonzero(np.abs(z - z1) <= tol)[0]
    planos = febio._caras_planos(malla)
    caras_base = planos[(2, 0)][0]
    caras_techo = malla["caras_techo"]
    # Anclas entre las ESQUINAS de la base: los motores que construyen su
    # espacio P2 solo reciben esquinas, y un GDL puntual en un nodo intermedio
    # no se puede nombrar igual en todos.
    esquinas = np.unique(elems[:, :4])
    a, b = _anclas(nodos, np.intersect1d(base, esquinas))
    from spinpy.escribe import area_caras
    A_osea = float(area_caras(nodos, caras_techo).sum())
    m = {"nombre": nombre, "tipo": "tet10", "apoyo": apoyo, "E": E_MPA,
         "nu": NU, "sigma_ref": SIGMA_REF_MPA, "A_bruta": malla["A_bruta"],
         "A_osea_techo": A_osea,
         "p": SIGMA_REF_MPA * malla["A_bruta"] / A_osea,
         "H": float(z1 - z0), "ancla_xy": a, "ancla_y": b,
         "forma": list(malla["forma"]),
         "spacing": np.asarray(malla["spacing"]).tolist(),
         "BVTV": float(malla["vol_elem"].sum()
                       / np.prod(np.asarray(malla["forma"])
                                 * malla["spacing"])),
         "n_nodos": int(nodos.shape[0]), "n_elems": int(elems.shape[0]),
         "n_gdl": int(3 * nodos.shape[0]), "analisis": "lineal",
         "rectificacion_max_h": mov / float(malla["spacing"][0])}
    m.update(meta)
    return {"nodos": nodos, "elems": elems, "caras_techo": caras_techo,
            "caras_base": caras_base, "base_nodos": base,
            "techo_nodos": techo, "superficie": malla["superficie"],
            "vol_elem": malla["vol_elem"], "meta": m}


def guardar(caso, ruta):
    ruta = Path(ruta)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    arr = {k: v for k, v in caso.items() if k != "meta"}
    np.savez(ruta, meta=json.dumps(caso["meta"]), **arr)
    return ruta


# ---------------------------------------------------------------------------
# Postproceso comun
# ---------------------------------------------------------------------------

def matriz_D(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[0, 0] = D[1, 1] = D[2, 2] = lam + 2 * mu
    D[3, 3] = D[4, 4] = D[5, 5] = mu
    return D


def _grad_hex_centro(nodos, elems):
    """Gradiente de las 8 funciones de forma en el centro de cada voxel."""
    P = nodos[elems]                               # (M, 8, 3)
    h = P[:, 6] - P[:, 0]                          # diagonal: (dx, dy, dz)
    s = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                  [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
    return 0.125 * s[None, :, :] * (2.0 / h)[:, None, :]    # (M, 8, 3)


def _grad_tet10_centroide(nodos, elems):
    """Gradiente de las 10 funciones de forma P2 en el centroide (orden C3D10).

    Con aristas rectas el campo de deformacion es LINEAL en el elemento, asi
    que su valor en el centroide es su media de volumen: la misma magnitud que
    FEBio escribe (media de los puntos de Gauss) y la que usa `febio.leer`.
    """
    L = np.full(4, 0.25)
    # dN/dL_i para N = L_i (2 L_i - 1) y N_ij = 4 L_i L_j
    dNdL = np.zeros((10, 4))
    for i in range(4):
        dNdL[i, i] = 4 * L[i] - 1
    for k, (i, j) in enumerate(((0, 1), (1, 2), (0, 2), (0, 3), (1, 3),
                                (2, 3))):
        dNdL[4 + k, i] = 4 * L[j]
        dNdL[4 + k, j] = 4 * L[i]
    # L0 = 1 - xi - eta - zeta
    dLdxi = np.array([[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    dNdxi = dNdL @ dLdxi                              # (10, 3)
    P = nodos[elems[:, :4]]
    J = np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]],
                 axis=2)                              # (M, 3, 3): dx/dxi
    Jinv = np.linalg.inv(J)
    return np.einsum("ak,mkj->maj", dNdxi, Jinv)      # (M, 10, 3)


def deformacion_elemental(nodos, elems, u):
    """Deformacion (Voigt, distorsiones de ingenieria) por elemento."""
    G = (_grad_hex_centro(nodos, elems) if elems.shape[1] == 8
         else _grad_tet10_centroide(nodos, elems))
    Ue = u[elems]                                     # (M, n, 3)
    H = np.einsum("mai,maj->mij", Ue, G)              # du_i/dx_j
    return np.stack([H[:, 0, 0], H[:, 1, 1], H[:, 2, 2],
                     H[:, 1, 2] + H[:, 2, 1], H[:, 0, 2] + H[:, 2, 0],
                     H[:, 0, 1] + H[:, 1, 0]], axis=1)


def von_mises(s):
    return np.sqrt(np.maximum(
        0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2
               + (s[:, 2] - s[:, 0]) ** 2)
        + 3.0 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2), 0.0))


def postproceso(caso, u):
    """Magnitudes del ensayo lineal a partir de `u` (N, 3) en mm."""
    m = caso["meta"]
    nodos, elems = caso["nodos"], caso["elems"]
    eps = deformacion_elemental(nodos, elems, u)
    sig = eps @ matriz_D(m["E"], m["nu"]).T           # MPa
    U = 0.5 * np.einsum("ij,ij->i", sig, eps)
    eps_eff = np.sqrt(np.maximum(2.0 * U / m["E"], 0.0))
    vm = von_mises(sig)
    vol = caso["vol_elem"]
    H = m["H"]
    techo = caso["techo_nodos"]
    # E_app, dos definiciones: media simple de los nodos del techo (la de la
    # app con voxeles) y media ponderada por area de las caras cargadas.
    uz_nodal = float(u[techo, 2].mean())
    from spinpy.escribe import area_caras
    caras = caso["caras_techo"]
    A = area_caras(nodos, caras)
    uzc = (u[caras, 2].mean(axis=1) if caras.shape[1] == 4
           else u[caras[:, 3:], 2].mean(axis=1))
    uz_area = float((A * uzc).sum() / A.sum())
    s_ref = m["sigma_ref"]
    F_ref = s_ref * m["A_bruta"]
    F_vol = float(-(sig[:, 2] * vol).sum() / H)       # N (MPa * mm^2)
    res = {"ok": True, "sigma_app": s_ref * 1e6, "F_total": F_ref * 1e6,
           "eps_eff_solido": eps_eff, "vm_solido": vm * 1e6,
           "superficie_solido": np.asarray(caso["superficie"], bool),
           "E_app": s_ref / (abs(uz_nodal) / H) * 1e6}
    if m["tipo"] == "tet10":
        res["vol_solido"] = vol
    p = criterio_pistoia(res)
    return {"E_app_nodal_MPa": s_ref / (abs(uz_nodal) / H),
            "E_app_area_MPa": s_ref / (abs(uz_area) / H),
            "F_vol_N": F_vol, "dF_rel": abs(F_vol - F_ref) / F_ref,
            "vm_p99_sup_MPa": p.get("vm_p99_superficie", np.nan) / 1e6,
            "vm_max_MPa": float(vm.max()),
            "sigma_fallo_MPa": p["sigma_fallo"] / 1e6,
            "factor_pistoia": p["factor"],
            "u_max_mm": float(np.abs(u).max()),
            "_vm": vm, "_sig": sig}


def dif_rel(a, b):
    """max|a - b| / max|b|."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.abs(a - b).max() / max(np.abs(b).max(), 1e-300))


# ---------------------------------------------------------------------------
# Soluciones cerradas del no lineal: traccion/compresion uniaxial homogenea
# ---------------------------------------------------------------------------

def uniaxial(material, lam_z, E=E_MPA, nu=NU):
    """Tension nominal P_zz (MPa) de un bloque con F = diag(a, a, lam_z).

    Bloque con base deslizante, laterales libres y plato sin friccion: la
    deformacion es homogenea y S_xx = S_yy = 0 fija el estiramiento lateral a.
      SVK:  S = lam tr(E) I + 2 mu E  ->  E_xx = -nu E_zz  (exacto),
            P_zz = lam_z * E_young * (lam_z^2 - 1) / 2.
      neo-Hookeano de FEBio: S = mu (I - C^-1) + lam ln(J) C^-1, J = a^2 lam_z;
            S_xx = 0  <=>  mu (a^2 - 1) + lam ln(a^2 lam_z) = 0 (raiz en a),
            P_zz = lam_z * [mu (1 - lam_z^-2) + lam ln(J) lam_z^-2].
    Devuelve (P_zz, a).
    """
    from scipy.optimize import brentq
    lmb = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    if material == "svk":
        Ezz = 0.5 * (lam_z ** 2 - 1)
        a = np.sqrt(1 - 2 * nu * Ezz)
        return lam_z * E * Ezz, a
    f = lambda a: mu * (a * a - 1) + lmb * np.log(a * a * lam_z)   # noqa: E731
    a = brentq(f, 0.5, 2.0, xtol=1e-15, rtol=1e-15)
    J = a * a * lam_z
    return lam_z * (mu * (1 - lam_z ** -2) + lmb * np.log(J) / lam_z ** 2), a
