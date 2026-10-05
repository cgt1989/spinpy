"""
correcciones.py: Evaluacion de las correcciones de A1, A2 y A3 frente a una
referencia embebida. Criterio de aceptacion: PRERREGISTRO.md (escrito antes).

    cd comparativa_motores/correcciones
    python correcciones.py referencia hex8 tet10    # por etapas o todo

Escribe resultados/*.json (una entrada por ensayo y region).

REFERENCIA. El mismo campo aleatorio del espinodoide se evalua en un bloque
que rodea al VOI de 5 mm con 1,25 mm (n/4 voxeles) por cada lado; los voxeles
interiores son identicos a los del VOI (se comprueba). El bloque se ensaya con
plato rigido y base deslizante y se mide SOLO en la region del VOI, que no
tiene ni techo cargado ni caras cortadas.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))
sys.path.insert(0, str(AQUI.parent.parent))
import casos                                                    # noqa: E402
from spinpy import fem, motores                                 # noqa: E402
from spinpy.grf import level_set, wave_directions, _evaluar_campo  # noqa: E402
from spinpy.resistencia import (_solo_portante, capa_superficie,  # noqa: E402
                                percentil_ponderado)

E_MPA, NU = 20000.0, 0.30
L = casos.LADO_MM
FRANJA = 0.25          # mm, espesor de cada franja del extensometro
NUCLEO = 0.625         # mm, margen de la ROI de F2 (n/8 voxeles)
RES = AQUI / "resultados"


# ---------------------------------------------------------------------------
# Geometria
# ---------------------------------------------------------------------------

def mascara(n, margen=0):
    """Mascara del espinodoide de referencia a resolucion n, ampliada con
    `margen` voxeles por cada lado evaluando el MISMO campo fuera de [0, 1]."""
    p = casos.ESPINODOIDE
    rng = np.random.default_rng(p["seed"])
    dirs = wave_directions(p["num_waves"], p["thetas"], rng=rng)
    fases = rng.uniform(0.0, 2.0 * np.pi, p["num_waves"])
    eje = np.arange(-margen, n + margen) / (n - 1.0)
    X, Y, Z = np.meshgrid(eje, eje, eje, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    s = _evaluar_campo(pts, dirs, fases, p["wave_number"], "numpy", "f64",
                       None)
    G = (np.sqrt(2.0 / p["num_waves"]) * s).reshape((n + 2 * margen,) * 3)
    return G <= level_set(p["rho"])


def comprobar_mascara(n):
    BW, _ = casos.espinodoide(n)
    m = n // 4
    Bb = mascara(n, m)
    ok = bool(np.array_equal(Bb[m:m + n, m:m + n, m:m + n], BW))
    return ok, float(BW.mean()), float(Bb.mean())


# ---------------------------------------------------------------------------
# Ensayo y magnitudes por region
# ---------------------------------------------------------------------------

def centroides_u(nodos, elems, u):
    """Centroide y u_z en el centroide de cada elemento."""
    if elems.shape[1] == 8:
        return nodos[elems].mean(1), u[elems, 2].mean(1)
    c = nodos[elems[:, :4]].mean(1)
    uz = -0.125 * u[elems[:, :4], 2].sum(1) + 0.25 * u[elems[:, 4:], 2].sum(1)
    return c, uz


def ensayo(malla, control, desp=0.0, voxeles=None):
    """Resuelve y devuelve los campos por elemento (coordenadas en el marco
    del VOI de hex8: [0, L]^3, con `desp` mm restados)."""
    carga = (1.0,) if control == "fuerza" else (0.001,)
    p = motores.problema_de_malla(malla, E=E_MPA, nu=NU, apoyo="deslizante",
                                  control=control, cargas=carga)
    t0 = time.perf_counter()
    out = motores.resolver(p, "ngsolve")
    t = time.perf_counter() - t0
    nodos, elems, u = p["nodos"], np.asarray(malla["elems"]), out["u"]
    sig = fem.tension_elemental(nodos, elems, u, E_MPA, NU)
    c, uz = centroides_u(nodos, elems, u)
    F = float(out["F_reac"][-1])
    H = malla["H"]
    if control == "fuerza":
        E_app = (F / malla["A_bruta"]) / (abs(fem._uz_medio(malla, u)) / H)
    else:
        E_app = (F / malla["A_bruta"]) / carga[0]
    return {"sig": sig, "vm": fem.von_mises(sig), "cen": c - desp, "uz": uz,
            "vol": np.asarray(malla["vol_elem"], float),
            "sup": np.asarray(malla["superficie"], bool), "F": F,
            "E_app": E_app, "t": t, "gdl": int(3 * nodos.shape[0]),
            "solver": out["meta"].get("solver"), "vox": voxeles}


def region(r, lo, hi):
    """Magnitudes en la caja [lo, hi]^3 (mm, marco del VOI)."""
    c = r["cen"]
    dentro = np.all((c >= lo) & (c <= hi), axis=1)
    V_caja = (hi - lo) ** 3
    w = r["vol"]
    sR = float(-(r["sig"][dentro, 2] * w[dentro]).sum() / V_caja)
    z = c[:, 2]
    arr = dentro & (z >= hi - FRANJA)
    aba = dentro & (z <= lo + FRANJA)
    zt, zb = (np.average(z[m], weights=w[m]) for m in (arr, aba))
    ut, ub = (np.average(r["uz"][m], weights=w[m]) for m in (arr, aba))
    eps = abs((ut - ub) / (zt - zb))
    s = dentro & r["sup"]
    p99 = float(percentil_ponderado(r["vm"][s], w[s], 99)) / sR
    return {"sigma": sR, "eps": float(eps), "E": sR / eps, "p99": p99,
            "n_sup": int(s.sum()), "_dentro": dentro}


def campo_vs_ref(r, ref, lo, hi, n, m):
    """Comparacion voxel a voxel (hex8) en la capa superficial de la region."""
    def tabla(x, off):
        d = {}
        reg = region(x, lo, hi)["_dentro"]
        for e in np.nonzero(reg)[0]:
            d[tuple(x["vox"][e] + off)] = e
        return d, reg
    sR = region(r, lo, hi)["sigma"]
    sRef = region(ref, lo, hi)["sigma"]
    dm, _ = tabla(r, 0)
    dr, _ = tabla(ref, -m)
    comunes = [(dr[k], dm[k]) for k in dr if k in dm and ref["sup"][dr[k]]]
    ie = np.array(comunes)
    vr = ref["vm"][ie[:, 0]] / sRef
    vm = r["vm"][ie[:, 1]] / sR
    err = np.abs(vm - vr) / np.maximum(vr, 1e-12)
    k = max(1, int(0.01 * len(vr)))
    top_r = set(np.argsort(vr)[-k:])
    top_m = set(np.argsort(vm)[-k:])
    return {"err_mediana": float(np.median(err)),
            "err_p90": float(np.percentile(err, 90)),
            "cola_recuperada": len(top_r & top_m) / k, "n": int(len(vr))}


def publico(d):
    return {k: v for k, v in d.items() if not k.startswith("_")}


def guardar(nombre, datos):
    RES.mkdir(exist_ok=True)
    (RES / f"{nombre}.json").write_text(json.dumps(datos, indent=1,
                                                   default=float))


# ---------------------------------------------------------------------------
# Etapas
# ---------------------------------------------------------------------------

_REF = {}


def referencia(n, control="plato"):
    if (n, control) in _REF:
        return _REF[(n, control)]
    m = n // 4
    h = L / n
    Bb = mascara(n, m)
    t0 = time.perf_counter()
    malla = fem.mallar(Bb, np.full(3, h), "hex8")
    B = _solo_portante(Bb)
    vox = np.stack(np.nonzero(B), axis=1)
    r = ensayo(malla, control, desp=m * h, voxeles=vox)
    r["t_malla"] = time.perf_counter() - t0
    _REF[(n, control)] = r
    return r


def etapa_referencia():
    out = {}
    for n in (32, 48, 64):
        ok, bv, bvb = comprobar_mascara(n)
        r = referencia(n)
        R = region(r, 0.0, L)
        Rn = region(r, NUCLEO, L - NUCLEO)
        out[f"n{n}"] = {"interior_identico": ok, "BVTV_voi": bv,
                        "BVTV_bloque": bvb, "gdl": r["gdl"], "t": r["t"],
                        "solver": r["solver"], "voi": publico(R),
                        "nucleo": publico(Rn)}
        print("ref", n, ok, r["gdl"], f"{r['t']:.1f}s",
              f"E_voi {R['E']:.2f} p99 {R['p99']:.2f}",
              f"E_nuc {Rn['E']:.2f} p99 {Rn['p99']:.2f}", flush=True)
    # Sensibilidad: traccion en el techo del bloque
    r = referencia(48, "fuerza")
    out["n48_traccion"] = {"voi": publico(region(r, 0.0, L)),
                           "nucleo": publico(region(r, NUCLEO, L - NUCLEO))}
    print("ref traccion 48", out["n48_traccion"]["voi"]["E"], flush=True)
    guardar("referencia", out)
    return out


def etapa_hex8():
    out = {}
    for n in (32, 48, 64):
        m = n // 4
        ref = referencia(n)
        BW, sp = casos.espinodoide(n)
        malla = fem.mallar(BW, sp, "hex8")
        vox = np.stack(np.nonzero(_solo_portante(BW)), axis=1)
        for control in ("fuerza", "plato"):
            r = ensayo(malla, control, voxeles=vox)
            R = region(r, 0.0, L)
            Rn = region(r, NUCLEO, L - NUCLEO)
            out[f"n{n}_{control}"] = {
                "E_app": r["E_app"], "gdl": r["gdl"], "t": r["t"],
                "voi": publico(R), "nucleo": publico(Rn),
                "campo_voi": campo_vs_ref(r, ref, 0.0, L, n, m),
                "campo_nucleo": campo_vs_ref(r, ref, NUCLEO, L - NUCLEO, n,
                                             m)}
            print("hex8", n, control, f"E_app {r['E_app']:.2f}",
                  f"E_voi {R['E']:.2f} E_nuc {Rn['E']:.2f}",
                  f"p99 {R['p99']:.2f}/{Rn['p99']:.2f}", flush=True)
    guardar("hex8", out)
    return out


#: Candidatos de malla suave para A3 (PRERREGISTRO.md)
MALLAS_TET = {"omision": {}, "taubin5": {"suavizado": 5, "decimado": 0},
              "sin_suavizado": {"suavizado": 0, "decimado": 0}}


def etapa_tet10(ns=(32, 48), candidatos=None):
    previo = RES / "tet10.json"
    out = json.loads(previo.read_text()) if previo.exists() else {}
    for n in ns:
        BW, sp = casos.espinodoide(n)
        h = L / n
        for nombre, op in MALLAS_TET.items():
            if candidatos and nombre not in candidatos:
                continue
            t0 = time.perf_counter()
            try:
                malla = fem.mallar(BW, sp, "tet10", **op)
            except Exception as e:                      # noqa: BLE001
                out[f"n{n}_{nombre}"] = {"error": f"{type(e).__name__}: {e}"}
                continue
            tm = time.perf_counter() - t0
            inf = malla["informe"]
            for control in ("fuerza", "plato"):
                r = ensayo(malla, control, desp=-0.5 * h)
                R = region(r, 0.0, L)
                Rn = region(r, NUCLEO, L - NUCLEO)
                out[f"n{n}_{nombre}_{control}"] = {
                    "E_app": r["E_app"], "gdl": r["gdl"], "t": r["t"],
                    "t_malla": tm, "perdida_pct": inf.get("perdida_pct"),
                    "n_elems": inf["n_elems"],
                    "voi": publico(R), "nucleo": publico(Rn)}
                print("tet10", n, nombre, control, r["gdl"],
                      f"E_app {r['E_app']:.2f} E_voi {R['E']:.2f}",
                      f"E_nuc {Rn['E']:.2f} p99 {R['p99']:.2f}", flush=True)
                guardar("tet10", out)
    return out


if __name__ == "__main__":
    etapas = sys.argv[1:] or ["referencia", "hex8", "tet10"]
    if "referencia" in etapas:
        etapa_referencia()
    if "hex8" in etapas:
        etapa_hex8()
    if "tet10" in etapas:
        etapa_tet10(ns=(32,))
    if "tet10_48" in etapas:
        etapa_tet10(ns=(48,), candidatos=sys.argv[sys.argv.index("tet10_48")
                                                    + 1:] or None)
