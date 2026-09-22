"""
validar_ajuste.py — Contrasta el ajuste portado contra el de la app de MATLAB.

QUE SE PUEDE VALIDAR Y QUE NO
------------------------------
El generador es estocastico, y MATLAB (Mersenne Twister) y NumPy (PCG64) tienen
flujos de numeros aleatorios sin ninguna relacion: la misma semilla entera NO
produce los mismos numeros. Cada plataforma evalua por tanto REALIZACIONES
DISTINTAS en los mismos puntos de la rejilla.

Eso parte la validacion en dos mitades de naturaleza opuesta, y mezclarlas
seria enganarse:

  DETERMINISTA — debe coincidir EXACTAMENTE:
      * rejilla de densidades de la etapa A          (formula sobre BV/TV)
      * rejilla de numeros de onda                   (rango fijo, G2)
      * presets de angulos conicos y su orden        (lista fija, G1+K2)
      * matriz de rotacion R_fit                     (eje principal del VOI)
      * numero de evaluaciones previstas
    Nada de esto depende del azar. Si difiere, el port esta mal.

  ESTOCASTICO — solo comparable de forma descriptiva:
      * el error alcanzado y los parametros ganadores

  Para lo segundo, la referencia util no es "¿coincide?" sino la banda de ruido
  del propio generador, ya medida en Validacion_Anexo. Dos ejecuciones de
  MATLAB con semillas distintas tampoco coinciden entre si.

Uso:  python validar_ajuste.py   (requiere resultados/ajuste_matlab.json)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from spinpy import ajustar_spinodoide, leer_voi, morfometria   # noqa: E402

RAIZ = Path(__file__).parent
SALIDA = RAIZ / "resultados"
REF = SALIDA / "ajuste_matlab.json"
VOIDIR = RAIZ.parent / "H4" / "Segmentadas"

THETAS0 = [15.0, 15.0, 15.0]
NUM_WAVES = 700
CAMPOS = ["BVTV", "DA", "BSBV", "TbTh", "TbSp", "TbN", "PoTot"]


def _apl(c):
    while isinstance(c, list) and len(c) == 1 and isinstance(c[0], list):
        c = c[0]
    return list(c)


def _arr(v):
    return np.array([np.nan if x is None else float(x)
                     for x in np.ravel(_apl(v) if isinstance(v, list) else [v])])


def main():
    if not REF.exists():
        print(f"Falta {REF}.\nEjecuta:  matlab -batch \"cd('{RAIZ}'); validar_ajuste_matlab\"")
        return 2

    ref = json.loads(REF.read_text(encoding="utf-8"))
    bloques = ref["resultados"]
    if isinstance(bloques, dict):
        bloques = [bloques]

    print("=" * 84)
    print(" AJUSTE: PORT A PYTHON  vs  APP DE MATLAB")
    print(f" modo={ref['modo']}   thetas0={_arr(ref['thetas0']).astype(int).tolist()}"
          f"   numWaves0={int(_arr(ref['numWaves0'])[0])}")
    print("=" * 84)

    fallos = 0
    salida = {"casos": []}

    for blk in bloques:
        nombre = str(_apl([blk["archivo"]])[0] if isinstance(blk["archivo"], list)
                     else blk["archivo"])
        print(f"\n{'-'*84}\n {nombre}\n{'-'*84}")

        VOI, spacing = leer_voi(VOIDIR / nombre)
        m_voi = morfometria(VOI, spacing)

        # alinear_fabrica=False A PROPOSITO: la correccion L1 no existe en
        # `AppFinal_V2.m`, asi que activarla aqui compararia dos cosas
        # distintas. Este script contrasta el port contra MATLAB; para usar la
        # correccion esta el ajuste normal, que la trae activada.
        r = ajustar_spinodoide(VOI, spacing, modo="completo", m_voi=m_voi,
                               num_waves=NUM_WAVES, precision="f32",
                               alinear_fabrica=False)
        rej = r["rejillas"]

        # ---------- parte DETERMINISTA ----------
        print("\n  DETERMINISTA (debe coincidir exactamente)")
        comprobaciones = []

        dm = _arr(blk["rejillaDensidad"]); dp = np.asarray(rej["densidad"])
        comprobaciones.append(("rejilla de densidad", dm, dp))

        wm = _arr(blk["rejillaOnda"]); wp = np.asarray(rej["onda"])
        comprobaciones.append(("rejilla de numero de onda", wm, wp))

        pm = np.asarray([_arr(p) for p in _apl(blk["presetsTheta"])], dtype=float)
        pp = np.asarray(rej["presets_theta"], dtype=float)
        comprobaciones.append(("presets de theta", pm.ravel(), pp.ravel()))

        Rm = np.asarray(_apl(blk["R"]), dtype=float).reshape(3, 3)
        Rp = np.asarray(r["parametros"]["R"], dtype=float)
        comprobaciones.append(("matriz R_fit", Rm.ravel(), Rp.ravel()))

        nm = float(_arr(blk["nEvalPrevistas"])[0])
        npy = float(rej["n_eval_previstas"])
        comprobaciones.append(("nº de evaluaciones previstas",
                               np.array([nm]), np.array([npy])))

        for etiqueta, a, b in comprobaciones:
            if a.shape != b.shape:
                print(f"    [FALLA] {etiqueta:32s} formas distintas "
                      f"{a.shape} vs {b.shape}")
                fallos += 1
                continue
            d = np.nanmax(np.abs(a - b)) if a.size else 0.0
            ok = d <= 1e-12
            fallos += (not ok)
            print(f"    [{'OK  ' if ok else 'FALLA'}] {etiqueta:32s} "
                  f"dif. maxima = {d:.3e}")

        # ---------- parte ESTOCASTICA ----------
        print("\n  ESTOCASTICO (realizaciones distintas; solo descriptivo)")
        em = float(_arr(blk["errorFinal"])[0])
        ep = float(r["error"])
        print(f"    error final          MATLAB {em:.6f}   Python {ep:.6f}")
        print(f"    densidad ganadora    MATLAB {float(_arr(blk['densidad'])[0]):.4f}"
              f"   Python {r['parametros']['densidad']:.4f}")
        print(f"    numero de onda       MATLAB {float(_arr(blk['waveNum'])[0]):g}"
              f"      Python {r['parametros']['wave_number_pi']:g}")
        print(f"    thetas ganadores     MATLAB {_arr(blk['thetas']).astype(int).tolist()}"
              f"   Python {[int(t) for t in r['parametros']['thetas']]}")

        mg = blk["metricasGanador"]
        print(f"\n    {'metrica':9s} {'VOI':>10s} {'MATLAB':>10s} {'Python':>10s}")
        for c in CAMPOS:
            v = float(m_voi[c])
            a = float(_arr([mg[c]])[0]) if c in mg else np.nan
            b = float(r["metricas_spin"][c])
            print(f"    {c:9s} {v:10.5f} {a:10.5f} {b:10.5f}")

        salida["casos"].append({
            "archivo": nombre,
            "matlab": {"error": em, "densidad": float(_arr(blk["densidad"])[0]),
                       "wave": float(_arr(blk["waveNum"])[0]),
                       "thetas": _arr(blk["thetas"]).tolist()},
            "python": {"error": ep, **{k: r["parametros"][k] for k in
                                       ("densidad", "wave_number_pi", "thetas")}},
            "k2": r["k2"],
            "n_eval": r["n_evaluaciones"], "tiempo_s": r["tiempo_s"],
        })

    print("\n" + "=" * 84)
    print(f"  Comprobaciones deterministas fuera de tolerancia: {fallos}")
    print("=" * 84)
    (SALIDA / "comparacion_ajuste.json").write_text(
        json.dumps(salida, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n  Resultados -> {SALIDA / 'comparacion_ajuste.json'}")
    return 1 if fallos else 0


if __name__ == "__main__":
    sys.exit(main())
