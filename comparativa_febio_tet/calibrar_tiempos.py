"""calibrar_tiempos.py — Mide FEBio y tetgen para recalibrar `spinpy.tiempos`.

    cd Port_Python
    python comparativa_febio_tet/calibrar_tiempos.py            # 24 32 40
    python comparativa_febio_tet/calibrar_tiempos.py --n 24 32 40 48

CORRER CON EL EQUIPO TRANQUILO. Los coeficientes de FEBio que hay hoy en
`tiempos.py` se midieron con la CPU al 100 % por procesos ajenos y son
PROVISIONALES (los tiempos por corrida se apartaban hasta un 70 % del ajuste).

Sobre el VOI proximal de H4 remuestreado a cada n, con hex8 y TET10: mallado
(tetgen, proceso hijo incluido) y una corrida lineal de FEBio (3 iteraciones)
con todos los hilos menos uno, como la GUI. Anota gdl, memoria de pico y
tiempos en `resultados/calibracion.jsonl` y ajusta al final las leyes de
potencia que usa `tiempos` (s por iteracion y MB frente a gdl; s de mallado
frente a TET10), para copiarlas a mano en `tiempos.py` junto con la fecha.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))
from spinpy import febio, leer_voi                  # noqa: E402
from spinpy.elastic import remuestrear_bw           # noqa: E402

VOI = AQUI.parent.parent / "H4" / "Segmentadas" / "VOI_proximal_cubico.vtk"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--n", nargs="*", type=int, default=[24, 32, 40])
    ap.add_argument("--salida", default=str(AQUI / "resultados"
                                            / "calibracion.jsonl"))
    a = ap.parse_args()
    BW0, sp0 = leer_voi(VOI)
    hilos = max(1, (os.cpu_count() or 2) - 1)
    carpeta = AQUI / "corridas" / "calibracion"
    filas = []
    for n in a.n:
        B, sp = remuestrear_bw(BW0, sp0, n)
        for tipo in ("hex8", "tet10"):
            t0 = time.perf_counter()
            try:
                m = febio.mallar_aislado(B, sp, tipo)
            except febio.ErrorMalla as e:
                print(f"n {n} {tipo}: {e}")
                continue
            t_malla = time.perf_counter() - t0
            feb = carpeta / f"{tipo}_{n}" / "c.feb"
            feb.parent.mkdir(parents=True, exist_ok=True)
            febio._escribir(m, feb, febio.protocolo("app"), sigma=1e3)
            r = febio.correr(feb, hilos=hilos)
            f = {"n": n, "tipo": tipo, "n_elems": int(m["elems"].shape[0]),
                 "gdl": int(3 * m["nodos"].shape[0]), "t_malla_s": t_malla,
                 "t_febio_s": r["tiempo_s"], "iteraciones": r["iteraciones"],
                 "memoria_MB": r["memoria_MB"], "hilos": hilos,
                 "fecha": time.strftime("%Y-%m-%d %H:%M")}
            filas.append(f)
            print(json.dumps(f), flush=True)
            with open(a.salida, "a", encoding="utf-8") as fh:
                fh.write(json.dumps(f) + "\n")
    for tipo, ref in (("hex8", 1e5), ("tet10", 3e5)):
        ff = [f for f in filas if f["tipo"] == tipo]
        if len(ff) < 2:
            continue
        g = np.array([f["gdl"] for f in ff], float) / ref
        it = np.array([f["t_febio_s"] / max(f["iteraciones"] or 3, 1)
                       for f in ff])
        mem = np.array([f["memoria_MB"] for f in ff])
        p1 = np.polyfit(np.log(g), np.log(it), 1)
        p2 = np.polyfit(np.log(g), np.log(mem), 1)
        print(f"{tipo}: s/iteracion = {np.exp(p1[1]):.3g} (gdl/{ref:g})^"
              f"{p1[0]:.3f};  memoria = {np.exp(p2[1]):.0f} MB (gdl/{ref:g})^"
              f"{p2[0]:.3f}")
    ft = [f for f in filas if f["tipo"] == "tet10"]
    if len(ft) >= 2:
        x = np.array([f["n_elems"] for f in ft], float) / 1.79e5
        y = np.array([f["t_malla_s"] for f in ft])
        p = np.polyfit(np.log(x), np.log(y), 1)
        print(f"mallado TET10 (con proceso hijo): {np.exp(p[1]):.3g} s "
              f"(n_tet/1.79e5)^{p[0]:.3f}")


if __name__ == "__main__":
    main()
