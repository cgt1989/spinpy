"""
p1_homogeneizacion.py: Espacio casi nulo del multigrid de la homogeneizacion.

    python comparativa_motores/libros/p1_homogeneizacion.py [n]

Espinodoides isotropos (conos 90, 90, 90; numero de onda 10 pi; semilla 1)
con rho de 0,15 a 0,40, homogeneizados con `elastic.homogeneizar` y las tres
variantes de `espacio_nulo`: 'escalar' (el comportamiento hasta la V2.1.1),
'traslaciones' y 'rigidos'. Se registran el residuo final, las iteraciones
del CG y el tiempo, y se compara el tensor con el del LU directo cuando el
tamano lo permite (n <= 24). Tolerancia de la app: 1e-8.

Escribe resultados/p1_homogeneizacion_n{n}.json.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parents[1]))
from spinpy import elastic, generar_mascara                     # noqa: E402

RHOS = (0.15, 0.20, 0.25, 0.30, 0.40)
VARIANTES = ("escalar", "traslaciones", "rigidos")


def main(n):
    out = {}
    archivo = AQUI / "resultados" / f"p1_homogeneizacion_n{n}.json"
    for rho in RHOS:
        BW, _, _ = generar_mascara(resolution=n, rho=rho, thetas=(90, 90, 90),
                                   wave_number=10 * np.pi, num_waves=700,
                                   seed=1)
        d = out.setdefault(str(rho), {"rho_real": float(BW.mean())})
        if n <= 24:
            umbral = elastic.UMBRAL_DIRECTO
            elastic.UMBRAL_DIRECTO = 10 ** 9      # fuerza el LU directo
            t0 = time.perf_counter()
            C_lu, inf = elastic.homogeneizar(BW, E_s=1.0, nu_s=0.3)
            elastic.UMBRAL_DIRECTO = umbral
            d["lu"] = {"C11": float(C_lu[0, 0]),
                       "t_s": time.perf_counter() - t0}
        for v in VARIANTES:
            umbral = elastic.UMBRAL_DIRECTO
            elastic.UMBRAL_DIRECTO = 0                # fuerza el iterativo
            elastic.UMBRAL_RESPALDO_LU = 0            # y sin respaldo LU
            t0 = time.perf_counter()
            C, inf = elastic.homogeneizar(BW, E_s=1.0, nu_s=0.3,
                                          espacio_nulo=v)
            elastic.UMBRAL_DIRECTO = umbral
            r = {"ok": bool(inf.get("ok")), "residuo": inf.get("residuo_rel"),
                 "iteraciones": inf.get("iteraciones"),
                 "t_s": time.perf_counter() - t0,
                 "C11": float(C[0, 0]) if np.isfinite(C[0, 0]) else None,
                 "msg": inf.get("msg", "")}
            if "lu" in d and r["C11"] is not None:
                r["dif_rel_lu"] = float(np.abs(C - C_lu).max()
                                        / np.abs(C_lu).max())
            d[v] = r
            print(n, rho, v, f"res={r['residuo']:.1e}",
                  f"it={max(r['iteraciones'] or [0])}",
                  f"t={r['t_s']:.1f}s", r.get("dif_rel_lu", ""), flush=True)
            archivo.parent.mkdir(exist_ok=True)
            archivo.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 16)
