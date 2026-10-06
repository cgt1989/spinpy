"""
p7_hex_orden2.py: Hexaedros cuadraticos (Q2) de NGSolve sobre la malla de
voxeles frente al hexaedro trilineal refinado.

    python comparativa_motores/libros/p7_hex_orden2.py [n_base ...]

Espinodoide de referencia de la comparativa generado a n_base; geometria
FIJA (cada voxel se divide en m^3 subvoxeles, como en p6_espinodoide.py).
Ensayo lineal con plato rigido (0,1 % de acortamiento) y base deslizante,
resuelto con NGSolve: trilineal (orden 1) con m = 1, 2, 3 y Q2 (orden 2) con
m = 1 y 2. Cada caso corre en un proceso aparte para medir su memoria de pico.
E_app = F / (A_bruta * eps). Escribe resultados/p7_hex_orden2.json.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parents[1]))
sys.path.insert(0, str(AQUI.parent))
ARCHIVO = AQUI / "resultados" / "p7_hex_orden2.json"
EPS = 1e-3
CASOS = ((1, 1), (1, 2), (1, 3), (2, 1), (2, 2))      # (orden, m)


def uno(n, orden, m, solver):
    from spinpy import fem, motores
    from spinpy.motores._comun import problema_de_malla
    import casos
    BW, sp = casos.espinodoide(n)
    if m > 1:
        BW = np.repeat(np.repeat(np.repeat(BW, m, 0), m, 1), m, 2)
    malla = fem.mallar(BW, sp / m, "hex8")
    p = problema_de_malla(malla, 20000.0, 0.3, control="plato", cargas=[EPS],
                          solver=solver)
    p["meta"]["orden_hex"] = orden
    t0 = time.perf_counter()
    r = motores.resolver(p, "ngsolve")
    F = abs(float(r["F_reac"][-1]))
    return {"E_app_MPa": F / (p["meta"]["A_bruta"] * EPS),
            "n_gdl": int(r["meta"]["n_gdl"]), "t_s": time.perf_counter() - t0,
            "rss_MB": motores.rss_pico_MB(),
            "solver": r["meta"].get("solver")}


def main(bases):
    ARCHIVO.parent.mkdir(exist_ok=True)
    out = json.loads(ARCHIVO.read_text()) if ARCHIVO.exists() else {}
    for n in bases:
        for orden, m in CASOS:
            clave = f"{n}_orden{orden}_m{m}"
            if clave in out:
                continue
            solver = "directo"
            cmd = [sys.executable, __file__, "--uno", str(n), str(orden),
                   str(m), solver]
            r = subprocess.run(cmd, capture_output=True, text=True)
            try:
                out[clave] = json.loads(r.stdout.strip().splitlines()[-1])
            except Exception:
                out[clave] = {"error": (r.stderr or r.stdout)[-800:]}
            print(clave, out[clave], flush=True)
            ARCHIVO.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--uno":
        n, orden, m, solver = sys.argv[2:6]
        print(json.dumps(uno(int(n), int(orden), int(m), solver)))
    else:
        main([int(x) for x in sys.argv[1:]] or [24])
