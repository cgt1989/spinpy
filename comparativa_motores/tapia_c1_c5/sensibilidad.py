"""
sensibilidad.py: la carga del techo, con plato rigido y en el nucleo.

    python sensibilidad.py

La tension uniforme sobre las caras oseas del techo corta trabeculas en la
cima y las hace trabajar como voladizos. El propio repositorio documenta que
eso sesga E_app (-43 a -55 %) y el p99 (+33 a +53 %) frente a la configuracion
embebida (`fem._corregido`). Este guion corre, a 40^3 con NGSolve directo, el
`lineal_plato` (plato rigido sin friccion en el techo) y el E corregido
(nucleo con plato), para ver si cambia el ORDEN entre estructuras.
Escribe `resultados/sensibilidad.jsonl`.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent.parent))
sys.path.insert(0, str(AQUI.parent / "tapia"))

from spinpy import fem                                          # noqa: E402
import replicar_tapia as RT                                     # noqa: E402

N = 40


def main():
    prot = fem.protocolo("tapia2026")
    sal = AQUI / "resultados" / "sensibilidad.jsonl"
    for v in ("C1", "C2", "C3", "C4", "C5"):
        for est in ("voi", "spinodoide", "dual-lattice"):
            BW, sp = RT.cargar(v, est)
            B, s2 = fem._remuestrear(BW, sp, N)
            malla = fem.mallar(B, s2, "hex8")
            reg = fem.ensayo(malla, prot, analisis=("lineal", "lineal_plato"),
                             motor="ngsolve", solver="directo", aislado=False)
            lin, pl, co = (reg.get("lineal") or {}, reg.get("lineal_plato") or {},
                           reg.get("corregido") or {})
            fila = {"voi": v, "estructura": est, "n": N,
                    "E_app_traccion_MPa": lin.get("E_app", np.nan) * 1e-6,
                    "E_app_plato_MPa": pl.get("E_app", np.nan) * 1e-6,
                    "cociente_plato_traccion": pl.get("cociente_plato_fuerza"),
                    "E_corregido_MPa": co.get("E_app", np.nan) * 1e-6,
                    "metodo_E": co.get("metodo_E"),
                    "fallos": reg["fallos"]}
            with open(sal, "a", encoding="utf-8") as fh:
                fh.write(json.dumps(fila, ensure_ascii=False, default=float)
                         + "\n")
            print(f"{v} {est:13s} traccion {fila['E_app_traccion_MPa']:8.1f}  "
                  f"plato {fila['E_app_plato_MPa']:8.1f}  corregido "
                  f"{fila['E_corregido_MPa']:8.1f} ({fila['metodo_E']})",
                  flush=True)


if __name__ == "__main__":
    main()
