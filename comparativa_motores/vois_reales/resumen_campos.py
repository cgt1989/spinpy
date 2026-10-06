"""
resumen_campos.py: Cifras del texto sobre los campos de tension.

    python resumen_campos.py

Lee resultados/campos.json y escribe resultados/resumen_campos.json con, por
VOI y por ensayo del cubo aislado (traccion, plato), los cocientes frente a
la referencia de la ecuacion (11): en la capa lateral exterior (d < 2h), su
media cuadratica en todas las capas laterales y el von Mises del hueso en el
5 % superior e inferior de la altura.
"""

import json
from pathlib import Path

import numpy as np

RES = Path(__file__).resolve().parent / "resultados"
cam = json.loads((RES / "campos.json").read_text())
out = {}
for nombre, d in cam.items():
    ref = d["referencia"]
    o = {"frac_portante": {c: d[c]["frac_portante"]
                           for c in ("traccion", "plato", "referencia")}}
    for caso in ("traccion", "plato"):
        q = d[caso]
        lat = np.divide(q["lateral"], ref["lateral"])
        ver = np.divide(q["vertical"], ref["vertical"])
        n5 = max(1, round(0.05 * len(ver)))
        o[caso] = {"lateral_exterior": float(lat[0]),
                   "lateral_rms": float(np.sqrt(np.mean((lat - 1) ** 2))),
                   "vm_techo": float(np.mean(ver[-n5:])),
                   "vm_base": float(np.mean(ver[:n5])),
                   "vm_centro": float(np.mean(ver[len(ver) // 4:
                                                  3 * len(ver) // 4]))}
    out[nombre] = o
(RES / "resumen_campos.json").write_text(json.dumps(out, indent=1))
for k, v in out.items():
    print(k, {c: {x: round(y, 3) for x, y in v[c].items()}
              for c in ("traccion", "plato")}, v["frac_portante"])
