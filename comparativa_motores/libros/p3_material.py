"""
p3_material.py: St. Venant-Kirchhoff frente a neo-Hookeano en compresion
uniaxial, con las soluciones cerradas de la comparativa (`comun.uniaxial`).

    python comparativa_motores/libros/p3_material.py

Para cada acortamiento se da la tension nominal y la rigidez tangente dP/dl
de los dos materiales (E = 20 GPa, nu = 0,3). El SVK alcanza un maximo de
carga en lambda = 1/sqrt(3) (42 % de acortamiento) y su tension vuelve a cero
cuando lambda -> 0: es la inestabilidad en compresion que describen Kamensky
(2022, sec. 4.3.4) y Bathe (1996, sec. 6.6.1). Escribe
resultados/p3_material.json.
"""

import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))
from comun import uniaxial                                       # noqa: E402


def main():
    filas = []
    for acort in (0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 1 - 1 / np.sqrt(3),
                  0.50, 0.70):
        lz = 1 - acort
        h = 1e-6
        f = {"acortamiento": acort}
        for mat in ("svk", "neohookeano"):
            P, _ = uniaxial(mat, lz)
            dP = (uniaxial(mat, lz + h)[0] - uniaxial(mat, lz - h)[0]) / (2 * h)
            f[mat] = {"P_MPa": float(P), "tangente_rel_E": float(dP / 20000.0)}
        f["cociente_P"] = f["svk"]["P_MPa"] / f["neohookeano"]["P_MPa"]
        filas.append(f)
        print(f"{acort:6.1%}  P svk {f['svk']['P_MPa']:10.1f}  "
              f"neo {f['neohookeano']['P_MPa']:10.1f}  "
              f"svk/neo {f['cociente_P']:.3f}  tangente svk "
              f"{f['svk']['tangente_rel_E']:.3f}  neo "
              f"{f['neohookeano']['tangente_rel_E']:.3f}")
    (AQUI / "resultados").mkdir(exist_ok=True)
    (AQUI / "resultados" / "p3_material.json").write_text(
        json.dumps(filas, indent=1))


if __name__ == "__main__":
    main()
