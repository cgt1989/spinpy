"""
p6_espinodoide.py: Error de discretizacion del hex8 y del hex8i en el
espinodoide de referencia, con la GEOMETRIA FIJA.

    python comparativa_motores/libros/p6_espinodoide.py [n_base ...]

Al subir la resolucion del generador cambia tambien la estructura (cada voxel
se clasifica de nuevo como hueso o vacio), de modo que la serie de E_app
frente a n mezcla dos efectos. Aqui la mascara se genera una sola vez a
n_base y cada voxel se divide en m^3 subvoxeles (m = 1, 2, 3, 4): la geometria
es identica y lo unico que cambia es la discretizacion. La referencia es la
extrapolacion de Richardson de las dos mallas mas finas de cada elemento.

Ensayo de la app (`resistencia.ensayo_compresion`): traccion uniforme en el
techo, base deslizante, E_s = 20 GPa, nu = 0,3. Escribe
resultados/p6_espinodoide.json.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parents[1]))
sys.path.insert(0, str(AQUI.parent))
from spinpy import morfometria                                  # noqa: E402
from spinpy.resistencia import ensayo_compresion                # noqa: E402
import casos                                                    # noqa: E402

MS = tuple(int(x) for x in os.environ.get("P6_MS", "1,2,3,4").split(","))


def subdividir(BW, m):
    return np.repeat(np.repeat(np.repeat(BW, m, 0), m, 1), m, 2)


def main(bases):
    out = {}
    archivo = AQUI / "resultados" / "p6_espinodoide.json"
    if archivo.exists():
        out = json.loads(archivo.read_text())
    for n in bases:
        BW, sp = casos.espinodoide(n)
        mo = morfometria(BW, sp)
        d = out.setdefault(str(n), {"BVTV": float(mo["BVTV"]),
                                    "TbTh_h": float(mo["TbTh"] / sp[0])})
        for el in ("hex8", "hex8i"):
            for m in MS:
                clave = f"{el}_m{m}"
                if clave in d:
                    continue
                t0 = time.perf_counter()
                r = ensayo_compresion(subdividir(BW, m), sp / m,
                                      E_s=20e9, elemento=el)
                d[clave] = {"E_app": float(r["E_app"]),
                            "t_s": time.perf_counter() - t0,
                            "n_dof": int(r.get("n_dof", 0) or 0)}
                print(n, clave, f"{r['E_app'] / 1e6:.3f} MPa",
                      f"{d[clave]['t_s']:.1f} s", flush=True)
                archivo.parent.mkdir(exist_ok=True)
                archivo.write_text(json.dumps(out, indent=1))
        for el in ("hex8", "hex8i"):
            a = d[f"{el}_m{MS[-2]}"]["E_app"]
            b = d[f"{el}_m{MS[-1]}"]["E_app"]
            r = MS[-1] / MS[-2]
            d[f"{el}_richardson"] = b + (b - a) / (r ** 2 - 1)
        archivo.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]] or [24])
