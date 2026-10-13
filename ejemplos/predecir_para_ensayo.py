"""
predecir_para_ensayo.py: calcula la rigidez aparente de cada VOI para contrastarla.

Escribe un CSV con las columnas de `spinpy-validar-exp` (id, magnitud,
predicho, medido). La columna `medido` queda VACIA: se rellena con el ensayo
fisico del mismo espécimen, y el CSV se pasa despues a

    spinpy-validar-exp predicciones.csv --ccc-min 0.9 --mape-max 0.25

Los criterios (`--ccc-min`, `--mape-max`, ...) se declaran ANTES de mirar los
datos medidos y se registran en el protocolo; fijarlos despues de ver el
resultado invalida el contraste.

La prediccion solo es comparable con el ensayo si se reproducen sus
condiciones: el mismo modulo del tejido (`--Es`, MPa), la misma condicion de
apoyo y la misma region del hueso. Los valores por omision son los de la
aplicacion (E_s = 20 GPa, plato deslizante) y NO equivalen a un ensayo real.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ))

from spinpy import ensayo_compresion, leer_voi  # noqa: E402

DATOS = RAIZ / "datos" / "porcino_koria2020"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--salida", default="predicciones.csv")
    ap.add_argument("--lado", type=int, default=96,
                    help="lado en voxeles del cubo central (187 = completo)")
    ap.add_argument("--Es", type=float, default=20000.0, help="MPa")
    ap.add_argument("--apoyo", default="deslizante")
    ap.add_argument("--vois", nargs="*",
                    default=[f"{g}{i}" for g in "CV" for i in range(1, 6)])
    a = ap.parse_args(argv)

    filas = []
    for nombre in a.vois:
        BW, sp = leer_voi(DATOS / f"VOI_{nombre}.mat")
        n = min(a.lado, *BW.shape)
        i0 = [(m - n) // 2 for m in BW.shape]
        sub = BW[i0[0]:i0[0] + n, i0[1]:i0[1] + n, i0[2]:i0[2] + n]
        r = ensayo_compresion(sub, sp, E_s=a.Es * 1e6, nu_s=0.3, sigma0=1e6,
                              apoyo=a.apoyo)
        E = r["E_app"] / 1e6 if r.get("ok") else ""
        print(f"VOI_{nombre}: lado {n}, E_app = {E if E == '' else round(E, 1)} MPa")
        filas.append({"id": f"VOI_{nombre}", "magnitud": "E_app_MPa",
                      "predicho": E, "medido": ""})
    with open(a.salida, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["id", "magnitud", "predicho", "medido"])
        w.writeheader()
        w.writerows(filas)
    print(f"Escrito {a.salida}. Complete la columna 'medido' con el ensayo fisico.")


if __name__ == "__main__":
    main()
