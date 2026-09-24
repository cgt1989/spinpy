"""barrido.py — Rigidez no lineal frente a la carga, con fuerza y con plato.

    cd Port_Python
    python comparativa_febio/barrido.py                   # proximal 32 y 48

Complementa a `comparar_febio.py` (una sola carga no lineal, 1 MPa) y a
`plato.py` (dos): recorre la carga aparente de 10 kPa a 3 MPa con los DOS
controles de carga y guarda E_nl / E_lin en cada punto. Es la figura que
dice a partir de que carga deja de valer el ensayo lineal, y si eso depende
de la estructura o de como se la carga.

* Fuerza: E_lin es el E_app de spinpy (coincide con FEBio extrapolado a carga
  nula a 3e-7, ver RESULTADOS.md); E_nl = sigma / (media de uz del techo / H).
* Plato: E_lin por extrapolacion de dos deformaciones pequenas; E_nl con la
  fuerza de la integral de volumen de sigma_zz. La deformacion de cada punto
  es la del ensayo lineal con fuerza a esa tension, para comparar a la misma
  carga.

Un punto que FEBio no converge se guarda como None: tambien es un resultado
(por encima de esa carga, el equilibrio con carga muerta deja de existir o el
Newton no lo encuentra).

Resultado: `resultados/barrido.jsonl`, una linea por caso.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
import comparar_febio as cf                                  # noqa: E402
from spinpy.escribe import escribir_febio                    # noqa: E402
from spinpy.resistencia import _solo_portante, ensayo_compresion  # noqa: E402
from spinpy.solido import malla_hex                          # noqa: E402

CARGAS = (1e4, 1e5, 3e5, 1e6, 3e6)          # Pa


def caso(clave, dir_corridas):
    BW, spc, apoyo, _f = cf.CASOS[clave]()
    nx, ny, _nz = BW.shape
    A = float(nx * spc[0] * ny * spc[1])
    res = ensayo_compresion(BW, spc, E_s=cf.E_S, nu_s=cf.NU_S,
                            sigma0=cf.SIGMA_SPINPY, apoyo=apoyo)
    E_lin = res["E_app"] / 1e6                                # MPa
    nodos, elems, _ = malla_hex(_solo_portante(BW), spc)
    z = nodos[:, 2]
    H = float(z.max() - z.min())
    techo = z >= z.max() - 1e-9 * max(H, 1.0)
    V_e = float(np.prod(spc))
    d = dir_corridas / f"{clave}_barrido"
    d.mkdir(parents=True, exist_ok=True)
    print(f"\n== {clave}: E_lin (fuerza) {E_lin:.1f} MPa")

    def corre(nombre, **kw):
        feb = d / f"{nombre}.feb"
        escribir_febio(nodos, elems, feb, E_s=cf.E_S, nu_s=cf.NU_S,
                       A_bruta=A, apoyo=apoyo, **kw)
        try:
            cf.correr_febio(feb)
        except RuntimeError:
            return None
        return feb

    def E_fuerza(sigma):
        feb = corre(f"fuerza_{sigma:.0f}", sigma_app=sigma)
        if feb is None:
            return None
        u = cf._ultimo_registro(feb.with_name(feb.stem + "_u.txt"), 3)
        return sigma / 1e6 / (abs(u[techo, 2].mean()) / H)

    def E_plato(eps):
        feb = corre(f"plato_{eps:.3e}", eps_plato=eps)
        if feb is None:
            return None
        s = cf._ultimo_registro(feb.with_name(feb.stem + "_s.txt"), 6)
        return -s[:, 2].sum() * V_e / H / A / eps

    eps0 = 1e3 / 1e6 / E_lin
    Ep_lin = 2 * E_plato(eps0) - E_plato(2 * eps0)
    fila = {"caso": clave, "BVTV": float(BW.mean()), "E_lin_fuerza_MPa": E_lin,
            "E_lin_plato_MPa": Ep_lin, "cargas_Pa": list(CARGAS),
            "fuerza": [], "plato": []}
    for sigma in CARGAS:
        Ef = E_fuerza(sigma)
        Ep = E_plato(sigma / 1e6 / E_lin)
        fila["fuerza"].append(None if Ef is None else Ef / E_lin)
        fila["plato"].append(None if Ep is None else Ep / Ep_lin)
        print(f"   {sigma / 1e6:5.2f} MPa   fuerza {fila['fuerza'][-1]}   "
              f"plato {fila['plato'][-1]}")
    return fila


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--casos", nargs="*", default=["proximal_32",
                                                   "proximal_48"])
    ap.add_argument("--salida",
                    default=str(AQUI / "resultados" / "barrido.jsonl"))
    ap.add_argument("--corridas", default=str(AQUI / "corridas"))
    a = ap.parse_args()
    for clave in a.casos:
        fila = caso(clave, Path(a.corridas))
        with open(a.salida, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(fila, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
