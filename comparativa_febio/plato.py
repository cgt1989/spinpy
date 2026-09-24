"""plato.py — ¿Es el ablandamiento no lineal cosa de los voladizos del techo?

    cd Port_Python
    python comparativa_febio/plato.py                     # proximal 32 y 48
    python comparativa_febio/plato.py --casos proximal_48

LA HIPOTESIS
------------
`comparar_febio.py` encontro que, a la carga de referencia de la app (1 MPa
con fuerza impuesta), FEBio no lineal da en el VOI proximal un E_app 16 %
(32^3) y 48 % (48^3) mas bajo que el lineal, y que el desvio se concentra en
trabeculas cortadas por la cara superior del VOI: reciben la carga en su
extremo libre y trabajan en voladizo.

Si es asi, un PLATO RIGIDO sin friccion —el mismo uz en todo el techo— debe
hacer desaparecer casi todo el ablandamiento, porque esos extremos ya no
pueden ceder por su cuenta.

PREDICCION, ESCRITA ANTES DE CORRER ESTE SCRIPT
-----------------------------------------------
Con plato, a la deformacion aparente equivalente a 1 MPa del ensayo con
fuerza, |E_nl / E_lin - 1| < 2 % en proximal 32 y 48 (frente a 16 y 48 %).
Si sale del mismo orden que con fuerza, la hipotesis de los voladizos es
FALSA y el ablandamiento es de la estructura, no de la condicion de carga.
Se corre tambien a la deformacion equivalente a la carga de fallo de Pistoia
(sin prediccion: es lo que decide si esa carga lineal tiene sentido).

COMO SE MIDE
------------
* E_lin con plato: dos deformaciones pequenas (las equivalentes a 1 y 2 kPa)
  y extrapolacion a cero, como en comparar_febio.py.
* La fuerza sale de la integral de volumen de sigma_zz (FEBio 4.5 escribe
  reacciones nulas). Con tension de Cauchy sobre el volumen de referencia el
  error es del orden de la deformacion del tejido: ~1e-3, frente a efectos
  del 10 %.
* E con plato y E con fuerza NO tienen por que coincidir ni en lineal: son
  condiciones de contorno distintas. Se reporta su cociente, que mide cuanto
  pesan los voladizos en el E_app lineal que publica la app.

Resultado: una linea por caso en `resultados/plato.jsonl`.
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
from spinpy import procedencia                               # noqa: E402
from spinpy.escribe import escribir_febio                    # noqa: E402
from spinpy.resistencia import (_solo_portante, criterio_pistoia,  # noqa: E402
                                ensayo_compresion)
from spinpy.solido import malla_hex                          # noqa: E402


def caso(clave, dir_corridas):
    BW, spc, apoyo, familia = cf.CASOS[clave]()
    nx, ny, nz = BW.shape
    A = float(nx * spc[0] * ny * spc[1])
    print(f"\n== {clave} con plato: {BW.shape}, BV/TV {BW.mean():.4f}")

    res = ensayo_compresion(BW, spc, E_s=cf.E_S, nu_s=cf.NU_S,
                            sigma0=cf.SIGMA_SPINPY, apoyo=apoyo)
    E_fuerza = res["E_app"]
    k_pistoia = criterio_pistoia(res)["factor"]
    nodos, elems, _ = malla_hex(_solo_portante(BW), spc)
    H = float(nodos[:, 2].max() - nodos[:, 2].min())
    V_e = float(np.prod(spc))
    d = dir_corridas / f"{clave}_plato"
    d.mkdir(parents=True, exist_ok=True)

    def E_plato(eps):
        feb = d / f"plato_{eps:.3e}.feb"
        escribir_febio(nodos, elems, feb, E_s=cf.E_S, nu_s=cf.NU_S,
                       A_bruta=A, apoyo=apoyo, eps_plato=eps)
        try:
            cf.correr_febio(feb)
        except RuntimeError as e:
            print(f"   eps {eps:.2e}: {e}")
            return None
        s = cf._ultimo_registro(feb.with_name(feb.stem + "_s.txt"), 6)
        F = -s[:, 2].sum() * V_e / H                  # N (MPa * mm^2)
        return F / A / eps                            # MPa

    # Deformaciones equivalentes: las del ensayo LINEAL con fuerza a esas
    # tensiones. Asi el plato se compara a la misma "carga" que la app.
    eps_de = lambda sigma: sigma / E_fuerza
    E1, E2 = E_plato(eps_de(1e3)), E_plato(eps_de(2e3))
    E_lin = 2 * E1 - E2
    fila = {"caso": clave, "BVTV": float(BW.mean()), "apoyo": apoyo,
            "E_fuerza_lineal_MPa": E_fuerza / 1e6,
            "E_plato_lineal_MPa": E_lin,
            "cociente_plato_fuerza_lineal": E_lin / (E_fuerza / 1e6),
            "sigma_fallo_pistoia_MPa": k_pistoia * cf.SIGMA_SPINPY / 1e6}
    for nombre, sigma in (("1MPa", cf.SIGMA_SPINPY),
                          ("fallo_pistoia", k_pistoia * cf.SIGMA_SPINPY)):
        Enl = E_plato(eps_de(sigma))
        fila[f"eps_{nombre}"] = eps_de(sigma)
        fila[f"E_plato_nl_{nombre}_MPa"] = Enl
        fila[f"dE_plato_nl_{nombre}"] = (None if Enl is None
                                         else Enl / E_lin - 1.0)

    # El dato del ensayo con fuerza, de la corrida principal, para la tabla.
    fuente = AQUI / "resultados" / "febio.jsonl"
    if fuente.exists():
        for l in fuente.read_text(encoding="utf-8").splitlines():
            r = json.loads(l)
            if r["caso"] == clave:
                fila["dE_fuerza_nl_1MPa"] = r.get("dE_nl_1MPa")
    fila["procedencia"] = procedencia.bloque(familia)
    fila["procedencia"]["febio"] = "4.5.0"

    print(f"   E lineal: fuerza {fila['E_fuerza_lineal_MPa']:.1f}  plato "
          f"{E_lin:.1f} MPa  (cociente {fila['cociente_plato_fuerza_lineal']:.3f})")
    print(f"   no lineal a 1 MPa: plato {fila['dE_plato_nl_1MPa']}  fuerza "
          f"{fila.get('dE_fuerza_nl_1MPa')}")
    print(f"   no lineal a la carga de Pistoia "
          f"({fila['sigma_fallo_pistoia_MPa']:.2f} MPa): plato "
          f"{fila['dE_plato_nl_fallo_pistoia']}")
    return fila


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--casos", nargs="*", default=["proximal_32",
                                                   "proximal_48"])
    ap.add_argument("--salida", default=str(AQUI / "resultados" / "plato.jsonl"))
    ap.add_argument("--corridas", default=str(AQUI / "corridas"))
    a = ap.parse_args()
    for clave in a.casos:
        fila = caso(clave, Path(a.corridas))
        with open(a.salida, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(fila, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
