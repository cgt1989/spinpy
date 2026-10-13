"""
inicio_rapido.py: del VOI a la morfometria y al ensayo de compresion, en un minuto.

Usa los VOIs porcinos de `datos/porcino_koria2020/` (Koria, Mengoni y Brockett
2020, CC BY 4.0), asi que corre en cualquier instalacion sin datos propios.

    python ejemplos/inicio_rapido.py                  # VOI_V1, cubo central de 64^3
    python ejemplos/inicio_rapido.py --voi C3 --lado 96
    python ejemplos/inicio_rapido.py --rapido         # lo usa la integracion continua

Que hace, en orden:
  1. Lee el VOI y comprueba su huella SHA-256 contra `SHA256SUMS.txt`.
  2. Mide la morfometria (BV/TV, Tb.Th, Tb.Sp, DA, ...).
  3. Ensaya a compresion en z y da el modulo aparente y la carga de fallo de
     Pistoia.
  4. Imprime la procedencia: version, condiciones y huella del VOI.

Sobre la lectura de los numeros mecanicos: son resultados de simulacion con la
condicion de carga declarada abajo, verificados numericamente (ver README); no
son predicciones validadas contra un ensayo fisico de este hueso. El recorte
central de 64^3 es una submuestra del VOI de 188^3: sirve para comprobar la
instalacion y para aprender el flujo, no para citar valores. Para el VOI
completo use `--lado 187`.
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import time
from pathlib import Path

import numpy as np

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ))

from spinpy import (__version__, criterio_pistoia, ensayo_compresion,  # noqa: E402
                    leer_voi, morfometria)

DATOS = RAIZ / "datos" / "porcino_koria2020"


def sha256(ruta):
    h = hashlib.sha256()
    with open(ruta, "rb") as f:
        for bloque in iter(lambda: f.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()


def huellas_esperadas():
    out = {}
    for linea in (DATOS / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
        if linea.strip():
            h, nombre = linea.split()
            out[nombre] = h
    return out


def recorte_central(BW, lado):
    lado = min(int(lado), *BW.shape)
    ini = [(n - lado) // 2 for n in BW.shape]
    return BW[ini[0]:ini[0] + lado, ini[1]:ini[1] + lado, ini[2]:ini[2] + lado]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--voi", default="V1",
                    help="C1..C5 o V1..V5 (o una ruta a un .mat/.vtk)")
    ap.add_argument("--lado", type=int, default=64,
                    help="lado en voxeles del cubo central (187 = completo)")
    ap.add_argument("--rapido", action="store_true",
                    help="cubo de 40^3 y sin la tabla extendida")
    a = ap.parse_args(argv)
    if a.rapido:
        a.lado = 40

    ruta = Path(a.voi)
    if not ruta.suffix:
        ruta = DATOS / f"VOI_{a.voi.upper()}.mat"
    if not ruta.exists():
        print(f"No existe {ruta}", file=sys.stderr)
        return 2

    huella = sha256(ruta)
    esperada = huellas_esperadas().get(ruta.name)
    estado = ("coincide con SHA256SUMS.txt" if esperada == huella
              else "no figura en SHA256SUMS.txt" if esperada is None
              else "NO COINCIDE con SHA256SUMS.txt")
    print(f"spinpy {__version__} | {ruta.name} | SHA-256 {huella[:16]}... {estado}")
    if esperada not in (None, huella):
        return 1

    BW, sp = leer_voi(ruta)
    print(f"VOI completo: {BW.shape} voxeles de {sp[0] * 1000:.1f} um, "
          f"BV/TV = {BW.mean():.3f}")
    sub = recorte_central(BW, a.lado)
    print(f"Submuestra central: {sub.shape}")

    t0 = time.time()
    m = morfometria(sub, sp)
    print(f"\nMorfometria ({time.time() - t0:.1f} s)")
    for k, u in (("BVTV", ""), ("TbTh", "mm"), ("TbSp", "mm"), ("TbN", "1/mm"),
                 ("DA", ""), ("ConnD", "1/mm3")):
        v = m.get(k)
        if v is not None and np.isfinite(v):
            print(f"  {k:6s} {v:10.4g} {u}")

    t0 = time.time()
    r = ensayo_compresion(sub, sp, E_s=20e9, nu_s=0.3, sigma0=1e6,
                          apoyo="deslizante")
    if not r.get("ok"):
        print(f"\nEl ensayo no se pudo resolver: {r.get('msg')}", file=sys.stderr)
        return 1
    print(f"\nCompresion en z, apoyo deslizante, E_s = 20 GPa, nu = 0.3, "
          f"1 MPa ({time.time() - t0:.1f} s)")
    print(f"  E_app            {r['E_app'] / 1e6:10.1f} MPa")
    p = criterio_pistoia(r)
    if p.get("ok"):
        print(f"  carga de fallo   {p['sigma_fallo'] / 1e6:10.2f} MPa (criterio "
              f"de Pistoia: 2 % del tejido a {p['eps_crit']:.4f})")
    print("\nLeer con cautela: simulacion con tejido homogeneo e isotropo; sin "
          "contraste fisico para este hueso (README, 'Validacion').")
    return 0


if __name__ == "__main__":
    sys.exit(main())
