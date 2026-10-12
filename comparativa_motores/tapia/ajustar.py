"""
ajustar.py: Ajusta el spinodoide y el dual-lattice a un VOI real y guarda las
mascaras generadas, para ensayarlas despues con cada motor.

    python ajustar.py C2 C4

Configuracion del informe automatico de la app: modo completo, resolucion de
comparacion por omision (96^3 para un VOI de 188^3), semilla 20260720, sin
termino mecanico ni replicas de incertidumbre. La mascara se regenera desde los
parametros por la misma traduccion que usa la app (`informe._regenerar`), de
modo que `estructuras/<VOI>_<familia>.npy` es la que reproduce el paquete de
reproduccion. Escribe tambien `estructuras/<VOI>_<familia>.json` con los
parametros, el error del ajuste, la morfometria y la huella SHA-256.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent.parent))

from spinpy import leer_voi, morfometria                       # noqa: E402
from spinpy.fit import ajustar_spinodoide                      # noqa: E402
from spinpy.fit_dual import ajustar_dual_lattice               # noqa: E402
from spinpy.informe import _regenerar, huella_mascara          # noqa: E402

SEMILLA = 20260720
SALIDA = AQUI / "estructuras"


def main(vois):
    SALIDA.mkdir(exist_ok=True)
    for v in vois:
        BW, sp = leer_voi(AQUI / "vois" / f"VOI_{v}.mat")
        BW, sp = np.asarray(BW, bool), np.asarray(sp, float)
        m_voi = morfometria(BW, sp)
        for fam, fn in (("spinodoide", ajustar_spinodoide),
                        ("dual-lattice", ajustar_dual_lattice)):
            dest = SALIDA / f"{v}_{fam}"
            if dest.with_suffix(".npy").exists():
                print("ya existe", dest.name, flush=True)
                continue
            t0 = time.perf_counter()
            kw = dict(modo="completo", m_voi=m_voi, seed=SEMILLA, replicas=0)
            if fam == "spinodoide":
                kw["precision"] = "f32"
            r = fn(BW, sp, **kw)
            t = time.perf_counter() - t0
            par = r["parametros"]
            M = np.asarray(_regenerar(par), bool)
            np.save(dest.with_suffix(".npy"), M)
            doc = {"voi": v, "familia": fam, "semilla": SEMILLA,
                   "tiempo_ajuste_s": t, "error": r.get("error"),
                   "parametros": par, "sha256_mascara": huella_mascara(M),
                   "forma": list(M.shape), "BVTV_generada": float(M.mean()),
                   "BVTV_voi": float(BW.mean()),
                   "spacing_mm": float(par["spacing_mm"])}
            json.dump(doc, open(dest.with_suffix(".json"), "w",
                                encoding="utf-8"), indent=1, default=float)
            print(f"{v} {fam}: error {doc['error']}, BV/TV "
                  f"{doc['BVTV_generada']:.3f} (VOI {doc['BVTV_voi']:.3f}), "
                  f"{t:.0f} s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or ["C2", "C4"])
