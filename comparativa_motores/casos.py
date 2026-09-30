"""
casos.py — Los problemas de la comparativa, todos construidos por spinpy.

Cada funcion devuelve un caso de `comun` (malla + apoyo + carga) listo para
los cinco motores. La geometria sale SIEMPRE del codigo de la app (mascara
del generador, `malla_hex`, `febio.mallar`), de modo que la comparativa mide
los motores sobre las mallas que la aplicacion produce de verdad.

Espinodoide de referencia: rho = 0.30, numero de onda 12 pi, 700 ondas, conos
(30, 30, 90) grados, semilla 1, lado 5 mm (el del VOI proximal de H4). Con la
misma semilla el campo aleatorio es el MISMO a cualquier resolucion: al subir
n cambia la discretizacion, no la estructura.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))
sys.path.insert(0, str(AQUI.parent / "comparativa_febio_tet"))

import comun                                                    # noqa: E402
from spinpy import febio, generar_mascara                       # noqa: E402

LADO_MM = 5.0
ESPINODOIDE = dict(wave_number=12 * np.pi, num_waves=700, thetas=(30, 30, 90),
                   rho=0.30, seed=1)


def espinodoide(n):
    BW, _, info = generar_mascara(resolution=int(n), **ESPINODOIDE)
    return BW, np.full(3, LADO_MM / n)


def bloque_hex(**kw):
    return comun.caso_hex(np.ones((8, 8, 12), bool), 0.05, "bloque_hex", **kw)


def bloque_tet():
    m = febio.mallar(np.ones((8, 8, 12), bool), np.full(3, 0.05), "tet10")
    return comun.caso_tet(m, "bloque_tet")


def cavidad_hex(n=32):
    from cavidad import geometria
    BW, spc, _, _ = geometria(n)
    return comun.caso_hex(BW, spc, f"cavidad_hex_n{n}", n=n)


def cavidad_tet(res_esfera=64, tam=0.06):
    """TET10 sobre la ESFERA ANALITICA, la referencia de FEBio de
    `comparativa_febio_tet/cavidad.py` (misma funcion, mismos argumentos)."""
    from cavidad import malla_analitica
    m = malla_analitica(32, res_esfera, tam)
    return comun.caso_tet(m, f"cavidad_tet_esfera{res_esfera}",
                          res_esfera=res_esfera, tam_max_mm=tam)


def espinodoide_hex(n, **kw):
    BW, sp = espinodoide(n)
    return comun.caso_hex(BW, sp, f"espinodoide_hex_n{n}", n=int(n), **kw)


def espinodoide_tet(n):
    BW, sp = espinodoide(n)
    t0 = time.perf_counter()
    m = febio.mallar(BW, sp, "tet10")
    t = time.perf_counter() - t0
    inf = m["informe"]
    return comun.caso_tet(m, f"espinodoide_tet_n{n}", n=int(n),
                          tiempo_malla_s=t,
                          perdida_volumen_pct=inf.get("perdida_pct"),
                          descartado_pct=inf.get("descartado_pct"))


def no_lineal(caso, material, eps):
    caso["meta"].update(analisis="nl", material=material,
                        eps_plato=[float(e) for e in eps])
    caso["meta"]["nombre"] += f"_nl_{material}"
    return caso
