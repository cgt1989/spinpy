"""
app.py — El resolvedor actual de spinpy (`resistencia`, `elastic`).

Es la LINEA BASE. Solo resuelve lo que la app sabe resolver hoy:
  * ensayo de compresion lineal sobre la malla de voxeles (hex8), con
    `resistencia.ensayo_compresion` tal cual lo llama el panel de analisis;
  * homogeneizacion periodica hex8 con `elastic.homogeneizar`.
No tiene elementos TET10 ni analisis no lineal: esos casos se declaran «no
disponible» en la tabla, que es justamente la brecha que hoy cubre FEBio.

Los tiempos por etapa salen de los mensajes de progreso de la propia funcion
(«Ensamblando», «Resolviendo», «Recuperando»), sin tocar su codigo.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from spinpy.elastic import homogeneizar                          # noqa: E402
from spinpy.resistencia import ensayo_compresion                # noqa: E402

from ._base import Cronometro, principal                        # noqa: E402

ETAPAS = {"Ensamblando": "montaje", "Resolviendo": "solucion",
          "Recuperando": "post", "Listo": "post"}


def resolver(p):
    m = p["meta"]
    if m["analisis"] == "homogeneizacion":
        return _homogeneizar(p)
    if m["tipo"] != "hex8" or m["analisis"] != "lineal":
        raise NotImplementedError("la app solo resuelve hex8 lineal")
    t = Cronometro()

    def progreso(_f, msg):
        for k, v in ETAPAS.items():
            if msg.startswith(k):
                t(v)

    t("montaje")
    r = ensayo_compresion(p["BW"], m["spacing"], E_s=m["E"], nu_s=m["nu"],
                          sigma0=m["sigma_ref"], apoyo=m["apoyo"],
                          progreso=progreso)
    tiempos = t.fin()
    if not r["ok"]:
        raise RuntimeError(r["msg"])
    u = np.asarray(r["u"]).reshape(-1, 3)
    if u.shape[0] != p["nodos"].shape[0]:
        raise RuntimeError("la malla de la app no coincide con la del caso")
    return u, {"motor": "app", "solver": r.get("solver"),
               "residuo_rel": r.get("residuo_rel"), "tiempos": tiempos,
               "n_gdl": int(r["n_dof"])}, None


def _homogeneizar(p):
    m = p["meta"]
    t = Cronometro()
    t("total")
    C, info = homogeneizar(p["BW"], E_s=m["E"], nu_s=m["nu"],
                           vox_size=m["spacing"], tol=m.get("tol", 1e-8))
    tiempos = t.fin()
    if not info["ok"]:
        raise RuntimeError(info["msg"])
    return np.zeros((1, 3)), {"motor": "app", "solver": info["solver"],
                              "residuo_rel": info.get("residuo_rel"),
                              "tiempos": tiempos,
                              "n_gdl": int(info["n_dof"])}, {"C": C}


if __name__ == "__main__":
    principal(resolver)
