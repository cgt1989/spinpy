"""
m_app.py — El resolvedor propio de spinpy (`resistencia.ensayo_compresion`).

Es el motor de siempre y el validado: hexaedros de un voxel, elasticidad
lineal, control por fuerza, LU directo por debajo de 6000 GDL y CG con
multigrid algebraico (pyamg, modos rigidos) por encima. No tiene elementos
TET10, ni plato rigido, ni analisis no lineal: esas combinaciones lanzan
`NoDisponible`.
"""

from __future__ import annotations

import numpy as np

from ._comun import Cronometro, ErrorMotor, NoDisponible, solucion

NOMBRE = "App (spinpy)"
LICENCIA = "MIT"
ETAPAS = {"Ensamblando": "montaje", "Resolviendo": "solucion",
          "Recuperando": "post", "Listo": "post"}


def version():
    from .. import __version__
    return __version__


def resolver(p):
    from ..resistencia import ensayo_compresion
    m = p["meta"]
    if m["tipo"] != "hex8" or m["analisis"] != "lineal" \
            or m["control"] != "fuerza" or "BW" not in p:
        raise NoDisponible("la app resuelve solo hex8, lineal, con fuerza")
    t = Cronometro()

    def progreso(_f, msg):
        for k, v in ETAPAS.items():
            if msg.startswith(k):
                t(v)

    t("montaje")
    # `cargas` ya viene como presion sobre el area osea; la app recibe la
    # tension aparente sobre la seccion bruta.
    sigma = m["cargas"][-1] * m["A_osea_techo"] / m["A_bruta"]
    r = ensayo_compresion(p["BW"], m["spacing"], E_s=m["E"], nu_s=m["nu"],
                          sigma0=sigma, apoyo=m["apoyo"], progreso=progreso)
    tiempos = t.fin()
    if not r["ok"]:
        raise ErrorMotor(r["msg"])
    u = np.asarray(r["u"]).reshape(-1, 3)
    if u.shape[0] != p["nodos"].shape[0]:
        raise ErrorMotor("la malla de la app no coincide con la del problema")
    return solucion(u, [float(r["F_total"])],
                    {"motor": "app", "solver": r.get("solver"),
                     "residuo_rel": r.get("residuo_rel"), "tiempos": tiempos,
                     "n_gdl": int(r["n_dof"]), "iteraciones_newton": None})
