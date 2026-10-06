"""
p2_newton.py: Newton globalizado (busqueda lineal sobre la energia y corte del
incremento de carga) frente al Newton completo de la V2.1.1, en NGSolve.

    python comparativa_motores/libros/p2_newton.py [bloque] [espinodoide]

1. BLOQUE. Bloque macizo de 4 x 4 x 6 voxeles (0,05 mm), plato rigido con
   5, 10, 15 y 20 % de acortamiento, SVK y neo-Hookeano: la globalizacion no
   debe cambiar la fuerza, que se compara con la solucion cerrada.
2. ESPINODOIDE. Espinodoide de referencia de la comparativa a 24^3 (hex8) y
   con malla suave TET10, fuerza impuesta sobre el techo en UN SOLO
   incremento (la opcion por omision de la app, pasos=1) con tensiones
   aparentes crecientes. Se registra si cada variante converge, cuantos
   cortes de carga y pasos recortados necesito y la fuerza final.

Escribe resultados/p2_newton.json.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parents[1]))
sys.path.insert(0, str(AQUI.parent))
from spinpy import fem, motores                                 # noqa: E402
from spinpy.motores._comun import ErrorMotor, problema_de_malla  # noqa: E402
import casos                                                    # noqa: E402
from comun import uniaxial                                      # noqa: E402

E_MPA, NU = 20000.0, 0.3
ARCHIVO = AQUI / "resultados" / "p2_newton.json"


def correr(malla, control, cargas, material, globalizar):
    p = problema_de_malla(malla, E_MPA, NU, analisis="nl", control=control,
                          cargas=cargas, material=material, solver="directo")
    p["meta"]["globalizar"] = globalizar
    t0 = time.perf_counter()
    try:
        out = motores.resolver(p, "ngsolve")
        m = out["meta"]
        u = np.asarray(out["u"])
        return {"ok": True, "F_N": [float(f) for f in out["F_reac"]],
                "uz_techo_mm": float(u[p["techo_nodos"], 2].mean()),
                "iter": m.get("iteraciones_newton"),
                "newton": m.get("newton"), "t_s": time.perf_counter() - t0,
                "A_bruta": p["meta"]["A_bruta"]}
    except ErrorMotor as e:
        return {"ok": False, "msg": str(e), "t_s": time.perf_counter() - t0}


def bloque(out):
    BW = np.ones((4, 4, 6), bool)
    malla = fem.mallar(BW, np.full(3, 0.05), "hex8")
    eps = [0.05, 0.10, 0.15, 0.20]
    for mat in ("svk", "neohookeano"):
        exacta = [-uniaxial(mat, 1 - e)[0] for e in eps]
        for g in (False, True):
            r = correr(malla, "plato", eps, mat, g)
            if r["ok"]:
                r["dif_rel_cerrada"] = float(max(
                    abs(f / r["A_bruta"] - x) / abs(x)
                    for f, x in zip(r["F_N"], exacta)))
            out[f"bloque_{mat}_{'glob' if g else 'v211'}"] = r
            print("bloque", mat, g, r.get("ok"), r.get("dif_rel_cerrada"),
                  r.get("iter"), flush=True)


def espinodoide(out, tipo, n, cargas_MPa, materiales=("svk", "neohookeano")):
    BW, sp = casos.espinodoide(n)
    malla = fem.mallar(BW, sp, tipo)
    for mat in materiales:
        for s in cargas_MPa:
            for g in (False, True):
                clave = f"esp_{tipo}_{n}_{mat}_{s:g}MPa_{'glob' if g else 'v211'}"
                if clave in out:
                    continue
                r = correr(malla, "fuerza", [s], mat, g)
                out[clave] = r
                print(clave, r.get("ok"), r.get("newton"),
                      f"{r['t_s']:.0f}s", flush=True)
                ARCHIVO.write_text(json.dumps(out, indent=1))


def main(que):
    ARCHIVO.parent.mkdir(exist_ok=True)
    out = json.loads(ARCHIVO.read_text()) if ARCHIVO.exists() else {}
    if "bloque" in que:
        bloque(out)
        ARCHIVO.write_text(json.dumps(out, indent=1))
    if "espinodoide" in que:
        espinodoide(out, "hex8", 24, (2, 5, 10, 20))
    if "tet10" in que:
        espinodoide(out, "tet10", 24, (2, 5, 10))
    ARCHIVO.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:] or ["bloque", "espinodoide"])
