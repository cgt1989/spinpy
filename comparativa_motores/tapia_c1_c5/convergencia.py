"""
convergencia.py: el protocolo de Tapia a varias resoluciones (NGSolve, directo).

    python convergencia.py [56 64 80]

La revision independiente mostro que los resultados a 40^3 NO estan convergidos
(C3 spinodoide: E_app -30 % de 40^3 a 64^3, p99 de von Mises +23 %). Este guion
repite las 15 estructuras a mas resoluciones, con la MISMA formulacion que
`correr.py` (base empotrada, 100 N, 18 GPa, nu = 0,30, hex8), y guarda una fila
por corrida en `resultados/convergencia.jsonl`.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent.parent))
sys.path.insert(0, str(AQUI.parent / "tapia"))

from spinpy import fem                                          # noqa: E402
import replicar_tapia as RT                                     # noqa: E402

VOIS = ("C1", "C2", "C3", "C4", "C5")
ESTR = ("voi", "spinodoide", "dual-lattice")


def main(res, solver="directo", solo=None):
    prot = fem.protocolo("tapia2026")
    sal = AQUI / "resultados" / "convergencia.jsonl"
    for n in res:
        for v in VOIS:
            for est in ESTR:
                if solo and (v, est) != solo:
                    continue
                BW, sp = RT.cargar(v, est)
                t0 = time.perf_counter()
                B, s2 = fem._remuestrear(BW, sp, n)
                malla = fem.mallar(B, s2, "hex8")
                reg = fem.ensayo(malla, prot, analisis=("lineal",),
                                 motor="ngsolve", solver=solver,
                                 aislado=False)
                lin = reg.get("lineal")
                fila = {"voi": v, "estructura": est, "n": n,
                        "h_mm": float(malla["spacing"][0]),
                        "n_gdl": int(3 * malla["nodos"].shape[0]),
                        "BVTV_malla": float(B.mean()),
                        "tiempo_s": time.perf_counter() - t0,
                        "fallos": reg["fallos"], "solver": solver,
                        "residuo_rel": reg.get("motor", {}).get("residuo_rel")}
                if lin and "_campos_lineal" in reg:
                    u = reg["_campos_lineal"]["u"]
                    vm, um, _ = RT.campos_elemento(malla, u, prot["E_s"],
                                                   prot["nu"])
                    fila.update(
                        E_app_MPa=lin["E_app"] * 1e-6,
                        E_app_nodal_MPa=lin["E_app_nodal"] * 1e-6,
                        dF_rel=lin["dF_rel"],
                        vm_p99_MPa=float(np.percentile(vm, 99)),
                        vm_p95_MPa=float(np.percentile(vm, 95)),
                        vm_max_MPa=float(vm.max()),
                        u_p99_mm=float(np.percentile(um, 99)),
                        u_max_mm=float(np.linalg.norm(u, axis=1).max()))
                with open(sal, "a", encoding="utf-8") as fh:
                    fh.write(json.dumps(fila, ensure_ascii=False,
                                        default=float) + "\n")
                print(f"n={n} {v} {est:13s} {fila['tiempo_s']:6.1f} s "
                      + (f"E_app {fila['E_app_MPa']:8.1f}  p99 "
                         f"{fila['vm_p99_MPa']:7.1f}" if lin else "FALLO"),
                      flush=True)


if __name__ == "__main__":
    # `python convergencia.py 80 iter`: CG + pyamg y un proceso por caso,
    # para acotar la memoria (el directo de n=80 agota 15 GB)
    a = sys.argv[1:]
    modo = "iterativo" if "iter" in a else "directo"
    ns = [int(x) for x in a if x.isdigit()] or [56, 64, 80]
    if modo == "iterativo" and "caso" not in a:
        import subprocess
        for n in ns:
            for v in VOIS:
                for e in ESTR:
                    subprocess.run([sys.executable, __file__, str(n), "iter",
                                    "caso", v, e], check=False)
    elif "caso" in a:
        i = a.index("caso")
        main(ns, modo, (a[i + 1], a[i + 2]))
    else:
        main(ns, modo)
