"""
replicar_tapia.py: Protocolo de Tapia et al. (2026) con cada motor FEM.

    python replicar_tapia.py [C2 C4]

Para cada VOI real y para su spinodoide y su dual-lattice ajustados
(`ajustar.py`): malla de ladrillos a 40^3 (la resolucion por omision de la
app), base empotrada, E_s = 18 GPa, nu = 0,3 y 100 N axiales sobre la seccion
bruta (`fem.PROTOCOLOS_FEM["tapia2026"]`), resuelta con la app y con los cuatro
motores externos. Cada motor corre en un proceso hijo con su resolvedor
automatico (lo que daria la app), y NGSolve ademas con el directo, que hace de
REFERENCIA de las diferencias.

Guarda por elemento la tension de von Mises (MPa) y la deformacion total |u|
(mm; media de los 8 nodos del elemento, la variable `campo_desp` de la app),
junto con el indice de voxel de cada elemento, en `campos/<VOI>_<estructura>.npz`
y una fila por corrida en `resultados/tapia.jsonl`.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent.parent))

from spinpy import fem, leer_voi, motores                       # noqa: E402
from spinpy.resistencia import _PERM, _solo_portante            # noqa: E402

N = 40
CORRIDAS = (("app", "auto"), ("ngsolve", "auto"), ("fenicsx", "auto"),
            ("skfem", "auto"), ("sfepy", "auto"), ("ngsolve", "directo"))
ESTRUCTURAS = ("voi", "spinodoide", "dual-lattice")


def clave(motor, solver):
    return motor if solver == "auto" else "referencia"


def cargar(v, est):
    if est == "voi":
        BW, sp = leer_voi(AQUI / "vois" / f"VOI_{v}.mat")
        return np.asarray(BW, bool), np.asarray(sp, float)
    M = np.load(AQUI / "estructuras" / f"{v}_{est}.npy")
    doc = json.load(open(AQUI / "estructuras" / f"{v}_{est}.json"))
    return M, np.full(3, doc["spacing_mm"])


def campos_elemento(malla, u, E_s, nu):
    s = fem.tension_elemental(malla["nodos"], malla["elems"], u, E_s, nu)
    vm = fem.von_mises(s) * 1e-6          # Pa -> MPa
    umag = np.linalg.norm(u, axis=1)[malla["elems"]].mean(axis=1)
    return vm, umag, s


def main(vois, claves=None):
    (AQUI / "campos").mkdir(exist_ok=True)
    (AQUI / "resultados").mkdir(exist_ok=True)
    prot = fem.protocolo("tapia2026")
    sal = AQUI / "resultados" / "tapia.jsonl"
    for v in vois:
        for est in ESTRUCTURAS:
            BW, sp = cargar(v, est)
            B, s2 = fem._remuestrear(BW, sp, N)
            malla = fem.mallar(B, s2, "hex8")
            BW_hex = _solo_portante(np.transpose(B, _PERM[2]))
            cen = malla["nodos"][malla["elems"]].mean(axis=1)
            ijk = np.floor(cen / malla["spacing"]).astype(np.int16)
            np.savez_compressed(
                AQUI / "campos" / f"{v}_{est}.npz", ijk=ijk,
                forma=np.array(malla["forma"]),
                spacing=np.asarray(malla["spacing"], float),
                huella=malla["huella"])
            for motor, solver in CORRIDAS:
                k = clave(motor, solver)
                if claves and k not in claves:
                    continue
                out = {}
                t0 = time.perf_counter()
                reg = fem.ensayo(malla, prot, analisis=("lineal",),
                                 motor=motor, solver=solver, aislado=True,
                                 BW=BW_hex if motor == "app" else None)
                t = time.perf_counter() - t0
                lin = reg.get("lineal")
                fila = {"voi": v, "estructura": est, "clave": k,
                        "motor": motor, "solver_pedido": solver,
                        "n": N, "n_elems": int(malla["elems"].shape[0]),
                        "n_gdl": int(3 * malla["nodos"].shape[0]),
                        "tiempo_pared_s": t,
                        "version": motores.version(motor),
                        "fallos": reg["fallos"],
                        "corridas": reg["corridas"],
                        "procedencia_malla": malla["huella"]}
                if lin and "_campos_lineal" in reg:
                    u = reg["_campos_lineal"]["u"]
                    vm, um, _ = campos_elemento(malla, u, prot["E_s"],
                                                prot["nu"])
                    out[f"vm_{k}"], out[f"u_{k}"] = vm, um
                    out[f"unod_{k}"] = u
                    fila.update(E_app_MPa=lin["E_app"] * 1e-6,
                                u_max_mm=lin["u_max_mm"], dF_rel=lin["dF_rel"],
                                **{kk: vv for kk, vv in lin["pistoia"].items()
                                   if kk in ("vm_p99_superficie", "vm_max",
                                             "factor_pistoia")})
                    fila["vm_max_MPa"] = float(vm.max())
                    fila["vm_p99_MPa"] = float(np.percentile(vm, 99))
                with open(sal, "a", encoding="utf-8") as fh:
                    fh.write(json.dumps(fila, ensure_ascii=False,
                                        default=float) + "\n")
                print(f"{v} {est:13s} {k:10s} {t:7.1f} s "
                      + (f"E_app {fila['E_app_MPa']:.1f} MPa  vm_p99 "
                         f"{fila['vm_p99_MPa']:.2f}" if lin else
                         f"FALLO {reg['fallos']}"), flush=True)
                if out:
                    np.savez_compressed(
                        AQUI / "campos" / f"{v}_{est}__{k}.npz", **out)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("vois", nargs="*", default=["C2", "C4"])
    ap.add_argument("--motores", nargs="+", default=None,
                    help="claves: app ngsolve fenicsx skfem sfepy referencia")
    a = ap.parse_args()
    main(a.vois, a.motores)
