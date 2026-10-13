"""
correr.py: protocolo de Tapia et al. (2026) sobre C1 a C5, sus spinodoides y
sus dual-lattice (15 estructuras).

    python correr.py [C1 C2 C3 C4 C5]

Malla de ladrillos a 40^3, base EMPOTRADA (ux = uy = uz = 0), cara superior con
x e y libres y 100 N repartidos en cargas nodales iguales en -z, E_s = 18 GPa,
nu = 0,3, elasticidad lineal (`fem.PROTOCOLOS_FEM["tapia2026"]`). Se resuelve
con la app y con NGSolve (resolvedor directo, la referencia).

Guarda por elemento la tension de von Mises (MPa) y la deformacion total |u|
(mm, media de los 8 nodos del elemento) en `campos/<VOI>_<estructura>.npz`, y
una fila por corrida en `resultados/tapia_c1_c5.jsonl`.

Las estructuras ajustadas las genera `../tapia/ajustar.py`; los VOIs son los de
`datos/porcino_koria2020/` (Koria, Mengoni y Brockett 2020, CC BY 4.0).
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

from spinpy import fem, motores                                 # noqa: E402
from spinpy.resistencia import _PERM, _solo_portante            # noqa: E402
import replicar_tapia as RT                                     # noqa: E402

N = 40
CORRIDAS = (("app", "auto", "app"), ("ngsolve", "directo", "referencia"))
ESTRUCTURAS = ("voi", "spinodoide", "dual-lattice")


def main(vois):
    (AQUI / "campos").mkdir(exist_ok=True)
    (AQUI / "resultados").mkdir(exist_ok=True)
    prot = fem.protocolo("tapia2026")
    sal = AQUI / "resultados" / "tapia_c1_c5.jsonl"
    for v in vois:
        for est in ESTRUCTURAS:
            BW, sp = RT.cargar(v, est)
            B, s2 = fem._remuestrear(BW, sp, N)
            malla = fem.mallar(B, s2, "hex8")
            BW_hex = _solo_portante(np.transpose(B, _PERM[2]))
            cen = malla["nodos"][malla["elems"]].mean(axis=1)
            ijk = np.floor(cen / malla["spacing"]).astype(np.int16)
            out = {"ijk": ijk, "forma": np.array(malla["forma"]),
                   "spacing": np.asarray(malla["spacing"], float)}
            for motor, solver, k in CORRIDAS:
                t0 = time.perf_counter()
                reg = fem.ensayo(malla, prot, analisis=("lineal",),
                                 motor=motor, solver=solver, aislado=True,
                                 BW=BW_hex if motor == "app" else None)
                t = time.perf_counter() - t0
                lin = reg.get("lineal")
                fila = {"voi": v, "estructura": est, "clave": k,
                        "motor": motor, "n": N,
                        "n_elems": int(malla["elems"].shape[0]),
                        "n_gdl": int(3 * malla["nodos"].shape[0]),
                        "BVTV_malla": float(B.mean()),
                        "lado_mm": float(malla["spacing"][0] * N),
                        "tiempo_s": t, "version": motores.version(motor),
                        "fallos": reg["fallos"],
                        "procedencia_malla": malla["huella"]}
                if lin and "_campos_lineal" in reg:
                    u = reg["_campos_lineal"]["u"]
                    vm, um, _ = RT.campos_elemento(malla, u, prot["E_s"],
                                                   prot["nu"])
                    out[f"vm_{k}"], out[f"u_{k}"] = vm, um
                    fila.update(E_app_MPa=lin["E_app"] * 1e-6,
                                dF_rel=lin["dF_rel"],
                                uz_max_mm=lin["u_max_mm"],
                                u_tot_max_mm=float(np.linalg.norm(
                                    u, axis=1).max()),
                                vm_max_MPa=float(vm.max()),
                                vm_p99_MPa=float(np.percentile(vm, 99)),
                                vm_media_MPa=float(vm.mean()))
                with open(sal, "a", encoding="utf-8") as fh:
                    fh.write(json.dumps(fila, ensure_ascii=False,
                                        default=float) + "\n")
                print(f"{v} {est:13s} {k:10s} {t:6.1f} s "
                      + (f"E_app {fila['E_app_MPa']:8.1f} MPa  vm_max "
                         f"{fila['vm_max_MPa']:7.1f}  |u|max "
                         f"{fila['u_tot_max_mm']:.4f} mm" if lin
                         else f"FALLO {reg['fallos']}"), flush=True)
            np.savez_compressed(AQUI / "campos" / f"{v}_{est}.npz", **out)


if __name__ == "__main__":
    main(sys.argv[1:] or ["C1", "C2", "C3", "C4", "C5"])
