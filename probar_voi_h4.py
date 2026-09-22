"""
probar_voi_h4.py — El port sobre los VOIs reales de H4.

POR QUE ESTE CONTRASTE ES MAS DURO QUE EL DEL GENERADOR
-------------------------------------------------------
`comparar.py` enfrenta dos implementaciones de un proceso ESTOCASTICO, asi que
el criterio tiene que ser estadistico (|z| frente a la banda de ruido). Aqui el
VOI es un dato fijo leido de disco: la morfometria es deterministica y las dos
implementaciones deben coincidir hasta la aritmetica. No hay varianza donde
esconder una discrepancia.

Cualquier diferencia en BV/TV, Tb.Sp, DA o la direccion principal es un error
real del port. Lo unico que se admite es el sesgo conocido del variante de
marching cubes en BS (sub-1%, creciente con la densidad), que se propaga a
BS/BV, Tb.Th y Tb.N por ser funciones de BS.

Uso:  python probar_voi_h4.py        (requiere resultados/referencia_voi_h4.json)
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from spinpy import leer_vtk_voi, morfometria   # noqa: E402

RAIZ = Path(__file__).parent
SALIDA = RAIZ / "resultados"
REF = SALIDA / "referencia_voi_h4.json"
VOIDIR = RAIZ.parent / "H4" / "Segmentadas"

# Tolerancias. Las metricas que NO dependen de la superficie deben coincidir
# practicamente al bit; las que si dependen arrastran el sesgo de marching
# cubes y se les concede 1.5%.
TOL_EXACTA = 1e-9
TOL_MC = 0.015
DEPENDEN_DE_BS = {"BS", "BSBV", "BSTV", "TbTh", "TbSp", "TbN"}


def _aplanar(c):
    while isinstance(c, list) and len(c) == 1 and isinstance(c[0], list):
        c = c[0]
    return list(c)


def main():
    if not REF.exists():
        print(f"Falta {REF}.\nEjecuta primero:  matlab -batch \"cd('{RAIZ}'); ref_voi_h4\"")
        return 2

    ref = json.loads(REF.read_text(encoding="utf-8"))
    campos = _aplanar(ref["campos"])
    bloques = ref["resultados"]
    if isinstance(bloques, dict):
        bloques = [bloques]

    print("=" * 88)
    print(" EL PORT SOBRE LOS VOIs REALES DE H4")
    print(" Morfometria deterministica: se espera coincidencia exacta salvo")
    print(" el sesgo del variante de marching cubes en BS y sus derivadas.")
    print("=" * 88)

    salida = {"campos": campos, "vois": []}
    peor_exacta = 0.0
    peor_mc = 0.0
    n_fallos = 0

    for blk in bloques:
        nombre = blk["archivo"]
        ruta = VOIDIR / nombre
        med_mat = np.array([np.nan if v is None else float(v)
                            for v in np.ravel(_aplanar(blk["metricas"]))])
        dir_mat = np.array(np.ravel(_aplanar(blk["dir_principal"])), dtype=float)

        t0 = time.time()
        VOI, spacing = leer_vtk_voi(ruta)
        m = morfometria(VOI, spacing)
        t = time.time() - t0

        print(f"\n--- {nombre} ---")
        print(f"    dims={VOI.shape}  spacing={spacing[0]:.6f} mm  "
              f"(python {t:.1f} s, matlab {float(np.ravel(blk['tiempo_s'])[0]):.1f} s)")
        print(f"    {'metrica':10s} {'MATLAB':>14s} {'Python':>14s} "
              f"{'dif rel':>11s}   tolerancia")

        fila = {"archivo": nombre, "dims": list(VOI.shape),
                "spacing": spacing.tolist(), "matlab": [], "python": [],
                "dif_rel": []}

        for q, c in enumerate(campos):
            vm = med_mat[q]
            vp = m.get(c, np.nan)
            vp = float(vp) if np.isscalar(vp) or isinstance(vp, float) else np.nan
            fila["matlab"].append(None if not np.isfinite(vm) else vm)
            fila["python"].append(None if not np.isfinite(vp) else vp)

            if not np.isfinite(vm) or not np.isfinite(vp):
                fila["dif_rel"].append(None)
                print(f"    {c:10s} {'—' if not np.isfinite(vm) else f'{vm:14.7g}':>14s} "
                      f"{vp:>14.7g} {'—':>11s}   (no en applib)")
                continue

            d = abs(vp - vm) / abs(vm) if vm != 0 else abs(vp)
            fila["dif_rel"].append(float(d))
            tol = TOL_MC if c in DEPENDEN_DE_BS else TOL_EXACTA
            if c in DEPENDEN_DE_BS:
                peor_mc = max(peor_mc, d)
            else:
                peor_exacta = max(peor_exacta, d)
            ok = d <= tol
            if not ok:
                n_fallos += 1
            marca = " " if ok else "*"
            etiq = "MC (1.5%)" if c in DEPENDEN_DE_BS else "exacta"
            print(f"  {marca} {c:10s} {vm:>14.7g} {vp:>14.7g} "
                  f"{100*d:>10.6f}%   {etiq}")

        # Direccion principal: se compara como eje, no como vector (el signo de
        # un autovector es arbitrario), de ahi el valor absoluto del producto.
        dir_py = np.asarray(m["dir_principal"], dtype=float)
        cos = abs(float(dir_py @ dir_mat))
        ang = np.degrees(np.arccos(np.clip(cos, -1, 1)))
        fila["dir_matlab"] = dir_mat.tolist()
        fila["dir_python"] = dir_py.tolist()
        fila["angulo_deg"] = float(ang)
        marca = " " if ang < 2.0 else "*"
        if ang >= 2.0:
            n_fallos += 1
        print(f"  {marca} dir.princ. MATLAB=[{dir_mat[0]:+.4f} {dir_mat[1]:+.4f} "
              f"{dir_mat[2]:+.4f}]  Python=[{dir_py[0]:+.4f} {dir_py[1]:+.4f} "
              f"{dir_py[2]:+.4f}]  angulo={ang:.3f} deg")

        salida["vois"].append(fila)

    print("\n" + "=" * 88)
    print(" RESUMEN")
    print("=" * 88)
    print(f"  Peor diferencia en metricas independientes de BS : {100*peor_exacta:.3e} %")
    print(f"  Peor diferencia en metricas que dependen de BS   : {100*peor_mc:.4f} %")
    print(f"  Comprobaciones fuera de tolerancia               : {n_fallos}")

    salida["resumen"] = {"peor_exacta_pct": 100 * peor_exacta,
                         "peor_mc_pct": 100 * peor_mc, "n_fallos": n_fallos}
    (SALIDA / "comparacion_voi_h4.json").write_text(
        json.dumps(salida, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n  Resultados -> {SALIDA / 'comparacion_voi_h4.json'}")
    return 1 if n_fallos else 0


if __name__ == "__main__":
    sys.exit(main())
