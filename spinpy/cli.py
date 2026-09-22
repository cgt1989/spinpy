"""cli.py — Entrada por linea de comandos, sin rutas propias de nadie.

POR QUE EXISTE. `correr_lote.py` sirve para trabajar, pero lleva el banco de
VOIs escrito dentro:

    BANCO = Path(r"C:\\Users\\carlo\\OneDrive\\Escritorio\\VOIs Caballos")

Eso esta bien para un script personal y es exactamente lo que impide que otra
persona use la herramienta: la primera ejecucion falla en una ruta que no
existe en su maquina, y el mensaje no le dice que la cambie. Aqui las rutas
son argumentos obligatorios.

    spinpy-lote VOI_*.vtk --replicas 20 --salida ./resultados

El script antiguo se conserva: tiene los valores por defecto del banco de este
proyecto y ahorra teclear. Los dos llaman al mismo `correr_lote`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main_lote(argv=None):
    p = argparse.ArgumentParser(
        prog="spinpy-lote",
        description="Ajusta un spinodoide a cada VOI, genera replicas "
                    "sinteticas y descompone la varianza.",
        epilog="AVISO: el numero de filas de la tabla NO es el N del estudio. "
               "Las replicas son tecnicas: reducen la incertidumbre de cada "
               "espécimen pero no anaden grados de libertad para comparar "
               "entre animales. Usa la columna n_efectivo.")
    p.add_argument("vois", nargs="+", type=Path,
                   help="archivos .vtk o .mat, o carpetas que los contengan")
    p.add_argument("--replicas", type=int, default=10)
    p.add_argument("--modo", default="completo", choices=["completo", "rapido"])
    p.add_argument("--mecanica", action="store_true",
                   help="anade homogeneizacion elastica (mucho mas lento)")
    p.add_argument("--res-homog", type=int, default=16)
    p.add_argument("--jobs", type=int, default=-1)
    p.add_argument("--salida", type=Path, default=Path("resultados"))
    p.add_argument("--familia", default="spinodoide",
                   choices=["spinodoide", "dual-lattice"])
    p.add_argument("--objetivo", type=Path, default=None,
                   help="JSON con opciones del ajuste: pesos, distancia, "
                        "replicas, seleccion, k_robusto, alineacion, "
                        "peso_mecanico. Ej.: {\"pesos\": {\"TbSp\": 0, "
                        "\"PoDm\": 1}, \"seleccion\": \"robusta\"}")
    a = p.parse_args(argv)

    opciones = None
    if a.objetivo is not None:
        import json
        opciones = json.loads(a.objetivo.read_text(encoding="utf-8"))
        if not isinstance(opciones, dict):
            p.error("--objetivo tiene que ser un objeto JSON")

    rutas = []
    for r in a.vois:
        if r.is_dir():
            rutas += sorted(list(r.glob("*.vtk")) + list(r.glob("*.mat")))
        elif r.exists():
            rutas.append(r)
        else:
            p.error(f"no existe: {r}")
    if not rutas:
        p.error("ningun VOI encontrado en lo indicado")

    from .estadistica import informe
    from .lote import METRICAS, correr_lote

    print(f"{len(rutas)} VOI(s), {a.replicas} replicas, modo {a.modo}")

    def avance(i, n, msg):
        print(f"  [{i}/{n}] {msg}", flush=True)

    df, aj, meta = correr_lote(
        rutas, n_replicas=a.replicas, modo=a.modo, mecanica=a.mecanica,
        res_homog=a.res_homog, n_jobs=a.jobs, progreso=avance,
        familia=a.familia, opciones_ajuste=opciones)
    inf = informe(df, [m for m in METRICAS if m in df.columns])

    a.salida.mkdir(parents=True, exist_ok=True)
    for nombre, tabla in (("especimenes", df), ("ajustes", aj),
                          ("varianza", inf)):
        destino = a.salida / f"lote_{nombre}.csv"
        tabla.to_csv(destino, index=False)
        print(f"  -> {destino}")

    n_ef = float(inf["n_efectivo"].min()) if len(inf) else float("nan")
    print(f"\n{len(df)} filas en la tabla, pero N EFECTIVO = {n_ef:.1f}. "
          f"Ese es el N que corresponde usar.")
    return 0


if __name__ == "__main__":
    sys.exit(main_lote())
