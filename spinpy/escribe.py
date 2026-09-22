"""
escribe.py — Escritores de malla para Abaqus, ANSYS, ParaView y STL.

SISTEMA DE UNIDADES: mm - N - MPa
----------------------------------
Ni Abaqus ni ANSYS llevan unidades; interpretan los numeros en el sistema
coherente que elija quien escribe el modelo. Las mallas de este proyecto estan
en MILIMETROS, porque el spacing del micro-CT lo esta, y el sistema coherente
con el milimetro es:

    longitud mm      fuerza N      tension y modulo MPa = N/mm^2

Las funciones reciben E en PASCALES —como lo maneja el resto del proyecto,
`localBaseMaterial` da 20e9— y lo escriben en MPa. Volcar 20e9 tal cual declara
un material un millon de veces mas rigido, y el analisis corre sin protestar:
solo se nota al ver desplazamientos absurdamente pequenos.

TIPOS DE ELEMENTO
-----------------
    hex8    Abaqus C3D8    ANSYS SOLID185
    tet10   Abaqus C3D10   ANSYS SOLID187
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pyvista as pv

TIPO_ABAQUS = {8: "C3D8", 10: "C3D10"}
TIPO_ANSYS = {8: "SOLID185", 10: "SOLID187"}
TIPO_VTK = {8: pv.CellType.HEXAHEDRON, 10: pv.CellType.QUADRATIC_TETRA}


def _caras_z(nodos, tol_rel=1e-9):
    """Nodos de la cara inferior y superior, como conjuntos con nombre."""
    z = np.asarray(nodos, float)[:, 2]
    zmin, zmax = float(z.min()), float(z.max())
    tol = max(tol_rel * (zmax - zmin), 1e-12)
    return (np.nonzero(z <= zmin + tol)[0], np.nonzero(z >= zmax - tol)[0])


def _lista(w, ids, por_linea=8):
    for i in range(0, len(ids), por_linea):
        w(", ".join(str(int(v)) for v in ids[i:i + por_linea]) + ",\n")


def escribir_stl(superficie, ruta):
    ruta = Path(ruta)
    superficie.save(str(ruta), binary=True)
    return {"ruta": str(ruta), "MB": ruta.stat().st_size / 1e6,
            "triangulos": int(superficie.n_cells)}


def escribir_vtu(nodos, elems, ruta, arrays=True):
    """Malla volumetrica en VTU, para inspeccionarla en ParaView o pyvista.

    Es el formato util para MIRAR la malla antes de mandarla a un solver: los
    .inp y .apdl son enormes y ningun visor los abre comodamente.

    POR QUE LLEVA ARRAYS. Una malla sin ningun campo escalar es geometria muda:
    ParaView no puede colorearla ni construir una tabla de color, y al pedirle
    la representacion "Volume" responde

        Failed to determine the LookupTable being used.
        Could not determine array range.

    que suena a archivo corrupto y solo significa que no hay nada que mapear.
    Se incluyen por tanto la calidad y el volumen de cada elemento, que ademas
    son lo que de verdad interesa mirar: permiten localizar en pantalla los
    elementos degenerados antes de mandar la malla a un solver.

    La calidad se evalua sobre la celda LINEAL —el hexaedro, o el tetraedro de
    las cuatro esquinas en el caso del TET10— porque el filtro de calidad de
    VTK no soporta celdas cuadraticas y devuelve -1 para todas, lo que aparenta
    una malla catastrofica cuando en realidad es que no sabe medirla. Como los
    TET10 son de lados rectos, la calidad de forma es la del tetraedro de
    esquinas y la medida es exacta, no una aproximacion.
    """
    ruta = Path(ruta)
    nodos = np.asarray(nodos, float)
    elems = np.asarray(elems, dtype=np.int64)
    n = elems.shape[1]
    celdas = np.hstack([np.full((elems.shape[0], 1), n, dtype=np.int64),
                        elems]).ravel()
    grid = pv.UnstructuredGrid(celdas,
                               np.full(elems.shape[0], TIPO_VTK[n]), nodos)

    if arrays:
        # Celda lineal equivalente para las medidas que VTK no sabe hacer
        # sobre celdas cuadraticas
        if n == 10:
            lin, tipo_lin, nl = elems[:, :4], pv.CellType.TETRA, 4
        else:
            lin, tipo_lin, nl = elems, pv.CellType.HEXAHEDRON, 8
        cl = np.hstack([np.full((lin.shape[0], 1), nl, dtype=np.int64),
                        lin]).ravel()
        g_lin = pv.UnstructuredGrid(cl, np.full(lin.shape[0], tipo_lin), nodos)

        try:
            if hasattr(g_lin, "cell_quality"):
                q = g_lin.cell_quality("scaled_jacobian")
            else:
                q = g_lin.compute_cell_quality(quality_measure="scaled_jacobian")
            clave = ("scaled_jacobian" if "scaled_jacobian" in q.cell_data
                     else "CellQuality")
            grid.cell_data["calidad"] = np.abs(
                np.asarray(q.cell_data[clave], dtype=float))
        except Exception:
            pass

        try:
            vol = g_lin.compute_cell_sizes(length=False, area=False, volume=True)
            grid.cell_data["volumen"] = np.abs(
                np.asarray(vol.cell_data["Volume"], dtype=float))
        except Exception:
            pass

        grid.cell_data["id_elemento"] = np.arange(elems.shape[0], dtype=np.int32)
        grid.point_data["z"] = nodos[:, 2].astype(np.float32)

    grid.save(str(ruta))
    return {"ruta": str(ruta), "MB": ruta.stat().st_size / 1e6,
            "volumen": float(grid.volume),
            "arrays": list(grid.cell_data.keys()) + list(grid.point_data.keys())}


def escribir_abaqus(nodos, elems, ruta, E_s=20e9, nu_s=0.30, nombre="SOLIDO",
                    conjuntos=True):
    """Malla volumetrica en formato Abaqus .inp.

    ANSYS Mechanical lo importa directamente (External Model), asi que el mismo
    archivo sirve para las dos plataformas. Se escriben los conjuntos de nodos
    BASE y TECHO para que aplicar el apoyo y la carga en el preprocesador sea
    un clic y no una seleccion manual sobre cientos de miles de nodos.
    """
    ruta = Path(ruta)
    nodos = np.asarray(nodos, float)
    elems = np.asarray(elems, dtype=np.int64)
    tipo = TIPO_ABAQUS.get(elems.shape[1])
    if tipo is None:
        raise ValueError(f"No se sabe escribir elementos de {elems.shape[1]} nodos")
    E_MPa = float(E_s) / 1e6

    with open(ruta, "w", encoding="utf-8") as fh:
        w = fh.write
        w(f"** {nombre} — malla {tipo} generada por spinpy/escribe.py\n")
        w("** Sistema de unidades coherente: mm, N, MPa\n")
        w(f"** E = {E_MPa:g} MPa ({E_s:g} Pa), nu = {nu_s:g}\n")
        w("*HEADING\n")
        w(f"{nombre}\n")
        w("*NODE\n")
        for i, p in enumerate(nodos, start=1):
            w(f"{i}, {p[0]:.9g}, {p[1]:.9g}, {p[2]:.9g}\n")
        w(f"*ELEMENT, TYPE={tipo}, ELSET=SOLIDO\n")
        for e, c in enumerate(elems, start=1):
            w(f"{e}, " + ", ".join(str(int(n) + 1) for n in c) + "\n")
        if conjuntos:
            base, techo = _caras_z(nodos)
            w("*NSET, NSET=BASE\n"); _lista(w, base + 1)
            w("*NSET, NSET=TECHO\n"); _lista(w, techo + 1)
        w("*SOLID SECTION, ELSET=SOLIDO, MATERIAL=MATRIZ\n")
        w("*MATERIAL, NAME=MATRIZ\n*ELASTIC\n")
        w(f"{E_MPa:g}, {nu_s:g}\n")

    return {"ruta": str(ruta), "MB": ruta.stat().st_size / 1e6, "tipo": tipo,
            "n_nodos": int(nodos.shape[0]), "n_elems": int(elems.shape[0]),
            "E_MPa": E_MPa, "unidades": "mm-N-MPa"}


def escribir_apdl(nodos, elems, ruta, E_s=20e9, nu_s=0.30):
    """Script de Mechanical APDL con NBLOCK/EBLOCK y conjuntos de nodos."""
    ruta = Path(ruta)
    nodos = np.asarray(nodos, float)
    elems = np.asarray(elems, dtype=np.int64)
    nn = elems.shape[1]
    tipo = TIPO_ANSYS.get(nn)
    if tipo is None:
        raise ValueError(f"No se sabe escribir elementos de {nn} nodos")
    E_MPa = float(E_s) / 1e6
    base, techo = _caras_z(nodos)

    with open(ruta, "w", encoding="utf-8") as fh:
        w = fh.write
        w("! Malla generada por spinpy/escribe.py\n")
        w("! Sistema de unidades coherente: mm, N, MPa\n")
        w("/PREP7\n")
        w(f"ET,1,{tipo}\n")
        w(f"MP,EX,1,{E_MPa:g}\nMP,PRXY,1,{nu_s:g}\n")
        w(f"NBLOCK,3,,{nodos.shape[0]}\n(1i9,3e20.9e3)\n")
        for i, p in enumerate(nodos, start=1):
            w(f"{i:9d}{p[0]:20.9e}{p[1]:20.9e}{p[2]:20.9e}\n")
        w("N,R5.3,LOC,-1,\n")
        w(f"EBLOCK,19,SOLID,,{elems.shape[0]}\n(19i9)\n")
        for e, c in enumerate(elems, start=1):
            cab = [1, 1, 1, 1, 0, 0, 0, 0, nn, 0, e]
            w("".join(f"{v:9d}" for v in cab)
              + "".join(f"{int(x)+1:9d}" for x in c) + "\n")
        w("-1\n")
        w("! Conjuntos de nodos para apoyo y carga\n")
        for ids, nom in ((base, "BASE"), (techo, "TECHO")):
            w("NSEL,NONE\n")
            for v in ids:
                w(f"NSEL,A,NODE,,{int(v)+1}\n")
            w(f"CM,{nom},NODE\n")
        w("ALLSEL,ALL\nFINISH\n")

    return {"ruta": str(ruta), "MB": ruta.stat().st_size / 1e6, "tipo": tipo,
            "n_base": int(base.size), "n_techo": int(techo.size),
            "E_MPa": E_MPa, "unidades": "mm-N-MPa"}
