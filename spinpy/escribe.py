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

from .resistencia import normalizar_apoyo

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


def escribir_febio(nodos, elems, ruta, E_s=20e9, nu_s=0.30, sigma_app=1e6,
                   A_bruta=None, apoyo="deslizante", datos=True,
                   eps_plato=None):
    """Ensayo de compresion en z listo para correr en FEBio 4 (`febio4 -i`).

    No es solo la malla: reproduce el MISMO problema que resuelve
    `resistencia.ensayo_compresion`, para que el resultado de FEBio pueda
    enfrentarse numero a numero con el de spinpy.

      material   'isotropic elastic' de FEBio (St. Venant-Kirchhoff): se reduce
                 a la elasticidad lineal de spinpy para deformaciones pequenas.
                 La diferencia es de orden de la deformacion del tejido; con
                 cargas pequenas es despreciable y, como el problema es lineal,
                 se puede cargar poco y reescalar.
      carga      presion NO seguidora (`linear` = 1) sobre las caras superiores
                 de los elementos que tocan el techo. La fuerza total es
                 sigma_app * A_bruta, repartida sobre el area osea del techo:
                 misma convencion de seccion BRUTA que spinpy. FEBio integra la
                 presion por su cuenta, asi que tambien pone a prueba el reparto
                 por area tributaria de spinpy (un cuarto por nodo y cara).
      apoyo      'deslizante': uz = 0 en la base, (ux, uy) en la esquina
                 (min x, min y) y uy en la (max x, min y), las mismas esquinas
                 geometricas que elige spinpy. 'empotrado': base fija entera.
      solver     un paso, Newton completo (max_ups = 0), Pardiso directo:
                 independiente del AMG + CG de spinpy.

    `A_bruta` es la seccion del VOI en mm^2. Si no se da, se toma la caja de
    los nodos, que COINCIDE con la del VOI salvo que el filtro de portantes
    haya vaciado una columna del borde: en ese caso hay que pasarla.

    Con `datos` se piden al logfile los desplazamientos nodales
    (`<ruta>_u.txt`) y las tensiones de Cauchy por elemento en orden Voigt del
    proyecto, xx yy zz yz xz xy (`<ruta>_s.txt`). FEBio da la media de los
    ocho puntos de Gauss, que en un hexaedro rectangular es exactamente el
    valor en el centro que usa spinpy. Las reacciones NO se piden: FEBio 4.5
    las escribe como cero en los GDL de un `zero displacement` (comprobado
    tambien con la variable del plotfile). El equilibrio se comprueba con la
    integral de volumen: sum(sigma_zz * V_e) = -F * H, exacta en el problema
    discreto.

    `eps_plato` cambia el control de carga: en lugar de la presion, un PLATO
    RIGIDO sin friccion impone el mismo uz = -eps_plato * H a todos los nodos
    del techo (ux, uy libres). No es el ensayo de la app —que controla la
    fuerza— sino su contrapunto: con fuerza impuesta, una trabecula cortada por
    la cara del VOI recibe carga en su extremo libre y trabaja en voladizo; con
    el plato, ese extremo se mueve con el resto del techo. `sigma_app` se
    ignora en ese caso; la fuerza se mide despues, por la misma integral de
    volumen de la tension.

    Unidades mm - N - MPa, como el resto de escritores: E y sigma se reciben en
    Pa y se escriben en MPa.
    """
    ruta = Path(ruta)
    nodos = np.asarray(nodos, float)
    elems = np.asarray(elems, dtype=np.int64)
    if elems.shape[1] != 8:
        raise ValueError("escribir_febio solo escribe hexaedros de 8 nodos; "
                         f"se recibieron elementos de {elems.shape[1]}")
    E_MPa = float(E_s) / 1e6
    base, techo = _caras_z(nodos)

    z = nodos[:, 2]
    tolz = 1e-9 * max(float(z.max() - z.min()), 1.0)
    caras = elems[np.all(z[elems[:, 4:8]] >= z.max() - tolz, axis=1), 4:8]
    if caras.shape[0] == 0:
        raise ValueError("Ningun elemento llega al techo: no hay donde cargar.")
    p0, p1, p3 = nodos[caras[:, 0]], nodos[caras[:, 1]], nodos[caras[:, 3]]
    A_osea = float(np.linalg.norm(np.cross(p1 - p0, p3 - p0), axis=1).sum())
    if A_bruta is None:
        A_bruta = float(np.ptp(nodos[:, 0]) * np.ptp(nodos[:, 1]))
    F_N = float(sigma_app) / 1e6 * float(A_bruta)
    p_MPa = F_N / A_osea

    cb = nodos[base]
    ancla_xy = int(base[int(np.argmin(cb[:, 0] + cb[:, 1]))])
    ancla_y = int(base[int(np.argmax(cb[:, 0] - cb[:, 1]))])
    empotrado = normalizar_apoyo(apoyo) == "empotrado"
    stem = ruta.stem

    with open(ruta, "w", encoding="utf-8") as fh:
        w = fh.write
        w('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
        w("<!-- Ensayo de compresion en z generado por spinpy/escribe.py.\n")
        w("     Unidades mm-N-MPa. E = %g MPa, nu = %g, sigma_app = %g MPa,\n"
          % (E_MPa, nu_s, float(sigma_app) / 1e6))
        w("     A_bruta = %.9g mm2, F = %.9g N, apoyo %s. -->\n"
          % (A_bruta, F_N, "empotrado" if empotrado else "deslizante"))
        w('<febio_spec version="4.0">\n')
        w('\t<Module type="solid">\n\t\t<units>mm-N-s</units>\n\t</Module>\n')
        w("\t<Control>\n\t\t<analysis>STATIC</analysis>\n")
        w("\t\t<time_steps>1</time_steps>\n\t\t<step_size>1</step_size>\n")
        w("\t\t<solver>\n\t\t\t<max_refs>50</max_refs>\n")
        w('\t\t\t<qn_method type="BFGS">\n\t\t\t\t<max_ups>0</max_ups>\n'
          "\t\t\t</qn_method>\n\t\t\t<dtol>1e-9</dtol>\n")
        w("\t\t\t<etol>1e-12</etol>\n\t\t\t<rtol>1e-12</rtol>\n")
        w("\t\t\t<lstol>0.9</lstol>\n")
        w('\t\t\t<linear_solver type="pardiso"/>\n\t\t</solver>\n')
        w("\t</Control>\n")
        w('\t<Material>\n\t\t<material id="1" name="hueso" '
          'type="isotropic elastic">\n')
        w(f"\t\t\t<E>{E_MPa:.9g}</E>\n\t\t\t<v>{nu_s:.9g}</v>\n")
        w("\t\t</material>\n\t</Material>\n")

        w('\t<Mesh>\n\t\t<Nodes name="todos">\n')
        w("".join(f'\t\t\t<node id="{i}">{p[0]:.12g},{p[1]:.12g},{p[2]:.12g}'
                  "</node>\n" for i, p in enumerate(nodos, start=1)))
        w('\t\t</Nodes>\n\t\t<Elements type="hex8" name="solido">\n')
        e1 = elems + 1
        w("".join(f'\t\t\t<elem id="{i}">' + ",".join(map(str, c))
                  + "</elem>\n" for i, c in enumerate(e1.tolist(), start=1)))
        w("\t\t</Elements>\n")
        for nom, ids in (("base", base), ("techo", techo),
                         ("ancla_xy", [ancla_xy]), ("ancla_y", [ancla_y])):
            w(f'\t\t<NodeSet name="{nom}">\n')
            ids = np.asarray(ids, dtype=np.int64) + 1
            for i in range(0, ids.size, 16):
                w("\t\t\t" + ", ".join(map(str, ids[i:i + 16].tolist()))
                  + ("," if i + 16 < ids.size else "") + "\n")
            w("\t\t</NodeSet>\n")
        w('\t\t<Surface name="techo_carga">\n')
        w("".join(f'\t\t\t<quad4 id="{i}">' + ",".join(map(str, c))
                  + "</quad4>\n"
                  for i, c in enumerate((caras + 1).tolist(), start=1)))
        w("\t\t</Surface>\n\t</Mesh>\n")
        w('\t<MeshDomains>\n\t\t<SolidDomain name="solido" mat="hueso"/>\n'
          "\t</MeshDomains>\n")

        w("\t<Boundary>\n")
        bcs = ([("base", 1, 1, 1)] if empotrado else
               [("base", 0, 0, 1), ("ancla_xy", 1, 1, 0), ("ancla_y", 0, 1, 0)])
        for ns, bx, by, bz in bcs:
            w(f'\t\t<bc name="fijo_{ns}" type="zero displacement" '
              f'node_set="{ns}">\n')
            w(f"\t\t\t<x_dof>{bx}</x_dof>\n\t\t\t<y_dof>{by}</y_dof>\n"
              f"\t\t\t<z_dof>{bz}</z_dof>\n\t\t</bc>\n")
        if eps_plato is not None:
            uz = -float(eps_plato) * float(z.max() - z.min())
            w('\t\t<bc name="plato" type="prescribed displacement" '
              'node_set="techo">\n')
            w(f'\t\t\t<dof>z</dof>\n\t\t\t<value lc="1">{uz:.17g}</value>\n'
              "\t\t\t<relative>0</relative>\n\t\t</bc>\n")
        w("\t</Boundary>\n")
        if eps_plato is None:
            w('\t<Loads>\n\t\t<surface_load name="compresion" '
              'type="pressure" surface="techo_carga">\n')
            w(f'\t\t\t<pressure lc="1">{p_MPa:.17g}</pressure>\n')
            w("\t\t\t<linear>1</linear>\n")
            w("\t\t\t<symmetric_stiffness>1</symmetric_stiffness>\n")
            w("\t\t</surface_load>\n\t</Loads>\n")
        w('\t<LoadData>\n\t\t<load_controller id="1" name="rampa" '
          'type="loadcurve">\n')
        w("\t\t\t<interpolate>LINEAR</interpolate>\n\t\t\t<points>\n")
        w("\t\t\t\t<point>0,0</point>\n\t\t\t\t<point>1,1</point>\n")
        w("\t\t\t</points>\n\t\t</load_controller>\n\t</LoadData>\n")
        w("\t<Output>\n")
        w('\t\t<plotfile type="febio">\n'
          '\t\t\t<var type="displacement"/>\n\t\t\t<var type="stress"/>\n'
          "\t\t</plotfile>\n")
        if datos:
            w("\t\t<logfile>\n")
            w(f'\t\t\t<node_data data="ux;uy;uz" delim=" " '
              f'file="{stem}_u.txt"/>\n')
            w(f'\t\t\t<element_data data="sx;sy;sz;syz;sxz;sxy" delim=" " '
              f'file="{stem}_s.txt"/>\n')
            w("\t\t</logfile>\n")
        w("\t</Output>\n</febio_spec>\n")

    return {"ruta": str(ruta), "MB": ruta.stat().st_size / 1e6,
            "tipo": "hex8", "n_nodos": int(nodos.shape[0]),
            "n_elems": int(elems.shape[0]), "n_caras_techo": int(caras.shape[0]),
            "A_bruta_mm2": float(A_bruta), "A_osea_techo_mm2": A_osea,
            "F_N": F_N, "presion_MPa": p_MPa, "E_MPa": E_MPa,
            "control": "fuerza" if eps_plato is None else "plato",
            "eps_plato": eps_plato,
            "apoyo": "empotrado" if empotrado else "deslizante",
            "ancla_xy": ancla_xy, "ancla_y": ancla_y,
            "unidades": "mm-N-MPa"}


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
