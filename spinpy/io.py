"""
io.py — Lectura de VOIs en VTK legacy ASCII (STRUCTURED_POINTS).

Port de `readVTKVOI` (AppFinal_V2.m:5147-5205) y de la binarizacion que aplica
`computeVOIMetrics` (AppFinal_V2.m:5229).

CONVENCION DE EJES — la invariante que hay que respetar
-------------------------------------------------------
El formato VTK escribe el indice X como el que varia mas rapido, asi que
`reshape(data, [nx ny nz])` de MATLAB (orden por columnas) deja la
DIMENSION 1 = X. Es lo que documenta la correccion del bug de ejes X/Y
intercambiados en `computeVOIMetrics`: la version anterior asignaba la
dimension 1 a Y y con voxeles anisotropos —lo normal en micro-CT— deformaba la
geometria medida y devolvia una direccion principal erronea, que es la que fija
la rotacion aplicada al spinodoide.

En numpy el equivalente exacto de ese reshape es `order='F'`. Usar el orden C
por defecto transpone el volumen en silencio.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np


def _validar_spacing(spacing, nombre):
    """El espaciado de vóxel tiene que ser finito y estrictamente positivo.

    Un valor nulo, negativo o no finito no produce un error aguas abajo: da
    volúmenes negativos, NaN o una escala equivocada en toda la morfometría y
    la rigidez. Se corta aquí, en la puerta de entrada.
    """
    sp = np.asarray(spacing, dtype=float).ravel()
    if sp.size != 3 or not np.all(np.isfinite(sp)) or np.any(sp <= 0):
        raise ValueError(f"{nombre}: espaciado de voxel invalido {sp.tolist()}; "
                         f"debe ser finito y mayor que 0 en los tres ejes.")
    return sp


def _sin_espaciado(nombre):
    warnings.warn(
        f"{nombre}: el archivo no declara el espaciado del voxel; se usa 1.0 "
        f"en los tres ejes, de modo que TODAS las magnitudes quedan en unidades "
        f"de voxel y no en mm. Indica la escala real antes de interpretar "
        f"morfometria o rigidez.", UserWarning, stacklevel=3)


def leer_vtk_voi(ruta, binarizar=True, umbral=0):
    """Lee un VOI de un VTK legacy ASCII STRUCTURED_POINTS.

    Devuelve (VOI, spacing) con VOI de forma (nx, ny, nz), dimension 1 = X.

    binarizar : aplica `VOI > umbral`, igual que `computeVOIMetrics`, que hace
        `if ~islogical(VOI), VOI = VOI > 0; end`. Los VOIs de H4 vienen como
        unsigned_char con valores 0 y 255.
    """
    ruta = Path(ruta)
    texto = ruta.read_text(encoding="latin-1")

    def _busca(patron, n):
        m = re.search(patron, texto, re.IGNORECASE)
        if not m:
            return None
        try:
            return [float(m.group(i + 1)) for i in range(n)]
        except ValueError:
            raise ValueError(f"{ruta.name}: valores no numericos en "
                             f"{m.group(0).split()[0]}") from None

    dims = _busca(r"DIMENSIONS\s+(\d+)\s+(\d+)\s+(\d+)", 3)
    if dims is None:
        raise ValueError(f"No se encontro DIMENSIONS en {ruta.name}")
    nx, ny, nz = (int(v) for v in dims)
    if min(nx, ny, nz) < 1:
        raise ValueError(f"{ruta.name}: DIMENSIONS {nx} {ny} {nz} no es un "
                         f"volumen valido.")

    # \S+ y no una clase numerica: "nan" o "inf" tienen que LEERSE para poder
    # rechazarlos, no saltarse en silencio y caer al valor por defecto.
    esp = _busca(r"SPACING\s+(\S+)\s+(\S+)\s+(\S+)", 3)
    if esp is None:
        # Formatos antiguos usan ASPECT_RATIO en lugar de SPACING
        esp = _busca(r"ASPECT_RATIO\s+(\S+)\s+(\S+)\s+(\S+)", 3)
    if esp:
        spacing = _validar_spacing(esp, ruta.name)
    else:
        _sin_espaciado(ruta.name)
        spacing = np.ones(3)

    m = re.search(r"POINT_DATA\s+(\d+)", texto, re.IGNORECASE)
    if not m:
        raise ValueError(f"No se encontro POINT_DATA en {ruta.name}")
    npts = int(m.group(1))

    # Tras POINT_DATA vienen SCALARS y LOOKUP_TABLE, y luego los valores
    resto = texto[m.end():]
    m2 = re.search(r"LOOKUP_TABLE\s+\S+\s*\n", resto, re.IGNORECASE)
    if not m2:
        raise ValueError(f"No se encontro LOOKUP_TABLE en {ruta.name}")
    datos = np.fromstring(resto[m2.end():], sep=" ", dtype=np.float64)

    if datos.size < npts:
        raise ValueError(
            f"{ruta.name} declara {npts} valores ({nx}x{ny}x{nz}) pero contiene "
            f"{datos.size} ({100*datos.size/npts:.0f}%). Esta incompleto o "
            f"truncado; vuelve a exportar el VOI.")
    datos = datos[:npts]
    if not np.all(np.isfinite(datos)):
        raise ValueError(f"{ruta.name}: contiene valores NaN o infinitos; "
                         f"se rechaza en vez de tratarlos como vacio.")
    if npts != nx * ny * nz:
        raise ValueError(f"{ruta.name}: POINT_DATA={npts} no cuadra con "
                         f"{nx}x{ny}x{nz}={nx*ny*nz}")

    # order='F' = reshape por columnas de MATLAB -> dimension 1 = X
    VOI = datos.reshape((nx, ny, nz), order="F")
    if binarizar:
        VOI = VOI > umbral
    return VOI, spacing


def leer_mat_voi(ruta, binarizar=True, umbral=0):
    """Lee un VOI guardado como .mat de MATLAB (variables `VOI` y `spacing`).

    Es el formato del subconjunto porcino publico (Koria, Mengoni & Brockett
    2020, Research Data Leeds, CC BY 4.0). El array ya viene con dimension 1 =
    X porque MATLAB lo guardo asi; scipy.io.loadmat lo devuelve con la misma
    disposicion, sin necesidad de transponer.
    """
    from scipy.io import loadmat

    d = loadmat(str(ruta))
    if "VOI" not in d:
        claves = [k for k in d if not k.startswith("__")]
        raise ValueError(f"{Path(ruta).name}: no contiene la variable 'VOI' "
                         f"(tiene {claves})")
    VOI = np.asarray(d["VOI"])
    if VOI.ndim != 3:
        raise ValueError(f"{Path(ruta).name}: 'VOI' no es un volumen 3D "
                         f"(forma {VOI.shape})")

    if "spacing" in d:
        spacing = np.asarray(d["spacing"], dtype=float).ravel()[:3]
        if spacing.size == 1:
            spacing = np.repeat(spacing, 3)
        spacing = _validar_spacing(spacing, Path(ruta).name)
    else:
        _sin_espaciado(Path(ruta).name)
        spacing = np.ones(3)

    if binarizar:
        VOI = VOI > umbral
    return VOI, spacing


def leer_voi(ruta, **kw):
    """Lee un VOI eligiendo el lector por lo que sea `ruta`.

    .vtk / .mat  volumen ya extraido, con su spacing dentro del archivo.
    carpeta      pila de rebanadas TIFF del micro-CT.
    .tif/.tiff   TIFF multipagina.

    Los dos ultimos casos van a `voi.leer_pila_tiff`, que ademas de apilar
    resuelve la ESCALA: la saca del `*_rec.log` de SkyScan o de un `tam_voxel`
    explicito. Un .vtk trae el spacing dentro; una pila de imagenes no trae
    nada, y suponer 1 mm/voxel escalaria en silencio toda la morfometria.

    Devuelve siempre `(BW, spacing)` en mm. `leer_pila_tiff` devuelve ademas un
    tercer valor con la trazabilidad (umbral aplicado, de donde salio la
    escala); si se necesita, hay que llamarla directamente.
    """
    p = Path(ruta)
    if p.is_dir():
        BW, spacing, _info = _leer_pila(p, **kw)
        return BW, spacing
    ext = p.suffix.lower()
    if ext == ".vtk":
        return leer_vtk_voi(ruta, **kw)
    if ext == ".mat":
        return leer_mat_voi(ruta, **kw)
    if ext in (".tif", ".tiff"):
        BW, spacing, _info = _leer_pila(p, **kw)
        return BW, spacing
    raise ValueError(f"Extension no soportada: {ext} "
                     f"(use .vtk, .mat, .tif/.tiff o una carpeta de rebanadas)")


def _leer_pila(ruta, **kw):
    # Import diferido: `voi` arrastra scikit-image y Pillow, y quien solo lee
    # un .vtk no tiene por que pagarlos.
    from .voi import leer_pila_tiff
    return leer_pila_tiff(ruta, **kw)
