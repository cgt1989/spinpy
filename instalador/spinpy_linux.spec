# -*- mode: python ; coding: utf-8 -*-
"""Receta de PyInstaller para spinpy en Linux (x86_64).

Es `spinpy_mac.spec` sin el paquete .app (BUNDLE) y con bibliotecas `.so` en
lugar de `.dylib`. Modo carpeta, por la misma razon que en Windows: arranque
rapido y sin descomprimir nada en cada ejecucion. Salida: `dist/spinpy/` con el
ejecutable `spinpy` y sus bibliotecas, que `construir_linux.sh` empaqueta en
`.tar.gz`.

VERSION DE PRUEBA: se construye y supera la autocomprobacion en un runner de
GitHub Actions (Ubuntu), pero la interfaz grafica no se ha probado a mano en un
escritorio Linux. Ver `LEEME_LINUX.txt`.

LICENCIA DEL BINARIO
--------------------
Incluye PyQt5 (GPL v3): la obra combinada se distribuye bajo GPL-3.0, igual
que los ejecutables de Windows y macOS. Ver `LICENCIA_BINARIO.txt`.
"""

import os
import sys
from pathlib import Path

from PyInstaller.utils.hooks import (collect_data_files, collect_dynamic_libs,
                                     collect_submodules)

RAIZ = Path(SPECPATH).parent          # noqa: F821
VERSION = os.environ.get("SPINPY_VERSION", "0.0.0")

ocultos = collect_submodules("vtkmodules")
ocultos += [
    "pyvista", "pyvistaqt", "pyamg", "pyamg.relaxation", "pyamg.krylov",
    "scipy._lib.array_api_compat.numpy.fft",
    "skimage.measure", "skimage.morphology",
    "tetgen", "pymeshfix",
    "matplotlib.backends.backend_qt5agg",
]

datos = [(str(Path(SPECPATH) / "spinpy_icono.png"), ".")]  # noqa: F821

_test = RAIZ / "Test"
if _test.is_dir():
    for _f in _test.rglob("*"):
        if _f.is_file() and "__pycache__" not in _f.parts \
                and _f.parent.name != "resultados":
            datos.append((str(_f),
                          str(Path("Test") / _f.relative_to(_test).parent)))

# Motores FEM internos: los mismos que viajan en Windows (NGSolve y
# scikit-fem). NGSolve carga sus bibliotecas desde su paquete en tiempo de
# ejecucion y PyInstaller no las ve en ningun import.
ocultos += collect_submodules("ngsolve") + collect_submodules("netgen")
ocultos += collect_submodules("skfem")
ocultos += ["spinpy.motores.m_app", "spinpy.motores.m_ngsolve",
            "spinpy.motores.m_skfem", "spinpy.motores.m_fenicsx",
            "spinpy.motores.m_sfepy"]
binarios = collect_dynamic_libs("ngsolve") + collect_dynamic_libs("netgen")

# Las ruedas manylinux llevan sus dependencias en `<paquete>.libs` junto al
# paquete (auditwheel): se copian tal cual a la misma ruta relativa, que es
# donde las busca el enlazador por rpath.
import importlib.util as _ilu  # noqa: E402
for _mod in ("netgen", "ngsolve", "netgen_occt"):
    _spec = _ilu.find_spec(_mod)
    if _spec is None or not _spec.submodule_search_locations:
        continue
    _pkg = Path(list(_spec.submodule_search_locations)[0])
    for _dir in (_pkg.parent / f"{_mod}.libs", _pkg / ".libs"):
        if _dir.is_dir():
            _rel = _dir.relative_to(_pkg.parent)
            binarios += [(str(_f), str(_rel)) for _f in _dir.glob("*.so*")]
            print(f"[spinpy] bibliotecas de {_mod}: {_dir}")

# OpenCASCADE (paquete netgen-occt): `netgen/libngpy.so` enlaza con sus
# bibliotecas (`libTKOffset.so.7.8.1`...) y el enlazador las busca en la raiz
# del paquete (misma disposicion que en macOS, medida en CI). En un entorno
# normal viven en `<entorno>/lib`, fuera de site-packages, asi que PyInstaller
# no las ve: se recogen de los archivos que declara la distribucion y van a la
# RAIZ.
# Sus metadatos NO viajan (filtro tras `Analysis`): sin ellos
# `netgen/__init__.py` no intenta precargarlas por ruta absoluta.
from importlib import metadata as _md  # noqa: E402
_occt = []
try:
    for _f in _md.files("netgen-occt") or []:
        if (_f.name.startswith("libTK") and ".so" in _f.name):
            _p = Path(_f.locate()).resolve()
            if _p.is_file():
                _occt.append((str(_p), "."))
except _md.PackageNotFoundError:
    pass
if sys.platform.startswith("linux") and not _occt:
    raise SystemExit("[spinpy] No se encontraron las bibliotecas de netgen-occt "
                     "(libTK*.so): NGSolve no funcionaria en el ejecutable.")
print(f"[spinpy] bibliotecas de OpenCASCADE para netgen: {len(_occt)}")
binarios += list(dict.fromkeys(_occt))
datos += collect_data_files("ngsolve") + collect_data_files("netgen")

datos += collect_data_files("pyvista")
datos += collect_data_files("vtkmodules")
datos += collect_data_files("matplotlib", subdir="mpl-data")
datos += collect_data_files("skimage", includes=["**/*.pyi"])

a = Analysis(                                        # noqa: F821
    [str(RAIZ / "visor.py")],
    pathex=[str(RAIZ)],
    binaries=binarios,
    datas=datos,
    hiddenimports=ocultos,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "tkinter", "_tkinter", "Tkinter",
        "PySide2", "PySide6", "PyQt6", "shiboken2", "shiboken6",
        "IPython", "jupyter", "notebook", "ipykernel", "jedi",
        "pytest", "sphinx", "docutils", "setuptools", "pip",
        "matplotlib.backends.backend_webagg",
        "matplotlib.backends.backend_tkagg",
    ],
    noarchive=False,
    optimize=0,
)

# Sin los metadatos de netgen-occt (ver arriba, junto a las bibliotecas OCC).
a.datas = [d for d in a.datas
           if not d[0].replace("\\", "/").lower().startswith("netgen_occt-")]

pyz = PYZ(a.pure)                                    # noqa: F821

exe = EXE(                                           # noqa: F821
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="spinpy",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
)

coll = COLLECT(                                      # noqa: F821
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="spinpy",
)
