# -*- mode: python ; coding: utf-8 -*-
"""Receta de PyInstaller para el ejecutable de Windows.

MODO CARPETA (--onedir), Y NO ARCHIVO UNICO
-------------------------------------------
`--onefile` produce un solo .exe, que suena mejor, pero descomprime el
contenido entero -del orden de 700 MB- en el directorio temporal EN CADA
ARRANQUE: entre 20 y 40 s mirando una pantalla vacia antes de ver la ventana,
todas las veces. En modo carpeta el arranque baja a unos pocos segundos, y el
usuario final no ve la carpeta de todos modos porque el instalador la deja en
su sitio y lo que el toca es un acceso directo.

Hay ademas una razon que no es de comodidad: en modo archivo unico la carpeta
temporal se BORRA al cerrar, y cualquier ruta derivada de `__file__` apunta
alli. Eso ya esta resuelto en `visor.py` (`DATOS` frente a `RAIZ`), pero el
modo carpeta quita el filo al problema en vez de depender solo de ese arreglo.

QUE SE EXCLUYE Y POR QUE
------------------------
La lista de `excludes` no es cosmetica: sin ella entran tkinter, IPython,
Jupyter y los demas enlaces de Qt, que no se usan y que anaden centenares de
megas. Dos de esas exclusiones importan de verdad:

  * PySide2/PySide6/PyQt6 - si alguno viaja dentro junto a PyQt5, Qt puede
    cargar los complementos de plataforma del enlace equivocado y la
    aplicacion muere con "could not find or load the Qt platform plugin".
  * tkinter - matplotlib lo detecta y arrastra todo Tcl/Tk para un backend
    que aqui no se usa nunca; el unico backend necesario es Qt5Agg.

LICENCIA DEL BINARIO
--------------------
Este ejecutable incluye PyQt5, que es GPL v3. El codigo fuente de spinpy es
MIT, pero la obra combinada que aqui se distribuye queda bajo GPL-3.0. Ver
`LICENCIA_BINARIO.txt`.
"""

import os
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

sys.path.insert(0, SPECPATH)          # noqa: F821
from dlls_conda import resolver_binarios      # noqa: E402

RAIZ = Path(SPECPATH).parent          # ...\Port_Python   # noqa: F821

# Con SPINPY_CONSOLA=1 el ejecutable se construye CON consola. Sin ella, un
# fallo de arranque solo ensena un cuadro que dice "Unhandled exception in
# script" y se traga la traza, que es justo lo que hace falta. Se deja como
# interruptor y no como una version aparte porque el mismo problema puede
# aparecer en la maquina de otro, donde no hay Python con que reproducirlo.
CONSOLA = bool(os.environ.get("SPINPY_CONSOLA"))
ICONO = Path(SPECPATH) / "spinpy.ico"  # noqa: F821

# vtkmodules se carga de forma perezosa: los modulos concretos no aparecen en
# ningun `import` que el analizador pueda ver, asi que hay que enumerarlos.
ocultos = collect_submodules("vtkmodules")
ocultos += [
    "pyvista",
    "pyvistaqt",
    "pyamg",
    "pyamg.relaxation",
    "pyamg.krylov",
    "scipy._lib.array_api_compat.numpy.fft",
    "skimage.measure",
    "skimage.morphology",
    "tetgen",
    "pymeshfix",
    "matplotlib.backends.backend_qt5agg",
]

# El icono viaja como DATO ademas de ir incrustado en el .exe: el que se
# incrusta lo usa el explorador para el archivo, pero la ventana y la barra de
# tareas lo piden en tiempo de ejecucion via `QIcon`, y para eso tiene que
# existir como archivo dentro del paquete.
datos = [(str(ICONO), ".")] if ICONO.exists() else []

# La carpeta `Test` viaja como DATO, no como codigo: sus guiones se importan
# por ruta en tiempo de ejecucion (ver dialogo_validacion.py) y sus figuras de
# referencia son mapas de bits. Sin esto, el menu Validacion del ejecutable no encuentra
# nada y la unica prueba de que el metodo es el publicado se queda fuera.
_test = RAIZ / "Test"
if _test.is_dir():
    for _f in _test.rglob("*"):
        if _f.is_file() and "__pycache__" not in _f.parts                 and _f.parent.name != "resultados":
            datos.append((str(_f),
                          str(Path("Test") / _f.relative_to(_test).parent)))

# MOTORES FEM INTERNOS (spinpy.motores). Desde V2.0.0 el ejecutable no lleva
# FEBio: los ensayos con malla suave y los no lineales los resuelven motores
# que corren DENTRO del proceso de spinpy. Viajan los que tienen rueda binaria
# para Windows en PyPI (comprobado 2026-09-30 con `pip download --platform
# win_amd64`): NGSolve (LGPL-2.1, el recomendado) y scikit-fem (BSD-3, Python
# puro). SfePy solo publica ruedas para Linux y FEniCSx no esta en PyPI
# (conda-forge, y compila C en tiempo de ejecucion): no viajan; si el usuario
# los instala en su Python, la app los detecta. NGSolve no necesita MKL: su
# Cholesky propio fue mas rapido que PARDISO en todas las medidas.
#
# NGSolve carga sus bibliotecas (libngsolve, netgen) desde su paquete en
# tiempo de ejecucion: PyInstaller no las ve en ningun import.
from PyInstaller.utils.hooks import collect_dynamic_libs  # noqa: E402
ocultos += collect_submodules("ngsolve") + collect_submodules("netgen")
ocultos += collect_submodules("skfem")
ocultos += ["spinpy.motores.m_app", "spinpy.motores.m_ngsolve",
            "spinpy.motores.m_skfem", "spinpy.motores.m_fenicsx",
            "spinpy.motores.m_sfepy"]
binarios = collect_dynamic_libs("ngsolve") + collect_dynamic_libs("netgen")
# Las carpetas `netgen_mesher.libs` y `ngsolve.libs` (delvewheel) van tal
# cual, con su msvcp140.dll 14.50: `spinpy/_msvc.py` la busca ahi y la carga
# antes que la 14.26 de PyQt5. Sin ella netgen muere al importarse (0xC0000005).
import importlib.util as _ilu  # noqa: E402
for _mod, _libs in (("netgen", "netgen_mesher.libs"), ("ngsolve", "ngsolve.libs")):
    _spec = _ilu.find_spec(_mod)
    if _spec is None or not _spec.submodule_search_locations:
        continue
    _dir = Path(list(_spec.submodule_search_locations)[0]).parent / _libs
    if _dir.is_dir():
        binarios += [(str(_f), _libs) for _f in _dir.glob("*.dll")]
datos += collect_data_files("ngsolve") + collect_data_files("netgen")

datos += collect_data_files("pyvista")
datos += collect_data_files("vtkmodules")
datos += collect_data_files("matplotlib", subdir="mpl-data")
datos += collect_data_files("skimage", includes=["**/*.pyi"])
# El paquete propio va como codigo, no como dato: `spinpy` es Python puro y
# PyInstaller lo sigue solo desde los imports de visor.py.

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

# --- DLL del interprete de Anaconda -----------------------------------------
# Sin esto el ejecutable se construye sin un aviso y muere al arrancar con
# "ImportError: DLL load failed while importing _ctypes", porque los .pyd de
# la biblioteca estandar de Anaconda enlazan contra DLL que viven en la
# carpeta Library/bin de la distribucion, que el analisis no mira. Las resuelve
# dlls_conda.py recorriendo la tabla de importaciones de cada binario.
_extra = resolver_binarios([src for _, src, _ in a.binaries])
if _extra:
    print(f"[spinpy] DLL de conda anadidas ({len(_extra)}): "
          + ", ".join(n for n, _, _ in _extra))
    a.binaries += _extra

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
    upx=False,          # UPX y las DLL de Qt/VTK se llevan mal; no se usa
    console=CONSOLA,    # ver SPINPY_CONSOLA arriba
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=str(ICONO) if ICONO.exists() else None,
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
