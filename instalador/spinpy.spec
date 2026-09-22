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

datos += collect_data_files("pyvista")
datos += collect_data_files("vtkmodules")
datos += collect_data_files("matplotlib", subdir="mpl-data")
datos += collect_data_files("skimage", includes=["**/*.pyi"])
# El paquete propio va como codigo, no como dato: `spinpy` es Python puro y
# PyInstaller lo sigue solo desde los imports de visor.py.

a = Analysis(                                        # noqa: F821
    [str(RAIZ / "visor.py")],
    pathex=[str(RAIZ)],
    binaries=[],
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
