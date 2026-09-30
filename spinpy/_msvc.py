"""
_msvc.py: Precarga en Windows la version mas reciente del runtime de C++
(msvcp140.dll) antes de que la cargue PyQt5.

EL PROBLEMA (medido en la compilacion de la V2.0.1)
---------------------------------------------------
PyQt5-Qt5 5.15.2 trae su propio msvcp140.dll, version 14.26 (2020). netgen y
NGSolve 6.2.2607 traen el suyo, 14.50, con el MISMO nombre (delvewheel no lo
renombra; VTK y numpy si). Windows resuelve una DLL por su nombre con la que
ya este cargada en el proceso: la interfaz importa PyQt5 al arrancar, asi que
cuando el motor NGSolve importa netgen, `pyngcore` se enlaza con el runtime
14.26 y muere al inicializarse con una violacion de acceso (codigo
0xC0000005). Las bibliotecas compiladas con Visual Studio 17.10 o posterior no
funcionan con un msvcp140 anterior; al reves si: el runtime nuevo es
compatible con lo compilado para el viejo.

LA SOLUCION
-----------
Cargar primero, por su ruta completa, el msvcp140.dll mas reciente que haya a
mano (el de netgen/NGSolve, el del sistema): las cargas posteriores por nombre,
las de Qt incluidas, reutilizan ese. Se llama desde `spinpy/__init__.py` y
desde lo primero de `visor.py`; los procesos hijo de los motores (spawn)
vuelven a importar `visor.py` y pasan por aqui antes de tocar PyQt5.

Fuera de Windows no hace nada.
"""

from __future__ import annotations

import os
import struct
import sys
from pathlib import Path

#: Lo que se cargo (ruta, version) o None; para la autocomprobacion.
CARGADO = None

#: DLL del runtime de Visual C++ que no deben mezclarse de versiones distintas.
RUNTIME = ("msvcp140.dll", "msvcp140_1.dll", "msvcp140_2.dll",
           "msvcp140_atomic_wait.dll", "msvcp140_codecvt_ids.dll",
           "vcruntime140.dll", "vcruntime140_1.dll", "concrt140.dll")


def version_dll(ruta):
    """Version de archivo (a, b, c, d) de una DLL, leida de su recurso
    VS_FIXEDFILEINFO; None si no se encuentra."""
    try:
        datos = Path(ruta).read_bytes()
    except OSError:
        return None
    i = datos.find("VS_VERSION_INFO".encode("utf-16-le"))
    if i < 0:
        return None
    j = datos.find(b"\xbd\x04\xef\xfe", i)
    if j < 0 or j + 16 > len(datos):
        return None
    ms, ls = struct.unpack_from("<II", datos, j + 8)
    return (ms >> 16, ms & 0xFFFF, ls >> 16, ls & 0xFFFF)


def _bases():
    """Directorios bajo los que buscar las carpetas `*.libs` de las ruedas."""
    bases = []
    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        bases.append(Path(meipass))
    try:
        import importlib.util
        for mod in ("netgen", "ngsolve"):
            spec = importlib.util.find_spec(mod)
            if spec is not None and spec.submodule_search_locations:
                for d in spec.submodule_search_locations:
                    bases.append(Path(d).parent)
    except Exception:
        pass
    return bases


def candidatos():
    """Rutas de msvcp140.dll que podrian cargarse, sin repetir."""
    rutas = []
    for b in _bases():
        for sub in ("netgen_mesher.libs", "ngsolve.libs", "."):
            rutas.append(b / sub / "msvcp140.dll")
    sistema = os.environ.get("SystemRoot")
    if sistema:
        rutas.append(Path(sistema) / "System32" / "msvcp140.dll")
    vistas, out = set(), []
    for r in rutas:
        try:
            k = str(r.resolve()).lower()
        except OSError:
            continue
        if k not in vistas and r.is_file():
            vistas.add(k)
            out.append(r)
    return out


def elegir(rutas):
    """(ruta, version) de la version mas alta; None si no hay ninguna."""
    mejor = None
    for r in rutas:
        v = version_dll(r)
        if v is not None and (mejor is None or v > mejor[1]):
            mejor = (r, v)
    return mejor


def precargar():
    """Carga el msvcp140.dll mas reciente disponible (solo Windows).

    No lanza nunca: si algo falla, la aplicacion sigue como antes.
    """
    global CARGADO
    if os.name != "nt" or CARGADO is not None:
        return CARGADO
    try:
        mejor = elegir(candidatos())
        if mejor is None:
            return None
        import ctypes
        ctypes.WinDLL(str(mejor[0]))
        CARGADO = (str(mejor[0]), ".".join(str(x) for x in mejor[1]))
    except Exception as e:                           # noqa: BLE001
        CARGADO = ("no se pudo precargar", f"{type(e).__name__}: {e}")
    return CARGADO


def cargados():
    """[(nombre, ruta, version)] de las DLL del runtime ya cargadas en este
    proceso (solo Windows); para diagnosticar mezclas de versiones."""
    if os.name != "nt":
        return []
    out = []
    try:
        import ctypes
        from ctypes import wintypes
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k32.GetModuleHandleW.restype = ctypes.c_void_p
        k32.GetModuleHandleW.argtypes = [wintypes.LPCWSTR]
        k32.GetModuleFileNameW.restype = wintypes.DWORD
        k32.GetModuleFileNameW.argtypes = [ctypes.c_void_p, wintypes.LPWSTR,
                                           wintypes.DWORD]
        for nombre in RUNTIME:
            h = k32.GetModuleHandleW(nombre)
            if not h:
                continue
            buf = ctypes.create_unicode_buffer(1024)
            k32.GetModuleFileNameW(h, buf, 1024)
            v = version_dll(buf.value)
            out.append((nombre, buf.value,
                        ".".join(map(str, v)) if v else "?"))
    except Exception:
        pass
    return out
