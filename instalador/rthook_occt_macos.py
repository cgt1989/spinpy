"""Precarga de OpenCASCADE en el .app de macOS con procesador Intel.

POR QUE HACE FALTA. En la rueda de `netgen-mesher` para macOS, la parte
x86_64 de `libngpy.so` NO enlaza ninguna biblioteca OpenCASCADE (cero entradas
`libTK*` en sus comandos de carga) y se apoya en que alguien las haya cargado
antes con alcance global; la parte arm64 si las enlaza por `@rpath`. Quien las
carga es `netgen/__init__.py` (`load_occ_libs`), con `ctypes.CDLL(...,
RTLD_GLOBAL)`, pero solo si existe el metadato del paquete `netgen-occt`.

La receta del .app elimina ese metadato a proposito (sus rutas apuntan a
`<entorno>/lib`, que en el .app no existe), asi que en Intel la precarga no
ocurre y `import netgen.libngpy` falla con «symbol not found in flat
namespace '__ZN16Standard_Failure14SetStackStringEPKc'». Medido en la
construccion de CI de la V2.3.0 (el arm64 pasa; el x86_64 no).

QUE HACE. Hace lo mismo que `load_occ_libs`, con las bibliotecas que la receta
copia a la raiz de `Frameworks`, y en el mismo orden. Solo actua en macOS
x86_64: en arm64 el enlace por `@rpath` ya funciona y no se toca.
"""
import ctypes
import glob
import os
import platform
import sys

# Orden de `netgen.load_occ_libs` (invertido: primero las de las que dependen
# las demas).
_ORDEN = ["TKernel", "TKDE", "TKMath", "TKG2d", "TKG3d", "TKGeomBase",
          "TKGeomAlgo", "TKTopAlgo", "TKShHealing", "TKService", "TKMesh",
          "TKHLR", "TKV3d", "TKCDF", "TKLCAF", "TKPrim", "TKBO", "TKCAF",
          "TKVCAF", "TKXCAF", "TKDESTL", "TKXSBase", "TKDESTEP", "TKBool",
          "TKDEIGES", "TKFillet", "TKOffset"]


def _precargar():
    if sys.platform != "darwin" or platform.machine() != "x86_64":
        return
    bases = [getattr(sys, "_MEIPASS", None),
             os.path.join(os.path.dirname(sys.executable), "..", "Frameworks")]
    base = next((b for b in bases if b and glob.glob(
        os.path.join(b, "libTKernel.*dylib"))), None)
    if base is None:
        return
    for nombre in _ORDEN:
        for ruta in sorted(glob.glob(os.path.join(base, f"lib{nombre}.*dylib"))):
            try:
                ctypes.CDLL(ruta, mode=ctypes.RTLD_GLOBAL)
            except OSError as e:
                sys.stderr.write(f"[spinpy] no se pudo precargar {ruta}: {e}\n")
            break


_precargar()
