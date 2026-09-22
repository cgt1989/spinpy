# -*- coding: utf-8 -*-
"""Resuelve las DLL del interprete de Anaconda que PyInstaller no sigue.

EL PROBLEMA
-----------
Este proyecto se desarrolla sobre Anaconda. Sus modulos de extension de la
biblioteca estandar -`_ctypes.pyd`, `_ssl.pyd`, `_sqlite3.pyd`...- no enlazan
contra las DLL del sistema, sino contra DLL propias de conda que viven en
`anaconda3\\Library\\bin`, una carpeta que no esta en la ruta de busqueda que
PyInstaller inspecciona. El resultado es un ejecutable que se construye sin un
solo aviso y muere al arrancar con

    ImportError: DLL load failed while importing _ctypes

que no dice cual es la DLL que falta. Medido aqui: `_ctypes.pyd` necesita
`ffi.dll` y el propio `python313.dll` necesita `zlib.dll`, ninguna de las dos
recogida por el analisis.

LO QUE HACE ESTE MODULO
-----------------------
Recorre las importaciones de cada .pyd y .dll ya recogidos, resuelve las que
falten contra las carpetas de conda, y repite hasta cerrar el grafo -las DLL de
conda dependen unas de otras-. Devuelve la lista en el formato de `binaries`
que espera PyInstaller.

Se ignoran a proposito las DLL del sistema (`KERNEL32`, `api-ms-win-*`,
`VCRUNTIME140`...): esas las pone Windows, y copiarlas dentro del paquete es
justamente lo que rompe una instalacion en otra maquina.

LA ALTERNATIVA, POR SI ALGUN DIA CONVIENE
-----------------------------------------
Todo esto desaparece construyendo desde un CPython de python.org en vez de
desde Anaconda. No se hizo asi para no tener que instalar un segundo Python en
la maquina de desarrollo; si el empaquetado se muda a un servidor de
integracion continua, alli lo correcto es el interprete oficial y este modulo
sobra.
"""
from __future__ import annotations

import sys
from pathlib import Path

try:
    import pefile
except ImportError:                                    # pragma: no cover
    pefile = None

# Prefijos de DLL que pone el sistema operativo y que NUNCA hay que copiar.
SISTEMA = (
    "api-ms-win-", "ext-ms-", "kernel32", "user32", "advapi32", "ole32",
    "oleaut32", "shell32", "shlwapi", "ws2_32", "wsock32", "gdi32", "comdlg32",
    "comctl32", "version", "winmm", "wintrust", "crypt32", "bcrypt", "ncrypt",
    "secur32", "setupapi", "iphlpapi", "dbghelp", "imm32", "rpcrt4", "msvcrt",
    "vcruntime", "ucrtbase", "msvcp", "concrt", "d3d", "dxgi", "opengl32",
    "glu32", "netapi32", "userenv", "powrprof", "propsys", "dwmapi",
    "uxtheme", "mfplat", "mf.dll", "authz", "cfgmgr32", "normaliz",
    "python3.dll",
)


def _es_del_sistema(nombre: str) -> bool:
    n = nombre.lower()
    return any(n.startswith(p) for p in SISTEMA)


def _importa(ruta: Path) -> set[str]:
    """Nombres de las DLL que importa un binario. Vacio si no se puede leer."""
    if pefile is None:
        return set()
    try:
        pe = pefile.PE(str(ruta), fast_load=True)
        pe.parse_data_directories(
            [pefile.DIRECTORY_ENTRY["IMAGE_DIRECTORY_ENTRY_IMPORT"]])
        out = {e.dll.decode("ascii", "ignore")
               for e in getattr(pe, "DIRECTORY_ENTRY_IMPORT", [])}
        pe.close()
        return out
    except Exception:
        return set()


def carpetas_conda(prefijo: Path | None = None) -> list[Path]:
    """Donde conda guarda sus DLL, en el orden en que Windows las buscaria."""
    base = Path(prefijo) if prefijo else Path(sys.base_prefix)
    cand = [base / "Library" / "bin", base / "DLLs", base,
            base / "Library" / "mingw-w64" / "bin",
            base / "Library" / "usr" / "bin"]
    return [c for c in cand if c.is_dir()]


def _cerrar_grafo(semillas: list[Path], presentes: set[str],
                  fuentes: list[Path]) -> dict[str, Path]:
    """Recorre importaciones hasta que no falte ninguna DLL de conda.

    Es un cierre transitivo y no una sola pasada porque las DLL de conda
    dependen unas de otras: `_ctypes.pyd` pide `ffi.dll`, y `libpq` arrastra
    `krb5`, `gssapi` y `comerr`. Con una sola pasada el ejecutable falla en el
    segundo nivel, que ademas es mas dificil de diagnosticar que el primero.
    """
    encontradas: dict[str, Path] = {}
    pendientes = list(semillas)
    vistos: set[str] = set()
    while pendientes:
        b = pendientes.pop()
        if b.name.lower() in vistos:
            continue
        vistos.add(b.name.lower())
        for dep in _importa(b):
            d = dep.lower()
            if d in presentes or d in encontradas or _es_del_sistema(d):
                continue
            for f in fuentes:
                cand = f / dep
                if cand.is_file():
                    encontradas[d] = cand
                    pendientes.append(cand)
                    break
    return encontradas


def resolver_binarios(rutas, prefijo: Path | None = None):
    """DLL de conda que faltan, partiendo de los binarios ya recogidos.

    `rutas` son las rutas ORIGEN de lo que PyInstaller ya ha decidido incluir
    (`[src for _, src, _ in a.binaries]`). Devuelve entradas listas para
    sumarse a `a.binaries`: (nombre_destino, ruta_origen, "BINARY").

    Se usa asi -y no con una lista escrita a mano- porque la lista cambia con
    la version de Anaconda y con las bibliotecas instaladas. Una lista fija
    caduca en silencio: el ejecutable se construye igual y falla al arrancar
    en la maquina de otro.
    """
    fuentes = carpetas_conda(prefijo)
    if not fuentes:
        return []
    semillas = [Path(r) for r in rutas
                if Path(r).suffix.lower() in (".dll", ".pyd")]
    presentes = {p.name.lower() for p in semillas}
    hallazgos = _cerrar_grafo(semillas, presentes, fuentes)
    return [(p.name, str(p), "BINARY") for p in sorted(hallazgos.values())]


def resolver(raices: list[Path], prefijo: Path | None = None):
    """Igual, pero mirando dentro de una carpeta ya construida.

    Sirve para diagnosticar un `dist/` que falla sin tener que reconstruirlo.
    """
    fuentes = carpetas_conda(prefijo)
    if not fuentes:
        return []
    semillas = [p for r in raices for p in Path(r).rglob("*")
                if p.suffix.lower() in (".dll", ".pyd")]
    presentes = {p.name.lower() for p in semillas}
    hallazgos = _cerrar_grafo(semillas, presentes, fuentes)
    return [(str(p), ".") for p in sorted(hallazgos.values())]


if __name__ == "__main__":
    # Uso:  python dlls_conda.py <carpeta_del_ejecutable>
    raiz = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    faltan = resolver([raiz])
    for origen, _ in faltan:
        print(origen)
    print(f"\n{len(faltan)} DLL de conda que faltaban.")
