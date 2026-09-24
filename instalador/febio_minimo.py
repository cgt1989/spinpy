"""febio_minimo.py — El conjunto MINIMO de archivos para correr febio4.exe.

    python instalador/febio_minimo.py                     # lista y prueba
    python instalador/febio_minimo.py --copiar DESTINO    # ademas lo copia

La carpeta bin de FEBio Studio pesa ~350 MB, casi todo interfaz (Qt, codecs
de video, OpenGL). Este guion lee la tabla de IMPORTACIONES de cada ejecutable
PE (sin dependencias: `pefile` no esta instalado) y sigue el cierre transitivo
desde febio4.exe y los plugins que declara febio.xml, quedandose con lo que
esta en esa misma carpeta (las DLL del sistema no se copian). Despues copia el
conjunto a una carpeta AISLADA y corre alli un modelo pequeno: si falta algo,
FEBio no arranca y el guion lo dice.

LICENCIA. Los binarios de FEBio Studio estan bajo la FEBio Software License
4.0 (doc/FEBio_EULA_4.pdf), que no permite redistribuirlos a terceros sin una
licencia aparte de la Universidad de Utah. Se incluyen en el Instalador y el
RAR por decision del responsable del laboratorio, con la licencia de
redistribucion EN TRAMITE; ver `febio/LICENCIA_FEBIO.txt` del paquete.
"""

from __future__ import annotations

import argparse
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

ORIGEN = Path(r"C:\Program Files\FEBioStudio2\bin")


def importaciones(ruta):
    """Nombres de DLL importadas (tabla normal y de carga diferida)."""
    b = Path(ruta).read_bytes()
    pe = struct.unpack_from("<I", b, 0x3C)[0]
    if b[pe:pe + 4] != b"PE\0\0":
        return set()
    n_sec = struct.unpack_from("<H", b, pe + 6)[0]
    t_opt = struct.unpack_from("<H", b, pe + 20)[0]
    opt = pe + 24
    magic = struct.unpack_from("<H", b, opt)[0]
    dd = opt + (112 if magic == 0x20B else 96)
    secs = []
    for i in range(n_sec):
        s = opt + t_opt + 40 * i
        vs, va, rs, rp = struct.unpack_from("<IIII", b, s + 8)
        secs.append((va, max(vs, rs), rp))

    def off(rva):
        for va, tam, rp in secs:
            if va <= rva < va + tam:
                return rva - va + rp
        return None

    def cadena(rva):
        o = off(rva)
        if o is None:
            return None
        return b[o:b.index(b"\0", o)].decode("ascii", "replace")

    out = set()
    for entrada, tam in ((1, 20), (13, 32)):          # import, delay import
        rva = struct.unpack_from("<I", b, dd + 8 * entrada)[0]
        o = off(rva) if rva else None
        while o is not None:
            campos = struct.unpack_from("<5I" if tam == 20 else "<8I", b, o)
            nombre_rva = campos[3] if tam == 20 else campos[1]
            if not any(campos):
                break
            n = cadena(nombre_rva)
            if n:
                out.add(n.lower())
            o += tam
    return out


def conjunto_minimo(origen=ORIGEN):
    presentes = {p.name.lower(): p for p in origen.iterdir() if p.is_file()}
    semillas = ["febio4.exe"]
    xml = origen / "febio.xml"
    if xml.exists():
        semillas += [m.lower() for m in
                     re.findall(r"<import>\s*([^<]+?)\s*</import>",
                                xml.read_text(errors="replace"))]
    hechos, cola = set(), list(semillas)
    while cola:
        n = cola.pop()
        if n in hechos or n not in presentes:
            continue
        hechos.add(n)
        cola += [d for d in importaciones(presentes[n]) if d in presentes]
    archivos = [presentes[n] for n in sorted(hechos)]
    if xml.exists():
        archivos.append(xml)
    return archivos


def probar(carpeta):
    """Corre un hexaedro en la carpeta aislada; devuelve la version o lanza."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from spinpy import febio
    exe = Path(carpeta) / "febio4.exe"
    febio._VERSIONES.pop(str(exe), None)
    v = febio.version(exe)
    if not v:
        raise RuntimeError(f"febio4.exe no arranca desde {carpeta}")
    return v


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--copiar", default=None)
    ap.add_argument("--origen", default=str(ORIGEN))
    a = ap.parse_args()
    archivos = conjunto_minimo(Path(a.origen))
    total = sum(p.stat().st_size for p in archivos)
    for p in archivos:
        print(f"  {p.name:28s} {p.stat().st_size / 2**20:7.2f} MB")
    print(f"  {len(archivos)} archivos, {total / 2**20:.1f} MB")
    destino = Path(a.copiar) if a.copiar else Path(tempfile.mkdtemp(
        prefix="febio_aislado_"))
    destino.mkdir(parents=True, exist_ok=True)
    for p in archivos:
        shutil.copy2(p, destino / p.name)
    print(f"copiado a {destino}; probando…")
    print("FEBio", probar(destino), "corre desde la copia aislada")


if __name__ == "__main__":
    main()
