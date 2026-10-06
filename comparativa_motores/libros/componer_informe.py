"""
componer_informe.py: Rellena el informe con las fuentes de `fuentes_libros.py`.

    python comparativa_motores/libros/componer_informe.py

Lee informe_libros/informe_libros.html, sustituye los marcadores
<!--FUENTES:Cn--> (bloque «De donde sale» de cada ficha) y <!--MATRIZ-->
(matriz de trazabilidad del anexo) y escribe
informe_libros/informe_libros_compuesto.html, que es el que se imprime:

    cd comparativa_motores/libros/informe_libros
    ENTRADA=informe_libros_compuesto.html SALIDA=INFORME_LIBROS.pdf sh imprimir.sh
"""

from __future__ import annotations

import html
import re
import sys
from pathlib import Path

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
from fuentes_libros import FUENTES, LIBROS, MEDIDAS               # noqa: E402

CARPETA = AQUI / "informe_libros"


def autor_corto(clave):
    ref = LIBROS[clave][0]
    autores, resto = ref.split(". ", 1)
    año = re.search(r"\b(?:19|20)\d{2}\b", ref).group(0)
    nombres = [a.strip().split(" ")[0] for a in autores.split(",")]
    if len(nombres) == 1:
        return f"{nombres[0]} ({año})"
    if len(nombres) == 2:
        return f"{nombres[0]} y {nombres[1]} ({año})"
    return f"{nombres[0]} et al. ({año})"


def bloque(cambio):
    filas = [f for f in FUENTES if f["cambio"] == cambio]
    e = html.escape
    partes = ['<div class="fuentes">', '<h4>De dónde sale</h4>']
    for f in filas:
        libro, archivo = LIBROS[f["libro"]]
        partes.append(
            '<div class="fuente">'
            f'<div class="fcab"><span class="codigo">{f["id"]}</span> '
            f'<b>{e(autor_corto(f["libro"]))}</b>, {e(f["ubicacion"])}</div>'
            f'<div class="floc">Archivo: <code>{e(archivo)}</code>, página '
            f'{e(f["archivo_pag"])} · {e(f["parrafo"])}</div>'
            f'<blockquote>{e(f["cita"])}</blockquote>'
            f'<div class="fresp"><b>Respalda:</b> {e(f["respalda"])}</div>'
            '</div>')
    if cambio in MEDIDAS:
        partes.append(f'<div class="fmed"><b>Lo que sale de las medidas, no '
                      f'de los libros:</b> {e(MEDIDAS[cambio])}</div>')
    partes.append("</div>")
    return "\n".join(partes)


def matriz():
    e = html.escape
    filas = ['<table class="matriz"><thead><tr><th>Fuente</th><th>Cambio</th>'
             '<th>Libro y ubicación</th><th>Archivo, página</th>'
             '<th>Respalda</th></tr></thead><tbody>']
    for f in FUENTES:
        libro, archivo = LIBROS[f["libro"]]
        filas.append(
            f'<tr><td><span class="codigo">{f["id"]}</span></td>'
            f'<td>{f["cambio"]}</td>'
            f'<td>{e(autor_corto(f["libro"]))}, {e(f["ubicacion"])}</td>'
            f'<td><code>{e(archivo)}</code>, p. {e(f["archivo_pag"])}</td>'
            f'<td>{e(f["respalda"])}</td></tr>')
    filas.append("</tbody></table>")
    return "\n".join(filas)


def main():
    t = (CARPETA / "informe_libros.html").read_text(encoding="utf-8")
    for c in sorted({f["cambio"] for f in FUENTES}):
        marca = f"<!--FUENTES:{c}-->"
        assert t.count(marca) == 1, marca
        t = t.replace(marca, bloque(c))
    assert t.count("<!--MATRIZ-->") == 1
    t = t.replace("<!--MATRIZ-->", matriz())
    (CARPETA / "informe_libros_compuesto.html").write_text(t, encoding="utf-8")
    print("escrito informe_libros_compuesto.html")


if __name__ == "__main__":
    main()
