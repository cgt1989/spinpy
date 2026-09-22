"""generar_manual.py — Compone el manual de referencia a partir del codigo.

    python docs/generar_manual.py

POR QUE SE GENERA Y NO SE ESCRIBE
----------------------------------
En este proyecto los docstrings NO son un resumen del codigo: son la
documentacion real. Llevan las decisiones de diseno, los sesgos medidos, las
trampas que costaron encontrar y las advertencias sobre que numeros no son
citables. Un manual escrito aparte seria una segunda version de esa
informacion, y las dos versiones se separan en cuanto alguien toca el codigo
sin acordarse de la otra.

Extrayendolos, el manual no puede quedar desactualizado sin que se note: si una
funcion cambia y su docstring no, el problema esta en el docstring, que es donde
tiene que estar.

Se usa `ast` y no `import`: importar `visor.py` exigiria PyQt5 y abriria una
ventana, y el manual tiene que poder componerse en una maquina sin interfaz.
"""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
# La copia vive DENTRO del repositorio a proposito: si apuntara a la de
# Validacion_Anexo, quien clone solo este repositorio no podria regenerar el
# manual. La unica diferencia con aquella es que el encabezado de pagina sale
# del titulo del documento en vez de estar fijo.
MD2PDF = Path(__file__).resolve().parent / "md2pdf.py"

# Orden deliberado: se lee de lo que sostiene todo lo demas hacia lo que se ve.
MODULOS = [
    ("spinpy/__init__.py", "El paquete"),
    ("spinpy/grf.py", "Generacion: el campo aleatorio gaussiano"),
    ("spinpy/io.py", "Lectura de volumenes de interes"),
    ("spinpy/morphometry.py", "Morfometria: BV/TV, superficie, MIL, Conn.D, SMI"),
    ("spinpy/espesor.py", "Espesor local por esferas inscritas"),
    ("spinpy/error.py", "Funcion de error entre VOI y candidato"),
    ("spinpy/fit.py", "Ajuste: busqueda escalonada y desempate mecanico"),
    ("spinpy/elastic.py", "Homogeneizacion elastica periodica"),
    ("spinpy/resistencia.py", "Ensayo de compresion, von Mises y Pistoia"),
    ("spinpy/solido.py", "Mallas solidas: hexaedrica y TET10"),
    ("spinpy/escribe.py", "Exportadores: Abaqus, APDL, VTU, STL"),
    ("spinpy/lote.py", "Lote de VOIs y dispersion del generador"),
    ("spinpy/estadistica.py", "Varianza, ICC, N efectivo y equivalencia TOST"),
    ("spinpy/cli.py", "Linea de comandos"),
    ("visor.py", "Interfaz grafica"),
]


def firma(nodo):
    """Firma legible de una funcion, con sus valores por defecto."""
    a = nodo.args
    partes = []
    pos = a.posonlyargs + a.args
    defs = [None] * (len(pos) - len(a.defaults)) + list(a.defaults)
    for arg, d in zip(pos, defs):
        s = arg.arg
        if d is not None:
            try:
                s += "=" + ast.unparse(d)
            except Exception:
                s += "=..."
        partes.append(s)
    if a.vararg:
        partes.append("*" + a.vararg.arg)
    for arg, d in zip(a.kwonlyargs, a.kw_defaults):
        s = arg.arg
        if d is not None:
            try:
                s += "=" + ast.unparse(d)
            except Exception:
                s += "=..."
        partes.append(s)
    if a.kwarg:
        partes.append("**" + a.kwarg.arg)
    return "%s(%s)" % (nodo.name, ", ".join(partes))


def limpiar(doc, sangria="", max_lineas=None):
    """Normaliza un docstring para que md2pdf lo componga bien.

    Dos cuidados. Las lineas de guiones que subrayan los titulos internos se
    convierten en regla horizontal, que es lo que uno querria de todos modos.
    Y las lineas que empiezan por `|` se escapan: md2pdf las tomaria por una
    tabla y fallaria al no encontrar la fila separadora.
    """
    if not doc:
        return ""
    lineas = [l.rstrip() for l in doc.strip("\n").split("\n")]
    if max_lineas:
        lineas = lineas[:max_lineas]
    salida = []
    for l in lineas:
        t = l.strip()
        if t.startswith("|") and not t.endswith("|"):
            l = l.replace("|", "/")
        salida.append(sangria + l)
    # md2pdf trata la sangria como texto corrido; se quita la comun para que
    # los parrafos no salgan con una indentacion aleatoria.
    return "\n".join(salida).strip() + "\n"


def publica(nombre):
    return not nombre.startswith("_") or nombre in ("__init__",)


def documentar(ruta, titulo):
    arbol = ast.parse(ruta.read_text(encoding="utf-8"))
    out = ["## %s" % titulo, "", "`%s`" % ruta.name, ""]

    doc_mod = ast.get_docstring(arbol)
    if doc_mod:
        out += [limpiar(doc_mod), ""]

    clases = [n for n in arbol.body if isinstance(n, ast.ClassDef)]
    funcs = [n for n in arbol.body if isinstance(n, ast.FunctionDef)
             and publica(n.name)]
    # Las funciones anidadas de visor.py son metodos: se recogen por clase.
    n_priv = len([n for n in arbol.body
                  if isinstance(n, ast.FunctionDef) and not publica(n.name)])

    if funcs:
        out += ["### Funciones", ""]
        for f in funcs:
            out += ["#### `%s`" % firma(f), ""]
            d = ast.get_docstring(f)
            out += [limpiar(d) if d else "_Sin documentar._", ""]

    for c in clases:
        out += ["### Clase `%s`" % c.name, ""]
        d = ast.get_docstring(c)
        if d:
            out += [limpiar(d), ""]
        metodos = [n for n in c.body if isinstance(n, ast.FunctionDef)
                   and publica(n.name) and ast.get_docstring(n)]
        for m in metodos:
            out += ["#### `%s.%s`" % (c.name, firma(m)), ""]
            out += [limpiar(ast.get_docstring(m)), ""]
        sin_doc = [n.name for n in c.body if isinstance(n, ast.FunctionDef)
                   and publica(n.name) and not ast.get_docstring(n)]
        if sin_doc:
            out += ["Metodos sin docstring propio (su comportamiento se "
                    "explica en la clase): " +
                    ", ".join("`%s`" % s for s in sin_doc), ""]

    if n_priv:
        out += ["", "_%d funcion(es) interna(s) de este modulo no se listan; "
                "su nombre empieza por guion bajo y no forman parte de la "
                "interfaz publica._" % n_priv, ""]
    return "\n".join(out)


def main():
    partes = [
        "# Manual de referencia de spinpy",
        "",
        "Ajuste de microestructuras espinodales a volumenes de hueso "
        "trabecular obtenidos por microtomografia.",
        "",
        "> Este manual esta **generado del codigo**. Cada apartado reproduce "
        "la documentacion que acompana a la funcion en el fuente, y no una "
        "reescritura: en este proyecto los docstrings llevan las decisiones "
        "de diseno, los sesgos medidos y las advertencias sobre que numeros "
        "no son citables. Si algo aqui esta desactualizado, lo esta tambien "
        "en el codigo, que es donde hay que corregirlo.",
        "",
        "---",
        "",
    ]
    n_mod = 0
    for rel, titulo in MODULOS:
        ruta = RAIZ / rel
        if not ruta.exists():
            print("  falta %s" % rel)
            continue
        partes.append(documentar(ruta, titulo))
        partes.append("")
        n_mod += 1

    md = RAIZ / "docs" / "MANUAL.md"
    md.write_text("\n".join(partes), encoding="utf-8")
    print("%d modulos -> %s (%d KB)" % (n_mod, md, md.stat().st_size / 1024))

    if MD2PDF.exists():
        pdf = RAIZ / "docs" / "MANUAL_spinpy.pdf"
        r = subprocess.run([sys.executable, str(MD2PDF), str(md), str(pdf)],
                           capture_output=True, text=True)
        if r.returncode == 0 and pdf.exists():
            print("-> %s (%.1f MB)" % (pdf, pdf.stat().st_size / 1e6))
        else:
            print("md2pdf fallo:", (r.stderr or r.stdout)[-1500:])
    else:
        print("md2pdf no encontrado en", MD2PDF)
    return 0


if __name__ == "__main__":
    sys.exit(main())
