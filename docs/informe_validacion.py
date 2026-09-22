"""
informe_validacion.py — Del registro JSON de la suite al informe Markdown.

Lee `resultados/validacion_literatura.json` (lo escribe `tests/conftest.py`
al terminar pytest) y produce `docs/validacion_literatura/INFORME.md`, con una
tabla por bloque y una lista aparte de lo que NO paso, porque eso es lo que
hay que leer primero. El PDF sale de `docs/md2pdf.py`.

Uso:
    python -m pytest tests            (o `-m "not lento"` para lo rapido)
    python docs/informe_validacion.py
    python docs/md2pdf.py docs/validacion_literatura/INFORME.md salida.pdf
"""
from __future__ import annotations

import json
from collections import OrderedDict
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
ENTRADA = RAIZ / "resultados" / "validacion_literatura.json"
SALIDA = RAIZ / "docs" / "validacion_literatura" / "INFORME.md"


def _num(v):
    if v is None:
        return "—"
    a = abs(v)
    if a == 0:
        return "0"
    if a < 1e-3 or a >= 1e5:
        return "%.2e" % v
    return ("%.6g" % v).replace(".", ",")


def _pct(v):
    return "—" if v is None else ("%+.2f %%" % (100 * v)).replace(".", ",")


def _celda(s):
    """Texto seguro dentro de una tabla Markdown: la barra vertical es el
    separador de columnas, y las tolerancias se escriben como |Th - t|."""
    return str(s).replace("|", "¦")


def main():
    d = json.loads(ENTRADA.read_text(encoding="utf-8"))
    filas = d["comprobaciones"]
    bloques = OrderedDict()
    for f in filas:
        bloques.setdefault(f["bloque"], []).append(f)
    fallos = [f for f in filas if not f["ok"]]

    L = []
    L.append("# Verificación de spinpy contra literatura y soluciones cerradas")
    L.append("")
    L.append("**Invariantes topológicos, definiciones publicadas, soluciones "
             "exactas de la elasticidad y cotas universales — con las "
             "tolerancias declaradas antes de medir**")
    L.append("")
    L.append("## Qué se verifica y con qué criterio")
    L.append("")
    L.append("Cada comprobación compara una función de `spinpy` con una "
             "respuesta que se conoce **de antemano** y que no depende del "
             "propio código: un número de Euler, el valor que una definición "
             "publicada asigna a una figura ideal, la solución cerrada de un "
             "problema elástico, o una desigualdad que toda microestructura "
             "cumple. La tolerancia de cada prueba está escrita en el módulo "
             "de pruebas, con su justificación, **antes** de ejecutarla; no se "
             "ajusta después para que pase.")
    L.append("")
    L.append("Una comprobación que no pasa no es un fallo de la suite: es un "
             "hallazgo sobre el código o sobre su documentación, y se lista "
             "primero.")
    L.append("")
    L.append("> **En palabras sencillas.** No comparamos el programa con otro "
             "programa, sino con cosas que se saben seguras: una rosquilla "
             "tiene exactamente un agujero, un bloque macizo se comprime "
             "exactamente como dice la ley de Hooke, y ningún material poroso "
             "puede ser más rígido que el sólido del que está hecho. Si el "
             "programa acierta esas cosas, sus números merecen confianza.")
    L.append("")
    L.append("Resumen: **%d de %d comprobaciones dentro de tolerancia**. "
             "Ejecutado el %s con Python %s, numpy %s, scipy %s, scikit-image "
             "%s." % (d["n_ok"], d["n_comprobaciones"], d["fecha"],
                      d["versiones"].get("python", "?"),
                      d["versiones"].get("numpy", "?"),
                      d["versiones"].get("scipy", "?"),
                      d["versiones"].get("scikit-image", "?")))
    L.append("")

    L.append("## Lo que NO pasó")
    L.append("")
    if not fallos:
        L.append("Ninguna comprobación quedó fuera de tolerancia.")
    else:
        L.append("| Bloque | Prueba | Esperado | Obtenido | Error | Tolerancia |")
        L.append("|---|---|---|---|---|---|")
        for f in fallos:
            L.append("| %s | %s | %s | %s | %s | %s |" % (
                _celda(f["bloque"]), _celda(f["prueba"]), _num(f["esperado"]),
                _num(f["obtenido"]), _pct(f["error_rel"]),
                _celda(f["tolerancia"])))
        L.append("")
        L.append("Cada uno de estos casos se discute en su bloque.")
    L.append("")

    for b, fs in bloques.items():
        n_ok = sum(f["ok"] for f in fs)
        L.append("## %s" % b)
        L.append("")
        refs = sorted(set(f["referencia"] for f in fs))
        L.append("Referencias: " + "; ".join(refs) + ".")
        L.append("")
        L.append("%d de %d dentro de tolerancia." % (n_ok, len(fs)))
        L.append("")
        L.append("| Prueba | Esperado | Obtenido | Error | Tolerancia | OK |")
        L.append("|---|---|---|---|---|---|")
        for f in fs:
            err = _pct(f["error_rel"]) if f["error_rel"] is not None else (
                _num(f["error_abs"]) if f["error_abs"] is not None else "—")
            L.append("| %s | %s | %s | %s | %s | %s |" % (
                _celda(f["prueba"]), _num(f["esperado"]), _num(f["obtenido"]),
                err, _celda(f["tolerancia"]), "sí" if f["ok"] else "**NO**"))
        notas = [f for f in fs if f.get("nota")]
        if notas:
            L.append("")
            for f in notas:
                L.append("- *%s*: %s" % (f["prueba"], f["nota"]))
        L.append("")

    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text("\n".join(L), encoding="utf-8")
    print("->", SALIDA, "(%d comprobaciones, %d fallos)" % (len(filas),
                                                            len(fallos)))


if __name__ == "__main__":
    main()
