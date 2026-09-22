# -*- coding: utf-8 -*-
"""
idioma_revisar.py — Comprueba que el diccionario de traduccion esta al dia.

POR QUE HACE FALTA
  El diccionario esta indexado por la cadena ESPANOLA. Es lo que mantiene el
  codigo legible y hace que un texto sin traducir salga en espanol en vez de
  romperse, pero tiene un filo: si alguien corrige una tilde o una coma en
  `visor.py` y no toca `idioma_textos.py`, la traduccion inglesa deja de
  encontrarse SIN NINGUN ERROR. El texto vuelve al espanol y nadie se entera
  hasta que un revisor de habla inglesa lo ve.

  Este guion construye la ventana de verdad, recorre su arbol de widgets y
  compara lo que hay con lo que el diccionario conoce.

QUE COMPRUEBA
  1. Textos de la interfaz que no estan traducidos.
  2. Entradas del diccionario que ya no usa nadie (texto renombrado).
  3. Marcadores {} que no cuadran entre las dos versiones: un {n} que se
     pierde en la traduccion es un KeyError en marcha, no una fealdad.

Uso:   python idioma_revisar.py            (resumen)
       python idioma_revisar.py -v         (lista completa, con repr)
Devuelve 1 si falta algo por traducir.
"""
from __future__ import annotations

import ast
import io
import os
import re
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
sys.path.insert(0, str(RAIZ))

MARCADOR = re.compile(r"\{([a-zA-Z_][a-zA-Z_0-9]*)")

# Cadenas que llegan a `_()` como VARIABLE y no como literal, de modo que el
# analisis del codigo no puede verlas. Son pocas y hay que declararlas aqui a
# mano; si no, el comprobador las daria por sobrantes y alguien acabaria
# borrandolas del diccionario.
DINAMICAS = {
    "malla", "voxel",        # `_(self._modo_medida())` en la barra de estado
    "inicial",               # `_(f['fase'])` en la tabla de DialogoSimulacion
}


def _textos_de_replicas():
    """Los textos de la tabla REPLICAS de la ventana de Validacion.

    Cabeceras, pies y lineas de parametros de cada replica llegan a `_()` como
    VARIABLE —salen de un diccionario, no de un literal— y el analisis del
    codigo no los ve. Se leen de la propia tabla en vez de mantener una lista
    a mano: asi, anadir una cuarta replica no deja sus textos sin traducir en
    silencio, que es exactamente el fallo que este guion existe para evitar.
    """
    try:
        from dialogo_validacion import REPLICAS
    except Exception:
        return set()
    fuera = set()
    for r in REPLICAS:
        for clave in ("cabecera", "params", "ayuda_rapido", "pie_izq"):
            v = r.get(clave)
            if isinstance(v, str) and len(v) > 3:
                fuera.add(v)
    return fuera


def _literales(nodo):
    """Cadenas literales de un argumento, incluidas las de un `a if c else b`.

    Se baja SOLO por condicionales -un criterio puede depender de la
    resolucion, como el umbral de ortotropia- y no por el arbol entero. Bajar
    por todo recogeria los trozos de los f-strings y las claves de los
    diccionarios que aparezcan dentro del argumento ("E_max_direccional" y
    demas), que no son textos de la interfaz y llenarian el diccionario de
    basura.
    """
    if isinstance(nodo, ast.Constant) and isinstance(nodo.value, str)             and len(nodo.value) > 3:
        return {nodo.value}
    if isinstance(nodo, ast.IfExp):
        return _literales(nodo.body) | _literales(nodo.orelse)
    return set()


def _marcadores(t):
    return set(MARCADOR.findall(t or ""))


def main():
    detalle = "-v" in sys.argv

    # La consola de Windows va en cp1252 y las lineas de parametros llevan rho,
    # beta y theta. Sin esto, el guion REVIENTA al imprimir una cadena sin
    # traducir que contenga una letra griega -es decir, justo cuando tiene algo
    # que decir- y `actualizar.ps1` lo interpreta como que el diccionario esta
    # roto por otra razon.
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    import pyvista as pv
    from PyQt5 import QtWidgets

    from spinpy.idioma import capturar, revisar
    from spinpy.idioma_textos import EN

    pv.set_plot_theme("document")
    # HAY QUE GUARDAR LA REFERENCIA, aunque no se use.
    #
    # Sin `app = ...` el objeto queda sin nadie que lo sujete, el recolector
    # se lo lleva, y construir la ventana despues mata el proceso de golpe:
    # sin traza, sin salida y con un codigo de salida sin sentido. Parece que
    # el guion "no imprime nada". No quitar por parecer una variable sin usar.
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    assert app is not None
    import visor as V

    v = V.Visor()

    # Se usa la captura QUE HIZO LA PROPIA VENTANA, no una nueva.
    #
    # `Visor.__init__` captura los textos y DESPUES aplica el idioma guardado
    # en QSettings. Volver a capturar aqui, con la ventana ya traducida,
    # grabaria las cadenas INGLESAS como si fueran los originales: ninguna
    # esta en el diccionario -que va indexado por el espanol- y el
    # comprobador las daria todas por sin traducir. Es decir, fallaria justo
    # cuando el usuario tiene la aplicacion en ingles.
    textos = getattr(v, "_textos", None) or capturar(v)
    usadas = {orig for _w, _p, orig, _e in textos}

    # Los deslizadores quedan fuera del recorrido -su etiqueta lleva el valor
    # dentro-, asi que su texto base se recoge aparte. Sin esto, renombrar un
    # deslizador dejaria su traduccion muerta sin que nadie lo notase, que es
    # exactamente el fallo que este guion existe para evitar.
    usadas |= {d._texto for d in v.findChildren(V.Deslizador)}

    # Las cadenas de marcha -mensajes, dialogos, informes- no cuelgan de
    # ningun widget, asi que se sacan del CODIGO: son las envueltas en _().
    #
    # Con AST y no con una expresion regular. La primera version usaba regex y
    # se dejaba en silencio las cadenas partidas en tres o mas trozos -que
    # aqui son casi todos los mensajes largos-: no aparecian como pendientes
    # porque el comprobador ni siquiera sabia que existian. Un comprobador con
    # puntos ciegos es peor que no tenerlo: da tranquilidad falsa.
    for ruta in (RAIZ / "visor.py", RAIZ / "dialogo_metodos.py",
                 RAIZ / "dialogo_validacion.py"):
        if not ruta.exists():
            continue
        arbol = ast.parse(ruta.read_text(encoding="utf-8"))
        for nodo in ast.walk(arbol):
            if (isinstance(nodo, ast.Call)
                    and isinstance(nodo.func, ast.Name)
                    and nodo.func.id == "_"
                    and len(nodo.args) == 1
                    and isinstance(nodo.args[0], ast.Constant)
                    and isinstance(nodo.args[0].value, str)):
                usadas.add(nodo.args[0].value)

    # Las predicciones y los criterios de la replica llegan a la interfaz
    # DESDE UN JSON, asi que `_()` los recibe como variable y el analisis del
    # codigo no los ve. Se sacan del guion que los define, leyendo los
    # literales que pasan por `anota(...)`. Asi no hay que mantener una lista
    # a mano que se quedaria atras en cuanto se anada una comprobacion.
    for rep in sorted((RAIZ / "Test").glob("replicar_*.py")):
        arbol = ast.parse(rep.read_text(encoding="utf-8"))
        for nodo in ast.walk(arbol):
            if (isinstance(nodo, ast.Call)
                    and isinstance(nodo.func, ast.Name)
                    and nodo.func.id == "anota"):
                for arg in nodo.args:
                    usadas |= _literales(arg)

    usadas |= DINAMICAS
    usadas |= _textos_de_replicas()

    sin, sobra = revisar(usadas)

    # Marcadores descuadrados
    descuadre = []
    for k, val in EN.items():
        a, b = _marcadores(k), _marcadores(val)
        if a != b:
            descuadre.append((k, sorted(a), sorted(b)))

    print("=" * 74)
    print(" REVISION DEL DICCIONARIO DE IDIOMA")
    print("=" * 74)
    print(f"  cadenas de la interfaz     : {len(usadas)}")
    print(f"  entradas del diccionario   : {len(EN)}")
    print(f"  SIN TRADUCIR               : {len(sin)}")
    print(f"  sobrantes (ya no se usan)  : {len(sobra)}")
    print(f"  marcadores descuadrados    : {len(descuadre)}")

    if descuadre:
        print("\n  MARCADORES DESCUADRADOS (esto revienta en marcha):")
        for k, a, b in descuadre:
            print(f"    es{a} != en{b}   {k[:60]!r}")

    if sin:
        print(f"\n  SIN TRADUCIR ({len(sin)}):")
        for c in sin:
            print(f"    {c!r}" if detalle else
                  f"    {c[:100]}{'…' if len(c) > 100 else ''}")

    if sobra and detalle:
        print(f"\n  SOBRANTES ({len(sobra)}):")
        for c in sobra:
            print(f"    {c[:100]}")

    if "--esqueleto" in sys.argv and sin:
        # Se vuelca un diccionario listo para pegar, con las claves EXACTAS.
        # Transcribirlas a mano es la forma segura de equivocarse en una
        # tilde y perder la traduccion sin enterarse.
        d = RAIZ / "idioma_pendiente.py"
        lineas = ["# Generado por idioma_revisar.py --esqueleto",
                  "# Rellena el lado ingles y pega en spinpy/idioma_textos.py",
                  "",
                  "PENDIENTE = {"]
        for c in sin:
            lineas.append("    %r:" % c)
            lineas.append('        "",')
        lineas.append("}")
        io.open(d, "w", encoding="utf-8").write("\n".join(lineas) + "\n")
        print("\n  esqueleto -> %s" % d)

    v.close()
    print("=" * 74)
    return 1 if (sin or descuadre) else 0


if __name__ == "__main__":
    rc = main()
    # SALIDA DURA, y no `sys.exit`.
    #
    # Construir la ventana deja vivos un panel de VTK y el hilo que genera la
    # primera vista. Al no entrar nunca en el bucle de eventos, el desmontaje
    # de Qt y VTK al terminar el interprete revienta el proceso
    # (STATUS_STACK_BUFFER_OVERRUN). Eso no afecta al resultado -ya esta
    # calculado- pero SI se lleva por delante la salida cuando esta
    # canalizada, porque Python la tiene en un buffer que nunca se vacia: el
    # guion parece no imprimir nada y su codigo de salida es basura.
    #
    # Se vacia a mano y se sale sin desmontar nada.
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(rc)
