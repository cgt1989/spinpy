# -*- coding: utf-8 -*-
"""
idioma.py — Traduccion de la interfaz: espanol (original) e ingles.

POR QUE ASI Y NO CON QTranslator
---------------------------------
Qt trae su propio sistema (`tr()`, `.ts`, `lrelease`). Se descarto por dos
razones concretas de este proyecto:

  * Obliga a una cadena de herramientas externa -`pylupdate5`, `lrelease`-
    para regenerar los `.qm` cada vez que cambia un texto. Este codigo se
    empaqueta con PyInstaller desde una maquina sin esas herramientas, y un
    paso de compilacion que se puede olvidar acaba distribuyendo traducciones
    viejas sin que nadie lo note.
  * Los textos de esta interfaz NO son etiquetas: varios tooltips son
    parrafos que explican por que una metrica no es citable o de donde sale
    un sesgo. Son documentacion, y conviene tenerlos donde se puedan leer y
    revisar juntos, no repartidos en un XML binario.

LA CLAVE ES EL TEXTO EN ESPANOL
--------------------------------
El diccionario va de espanol a ingles, y la clave es la propia cadena
espanola. Tres consecuencias, dos buenas y una que hay que vigilar:

  + El codigo fuente sigue leyendose en espanol, que es como esta escrito
    todo el proyecto y como lo lee su autor.
  + Una cadena sin traducir sale EN ESPANOL, no como `MISSING_KEY_42`. La
    interfaz queda mezclada, que es feo, pero nunca rota.
  - Si alguien edita un texto espanol en `visor.py` y no toca el
    diccionario, la traduccion inglesa se pierde en silencio. Para eso esta
    `revisar()`, que compara ambos lados y lista los descuadres; el guion
    `idioma_revisar.py` lo ejecuta.

COMO SE APLICA UN CAMBIO DE IDIOMA
-----------------------------------
No se reconstruye la ventana. Perder el VOI cargado, las metricas medidas y
los campos calculados por cambiar de idioma seria inaceptable: son minutos u
horas de calculo.

En vez de eso, `capturar()` recorre el arbol de widgets UNA VEZ, al arrancar,
y anota el texto original de cada uno. `aplicar()` vuelve a ponerlos
traducidos. Como el original queda guardado, cambiar de idioma es idempotente
y se puede ir y volver sin degradar nada.

Lo que ese recorrido NO ve son los textos que se generan en marcha -mensajes
de la barra de estado, cuadros de dialogo, informes en HTML-. Esos llevan
`_()` explicito en el punto donde se construyen.
"""
from __future__ import annotations

IDIOMAS = {"es": "Español", "en": "English"}

_actual = "es"

# Se rellena desde `idioma_textos.py`, que es solo datos. Se separan para que
# este modulo -la logica- se pueda leer entero sin desplazarse por 5000
# palabras de traduccion.
from .idioma_textos import EN                                  # noqa: E402


def idioma() -> str:
    return _actual


def fijar_idioma(cod: str) -> str:
    global _actual
    _actual = cod if cod in IDIOMAS else "es"
    return _actual


def _(texto: str) -> str:
    """Traduce `texto` al idioma activo. En espanol devuelve el original."""
    if _actual == "es" or not texto:
        return texto
    return EN.get(texto, texto)


# ---------------------------------------------------------------------------
# Recorrido del arbol de widgets
# ---------------------------------------------------------------------------
#
# Cada entrada del registro es (widget, propiedad, original). `propiedad`
# dice como volver a ponerlo; los combos y las tablas necesitan indice, asi
# que llevan una tupla.

def capturar(raiz):
    """Anota el texto ORIGINAL de cada widget que cuelgue de `raiz`.

    Llamar UNA sola vez, con la interfaz recien construida y todavia en
    espanol. Si se llamase despues de traducir, el ingles quedaria grabado
    como original y volver al espanol seria imposible.
    """
    from PyQt5 import QtWidgets

    reg = []

    def anota(w, prop, valor, extra=None):
        if valor and isinstance(valor, str) and valor.strip():
            reg.append((w, prop, valor, extra))

    hijos = raiz.findChildren(QtWidgets.QWidget)
    for w in [raiz] + hijos:
        # Los widgets cuyo texto se construye en marcha se excluyen: la
        # etiqueta de un deslizador dice "Densidad relativa: 35 %", y grabar
        # eso como original significaria buscar en el diccionario una cadena
        # que cambia cada vez que se mueve la barra. Esos se retraducen ellos
        # solos; ver `Deslizador.retraducir` en visor.py.
        if w.property("i18n_dinamico"):
            continue
        # Texto principal
        if isinstance(w, (QtWidgets.QLabel, QtWidgets.QPushButton,
                          QtWidgets.QCheckBox, QtWidgets.QRadioButton)):
            anota(w, "texto", w.text())
        if isinstance(w, QtWidgets.QGroupBox):
            anota(w, "titulo", w.title())
        if isinstance(w, QtWidgets.QComboBox):
            for i in range(w.count()):
                anota(w, "combo", w.itemText(i), i)
        if isinstance(w, QtWidgets.QTableWidget):
            for c in range(w.columnCount()):
                it = w.horizontalHeaderItem(c)
                if it is not None:
                    anota(w, "cabecera", it.text(), c)
        # Ayuda emergente: la lleva cualquier widget
        anota(w, "ayuda", w.toolTip())
        # Titulo de ventana: solo lo tienen la principal y los dialogos, pero
        # es el texto mas visible de todos y se olvidaba con facilidad.
        if w.isWindow():
            anota(w, "ventana", w.windowTitle())

    # Menus y acciones: no son QWidget, hay que pedirlos aparte
    for m in raiz.findChildren(QtWidgets.QMenu):
        anota(m, "titulo_menu", m.title())
    for a in raiz.findChildren(QtWidgets.QAction):
        anota(a, "accion", a.text())
        anota(a, "ayuda", a.toolTip())

    return reg


def aplicar(registro):
    """Reescribe cada widget del registro en el idioma activo."""
    for w, prop, original, extra in registro:
        txt = _(original)
        try:
            if prop == "texto":
                w.setText(txt)
            elif prop in ("titulo", "titulo_menu"):
                w.setTitle(txt)
            elif prop == "accion":
                w.setText(txt)
            elif prop == "ayuda":
                w.setToolTip(txt)
            elif prop == "ventana":
                w.setWindowTitle(txt)
            elif prop == "combo":
                # setItemText y NO clear()+addItems(): vaciar el combo
                # dispara currentIndexChanged y, en esta interfaz, eso relanza
                # un calculo o cambia el modo de medida. Cambiar de idioma no
                # puede tener efectos secundarios.
                w.setItemText(extra, txt)
            elif prop == "cabecera":
                it = w.horizontalHeaderItem(extra)
                if it is not None:
                    it.setText(txt)
        except Exception:
            # Un widget destruido no puede impedir que se traduzca el resto.
            pass


# ---------------------------------------------------------------------------
# Revision del diccionario
# ---------------------------------------------------------------------------

def revisar(cadenas_usadas):
    """Compara las cadenas que usa la interfaz con las del diccionario.

    Devuelve (sin_traducir, sobrantes). La primera lista es la que importa:
    son textos que saldrian en espanol dentro de la interfaz inglesa.
    """
    usadas = {c for c in cadenas_usadas if c and c.strip()}
    sin = sorted(c for c in usadas if c not in EN)
    sobra = sorted(k for k in EN if k not in usadas)
    return sin, sobra
