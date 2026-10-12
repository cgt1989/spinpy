"""
pantalla.py: dialogos que se ajustan a la pantalla del equipo.

Un dialogo con muchos controles (Informe FEM, Informe automatico, Ajustar con
todos los metodos...) mide mas de lo que cabe en una pantalla de 1366x768 o en
una de 1920x1080 con escala del 150 %. Como el minimo del layout supera la
pantalla, Qt no deja encogerlo y la ventana queda cortada, sin manera de
verla entera.

`DialogoAjustable` arregla eso para cualquier QDialog sin tocar su contenido:

  * la primera vez que se muestra, mete el layout en una QScrollArea, de modo
    que el minimo del dialogo ya no depende de lo que contiene. Una fila final
    de botones (QDialogButtonBox o una fila de QPushButton) se deja FUERA del
    desplazamiento, para que Empezar/Cancelar no se pierdan;
  * limita el tamano inicial al area util de la pantalla donde aparece (sin
    la barra de tareas) y lo centra;
  * habilita maximizar y minimizar, ademas de redimensionar con el raton.

El tamano que el dialogo pedia (`resize`) se respeta mientras quepa.
"""
from __future__ import annotations

from PyQt5 import QtCore, QtGui, QtWidgets, sip

# Fraccion maxima del area util que puede ocupar un dialogo al abrirse.
FRACCION = 0.94
# Margen que se suma al contenido para la barra de titulo y los bordes.
_MARGEN = 40


def area_util(widget=None) -> QtCore.QRect:
    """Rectangulo disponible (sin barra de tareas) de la pantalla del widget,
    o de la pantalla principal si aun no esta en ninguna."""
    app = QtWidgets.QApplication.instance()
    pantalla = None
    if widget is not None:
        wh = widget.windowHandle()
        if wh is not None:
            pantalla = wh.screen()
        if pantalla is None and widget.parentWidget() is not None:
            wh = widget.parentWidget().window().windowHandle()
            if wh is not None:
                pantalla = wh.screen()
    if pantalla is None and app is not None:
        pantalla = app.primaryScreen()
    if pantalla is None:
        return QtCore.QRect(0, 0, 1280, 720)
    return pantalla.availableGeometry()


def _es_fila_de_botones(item) -> bool:
    w = item.widget()
    if w is not None:
        return isinstance(w, QtWidgets.QDialogButtonBox)
    lay = item.layout()
    if lay is None or lay.count() == 0:
        return False
    hay_boton = False
    for i in range(lay.count()):
        sub = lay.itemAt(i)
        sw = sub.widget()
        if sw is None:
            continue                       # espaciadores y stretch
        if isinstance(sw, (QtWidgets.QPushButton,
                           QtWidgets.QDialogButtonBox)):
            hay_boton = True
        else:
            return False
    return hay_boton


def envolver_en_scroll(dlg: QtWidgets.QDialog) -> bool:
    """Mueve el layout de `dlg` a una QScrollArea. Devuelve False si no hay
    nada que envolver (sin layout, o ya envuelto)."""
    lay = dlg.layout()
    if (lay is None or dlg.property("_pantalla_envuelto")
            or not isinstance(lay, QtWidgets.QBoxLayout)):
        return False
    pie = None
    if lay.count() >= 2:
        if _es_fila_de_botones(lay.itemAt(lay.count() - 1)):
            pie = lay.takeAt(lay.count() - 1)

    margenes = lay.contentsMargins()
    # Un layout no puede cambiar de propietario: se rehace uno equivalente
    # para el contenido y se le pasan los elementos con sus estiramientos.
    nuevo = (QtWidgets.QHBoxLayout()
             if isinstance(lay, QtWidgets.QHBoxLayout)
             else QtWidgets.QVBoxLayout())
    nuevo.setSpacing(lay.spacing())
    nuevo.setContentsMargins(margenes)
    estira = [lay.stretch(i) for i in range(lay.count())]
    i = 0
    while lay.count():
        nuevo.addItem(lay.takeAt(0))
        nuevo.setStretch(i, estira[i])
        i += 1
    contenido = QtWidgets.QWidget()
    contenido.setLayout(nuevo)
    # El layout vacio original se retira del dialogo para poder poner el raiz.
    sip.delete(lay)

    area = QtWidgets.QScrollArea()
    area.setWidgetResizable(True)
    area.setFrameShape(QtWidgets.QFrame.NoFrame)
    area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
    area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
    area.setWidget(contenido)

    raiz = QtWidgets.QVBoxLayout(dlg)
    raiz.setContentsMargins(0, 0, 0, 0)
    raiz.setSpacing(0)
    raiz.addWidget(area, 1)
    if pie is not None:
        marco = QtWidgets.QWidget()
        pl = QtWidgets.QVBoxLayout(marco)
        pl.setContentsMargins(margenes.left(), 4, margenes.right(),
                              margenes.bottom())
        if pie.widget() is not None:
            pl.addWidget(pie.widget())
        else:
            pl.addLayout(pie.layout())
        raiz.addWidget(marco, 0)
    dlg.setProperty("_pantalla_envuelto", True)
    dlg._area_scroll = area
    return True


def ajustar_tamano(dlg: QtWidgets.QDialog) -> None:
    """Fija el tamano inicial: lo que pide el contenido o el `resize` previo,
    sin pasar de la fraccion permitida del area util; y centra."""
    r = area_util(dlg)
    max_w = int(r.width() * FRACCION)
    max_h = int(r.height() * FRACCION)
    dlg.setMaximumSize(r.width(), r.height())
    pedido = dlg.size()
    hint = dlg.sizeHint()
    ancho = max(pedido.width(), hint.width() + _MARGEN)
    alto = max(pedido.height(), hint.height() + _MARGEN)
    ancho, alto = min(ancho, max_w), min(alto, max_h)
    dlg.resize(ancho, alto)
    marco = dlg.frameGeometry()
    marco.moveCenter(r.center())
    x = max(r.left(), min(marco.left(), r.right() - ancho))
    y = max(r.top(), min(marco.top(), r.bottom() - alto))
    dlg.move(x, y)


class DialogoAjustable(QtWidgets.QDialog):
    """QDialog que nunca excede la pantalla y se puede redimensionar."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.setWindowFlags(self.windowFlags()
                            | QtCore.Qt.WindowMaximizeButtonHint
                            | QtCore.Qt.WindowMinimizeButtonHint)
        self.setSizeGripEnabled(True)
        self._ajustado = False

    def setVisible(self, visible):
        # Antes de mostrarse por primera vez (exec_, show, open): asi la
        # ventana aparece ya con su tamano final, sin parpadeo.
        if visible and not self._ajustado:
            self._ajustado = True
            envolver_en_scroll(self)
            ajustar_tamano(self)
        super().setVisible(visible)
