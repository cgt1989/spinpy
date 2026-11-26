"""
dialogo_cargas.py: Análisis mecánico con cargas elegidas con el ratón.

POR QUÉ ESTA VENTANA
--------------------
«Resolver y estimar el fallo» hace una pregunta fija: compresión axial entre dos
platos. Aquí la pregunta la formula quien mira el cubo: dónde se apoya, dónde
se empuja o se tira, con qué fuerza y en qué dirección. Se elige sobre el VOI
del hueso, el espinodoide, el dual-lattice o el level-set foam, con el mismo cubo y el mismo caso
de carga, para poder comparar las tres estructuras bajo la misma solicitación.

DISTRIBUCIÓN
------------
Es la de «Recortar VOI»: el 3D a la izquierda, tres cortes debajo y los mandos
agrupados a la derecha.

  3D        la estructura, el cubo con el nombre de sus caras, cada carga como
            una flecha de su color sobre su región, los apoyos en azul y, al
            resolver, el campo elegido sobre la superficie (deformada
            amplificada y corte opcionales).
  Cortes    tres planos por el cubo, del campo si hay resultado y de la
            estructura si no, cada uno con su deslizador.
  Mandos    Estructura, Plantilla, Cargas, Apoyos, Material y solver,
            Resultados, Barrido de ángulo y Exportar.

CON EL RATÓN (herramienta «Elegir con el ratón»)
------------------------------------------------
  Cara del cubo    un clic sobre el cubo o sobre el hueso elige la cara del
                   cubo más cercana (+x … -z).
  Parche           un clic sobre la superficie del hueso elige las caras de
                   voxel de su entorno (radio ajustable) que miran al punto.
  Nodos            un clic añade el nodo de superficie más cercano; con Mayús
                   lo quita. Es la carga o el apoyo puntual exacto.
  Punta de flecha  una esfera en la punta de la flecha de la carga activa:
                   arrastrarla fija la dirección.
La región elegida se asigna a la carga o apoyo activo (el último que se tocó en
las listas). Todo cambio se puede deshacer y rehacer (Ctrl+Z, Ctrl+Mayús+Z).

NADA CALCULA AQUÍ
-----------------
Las regiones, las fuerzas nodales, el solver, las reacciones y el barrido son de
`spinpy.cargas`; el dibujo, de `spinpy.escena_cargas`. Esta ventana construye
el caso y lo manda al solver en un hilo. Los textos van a `_()` como LITERALES,
porque `idioma_revisar.py` los recoge del árbol sintáctico de este archivo.
"""

from __future__ import annotations

import copy
import json
import os
import time
import traceback
from pathlib import Path

import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets

from spinpy import cargas as C
from spinpy.elastic import remuestrear_bw
from spinpy.idioma import _
from spinpy.pantalla import DialogoAjustable

try:                                   # el 3D es opcional, como en el recorte
    from pyvistaqt import QtInteractor
except Exception:                      # pragma: no cover
    QtInteractor = None
# Sin OpenGL, VTK no lanza una excepcion: termina el proceso. Las pruebas que
# solo ejercen los mandos (no el dibujo) arrancan la ventana con esta variable
# y la vista 3D desactivada.
if os.environ.get("SPINPY_SIN_3D"):
    QtInteractor = None

AZUL = "#1f6fb5"
LIMITE_HISTORIAL = 100


# ---------------------------------------------------------------------------
# Piezas pequeñas
# ---------------------------------------------------------------------------

class _Hilo(QtCore.QThread):
    """Una función lenta fuera del hilo de la interfaz."""
    listo = QtCore.pyqtSignal(object)
    fallo = QtCore.pyqtSignal(str)
    avance = QtCore.pyqtSignal(int, str)

    def __init__(self, fn):
        super().__init__()
        self._fn = fn

    def informar(self, f, m):
        self.avance.emit(int(round(100 * f)), str(m))

    def run(self):
        try:
            self.listo.emit(self._fn(self.informar))
        except Exception:
            self.fallo.emit(traceback.format_exc())


def _spin(lo, hi, val, dec=3, suf="", paso=None, ancho=84):
    s = QtWidgets.QDoubleSpinBox()
    s.setDecimals(dec)
    s.setRange(lo, hi)
    s.setValue(val)
    if suf:
        s.setSuffix(" " + suf)
    if paso:
        s.setSingleStep(paso)
    s.setKeyboardTracking(False)
    s.setMinimumWidth(ancho)
    return s


def _fila(*ws):
    w = QtWidgets.QWidget()
    h = QtWidgets.QHBoxLayout(w)
    h.setContentsMargins(0, 0, 0, 0)
    for x in ws:
        h.addWidget(x) if not isinstance(x, str) else h.addWidget(
            QtWidgets.QLabel(x))
    return w


def _lut_turbo():
    try:
        import matplotlib
        cm = matplotlib.colormaps["turbo"]
        return (cm(np.linspace(0, 1, 256))[:, :3] * 255).astype(np.uint8)
    except Exception:
        t = np.linspace(0, 1, 256)
        return np.stack([255 * t, 255 * (1 - abs(2 * t - 1)),
                         255 * (1 - t)], axis=1).astype(np.uint8)


_LUT = None


class _CorteCampo(QtWidgets.QWidget):
    """Un corte plano del cubo: el campo con su escala de color, o el hueso."""
    FONDO = (28, 34, 43)
    HUESO = (217, 199, 160)
    NODATO = (120, 120, 120)

    def __init__(self, titulo):
        super().__init__()
        self._titulo = titulo
        self._qi = None
        self._rango = None
        self._pie = ""
        self.setMinimumSize(110, 90)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                           QtWidgets.QSizePolicy.Expanding)

    def poner(self, rgb, pie="", rango=None):
        """`rgb` (nb, na, 3) uint8 con la fila 0 abajo."""
        rgb = np.ascontiguousarray(rgb[::-1])
        nb, na = rgb.shape[:2]
        self._qi = QtGui.QImage(rgb.data, na, nb, 3 * na,
                                QtGui.QImage.Format_RGB888).copy()
        self._pie, self._rango = pie, rango
        self.update()

    def paintEvent(self, _ev):
        p = QtGui.QPainter(self)
        p.fillRect(self.rect(), QtGui.QColor(*self.FONDO))
        if self._qi is None:
            return
        franja = 18
        w, h = self.width(), self.height() - franja
        s = min(w / self._qi.width(), h / self._qi.height())
        ox = (w - s * self._qi.width()) / 2.0
        oy = franja + (h - s * self._qi.height()) / 2.0
        p.drawImage(QtCore.QRectF(ox, oy, s * self._qi.width(),
                                  s * self._qi.height()), self._qi)
        p.setPen(QtGui.QColor(235, 235, 235))
        f = p.font()
        f.setPointSizeF(8.5)
        p.setFont(f)
        p.drawText(QtCore.QPointF(4, 13), self._titulo + "  " + self._pie)
        if self._rango:
            p.drawText(QtCore.QPointF(ox + 4, oy + s * self._qi.height() - 4),
                       "%.3g … %.3g" % self._rango)


class _Rosa(QtWidgets.QWidget):
    """Gráfico polar de una magnitud frente al ángulo de la carga."""

    def __init__(self):
        super().__init__()
        self._d = None
        self.setMinimumSize(300, 300)

    def poner(self, theta_deg, valores, rotulo):
        self._d = (np.asarray(theta_deg, float), np.asarray(valores, float),
                   rotulo)
        self.update()

    def paintEvent(self, _ev):
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        p.fillRect(self.rect(), QtGui.QColor("#ffffff"))
        if self._d is None:
            return
        th, v, rot = self._d
        ok = np.isfinite(v)
        if not ok.any():
            p.drawText(self.rect(), QtCore.Qt.AlignCenter, _("Sin datos"))
            return
        vmax = float(v[ok].max()) or 1.0
        cx, cy = self.width() / 2.0, self.height() / 2.0 + 8
        R = min(cx, cy) - 28
        p.setPen(QtGui.QPen(QtGui.QColor("#cccccc"), 1))
        for f in (0.25, 0.5, 0.75, 1.0):
            p.drawEllipse(QtCore.QPointF(cx, cy), R * f, R * f)
        for a in range(0, 360, 45):
            t = np.radians(a)
            p.drawLine(QtCore.QPointF(cx, cy), QtCore.QPointF(
                cx + R * np.cos(t), cy - R * np.sin(t)))
            p.drawText(QtCore.QPointF(cx + (R + 6) * np.cos(t) - 8,
                                      cy - (R + 6) * np.sin(t) + 4),
                       "%d°" % a)
        pts = [QtCore.QPointF(cx + R * x / vmax * np.cos(np.radians(t)),
                              cy - R * x / vmax * np.sin(np.radians(t)))
               for t, x in zip(th, v) if np.isfinite(x)]
        if pts:
            pts.append(pts[0])
            p.setPen(QtGui.QPen(QtGui.QColor(AZUL), 2))
            p.drawPolyline(QtGui.QPolygonF(pts))
        p.setPen(QtGui.QColor("#222222"))
        p.drawText(QtCore.QPointF(8, 14),
                   "%s (máx. %.4g)" % (rot, vmax))


class _EditorRegion(QtWidgets.QWidget):
    """Los parámetros de una región, como los guarda `spinpy.cargas`.

    Emite `cambiado` solo ante gestos del usuario, nunca al poner una región
    desde fuera (`poner`), para que elegir con el ratón no dispare un bucle.
    """
    cambiado = QtCore.pyqtSignal()

    def __init__(self, ext, paso):
        super().__init__()
        self._ext = np.asarray(ext, float)
        self._silencio = False
        self._ids = []
        self._xyz = []
        L = float(self._ext.max())
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(2)

        self.cmb_forma = QtWidgets.QComboBox()
        for k, t in (("cara", _("Cara del cubo")),
                     ("esfera", _("Parche esférico")),
                     ("caja", _("Caja")), ("cilindro", _("Cilindro")),
                     ("nodos", _("Nodos elegidos"))):
            self.cmb_forma.addItem(t, k)
        v.addWidget(_fila(QtWidgets.QLabel(_("Región")), self.cmb_forma))

        self.pila = QtWidgets.QStackedWidget()
        v.addWidget(self.pila)

        def xyz(val):
            return [_spin(-L, 2 * L, val[k], 3, "mm", paso) for k in range(3)]

        # cara
        w = QtWidgets.QWidget()
        h = QtWidgets.QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        self.cmb_cara = QtWidgets.QComboBox()
        self.cmb_cara.addItems(list(C.CARAS))
        self.sp_capa = QtWidgets.QSpinBox()
        self.sp_capa.setRange(0, 50)
        self.sp_capa.setToolTip(_(
            "Capas de voxel desde la cara. Solo cuenta con ámbito «volumen»: "
            "los nodos a menos de esa profundidad."))
        h.addWidget(QtWidgets.QLabel(_("Cara")))
        h.addWidget(self.cmb_cara)
        h.addWidget(QtWidgets.QLabel(_("Capas")))
        h.addWidget(self.sp_capa)
        h.addStretch(1)
        self.pila.addWidget(w)
        # esfera
        w = QtWidgets.QWidget()
        g = QtWidgets.QGridLayout(w)
        g.setContentsMargins(0, 0, 0, 0)
        self.e_centro = xyz(self._ext / 2)
        self.sp_radio = _spin(1e-4, 10 * L, 0.15 * L, 3, "mm", paso)
        g.addWidget(QtWidgets.QLabel(_("Centro")), 0, 0)
        for k, s in enumerate(self.e_centro):
            g.addWidget(s, 0, 1 + k)
        g.addWidget(QtWidgets.QLabel(_("Radio")), 1, 0)
        g.addWidget(self.sp_radio, 1, 1)
        self.pila.addWidget(w)
        # caja
        w = QtWidgets.QWidget()
        g = QtWidgets.QGridLayout(w)
        g.setContentsMargins(0, 0, 0, 0)
        self.e_lo = xyz(np.zeros(3))
        self.e_hi = xyz(self._ext)
        g.addWidget(QtWidgets.QLabel(_("Mínimo")), 0, 0)
        g.addWidget(QtWidgets.QLabel(_("Máximo")), 1, 0)
        for k in range(3):
            g.addWidget(self.e_lo[k], 0, 1 + k)
            g.addWidget(self.e_hi[k], 1, 1 + k)
        self.pila.addWidget(w)
        # cilindro
        w = QtWidgets.QWidget()
        g = QtWidgets.QGridLayout(w)
        g.setContentsMargins(0, 0, 0, 0)
        self.c_centro = xyz(self._ext / 2)
        self.c_radio = _spin(1e-4, 10 * L, 0.2 * L, 3, "mm", paso)
        self.c_alto = _spin(1e-4, 10 * L, 0.5 * L, 3, "mm", paso)
        self.c_eje = QtWidgets.QComboBox()
        self.c_eje.addItems(["x", "y", "z"])
        self.c_eje.setCurrentIndex(2)
        g.addWidget(QtWidgets.QLabel(_("Centro")), 0, 0)
        for k, s in enumerate(self.c_centro):
            g.addWidget(s, 0, 1 + k)
        g.addWidget(QtWidgets.QLabel(_("Radio")), 1, 0)
        g.addWidget(self.c_radio, 1, 1)
        g.addWidget(QtWidgets.QLabel(_("Altura")), 1, 2)
        g.addWidget(self.c_alto, 1, 3)
        g.addWidget(QtWidgets.QLabel(_("Eje")), 2, 0)
        g.addWidget(self.c_eje, 2, 1)
        self.pila.addWidget(w)
        # nodos
        w = QtWidgets.QWidget()
        h = QtWidgets.QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        self.lab_nodos = QtWidgets.QLabel()
        self.btn_vaciar = QtWidgets.QPushButton(_("Vaciar"))
        h.addWidget(self.lab_nodos, 1)
        h.addWidget(self.btn_vaciar)
        self.pila.addWidget(w)

        # comunes
        self.cmb_ambito = QtWidgets.QComboBox()
        self.cmb_ambito.addItem(_("Superficie (caras expuestas)"), "superficie")
        self.cmb_ambito.addItem(_("Volumen (todos los nodos)"), "volumen")
        self.cmb_ambito.setToolTip(_(
            "Superficie: las caras de voxel con hueso a un lado y vacío al "
            "otro; una carga se reparte por área. Volumen: todos los nodos "
            "de hueso dentro de la región; una fuerza se reparte por nodo."))
        v.addWidget(_fila(QtWidgets.QLabel(_("Ámbito")), self.cmb_ambito))
        self.chk_filtro = QtWidgets.QCheckBox(_("Solo caras que miran a…"))
        self.chk_filtro.setToolTip(_(
            "Filtra por la normal: solo las caras cuya normal saliente forma "
            "menos del ángulo indicado con la dirección. Sobre una trabécula "
            "deja las del lado pinchado."))
        self.f_dir = [_spin(-1, 1, k == 2, 3, "", 0.1, 64) for k in range(3)]
        self.sp_ang = _spin(1, 179, 60, 0, "°", 5, 70)
        v.addWidget(_fila(self.chk_filtro, *self.f_dir, self.sp_ang))

        for w in (self.cmb_forma, self.cmb_cara, self.c_eje, self.cmb_ambito):
            w.currentIndexChanged.connect(self._gesto)
        for w in (self.sp_capa, self.sp_radio, self.c_radio, self.c_alto,
                  self.sp_ang, *self.e_centro, *self.e_lo, *self.e_hi,
                  *self.c_centro, *self.f_dir):
            w.valueChanged.connect(self._gesto)
        self.chk_filtro.toggled.connect(self._gesto)
        self.btn_vaciar.clicked.connect(self._vaciar)
        self.cmb_forma.currentIndexChanged.connect(self._ver)
        self._ver()

    # --
    def _ver(self, *_a):
        k = self.cmb_forma.currentIndex()
        self.pila.setCurrentIndex(k)
        # Un QStackedWidget reserva la altura de su pagina MAS ALTA: las que
        # no se ven se declaran sin tamano para que no dejen huecos.
        for i in range(self.pila.count()):
            w = self.pila.widget(i)
            pol = (QtWidgets.QSizePolicy.Preferred if i == k
                   else QtWidgets.QSizePolicy.Ignored)
            w.setSizePolicy(pol, pol)
        self.pila.adjustSize()

    def _gesto(self, *_a):
        if not self._silencio:
            self.cambiado.emit()

    def _vaciar(self):
        self._ids, self._xyz = [], []
        self._rotulo_nodos()
        self._gesto()

    def _rotulo_nodos(self):
        self.lab_nodos.setText(_("{n} nodo(s) elegidos con el ratón").format(
            n=len(self._ids)))

    def poner(self, region):
        self._silencio = True
        try:
            f = region.get("forma", "cara")
            self.cmb_forma.setCurrentIndex(max(0, self.cmb_forma.findData(f)))
            self._ver()
            if f == "cara":
                self.cmb_cara.setCurrentText(str(region.get("cara", "+z")))
                self.sp_capa.setValue(int(region.get("capa", 0)))
            elif f == "esfera":
                for s, x in zip(self.e_centro, region["centro"]):
                    s.setValue(float(x))
                self.sp_radio.setValue(float(region["radio"]))
            elif f == "caja":
                for s, x in zip(self.e_lo, region["lo"]):
                    s.setValue(float(x))
                for s, x in zip(self.e_hi, region["hi"]):
                    s.setValue(float(x))
            elif f == "cilindro":
                for s, x in zip(self.c_centro, region["centro"]):
                    s.setValue(float(x))
                self.c_radio.setValue(float(region["radio"]))
                self.c_alto.setValue(float(region["altura"]))
                self.c_eje.setCurrentIndex(int(region.get("eje", 2)))
            elif f == "nodos":
                self._ids = [int(i) for i in region.get("ids", [])]
                self._xyz = [list(map(float, x)) for x in region.get("xyz", [])]
            self._rotulo_nodos()
            self.cmb_ambito.setCurrentIndex(max(0, self.cmb_ambito.findData(
                region.get("ambito", "superficie"))))
            fn = region.get("filtro_normal")
            self.chk_filtro.setChecked(bool(fn))
            if fn:
                for s, x in zip(self.f_dir, fn["dir"]):
                    s.setValue(float(x))
                self.sp_ang.setValue(float(fn.get("angulo", 60)))
        finally:
            self._silencio = False

    def region(self):
        f = self.cmb_forma.currentData()
        r = {"forma": f}
        if f == "cara":
            r["cara"] = self.cmb_cara.currentText()
            if self.sp_capa.value():
                r["capa"] = int(self.sp_capa.value())
        elif f == "esfera":
            r["centro"] = [s.value() for s in self.e_centro]
            r["radio"] = self.sp_radio.value()
        elif f == "caja":
            r["lo"] = [s.value() for s in self.e_lo]
            r["hi"] = [s.value() for s in self.e_hi]
        elif f == "cilindro":
            r["centro"] = [s.value() for s in self.c_centro]
            r["radio"] = self.c_radio.value()
            r["altura"] = self.c_alto.value()
            r["eje"] = self.c_eje.currentIndex()
        else:
            r["ids"] = list(self._ids)
            r["xyz"] = [list(x) for x in self._xyz]
        r["ambito"] = self.cmb_ambito.currentData()
        if r["ambito"] == "superficie":
            del r["ambito"]
        if f != "nodos" and self.chk_filtro.isChecked():
            r["filtro_normal"] = {"dir": [s.value() for s in self.f_dir],
                                  "angulo": self.sp_ang.value()}
        return r


# ---------------------------------------------------------------------------
# La ventana
# ---------------------------------------------------------------------------


def direcciones():
    return (("vector", _("Vector")), ("+x", "+x"), ("-x", "-x"),
            ("+y", "+y"), ("-y", "-y"), ("+z", "+z"), ("-z", "-z"),
            ("normal_entrante", _("Normal entrante (empuja)")),
            ("normal_saliente", _("Normal saliente (tira)")))


def tipos_carga_txt():
    return (("fuerza", _("Fuerza (N)")), ("presion", _("Presión (MPa)")),
            ("momento", _("Momento (N·mm)")),
            ("desplazamiento", _("Desplazamiento impuesto (mm)")))


def tipos_apoyo_txt():
    return (("empotrado", _("Empotrado (x, y, z)")),
            ("deslizante", _("Deslizante (solo el normal)")),
            ("personalizado", _("Personalizado (grados de libertad)")))


def plantillas_txt():
    return (("compresion", _("Compresión axial")),
            ("traccion", _("Tracción axial")),
            ("cizalla", _("Cizalla")), ("flexion", _("Flexión")),
            ("torsion", _("Torsión")),
            ("puntual", _("Carga puntual sobre la superficie")))


class DialogoCargas(DialogoAjustable):
    """Ver, elegir y resolver un caso de carga sobre el VOI o un candidato."""

    registrado = QtCore.pyqtSignal(object)

    def __init__(self, padre, estructuras, clave=None, caso=None, n_calculo=40,
                 E_s=C.E_S_DEF, nu_s=C.NU_S_DEF, carpeta=None):
        super().__init__(padre)
        self.setWindowTitle(_("Análisis mecánico con cargas"))
        self.resize(1380, 880)
        self._n0 = int(n_calculo)
        self.estructuras = dict(estructuras)
        if not self.estructuras:
            raise ValueError("no hay estructuras")
        self.clave = (clave if clave in self.estructuras
                      else next(iter(self.estructuras)))
        self.carpeta = Path(carpeta) if carpeta else Path.home()
        self.escena = None
        self.res = None
        self.resultados = {}                 # clave -> (res, caso)
        self.registros = []
        self.hilo = None
        self._silencio = False
        self._foco = ("carga", 0)
        self._historial, self._pos = [], -1
        self._ultimo_snap = (None, 0.0)
        self._modo = "orbitar"
        self._nueva = False
        self._cache_sel = {}
        self._BW = self._sp = None
        self.caso = caso if caso else C.caso_vacio()
        self.caso.setdefault("material", {})
        self.caso["material"].setdefault("E_s", E_s)
        self.caso["material"].setdefault("nu", nu_s)
        self._cortes_img = {}

        raiz = QtWidgets.QHBoxLayout(self)

        # ---- izquierda: 3D + cortes ---------------------------------------
        izq = QtWidgets.QVBoxLayout()
        raiz.addLayout(izq, 1)
        izq.addLayout(self._barra_herramientas())
        self.vista3d = None
        if QtInteractor is not None:
            try:
                self.vista3d = QtInteractor(self)
                izq.addWidget(self.vista3d.interactor, 5)
            except Exception:
                self.vista3d = None
        if self.vista3d is None:
            aus = QtWidgets.QLabel(_(
                "La vista 3D no está disponible en este equipo. Los mandos y "
                "los cortes siguen funcionando."))
            aus.setWordWrap(True)
            aus.setAlignment(QtCore.Qt.AlignCenter)
            izq.addWidget(aus, 5)
        self.lab_estado = QtWidgets.QLabel()
        self.lab_estado.setStyleSheet("color:#555; font-size:11px;")
        izq.addWidget(self.lab_estado)
        fila = QtWidgets.QHBoxLayout()
        izq.addLayout(fila, 2)
        self.cortes = []
        for titulo, eje in ((_("Corte x"), 0), (_("Corte y"), 1),
                            (_("Corte z"), 2)):
            col = QtWidgets.QVBoxLayout()
            v = _CorteCampo(titulo)
            sl = QtWidgets.QSlider(QtCore.Qt.Horizontal)
            sl.setRange(0, 100)
            sl.setValue(50)
            sl.valueChanged.connect(lambda _v, e=eje: self._pintar_corte(e))
            col.addWidget(v, 1)
            col.addWidget(sl)
            fila.addLayout(col, 1)
            self.cortes.append((eje, v, sl))

        # ---- derecha: mandos -------------------------------------------------
        cont = QtWidgets.QWidget()
        der = QtWidgets.QVBoxLayout(cont)
        der.setContentsMargins(0, 0, 4, 0)
        der.addWidget(self._grupo_estructura())
        der.addWidget(self._grupo_plantilla())
        self.g_cargas = self._grupo_cargas()
        der.addWidget(self.g_cargas)
        self.g_apoyos = self._grupo_apoyos()
        der.addWidget(self.g_apoyos)
        der.addWidget(self._grupo_material())
        der.addWidget(self._grupo_resultados())
        der.addWidget(self._grupo_barrido())
        der.addWidget(self._grupo_exportar())
        der.addStretch(1)
        sc = QtWidgets.QScrollArea()
        sc.setWidget(cont)
        sc.setWidgetResizable(True)
        sc.setFixedWidth(520)
        sc.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        raiz.addWidget(sc)

        for tecla, fn in (("Ctrl+Z", self.deshacer),
                          ("Ctrl+Shift+Z", self.rehacer),
                          ("Ctrl+Y", self.rehacer)):
            a = QtWidgets.QAction(self)
            a.setShortcut(tecla)
            a.triggered.connect(fn)
            self.addAction(a)

        self._preparar_estructura()
        if not self.caso["cargas"] and not self.caso["apoyos"]:
            self._aplicar_plantilla("compresion", registrar=False)
        self._snap(_("Estado inicial"), forzar=True)
        self._refrescar_listas()
        self._sel_carga(0 if self.caso["cargas"] else -1)

    # ------------------------------------------------------------------
    # Construcción de los grupos
    # ------------------------------------------------------------------

    def _barra_herramientas(self):
        h = QtWidgets.QHBoxLayout()
        self._botones_vista = []
        for clave, texto in (("iso", _("Perspectiva")), ("xy", _("Plano x-y")),
                             ("xz", _("Plano x-z")), ("yz", _("Plano y-z"))):
            b = QtWidgets.QPushButton(texto)
            b.clicked.connect(lambda _c=False, k=clave: self._vista(k))
            h.addWidget(b)
            self._botones_vista.append(b)
        h.addSpacing(10)
        h.addWidget(QtWidgets.QLabel(_("Elegir con el ratón:")))
        self.cmb_modo = QtWidgets.QComboBox()
        for k, t in (("orbitar", _("Solo girar la vista")),
                     ("cara", _("Cara del cubo")),
                     ("parche", _("Parche en la superficie")),
                     ("nodos", _("Nodos (Mayús quita)"))):
            self.cmb_modo.addItem(t, k)
        self.cmb_modo.setToolTip(_(
            "Con una herramienta activa, un clic (sin arrastrar) sobre el "
            "hueso o el cubo asigna la región a la carga o al apoyo activo."))
        self.cmb_modo.currentIndexChanged.connect(self._cambia_modo)
        h.addWidget(self.cmb_modo)
        self.sp_parche = _spin(0.001, 100, 0.2, 3, "mm", 0.05, 90)
        self.sp_parche.setToolTip(_("Radio del parche que crea un clic."))
        h.addWidget(self.sp_parche)
        self.btn_nueva_pt = QtWidgets.QPushButton(_("Nueva carga con un clic"))
        self.btn_nueva_pt.setToolTip(_(
            "El próximo clic sobre la superficie crea una fuerza normal "
            "entrante sobre un parche."))
        self.btn_nueva_pt.clicked.connect(self._nueva_con_clic)
        h.addWidget(self.btn_nueva_pt)
        h.addStretch(1)
        self.chk_tras = QtWidgets.QCheckBox(_("Traslúcido"))
        self.chk_tras.setChecked(True)
        self.chk_tras.toggled.connect(self._traslucido)
        h.addWidget(self.chk_tras)
        self.btn_atras = QtWidgets.QToolButton()
        self.btn_atras.setText("↶")
        self.btn_atras.setToolTip(_("Deshacer (Ctrl+Z)"))
        self.btn_atras.clicked.connect(self.deshacer)
        self.btn_adel = QtWidgets.QToolButton()
        self.btn_adel.setText("↷")
        self.btn_adel.setToolTip(_("Rehacer (Ctrl+Mayús+Z)"))
        self.btn_adel.clicked.connect(self.rehacer)
        h.addWidget(self.btn_atras)
        h.addWidget(self.btn_adel)
        return h

    def _grupo_estructura(self):
        g = QtWidgets.QGroupBox(_("Estructura"))
        v = QtWidgets.QVBoxLayout(g)
        self.cmb_estructura = QtWidgets.QComboBox()
        for k, d in self.estructuras.items():
            self.cmb_estructura.addItem(d["etq"], k)
        self.cmb_estructura.setCurrentIndex(
            self.cmb_estructura.findData(self.clave))
        self.cmb_estructura.currentIndexChanged.connect(self._cambia_estructura)
        v.addWidget(_fila(QtWidgets.QLabel(_("Aplicar a")),
                          self.cmb_estructura))
        self.sp_n = QtWidgets.QSpinBox()
        self.sp_n.setRange(8, 160)
        self.sp_n.setValue(int(max(8, min(160, self._n_defecto))))
        self.sp_n.setToolTip(_(
            "Vóxeles por lado con que se resuelve. Más es más fino y mucho "
            "más lento. Las regiones se dan en mm y no cambian al variarlo."))
        self.sp_n.editingFinished.connect(self._cambia_n)
        v.addWidget(_fila(QtWidgets.QLabel(_("Resolución")), self.sp_n,
                          QtWidgets.QLabel(_("vóxeles por lado"))))
        self.lab_info = QtWidgets.QLabel()
        self.lab_info.setWordWrap(True)
        self.lab_info.setStyleSheet("color:#555; font-size:11px;")
        v.addWidget(self.lab_info)
        return g

    @property
    def _n_defecto(self):
        return getattr(self, "_n0", 40)

    def _grupo_plantilla(self):
        g = QtWidgets.QGroupBox(_("Plantilla de partida"))
        v = QtWidgets.QVBoxLayout(g)
        self.cmb_plantilla = QtWidgets.QComboBox()
        for k, t in plantillas_txt():
            self.cmb_plantilla.addItem(t, k)
        self.cmb_eje_pl = QtWidgets.QComboBox()
        self.cmb_eje_pl.addItems(["x", "y", "z"])
        self.cmb_eje_pl.setCurrentIndex(2)
        self.sp_pl = _spin(0.0, 1e9, 1.0, 3, "", 0.5)
        v.addWidget(_fila(self.cmb_plantilla, QtWidgets.QLabel(_("eje")),
                          self.cmb_eje_pl))
        self.lab_pl = QtWidgets.QLabel()
        b = QtWidgets.QPushButton(_("Aplicar plantilla"))
        b.setToolTip(_(
            "Sustituye las cargas y los apoyos por los de la plantilla. Se "
            "puede deshacer."))
        b.clicked.connect(lambda: self._aplicar_plantilla(
            self.cmb_plantilla.currentData()))
        v.addWidget(_fila(self.lab_pl, self.sp_pl, b))
        self.cmb_plantilla.currentIndexChanged.connect(self._rotulo_plantilla)
        self._rotulo_plantilla()
        v.addWidget(_fila(self._boton(_("Cargar caso…"), self.cargar_caso),
                          self._boton(_("Guardar caso…"), self.guardar_caso)))
        return g

    def _boton(self, texto, fn, tip=None):
        b = QtWidgets.QPushButton(texto)
        b.clicked.connect(lambda _c=False: fn())
        if tip:
            b.setToolTip(tip)
        return b

    def _rotulo_plantilla(self, *_a):
        k = self.cmb_plantilla.currentData()
        self.lab_pl.setText(_("Tensión (MPa)") if k in ("compresion", "traccion")
                            else _("Fuerza (N)"))

    def _grupo_cargas(self):
        g = QtWidgets.QGroupBox(_("Cargas"))
        v = QtWidgets.QVBoxLayout(g)
        self.lista_c = QtWidgets.QListWidget()
        self.lista_c.setMaximumHeight(96)
        self.lista_c.currentRowChanged.connect(self._sel_carga)
        v.addWidget(self.lista_c)
        h = QtWidgets.QHBoxLayout()
        self.btn_add_c = QtWidgets.QToolButton()
        self.btn_add_c.setText(_("Añadir ▾"))
        self.btn_add_c.setPopupMode(QtWidgets.QToolButton.InstantPopup)
        m = QtWidgets.QMenu(self.btn_add_c)
        for k, t in tipos_carga_txt():
            m.addAction(t).triggered.connect(
                lambda _c=False, k=k: self._anadir_carga(k))
        self.btn_add_c.setMenu(m)
        h.addWidget(self.btn_add_c)
        h.addWidget(self._boton(_("Duplicar"), self._duplicar_carga))
        h.addWidget(self._boton(_("Quitar"), self._quitar_carga))
        h.addStretch(1)
        v.addLayout(h)

        self.f_c = QtWidgets.QWidget()
        f = QtWidgets.QVBoxLayout(self.f_c)
        f.setContentsMargins(0, 0, 0, 0)
        self.cmb_tipo_c = QtWidgets.QComboBox()
        for k, t in tipos_carga_txt():
            self.cmb_tipo_c.addItem(t, k)
        self.sp_mag = _spin(-1e12, 1e12, 100.0, 4, "", 10.0, 120)
        f.addWidget(_fila(self.cmb_tipo_c))
        f.addWidget(_fila(QtWidgets.QLabel(_("Magnitud")), self.sp_mag))

        # dirección (fuerza, presión, desplazamiento con magnitud)
        self.w_dir = QtWidgets.QWidget()
        d = QtWidgets.QVBoxLayout(self.w_dir)
        d.setContentsMargins(0, 0, 0, 0)
        self.cmb_dir = QtWidgets.QComboBox()
        for k, t in direcciones():
            self.cmb_dir.addItem(t, k)
        self.v_dir = [_spin(-1e6, 1e6, k == 2, 3, "", 0.1, 64)
                      for k in range(3)]
        self.sp_azi = _spin(-360, 360, 0, 1, "°", 5, 78)
        self.sp_ele = _spin(-90, 90, -90, 1, "°", 5, 78)
        d.addWidget(_fila(QtWidgets.QLabel(_("Dirección")), self.cmb_dir))
        d.addWidget(_fila(QtWidgets.QLabel(_("Vector")), *self.v_dir))
        d.addWidget(_fila(QtWidgets.QLabel(_("Azimut")), self.sp_azi,
                          QtWidgets.QLabel(_("Elevación")), self.sp_ele))
        self.chk_mango = QtWidgets.QCheckBox(_(
            "Arrastrar la punta de la flecha con el ratón"))
        self.chk_mango.setToolTip(_(
            "Muestra una esfera en la punta de la flecha de la carga activa: "
            "al arrastrarla se fija la dirección."))
        d.addWidget(self.chk_mango)
        f.addWidget(self.w_dir)

        # presión
        self.w_pres = QtWidgets.QWidget()
        self.cmb_area = QtWidgets.QComboBox()
        self.cmb_area.addItem(_("Sección bruta de la cara"), "bruta")
        self.cmb_area.addItem(_("Área de hueso expuesta"), "expuesta")
        self.cmb_area.setToolTip(_(
            "Sobre qué área se aplica la presión. Bruta: la sección del cubo, "
            "como la tensión aparente del ensayo de siempre."))
        h2 = QtWidgets.QHBoxLayout(self.w_pres)
        h2.setContentsMargins(0, 0, 0, 0)
        h2.addWidget(QtWidgets.QLabel(_("Área de referencia")))
        h2.addWidget(self.cmb_area)
        f.addWidget(self.w_pres)

        # momento
        self.w_mom = QtWidgets.QWidget()
        m2 = QtWidgets.QVBoxLayout(self.w_mom)
        m2.setContentsMargins(0, 0, 0, 0)
        self.v_eje = [_spin(-1e6, 1e6, k == 2, 3, "", 0.1, 64)
                      for k in range(3)]
        self.chk_centro_m = QtWidgets.QCheckBox(_("Centro del par"))
        self.chk_centro_m.setToolTip(_(
            "Por omisión, el baricentro de la región. Fíjalo para girar "
            "alrededor de otro punto."))
        self.v_cm = [_spin(-1e6, 1e6, 0, 3, "mm", 0.1, 64) for _k in range(3)]
        m2.addWidget(_fila(QtWidgets.QLabel(_("Eje")), *self.v_eje))
        m2.addWidget(_fila(self.chk_centro_m, *self.v_cm))
        f.addWidget(self.w_mom)

        # desplazamiento
        self.w_des = QtWidgets.QWidget()
        g2 = QtWidgets.QGridLayout(self.w_des)
        g2.setContentsMargins(0, 0, 0, 0)
        self.k_des = [QtWidgets.QCheckBox(e) for e in ("ux", "uy", "uz")]
        self.v_des = [_spin(-1e6, 1e6, 0.0, 5, "mm", 0.001, 110)
                      for _k in range(3)]
        for k in range(3):
            g2.addWidget(self.k_des[k], k, 0)
            g2.addWidget(self.v_des[k], k, 1)
        self.k_des[2].setChecked(True)
        self.v_des[2].setValue(-0.01)
        f.addWidget(self.w_des)

        self.reg_c = _EditorRegion(self._ext, self._paso)
        f.addWidget(self.reg_c)
        v.addWidget(self.f_c)

        for w in (self.cmb_tipo_c, self.cmb_dir, self.cmb_area):
            w.currentIndexChanged.connect(self._gesto_carga)
        for w in (self.sp_mag, *self.v_dir, *self.v_eje, *self.v_cm,
                  *self.v_des):
            w.valueChanged.connect(self._gesto_carga)
        self.sp_azi.valueChanged.connect(self._de_angulos)
        self.sp_ele.valueChanged.connect(self._de_angulos)
        for w in (self.chk_centro_m, *self.k_des):
            w.toggled.connect(self._gesto_carga)
        self.chk_mango.toggled.connect(self._mango)
        self.reg_c.cambiado.connect(self._gesto_carga)
        return g

    @property
    def _ext(self):
        return C.extension(self._BW.shape, self._sp) if self._BW is not None \
            else np.ones(3)

    @property
    def _paso(self):
        return float(np.min(self._sp)) if self._sp is not None else 0.01

    def _grupo_apoyos(self):
        g = QtWidgets.QGroupBox(_("Apoyos"))
        v = QtWidgets.QVBoxLayout(g)
        self.lista_a = QtWidgets.QListWidget()
        self.lista_a.setMaximumHeight(80)
        self.lista_a.currentRowChanged.connect(self._sel_apoyo)
        v.addWidget(self.lista_a)
        h = QtWidgets.QHBoxLayout()
        self.btn_add_a = QtWidgets.QToolButton()
        self.btn_add_a.setText(_("Añadir ▾"))
        self.btn_add_a.setPopupMode(QtWidgets.QToolButton.InstantPopup)
        m = QtWidgets.QMenu(self.btn_add_a)
        for k, t in tipos_apoyo_txt():
            m.addAction(t).triggered.connect(
                lambda _c=False, k=k: self._anadir_apoyo(k))
        self.btn_add_a.setMenu(m)
        h.addWidget(self.btn_add_a)
        h.addWidget(self._boton(_("Duplicar"), self._duplicar_apoyo))
        h.addWidget(self._boton(_("Quitar"), self._quitar_apoyo))
        h.addStretch(1)
        v.addLayout(h)
        self.f_a = QtWidgets.QWidget()
        f = QtWidgets.QVBoxLayout(self.f_a)
        f.setContentsMargins(0, 0, 0, 0)
        self.cmb_tipo_a = QtWidgets.QComboBox()
        for k, t in tipos_apoyo_txt():
            self.cmb_tipo_a.addItem(t, k)
        f.addWidget(self.cmb_tipo_a)
        self.k_gdl = [QtWidgets.QCheckBox(e) for e in ("x", "y", "z")]
        f.addWidget(_fila(QtWidgets.QLabel(_("Grados de libertad fijos")),
                          *self.k_gdl))
        self.reg_a = _EditorRegion(self._ext, self._paso)
        f.addWidget(self.reg_a)
        v.addWidget(self.f_a)
        self.cmb_tipo_a.currentIndexChanged.connect(self._gesto_apoyo)
        for k in self.k_gdl:
            k.toggled.connect(self._gesto_apoyo)
        self.reg_a.cambiado.connect(self._gesto_apoyo)
        return g

    def _grupo_material(self):
        g = QtWidgets.QGroupBox(_("Material y solver"))
        v = QtWidgets.QVBoxLayout(g)
        mat = self.caso["material"]
        self.sp_E = _spin(0.001, 1e4, mat["E_s"] / 1e9, 3, "GPa", 1.0, 110)
        self.sp_nu = _spin(-0.99, 0.499, mat["nu"], 3, "", 0.01, 80)
        v.addWidget(_fila(QtWidgets.QLabel(_("Tejido: E")), self.sp_E,
                          QtWidgets.QLabel("ν"), self.sp_nu))
        self.chk_estab = QtWidgets.QCheckBox(_("Estabilizar automáticamente"))
        self.chk_estab.setChecked(bool(self.caso.get("estabilizar", True)))
        self.chk_estab.setToolTip(_(
            "Si los apoyos no impiden el movimiento de sólido rígido, añade "
            "las restricciones mínimas y comprueba que sus reacciones son "
            "nulas (es decir, que las cargas están equilibradas)."))
        self.chk_portante = QtWidgets.QCheckBox(_(
            "Solo los fragmentos que unen cargas y apoyos"))
        self.chk_portante.setChecked(bool(self.caso.get("solo_portante", True)))
        v.addWidget(self.chk_estab)
        v.addWidget(self.chk_portante)
        for w in (self.sp_E, self.sp_nu):
            w.valueChanged.connect(self._gesto_material)
        for w in (self.chk_estab, self.chk_portante):
            w.toggled.connect(self._gesto_material)
        self.btn_resolver = QtWidgets.QPushButton(_("Resolver"))
        self.btn_resolver.setStyleSheet("font-weight:bold;")
        self.btn_resolver.clicked.connect(self.resolver)
        self.btn_comparar = QtWidgets.QPushButton(_("Resolver en todas las estructuras"))
        self.btn_comparar.setToolTip(_(
            "El mismo caso sobre cada estructura disponible, con una tabla "
            "comparativa."))
        self.btn_comparar.clicked.connect(self.comparar)
        v.addWidget(_fila(self.btn_resolver, self.btn_comparar))
        self.barra = QtWidgets.QProgressBar()
        self.barra.setRange(0, 100)
        self.barra.setVisible(False)
        v.addWidget(self.barra)
        self.lab_aviso = QtWidgets.QLabel()
        self.lab_aviso.setWordWrap(True)
        v.addWidget(self.lab_aviso)
        return g

    def _grupo_resultados(self):
        g = QtWidgets.QGroupBox(_("Resultados"))
        v = QtWidgets.QVBoxLayout(g)
        self.cmb_campo = QtWidgets.QComboBox()
        self.cmb_campo.addItem(_("(sin campo)"), None)
        for k in C.CAMPOS:
            self.cmb_campo.addItem(_(C.ETIQUETA_CAMPO[k]), k)
        self.cmb_campo.setCurrentIndex(1)
        v.addWidget(_fila(QtWidgets.QLabel(_("Campo")), self.cmb_campo))
        self.chk_def = QtWidgets.QCheckBox(_("Deformada ×"))
        self.sp_def = _spin(0.0, 1e9, 1.0, 1, "", 10.0, 90)
        self.btn_def_auto = self._boton(_("Auto"), self._def_auto)
        v.addWidget(_fila(self.chk_def, self.sp_def, self.btn_def_auto))
        self.chk_rango = QtWidgets.QCheckBox(_("Rango de color fijo"))
        self.sp_cmin = _spin(-1e12, 1e12, 0.0, 4, "", 1.0, 90)
        self.sp_cmax = _spin(-1e12, 1e12, 1.0, 4, "", 1.0, 90)
        v.addWidget(_fila(self.chk_rango, self.sp_cmin, self.sp_cmax))
        self.cmb_modo_c = QtWidgets.QComboBox()
        self.cmb_modo_c.addItem(_("Superficie suave"), "suave")
        self.cmb_modo_c.addItem(_("Caras de voxel (valor exacto)"), "voxeles")
        v.addWidget(_fila(QtWidgets.QLabel(_("Dibujo")), self.cmb_modo_c))
        self.chk_corte = QtWidgets.QCheckBox(_("Cortar en"))
        self.cmb_eje_corte = QtWidgets.QComboBox()
        self.cmb_eje_corte.addItems(["x", "y", "z"])
        self.cmb_eje_corte.setCurrentIndex(1)
        self.sl_corte = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.sl_corte.setRange(0, 100)
        self.sl_corte.setValue(50)
        v.addWidget(_fila(self.chk_corte, self.cmb_eje_corte, self.sl_corte))
        for w in (self.cmb_campo, self.cmb_modo_c, self.cmb_eje_corte):
            w.currentIndexChanged.connect(self._pintar_resultado)
        for w in (self.chk_def, self.chk_rango, self.chk_corte):
            w.toggled.connect(self._pintar_resultado)
        for w in (self.sp_def, self.sp_cmin, self.sp_cmax):
            w.valueChanged.connect(self._pintar_resultado)
        self.sl_corte.valueChanged.connect(self._pintar_resultado)

        self.tabla = QtWidgets.QTableWidget(0, 4)
        self.tabla.setHorizontalHeaderLabels(
            [_("Origen"), _("Fuerza (N)"), _("Módulo (N)"), _("Momento (N·mm)")])
        self.tabla.horizontalHeader().setStretchLastSection(True)
        self.tabla.setMaximumHeight(120)
        self.tabla.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        v.addWidget(self.tabla)
        self.lab_res = QtWidgets.QLabel(_("Todavía no se ha resuelto."))
        self.lab_res.setWordWrap(True)
        self.lab_res.setTextFormat(QtCore.Qt.RichText)
        v.addWidget(self.lab_res)
        return g

    def _grupo_barrido(self):
        g = QtWidgets.QGroupBox(_("Barrido de ángulo"))
        v = QtWidgets.QVBoxLayout(g)
        self.cmb_bar_carga = QtWidgets.QComboBox()
        self.cmb_bar_plano = QtWidgets.QComboBox()
        for t, k in ((_("plano x-y"), (0, 1)), (_("plano x-z"), (0, 2)),
                     (_("plano y-z"), (1, 2))):
            self.cmb_bar_plano.addItem(t, k)
        self.sp_bar_n = QtWidgets.QSpinBox()
        self.sp_bar_n.setRange(4, 180)
        self.sp_bar_n.setValue(24)
        v.addWidget(_fila(QtWidgets.QLabel(_("Girar")), self.cmb_bar_carga))
        v.addWidget(_fila(self.cmb_bar_plano, QtWidgets.QLabel(_("pasos")),
                          self.sp_bar_n))
        lab = QtWidgets.QLabel(_(
            "Con apoyo empotrado. Resuelve dos veces, no una por ángulo: el "
            "problema es lineal."))
        lab.setWordWrap(True)
        lab.setStyleSheet("color:#555; font-size:11px;")
        v.addWidget(lab)
        v.addWidget(self._boton(_("Barrer…"), self.barrer))
        return g

    def _grupo_exportar(self):
        g = QtWidgets.QGroupBox(_("Exportar"))
        v = QtWidgets.QVBoxLayout(g)
        v.addWidget(_fila(self._boton(_("Imagen 3D…"), self.exportar_imagen),
                          self._boton(_("VTK (ParaView)…"), self.exportar_vtk)))
        v.addWidget(_fila(self._boton(_("CSV…"), self.exportar_csv),
                          self._boton(_("JSON…"), self.exportar_json)))
        b = self._boton(_("Añadir al informe de publicación"), self.registrar)
        b.setToolTip(_(
            "Guarda este caso y su resultado en la sesión: el informe de "
            "publicación los declara en una sección propia."))
        v.addWidget(b)
        return g

    # ------------------------------------------------------------------
    # Estructura y escena
    # ------------------------------------------------------------------

    def _estructura_calculo(self, clave, n):
        d = self.estructuras[clave]
        # `remuestrear_bw` toma un escalar por un espaciado inválido y lo
        # cambia por 1 mm: se le pasa siempre el vector de tres.
        return remuestrear_bw(d["BW"], C._spacing3(d["spacing"]), n)

    def _preparar_estructura(self):
        from spinpy.escena_cargas import EscenaCargas
        n = int(self.sp_n.value()) if hasattr(self, "sp_n") else self._n_defecto
        self._BW, self._sp = self._estructura_calculo(self.clave, n)
        self._sp = C._spacing3(self._sp)
        self._cache_sel = {}
        # Las caras de voxel expuestas no dependen del 3D: con ellas se
        # resuelven los clics y los nodos aunque la vista no este disponible.
        self._caras = C.caras_expuestas(self._BW, self._sp)
        self.escena = None
        if self.vista3d is not None:
            try:
                self.vista3d.clear()
                self.escena = EscenaCargas(self.vista3d, self._BW, self._sp,
                                           caras=self._caras)
                self.escena.fijar_traslucido(self.chk_tras.isChecked()
                                             if hasattr(self, "chk_tras")
                                             else True)
                self.escena.vista("iso")
                if hasattr(self, "cmb_modo"):
                    self._cambia_modo()
            except Exception:
                self.escena = None
        ext = C.extension(self._BW.shape, self._sp)
        self.lab_info.setText(_(
            "{n}×{m}×{k} vóxeles, {dx:.3g} × {dy:.3g} × {dz:.3g} mm, "
            "BV/TV {r:.3f}").format(
                n=self._BW.shape[0], m=self._BW.shape[1], k=self._BW.shape[2],
                dx=ext[0], dy=ext[1], dz=ext[2], r=float(self._BW.mean())))
        self._remapear_nodos()
        for ed in ("reg_c", "reg_a"):
            if hasattr(self, ed):
                self._ajustar_editor(getattr(self, ed))
        self._pintar_caso()
        self._pintar_cortes()

    def _ajustar_editor(self, ed):
        ext = self._ext
        L = float(ext.max())
        spins = (*ed.e_centro, *ed.e_lo, *ed.e_hi, *ed.c_centro)
        for s in spins:
            s.setRange(-L, 2 * L)
            s.setSingleStep(self._paso)

    def _remapear_nodos(self):
        """Los nodos elegidos se guardan con sus coordenadas: al cambiar de
        resolución o de estructura se vuelven a buscar en la malla nueva."""
        for lista in (self.caso["cargas"], self.caso["apoyos"]):
            for it in lista:
                r = it["region"]
                if r.get("forma") != "nodos" or not r.get("xyz"):
                    continue
                ids = []
                for x in r["xyz"]:
                    i = C.nodo_cercano(self._caras, x, self._BW.shape,
                                       self._sp)
                    if i is not None and i not in ids:
                        ids.append(i)
                r["ids"] = ids

    def _cambia_estructura(self, *_a):
        k = self.cmb_estructura.currentData()
        if k == self.clave:
            return
        self.resultados[self.clave] = (self.res, copy.deepcopy(self.caso))
        self.clave = k
        self._preparar_estructura()
        prev = self.resultados.get(k)
        self.res = prev[0] if prev else None
        self._mostrar_res_o_limpiar()

    def _cambia_n(self):
        if self._BW is None:
            return
        n = int(self.sp_n.value())
        if max(self._BW.shape) == min(n, max(self.estructuras[self.clave][
                "BW"].shape)):
            return
        self.res = None
        self._preparar_estructura()
        self._mostrar_res_o_limpiar()

    def _vista(self, k):
        if self.escena is not None:
            self.escena.vista(k)

    def _traslucido(self, v):
        if self.escena is not None:
            self.escena.fijar_traslucido(v)

    # ------------------------------------------------------------------
    # Historial
    # ------------------------------------------------------------------

    def _snap(self, etiqueta="", forzar=False):
        """Guarda el estado del caso. Los cambios seguidos del MISMO mando
        (arrastrar un deslizador) se funden en un solo paso."""
        txt = json.dumps(C._json_seguro(self.caso), sort_keys=True)
        if self._historial and self._historial[self._pos] == txt:
            return
        ahora = time.time()
        if (not forzar and etiqueta and etiqueta == self._ultimo_snap[0]
                and ahora - self._ultimo_snap[1] < 0.8
                and self._pos == len(self._historial) - 1):
            self._historial[self._pos] = txt
        else:
            del self._historial[self._pos + 1:]
            self._historial.append(txt)
            if len(self._historial) > LIMITE_HISTORIAL:
                self._historial.pop(0)
            self._pos = len(self._historial) - 1
        self._ultimo_snap = (etiqueta, ahora)
        self._botones_historial()

    def _botones_historial(self):
        self.btn_atras.setEnabled(self._pos > 0)
        self.btn_adel.setEnabled(self._pos < len(self._historial) - 1)

    def _restaurar(self):
        self.caso = json.loads(self._historial[self._pos])
        self._sincronizar_material()
        self._remapear_nodos()
        self._refrescar_listas()
        self._pintar_caso()
        self._botones_historial()

    def deshacer(self):
        if self._pos > 0:
            self._pos -= 1
            self._restaurar()

    def rehacer(self):
        if self._pos < len(self._historial) - 1:
            self._pos += 1
            self._restaurar()

    # ------------------------------------------------------------------
    # Cargas y apoyos: listas y formularios
    # ------------------------------------------------------------------

    def _refrescar_listas(self, c=None, a=None):
        self._silencio = True
        try:
            for lista, items, f in ((self.lista_c, self.caso["cargas"],
                                     C.describir_carga),
                                    (self.lista_a, self.caso["apoyos"],
                                     C.describir_apoyo)):
                actual = lista.currentRow()
                lista.clear()
                for i, it in enumerate(items):
                    lista.addItem(f(it, i))
                lista.setCurrentRow(min(max(actual, 0), len(items) - 1)
                                    if items else -1)
            self.cmb_bar_carga.clear()
            for i, it in enumerate(self.caso["cargas"]):
                if it["tipo"] in ("fuerza", "presion", "momento"):
                    self.cmb_bar_carga.addItem(C.describir_carga(it, i), i)
        finally:
            self._silencio = False
        self.f_c.setVisible(bool(self.caso["cargas"]))
        self.f_a.setVisible(bool(self.caso["apoyos"]))

    def _sel_carga(self, i):
        if self._silencio:
            return
        self._foco = ("carga", i)
        if 0 <= i < len(self.caso["cargas"]):
            self._cargar_carga(self.caso["cargas"][i])
            self.lista_c.setCurrentRow(i)
        self._pintar_caso()
        self._mango(self.chk_mango.isChecked())

    def _sel_apoyo(self, i):
        if self._silencio:
            return
        self._foco = ("apoyo", i)
        if 0 <= i < len(self.caso["apoyos"]):
            self._cargar_apoyo(self.caso["apoyos"][i])
        self._pintar_caso()

    def _cargar_carga(self, c):
        self._silencio = True
        try:
            t = c["tipo"]
            self.cmb_tipo_c.setCurrentIndex(self.cmb_tipo_c.findData(t))
            self.sp_mag.setValue(float(c.get("magnitud", 0.0)))
            d = c.get("direccion", "normal_entrante" if t != "fuerza"
                      else "-z")
            if isinstance(d, str):
                self.cmb_dir.setCurrentIndex(self.cmb_dir.findData(d))
                v = C.direccion_vector(d)
                if v is not None:
                    for s, x in zip(self.v_dir, v):
                        s.setValue(float(x))
            else:
                self.cmb_dir.setCurrentIndex(self.cmb_dir.findData("vector"))
                for s, x in zip(self.v_dir, d):
                    s.setValue(float(x))
            self._angulos_desde_vector()
            self.cmb_area.setCurrentIndex(max(0, self.cmb_area.findData(
                c.get("area_ref", "bruta"))))
            for s, x in zip(self.v_eje, c.get("eje", [0, 0, 1])):
                s.setValue(float(x))
            cm = c.get("centro")
            self.chk_centro_m.setChecked(cm is not None)
            for s, x in zip(self.v_cm, cm or [0, 0, 0]):
                s.setValue(float(x))
            u = c.get("u")
            if u is None and t == "desplazamiento":
                dv = C.direccion_vector(c.get("direccion", "+z"))
                m = float(c.get("magnitud", 0.0))
                u = ([None] * 3 if dv is None else
                     [m * dv[k] if abs(dv[k]) > 1e-12 else None
                      for k in range(3)])
            for k in range(3):
                x = None if u is None else u[k]
                self.k_des[k].setChecked(x is not None)
                self.v_des[k].setValue(0.0 if x is None else float(x))
            self.reg_c.poner(c["region"])
            self._ver_tipo(t)
        finally:
            self._silencio = False

    def _ver_tipo(self, t):
        self.w_dir.setVisible(t in ("fuerza", "presion"))
        self.w_pres.setVisible(t == "presion")
        self.w_mom.setVisible(t == "momento")
        self.w_des.setVisible(t == "desplazamiento")
        self.sp_mag.setVisible(t != "desplazamiento")
        self.sp_mag.setSuffix({"fuerza": " N", "presion": " MPa",
                               "momento": " N·mm"}.get(t, ""))

    def _leer_carga(self):
        t = self.cmb_tipo_c.currentData()
        c = {"tipo": t, "region": self.reg_c.region()}
        if t in ("fuerza", "presion", "momento"):
            c["magnitud"] = self.sp_mag.value()
        if t in ("fuerza", "presion"):
            m = self.cmb_dir.currentData()
            c["direccion"] = ([s.value() for s in self.v_dir]
                              if m == "vector" else m)
            if t == "presion":
                c["area_ref"] = self.cmb_area.currentData()
        if t == "momento":
            c["eje"] = [s.value() for s in self.v_eje]
            if self.chk_centro_m.isChecked():
                c["centro"] = [s.value() for s in self.v_cm]
        if t == "desplazamiento":
            c["u"] = [self.v_des[k].value() if self.k_des[k].isChecked()
                      else None for k in range(3)]
        return c

    def _gesto_carga(self, *_a):
        if self._silencio:
            return
        i = self.lista_c.currentRow()
        if not 0 <= i < len(self.caso["cargas"]):
            return
        viejo = self.caso["cargas"][i]
        nuevo = self._leer_carga()
        if nuevo["tipo"] != viejo["tipo"]:
            self._ver_tipo(nuevo["tipo"])
            nuevo = self._con_defecto_tipo(nuevo)
        sel_dir = self.cmb_dir.currentData()
        if sel_dir not in (None, "vector") and nuevo["tipo"] in (
                "fuerza", "presion"):
            v = C.direccion_vector(sel_dir)
            if v is not None:
                self._silencio = True
                for s, x in zip(self.v_dir, v):
                    s.setValue(float(x))
                self._silencio = False
        self.caso["cargas"][i] = nuevo
        self._angulos_desde_vector()
        self._tras_cambio(_("carga"), carga=i)

    def _con_defecto_tipo(self, c):
        t = c["tipo"]
        c.setdefault("region", {"forma": "cara", "cara": "+z"})
        if t == "fuerza":
            c.setdefault("magnitud", 100.0)
            c["direccion"] = c.get("direccion", "normal_entrante")
        elif t == "presion":
            c["magnitud"] = 1.0
            c["direccion"] = "normal_entrante"
            c["area_ref"] = "bruta"
        elif t == "momento":
            c["magnitud"] = 100.0
            c["eje"] = [0, 0, 1]
        else:
            c["u"] = [None, None, -0.01]
        self._cargar_carga(c)
        return c

    def _de_angulos(self, *_a):
        if self._silencio:
            return
        a, e = np.radians(self.sp_azi.value()), np.radians(self.sp_ele.value())
        v = [np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)]
        self._silencio = True
        self.cmb_dir.setCurrentIndex(self.cmb_dir.findData("vector"))
        for s, x in zip(self.v_dir, v):
            s.setValue(float(x))
        self._silencio = False
        self._gesto_carga()

    def _angulos_desde_vector(self):
        v = np.array([s.value() for s in self.v_dir])
        n = np.linalg.norm(v)
        if n < 1e-12:
            return
        v = v / n
        bloqueado = self._silencio
        self._silencio = True
        self.sp_azi.setValue(float(np.degrees(np.arctan2(v[1], v[0]))))
        self.sp_ele.setValue(float(np.degrees(np.arcsin(np.clip(v[2], -1, 1)))))
        self._silencio = bloqueado

    def _cargar_apoyo(self, a):
        self._silencio = True
        try:
            self.cmb_tipo_a.setCurrentIndex(self.cmb_tipo_a.findData(a["tipo"]))
            g = C._gdl(a.get("gdl")) or {"empotrado": [0, 1, 2],
                                         "deslizante": [], "personalizado": []
                                         }[a["tipo"]]
            if a["tipo"] == "deslizante" and not g:
                r = a["region"]
                g = ([C.CARAS[r["cara"]][0]] if r.get("forma") == "cara"
                     else [2])
            for k in range(3):
                self.k_gdl[k].setChecked(k in g)
            self.reg_a.poner(a["region"])
        finally:
            self._silencio = False

    def _gesto_apoyo(self, *_a):
        if self._silencio:
            return
        i = self.lista_a.currentRow()
        if not 0 <= i < len(self.caso["apoyos"]):
            return
        t = self.cmb_tipo_a.currentData()
        a = {"tipo": t, "region": self.reg_a.region()}
        if t == "personalizado":
            a["gdl"] = [e for k, e in enumerate("xyz")
                        if self.k_gdl[k].isChecked()]
        elif t == "empotrado":
            self._silencio = True
            for k in self.k_gdl:
                k.setChecked(True)
            self._silencio = False
        if t == "deslizante":
            r = a["region"]
            g = ([C.CARAS[r["cara"]][0]] if r.get("forma") == "cara" else None)
            if g:
                self._silencio = True
                for k in range(3):
                    self.k_gdl[k].setChecked(k in g)
                self._silencio = False
        for k in self.k_gdl:
            k.setEnabled(t == "personalizado")
        self.caso["apoyos"][i] = a
        self._tras_cambio(_("apoyo"), apoyo=i)

    def _tras_cambio(self, etiqueta, carga=None, apoyo=None):
        self._snap(etiqueta)
        self._refrescar_item(carga, apoyo)
        self._pintar_caso()
        self._invalidar_resultado()
        if self.chk_mango.isChecked():
            self._mango(True)

    def _refrescar_item(self, carga, apoyo):
        self._silencio = True
        try:
            if carga is not None and carga < self.lista_c.count():
                self.lista_c.item(carga).setText(
                    C.describir_carga(self.caso["cargas"][carga], carga))
            if apoyo is not None and apoyo < self.lista_a.count():
                self.lista_a.item(apoyo).setText(
                    C.describir_apoyo(self.caso["apoyos"][apoyo], apoyo))
        finally:
            self._silencio = False

    # añadir / duplicar / quitar ------------------------------------------------

    def _nueva_carga(self, tipo):
        reg = {"forma": "cara", "cara": "+z"}
        return {
            "fuerza": {"tipo": tipo, "region": reg, "magnitud": 100.0,
                       "direccion": "normal_entrante"},
            "presion": {"tipo": tipo, "region": reg, "magnitud": 1.0,
                        "direccion": "normal_entrante", "area_ref": "bruta"},
            "momento": {"tipo": tipo, "region": reg, "magnitud": 100.0,
                        "eje": [0, 0, 1]},
            "desplazamiento": {"tipo": tipo, "region": reg,
                               "u": [None, None, -0.01]}}[tipo]

    def _anadir_carga(self, tipo):
        self.caso["cargas"].append(self._nueva_carga(tipo))
        self._snap(_("añadir carga"), forzar=True)
        self._refrescar_listas()
        self.lista_c.setCurrentRow(len(self.caso["cargas"]) - 1)
        self._invalidar_resultado()

    def _duplicar_carga(self):
        i = self.lista_c.currentRow()
        if 0 <= i < len(self.caso["cargas"]):
            self.caso["cargas"].insert(i + 1, copy.deepcopy(
                self.caso["cargas"][i]))
            self._snap(_("duplicar carga"), forzar=True)
            self._refrescar_listas()
            self.lista_c.setCurrentRow(i + 1)
            self._invalidar_resultado()

    def _quitar_carga(self):
        i = self.lista_c.currentRow()
        if 0 <= i < len(self.caso["cargas"]):
            del self.caso["cargas"][i]
            self._snap(_("quitar carga"), forzar=True)
            self._refrescar_listas()
            self._sel_carga(self.lista_c.currentRow())
            self._invalidar_resultado()

    def _anadir_apoyo(self, tipo):
        a = {"tipo": tipo, "region": {"forma": "cara", "cara": "-z"}}
        if tipo == "personalizado":
            a["gdl"] = ["z"]
        self.caso["apoyos"].append(a)
        self._snap(_("añadir apoyo"), forzar=True)
        self._refrescar_listas()
        self.lista_a.setCurrentRow(len(self.caso["apoyos"]) - 1)
        self._invalidar_resultado()

    def _duplicar_apoyo(self):
        i = self.lista_a.currentRow()
        if 0 <= i < len(self.caso["apoyos"]):
            self.caso["apoyos"].insert(i + 1, copy.deepcopy(
                self.caso["apoyos"][i]))
            self._snap(_("duplicar apoyo"), forzar=True)
            self._refrescar_listas()
            self.lista_a.setCurrentRow(i + 1)
            self._invalidar_resultado()

    def _quitar_apoyo(self):
        i = self.lista_a.currentRow()
        if 0 <= i < len(self.caso["apoyos"]):
            del self.caso["apoyos"][i]
            self._snap(_("quitar apoyo"), forzar=True)
            self._refrescar_listas()
            self._sel_apoyo(self.lista_a.currentRow())
            self._invalidar_resultado()

    def _aplicar_plantilla(self, nombre, registrar=True):
        eje = self.cmb_eje_pl.currentIndex() if hasattr(self, "cmb_eje_pl") else 2
        val = self.sp_pl.value() if hasattr(self, "sp_pl") else 1.0
        kw = ({"sigma_MPa": val} if nombre in ("compresion", "traccion")
              else {"F_N": val})
        mat = self.caso["material"]
        est, por = self.caso.get("estabilizar", True), self.caso.get(
            "solo_portante", True)
        nuevo = C.plantilla(nombre, self._BW.shape, self._sp, eje=eje, **kw)
        nuevo["material"] = mat
        nuevo["estabilizar"], nuevo["solo_portante"] = est, por
        self.caso = nuevo
        if registrar:
            self._snap(_("plantilla"), forzar=True)
            self._refrescar_listas()
            self._sel_carga(0)
            self._pintar_caso()
            self._invalidar_resultado()

    # material -----------------------------------------------------------------

    def _sincronizar_material(self):
        m = self.caso["material"]
        self._silencio = True
        try:
            self.sp_E.setValue(float(m.get("E_s", C.E_S_DEF)) / 1e9)
            self.sp_nu.setValue(float(m.get("nu", C.NU_S_DEF)))
            self.chk_estab.setChecked(bool(self.caso.get("estabilizar", True)))
            self.chk_portante.setChecked(bool(self.caso.get("solo_portante",
                                                            True)))
        finally:
            self._silencio = False

    def _gesto_material(self, *_a):
        if self._silencio:
            return
        self.caso["material"] = {"E_s": self.sp_E.value() * 1e9,
                                 "nu": self.sp_nu.value()}
        self.caso["estabilizar"] = self.chk_estab.isChecked()
        self.caso["solo_portante"] = self.chk_portante.isChecked()
        self._snap(_("material"))
        self._invalidar_resultado()

    # ------------------------------------------------------------------
    # Ratón
    # ------------------------------------------------------------------

    def _cambia_modo(self, *_a):
        self._modo = self.cmb_modo.currentData()
        if self.escena is None:
            return
        if self._modo == "orbitar" and not self._nueva:
            self.escena.desactivar_pinchazo()
        else:
            self.escena.activar_pinchazo(self._al_pinchar)

    def _nueva_con_clic(self):
        self._nueva = True
        if self.escena is not None:
            self.escena.activar_pinchazo(self._al_pinchar)
        self.lab_estado.setText(_(
            "Haz clic sobre la superficie del hueso para crear la carga."))

    def _al_pinchar(self, punto, shift=False):
        """Un clic sobre la escena: lo convierte en región (ver cabecera)."""
        caras = self._caras
        shp = self._BW.shape
        if self._nueva:
            self._nueva = False
            try:
                reg = C.region_parche(punto, self.sp_parche.value(), caras)
            except C.ErrorCaso as e:
                self.lab_estado.setText(str(e))
                return
            c = self._nueva_carga("fuerza")
            c["region"] = reg
            self.caso["cargas"].append(c)
            self._snap(_("carga con un clic"), forzar=True)
            self._refrescar_listas()
            self.lista_c.setCurrentRow(len(self.caso["cargas"]) - 1)
            self._invalidar_resultado()
            if self._modo == "orbitar" and self.escena is not None:
                self.escena.desactivar_pinchazo()
            self.lab_estado.setText(_("Carga creada."))
            return
        tipo, i = self._foco
        lista = self.caso["cargas" if tipo == "carga" else "apoyos"]
        if not 0 <= i < len(lista):
            self.lab_estado.setText(_(
                "Elige o crea primero una carga o un apoyo."))
            return
        reg = lista[i]["region"]
        amb = reg.get("ambito")
        if self._modo == "cara":
            nombre, _d = C.cara_del_cubo(punto, shp, self._sp)
            nuevo = {"forma": "cara", "cara": nombre}
        elif self._modo == "parche":
            try:
                nuevo = C.region_parche(punto, self.sp_parche.value(), caras)
            except C.ErrorCaso as e:
                self.lab_estado.setText(str(e))
                return
        elif self._modo == "nodos":
            nid = C.nodo_cercano(caras, punto, shp, self._sp)
            if nid is None:
                return
            xyz = C._coord_nodos(np.array([nid]), shp, self._sp)[0].tolist()
            nuevo = {"forma": "nodos",
                     "ids": list(reg.get("ids", [])) if reg.get(
                         "forma") == "nodos" else [],
                     "xyz": [list(x) for x in reg.get("xyz", [])]
                     if reg.get("forma") == "nodos" else []}
            if shift:
                if nid in nuevo["ids"]:
                    k = nuevo["ids"].index(nid)
                    del nuevo["ids"][k]
                    if k < len(nuevo["xyz"]):
                        del nuevo["xyz"][k]
            elif nid not in nuevo["ids"]:
                nuevo["ids"].append(nid)
                nuevo["xyz"].append(xyz)
        else:
            return
        if amb and nuevo.get("forma") != "nodos":
            nuevo["ambito"] = amb
        lista[i]["region"] = nuevo
        if tipo == "carga":
            self._silencio = True
            self.reg_c.poner(nuevo)
            self._silencio = False
        else:
            self._silencio = True
            self.reg_a.poner(nuevo)
            self._silencio = False
        self._tras_cambio(_("región con el ratón"),
                          carga=i if tipo == "carga" else None,
                          apoyo=i if tipo == "apoyo" else None)
        self.lab_estado.setText(_("Región asignada: {r}").format(
            r=C.describir_region(nuevo)))

    def _mango(self, activo):
        if self.escena is None:
            return
        self.escena.quitar_mango()
        i = self.lista_c.currentRow()
        if not activo or not 0 <= i < len(self.caso["cargas"]):
            return
        c = self.caso["cargas"][i]
        if c["tipo"] not in ("fuerza", "presion", "momento"):
            return
        fl = self.escena.flechas.get(i)
        if fl is None:
            return
        o, d = fl
        L = 0.22 * float(np.linalg.norm(self._ext)) * 1.35
        # la esfera va en el extremo que apunta la carga
        if c["tipo"] == "momento":
            origen, largo = o, 0.45 * L
        else:
            origen, largo = o - d * L, L
        self.escena.activar_mango(origen, d, largo, self._mango_movido)

    def _mango_movido(self, direccion, _largo):
        i = self.lista_c.currentRow()
        if not 0 <= i < len(self.caso["cargas"]):
            return
        c = self.caso["cargas"][i]
        clave = "eje" if c["tipo"] == "momento" else "direccion"
        if c["tipo"] == "momento":
            c["eje"] = [float(x) for x in direccion]
        else:
            c["direccion"] = [float(x) for x in direccion]
        self._silencio = True
        for s, x in zip(self.v_eje if clave == "eje" else self.v_dir,
                        direccion):
            s.setValue(float(x))
        if clave == "direccion":
            self.cmb_dir.setCurrentIndex(self.cmb_dir.findData("vector"))
            self._angulos_desde_vector()
        self._silencio = False
        self._snap(_("flecha"))
        self._refrescar_item(i, None)
        if self.escena is not None:
            self.escena.poner_caso(self.caso, activa=i, render=True)
        self._invalidar_resultado()

    # ------------------------------------------------------------------
    # Dibujo
    # ------------------------------------------------------------------

    def _pintar_caso(self):
        if self.escena is None:
            return
        tipo, i = self._foco
        try:
            self.escena.poner_caso(
                self.caso, activa=i if tipo == "carga" else None,
                activo_apoyo=i if tipo == "apoyo" else None)
        except C.ErrorCaso as e:
            self.lab_estado.setText(str(e))

    def _imagen_corte(self, eje, frac):
        """RGB del corte `eje` y el pie con su posición; campo si lo hay."""
        BW = self._BW
        k = int(round(frac * (BW.shape[eje] - 1)))
        ejes = [a for a in range(3) if a != eje]
        sl = [slice(None)] * 3
        sl[eje] = k
        hueso = BW[tuple(sl)]
        nm = hueso.T                              # (fila = 2.o eje)
        campo = self._campo_actual()
        rgb = np.empty(nm.shape + (3,), np.uint8)
        rgb[:] = _CorteCampo.FONDO
        rango = None
        if campo is not None:
            global _LUT
            if _LUT is None:
                _LUT = _lut_turbo()
            vals = np.asarray(campo[tuple(sl)], float).T
            fin = np.isfinite(vals)
            rango = self._rango_color()
            lo, hi = rango
            t = np.clip(np.nan_to_num((vals - lo) / ((hi - lo) or 1.0)), 0, 1)
            col = _LUT[(t * 255).astype(int)]
            rgb[nm] = _CorteCampo.NODATO
            rgb[fin] = col[fin]
        else:
            rgb[nm] = _CorteCampo.HUESO
        pos = (k + 0.5) * self._sp[eje]
        return rgb, "= %.3g mm" % pos, rango, ejes

    def _pintar_corte(self, eje):
        for e, v, sl in self.cortes:
            if e != eje:
                continue
            rgb, pie, rango, _e = self._imagen_corte(eje, sl.value() / 100.0)
            v.poner(rgb, pie, rango)

    def _pintar_cortes(self):
        for e, _v, _s in self.cortes:
            self._pintar_corte(e)

    def _campo_actual(self):
        k = self.cmb_campo.currentData()
        if self.res is None or not self.res.get("ok") or k is None:
            return None
        return self.res["campos"].get(k)

    def _rango_color(self):
        if self.chk_rango.isChecked():
            return (self.sp_cmin.value(), self.sp_cmax.value())
        st = self.res["stats"][self.cmb_campo.currentData()]
        lo, hi = st["min"], st["max"]
        return (lo, hi if hi > lo else lo + 1e-12)

    def _def_auto(self):
        if self.res is not None and self.res.get("ok"):
            self.sp_def.setValue(C.escala_deformada(self.res))

    def _pintar_resultado(self, *_a):
        if self._silencio:
            return
        if self.escena is None:
            self._pintar_cortes()
            return
        k = self.cmb_campo.currentData()
        if self.res is None or not self.res.get("ok") or k is None:
            self.escena.quitar_resultado()
            self._pintar_caso()
            self._pintar_cortes()
            return
        self.escena.modo_campo = self.cmb_modo_c.currentData()
        corte = None
        if self.chk_corte.isChecked():
            ax = self.cmb_eje_corte.currentIndex()
            corte = (ax, self.sl_corte.value() / 100.0 * float(self._ext[ax]),
                     1)
        clim = ((self.sp_cmin.value(), self.sp_cmax.value())
                if self.chk_rango.isChecked() else None)
        self.escena.mostrar_resultado(
            self.res, k, clim=clim, deformada=self.chk_def.isChecked(),
            escala=self.sp_def.value(), corte=corte)
        self._pintar_cortes()

    # ------------------------------------------------------------------
    # Resolver
    # ------------------------------------------------------------------

    def _invalidar_resultado(self):
        if self.res is not None:
            self.res = None
            self.lab_res.setText(_(
                "El caso ha cambiado: resuelve de nuevo para ver el campo."))
            self.tabla.setRowCount(0)
            if self.escena is not None:
                self.escena.quitar_resultado()
                self._pintar_caso()
            self._pintar_cortes()

    def _ocupado(self, si):
        for w in (self.btn_resolver, self.btn_comparar):
            w.setEnabled(not si)
        self.barra.setVisible(si)
        self.barra.setValue(0)

    def resolver(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        try:
            C.validar_caso(self.caso)
        except C.ErrorCaso as e:
            self.lab_aviso.setText("<span style='color:#b00020'>%s</span>" % e)
            return
        BW, sp, caso = self._BW, self._sp, copy.deepcopy(self.caso)
        self._ocupado(True)
        self.lab_aviso.setText("")
        self._t0 = time.time()
        self.hilo = _Hilo(lambda inf: C.resolver_caso(
            BW, sp, caso, progreso=inf))
        self.hilo.avance.connect(lambda p, m: (self.barra.setValue(p),
                                               self.barra.setFormat(m)))
        self.hilo.listo.connect(self._listo)
        self.hilo.fallo.connect(self._fallo)
        self.hilo.start()

    def _fallo(self, tb):
        self._ocupado(False)
        self.lab_aviso.setText("<span style='color:#b00020'>%s</span>" % (
            tb.strip().splitlines()[-1] if tb.strip() else "error"))

    def _listo(self, res):
        self._ocupado(False)
        if not res.get("ok"):
            self.res = None
            self.lab_aviso.setText(
                "<span style='color:#b00020'>%s</span>" % res.get("msg", ""))
            return
        self.res = res
        self.resultados[self.clave] = (res, copy.deepcopy(self.caso))
        self.sp_def.setValue(C.escala_deformada(res))
        k = self.cmb_campo.currentData()
        if k is not None:
            st = res["stats"][k]
            self._silencio = True
            self.sp_cmin.setValue(st["min"])
            self.sp_cmax.setValue(st["max"])
            self._silencio = False
        av = res.get("avisos") or []
        self.lab_aviso.setText(
            "<span style='color:#9a6700'>%s</span>" % " ".join(av)
            if av else "")
        self._rellenar_resultados()
        self._pintar_resultado()

    def _mostrar_res_o_limpiar(self):
        if self.res is not None and self.res.get("ok"):
            self._rellenar_resultados()
            self._pintar_resultado()
        else:
            self.tabla.setRowCount(0)
            self.lab_res.setText(_("Todavía no se ha resuelto."))
            if self.escena is not None:
                self.escena.quitar_resultado()
            self._pintar_cortes()

    def _rellenar_resultados(self):
        r = self.res
        self.tabla.setRowCount(0)
        for x in r.get("reacciones", []):
            n = x["indice"]
            etq = ("A%d" % (n + 1) if x["origen"] == "apoyo"
                   else "C%d" % (n + 1))
            k = self.tabla.rowCount()
            self.tabla.insertRow(k)
            vals = (etq,
                    "(%s)" % ", ".join("%.4g" % a for a in x["F_N"]),
                    "%.4g" % x["modulo_N"],
                    "(%s)" % ", ".join("%.4g" % a for a in x["M_Nmm"]))
            for j, t in enumerate(vals):
                self.tabla.setItem(k, j, QtWidgets.QTableWidgetItem(t))
        eq = r["equilibrio"]
        st = r["stats"]
        partes = [
            _("<b>Equilibrio</b>: residuo de fuerzas {f:.1e}, de momentos "
              "{m:.1e}.").format(f=eq["residuo_F"], m=eq["residuo_M"]),
            _("<b>von Mises</b>: máx. {a:.4g} MPa, p99 {b:.4g} MPa.").format(
                a=st["vm"]["max"], b=st["vm"]["p99"]) if "vm" in st else "",
            _("{n} elementos, {d} grados de libertad, {s}.").format(
                n=r["n_elem"], d=r["n_dof"], s=r.get("solver", ""))]
        for g in r.get("rigidez", []):
            if "rigidez_N_mm" in g and np.isfinite(g["rigidez_N_mm"]):
                partes.append(_("<b>Rigidez de C{n}</b>: {k:.4g} N/mm.").format(
                    n=g["indice"] + 1, k=g["rigidez_N_mm"]))
        p = r.get("pistoia")
        if p:
            partes.append(_(
                "<b>Pistoia</b>: las cargas llegarían al criterio con un "
                "factor ×{f:.3g}. {nota}").format(f=p["factor"], nota=p["nota"]))
        self.lab_res.setText("<br>".join(t for t in partes if t))

    # ------------------------------------------------------------------
    # Comparar estructuras
    # ------------------------------------------------------------------

    def comparar(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        n = int(self.sp_n.value())
        caso = copy.deepcopy(self.caso)
        ests = {k: (d["etq"], d["BW"], d["spacing"])
                for k, d in self.estructuras.items()}

        def tarea(inf):
            filas = []
            for j, (k, (etq, BW, sp)) in enumerate(ests.items()):
                bw, spr = remuestrear_bw(BW, C._spacing3(sp), n)
                c = copy.deepcopy(caso)
                _remapear_en(c, bw, spr)
                r = C.resolver_caso(
                    bw, spr, c, campos=("vm",),
                    progreso=lambda f, m, j=j: inf((j + f) / len(ests), m))
                filas.append((k, etq, r))
            return filas

        self._ocupado(True)
        self.hilo = _Hilo(tarea)
        self.hilo.avance.connect(lambda p, m: (self.barra.setValue(p),
                                               self.barra.setFormat(m)))
        self.hilo.listo.connect(self._comparado)
        self.hilo.fallo.connect(self._fallo)
        self.hilo.start()

    def _comparado(self, filas):
        self._ocupado(False)
        d = QtWidgets.QDialog(self)
        d.setWindowTitle(_("Comparación de estructuras"))
        v = QtWidgets.QVBoxLayout(d)
        t = QtWidgets.QTableWidget(len(filas), 7)
        t.setHorizontalHeaderLabels([
            _("Estructura"), _("BV/TV"), _("Rigidez (N/mm)"),
            _("von Mises máx. (MPa)"), _("von Mises p99 (MPa)"),
            _("Factor de Pistoia"), _("Equilibrio")])
        for i, (k, etq, r) in enumerate(filas):
            bw = self.estructuras[k]["BW"]
            if r.get("ok"):
                rig = next((g["rigidez_N_mm"] for g in r.get("rigidez", [])
                            if "rigidez_N_mm" in g), float("nan"))
                vals = (etq, "%.3f" % float(bw.mean()), "%.4g" % rig,
                        "%.4g" % r["stats"]["vm"]["max"],
                        "%.4g" % r["stats"]["vm"]["p99"],
                        "%.3g" % (r["pistoia"]["factor"] if r.get("pistoia")
                                  else float("nan")),
                        "%.1e" % r["equilibrio"]["residuo_F"])
            else:
                vals = (etq, "%.3f" % float(bw.mean()), r.get("msg", ""),
                        "", "", "", "")
            for j, x in enumerate(vals):
                t.setItem(i, j, QtWidgets.QTableWidgetItem(x))
        t.resizeColumnsToContents()
        v.addWidget(t)
        v.addWidget(QtWidgets.QLabel(_(
            "El mismo caso de carga sobre cada estructura, con las regiones "
            "dadas en mm. Las diferencias son de estructura, no de carga.")))
        d.resize(900, 220)
        d.exec_()

    # ------------------------------------------------------------------
    # Barrido
    # ------------------------------------------------------------------

    def barrer(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        i = self.cmb_bar_carga.currentData()
        if i is None:
            QtWidgets.QMessageBox.information(
                self, _("Nada que girar"),
                _("El caso no tiene ninguna fuerza, presión ni momento."))
            return
        plano = self.cmb_bar_plano.currentData()
        n = int(self.sp_bar_n.value())
        BW, sp, caso = self._BW, self._sp, copy.deepcopy(self.caso)
        self._ocupado(True)
        self.hilo = _Hilo(lambda inf: C.barrido_direccion(
            BW, sp, caso, indice=i, plano=plano, n=n, progreso=inf))
        self.hilo.avance.connect(lambda p, m: (self.barra.setValue(p),
                                               self.barra.setFormat(m)))
        self.hilo.listo.connect(self._barrido_listo)
        self.hilo.fallo.connect(self._fallo)
        self.hilo.start()

    def _barrido_listo(self, b):
        self._ocupado(False)
        if not b.get("ok"):
            QtWidgets.QMessageBox.warning(self, _("Barrido"), b.get("msg", ""))
            return
        self.barrido = b
        d = QtWidgets.QDialog(self)
        d.setWindowTitle(_("Barrido de ángulo"))
        v = QtWidgets.QVBoxLayout(d)
        cmb = QtWidgets.QComboBox()
        opciones = [("rigidez_N_mm", _("Rigidez (N/mm)")),
                    ("vm_max", _("von Mises máximo (MPa)")),
                    ("vm_p99", _("von Mises p99 (MPa)")),
                    ("factor", _("Factor de Pistoia"))]
        if b["tipo"] == "momento":
            opciones = opciones[1:]
        for k, t in opciones:
            cmb.addItem(t, k)
        rosa = _Rosa()

        def pinta(*_a):
            k = cmb.currentData()
            vals = np.where(np.isfinite(b[k]), b[k], np.nan)
            rosa.poner(b["theta_deg"], vals, cmb.currentText())
        cmb.currentIndexChanged.connect(pinta)
        v.addWidget(cmb)
        v.addWidget(rosa, 1)
        v.addWidget(QtWidgets.QLabel(_(
            "Ángulo medido desde el primer eje del plano hacia el segundo.")))
        v.addWidget(self._boton(_("Guardar CSV…"),
                                lambda: self._guardar_barrido(b)))
        pinta()
        d.resize(420, 520)
        d.exec_()

    def _guardar_barrido(self, b):
        f, _x = QtWidgets.QFileDialog.getSaveFileName(
            self, _("Guardar el barrido"), str(self.carpeta / "barrido.csv"),
            "CSV (*.csv)")
        if not f:
            return
        import csv
        with open(f, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            cols = ("theta_deg", "rigidez_N_mm", "u_dir_mm", "vm_max",
                    "vm_p99", "factor")
            w.writerow(cols)
            for i in range(len(b["theta_deg"])):
                w.writerow([b[c][i] for c in cols])

    # ------------------------------------------------------------------
    # Archivos
    # ------------------------------------------------------------------

    def _ruta_guardar(self, titulo, nombre, filtro):
        f, _x = QtWidgets.QFileDialog.getSaveFileName(
            self, titulo, str(self.carpeta / nombre), filtro)
        return Path(f) if f else None

    def _necesita_resultado(self):
        if self.res is None or not self.res.get("ok"):
            QtWidgets.QMessageBox.information(
                self, _("Sin resultado"), _("Resuelve el caso primero."))
            return False
        return True

    def _estructura_info(self):
        return {"clave": self.clave,
                "etiqueta": self.estructuras[self.clave]["etq"],
                "forma": [int(x) for x in self._BW.shape],
                "spacing_mm": [float(x) for x in self._sp],
                "BV_TV": float(self._BW.mean())}

    def exportar_imagen(self):
        if self.escena is None:
            return
        r = self._ruta_guardar(_("Guardar la imagen 3D"), "cargas.png",
                               "PNG (*.png)")
        if r:
            self.escena.captura(r)

    def exportar_vtk(self):
        if not self._necesita_resultado():
            return
        r = self._ruta_guardar(_("Guardar para ParaView"), "cargas.vti",
                               "VTK ImageData (*.vti)")
        if r:
            C.exportar_vtk(self.res, r)

    def exportar_csv(self):
        if not self._necesita_resultado():
            return
        r = self._ruta_guardar(_("Guardar el CSV"), "cargas.csv",
                               "CSV (*.csv)")
        if r:
            C.exportar_csv(self.res, self.caso, r)

    def exportar_json(self):
        if not self._necesita_resultado():
            return
        r = self._ruta_guardar(_("Guardar el resultado"), "cargas.json",
                               "JSON (*.json)")
        if r:
            C.exportar_json(self.res, self.caso, r, self._estructura_info())

    def guardar_caso(self):
        r = self._ruta_guardar(_("Guardar el caso de carga"), "caso.json",
                               "JSON (*.json)")
        if r:
            C.guardar_caso(self.caso, r, self._estructura_info())

    def cargar_caso(self):
        f, _x = QtWidgets.QFileDialog.getOpenFileName(
            self, _("Cargar un caso de carga"), str(self.carpeta),
            "JSON (*.json)")
        if not f:
            return
        try:
            caso = C.cargar_caso(f)
        except C.ErrorCaso as e:
            QtWidgets.QMessageBox.warning(self, _("Caso no válido"), str(e))
            return
        self.poner_caso(caso)

    def poner_caso(self, caso):
        """Sustituye el caso (se puede deshacer)."""
        self.caso = copy.deepcopy(caso)
        self.caso.setdefault("material", {"E_s": C.E_S_DEF, "nu": C.NU_S_DEF})
        self._sincronizar_material()
        self._remapear_nodos()
        self._snap(_("cargar caso"), forzar=True)
        self._refrescar_listas()
        self._sel_carga(0 if self.caso["cargas"] else -1)
        self._invalidar_resultado()

    def registrar(self):
        if not self._necesita_resultado():
            return
        reg = C.resumen(self.res, self.caso, self._estructura_info())
        self.registros.append(reg)
        self.registrado.emit(reg)
        self.lab_estado.setText(_(
            "Caso añadido al informe de publicación ({n} en la sesión)."
        ).format(n=len(self.registros)))

    # ------------------------------------------------------------------
    def done(self, r):
        if self.escena is not None:
            self.escena.quitar_mango()
            self.escena.desactivar_pinchazo()
        if self.hilo is not None and self.hilo.isRunning():
            self.hilo.wait(30000)
        if self.vista3d is not None:
            try:
                self.vista3d.close()
            except Exception:
                pass
        super().done(r)


def _remapear_en(caso, BW, sp):
    """Vuelve a buscar los nodos elegidos (por sus coordenadas) en otra malla."""
    from spinpy import cargas as c
    caras = c.caras_expuestas(BW, sp)
    for lista in (caso["cargas"], caso["apoyos"]):
        for it in lista:
            r = it["region"]
            if r.get("forma") == "nodos" and r.get("xyz"):
                ids = []
                for x in r["xyz"]:
                    i = c.nodo_cercano(caras, x, BW.shape, sp)
                    if i is not None and i not in ids:
                        ids.append(i)
                r["ids"] = ids
