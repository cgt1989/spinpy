"""
dialogo_febio.py: Las ventanas de elementos finitos del visor.

(El nombre del archivo es el de antes de los motores internos; FEBio ya no se
ejecuta.) Dos puertas, un solo camino de calculo (`spinpy.fem.analizar` y
`fem.homogeneizar`, que el visor encadena en `Visor._febio_una`):

  DialogoFEMAuto   «Informe FEM (Auto)…», en la barra superior junto al
                   informe automatico. Corrida desatendida: estructuras x
                   protocolos x mallas x MOTORES, con tiempo y memoria
                   estimados y tabla comparativa entre motores.
  DialogoFEBio     «Comparar motores…», en la seccion Analisis mecanico. La
                   version corta: la estructura activa y los controles del
                   panel (resolucion, direccion, apoyo, motor).

Y las de apoyo: DialogoParametrosProtocolo (valores editables, en naranja los
que difieren del protocolo publicado), DialogoResultadosFEBio (tabla app vs
motores con hex8 y TET10, con citabilidad) y DialogoMapasFEBio (von Mises con
UNA escala para todos los paneles).

Nada aqui calcula: construye «tareas» (dict) que el visor ejecuta. Los textos
van a `_()` como LITERALES, porque `idioma_revisar.py` los recoge del arbol
sintactico de este archivo.
"""

from __future__ import annotations

import html
import json
import os
from pathlib import Path

import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets

from spinpy import fem as febio
from spinpy import motores, tiempos
from spinpy.idioma import _

NARANJA = "#b06000"
AJUSTES = ("spinpy", "visor")


def textos_analisis():
    return [("lineal", _("Lineal (extrapolado a carga nula; comparable 1:1 "
                         "con la app)")),
            ("nl_fuerza", _("No lineal, fuerza impuesta a la carga del "
                            "protocolo")),
            ("nl_plato", _("No lineal, plato rígido")),
            ("nl_pistoia", _("No lineal a la carga de fallo de Pistoia"))]


def etiqueta_estructura(codigo):
    return {"voi": _("VOI de referencia"),
            "spinodoide": _("Spinodoide ajustado"),
            "dual-lattice": _("Dual-lattice ajustado")}[codigo]


def etiqueta_malla(tipo):
    return {"hex8": _("Ladrillos (hex8)"),
            "tet10": _("Mallado suave (TET10)")}[tipo]


def descripcion_protocolo(p):
    """Una linea gris bajo cada protocolo, con los valores que se usaran."""
    ap = _("deslizante") if p.get("apoyo") == "deslizante" else _("empotrado")
    if p["tipo"] == "homogeneizacion":
        return _("periódica con ladrillos · cotas KUBC/SUBC con tetraedros · "
                 "E_s {E} GPa · ν {nu}").format(E=f"{p['E_s'] / 1e9:g}",
                                                nu=f"{p['nu']:g}")
    if p.get("carga_N") is not None:
        carga = f"{p['carga_N']:g} N"
    else:
        carga = _("{s} MPa sobre la sección bruta").format(
            s=f"{p['sigma_app'] / 1e6:g}")
    return _("{carga} · apoyo {apoyo} · E_s {E} GPa · ν {nu}").format(
        carga=carga, apoyo=ap, E=f"{p['E_s'] / 1e9:g}", nu=f"{p['nu']:g}")


def factores_febio():
    """Factores real/estimado de las etapas FEM en ESTE equipo."""
    try:
        crudo = QtCore.QSettings(*AJUSTES).value("tiempos/factores", "")
        f = json.loads(crudo) if crudo else {}
    except (TypeError, ValueError):
        f = {}
    return {k: v for k, v in f.items() if str(k).startswith("febio")}


# ---------------------------------------------------------------------------
# Parametros de un protocolo
# ---------------------------------------------------------------------------

class DialogoParametrosProtocolo(QtWidgets.QDialog):
    """Valores editables de un protocolo. En naranja los que difieren del
    publicado: si se cambia alguno del de Tapia, el registro deja de llamarse
    «Tapia» y pasa a «Tapia (modificado)» (`febio.protocolo`)."""

    def __init__(self, padre, clave, cambios=None):
        super().__init__(padre)
        self.clave = clave
        self.base = febio.PROTOCOLOS_FEM[clave]
        actual = dict(self.base)
        actual.update(cambios or {})
        self.setWindowTitle(_("Parámetros del protocolo"))
        lay = QtWidgets.QFormLayout(self)
        self.campos = {}

        def num(clave_p, etq, esc, dec, lo, hi, suf):
            w = QtWidgets.QDoubleSpinBox()
            w.setDecimals(dec)
            w.setRange(lo, hi)
            w.setSuffix(suf)
            w.setValue(float(actual[clave_p]) / esc)
            w.valueChanged.connect(self._marcar)
            self.campos[clave_p] = (w, esc)
            lay.addRow(etq, w)

        num("E_s", _("Módulo del tejido E_s:"), 1e9, 3, 0.1, 100.0, " GPa")
        num("nu", _("Coeficiente de Poisson ν:"), 1.0, 3, 0.0, 0.499, "")
        if self.base["tipo"] == "compresion":
            if self.base.get("carga_N") is not None:
                num("carga_N", _("Carga axial:"), 1.0, 2, 0.01, 1e5, " N")
            else:
                num("sigma_app", _("Tensión aparente:"), 1e6, 4, 1e-4, 1e3,
                    " MPa")
            self.cmb_apoyo = QtWidgets.QComboBox()
            self.cmb_apoyo.addItems([_("deslizante"), _("empotrado")])
            self.cmb_apoyo.setCurrentIndex(
                0 if actual["apoyo"] == "deslizante" else 1)
            self.cmb_apoyo.currentIndexChanged.connect(self._marcar)
            lay.addRow(_("Apoyo:"), self.cmb_apoyo)
        else:
            self.cmb_apoyo = None
            num("amplitud", _("Amplitud de la deformación:"), 1.0, 7, 1e-7,
                1e-2, "")
        self.lab = QtWidgets.QLabel()
        self.lab.setWordWrap(True)
        lay.addRow(self.lab)
        bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok
                                        | QtWidgets.QDialogButtonBox.Cancel)
        b = bb.addButton(_("Valores publicados"),
                         QtWidgets.QDialogButtonBox.ResetRole)
        b.clicked.connect(self._restaurar)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay.addRow(bb)
        self._marcar()

    def cambios(self):
        c = {k: w.value() * esc for k, (w, esc) in self.campos.items()}
        if self.cmb_apoyo is not None:
            c["apoyo"] = ("deslizante" if self.cmb_apoyo.currentIndex() == 0
                          else "empotrado")
        def igual(v, b):
            # El spinbox devuelve valor/escala*escala: tolerancia relativa,
            # o un valor sin tocar podria salir «modificado» por redondeo.
            if isinstance(v, float) and isinstance(b, (int, float)):
                return abs(v - b) <= 1e-9 * max(abs(b), 1e-300)
            return v == b
        return {k: v for k, v in c.items() if not igual(v, self.base.get(k))}

    def _restaurar(self):
        for k, (w, esc) in self.campos.items():
            w.setValue(float(self.base[k]) / esc)
        if self.cmb_apoyo is not None:
            self.cmb_apoyo.setCurrentIndex(
                0 if self.base["apoyo"] == "deslizante" else 1)

    def _marcar(self, *_a):
        ch = self.cambios()
        for k, (w, _esc) in self.campos.items():
            w.setStyleSheet(f"color:{NARANJA}; font-weight:bold;"
                            if k in ch else "")
        if self.cmb_apoyo is not None:
            self.cmb_apoyo.setStyleSheet(
                f"color:{NARANJA}; font-weight:bold;" if "apoyo" in ch
                else "")
        self.lab.setText(
            _("Registro: «{n}». Los valores en naranja difieren del protocolo "
              "publicado.").format(n=febio.nombre_protocolo(self.clave,
                                                             bool(ch)))
            if ch else _("Valores publicados del protocolo."))


# ---------------------------------------------------------------------------
# Piezas comunes de las dos ventanas
# ---------------------------------------------------------------------------

class _BaseFEBio(QtWidgets.QDialog):
    """Malla, analisis, material y solver: lo que comparten las dos puertas."""

    def _construir_malla(self, v, n_hex_def):
        g = QtWidgets.QGroupBox(_("Malla"))
        gl = QtWidgets.QVBoxLayout(g)
        self.rb_ambas = QtWidgets.QRadioButton(_("Ambas (recomendada)"))
        self.rb_hex = QtWidgets.QRadioButton(
            _("Ladrillos (hex8, un vóxel = un elemento)"))
        self.rb_tet = QtWidgets.QRadioButton(
            _("Mallado suave (tetraedros TET10)"))
        self.rb_ambas.setChecked(True)
        for rb in (self.rb_ambas, self.rb_hex, self.rb_tet):
            gl.addWidget(rb)
        self.lab_malla = QtWidgets.QLabel()
        self.lab_malla.setWordWrap(True)
        self.lab_malla.setStyleSheet("color:#666; font-size:10px;")
        gl.addWidget(self.lab_malla)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(_("Ladrillos, resolución:")))
        self.spin_nhex = QtWidgets.QSpinBox()
        self.spin_nhex.setRange(12, 96)
        self.spin_nhex.setValue(int(n_hex_def))
        self.spin_nhex.setSuffix(" vox")
        f.addWidget(self.spin_nhex)
        f.addStretch(1)
        gl.addLayout(f)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(_("Suave, resolución de la máscara:")))
        self.spin_ntet = QtWidgets.QSpinBox()
        self.spin_ntet.setRange(16, 96)
        self.spin_ntet.setValue(febio.N_TET_DEF)
        self.spin_ntet.setSuffix(" vox")
        self.spin_ntet.setToolTip(_(
            "Por debajo de 48³ los puntales tienen ~2 vóxeles y el suavizado "
            "los estrecha: en el VOI proximal de H4 a 32³ la malla suave salió "
            "2,3 veces más blanda que la de ladrillos. A 48³ hacen falta ~7 GB "
            "de memoria."))
        f.addWidget(self.spin_ntet)
        f.addStretch(1)
        gl.addLayout(f)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(_("Tamaño máx. de tetraedro:")))
        self.spin_tam = QtWidgets.QDoubleSpinBox()
        self.spin_tam.setDecimals(3)
        self.spin_tam.setRange(0.0, 5.0)
        self.spin_tam.setSingleStep(0.005)
        self.spin_tam.setSuffix(" mm")
        self.spin_tam.setSpecialValueText(_("el de la superficie"))
        tbth = (v.m_voi or {}).get("TbTh") if v.m_voi else None
        if tbth and np.isfinite(tbth):
            self.spin_tam.setToolTip(_("Sugerido ≈ Tb.Th/3 = {t} mm").format(
                t=f"{tbth / 3:.3f}"))
        f.addWidget(self.spin_tam)
        b = QtWidgets.QPushButton(_("Tapia: 0,05 mm"))
        b.clicked.connect(lambda: self.spin_tam.setValue(0.05))
        f.addWidget(b)
        f.addStretch(1)
        gl.addLayout(f)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(_("Suavizado Taubin:")))
        self.spin_taubin = QtWidgets.QSpinBox()
        self.spin_taubin.setRange(0, 400)
        self.spin_taubin.setValue(20)
        self.spin_taubin.setSuffix(" it.")
        f.addWidget(self.spin_taubin)
        f.addWidget(QtWidgets.QLabel(_("Decimado:")))
        self.spin_dec = QtWidgets.QDoubleSpinBox()
        self.spin_dec.setRange(0.0, 0.9)
        self.spin_dec.setSingleStep(0.1)
        self.spin_dec.setValue(0.5)
        f.addWidget(self.spin_dec)
        f.addStretch(1)
        gl.addLayout(f)
        self.chk_corr = QtWidgets.QCheckBox(
            _("Corregir la pérdida de volumen del suavizado"))
        self.chk_corr.setChecked(True)
        self.chk_corr.setToolTip(_(
            "Desplaza la superficie por su normal hasta el volumen de los "
            "vóxeles. Medido en el VOI proximal de H4 a 48³: −6,8 % sin "
            "corregir, −0,15 % corregido."))
        gl.addWidget(self.chk_corr)
        self.chk_conv = QtWidgets.QCheckBox(_(
            "Verificar convergencia con dos tamaños de malla (repite el "
            "lineal con la mitad del tamaño)"))
        gl.addWidget(self.chk_conv)
        for rb in (self.rb_ambas, self.rb_hex, self.rb_tet):
            rb.toggled.connect(self._malla_cambiada)
        self._malla_cambiada()
        return g

    def _malla_cambiada(self, *_a):
        suave = not self.rb_hex.isChecked()
        for w in (self.spin_ntet, self.spin_tam, self.spin_taubin,
                  self.spin_dec, self.chk_corr, self.chk_conv):
            w.setEnabled(suave)
        self.spin_nhex.setEnabled(not self.rb_tet.isChecked())
        self.lab_malla.setText(
            _("Ambas separa el efecto del programa (app frente a los motores "
              "con ladrillos, debe dar ~0) del efecto de la malla (ladrillos "
              "frente a tetraedros, mismo motor y carga).")
            if self.rb_ambas.isChecked() else
            _("Ladrillos: la malla que resuelve la app; los motores la "
              "reproducen a ~1e-10. Conserva los escalones.")
            if self.rb_hex.isChecked() else
            _("Mallado suave: quita los escalones, pero su borde es una de "
              "muchas superficies compatibles con la imagen; el pico de von "
              "Mises depende del suavizado."))

    def mallas(self):
        if self.rb_hex.isChecked():
            return ["hex8"]
        if self.rb_tet.isChecked():
            return ["tet10"]
        return ["hex8", "tet10"]

    def opciones_malla(self):
        return {"suavizado": int(self.spin_taubin.value()),
                "decimado": float(self.spin_dec.value()),
                "tam_max_mm": (float(self.spin_tam.value())
                               if self.spin_tam.value() > 0 else None),
                "corregir_volumen": self.chk_corr.isChecked()}

    def _construir_analisis(self, con_direccion=True):
        g = QtWidgets.QGroupBox(_("Tipo de análisis (protocolos de "
                                  "compresión)"))
        gl = QtWidgets.QVBoxLayout(g)
        self.chk_an = {}
        for k, t in textos_analisis():
            c = QtWidgets.QCheckBox(t)
            c.setChecked(k in ("lineal", "nl_fuerza", "nl_plato"))
            self.chk_an[k] = c
            gl.addWidget(c)
        nota = QtWidgets.QLabel(_(
            "El plato rígido y la carga de Pistoia necesitan el lineal: si no "
            "se marca, se corre igual."))
        nota.setWordWrap(True)
        nota.setStyleSheet("color:#666; font-size:10px;")
        gl.addWidget(nota)
        self.cmb_dir = None
        if con_direccion:
            f = QtWidgets.QHBoxLayout()
            f.addWidget(QtWidgets.QLabel(_("Dirección de carga:")))
            self.cmb_dir = QtWidgets.QComboBox()
            self.cmb_dir.addItems([_("Z (axial)"), _("X, Y y Z")])
            f.addWidget(self.cmb_dir, 1)
            gl.addLayout(f)
        return g

    def analisis(self):
        return [k for k, c in self.chk_an.items() if c.isChecked()]

    def _construir_solver(self):
        g = QtWidgets.QGroupBox(_("Material y solver"))
        gl = QtWidgets.QFormLayout(g)
        self.cmb_mat = QtWidgets.QComboBox()
        self.cmb_mat.addItem(_("St. Venant-Kirchhoff (validado)"), "svk")
        self.cmb_mat.addItem(_("Neo-Hookeano"), "neohookeano")
        gl.addRow(_("Material no lineal:"), self.cmb_mat)
        self.spin_hilos = QtWidgets.QSpinBox()
        n = os.cpu_count() or 2
        self.spin_hilos.setRange(1, n)
        self.spin_hilos.setValue(max(1, n - 1))
        self.spin_hilos.setToolTip(_("Núcleos − 1 por omisión, para que la "
                                     "ventana siga respondiendo."))
        gl.addRow(_("Hilos de CPU:"), self.spin_hilos)
        self.spin_pasos = QtWidgets.QSpinBox()
        self.spin_pasos.setRange(1, 50)
        self.spin_pasos.setValue(1)
        self.spin_pasos.setToolTip(_(
            "Incrementos de carga del análisis no lineal. Con cargas "
            "moderadas basta uno; si Newton no converge, subirlos."))
        gl.addRow(_("Pasos de carga (no lineal):"), self.spin_pasos)
        self.cmb_solver = QtWidgets.QComboBox()
        self.cmb_solver.addItem(_("Automático (según el tamaño)"), "auto")
        self.cmb_solver.addItem(_("Directo"), "directo")
        self.cmb_solver.addItem(_("Iterativo (multigrid)"), "iterativo")
        gl.addRow(_("Resolvedor lineal:"), self.cmb_solver)
        nota = QtWidgets.QLabel(_(
            "Todo se resuelve dentro de spinpy, sin programas externos. "
            "Newton completo con tangente consistente; tolerancia 1e-10."))
        nota.setWordWrap(True)
        nota.setStyleSheet("color:#666; font-size:10px;")
        gl.addRow(nota)
        return g

    def _construir_motores(self, varios=True):
        """Casillas de los motores FEM: los instalados, habilitados; los
        demas, en gris con «(no instalado)». La app va siempre."""
        g = QtWidgets.QGroupBox(_("Motores FEM"))
        gl = QtWidgets.QVBoxLayout(g)
        self.chk_mot = {}
        disp = motores.disponibles()
        for k in motores.MOTORES:
            ver = motores.version(k) if k in disp else None
            t = motores.ETIQUETAS[k] + (f" {ver}" if ver else "")
            if k == motores.RECOMENDADO:
                t += "  " + _("(recomendado para TET10 y no lineal)")
            if k not in disp:
                t += "  " + _("(no instalado)")
            c = QtWidgets.QCheckBox(t)
            c.setEnabled(k in disp)
            c.setChecked(k in ("app", motores.RECOMENDADO) and k in disp)
            self.chk_mot[k] = c
            gl.addWidget(c)
        nota = QtWidgets.QLabel(_(
            "Cada motor resuelve la MISMA malla: la tabla compara los motores "
            "entre sí y con la app. La app solo resuelve ladrillos en lineal; "
            "en lo demás su columna queda «no disponible»."))
        nota.setWordWrap(True)
        nota.setStyleSheet("color:#666; font-size:10px;")
        gl.addWidget(nota)
        return g

    def motores_sel(self):
        return [k for k, c in self.chk_mot.items()
                if c.isChecked() and c.isEnabled()]

    def _cabecera_febio(self):
        ms = self.motores_sel()
        if not ms:
            return ("<span style='color:#b62324'><b>"
                    + _("Ningún motor elegido.") + "</b></span>")
        return _("Resuelve dentro de spinpy con <b>{m}</b>").format(
            m=html.escape(", ".join(motores.ETIQUETAS[k] for k in ms)))

    def comunes(self):
        return {"mallas": self.mallas(), "analisis": self.analisis(),
                "n_hex": int(self.spin_nhex.value()),
                "n_tet": int(self.spin_ntet.value()),
                "opciones_malla": self.opciones_malla(),
                "conv_malla": self.chk_conv.isChecked(),
                "material": self.cmb_mat.currentData(),
                "hilos": int(self.spin_hilos.value()),
                "pasos": int(self.spin_pasos.value()),
                "motores": self.motores_sel(),
                "solver": self.cmb_solver.currentData()}


# ---------------------------------------------------------------------------
# Estimacion
# ---------------------------------------------------------------------------

def estimar_tareas(tareas, tamanos, factores=None):
    """[(segundos, MB de pico)] de cada tarea, con `spinpy.tiempos`.

    `tamanos[(estructura, malla, n)]` es `fem.tamano_previsto`. El modelo es
    el medido con FEBio (PARDISO); cada motor adicional suma una resolucion.
    """
    f = factores or {}
    out = []
    for t in tareas:
        tam = tamanos.get((t["estructura"], t["malla"], t["n"]))
        if tam is None:
            out.append((None, None))
            continue
        g = tam["gdl"]
        if t["protocolo"]["tipo"] == "homogeneizacion":
            s = tiempos.febio_homog(t["malla"], g if t["malla"] == "tet10"
                                    else 3 * (t["n"] + 1) ** 3,
                                    n_tet=tam["n_elems"])
        else:
            s = tiempos.febio_ensayo(t["malla"], g, t["analisis"],
                                     n_tet=tam["n_elems"],
                                     n_ejes=len(t["ejes"]),
                                     conv_malla=t.get("conv_malla"))
        s *= float(f.get(f"febio_{t['malla']}", 1.0))
        s *= max(1, len(t.get("motores") or [1]))
        out.append((s, tam["memoria_MB"] * (2 if t.get("conv_malla")
                                            and t["malla"] == "tet10"
                                            else 1)))
    return out


# ---------------------------------------------------------------------------
# FEM automatico
# ---------------------------------------------------------------------------

class DialogoFEMAuto(_BaseFEBio):
    """«Informe FEM (Auto)»: la maqueta aprobada, `Maqueta_FEM_automatico
    .html`. Mismo formato que `DialogoInformeAuto`: dos columnas, tiempo
    estimado abajo y Empezar/Cancelar."""

    def __init__(self, v, carpeta_def=""):
        super().__init__(v)
        self.v = v
        self._cambios = {k: {} for k in febio.PROTOCOLOS_FEM}
        self._tamanos = {}
        self._ocupado_act = False
        self.setWindowTitle(_("Informe FEM (Auto)"))
        raiz = QtWidgets.QVBoxLayout(self)
        self.intro = QtWidgets.QLabel()
        self.intro.setWordWrap(True)
        self.intro.setTextFormat(QtCore.Qt.RichText)
        raiz.addWidget(self.intro)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(_("Carpeta de resultados:")))
        self.ed_carpeta = QtWidgets.QLineEdit(carpeta_def)
        self.ed_carpeta.setReadOnly(True)
        f.addWidget(self.ed_carpeta, 1)
        b = QtWidgets.QPushButton(_("Carpeta…"))
        b.clicked.connect(self._elegir_carpeta)
        f.addWidget(b)
        raiz.addLayout(f)

        cols = QtWidgets.QHBoxLayout()
        izq, der = QtWidgets.QVBoxLayout(), QtWidgets.QVBoxLayout()
        cols.addLayout(izq, 1)
        cols.addLayout(der, 1)
        raiz.addLayout(cols)

        # -- estructuras --
        g = QtWidgets.QGroupBox(_("Estructuras"))
        gl = QtWidgets.QVBoxLayout(g)
        self.chk_est = {}
        for cod in ("voi", "spinodoide", "dual-lattice"):
            c = QtWidgets.QCheckBox(etiqueta_estructura(cod))
            self.chk_est[cod] = c
            gl.addWidget(c)
        self.chk_ajustar = QtWidgets.QCheckBox(_(
            "Ajustar antes si falta (mejor ajuste: rápido, completo, "
            "equitativo)"))
        self.chk_ajustar.setChecked(True)
        gl.addWidget(self.chk_ajustar)
        izq.addWidget(g)

        # -- protocolos --
        g = QtWidgets.QGroupBox(_("Protocolos"))
        gl = QtWidgets.QGridLayout(g)
        self.chk_prot, self.lab_prot = {}, {}
        for i, (clave, p) in enumerate(febio.PROTOCOLOS_FEM.items()):
            c = QtWidgets.QCheckBox(self._etq_protocolo(clave))
            c.setChecked(clave != "homogeneizacion")
            self.chk_prot[clave] = c
            gl.addWidget(c, 2 * i, 0)
            b = QtWidgets.QPushButton(_("Parámetros…"))
            b.clicked.connect(lambda _c=False, k=clave: self._parametros(k))
            gl.addWidget(b, 2 * i, 1)
            lab = QtWidgets.QLabel()
            lab.setStyleSheet("color:#666; font-size:10px; margin-left:22px;")
            lab.setWordWrap(True)
            self.lab_prot[clave] = lab
            gl.addWidget(lab, 2 * i + 1, 0, 1, 2)
        izq.addWidget(g)
        izq.addWidget(self._construir_analisis(True))
        izq.addStretch(1)

        der.addWidget(self._construir_motores())
        der.addWidget(self._construir_malla(v, febio.N_HEX_DEF))
        der.addWidget(self._construir_solver())

        # -- salidas --
        g = QtWidgets.QGroupBox(_("Salidas"))
        gl = QtWidgets.QVBoxLayout(g)
        self.chk_tabla = QtWidgets.QCheckBox(_(
            "Tabla comparativa entre motores y frente a la app, con "
            "citabilidad"))
        self.chk_tabla.setChecked(True)
        self.chk_tabla.setEnabled(False)
        self.chk_fig = QtWidgets.QCheckBox(_(
            "Figura de validación (5 columnas: E_app, p99 de superficie, "
            "fallo, desvío no lineal, BV/TV de la malla)"))
        self.chk_fig.setChecked(True)
        self.chk_mapas = QtWidgets.QCheckBox(_(
            "Mapas 3D de von Mises con una sola escala"))
        self.chk_mapas.setChecked(True)
        self.chk_conservar = QtWidgets.QCheckBox(_(
            "Guardar los campos (desplazamientos) de cada motor (pesan)"))
        self.chk_sesion = QtWidgets.QCheckBox(_(
            "Incluir los resultados en el informe de publicación"))
        self.chk_sesion.setChecked(True)
        for c in (self.chk_tabla, self.chk_fig, self.chk_mapas,
                  self.chk_conservar, self.chk_sesion):
            gl.addWidget(c)
        der.addWidget(g)
        der.addStretch(1)

        self.lab_avisos = QtWidgets.QLabel()
        self.lab_avisos.setWordWrap(True)
        self.lab_avisos.setStyleSheet(
            "background:#fff4ce; color:#7a5700; border:1px solid #e8d38a;"
            "border-radius:4px; padding:6px;")
        raiz.addWidget(self.lab_avisos)

        marco = QtWidgets.QFrame()
        marco.setStyleSheet("QFrame { background:#f3f2ee; border-radius:4px; }")
        fl = QtWidgets.QVBoxLayout(marco)
        self.lab_total = QtWidgets.QLabel()
        self.lab_total.setTextFormat(QtCore.Qt.RichText)
        fl.addWidget(self.lab_total)
        self.lab_mem = QtWidgets.QLabel()
        self.lab_mem.setWordWrap(True)
        fl.addWidget(self.lab_mem)
        self.lab_nota = QtWidgets.QLabel(_(
            "Modelo PROVISIONAL: medido con el equipo de desarrollo cargado; "
            "puede errar un factor 2 hasta recalibrarlo. Se corrige solo con "
            "cada etapa terminada en este equipo."))
        self.lab_nota.setWordWrap(True)
        self.lab_nota.setStyleSheet("color:#666; font-size:10px;")
        fl.addWidget(self.lab_nota)
        raiz.addWidget(marco)

        bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Cancel)
        self.b_ok = bb.addButton(_("Empezar"),
                                 QtWidgets.QDialogButtonBox.AcceptRole)
        bb.accepted.connect(self._aceptar)
        bb.rejected.connect(self.reject)
        raiz.addWidget(bb)

        self._restaurar()
        for c in self.findChildren(QtWidgets.QCheckBox):
            c.toggled.connect(self._actualizar)
        for c in self.findChildren(QtWidgets.QRadioButton):
            c.toggled.connect(self._actualizar)
        for s in self.findChildren(QtWidgets.QAbstractSpinBox):
            s.editingFinished.connect(self._actualizar)
        for c in self.findChildren(QtWidgets.QComboBox):
            c.currentIndexChanged.connect(self._actualizar)
        self._actualizar()

    # -- protocolos --

    def _etq_protocolo(self, clave):
        return {"app": _("Ensayo de la app"),
                "tapia2026": _("Protocolo de Tapia et al. (2026)"),
                "homogeneizacion": _("Homogeneización (tensor elástico)")
                }.get(clave, clave)

    def _parametros(self, clave):
        d = DialogoParametrosProtocolo(self, clave, self._cambios[clave])
        if d.exec_() == QtWidgets.QDialog.Accepted:
            self._cambios[clave] = d.cambios()
            self._actualizar()

    def protocolos(self):
        return [febio.protocolo(k, **self._cambios[k])
                for k, c in self.chk_prot.items() if c.isChecked()]

    # -- estructuras --

    def _disponible(self, cod):
        if cod == "voi":
            return self.v.VOI is not None
        return self.v._de(cod, "BW") is not None

    def estructuras(self):
        return [c for c, w in self.chk_est.items()
                if w.isChecked() and w.isEnabled()]

    def por_ajustar(self):
        return [c for c in self.estructuras()
                if c != "voi" and not self._disponible(c)]

    # -- plan --

    def tareas(self):
        com = self.comunes()
        ejes = [2] if (self.cmb_dir is None or self.cmb_dir.currentIndex()
                       == 0) else [0, 1, 2]
        carpeta = self.ed_carpeta.text()
        out = []
        for est in self.estructuras():
            for p in self.protocolos():
                for m in com["mallas"]:
                    out.append({
                        "estructura": est, "protocolo": p, "malla": m,
                        "analisis": com["analisis"] or ["lineal"],
                        "ejes": ejes if p["tipo"] == "compresion" else [2],
                        "n": com["n_hex"] if m == "hex8" else com["n_tet"],
                        "opciones_malla": com["opciones_malla"],
                        "conv_malla": com["conv_malla"] and m == "tet10",
                        "material": com["material"], "pasos": com["pasos"],
                        "hilos": com["hilos"], "motores": com["motores"],
                        "solver": com["solver"],
                        "carpeta": carpeta,
                        "conservar": self.chk_conservar.isChecked(),
                        "comparar_app": True})
        return out

    def _tamano(self, est, malla, n):
        clave = (est, malla, n)
        if clave not in self._tamanos:
            if est == "voi" or not self._disponible(est):
                # Un candidato por ajustar se estima con el VOI: tiene su
                # BV/TV y su Tb.Th, que es lo que fija el tamano.
                if self.v.VOI is None:
                    return None
                BW, sp = self.v.VOI, self.v.VOI_spacing
            else:
                BW = self.v._de(est, "BW")
                sp = self.v._spacing(BW.shape[0])
            self._tamanos[clave] = febio.tamano_previsto(BW, sp, malla, n)
        return self._tamanos[clave]

    def _actualizar(self, *_a):
        if self._ocupado_act:
            return
        self._ocupado_act = True
        try:
            self._actualizar_una()
        finally:
            self._ocupado_act = False

    def _actualizar_una(self):
        v = self.v
        for cod, w in self.chk_est.items():
            ok = self._disponible(cod)
            puede = ok or (cod != "voi" and self.chk_ajustar.isChecked()
                           and v.VOI is not None)
            w.setEnabled(puede)
            w.setText(etiqueta_estructura(cod) + (
                f"  —  {v.VOI_nombre}" if cod == "voi" and ok else
                "" if ok else "  " + _("(no generado)")))
            w.setStyleSheet("" if ok else "color:#999;")
        for clave, lab in self.lab_prot.items():
            p = febio.protocolo(clave, **self._cambios[clave])
            txt = descripcion_protocolo(p)
            if p["modificado"]:
                txt += "  —  " + _("modificado: «{n}»").format(n=p["nombre"])
            lab.setText(txt)
            lab.setStyleSheet(("color:%s;" % NARANJA if p["modificado"]
                               else "color:#666;")
                              + " font-size:10px; margin-left:22px;")
        self.intro.setText(
            self._cabecera_febio() + " " + _(
                "las estructuras, protocolos y mallas marcados, y los compara "
                "con el ensayo de la app. Los diálogos de resultados no se "
                "abren: todo queda en la sesión y en la carpeta. Puede tardar "
                "<b>horas</b>; la ventana sigue respondiendo y se puede "
                "detener tras la etapa en curso."))

        tareas = self.tareas()
        for t in tareas:
            self._tamano(t["estructura"], t["malla"], t["n"])
        est = estimar_tareas(tareas, self._tamanos, factores_febio())
        self._estimados = [s for s, _m in est]
        s = sum(x for x in self._estimados if x)
        s += 1100.0 * len(self.por_ajustar())       # ajuste: ~18 min medido
        fin = QtCore.QDateTime.currentDateTime().addSecs(int(s))
        self.lab_total.setText(
            _("<b>Tiempo estimado: {t}</b> · terminaría hacia las {h} · "
              "{n} etapa(s)").format(t=tiempos.texto(s),
                                     h=fin.toString("HH:mm"), n=len(tareas)))
        mem = max([m for _s, m in est if m] or [0.0])
        tot = febio.memoria_equipo_MB()
        self._memoria = (mem, tot)
        if tot:
            self.lab_mem.setText(
                _("Memoria máx. estimada: {m} GB de {t} GB").format(
                    m=f"{mem / 1024:.1f}", t=f"{tot / 1024:.0f}"))
            self.lab_mem.setStyleSheet(
                "color:#b62324; font-weight:bold;" if mem > 0.7 * tot
                else "color:#52514e;")

        av = self.avisos(tareas)
        self.lab_avisos.setText("\n".join("⚠ " + a for a in av))
        self.lab_avisos.setVisible(bool(av))
        self.b_ok.setEnabled(bool(tareas) and bool(self.motores_sel())
                             and bool(self.ed_carpeta.text())
                             and bool(self.analisis() or all(
                                 t["protocolo"]["tipo"] == "homogeneizacion"
                                 for t in tareas)))

    def avisos(self, tareas):
        from spinpy.avisos import desalineacion, voi_no_trabecular
        v = self.v
        a = []
        if not self.motores_sel():
            a.append(_("Elige al menos un motor FEM."))
        if "tet10" in self.mallas() and self.motores_sel() == ["app"]:
            a.append(_("La app no resuelve la malla suave: marca otro motor "
                       "para TET10."))
        if v.m_voi:
            nt = voi_no_trabecular(v.m_voi)
            if nt.get("aviso"):
                a.append(_("VOI no trabecular (BV/TV {b}): ninguna familia es "
                           "bicontinua ahí.").format(b=f"{nt['BVTV']:.2f}"))
        for est in self.estructuras():
            if est == "voi":
                continue
            m = v._de(est, "m")
            if v.m_voi and m:
                try:
                    d = desalineacion(v.m_voi, m)
                    if d.get("aviso"):
                        a.append(_("{e}: eje a {g}° del VOI.").format(
                            e=etiqueta_estructura(est),
                            g=f"{d['angulo_deg']:.0f}"))
                except Exception:
                    pass
        if "tet10" in self.mallas() and self.spin_ntet.value() < 48:
            a.append(_("Mallado suave a {n}³, por debajo de 48³: en el VOI "
                       "proximal de H4 a 32³ el suavizado estrecha los "
                       "puntales y la malla sale 2,3 veces más blanda que la de "
                       "ladrillos.").format(n=self.spin_ntet.value()))
        if "hex8" in self.mallas() and self.spin_nhex.value() < 40:
            a.append(_("Ladrillos a {n}³, por debajo de 40³: la rigidez del "
                       "VOI proximal de H4 cae ~10 % a 32³.").format(
                n=self.spin_nhex.value()))
        vistos = set()
        for t in tareas:
            tam = self._tamanos.get((t["estructura"], t["malla"], t["n"]))
            # Solo de lo que existe: un candidato por ajustar se ESTIMA con
            # el VOI, y su desconexion no es la del VOI.
            if tam and tam["frac_portante"] < 0.99 and t["estructura"] \
                    not in vistos and self._disponible(t["estructura"]):
                vistos.add(t["estructura"])
                a.append(_("{e}: {p} % del hueso desconectado (se filtra "
                           "antes de mallar).").format(
                    e=etiqueta_estructura(t["estructura"]),
                    p=f"{100 * (1 - tam['frac_portante']):.1f}"))
        for est in self.por_ajustar():
            a.append(_("{e} no generado: se ajustará antes porque está marcado "
                       "«Ajustar antes si falta».").format(
                e=etiqueta_estructura(est)))
        mem, tot = getattr(self, "_memoria", (0, None))
        if tot and mem > 0.7 * tot:
            a.append(_("La memoria estimada supera el 70 % de la del equipo: "
                       "baja la resolución de la malla suave."))
        return a

    def _elegir_carpeta(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(
            self, _("Carpeta de resultados FEM"), self.ed_carpeta.text())
        if d:
            self.ed_carpeta.setText(d)
            self._actualizar()

    # -- recordar la ultima configuracion --

    def _estado(self):
        return {"estructuras": {k: c.isChecked()
                                for k, c in self.chk_est.items()},
                "ajustar": self.chk_ajustar.isChecked(),
                "protocolos": {k: c.isChecked()
                               for k, c in self.chk_prot.items()},
                "cambios": self._cambios,
                "analisis": {k: c.isChecked() for k, c in self.chk_an.items()},
                "dir": self.cmb_dir.currentIndex(),
                "malla": ("hex8" if self.rb_hex.isChecked() else "tet10"
                          if self.rb_tet.isChecked() else "ambas"),
                "n_hex": self.spin_nhex.value(),
                "n_tet": self.spin_ntet.value(),
                "tam": self.spin_tam.value(),
                "taubin": self.spin_taubin.value(),
                "dec": self.spin_dec.value(),
                "corr": self.chk_corr.isChecked(),
                "conv": self.chk_conv.isChecked(),
                "mat": self.cmb_mat.currentIndex(),
                "motores": self.motores_sel(),
                "solver": self.cmb_solver.currentIndex(),
                "pasos": self.spin_pasos.value(),
                "fig": self.chk_fig.isChecked(),
                "mapas": self.chk_mapas.isChecked(),
                "conservar": self.chk_conservar.isChecked(),
                "sesion": self.chk_sesion.isChecked()}

    def _restaurar(self):
        try:
            crudo = QtCore.QSettings(*AJUSTES).value("fem_auto", "")
            e = json.loads(crudo) if crudo else {}
        except (TypeError, ValueError):
            e = {}
        if not e:
            for k, c in self.chk_est.items():
                c.setChecked(True)
            return
        try:
            for k, c in self.chk_est.items():
                c.setChecked(bool(e["estructuras"].get(k, True)))
            self.chk_ajustar.setChecked(bool(e["ajustar"]))
            for k, c in self.chk_prot.items():
                c.setChecked(bool(e["protocolos"].get(k, False)))
            self._cambios.update({k: dict(v) for k, v in e["cambios"].items()
                                  if k in self._cambios})
            for k, c in self.chk_an.items():
                c.setChecked(bool(e["analisis"].get(k, False)))
            self.cmb_dir.setCurrentIndex(int(e["dir"]))
            {"hex8": self.rb_hex, "tet10": self.rb_tet}.get(
                e["malla"], self.rb_ambas).setChecked(True)
            self.spin_nhex.setValue(int(e["n_hex"]))
            self.spin_ntet.setValue(int(e["n_tet"]))
            self.spin_tam.setValue(float(e["tam"]))
            self.spin_taubin.setValue(int(e["taubin"]))
            self.spin_dec.setValue(float(e["dec"]))
            self.chk_corr.setChecked(bool(e["corr"]))
            self.chk_conv.setChecked(bool(e["conv"]))
            self.cmb_mat.setCurrentIndex(int(e["mat"]))
            for k, c in self.chk_mot.items():
                if c.isEnabled() and "motores" in e:
                    c.setChecked(k in e["motores"])
            self.cmb_solver.setCurrentIndex(int(e.get("solver", 0)))
            self.spin_pasos.setValue(int(e["pasos"]))
            self.chk_fig.setChecked(bool(e["fig"]))
            self.chk_mapas.setChecked(bool(e["mapas"]))
            self.chk_conservar.setChecked(bool(e["conservar"]))
            self.chk_sesion.setChecked(bool(e["sesion"]))
        except (KeyError, TypeError, ValueError):
            pass                  # una clave vieja no rompe la ventana

    def _aceptar(self):
        QtCore.QSettings(*AJUSTES).setValue("fem_auto",
                                            json.dumps(self._estado()))
        self.accept()

    def opciones(self):
        return {"carpeta": self.ed_carpeta.text(),
                "tareas": self.tareas(), "estimados": list(self._estimados),
                "por_ajustar": self.por_ajustar(),
                "figura": self.chk_fig.isChecked(),
                "mapas": self.chk_mapas.isChecked(),
                "sesion": self.chk_sesion.isChecked(),
                "configuracion": self._estado()}


# ---------------------------------------------------------------------------
# Comparar motores (panel)
# ---------------------------------------------------------------------------

class DialogoFEBio(_BaseFEBio):
    """La version corta del panel: la estructura activa (y el VOI) con los
    controles del ensayo del panel. Protocolo, malla y analisis aqui."""

    def __init__(self, v, carpeta_def=""):
        super().__init__(v)
        self.v = v
        self._cambios = {}
        self.setWindowTitle(_("Comparar motores FEM"))
        raiz = QtWidgets.QVBoxLayout(self)
        self.intro = QtWidgets.QLabel()
        self.intro.setWordWrap(True)
        self.intro.setTextFormat(QtCore.Qt.RichText)
        raiz.addWidget(self.intro)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(_("Protocolo:")))
        self.cmb_prot = QtWidgets.QComboBox()
        for k in febio.PROTOCOLOS_FEM:
            self.cmb_prot.addItem({"app": _("Ensayo de la app"),
                                   "tapia2026": _("Protocolo de Tapia et al. "
                                                  "(2026)"),
                                   "homogeneizacion": _("Homogeneización "
                                                        "(tensor elástico)")
                                   }[k], k)
        f.addWidget(self.cmb_prot, 1)
        b = QtWidgets.QPushButton(_("Parámetros…"))
        b.clicked.connect(self._parametros)
        f.addWidget(b)
        raiz.addLayout(f)
        self.lab_prot = QtWidgets.QLabel()
        self.lab_prot.setStyleSheet("color:#666; font-size:10px;")
        raiz.addWidget(self.lab_prot)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(_("Carpeta de resultados:")))
        self.ed_carpeta = QtWidgets.QLineEdit(carpeta_def)
        self.ed_carpeta.setReadOnly(True)
        f.addWidget(self.ed_carpeta, 1)
        b = QtWidgets.QPushButton(_("Carpeta…"))
        b.clicked.connect(self._elegir_carpeta)
        f.addWidget(b)
        raiz.addLayout(f)

        cols = QtWidgets.QHBoxLayout()
        izq, der = QtWidgets.QVBoxLayout(), QtWidgets.QVBoxLayout()
        cols.addLayout(izq, 1)
        cols.addLayout(der, 1)
        raiz.addLayout(cols)
        izq.addWidget(self._construir_analisis(False))
        izq.addWidget(self._construir_solver())
        der.addWidget(self._construir_motores())
        der.addWidget(self._construir_malla(v, max(
            int(v.spin_res_fe.value()), 12)))
        self.lab_total = QtWidgets.QLabel()
        raiz.addWidget(self.lab_total)
        bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Cancel)
        self.b_ok = bb.addButton(_("Empezar"),
                                 QtWidgets.QDialogButtonBox.AcceptRole)
        bb.accepted.connect(self._aceptar)
        bb.rejected.connect(self.reject)
        raiz.addWidget(bb)
        self.cmb_prot.currentIndexChanged.connect(self._cambiar_protocolo)
        for c in self.findChildren(QtWidgets.QAbstractButton):
            if isinstance(c, (QtWidgets.QCheckBox, QtWidgets.QRadioButton)):
                c.toggled.connect(self._actualizar)
        # El motor del panel viene marcado.
        mp = getattr(v, "motor_fe", lambda: None)()
        if mp in self.chk_mot and self.chk_mot[mp].isEnabled():
            self.chk_mot[mp].setChecked(True)
        self._cambiar_protocolo()

    def _cambiar_protocolo(self, *_a):
        clave = self.cmb_prot.currentData()
        # El ensayo de la app toma el apoyo del panel: si no es el publicado,
        # el registro lo declara como modificado.
        from spinpy.resistencia import APOYOS
        self._cambios = {}
        if clave == "app":
            ap = APOYOS[max(0, self.v.cmb_apoyo.currentIndex())]
            if ap != febio.PROTOCOLOS_FEM["app"]["apoyo"]:
                self._cambios["apoyo"] = ap
        self._actualizar()

    def _parametros(self):
        clave = self.cmb_prot.currentData()
        d = DialogoParametrosProtocolo(self, clave, self._cambios)
        if d.exec_() == QtWidgets.QDialog.Accepted:
            self._cambios = d.cambios()
            self._actualizar()

    def _elegir_carpeta(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(
            self, _("Carpeta de resultados FEM"), self.ed_carpeta.text())
        if d:
            self.ed_carpeta.setText(d)
            self._actualizar()

    def protocolo(self):
        return febio.protocolo(self.cmb_prot.currentData(), **self._cambios)

    def _actualizar(self, *_a):
        p = self.protocolo()
        self.lab_prot.setText(descripcion_protocolo(p) + (
            "  —  " + _("modificado: «{n}»").format(n=p["nombre"])
            if p["modificado"] else ""))
        i = self.v.cmb_eje_fe.currentIndex()
        self.intro.setText(self._cabecera_febio() + " " + _(
            "la estructura activa y el VOI, con la dirección del panel "
            "(<b>{eje}</b>).").format(eje=html.escape(
                self.v.cmb_eje_fe.itemText(i))))
        for c in self.chk_an.values():
            c.setEnabled(p["tipo"] == "compresion")
        self.b_ok.setEnabled(bool(self.motores_sel())
                             and bool(self.ed_carpeta.text())
                             and (p["tipo"] != "compresion"
                                  or bool(self.analisis())))

    def _aceptar(self):
        self.accept()

    def opciones(self):
        i = self.v.cmb_eje_fe.currentIndex()
        ejes = [0, 1, 2] if i == 3 else [(2,), (0,), (1,)][i]
        return {**self.comunes(), "protocolo": self.protocolo(),
                "ejes": list(ejes), "carpeta": self.ed_carpeta.text()}


# ---------------------------------------------------------------------------
# Resultados
# ---------------------------------------------------------------------------

def _estado_registro(reg):
    """Peor estado de citabilidad de un registro (informe.comprobar)."""
    from spinpy import informe
    items = informe.comprobar({"resultados": {"febio": {"registros": [reg]}},
                               "procedencia": {"version_formato": 99}})
    peor = "citable"
    orden = {"citable": 0, "reservas": 1, "no_citable": 2}
    for it in items:
        if it["magnitud"] in ("vm_max", "vm_p99"):
            continue                    # no se citan nunca: no marcan la fila
        if orden[it["estado"]] > orden[peor]:
            peor = it["estado"]
    return peor


def tabla_html(registros):
    marcas = {"reservas": "*", "no_citable": "†", "citable": ""}

    def c(x, fmt="{:.4g}", pct=False):
        if x is None or (isinstance(x, float) and not np.isfinite(x)):
            return "—"
        return f"{100 * x:+.2f} %" if pct else fmt.format(x)

    t = ["<table cellpadding=4 cellspacing=0 border=1 "
         "style='border-collapse:collapse'>",
         "<tr><th>" + _("estructura") + "</th><th>" + _("protocolo")
         + "</th><th>" + _("malla") + "</th><th>" + _("eje") + "</th><th>"
         + _("motor") + "</th>"
         "<th>E<sub>app</sub> app [MPa]</th><th>E<sub>app</sub> "
         + _("motor") + " [MPa]</th><th>" + _("Δ implementación") + "</th>"
         "<th>p99 " + _("superficie") + " [MPa]</th><th>σ<sub>fallo</sub> "
         "[MPa]</th><th>" + _("NL fuerza") + "</th><th>" + _("NL plato")
         + "</th><th>" + _("NL Pistoia") + "</th><th>BV/TV "
         + _("malla") + "</th><th>" + _("Δ volumen") + "</th></tr>"]
    for r in registros:
        if r.get("tipo") == "homogeneizacion":
            C = r.get("C_periodico", r.get("C_KUBC"))
            txt = (_("periódico") if r.get("malla") == "hex8"
                   else _("cotas KUBC/SUBC"))
            t.append(f"<tr><td>{html.escape(str(r.get('estructura')))}</td>"
                     f"<td>{html.escape(str(r.get('nombre')))}</td>"
                     f"<td>{r.get('malla')}</td><td colspan=12>{txt}"
                     + (f"; C<sub>33</sub> = {np.asarray(C)[2][2] / 1e6:.4g} "
                        "MPa" if C is not None else "")
                     + (f"; Δ app {r['dC_app_rel']:.1e}"
                        if r.get("dC_app_rel") is not None else "")
                     + "</td></tr>")
            continue
        f = febio.fila_tabla(r)
        m = marcas[_estado_registro(r)]
        t.append(
            f"<tr><td>{html.escape(str(f['estructura']))}</td>"
            f"<td>{html.escape(str(f['protocolo']))}</td>"
            f"<td>{f['malla']}</td><td>{f['eje']}</td>"
            f"<td>{html.escape(str(f.get('motor') or '—'))}</td>"
            f"<td>{c(f.get('E_app_app_MPa'))}</td>"
            f"<td>{c(f.get('E_app_MPa'))}{m}</td>"
            f"<td>{c(f.get('dE_app_implementacion'), pct=True)}</td>"
            f"<td>{c(f.get('vm_p99_sup_MPa'))}</td>"
            f"<td>{c(f.get('sigma_fallo_MPa'))}</td>"
            f"<td>{c(f.get('dE_nl_fuerza'), pct=True)}</td>"
            f"<td>{c(f.get('dE_nl_plato'), pct=True)}</td>"
            f"<td>{c(f.get('dE_nl_pistoia'), pct=True)}</td>"
            f"<td>{c(f.get('BVTV_malla'), '{:.4f}')}</td>"
            f"<td>{c(None if f.get('perdida_volumen_pct') is None else -f['perdida_volumen_pct'] / 100, pct=True)}</td></tr>")
    t.append("</table>")
    t.append("<p style='color:#666'>" + _(
        "* con reservas; † no citable (ver el informe de publicación). "
        "E<sub>app</sub> del motor con el desplazamiento del techo ponderado "
        "por área; «Δ implementación» compara el motor con ladrillos y la app "
        "con la definición de la app y debe ser ~0. Los desvíos no lineales "
        "son frente al lineal de la misma malla.") + "</p>")
    t.append(tabla_motores_html(registros))
    return "".join(t)


def tabla_motores_html(registros):
    """Tabla comparativa ENTRE MOTORES sobre la misma malla
    (`fem.tabla_motores`). Vacia si solo hubo un motor por malla."""
    filas = febio.tabla_motores(registros)
    grupos = {}
    for f in filas:
        grupos.setdefault((f["estructura"], f["protocolo"], f["malla"],
                           f["eje"]), []).append(f)
    if not any(len(g) > 1 for g in grupos.values()):
        return ""

    def c(x, fmt="{:.4g}"):
        if x is None or (isinstance(x, float) and not np.isfinite(x)):
            return "—"
        return fmt.format(x)

    t = ["<p><b>" + _("Comparación entre motores (misma malla)") + "</b></p>",
         "<table cellpadding=4 cellspacing=0 border=1 "
         "style='border-collapse:collapse'><tr><th>" + _("estructura")
         + "</th><th>" + _("malla") + "</th><th>" + _("eje") + "</th><th>"
         + _("motor") + "</th><th>" + _("resolvedor") + "</th>"
         "<th>E<sub>app</sub> [MPa]</th><th>p99 " + _("superficie")
         + " [MPa]</th><th>σ<sub>fallo</sub> [MPa]</th><th>ΔE<sub>app</sub>"
         "</th><th>Δp99</th><th>Δσ<sub>fallo</sub></th><th>"
         + _("tiempo") + " [s]</th><th>" + _("memoria") + " [MB]</th></tr>"]
    for g in grupos.values():
        if len(g) < 2:
            continue
        for f in g:
            if f.get("no_disponible"):
                t.append(f"<tr><td>{html.escape(str(f['estructura']))}</td>"
                         f"<td>{f['malla']}</td><td>{f['eje']}</td>"
                         f"<td>{html.escape(str(f['motor']))}</td>"
                         "<td colspan=9 style='color:#888'>"
                         + _("no disponible con esta malla o análisis")
                         + "</td></tr>")
                continue
            t.append(
                f"<tr><td>{html.escape(str(f['estructura']))}</td>"
                f"<td>{f['malla']}</td><td>{f['eje']}</td>"
                f"<td>{html.escape(str(f['motor']))}</td>"
                f"<td>{html.escape(str(f.get('solver') or '—'))}</td>"
                f"<td>{c(f.get('E_app_MPa'))}</td>"
                f"<td>{c(f.get('vm_p99_sup_MPa'))}</td>"
                f"<td>{c(f.get('sigma_fallo_MPa'))}</td>"
                f"<td>{c(f.get('dE_rel'), '{:.1e}')}</td>"
                f"<td>{c(f.get('dp99_rel'), '{:.1e}')}</td>"
                f"<td>{c(f.get('dfallo_rel'), '{:.1e}')}</td>"
                f"<td>{c(f.get('tiempo_s'), '{:.1f}')}</td>"
                f"<td>{c(f.get('memoria_MB'), '{:.0f}')}</td></tr>")
    t.append("</table><p style='color:#666'>" + _(
        "Diferencias relativas frente al primer motor de cada grupo. Con la "
        "misma malla, cargas y postproceso miden solo la implementación: "
        "lo esperable es del orden de la tolerancia del resolvedor (1e-10 a "
        "1e-8).") + "</p>")
    return "".join(t)


class DialogoResultadosFEBio(QtWidgets.QDialog):
    def __init__(self, v, registros, carpeta, pendiente=False, mapas=None,
                 figura=None):
        super().__init__(v)
        self.v = v
        self.registros = registros
        self.carpeta = Path(carpeta) if carpeta else None
        self.mapas = mapas or []
        self.setWindowTitle(_("Resultados FEM"))
        self.resize(1200, 620)
        lay = QtWidgets.QVBoxLayout(self)
        vistos = []
        for r in registros:
            m = r.get("motor") or {}
            et = " ".join(str(x) for x in (m.get("nombre"), m.get("version"))
                          if x)
            if et and et not in vistos:
                vistos.append(et)
        cab = QtWidgets.QLabel(_("{m} · {n} registro(s) · carpeta "
                                 "{c}").format(m=", ".join(vistos) or "—",
                                               n=len(registros),
                                               c=self.carpeta or "—"))
        cab.setWordWrap(True)
        lay.addWidget(cab)
        tb = QtWidgets.QTextBrowser()
        txt = tabla_html(registros)
        fallos = [(r.get("estructura"), (r.get("motor") or {}).get("nombre"),
                   r.get("nombre"), r.get("malla"), f.get("corrida"),
                   f.get("msg")) for r in registros
                  for f in r.get("fallos", [])]
        if fallos:
            txt += "<p><b>" + _("Resoluciones que fallaron o no están "
                                "disponibles") + "</b><br>" + \
                "<br>".join(html.escape(" · ".join(str(x) for x in f))
                            for f in fallos) + "</p>"
        if figura:
            txt += f"<p><img src='{Path(figura).as_uri()}' width=900></p>"
        tb.setHtml(txt)
        lay.addWidget(tb, 1)
        f = QtWidgets.QHBoxLayout()
        b = QtWidgets.QPushButton(_("Mapas 3D de von Mises…"))
        b.setEnabled(bool(self.mapas))
        b.clicked.connect(self._mapas)
        f.addWidget(b)
        b = QtWidgets.QPushButton(_("Exportar CSV…"))
        b.clicked.connect(self._csv)
        f.addWidget(b)
        b = QtWidgets.QPushButton(_("Abrir carpeta"))
        b.setEnabled(self.carpeta is not None)
        b.clicked.connect(lambda: QtGui.QDesktopServices.openUrl(
            QtCore.QUrl.fromLocalFile(str(self.carpeta))))
        f.addWidget(b)
        self.b_incluir = QtWidgets.QPushButton(
            _("Incluir en el informe de publicación"))
        self.b_incluir.setEnabled(pendiente)
        self.b_incluir.clicked.connect(self._incluir)
        f.addWidget(self.b_incluir)
        f.addStretch(1)
        b = QtWidgets.QPushButton(_("Cerrar"))
        b.clicked.connect(self.accept)
        f.addWidget(b)
        lay.addLayout(f)

    def _incluir(self):
        self.v._febio_registrar(self.registros)
        self.b_incluir.setEnabled(False)

    def _csv(self):
        r, _f = QtWidgets.QFileDialog.getSaveFileName(
            self, _("Exportar tabla CSV"),
            str((self.carpeta or Path.home()) / "tabla_fem.csv"), "CSV (*.csv)")
        if r:
            escribir_csv(self.registros, r)

    def _mapas(self):
        DialogoMapasFEBio(self, self.mapas).exec_()


def escribir_csv(registros, ruta):
    """Tabla de `fila_tabla`, una fila por registro, con procedencia."""
    import csv

    from spinpy import procedencia
    filas = [febio.fila_tabla(r) for r in registros
             if r.get("tipo") != "homogeneizacion"]
    claves = []
    for f in filas:
        for k in f:
            if k not in claves:
                claves.append(k)
    cols = list(procedencia.COLUMNAS) + ["motor_version", "solver",
                                         "huella_malla"] + claves
    with open(ruta, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r, f in zip([r for r in registros
                         if r.get("tipo") != "homogeneizacion"], filas):
            p = r.get("procedencia") or {}
            fila = {k: p.get(k) for k in procedencia.COLUMNAS}
            mt = r.get("motor") or {}
            fila.update(f, motor_version=mt.get("version")
                        or (r.get("febio") or {}).get("version"),
                        solver=mt.get("solver"),
                        huella_malla=r.get("huella_malla"))
            w.writerow(fila)


class DialogoMapasFEBio(QtWidgets.QDialog):
    """von Mises por elemento (lineal, a la carga del protocolo), con UNA
    escala de color para todos los paneles; la de p99 del conjunto, porque el
    maximo no converge y aplastaria el resto."""

    def __init__(self, padre, mapas):
        super().__init__(padre)
        import pyvista as pv
        from pyvistaqt import QtInteractor
        self.setWindowTitle(_("Mapas 3D de von Mises"))
        self.resize(1200, 700)
        lay = QtWidgets.QGridLayout(self)
        todos = np.concatenate([m["vm"] for m in mapas]) if mapas else [1.0]
        lim = (0.0, float(np.percentile(todos, 99)))
        self._vistas = []
        cols = min(3, max(1, len(mapas)))
        for i, m in enumerate(mapas[:9]):
            caja = QtWidgets.QWidget()
            cl = QtWidgets.QVBoxLayout(caja)
            cl.setContentsMargins(0, 0, 0, 0)
            cl.addWidget(QtWidgets.QLabel(m["titulo"]))
            vista = QtInteractor(caja)
            vista.set_background("white")
            cl.addWidget(vista.interactor, 1)
            e = np.asarray(m["elems"])
            lin = e[:, :4] if e.shape[1] == 10 else e
            tipo = pv.CellType.TETRA if e.shape[1] == 10 else \
                pv.CellType.HEXAHEDRON
            celdas = np.hstack([np.full((lin.shape[0], 1), lin.shape[1]),
                                lin]).ravel()
            g = pv.UnstructuredGrid(celdas, np.full(lin.shape[0], tipo),
                                    np.asarray(m["nodos"], float))
            g.cell_data["von Mises [MPa]"] = np.asarray(m["vm"], float)
            vista.add_mesh(g.extract_surface(), scalars="von Mises [MPa]",
                           cmap="inferno", clim=lim, show_scalar_bar=(i == 0))
            vista.view_isometric()
            vista.enable_parallel_projection()
            self._vistas.append(vista)
            lay.addWidget(caja, i // cols, i % cols)

    def closeEvent(self, ev):
        for v in self._vistas:
            v.close()
        super().closeEvent(ev)


def figura_validacion(registros, ruta):
    """Cinco columnas: E_app, p99 de superficie, tension de fallo, desvio no
    lineal (fuerza) y BV/TV de la malla; barras app y una por cada (malla,
    motor), por estructura y protocolo. Devuelve la ruta o None."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    filas = [febio.fila_tabla(r) for r in registros
             if r.get("tipo") != "homogeneizacion" and r.get("lineal")]
    if not filas:
        return None
    grupos = []
    for f in filas:
        g = (f["estructura"], f["protocolo"], f["eje"])
        if g not in grupos:
            grupos.append(g)
    paneles = [("E_app_MPa", "E_app_app_MPa", "E_app [MPa]"),
               ("vm_p99_sup_MPa", "vm_p99_sup_app_MPa", "p99 σ_vM sup. [MPa]"),
               ("sigma_fallo_MPa", "sigma_fallo_app_MPa", "σ_fallo [MPa]"),
               ("dE_nl_fuerza", None, "ΔE no lineal [%]"),
               ("BVTV_malla", None, "BV/TV malla")]
    combos = []
    for f in filas:
        k = (f["malla"], f.get("motor"))
        if k not in combos:
            combos.append(k)
    paleta = ["#1f4e9c", "#b06000", "#2e7d32", "#7b1fa2", "#00838f",
              "#c62828", "#5d4037", "#455a64", "#9e9d24"]
    series = [("app", None, "#8c8c8c")] + [
        (f"{mo} {ma}", (ma, mo), paleta[i % len(paleta)])
        for i, (ma, mo) in enumerate(combos)]
    ancho = 0.8 / len(series)
    fig, axs = plt.subplots(1, 5, figsize=(17, 3.6))
    x = np.arange(len(grupos))
    for ax, (k, k_app, titulo) in zip(axs, paneles):
        for j, (etq, combo, col) in enumerate(series):
            vals = []
            for g in grupos:
                ff = [f for f in filas if (f["estructura"], f["protocolo"],
                                           f["eje"]) == g]
                if combo is None:
                    v_ = next((f.get(k_app) for f in ff if k_app
                               and f.get(k_app) is not None), None)
                else:
                    v_ = next((f.get(k) for f in ff
                               if (f["malla"], f.get("motor")) == combo), None)
                if v_ is not None and k == "dE_nl_fuerza":
                    v_ = 100 * v_
                vals.append(np.nan if v_ is None else v_)
            ax.bar(x + (j - (len(series) - 1) / 2) * ancho, vals, ancho,
                   color=col, label=etq)
        ax.set_title(titulo, fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{g[0]}\n{g[1]} {g[2]}" for g in grupos],
                           fontsize=7)
    axs[0].legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(ruta, dpi=200)
    plt.close(fig)
    return ruta
