"""dialogo_metodos.py — Ajustar con todos los metodos, comparar y elegir.

FLUJO
-----
    1. una fila por metodo, con su barra de progreso, corriendo a la vez
    2. al terminar, una columna por metodo: vista previa, metricas frente al
       VOI y el error
    3. el usuario elige UNO; el resto se descarta
    4. con el elegido: seguir trabajando, o ademas generar N replicas
       —tecnicas o con variacion en las metricas que se marquen

POR QUE HILOS Y NO PROCESOS
----------------------------
Cada metodo corre en su propio QThread. Con procesos habria paralelismo real,
pero cada uno tendria que recibir el VOI por copia —decenas de MB— y devolver
la mascara ganadora por el mismo camino, y el progreso por metodo exigiria una
tuberia aparte. Con hilos el reparto es peor (numpy libera el GIL solo a
ratos) pero el progreso llega solo y la memoria se comparte, que en esta
maquina importa mas que los segundos.

Consecuencia honesta: **los cuatro metodos juntos no tardan la cuarta parte de
lo que tardarian en serie**. Tardan algo menos, y sobre todo se ven avanzar.
"""

from __future__ import annotations

import sys
import traceback
from pathlib import Path

import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spinpy.idioma import _                                 # noqa: E402
from spinpy.metodos import (COMPARADAS, REGISTRO,  # noqa: E402
                            comparar_resultados, generar_replicas,
                            gobierna_de)


class HiloMetodo(QtCore.QThread):
    """Un metodo de ajuste corriendo aparte, informando de su avance."""

    avance = QtCore.pyqtSignal(str, int, int, str)   # metodo, i, n, etapa
    listo = QtCore.pyqtSignal(str, object)
    fallo = QtCore.pyqtSignal(str, str)

    def __init__(self, clave, VOI, spacing, m_voi, kw, registro=None):
        super().__init__()
        self.clave = clave
        self.registro = registro if registro is not None else REGISTRO
        self._args = (VOI, spacing, m_voi, kw)

    def run(self):
        VOI, spacing, m_voi, kw = self._args
        try:
            def pr(i, n, etapa):
                self.avance.emit(self.clave, int(i), int(n), str(etapa))
            r = self.registro[self.clave]["correr"](
                VOI, spacing, m_voi=m_voi, progreso=pr, **kw)
            self.listo.emit(self.clave, r)
        except Exception:
            self.fallo.emit(self.clave, traceback.format_exc())


class DialogoMetodos(QtWidgets.QDialog):
    """Ventana de ajuste comparado. Devuelve el resultado elegido."""

    def __init__(self, VOI, spacing, m_voi, parent=None, semilla=20260720,
                 num_waves=700, claves=None, registro=None):
        super().__init__(parent)
        self.setWindowTitle(_("Ajustar con todos los metodos"))
        self.resize(1180, 780)
        self.VOI, self.spacing, self.m_voi = VOI, spacing, m_voi
        # El registro de la FAMILIA que se ajusta. Por omision el del
        # spinodoide, que es lo que esta ventana hizo siempre.
        self.registro = registro if registro is not None else REGISTRO
        self.claves = list(claves or self.registro.keys())
        self.resultados = {}
        self.elegido = None
        self.replicas = None
        self._hilos = {}

        self.pila = QtWidgets.QStackedWidget()
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(self.pila)

        self.pila.addWidget(self._pagina_progreso())
        self.pila.addWidget(QtWidgets.QWidget())     # la comparativa, al final

        for c in self.claves:
            h = HiloMetodo(c, VOI, spacing, m_voi,
                           {"semilla": semilla, "num_waves": num_waves},
                           registro=self.registro)
            h.avance.connect(self._avance)
            h.listo.connect(self._listo)
            h.fallo.connect(self._fallo)
            self._hilos[c] = h
            h.start()

    # -- fase 1: progreso ---------------------------------------------------

    def _pagina_progreso(self):
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        t = QtWidgets.QLabel(
            _("<b>Ajustando el VOI con {n} metodos a la vez.</b>").format(
                n=len(self.claves)))
        lay.addWidget(t)
        nota = QtWidgets.QLabel(_(
            "Corren en hilos y comparten memoria, asi que el reparto de CPU "
            "no es perfecto: juntos tardan algo menos que en serie, no la "
            "cuarta parte. El desempate mecanico es el mas lento con "
            "diferencia, porque homogeneiza cada finalista."))
        nota.setWordWrap(True)
        nota.setStyleSheet("color:#666; font-size:10px;")
        lay.addWidget(nota)

        self.barras, self.etiquetas = {}, {}
        rej = QtWidgets.QGridLayout()
        for i, c in enumerate(self.claves):
            m = self.registro[c]
            n = QtWidgets.QLabel("<b>%s</b>" % m["etiqueta"])
            n.setToolTip(m["descripcion"])
            b = QtWidgets.QProgressBar()
            b.setRange(0, max(1, m["n_eval"]))
            b.setFormat("%p%")
            e = QtWidgets.QLabel(_("en cola"))
            e.setStyleSheet("color:#666;")
            rej.addWidget(n, i, 0)
            rej.addWidget(b, i, 1)
            rej.addWidget(e, i, 2)
            rej.setColumnStretch(1, 1)
            self.barras[c], self.etiquetas[c] = b, e
        lay.addLayout(rej)
        lay.addStretch(1)

        self.btn_cancelar = QtWidgets.QPushButton(_("Cancelar"))
        self.btn_cancelar.clicked.connect(self.reject)
        f = QtWidgets.QHBoxLayout()
        f.addStretch(1)
        f.addWidget(self.btn_cancelar)
        lay.addLayout(f)
        return w

    def _avance(self, clave, i, n, etapa):
        b = self.barras[clave]
        if n and b.maximum() != n:
            b.setRange(0, n)
        b.setValue(i)
        self.etiquetas[clave].setText("etapa %s — %d/%d" % (etapa, i, n))

    def _listo(self, clave, r):
        self.resultados[clave] = r
        b = self.barras[clave]
        b.setValue(b.maximum())
        self.etiquetas[clave].setText(
            "<span style='color:#1a7f37'>listo — error %.5f, %.0f s</span>"
            % (r["error"], r["tiempo_s"]))
        self._quizas_terminar()

    def _fallo(self, clave, tb):
        self.resultados[clave] = {"metodo": clave, "fallo": tb.strip().split("\n")[-1]}
        self.etiquetas[clave].setText(
            "<span style='color:#b62324'>fallo</span>")
        self.barras[clave].setValue(0)
        self._quizas_terminar()

    def _quizas_terminar(self):
        if len(self.resultados) < len(self.claves):
            return
        ok = [r for r in self.resultados.values() if not r.get("fallo")]
        if not ok:
            QtWidgets.QMessageBox.critical(
                self, _("Ningun metodo termino"),
                _("Los {n} metodos fallaron. El primero dice:\n\n{msg}"
                  ).format(n=len(self.claves),
                           msg=next(iter(self.resultados.values())).get(
                               "fallo", "")))
            self.reject()
            return
        self.pila.removeWidget(self.pila.widget(1))
        self.pila.addWidget(self._pagina_comparativa())
        self.pila.setCurrentIndex(1)

    # -- fase 2: comparativa y seleccion ------------------------------------

    def _vista_previa(self, BW, spacing, lado=210):
        """Miniatura de la estructura. Render 3D si se puede, corte si no.

        El render fuera de pantalla necesita un contexto OpenGL, que no existe
        en todas las maquinas ni en modo headless. Cuando falta, un corte
        central sigue diciendo lo esencial —si la estructura es fina o gruesa,
        abierta o cerrada— y es preferible a un hueco.
        """
        try:
            import pyvista as pv
            from visor import malla_de_mascara
            malla = malla_de_mascara(BW, spacing)
            p = pv.Plotter(off_screen=True, window_size=(lado, lado))
            p.set_background("white")
            p.add_mesh(malla, color="#d9c7a0", smooth_shading=True)
            p.camera_position = "iso"
            img = p.screenshot(return_img=True)
            p.close()
            h, w, _ = img.shape
            qi = QtGui.QImage(np.ascontiguousarray(img).data, w, h, 3 * w,
                              QtGui.QImage.Format_RGB888).copy()
            return QtGui.QPixmap.fromImage(qi)
        except Exception:
            z = BW.shape[2] // 2
            corte = np.ascontiguousarray(
                np.where(BW[:, :, z].T, 40, 245).astype(np.uint8))
            h, w = corte.shape
            qi = QtGui.QImage(corte.data, w, h, w,
                              QtGui.QImage.Format_Grayscale8).copy()
            return QtGui.QPixmap.fromImage(qi).scaled(
                lado, lado, QtCore.Qt.KeepAspectRatio,
                QtCore.Qt.SmoothTransformation)

    def _pagina_comparativa(self):
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        lay.addWidget(QtWidgets.QLabel(
            _("<b>Elige con que ajuste seguir.</b> El resto se descarta.")))

        filas = comparar_resultados(
            [self.resultados[c] for c in self.claves if c in self.resultados],
            self.m_voi)
        por_metodo = {f["metodo"]: f for f in filas}
        mejor = min((f for f in filas if "error" in f),
                    key=lambda f: f["error"], default=None)

        cont = QtWidgets.QWidget()
        rej = QtWidgets.QHBoxLayout(cont)
        self.grupo = QtWidgets.QButtonGroup(self)
        for c in self.claves:
            r = self.resultados.get(c)
            col = QtWidgets.QGroupBox(self.registro[c]["etiqueta"])
            cl = QtWidgets.QVBoxLayout(col)
            if not r or r.get("fallo"):
                cl.addWidget(QtWidgets.QLabel(
                    "<span style='color:#b62324'>fallo</span><br>"
                    + (r or {}).get("fallo", "")[:120]))
                rej.addWidget(col)
                continue

            img = QtWidgets.QLabel()
            img.setPixmap(self._vista_previa(r["mascara"],
                                             np.full(3, r["parametros"]["spacing_mm"])))
            img.setAlignment(QtCore.Qt.AlignCenter)
            cl.addWidget(img)

            f = por_metodo.get(c, {})
            t = QtWidgets.QTableWidget(len(COMPARADAS) + 1, 3)
            t.setHorizontalHeaderLabels(["", "valor", "vs VOI"])
            t.verticalHeader().setVisible(False)
            t.horizontalHeader().setSectionResizeMode(
                QtWidgets.QHeaderView.Stretch)
            t.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
            t.setItem(0, 0, QtWidgets.QTableWidgetItem("error"))
            t.setItem(0, 1, QtWidgets.QTableWidgetItem("%.5f" % r["error"]))
            t.setItem(0, 2, QtWidgets.QTableWidgetItem(
                "%d term." % f.get("n_terminos", 0) if f.get("mecanico")
                else ""))
            for i, k in enumerate(COMPARADAS, start=1):
                t.setItem(i, 0, QtWidgets.QTableWidgetItem(k))
                v = f.get(k)
                t.setItem(i, 1, QtWidgets.QTableWidgetItem(
                    "%.4g" % v if v is not None and np.isfinite(v) else "—"))
                d = f.get("d_" + k)
                it = QtWidgets.QTableWidgetItem(
                    "%+.1f %%" % d if d is not None else "—")
                if d is not None:
                    a = abs(d)
                    it.setForeground(QtGui.QBrush(QtGui.QColor(
                        "#1a7f37" if a < 5 else "#9a6700" if a < 15
                        else "#b62324")))
                t.setItem(i, 2, it)
            t.setFixedHeight(26 * (len(COMPARADAS) + 1) + 26)
            cl.addWidget(t)

            pie = [_("{t} s · {n} evaluaciones").format(
                t="%.0f" % r["tiempo_s"], n=r.get("n_evaluaciones", 0))]
            if not r["diagnostico"]["bicontinuo"]:
                pie.append("<span style='color:#b62324'>"
                           + _("fase poro fragmentada: no es trabecula")
                           + "</span>")
            if (r.get("mecanico") or {}).get("traza"):
                pie.append("desempate mecanico: "
                           + ("reordeno" if r["mecanico"]["reordeno"]
                              else "confirmo"))
            n = QtWidgets.QLabel("<br>".join(pie))
            n.setWordWrap(True)
            n.setStyleSheet("color:#666; font-size:10px;")
            cl.addWidget(n)

            rb = QtWidgets.QRadioButton("Usar este")
            if mejor and c == mejor["metodo"]:
                rb.setChecked(True)
                rb.setText("Usar este (menor error)")
            rb.clave = c
            self.grupo.addButton(rb)
            cl.addWidget(rb)
            rej.addWidget(col)

        area = QtWidgets.QScrollArea()
        area.setWidget(cont)
        area.setWidgetResizable(True)
        lay.addWidget(area, 1)

        aviso = QtWidgets.QLabel(_(
            "El <b>error no es comparable entre metodos que usan distinto "
            "numero de terminos</b>: el del desempate mecanico incluye "
            "E<sub>z</sub>/E<sub>s</sub> y E<sub>z</sub>/E<sub>x</sub> y los "
            "demas no. Para comparar de igual a igual, mira las diferencias "
            "por metrica."))
        aviso.setWordWrap(True)
        aviso.setTextFormat(QtCore.Qt.RichText)
        aviso.setStyleSheet("color:#8a6d00; font-size:10px;")
        lay.addWidget(aviso)

        f = QtWidgets.QHBoxLayout()
        b = QtWidgets.QPushButton(_("Cancelar"))
        b.clicked.connect(self.reject)
        f.addWidget(b)
        f.addStretch(1)
        b = QtWidgets.QPushButton(_("Usar el seleccionado"))
        b.clicked.connect(self._aceptar_solo)
        f.addWidget(b)
        b = QtWidgets.QPushButton(_("Usar y generar replicas…"))
        b.setStyleSheet("font-weight:bold; padding:6px;")
        b.clicked.connect(self._aceptar_con_replicas)
        f.addWidget(b)
        lay.addLayout(f)
        return w

    def _seleccion(self):
        b = self.grupo.checkedButton()
        if b is None:
            QtWidgets.QMessageBox.information(
                self, _("Sin seleccion"),
                _("Marca con que ajuste quieres seguir."))
            return None
        return self.resultados[b.clave]

    def _aceptar_solo(self):
        r = self._seleccion()
        if r is not None:
            self.elegido = r
            self.accept()

    def _aceptar_con_replicas(self):
        r = self._seleccion()
        if r is None:
            return
        d = DialogoReplicas(r, self.spacing, self)
        if d.exec_() != QtWidgets.QDialog.Accepted:
            return
        self.elegido = r
        self.replicas = d.salida
        self.accept()

    def closeEvent(self, ev):
        for h in self._hilos.values():
            if h.isRunning():
                h.wait(50)
        super().closeEvent(ev)


class DialogoReplicas(QtWidgets.QDialog):
    """Cuantas replicas, y si varian en algo.

    LA DISTINCION QUE ESTA VENTANA TIENE QUE DEJAR CLARA. Replicas identicas
    (misma parametrizacion, distinta semilla) miden el ruido del generador y
    **no son especimenes**: no anaden grados de libertad a ninguna comparacion
    biologica. Replicas con variacion son puntos de diseno de una familia
    sintetica, y tampoco lo son. Confundir cualquiera de las dos con N es
    pseudorreplicacion, y es el error que `spinpy.estadistica` existe para
    impedir.
    """

    def __init__(self, resultado, spacing, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Generar replicas")
        self.resize(560, 520)
        self.resultado, self.spacing = resultado, spacing
        self.salida = None

        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(QtWidgets.QLabel(
            _("<b>{metodo}</b> — error {err}").format(
                metodo=resultado.get("etiqueta", ""),
                err="%.5f" % resultado["error"])))

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(_("Numero de replicas:")))
        self.spin_n = QtWidgets.QSpinBox()
        self.spin_n.setRange(1, 200)
        self.spin_n.setValue(10)
        f.addWidget(self.spin_n, 1)
        lay.addLayout(f)

        self.rb_iguales = QtWidgets.QRadioButton(
            _("Identicas (misma parametrizacion, distinta semilla)"))
        self.rb_iguales.setChecked(True)
        self.rb_variar = QtWidgets.QRadioButton(
            _("Con variacion en las metricas que elija"))
        lay.addWidget(self.rb_iguales)
        lay.addWidget(self.rb_variar)

        self.caja = QtWidgets.QGroupBox(_("Metricas que deben variar"))
        cl = QtWidgets.QVBoxLayout(self.caja)
        self.checks = {}
        familia = (resultado.get("parametros") or {}).get("familia")
        for k, (par, expl) in gobierna_de(familia).items():
            c = QtWidgets.QCheckBox("%s  —  %s" % (k, expl))
            c.setEnabled(False)
            self.checks[k] = c
            cl.addWidget(c)
        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(
            _("Amplitud (desviacion relativa):")))
        self.spin_amp = QtWidgets.QDoubleSpinBox()
        self.spin_amp.setRange(0.01, 0.50)
        self.spin_amp.setSingleStep(0.01)
        self.spin_amp.setValue(0.05)
        self.spin_amp.setEnabled(False)
        f.addWidget(self.spin_amp, 1)
        cl.addLayout(f)
        lay.addWidget(self.caja)

        self.rb_variar.toggled.connect(self._alternar)

        nota = QtWidgets.QLabel(_(
            "<b>Ni las identicas ni las variadas son especimenes.</b> Las "
            "identicas miden el ruido del generador; las variadas son puntos "
            "de diseno de una familia sintetica. Tratar N replicas como N "
            "especimenes infla los grados de libertad y estrecha los "
            "intervalos de confianza: es pseudorreplicacion.<br><br>"
            "La amplitud se aplica al PARAMETRO que gobierna cada metrica, no "
            "a la metrica: el generador toma densidad, numero de onda y "
            "angulos, y las metricas salen de ahi. La dispersion realmente "
            "lograda se reporta al terminar, y no tiene por que coincidir con "
            "la pedida."))
        nota.setWordWrap(True)
        nota.setTextFormat(QtCore.Qt.RichText)
        nota.setStyleSheet("color:#666; font-size:10px;")
        lay.addWidget(nota)
        lay.addStretch(1)

        self.barra = QtWidgets.QProgressBar()
        self.barra.hide()
        lay.addWidget(self.barra)

        f = QtWidgets.QHBoxLayout()
        b = QtWidgets.QPushButton(_("Cancelar"))
        b.clicked.connect(self.reject)
        f.addWidget(b)
        f.addStretch(1)
        self.btn_ok = QtWidgets.QPushButton("Generar")
        self.btn_ok.setStyleSheet("font-weight:bold; padding:6px;")
        self.btn_ok.clicked.connect(self._generar)
        f.addWidget(self.btn_ok)
        lay.addLayout(f)

    def _alternar(self, v):
        for c in self.checks.values():
            c.setEnabled(v)
        self.spin_amp.setEnabled(v)

    def _generar(self):
        variar = ([k for k, c in self.checks.items() if c.isChecked()]
                  if self.rb_variar.isChecked() else [])
        if self.rb_variar.isChecked() and not variar:
            QtWidgets.QMessageBox.information(
                self, _("Sin metricas"),
                _("Marca al menos una metrica que variar, o elige replicas "
                  "identicas."))
            return
        n = int(self.spin_n.value())
        self.btn_ok.setEnabled(False)
        self.barra.setRange(0, n)
        self.barra.show()

        def pr(i, tot, msg):
            self.barra.setValue(i)
            QtWidgets.QApplication.processEvents()

        try:
            self.salida = generar_replicas(
                self.resultado["parametros"], self.spacing, n=n,
                variar=variar, amplitud=float(self.spin_amp.value()),
                progreso=pr)
        except Exception:
            QtWidgets.QMessageBox.critical(self, _("Error al generar"),
                                           traceback.format_exc())
            self.btn_ok.setEnabled(True)
            self.barra.hide()
            return
        self.accept()
