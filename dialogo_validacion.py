# -*- coding: utf-8 -*-
"""
dialogo_validacion.py — Ventana «Validación»: las réplicas de los resultados
                        publicados, cada una al lado de su original.

POR QUE ESTA VENTANA EXISTE
---------------------------
El método que implementa esta aplicación no es nuestro. Es el de

    Kumar, S., Tan, S., Zheng, L., & Kochmann, D. M. (2020). Inverse-designed
        spinodoid metamaterials. npj Computational Materials, 6(1), 1-10.
        Articulo 73. https://doi.org/10.1038/s41524-020-0341-6

y la pregunta que cualquiera hace la primera vez —«¿de donde sale esto y como
se que esta bien implementado?»— merece una respuesta que se pueda mirar, no
una lista de pruebas internas. Aqui se pone la figura publicada a la izquierda
y la que genera este programa a la derecha, con las comprobaciones numericas
que tendrian que fallar si la implementacion se hubiera desviado del articulo.

UNA PESTANA POR ARTICULO
------------------------
Son tres replicas y cada una comprueba algo distinto, asi que no se mezclan en
una sola tabla: quien mire tiene que poder decir de que articulo sale cada
afirmacion.

  Kumar et al. (2020)   el metodo en si: las cuatro clases de anisotropia y el
                        muestreo por rechazo de la ecuacion (2).
  Zheng et al. (2021)   las COTAS: cada superficie elastica dentro de Voigt y
                        de Hashin-Shtrikman, la ortotropia del tensor, y por
                        que el articulo impone rho >= 0.3.
  Guo et al. (2024)     las CURVATURAS: el perfil (k1, k2) de un spinodoide de
                        parametros publicados y el de una superficie nodal
                        periodica, cuya respuesta es cerrada.

La logica de cada replica NO esta aqui: esta en `Test/replicar_*.py`, que se
pueden ejecutar solos, sin interfaz, y dejan su informe en JSON. Esta ventana
unicamente los lanza y muestra lo que produjeron. Anadir una cuarta replica es
anadir una entrada a `REPLICAS` y su guion; no hay que tocar la ventana.

LICENCIA DE LAS FIGURAS REPRODUCIDAS
------------------------------------
Los tres articulos son de acceso abierto bajo Creative Commons Attribution
(CC BY 4.0), https://creativecommons.org/licenses/by/4.0/. Esa licencia es la
que permite incluir sus figuras aqui, sin modificar y con atribucion. Ver
`Test/referencia/FUENTES.md`, donde consta la procedencia exacta de cada una.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from PyQt5 import QtCore, QtGui, QtWidgets

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spinpy.idioma import _                                   # noqa: E402

# ---------------------------------------------------------------------------
# Las replicas disponibles
# ---------------------------------------------------------------------------
# Cada entrada se basta a si misma: quien lea SOLO este archivo tiene que
# encontrar la cita completa, no una referencia a otro sitio. Por eso las citas
# estan repetidas aqui aunque tambien esten en la cabecera de cada guion.
REPLICAS = [
    {
        "clave": "kumar2020",
        "pestana": "Kumar 2020",
        "opcion": "--replicar-kumar2020",
        "guion": "replicar_kumar2020.py",
        "figura_ref": "Kumar2020_Fig2.png",
        "pie_izq": "Original publicado (Fig. 2)",
        "json": "replica_kumar2020.json",
        "figura": "replica_kumar2020_fig2.png",
        "cita": (
            "Kumar, S., Tan, S., Zheng, L., &amp; Kochmann, D. M. (2020). "
            "Inverse-designed spinodoid metamaterials. <i>npj Computational "
            "Materials, 6</i>(1), 1-10. Articulo 73. "
            "https://doi.org/10.1038/s41524-020-0341-6"),
        "licencia": (
            "Figura original bajo "
            "<a href='https://creativecommons.org/licenses/by/4.0/'>CC BY "
            "4.0</a> &mdash; &copy; The Author(s) 2020. Reproducida sin "
            "modificaciones."),
        "cabecera": (
            "<b>Este método no es nuestro: es el de Kumar et al. (2020).</b> "
            "Lo que sigue es su figura publicada al lado de la misma figura "
            "regenerada con este programa, y las comprobaciones numéricas que "
            "fallarían si nuestra implementación se hubiera desviado del "
            "artículo."),
        "params": ("ρ = {rho} · β = {beta_pi}π · N = {num_waves} · vista "
                   "{resolucion_vista}³ · homogeneización "
                   "{resolucion_homogeneizacion}³ · semilla {semilla} · "
                   "muestreo: {esquema}"),
        "ayuda_rapido": (
            "Baja la malla de vista a 96³ y la de homogeneización a 16³. La "
            "réplica sale en algo más de un minuto en vez de cuatro, y las "
            "seis comprobaciones se siguen pasando; lo que se pierde es "
            "definición en la figura."),
    },
    {
        "clave": "zheng2021",
        "pestana": "Zheng 2021",
        "opcion": "--replicar-zheng2021",
        "guion": "replicar_zheng2021.py",
        "figura_ref": "Zheng2021_Fig1.png",
        "pie_izq": "Original publicado (Fig. 1)",
        "json": "replica_zheng2021.json",
        "figura": "replica_zheng2021_fig1.png",
        "cita": (
            "Zheng, L., Kumar, S., &amp; Kochmann, D. M. (2021). Data-driven "
            "topology optimization of spinodoid metamaterials with seamlessly "
            "tunable anisotropy. <i>Computer Methods in Applied Mechanics and "
            "Engineering, 383</i>, 113894. "
            "https://doi.org/10.1016/j.cma.2021.113894"),
        "licencia": (
            "Figura original bajo "
            "<a href='http://creativecommons.org/licenses/by/4.0/'>CC BY</a> "
            "&mdash; &copy; 2021 The Author(s), Elsevier B.V. Reproducida sin "
            "modificaciones."),
        "cabecera": (
            "<b>Aquí se comprueban COTAS, no proporciones.</b> Su Figura 1 "
            "dibuja cada superficie elástica dentro de la esfera de Voigt y, "
            "en la isótropa, de la de Hashin-Shtrikman. Son cotas con fórmula "
            "cerrada: un error de escala en nuestra homogeneización no las "
            "pasa. Se comprueba además la ortotropía del tensor y la razón "
            "por la que el artículo impone ρ ≥ 0.3."),
        "params": ("ρ = {rho} · β = {beta_pi}π · ν_s = {nu_s} · N = "
                   "{num_waves} · vista {resolucion_vista}³ · "
                   "homogeneización {resolucion_homogeneizacion}³ · semilla "
                   "{semilla} · contorno: {condiciones_contorno}"),
        "ayuda_rapido": (
            "Baja la malla de vista a 96³ y la de homogeneización a 24³. Las "
            "cotas se cumplen igual —una malla gruesa da rigideces más bajas, "
            "o sea más holgura—, pero las cifras que se reporten deberían "
            "salir de la ejecución completa, que homogeneiza a 40³ como pide "
            "el estudio de convergencia."),
    },
    {
        "clave": "guo2024",
        "pestana": "Guo 2024",
        "opcion": "--replicar-guo2024",
        "guion": "replicar_guo2024.py",
        "figura_ref": "Guo2024_Fig7.png",
        "pie_izq": "Original publicado (Fig. 7)",
        "json": "replica_guo2024.json",
        "figura": "replica_guo2024_fig7.png",
        "cita": (
            "Guo, Y., Sharma, S., &amp; Kumar, S. (2024). Inverse designing "
            "surface curvatures by deep learning. <i>Advanced Intelligent "
            "Systems, 6</i>(6), 2300789. "
            "https://doi.org/10.1002/aisy.202300789"),
        "licencia": (
            "Figura original bajo "
            "<a href='https://creativecommons.org/licenses/by/4.0/'>CC BY "
            "4.0</a> &mdash; &copy; 2024 The Authors, Wiley-VCH GmbH. "
            "Reproducida sin modificaciones."),
        "cabecera": (
            "<b>Aquí se replican las dos columnas «Target» de su Figura 7.</b> "
            "El spinodoide lleva parámetros publicados verbatim, y la "
            "superficie nodal periódica tiene curvaturas con fórmula cerrada: "
            "es la que valida nuestro estimador contra una respuesta exacta "
            "que no hemos calculado nosotros. No se replica su diseño inverso "
            "por redes neuronales, ni su muestra de hueso, que no tenemos."),
        "params": ("ρ = {rho} · β = {beta_pi}π · θ = {thetas} · N = "
                   "{num_waves} · spinodoide {resolucion_spinodoide}³ · PNS "
                   "{resolucion_pns}³ · semilla {semilla} · {unidad}"),
        "ayuda_rapido": (
            "Baja las rejillas a 96³ y 120³. Sirve para comprobar que la "
            "tubería funciona: el error del estimador discreto sube al 10 %, "
            "que es justo el criterio, y pasar por una centésima es pasar por "
            "suerte. Para dar por buena la implementación, ejecución "
            "completa."),
    },
]


def ruta_test():
    """La carpeta `Test`, este el codigo suelto o empaquetado.

    Empaquetada, PyInstaller deja los datos en `sys._MEIPASS`. Se prueban los
    dos sitios y se devuelve None si no esta: la validacion es un extra, y su
    ausencia no puede impedir que la aplicacion funcione.
    """
    cand = []
    emp = getattr(sys, "_MEIPASS", None)
    if emp:
        cand.append(Path(emp) / "Test")
    cand.append(Path(__file__).resolve().parent / "Test")
    for c in cand:
        try:
            if (c / "replicar_kumar2020.py").is_file():
                return c
        except Exception:
            pass
    return None


def replica_por_clave(clave):
    for r in REPLICAS:
        if r["clave"] == clave:
            return r
    return None


def orden_replica(rapido, clave="kumar2020"):
    """Programa y argumentos para lanzar una replica COMO PROCESO APARTE.

    POR QUE UN PROCESO Y NO UN HILO
      Las replicas terminan componiendo su figura con
      `pyvista.Plotter(off_screen=True)`, y VTK no puede crear su contexto
      OpenGL fuera del hilo principal. Lanzadas desde un QThread, la
      aplicacion MUERE de golpe -codigo 0xC0000409, sin excepcion que
      capturar- justo al llegar a la figura, despues de haber calculado
      correctamente todas las comprobaciones. Se reprodujo tal cual desde el
      ejecutable empaquetado y desde el codigo.

      En un proceso nuevo la figura se compone en su propio hilo principal y
      no hay conflicto. De regalo, un fallo en la replica ya no puede
      llevarse por delante la aplicacion, y la ventana sigue respondiendo.

    Empaquetada, `sys.executable` es el propio spinpy.exe y entiende la
    opcion; suelta, hay que pasarle ademas la ruta de visor.py.
    """
    rep = replica_por_clave(clave) or REPLICAS[0]
    args = [rep["opcion"]] + (["--rapido"] if rapido else [])
    if getattr(sys, "frozen", False):
        return sys.executable, args
    visor = Path(__file__).resolve().parent / "visor.py"
    return sys.executable, [str(visor)] + args


class _Imagen(QtWidgets.QLabel):
    """Imagen que se reescala con el hueco disponible, sin deformarse."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._orig = None
        self.setAlignment(QtCore.Qt.AlignCenter)
        self.setMinimumSize(220, 260)
        self.setStyleSheet("background:#fafafa; border:1px solid #ddd;")

    def poner(self, ruta):
        p = QtGui.QPixmap(str(ruta)) if ruta and Path(ruta).is_file() else None
        self._orig = p if (p and not p.isNull()) else None
        self._pintar()

    def _pintar(self):
        if self._orig is None:
            self.setText(_("(sin imagen)"))
            return
        self.setPixmap(self._orig.scaled(
            self.size(), QtCore.Qt.KeepAspectRatio,
            QtCore.Qt.SmoothTransformation))

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        self._pintar()


class PanelReplica(QtWidgets.QWidget):
    """Una replica: original a la izquierda, la nuestra a la derecha."""

    def __init__(self, rep, ruta, parent=None):
        super().__init__(parent)
        self.rep = rep
        self.ruta = ruta
        self.proc = None
        self._salida = []

        lay = QtWidgets.QVBoxLayout(self)

        cab = QtWidgets.QLabel(_(rep["cabecera"]))
        cab.setWordWrap(True)
        lay.addWidget(cab)

        cita = QtWidgets.QLabel(rep["cita"])
        cita.setWordWrap(True)
        cita.setOpenExternalLinks(True)
        cita.setStyleSheet("color:#1b4f72; font-size:11px;")
        lay.addWidget(cita)

        # --- las dos imagenes ---
        fila = QtWidgets.QHBoxLayout()
        for etq, attr in ((_(rep["pie_izq"]), "img_ref"),
                          (_("Regenerado con esta aplicación"), "img_rep")):
            caja = QtWidgets.QVBoxLayout()
            t = QtWidgets.QLabel("<b>%s</b>" % etq)
            t.setAlignment(QtCore.Qt.AlignCenter)
            caja.addWidget(t)
            im = _Imagen()
            setattr(self, attr, im)
            caja.addWidget(im, 1)
            fila.addLayout(caja, 1)
        lay.addLayout(fila, 1)

        pie = QtWidgets.QLabel(rep["licencia"])
        pie.setOpenExternalLinks(True)
        pie.setStyleSheet("color:#777; font-size:10px;")
        lay.addWidget(pie)

        # --- comprobaciones ---
        self.tabla = QtWidgets.QTableWidget(0, 4)
        self.tabla.setHorizontalHeaderLabels(
            [_("Código"), _("Predicción del artículo"), _("Obtenido"),
             _("Veredicto")])
        # Anchos FIJOS por modo, no `resizeColumnsToContents`: la columna
        # "Obtenido" lleva frases largas y se comia todo el ancho, dejando la
        # prediccion -que es lo que hay que leer- en tres palabras cortadas.
        c = self.tabla.horizontalHeader()
        c.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
        c.setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
        c.setSectionResizeMode(2, QtWidgets.QHeaderView.Stretch)
        c.setSectionResizeMode(3, QtWidgets.QHeaderView.ResizeToContents)
        self.tabla.setWordWrap(True)
        self.tabla.verticalHeader().setVisible(False)
        self.tabla.setMinimumHeight(190)
        self.tabla.setMaximumHeight(260)
        lay.addWidget(self.tabla)

        self.estado = QtWidgets.QLabel("")
        self.estado.setWordWrap(True)
        lay.addWidget(self.estado)

        self.barra = QtWidgets.QProgressBar()
        self.barra.setRange(0, 0)
        self.barra.hide()
        lay.addWidget(self.barra)

        # --- botones ---
        f = QtWidgets.QHBoxLayout()
        self.chk_rapido = QtWidgets.QCheckBox(_("Modo rápido (mallas menores)"))
        self.chk_rapido.setToolTip(_(rep["ayuda_rapido"]))
        f.addWidget(self.chk_rapido)
        f.addStretch(1)
        b = QtWidgets.QPushButton(_("Abrir la carpeta Test"))
        b.clicked.connect(self._abrir_carpeta)
        f.addWidget(b)
        self.btn = QtWidgets.QPushButton(_("Ejecutar la réplica"))
        self.btn.setStyleSheet("font-weight:bold; padding:6px;")
        self.btn.clicked.connect(self._correr)
        f.addWidget(self.btn)
        lay.addLayout(f)

        self._cargar_existente()

    # -- carga de lo que ya haya en disco ---------------------------------
    def _cargar_existente(self):
        if self.ruta is None:
            self.estado.setText(
                "<span style='color:#b62324'>" +
                _("No se encuentra la carpeta <b>Test</b>. La validación "
                  "necesita los guiones de réplica y "
                  "<code>Test/referencia/</code>.") + "</span>")
            self.btn.setEnabled(False)
            return
        self.img_ref.poner(self.ruta / "referencia" / self.rep["figura_ref"])
        self.img_rep.poner(self.ruta / "resultados" / self.rep["figura"])
        j = self.ruta / "resultados" / self.rep["json"]
        if j.is_file():
            try:
                self._mostrar(json.loads(j.read_text(encoding="utf-8")))
                return
            except Exception:
                pass
        self.estado.setText(_(
            "Todavía no se ha ejecutado esta réplica en este equipo. Pulsa "
            "«Ejecutar la réplica»."))

    def _mostrar(self, doc):
        ch = doc.get("comprobaciones", [])
        self.tabla.setRowCount(len(ch))
        for i, c in enumerate(ch):
            pasa = c.get("pasa")
            # La prediccion y el criterio son cadenas FIJAS del guion, asi
            # que se traducen. Lo "obtenido" no: son cifras y simbolos, que no
            # tienen idioma -y por eso el guion lo emite ya neutro-.
            for j, txt in enumerate((c.get("codigo", ""),
                                     _(c.get("prediccion", "")),
                                     c.get("obtenido", ""),
                                     _("cumple") if pasa else _("NO cumple"))):
                it = QtWidgets.QTableWidgetItem(str(txt))
                it.setFlags(QtCore.Qt.ItemIsEnabled)
                if j == 1 and c.get("criterio"):
                    it.setToolTip(_("Criterio: {c}").format(
                        c=_(c["criterio"])))
                if j == 3:
                    it.setForeground(QtGui.QBrush(QtGui.QColor(
                        "#1a7f37" if pasa else "#b62324")))
                self.tabla.setItem(i, j, it)
        self.tabla.resizeRowsToContents()

        p = dict(doc.get("parametros", {}))
        n = doc.get("fallos", 0)
        col = "#1a7f37" if n == 0 else "#b62324"
        try:
            det = _(self.rep["params"]).format(**p)
        except Exception:
            # Un informe viejo puede no llevar alguna clave. Se ensena lo que
            # haya en vez de romper la ventana entera por una linea de pie.
            det = " · ".join(f"{k} = {v}" for k, v in p.items())
        self.estado.setText(
            "<span style='color:%s'><b>%s</b></span><br>"
            "<span style='color:#666;font-size:10px'>%s</span>" % (
                col,
                (_("Las {n} comprobaciones se cumplen: la implementación "
                   "reproduce el artículo.").format(n=len(ch)) if n == 0
                 else _("{n} comprobación/es no se cumplen.").format(n=n)),
                det))

    # -- ejecucion ---------------------------------------------------------
    def _correr(self):
        if self.proc is not None and \
                self.proc.state() != QtCore.QProcess.NotRunning:
            return
        prog, args = orden_replica(self.chk_rapido.isChecked(),
                                   self.rep["clave"])
        self.btn.setEnabled(False)
        self.barra.show()
        self.estado.setText(_(
            "Calculando. Tarda unos minutos; la ventana responde igual."))
        self._salida = []
        self.proc = QtCore.QProcess(self)
        self.proc.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        self.proc.setWorkingDirectory(str(self.ruta))
        self.proc.readyReadStandardOutput.connect(self._leer)
        self.proc.finished.connect(self._termino)
        self.proc.errorOccurred.connect(self._no_arranco)
        self.proc.start(prog, args)

    def _leer(self):
        """Ultima linea util del proceso, como senal de vida."""
        txt = bytes(self.proc.readAllStandardOutput()).decode(
            "utf-8", "replace")
        self._salida.append(txt)
        util = [l.strip() for l in txt.splitlines() if l.strip()]
        if util:
            self.estado.setText(util[-1][:160])

    def _no_arranco(self, _err):
        self.barra.hide()
        self.btn.setEnabled(True)
        QtWidgets.QMessageBox.critical(
            self, _("La réplica falló"),
            _("No se pudo lanzar el proceso de réplica:\n{p}").format(
                p=self.proc.program()))

    def _termino(self, codigo, estado):
        self.barra.hide()
        self.btn.setEnabled(True)
        self.img_rep.poner(self.ruta / "resultados" / self.rep["figura"])
        j = self.ruta / "resultados" / self.rep["json"]
        if estado == QtCore.QProcess.NormalExit and j.is_file():
            try:
                self._mostrar(json.loads(j.read_text(encoding="utf-8")))
                return
            except Exception:
                pass
        # El proceso se cayo, o no dejo informe. Se ensena su salida: es lo
        # unico que permite entender que paso en la maquina de otro.
        QtWidgets.QMessageBox.critical(
            self, _("La réplica falló"),
            _("El proceso terminó con código {c}.").format(c=codigo)
            + "\n\n" + "".join(self._salida)[-4000:])

    def _abrir_carpeta(self):
        if self.ruta is None:
            return
        QtGui.QDesktopServices.openUrl(
            QtCore.QUrl.fromLocalFile(str(self.ruta)))


class DialogoValidacion(QtWidgets.QDialog):
    """Una pestana por articulo replicado."""

    def __init__(self, parent=None, clave=None):
        super().__init__(parent)
        self.setWindowTitle(_("Validación contra los resultados publicados"))
        self.resize(1280, 900)
        self.ruta = ruta_test()

        lay = QtWidgets.QVBoxLayout(self)
        self.tabs = QtWidgets.QTabWidget()
        lay.addWidget(self.tabs, 1)

        self.paneles = {}
        for rep in REPLICAS:
            # Una replica cuyo guion no esta -por ejemplo en un empaquetado
            # antiguo- no se ensena vacia: se omite y las demas siguen.
            if self.ruta is not None and \
                    not (self.ruta / rep["guion"]).is_file():
                continue
            pan = PanelReplica(rep, self.ruta, self)
            self.paneles[rep["clave"]] = pan
            self.tabs.addTab(pan, rep["pestana"])

        if not self.paneles:
            aviso = QtWidgets.QLabel(
                "<span style='color:#b62324'>" +
                _("No se encuentra la carpeta <b>Test</b>. La validación "
                  "necesita los guiones de réplica y "
                  "<code>Test/referencia/</code>.") + "</span>")
            aviso.setWordWrap(True)
            self.tabs.addTab(aviso, _("Validación"))

        if clave and clave in self.paneles:
            self.tabs.setCurrentWidget(self.paneles[clave])

        f = QtWidgets.QHBoxLayout()
        f.addStretch(1)
        c = QtWidgets.QPushButton(_("Cerrar"))
        c.clicked.connect(self.accept)
        f.addWidget(c)
        lay.addLayout(f)
