# -*- coding: utf-8 -*-
"""Genera el banner y los GIF animados del README de GitHub.

LO QUE SE VE ES LO QUE EL PROGRAMA HACE
  Igual que el icono (`instalador/hacer_icono.py`), nada de esto es un dibujo
  hecho aparte. El banner se renderiza con el mismo generador y la misma
  `malla_de_mascara` que la aplicacion, sobre el VOI equino real, y los GIF
  son capturas de la ventana del visor pulsando sus propios botones: si un
  dia la interfaz cambia, se vuelven a generar y siguen diciendo la verdad.

DOS DECISIONES
  1. CAPTURA DE PANTALLA, no `QWidget.grab()`. Los paneles 3D son OpenGL y
     `grab()` los devuelve en negro; `QScreen.grabWindow` lee los pixeles ya
     compuestos. Por eso la ventana se muestra de verdad mientras corre: no
     hay que tocar el raton ni taparla hasta que termine.
  2. GIF y no video. GitHub no reproduce un .mp4 del propio repositorio
     dentro del README; un GIF si. Se reduce a 960 px de ancho y paleta
     adaptativa, que deja cada animacion en unos pocos MB.

Uso (desde la raiz del repositorio):

    python docs/media/hacer_medios.py            # banner + GIF
    python docs/media/hacer_medios.py --banner   # solo las imagenes fijas
    python docs/media/hacer_medios.py --gif      # solo las animaciones
    python docs/media/hacer_medios.py --idioma en   # lo mismo en ingles

`--idioma en` escribe `banner_en.png`, `social_preview_en.png` y
`0X_*_en.gif` para `README.en.md`, con la interfaz del visor cambiada a
ingles por su propio `cambiar_idioma` (sin guardarlo en QSettings). El giro
de la cabecera no lleva texto y es el mismo para los dos README.

El GIF de ajuste necesita un VOI; por defecto busca
`../H4/Segmentadas/VOI_proximal_cubico.vtk` y se puede pasar otro con
`--voi ruta.vtk`.

El GIF de recorte (`00_recorte.gif`) parte de una CARPETA DE TIFF. Si no se
pasa `--pila carpeta`, la escribe el propio script con las rebanadas del VOI
porcino real de `datos/porcino_koria2020/VOI_C1.mat` (CC BY 4.0, 16 um/voxel),
que esta en el repositorio: el GIF se puede regenerar sin datos externos.
`--solo recorte` regenera solo esa animacion.
"""
from pathlib import Path
import argparse
import sys
import time

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[1]
sys.path.insert(0, str(RAIZ))

from spinpy import (generar_dual_lattice, generar_level_set,  # noqa: E402
                    generar_mascara, leer_voi)
from visor import malla_de_mascara                                  # noqa: E402

VOI_DEFECTO = RAIZ.parent / "H4" / "Segmentadas" / "VOI_proximal_cubico.vtk"

FONDO = (22, 24, 28)          # el negro de la tesela del icono, algo azulado
TINTA = (236, 236, 236)
SUAVE = (150, 156, 166)
ACENTO = (232, 163, 61)       # ambar: resalta sin competir con el hueso
HUESO = "#e8dcc4"             # marfil, el mismo tono de la vista del visor
ANCHO_GIF = 960

_FUENTES = Path("C:/Windows/Fonts")

# Todo texto que se dibuja sobre una imagen, por idioma.
TEXTOS = {
    "es": {
        "lema1": "Microestructuras sintéticas ajustadas",
        "lema2": "a hueso trabecular de micro-CT",
        "etapas": "morfometría  ·  ajuste  ·  homogeneización  ·  ensayo FE",
        "voi": "VOI micro-CT", "spin": "Spinodoide", "dual": "Dual-lattice",
        "ls": "Level-set foam",
        "isotropo": "Isótropo", "columnar": "Columnar", "lamelar": "Lamelar",
        "cubico": "Cúbico", "ajustado": "Ajustado al hueso",
        "menos": "Menos densidad",
        "cargado": "1 · VOI micro-CT cargado", "ajustar": "2 · Ajustar al VOI",
        "ajustando": "2 · Ajustando…",
        "medido": "3 · Candidato ajustado y medido",
        "ensayo": "4 · Ensayo de compresión en Z",
        "vm": "5 · Tensión de von Mises",
        "r_carpeta": "1 · Carpeta de rebanadas TIFF",
        "r_escala": "2 · Tamaño de voxel",
        "r_cargada": "3 · Pila cargada",
        "r_cubo": "4 · Cubo sobre el hueso",
        "r_mover": "5 · Se desplaza por el hueso",
        "r_lado": "6 · Tamaño por lado",
        "r_vol": "7 · Tamaño por volumen absoluto",
        "r_cil": "8 · Cilindro",
        "r_eje": "9 · Eje del cilindro: {eje}",
        "r_fin": "10 · VOI recortado",
    },
    "en": {
        "lema1": "Synthetic microstructures fitted",
        "lema2": "to micro-CT trabecular bone",
        "etapas": "morphometry  ·  fitting  ·  homogenization  ·  FE testing",
        "voi": "micro-CT VOI", "spin": "Spinodoid", "dual": "Dual-lattice",
        "ls": "Level-set foam",
        "isotropo": "Isotropic", "columnar": "Columnar", "lamelar": "Lamellar",
        "cubico": "Cubic", "ajustado": "Fitted to bone",
        "menos": "Lower density",
        "cargado": "1 · micro-CT VOI loaded", "ajustar": "2 · Fit to the VOI",
        "ajustando": "2 · Fitting…",
        "medido": "3 · Fitted candidate, measured",
        "ensayo": "4 · Compression test along Z",
        "vm": "5 · von Mises stress",
        "r_carpeta": "1 · Folder of TIFF slices",
        "r_escala": "2 · Voxel size",
        "r_cargada": "3 · Stack loaded",
        "r_cubo": "4 · Cube over the bone",
        "r_mover": "5 · Moved through the bone",
        "r_lado": "6 · Size by side",
        "r_vol": "7 · Size by absolute volume",
        "r_cil": "8 · Cylinder",
        "r_eje": "9 · Cylinder axis: {eje}",
        "r_fin": "10 · Cropped VOI",
    },
}


def _sufijo(idioma):
    """'' en espanol (los nombres de siempre), '_en' en ingles."""
    return "" if idioma == "es" else f"_{idioma}"


def _fuente(n, negrita=False):
    for nombre in (("segoeuib.ttf", "arialbd.ttf") if negrita
                   else ("segoeui.ttf", "arial.ttf")):
        f = _FUENTES / nombre
        if f.exists():
            return ImageFont.truetype(str(f), n)
    # Linux (donde se corre en un servidor virtual): DejaVu o Liberation.
    for ruta in (("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                  "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf")
                 if negrita else
                 ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
                  "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf")):
        if Path(ruta).exists():
            return ImageFont.truetype(ruta, n)
    return ImageFont.load_default()


# ---------------------------------------------------------------------------
# Render de un cubo, con la tuberia del icono
# ---------------------------------------------------------------------------

def _render_cubo(BW, spacing, color=HUESO, lado=900, azimut=0.0):
    malla = malla_de_mascara(BW, spacing)
    # Suavizado de Taubin SOLO para la imagen, como en las figuras del
    # informe: a tamano de banner la escalera de voxeles se ve como ruido.
    # Nunca se mide sobre una malla suavizada.
    malla = malla.smooth_taubin(n_iter=30, pass_band=0.1).compute_normals(
        auto_orient_normals=True, consistent_normals=False)
    p = pv.Plotter(off_screen=True, window_size=(lado, lado))
    p.enable_anti_aliasing("ssaa")
    p.add_mesh(malla, color=color, smooth_shading=True, ambient=0.18,
               diffuse=0.85, specular=0.18, specular_power=24)
    b = np.asarray(malla.bounds, float)
    d = float(np.linalg.norm([b[1] - b[0], b[3] - b[2], b[5] - b[4]]))
    p.enable_ssao(radius=0.05 * d, bias=0.001 * d, kernel_size=128, blur=True)
    p.camera_position = "iso"
    if azimut:
        p.camera.azimuth = azimut
    p.camera.zoom(0.95)
    img = p.screenshot(return_img=True, transparent_background=True)
    p.close()
    return _recortar(Image.fromarray(img).convert("RGBA"))


def _recortar(im):
    # La oclusion ambiental deja un velo casi transparente en todo el fondo;
    # sobre un lienzo oscuro se ve como un rectangulo. Se corta por debajo de
    # un umbral y se reescala el resto para conservar el borde suave.
    a = np.asarray(im.getchannel("A"), float)
    a = np.clip((a - 40) * 255 / 215, 0, 255).astype(np.uint8)
    im.putalpha(Image.fromarray(a))
    alfa = a
    ys, xs = np.where(alfa > 8)
    if not ys.size:
        return im
    return im.crop((xs.min(), ys.min(), xs.max() + 1, ys.max() + 1))


def _encajar(im, alto):
    return im.resize((max(1, round(im.width * alto / im.height)), alto),
                     Image.LANCZOS)


def _cubos(voi, T):
    """Las cuatro estructuras del banner: el hueso y las tres familias."""
    out = []
    if voi is not None and Path(voi).exists():
        VOI, sp = leer_voi(str(voi))
        out.append((T["voi"], _render_cubo(VOI, sp)))
    BW, _c, _i = generar_mascara(resolution=96, wave_number=15 * np.pi,
                                 num_waves=700, thetas=[15., 15., 60.],
                                 rho=0.30, seed=7)
    out.append((T["spin"], _render_cubo(BW, np.full(3, 1 / 96),
                                           color="#d9d9d9")))
    BW, _c, _i = generar_dual_lattice(96, celdas=5, rho=0.28,
                                      estiramiento=(1.0, 1.0, 1.6), seed=7)
    out.append((T["dual"], _render_cubo(BW, np.full(3, 1 / 96),
                                             color="#b9c7d6")))
    # Level-set foam (morphoGen): la forma que el ajuste encontro para los VOI
    # porcinos (celdas moderadas, un eje alargado, secciones convexas y nodos
    # algo engrosados).
    BW, _c, _i = generar_level_set(96, 134, rho=0.32,
                                   estiramiento=(1.0, 1.0, 1.5),
                                   curvatura=-0.5, c1=10.0, c2=0.5, seed=7)
    out.append((T["ls"], _render_cubo(BW, np.full(3, 1 / 96),
                                          color="#e6cf9a")))
    return out


def _banner(cubos, ancho, alto, titulo_px, destino, T, vertical=False):
    """Titulo y cubos. `vertical`: titulo arriba y cubos debajo (1280x640)."""
    g = np.linspace(0, 1, ancho)[None, :, None]
    base = np.zeros((alto, ancho, 4))
    base[..., 3] = 255
    base[..., :3] = np.array(FONDO) * (1 - 0.4 * g) + np.array([40, 46, 58]) * 0.4 * g
    lienzo = Image.fromarray(base.astype(np.uint8), "RGBA")
    d = ImageDraw.Draw(lienzo)

    ft = _fuente(titulo_px, negrita=True)
    fs = _fuente(int(titulo_px * 0.27))
    fp = _fuente(int(titulo_px * 0.19))
    x0 = int(titulo_px * 0.45)
    y = int(titulo_px * (0.25 if vertical else 0.55))
    d.text((x0, y), "spin", font=ft, fill=TINTA)
    w = d.textlength("spin", font=ft)
    d.text((x0 + w, y), "py", font=ft, fill=ACENTO)
    y += int(titulo_px * 1.34)
    d.text((x0, y), T["lema1"], font=fs, fill=TINTA)
    y += int(titulo_px * 0.36)
    d.text((x0, y), T["lema2"], font=fs, fill=TINTA)
    y += int(titulo_px * 0.48)
    d.text((x0, y), T["etapas"], font=fp, fill=SUAVE)
    y_fin_texto = y + int(titulo_px * 0.3)

    # Zona de los cubos: a la derecha del texto, o debajo.
    if vertical:
        zx0, zx1 = int(ancho * 0.04), int(ancho * 0.96)
        zy0, zy1 = y_fin_texto + 10, alto - int(alto * 0.07)
    else:
        zx0, zx1 = int(ancho * 0.50), ancho - int(alto * 0.06)
        zy0, zy1 = int(alto * 0.10), alto - int(alto * 0.16)
    n = len(cubos)
    sol = 0.06 if n <= 3 else -0.03      # con cuatro, un hueco en vez de solape
    hueco = (zx1 - zx0) / (n - (n - 1) * sol)
    fe = _fuente(max(14, int(titulo_px * 0.15)))
    ims = []
    for nombre, im in cubos:
        h = zy1 - zy0
        im = _encajar(im, h)
        if im.width > hueco:
            im = im.resize((int(hueco), round(im.height * hueco / im.width)),
                           Image.LANCZOS)
        ims.append((nombre, im))
    x = zx0
    for nombre, im in ims:
        cx = x + (hueco - im.width) / 2
        cy = zy0 + (zy1 - zy0 - im.height) / 2
        lienzo.alpha_composite(im, (int(cx), int(cy)))
        tw = d.textlength(nombre, font=fe)
        ImageDraw.Draw(lienzo).text((x + (hueco - tw) / 2, cy + im.height + 6), nombre,
                                    font=fe, fill=SUAVE)
        x += hueco * (1 - sol)
    lienzo.convert("RGB").save(destino, optimize=True)
    print("escrito", destino)


def hacer_banner(voi, idioma="es"):
    T, suf = TEXTOS[idioma], _sufijo(idioma)
    cubos = _cubos(voi, T)
    _banner(cubos, 1600, 440, 124, AQUI / f"banner{suf}.png", T)
    # 1280x640 es el tamano que pide GitHub para la vista previa social
    # (Settings -> Social preview). Se sube a mano: no vive en el repo.
    _banner(cubos, 1280, 640, 96, AQUI / f"social_preview{suf}.png", T,
            vertical=True)
    if idioma != "es":
        return              # el giro no lleva texto: sirve para los dos
    # Giro de 360 grados de la estructura ajustada, para la cabecera.
    BW, _c, _i = generar_mascara(resolution=72, wave_number=15 * np.pi,
                                 num_waves=700, thetas=[15., 15., 60.],
                                 rho=0.30, seed=7)
    cuadros = []
    for k in range(36):
        im = _render_cubo(BW, np.full(3, 1 / 72), color=HUESO, lado=520,
                          azimut=10.0 * k)
        c = Image.new("RGBA", (440, 440), FONDO + (255,))
        im = _encajar(im, 400)
        c.alpha_composite(im, ((440 - im.width) // 2, 20))
        cuadros.append(c.convert("RGB"))
    _guardar_gif(cuadros, AQUI / "giro.gif", ms=70, ancho=440)


# ---------------------------------------------------------------------------
# GIF de la aplicacion
# ---------------------------------------------------------------------------

def _guardar_gif(cuadros, destino, ms=90, ancho=ANCHO_GIF, duraciones=None):
    red = []
    for c in cuadros:
        if c.width != ancho:
            c = c.resize((ancho, round(c.height * ancho / c.width)),
                         Image.LANCZOS)
        red.append(c.convert("RGB"))
    # Una sola paleta para toda la animacion: con una por cuadro el fondo
    # parpadea. Se saca de un mosaico de varios cuadros repartidos.
    muestra = Image.new("RGB", (red[0].width, red[0].height * 4))
    for i, j in enumerate(np.linspace(0, len(red) - 1, 4).astype(int)):
        muestra.paste(red[j], (0, red[0].height * i))
    pal = muestra.quantize(colors=224, method=Image.MEDIANCUT)
    q = [c.quantize(palette=pal, dither=Image.Dither.NONE) for c in red]
    q[0].save(destino, save_all=True, append_images=q[1:], loop=0,
              duration=duraciones or ms, optimize=True, disposal=1)
    print("escrito", destino, f"{destino.stat().st_size / 1e6:.1f} MB",
          f"{len(q)} cuadros")


def _a_pil(pix):
    """QPixmap -> PIL. `Image.fromqpixmap` no reconoce PyQt5 en este entorno."""
    from PyQt5.QtGui import QImage
    q = pix.toImage().convertToFormat(QImage.Format_RGB888)
    w, h, b = q.width(), q.height(), q.bytesPerLine()
    ptr = q.bits()
    ptr.setsize(h * b)
    a = np.frombuffer(ptr, np.uint8).reshape(h, b)[:, :3 * w]
    return Image.fromarray(a.reshape(h, w, 3).copy(), "RGB")


class Grabador:
    """Conduce el visor y va sacando cuadros de la ventana."""

    def __init__(self, voi, idioma="es"):
        from PyQt5 import QtCore, QtWidgets
        import visor
        self.QtCore, self.QtWidgets = QtCore, QtWidgets
        pv.set_plot_theme("document")
        self.app = QtWidgets.QApplication.instance() or \
            QtWidgets.QApplication(sys.argv)
        v = visor.Visor()
        # El Visor se construye siempre en espanol (captura ahi sus textos
        # originales) y despues se traduce; sin guardar, para no cambiar el
        # idioma con que el usuario abre la aplicacion.
        v.cambiar_idioma(idioma, guardar=False)
        # Sin nadie delante: los dialogos de resultados no se abren y la
        # pregunta de orientacion se da por contestada (como en el informe
        # automatico, que hace exactamente esto).
        v._mostrar = lambda d: (None if callable(d) else d.deleteLater())
        v._confirmar_orientacion = lambda: True
        v.resize(1400, 820)
        v.move(20, 20)
        v.show()
        v.raise_()
        v.activateWindow()
        self.v, self.voi = v, voi
        self.T, self.suf = TEXTOS[idioma], _sufijo(idioma)
        self.scroll = v.findChild(QtWidgets.QScrollArea)
        self.cuadros, self.duraciones = [], []
        self.esperar(2.0)

    def esperar(self, minimo=0.2):
        t0 = time.time()
        while ((self.v.hilo is not None and self.v.hilo.isRunning())
               or time.time() - t0 < minimo):
            self.app.processEvents()
            time.sleep(0.02)
        for _ in range(5):
            self.app.processEvents()

    def ver(self, w):
        self.scroll.ensureWidgetVisible(w, 0, 60)
        self.esperar(0.1)

    def foto(self, rotulo, resalta=None, ms=90, n=1, region=False,
             pos="abajo-der"):
        """Un cuadro. `region=True` captura el RECTANGULO de la ventana
        principal desde la pantalla, de modo que un dialogo colocado encima
        sale en la imagen (con `grabWindow(winId)` solo saldria la ventana)."""
        if region:
            g = self.v.geometry()
            pix = self.app.primaryScreen().grabWindow(
                0, g.x(), g.y(), g.width(), g.height())
        else:
            pix = self.app.primaryScreen().grabWindow(int(self.v.winId()))
        im = _a_pil(pix)
        d = ImageDraw.Draw(im, "RGBA")
        if resalta is not None:
            if region:
                p = (resalta.mapToGlobal(self.QtCore.QPoint(0, 0))
                     - self.v.mapToGlobal(self.QtCore.QPoint(0, 0)))
            else:
                p = resalta.mapTo(self.v, self.QtCore.QPoint(0, 0))
            r = [p.x() - 5, p.y() - 5, p.x() + resalta.width() + 5,
                 p.y() + resalta.height() + 5]
            d.rounded_rectangle(r, radius=8, outline=ACENTO + (255,), width=4)
        if rotulo:
            f = _fuente(30, negrita=True)
            tw = d.textlength(rotulo, font=f)
            x, y = im.width - tw - 60, im.height - 110
            if pos == "arriba":
                x, y = (im.width - tw) / 2, 64
            d.rounded_rectangle([x - 22, y - 12, x + tw + 22, y + 50],
                                radius=12, fill=FONDO + (228,))
            d.text((x, y), rotulo, font=f, fill=TINTA)
        for _ in range(n):
            self.cuadros.append(im)
            self.duraciones.append(ms)

    def girar(self, rotulo, pasos=24, grados=7.5, ms=80):
        vs = [self.v.vis_voi] + [w for w in self.v._vis.values()
                                 if w.interactor.isVisible()]
        for _ in range(pasos):
            for w in vs:
                w.camera.Azimuth(grados)
                w.render()
            self.esperar(0.03)
            self.foto(rotulo, ms=ms)

    def guardar(self, nombre):
        base, ext = nombre.rsplit(".", 1)
        _guardar_gif(self.cuadros, AQUI / f"{base}{self.suf}.{ext}",
                     duraciones=self.duraciones)
        self.cuadros, self.duraciones = [], []

    # -- escenas --------------------------------------------------------------

    def explorar(self):
        v, sl = self.v, self.v.sl
        sl["resv"].fijar(56)
        self.esperar(0.5)
        self.ver(sl["dens"])
        T = self.T
        clases = [(T["isotropo"], 35, (90, 90, 90)),
                  (T["columnar"], 35, (15, 15, 0)),
                  (T["lamelar"], 40, (0, 0, 15)),
                  (T["cubico"], 35, (15, 15, 15)),
                  (T["ajustado"], 30, (15, 15, 60)),
                  (T["menos"], 22, (15, 15, 60))]
        for nombre, dens, (tx, ty, tz) in clases:
            v._aplicando = True
            sl["dens"].fijar(dens)
            sl["thx"].fijar(tx)
            sl["thy"].fijar(ty)
            sl["thz"].fijar(tz)
            v._aplicando = False
            v._generar_vista()
            self.esperar(0.3)
            self.foto(nombre, resalta=v.btn_gen, ms=700)
            self.girar(nombre, pasos=10, grados=6, ms=70)
        self.guardar("01_explorar.gif")

    def ajustar(self):
        v, T = self.v, self.T
        if not Path(self.voi).exists():
            print("sin VOI, se omite el GIF de ajuste:", self.voi)
            return False
        v.cargar_voi(str(self.voi))
        self.esperar(1.0)
        self.foto(T["cargado"], ms=1400)
        self.ver(v.btn_fit)
        self.foto(T["ajustar"], resalta=v.btn_fit, ms=1000)
        v.ajustar()
        t0 = time.time()
        while v.hilo is not None and v.hilo.isRunning():
            self.esperar(0.1)
            if time.time() - t0 > 4:
                self.foto(T["ajustando"], ms=260)
                t0 = time.time()
        self.esperar(1.0)
        self.ver(v.btn_med)
        v.medir()
        self.esperar(0.5)
        self.foto(T["medido"], resalta=v.btn_med, ms=1800)
        self.girar(T["medido"], pasos=36, grados=10)
        self.guardar("02_ajuste.gif")
        return True

    def mecanica(self):
        v, T = self.v, self.T
        v.spin_res_fe.setValue(40)
        v.cmb_eje_fe.setCurrentIndex(0)
        self.ver(v.btn_fe)
        self.foto(T["ensayo"], resalta=v.btn_fe, ms=1000)
        v.ensayo_fe()
        self.esperar(0.5)
        v.cmb_color.setCurrentIndex(3)       # tension de von Mises
        self.esperar(1.5)
        self.foto(T["vm"], resalta=v.cmb_color, ms=1600)
        self.girar(T["vm"], pasos=36, grados=10)
        self.guardar("03_mecanica.gif")

    # -- recorte: carpeta de TIFF -> cubo / cilindro ---------------------------

    def _centrar(self, d):
        """Pone un dialogo en el centro de la ventana principal."""
        g = self.v.geometry()
        d.move(g.x() + (g.width() - d.width()) // 2,
               g.y() + (g.height() - d.height()) // 2)

    def _exec_pila(self, d):
        """Sustituye a `DialogoPilaTiff.exec_`: rellena el dialogo como lo haria
        quien carga, y lo acepta."""
        T = self.T
        d.show()
        d.adjustSize()
        self._centrar(d)
        d.raise_()
        self.esperar(0.6)
        self.foto(T["r_carpeta"], resalta=d.ed_origen, ms=900, region=True,
                  pos="arriba")
        d.ed_origen.setText(str(self.pila))
        d._buscar_escala(str(self.pila))     # sin *_rec.log: avisa y pide la escala
        d.sp_tam.setValue(16.0)
        self.esperar(0.4)
        self.foto(T["r_escala"], resalta=d.sp_tam, ms=1800, region=True,
                  pos="arriba")
        d.accept()
        return self.QtWidgets.QDialog.Accepted

    def _exec_recorte(self, d):
        """Sustituye a `DialogoRecorte.exec_`: usa los mandos reales del dialogo."""
        T = self.T
        v = self.v
        d.show()
        d.raise_()
        g = v.geometry()
        d.resize(g.width(), g.height())      # tras show(): antes lo recoloca
        d.move(g.x(), g.y())
        self.esperar(3.0)                    # malla 3D, cubo y perfil
        print("dialogo", d.geometry(), "ventana", g)

        def f(rot, res=None, ms=700, esp=0.45):
            self.esperar(esp)
            self.foto(rot, resalta=res, ms=ms, region=True, pos="arriba")

        # 1. cubo por defecto
        f(T["r_cubo"], d.cmb_forma, ms=1600, esp=0.6)
        # 2. se desplaza por el hueso con el mando de posicion
        d.m_t.fijar(12.0)
        d._posicion_editada(0)
        f(T["r_mover"], d.m_t, ms=900, esp=0.7)
        for x in (24, 36, 48, 60, 72, 84):
            d.m_t.fijar(float(x))
            d._posicion_editada(0)
            f(T["r_mover"], d.m_t, ms=380, esp=0.5)
        d.m_o2.fijar(d.m_o2.sp.maximum() * 0.25)
        d._posicion_editada(0)
        f(T["r_mover"], d.m_o2, ms=700, esp=0.5)
        # de vuelta al centro del hueso: los tamanos se ensenan sobre hueso
        # y no sobre el borde, donde el BV/TV saldria subestimado
        d.m_t.fijar(50.0)
        d.m_o2.fijar(0.0)
        d._posicion_editada(0)
        f(T["r_mover"], d.m_t, ms=900, esp=0.7)
        # 3. tamano por lado
        d.sp_lado.setValue(1200.0)
        f(T["r_lado"], d.sp_lado, ms=1200, esp=0.6)
        # 4. tamano por volumen absoluto
        d.rb_vol.setChecked(True)
        d.ed_vol.setText("1.2e9")
        d._tamano_editado()
        f(T["r_vol"], d.ed_vol, ms=1600, esp=0.7)
        # 5. cilindro
        d.cmb_forma.setCurrentIndex(1)
        f(T["r_cil"], d.cmb_forma, ms=1600, esp=2.2)
        d.sp_radio.setValue(500.0)
        d.sp_alto.setValue(1400.0)
        f(T["r_cil"], d.sp_radio, ms=1000, esp=1.6)
        for k, eje in enumerate(("PC2", "PC3", "PC1")):
            d.cmb_eje_cil.setCurrentIndex((k + 1) % 3)
            f(T["r_eje"].format(eje=eje), d.cmb_eje_cil, ms=1300, esp=1.4)
        # 6. aceptar: el VOI recortado pasa a la ventana principal
        self.foto(T["r_fin"], resalta=d.findChild(
            self.QtWidgets.QDialogButtonBox), ms=1200, region=True,
            pos="arriba")
        d._aceptar()
        self._recorte_listo = True
        return self.QtWidgets.QDialog.Accepted

    def recortar(self, pila):
        import visor
        v, T = self.v, self.T
        # El dialogo de recorte necesita mas sitio que la ventana de los otros
        # GIF (su 3D, tres cortes y los mandos): se agranda la ventana.
        v.resize(1500, 1080)
        v.move(20, 20)
        self.esperar(0.5)
        self.pila = Path(pila)
        visor.DialogoPilaTiff.exec_ = lambda d: self._exec_pila(d)
        visor.DialogoRecorte.exec_ = lambda d: self._exec_recorte(d)
        self._recorte_listo = False
        v.cargar_pila_tiff()
        t0 = time.time()
        while v.VOI is None and time.time() - t0 < 60:
            self.esperar(0.2)
        self.esperar(1.0)
        self.foto(T["r_cargada"], resalta=v.b_recorte, ms=1800)
        v.recortar_voi()
        t0 = time.time()
        while not self._recorte_listo and time.time() - t0 < 180:
            self.app.processEvents()
            time.sleep(0.05)
        self.esperar(1.5)
        self.foto(T["r_fin"], ms=2200, pos="arriba")
        self.guardar("00_recorte.gif")

    def cerrar(self):
        self.v.close()


def _pila_de_ejemplo(destino):
    """Escribe como rebanadas TIFF el VOI porcino real del repositorio.

    Se guarda como (fila, columna) = (y, x), que es lo que escribe un escaner
    (el lector traspone al leer; ver tests/test_08_pila_tiff.py).
    """
    import tifffile
    VOI, _sp = leer_voi(str(RAIZ / "datos" / "porcino_koria2020" / "VOI_C1.mat"))
    destino.mkdir(parents=True, exist_ok=True)
    for k in range(VOI.shape[2]):
        tifffile.imwrite(destino / f"C1_{k:05d}.tif",
                         VOI[:, :, k].T.astype(np.uint8) * 255)
    return destino


def hacer_gif(voi, idioma="es", pila=None, solo=None):
    import tempfile
    g = Grabador(voi, idioma)
    try:
        if solo in (None, "recorte"):
            with tempfile.TemporaryDirectory() as tmp:
                g.recortar(pila or _pila_de_ejemplo(Path(tmp) / "pila"))
        if solo == "recorte":
            return
        g.explorar()
        if g.ajustar():
            g.mecanica()
    finally:
        g.cerrar()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--banner", action="store_true")
    ap.add_argument("--gif", action="store_true")
    ap.add_argument("--voi", default=str(VOI_DEFECTO))
    ap.add_argument("--pila", default=None,
                    help="carpeta de TIFF para el GIF de recorte")
    ap.add_argument("--solo", choices=["recorte"], default=None,
                    help="regenera solo esa animacion")
    ap.add_argument("--idioma", choices=sorted(TEXTOS), default="es")
    a = ap.parse_args()
    if a.solo:
        hacer_gif(a.voi, a.idioma, a.pila, a.solo)
        return
    todo = not (a.banner or a.gif)
    if a.banner or todo:
        hacer_banner(a.voi, a.idioma)
    if a.gif or todo:
        hacer_gif(a.voi, a.idioma, a.pila)


if __name__ == "__main__":
    main()
