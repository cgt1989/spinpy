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

from spinpy import generar_mascara, generar_dual_lattice, leer_voi  # noqa: E402
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
        "lema1": "Microestructuras espinodales ajustadas",
        "lema2": "a hueso trabecular de micro-CT",
        "etapas": "morfometría  ·  ajuste  ·  homogeneización  ·  ensayo FE",
        "voi": "VOI micro-CT", "spin": "Spinodoide", "dual": "Dual-lattice",
        "isotropo": "Isótropo", "columnar": "Columnar", "lamelar": "Lamelar",
        "cubico": "Cúbico", "ajustado": "Ajustado al hueso",
        "menos": "Menos densidad",
        "cargado": "1 · VOI micro-CT cargado", "ajustar": "2 · Ajustar al VOI",
        "ajustando": "2 · Ajustando…",
        "medido": "3 · Candidato ajustado y medido",
        "ensayo": "4 · Ensayo de compresión en Z",
        "vm": "5 · Tensión de von Mises",
    },
    "en": {
        "lema1": "Spinodal microstructures fitted",
        "lema2": "to micro-CT trabecular bone",
        "etapas": "morphometry  ·  fitting  ·  homogenization  ·  FE testing",
        "voi": "micro-CT VOI", "spin": "Spinodoid", "dual": "Dual-lattice",
        "isotropo": "Isotropic", "columnar": "Columnar", "lamelar": "Lamellar",
        "cubico": "Cubic", "ajustado": "Fitted to bone",
        "menos": "Lower density",
        "cargado": "1 · micro-CT VOI loaded", "ajustar": "2 · Fit to the VOI",
        "ajustando": "2 · Fitting…",
        "medido": "3 · Fitted candidate, measured",
        "ensayo": "4 · Compression test along Z",
        "vm": "5 · von Mises stress",
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
    """Las tres estructuras del banner: el hueso y las dos familias."""
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
    n, sol = len(cubos), 0.06
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

    def foto(self, rotulo, resalta=None, ms=90, n=1):
        pix = self.app.primaryScreen().grabWindow(int(self.v.winId()))
        im = _a_pil(pix)
        d = ImageDraw.Draw(im, "RGBA")
        if resalta is not None:
            p = resalta.mapTo(self.v, self.QtCore.QPoint(0, 0))
            r = [p.x() - 5, p.y() - 5, p.x() + resalta.width() + 5,
                 p.y() + resalta.height() + 5]
            d.rounded_rectangle(r, radius=8, outline=ACENTO + (255,), width=4)
        if rotulo:
            f = _fuente(30, negrita=True)
            tw = d.textlength(rotulo, font=f)
            x, y = im.width - tw - 60, im.height - 110
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

    def cerrar(self):
        self.v.close()


def hacer_gif(voi, idioma="es"):
    g = Grabador(voi, idioma)
    try:
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
    ap.add_argument("--idioma", choices=sorted(TEXTOS), default="es")
    a = ap.parse_args()
    todo = not (a.banner or a.gif)
    if a.banner or todo:
        hacer_banner(a.voi, a.idioma)
    if a.gif or todo:
        hacer_gif(a.voi, a.idioma)


if __name__ == "__main__":
    main()
