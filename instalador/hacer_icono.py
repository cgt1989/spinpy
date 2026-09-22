# -*- coding: utf-8 -*-
"""Genera `spinpy.ico`: un spinodoide poroso en blanco y negro.

EL ICONO NO ES DECORACION
  Es como se encuentra la ventana en la barra de tareas, entre otras diez.
  Se renderiza con la MISMA tuberia que dibuja la aplicacion -mismo generador,
  mismo material, misma oclusion ambiental- de modo que lo que se ve en el
  icono es literalmente lo que el programa produce.

TRES DECISIONES, Y LAS TRES SE TOMARON MIRANDO EL RESULTADO A 16 PIXELES
  1. TESELA NEGRA. Se probaron cuatro variantes: estructura clara sobre fondo
     transparente, estructura oscura sobre transparente, y las dos sobre
     tesela. Sobre transparente NINGUNA funciona en los dos temas de Windows:
     la clara se lava sobre fondo claro y la oscura desaparece sobre fondo
     oscuro. La tesela negra con la estructura blanca se lee en ambos, que es
     el unico requisito que no se puede negociar.

  2. TRABECULAS GRUESAS. El numero de onda del icono es 7*pi y no el 13*pi
     que usa la aplicacion por defecto. Con 13*pi la estructura es fina y
     bonita a 256 pixeles, pero a 32 se convierte en ruido gris y a 16 en una
     mancha. Con 7*pi se distinguen los poros hasta el tamano mas pequeno.
     Es una eleccion de LEGIBILIDAD, no de fidelidad: el icono no tiene que
     representar una densidad concreta.

  3. ESCALA DE GRISES DE VERDAD, no un material gris. Se renderiza claro y se
     convierte a luminancia despues, para que el sombreado de la oclusion
     ambiental -que es lo que hace legible el relieve- sobreviva intacto.

Uso:  python hacer_icono.py
"""
from pathlib import Path
import sys

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))

from spinpy import generar_mascara            # noqa: E402
from visor import malla_de_mascara            # noqa: E402

LADO = 512          # se renderiza grande y se reduce: el remuestreo suaviza
RES = 72            # voxeles del cubo del icono
ONDA = 7            # x pi. Ver decision 2 en la cabecera.
DENSIDAD = 0.34
SEMILLA = 7
TESELA = (24, 24, 24, 255)     # negro, no negro puro: el puro se ve duro


def _render():
    BW, _, info = generar_mascara(resolution=RES, wave_number=ONDA * np.pi,
                                  num_waves=700, thetas=[15., 15., 45.],
                                  rho=DENSIDAD, seed=SEMILLA)
    malla = malla_de_mascara(BW, np.full(3, 1.0 / RES))

    p = pv.Plotter(off_screen=True, window_size=(LADO, LADO))
    p.enable_anti_aliasing("ssaa")
    p.add_mesh(malla, color="#ededed", smooth_shading=True, ambient=0.16,
               diffuse=0.86, specular=0.15, specular_power=20)
    b = np.asarray(malla.bounds, float)
    d = float(np.linalg.norm([b[1] - b[0], b[3] - b[2], b[5] - b[4]]))
    p.enable_ssao(radius=0.05 * d, bias=0.001 * d, kernel_size=192, blur=True)
    p.camera_position = "iso"
    p.camera.zoom(1.42)
    img = p.screenshot(return_img=True, transparent_background=True)
    p.close()
    return Image.fromarray(img).convert("RGBA"), info


def _recortar(im):
    """Recorta al contenido y cuadra. Sin esto el cubo sale diminuto."""
    alfa = np.asarray(im.getchannel("A"))
    ys, xs = np.where(alfa > 8)
    if not ys.size:
        return im
    m = 4
    y0, y1 = max(0, ys.min() - m), min(im.height, ys.max() + m)
    x0, x1 = max(0, xs.min() - m), min(im.width, xs.max() + m)
    lado = max(y1 - y0, x1 - x0)
    cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
    x0 = max(0, cx - lado // 2)
    y0 = max(0, cy - lado // 2)
    return im.crop((x0, y0, x0 + lado, y0 + lado))


def _a_grises(im):
    """Luminancia, conservando el canal alfa."""
    g = im.convert("L").convert("RGBA")
    g.putalpha(im.getchannel("A"))
    return g


def _tesela(contenido, n=256, radio=0.20, margen=0.10):
    base = Image.new("RGBA", (n, n), (0, 0, 0, 0))
    mascara = Image.new("L", (n, n), 0)
    ImageDraw.Draw(mascara).rounded_rectangle(
        [0, 0, n - 1, n - 1], radius=int(n * radio), fill=255)
    base.paste(Image.new("RGBA", (n, n), TESELA), (0, 0), mascara)
    c = int(n * (1 - 2 * margen))
    o = (n - c) // 2
    base.alpha_composite(contenido.resize((c, c), Image.LANCZOS), (o, o))
    return base


def main():
    im, info = _render()
    im = _a_grises(_recortar(im))
    icono = _tesela(im.resize((256, 256), Image.LANCZOS))

    destino = AQUI / "spinpy.ico"
    icono.save(destino, format="ICO",
               sizes=[(256, 256), (128, 128), (64, 64), (48, 48),
                      (32, 32), (16, 16)])
    icono.save(AQUI / "spinpy_icono.png")        # para el README y la web
    print(f"icono -> {destino}  ({destino.stat().st_size} bytes, "
          f"BV/TV {info['rho_obtenida']:.3f})")


if __name__ == "__main__":
    main()
