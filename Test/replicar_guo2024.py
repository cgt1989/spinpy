# -*- coding: utf-8 -*-
"""
replicar_guo2024.py — Reproduce los objetivos de la Figura 7 de Guo et al.
                      (2024): el perfil de curvaturas de un spinodoide de
                      parametros publicados y el de una superficie nodal
                      periodica de ecuacion cerrada.

REFERENCIA REPLICADA (APA 7)
    Guo, Y., Sharma, S., & Kumar, S. (2024). Inverse designing surface
        curvatures by deep learning. Advanced Intelligent Systems, 6(6),
        2300789. https://doi.org/10.1002/aisy.202300789

    Licencia del articulo: Creative Commons Attribution (CC BY 4.0),
    https://creativecommons.org/licenses/by/4.0/
    (c) 2024 The Authors. Advanced Intelligent Systems published by Wiley-VCH.

    La figura original reproducida en `referencia/Guo2024_Fig7.png` se incluye
    al amparo de esa licencia, sin modificaciones y con atribucion. Ver
    `referencia/FUENTES.md`.


QUE SE REPLICA, Y QUE NO
------------------------
El articulo entrena dos redes neuronales para DISENAR una topologia con un
perfil de curvaturas dado. Eso no se replica: no tenemos su marco de campo de
fase ni sus 18.000 casos de entrenamiento, y no hace falta para lo que aqui
interesa. Se replican las dos columnas "Target" de su Figura 7, que son
geometrias completamente especificadas y medibles:

  (b) SPINODOIDE. beta = 15*pi, Q = 1000 ondas, rho = 0.3, thetas =
      (60, 30, 10) grados. VERBATIM del recuadro de la Fig. 7b. Es una
      especificacion MEJOR que la de la Fig. 2 de Kumar et al. (2020), donde
      las ternas de los paneles no se publicaron y hubo que deducir una.

  (c) SUPERFICIE NODAL PERIODICA (PNS). La ecuacion (10) del articulo,
          sin(x) sin(1.8 y) + sin(y) sin(1.8 z) + sin(z) sin(1.8 x) = 0.5
      VERBATIM. Al ser implicita y analitica, sus curvaturas tienen FORMULA
      CERRADA: es la unica de las tres geometrias del articulo contra la que se
      puede validar el estimador discreto sin discutir. El propio articulo la
      elige por eso —dice explicitamente que descarta las superficies minimas
      porque su distribucion (k2 = -k1) seria un objetivo trivial—, de modo que
      es ademas NO minima y sirve para comprobar el signo de la curvatura
      media.

  (a) HUESO TRABECULAR. NO se replica. Su muestra es de Tozzi et al. (2017) y
      no la tenemos; medir nuestros VOIs de sesamoideo y ponerlos al lado seria
      una comparacion entre especimenes distintos disfrazada de replica. El
      modulo `spinpy.curvatura` sirve igual para eso, pero es otro trabajo.


LA UNIDAD DE LA CURVATURA — ESTO ES DEDUCIDO, NO CITADO
--------------------------------------------------------
Los ejes de sus perfiles van de -100 a +100 y sus manchas caen hacia k1 ~ 30.
El articulo dice que el dominio es [0, 100]^3 y que "todas las dimensiones
estan normalizadas respecto del tamano del dominio", lo que admite dos
lecturas. La que cuadra con la fisica es que la curvatura esta en unidades de
1 / (lado del dominio):

  * spinodoide: la longitud de onda es 2*pi/beta = 2/15 del lado, y el radio de
    una trabecula es del orden de un cuarto de eso, 0.033 lados -> k ~ 30.
  * PNS: su celda de periodicidad es 10*pi (el minimo comun de 2*pi y 2*pi/1.8)
    y los radios son del orden de 1 en las unidades de la ecuacion, es decir
    1/31 del lado -> k ~ 30.
  * el limite del eje, k = 100, es exactamente un radio de 1/100 del lado, o
    sea UN VOXEL en una rejilla de 100^3: el borde del grafico cae justo en el
    limite de resolucion, que es donde tiene sentido ponerlo.

Las tres cosas cuadran a la vez, asi que aqui la curvatura se reporta en
1/(lado del cubo) y es directamente comparable con sus ejes. Si la lectura
fuera la otra, nuestras cifras estarian todas multiplicadas por 100 y la
comprobacion R1 fallaria a gritos, no en silencio.


UNA DIFERENCIA VISUAL QUE CONVIENE DECLARAR
-------------------------------------------
Nuestras manchas salen MAS ANCHAS que las suyas. No es un error de sitio —el
centro cae donde tiene que caer, que es lo que comprueba R1— sino de anchura, y
tiene dos causas que conviene separar:

  * ellos miden sobre un campo de fase con interfaz DIFUSA de espesor epsilon,
    que es una superficie intrinsecamente mas lisa que nuestro conjunto de
    nivel de un GRF muestreado en una rejilla;
  * el estimador discreto anade su propia dispersion.

La segunda se puede medir y se mide: el informe lleva la dispersion de H por
los DOS caminos, el discreto y el exacto, sobre los mismos vertices. Si las dos
salen parecidas, la anchura es de la geometria y no del estimador.


LO QUE SE COMPRUEBA
-------------------
  R1  El perfil del spinodoide cae donde cae el suyo: la mancha esta en el
      cuadrante de SILLA (k1 > 0 > k2) y su centro tiene k1 y k2 dentro de la
      caja leida de su panel (b). Los limites de esa caja se leyeron de la
      figura publicada ANTES de medir nada, y se anotan como lo que son:
      lectura a ojo de un mapa de densidad, no una cifra tabulada.
  R2  La curvatura escala con beta. Es exacto: el GRF con beta' = 2*beta es el
      mismo campo con las longitudes a la mitad, luego TODAS las curvaturas se
      duplican. Comprueba de una vez la normalizacion de longitudes, que es lo
      que mas facil se equivoca, y no depende de ninguna lectura de figura.
  R3  A rho = 0.5 la curvatura media es CERO. El campo y su opuesto tienen la
      misma distribucion, asi que al 50 % solido y vacio son intercambiables y
      H tiene que salir simetrico. Es una prediccion exacta del ensemble; una
      realizacion suelta la cumple con dispersion, y por eso el criterio se
      mide contra la anchura de la propia distribucion.
  R4  En el SPINODOIDE, el estimador discreto coincide con la formula cerrada.
      La esfera del bloque 09 de `tests/` ya lo comprueba, pero una esfera no
      tiene sillas ni cuellos: esta es la comprobacion sobre la geometria de
      verdad.
  R5  En la PNS, lo mismo contra su formula cerrada, que es independiente de
      la nuestra: la PNS no la genera `spinpy`, viene de la ecuacion (10).

      EL CRITERIO DEL 10 % ES PARA LA EJECUCION COMPLETA. Medido: a resolucion
      completa el error mediano sale del 5.3 % (spinodoide) y 4.8 % (PNS), con
      holgura; en `--rapido` sube al 10.0 % y 9.4 %, es decir JUSTO en el
      limite. Pasar por una centesima es pasar por suerte, asi que el modo
      rapido vale para comprobar que la tuberia funciona y no para dar por
      buena la implementacion.
  R6  La PNS NO es una superficie minima: su curvatura media no es cero.
  R7  El area de la isosuperficie del campo CONTINUO es menor que la que da
      marching cubes sobre la mascara BINARIA. Es la correccion C4 de este
      proyecto vista desde otro lado: el sesgo del binario se midio en +8.5 %
      sobre una esfera, y aqui tiene que reaparecer con el mismo signo.

Uso:
    python replicar_guo2024.py                # completo (~4 min)
    python replicar_guo2024.py --rapido       # rejillas menores (~1 min)
    python replicar_guo2024.py --solo-figura  # sin las comprobaciones caras
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))

from spinpy.curvatura import (curvaturas_de_campo,  # noqa: E402
                              mediana_ponderada, perfil, resumen)
from spinpy.grf import campo_grf, derivadas_grf, level_set  # noqa: E402
from spinpy.morphometry import area_superficie  # noqa: E402

CITA_APA = (
    "Guo, Y., Sharma, S., & Kumar, S. (2024). Inverse designing surface "
    "curvatures by deep learning. Advanced Intelligent Systems, 6(6), "
    "2300789. https://doi.org/10.1002/aisy.202300789"
)
DOI = "10.1002/aisy.202300789"
LICENCIA = "CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)"

# --- parametros del panel (b), VERBATIM del recuadro de la figura ----------
BETA_PI = 15.0
Q_ONDAS = 1000
RHO = 0.3
THETAS = [60.0, 30.0, 10.0]
SEMILLA = 20260720            # la del resto del proyecto; no sale del articulo

# --- caja leida a ojo del panel (b) de la figura publicada, ANTES de medir --
# La mancha de densidad de su objetivo esta centrada en el cuadrante de silla,
# alargada a lo largo de k1. Se anota el rectangulo con generosidad: lo que se
# comprueba es que caemos en su sitio, no una cifra que ellos no tabulan.
CAJA_B = {"k1": (10.0, 60.0), "k2": (-40.0, 10.0)}

# --- PNS, ecuacion (10) ----------------------------------------------------
PNS_NIVEL = 0.5
PNS_W = 1.8
# Celda de periodicidad: 2*pi y 2*pi/1.8 se repiten juntas cada 10*pi.
PNS_L = 10.0 * np.pi


# ---------------------------------------------------------------------------
# Campos y sus derivadas analiticas
# ---------------------------------------------------------------------------
def campo_pns(res):
    """Muestrea la PNS en una celda de periodicidad completa.

    Devuelve (campo, paso_normalizado, exacta) con el campo en las unidades de
    la ecuacion y el paso ya NORMALIZADO al lado del cubo, de modo que las
    curvaturas salen en 1/(lado) sin conversion posterior. `exacta` espera los
    vertices en esas mismas coordenadas normalizadas.
    """
    u = np.linspace(0.0, PNS_L, res)
    X, Y, Z = np.meshgrid(u, u, u, indexing="ij")
    F = (np.sin(X) * np.sin(PNS_W * Y)
         + np.sin(Y) * np.sin(PNS_W * Z)
         + np.sin(Z) * np.sin(PNS_W * X))

    def exacta(p):
        """grad y hess respecto de la coordenada NORMALIZADA (x/L).

        Por la regla de la cadena, derivar respecto de u = x/L multiplica el
        gradiente por L y el hessiano por L^2. Hacerlo aqui y no fuera evita
        que la conversion se pierda por el camino, que es como se cuelan los
        factores de escala.
        """
        q = np.asarray(p, float) * PNS_L
        x, y, z = q[:, 0], q[:, 1], q[:, 2]
        sx, sy, sz = np.sin(x), np.sin(y), np.sin(z)
        cx, cy, cz = np.cos(x), np.cos(y), np.cos(z)
        sX, sY, sZ = np.sin(PNS_W * x), np.sin(PNS_W * y), np.sin(PNS_W * z)
        cX, cY, cZ = np.cos(PNS_W * x), np.cos(PNS_W * y), np.cos(PNS_W * z)
        w = PNS_W

        g = np.stack([cx * sY + w * sz * cX,
                      w * sx * cY + cy * sZ,
                      w * sy * cZ + cz * sX], axis=1) * PNS_L

        H = np.zeros((len(q), 3, 3))
        H[:, 0, 0] = -sx * sY - w * w * sz * sX
        H[:, 1, 1] = -w * w * sx * sY - sy * sZ
        H[:, 2, 2] = -w * w * sy * sZ - sz * sX
        H[:, 0, 1] = H[:, 1, 0] = w * cx * cY
        H[:, 1, 2] = H[:, 2, 1] = w * cy * cZ
        H[:, 0, 2] = H[:, 2, 0] = w * cz * cX
        return g, H * PNS_L * PNS_L

    return F, 1.0 / (res - 1.0), exacta


def campo_spinodoide(res, beta_pi=BETA_PI, q=Q_ONDAS, thetas=None,
                     semilla=SEMILLA):
    """GRF del spinodoide y sus derivadas analiticas, en el cubo unidad.

    El paso de rejilla es 1/(res-1) y NO 1/res: `campo_grf` muestrea en
    linspace(0, 1, res), que incluye los dos extremos. Con 1/res la curvatura
    saldria un 1/res demasiado pequena — un 0.6 % a 160^3—, un error sin
    sintoma visible y del mismo tipo que el del tamano de voxel de las pilas.
    """
    thetas = THETAS if thetas is None else thetas
    beta = beta_pi * np.pi
    GRF, dirs, fases = campo_grf(res, beta, q, thetas, seed=semilla)

    def exacta(p):
        return derivadas_grf(p, dirs, fases, beta, q)

    return GRF, 1.0 / (res - 1.0), exacta


# ---------------------------------------------------------------------------
# Medida
# ---------------------------------------------------------------------------
def medir(campo, paso, nivel, exacta, margen=2, etiqueta=""):
    """Curvaturas y resumen de un conjunto de nivel, por los dos caminos."""
    t0 = time.time()
    d = curvaturas_de_campo(campo, paso, nivel, exacta=exacta, margen=margen)
    a = d["areas"]
    out = {
        "etiqueta": etiqueta,
        "n_vertices": d["n_vertices"],
        "n_descartados": d["n_descartados"],
        "discreta": resumen(d["k1"], d["k2"], a),
        "exacta": resumen(d["k1_exacta"], d["k2_exacta"], a),
        "segundos": round(time.time() - t0, 1),
    }
    dentro = a > 0
    out["error_discreta_vs_exacta"] = _error_relativo(
        d["k1"][dentro], d["k1_exacta"][dentro],
        d["k2"][dentro], d["k2_exacta"][dentro])
    return d, out


def _error_relativo(k1d, k1e, k2d, k2e):
    """Mediana de |k_discreta - k_exacta| en unidades de la escala tipica.

    Se normaliza por la mediana de |k_exacta| de las DOS curvaturas juntas, y
    no punto a punto: k2 pasa por cero en cualquier superficie con sillas, y un
    error relativo punto a punto explotaria ahi sin que pase nada raro.
    """
    esc = float(np.median(np.abs(np.r_[k1e, k2e])))
    esc = esc if esc > 0 else 1.0
    return {
        "escala": esc,
        "k1": float(np.median(np.abs(k1d - k1e)) / esc),
        "k2": float(np.median(np.abs(k2d - k2e)) / esc),
    }


def centro_del_perfil(k1, k2, areas):
    """Centro de la mancha de densidad: la mediana ponderada por area.

    Se usa la MEDIANA y no la media porque las colas del perfil llegan hasta el
    limite de resolucion —las esquinas de la malla tienen curvaturas enormes—
    y la media las sigue. Su figura ensena una mancha; el centro de la mancha
    es la mediana.
    """
    return (mediana_ponderada(np.asarray(k1, float),
                              np.asarray(areas, float)),
            mediana_ponderada(np.asarray(k2, float),
                              np.asarray(areas, float)))


# ---------------------------------------------------------------------------
# Comprobaciones
# ---------------------------------------------------------------------------
def comprobar(spin, pns, escala, simetria, areas):
    """Las siete predicciones, con sus criterios declarados de antemano."""
    C = []

    def anota(cod, texto, cond, obtenido, criterio):
        C.append({"codigo": cod, "prediccion": texto, "criterio": criterio,
                  "obtenido": obtenido, "pasa": bool(cond)})

    if spin:
        c1, c2 = spin["centro"]
        b1, b2 = CAJA_B["k1"], CAJA_B["k2"]
        dentro = (b1[0] <= c1 <= b1[1]) and (b2[0] <= c2 <= b2[1])
        anota("R1",
              "El perfil del spinodoide cae en la mancha de su panel (b): "
              "silla, con el centro dentro de la caja leida de la figura",
              dentro,
              f"centro (k1, k2) = ({c1:.1f}, {c2:.1f}) frente a la caja "
              f"[{b1[0]:.0f}, {b1[1]:.0f}] x [{b2[0]:.0f}, {b2[1]:.0f}]; "
              f"area de silla {100*spin['discreta']['silla']:.0f} %",
              # El criterio va como cadena FIJA y los limites en lo obtenido:
              # una cadena compuesta no se puede traducir, y este texto sale
              # en la ventana de Validacion tambien en ingles.
              "el centro cae dentro de la caja leida de su panel (b), "
              "en 1/(lado del cubo)")

        e = spin["error_discreta_vs_exacta"]
        peor = max(e["k1"], e["k2"])
        anota("R4",
              "En el spinodoide, el estimador discreto coincide con la "
              "formula cerrada del propio campo",
              peor <= 0.10,
              f"mediana del error: k1 {100*e['k1']:.1f} %, "
              f"k2 {100*e['k2']:.1f} % de la escala",
              "error mediano <= 10 % de la escala de curvatura")

    if escala:
        r, esperado = escala["razon"], 2.0
        anota("R2",
              "Doblar beta dobla todas las curvaturas: es el mismo campo con "
              "las longitudes a la mitad",
              abs(r - esperado) / esperado <= 0.05,
              f"k1(2 beta) / k1(beta) = {r:.3f}",
              "razon = 2.00 con 5 % de holgura")

    if simetria:
        H, s = simetria["H_medio"], simetria["H_desv"]
        anota("R3",
              "A rho = 0.5 solido y vacio son intercambiables: la curvatura "
              "media del ensemble es cero",
              abs(H) <= 0.15 * s,
              f"H medio = {H:+.2f} frente a una dispersion de {s:.1f}",
              "|H medio| <= 15 % de la desviacion tipica de H")

    if pns:
        e = pns["error_discreta_vs_exacta"]
        peor = max(e["k1"], e["k2"])
        anota("R5",
              "En la superficie nodal periodica —cuya formula NO es nuestra— "
              "el estimador discreto coincide con la respuesta cerrada",
              peor <= 0.10,
              f"mediana del error: k1 {100*e['k1']:.1f} %, "
              f"k2 {100*e['k2']:.1f} % de la escala",
              "error mediano <= 10 % de la escala de curvatura")

        H = pns["exacta"]["H_medio"]
        esc = e["escala"]
        anota("R6",
              "La PNS no es una superficie minima; el articulo la elige por "
              "eso, porque k2 = -k1 seria un objetivo trivial",
              abs(H) > 0.10 * esc,
              f"H medio = {H:+.1f}, escala de curvatura {esc:.1f}",
              "|H medio| > 10 % de la escala de curvatura")

    if areas:
        r = areas["razon"]
        anota("R7",
              "Marching cubes sobre la mascara BINARIA sobreestima el area "
              "frente a la isosuperficie del campo continuo (correccion C4)",
              1.02 <= r <= 1.20,
              f"BS_binaria / BS_campo = {r:.3f}",
              "entre 1.02 y 1.20; medido +8.5 % sobre una esfera")
    return C


# ---------------------------------------------------------------------------
# Figura
# ---------------------------------------------------------------------------
def figura(paneles, destino):
    """Dos filas: la estructura y su perfil, con los ejes de la Fig. 7."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pyvista as pv

    n = len(paneles)
    fig = plt.figure(figsize=(9.2, 4.0 * n))

    for i, p in enumerate(paneles):
        # --- estructura
        ax = fig.add_subplot(n, 2, 2 * i + 1)
        try:
            pl = pv.Plotter(off_screen=True, window_size=(700, 700))
            pl.set_background("white")
            try:
                pl.enable_anti_aliasing("ssaa")
            except Exception:
                pass
            malla = pv.PolyData(
                p["verts"],
                np.hstack([np.full((len(p["faces"]), 1), 3), p["faces"]]).ravel())
            malla = malla.compute_normals(auto_orient_normals=True,
                                          consistent_normals=False)
            pl.add_mesh(malla, color=p["color"], smooth_shading=True,
                        ambient=0.18, diffuse=0.82, specular=0.2,
                        specular_power=18)
            pl.camera_position = "iso"
            pl.camera.zoom(1.1)
            img = pl.screenshot(return_img=True)
            pl.close()
            ax.imshow(img)
        except Exception as e:                       # pragma: no cover
            ax.text(0.5, 0.5, f"sin render: {type(e).__name__}", ha="center")
        ax.axis("off")
        ax.set_title(p["titulo"], fontsize=11, loc="left")

        # --- perfil de curvaturas, con los ejes de su figura
        ax = fig.add_subplot(n, 2, 2 * i + 2)
        P, bx, by = perfil(p["k1"], p["k2"], p["areas"], limite=100.0,
                           nbins=200)
        # Se dibuja la transpuesta: histogram2d deja k1 en las FILAS, y en su
        # figura k1 es el eje horizontal.
        #
        # Se enmascara por debajo de la milesima del pico, y no solo el cero
        # exacto. Las colas del perfil llegan a las esquinas de la malla con
        # densidades de 1e-7, y con un mapa de color que empieza en negro eso
        # pinta de negro TODO el cuadrante y esconde la mancha, que es lo unico
        # que hay que mirar. Es una decision de dibujo y por eso se escribe: no
        # se tira ningun dato, el JSON los lleva todos.
        #
        # El mapa de color es `inferno` sobre fondo blanco, el mismo sentido
        # que el de su figura -oscuro para la densidad baja, amarillo en el
        # pico-, porque estos dos paneles se miran AL LADO de los suyos y dos
        # escalas invertidas entre si se comparan mal.
        M = np.ma.masked_where(P.T < 1e-3 * P.max(), P.T)
        im = ax.imshow(M, origin="lower",
                       extent=(bx[0], bx[-1], by[0], by[-1]), aspect="equal",
                       cmap="inferno")
        ax.set_facecolor("white")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).ax.tick_params(
            labelsize=7)
        ax.plot([-100, 100], [-100, 100], color="#888", lw=0.8)
        ax.axhline(0, color="#888", lw=0.6)
        ax.axvline(0, color="#888", lw=0.6)
        ax.set_xlim(-100, 100)
        ax.set_ylim(-100, 100)
        ax.set_xlabel("$\\kappa_1$")
        ax.set_ylabel("$\\kappa_2$")
        ax.set_title(p["titulo_perfil"], fontsize=10)

    fig.suptitle(
        "Replica de los objetivos de Guo et al. (2024), Fig. 7 — con spinpy\n"
        "curvaturas en 1/(lado del cubo), mismos ejes que el articulo",
        fontsize=12)
    fig.text(0.5, 0.014, "Replica de: " + CITA_APA.replace(". Advanced",
                                                           ".\nAdvanced"),
             ha="center", va="bottom", fontsize=7, color="#444")
    fig.text(0.5, 0.003, "Figura original bajo " + LICENCIA
             + " · (c) 2024 The Authors", ha="center", va="bottom",
             fontsize=7, color="#777")
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(destino, dpi=130)
    plt.close(fig)
    return destino


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--rapido", action="store_true",
                    help="rejillas menores; sirve para comprobar la tuberia")
    ap.add_argument("--solo-figura", action="store_true",
                    help="omite las comprobaciones que piden mallas extra")
    a = ap.parse_args()

    res = 96 if a.rapido else 160
    res_pns = 120 if a.rapido else 200
    q = 400 if a.rapido else Q_ONDAS
    extras = not a.solo_figura

    salida = AQUI / "resultados"
    salida.mkdir(exist_ok=True)

    print("=" * 78)
    print(" REPLICA DE Guo et al. (2024), objetivos de la Fig. 7")
    print(" " + CITA_APA)
    print(" licencia del articulo: " + LICENCIA)
    print("=" * 78)
    print(f"  spinodoide: beta={BETA_PI}pi  Q={q}  rho={RHO}  "
          f"thetas={THETAS}  rejilla {res}^3  semilla {SEMILLA}")
    print(f"  PNS: sin(x)sin(1.8y)+sin(y)sin(1.8z)+sin(z)sin(1.8x) = "
          f"{PNS_NIVEL}  ·  celda 10pi  ·  rejilla {res_pns}^3")
    print()

    t0 = time.time()
    paneles = []

    # --- (b) spinodoide -----------------------------------------------------
    print("  spinodoide: campo ...", end="", flush=True)
    GRF, paso, exacta = campo_spinodoide(res, q=q)
    phi0 = level_set(RHO)
    print(f" listo ({time.time()-t0:.0f} s), curvaturas ...", end="",
          flush=True)
    d_sp, r_sp = medir(GRF, paso, phi0, exacta, etiqueta="spinodoide")
    r_sp["centro"] = centro_del_perfil(d_sp["k1"], d_sp["k2"], d_sp["areas"])
    print(f" {r_sp['segundos']:.0f} s  ({r_sp['n_vertices']} vertices, "
          f"{r_sp['n_descartados']} de borde descartados)")
    print(f"    centro del perfil (k1, k2) = "
          f"({r_sp['centro'][0]:.1f}, {r_sp['centro'][1]:.1f})   "
          f"silla {100*r_sp['discreta']['silla']:.0f} %   "
          f"H = {r_sp['discreta']['H_medio']:+.1f}")
    paneles.append({
        "verts": d_sp["verts"], "faces": d_sp["faces"],
        "k1": d_sp["k1"], "k2": d_sp["k2"], "areas": d_sp["areas"],
        "color": "#7fa8d9",
        "titulo": (f"Spinodoide   $\\beta$ = {BETA_PI:g}$\\pi$ · "
                   f"$\\rho$ = {RHO} · $\\theta$ = "
                   f"({THETAS[0]:.0f}, {THETAS[1]:.0f}, {THETAS[2]:.0f})"
                   "$^\\circ$"),
        "titulo_perfil": "perfil de curvaturas (objetivo de su Fig. 7b)"})

    # --- (c) superficie nodal periodica ------------------------------------
    print("  PNS: campo ...", end="", flush=True)
    F, paso_p, exacta_p = campo_pns(res_pns)
    # El solido es el lado en el que la PNS forma tubos: {F >= nivel}. Como
    # `curvaturas_de_campo` toma siempre {campo <= nivel} como solido, se le
    # pasa el campo cambiado de signo. No es un apano: es la unica forma de
    # que la normal apunte al vacio, y de eso depende el SIGNO de todo.
    d_pns, r_pns = medir(-F, paso_p, -PNS_NIVEL,
                         lambda p: tuple(-x for x in exacta_p(p)),
                         etiqueta="PNS")
    r_pns["centro"] = centro_del_perfil(d_pns["k1"], d_pns["k2"],
                                        d_pns["areas"])
    print(f" {r_pns['segundos']:.0f} s  ({r_pns['n_vertices']} vertices)")
    print(f"    centro del perfil (k1, k2) = "
          f"({r_pns['centro'][0]:.1f}, {r_pns['centro'][1]:.1f})   "
          f"H = {r_pns['exacta']['H_medio']:+.1f}")
    paneles.append({
        "verts": d_pns["verts"], "faces": d_pns["faces"],
        "k1": d_pns["k1"], "k2": d_pns["k2"], "areas": d_pns["areas"],
        "color": "#8fbf8f",
        "titulo": ("Superficie nodal periodica\n"
                   "sin$x\\,$sin$1.8y$ + sin$y\\,$sin$1.8z$ + "
                   f"sin$z\\,$sin$1.8x$ = {PNS_NIVEL}"),
        "titulo_perfil": "perfil de curvaturas (objetivo de su Fig. 7c)"})

    # --- comprobaciones que piden mallas adicionales -----------------------
    escala = simetria = areas = None
    if extras:
        print("\n  R2  escalado con beta ...", end="", flush=True)
        n2 = 64 if a.rapido else 96
        med = []
        for factor in (1.0, 2.0):
            G2, p2, e2 = campo_spinodoide(int(n2 * factor),
                                          beta_pi=BETA_PI * factor,
                                          q=min(q, 400))
            d2 = curvaturas_de_campo(G2, p2, phi0, exacta=e2, margen=2)
            m = d2["areas"] > 0
            med.append(float(np.median(d2["k1_exacta"][m])))
        escala = {"k1_beta": med[0], "k1_2beta": med[1],
                  "razon": med[1] / med[0] if med[0] else np.nan}
        print(f" k1: {med[0]:.1f} -> {med[1]:.1f}   "
              f"razon {escala['razon']:.3f}")

        print("  R3  simetria a rho = 0.5 ...", end="", flush=True)
        n3 = 96 if a.rapido else 128
        G3, p3, e3 = campo_spinodoide(n3, q=min(q, 400))
        d3 = curvaturas_de_campo(G3, p3, level_set(0.5), exacta=e3, margen=2)
        m3 = d3["areas"] > 0
        H3 = (d3["k1_exacta"][m3] + d3["k2_exacta"][m3]) / 2.0
        w3 = d3["areas"][m3]
        Hm = float((w3 * H3).sum() / w3.sum())
        Hs = float(np.sqrt((w3 * (H3 - Hm) ** 2).sum() / w3.sum()))
        simetria = {"H_medio": Hm, "H_desv": Hs, "resolucion": n3}
        print(f" H = {Hm:+.2f} +- {Hs:.1f}")

        print("  R7  area binaria frente a area del campo ...", end="",
              flush=True)
        BW = GRF <= phi0
        BS_bin = area_superficie(BW, np.full(3, paso))
        BS_campo = float(np.sum(np.linalg.norm(np.cross(
            d_sp["verts"][d_sp["faces"][:, 1]] - d_sp["verts"][d_sp["faces"][:, 0]],
            d_sp["verts"][d_sp["faces"][:, 2]] - d_sp["verts"][d_sp["faces"][:, 0]]),
            axis=1)) / 2.0)
        areas = {"BS_binaria": float(BS_bin), "BS_campo": BS_campo,
                 "razon": float(BS_bin / BS_campo) if BS_campo else np.nan}
        print(f" {BS_bin:.3f} / {BS_campo:.3f} = {areas['razon']:.3f}")

    checks = comprobar(r_sp, r_pns, escala, simetria, areas)
    print("\n  PREDICCIONES")
    fallos = 0
    for c in checks:
        ok = "OK   " if c["pasa"] else "FALLA"
        fallos += (not c["pasa"])
        print(f"    [{ok}] {c['codigo']}  {c['prediccion']}")
        print(f"             {c['obtenido']}   (criterio: {c['criterio']})")

    fig = None
    print("\n  componiendo la figura...", end="", flush=True)
    try:
        fig = figura(paneles, salida / "replica_guo2024_fig7.png")
        print(f" -> {fig.name}")
    except Exception as e:
        print(f" fallo: {type(e).__name__}: {e}")

    doc = {
        "replica_de": {"cita_apa": CITA_APA, "doi": DOI, "licencia": LICENCIA,
                       "figura_original": "referencia/Guo2024_Fig7.png"},
        "generado": time.strftime("%Y-%m-%d %H:%M:%S"),
        "parametros": {"rho": RHO, "beta_pi": BETA_PI, "thetas": THETAS,
                       "num_waves": q, "semilla": SEMILLA,
                       "resolucion_spinodoide": res, "resolucion_pns": res_pns,
                       "nivel_pns": PNS_NIVEL, "celda_pns": "10 pi",
                       "unidad": "curvatura en 1/(lado del cubo)",
                       "caja_leida_panel_b": CAJA_B},
        "spinodoide": {k: v for k, v in r_sp.items()},
        "pns": {k: v for k, v in r_pns.items()},
        "escalado_beta": escala,
        "simetria_rho05": simetria,
        "areas": areas,
        "comprobaciones": checks,
        "fallos": fallos,
        "figura": fig.name if fig else None,
    }
    j = salida / "replica_guo2024.json"
    j.write_text(json.dumps(doc, indent=2, ensure_ascii=False,
                            default=float), encoding="utf-8")

    print(f"\n  total {time.time()-t0:.0f} s   ·   informe -> {j.name}")
    print("=" * 78)
    print(f"  COMPROBACIONES QUE FALLAN: {fallos}")
    print("=" * 78)
    return 1 if fallos else 0


if __name__ == "__main__":
    sys.exit(main())
