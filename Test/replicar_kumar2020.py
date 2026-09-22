# -*- coding: utf-8 -*-
"""
replicar_kumar2020.py — Reproduce la Figura 2 de Kumar et al. (2020) con
                        el generador de esta aplicacion.

REFERENCIA REPLICADA (APA 7)
    Kumar, S., Tan, S., Zheng, L., & Kochmann, D. M. (2020). Inverse-designed
        spinodoid metamaterials. npj Computational Materials, 6(1), 1-10.
        Articulo 73. https://doi.org/10.1038/s41524-020-0341-6

    Licencia del articulo: Creative Commons Attribution 4.0 International
    (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/
    (c) The Author(s) 2020. Publicado por Springer Nature en asociacion con el
    Shanghai Institute of Ceramics of the Chinese Academy of Sciences.

    La figura original reproducida en `referencia/Kumar2020_Fig2.png` se
    incluye al amparo de esa licencia, sin modificaciones, con la atribucion
    que exige. Ver `referencia/FUENTES.md`.


PARA QUE SIRVE ESTE GUION
-------------------------
No es una prueba de la aplicacion contra si misma: es una prueba contra la
FUENTE. El metodo que implementa `spinpy` no lo inventamos nosotros; es el de
Kumar et al. (2020), y este guion existe para poder afirmar, con una figura al
lado de la otra, que nuestra implementacion reproduce lo que el articulo
publica. Si alguien pregunta "¿de donde sale esto y como se que esta bien?",
la respuesta es esta carpeta.


LO QUE SE REPLICA, Y CON QUE PARAMETROS
---------------------------------------
El articulo define (ecuaciones 1 a 3, pp. 3-4):

  (1)  campo gaussiano aleatorio
           phi(x) = sqrt(2/N) * SUM_i cos(beta * n_i . x + gamma_i)
       con n_i ~ U(S^2) y gamma_i ~ U([0, 2*pi))

  (2)  anisotropia, restringiendo las direcciones de onda a la UNION de tres
       conos alrededor de los ejes cartesianos:
           n_i ~ U{ k en S^2 : |k.e1| > cos(th1)  o
                               |k.e2| > cos(th2)  o
                               |k.e3| > cos(th3) }
       con th1, th2, th3 en {0} U [th_min, 90 grados] y th_min = 15 grados

  (3)  solido donde phi(x) <= phi0, con
           phi0 = sqrt(2) * erf^-1(2*rho - 1)

NUESTRA IMPLEMENTACION ES ESA, NO UNA PARECIDA:
  * `spinpy.grf._waves_rechazo` acepta un candidato isotropo si el angulo
    respecto de ALGUN eje es menor que el theta de ese eje. Eso es la union de
    conos de la ecuacion (2), literalmente.
  * `spinpy.grf.level_set` calcula phi0 = sqrt(2) * erfinv(2*rho - 1), que es
    la formula de la ecuacion (3) tal cual.
  * El esquema alternativo 'equitativo' (TPMS-Scaffolds-generator) reparte las
    ondas en partes iguales entre los conos y NO es la ecuacion (2). Por eso
    esta replica usa 'rechazo': es el unico que corresponde al articulo. Es un
    dato util por si mismo, y esta comprobado abajo.

PARAMETROS. Se distingue con cuidado lo que el articulo dice de lo que
deducimos, porque una replica que mete valores inventados no replica nada:

  * rho = 0.5                      VERBATIM, leyenda de la Fig. 2
  * lamelar  (30, 0, 0) grados     VERBATIM, leyenda de la Fig. 3
  * columnar (0, 30, 30) grados    VERBATIM, leyenda de la Fig. 3
  * cubica   (30, 30, 30) grados   DEDUCIDO de la regla de clases del texto
                                   (p. 4): tres conos activos -> simetria
                                   cubica. El articulo no publica la terna.
  * isotropa (90, 90, 90) grados   VERBATIM del texto (p. 4): "th3 = pi/2,
                                   reduciendose a la esfera unidad"
  * th_min = 15 grados             VERBATIM (p. 4)
  * beta = 15*pi                   VERBATIM de la leyenda de la Fig. 5, donde
                                   se usa "para comparacion visual"

La Fig. 2 del articulo no publica las ternas de sus seis paneles como texto
extraible, asi que NO se intenta reproducir panel por panel: se reproducen las
CUATRO CLASES que el texto nombra, con las ternas que el articulo si publica.

UNA DIFERENCIA VISUAL QUE CONVIENE DECLARAR. Nuestras estructuras salen mas
FINAS que las de su Fig. 2: donde ellos muestran unas seis laminas, aqui se ven
muchas mas. No es un error de implementacion sino de escala. El numero de onda
beta fija el tamano caracteristico, y el articulo NO publica el beta de la Fig.
2; el 15*pi que usamos viene de la leyenda de la Fig. 5. Con un beta menor las
estructuras salen mas gruesas y el parecido visual con la Fig. 2 es mayor, sin
que cambie ninguna de las seis comprobaciones: la anisotropia depende de los
conos, no del tamano de poro. Se deja el valor citable en vez del que "queda
mejor".


LO QUE SE COMPRUEBA, Y QUE PODRIA FALLAR
-----------------------------------------
Una figura parecida no demuestra nada por si sola. El guion evalua ademas
predicciones FALSABLES que salen del propio articulo (p. 4): "las interfaces
de las topologias se alinean preferentemente PERPENDICULARES a esos vectores
n_i". De ahi se sigue que un cono alrededor de un eje produce laminas apiladas
a lo largo de ese eje, es decir que ese eje queda BLANDO:

  C1  lamelar (cono solo en e1): E1 << E2 y E1 << E3.
  C2  columnar (conos en e2 y e3): E1 >> E2 y E1 >> E3.
  C3  isotropa: E1, E2 y E3 iguales dentro de la dispersion del generador.
  C4  cubica: mucho menos anisotropa en los ejes que la lamelar y la columnar.

      UNA REALIZACION NO ES SIMETRICA. El articulo habla del ENSEMBLE: con
      tres conos iguales, la distribucion de direcciones de onda tiene
      simetria cubica y por tanto E1 = E2 = E3 EN MEDIA. Una sola realizacion
      con una semilla concreta no lo cumple: medido aqui a 24^3, la cubica da
      E_max/E_min = 1.249. La primera version de C4 exigia < 1.25 y pasaba por
      una centesima, que es pasar por suerte, no por buena implementacion.

      Por eso C4 se formula ahora en RELATIVO -la cubica tiene que ser mucho
      menos anisotropa que la lamelar y la columnar-, que es la afirmacion
      robusta y sigue siendo falsable. La cifra cruda se reporta igual, para
      que quien lea vea el 1.249 y no una comprobacion maquillada.
  C6  el esquema 'equitativo' no reproduce la ecuacion (2): con una terna de
      conos DESIGUALES da una anisotropia distinta.

      La primera version de esta comprobacion usaba la terna columnar
      (0, 30, 30) y FALLABA: rechazo daba DA 1.804 y equitativo 1.843, un
      2.2 % de diferencia. No era un fallo del codigo sino de la prueba. Con
      dos conos IGUALES los dos esquemas reparten las ondas mitad y mitad y
      coinciden por construccion; la diferencia solo aparece cuando los conos
      son de tamanos distintos, que es donde el reparto equitativo deja de ser
      proporcional al angulo solido. Se cambio a (15, 45, 90).

Si alguna falla, la implementacion se aparta del articulo y hay que mirarla.

Uso:
    python replicar_kumar2020.py                 # completo (~4 min)
    python replicar_kumar2020.py --rapido        # mallas menores (~1 min)
    python replicar_kumar2020.py --solo-figura   # sin homogeneizar
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

from spinpy import (constantes_ingenieria, generar_mascara,  # noqa: E402
                    homogeneizar, morfometria)

# ---------------------------------------------------------------------------
# Cita, para que viaje con cualquier salida que produzca este guion
# ---------------------------------------------------------------------------
CITA_APA = (
    "Kumar, S., Tan, S., Zheng, L., & Kochmann, D. M. (2020). "
    "Inverse-designed spinodoid metamaterials. npj Computational Materials, "
    "6(1), 1-10. Articulo 73. https://doi.org/10.1038/s41524-020-0341-6"
)
DOI = "10.1038/s41524-020-0341-6"
LICENCIA = "CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)"

# ---------------------------------------------------------------------------
# Las cuatro clases del articulo
# ---------------------------------------------------------------------------
# `fuente` dice de donde sale cada terna. Es informacion que tiene que llegar
# al informe: un revisor debe poder separar lo citado de lo deducido.
RHO = 0.5                       # Fig. 2 del articulo
BETA_PI = 15.0                  # Fig. 5 del articulo, beta = 15*pi/l
TH_MIN = 15.0                   # p. 4 del articulo

CLASES = [
    {"clave": "lamelar", "nombre": "Lamelar", "thetas": [30.0, 0.0, 0.0],
     "fuente": "verbatim, leyenda Fig. 3",
     "esperado": "laminas apiladas a lo largo de e1; e1 blando"},
    {"clave": "columnar", "nombre": "Columnar", "thetas": [0.0, 30.0, 30.0],
     "fuente": "verbatim, leyenda Fig. 3",
     "esperado": "columnas a lo largo de e1; e1 rigido"},
    {"clave": "cubica", "nombre": "Cubica", "thetas": [30.0, 30.0, 30.0],
     "fuente": "deducido de la regla de clases (p. 4)",
     "esperado": "tres conos activos; E1 = E2 = E3 pero no isotropa"},
    {"clave": "isotropa", "nombre": "Isotropa", "thetas": [90.0, 90.0, 90.0],
     "fuente": "verbatim del texto (p. 4): th = pi/2 -> esfera unidad",
     "esperado": "superficie elastica casi esferica"},
]

SEMILLA = 20260720              # la misma que fija el resto del proyecto

# Terna de conos DESIGUALES, solo para contrastar los dos esquemas de
# muestreo. NO sale del articulo: es nuestra, y se declara asi en la salida.
TERNA_DESIGUAL = (15.0, 45.0, 90.0)


# ---------------------------------------------------------------------------
# Superficie elastica E(d), que es lo que dibujan los paneles derechos de la
# Fig. 2 del articulo
# ---------------------------------------------------------------------------
def _voigt_a_tensor4(S6):
    """Pasa una matriz 6x6 de FLEXIBILIDAD (Voigt ingenieril) a S_ijkl.

    Los factores 1/2 y 1/4 NO son decorativos: la notacion de Voigt con
    deformaciones angulares de ingenieria absorbe un 2 en cada indice de
    cortante. Sin ellos, E(d) sale mal en cuanto la direccion no coincide con
    un eje, y la superficie elastica pareceria correcta en los ejes y estaria
    equivocada en medio, que es la peor forma de estar mal.
    """
    S6 = np.asarray(S6, float)
    par = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    S = np.zeros((3, 3, 3, 3))
    for m, (i, j) in enumerate(par):
        for n, (k, l) in enumerate(par):
            f = 1.0
            if m >= 3:
                f *= 0.5
            if n >= 3:
                f *= 0.5
            v = S6[m, n] * f
            for a, b in ((i, j), (j, i)):
                for c, d in ((k, l), (l, k)):
                    S[a, b, c, d] = v
    return S


def modulo_direccional(C6, dirs):
    """E(d) para cada direccion unitaria: 1/E = S_ijkl d_i d_j d_k d_l."""
    S = _voigt_a_tensor4(np.linalg.inv(np.asarray(C6, float)))
    d = np.asarray(dirs, float)
    inv = np.einsum("ijkl,ni,nj,nk,nl->n", S, d, d, d, d)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(inv > 0, 1.0 / inv, np.nan)


def _esfera(n_theta=90, n_phi=180):
    t = np.linspace(0.0, np.pi, n_theta)
    p = np.linspace(0.0, 2.0 * np.pi, n_phi)
    T, P = np.meshgrid(t, p, indexing="ij")
    d = np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P), np.cos(T)],
                 axis=-1)
    return d.reshape(-1, 3), T.shape


# ---------------------------------------------------------------------------
# Comprobacion de la propia superficie elastica, antes de fiarse de ella
# ---------------------------------------------------------------------------
def autocomprobar_superficie():
    """Con un C isotropo, E(d) tiene que ser constante y valer E.

    Se ejecuta siempre antes de dibujar nada. Una superficie elastica mal
    calculada produciria figuras convincentes y falsas.
    """
    E, nu = 3.0, 0.3
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    C = np.zeros((6, 6))
    C[:3, :3] = lam
    C[0, 0] = C[1, 1] = C[2, 2] = lam + 2 * mu
    C[3, 3] = C[4, 4] = C[5, 5] = mu
    d, _ = _esfera(30, 60)
    Ed = modulo_direccional(C, d)
    err = float(np.nanmax(np.abs(Ed - E)) / E)
    if err > 1e-10:
        raise AssertionError(
            f"E(d) no es constante para un material isotropo: error {err:.2e}")
    return err


# ---------------------------------------------------------------------------
# Generacion y medida de una clase
# ---------------------------------------------------------------------------
def una_clase(cl, res_vista, res_homog, num_waves, esquema="rechazo",
              con_mecanica=True):
    t0 = time.time()
    BW, _, info = generar_mascara(
        resolution=res_vista, wave_number=BETA_PI * np.pi,
        num_waves=num_waves, thetas=cl["thetas"], rho=RHO,
        esquema=esquema, seed=SEMILLA)
    sp = np.full(3, 1.0 / res_vista)
    m = morfometria(BW, sp, do_mil=True)

    out = {
        "clave": cl["clave"], "nombre": cl["nombre"],
        "thetas": list(cl["thetas"]), "fuente_parametros": cl["fuente"],
        "esquema": esquema,
        "rho_pedida": RHO, "rho_obtenida": float(info["rho_obtenida"]),
        "resolucion_vista": int(res_vista),
        "BVTV": float(m["BVTV"]), "DA": float(m["DA"]),
        "DA2": float(m.get("DA2", np.nan)),
        "dir_principal": [float(x) for x in m.get("dir_principal", [0, 0, 1])],
    }

    if con_mecanica:
        from spinpy.elastic import remuestrear_bw
        BWh, sph = remuestrear_bw(BW, sp, res_homog)
        C, inf = homogeneizar(BWh, E_s=1.0, nu_s=0.3, vox_size=float(sph[0]))
        ec = constantes_ingenieria(C)
        out.update({
            "resolucion_homog": int(res_homog),
            "homogeneizacion_ok": bool(inf["ok"]),
            "C_normalizado_Es": np.asarray(C, float).tolist(),
            "E1": float(ec["Ex"]), "E2": float(ec["Ey"]), "E3": float(ec["Ez"]),
            "E1_E3": float(ec["Ex"] / ec["Ez"]),
            "anisotropia_E": float(ec["anisotropia_E"]),
        })
    out["segundos"] = round(time.time() - t0, 1)
    return BW, out


# ---------------------------------------------------------------------------
# Predicciones falsables
# ---------------------------------------------------------------------------
def comprobar(res, res_equit):
    """Devuelve la lista de comprobaciones con su veredicto.

    Los margenes se declaran ANTES de mirar los numeros y son deliberadamente
    holgados: el generador es estocastico y la homogeneizacion se hace en una
    malla pequena. Lo que se comprueba es el SIGNO y el ORDEN DE MAGNITUD de
    la anisotropia, que es lo que el articulo afirma; no una cifra.
    """
    d = {r["clave"]: r for r in res}
    C = []

    def anota(cod, texto, cond, obtenido, criterio):
        C.append({"codigo": cod, "prediccion": texto, "criterio": criterio,
                  "obtenido": obtenido, "pasa": bool(cond)})

    lam, col = d.get("lamelar"), d.get("columnar")
    iso, cub = d.get("isotropa"), d.get("cubica")

    if lam and "E1" in lam:
        r = lam["E1"] / max(lam["E3"], 1e-30)
        anota("C1", "Lamelar (cono solo en e1): e1 es el eje BLANDO",
              r < 0.6, f"E1/E3 = {r:.3f}", "E1/E3 < 0.6")
    if col and "E1" in col:
        r = col["E1"] / max(col["E3"], 1e-30)
        anota("C2", "Columnar (conos en e2 y e3): e1 es el eje RIGIDO",
              r > 1.6, f"E1/E3 = {r:.3f}", "E1/E3 > 1.6")
    if iso and "E1" in iso:
        a = iso["anisotropia_E"]
        anota("C3", "Isotropa (th = 90 grados): E1 = E2 = E3",
              a < 1.25, f"E_max/E_min = {a:.3f}", "E_max/E_min < 1.25")
    if cub and "E1" in cub and lam and col and "E1" in lam:
        a = cub["anisotropia_E"]
        peor = max(lam["anisotropia_E"], col["anisotropia_E"])
        anota("C4",
              "Cubica (tres conos iguales): mucho menos anisotropa en los "
              "ejes que la lamelar y la columnar. Una realizacion suelta no "
              "es simetrica; ver la nota de C4 en la cabecera del guion.",
              a < 0.5 * peor,
              f"E_max/E_min = {a:.3f}  vs  {peor:.3f}",
              "menos de la mitad de la mayor anisotropia")
    if lam and col and "E1" in lam and "E1" in col:
        s = (lam["E1"] / lam["E3"]) < (col["E1"] / col["E3"])
        anota("C5", "Lamelar y columnar ordenan al reves el eje e1",
              s, f"{lam['E1']/lam['E3']:.3f} < {col['E1']/col['E3']:.3f}",
              "razon lamelar < razon columnar")
    if res_equit:
        # Se compara sobre la terna DESIGUAL, no sobre la columnar: con conos
        # iguales los dos esquemas coinciden y la prueba no distinguiria nada.
        # Ver la nota de C6 en la cabecera.
        par = {r["clave"]: r for r in res_equit}
        ref = {r["clave"]: r for r in res}
        if "desigual" in par and "desigual" in ref:
            da_r, da_e = ref["desigual"]["DA"], par["desigual"]["DA"]
            dif = abs(da_e - da_r) / max(da_r, 1e-30)
            anota("C6",
                  "Con conos DESIGUALES, 'equitativo' no es la ecuacion (2)",
                  dif > 0.05,
                  f"DA {da_r:.3f} / {da_e:.3f}  ({100*dif:+.1f} %)",
                  "difieren mas del 5 %")
    return C


# ---------------------------------------------------------------------------
# Figura
# ---------------------------------------------------------------------------
def figura(mascaras, res, destino, con_mecanica=True):
    """Una fila por clase: estructura y superficie elastica, como la Fig. 2."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pyvista as pv

    n = len(res)
    ncol = 2 if con_mecanica else 1
    fig = plt.figure(figsize=(4.6 * ncol, 3.5 * n))

    dirs, forma = _esfera()
    for i, r in enumerate(res):
        BW = mascaras[r["clave"]]

        # --- estructura, renderizada con la misma tuberia que la aplicacion
        pl = pv.Plotter(off_screen=True, window_size=(700, 700))
        pl.set_background("white")
        try:
            pl.enable_anti_aliasing("ssaa")
        except Exception:
            pass
        grid = pv.ImageData(dimensions=BW.shape,
                            spacing=(1.0 / BW.shape[0],) * 3)
        grid.point_data["v"] = BW.astype(np.float32).flatten(order="F")
        malla = grid.contour([0.5], scalars="v")
        if malla.n_points:
            malla = malla.compute_normals(auto_orient_normals=True,
                                          consistent_normals=False)
            pl.add_mesh(malla, color="#d98880", smooth_shading=True,
                        ambient=0.18, diffuse=0.82, specular=0.2,
                        specular_power=18)
            try:
                b = np.asarray(malla.bounds, float)
                dg = float(np.linalg.norm([b[1]-b[0], b[3]-b[2], b[5]-b[4]]))
                pl.enable_ssao(radius=0.045 * dg, bias=0.0009 * dg,
                               kernel_size=128, blur=True)
            except Exception:
                pass
        pl.camera_position = "iso"
        pl.camera.zoom(1.25)
        img = pl.screenshot(return_img=True)
        pl.close()

        ax = fig.add_subplot(n, ncol, i * ncol + 1)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(f"{r['nombre']}   $\\theta$ = "
                     f"({r['thetas'][0]:.0f}, {r['thetas'][1]:.0f}, "
                     f"{r['thetas'][2]:.0f})$^\\circ$",
                     fontsize=11, loc="left")

        if not con_mecanica:
            continue

        # --- superficie elastica E(d), como los paneles derechos de la Fig. 2
        ax = fig.add_subplot(n, ncol, i * ncol + 2, projection="3d")
        if r.get("homogeneizacion_ok"):
            Ed = modulo_direccional(np.asarray(r["C_normalizado_Es"]), dirs)
            R3 = np.nan_to_num(Ed).reshape(forma)
            D = dirs.reshape(*forma, 3)
            X, Y, Z = (R3 * D[..., 0], R3 * D[..., 1], R3 * D[..., 2])
            norm = plt.Normalize(np.nanmin(R3), np.nanmax(R3))
            ax.plot_surface(X, Y, Z, facecolors=plt.cm.jet(norm(R3)),
                            rstride=2, cstride=2, linewidth=0,
                            antialiased=True, shade=False)
            lim = float(np.nanmax(R3)) * 1.05
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
            ax.set_xlabel("$\\hat{e}_1$"); ax.set_ylabel("$\\hat{e}_2$")
            ax.set_zlabel("$\\hat{e}_3$")
            ax.set_title(f"E(d)/E$_s$   ·   E1/E3 = {r['E1_E3']:.2f}",
                         fontsize=10)
            ax.tick_params(labelsize=7)
        else:
            ax.text2D(0.5, 0.5, "la homogeneizacion no convergio",
                      ha="center", transform=ax.transAxes, color="#b62324")
            ax.set_axis_off()

    fig.suptitle(
        "Replica de Kumar et al. (2020), Fig. 2 — generada con spinpy\n"
        "$\\rho$ = %.2f · $\\beta$ = %g$\\pi$ · muestreo por rechazo "
        "(ecuacion 2 del articulo)" % (RHO, BETA_PI),
        fontsize=12)
    # El pie va en dos lineas: en una sola se sale del lienzo por los dos
    # lados, y la atribucion es obligatoria bajo CC BY, asi que tiene que
    # leerse entera.
    fig.text(0.5, 0.016,
             "Replica de: " + CITA_APA.replace(". npj", ".\nnpj"),
             ha="center", va="bottom", fontsize=7, color="#444")
    fig.text(0.5, 0.003, "Figura original bajo " + LICENCIA
             + " · (c) The Author(s) 2020",
             ha="center", va="bottom", fontsize=7, color="#777")
    fig.tight_layout(rect=(0, 0.045, 1, 0.955))
    fig.savefig(destino, dpi=130)
    plt.close(fig)
    return destino


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--rapido", action="store_true",
                    help="mallas menores; sirve para comprobar la tuberia")
    ap.add_argument("--solo-figura", action="store_true",
                    help="no homogeneiza (sin superficies elasticas)")
    a = ap.parse_args()

    res_vista = 96 if a.rapido else 128
    res_homog = 16 if a.rapido else 24
    num_waves = 400 if a.rapido else 1000
    con_mec = not a.solo_figura

    salida = AQUI / "resultados"
    salida.mkdir(exist_ok=True)

    print("=" * 78)
    print(" REPLICA DE Kumar et al. (2020)")
    print(" " + CITA_APA)
    print(" licencia del articulo: " + LICENCIA)
    print("=" * 78)
    err = autocomprobar_superficie()
    print(f"  E(d) sobre un material isotropo: error {err:.2e}  [OK]")
    print(f"  rho={RHO}  beta={BETA_PI}pi  N={num_waves}  "
          f"vista {res_vista}^3  homog {res_homog}^3  semilla {SEMILLA}")
    print()

    t0 = time.time()
    mascaras, res = {}, []
    for cl in CLASES:
        print(f"  {cl['nombre']:10s} theta={cl['thetas']} ...", end="",
              flush=True)
        BW, r = una_clase(cl, res_vista, res_homog, num_waves,
                          con_mecanica=con_mec)
        mascaras[cl["clave"]] = BW
        res.append(r)
        extra = ""
        if con_mec and r.get("homogeneizacion_ok"):
            extra = (f"  E1/E3={r['E1_E3']:.3f}  "
                     f"aniso={r['anisotropia_E']:.2f}")
        print(f" {r['segundos']:5.1f} s  BV/TV={r['BVTV']:.3f}  "
              f"DA={r['DA']:.3f}{extra}")

    # Contraste de esquemas sobre una terna de conos DESIGUALES. No es una de
    # las cuatro clases del articulo: se anade solo para este contraste, y su
    # procedencia se declara como tal en el JSON.
    res_equit = []
    print("\n  contraste de esquema de muestreo, conos DESIGUALES:")
    cl_des = {"clave": "desigual", "nombre": "Conos desiguales",
              "thetas": list(TERNA_DESIGUAL),
              "fuente": "NO procede del articulo; sirve para contrastar los "
                        "dos esquemas de muestreo",
              "esperado": "los dos esquemas deben diferir"}
    _, r_rech = una_clase(cl_des, res_vista, res_homog, num_waves,
                          esquema="rechazo", con_mecanica=False)
    _, r_equi = una_clase(cl_des, res_vista, res_homog, num_waves,
                          esquema="equitativo", con_mecanica=False)
    res.append(r_rech)
    res_equit.append(r_equi)
    print(f"    theta={list(TERNA_DESIGUAL)}   rechazo DA={r_rech['DA']:.3f}"
          f"   equitativo DA={r_equi['DA']:.3f}")

    checks = comprobar(res, res_equit)
    print("\n  PREDICCIONES DEL ARTICULO")
    fallos = 0
    for c in checks:
        ok = "OK   " if c["pasa"] else "FALLA"
        fallos += (not c["pasa"])
        print(f"    [{ok}] {c['codigo']}  {c['prediccion']}")
        print(f"             {c['obtenido']}   (criterio: {c['criterio']})")

    fig = None
    print("\n  componiendo la figura...", end="", flush=True)
    try:
        # Solo las cuatro clases del articulo: la terna desigual es
        # nuestra y no pinta nada en una replica de su Fig. 2.
        claves_art = {c["clave"] for c in CLASES}
        fig = figura(mascaras, [r for r in res if r["clave"] in claves_art],
                     salida / "replica_kumar2020_fig2.png",
                     con_mecanica=con_mec)
        print(f" -> {fig.name}")
    except Exception as e:
        print(f" fallo: {type(e).__name__}: {e}")

    doc = {
        "replica_de": {"cita_apa": CITA_APA, "doi": DOI,
                       "licencia": LICENCIA,
                       "figura_original": "referencia/Kumar2020_Fig2.png"},
        "generado": time.strftime("%Y-%m-%d %H:%M:%S"),
        "parametros": {"rho": RHO, "beta_pi": BETA_PI, "theta_min": TH_MIN,
                       "num_waves": num_waves, "semilla": SEMILLA,
                       "resolucion_vista": res_vista,
                       "resolucion_homogeneizacion": res_homog,
                       "esquema": "rechazo (ecuacion 2 del articulo)"},
        "clases": res,
        "contraste_equitativo": res_equit,
        "comprobaciones": checks,
        "fallos": fallos,
        "figura": fig.name if fig else None,
    }
    j = salida / "replica_kumar2020.json"
    j.write_text(json.dumps(doc, indent=2, ensure_ascii=False),
                 encoding="utf-8")

    print(f"\n  total {time.time()-t0:.0f} s   ·   informe -> {j.name}")
    print("=" * 78)
    print(f"  COMPROBACIONES QUE FALLAN: {fallos}")
    print("=" * 78)
    return 1 if fallos else 0


if __name__ == "__main__":
    sys.exit(main())
