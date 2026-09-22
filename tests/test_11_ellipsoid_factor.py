"""
Bloque 11 — Ellipsoid Factor (placas frente a barras).

REFERENCIAS
  Doube M. The ellipsoid factor for quantification of rods, plates, and
  intermediate forms in 3D geometries. Front Endocrinol 2015;6:15.
    EF = a/b - b/c con a <= b <= c los semiejes del mayor elipsoide inscrito
    que contiene el punto. Placa -1, esfera 0, barra +1.
  Vafaeefar M et al. Gyroid, spinodoid and dual-lattice algorithms as
  structural models of trabecular bone. JMBBM 2022.
    Criterio de validez: relleno > 95 % del solido.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  La discretizacion es por centros de voxel y el exterior del cubo cuenta como
  fondo (ver la cabecera de `spinpy/elipsoide.py`). Con esa regla los valores
  esperados NO son los ideales -1, 0, +1 sino los de la figura finita:

    losa  40^3, espesor 6    el mayor elipsoide es el de la caja 7 x 41 x 41
                             (medio voxel mas por cada lado hasta el centro de
                             fondo): EF = 7/41 - 1 = -0.829.
                             Mediana en [-1.00, -0.70].
    cilindro d = 8, largo 60 radio hasta el primer centro de fondo ~4.1,
                             semilongitud 30.5: EF = 1 - 4.1/30.5 = +0.866.
                             Mediana en [+0.65, +1.00].
    cilindro OBLICUO d = 8   el mismo, orientado en (1,1,0) dentro de 48^3.
                             Semilongitud menor y cortado por las caras, asi
                             que se pide menos: mediana >= +0.60. Es la prueba
                             de la busqueda de orientacion: un metodo con los
                             elipsoides alineados a los ejes daria ~0.
    esfera r = 12            |mediana| <= 0.15 (esfera digital, no perfecta).
    relleno                  >= 0.95 en las cuatro (el criterio de Vafaeefar).
    reproducible             la misma mascara da exactamente el mismo EF.
    clases de espinodoide    lamelar (0,0,15) con EF MENOR que columnar
                             (15,15,0) a igual densidad: es lo que las dos
                             clases son por construccion (Kumar 2020). Solo se
                             declara el orden, no los valores.

  Las bandas son anchas a proposito: lo que se verifica es que el metodo
  distingue placa, esfera y barra con el signo correcto y cerca del valor de la
  figura finita. Un EF de placa POSITIVO seria el mismo tipo de error que el
  signo de dilatacion del SMI (bloque 03): no salta a la vista y lo invierte
  todo.
"""
import numpy as np

from conftest import bola, cilindro, losa
from spinpy.elipsoide import factor_elipsoide

BLOQUE = "11 Ellipsoid Factor (Doube)"
REF = "Doube 2015; Vafaeefar et al. 2022"
SP = [1.0, 1.0, 1.0]


def cilindro_oblicuo(n=48, diametro=8.0):
    """Cilindro de eje (1,1,0)/sqrt(2) que pasa por el centro del cubo."""
    g = np.indices((n, n, n), dtype=float) - (n - 1) / 2.0
    e = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    t = g[0] * e[0] + g[1] * e[1]
    d2 = (g ** 2).sum(0) - t ** 2
    return d2 <= (diametro / 2.0) ** 2


def _relleno(r, etiqueta, registro):
    ok = r["EF_relleno"] >= 0.95
    registro.anotar(BLOQUE, f"{etiqueta}: relleno", REF, None,
                    r["EF_relleno"], ">= 0.95", "criterio de Vafaeefar", ok)
    return ok


def test_losa_ef_placa(registro):
    r = factor_elipsoide(losa(n=40, espesor=6), SP)
    ok = -1.0 <= r["EF_mediana"] <= -0.70
    registro.anotar(BLOQUE, "losa e=6: EF mediana", REF, 7 / 41 - 1,
                    r["EF_mediana"], "[-1.00, -0.70]",
                    "placa; caja finita 7x41x41", ok)
    assert _relleno(r, "losa", registro), r
    assert ok, r


def test_cilindro_ef_barra(registro):
    r = factor_elipsoide(cilindro(n_xy=20, largo=60, diametro=8.0), SP)
    ok = 0.65 <= r["EF_mediana"] <= 1.0
    registro.anotar(BLOQUE, "cilindro d=8 L=60: EF mediana", REF,
                    1 - 4.1 / 30.5, r["EF_mediana"], "[+0.65, +1.00]",
                    "barra; semilongitud acotada por el cubo", ok)
    assert _relleno(r, "cilindro", registro), r
    assert ok, r


def test_cilindro_oblicuo_orientacion(registro):
    r = factor_elipsoide(cilindro_oblicuo(n=48, diametro=8.0), SP)
    ok = r["EF_mediana"] >= 0.60
    registro.anotar(BLOQUE, "cilindro oblicuo (1,1,0): EF mediana", REF, None,
                    r["EF_mediana"], ">= +0.60",
                    "busqueda de orientacion; alineado a ejes daria ~0", ok)
    assert _relleno(r, "cilindro oblicuo", registro), r
    assert ok, r


def test_esfera_ef_cero(registro):
    r = factor_elipsoide(bola(n=32, radio=12.0), SP)
    ok = abs(r["EF_mediana"]) <= 0.15
    registro.anotar(BLOQUE, "esfera r=12: EF mediana", REF, 0.0,
                    r["EF_mediana"], "|EF| <= 0.15", "esfera ideal = 0", ok)
    assert _relleno(r, "esfera", registro), r
    assert ok, r


def test_reproducible(registro):
    BW = cilindro(n_xy=16, largo=32, diametro=6.0)
    a = factor_elipsoide(BW, SP)
    b = factor_elipsoide(BW, SP)
    ok = a["EF"] == b["EF"] and a["EF_relleno"] == b["EF_relleno"]
    registro.anotar(BLOQUE, "misma mascara, mismo EF", REF, a["EF"], b["EF"],
                    "identico", "semilla fija 20260720", ok)
    assert ok, (a, b)


def test_clases_espinodoide(registro):
    from spinpy.grf import campo_grf, level_set

    def ef(thetas):
        G, _, _ = campo_grf(48, 8.0 * np.pi, 1000, list(thetas), seed=20260720)
        return factor_elipsoide(G <= level_set(0.35), np.full(3, 1.0 / 48))

    lam = ef((0.0, 0.0, 15.0))
    col = ef((15.0, 15.0, 0.0))
    ok = lam["EF_mediana"] < col["EF_mediana"]
    registro.anotar(BLOQUE, "lamelar (0,0,15): EF mediana", REF, None,
                    lam["EF_mediana"], "< columnar", "Kumar 2020, clases", ok)
    registro.anotar(BLOQUE, "columnar (15,15,0): EF mediana", REF, None,
                    col["EF_mediana"], "> lamelar", "Kumar 2020, clases", ok)
    assert ok, (lam, col)
