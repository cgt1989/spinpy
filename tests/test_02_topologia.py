"""
Bloque 2 — Topologia: numero de Euler y densidad de conectividad.

REFERENCIAS
  Odgaard A, Gundersen HJG. Quantification of connectivity in cancellous bone,
  with special emphasis on 3-D reconstructions. Bone 1993;14(2):173-182.
  Formula de Euler-Poincare: chi = b0 - b1 + b2 (componentes - tuneles +
  cavidades).

POR QUE ES LA VERIFICACION MAS LIMPIA DE TODAS
  El numero de Euler es un INVARIANTE TOPOLOGICO ENTERO. No hay interpolacion,
  ni azar, ni sesgo de discretizacion que valga: una bola tiene chi = 1, un
  toro chi = 0, y un reticulo con V nodos y E barras tiene chi = V - E. O sale
  el entero exacto o la implementacion esta mal.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  Exactas. Cero error en chi y en Conn = 1 - chi.

LO QUE LA FORMULA DE ODGAARD PRESUPONE, Y SE COMPRUEBA A PROPOSITO
  Conn = 1 - chi solo cuenta conexiones redundantes si hay UNA componente y
  NINGUNA cavidad cerrada (b0 = 1, b2 = 0). Con dos componentes o con una
  burbuja, chi sube y Conn sale negativo. No es un fallo: es la hipotesis de
  la definicion, y conviene que quede demostrada en una prueba para que nadie
  lea un Conn.D negativo como "menos conectado".
"""
import numpy as np

from conftest import bola, cascara, dos_bolas, reticulo, toro
from spinpy.morphometry import conectividad

BLOQUE = "02 Topologia (Euler, Conn.D)"
REF = "Odgaard & Gundersen 1993, Bone 14:173; Euler-Poincare chi = b0-b1+b2"


def _caso(registro, nombre, BW, chi_esp, nota=""):
    r = conectividad(BW, [1.0, 1.0, 1.0])
    conn_esp = 1.0 - chi_esp
    registro.anotar(BLOQUE, f"{nombre}: chi", REF, chi_esp, r["euler"],
                    "exacta", "entero exacto", r["euler"] == chi_esp, nota)
    registro.anotar(BLOQUE, f"{nombre}: Conn = 1 - chi", REF, conn_esp,
                    r["Conn"], "exacta", "entero exacto",
                    r["Conn"] == conn_esp, nota)
    assert r["euler"] == chi_esp, f"{nombre}: chi={r['euler']} != {chi_esp}"
    assert r["Conn"] == conn_esp
    return r


def test_bola_maciza(registro):
    """b0=1, b1=0, b2=0 -> chi=1, ninguna conexion redundante."""
    _caso(registro, "bola maciza", bola(), 1)


def test_toro_macizo(registro):
    """b1=1 -> chi=0, exactamente UNA conexion redundante."""
    _caso(registro, "toro macizo", toro(), 0)


def test_reticulo_de_barras(registro):
    """V nodos y E barras: chi = V - E, redundantes b1 = E - V + 1."""
    BW, t = reticulo(nodos=3, paso=10, grosor=3)
    r = _caso(registro, f"reticulo {t['V']} nodos / {t['E']} barras", BW,
              t["euler"],
              nota=f"b1 = E - V + 1 = {t['b1']} conexiones redundantes")
    assert r["Conn"] == t["b1"]


def test_dos_componentes_delatan_la_hipotesis(registro):
    """Dos bolas: chi=2 y Conn=-1. Negativo = mas de una componente."""
    _caso(registro, "dos bolas separadas", dos_bolas(), 2,
          nota="Conn negativo: la formula presupone b0 = 1")


def test_cavidad_delata_la_hipotesis(registro):
    """Cascara hueca: b2=1 -> chi=2, Conn=-1. Negativo = cavidad cerrada."""
    _caso(registro, "cascara esferica hueca", cascara(), 2,
          nota="Conn negativo: la formula presupone b2 = 0")


def test_conn_d_escala_con_el_volumen(registro):
    """Conn.D = Conn / TV: con la mitad de voxel sale ocho veces mas."""
    BW = toro()
    a = conectividad(BW, [1.0, 1.0, 1.0])["ConnD"]
    b = conectividad(BW, [0.5, 0.5, 0.5])["ConnD"]
    registro.anotar(BLOQUE, "Conn.D con spacing/2 = 8 x Conn.D", REF, 8.0,
                    b / a, 1e-12, "cociente exacto", abs(b / a - 8.0) < 1e-12)
    assert abs(b / a - 8.0) < 1e-12
