"""
Bloque 3 — Structure Model Index.

REFERENCIAS
  Hildebrand T, Ruegsegger P. Quantification of bone microarchitecture with
  the structure model index. Comput Methods Biomech Biomed Engin 1997;1:15-23.
    Valores ideales: placa 0, cilindro 3, esfera 4.
  Salmon PL, Ohlsson C, Shefelbine SJ, Doube M. Structure Model Index does not
  measure rods and plates in trabecular bone. Front Endocrinol 2015;6:162.
    Las superficies CONCAVAS aportan area negativa al dilatar: SMI < 0.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  El calculo va sobre la isosuperficie de marching cubes, que sobreestima el
  area de una figura binaria un ~8.5 % (medido en validar_analitico.py). Como
  SMI ~ 1/BS^2, cabe esperar un deficit del orden del 15 % en las figuras
  curvas, y ninguno en la placa (su dBS/dr es cero). Por tanto:

    placa      |SMI| <= 0.05
    cilindro   2.40 <= SMI <= 3.00     (hasta -20 %, nunca por encima)
    esfera     3.20 <= SMI <= 4.00     (idem)
    cavidad    SMI < 0                 (Salmon 2015: concavidad)
    signo      dBS/dr > 0 para una esfera (la dilatacion es hacia fuera)

  La cota superior es tan informativa como la inferior: si el SMI saliera POR
  ENCIMA del ideal, el sesgo de marching cubes estaria entrando con el signo
  cambiado, que es exactamente el error que el codigo dice haber corregido.
"""
from conftest import bola, cavidad, cilindro, losa
from spinpy.morphometry import indice_smi

BLOQUE = "03 SMI (Hildebrand & Ruegsegger)"
REF = "Hildebrand & Ruegsegger 1997; Salmon et al. 2015"
SP = [1.0, 1.0, 1.0]


def test_placa_smi_cero(registro):
    r = indice_smi(losa(n=40, espesor=10), SP)
    ok = abs(r["SMI"]) <= 0.05
    registro.anotar(BLOQUE, "losa: SMI", REF, 0.0, r["SMI"], "|SMI| <= 0.05",
                    "placa ideal = 0; dBS/dr es nulo", ok)
    assert ok, r


def test_cilindro_smi_tres(registro):
    r = indice_smi(cilindro(n_xy=40, largo=40, diametro=16.0), SP)
    ok = 2.40 <= r["SMI"] <= 3.00
    registro.anotar(BLOQUE, "cilindro d=16: SMI", REF, 3.0, r["SMI"],
                    "2.40 <= SMI <= 3.00", "ideal 3; deficit <= 20 % por MC",
                    ok)
    assert ok, r


def test_esfera_smi_cuatro(registro):
    r = indice_smi(bola(n=48, radio=14.0), SP)
    ok = 3.20 <= r["SMI"] <= 4.00
    registro.anotar(BLOQUE, "esfera d=28: SMI", REF, 4.0, r["SMI"],
                    "3.20 <= SMI <= 4.00", "ideal 4; deficit <= 20 % por MC",
                    ok)
    assert ok, r
    registro.anotar(BLOQUE, "esfera: dBS/dr > 0 (dilata hacia fuera)", REF,
                    None, r["dBSdr"], "> 0", "signo de la dilatacion",
                    r["dBSdr"] > 0)
    assert r["dBSdr"] > 0


def test_cavidad_smi_negativo(registro):
    """Bloque con una burbuja: la superficie es concava y el SMI sale < 0."""
    r = indice_smi(cavidad(n=40, radio=10.0), SP)
    ok = r["SMI"] < 0
    registro.anotar(BLOQUE, "cavidad esferica: SMI < 0 (concavidad)", REF,
                    None, r["SMI"], "< 0", "Salmon 2015", ok,
                    nota="es la razon por la que el SMI no es fiable a BV/TV alto")
    assert ok, r
