"""
Bloque 4 — Espesor local por esferas maximas inscritas.

REFERENCIA
  Hildebrand T, Ruegsegger P. A new method for the model-independent
  assessment of thickness in three-dimensional images. J Microsc
  1997;185:67-75.
    El espesor local en un punto es el DIAMETRO de la mayor esfera que lo
    contiene y cabe en la estructura. Para una losa de grosor t vale t en
    todo punto; para un cilindro de diametro d vale d; para una bola, d.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  El docstring de `espesor.py` afirma dos cosas medibles:
    (a) "una losa plana sale EXACTA (0.00 %)"
    (b) para un cilindro, el error por diametro en voxeles es
          d = 2: 0 %   3: -33 %   4: -29 %   6: -5.7 %   10: -11 %   20: -5.7 %
  Se toman esas afirmaciones como prediccion:
    losa       |Th - t| <= 0.5 voxel   (la afirmacion (a), con medio voxel de
                                         margen por la discretizacion)
    cilindro   d >= 6:  |err| <= 15 %    (la tabla (b) con margen)
               d <= 4:  |err| <= 40 %    (idem; la zona donde el sesgo es grande)
    bola       |err| <= 15 %
  Y ademas la condicion de que el metodo esta bien planteado:
    el espesor NUNCA supera el diametro de la figura + 1 voxel, porque no cabe
    una esfera mayor que la propia figura.

  Si (a) falla, es un hallazgo sobre la propia documentacion del modulo, no
  un fallo de la suite, y asi se reporta.

RESULTADO DE LA PRIMERA EJECUCION, DEJADO AQUI A PROPOSITO
  (a) FALLA: la losa de t voxeles devuelve t+1 en todas las capas, para todo
      t. La transformada de distancia mide entre CENTROS de voxel, y el
      centro del primer voxel de fondo esta a un voxel entero del ultimo
      centro solido, no a medio: r_max = (t+1)/2 y Th = t+1. El docstring
      dice lo contrario ("sale corto"): para una figura plana sale LARGO.
  (b) PASA, pero la tabla del docstring no es una tabla de sesgo: es una
      tabla de ALINEACION. Con el eje del cilindro ENTRE voxeles (n_xy par)
      la mayor esfera inscrita queda corta; con el eje SOBRE un voxel (n_xy
      impar) hereda el +1 de la losa y sale larga. Lo mide la ultima prueba.
  Tolerancia del escalado con spacing: el campo se guarda en float32 por
  diseno, asi que el criterio coherente es 1e-6, no 1e-9 (corregido).
"""
import numpy as np

from conftest import bola, cilindro, losa
from spinpy.espesor import espesor_local

BLOQUE = "04 Espesor local (esferas inscritas)"
REF = "Hildebrand & Ruegsegger 1997, J Microsc 185:67"
SP = [1.0, 1.0, 1.0]


def _mediana_solido(esp, BW):
    v = esp[BW]
    return float(np.median(v)), float(v.max())


def test_losa_es_exacta(registro):
    fallos = []
    for t in (3, 5, 9, 15):
        BW = losa(n=32, espesor=t)
        med, mx = _mediana_solido(espesor_local(BW, SP), BW)
        ok = abs(med - t) <= 0.5
        registro.anotar(BLOQUE, f"losa t={t}: mediana Th", REF, t, med,
                        "|Th - t| <= 0.5 vox", "afirmacion 'losa exacta'", ok,
                        nota=f"max = {mx:.2f}")
        if not ok:
            fallos.append((t, med))
    assert not fallos, f"la losa NO sale exacta: (t, mediana) = {fallos}"


def test_cilindro_sesgo_acotado(registro):
    fallos = []
    for d, tol in ((4, 0.40), (6, 0.15), (10, 0.15), (20, 0.15)):
        BW = cilindro(n_xy=d + 12, largo=24, diametro=float(d))
        med, mx = _mediana_solido(espesor_local(BW, SP), BW)
        err = (med - d) / d
        ok = abs(err) <= tol
        registro.anotar(BLOQUE, f"cilindro d={d}: mediana Th", REF, d, med,
                        f"|err| <= {tol:.0%}", "tabla del docstring", ok,
                        nota=f"max = {mx:.2f}")
        if not ok:
            fallos.append((d, med, err))
        # nunca mayor que la figura
        assert mx <= d + 1.0 + 1e-9, f"Th max {mx} > d + 1 en cilindro d={d}"
    assert not fallos, fallos


def test_bola(registro):
    d = 20.0
    BW = bola(n=32, radio=d / 2)
    med, mx = _mediana_solido(espesor_local(BW, SP), BW)
    err = (med - d) / d
    ok = abs(err) <= 0.15
    registro.anotar(BLOQUE, "bola d=20: mediana Th", REF, d, med,
                    "|err| <= 15 %", "diametro de la bola", ok,
                    nota=f"max = {mx:.2f}")
    assert mx <= d + 1.0 + 1e-9
    assert ok, (med, err)


def test_espesor_escala_con_spacing(registro):
    """El campo esta en unidades fisicas: con spacing 0.05 sale 0.05 veces."""
    BW = cilindro(n_xy=24, largo=16, diametro=10.0)
    a = float(np.median(espesor_local(BW, [1.0] * 3)[BW]))
    b = float(np.median(espesor_local(BW, [0.05] * 3)[BW]))
    # El campo es float32 por diseno: la precision relativa es ~1.2e-7.
    ok = abs(b / a - 0.05) / 0.05 < 1e-6
    registro.anotar(BLOQUE, "Th escala linealmente con el spacing", REF, 0.05,
                    b / a, 1e-6, "cociente; campo float32", ok)
    assert ok


def test_alineacion_con_la_rejilla_cambia_el_signo(registro):
    """El mismo cilindro, eje entre voxeles o sobre un voxel: signo opuesto.

    No es una prueba de exactitud sino de DIAGNOSTICO: deja registrado que el
    error del espesor no es una funcion del diametro sino de donde cae la
    figura en la rejilla. Con eje sobre voxel se hereda el +1 de la losa.
    """
    d = 10.0
    filas = []
    for n_xy, donde in ((22, "eje ENTRE voxeles"), (23, "eje SOBRE un voxel")):
        BW = cilindro(n_xy=n_xy, largo=16, diametro=d)
        med, mx = _mediana_solido(espesor_local(BW, SP), BW)
        err = (med - d) / d
        registro.anotar(BLOQUE, f"cilindro d=10, {donde}: mediana Th", REF, d,
                        med, "|err| <= 40 %", "diagnostico de alineacion",
                        abs(err) <= 0.40, nota=f"error {err:+.1%}")
        filas.append(err)
        assert abs(err) <= 0.40
    # Lo que se documenta: entre voxeles queda corto, sobre voxel queda largo.
    registro.anotar(BLOQUE, "signo del error: entre voxeles < 0 < sobre voxel",
                    REF, None, filas[1] - filas[0], "> 0",
                    "el sesgo depende de la alineacion",
                    filas[0] < 0 < filas[1],
                    nota="entre %+.1f%%, sobre %+.1f%%" % (100 * filas[0],
                                                          100 * filas[1]))
    assert filas[0] < 0 < filas[1], filas
