"""
Bloque 32: contraste con ensayos fisicos (`spinpy.validacion_exp`).

Se verifica que las estadisticas dan la respuesta conocida en casos cuya
solucion se calcula a mano, y que el modulo NO inventa un veredicto cuando no
se declararon criterios.

CASOS DE RESPUESTA CONOCIDA
  (1) Acuerdo perfecto: CCC = 1, pendiente = 1, intercepto = 0, RMSE = 0.
  (2) CCC de Lin analitico: p = (1, 2, 3), m = (2, 3, 4) da 2*(2/3) /
      (2/3 + 2/3 + 1) = 4/7 exactamente (varianzas poblacionales).
  (3) Error de escala puro, m = 1,3 p: pendiente 1,3, sesgo relativo
      (p - m)/m = -3/13, factor de Bland-Altman 1,3 con desviacion nula,
      y la calibracion de E_s recupera k = 1,3 con error de validacion
      cruzada nulo.
  (4) Cobertura: con ruido gaussiano multiplicativo y semilla fija, el IC 95 %
      de la pendiente contiene a la pendiente verdadera en una fraccion de
      2000 repeticiones entre 0,93 y 0,97.
  (5) Veredicto: sin criterios es "no evaluado"; con criterios cumple o no
      cumple segun corresponda; una clave desconocida es un error.
  (6) Entrada: menos de 3 parejas, valores no finitos o largos distintos
      fallan; con n < 8 se avisa; el CSV sin 'medido' (la plantilla) no da
      ningun grupo y la linea de ordenes devuelve 2.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  (1) a (3): 1e-12 (aritmetica exacta salvo redondeo).
  (4): cobertura entre 0,93 y 0,97 (nominal 0,95; con 2000 repeticiones el
       error tipico de la cobertura es 0,005).
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from spinpy import validacion_exp as ve

RAIZ = Path(__file__).resolve().parent.parent
PLANTILLA = RAIZ / "datos" / "porcino_koria2020" / \
    "ensayos_experimentales_PLANTILLA.csv"


def test_acuerdo_perfecto():
    p = np.array([1.0, 2.0, 4.0, 7.0, 11.0])
    r = ve.acuerdo(p, p.copy())
    assert r["ccc"] == pytest.approx(1.0, abs=1e-12)
    assert r["regresion"]["pendiente"] == pytest.approx(1.0, abs=1e-12)
    assert r["regresion"]["intercepto"] == pytest.approx(0.0, abs=1e-12)
    assert r["rmse"] == pytest.approx(0.0, abs=1e-12)


def test_ccc_de_lin_analitico():
    assert ve.ccc_lin(np.array([1.0, 2, 3]), np.array([2.0, 3, 4])) == \
        pytest.approx(4 / 7, abs=1e-12)


def test_error_de_escala_puro():
    p = np.array([500.0, 900.0, 1500.0, 2400.0, 3100.0])
    r = ve.acuerdo(p, 1.3 * p)
    assert r["regresion"]["pendiente"] == pytest.approx(1.3, abs=1e-12)
    assert r["sesgo_rel"] == pytest.approx(-3 / 13, abs=1e-12)
    ba = r["bland_altman_log"]
    assert ba["factor_sesgo"] == pytest.approx(1.3, abs=1e-12)
    assert ba["sd_log"] == pytest.approx(0.0, abs=1e-12)
    cal = ve.calibrar_Es(p, 1.3 * p, Es_usado=20e9)
    assert cal["k"] == pytest.approx(1.3, rel=1e-12)
    assert cal["Es_calibrado"] == pytest.approx(26e9, rel=1e-12)
    assert cal["error_rel_loo_medio_abs"] == pytest.approx(0.0, abs=1e-12)


def test_calibracion_distingue_error_de_forma():
    """Si lo medido NO es proporcional a lo predicho, la validacion cruzada lo delata."""
    p = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    m = np.array([3.0, 1.0, 5.0, 2.0, 6.0, 4.0])      # sin relacion lineal
    cal = ve.calibrar_Es(p, m)
    assert cal["error_rel_loo_medio_abs"] > 0.3


def test_cobertura_del_ic_de_la_pendiente():
    rng = np.random.default_rng(20261013)
    verdadera, n, reps, dentro = 1.2, 12, 2000, 0
    for _ in range(reps):
        p = rng.uniform(500, 3000, n)
        m = verdadera * p + rng.normal(0, 150, n)
        lo, hi = ve.acuerdo(p, m)["regresion"]["pendiente_ic"]
        dentro += lo <= verdadera <= hi
    assert 0.93 <= dentro / reps <= 0.97


def test_veredicto_sin_criterios_es_no_evaluado():
    p = np.array([1.0, 2, 3, 4, 5, 6, 7, 8])
    r = ve.acuerdo(p, p * 1.05)
    assert r["veredicto"] == "no evaluado" and r["criterios"] is None


def test_veredicto_con_criterios():
    p = np.array([1.0, 2, 3, 4, 5, 6, 7, 8])
    m = p * 1.05 + np.array([0.01, -0.02, 0.0, 0.02, -0.01, 0.01, 0.0, -0.01])
    ok = ve.acuerdo(p, m, criterios={"ccc_min": 0.95, "mape_max": 0.10,
                                     "pendiente_ic_incluye_1": False})
    assert ok["veredicto"] == "cumple"
    mal = ve.acuerdo(p, 2 * p + 0.3, criterios={"sesgo_rel_max": 0.10,
                                                "pendiente_ic_incluye_1": True})
    assert mal["veredicto"] == "no cumple"
    assert any(d.startswith("NO cumple") for d in mal["detalle_veredicto"])


def test_criterio_desconocido_falla():
    p = np.arange(1.0, 9.0)
    with pytest.raises(ValueError, match="criterios desconocidos"):
        ve.acuerdo(p, p, criterios={"r2_min": 0.9})


def test_entradas_invalidas():
    with pytest.raises(ValueError):
        ve.acuerdo([1, 2], [1, 2])
    with pytest.raises(ValueError):
        ve.acuerdo([1, 2, 3], [1, 2])
    with pytest.raises(ValueError):
        ve.acuerdo([1, 2, np.nan], [1, 2, 3])


def test_aviso_con_pocas_parejas():
    p = np.array([1.0, 2, 3, 4, 5])
    r = ve.acuerdo(p, p * 1.1)
    assert any("n = 5" in a for a in r["avisos"])


def test_plantilla_sin_rellenar_no_da_grupos(capsys):
    assert PLANTILLA.exists()
    assert ve.leer_csv(PLANTILLA) == {}
    assert ve.main([str(PLANTILLA)]) == 2


def test_csv_rellenado(tmp_path):
    f = tmp_path / "e.csv"
    f.write_text("id,magnitud,predicho,medido\n"
                 "a,E_app_MPa,1000,1100\nb,E_app_MPa,2000,2100\n"
                 "c,E_app_MPa,3000,3300\nd,E_app_MPa,4000,\n",
                 encoding="utf-8")
    g = ve.leer_csv(f)
    assert list(g) == ["E_app_MPa"] and len(g["E_app_MPa"][0]) == 3
    assert ve.main([str(f)]) == 0
