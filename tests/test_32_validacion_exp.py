"""
Bloque 32: contraste con ensayos fisicos (`spinpy.validacion_exp`).

Se verifica que las estadisticas dan la respuesta conocida en casos cuya
solucion se calcula a mano, y que el modulo NO inventa un veredicto cuando no
se declararon criterios o cuando un criterio no se puede calcular.

CASOS DE RESPUESTA CONOCIDA
  (1) Acuerdo perfecto: CCC = 1, pendiente = 1, intercepto = 0, RMSE = 0.
  (2) CCC de Lin analitico: p = (1, 2, 3), m = (2, 3, 4) da 2*(2/3) /
      (2/3 + 2/3 + 1) = 4/7 exactamente (varianzas poblacionales).
  (3) Error de escala puro, m = 1,3 p: pendiente 1,3, sesgo relativo
      (p - m)/m = -3/13, factor de Bland-Altman 1,3 con desviacion nula,
      y la calibracion de E_s recupera k = 1,3 con error de validacion
      cruzada nulo.
  (4) Bland-Altman a mano: m/p = (1, 2, 4) da d = ln(m/p) = (0, ln 2, ln 4),
      media ln 2 y desviacion muestral ln 2; limites 2^(1 -+ z) y IC con
      t(2 g.l.) y el error estandar de Bland y Altman (1999).
  (5) MAPE, sesgo relativo, Spearman e intercepto con valores de mano.
  (6) Cobertura simulada (semilla fija, 2000 repeticiones) del IC 95 % de la
      pendiente, y del IC del CCC.
  (7) Veredicto: sin criterios "no evaluado"; con criterios cumple o no
      cumple; la equivalencia de la pendiente NO se cumple con un IC ancho;
      un criterio que da NaN es "no evaluable"; una clave desconocida falla.
  (8) Entrada: pocas parejas, no finitos, largos distintos; avisos con n
      pequeño y con medido <= 0; CSV con BOM, punto y coma, coma decimal y
      errores con numero de fila; JSON sin NaN.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  (1) a (5): 1e-12 (aritmetica exacta salvo redondeo).
  (6): cobertura de la pendiente entre 0,93 y 0,97 (nominal 0,95; error tipico
       con 2000 repeticiones 0,005). Cobertura del IC del CCC entre 0,90 y
       0,99 con n = 40: es una aproximacion asintotica de Lin (1989).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from spinpy import validacion_exp as ve

RAIZ = Path(__file__).resolve().parent.parent
PLANTILLA = RAIZ / "datos" / "porcino_koria2020" / \
    "ensayos_experimentales_PLANTILLA.csv"
Z = stats.norm.ppf(0.975)


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
    assert cal["k_geometrico"] == pytest.approx(1.3, rel=1e-12)
    assert cal["Es_calibrado"] == pytest.approx(26e9, rel=1e-12)
    assert cal["error_rel_loo_medio_abs"] == pytest.approx(0.0, abs=1e-12)


def test_calibracion_distingue_error_de_forma():
    """Si lo medido NO es proporcional a lo predicho, la validacion cruzada lo delata."""
    p = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    m = np.array([3.0, 1.0, 5.0, 2.0, 6.0, 4.0])      # sin relacion lineal
    cal = ve.calibrar_Es(p, m)
    assert cal["error_rel_loo_medio_abs"] > 0.3


def test_calibracion_exige_positivos():
    with pytest.raises(ValueError, match="positivos"):
        ve.calibrar_Es([1.0, 2.0, 3.0], [1.0, 0.0, 3.0])


def test_bland_altman_log_a_mano():
    p = np.array([1.0, 1.0, 1.0])
    m = np.array([1.0, 2.0, 4.0])
    ba = ve.acuerdo(p, m)["bland_altman_log"]
    ln2, n = math.log(2.0), 3
    assert ba["factor_sesgo"] == pytest.approx(2.0, abs=1e-12)
    assert ba["sd_log"] == pytest.approx(ln2, abs=1e-12)
    assert ba["limite_inf"] == pytest.approx(math.exp(ln2 - Z * ln2), abs=1e-12)
    assert ba["limite_sup"] == pytest.approx(math.exp(ln2 + Z * ln2), abs=1e-12)
    t = stats.t.ppf(0.975, n - 1)
    ee_sesgo = ln2 / math.sqrt(n)
    assert ba["factor_sesgo_ic"][0] == pytest.approx(
        math.exp(ln2 - t * ee_sesgo), abs=1e-12)
    ee_lim = ln2 * math.sqrt(1 / n + Z ** 2 / (2 * (n - 1)))
    assert ba["limite_sup_ic"][1] == pytest.approx(
        math.exp(ln2 + Z * ln2 + t * ee_lim), abs=1e-12)
    assert ba["limite_inf_ic"][0] == pytest.approx(
        math.exp(ln2 - Z * ln2 - t * ee_lim), abs=1e-12)


def test_mape_sesgo_spearman_e_intercepto_a_mano():
    r = ve.acuerdo([110.0, 180.0, 330.0], [100.0, 200.0, 300.0])
    assert r["mape"] == pytest.approx(0.1, abs=1e-12)
    assert r["sesgo_rel"] == pytest.approx(0.1 / 3, abs=1e-12)
    r2 = ve.acuerdo([1.0, 2.0, 3.0, 4.0], [10.0, 30.0, 20.0, 40.0])
    assert r2["spearman"] == pytest.approx(0.8, abs=1e-12)
    p = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    r3 = ve.acuerdo(p, 2 * p + 5)
    assert r3["regresion"]["intercepto"] == pytest.approx(5.0, abs=1e-9)


def test_cobertura_del_ic_de_la_pendiente():
    rng = np.random.default_rng(20261013)
    verdadera, n, reps, dentro = 1.2, 12, 2000, 0
    for _ in range(reps):
        p = rng.uniform(500, 3000, n)
        m = verdadera * p + rng.normal(0, 150, n)
        lo, hi = ve.acuerdo(p, m)["regresion"]["pendiente_ic"]
        dentro += lo <= verdadera <= hi
    assert 0.93 <= dentro / reps <= 0.97


def test_cobertura_del_ic_del_ccc():
    """El CCC poblacional se calcula con los momentos verdaderos del modelo."""
    rng = np.random.default_rng(7)
    n, reps, dentro = 40, 1500, 0
    # p ~ U(500, 3000); m = 1,2 p + e, e ~ N(0, 250^2)
    vp = (3000 - 500) ** 2 / 12
    mp = 1750.0
    mm, vm = 1.2 * mp, 1.44 * vp + 250 ** 2
    cov = 1.2 * vp
    ccc_pob = 2 * cov / (vp + vm + (mp - mm) ** 2)
    for _ in range(reps):
        p = rng.uniform(500, 3000, n)
        m = 1.2 * p + rng.normal(0, 250, n)
        lo, hi = ve.acuerdo(p, m)["ccc_ic"]
        dentro += lo <= ccc_pob <= hi
    assert 0.90 <= dentro / reps <= 0.99


def test_veredicto_sin_criterios_es_no_evaluado():
    p = np.array([1.0, 2, 3, 4, 5, 6, 7, 8])
    r = ve.acuerdo(p, p * 1.05)
    assert r["veredicto"] == "no evaluado" and r["criterios"] is None


def test_veredicto_con_criterios():
    p = np.linspace(1000, 4000, 20)
    ruido = np.array([10, -20, 0, 20, -10, 10, 0, -10, 15, -15,
                      5, -5, 20, -20, 10, -10, 0, 5, -5, 0.0])
    m = p * 1.02 + ruido
    ok = ve.acuerdo(p, m, criterios={"ccc_min": 0.95, "mape_max": 0.10,
                                     "pendiente_equiv_delta": 0.10})
    assert ok["veredicto"] == "cumple", ok["detalle_veredicto"]
    mal = ve.acuerdo(p, 2 * p + 300, criterios={"sesgo_rel_max": 0.10,
                                                "pendiente_equiv_delta": 0.10})
    assert mal["veredicto"] == "no cumple"
    assert any(d.startswith("NO cumple") for d in mal["detalle_veredicto"])


def test_equivalencia_no_se_cumple_con_ic_ancho():
    """Pocos datos y mucho ruido: el IC de la pendiente es ancho y no prueba equivalencia."""
    p = np.array([1000.0, 1500.0, 2200.0, 3000.0])
    m = np.array([1300.0, 1100.0, 2600.0, 2700.0])
    r = ve.acuerdo(p, m, criterios={"pendiente_equiv_delta": 0.10})
    lo, hi = r["regresion"]["pendiente_ic"]
    assert hi - lo > 0.2
    assert r["veredicto"] == "no cumple"


def test_criterio_que_da_nan_es_no_evaluable():
    r = ve.acuerdo([1.0, 2.0, 3.0], [0.0, 2.0, 3.0], criterios={"mape_max": 0.1})
    assert r["veredicto"] == "no evaluable"
    assert any("<= 0" in a for a in r["avisos"])
    c = ve.acuerdo([5.0, 5.0, 5.0], [5.0, 5.0, 5.0], criterios={"ccc_min": 0.9})
    assert c["veredicto"] == "no evaluable"
    assert any("CCC no definido" in a for a in c["avisos"])


def test_criterio_desconocido_falla():
    p = np.arange(1.0, 9.0)
    with pytest.raises(ValueError, match="criterios desconocidos"):
        ve.acuerdo(p, p, criterios={"r2_min": 0.9})
    with pytest.raises(ValueError, match="criterios desconocidos"):
        ve.acuerdo(p, p, criterios={"pendiente_ic_incluye_1": True})


def test_entradas_invalidas():
    with pytest.raises(ValueError):
        ve.acuerdo([1, 2], [1, 2])
    with pytest.raises(ValueError):
        ve.acuerdo([1, 2, 3], [1, 2])
    with pytest.raises(ValueError):
        ve.acuerdo([1, 2, np.nan], [1, 2, 3])


def test_avisos_con_pocas_parejas():
    p = np.array([1.0, 2, 3, 4, 5])
    r = ve.acuerdo(p, p * 1.1 + np.array([0.1, -0.1, 0.05, 0.0, -0.05]))
    assert any("n = 5" in a and "no permite concluir" in a for a in r["avisos"])
    assert any("anticonservador" in a for a in r["avisos"])


def test_plantilla_sin_rellenar_no_da_grupos():
    assert PLANTILLA.exists()
    assert ve.leer_csv(PLANTILLA) == {}
    assert ve.main([str(PLANTILLA)]) == 2


def test_csv_rellenado_y_json(tmp_path):
    f = tmp_path / "e.csv"
    f.write_text("id,magnitud,predicho,medido\n"
                 "a,E_app_MPa,1000,1100\nb,E_app_MPa,2000,2100\n"
                 "c,E_app_MPa,3000,3300\nd,E_app_MPa,4000,\n",
                 encoding="utf-8")
    g = ve.leer_csv(f)
    assert list(g) == ["E_app_MPa"] and g["E_app_MPa"][0] == ["a", "b", "c"]
    salida = tmp_path / "r.json"
    assert ve.main([str(f), "--json", str(salida), "--Es-usado", "20000",
                    "--mape-max", "0.2"]) == 0
    datos = json.loads(salida.read_text(encoding="utf-8"))
    assert datos[0]["n"] == 3 and datos[0]["calibracion"]["Es_usado"] == 20000.0
    assert datos[0]["veredicto"] in ("cumple", "no cumple")


def test_csv_excel_en_español_con_bom(tmp_path):
    f = tmp_path / "x.csv"
    f.write_bytes("id;magnitud;predicho;medido\r\n"
                  "a;E_app_MPa;1000,5;1100,25\r\n"
                  "b;E_app_MPa;2000;2100\r\n".encode("utf-8-sig"))
    g = ve.leer_csv(f)
    assert g["E_app_MPa"][1] == [1000.5, 2000.0]
    assert g["E_app_MPa"][2] == [1100.25, 2100.0]


def test_csv_error_indica_la_fila(tmp_path):
    f = tmp_path / "x.csv"
    f.write_text("id,magnitud,predicho,medido\n"
                 "a,E_app_MPa,1000,1100\nb,E_app_MPa,,2100\n", encoding="utf-8")
    with pytest.raises(ValueError, match=r"fila 3 \(id 'b'\)"):
        ve.leer_csv(f)


def test_json_sin_nan(tmp_path):
    f = tmp_path / "x.csv"
    f.write_text("id,magnitud,predicho,medido\n"
                 "a,M,1,0\nb,M,2,2\nc,M,3,3\n", encoding="utf-8")
    salida = tmp_path / "r.json"
    assert ve.main([str(f), "--json", str(salida)]) == 0
    json.loads(salida.read_text(encoding="utf-8"),
               parse_constant=lambda c: pytest.fail(f"JSON invalido: {c}"))
