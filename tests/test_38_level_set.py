"""
Bloque 38: la tercera familia, level-set foam (morphoGen).

Se salta entera si morphoGen no esta instalado: es una dependencia OPCIONAL
(GPL-3.0-or-later, Python >= 3.12) y el nucleo de spinpy es MIT y funciona sin
ella. Lo que SI se comprueba siempre, con o sin morphoGen, es que el nucleo no
la importa al cargar y que la fuente unica de familias es coherente.

Tolerancias declaradas antes de medir:
  * densidad: morphoGen resuelve el espesor para la fraccion pedida sobre su
    rejilla de generacion; en la del candidato (res^3) la diferencia es de
    discretizacion. Se exige < 2 % relativo a 48^3.
  * reproducibilidad: misma semilla, misma mascara, bit a bit.
  * Conn.D: `metricas_forma(conectividad=True)` y `conectividad()` deben dar
    EXACTAMENTE el mismo numero (es la misma funcion).
"""

import subprocess
import sys

import numpy as np
import pytest

from spinpy import familias, procedencia
from spinpy.error import (PESOS_LEVEL_SET, TERMINOS_OPCIONALES, medidas_necesarias,
                          validar_pesos)
from spinpy.morphometry import conectividad, metricas_forma
from spinpy.level_set import disponible, es_viable

# ---------------------------------------------------------------------------
# Siempre: no necesitan morphoGen
# ---------------------------------------------------------------------------


def test_nucleo_no_importa_morphogen_al_cargar():
    """`import spinpy` no puede arrastrar una libreria GPL ni exigir Python 3.12."""
    codigo = ("import sys, spinpy; "
              "assert 'morphogen' not in sys.modules, 'spinpy importa morphogen'")
    r = subprocess.run([sys.executable, "-c", codigo], capture_output=True,
                       text=True)
    assert r.returncode == 0, r.stderr


def test_fuente_unica_de_familias():
    assert familias.FAMILIAS == ("spinodoide", "dual-lattice", "level-set-foam")
    assert familias.MODOS_FAMILIA == familias.FAMILIAS + ("todas",)
    # el spinodoide conserva sus claves de siempre; las otras llevan sufijo
    assert familias.clave("ajuste", "spinodoide") == "ajuste"
    assert familias.clave("ajuste", "dual-lattice") == "ajuste_dual"
    assert familias.clave("ajuste", "level-set-foam") == "ajuste_ls"
    sufijos = [familias.sufijo(f) for f in familias.FAMILIAS]
    assert len(set(sufijos)) == len(sufijos), "dos familias con el mismo sufijo"
    for f in familias.FAMILIAS:
        assert familias.clave_morfometria(f).startswith("morfometria_")


def test_es_viable_regla_de_morphogen():
    """Con pocas semillas muy alargadas el empaquetado denso no cabe."""
    assert es_viable(134, (1.0, 1.5, 1.0))
    assert es_viable(800, (1.0, 3.0, 4.35))
    assert not es_viable(50, (1.0, 3.0, 4.35))
    assert not es_viable(100, (1.0, -1.0, 1.0))      # estiramiento no positivo


def test_termino_connd_registrado_y_pesos_level_set():
    assert "ConnD" in TERMINOS_OPCIONALES
    assert TERMINOS_OPCIONALES["ConnD"]["medida"] == "conectividad"
    w = validar_pesos(PESOS_LEVEL_SET)
    # los tres discriminantes de spinpy, mas Conn.D, y TbSp apagado
    assert w["TbSp"] == 0.0
    assert all(w[k] > 0 for k in ("PoDm", "TbTh_CV", "H_desv", "ConnD"))
    assert medidas_necesarias(w)["forma"] == {"poro", "tbth", "curvatura",
                                              "conectividad"}


def test_metricas_forma_conectividad_es_la_funcion_de_conectividad():
    rng = np.random.default_rng(3)
    BW = rng.random((24, 24, 24)) < 0.35
    sp = np.full(3, 0.1)
    f = metricas_forma(BW, sp, conectividad=True)
    assert f["ConnD"] == conectividad(BW, sp)["ConnD"]
    # la compuerta apagada no calcula nada nuevo
    assert "ConnD" not in metricas_forma(BW, sp)


def test_procedencia_sin_onda_ni_esquema():
    b = procedencia.bloque("level-set-foam", semilla=3, morphogen="1.1.0")
    assert b["familia"] == "level-set-foam"
    assert b["esquema"] == "no_aplica"
    assert b["wave_number_pi"] is None and b["wave_number_rad"] is None
    assert b["morphogen"] == "1.1.0"


def test_kwargs_generador_level_set_sin_R():
    par = {"familia": "level-set-foam", "resolucion": 32, "n_celdas": 120.0,
           "densidad": 0.3, "estiramiento": [1.0, 1.5, 1.0], "curvatura": -0.5,
           "c1": 10.0, "c2": 0.5, "semilla": 7}
    kw = procedencia.kwargs_generador(par)
    assert "R" not in kw, "level-set foam no admite rotacion"
    assert kw["seed"] == 7 and kw["n_celdas"] == 120.0 and kw["rho"] == 0.3
    assert procedencia.kwargs_generador(par, semilla=9)["seed"] == 9


def test_figuras_tienen_color_orden_y_nombre_de_cada_familia():
    from spinpy import figuras
    for f in familias.FAMILIAS:
        assert f in figuras.COLOR and f in figuras.ORDEN and f in figuras.NOMBRES
    colores = [figuras.COLOR[e] for e in figuras.ORDEN]
    assert len(set(colores)) == len(colores), "dos estructuras con el mismo color"


# ---------------------------------------------------------------------------
# Con morphoGen
# ---------------------------------------------------------------------------

pytestmark_ls = pytest.mark.skipif(
    not disponible(), reason="morphoGen no instalado (extra spinpy[level-set])")


@pytestmark_ls
def test_generador_contrato_densidad_y_reproducibilidad():
    from spinpy.level_set import generar_level_set
    BW, campo, info = generar_level_set(48, 100, 0.30, (1.0, 1.5, 1.0),
                                        curvatura=-0.5, c1=10.0, c2=0.5, seed=5)
    assert BW.shape == (48, 48, 48) and BW.dtype == bool
    assert campo.shape == BW.shape
    assert np.array_equal(BW, campo <= 0.0)
    assert abs(BW.mean() - 0.30) / 0.30 < 0.02, BW.mean()
    assert info["familia"] == "level-set-foam"
    assert info["rho_objetivo"] == 0.30
    assert info["esquema"] == "no_aplica"
    BW2, _, _ = generar_level_set(48, 100, 0.30, (1.0, 1.5, 1.0),
                                  curvatura=-0.5, c1=10.0, c2=0.5, seed=5)
    assert np.array_equal(BW, BW2), "misma semilla, distinta mascara"
    BW3, _, _ = generar_level_set(48, 100, 0.30, (1.0, 1.5, 1.0),
                                  curvatura=-0.5, c1=10.0, c2=0.5, seed=6)
    assert not np.array_equal(BW, BW3), "la semilla no cambia nada"


@pytestmark_ls
def test_generador_rechaza_rotacion_e_inviables():
    from spinpy.level_set import generar_level_set
    giro = np.eye(3)[[1, 0, 2]]
    with pytest.raises(ValueError, match="rotacion"):
        generar_level_set(24, 100, 0.3, R=giro)
    with pytest.raises(ValueError, match="no cabe"):
        generar_level_set(24, 50, 0.3, (1.0, 3.0, 4.35))
    with pytest.raises(ValueError):
        generar_level_set(24, 100, 1.2)                # rho fuera de (0, 1)
    with pytest.raises(ValueError):
        generar_level_set(24, 100, 0.3, curvatura=1.0)  # K_c < 1


@pytestmark_ls
def test_el_estiramiento_alarga_el_eje_pedido():
    """El DA crece con el estiramiento: es lo que el ajuste usa para orientar."""
    from spinpy.level_set import generar_level_set
    from spinpy.morphometry import morfometria
    sp = np.full(3, 3.0 / 48)
    da = {}
    assert es_viable(200, (1.0, 1.0, 2.0))
    for nombre, e in (("iso", (1.0, 1.0, 1.0)), ("z", (1.0, 1.0, 2.0))):
        BW, _, _ = generar_level_set(48, 200, 0.35, e, seed=2)
        da[nombre] = morfometria(BW, sp)["DA"]
    assert da["z"] > da["iso"] + 0.1, da


@pytestmark_ls
def test_desde_parametros_coincide_con_el_generador():
    from spinpy.level_set import generar_level_set
    from spinpy.lote import generar_desde_parametros
    par = {"familia": "level-set-foam", "resolucion": 32, "n_celdas": 90.0,
           "densidad": 0.3, "estiramiento": [1.0, 1.0, 1.5], "curvatura": 0.0,
           "c1": 0.0, "c2": 0.0, "semilla": 4}
    BW = generar_desde_parametros(par, 4)
    ref, _, _ = generar_level_set(32, 90.0, 0.3, (1.0, 1.0, 1.5), seed=4)
    assert np.array_equal(BW, ref)


@pytestmark_ls
def test_registro_de_metodos_y_replicas():
    from spinpy import metodos
    assert set(metodos.REGISTROS) == set(familias.FAMILIAS)
    reg = metodos.REGISTROS["level-set-foam"]
    assert {"rapido", "completo"} <= set(reg)
    assert reg["completo"]["n_eval"] == 60 and reg["rapido"]["n_eval"] == 15
    assert metodos.gobierna_de("level-set-foam") is metodos.GOBIERNA_LS
    par = {"familia": "level-set-foam", "resolucion": 24, "n_celdas": 80.0,
           "densidad": 0.3, "estiramiento": [1.0, 1.0, 1.0], "curvatura": 0.0,
           "c1": 0.0, "c2": 0.0, "semilla": 11, "R": np.eye(3).tolist()}
    rep = metodos.generar_replicas(par, np.full(3, 0.1), n=2, medir=False)
    assert rep["familia"] == "level-set-foam" and len(rep["filas"]) == 2
    assert [f["semilla"] for f in rep["filas"]] == [11, 12]
    assert all(abs(f["rho_obtenida"] - 0.3) < 0.03 for f in rep["filas"])


@pytestmark_ls
def test_paquete_de_reproduccion_regenera_la_estructura():
    """El informe regenera el level-set foam por la unica traduccion."""
    from spinpy import informe
    par = {"familia": "level-set-foam", "resolucion": 24, "n_celdas": 80.0,
           "densidad": 0.3, "estiramiento": [1.0, 1.0, 1.0], "curvatura": 0.0,
           "c1": 0.0, "c2": 0.0, "semilla": 4}
    doc = {"voi": {"nombre": "x", "forma": [24, 24, 24],
                   "spacing_mm": [0.1, 0.1, 0.1]},
           "resultados": {"ajuste_ls": {"parametros": par,
                                        "procedencia": procedencia.desde_parametros(par)}}}
    paq = informe.paquete_reproduccion(doc)
    a = paq["ajustes"]["level-set-foam"]
    assert a["forma"] == [24, 24, 24]
    assert abs(a["BVTV"] - 0.3) < 0.03 and a["sha256_mascara"]
    assert "level-set-foam" in informe.ESTRUCTURAS


@pytestmark_ls
@pytest.mark.lento
def test_ajuste_rapido_contrato_y_densidad_fija():
    """El ajuste devuelve el contrato de las otras familias y NO busca la densidad."""
    from spinpy.fit_levelset import ajustar_level_set
    from spinpy.level_set import generar_level_set
    # VOI sintetico: un level-set foam alargado en z, a 48^3 de 2.4 mm de lado
    VOI, _, _ = generar_level_set(48, 110, 0.33, (1.0, 1.0, 1.8), seed=21)
    sp = np.full(3, 2.4 / 48)
    r = ajustar_level_set(VOI, sp, modo="rapido", replicas=0, seed=3)
    for clave in ("familia", "procedencia", "parametros", "objetivo",
                  "diagnostico", "alineacion", "error", "incertidumbre",
                  "metricas_voi", "metricas_spin", "mascara", "traza",
                  "n_evaluaciones", "tiempo_s", "k2", "mecanico", "rejillas"):
        assert clave in r, clave
    p = r["parametros"]
    assert p["familia"] == "level-set-foam"
    assert p["densidad"] == pytest.approx(VOI.mean(), rel=1e-9)   # fijada, no buscada
    assert p["R"] == np.eye(3).tolist()
    assert np.isfinite(r["error"]) and r["error"] < 0.2
    assert r["n_evaluaciones"] == r["rejillas"]["n_eval_previstas"] == 15
    assert abs(r["mascara"].mean() - VOI.mean()) / VOI.mean() < 0.03
    # el objetivo es el ampliado y declara sus medidas de forma
    assert set(r["objetivo"]["medidas_por_evaluacion"]) == {
        "poro", "tbth", "curvatura", "conectividad"}
    # lo que devuelve se regenera igual desde sus parametros
    from spinpy.lote import generar_desde_parametros
    BW = generar_desde_parametros(p, p["semilla"])
    assert np.array_equal(BW, r["mascara"])
    # la alineacion declara su metodo
    assert r["alineacion"]["metodo"] == "permutacion de ejes"
