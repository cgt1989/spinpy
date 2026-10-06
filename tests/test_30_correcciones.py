"""
Bloque 30 — Correcciones de los artefactos de borde (plato rigido y nucleo).

REFERENCIA
  `comparativa_motores/correcciones/`: el espinodoide de referencia ensayado
  dentro de un bloque mayor del mismo hueso (configuracion embebida), sin
  techo cargado ni caras cortadas en la region del VOI. Frente a esa
  referencia, la E aparente de la app con traccion en el techo erraba un
  -43 a -55 %; con plato rigido y medida en el nucleo, +0,8 a +1,5 %.

LO QUE SE VERIFICA
  (1) Bloque macizo con traccion y con plato: la E del nucleo
      (`fem.magnitudes_nucleo`, extensometro virtual) es E_s. La deformacion
      es uniforme, asi que es exacta en el problema discreto (hex8 y TET10).
  (2) `fem.ensayo` con 'lineal_plato' en el bloque: E_app del plato = E_s,
      cociente plato/traccion = 1 y `corregido` usa el nucleo con plato.
  (3) Con el motor de la app (no resuelve el plato): 'lineal_plato' se
      declara no disponible, no es un fallo, y `corregido` usa el nucleo con
      traccion.
  (4) VOI demasiado pequeno para el margen: el nucleo devuelve ok = False.
  (6) Por `fem.analizar` (el camino de la GUI y la CLI), con el motor de la
      app y con NGSolve: 'lineal_plato' no se cuenta como fallo y el registro
      queda ok (regresion de la V2.1.0, corregida en la V2.1.1).
  (7) VOI demasiado pequeno para el nucleo, con plato: el E corregido es el
      del plato sobre el VOI completo ('plato_voi'), no falta (V2.1.1).
  (5) Espinodoide de referencia a 32^3 (hex8): la E corregida queda a menos
      del 3 % de la referencia embebida medida en
      `comparativa_motores/correcciones/resultados/referencia.json`
      (nucleo, 678,6 MPa) y la E de la linea base, a mas del 40 %.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  (1), (2) 1e-8 relativa (resolvedor directo).
  (5) 3 % (lo medido fue +0,8 %); la linea base, error > 40 % (medido -55 %).
"""
from __future__ import annotations

import numpy as np
import pytest

from spinpy import fem, generar_mascara, motores

BLOQUE = "30 Correcciones de borde"
REF = "comparativa_motores/correcciones/"
E_MPA, NU = 20000.0, 0.30
#: E del nucleo en la configuracion embebida, 32^3 (referencia.json).
E_REF_NUCLEO_32 = 678.6
NG = "ngsolve" in motores.disponibles()


def _anotar(registro, prueba, esp, obt, err, crit, ok, nota=""):
    registro.anotar(BLOQUE, prueba, REF, esp, obt, err, crit, bool(ok),
                    nota=nota)


def _bloque(tipo):
    """Bloque macizo. TET10: el mismo del bloque 29 (8 x 8 x 12 a 0,05 mm;
    un macizo de 12^3 a 0,25 mm deja astillas planas en tetgen), con margen
    y franjas a su escala."""
    if tipo == "hex8":
        return (fem.mallar(np.ones((12, 12, 12), bool), np.full(3, 0.25),
                           tipo), {})
    return (fem.mallar(np.ones((8, 8, 12), bool), np.full(3, 0.05), tipo),
            {"margen": 0.1, "franja": 0.05})


def _resolver(malla, control, motor="ngsolve"):
    carga = (1.0,) if control == "fuerza" else (1e-3,)
    p = motores.problema_de_malla(malla, E=E_MPA, nu=NU, control=control,
                                  cargas=carga, solver="directo")
    out = motores.resolver(p, motor)
    s = fem.tension_elemental(p["nodos"], malla["elems"], out["u"], E_MPA * 1e6,
                              NU)
    return out, s


@pytest.mark.skipif(not NG, reason="NGSolve no instalado")
@pytest.mark.parametrize("tipo", ["hex8", "tet10"])
@pytest.mark.parametrize("control", ["fuerza", "plato"])
def test_nucleo_bloque_macizo(registro, tipo, control):
    malla, kw = _bloque(tipo)
    out, s = _resolver(malla, control)
    nuc = fem.magnitudes_nucleo(malla, out["u"], s, **kw)
    err = abs(nuc["E_app"] / (E_MPA * 1e6) - 1)
    ok = nuc["ok"] and err < 1e-8
    _anotar(registro, f"nucleo {tipo} {control}", E_MPA, nuc["E_app"] / 1e6,
            err, "< 1e-8", ok)
    assert ok, nuc


@pytest.mark.skipif(not NG, reason="NGSolve no instalado")
def test_ensayo_lineal_plato(registro):
    malla, _ = _bloque("hex8")
    prot = fem.protocolo("app")
    reg = fem.ensayo(malla, prot, analisis=["lineal", "lineal_plato"],
                     motor="ngsolve", aislado=False, solver="directo")
    lp, cor = reg["lineal_plato"], reg["corregido"]
    err = max(abs(lp["E_app"] / prot["E_s"] - 1),
              abs(lp["cociente_plato_fuerza"] - 1),
              abs(cor["E_app"] / prot["E_s"] - 1))
    ok = (reg["ok"] and err < 1e-8 and cor["metodo_E"] == "plato_nucleo"
          and cor["metodo_p99"] == "plato_voi")
    _anotar(registro, "ensayo con lineal_plato", 1.0,
            lp["cociente_plato_fuerza"], err, "< 1e-8", ok)
    assert ok, (reg.get("fallos"), lp, cor)


def test_app_sin_plato(registro):
    malla, _ = _bloque("hex8")
    BW = np.ones((12, 12, 12), bool)
    reg = fem.ensayo(malla, fem.protocolo("app"),
                     analisis=["lineal", "lineal_plato"], motor="app",
                     aislado=False, BW=BW)
    ok = (reg["ok"] and reg["lineal_plato"].get("no_disponible")
          and reg["corregido"]["metodo_E"] == "traccion_nucleo")
    _anotar(registro, "motor app sin plato (1 = no disponible)", 1.0,
            float(bool(reg["lineal_plato"].get("no_disponible"))), 0.0,
            "exacto", ok)
    assert ok, reg.get("fallos")


def test_nucleo_voi_pequeno(registro):
    malla = fem.mallar(np.ones((8, 8, 8), bool), np.full(3, 0.1), "hex8")
    u = np.zeros((malla["nodos"].shape[0], 3))
    s = np.zeros((malla["elems"].shape[0], 6))
    nuc = fem.magnitudes_nucleo(malla, u, s)
    ok = nuc["ok"] is False
    _anotar(registro, "VOI pequeno (0 = sin nucleo)", 0.0,
            float(bool(nuc["ok"])), 0.0, "exacto", ok)
    assert ok


@pytest.mark.skipif(not NG, reason="NGSolve no instalado")
def test_espinodoide_frente_a_embebido(registro):
    n = 32
    BW, _, _ = generar_mascara(resolution=n, wave_number=12 * np.pi,
                               num_waves=700, thetas=(30, 30, 90), rho=0.30,
                               seed=1)
    malla = fem.mallar(BW, np.full(3, 5.0 / n), "hex8")
    reg = fem.ensayo(malla, fem.protocolo("app", E_s=E_MPA * 1e6),
                     analisis=["lineal", "lineal_plato"], motor="ngsolve",
                     aislado=False)
    E_cor = reg["corregido"]["E_app"] / 1e6
    E_base = reg["lineal"]["E_app"] / 1e6
    err = abs(E_cor / E_REF_NUCLEO_32 - 1)
    err_base = abs(E_base / E_REF_NUCLEO_32 - 1)
    ok = err < 0.03 and err_base > 0.40
    _anotar(registro, "espinodoide 32^3 frente a embebido", E_REF_NUCLEO_32,
            E_cor, err, "< 3 % (linea base > 40 %)", ok,
            nota=f"linea base {E_base:.1f} MPa ({err_base:.1%})")
    assert ok, (E_cor, E_base)


@pytest.mark.parametrize("motor", ["app", "ngsolve"])
def test_analizar_lineal_plato_no_es_fallo(registro, motor):
    if motor not in motores.disponibles():
        pytest.skip(f"{motor} no instalado")
    BW, _, _ = generar_mascara(resolution=16, wave_number=12 * np.pi,
                               num_waves=700, thetas=(30, 30, 90), rho=0.30,
                               seed=1)
    reg = fem.analizar(BW, np.full(3, 5.0 / 16), fem.protocolo("app"),
                       malla="hex8", analisis=["lineal", "lineal_plato"],
                       n=16, motores_fem=[motor], aislado=False,
                       comparar_app=False)[0]
    esperado = "traccion_nucleo" if motor == "app" else "plato_nucleo"
    ok = (reg["ok"] and not reg["fallos"]
          and reg["corregido"].get("metodo_E") == esperado)
    _anotar(registro, f"analizar con lineal_plato ({motor}), 1 = ok", 1.0,
            float(ok), 0.0, "exacto", ok, nota=str(reg["fallos"]))
    assert ok, reg["fallos"]


@pytest.mark.skipif(not NG, reason="NGSolve no instalado")
def test_corregido_sin_nucleo_usa_plato(registro):
    malla = fem.mallar(np.ones((8, 8, 8), bool), np.full(3, 0.2), "hex8")
    prot = fem.protocolo("app")
    reg = fem.ensayo(malla, prot, analisis=["lineal", "lineal_plato"],
                     motor="ngsolve", aislado=False, solver="directo")
    cor = reg["corregido"]
    err = abs(cor.get("E_app", 0.0) / prot["E_s"] - 1)
    ok = (not reg["lineal"]["nucleo"]["ok"] and cor.get("metodo_E") ==
          "plato_voi" and err < 1e-8)
    _anotar(registro, "E corregido sin nucleo (plato sobre el VOI)",
            prot["E_s"] / 1e6, cor.get("E_app", np.nan) / 1e6, err, "< 1e-8",
            ok)
    assert ok, cor
