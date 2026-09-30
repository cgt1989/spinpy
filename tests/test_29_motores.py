"""
Bloque 29 — Motores FEM internos (spinpy.motores) contra soluciones cerradas.

REFERENCIA
  `comparativa_motores/INFORME.md`: la comparativa de la app, NGSolve,
  FEniCSx, scikit-fem y SfePy sobre las mismas mallas. Esta prueba es su
  version pequena y repetible, sobre los motores INSTALADOS en el entorno
  (los que falten se saltan y se anotan).

LO QUE SE VERIFICA
  (1) Bloque macizo, lineal, con cada motor y cada malla que resuelve (hex8,
      TET10): E_app = E_s y la fuerza de equilibrio igual a la aplicada. Es
      solucion exacta del problema discreto: deformacion uniforme.
  (2) Plato rigido lineal en el bloque, hex8 y TET10: la reaccion del plato
      es E_s A eps (exacta).
  (3) No lineal, bloque en compresion uniaxial hasta el 20 %, con plato y
      laterales libres: la fuerza de reaccion por paso frente a la solucion
      cerrada (`uniaxial`) de St. Venant-Kirchhoff y del neo-Hookeano de
      FEBio. La deformacion es homogenea, asi que el problema discreto la
      reproduce exactamente.
  (4) Entre motores: el mismo espinodoide 16^3 (hex8, lineal) con todos los
      motores instalados da el mismo campo de desplazamientos que la app.
  (5) Combinaciones no disponibles: la app con TET10 o no lineal, y SfePy
      con el neo-Hookeano, lanzan `NoDisponible` (nunca un resultado de otro
      modelo).
  (6) Proceso hijo (`resolver_aislado`, el de la GUI): da lo mismo que en el
      proceso actual, y cancelar lo detiene con `Cancelado`.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  (1), (2) 1e-9 relativa (resolvedor directo; iterativo a 1e-10).
  (3) 1e-8 relativa a la fuerza cerrada (Newton a 1e-10 del residuo).
  (4) 1e-6 del maximo de |u| (la app resuelve con CG a 1e-8).
  (5), (6) exacto.
"""
from __future__ import annotations

import threading
import time

import numpy as np
import pytest

from spinpy import fem, motores
from spinpy.motores import NoDisponible

BLOQUE = "29 Motores FEM"
REF = "comparativa_motores/INFORME.md"
DISP = motores.disponibles()
E_MPA, NU = 20000.0, 0.30


def _anotar(registro, prueba, esp, obt, err, crit, ok, nota=""):
    registro.anotar(BLOQUE, prueba, REF, esp, obt, err, crit, bool(ok),
                    nota=nota)


def _uniaxial(material, lam_z, E=E_MPA, nu=NU):
    """Tension nominal P_zz (MPa) con F = diag(a, a, lam_z) y S_xx = 0."""
    from scipy.optimize import brentq
    lmb = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    if material == "svk":
        return lam_z * E * 0.5 * (lam_z ** 2 - 1)
    a = brentq(lambda a: mu * (a * a - 1) + lmb * np.log(a * a * lam_z),
               0.5, 2.0, xtol=1e-15, rtol=1e-15)
    J = a * a * lam_z
    return lam_z * (mu * (1 - lam_z ** -2) + lmb * np.log(J) / lam_z ** 2)


def _malla_bloque(tipo, forma=(4, 4, 6), h=0.05):
    return fem.mallar(np.ones(forma, bool), np.full(3, h), tipo)


def _problema(m, motor, **kw):
    BW = np.ones(m["forma"], bool) if motor == "app" else None
    return motores.problema_de_malla(m, E=E_MPA, nu=NU, solver="directo",
                                     BW=BW, **kw)


def _casos_lineal():
    for mot in motores.MOTORES:
        for tipo in ("hex8", "tet10"):
            yield pytest.param(mot, tipo, id=f"{mot}-{tipo}")


@pytest.mark.parametrize("motor,tipo", list(_casos_lineal()))
def test_bloque_macizo(registro, motor, tipo):
    if motor not in DISP:
        pytest.skip(f"{motor} no instalado")
    m = _malla_bloque(tipo)
    if not motores.puede(motor, tipo):
        with pytest.raises(NoDisponible):
            motores.resolver(_problema(m, motor), motor)
        return
    out = motores.resolver(_problema(m, motor, cargas=[1.0]), motor)
    u = out["u"]
    s = fem.tension_elemental(m["nodos"], m["elems"], u, 20e9, NU)
    E = 1e6 / (abs(fem._uz_medio(m, u)) / m["H"])
    dE = abs(E / 20e9 - 1)
    dF = abs(fem._fuerza(m, s) / (1e-6 * 1e6 * m["A_bruta"]) - 1)
    dR = abs(out["F_reac"][-1] / m["A_bruta"] - 1)       # 1 MPa * A (N)
    ok = dE < 1e-9 and dF < 1e-9 and dR < 1e-9
    _anotar(registro, f"bloque macizo {tipo}, {motor}: E_app = E_s, "
            "equilibrio", 0.0, max(dE, dF, dR), max(dE, dF, dR), "1e-9", ok)
    assert ok


@pytest.mark.parametrize("motor,tipo", [
    pytest.param(m, t, id=f"{m}-{t}") for m in motores.MOTORES if m != "app"
    for t in ("hex8", "tet10")])
def test_plato_lineal(registro, motor, tipo):
    """Con TET10 detecta el espacio jerarquico de NGSolve: si el plato se
    impusiera tambien en los GDL de burbuja, o la reaccion los sumara, la
    fuerza saldria mal (medido: 0,875 de la exacta)."""
    if motor not in DISP:
        pytest.skip(f"{motor} no instalado")
    m = _malla_bloque(tipo)
    eps = 1e-3
    out = motores.resolver(_problema(m, motor, control="plato",
                                     cargas=[eps]), motor)
    F = out["F_reac"][-1]
    exacta = E_MPA * m["A_bruta"] * eps
    err = abs(F / exacta - 1)
    _anotar(registro, f"plato lineal {tipo}, {motor}: F = E A eps", exacta, F,
            err, "1e-9", err < 1e-9)
    assert err < 1e-9


def _casos_nl():
    for mot in motores.MOTORES:
        for mat in ("svk", "neohookeano"):
            yield pytest.param(mot, mat, id=f"{mot}-{mat}")


@pytest.mark.parametrize("motor,material", list(_casos_nl()))
def test_uniaxial_no_lineal(registro, motor, material):
    if motor not in DISP:
        pytest.skip(f"{motor} no instalado")
    m = _malla_bloque("hex8", forma=(2, 2, 3))
    eps = [0.05, 0.10, 0.15, 0.20]
    p = _problema(m, motor, analisis="nl", control="plato", cargas=eps,
                  material=material)
    if not motores.puede(motor, "hex8", "nl", material, "plato"):
        with pytest.raises(NoDisponible):
            motores.resolver(p, motor)
        return
    out = motores.resolver(p, motor)
    exacta = np.array([-_uniaxial(material, 1 - e) * m["A_bruta"]
                       for e in eps])
    err = float(np.abs(out["F_reac"] / exacta - 1).max())
    _anotar(registro, f"uniaxial no lineal {material}, {motor}, hasta 20 %",
            0.0, err, err, "1e-8", err < 1e-8,
            nota=f"iteraciones {out['meta'].get('iteraciones_newton')}")
    assert err < 1e-8


def test_entre_motores(registro):
    from spinpy import generar_mascara
    BW, _, _ = generar_mascara(resolution=16, wave_number=8 * np.pi,
                               num_waves=300, thetas=(90, 90, 90), rho=0.35,
                               seed=4)
    regs = fem.analizar(BW, np.full(3, 5.0 / 16), fem.protocolo("app"),
                        malla="hex8", n=16, motores_fem=DISP, aislado=False,
                        comparar_app=False)
    ref = next(r for r in regs if r["motor"]["clave"] == "app")
    u0 = ref["_campos_lineal"]["u"]
    peor = 0.0
    for r in regs:
        u = r["_campos_lineal"]["u"]
        peor = max(peor, float(np.abs(u - u0).max() / np.abs(u0).max()))
    filas = fem.tabla_motores(regs, referencia="app")
    ok = peor < 1e-6 and all(abs(f.get("dE_rel", 0.0)) < 1e-6 for f in filas)
    _anotar(registro, "espinodoide 16^3: todos los motores = app", 0.0, peor,
            peor, "1e-6 de max|u|", ok, nota=", ".join(DISP))
    assert ok


def test_no_disponible(registro):
    m = _malla_bloque("tet10", forma=(2, 2, 3))
    casos = [("app", _problema(m, "app"))]
    mh = _malla_bloque("hex8", forma=(2, 2, 3))
    casos.append(("app", _problema(mh, "app", analisis="nl",
                                   control="plato", cargas=[0.01])))
    if "sfepy" in DISP:
        casos.append(("sfepy", _problema(mh, "sfepy", analisis="nl",
                                         control="plato", cargas=[0.01],
                                         material="neohookeano")))
    ok = True
    for mot, p in casos:
        try:
            motores.resolver(p, mot)
            ok = False
        except NoDisponible:
            pass
    _anotar(registro, "combinaciones sin soporte: NoDisponible", 1, int(ok),
            0.0, "exacto", ok)
    assert ok


def test_proceso_hijo_y_cancelar(registro):
    mot = motores.RECOMENDADO if motores.RECOMENDADO in DISP else "app"
    m = _malla_bloque("hex8")
    p = _problema(m, mot, cargas=[1.0])
    a = motores.resolver(p, mot)["u"]
    b = motores.resolver_aislado(p, mot)["u"]
    igual = float(np.abs(a - b).max()) <= 1e-12 * float(np.abs(a).max())
    ev = threading.Event()
    ev.set()
    t0 = time.perf_counter()
    try:
        motores.resolver_aislado(p, mot, cancelar=ev)
        cancelado = False
    except motores.Cancelado:
        cancelado = True
    ok = igual and cancelado and time.perf_counter() - t0 < 30
    _anotar(registro, f"proceso hijo = proceso actual; cancelar detiene "
            f"({mot})", 1, int(ok), 0.0, "exacto", ok)
    assert ok
