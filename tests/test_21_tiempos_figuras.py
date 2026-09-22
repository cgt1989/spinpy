"""
Bloque 21 — Estimacion de tiempos del informe automatico y opciones de figuras.

QUE SE VERIFICA
  puntos        el modelo (`spinpy.tiempos`) reproduce los tiempos MEDIDOS con
                los que se calibro, dentro de ±20 % cada uno. Si alguien toca
                un exponente o un coeficiente sin volver a medir, falla aqui.
  validacion    el total de la ejecucion completa del cerdo (1755 s medidos)
                queda dentro de ±10 %, y el ajuste dual de 32^3 (71 s), que no
                se uso para calibrar, dentro de ±20 %.
  coherencia    marcar una etapa o subir una resolucion nunca baja el total;
                la convergencia con menos de tres mallas cuesta 0 (no se
                lanza); el factor de equipo multiplica, y su actualizacion
                queda acotada a [0.2, 5].
  figuras       `opciones_3d` descarta estilos y vistas desconocidos, rellena
                lo que falta y mete siempre el estilo y la vista principales
                en la galeria.

Tolerancias declaradas antes de medir. Los tiempos de referencia son del
equipo de desarrollo (2026-09-16); el test no cronometra nada, compara el
modelo con esas cifras, asi que no depende de la maquina que lo corre.
"""
import numpy as np
import pytest

from spinpy import figuras as F
from spinpy import tiempos as T

BLOQUE = "21 Tiempos y figuras"
REF = "tiempos.py (calibracion 2026-09-16, VOI_V1 porcino)"

# (descripcion, funcion, medido en s)
PUNTOS = [
    ("generar spinodoide 48^3", lambda: T.gen("spinodoide", 48), 1.25),
    ("generar spinodoide 96^3", lambda: T.gen("spinodoide", 96), 8.75),
    ("generar dual-lattice 48^3", lambda: T.gen("dual-lattice", 48), 0.60),
    ("generar dual-lattice 96^3", lambda: T.gen("dual-lattice", 96), 1.30),
    ("morfometria 96^3", lambda: T.morfo(96 ** 3), 1.90),
    ("morfometria 188^3", lambda: T.morfo(188 ** 3), 4.46),
    ("morfometria + extra 188^3", lambda: T.morfo(188 ** 3, extra=True), 8.45),
    ("morfometria + Po.Dm 188^3",
     lambda: T.morfo(188 ** 3, extra=True, poro=True), 58.6),
    ("EF 32^3", lambda: T.morfo(32 ** 3, ef=True) - T.morfo(32 ** 3), 184.0),
    ("homogeneizacion 16^3", lambda: T.homog(16), 6.4),
    ("homogeneizacion 24^3", lambda: T.homog(24), 31.6),
    ("ensayo un eje 24^3", lambda: T.ensayo(24), 2.6),
    ("ensayo un eje 32^3", lambda: T.ensayo(32), 8.0),
    ("ensayo un eje 48^3", lambda: T.ensayo(48), 26.6),
    ("distribuciones 188^3", lambda: T.distribuciones(188 ** 3), 73.0),
    ("distribuciones 96^3", lambda: T.distribuciones(96 ** 3), 4.5),
]


@pytest.mark.parametrize("nombre,fn,medido", PUNTOS,
                         ids=[p[0] for p in PUNTOS])
def test_puntos_de_calibracion(registro, nombre, fn, medido):
    obtenido = fn()
    ok = abs(obtenido / medido - 1.0) <= 0.20
    registro.anotar(BLOQUE, f"calibracion: {nombre}", REF, medido, obtenido,
                    "±20 %", "modelo frente a tiempo medido", ok)
    assert ok


def _plan_cerdo(**cambios):
    plan = dict(familias=["spinodoide"], vox_voi=188 ** 3, n_fit=96,
                ajuste=True, morfometria=True, extra=True, poro=True,
                ef=False, n_cand=96, elastico=True, n_homog=32, ensayo=True,
                n_fe=40, n_ejes=1, comparado=True, n_estilos=1, n_vistas=1,
                suavizar=False, K=8, n_resm=64, pasos=8,
                protocolo="adelgazamiento", n_sim=188)
    plan.update(cambios)
    return plan


def test_total_ejecucion_cerdo(registro):
    obtenido = T.total(T.estimar(_plan_cerdo()))
    ok = abs(obtenido / 1755.0 - 1.0) <= 0.10
    registro.anotar(BLOQUE, "total de la ejecucion completa (cerdo)", REF,
                    1755.0, obtenido, "±10 %",
                    "informe automatico por defecto sobre VOI_V1", ok)
    assert ok


def test_ajuste_dual_no_calibrado(registro):
    obtenido = T.ajuste("dual-lattice", 32, 32 ** 3)
    ok = abs(obtenido / 71.3 - 1.0) <= 0.20
    registro.anotar(BLOQUE, "ajuste dual-lattice 32^3 (no calibrado)", REF,
                    71.3, obtenido, "±20 %",
                    "coeficiente derivado de contar evaluaciones", ok)
    assert ok


def test_coherencia(registro):
    base = T.total(T.estimar(_plan_cerdo()))
    mas = [T.total(T.estimar(_plan_cerdo(**c))) for c in (
        {"ef": True}, {"n_ejes": 3}, {"n_fe": 48}, {"dispersion": True},
        {"perdida": True}, {"fallo": True}, {"familias": ["spinodoide",
                                                          "dual-lattice"]},
        {"n_estilos": 3, "n_vistas": 7})]
    menos = T.total(T.estimar(_plan_cerdo(ajuste=False, elastico=False)))
    conv_baja = T.estimar(_plan_cerdo(convergencia=True, n_fe=16))
    conv_ok = T.estimar(_plan_cerdo(convergencia=True, n_fe=40))
    doble = T.total(T.estimar(_plan_cerdo(),
                              {k: 2.0 for k in T.estimar(_plan_cerdo())}))
    acot = (T.actualizar_factor(None, 100.0, 1.0),
            T.actualizar_factor(None, 1.0, 100.0),
            T.actualizar_factor(1.0, 30.0, 10.0))
    ok = (all(m > base for m in mas) and menos < base
          and conv_baja["convergencia"][None] == 0.0
          and conv_ok["convergencia"][None] > 0.0
          and abs(doble / base - 2.0) < 1e-9
          and acot[0] == 5.0 and acot[1] == 0.2 and acot[2] == 2.0)
    registro.anotar(BLOQUE, "coherencia de la estimacion", REF, None, None,
                    "exacto", "mas trabajo nunca cuesta menos; factor "
                    "multiplica; actualizacion acotada", ok,
                    nota=f"base {base:.0f} s, EF {mas[0]:.0f} s")
    assert ok


def test_opciones_3d(registro):
    o = F.opciones_3d({"estilos": ["hueso", "inventado"],
                       "vistas": ["lateral", "diagonal"],
                       "estilo_principal": "altura",
                       "vista_principal": "superior", "otra": 1})
    vacio = F.opciones_3d(None)
    malo = F.opciones_3d({"estilo_principal": "x", "vista_principal": "y"})
    ok = (o["estilos"] == ["altura", "hueso"]
          and o["vistas"] == ["superior", "lateral"]
          and "otra" not in o
          and vacio == F.opciones_3d(F.OPCIONES_3D)
          and malo["estilo_principal"] == F.OPCIONES_3D["estilo_principal"]
          and malo["vista_principal"] == F.OPCIONES_3D["vista_principal"]
          and all(e in F.ESTILOS_3D for e in vacio["estilos"])
          and all(v in F.VISTAS_3D for v in vacio["vistas"]))
    registro.anotar(BLOQUE, "validacion de opciones de figuras", REF, None,
                    None, "exacto", "desconocidos fuera, principales dentro",
                    ok, nota=str(o))
    assert ok


def _C_isotropo(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    C = np.zeros((6, 6))
    C[:3, :3] = lam
    C[np.arange(3), np.arange(3)] = lam + 2 * mu
    C[np.arange(3, 6), np.arange(3, 6)] = mu     # cizalla ingenieril
    return C


def test_modulo_direccional(registro):
    """E(n) de la figura 6 contra soluciones cerradas.

    Isotropo: E(n) = E en toda direccion. Ortotropo: sobre los ejes, E = 1/S_ii;
    y la convencion de cizalla ingenieril se comprueba en la diagonal a 45° del
    plano XY, donde 1/E = (S11 + S22 + 2 S12 + S66)/4.
    """
    rng = np.random.default_rng(3)
    n = rng.normal(size=(200, 3))
    n /= np.linalg.norm(n, axis=1)[:, None]
    iso = F.modulo_direccional(_C_isotropo(18e3, 0.3), n)
    S = np.diag([1 / 100., 1 / 200., 1 / 400., 1 / 30., 1 / 40., 1 / 50.])
    S[0, 1] = S[1, 0] = -0.3 / 200.
    S[0, 2] = S[2, 0] = -0.2 / 400.
    S[1, 2] = S[2, 1] = -0.25 / 400.
    C = np.linalg.inv(S)
    ejes = F.modulo_direccional(C, np.eye(3))
    d = np.array([[1, 1, 0]]) / np.sqrt(2)
    diag = F.modulo_direccional(C, d)[0]
    diag_ref = 4.0 / (S[0, 0] + S[1, 1] + 2 * S[0, 1] + S[5, 5])
    ok = (np.allclose(iso, 18e3, rtol=1e-9)
          and np.allclose(ejes, [100., 200., 400.], rtol=1e-9)
          and abs(diag / diag_ref - 1) < 1e-9)
    registro.anotar(BLOQUE, "modulo de Young direccional E(n)", REF,
                    diag_ref, diag, "1e-9", "isotropo, ejes y diagonal a 45°",
                    ok)
    assert ok
