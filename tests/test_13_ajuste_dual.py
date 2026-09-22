"""
Bloque 13 — Ajuste de un dual-lattice: recuperar parametros conocidos.

QUE SE VERIFICA
  Un "VOI" que ES un dual-lattice, de parametros conocidos, y el ajuste rapido
  lanzado sobre el con OTRA semilla: la busqueda no ve la misma realizacion
  que genero el objetivo, igual que con un hueso real no hay semilla que
  acertar. Si el ajuste no recupera los parametros de algo que su propia
  familia puede producir, no hay motivo para creer lo que devuelva sobre un
  hueso.

EL OBJETIVO
  64^3, voxel 0.075 mm (lado 4.8 mm, del orden del VOI equino),
  rho 0.30, 6 celdas, estiramiento (1, 1, 1.5), semilla 20260727.

TOLERANCIAS DECLARADAS ANTES DE MEDIR, Y DE DONDE SALEN
  densidad     |d - 0.30| / 0.30 <= 0.16. El modo rapido prueba
               BV/TV x [0.85, 1.00, 1.15]; ganar el punto central es lo
               esperado y el vecino esta a 0.15.
  celdas       |c - 6| / 6 <= 0.45. La rejilla rapida es la estimacion desde
               Tb.N por [0.7, 1.0, 1.4]; con la calibracion (celdas/Tb.N 0.54
               a rho 0.30) la estimacion cae cerca de 6, y el vecino mas
               lejano esta a un factor 1.4.
  estiramiento el preset ganador alarga UN eje (no el isotropo). Cual de los
               tres depende de la rotacion de partida, que lleva z sobre la
               direccion principal del VOI, asi que se pide ademas la
               siguiente condicion y no un preset concreto.
  orientacion  direccion principal del candidato final a <= 20 grados de la
               del VOI (correccion L1 incluida).
  error        error morfometrico final <= 0.02. Referencia: los ajustes de
               spinodoide al VOI proximal de H4 dan 0.003-0.015.
  estructura   el resultado trae las mismas claves que el del spinodoide
               (`parametros`, `metricas_spin`, `diagnostico`, `traza`,
               `alineacion`), que es lo que la interfaz lee.
"""
import numpy as np

from spinpy.dual_lattice import generar_dual_lattice
from spinpy.metodos import REGISTRO, REGISTROS
from spinpy.morphometry import morfometria

BLOQUE = "13 Ajuste dual-lattice (recuperacion)"
REF = "construccion: objetivo de la propia familia"


def _angulo(u, v):
    u = np.asarray(u, float) / np.linalg.norm(u)
    v = np.asarray(v, float) / np.linalg.norm(v)
    return float(np.degrees(np.arccos(min(1.0, abs(float(u @ v))))))


def test_recupera_parametros(registro):
    VOI, _, _ = generar_dual_lattice(64, 6.0, 0.30, estiramiento=(1, 1, 1.5),
                                     seed=20260727)
    sp = np.full(3, 0.075)
    m_voi = morfometria(VOI, sp)

    r = REGISTROS["dual-lattice"]["rapido"]["correr"](
        VOI, sp, m_voi=m_voi, num_waves=700, esquema="rechazo")
    p = r["parametros"]

    for clave in ("parametros", "metricas_spin", "diagnostico", "traza",
                  "alineacion", "error"):
        assert clave in r, clave

    d_rel = abs(p["densidad"] - 0.30) / 0.30
    ok_d = d_rel <= 0.16
    registro.anotar(BLOQUE, "densidad recuperada", REF, 0.30, p["densidad"],
                    "rel <= 0.16", "rejilla rapida x[0.85,1,1.15]", ok_d)

    c_rel = abs(p["celdas"] - 6.0) / 6.0
    ok_c = c_rel <= 0.45
    registro.anotar(BLOQUE, "celdas recuperadas", REF, 6.0, p["celdas"],
                    "rel <= 0.45", "estimacion desde Tb.N x[0.7,1,1.4]", ok_c,
                    nota=f"estimadas {r['diagnostico']['celdas_estimadas']:.2f}")

    e = np.asarray(p["estiramiento"], float)
    ok_e = bool(e.max() > 1.0 + 1e-9)
    registro.anotar(BLOQUE, "preset ganador alarga un eje", REF, 1.5,
                    float(e.max()), "> 1", "no gana el isotropo", ok_e,
                    nota=str(p["estiramiento"]))

    ang = _angulo(r["metricas_spin"]["dir_principal"], m_voi["dir_principal"])
    ok_a = ang <= 20.0
    registro.anotar(BLOQUE, "angulo candidato-VOI", REF, 0.0, ang,
                    "<= 20 grados", "L1", ok_a)

    ok_err = r["error"] <= 0.02
    registro.anotar(BLOQUE, "error morfometrico final", REF, None, r["error"],
                    "<= 0.02", "spinodoide en H4: 0.003-0.015", ok_err)

    assert ok_d and ok_c and ok_e and ok_a and ok_err, {
        k: v for k, v in r.items()
        if k in ("parametros", "error", "diagnostico", "alineacion")}


def test_registro_spinodoide_intacto(registro):
    """Anadir la familia no puede haber cambiado los metodos del spinodoide."""
    ok = (list(REGISTRO) == ["rapido", "completo", "mecanico", "equitativo"]
          and REGISTROS["spinodoide"] is REGISTRO
          and list(REGISTROS["dual-lattice"]) == ["rapido", "completo"])
    registro.anotar(BLOQUE, "registros por familia", REF, None, None,
                    "spinodoide 4 metodos, dual-lattice 2", "sin regresion", ok)
    assert ok
