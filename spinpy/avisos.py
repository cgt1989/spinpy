"""
avisos.py — Condiciones en las que un resultado deja de ser interpretable.

Dos comprobaciones que el estudio de familias (`Estudio_Familias/INFORME.md`,
apartados 4.5 y 4.6) mostro necesarias y que ningun calculo detecta por si
solo: los dos casos producen numeros perfectamente finitos, con buen error de
ajuste y residuo bajo, y aun asi no significan lo que parecen.

1. EL VOI NO ES TRABECULAR (BV/TV >= 0.83)
------------------------------------------
Los cuatro VOIs equinos con BV/TV de 0.83 o mas —H1 medio 0.89, H1 distal
0.92, H2 distal 0.83, H3 distal 0.87— tienen la fase poro fragmentada, y
ninguna de las dos familias ajustadas a ellos resulta bicontinua. A esa
densidad no queda trabecula que imitar: es hueso compacto con poros aislados.
El ajuste converge y da un error pequeno, que se lee como "buena replica de
hueso trabecular" cuando no lo es.

El umbral es el del banco medido, no uno teorico. `fit.RHO_MAX` explica por
que el ajuste no se acota ahi (existen especimenes asi y hay que poder
ajustarlos); este aviso es lo complementario: dejarlos ajustar, pero marcarlos.

2. EL CANDIDATO ESTA MAL ORIENTADO PARA UN ENSAYO MECANICO (> 30 grados)
-----------------------------------------------------------------------
El error morfometrico usa el DA como ESCALAR, asi que no ve la orientacion. En
10 de los 12 VOIs equinos el espinodoide ganador quedo con su eje principal a
mas de 30 grados del eje del hueso, y en siete casi perpendicular. Para la
morfometria da igual; para un ensayo no, porque se carga el candidato en una
direccion que no es la del hueso: en el VOI proximal de H4, un eje a 86 grados
dio un modulo aparente 57 veces menor que el del VOI (ver `fit.alinear_por_fabrica`).

El angulo solo tiene sentido si las dos estructuras TIENEN eje. Por debajo de
DA 1.10 —el mismo `da_min` de la correccion L1, justo por encima del suelo de
ruido del MIL (~1.07)— la direccion principal es ruido y no se avisa: una
estructura isotropa no esta mal orientada en ninguna direccion.

Las funciones devuelven datos, no frases: la interfaz las traduce con `_()` y
un guion puede usarlas sin arrastrar textos.
"""

from __future__ import annotations

import numpy as np

BVTV_NO_TRABECULAR = 0.83
ANGULO_MECANICO_MAX_DEG = 30.0
DA_MIN_EJE = 1.10
# Por debajo de este DA2 el candidato es PLANO (los dos autovalores menores del
# MIL iguales): no tiene direccion principal, tiene normal. Mismo umbral que la
# alineacion impuesta L2 (`fit.alinear_marco_fabrica`).
DA2_PLANO = 1.06


def _escalar(m, clave):
    if not isinstance(m, dict):
        return None
    v = m.get(clave)
    try:
        x = float(np.asarray(v, float).reshape(()))
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def voi_no_trabecular(m_o_bvtv, umbral=BVTV_NO_TRABECULAR):
    """¿El VOI es demasiado denso para ser hueso trabecular?

    Acepta el dict de morfometria o directamente un BV/TV. Devuelve
    {"aviso": bool, "BVTV": float|None, "umbral": float}. Sin BV/TV medible,
    no se avisa: no hay nada que afirmar.
    """
    if isinstance(m_o_bvtv, dict):
        bv = _escalar(m_o_bvtv, "BVTV")
    else:
        try:
            bv = float(m_o_bvtv)
        except (TypeError, ValueError):
            bv = None
        if bv is not None and not np.isfinite(bv):
            bv = None
    return {"aviso": bool(bv is not None and bv >= umbral), "BVTV": bv,
            "umbral": float(umbral)}


def _eje(m):
    d = m.get("dir_principal") if isinstance(m, dict) else None
    if d is None:
        return None
    d = np.asarray(d, float).ravel()
    if d.size != 3 or not np.all(np.isfinite(d)) or np.linalg.norm(d) <= 0:
        return None
    return d / np.linalg.norm(d)


def _marco(m):
    V = m.get("eigenvectors") if isinstance(m, dict) else None
    if V is None:
        return None
    V = np.asarray(V, float)
    if V.shape != (3, 3) or not np.all(np.isfinite(V)):
        return None
    return V


def desalineacion(m_voi, m_cand, umbral_deg=ANGULO_MECANICO_MAX_DEG,
                  da_min=DA_MIN_EJE, da2_plano=DA2_PLANO):
    """Angulo entre los ejes principales del VOI y del candidato.

    CANDIDATO PLANO. Si el DA2 del candidato es menor que `da2_plano`, su
    direccion principal es una cualquiera dentro del plano y compararla no
    dice nada. Se compara entonces la NORMAL del plano (tercer autovector del
    MIL) con el eje menor del VOI, que es el criterio de la alineacion
    impuesta L2. Medido sobre el VOI proximal de H4: tras L2 la normal quedo a
    5.5 grados, y la comparacion de direcciones principales daba 80 — un aviso
    falso sobre un candidato bien colocado.

    Devuelve {"aviso", "angulo_deg", "umbral_deg", "evaluable", "motivo",
    "eje"}. `eje` es "direccion principal" o "normal del plano". `motivo` es
    un codigo, no un texto:
        ""                  evaluado
        "sin_datos"         falta la direccion principal de alguno
        "voi_isotropo"      DA del VOI < da_min: su eje es ruido
        "candidato_isotropo" DA del candidato < da_min
    Solo `aviso=True` cuando es evaluable y el angulo supera el umbral.
    """
    out = {"aviso": False, "angulo_deg": float("nan"),
           "umbral_deg": float(umbral_deg), "evaluable": False, "motivo": "",
           "eje": "direccion principal"}
    a, b = _eje(m_voi), _eje(m_cand)
    if a is None or b is None:
        out["motivo"] = "sin_datos"
        return out
    # Como eje y no como vector: el signo de un autovector es arbitrario.
    out["angulo_deg"] = float(np.degrees(np.arccos(min(1.0, abs(float(a @ b))))))
    da_v, da_c = _escalar(m_voi, "DA"), _escalar(m_cand, "DA")
    if da_v is None or da_v < da_min:
        out["motivo"] = "voi_isotropo"
        return out
    if da_c is None or da_c < da_min:
        out["motivo"] = "candidato_isotropo"
        return out
    da2_c = _escalar(m_cand, "DA2")
    Vv, Vc = _marco(m_voi), _marco(m_cand)
    if (da2_c is not None and da2_c < da2_plano and Vv is not None
            and Vc is not None):
        nv, nc = Vv[:, 2], Vc[:, 2]
        c = abs(float(nv @ nc)) / (np.linalg.norm(nv) * np.linalg.norm(nc))
        out["angulo_deg"] = float(np.degrees(np.arccos(min(1.0, c))))
        out["eje"] = "normal del plano"
    out["evaluable"] = True
    out["aviso"] = bool(out["angulo_deg"] > umbral_deg)
    return out
