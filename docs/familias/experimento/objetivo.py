"""Objetivo ampliado, aplicado EN MEMORIA sobre spinpy (el repositorio no se modifica).

* Terminos nuevos en el error: los tres que spinpy ya admite en `pesos` (PoDm, TbTh_CV, H_desv; son su propia
  sugerencia `PESOS_DISCRIMINANTES`, con TbSp a 0) mas Conn.D, que se registra aqui como termino opcional.
* Los terminos de forma del VOI se miden al MISMO tamano de voxel que los candidatos (el VOI se submuestrea
  al lado del candidato, `remuestrear_bw`): curvatura, espesor local, poro y conectividad dependen del voxel.
  Los terminos base (BVTV, DA, BSBV, TbTh, TbN, PoTot) siguen midiendose como en spinpy (VOI a resolucion nativa).
* peso_da2 = 1 en TODAS las familias (spinpy lo aplica por omision solo en el dual-lattice).
"""
import sys
sys.path.insert(0, '../../..')
import numpy as np
import spinpy.error as E
import spinpy.fit as F
import spinpy.fit_dual as FD
from spinpy import morphometry as M
from spinpy.elastic import remuestrear_bw

RES_EQ = 96
PESOS_EXT = {"TbSp": 0.0, "PoDm": 1.0, "TbTh_CV": 1.0, "H_desv": 1.0, "ConnD": 1.0}   # PESOS_DISCRIMINANTES + ConnD
_orig_forma = M.metricas_forma
_orig_err = E.error_morfometrico


def metricas_forma_ext(BW, spacing, tbth=False, poro=False, curvatura=False, ef=False,
                       conectividad=False, **kw):
    BW = np.asarray(BW, dtype=bool)
    sp = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if sp.size == 1:
        sp = np.repeat(sp, 3)
    if BW.shape[0] > RES_EQ:                       # el VOI: al voxel del candidato
        BW, sp = remuestrear_bw(BW, sp, RES_EQ)
    out = _orig_forma(BW, sp, tbth=tbth, poro=poro, curvatura=curvatura, ef=ef, **kw)
    if conectividad:
        out['ConnD'] = M.conectividad(BW, sp)['ConnD']
    return out


def error_ext(m_a, m_b, *a, **kw):
    kw.setdefault('peso_da2', 1.0)
    return _orig_err(m_a, m_b, *a, **kw)


def activar():
    E.TERMINOS_OPCIONALES["ConnD"] = {"tipo": "escalar", "medida": "conectividad", "coste": "forma"}
    M.metricas_forma = metricas_forma_ext
    F.metricas_forma = metricas_forma_ext
    FD.metricas_forma = metricas_forma_ext
    F.error_morfometrico = error_ext                 # el spinodoide tambien con DA2
    E.error_morfometrico = error_ext
    return PESOS_EXT


def medir_ext(BW, sp):
    """Morfometria base + terminos de forma ampliados (mismo camino para VOI y candidatos)."""
    m = M.morfometria(BW, sp)
    m.update(metricas_forma_ext(BW, sp, tbth=True, poro=True, curvatura=True, conectividad=True))
    return m
