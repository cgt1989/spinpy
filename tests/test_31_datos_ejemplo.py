"""
Bloque 31: datos de ejemplo abiertos y ejemplo de inicio rapido.

QUE SE VERIFICA
  (1) Los 10 VOIs de `datos/porcino_koria2020/` (Koria, Mengoni y Brockett
      2020, CC BY 4.0) coinciden con las huellas de `SHA256SUMS.txt`: un
      archivo alterado o truncado se detecta antes de producir resultados.
  (2) Cada uno se lee con `leer_voi` como un volumen booleano de 188^3 voxeles
      de 0,016 mm, con una fraccion de hueso plausible para hueso subcondral
      (0,30 a 0,45; medido 0,333 a 0,398).
  (3) La atribucion exigida por la licencia esta en FUENTE.txt (DOI y CC BY).
  (4) `ejemplos/inicio_rapido.py --rapido` corre de punta a punta y devuelve 0.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  (2) BV/TV entre 0,30 y 0,45: el README documenta 0,40 para V1 y el resto del
      conjunto es hueso subcondral del mismo tipo.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pytest

from spinpy import leer_voi

RAIZ = Path(__file__).resolve().parent.parent
DATOS = RAIZ / "datos" / "porcino_koria2020"
NOMBRES = [f"VOI_{g}{i}.mat" for g in "CV" for i in range(1, 6)]


def _huellas():
    out = {}
    for linea in (DATOS / "SHA256SUMS.txt").read_text().splitlines():
        if linea.strip():
            h, n = linea.split()
            out[n] = h
    return out


def test_manifiesto_cubre_los_diez_vois():
    assert sorted(_huellas()) == sorted(NOMBRES)


@pytest.mark.parametrize("nombre", NOMBRES)
def test_huella_sha256(nombre):
    h = hashlib.sha256((DATOS / nombre).read_bytes()).hexdigest()
    assert h == _huellas()[nombre]


@pytest.mark.parametrize("nombre", NOMBRES)
def test_voi_se_lee_con_forma_escala_y_bvtv_plausibles(nombre):
    BW, sp = leer_voi(DATOS / nombre)
    assert BW.dtype == bool and BW.shape == (188, 188, 188)
    np.testing.assert_allclose(sp, 0.016, rtol=1e-9)
    assert 0.30 <= BW.mean() <= 0.45


def test_atribucion_de_la_licencia():
    t = (DATOS / "FUENTE.txt").read_text()
    assert "10.5518/787" in t and "CC BY 4.0" in t and "Koria" in t


@pytest.mark.lento
def test_inicio_rapido_corre_de_punta_a_punta(capsys):
    sys.path.insert(0, str(RAIZ / "ejemplos"))
    import inicio_rapido
    assert inicio_rapido.main(["--rapido"]) == 0
    salida = capsys.readouterr().out
    assert "coincide con SHA256SUMS.txt" in salida
    assert "E_app" in salida
