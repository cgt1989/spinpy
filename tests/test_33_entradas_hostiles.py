"""
Bloque 33 - Entradas hostiles o fuera de rango (verificacion Red Team).

LO QUE SE VERIFICA
  (1) Lector VTK: espaciado nulo, negativo, NaN o infinito se rechaza; un
      archivo sin espaciado avisa en vez de suponer 1.0 en silencio; los
      valores NaN y las dimensiones nulas se rechazan.
  (2) Lector .mat: espaciado invalido se rechaza.
  (3) `elastic.homogeneizar`: nu fuera de (-1, 0.5), E <= 0, voxel <= 0 y
      escala_vacio fuera de [0, 1] devuelven NaN con ok = False y el motivo,
      sin excepcion ni tensor aparentemente valido.
  (4) `escribir_febio_ensayo`: un nombre de archivo con &, comillas o < y un
      comentario con guiones o caracteres no ASCII dan un XML bien formado.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  Son comprobaciones logicas (rechaza / no rechaza); no hay tolerancia numerica.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from spinpy import elastic, escribe, leer_vtk_voi


def _vtk(ruta, dims="4 4 4", spacing="1 1 1", n=64, valores=None):
    v = valores if valores is not None else " ".join(["1", "0"] * (n // 2))
    lineas = ["# vtk DataFile Version 3.0", "x", "ASCII",
              "DATASET STRUCTURED_POINTS", f"DIMENSIONS {dims}"]
    if spacing is not None:
        lineas.append(f"SPACING {spacing}")
    lineas += ["ORIGIN 0 0 0", f"POINT_DATA {n}",
               "SCALARS s unsigned_char 1", "LOOKUP_TABLE default", v]
    ruta.write_text("\n".join(lineas) + "\n")
    return ruta


def test_vtk_valido_se_lee(tmp_path):
    vol, sp = leer_vtk_voi(_vtk(tmp_path / "ok.vtk", spacing="0.05 0.05 0.05"))
    assert vol.shape == (4, 4, 4)
    assert np.allclose(sp, 0.05)


@pytest.mark.parametrize("spacing", ["0 0 0", "-1 1 1", "nan 1 1", "inf 1 1",
                                     "abc 1 1"])
def test_vtk_espaciado_invalido_se_rechaza(tmp_path, spacing):
    with pytest.raises(ValueError):
        leer_vtk_voi(_vtk(tmp_path / "s.vtk", spacing=spacing))


def test_vtk_sin_espaciado_avisa(tmp_path):
    with pytest.warns(UserWarning, match="espaciado"):
        _, sp = leer_vtk_voi(_vtk(tmp_path / "sin.vtk", spacing=None))
    assert np.allclose(sp, 1.0)


def test_vtk_nan_en_los_datos_se_rechaza(tmp_path):
    with pytest.raises(ValueError, match="NaN"):
        leer_vtk_voi(_vtk(tmp_path / "n.vtk", valores=" ".join(["nan"] * 64)))


def test_vtk_dimension_nula_se_rechaza(tmp_path):
    with pytest.raises(ValueError):
        leer_vtk_voi(_vtk(tmp_path / "d.vtk", dims="0 4 4", n=0, valores=""))


def test_mat_espaciado_invalido_se_rechaza(tmp_path):
    from scipy.io import savemat

    from spinpy import leer_mat_voi
    ruta = tmp_path / "v.mat"
    savemat(str(ruta), {"VOI": np.ones((4, 4, 4)), "spacing": [0.0, 1.0, 1.0]})
    with pytest.raises(ValueError):
        leer_mat_voi(ruta)


@pytest.mark.parametrize("kw", [
    {"nu_s": 0.5}, {"nu_s": 0.6}, {"nu_s": -1.0}, {"nu_s": -2.0},
    {"nu_s": float("nan")}, {"E_s": -1.0}, {"E_s": 0.0},
    {"vox_size": 0.0}, {"vox_size": -1.0}, {"escala_vacio": -0.1},
    {"escala_vacio": 1.5},
])
def test_homogeneizar_parametros_invalidos(kw):
    bw = np.random.default_rng(0).random((4, 4, 4)) < 0.5
    C, info = elastic.homogeneizar(bw, **kw)
    assert info["ok"] is False
    assert info["msg"]
    assert np.isnan(C).all()


def test_homogeneizar_parametros_validos_siguen_resolviendo():
    bw = np.random.default_rng(0).random((4, 4, 4)) < 0.5
    C, info = elastic.homogeneizar(bw, E_s=1.0, nu_s=0.3)
    assert info["ok"] is True
    assert np.isfinite(C).all()


@pytest.mark.parametrize("nombre", ["normal", "a&b", 'a"onda', "x<y",
                                    "muestra_ñandú"])
def test_febio_xml_bien_formado(tmp_path, nombre):
    nodos = np.array([[x, y, z] for z in (0, 1) for y in (0, 1) for x in (0, 1)],
                     float)
    elems = np.array([[0, 1, 3, 2, 4, 5, 7, 6]])
    ruta = tmp_path / f"{nombre}.feb"
    escribe.escribir_febio_ensayo(nodos, elems, ruta,
                                  comentario="muéstra µm --- x--")
    ET.parse(ruta)
    assert "µm" in ruta.read_text(encoding="utf-8")
