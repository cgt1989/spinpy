"""
Bloque 18 — Muestreo del MIL (correccion M1, opcional).

Sobre voxeles, la direccion principal del MIL se tuerce hacia la rejilla cuando
la estructura esta oblicua (Estudio_MIL). M1 anade `muestreo="suavizado"` a
`tensor_mil` SIN tocar el camino por omision, que sigue siendo el de MATLAB.

Tolerancias declaradas ANTES de medir en este bloque; las cifras de referencia
son las de Estudio_MIL/resultados/rotacion.jsonl y da.jsonl (48^3, R0 =
euler(30,50,20), estiramiento 2.5, semillas 1-5):
    voxel      13.3 +- 0.2 a 13.9 +- 0.5 grados
    suavizado   1.5 +- 0.6 grados
"""

import numpy as np
import pytest

from spinpy.dual_lattice import generar_dual_lattice
from spinpy.grf import euler_R
from spinpy.morphometry import morfometria, tensor_mil

BLOQUE = "18 muestreo del MIL"


def _ang(u, v):
    return float(np.degrees(np.arccos(min(1.0, abs(float(u @ v))
                                          / np.linalg.norm(u) / np.linalg.norm(v)))))


def test_omision_intacta(registro):
    """El camino por omision no cambia: 'voxel' explicito == sin argumento."""
    BW = generar_dual_lattice(24, 3.0, 0.3, (1, 1, 2.0), seed=1)[0]
    sp = np.full(3, 0.1)
    a = tensor_mil(BW, sp)
    b = tensor_mil(BW, sp, muestreo="voxel")
    ok = (np.array_equal(a["MIL"], b["MIL"]) and a["DA"] == b["DA"]
          and np.array_equal(a["eigenvectors"], b["eigenvectors"])
          and a["muestreo"] == "voxel")
    m = morfometria(BW, sp)
    ok &= m["DA"] == a["DA"]
    registro.anotar(BLOQUE, "omision bit a bit", "construccion", 0.0,
                    0.0 if ok else 1.0, "exacto", "voxel == por omision", ok)
    assert ok


def test_muestreo_invalido():
    BW = np.zeros((8, 8, 8), bool)
    BW[2:5] = True
    with pytest.raises(ValueError):
        tensor_mil(BW, np.full(3, 0.1), muestreo="trilineal")


@pytest.mark.lento
def test_sesgo_orientacion(registro):
    """Estructura girada de forma exacta: el suavizado quita el sesgo."""
    sp = np.full(3, 0.1)
    R0 = euler_R(30, 50, 20)
    eje = R0[:, 2]
    ang_v, ang_s, da_s = [], [], []
    for s in (1, 2, 3):
        BW = generar_dual_lattice(48, 5.0, 0.3, (1, 1, 2.5), R=R0, seed=s)[0]
        ang_v.append(_ang(tensor_mil(BW, sp)["dir_principal"], eje))
        f = tensor_mil(BW, sp, muestreo="suavizado")
        ang_s.append(_ang(f["dir_principal"], eje))
        da_s.append(f["DA"])
        BWi = generar_dual_lattice(48, 5.0, 0.3, (1, 1, 2.5), seed=s)[0]
        da_s.append(tensor_mil(BWi, sp, muestreo="suavizado")["DA"])
    ok_v = np.mean(ang_v) > 8.0          # el sesgo existe (si no, el bloque no prueba nada)
    ok_s = np.mean(ang_s) <= 4.0
    registro.anotar(BLOQUE, "sesgo con voxel (girado 30/50/20)", "Estudio_MIL",
                    13.3, float(np.mean(ang_v)), "> 8 grados", "media 3 semillas", ok_v)
    registro.anotar(BLOQUE, "sesgo con suavizado (girado 30/50/20)", "Estudio_MIL",
                    1.5, float(np.mean(ang_s)), "<= 4 grados", "media 3 semillas", ok_s,
                    nota=f"DA suavizado girado/alineado {np.round(da_s, 3).tolist()}")
    assert ok_v and ok_s
