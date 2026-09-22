"""
Bloque 7 — Medida sobre malla suavizada y reparada, y estanqueidad.

REFERENCIAS
  Taubin G. A signal processing approach to fair surface design. SIGGRAPH
  1995. Suavizado que no encoge: conserva el volumen a diferencia del
  laplaciano.
  Attene M. A lightweight approach to repairing digitized polygon meshes.
  Vis Comput 2010;26:1393. Es el algoritmo de PyMeshFix.
  Parfitt 1987 para las derivadas (Tb.Th = 2BV/BS, etc.).

QUE SE VERIFICA
  `morfometria_malla` mide sobre la superficie cerrada + Taubin + PyMeshFix,
  recortada a la caja del VOI y sin las seis tapas. Contra figuras de area y
  volumen conocidos debe estar mucho mas cerca que el modo de voxeles, cuyo
  +8.5 % de area esta medido en validar_analitico.py.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  esfera r=14:   BV a 1 % del exacto;  BS a 2 % del exacto  (medido a mano el
                 mismo dia: +0.1 % y +0.6 %; se deja margen)
  losa t=10:     BS = 2 * L^2 a 1 % con L = (n-1)*sp, la caja de centros de
                 voxel en la que mide el modo malla: las dos caras grandes son
                 interfaz y las cuatro laterales son tapas, excluidas.
                 BV = L^2 * t a 1 %.
  TV:            TV_malla = ((n-1)/n)^3 * TV_voxel a 1e-12 (por diseno).
  spinodoide:    estanca; perdida_vs_mc (suavizar + reparar, misma caja)
                 dentro del 3 %; BS_malla / BS_voxel en [0.85, 0.97]: el
                 suavizado solo puede QUITAR area, y lo que quita es el sesgo
                 de voxel de ~8.5 % medido en la esfera. BV/TV malla frente a
                 voxel se REPORTA con su descomposicion y solo se acota por
                 sanidad en el 8 %: es diferencia de definicion, no perdida.
  malla_tet10:   estanca = True y 0 bordes abiertos; vol_pct_MC >= 95 a 40^3
                 (medido -1.9 % a 48^3 con decimate 0.5; a 40^3 las trabeculas
                 son mas finas, se deja margen).

RESULTADO DE LA PRIMERA EJECUCION, DEJADO AQUI A PROPOSITO
  La primera version recortaba en la caja completa y la losa dio BS +27 %:
  el suavizado redondea las aristas contra las caras del cubo y las puntas
  redondeadas contaban como interfaz. Se paso a recortar medio voxel por
  dentro (ver `morfometria_malla`). La primera version tambien exigia
  |BV_malla/BV_voxel - 1| <= 6 % y a 40^3 salio -7.1 %, de los que -4.6 %
  eran definicion (MC a nivel 0.5 frente a conteo de cubos) y -2.6 % la
  reparacion de puntas redondeadas ya inexistentes en la version final. Se
  separo en dos numeros porque acotar la suma mezclaba una definicion con
  una perdida.
"""
import numpy as np
import pytest

from conftest import bola, losa
from spinpy import generar_mascara
from spinpy.morphometry import morfometria, morfometria_malla

BLOQUE = "07 Medida sobre malla y estanqueidad"
REF = "Taubin 1995; Attene 2010 (PyMeshFix); Parfitt 1987"


def test_esfera_area_y_volumen(registro):
    r, n = 14.0, 40
    BW = bola(n=n, radio=r)
    m = morfometria_malla(BW, [1.0] * 3)
    A_ex, V_ex = 4 * np.pi * r ** 2, 4 / 3 * np.pi * r ** 3
    eA = (m["BS"] - A_ex) / A_ex
    eV = (m["BV"] - V_ex) / V_ex
    mv = morfometria(BW, [1.0] * 3)
    eA_vox = (mv["BS"] - A_ex) / A_ex
    registro.anotar(BLOQUE, "esfera r=14: BS malla", REF, A_ex, m["BS"], "2 %",
                    "relativo al exacto", abs(eA) <= 0.02,
                    nota="voxeles: %+.1f %%" % (100 * eA_vox))
    registro.anotar(BLOQUE, "esfera r=14: BV malla", REF, V_ex, m["BV"], "1 %",
                    "relativo al exacto", abs(eV) <= 0.01)
    registro.anotar(BLOQUE, "esfera: la malla mejora el area de voxeles", REF,
                    None, abs(eA) / abs(eA_vox), "< 0.5",
                    "cociente de errores", abs(eA) < 0.5 * abs(eA_vox))
    assert abs(eA) <= 0.02, (m["BS"], A_ex)
    assert abs(eV) <= 0.01, (m["BV"], V_ex)
    assert abs(eA) < 0.5 * abs(eA_vox)
    assert m["estanca"] and m["BS_tapas"] == 0.0


def test_losa_tapas_excluidas(registro):
    n, t = 32, 10
    BW = losa(n=n, espesor=t)
    m = morfometria_malla(BW, [1.0] * 3)
    L = float(n - 1)                    # caja de centros de voxel
    A_ex, V_ex = 2 * L * L, L * L * t
    eA = (m["BS"] - A_ex) / A_ex
    eV = (m["BV"] - V_ex) / V_ex
    registro.anotar(BLOQUE, "losa t=10: BS = 2 L^2 (tapas excluidas)", REF,
                    A_ex, m["BS"], "1 %", "relativo", abs(eA) <= 0.01,
                    nota="tapas excluidas: %.1f" % m["BS_tapas"])
    registro.anotar(BLOQUE, "losa t=10: BV = L^2 t", REF, V_ex, m["BV"], "1 %",
                    "relativo", abs(eV) <= 0.01)
    assert abs(eA) <= 0.01, (m["BS"], A_ex, m["BS_tapas"])
    assert abs(eV) <= 0.01, (m["BV"], V_ex)
    assert m["BS_tapas"] > 0          # las cuatro caras laterales SON tapas


def test_tv_es_la_caja_de_centros(registro):
    n = 24
    BW = bola(n=n, radio=8.0)
    sp = [0.05, 0.05, 0.05]
    a = morfometria(BW, sp)["TV"]
    m = morfometria_malla(BW, sp)
    esp = a * ((n - 1) / n) ** 3
    ok = abs(m["TV"] - esp) / esp <= 1e-12 and abs(m["TV_voxel"] - a) <= 1e-12
    registro.anotar(BLOQUE, "TV malla = ((n-1)/n)^3 TV voxeles", REF, esp,
                    m["TV"], 1e-12, "relativo; caja de centros", ok)
    assert ok


@pytest.mark.lento
def test_spinodoide_coherente_y_estanco(registro):
    N = 40
    sp = [4.994 / N] * 3
    BW, _, _ = generar_mascara(rho=0.35, wave_number=15 * np.pi, num_waves=700,
                               thetas=[15, 15, 45], resolution=N, seed=20260720)
    m = morfometria_malla(BW, sp)
    mv = morfometria(BW, sp)
    dBVTV = m["BVTV"] / mv["BVTV"] - 1.0        # densidades, cada una en su caja
    rBS = m["BS"] / mv["BS"]
    registro.anotar(BLOQUE, "spinodoide 40^3: estanca", REF, None,
                    float(m["bordes_abiertos"]), "0 bordes", "manifold",
                    m["estanca"])
    registro.anotar(BLOQUE, "spinodoide 40^3: perdida vs superficie cruda "
                    "(misma caja)", REF, 0.0, m["perdida_vs_mc_pct"] / 100.0,
                    "|.| <= 3 %", "suavizar + reparar",
                    abs(m["perdida_vs_mc_pct"]) <= 3.0)
    registro.anotar(BLOQUE, "spinodoide 40^3: BS malla / BS voxel", REF, 0.915,
                    rBS, "[0.85, 0.97]", "el suavizado quita el sesgo de voxel",
                    0.85 <= rBS <= 0.97)
    registro.anotar(BLOQUE, "spinodoide 40^3: BV/TV malla / BV/TV voxel - 1",
                    REF, 0.0, dBVTV, "|.| <= 8 % (sanidad)",
                    "diferencia de definicion, se reporta",
                    abs(dBVTV) <= 0.08,
                    nota="definicion MC vs conteo %+.1f %%; perdida vs MC "
                         "%+.1f %%" % (m["dif_BV_vs_voxel_pct"],
                                       m["perdida_vs_mc_pct"]))
    assert m["estanca"]
    assert abs(m["perdida_vs_mc_pct"]) <= 3.0, m["perdida_vs_mc_pct"]
    assert 0.85 <= rBS <= 0.97, rBS
    assert abs(dBVTV) <= 0.08, dBVTV


@pytest.mark.lento
def test_tet10_estanqueidad_y_volumen(registro):
    from spinpy.solido import malla_tet10
    N = 40
    sp = [4.994 / N] * 3
    BW, _, _ = generar_mascara(rho=0.35, wave_number=15 * np.pi, num_waves=700,
                               thetas=[15, 15, 45], resolution=N, seed=20260720)
    nodos, elems, sup, inf = malla_tet10(BW, sp)
    registro.anotar(BLOQUE, "TET10 40^3: superficie estanca antes de tetgen",
                    REF, None, float(inf["bordes_abiertos"]), "0 bordes",
                    "manifold", inf["estanca"],
                    nota="%d componente(s)" % inf["n_componentes"])
    registro.anotar(BLOQUE, "TET10 40^3: volumen frente a superficie cruda",
                    REF, 100.0, inf["vol_pct_MC"], ">= 95 %",
                    "decimate 0.5 + PyMeshFix", inf["vol_pct_MC"] >= 95.0,
                    nota="frente al conteo de voxeles: %.1f %%" % inf["vol_pct_BV"])
    assert inf["estanca"] and inf["bordes_abiertos"] == 0
    assert inf["orden_ok"]
    assert inf["vol_pct_MC"] >= 95.0, inf["vol_pct_MC"]
