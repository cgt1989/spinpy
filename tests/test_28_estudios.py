"""
Bloque 28: estudios complementarios del informe (`spinpy.estudios`) y sus
figuras (`spinpy.figuras_estudios`).

REFERENCIAS
  Hara T et al. Bone 2002;31:107. Un cambio de umbral mueve BV/TV y rigidez.
  Pahr DH, Zysset PK. Biomech Model Mechanobiol 2008;7:463. La rigidez
  calculada depende de la condicion de contorno.
  Guo Y et al. Adv Intell Syst 2024;6:2300789. Perfil de curvaturas.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  umbral equivalente    t(0) = 0,5 exacto; t(-d) = 1 - t(d) a 1e-12.
  desplazamiento        losa de t voxeles alineada con la rejilla:
                        delta = +1 -> t + 2 voxeles; -1 -> t - 2; 0 -> t.
                        Exacto: la interfaz plana es el caso en que la
                        equivalencia umbral-desplazamiento es una identidad.
                        (Con delta = +-0,5 una cara alineada queda justo en el
                        umbral y el resultado es un empate; no se prueba.)
  contorno, macizo      deslizante: E/E_s = 1 a 1e-6 (estado uniaxial);
                        periodica: E_z/E_s = 1 a 1e-6;
                        empotrado: 1 <= E/E_s <= M_oed/E_s (+1e-9).
  curvatura de una bola fraccion convexa >= 0,9 y silla <= 0,05: una esfera
                        es convexa en todas partes; el margen cubre la
                        discretizacion de la malla.
  figuras               las ocho salen de un documento con sus datos, en ES
                        y EN, y ninguna sale de un documento sin ellos.

RESULTADO DE LA PRIMERA EJECUCION, DEJADO AQUI A PROPOSITO
  La primera version daba la bola CONCAVA (H = -1/R): `curvaturas_de_campo`
  orienta la normal hacia donde el campo crece, y la mascara suavizada crece
  hacia el hueso. Corregido pasando su complemento. Con eso H = 1/R a un
  0,2 % (R = 12) y a un 0,1 % (R = 20), y se comprueba en su propia prueba.

  Las FRACCIONES por signo no cumplen lo declarado: el 9,4 % del area de la
  bola de R = 12 voxeles sale como silla, y el 18,9 % con R = 20. No es el
  signo: es que la escalera de voxeles, que la gaussiana de sigma = 1 no borra
  del todo, deja k2 cerca de cero con signo aleatorio, y cuanto mayor el radio
  menor es 1/R frente a ese ruido. La cota no se mueve; la prueba queda como
  xfail estricto y el pie de la figura 15 declara el sesgo medido. La media H,
  que promedia el ruido, si es fiable.
"""
import json

import numpy as np
import pytest

from conftest import bola, losa
from spinpy import estudios as S
from spinpy import figuras_estudios as FE
from spinpy import informe

BLOQUE = "28 Estudios complementarios del informe"
REF = "estudios.py; Hara 2002; Pahr y Zysset 2008; Guo 2024"
NU = 0.30
M_OED = (1 - NU) / ((1 + NU) * (1 - 2 * NU))


def test_umbral_equivalente(registro):
    t0 = S.nivel_de_desplazamiento(0.0)
    sim = max(abs(S.nivel_de_desplazamiento(-d) - (1 - S.nivel_de_desplazamiento(d)))
              for d in (0.25, 0.5, 1.0, 2.0))
    ok = t0 == 0.5 and sim < 1e-12 and S.nivel_de_desplazamiento(1) < 0.5
    registro.anotar(BLOQUE, "t = Phi(-delta/sigma): t(0) = 0,5 y simetria",
                    REF, 0.5, t0, 1e-12, "identidad", ok,
                    nota=f"asimetria max {sim:.1e}")
    assert ok


def test_desplazamiento_de_una_losa(registro):
    t = 8
    BW = losa(n=32, espesor=t)
    medidos = {}
    for d in (-1.0, 0.0, 1.0):
        bw = S.desplazar_superficie(BW, d)
        medidos[d] = int(bw[16, 16, :].sum())
    esperado = {-1.0: t - 2, 0.0: t, 1.0: t + 2}
    ok = medidos == esperado
    registro.anotar(BLOQUE, "losa: delta = -1, 0, +1 voxel -> t-2, t, t+2",
                    REF, None, None, "exacto", "interfaz plana", ok,
                    nota=str(medidos))
    assert ok, medidos


def test_condiciones_contorno_en_un_cubo_macizo(registro):
    r = S.condiciones_contorno(np.ones((10, 10, 10), bool), [1.0] * 3,
                               n_mec=10, E_s=1.0, nu_s=NU)
    desl = r["deslizante"]["E_rel"]
    emp = r["empotrado"]["E_rel"]
    per = r["periodica"]["E_rel"]
    ok = (abs(desl - 1) < 1e-6 and abs(per - 1) < 1e-6
          and 1 - 1e-9 <= emp <= M_OED + 1e-9)
    registro.anotar(BLOQUE, "macizo: deslizante = periodica = E_s; "
                    "empotrado entre E_s y M_oed", "Hooke; modulo edometrico",
                    1.0, desl, 1e-6, "relativo y desigualdad", ok,
                    nota=f"deslizante {desl:.9f}, periodica {per:.9f}, "
                         f"empotrado {emp:.6f} (M_oed {M_OED:.4f})")
    assert ok


def test_tensor_de_sesion_evita_homogeneizar(registro):
    """Con el tensor de la sesion la periodica sale de el, no de un calculo."""
    C = np.eye(6) * 7.0
    r = S.condiciones_contorno(np.ones((6, 6, 6), bool), [1.0] * 3, n_mec=6,
                               E_s=1.0, nu_s=NU, C_periodica=C)
    ok = r["periodica"]["origen"] == "sesion" and abs(
        r["periodica"]["E_rel"] - 7.0) < 1e-9
    registro.anotar(BLOQUE, "periodica: se usa el tensor de la sesion", REF,
                    7.0, r["periodica"]["E_rel"], 1e-9, "exacto", ok)
    assert ok


def test_curvatura_media_de_una_bola(registro):
    """H = 1/R: la media promedia el ruido de la escalera y es fiable."""
    R = 12.0
    c = S.curvaturas(bola(n=40, radio=R), [1.0] * 3)
    H = S.perfiles_curvatura({"bola": c})["bola"]["resumen"]["H_medio"]
    ok = abs(H * R - 1.0) <= 0.05
    registro.anotar(BLOQUE, "bola R = 12: H = 1/R (signo y valor)",
                    "Guo 2024; do Carmo", 1.0 / R, H, "5 %", "relativo", ok)
    assert ok, H


@pytest.mark.xfail(strict=True, reason=(
    "hallazgo documentado: sobre una esfera digitalizada el 9 % del area sale "
    "como silla por el ruido de la escalera; ver el docstring del modulo"))
def test_curvatura_de_una_bola(registro):
    c = S.curvaturas(bola(n=40, radio=12.0), [1.0] * 3)
    p = S.perfiles_curvatura({"bola": c})["bola"]
    r = p["resumen"]
    ok = r["convexa"] >= 0.9 and r["silla"] <= 0.05 and len(p["perfil"]) == \
        S.NBINS_CURVATURA
    registro.anotar(BLOQUE, "bola: superficie convexa", "Guo 2024; do Carmo",
                    1.0, r["convexa"], ">= 0,9", "fraccion de area", ok,
                    nota=f"silla {r['silla']:.3f}")
    assert ok


def _doc_con_estudios():
    """Documento minimo con un dato de cada figura 12-19."""
    rng = np.random.default_rng(0)
    pasos = [{"paso": k, "fase": "inicial" if k == 0 else "perdida",
              "hueso_rel": 1 - 0.05 * k, "BVTV": 0.4 * (1 - 0.05 * k),
              "TbTh": 0.14, "TbN": 2.9, "ConnD": 10.0 * (1 - 0.1 * k),
              "E_rel": (1 - 0.05 * k) ** 3,
              "sigma_fallo_rel": (1 - 0.05 * k) ** 2} for k in range(5)]
    fallo = [{"paso": k, "F_fallo": 100.0 + 10 * k - 3 * k * k,
              "E_rel": 1 - 0.08 * k, "dano_acumulado": 0.02 * k}
             for k in range(6)]
    P = np.zeros((S.NBINS_CURVATURA,) * 2)
    P[40, 20] = 0.7
    P[45, 45] = 0.3
    hist = np.histogram(rng.normal(0, 0.4, 2000), bins=40, range=(-1, 1))[0]
    forma = {"EF": {"EF": 0.0, "EF_frac_placa": 0.3, "EF_frac_barra": 0.3,
                    "EF_relleno": 0.97,
                    "histograma": (hist / hist.sum()).tolist(),
                    "bordes": np.linspace(-1, 1, 41).tolist()},
             "curvatura": {"perfil": P.tolist(), "limite": 30.0,
                           "resumen": {"H_medio": -1.0, "silla": 0.7,
                                       "convexa": 0.1, "concava": 0.2}}}
    C = (np.diag([3e9, 3e9, 5e9, 1e9, 1e9, 1e9])
         + np.pad(np.full((3, 3), 1e9) - np.diag([1e9] * 3), ((0, 3), (0, 3))))
    sens = {"filas": [{"desplazamiento_vox": d, "desplazamiento_um": 16 * d,
                       "BVTV": 0.4 + 0.05 * d, "TbTh": 0.14 + 0.01 * d,
                       "DA": 1.3, "E_app": 2e9 * (1 + 0.3 * d)}
                      for d in (-1.0, 0.0, 1.0)]}
    reg = {"estructura": "VOI", "malla": "hex8", "nombre": "app",
           "lineal": {"E_app": 4e9, "pistoia": {"vm_p99_superficie": 3e7,
                                                "sigma_fallo": 2e7}},
           "app": {"ok": True, "E_app": 4e9 * (1 + 1e-8),
                   "pistoia": {"vm_p99_superficie": 3e7,
                               "sigma_fallo": 2e7}}}
    return {
        "voi": {"nombre": "sintetico", "forma": [64] * 3,
                "spacing_mm": [0.016] * 3},
        "morfometria_voi": {"BVTV": 0.4, "TbTh": 0.14, "TbSp": 0.2,
                            "TbN": 2.9, "DA": 1.3},
        "resultados": {
            "simulacion_perdida": {"protocolo": "adelgazamiento",
                                   "pasos": pasos, "estructura_codigo": "voi",
                                   "parametros": {"pasos": 4,
                                                  "perdida_paso": 0.05,
                                                  "n_mec": 32}},
            "fallo_progresivo": {"pasos": fallo, "estructura_codigo": "voi",
                                 "parametros": {"pasos": 5, "n_mec": 32,
                                                "rigidez_danada": 0.05},
                                 "resumen": {"F_max": 108.0, "paso_F_max": 2,
                                             "colapso_en_paso": None,
                                             "E_rel_colapso": 0.5}},
            "elastico": {"resolucion": 32, "E_s_Pa": 2e10, "nu_s": 0.3,
                         "por_estructura": {"voi": {"C_Pa": C.tolist(),
                                                    "Ex": 1e9}}},
            "forma": {"por_estructura": {"voi": forma},
                      "parametros": {"sigma_vox": 1.0}},
            "sensibilidad_superficie": {"por_estructura": {"voi": sens},
                                        "parametros": {"n_mec": 40,
                                                       "apoyo": "deslizante",
                                                       "sigma_vox": 1.0}},
            "condiciones_contorno": {"por_estructura": {"voi": {
                "n": 32, "deslizante": {"E_rel": 0.1},
                "empotrado": {"E_rel": 0.12}, "periodica": {"E_rel": 0.14}}}},
            "febio": {"registros": [reg]},
        }}


def test_figuras_salen_solo_con_su_dato(registro, tmp_path):
    doc = _doc_con_estudios()
    hechas, vacias = [], []
    for clave, nombres, hacer in FE.FIGURAS:
        for i, idioma in enumerate(("es", "en")):
            if hacer(doc, tmp_path / nombres[i], idioma):
                hechas.append(f"{clave}/{idioma}")
                informe.PIES[clave][i].format(
                    **informe.datos_pies(doc, idioma))
        if hacer({"resultados": {}}, tmp_path / "vacia", "es"):
            vacias.append(clave)
    ok = len(hechas) == 2 * len(FE.FIGURAS) and not vacias
    registro.anotar(BLOQUE, "figuras 12-19: con dato salen, sin dato no",
                    REF, float(2 * len(FE.FIGURAS)), float(len(hechas)),
                    "exacto", "ES y EN", ok, nota=f"sin dato: {vacias}")
    assert ok, (hechas, vacias)


def test_metodos_y_referencias(registro):
    doc = _doc_con_estudios()
    md = informe.informe_markdown(doc, informe.comprobar(doc), {}, "es")
    ok = all(doi in md for doi in ("10.1016/s8756-3282(02)00782-2",
                                   "10.1007/s10237-007-0109-7",
                                   "10.1002/aisy.202300789"))
    ok = ok and "sustituto de la incertidumbre de segmentación" in md
    filas = __import__("spinpy.lista_chequeo", fromlist=["x"]).evaluar(doc)
    seg = next(f for f in filas if f["item"] == "Segmentación")
    ok = ok and seg["estado"] == "parcial"
    registro.anotar(BLOQUE, "metodos citan Hara, Pahr y Guo; la lista de "
                    "chequeo reconoce la sensibilidad", REF, None, None,
                    "exacto", "presentes", ok, nota=seg["estado"])
    assert ok


def test_completar_guarda_en_el_documento(registro):
    """forma + contorno sobre una estructura pequena, sin tocar lo que hay."""
    from spinpy.grf import generar_mascara
    BW, _, _ = generar_mascara(rho=0.4, wave_number=10 * np.pi,
                               num_waves=200, thetas=[90, 90, 90],
                               resolution=24, seed=5)
    doc = {"resultados": {"sensibilidad_superficie": {"marca": 1}}}
    orig = json.dumps(doc["resultados"]["sensibilidad_superficie"])
    hechos = S.completar(doc, [("spinodoide", BW, [0.05] * 3)],
                         que=("superficie", "contorno"), n_mec=12)
    r = doc["resultados"]
    ok = (hechos == ["contorno"]
          and json.dumps(r["sensibilidad_superficie"]) == orig
          and r["condiciones_contorno"]["por_estructura"]["spinodoide"][
              "periodica"]["origen"] == "calculada")
    json.dumps(r)                              # serializable
    registro.anotar(BLOQUE, "completar: calcula lo que falta, respeta lo que "
                    "hay", REF, None, None, "exacto", "claves del documento",
                    ok, nota=str(hechos))
    assert ok
