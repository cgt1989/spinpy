"""
Bloque 22 — La seccion «Modelos matematicos» del informe dice lo que hace el
codigo.

QUE SE VERIFICA
  Cada ecuacion de `spinpy/informe_modelos.py` se implementa aqui a mano, tal
  como esta escrita en el informe, y se compara con lo que devuelve el codigo
  que el informe describe. Si alguien cambia el calculo sin cambiar el texto
  (o al reves), falla aqui.

  convenciones  cuantiles Q5 y Q7 de la seccion 1 == numpy "hazen" y "linear";
                sigma_0 = 1e6 F / A con A en mm^2.
  morfometria   BV/TV como media; Tb.Th = 2 BV/BS, Tb.Sp, Tb.N; DA y DA2 de
                los autovalores del tensor MIL.
  generadores   umbral phi0 = sqrt(2) erfinv(2 rho - 1) y fase solida
                {phi <= phi0}; en el dual-lattice k = round(rho n^3) voxeles.
  ajuste        los numeros de onda que cita el texto son los que recorre
                `fit.ajustar_spinodoide`.
  elastico      remuestreo i = round(i'(n-1)/(n'-1)), h' = h n/n'; tensor de un
                solido macizo == D(E, nu) de la ecuacion 2; E_x = 1/S11,
                G_yz = 1/S44, nu_xy = -S12/S11.
  ensayo        sobre una estructura porosa: sigma_vM de la formula del
                informe (y sqrt(3 J2) del desviador) == campo del codigo;
                eps_eff = sqrt(2U/E_s) con U = sigma.eps/2; en un bloque macizo
                E_app = E_s y sigma_vM = sigma_0 (solucion exacta uniaxial).
  estadisticos  capa superficial segun la ecuacion (vecino por cara en el
                vacio, B = 1 fuera) == `capa_superficie`; p99 de la capa con
                Q5 == `vm_p99_superficie`; Pistoia con Q7 == `criterio_pistoia`.
  composicion   todas las ecuaciones de un informe completo, en ES y EN, se
                componen con mathtext (ninguna cae a bloque de codigo).

Tolerancia 1e-9 relativa salvo donde se indica: son identidades, no medidas.
"""
import inspect
import re

import numpy as np
import pytest
from scipy.special import erfinv

from spinpy import generar_mascara, informe, morfometria
from spinpy import informe_modelos as IM
from spinpy.dual_lattice import generar_dual_lattice
from spinpy.elastic import (constantes_ingenieria, homogeneizar,
                            remuestrear_bw)
from spinpy.fit import ajustar_spinodoide
from spinpy.grf import campo_grf, level_set
from spinpy.resistencia import (EPS_CRITICA, FRAC_CRITICA, capa_superficie,
                                criterio_pistoia, ensayo_compresion,
                                estadisticos_vm)

BLOQUE = "22 Modelos del informe"
REF = "informe_modelos.py frente al codigo que describe"
TOL = 1e-9


def _cerca(a, b, tol=TOL):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return bool(np.all(np.abs(a - b) <= tol * np.maximum(1.0, np.abs(b))))


def _q_manual(x, p, tipo):
    """Ecuaciones 4 y 5 del informe, literalmente."""
    x = np.sort(np.asarray(x, float))
    n = x.size
    h = (n - 1) * p + 1 if tipo == 7 else n * p + 0.5
    h = min(max(h, 1.0), float(n))
    lo = int(np.floor(h))
    hi = min(lo + 1, n)
    return x[lo - 1] + (h - lo) * (x[hi - 1] - x[lo - 1])


@pytest.fixture(scope="module")
def poroso():
    BW = generar_mascara(16, 6 * np.pi, 200, [30, 30, 60], 0.45, seed=5)[0]
    sp = np.full(3, 0.1)
    r = ensayo_compresion(BW, sp, E_s=20e9, nu_s=0.3, apoyo="deslizante")
    assert r["ok"], r["msg"]
    return BW, sp, r


def test_cuantiles(registro):
    rng = np.random.default_rng(1)
    x = rng.lognormal(size=257)
    ok = all(_cerca(_q_manual(x, p, 7), np.percentile(x, 100 * p))
             and _cerca(_q_manual(x, p, 5),
                        np.percentile(x, 100 * p, method="hazen"))
             for p in (0.5, 0.9, 0.98, 0.99))
    registro.anotar(BLOQUE, "Q7 y Q5 del informe == numpy linear y hazen",
                    REF, None, None, "1e-9", "identidad", ok)
    assert ok


def test_tension_desde_fuerza(registro):
    BW = np.ones((4, 4, 4), bool)
    sp = np.full(3, 0.25)                         # A = 1 mm^2
    r = ensayo_compresion(BW, sp, carga_N=100.0, unidad="mm")
    ok = r["ok"] and _cerca(r["sigma_app"], 1e6 * 100.0 / 1.0)
    registro.anotar(BLOQUE, "sigma0 = 1e6 F / A[mm2]", REF, 1e8,
                    r["sigma_app"], "1e-9", "ecuacion 3", ok)
    assert ok


def test_parfitt_y_da(registro):
    BW = generar_mascara(40, 8 * np.pi, 300, [15, 15, 60], 0.35, seed=2)[0]
    sp = np.full(3, 0.05)
    m = morfometria(BW, sp)
    lam = np.asarray(m["eigenvalues"], float)
    ok = (_cerca(m["BVTV"], BW.mean())
          and _cerca(m["TbTh"], 2 * m["BV"] / m["BS"])
          and _cerca(m["TbSp"], m["TbTh"] * (m["TV"] / m["BV"] - 1))
          and _cerca(m["TbN"], m["BVTV"] / m["TbTh"])
          and _cerca(m["DA"], np.sqrt(lam[2] / lam[0]))
          and _cerca(m["DA2"], np.sqrt(lam[1] / lam[0])))
    registro.anotar(BLOQUE, "Parfitt, BV/TV, DA y DA2 (ecuaciones 6-8)", REF,
                    None, None, "1e-9", "identidad", ok)
    assert ok


def test_umbral_espinodoide(registro):
    rho = 0.37
    phi0 = np.sqrt(2) * erfinv(2 * rho - 1)
    GRF = campo_grf(20, 7 * np.pi, 150, [20, 40, 90], seed=9)[0]
    BW = generar_mascara(20, 7 * np.pi, 150, [20, 40, 90], rho, seed=9)[0]
    ok = _cerca(level_set(rho), phi0) and np.array_equal(BW, GRF <= phi0)
    registro.anotar(BLOQUE, "phi0 y fase solida del espinodoide (ec. 12)",
                    REF, phi0, level_set(rho), "1e-9", "identidad", ok)
    assert ok


def test_densidad_dual(registro):
    n, rho = 24, 0.3
    BW = generar_dual_lattice(n, 3.0, rho, seed=4)[0]
    k = round(rho * n ** 3)
    # Solo los empates exactos de distancia pueden sumar voxeles.
    ok = k <= int(BW.sum()) <= k + 50
    registro.anotar(BLOQUE, "dual-lattice: k = round(rho n^3) (ec. 13)", REF,
                    k, int(BW.sum()), "+50 voxeles (empates)", "cuantil", ok)
    assert ok


def test_ondas_del_ajuste(registro):
    fuente = inspect.getsource(ajustar_spinodoide)
    ok = ("np.clip([10.0, 15.0, 20.0], 8, 25)" in fuente
          and "np.round(np.linspace(8, 25, 6))" in fuente
          and list(IM.ONDAS_COMPLETO) == [8, 11, 15, 18, 22, 25]
          and list(IM.ONDAS_RAPIDO) == [10, 15, 20])
    registro.anotar(BLOQUE, "numeros de onda del texto == los del ajuste",
                    REF, None, None, "exacto", "fuente de fit.py", ok)
    assert ok


def test_remuestreo(registro):
    rng = np.random.default_rng(3)
    BW = rng.random((30, 30, 30)) < 0.4
    sp = np.full(3, 0.02)
    bw, spr = remuestrear_bw(BW, sp, 12)
    n, n2 = 30, bw.shape[0]
    idx = [int(round(i * (n - 1) / (n2 - 1))) for i in range(n2)]
    ok = (np.array_equal(bw, BW[np.ix_(idx, idx, idx)])
          and _cerca(spr, sp * n / n2))
    registro.anotar(BLOQUE, "remuestreo i = round(i'(n-1)/(n'-1)), h n/n'",
                    REF, None, None, "exacto", "identidad", ok)
    assert ok


def test_tensor_macizo_y_constantes(registro):
    E, nu = 20e9, 0.3
    C, info = homogeneizar(np.ones((4, 4, 4), bool), E, nu, vox_size=0.1)
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[np.arange(3), np.arange(3)] = lam + 2 * mu
    D[np.arange(3, 6), np.arange(3, 6)] = mu
    S = np.linalg.inv(C)
    ec = constantes_ingenieria(C)
    # Adimensional: los ceros fuera de la diagonal salen ~1e-6 Pa por redondeo
    # de E_s = 2e10, y en Pa absolutos ninguna tolerancia relativa los admite.
    ok = (info["ok"] and _cerca(C / E, D / E, 1e-8)
          and _cerca(ec["Ex"], 1 / S[0, 0]) and _cerca(ec["Gyz"], 1 / S[3, 3])
          and _cerca(ec["nu_xy"], -S[0, 1] / S[0, 0]))
    registro.anotar(BLOQUE, "C de un solido macizo == D (ec. 2); E, G, nu "
                    "de S (ec. 19)", REF, E, ec["Ex"], "1e-8", "identidad", ok)
    assert ok


def test_von_mises_y_deformacion_efectiva(registro, poroso):
    _BW, _sp, r = poroso
    s, e = np.asarray(r["sig"]), np.asarray(r["eps"])
    formula = np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2
                             + (s[:, 2] - s[:, 0]) ** 2)
                      + 3 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2))
    # sqrt(3 J2) desde el tensor desviador, sin pasar por la formula
    T = np.zeros((len(s), 3, 3))
    T[:, 0, 0], T[:, 1, 1], T[:, 2, 2] = s[:, 0], s[:, 1], s[:, 2]
    T[:, 1, 2] = T[:, 2, 1] = s[:, 3]
    T[:, 0, 2] = T[:, 2, 0] = s[:, 4]
    T[:, 0, 1] = T[:, 1, 0] = s[:, 5]
    dev = T - np.trace(T, axis1=1, axis2=2)[:, None, None] / 3 * np.eye(3)
    j2 = 0.5 * np.einsum("nij,nij->n", dev, dev)
    U = 0.5 * np.einsum("ij,ij->i", s, e)
    ok = (_cerca(r["vm"], formula, 1e-9) and _cerca(formula, np.sqrt(3 * j2),
                                                     1e-9)
          and _cerca(r["eps_eff"], np.sqrt(2 * U / 20e9), 1e-9))
    registro.anotar(BLOQUE, "sigma_vM (ec. 22) == sqrt(3 J2) == codigo; "
                    "eps_eff (ec. 21)", REF, None, None, "1e-9", "identidad",
                    ok)
    assert ok


def test_bloque_uniaxial(registro):
    r = ensayo_compresion(np.ones((6, 6, 6), bool), np.full(3, 0.1),
                          E_s=20e9, nu_s=0.3, apoyo="deslizante")
    ok = (r["ok"] and _cerca(r["E_app"], 20e9, 1e-9)
          and _cerca(r["vm_solido"], r["sigma_app"], 1e-9)
          and _cerca(r["eps_eff_solido"], r["sigma_app"] / 20e9, 1e-9))
    registro.anotar(BLOQUE, "bloque macizo: E_app = E_s, sigma_vM = sigma0 "
                    "(ec. 23)", REF, 20e9, r["E_app"], "1e-9",
                    "solucion exacta", ok)
    assert ok


def test_capa_y_percentil_superficie(registro, poroso):
    BW, _sp, r = poroso
    # Ecuacion 24 escrita a mano: vecino por cara en el vacio, B = 1 fuera.
    Bp = np.pad(BW, 1, constant_values=True)
    vac = np.zeros(BW.shape, bool)
    for eje in range(3):
        for d in (1, -1):
            vac |= ~np.roll(Bp, d, axis=eje)[1:-1, 1:-1, 1:-1]
    capa = BW & vac
    campo = np.asarray(r["campo_vm"], float)
    dentro = np.isfinite(campo)
    vs = np.asarray(r["vm_solido"], float)[
        np.asarray(r["superficie_solido"], bool)]
    est = estadisticos_vm(r)
    ok = (np.array_equal(capa_superficie(BW), capa)
          and _cerca(est["vm_p99_superficie"], _q_manual(vs, 0.99, 5))
          and est["vm_n_superficie"] == vs.size and dentro.any())
    registro.anotar(BLOQUE, "capa superficial (ec. 24) y p99 con Q5 "
                    "(ec. 25)", REF, _q_manual(vs, 0.99, 5),
                    est["vm_p99_superficie"], "1e-9", "identidad", ok)
    assert ok


def test_pistoia(registro, poroso):
    _BW, _sp, r = poroso
    p = criterio_pistoia(r)
    e = np.asarray(r["eps_eff_solido"], float)
    e = e[np.isfinite(e)]
    k = EPS_CRITICA / _q_manual(e, 1 - FRAC_CRITICA, 7)
    ok = (_cerca(p["factor"], k) and _cerca(p["sigma_fallo"],
                                            k * r["sigma_app"])
          and _cerca(p["F_fallo"], k * r["F_total"]))
    registro.anotar(BLOQUE, "Pistoia con Q7 (ec. 26)", REF, k, p["factor"],
                    "1e-9", "identidad", ok)
    assert ok


def test_ecuaciones_componen(registro):
    md2pdf = pytest.importorskip("spinpy.md2pdf")
    doc = _doc_completo()
    malos, total = [], 0
    for idioma in ("es", "en"):
        md = informe.informe_markdown(doc, informe.comprobar(doc), {},
                                      idioma, "reproduccion.json",
                                      [("figuras/x.png", "x")])
        for b in re.findall(r"```math\n(.*?)\n```", md, re.S):
            total += 1
            tex, _num = md2pdf.texto_ecuacion(b.split("\n"))
            if md2pdf.imagen_ecuacion(tex) is None:
                malos.append(tex[:60])
    ok = total >= 40 and not malos
    registro.anotar(BLOQUE, "todas las ecuaciones componen (ES y EN)", REF,
                    0, len(malos), "exacto", f"{total} ecuaciones", ok,
                    nota="; ".join(malos))
    assert ok


def _doc_completo():
    """Documento con todos los calculos, para que salgan todas las secciones."""
    m = {"BVTV": 0.3, "TbTh": 0.19, "DA": 1.4, "BSBV": 10.0, "ConnD": 3.0,
         "PoDm": 0.5, "EF": 0.1, "modo": "voxel"}
    par = {"familia": "spinodoide", "densidad": 0.3, "wave_number_pi": 15.0,
           "wave_number_rad": 15 * np.pi}
    ens = {"resolucion": 40, "E_s_Pa": 20e9, "nu_s": 0.3,
           "apoyo": "deslizante", "por_estructura": {}}
    return {"morfometria_voi": m, "morfometria_spin": m,
            "morfometria_dual": m,
            "resultados": {
                "morfometria": {"modo": "voxel"},
                "ajuste": {"parametros": par, "objetivo": {},
                           "incertidumbre": {"K": 5}},
                "ajuste_dual": {"parametros": dict(par, familia="dual-lattice"),
                                "objetivo": {}},
                "ajuste_comparado": {"elegido": "completo"},
                "elastico": dict(ens), "resistencia": dict(ens),
                "analisis_comparado": {"resolucion": 40, "E_s_Pa": 18e9,
                                       "nu_s": 0.3, "carga_N": 100.0,
                                       "apoyo": "empotrado",
                                       "protocolo": "Tapia et al.",
                                       "por_estructura": {}}}}
