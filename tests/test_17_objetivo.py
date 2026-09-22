"""
Bloque 17 — Funcion objetivo, tensores, region degenerada, incertidumbre y
alineacion impuesta.

QUE SE VERIFICA, Y CONTRA QUE
  identidad    con `pesos=None` —y con los pesos base escritos a mano— el
               error es bit a bit el de antes, y dentro de 1e-12 del de MATLAB
               en los 1173 pares de `validar_error.py` (se salta si faltan).
  C6           un termino opcional ausente se salta y renormaliza.
  distancias   W1 entre dos uniformes desplazadas = desplazamiento (exacto);
               escalar una distribucion por s da el termino (s - 1)^2;
               Hellinger 0 si iguales, 1 si disjuntas, y entre N(0,1) y N(1,1)
               H^2 = 1 - exp(-1/8) = 0.1175 (tol 0.01: 32 clases).
  tensores     d_LE(aC, C) = sqrt(6)|ln a| (1e-12); invariante si se giran los
               dos; un laminado de Backus girado conserva E1..E3 y G (1e-9) y
               sale transversal; isotropo, cubico y ortotropo a mano.
  degenerada   (60,60,60), (70,70,30), (55,55,55), (15,45,90) generan la MISMA
               mascara que (90,90,90); (15,15,60) no; en equitativo no hay
               degeneracion; junto al borde, fraccion sin cubrir > 0 y pequena.
  incertidumbre K < 5 se rechaza; la seleccion robusta elige el de menor
               media + k sd con errores fabricados cuya respuesta se conoce.
  alineacion   un dual-lattice girado 30/50/20 grados como VOI: L2 deja el eje
               a <= 10 grados, nunca peor que la rotacion clasica; y un
               candidato PLANO se alinea por su normal.
"""
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import norm

from spinpy import incertidumbre
from spinpy.dual_lattice import generar_dual_lattice
from spinpy.elastic import (backus_laminado, clase_anisotropia,
                            constantes_principales, distancia_log_euclidea,
                            en_ejes)
from spinpy.error import (CAMPOS, PESOS_BASE, error_morfometrico,
                          hellinger_cuantiles, wasserstein_cuantiles)
from spinpy.fit import alinear_marco_fabrica
from spinpy.grf import euler_R, generar_mascara, region_degenerada
from spinpy.morphometry import metricas_forma, morfometria

BLOQUE = "17 Funcion objetivo y comparacion"
RAIZ = Path(__file__).resolve().parents[1]


def _niveles(n):
    return (np.arange(n) + 0.5) / n


# ---------------------------------------------------------------------------

def test_identidad_con_matlab(registro):
    pares = RAIZ / "resultados" / "pares_error.json"
    mat = RAIZ / "resultados" / "error_matlab.json"
    if not (pares.exists() and mat.exists()):
        pytest.skip("sin los pares de validar_error.py")
    d = json.loads(pares.read_text(encoding="utf-8"))
    m = json.loads(mat.read_text(encoding="utf-8"))

    def plano(c):
        while isinstance(c, list) and len(c) == 1 and isinstance(c[0], list):
            c = c[0]
        return c

    err_mat = np.array([np.nan if v is None else float(v)
                        for v in np.ravel(plano([m["errores"]]))])
    base = dict(zip(CAMPOS, PESOS_BASE))
    peor, identicos = 0.0, True
    for i, p in enumerate(d["pares"]):
        a = {c: v for c, v in zip(d["campos"], p["a"]) if v is not None}
        b = {c: v for c, v in zip(d["campos"], p["b"]) if v is not None}
        e0 = error_morfometrico(a, b)
        e1 = error_morfometrico(a, b, pesos=base)
        e2 = error_morfometrico(a, b, pesos={})
        identicos &= (e0 == e1 == e2)
        if np.isfinite(err_mat[i]):
            rel = abs(e0[0] - err_mat[i]) / (abs(err_mat[i]) or 1.0)
            peor = max(peor, rel)
    ok = identicos and peor <= 1e-12
    registro.anotar(BLOQUE, "error por omision = MATLAB, bit a bit",
                    "validar_error.py (1173 pares)", 0.0, peor, "<= 1e-12",
                    "pesos None == base explicitos == {}", ok)
    assert ok


def test_c6_terminos_opcionales(registro):
    a = {"BVTV": 0.3, "TbTh": 0.2, "PoDm": 0.6}
    b = {"BVTV": 0.33, "TbTh": 0.2}
    e_sin, n_sin = error_morfometrico(a, b)
    e_con, n_con = error_morfometrico(a, b, pesos={"PoDm": 5.0})
    b2 = dict(b, PoDm=0.3)
    e2, n2 = error_morfometrico(a, b2, pesos={"PoDm": 5.0})
    esperado = (3 * 0.1 ** 2 + 5 * 0.5 ** 2) / (3 + 1 + 5)
    ok = (e_sin == e_con and n_sin == n_con and n2 == n_sin + 1
          and abs(e2 - esperado) < 1e-15)
    with pytest.raises(ValueError):
        error_morfometrico(a, b, pesos={"PoDM": 1})
    with pytest.raises(ValueError):
        error_morfometrico(a, b, pesos={"DA2": 1})
    registro.anotar(BLOQUE, "termino opcional ausente se salta", "error.py C6",
                    esperado, e2, "1e-15", "renormaliza por el peso usado", ok)
    assert ok


def test_distancias_entre_distribuciones(registro):
    u = np.linspace(0.0, 1.0, 101)
    w = wasserstein_cuantiles(u, u + 0.25)
    q = norm.ppf(_niveles(1001), loc=2.0, scale=0.5)
    a = {"cuantiles": {"PoDm": q.tolist()}}
    b = {"cuantiles": {"PoDm": (1.3 * q).tolist()}}
    t_scale, _ = error_morfometrico(a, b, pesos={"dist_PoDm": 1.0,
                                                 "BVTV": 0, "DA": 0,
                                                 "BSBV": 0, "TbTh": 0,
                                                 "TbSp": 0, "TbN": 0,
                                                 "PoTot": 0})
    h0 = hellinger_cuantiles(q, q)
    h1 = hellinger_cuantiles(u, u + 5.0)
    hn = hellinger_cuantiles(norm.ppf(_niveles(1001)),
                             norm.ppf(_niveles(1001), loc=1.0))
    h2_exacta = 1.0 - np.exp(-1.0 / 8.0)
    # Hellinger de una distribucion consigo misma: sqrt(1 - bc) con bc = 1 -
    # redondeo, del orden de 1e-8. Cero en la practica, no bit a bit.
    ok = (abs(w - 0.25) < 1e-12 and abs(t_scale - 0.09) < 1e-3
          and h0 < 1e-6 and abs(h1 - 1.0) < 1e-9
          and abs(hn ** 2 - h2_exacta) < 0.01)
    registro.anotar(BLOQUE, "Hellinger N(0,1) vs N(1,1)",
                    "1 - exp(-1/8)", h2_exacta, hn ** 2, "0.01",
                    "32 clases desde cuantiles", ok,
                    nota=f"W1={w:.3g} termino escala={t_scale:.4g}")
    assert ok


# ---------------------------------------------------------------------------

def _iso(E=1.0, nu=0.3):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    C = np.zeros((6, 6))
    C[:3, :3] = lam
    C[[0, 1, 2], [0, 1, 2]] = lam + 2 * mu
    C[[3, 4, 5], [3, 4, 5]] = mu
    return C


def test_tensores(registro):
    C = backus_laminado([0.5, 0.5], [1.0, 0.3], [0.3, 0.3])
    C2 = backus_laminado([0.3, 0.7], [1.0, 0.1], [0.3, 0.25])
    R = euler_R(20, 35, -50)
    d_esc = distancia_log_euclidea(2.0 * C, C)
    d_rot = distancia_log_euclidea(en_ejes(C, R), en_ejes(C2, R))
    d = distancia_log_euclidea(C, C2)
    p = constantes_principales(C)
    pr = constantes_principales(en_ejes(C, R))
    claves = ("E1", "E2", "E3", "G23", "G13", "G12")
    dif = max(abs(p[k] - pr[k]) / abs(p[k]) for k in claves)
    cubica = _iso()
    cubica[[3, 4, 5], [3, 4, 5]] = 0.2
    orto = C.copy()
    orto[0, 0] *= 1.3
    clases = (clase_anisotropia(C)["clase"],
              clase_anisotropia(en_ejes(C, R))["clase"],
              clase_anisotropia(_iso())["clase"],
              clase_anisotropia(cubica)["clase"],
              clase_anisotropia(orto)["clase"])
    ok = (abs(d_esc - np.sqrt(6) * np.log(2)) < 1e-12
          and abs(d_rot - d) < 1e-9 and dif < 1e-9
          and clases == ("transversal", "transversal", "isotropa", "cubica",
                         "ortotropa"))
    registro.anotar(BLOQUE, "laminado girado: constantes principales",
                    "Backus 1962", 0.0, dif, "1e-9",
                    "d_LE(2C,C)=sqrt6 ln2; clases a mano", ok,
                    nota=str(clases))
    assert ok, (d_esc, d_rot, d, dif, clases)


def test_region_degenerada(registro):
    kw = dict(resolution=24, wave_number=10 * np.pi, num_waves=300, rho=0.35,
              seed=7)
    iso, _, _ = generar_mascara(thetas=(90, 90, 90), **kw)
    iguales = all(np.array_equal(generar_mascara(thetas=t, **kw)[0], iso)
                  for t in [(60, 60, 60), (70, 70, 30), (55, 55, 55),
                            (15, 45, 90)])
    distinta = not np.array_equal(generar_mascara(thetas=(15, 15, 60),
                                                  **kw)[0], iso)
    flags = (region_degenerada((70, 70, 30))["degenerado"],
             region_degenerada((15, 15, 60))["degenerado"],
             region_degenerada((60, 60, 60), "equitativo")["degenerado"])
    borde = region_degenerada((60, 60, 40), num_waves=700)
    ok = (iguales and distinta and flags == (True, False, False)
          and 0.0 < borde["fraccion_sin_cubrir"] < 0.01
          and borde["prob_realizacion_isotropa"] is not None)
    registro.anotar(BLOQUE, "region degenerada: misma mascara",
                    "sum cos^2 <= 1 o theta >= 90", None,
                    borde["fraccion_sin_cubrir"], "identica bit a bit",
                    "(70,70,30) incluido: mas ancha que 54.74 grados", ok,
                    nota=f"P(iso | 60,60,40, 700 ondas) = "
                         f"{borde['prob_realizacion_isotropa']:.3f}")
    assert ok


def test_metricas_forma_losa(registro):
    n, e = 40, 10
    BW = np.zeros((n, n, n), bool)
    BW[:, :, 15:15 + e] = True
    f = metricas_forma(BW, 1.0, tbth=True, poro=True, curvatura=True)
    q = np.asarray(f["cuantiles"]["TbTh"])
    ok = (abs(f["TbTh_local"] - e) < 1e-9 and f["TbTh_CV"] < 1e-9
          and len(q) == 101 and np.all(np.diff(q) >= 0)
          and "PoDm" in f["cuantiles"] and "H" in f["cuantiles"])
    registro.anotar(BLOQUE, "espesor local de una losa", "espesor.py (exacto)",
                    float(e), f["TbTh_local"], "1e-9", "CV = 0", ok)
    assert ok


# ---------------------------------------------------------------------------

def test_incertidumbre_y_seleccion_robusta(registro):
    with pytest.raises(ValueError):
        incertidumbre.validar_replicas(3)
    assert incertidumbre.validar_replicas(0) == 0
    tabla = {"A": [0.002, 0.018, 0.004, 0.016, 0.010],   # media 0.010, sd 0.0071
             "B": [0.011, 0.013, 0.012, 0.012, 0.012]}   # media 0.012, sd 0.0007
    sems = incertidumbre.semillas(100, 5)
    fins = [{"id": "A"}, {"id": "B"}]

    def generar_de(t, s):
        return (t["id"], sems.index(s))

    def medir(x):
        return {"e": tabla[x[0]][x[1]]}

    def error(_a, m):
        return m["e"], 1

    i1, t1 = incertidumbre.seleccion_robusta(fins, generar_de, medir, error,
                                             {}, 100, 5, k=1.0)
    i0, _ = incertidumbre.seleccion_robusta(fins, generar_de, medir, error,
                                            {}, 100, 5, k=0.0)
    ok = (i1 == 1 and i0 == 0 and 100 not in sems
          and abs(t1[0]["error"]["media"] - 0.010) < 1e-12)
    registro.anotar(BLOQUE, "seleccion robusta", "construccion", 1.0,
                    float(i1), "exacto", "k=1 elige B, k=0 elige A", ok)
    assert ok


def _eje(m, col):
    return np.asarray(m["eigenvectors"], float)[:, col]


def _ang(u, v):
    return float(np.degrees(np.arccos(min(1.0, abs(float(
        u @ v) / np.linalg.norm(u) / np.linalg.norm(v))))))


@pytest.mark.lento
@pytest.mark.parametrize("estir,col,nombre", [((1.0, 1.0, 2.5), 0, "axial"),
                                               ((2.5, 2.5, 1.0), 2, "plano")])
def test_alineacion_impuesta(registro, estir, col, nombre):
    n, sp = 48, np.full(3, 0.1)
    R0 = euler_R(30, 50, 20)
    VOI, _, _ = generar_dual_lattice(n, 5.0, 0.3, estiramiento=estir, R=R0,
                                     seed=1)
    m_voi = morfometria(VOI, sp)

    def generar(R):
        return generar_dual_lattice(n, 5.0, 0.3, estiramiento=estir, R=R,
                                    seed=2)[0]

    def medir(BW):
        return morfometria(BW, sp)

    R, info = alinear_marco_fabrica(m_voi, generar, medir, np.eye(3))
    m = medir(generar(R))
    ang = _ang(_eje(m, col), _eje(m_voi, col))
    ok = (ang <= 10.0 and info["angulo_despues_deg"]
          <= info["angulo_antes_deg"] + 1e-9)
    registro.anotar(BLOQUE, f"alineacion impuesta ({nombre})",
                    "construccion: VOI girado 30/50/20", 0.0, ang,
                    "<= 10 grados", info.get("eje_comparado", ""), ok,
                    nota=f"historial {info.get('historial_deg')}")
    assert ok, info
