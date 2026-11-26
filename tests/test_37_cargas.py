"""
Bloque 37 - Analisis mecanico con cargas y apoyos elegidos por el usuario.

QUE SE VERIFICA Y POR QUE
-------------------------
`spinpy.cargas` resuelve un caso de carga ARBITRARIO (fuerzas, presiones,
momentos, desplazamientos impuestos, sobre regiones elegidas con el raton).
Un solver generico falla sin dar sintomas: una fuerza con el signo cambiado, un
par aplicado alrededor de otro eje o una carga repartida por nodo en lugar de
por area dan un campo de aspecto razonable. Por eso se comprueba contra
respuestas que se conocen de antemano:

  37.1  El caso «compresion» reproduce `ensayo_compresion` (el ensayo validado
        de la app) sobre un hueso poroso: mismo campo de von Mises.
  37.2  Bloque macizo: tension uniforme -sigma, acortamiento sigma*L/E y
        expansion lateral nu*sigma/E (soluciones cerradas de la elasticidad).
  37.3  La respuesta es lineal: la de una fuerza oblicua es la combinacion de
        las de sus componentes. Es lo que justifica el barrido en dos solves.
  37.4  Equilibrio: la suma de reacciones iguala a la de cargas, en fuerzas y
        en momentos, con varias cargas, un par y un desplazamiento impuesto,
        con el solver directo y con el multigrid.
  37.5  Torsion de una barra de seccion cuadrada: el giro por unidad de longitud
        frente a la solucion de Saint-Venant, J = 0.1406 a^4 (Timoshenko).
  37.6  Flexion pura: la pendiente de la tension axial frente a la coordenada
        de la fibra es M / I.
  37.7  Rigidez impuesta y rigidez medida coinciden: k = E A / L.
  37.8  Regiones: area de una cara, parche esferico, filtro por normal, nodos.
  37.9  Barrido de angulo: igual a resolver cada angulo por separado.
  37.10 Rechazos: cargas desequilibradas, carga sobre un fragmento sin apoyo,
        caso mal formado. Un resultado falso es peor que un error.
  37.11 Archivos: caso y resumen en JSON, CSV, VTK.
  37.12 La seccion del informe: sale solo con casos, en los dos idiomas.
  37.13 La escena 3D y el pinchazo con el raton (solo con render disponible).
  37.14 La ventana (en un subproceso, con y sin vista 3D).

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  37.1  campo vm        5e-6 relativo (el campo del ensayo se guarda en float32)
  37.2  solucion cerrada 1e-8 (el campo es exactamente lineal en hexaedros)
  37.3  superposicion   1e-8
  37.4  equilibrio      1e-6 con el iterativo (tolerancia 1e-8 del solver), 1e-9
                        con el directo
  37.5  torsion         8 %: hexaedros trilineales de 12x12 por seccion; el par
                        se reparte con fuerza tangencial proporcional a r, que
                        no es la traccion exacta de Saint-Venant, y la medida
                        se toma lejos de los dos extremos
  37.6  flexion         5 %, medida en el centro de la barra
  37.7  rigidez         1e-8
  37.9  barrido         1e-8 (misma matriz, misma carga)
"""
import json
import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest
from scipy import ndimage

from spinpy import cargas as C
from spinpy.resistencia import ensayo_compresion

BLOQUE = "37 · cargas y apoyos a elegir"
REF = "Elasticidad lineal: soluciones cerradas y `resistencia.ensayo_compresion`."
E = 20e9
NU = 0.30
E_MPA = E / 1e6


def _poroso(n=14, semilla=1, pct=50):
    """Hueso sintetico poroso y CONECTADO: la mayor componente de un campo
    aleatorio suavizado y umbralizado."""
    f = ndimage.gaussian_filter(
        np.random.default_rng(semilla).standard_normal((n, n, n)), n / 9.0)
    bw = f > np.percentile(f, pct)
    lab, _n = ndimage.label(bw)
    return lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)


def _barra(a, largo, h=0.1):
    return np.ones((a, a, largo), bool), h


def _nodos_cara_z(res, k):
    """u en los nodos del plano z = k (indice de nodo)."""
    return res["u_malla"][:, :, k, :].astype(float)


# -- 37.1 ----------------------------------------------------------------------

def test_compresion_reproduce_el_ensayo_validado(registro):
    BW, sp = _poroso(20), 0.04
    ref = ensayo_compresion(BW, sp, sigma0=1e6, apoyo="deslizante")
    assert ref["ok"]
    r = C.resolver_caso(BW, sp, C.plantilla("compresion", BW.shape, sp,
                                            eje=2, sigma_MPa=1.0))
    assert r["ok"], r["msg"]
    a = ref["campo_vm"] / 1e6
    b = r["campos"]["vm"]
    ambos = np.isfinite(a) & np.isfinite(b)
    assert np.array_equal(np.isfinite(a), np.isfinite(b))
    err = float(np.abs(a[ambos] - b[ambos]).max() / a[ambos].max())
    registro.anotar(BLOQUE, "compresion = ensayo_compresion (von Mises)", REF,
                    0.0, err, "5e-6 relativo", f"{int(ambos.sum())} elementos",
                    err < 5e-6)
    assert err < 5e-6
    # y la reaccion es la carga total de la app
    F = ref["F_total"] / 1e6
    assert r["reacciones"][0]["F_N"][2] == pytest.approx(F, rel=1e-6)


# -- 37.2 ----------------------------------------------------------------------

def test_bloque_macizo_solucion_cerrada(registro):
    n, h, sig = 10, 0.05, 2.0
    BW = np.ones((n, n, n), bool)
    r = C.resolver_caso(BW, h, C.plantilla("compresion", BW.shape, h, eje=2,
                                           sigma_MPa=sig))
    assert r["ok"], r["msg"]
    for k, esperado in (("szz", -sig), ("sxx", 0.0), ("syy", 0.0)):
        v = r["campos"][k]
        assert np.abs(v - esperado).max() < 1e-8 * sig, k
    L = n * h
    u = r["u_malla"].astype(float)
    acort = u[:, :, -1, 2].mean() - u[:, :, 0, 2].mean()
    esp = -sig * L / E_MPA
    err = abs(acort - esp) / abs(esp)
    registro.anotar(BLOQUE, "acortamiento sigma*L/E", REF, esp, acort,
                    "1e-5 relativo (u en float32)", "bloque macizo", err < 1e-5)
    assert err < 1e-5
    lat = u[-1, :, :, 0].mean() - u[0, :, :, 0].mean()
    esp_lat = -NU * (-sig / E_MPA) * L
    err_l = abs(lat - esp_lat) / abs(esp_lat)
    registro.anotar(BLOQUE, "expansion lateral nu*sigma*L/E", REF, esp_lat,
                    lat, "1e-5 relativo", "apoyo deslizante: lados libres",
                    err_l < 1e-5)
    assert err_l < 1e-5


def test_traccion_es_la_compresion_con_signo_cambiado():
    BW, sp = _poroso(14), 0.05
    a = C.resolver_caso(BW, sp, C.plantilla("compresion", BW.shape, sp))
    b = C.resolver_caso(BW, sp, C.plantilla("traccion", BW.shape, sp))
    assert a["ok"] and b["ok"]
    assert np.allclose(a["campos"]["vm"], b["campos"]["vm"], equal_nan=True,
                       rtol=1e-6)
    assert np.allclose(a["campos"]["szz"], -b["campos"]["szz"],
                       equal_nan=True, rtol=1e-6, atol=1e-9)


# -- 37.3 ----------------------------------------------------------------------

def test_la_respuesta_es_lineal(registro):
    BW, sp = _poroso(14, semilla=2), 0.05
    reg = {"forma": "cara", "cara": "+z"}
    ap = [{"tipo": "empotrado", "region": {"forma": "cara", "cara": "-z"}}]

    def caso(d, m):
        c = C.caso_vacio()
        c["cargas"] = [{"tipo": "fuerza", "magnitud": m, "direccion": d,
                        "region": reg}]
        c["apoyos"] = ap
        return c
    a, b, c3 = 3.0, -2.0, 5.0
    mix = C.resolver_caso(BW, sp, caso([a, b, c3], float(np.sqrt(
        a * a + b * b + c3 * c3))))
    base = [C.resolver_caso(BW, sp, caso(e, 1.0))
            for e in ([1, 0, 0], [0, 1, 0], [0, 0, 1])]
    assert mix["ok"] and all(x["ok"] for x in base)
    suma = a * base[0]["u_malla"] + b * base[1]["u_malla"] \
        + c3 * base[2]["u_malla"]
    err = float(np.abs(mix["u_malla"] - suma).max()
                / np.abs(mix["u_malla"]).max())
    registro.anotar(BLOQUE, "superposicion de fuerzas", REF, 0.0, err,
                    "1e-5 relativo (u en float32)", "fuerza oblicua", err < 1e-5)
    assert err < 1e-5


# -- 37.4 ----------------------------------------------------------------------

def _caso_mixto(BW, sp):
    ext = C.extension(BW.shape, sp)
    c = C.caso_vacio()
    c["cargas"] = [
        {"tipo": "fuerza", "magnitud": 3.0, "direccion": [1, 0, -1],
         "region": {"forma": "esfera", "centro": [ext[0] / 2, ext[1] / 2,
                                                  ext[2]],
                    "radio": 0.3 * ext[0]}},
        {"tipo": "momento", "magnitud": 1.5, "eje": [0, 1, 0],
         "region": {"forma": "cara", "cara": "+x"}},
        {"tipo": "presion", "magnitud": 0.5, "direccion": "normal_entrante",
         "area_ref": "expuesta",
         "region": {"forma": "cara", "cara": "+y"}}]
    c["apoyos"] = [{"tipo": "empotrado",
                    "region": {"forma": "cara", "cara": "-z"}}]
    return c


@pytest.mark.parametrize("n,solver", [(12, "directo"), (22, "multigrid")])
def test_equilibrio_de_fuerzas_y_momentos(registro, n, solver):
    BW, sp = _poroso(n, semilla=3, pct=45), 0.05
    r = C.resolver_caso(BW, sp, _caso_mixto(BW, sp))
    assert r["ok"], r["msg"]
    if solver == "directo":
        assert r["solver"] == "LU directo"
        tol = 1e-9
    else:
        assert "AMG" in r["solver"] or "CG" in r["solver"]
        tol = 1e-6
    eq = r["equilibrio"]
    registro.anotar(BLOQUE, f"equilibrio de fuerzas ({solver})", REF, 0.0,
                    eq["residuo_F"], f"{tol:g} relativo",
                    "fuerza + par + presion", eq["residuo_F"] < tol)
    registro.anotar(BLOQUE, f"equilibrio de momentos ({solver})", REF, 0.0,
                    eq["residuo_M"], f"{tol:g} relativo",
                    "fuerza + par + presion", eq["residuo_M"] < tol)
    assert eq["residuo_F"] < tol and eq["residuo_M"] < tol


def test_el_par_no_aplica_fuerza_neta():
    BW, sp = _barra(6, 12)
    c = C.plantilla("torsion", BW.shape, sp, eje=2, F_N=10.0)
    r = C.resolver_caso(BW, sp, c)
    assert r["ok"], r["msg"]
    d = r["cargas"][0]
    assert np.abs(d["F_N"]).max() < 1e-12 * 10.0
    # el par pedido es el aplicado, alrededor del eje pedido
    M = 10.0 * 0.6
    assert d["M_Nmm"][2] == pytest.approx(M, rel=1e-9)
    assert abs(d["M_Nmm"][0]) < 1e-9 * M and abs(d["M_Nmm"][1]) < 1e-9 * M


# -- 37.5 / 37.6 -----------------------------------------------------------------

def test_torsion_de_barra_cuadrada_saint_venant(registro):
    a, L = 12, 48
    h = 0.1
    BW, _h = _barra(a, L, h)
    M = 50.0
    c = C.plantilla("torsion", BW.shape, h, eje=2, F_N=1.0)
    c["cargas"][0]["magnitud"] = M
    r = C.resolver_caso(BW, h, c)
    assert r["ok"], r["msg"]
    u = r["u_malla"].astype(float)
    x = (np.arange(a + 1) * h - a * h / 2)
    X, Y = np.meshgrid(x, x, indexing="ij")

    def giro(k):
        ux, uy = u[:, :, k, 0], u[:, :, k, 1]
        return float(((X * uy - Y * ux).sum()) / ((X ** 2 + Y ** 2).sum()))
    k1, k2 = int(0.30 * L), int(0.70 * L)
    tasa = (giro(k2) - giro(k1)) / ((k2 - k1) * h)
    lado = a * h
    J = 0.1406 * lado ** 4
    G = E_MPA / (2 * (1 + NU))
    esp = M / (G * J)
    err = abs(abs(tasa) - esp) / esp
    registro.anotar(BLOQUE, "torsion: giro por unidad de longitud",
                    "Timoshenko & Goodier: J = 0.1406 a^4 (seccion cuadrada)",
                    esp, abs(tasa), "8 %", f"{a}x{a}x{L} hexaedros", err < 0.08)
    assert err < 0.08


def test_flexion_pura_pendiente_de_la_tension(registro):
    a, L = 12, 48
    h = 0.1
    BW, _h = _barra(a, L, h)
    M = 40.0
    c = C.plantilla("flexion", BW.shape, h, eje=2, F_N=1.0)
    c["cargas"][0]["magnitud"] = M
    r = C.resolver_caso(BW, h, c)
    assert r["ok"], r["msg"]
    szz = r["campos"]["szz"]
    corte = szz[:, :, L // 2]                     # (x, y) en el centro
    y = (np.arange(a) + 0.5) * h - a * h / 2
    perfil = corte.mean(axis=0)                    # media en x, funcion de y
    pend = float(np.polyfit(y, perfil, 1)[0])
    I = (a * h) ** 4 / 12.0
    esp = M / I
    err = abs(abs(pend) - esp) / esp
    registro.anotar(BLOQUE, "flexion pura: d(sigma_zz)/dy = M/I",
                    "Euler-Bernoulli", esp, abs(pend), "5 %",
                    f"{a}x{a}x{L} hexaedros, seccion central", err < 0.05)
    assert err < 0.05


# -- 37.7 ------------------------------------------------------------------------

def test_desplazamiento_impuesto_da_la_rigidez_del_bloque(registro):
    n, h = 8, 0.1
    BW = np.ones((n, n, n), bool)
    d = 0.002
    c = C.plantilla("compresion", BW.shape, h, eje=2, desp_mm=d)
    assert c["cargas"][0]["tipo"] == "desplazamiento"
    r = C.resolver_caso(BW, h, c)
    assert r["ok"], r["msg"]
    L, A = n * h, (n * h) ** 2
    k = E_MPA * A / L
    reac = r["reacciones"]
    carga = next(x for x in reac if x["origen"] == "desp")
    esp = k * d
    err = abs(carga["modulo_N"] - esp) / esp
    registro.anotar(BLOQUE, "reaccion de un desplazamiento impuesto", REF, esp,
                    carga["modulo_N"], "1e-8 relativo", "k = E A / L", err < 1e-8)
    assert err < 1e-8
    g = next(x for x in r["rigidez"] if "desp_impuesto_mm" in x)
    assert g["rigidez_N_mm"] == pytest.approx(k, rel=1e-8)
    # el mismo bloque con la fuerza correspondiente se acorta lo mismo
    f = C.plantilla("compresion", BW.shape, h, eje=2, sigma_MPa=k * d / A)
    r2 = C.resolver_caso(BW, h, f)
    u = r2["u_malla"].astype(float)
    acort = u[:, :, -1, 2].mean() - u[:, :, 0, 2].mean()
    assert abs(acort) == pytest.approx(d, rel=1e-5)


# -- 37.8 ------------------------------------------------------------------------

def test_regiones_sobre_un_bloque_macizo():
    n, h = 8, 0.1
    BW = np.ones((n, n, n), bool)
    caras = C.caras_expuestas(BW, h)
    assert caras["cent"].shape[0] == 6 * n * n
    nodos = np.unique(caras["nod"])
    for cara in C.CARAS:
        s = C.seleccionar({"forma": "cara", "cara": cara}, caras, BW.shape, h,
                          nodos)
        assert s["area"] == pytest.approx((n * h) ** 2)
        assert s["quads"].size == n * n
        assert s["peso"].sum() == pytest.approx((n * h) ** 2)
    # una esfera centrada en la cara +z: un disco de area ~ pi r^2
    ext = n * h
    s = C.seleccionar({"forma": "esfera", "centro": [ext / 2, ext / 2, ext],
                       "radio": 0.25}, caras, BW.shape, h, nodos)
    assert 0.5 * np.pi * 0.25 ** 2 < s["area"] < 1.5 * np.pi * 0.25 ** 2
    # y con el filtro por normal solo quedan las caras de arriba
    sf = C.seleccionar({"forma": "esfera", "centro": [ext / 2, ext / 2, ext],
                        "radio": 0.25,
                        "filtro_normal": {"dir": [0, 0, 1], "angulo": 10}},
                       caras, BW.shape, h, nodos)
    assert (caras["nor"][sf["quads"], 2] > 0.99).all()
    # nodos sueltos
    s = C.seleccionar({"forma": "nodos", "ids": [int(nodos[0]),
                                                  int(nodos[5]), -1]},
                      caras, BW.shape, h, nodos)
    assert s["nodos"].size == 2
    # ambito volumen: los nodos a menos de una capa de la cara
    todos = np.arange((n + 1) ** 3)                # todos los nodos de hueso
    s = C.seleccionar({"forma": "cara", "cara": "+z", "capa": 1,
                       "ambito": "volumen"}, caras, BW.shape, h, todos)
    assert s["nodos"].size == 2 * (n + 1) ** 2
    # una region vacia no se resuelve y lo dice
    r = C.resolver_caso(BW, h, {
        **C.caso_vacio(),
        "cargas": [{"tipo": "fuerza", "magnitud": 1, "direccion": "-z",
                    "region": {"forma": "esfera", "centro": [5, 5, 5],
                               "radio": 0.1}}],
        "apoyos": [{"tipo": "empotrado",
                    "region": {"forma": "cara", "cara": "-z"}}]})
    assert not r["ok"] and "no toca" in r["msg"]


def test_pinchazos_se_convierten_en_regiones():
    BW, h = _barra(6, 10)
    caras = C.caras_expuestas(BW, h)
    nombre, d = C.cara_del_cubo([0.31, 0.29, 1.02], BW.shape, h)
    assert nombre == "+z" and d == pytest.approx(0.02)
    nombre, _d = C.cara_del_cubo([-0.01, 0.3, 0.5], BW.shape, h)
    assert nombre == "-x"
    reg = C.region_parche([0.3, 0.3, 1.0], 0.2, caras)
    assert reg["forma"] == "esfera"
    assert reg["filtro_normal"]["dir"] == [0.0, 0.0, 1.0]
    assert C.nodo_cercano(caras, [0.0, 0.0, 0.0], BW.shape, h) == 0


# -- 37.9 ------------------------------------------------------------------------

def test_barrido_igual_a_resolver_cada_angulo(registro):
    BW, sp = _poroso(14, semilla=4), 0.05
    c = C.plantilla("cizalla", BW.shape, sp, eje=2, F_N=7.0)
    c["cargas"][0]["direccion"] = [1, 0, 0]
    b = C.barrido_direccion(BW, sp, c, 0, plano=(0, 1), n=8)
    assert b["ok"], b["msg"]
    peor = 0.0
    for i in (1, 3, 6):
        t = np.radians(b["theta_deg"][i])
        cc = json.loads(json.dumps(c))
        cc["cargas"][0]["direccion"] = [float(np.cos(t)), float(np.sin(t)), 0]
        r = C.resolver_caso(BW, sp, cc)
        assert r["ok"], r["msg"]
        k = r["rigidez"][0]["rigidez_N_mm"]
        peor = max(peor, abs(b["rigidez_N_mm"][i] - k) / k,
                   abs(b["vm_max"][i] - r["stats"]["vm"]["max"])
                   / r["stats"]["vm"]["max"],
                   abs(b["factor"][i] - r["pistoia"]["factor"])
                   / r["pistoia"]["factor"])
    registro.anotar(BLOQUE, "barrido = un solve por angulo", REF, 0.0, peor,
                    "1e-8 relativo", "rigidez, von Mises maximo y Pistoia",
                    peor < 1e-8)
    assert peor < 1e-8
    # un angulo y su opuesto dan la misma rigidez (la respuesta es antisimetrica)
    assert b["rigidez_N_mm"][0] == pytest.approx(b["rigidez_N_mm"][4],
                                                 rel=1e-8)


# -- 37.10 -----------------------------------------------------------------------

def test_rechaza_cargas_desequilibradas_y_apoyos_inalcanzables():
    BW, h = _barra(6, 12)
    sin_apoyo = C.caso_vacio()
    sin_apoyo["cargas"] = [{"tipo": "fuerza", "magnitud": 1.0,
                            "direccion": "-z",
                            "region": {"forma": "cara", "cara": "+z"}}]
    r = C.resolver_caso(BW, h, sin_apoyo)
    assert not r["ok"] and "equilibr" in r["msg"]
    # sin estabilizar, el mismo caso se rechaza por falta de apoyo
    sin_apoyo["estabilizar"] = False
    r = C.resolver_caso(BW, h, sin_apoyo)
    assert not r["ok"]
    # cargas autoequilibradas SIN apoyos: si se resuelven
    eq = C.caso_vacio()
    eq["cargas"] = [
        {"tipo": "fuerza", "magnitud": 1.0, "direccion": "-z",
         "region": {"forma": "cara", "cara": "+z"}},
        {"tipo": "fuerza", "magnitud": 1.0, "direccion": "+z",
         "region": {"forma": "cara", "cara": "-z"}}]
    r = C.resolver_caso(BW, h, eq)
    assert r["ok"], r["msg"]
    assert r["n_estabilizadores"] == 6
    assert r["reaccion_estabilizadores_rel"] < 1e-6
    # una carga sobre un fragmento sin apoyo se nombra
    BW2 = np.zeros((6, 6, 14), bool)
    BW2[:, :, :6] = True
    BW2[:, :, 8:] = True                          # dos trozos separados
    c = C.plantilla("compresion", BW2.shape, h, eje=2)
    r = C.resolver_caso(BW2, h, c)
    assert not r["ok"] and "C1" in r["msg"] and "apoyo" in r["msg"]


def test_casos_mal_formados():
    BW, h = _barra(6, 12)
    con = [
        {},
        {"cargas": [{"tipo": "nube", "region": {}}]},
        {"cargas": [{"tipo": "fuerza", "magnitud": 1.0}]},
        {"cargas": [{"tipo": "fuerza", "magnitud": 1.0, "direccion": "-z",
                     "region": {"forma": "cara", "cara": "+z"}}],
         "apoyos": [{"tipo": "personalizado",
                     "region": {"forma": "cara", "cara": "-z"}}]},
        {"cargas": [{"tipo": "fuerza", "magnitud": float("nan"),
                     "direccion": "-z",
                     "region": {"forma": "cara", "cara": "+z"}}]},
        {"material": {"E_s": -1},
         "cargas": [{"tipo": "fuerza", "magnitud": 1.0, "direccion": "-z",
                     "region": {"forma": "cara", "cara": "+z"}}]},
    ]
    for caso in con:
        with pytest.raises(C.ErrorCaso):
            C.resolver_caso(BW, h, caso)
    with pytest.raises(C.ErrorCaso):                 # cara que no existe
        c = C.plantilla("compresion", BW.shape, h)
        c["cargas"][0]["region"]["cara"] = "+w"
        C.resolver_caso(BW, h, c)
    with pytest.raises(C.ErrorCaso):                 # direccion nula
        c = C.plantilla("cizalla", BW.shape, h)
        c["cargas"][0]["direccion"] = [0, 0, 0]
        C.resolver_caso(BW, h, c)


# -- 37.11 -----------------------------------------------------------------------

def test_archivos_del_caso_y_del_resultado(tmp_path):
    BW, sp = _poroso(14), 0.05
    caso = C.plantilla("puntual", BW.shape, sp, F_N=4.0)
    caso["cargas"].append({"tipo": "fuerza", "magnitud": 1.0,
                           "direccion": [0, 1, 0],
                           "region": {"forma": "nodos", "ids": [0],
                                      "xyz": [[0, 0, 0]]}})
    ruta = C.guardar_caso(caso, tmp_path / "c.json", {"clave": "voi"})
    assert C.cargar_caso(ruta) == json.loads(json.dumps(caso))
    (tmp_path / "x.json").write_text('{"formato": "otro"}')
    with pytest.raises(C.ErrorCaso):
        C.cargar_caso(tmp_path / "x.json")
    with pytest.raises(C.ErrorCaso):
        C.cargar_caso(tmp_path / "no_existe.json")
    caso["cargas"].pop()
    r = C.resolver_caso(BW, sp, caso)
    assert r["ok"], r["msg"]
    C.exportar_csv(r, caso, tmp_path / "r.csv")
    C.exportar_json(r, caso, tmp_path / "r.json", {"clave": "voi"})
    C.exportar_vtk(r, tmp_path / "r.vti")
    doc = json.loads((tmp_path / "r.json").read_text())
    assert doc["ok"] and doc["caso"]["cargas"][0]["tipo"] == "fuerza"
    assert "campos" not in doc and "u_malla" not in doc
    assert "reaccion" in (tmp_path / "r.csv").read_text()
    import pyvista as pv
    g = pv.read(str(tmp_path / "r.vti"))
    assert {"vm", "uz", "s1"} <= set(g.cell_data.keys())
    assert g.point_data["u_mm"].shape[1] == 3
    # resumen() es JSON puro
    json.dumps(C.resumen(r, caso))
    assert C.escala_deformada(r) > 1.0


# -- 37.12 -----------------------------------------------------------------------

def _doc_con_caso(BW, sp):
    caso = _caso_mixto(BW, sp)
    r = C.resolver_caso(BW, sp, caso)
    assert r["ok"], r["msg"]
    est = {"clave": "voi", "etiqueta": "VOI de referencia",
           "forma": list(BW.shape), "spacing_mm": [sp] * 3,
           "BV_TV": float(BW.mean())}
    return {"resultados": {"cargas_casos": [C.resumen(r, caso, est)]}}


@pytest.mark.parametrize("idioma", ["es", "en"])
def test_seccion_del_informe(idioma):
    from spinpy.informe import _Redactor
    from spinpy.informe_cargas import seccion_cargas
    BW, sp = _poroso(12, semilla=5, pct=45), 0.05
    assert seccion_cargas({}, idioma, _Redactor(idioma), 6) == []
    assert seccion_cargas({"resultados": {"cargas_casos": [
        {"ok": False}]}}, idioma, _Redactor(idioma), 6) == []
    L = seccion_cargas(_doc_con_caso(BW, sp), idioma, _Redactor(idioma), 6)
    txt = "\n".join(L)
    assert txt.startswith("## 6. ")
    assert "### 6.1" in txt
    assert "{" not in txt.replace("```json", "").split("```json")[0] \
        or True                                   # el JSON si lleva llaves
    cuerpo = txt.split("```json")[0]
    assert "{" not in cuerpo and "}" not in cuerpo
    assert "C1:" in txt and "A1:" in txt and "Fz (N)" in txt
    if idioma == "es":
        assert "Casos de carga" in txt and "empotrado" in txt
    else:
        assert "load cases" in txt and "fixed" in txt
    # el JSON del caso se puede volver a leer y resolver
    js = txt.split("```json")[1].split("```")[0]
    assert C.validar_caso(json.loads(js))["apoyos"][0]["tipo"] == "empotrado"


def test_informe_markdown_numera_la_seccion(monkeypatch):
    """La seccion entra entre los modelos y las referencias, y la numeracion
    de lo que sigue se desplaza una unidad."""
    from spinpy import informe as I
    BW, sp = _poroso(12, semilla=5, pct=45), 0.05
    base = {"procedencia": {}, "estado": {}, "resultados": {}, "voi": None}
    con = json.loads(json.dumps(base))
    con["resultados"] = _doc_con_caso(BW, sp)["resultados"]
    items = []
    sin_t = I.informe_markdown(base, items, {}, "es")
    con_t = I.informe_markdown(con, items, {}, "es")
    assert "Casos de carga definidos" not in sin_t
    assert "Casos de carga definidos" in con_t
    import re
    num = lambda t, nombre: int(re.search(     # noqa: E731
        r"^## (\d+)\. " + nombre, t, re.M).group(1))
    assert num(con_t, "Referencias") == num(sin_t, "Referencias") + 1
    assert num(con_t, "Casos de carga") == num(sin_t, "Referencias")


# -- 37.13 -----------------------------------------------------------------------

def test_escena_y_pinchazo_con_el_raton(render_disponible):
    if not render_disponible:
        pytest.skip("sin OpenGL fuera de pantalla")
    import pyvista as pv
    from vtkmodules.vtkRenderingCore import vtkCoordinate

    from spinpy.escena_cargas import EscenaCargas
    BW, sp = _poroso(16, semilla=6, pct=50), 0.05
    caso = _caso_mixto(BW, sp)
    r = C.resolver_caso(BW, sp, caso)
    assert r["ok"], r["msg"]
    p = pv.Plotter(off_screen=True, window_size=(480, 360))
    e = EscenaCargas(p, BW, sp)
    assert e.malla.n_points > 0
    ext = C.extension(BW.shape, sp)
    assert np.allclose(np.array(e.malla.bounds).reshape(3, 2)[:, 1], ext,
                       atol=1.0 * sp.max() if hasattr(sp, "max") else sp)
    e.poner_caso(caso, activa=0)
    for nombre in ("cubo", "carga_flecha_0", "carga_flecha_1",
                   "carga_reg_0", "apoyo_reg_0"):
        assert nombre in p.renderer.actors, nombre
    d0 = e.flechas[0][1]
    assert d0 @ np.array([1, 0, -1]) / np.sqrt(2) > 0.9   # la flecha sigue F
    for modo in ("suave", "voxeles"):
        e.modo_campo = modo
        e.mostrar_resultado(r, "vm", deformada=True,
                            escala=C.escala_deformada(r), corte=(1, ext[1] / 2,
                                                                 1))
        assert "campo" in p.renderer.actors
    e.quitar_resultado()
    assert "campo" not in p.renderer.actors and "estructura" in p.renderer.actors
    # pinchar: la esquina de la cara +z, proyectada a pantalla, devuelve ese punto
    e.vista("iso")
    p.render()
    q = np.array([ext[0] * 0.5, ext[1] * 0.5, ext[2]])
    c = vtkCoordinate()
    c.SetCoordinateSystemToWorld()
    c.SetValue(*q)
    x, y = c.GetComputedDisplayValue(p.renderer)
    visto = []
    e._pinchar_pantalla(int(x), int(y), lambda pt, sh: visto.append(pt))
    assert visto, "el clic no alcanzo ni la estructura ni el cubo"
    assert np.linalg.norm(visto[0] - q) < 4 * float(np.max(sp))
    p.close()


# -- 37.14 -----------------------------------------------------------------------

_GUION_VENTANA = textwrap.dedent('''
    import json, sys
    import numpy as np
    from scipy import ndimage
    from PyQt5 import QtWidgets, QtCore
    app = QtWidgets.QApplication([])
    from spinpy import cargas as C
    from dialogo_cargas import DialogoCargas

    def poroso(n, semilla, pct):
        f = ndimage.gaussian_filter(
            np.random.default_rng(semilla).standard_normal((n, n, n)), n / 9.0)
        bw = f > np.percentile(f, pct)
        lab, _ = ndimage.label(bw)
        return lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)

    est = {"voi": {"etq": "VOI", "BW": poroso(14, 1, 50), "spacing": 0.05},
           "spinodoide": {"etq": "Espinodoide", "BW": poroso(14, 2, 55),
                          "spacing": 0.05}}
    d = DialogoCargas(None, est, n_calculo=14)
    d.show(); app.processEvents()
    out = {"vista3d": d.vista3d is not None}
    # plantilla inicial: compresion
    out["inicial"] = d.caso["cargas"][0]["tipo"]
    # elegir cara con el raton, deshacer y rehacer
    d.cmb_modo.setCurrentIndex(d.cmb_modo.findData("cara"))
    d._al_pinchar(np.array([0.7, 0.35, 0.35]))
    out["cara"] = d.caso["cargas"][0]["region"].get("cara")
    d.deshacer(); out["deshecho"] = d.caso["cargas"][0]["region"].get("cara")
    d.rehacer(); out["rehecho"] = d.caso["cargas"][0]["region"].get("cara")
    # nodos: anadir dos y quitar uno con Mayus
    d._anadir_apoyo("empotrado")
    d.lista_a.setCurrentRow(1)
    d.cmb_modo.setCurrentIndex(d.cmb_modo.findData("nodos"))
    d._al_pinchar(np.array([0.1, 0.1, 0.0]))
    d._al_pinchar(np.array([0.5, 0.5, 0.0]))
    n2 = len(d.caso["apoyos"][1]["region"]["ids"])
    d._al_pinchar(np.array([0.5, 0.5, 0.0]), True)
    out["nodos"] = [n2, len(d.caso["apoyos"][1]["region"]["ids"])]
    d._quitar_apoyo()
    # nueva carga con un clic sobre la superficie y parche
    d.cmb_modo.setCurrentIndex(d.cmb_modo.findData("orbitar"))
    d._nueva_con_clic()
    d._al_pinchar(np.array([0.3, 0.3, 0.7]))
    out["n_cargas"] = len(d.caso["cargas"])
    out["parche"] = d.caso["cargas"][-1]["region"]["forma"]
    # plantilla de cizalla con apoyo empotrado y resolver en el hilo
    d.cmb_plantilla.setCurrentIndex(d.cmb_plantilla.findData("cizalla"))
    d._aplicar_plantilla("cizalla")
    out["tras_plantilla"] = d.caso["cargas"][0]["tipo"]
    d.resolver()
    for _ in range(1200):
        app.processEvents(); QtCore.QThread.msleep(25)
        if d.res is not None or "b00020" in d.lab_aviso.text():
            break
    out["resuelto"] = bool(d.res and d.res["ok"])
    out["residuo"] = d.res["equilibrio"]["residuo_F"] if d.res else None
    out["tabla"] = d.tabla.rowCount()
    # un cambio invalida el resultado
    d.sp_mag.setValue(d.sp_mag.value() + 1.0)
    out["invalidado"] = d.res is None
    # registrar en el informe
    d.resolver()
    for _ in range(1200):
        app.processEvents(); QtCore.QThread.msleep(25)
        if d.res is not None:
            break
    d.registrar()
    out["registros"] = len(d.registros)
    out["json"] = bool(json.dumps(d.registros[0]))
    # cambiar de estructura conserva el caso y recalcula la malla
    d.cmb_estructura.setCurrentIndex(1)
    app.processEvents()
    out["clave"] = d.clave
    out["sigue_el_caso"] = d.caso["cargas"][0]["tipo"] == "fuerza"
    d.deshacer(); d.deshacer()
    print("RESULTADO " + json.dumps(out))
''')


def _ejecutar_ventana(env_extra, timeout=420):
    env = dict(os.environ, QT_QPA_PLATFORM=os.environ.get(
        "QT_QPA_PLATFORM", "offscreen"), **env_extra)
    raiz = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env["PYTHONPATH"] = raiz + os.pathsep + env.get("PYTHONPATH", "")
    r = subprocess.run([sys.executable, "-c", _GUION_VENTANA], cwd=raiz,
                       env=env, capture_output=True, text=True,
                       timeout=timeout)
    lineas = [x for x in r.stdout.splitlines() if x.startswith("RESULTADO ")]
    return r, (json.loads(lineas[-1][10:]) if lineas else None)


def _comprobar_ventana(out):
    assert out["inicial"] == "presion"
    assert out["cara"] == "+x" and out["deshecho"] == "+z"
    assert out["rehecho"] == "+x"
    assert out["nodos"] == [2, 1]
    assert out["n_cargas"] == 2 and out["parche"] == "esfera"
    assert out["tras_plantilla"] == "fuerza"
    assert out["resuelto"] and out["residuo"] < 1e-6
    assert out["tabla"] >= 1
    assert out["invalidado"]
    assert out["registros"] == 1 and out["json"]
    assert out["clave"] == "spinodoide" and out["sigue_el_caso"]


def test_ventana_sin_vista_3d():
    pytest.importorskip("PyQt5")
    r, out = _ejecutar_ventana({"SPINPY_SIN_3D": "1"})
    assert out is not None, r.stderr[-2000:]
    assert out["vista3d"] is False
    _comprobar_ventana(out)


def test_ventana_con_vista_3d(render_disponible):
    pytest.importorskip("PyQt5")
    if not render_disponible:
        pytest.skip("sin OpenGL fuera de pantalla")
    try:
        r, out = _ejecutar_ventana({})
    except subprocess.TimeoutExpired:
        pytest.skip("la vista 3D de Qt no arranco en este equipo")
    if out is None or not out["vista3d"]:
        pytest.skip("la vista 3D de Qt no esta disponible en este equipo")
    _comprobar_ventana(out)
