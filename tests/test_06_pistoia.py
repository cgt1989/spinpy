"""
Bloque 6 — Criterio de fallo de Pistoia.

REFERENCIA
  Pistoia W, van Rietbergen B, Lochmuller EM, Lill CA, Eckstein F, Ruegsegger P.
  Estimation of distal radius failure load with micro-finite element analysis
  models based on three-dimensional peripheral quantitative computed tomography
  images. Bone 2002;30:842-848.
    El hueso falla cuando una fraccion `frac` (2 %) del tejido supera una
    deformacion efectiva `eps_crit` (0.7 %). Como el problema es lineal, la
    carga de fallo es la aplicada escalada por eps_crit / percentil.

LO QUE SE VERIFICA Y LO QUE NO
  Los dos parametros son CONVENCIONES calibradas en radio distal humano; no
  hay solucion cerrada que los valide. Lo que si se puede verificar de forma
  exacta es la IMPLEMENTACION:
    (1) en un cubo macizo la deformacion es uniforme, asi que el percentil es
        sigma0/E_s y la tension de fallo vale exactamente eps_crit * E_s;
    (2) la carga de fallo no depende de la carga aplicada (linealidad);
    (3) la semantica de "fraccion del tejido": con una distribucion construida
        a mano en la que el 5 % de los elementos tiene el doble de deformacion,
        frac = 2 % debe seleccionar el grupo alto y frac = 10 % el bajo.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  (1) y (2): 1e-9 relativo.   (3): 1e-12 relativo (aritmetica pura).
"""
import numpy as np

from spinpy.resistencia import criterio_pistoia, ensayo_compresion

BLOQUE = "06 Criterio de Pistoia"
REF = "Pistoia et al. 2002, Bone 30:842"
E_S, NU = 1.0, 0.30


def _macizo(s0):
    r = ensayo_compresion(np.ones((10, 10, 10), bool), [1.0] * 3, E_s=E_S,
                          nu_s=NU, sigma0=s0, apoyo="deslizante")
    assert r["ok"], r["msg"]
    return r


def test_cubo_macizo_tension_de_fallo_exacta(registro):
    s0 = 1e-3
    r = _macizo(s0)
    p = criterio_pistoia(r, frac=0.02, eps_crit=0.007)
    assert p["ok"], p
    esp_sigma = 0.007 * E_S
    esp_F = esp_sigma * r["A_bruta"]
    esp_k = 0.007 * E_S / s0
    for nombre, esp, obt in (("sigma_fallo = eps_crit * E_s", esp_sigma,
                              p["sigma_fallo"]),
                             ("F_fallo = sigma_fallo * A_bruta", esp_F,
                              p["F_fallo"]),
                             ("factor = eps_crit * E_s / sigma0", esp_k,
                              p["factor"])):
        e = abs(obt - esp) / abs(esp)
        registro.anotar(BLOQUE, f"cubo macizo: {nombre}", REF, esp, obt, 1e-9,
                        "relativo; deformacion uniforme", e <= 1e-9)
        assert e <= 1e-9, (nombre, obt, esp)


def test_linealidad_carga_de_fallo(registro):
    pa = criterio_pistoia(_macizo(1e-3))
    pb = criterio_pistoia(_macizo(2e-3))
    e = abs(pa["sigma_fallo"] - pb["sigma_fallo"]) / pa["sigma_fallo"]
    registro.anotar(BLOQUE, "sigma_fallo independiente de sigma0 (linealidad)",
                    REF, pa["sigma_fallo"], pb["sigma_fallo"], 1e-9,
                    "relativo", e <= 1e-9)
    assert e <= 1e-9
    ok = abs(pa["factor"] / pb["factor"] - 2.0) < 1e-9
    registro.anotar(BLOQUE, "factor se divide por 2 al duplicar sigma0", REF,
                    2.0, pa["factor"] / pb["factor"], 1e-9, "relativo", ok)
    assert ok


def test_deformacion_total_cubo_macizo(registro):
    """|u| en el cubo macizo, con la expansion de Poisson incluida.

    La "Total Deformation" de ANSYS es el MODULO del desplazamiento, no su
    componente axial, asi que en un ensayo deslizante incluye la dilatacion
    lateral. Con uz = 0 en toda la base y los movimientos de solido rigido
    eliminados en DOS ESQUINAS, esa dilatacion se mide desde una esquina y
    vale nu * eps_z * L en cada direccion transversal, de modo que

        |u|_max = sqrt( uz^2 + 2 (nu eps_z L)^2 )     con uz = eps_z L

    Es una solucion cerrada, y verifica de paso que `carga_N` reparte una
    fuerza total correcta: sigma_app tiene que salir carga_N / A_bruta.

    TOLERANCIA DECLARADA ANTES DE MEDIR: 1e-6 relativo (solver iterativo con
    tol 1e-8; aqui no se compara elemento a elemento sino un maximo).

    CONSECUENCIA QUE CONVIENE TENER PRESENTE: bajo apoyo DESLIZANTE el maximo
    de |u| depende de donde se anclen los modos de solido rigido, porque la
    dilatacion lateral se mide desde ese ancla. Bajo apoyo EMPOTRADO -el del
    articulo- la base esta toda fija y esa ambiguedad no existe.
    """
    E_s, nu, n = 20e9, 0.30, 10
    sp = 1e-4                      # 0.1 mm -> cubo de 1 mm
    L = n * sp
    F = 100.0
    BW = np.ones((n, n, n), bool)
    # unidad='m': aqui el sistema es SI puro (spacing en metros, E en Pa), a
    # diferencia de la aplicacion, que trabaja en mm con E en Pa.
    r = ensayo_compresion(BW, [sp] * 3, E_s=E_s, nu_s=nu,
                          apoyo="deslizante", carga_N=F, unidad="m")
    assert r["ok"], r["msg"]

    sigma = F / (L * L)
    e_sig = abs(r["sigma_app"] - sigma) / sigma
    registro.anotar(BLOQUE, "carga_N: sigma_app = carga_N / A_bruta",
                    "convencion de seccion bruta", sigma, r["sigma_app"],
                    1e-12, "relativo", e_sig <= 1e-12)
    assert e_sig <= 1e-12

    eps_z = sigma / E_s
    esperado = np.sqrt((eps_z * L) ** 2 + 2 * (nu * eps_z * L) ** 2)
    e = abs(r["desp_max"] - esperado) / esperado
    registro.anotar(BLOQUE,
                    "deformacion total max = sqrt(uz^2 + 2(nu eps L)^2)",
                    "modulo del desplazamiento con Poisson", esperado,
                    r["desp_max"], 1e-6, "relativo; solver tol 1e-8",
                    e <= 1e-6,
                    nota="la componente axial sola valdria %.4e" % (eps_z * L))
    assert e <= 1e-6, (r["desp_max"], esperado)

    # El campo por elemento existe, es finito donde hay hueso y no supera el
    # maximo nodal (es un promedio de ocho nodos).
    c = r["campo_desp"]
    assert c.shape == BW.shape
    assert np.isfinite(c).sum() == BW.size
    registro.anotar(BLOQUE, "campo_desp <= maximo nodal (es un promedio)",
                    "coherencia del promediado", None,
                    float(np.nanmax(c)) / r["desp_max"], "<= 1",
                    "cociente", np.nanmax(c) <= r["desp_max"] * (1 + 1e-12))
    assert np.nanmax(c) <= r["desp_max"] * (1 + 1e-12)


def test_unidad_de_la_carga(registro):
    """En mm con E en Pa, 100 N exigen el factor 1e6. Sin el, todo es 1e-6.

    Es el error que tuvo la primera version del analisis comparado: pedia
    100 N con el spacing en milimetros y aplicaba de hecho 1e-4 N, con lo que
    la von Mises salia 0.000 MPa y la deformacion total 0.00000 mm. La prueba
    fija la equivalencia entre los dos sistemas.
    """
    E_s, nu, n = 20e9, 0.30, 8
    BW = np.ones((n, n, n), bool)
    F = 100.0
    r_mm = ensayo_compresion(BW, [0.1] * 3, E_s=E_s, nu_s=nu,          # mm
                             apoyo="deslizante", carga_N=F, unidad="mm")
    r_m = ensayo_compresion(BW, [1e-4] * 3, E_s=E_s, nu_s=nu,          # m
                            apoyo="deslizante", carga_N=F, unidad="m")
    assert r_mm["ok"] and r_m["ok"]
    # La misma probeta fisica: las tensiones deben coincidir y los
    # desplazamientos diferir exactamente en el factor 1000 (mm frente a m).
    e_sig = abs(r_mm["sigma_app"] - r_m["sigma_app"]) / r_m["sigma_app"]
    e_u = abs(r_mm["desp_max"] / r_m["desp_max"] - 1000.0) / 1000.0
    registro.anotar(BLOQUE, "carga en mm-Pa = carga en m-Pa (tension)",
                    "coherencia de unidades", r_m["sigma_app"],
                    r_mm["sigma_app"], 1e-9, "relativo", e_sig <= 1e-9)
    registro.anotar(BLOQUE, "desplazamiento en mm = 1000 x el de metros",
                    "coherencia de unidades", 1000.0,
                    r_mm["desp_max"] / r_m["desp_max"], 1e-6, "relativo",
                    e_u <= 1e-6)
    assert e_sig <= 1e-9, (r_mm["sigma_app"], r_m["sigma_app"])
    assert e_u <= 1e-6, (r_mm["desp_max"], r_m["desp_max"])


def test_semantica_de_la_fraccion(registro):
    """5 % de los elementos a 2e, el resto a e. frac=2 % -> 2e; frac=10 % -> e."""
    e0 = 1e-3
    eps = np.full(1000, e0)
    eps[:50] = 2 * e0                     # el 5 % mas deformado
    res = {"ok": True, "eps_eff_solido": eps, "sigma_app": 1.0,
           "F_total": 1.0, "vm_solido": np.array([])}
    p2 = criterio_pistoia(res, frac=0.02, eps_crit=0.007)
    p10 = criterio_pistoia(res, frac=0.10, eps_crit=0.007)
    ok2 = abs(p2["eps_eff_p"] - 2 * e0) / (2 * e0) <= 1e-12
    ok10 = abs(p10["eps_eff_p"] - e0) / e0 <= 1e-12
    registro.anotar(BLOQUE, "frac = 2 % selecciona el grupo mas deformado", REF,
                    2 * e0, p2["eps_eff_p"], 1e-12, "percentil por conteo", ok2)
    registro.anotar(BLOQUE, "frac = 10 % cae en el grupo base", REF, e0,
                    p10["eps_eff_p"], 1e-12, "percentil por conteo", ok10)
    assert ok2 and ok10, (p2["eps_eff_p"], p10["eps_eff_p"])
    ok = abs(p2["factor"] * 2 - p10["factor"]) / p10["factor"] <= 1e-12
    registro.anotar(BLOQUE, "factor(2 %) = factor(10 %) / 2", REF,
                    p10["factor"] / 2, p2["factor"], 1e-12, "relativo", ok)
    assert ok
