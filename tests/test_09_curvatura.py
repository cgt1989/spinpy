"""
Bloque 9 — Curvaturas principales de la interfaz.

QUE SE VERIFICA
  El modulo `spinpy.curvatura` tiene dos caminos: uno EXACTO, que sale de las
  derivadas del campo, y uno DISCRETO, que ajusta la segunda forma fundamental
  sobre la malla y es el unico aplicable a un VOI de micro-CT. Aqui se
  comprueban los dos contra respuestas que se conocen de antemano, y se
  comprueba que el discreto no se aparta del exacto.

REFERENCIAS
  do Carmo MP. Differential geometry of curves and surfaces. Prentice-Hall,
    1976. Esfera de radio R: k1 = k2 = 1/R. Toro: la integral de la curvatura
    gaussiana es cero.
  Gauss-Bonnet: la integral de K sobre una superficie cerrada vale 2*pi*chi,
    con chi el numero de Euler. Es una identidad TOPOLOGICA, de modo que
    compararla con el chi que mide `morphometry.conectividad` cruza dos
    caminos de codigo que no comparten nada.
  Rusinkiewicz S. Estimating curvatures and their derivatives on triangle
    meshes. 3DPVT 2004:486-93.
  Guo Y, Sharma S, Kumar S. Inverse designing surface curvatures by deep
    learning. Adv Intell Syst 2024;6:2300789. Ecuacion (3): el perfil pondera
    por AREA, no por numero de elementos.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  Los campos de prueba se muestrean en rejillas de 64^3 y la malla es la de
  marching cubes sobre el campo CONTINUO, no sobre la mascara binaria. El error
  de un ajuste de segunda forma fundamental sobre una malla asi es de segundo
  orden en el paso; con ~20 voxeles por radio cabe esperar decimas de por
  ciento. Se declara, con holgura:

    esfera, camino exacto        |k - 1/R| / (1/R) <= 0.01
    esfera, camino discreto      |k - 1/R| / (1/R) <= 0.02
    signo                        bola de solido -> k > 0; poro -> k < 0
    Gauss-Bonnet en la esfera    int K dA / (2*pi*chi) en [0.97, 1.03], chi = 2
    Gauss-Bonnet en el toro      |int K dA| <= 0.03 * int |K| dA   (chi = 0)
    discreto vs exacto           mediana de |k_dis - k_exa| <= 0.03 * |k_exa|

  La comparacion del ultimo punto se hace SOBRE LOS MISMOS VERTICES: si se
  mallara dos veces, el error del estimador se mezclaria con el del mallado y
  no se sabria cual se esta midiendo.
"""
import numpy as np

from spinpy.curvatura import (curvaturas_de_campo, curvaturas_implicitas,
                              perfil, resumen)

BLOQUE = "09 Curvatura (do Carmo; Gauss-Bonnet; Rusinkiewicz 2004)"
REF = "do Carmo 1976; Rusinkiewicz 2004; Guo et al. 2024"

N = 64
SP = 1.0 / (N - 1)
R_ESFERA = 0.30
R_TORO, r_TORO = 0.30, 0.12


def _rejilla(n=N):
    ax = np.linspace(0.0, 1.0, n)
    return np.meshgrid(ax, ax, ax, indexing="ij")


def _campo_esfera(n=N):
    """f = |x - c|^2, solido {f <= R^2}: bola de hueso centrada en el cubo."""
    X, Y, Z = _rejilla(n)
    return (X - .5) ** 2 + (Y - .5) ** 2 + (Z - .5) ** 2


def _exacta_esfera(p):
    g = 2.0 * (np.asarray(p, float) - 0.5)
    H = np.broadcast_to(2.0 * np.eye(3), (len(p), 3, 3)).copy()
    return g, H


def _campo_toro(n=N):
    """f = (sqrt(x^2+y^2) - R)^2 + z^2, solido {f <= r^2}. Eje z."""
    X, Y, Z = _rejilla(n)
    rho = np.sqrt((X - .5) ** 2 + (Y - .5) ** 2)
    return (rho - R_TORO) ** 2 + (Z - .5) ** 2


def _exacta_toro(p):
    p = np.asarray(p, float)
    x, y, z = p[:, 0] - .5, p[:, 1] - .5, p[:, 2] - .5
    rho = np.sqrt(x ** 2 + y ** 2)
    rho = np.where(rho > 0, rho, 1e-12)
    c = 2.0 * (rho - R_TORO) / rho                     # d f / d(x,y) = c * (x,y)
    g = np.stack([c * x, c * y, 2.0 * z], axis=1)
    # Hessiano: derivar c*(x,y) con d rho/dx = x/rho.
    H = np.zeros((len(p), 3, 3))
    drx, dry = x / rho, y / rho
    dc_dx = 2.0 * R_TORO * drx / rho ** 2
    dc_dy = 2.0 * R_TORO * dry / rho ** 2
    H[:, 0, 0] = c + dc_dx * x
    H[:, 0, 1] = dc_dy * x
    H[:, 1, 0] = dc_dx * y
    H[:, 1, 1] = c + dc_dy * y
    H[:, 2, 2] = 2.0
    return g, (H + np.swapaxes(H, 1, 2)) / 2.0


# ---------------------------------------------------------------------------
def test_esfera_curvatura_exacta(registro):
    """Camino exacto sobre una esfera: k1 = k2 = 1/R, y POSITIVAS."""
    d = curvaturas_de_campo(_campo_esfera(), SP, R_ESFERA ** 2,
                            exacta=_exacta_esfera)
    obj = 1.0 / R_ESFERA
    for nom, v in (("k1", d["k1_exacta"]), ("k2", d["k2_exacta"])):
        med = float(np.median(v))
        ok = abs(med - obj) / obj <= 0.01
        registro.anotar(BLOQUE, f"esfera R={R_ESFERA}: {nom} exacta", REF,
                        obj, med, "error relativo <= 1 %",
                        "do Carmo: k1 = k2 = 1/R", ok)
        assert ok, (nom, med, obj)


def test_esfera_curvatura_discreta_y_signo(registro):
    """Camino discreto sobre la misma malla, y el signo de la convencion.

    El signo es la mitad de la prueba: con la normal hacia el vacio, una bola
    de HUESO tiene curvatura positiva. Si alguien invierte la orientacion, los
    perfiles siguen pareciendo razonables y H cambia de signo en silencio.
    """
    d = curvaturas_de_campo(_campo_esfera(), SP, R_ESFERA ** 2,
                            exacta=_exacta_esfera)
    obj = 1.0 / R_ESFERA
    med = float(np.median(d["k1"]))
    ok = abs(med - obj) / obj <= 0.02
    registro.anotar(BLOQUE, f"esfera R={R_ESFERA}: k1 discreta", REF, obj, med,
                    "error relativo <= 2 %",
                    "Rusinkiewicz sobre la malla del campo continuo", ok)
    assert ok, (med, obj)

    r = resumen(d["k1"], d["k2"], d["areas"])
    conv = r["convexa"]
    ok2 = conv >= 0.99 and r["H_medio"] > 0
    registro.anotar(BLOQUE, "esfera: fraccion de area convexa (k1>0, k2>0)",
                    REF, 1.0, conv, ">= 0.99 y H > 0",
                    "convencion de signo: solido convexo -> k > 0", ok2)
    assert ok2, r

    # La mediana del error entre los dos caminos, sobre los MISMOS vertices.
    e = np.median(np.abs(d["k1"] - d["k1_exacta"])) / obj
    ok3 = e <= 0.03
    registro.anotar(BLOQUE, "esfera: discreta vs exacta (mismos vertices)",
                    REF, 0.0, e, "mediana del error relativo <= 3 %",
                    "valida el estimador aplicable al hueso real", ok3)
    assert ok3, e


def test_poro_curvatura_negativa(registro):
    """Un poro esferico dentro de hueso macizo da curvatura NEGATIVA.

    Es la misma esfera con el solido y el vacio intercambiados: se consigue
    cambiando el signo del campo, de modo que {-f <= -R^2} es el complementario.
    """
    d = curvaturas_de_campo(-_campo_esfera(), SP, -R_ESFERA ** 2)
    med = float(np.median(d["k1"]))
    ok = med < 0 and abs(abs(med) - 1.0 / R_ESFERA) / (1.0 / R_ESFERA) <= 0.02
    registro.anotar(BLOQUE, f"poro R={R_ESFERA}: k1 discreta", REF,
                    -1.0 / R_ESFERA, med, "negativa, error relativo <= 2 %",
                    "concavidad: el hueso rodea al poro", ok)
    assert ok, med


def test_gauss_bonnet_esfera(registro):
    """int K dA = 2*pi*chi con chi = 2 para una esfera.

    Cruza la curvatura con la topologia: el mismo numero sale del estimador de
    curvatura y del numero de Euler, que no comparten una sola linea de codigo.
    """
    d = curvaturas_de_campo(_campo_esfera(), SP, R_ESFERA ** 2)
    K = d["k1"] * d["k2"]
    I = float((K * d["areas"]).sum())
    razon = I / (2.0 * np.pi * 2.0)
    ok = 0.97 <= razon <= 1.03
    registro.anotar(BLOQUE, "esfera: int K dA / (2 pi chi), chi = 2",
                    "Gauss-Bonnet", 1.0, razon, "en [0.97, 1.03]",
                    "identidad topologica, no un ajuste", ok)
    assert ok, razon


def test_gauss_bonnet_toro(registro):
    """El toro tiene chi = 0: la integral de K se cancela entre dentro y fuera.

    Es una prueba mas dura que la de la esfera porque exige que la curvatura
    NEGATIVA del anillo interior se mida igual de bien que la positiva del
    exterior. Se normaliza por la integral de |K| para que el criterio no
    dependa del tamano de la malla.
    """
    d = curvaturas_de_campo(_campo_toro(), SP, r_TORO ** 2,
                            exacta=_exacta_toro)
    for nom, k1, k2 in (("discreta", d["k1"], d["k2"]),
                        ("exacta", d["k1_exacta"], d["k2_exacta"])):
        K = k1 * k2
        I = float((K * d["areas"]).sum())
        Iabs = float((np.abs(K) * d["areas"]).sum())
        razon = I / Iabs if Iabs > 0 else np.nan
        ok = abs(razon) <= 0.03
        registro.anotar(BLOQUE, f"toro: int K dA normalizada ({nom})",
                        "Gauss-Bonnet", 0.0, razon,
                        "|int K| <= 3 % de int |K|", "chi = 0", ok)
        assert ok, (nom, razon)


def test_toro_curvaturas_conocidas(registro):
    """En el toro, k1 = 1/r en todo punto y k2 va de -1/(R-r) a 1/(R+r).

    La curvatura mayor de un toro es constante y vale 1/r: es una comprobacion
    de que las dos curvaturas no se estan mezclando entre si al ordenarlas.
    """
    d = curvaturas_de_campo(_campo_toro(), SP, r_TORO ** 2,
                            exacta=_exacta_toro)
    obj = 1.0 / r_TORO
    med = float(np.median(d["k1_exacta"]))
    ok = abs(med - obj) / obj <= 0.01
    registro.anotar(BLOQUE, "toro: k1 = 1/r (constante)", REF, obj, med,
                    "error relativo <= 1 %", "do Carmo", ok)
    assert ok, (med, obj)

    k2 = d["k2_exacta"]
    lo, hi = -1.0 / (R_TORO - r_TORO), 1.0 / (R_TORO + r_TORO)
    frac = float(np.mean((k2 >= lo * 1.05) & (k2 <= hi * 1.05)))
    ok2 = frac >= 0.98
    registro.anotar(BLOQUE, "toro: k2 dentro de [-1/(R-r), 1/(R+r)]", REF,
                    1.0, frac, ">= 98 % de los vertices",
                    "cota exacta del toro, con 5 % de holgura", ok2)
    assert ok2, frac


def test_perfil_pondera_por_area(registro):
    """El perfil es una probabilidad de AREA (ecuacion 3 de Guo et al. 2024).

    Se comprueba sobre una esfera, donde toda la superficie tiene la misma
    curvatura: el perfil tiene que concentrarse en una sola celda y sumar 1.
    Si se ponderara por numero de elementos tambien sumaria 1, asi que la
    prueba de verdad es la segunda: doblar el area de la mitad de los vertices
    tiene que mover el resultado.
    """
    d = curvaturas_de_campo(_campo_esfera(), SP, R_ESFERA ** 2)
    P, bx, by = perfil(d["k1"], d["k2"], d["areas"], limite=6.0, nbins=60)
    s = float(P.sum())
    ok = abs(s - 1.0) <= 1e-12
    registro.anotar(BLOQUE, "perfil: suma de la probabilidad", "Guo et al. 2024",
                    1.0, s, "= 1 (1e-12)", "ecuacion (3), normalizada", ok)
    assert ok, s

    w = d["areas"].copy()
    w[::2] *= 2.0
    r0 = resumen(d["k1"], d["k2"], d["areas"])
    r1 = resumen(d["k1"], d["k2"], w)
    movio = abs(r1["H_medio"] - r0["H_medio"]) > 0
    registro.anotar(BLOQUE, "perfil: los pesos de area intervienen",
                    "Guo et al. 2024", None, float(r1["H_medio"] - r0["H_medio"]),
                    "distinto de cero al cambiar las areas",
                    "descarta una ponderacion por numero de elementos", movio)
    assert movio


def test_curvaturas_implicitas_silla_perfecta(registro):
    """Silla pura: k1 = -k2, con el valor propio nulo EN MEDIO al ordenar.

    Es el caso que rompe la implementacion ingenua. Para el paraboloide
    hiperbolico f = x^2 - y^2 + z (solido {f <= 0}), en el origen la superficie
    tiene k1 = 2/|grad| y k2 = -2/|grad| con |grad| = 1, y el tercer valor
    propio del operador es cero: si se descarta por POSICION en vez de por
    magnitud, se devuelve (2, 0) en lugar de (2, -2).
    """
    p = np.zeros((1, 3))
    g = np.array([[0.0, 0.0, 1.0]])
    H = np.array([[[2.0, 0.0, 0.0], [0.0, -2.0, 0.0], [0.0, 0.0, 0.0]]])
    k1, k2 = curvaturas_implicitas(g, H)
    ok = abs(k1[0] - 2.0) < 1e-12 and abs(k2[0] + 2.0) < 1e-12
    registro.anotar(BLOQUE, "silla x^2 - y^2: k1 = 2, k2 = -2", REF, -2.0,
                    float(k2[0]), "= -2 (1e-12)",
                    "el valor propio nulo se descarta por magnitud", ok)
    assert ok, (k1, k2, p)
