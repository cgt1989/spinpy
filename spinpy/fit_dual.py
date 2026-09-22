"""
fit_dual.py — Ajuste de un dual-lattice a un VOI por busqueda escalonada.

Es la misma busqueda que `fit.ajustar_spinodoide` —etapas A, B, C y D, el
desempate K2 por orientacion, la alineacion impuesta L2, el objetivo O1 y la
incertidumbre U1— aplicada a la otra familia. Se escribe aparte: el ajuste del
spinodoide esta contrastado contra `AppFinal_V2.m` (`validar_ajuste.py`), y lo
que aqui se reutiliza de `fit.py` se importa tal cual.

LAS VARIABLES DE BUSQUEDA, UNA A UNA FRENTE AL SPINODOIDE
--------------------------------------------------------
    etapa   spinodoide                     dual-lattice
    A       densidad x numero de onda      densidad x celdas
    B       presets de thetas (conos)      presets de estiramiento
    C       densidad x numero de ondas     densidad x celdas (refinado)
    D       desempate entre finalistas     desempate entre finalistas

La funcion de error es la MISMA (`error.error_morfometrico`) y el candidato se
mide por el MISMO camino que el VOI (`morphometry.morfometria`) al MISMO tamano
fisico (correccion F1). Sin eso, comparar el error de las dos familias no
significaria nada.

LA REJILLA DE CELDAS SALE DEL VOI, NO DE UN RANGO FIJO
------------------------------------------------------
El spinodoide recorre un rango fijo de numero de onda, [8, 25] pi, porque la
correccion G2 mostro que una ventana estrecha alrededor del deslizador no
alcanzaba las trabeculas finas. Aqui no hay deslizador del que partir, y un
rango fijo tendria el mismo problema al reves: el numero util de celdas depende
de cuantas trabeculas caben en el VOI.

Se estima a partir de lo que el VOI ya dice: Tb.N (trabeculas por mm) por el
lado fisico da trabeculas por lado del cubo, y la relacion entre celdas y Tb.N
del dual-lattice se MIDIO (`CELDAS_POR_TBN`, ver su comentario). La rejilla se
abre alrededor de esa estimacion lo bastante —un factor 2 a cada lado en el
modo completo— para que un error de la calibracion no deje la respuesta fuera.
Y se acota por arriba para que cada celda tenga al menos `VOX_POR_CELDA_MIN`
voxeles: por debajo, los puntales tienen uno o dos voxeles de grosor y la
morfometria mide la discretizacion, no la estructura.
"""

from __future__ import annotations

import time

import numpy as np

from . import incertidumbre, procedencia
from .dual_lattice import generar_dual_lattice
from .error import error_morfometrico, medidas_necesarias
from .fit import (CLAVES_MECANICAS, E_S_DEF, NU_S_DEF, REPLICAS_DEF, RHO_MAX,
                  RHO_MIN, SEMILLA_POR_DEFECTO, _angulo_ejes,
                  _desempatar_theta, _rot_entre_ejes, alinear_marco_fabrica,
                  claves_incertidumbre, dir_a_euler, fraccion_poro_conexa,
                  longitud_caracteristica, metricas_mecanicas,
                  resolucion_comparacion, validar_opciones)
from .grf import euler_R
from .morphometry import metricas_forma, morfometria

# Celdas por cada trabecula por lado del cubo, MEDIDO sobre el generador
# (cubo unidad, 64^3, semilla del proyecto, Tb.N = BV/TV / Tb.Th de Parfitt):
#
#     rho          0.20    0.30    0.40
#     celdas/Tb.N  0.61    0.54    0.52      (celdas 4, 6 y 8; estiramiento
#                                              1 y 1.5 dan lo mismo a +-0.02)
#
# No es constante: a mas densidad los puntales engordan y Parfitt cuenta menos
# trabeculas por el mismo numero de celdas. Se toma el valor de rho 0.30, que es
# el de los VOIs proximales, y la rejilla se abre un factor 2 a cada lado (1.4
# en el modo rapido), mucho mas que el +-12 % que varia la relacion.
CELDAS_POR_TBN = 0.55

VOX_POR_CELDA_MIN = 5.0
CELDAS_MIN = 2.0

# Los presets de estiramiento salen de lo que el estiramiento HACE al DA, que
# se midio (VOI proximal de H4: 4.55 celdas, rho 0.28, 96^3, dos semillas):
#
#     estiramiento z   1.0    1.25   1.5    2.0    2.5    3.0    4.0
#     DA               1.04   1.09   1.18   1.38   1.56   1.69   1.85
#     DA2              1.02   1.07   1.18   1.37   1.53   1.67   1.82
#
# El VOI tiene DA 1.51: cae entre 2.0 y 2.5. Una primera lista con tope en 2.0
# —puesta a ojo, antes de medir— dejaba la anisotropia del hueso fuera del
# espacio de busqueda, el mismo defecto que la correccion G1 arreglo en el
# spinodoide. Por eso las intensidades llegan a 3.0 (DA 1.69): hueso algo mas
# anisotropo que este sigue dentro.
#
# Alargar un solo eje da DA2 ~ DA: estructura AXISIMETRICA. El VOI no lo es
# (DA2 1.18 frente a DA 1.51), asi que se anaden presets TRIAXIALES, con el eje
# largo en z porque la rotacion de partida lleva z sobre la direccion principal
# del VOI.
#
# Las permutaciones de eje (K2) se mantienen: el eje alargado es el eje RIGIDO,
# y sin permutarlo el desempate por orientacion no tendria entre que elegir.
#
# LOS TRIAXIALES, MEDIDOS (Estudio_Familias/code/presets_triaxiales.py, mismas
# condiciones que la tabla de arriba, dos semillas; DA2 = autovalor mayor
# entre el intermedio):
#
#     (1, ey, ez)      DA     DA2        (1, ey, ez)      DA     DA2
#     (1, 1.0, 2.5)   1.56   1.53        (1, 2.0, 2.5)   1.30   1.14
#     (1, 1.5, 2.5)   1.38   1.29        (1, 2.5, 3.25)  1.37   1.18
#     (1, 1.5, 3.0)   1.54   1.46        (1, 3.0, 3.9)   1.40   1.19
#     (1, 2.0, 3.0)   1.40   1.24        (1, 3.0, 4.35)  1.47   1.27
#
# DA2 lo fija sobre todo el cociente ez/ey, y DA el alargamiento de esos dos
# ejes frente a x. Pero alargar DOS ejes satura el DA: con ez/ey 1.3 no pasa
# de ~1.40 aunque ey llegue a 4 (DA 1.415, DA2 1.19). El hueso equino (DA
# 1.45-1.51, DA2 1.12-1.18) queda justo fuera de lo que la familia alcanza;
# lo mas cercano es DA ~1.40 con DA2 ~1.18. Los presets anteriores solo
# ofrecian DA2 >= 1.29 a ese DA, que es el sesgo que midio el informe.
#
# Se anaden en las dos permutaciones del eje INTERMEDIO (y o x), con el largo
# siempre en z: la rotacion de partida lleva z sobre la direccion principal del
# VOI, pero el eje secundario del hueso puede caer en cualquiera de los dos.
PRESETS_RAPIDO = [(1.0, 1.0, 1.0), (1.0, 1.0, 1.5), (1.0, 1.0, 2.0),
                  (1.0, 1.0, 2.5), (1.0, 2.0, 1.0), (2.0, 1.0, 1.0),
                  (1.0, 2.5, 3.25)]

PRESETS_COMPLETO = [(1.0, 1.0, 1.0),
                    (1.0, 1.0, 1.5), (1.0, 1.0, 2.0), (1.0, 1.0, 2.5),
                    (1.0, 1.0, 3.0),
                    (1.0, 1.5, 1.0), (1.0, 2.0, 1.0), (1.0, 2.5, 1.0),
                    (1.0, 3.0, 1.0),
                    (1.5, 1.0, 1.0), (2.0, 1.0, 1.0), (2.5, 1.0, 1.0),
                    (3.0, 1.0, 1.0),
                    (1.0, 1.5, 2.5), (1.5, 1.0, 2.5),
                    (1.0, 2.0, 2.5), (2.0, 1.0, 2.5),
                    (1.0, 2.5, 3.25), (2.5, 1.0, 3.25),
                    (1.0, 3.0, 3.9), (3.0, 1.0, 3.9),
                    (1.0, 3.0, 4.35), (3.0, 1.0, 4.35)]


# Peso del termino DA2 en el error de ESTA familia (en el spinodoide es 0, ver
# `error.error_morfometrico`). Sin el, un preset triaxial que iguala el DA solo
# empata con el axial del mismo DA, y el desempate K2 decide por orientacion,
# no por la forma de la fabrica. 1.0 es el peso de BS/BV, Tb.Th o Tb.N: un
# termino mas, no uno que mande.
PESO_DA2 = 1.0


def celdas_estimadas(m_voi, lc):
    """Numero de celdas que corresponde al Tb.N del VOI en un cubo de lado lc."""
    tbn = float(m_voi.get("TbN", np.nan))
    if not np.isfinite(tbn) or tbn <= 0:
        return 6.0
    return CELDAS_POR_TBN * tbn * float(lc)


def alinear_dual_por_fabrica(m_voi, generar, sp_cand, R_clasica, da_min=1.10,
                             da2_max=1.06):
    """Correccion L1 para el dual-lattice: misma logica que
    `fit.alinear_por_fabrica`, con el generador pasado como funcion R -> BW.

    Se conserva como `alineacion="eje"`; por omision se usa L2
    (`fit.alinear_marco_fabrica`).

    En el dual-lattice la rotacion es EXACTA: la red se genera en su marco con
    la misma semilla y despues se gira, asi que girar no remuestrea nada. Aun
    asi se vuelve a medir tras rotar, porque la voxelizacion de la red girada
    no es la voxelizacion girada.
    """
    u_voi = np.asarray(m_voi.get("dir_principal", [0.0, 0.0, 1.0]), float)
    info = {"metodo": "clasica (eje z del candidato sobre la del VOI)",
            "aplicada": False, "motivo": "", "angulo_antes_deg": float("nan"),
            "angulo_despues_deg": float("nan"),
            "DA_candidato": float("nan"), "DA2_candidato": float("nan")}

    da_voi = float(m_voi.get("DA", np.nan))
    if not np.isfinite(da_voi) or da_voi < da_min:
        info["motivo"] = (f"el VOI es casi isotropo (DA {da_voi:.3f} < "
                          f"{da_min}): su direccion principal es ruido")
        return R_clasica, info

    m0 = morfometria(generar(np.eye(3)), sp_cand)
    u0 = np.asarray(m0.get("dir_principal", [0.0, 0.0, 1.0]), float)
    da0 = float(m0.get("DA", np.nan))
    da2_0 = float(m0.get("DA2", np.nan))
    info["DA_candidato"] = da0
    info["DA2_candidato"] = da2_0
    if not np.isfinite(da0) or da0 < da_min:
        info["motivo"] = (f"el candidato es casi isotropo (DA {da0:.3f}): no "
                          f"hay eje que alinear")
        return R_clasica, info
    if np.isfinite(da2_0) and da2_0 < da2_max:
        info["motivo"] = (
            f"la fabrica del candidato es degenerada (DA2 {da2_0:.3f} < "
            f"{da2_max}): es rigido en un PLANO, no en un eje")
        return R_clasica, info

    R = _rot_entre_ejes(u0, u_voi)
    m1 = morfometria(generar(R), sp_cand)
    u1 = np.asarray(m1.get("dir_principal", [0.0, 0.0, 1.0]), float)
    info.update({
        "metodo": "fabrica medida (correccion L1)", "aplicada": True,
        "angulo_antes_deg": _angulo_ejes(u0, u_voi),
        "angulo_despues_deg": _angulo_ejes(u1, u_voi),
        "dir_candidato_sin_rotar": [float(x) for x in u0],
        "dir_candidato_rotado": [float(x) for x in u1],
        "dir_voi": [float(x) for x in u_voi],
    })
    if info["angulo_despues_deg"] > info["angulo_antes_deg"]:
        info["aplicada"] = False
        info["motivo"] = ("tras rotar, el eje quedo mas desalineado que antes; "
                          "se conserva la rotacion clasica")
        return R_clasica, info
    return R, info


def ajustar_dual_lattice(VOI, spacing, modo="completo", m_voi=None,
                         seed=SEMILLA_POR_DEFECTO, irregularidad=0.5,
                         resolucion=None, progreso=None, rho_max=RHO_MAX,
                         peso_mecanico=0.0, res_mec=16, n_finalistas=5,
                         E_s=E_S_DEF, nu_s=NU_S_DEF, alinear_fabrica=True,
                         peso_da2=PESO_DA2, pesos=None,
                         distancia="wasserstein", replicas=REPLICAS_DEF,
                         seleccion="minimo", k_robusto=1.0,
                         alineacion="impuesta"):
    """Busca los parametros de dual-lattice que mejor reproducen el VOI.

    modo : 'completo' (A+B+C, 23 presets) | 'rapido' (A+B, 7 presets)

    El resto de argumentos significan lo mismo que en `ajustar_spinodoide`.
    Devuelve un dict con la MISMA estructura que aquel, para que la interfaz
    y `metodos.comparar_resultados` lo lean igual; los parametros propios de
    la familia van en `parametros` y `familia` vale "dual-lattice".
    """
    pesos, K = validar_opciones(pesos, distancia, replicas, seleccion,
                                alineacion)
    nec = medidas_necesarias(pesos)
    forma = {k: True for k in nec["forma"]}
    caro = {k: True for k in nec["caro"]}
    usa_mec = bool(peso_mecanico > 0 or nec["mecanico"])
    etapa_d = bool(usa_mec or caro)

    VOI = np.asarray(VOI, dtype=bool)
    spacing = np.asarray(spacing, dtype=float).ravel()
    if m_voi is None:
        m_voi = morfometria(VOI, spacing)
    if forma:
        m_voi = dict(m_voi)
        m_voi.update(metricas_forma(VOI, spacing, **forma))

    lc = longitud_caracteristica(VOI.shape, spacing)
    res = int(resolucion) if resolucion else resolucion_comparacion(VOI.shape)
    sp_cand = np.full(3, lc / res)

    bvtv = float(m_voi["BVTV"])
    rx, ry, rz = dir_a_euler(m_voi["dir_principal"])
    R_fit = euler_R(rx, ry, rz)

    c_est = celdas_estimadas(m_voi, lc)
    c_max = max(CELDAS_MIN, res / VOX_POR_CELDA_MIN)
    e0 = PRESETS_RAPIDO[0]

    if modo == "rapido":
        dens = np.linspace(max(RHO_MIN, bvtv * 0.85), min(rho_max, bvtv * 1.15), 3)
        celdas = c_est * np.array([0.7, 1.0, 1.4])
        presets = PRESETS_RAPIDO
        refinar = False
    else:
        dens = np.linspace(max(RHO_MIN, bvtv * 0.7), min(rho_max, bvtv * 1.3), 5)
        celdas = c_est * np.geomspace(0.5, 2.0, 6)
        presets = PRESETS_COMPLETO
        refinar = True
    celdas = np.unique(np.round(np.clip(celdas, CELDAS_MIN, c_max), 2))

    n_total = len(dens) * len(celdas) + len(presets) + (9 if refinar else 0)
    traza = []
    t0 = time.time()
    contador = {"i": 0}

    def generar(d, c, e, R, s=None):
        BW, _, _ = generar_dual_lattice(res, float(c), float(d),
                                        estiramiento=e,
                                        irregularidad=irregularidad, R=R,
                                        seed=seed if s is None else int(s))
        return BW

    def medir(BW, completo=False):
        m = morfometria(BW, sp_cand)
        if forma:
            m.update(metricas_forma(BW, sp_cand, **forma))
        if completo:
            if caro:
                m.update(metricas_forma(BW, sp_cand, **caro))
            if usa_mec:
                m.update(metricas_mecanicas(BW, sp_cand, res_mec, E_s, nu_s))
        return m

    def err(ma, mb, con_mec=False):
        return error_morfometrico(
            ma, mb, peso_mecanico=peso_mecanico if con_mec else 0.0,
            peso_da2=peso_da2, pesos=pesos, distancia=distancia)

    def pr_extra(i, n, etapa):
        if progreso:
            progreso(contador["i"], n_total, f"{etapa} ({i + 1}/{n})")

    def evaluar(d, c, e, etapa):
        contador["i"] += 1
        if progreso:
            progreso(contador["i"], n_total, etapa)
        m = medir(generar(d, c, e, R_fit))
        er, n_us = err(m_voi, m)
        traza.append({"etapa": etapa, "dens": float(d), "celdas": float(c),
                      "estiramiento": [float(x) for x in e],
                      "err": float(er), "n_terminos": int(n_us)})
        return m, er

    def finalistas_de_traza(n):
        vistos, fin = set(), []
        for t in sorted(traza, key=lambda t: t["err"]):
            clave = (round(t["dens"], 9), round(t["celdas"], 9),
                     tuple(t["estiramiento"]))
            if clave in vistos:
                continue
            vistos.add(clave)
            fin.append(t)
            if len(fin) >= int(n):
                break
        return fin

    # ---------- ETAPA A: densidad x celdas ----------
    best = {"err": np.inf, "d": dens[0], "c": celdas[0], "e": e0, "m": None}
    for d in dens:
        for c in celdas:
            m, er = evaluar(d, c, e0, "A")
            if er < best["err"]:
                best.update(err=er, d=d, c=c, m=m)

    # ---------- ETAPA B: presets de estiramiento ----------
    eB, errB, mB = [], [], []
    for e in presets:
        m, er = evaluar(best["d"], best["c"], e, "B")
        eB.append(tuple(e)); errB.append(er); mB.append(m)
    e_best, err_best, m_best, info_k2 = _desempatar_theta(
        eB, errB, mB, best["err"], best["e"], best["m"], m_voi)
    best.update(err=err_best, e=tuple(e_best), m=m_best)

    # ---------- ETAPA C: refinado ----------
    if refinar:
        d_ref = np.unique(np.clip(best["d"] * np.array([0.92, 1.0, 1.08]),
                                  RHO_MIN, rho_max))
        c_ref = np.unique(np.round(np.clip(
            best["c"] * np.array([0.9, 1.0, 1.1]), CELDAS_MIN, c_max), 2))
        for d in d_ref:
            for c in c_ref:
                m, er = evaluar(d, c, best["e"], "C")
                if er < best["err"]:
                    best.update(err=er, d=d, c=c, m=m)

    # ---------- ETAPA D: desempate entre finalistas ----------
    traza_mec = []
    m_voi_mec = {}
    finalistas = []
    if etapa_d:
        if progreso:
            progreso(contador["i"], n_total, "D")
        m_voi = dict(m_voi)
        if usa_mec:
            m_voi_mec = metricas_mecanicas(VOI, spacing, res_mec, E_s, nu_s)
            m_voi.update({k: v for k, v in m_voi_mec.items()
                          if k in CLAVES_MECANICAS})
        if caro:
            m_voi.update(metricas_forma(VOI, spacing, **caro))
        finalistas = finalistas_de_traza(n_finalistas)
        mejor_mec = None
        for k, t in enumerate(finalistas):
            if progreso:
                progreso(contador["i"], n_total, f"D ({k+1}/{len(finalistas)})")
            mc = medir(generar(t["dens"], t["celdas"], t["estiramiento"], R_fit),
                       completo=True)
            err_c, n_us = err(m_voi, mc, con_mec=True)
            traza_mec.append({
                "orden_morfo": k, "err_morfo": t["err"],
                "err_con_mec": float(err_c), "n_terminos": int(n_us),
                "dens": t["dens"], "celdas": t["celdas"],
                "estiramiento": t["estiramiento"],
                "Ez_rel": mc.get("Ez_rel"), "Ez_Ex": mc.get("Ez_Ex"),
                "homogeneizo": bool(usa_mec and "Ez_rel" in mc)})
            if mejor_mec is None or err_c < mejor_mec[0]:
                mejor_mec = (err_c, t, mc)
        if mejor_mec is not None:
            _, t, mc = mejor_mec
            best.update(err=mejor_mec[0], d=t["dens"], c=t["celdas"],
                        e=tuple(t["estiramiento"]), m=mc)

    # ---------- seleccion robusta (U1) ----------
    info_sel = {"criterio": seleccion, "k": float(k_robusto), "K": int(K)}
    if seleccion == "robusta":
        base_fin = finalistas or finalistas_de_traza(n_finalistas)
        idx, tabla = incertidumbre.seleccion_robusta(
            base_fin,
            lambda t, s: generar(t["dens"], t["celdas"], t["estiramiento"],
                                 R_fit, s),
            lambda BW: medir(BW, completo=etapa_d),
            lambda a, b: err(a, b, con_mec=etapa_d),
            m_voi, seed, K, k_robusto, progreso=pr_extra)
        for f, t in zip(tabla, base_fin):
            f.update(dens=t["dens"], celdas=t["celdas"],
                     estiramiento=t["estiramiento"])
        t = base_fin[idx]
        antes = (float(best["d"]), float(best["c"]),
                 tuple(float(x) for x in best["e"]))
        despues = (float(t["dens"]), float(t["celdas"]),
                   tuple(float(x) for x in t["estiramiento"]))
        info_sel.update(tabla=tabla, elegido=idx, cambio=antes != despues)
        best.update(err=tabla[idx]["error"]["media"], d=t["dens"],
                    c=t["celdas"], e=tuple(t["estiramiento"]),
                    m=None if despues != antes else best["m"])

    # ---------- orientacion: L2 (impuesta) o L1 (eje) ----------
    R_final = R_fit
    info_alin = {"aplicada": False, "metodo": "clasica (desactivada)",
                 "motivo": "alinear_fabrica=False"}
    if alinear_fabrica:
        gen_R = lambda R: generar(best["d"], best["c"], best["e"], R)  # noqa: E731
        if alineacion == "eje":
            R_final, info_alin = alinear_dual_por_fabrica(
                m_voi, gen_R, sp_cand, R_fit)
        else:
            R_final, info_alin = alinear_marco_fabrica(
                m_voi, gen_R, lambda BW: morfometria(BW, sp_cand), R_fit,
                m_clasica=best["m"])

    # ---------- medida final ----------
    BW, _, info_gen = generar_dual_lattice(
        res, float(best["c"]), float(best["d"]), estiramiento=best["e"],
        irregularidad=irregularidad, R=R_final, seed=seed)
    m_final = medir(BW, completo=etapa_d)
    err_final, _ = err(m_voi, m_final, con_mec=True)

    # ---------- incertidumbre del ganador (U1) ----------
    info_inc = {"K": 0}
    if K > 0:
        info_inc = incertidumbre.informe(
            lambda s: generar(best["d"], best["c"], best["e"], R_final, s),
            lambda BWr: medir(BWr, completo=etapa_d),
            lambda a, b: err(a, b, con_mec=True),
            m_voi, seed, K, claves_incertidumbre(pesos, usa_mec),
            progreso=pr_extra)
        info_inc["error_busqueda"] = float(best["err"])
        info_inc["error_realizacion_devuelta"] = float(err_final)

    fp, n_poros = fraccion_poro_conexa(BW)
    fp_voi, n_poros_voi = fraccion_poro_conexa(VOI)
    en_tope_d = bool(abs(best["d"] - min(rho_max, bvtv * 1.3)) < 1e-9)
    en_tope_c = bool(best["c"] >= c_max - 1e-9 or best["c"] <= CELDAS_MIN + 1e-9)

    parametros = {
        "familia": "dual-lattice",
        "densidad": float(best["d"]),
        "celdas": float(best["c"]),
        "estiramiento": [float(x) for x in best["e"]],
        "irregularidad": float(irregularidad),
        "R": R_final.tolist(),
        "R_clasica": R_fit.tolist(),
        "euler_deg": [rx, ry, rz],
        "resolucion": res,
        "spacing_mm": float(sp_cand[0]),
        "semilla": int(seed),
        "rho_max": float(rho_max),
        "peso_da2": float(peso_da2),
        "espesor_puntal_mm": float(info_gen["espesor_puntal"] * lc),
    }

    return {
        "familia": "dual-lattice",
        "procedencia": procedencia.desde_parametros(parametros),
        "parametros": parametros,
        "objetivo": {
            "pesos": pesos, "distancia": distancia,
            "peso_mecanico": float(peso_mecanico), "peso_da2": float(peso_da2),
            "medidas_por_evaluacion": sorted(forma),
            "medidas_finalistas": sorted(caro) + (["mecanica"] if usa_mec
                                                  else []),
            "etapa_d": etapa_d,
        },
        "diagnostico": {
            "en_tope_densidad": en_tope_d,
            # Especifico de esta familia: si el ganador esta en un extremo de
            # la rejilla de celdas, la estimacion desde Tb.N se quedo corta o
            # la resolucion no permite celdas mas finas.
            "en_tope_celdas": en_tope_c,
            "celdas_estimadas": float(c_est),
            "poro_conexo_spin": fp,
            "n_poros_spin": n_poros,
            "poro_conexo_voi": fp_voi,
            "n_poros_voi": n_poros_voi,
            "bicontinuo": bool(fp >= 0.9),
            "aviso": ("El candidato NO es bicontinuo: la fase poro esta "
                      "fragmentada." if fp < 0.9 else ""),
        },
        "alineacion": info_alin,
        "error": float(err_final),
        "error_busqueda": float(best["err"]),
        "incertidumbre": info_inc,
        "seleccion": info_sel,
        "metricas_voi": m_voi,
        # La clave conserva el nombre de la del spinodoide porque la leen la
        # interfaz y `metodos.comparar_resultados`; `metricas_candidato` es el
        # mismo dict con un nombre que no miente.
        "metricas_spin": m_final,
        "metricas_candidato": m_final,
        "mascara": BW,
        "rejillas": {
            "densidad": [float(x) for x in dens],
            "celdas": [float(x) for x in celdas],
            "presets_estiramiento": [[float(x) for x in p] for p in presets],
            "n_eval_previstas": int(n_total),
        },
        "k2": info_k2,
        "mecanico": {
            "peso": float(peso_mecanico), "res_mec": int(res_mec),
            "voi": m_voi_mec, "traza": traza_mec,
            "reordeno": bool(traza_mec and min(
                range(len(traza_mec)),
                key=lambda i: traza_mec[i]["err_con_mec"]) != 0),
        },
        "traza": traza,
        "n_evaluaciones": contador["i"],
        "tiempo_s": time.time() - t0,
    }
