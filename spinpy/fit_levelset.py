"""
fit_levelset.py — Ajuste de level-set-foam a un VOI por busqueda escalonada.

Misma estructura que `fit_dual.ajustar_dual_lattice` (el contrato del resultado,
el objetivo O1, la incertidumbre U1, la alineacion) con la busqueda validada
sobre cuatro VOI porcinos (`docs/familias/LEVEL_SET_FOAM.md`). Se escribe aparte
para no tocar los ajustes de las otras dos familias, contrastados contra
`AppFinal_V2.m` y entre si.

LAS ETAPAS, FRENTE A LAS DE LAS OTRAS DOS FAMILIAS
--------------------------------------------------
    etapa   spinodoide              dual-lattice            level-set-foam
    A       densidad x onda         densidad x celdas       tamano (celdas)
    B       presets de conos        presets estiramiento    presets estiramiento
    C       refinado                refinado                refinado del tamano
    D       desempate               desempate               FORMA: K_c x nodos
    E       -                       -                       reajuste del tamano
    F       -                       -                       afinado de forma

LA DENSIDAD NO SE BUSCA. morphoGen resuelve el espesor del puntal para la
fraccion volumetrica pedida, asi que se fija al BV/TV medido del VOI. Dejarla
libre (se probo) hace que el ajuste la sacrifique hasta un 15 % para mejorar
los terminos de forma, y la rigidez empeora (docs/familias/LEVEL_SET_FOAM.md,
ronda con densidad libre). Es la desviacion principal frente a las otras dos
familias y es deliberada.

EL OBJETIVO ES EL AMPLIADO
--------------------------
`error.PESOS_LEVEL_SET`: los terminos base mas PoDm, TbTh_CV, H_desv y Conn.D
(TbSp a 0), con DA2. Los terminos de forma del VOI se miden al MISMO tamano de
voxel que el candidato (el VOI se submuestrea al lado del candidato): la
curvatura, el espesor local, el poro y la conectividad dependen del voxel, y
compararlos a resoluciones distintas mide la resolucion, no la estructura.

LA ORIENTACION SE RESUELVE PERMUTANDO EJES
------------------------------------------
Esta familia no admite rotacion. El objetivo no ve la orientacion (DA es un
escalar), de modo que la busqueda elige un eje alargado cualquiera entre los
equivalentes. Al terminar, se prueban las permutaciones del estiramiento y se
queda la cuyo tensor de fabrica (MIL) es mas parecido al del VOI. Es exacto, no
interpola, y solo se aplica si el VOI tiene anisotropia medible (DA > 1.10).
"""

from __future__ import annotations

import itertools
import time

import numpy as np

from . import incertidumbre, procedencia
from .elastic import remuestrear_bw
from .error import (PESOS_LEVEL_SET, error_morfometrico, fabric_cos,
                    medidas_necesarias)
from .fit import (CLAVES_MECANICAS, E_S_DEF, NU_S_DEF, REPLICAS_DEF,
                  SEMILLA_POR_DEFECTO, claves_incertidumbre, fraccion_poro_conexa,
                  longitud_caracteristica, metricas_mecanicas,
                  resolucion_comparacion, validar_opciones)
from .fit_dual import PESO_DA2, PRESETS_COMPLETO, PRESETS_RAPIDO
from .level_set import es_viable, generar_level_set
from .morphometry import metricas_forma, morfometria

# Rejillas de la busqueda (las de la validacion).
KC = (-1.0, -0.7, -0.35, 0.0, 0.3)
NODOS = ((0.0, 0.0), (10.0, 0.5), (50.0, 1.0), (150.0, 2.0))
KC_RAPIDO = (-0.7, 0.0)
NODOS_RAPIDO = ((0.0, 0.0), (10.0, 0.5))

# Celdas por (trabeculas por lado)^3, MEDIDO sobre cuatro VOI porcinos: el
# optimo estuvo entre 122 y 211 semillas con Tb.N*lc de ~9 a ~11. Solo ancla la
# rejilla; la busqueda la abre de 0.37 a 6 veces ese valor.
CELDAS_POR_TRAB3 = 0.18
CELDAS_MIN = 24.0
FACTORES_CELDAS = (0.37, 6.0)
VOX_POR_CELDA_MIN = 6.0

INVIABLE = 1e6


def celdas_estimadas(m_voi, lc):
    """Numero de celdas que corresponde al Tb.N del VOI en un cubo de lado lc."""
    tbn = float(m_voi.get("TbN", np.nan))
    if not np.isfinite(tbn) or tbn <= 0:
        return 134.0
    return CELDAS_POR_TRAB3 * (tbn * float(lc)) ** 3


def _tensor_fabrica(m):
    """Tensor MIL normalizado a traza 1, o None si no es utilizable."""
    if not m.get("MIL_valid", False):
        return None
    w = np.asarray(m["eigenvalues"], float)
    V = np.asarray(m["eigenvectors"], float)
    T = V @ np.diag(w) @ V.T
    tr = float(np.trace(T))
    return T / tr if tr > 0 else None


def alinear_por_permutacion(m_voi, generar, medir, e_mejor, da_min=1.10):
    """Permutacion del estiramiento cuyo tensor de fabrica se parece al del VOI.

    `generar(e)` da la mascara y `medir(BW)` su morfometria. Devuelve
    (estiramiento, info). No aplica nada si el VOI no tiene anisotropia
    medible: por debajo de DA 1.10 su direccion principal es ruido (suelo del
    MIL ~1.07).
    """
    e0 = tuple(float(x) for x in e_mejor)
    info = {"aplicada": False, "metodo": "permutacion de ejes"}
    T_voi = _tensor_fabrica(m_voi)
    da = float(m_voi.get("DA", np.nan))
    if T_voi is None or not np.isfinite(da) or da < da_min:
        info["motivo"] = ("el VOI no tiene anisotropia medible (DA < "
                          f"{da_min}): la direccion principal no significa nada")
        return e0, info
    perms = sorted(set(itertools.permutations(e0)))
    if len(perms) == 1:
        info["motivo"] = "estiramiento isotropo: todas las permutaciones coinciden"
        return e0, info
    filas = []
    for p in perms:
        m = medir(generar(p))
        T = _tensor_fabrica(m)
        if T is None:
            continue
        filas.append((float(np.linalg.norm(T - T_voi)), p,
                      float(np.degrees(np.arccos(
                          min(1.0, max(0.0, fabric_cos(m_voi, m))))))))
    if not filas:
        info["motivo"] = "ninguna permutacion dio un MIL valido"
        return e0, info
    d_ant = next((d for d, p, _a in filas if p == e0), None)
    ang_ant = next((a for d, p, a in filas if p == e0), None)
    d, p, ang = min(filas)
    info.update(aplicada=bool(p != e0), permutacion=list(p),
                distancia_antes=d_ant, distancia_despues=d,
                angulo_antes_deg=ang_ant, angulo_despues_deg=ang,
                n_probadas=len(filas))
    return tuple(p), info


def ajustar_level_set(VOI, spacing, modo="completo", m_voi=None,
                      seed=SEMILLA_POR_DEFECTO, resolucion=None, progreso=None,
                      peso_mecanico=0.0, res_mec=16, n_finalistas=5,
                      E_s=E_S_DEF, nu_s=NU_S_DEF, alinear_fabrica=True,
                      peso_da2=PESO_DA2, pesos=None, distancia="wasserstein",
                      replicas=REPLICAS_DEF, seleccion="minimo",
                      k_robusto=1.0, alineacion="impuesta"):
    """Busca los parametros de level-set-foam que mejor reproducen el VOI.

    modo : 'completo' (A+B+C+D+E+F, ~60 evaluaciones) | 'rapido' (A+B+D
           reducida, ~15)

    pesos : por omision `error.PESOS_LEVEL_SET`. Se puede pasar `{}` para el
        objetivo base de las otras familias (no recomendado: el error queda
        casi saturado y no separa estructuras).

    El resto de argumentos significan lo mismo que en `ajustar_dual_lattice`.
    `alineacion` se acepta por simetria con las otras familias; aqui la unica
    alineacion es la permutacion de ejes.
    """
    if modo not in ("completo", "rapido"):
        raise ValueError(f"modo: 'completo' o 'rapido', no {modo!r}")
    if pesos is None:
        pesos = dict(PESOS_LEVEL_SET)
    pesos, K = validar_opciones(pesos, distancia, replicas, seleccion, "impuesta")
    nec = medidas_necesarias(pesos)
    forma = {k: True for k in nec["forma"]}
    caro = {k: True for k in nec["caro"]}
    usa_mec = bool(peso_mecanico > 0 or nec["mecanico"])
    etapa_d = bool(usa_mec or caro)

    VOI = np.asarray(VOI, dtype=bool)
    spacing = np.asarray(spacing, dtype=float).ravel()
    if m_voi is None:
        m_voi = morfometria(VOI, spacing)
    lc = longitud_caracteristica(VOI.shape, spacing)
    res = int(resolucion) if resolucion else resolucion_comparacion(VOI.shape)
    sp_cand = np.full(3, lc / res)
    bvtv = float(m_voi["BVTV"])

    # Terminos de forma del VOI al voxel del candidato.
    VOI_eq, sp_eq = remuestrear_bw(VOI, spacing, res)
    m_voi = dict(m_voi)
    if forma:
        m_voi.update(metricas_forma(VOI_eq, sp_eq, **forma))

    n0 = celdas_estimadas(m_voi, lc)
    n_max = max(CELDAS_MIN * 2, (res / VOX_POR_CELDA_MIN) ** 3)
    e0 = (1.0, 1.0, 1.0)
    f_lo, f_hi = FACTORES_CELDAS
    if modo == "rapido":
        tamanos = n0 * np.geomspace(0.6, 3.0, 4)
        presets, kcs, nodos = PRESETS_RAPIDO, KC_RAPIDO, NODOS_RAPIDO
    else:
        tamanos = n0 * np.geomspace(f_lo, f_hi, 6)
        presets, kcs, nodos = PRESETS_COMPLETO, KC, NODOS
    tamanos = np.unique(np.round(np.clip(tamanos, CELDAS_MIN, n_max)))

    n_total = (len(tamanos) + len(presets) + len(kcs) * len(nodos)
               + (0 if modo == "rapido" else 3 + 4 + 4))
    traza, cache = [], {}
    t0 = time.time()
    contador = {"i": 0}

    def generar(N, e, kc, c1, c2, s=None):
        BW, _, _ = generar_level_set(
            res, float(N), bvtv, estiramiento=tuple(e), curvatura=float(kc),
            c1=float(c1), c2=float(c2), seed=seed if s is None else int(s))
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

    best = {"err": np.inf, "N": tamanos[0], "e": e0, "kc": 0.0, "c1": 0.0,
            "c2": 0.0, "m": None}

    def evaluar(N, e, kc, c1, c2, etapa):
        contador["i"] += 1
        if progreso:
            progreso(contador["i"], n_total, etapa)
        clave = (round(float(N), 6), tuple(round(float(x), 9) for x in e),
                 round(float(kc), 9), round(float(c1), 9), round(float(c2), 9))
        if clave in cache:
            return cache[clave]
        par = {"etapa": etapa, "N": float(N), "estiramiento": [float(x) for x in e],
               "curvatura": float(kc), "c1": float(c1), "c2": float(c2)}
        if not es_viable(N, e):
            traza.append(dict(par, err=INVIABLE, n_terminos=0, inviable=True))
            cache[clave] = (None, INVIABLE)
            return cache[clave]
        try:
            m = medir(generar(N, e, kc, c1, c2))
        except (ValueError, ArithmeticError):
            traza.append(dict(par, err=INVIABLE, n_terminos=0, inviable=True))
            cache[clave] = (None, INVIABLE)
            return cache[clave]
        er, n_us = err(m_voi, m)
        traza.append(dict(par, err=float(er), n_terminos=int(n_us)))
        cache[clave] = (m, float(er))
        if er < best["err"]:
            best.update(err=float(er), N=N, e=tuple(e), kc=kc, c1=c1, c2=c2, m=m)
        return cache[clave]

    def estado():
        return (best["N"], best["e"], best["kc"], best["c1"], best["c2"])

    # ---------- ETAPA A: tamano ----------
    for N in tamanos:
        evaluar(N, e0, 0.0, 0.0, 0.0, "A")
    # ---------- ETAPA B: presets de estiramiento ----------
    for e in presets:
        evaluar(best["N"], e, 0.0, 0.0, 0.0, "B")
    if modo == "completo":
        # ---------- ETAPA C: refinado del tamano ----------
        for f in (0.8, 1.0, 1.25):
            evaluar(round(best["N"] * f), best["e"], 0.0, 0.0, 0.0, "C")
    # ---------- ETAPA D: forma (curvatura x engrosamiento de nodos) ----------
    N_b, e_b = best["N"], best["e"]
    for kc in kcs:
        for c1, c2 in nodos:
            evaluar(N_b, e_b, kc, c1, c2, "D")
    if modo == "completo":
        # ---------- ETAPA E: reajuste del tamano con la forma elegida ----------
        for f in (0.75, 0.9, 1.1, 1.3):
            evaluar(round(best["N"] * f), best["e"], best["kc"], best["c1"],
                    best["c2"], "E")
        # ---------- ETAPA F: afinado de la forma ----------
        kc0, c10, c20 = best["kc"], best["c1"], best["c2"]
        for dk in (-0.15, 0.15):
            if kc0 + dk < 0.95:
                evaluar(best["N"], best["e"], kc0 + dk, c10, c20, "F")
        for fn in (0.6, 1.6):
            if c10 > 0:
                evaluar(best["N"], best["e"], best["kc"], c10 * fn, c20, "F")
            else:
                evaluar(best["N"], best["e"], best["kc"], 5.0 * fn, 0.5, "F")

    if best["m"] is None:
        raise RuntimeError("level-set-foam: ninguna combinacion fue viable para "
                           "este VOI (revise la resolucion y el tamano del VOI)")

    def finalistas_de_traza(n):
        vistos, fin = set(), []
        for t in sorted((t for t in traza if not t.get("inviable")),
                        key=lambda t: t["err"]):
            clave = (round(t["N"], 6), tuple(round(x, 9) for x in t["estiramiento"]),
                     round(t["curvatura"], 9), round(t["c1"], 9), round(t["c2"], 9))
            if clave in vistos:
                continue
            vistos.add(clave)
            fin.append(t)
            if len(fin) >= int(n):
                break
        return fin

    # ---------- desempate mecanico / medidas caras entre finalistas ----------
    traza_mec, m_voi_mec = [], {}
    if etapa_d:
        if progreso:
            progreso(contador["i"], n_total, "D")
        m_voi = dict(m_voi)
        if usa_mec:
            m_voi_mec = metricas_mecanicas(VOI, spacing, res_mec, E_s, nu_s)
            m_voi.update({k: v for k, v in m_voi_mec.items()
                          if k in CLAVES_MECANICAS})
        if caro:
            m_voi.update(metricas_forma(VOI_eq, sp_eq, **caro))
        mejor_mec = None
        for k, t in enumerate(finalistas_de_traza(n_finalistas)):
            mc = medir(generar(t["N"], t["estiramiento"], t["curvatura"],
                               t["c1"], t["c2"]), completo=True)
            err_c, n_us = err(m_voi, mc, con_mec=True)
            traza_mec.append({
                "orden_morfo": k, "err_morfo": t["err"],
                "err_con_mec": float(err_c), "n_terminos": int(n_us),
                "N": t["N"], "estiramiento": t["estiramiento"],
                "Ez_rel": mc.get("Ez_rel"), "Ez_Ex": mc.get("Ez_Ex"),
                "homogeneizo": bool(usa_mec and "Ez_rel" in mc)})
            if mejor_mec is None or err_c < mejor_mec[0]:
                mejor_mec = (err_c, t, mc)
        if mejor_mec is not None:
            _, t, mc = mejor_mec
            best.update(err=mejor_mec[0], N=t["N"], e=tuple(t["estiramiento"]),
                        kc=t["curvatura"], c1=t["c1"], c2=t["c2"], m=mc)

    # ---------- seleccion robusta (U1) ----------
    info_sel = {"criterio": seleccion, "k": float(k_robusto), "K": int(K)}
    if seleccion == "robusta":
        base_fin = finalistas_de_traza(n_finalistas)
        idx, tabla = incertidumbre.seleccion_robusta(
            base_fin,
            lambda t, s: generar(t["N"], t["estiramiento"], t["curvatura"],
                                 t["c1"], t["c2"], s),
            lambda BW: medir(BW, completo=etapa_d),
            lambda a, b: err(a, b, con_mec=etapa_d),
            m_voi, seed, K, k_robusto, progreso=pr_extra)
        for f, t in zip(tabla, base_fin):
            f.update(N=t["N"], estiramiento=t["estiramiento"],
                     curvatura=t["curvatura"], c1=t["c1"], c2=t["c2"])
        t = base_fin[idx]
        cambio = (t["N"], tuple(t["estiramiento"]), t["curvatura"], t["c1"],
                  t["c2"]) != estado()
        info_sel.update(tabla=tabla, elegido=idx, cambio=cambio)
        best.update(err=tabla[idx]["error"]["media"], N=t["N"],
                    e=tuple(t["estiramiento"]), kc=t["curvatura"], c1=t["c1"],
                    c2=t["c2"], m=None if cambio else best["m"])

    # ---------- orientacion: permutacion de ejes ----------
    e_final = best["e"]
    info_alin = {"aplicada": False, "metodo": "permutacion de ejes",
                 "motivo": "alinear_fabrica=False"}
    if alinear_fabrica:
        e_final, info_alin = alinear_por_permutacion(
            m_voi,
            lambda e: generar(best["N"], e, best["kc"], best["c1"], best["c2"]),
            lambda BW: morfometria(BW, sp_cand), best["e"])

    # ---------- medida final ----------
    BW, _, info_gen = generar_level_set(
        res, float(best["N"]), bvtv, estiramiento=e_final,
        curvatura=best["kc"], c1=best["c1"], c2=best["c2"], seed=seed)
    m_final = medir(BW, completo=etapa_d)
    err_final, _ = err(m_voi, m_final, con_mec=True)

    # ---------- incertidumbre del ganador (U1) ----------
    info_inc = {"K": 0}
    if K > 0:
        info_inc = incertidumbre.informe(
            lambda s: generar(best["N"], e_final, best["kc"], best["c1"],
                              best["c2"], s),
            lambda BWr: medir(BWr, completo=etapa_d),
            lambda a, b: err(a, b, con_mec=True),
            m_voi, seed, K, claves_incertidumbre(pesos, usa_mec),
            progreso=pr_extra)
        info_inc["error_busqueda"] = float(best["err"])
        info_inc["error_realizacion_devuelta"] = float(err_final)

    fp, n_poros = fraccion_poro_conexa(BW)
    fp_voi, n_poros_voi = fraccion_poro_conexa(VOI)
    n_inv = sum(1 for t in traza if t.get("inviable"))
    c_extremo = bool(best["N"] >= tamanos.max() - 1e-9
                     or best["N"] <= tamanos.min() + 1e-9)

    I3 = np.eye(3).tolist()
    parametros = {
        "familia": "level-set-foam",
        "densidad": bvtv,
        "n_celdas": float(best["N"]),
        "estiramiento": [float(x) for x in e_final],
        "curvatura": float(best["kc"]),
        "c1": float(best["c1"]),
        "c2": float(best["c2"]),
        # Sin rotacion (ver cabecera): la identidad, para que el codigo que lee
        # `R` de cualquier familia no tenga que distinguir.
        "R": I3, "R_clasica": I3, "euler_deg": [0.0, 0.0, 0.0],
        "resolucion": res,
        "spacing_mm": float(sp_cand[0]),
        "semilla": int(seed),
        "morphogen": info_gen.get("morphogen"),
    }

    return {
        "familia": "level-set-foam",
        "procedencia": procedencia.desde_parametros(parametros),
        "parametros": parametros,
        "objetivo": {
            "pesos": pesos, "distancia": distancia,
            "peso_mecanico": float(peso_mecanico), "peso_da2": float(peso_da2),
            "medidas_por_evaluacion": sorted(forma),
            "medidas_finalistas": sorted(caro) + (["mecanica"] if usa_mec else []),
            "etapa_d": etapa_d,
        },
        "diagnostico": {
            "en_tope_densidad": False,        # la densidad no se busca
            "en_tope_celdas": c_extremo,
            "celdas_estimadas": float(n0),
            "n_inviables": int(n_inv),
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
        "metricas_spin": m_final,
        "metricas_candidato": m_final,
        "mascara": BW,
        "rejillas": {
            "celdas": [float(x) for x in tamanos],
            "presets_estiramiento": [[float(x) for x in p] for p in presets],
            "curvatura": [float(x) for x in kcs],
            "nodos": [[float(a), float(b)] for a, b in nodos],
            "n_eval_previstas": int(n_total),
        },
        "k2": {},
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
