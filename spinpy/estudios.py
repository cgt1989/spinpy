"""
estudios.py: estudios complementarios del informe para publicacion.

POR QUE EXISTE
--------------
El informe describe lo que la sesion calculo. Tres preguntas que un revisor
hace sobre un modelo de micro-EF no tenian respuesta en la sesion:

  1. FORMA. BV/TV, Tb.Th y DA no distinguen placas de barras ni una red de un
     monton de islas. El Ellipsoid Factor (Doube 2015) y el perfil de
     curvaturas de la interfaz (Guo et al. 2024) si.
  2. SEGMENTACION. Hara et al. (2002, Bone 31:107) midieron que un 0,5 % de
     cambio de umbral mueve el BV/TV un 5 % y la rigidez un 9 % cuando
     BV/TV < 0,15, y un 2 % y un 3 % por encima de 0,2. El VOI llega casi
     siempre ya binario, asi que el umbral no se puede volver a elegir; lo que
     si se puede es mover la superficie como la moveria un cambio de umbral.
  3. CONDICIONES DE CONTORNO. Pahr y Zysset (2008, Biomech Model Mechanobiol
     7:463) mostraron que en muestras de menos de ~5 mm de lado las
     condiciones de desplazamiento uniforme pueden sobrestimar la rigidez
     efectiva. spinpy tiene tres: el apoyo deslizante y el empotrado del
     ensayo, y la condicion periodica de la homogeneizacion.

Cada funcion devuelve un diccionario serializable en JSON, y `completar` los
guarda en `doc["resultados"]` con las claves que lee el informe:

  forma                     EF y curvaturas por estructura
  sensibilidad_superficie   metricas y E_app frente al desplazamiento
  condiciones_contorno      E/E_s con cada condicion, misma rejilla

EL SUSTITUTO DEL UMBRAL, CON SUS LIMITES
----------------------------------------
Sobre la mascara binaria se aplica una gaussiana de sigma = 1 voxel, que
emula el emborronamiento del escaner, y se umbraliza a t. Para una interfaz
plana, umbralizar a t equivale a desplazarla

    delta = -sigma * Phi^-1(t)        (Phi: normal estandar acumulada)

hacia el vacio si delta > 0 (engrosa) o hacia el hueso si delta < 0. Por eso
los niveles se piden como desplazamientos en voxeles y se convierten a t. En
puntales de pocos voxeles la interfaz no es plana y la equivalencia es solo
aproximada; el informe lo declara. Con t = 0,5 se recupera la estructura
original salvo el redondeo de las aristas.

Todos los caminos de medida son los de la app: `morfometria`,
`remuestrear_bw` + `ensayo_compresion_eje` y `homogeneizar`. Nada aqui mide
por su cuenta.
"""

from __future__ import annotations

import time

import numpy as np

from .resistencia import E_S_DEF, NU_S_DEF

# Desplazamientos de la superficie, en voxeles de la imagen original.
DESPLAZAMIENTOS = (-1.0, -0.5, 0.0, 0.5, 1.0)
SIGMA_VOX = 1.0
N_MEC_DEF = 40
NBINS_CURVATURA = 64
NBINS_EF = 40
# Ellipsoid Factor: el coste crece con el volumen (un VOI de 188^3 serian
# horas en 4 nucleos). Es una medida LOCAL de forma, asi que se calcula en un
# subcubo central a la resolucion original y no en el volumen reducido, que
# adelgazaria los puntales a pocos voxeles. El lado del subcubo es el MISMO en
# milimetros para todas las estructuras: EF_SUBCUBO_VOX voxeles de la rejilla
# mas gruesa. Medido sobre el spinodoide porcino (96^3): 64^3, 3000 semillas,
# 100 evaluaciones y 4 rondas cubren el 96,5 % del solido en 154 s.
EF_SUBCUBO_VOX = 64
EF_PRESUPUESTO = {"max_semillas": 3000, "max_evaluaciones": 100, "rondas": 4,
                  "relleno_objetivo": 0.95}
# Vertices a menos de este numero de voxeles de las caras del cubo no
# cuentan en la curvatura: alli la malla esta cortada (ver curvatura.py).
MARGEN_CURVATURA = 2
APOYOS_ENSAYO = ("deslizante", "empotrado")


def _num(x):
    """float finito o None: lo que se puede escribir en JSON sin NaN."""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def _sp3(spacing):
    sp = np.atleast_1d(np.asarray(spacing, float)).ravel()
    return np.repeat(sp, 3) if sp.size == 1 else sp


def _progreso(cb, frac, msg):
    if cb is not None:
        cb(float(frac), msg)


# ---------------------------------------------------------------------------
# 1. Forma: Ellipsoid Factor y curvaturas
# ---------------------------------------------------------------------------

def subcubo_central(BW, spacing, lado_mm):
    """Cubo central de `lado_mm` (o el volumen entero si es menor)."""
    BW = np.asarray(BW, bool)
    sp = _sp3(spacing)
    m = [min(n, max(2, int(round(lado_mm / h)))) for n, h in zip(BW.shape, sp)]
    o = [(n - k) // 2 for n, k in zip(BW.shape, m)]
    return BW[o[0]:o[0] + m[0], o[1]:o[1] + m[1], o[2]:o[2] + m[2]]


def ellipsoid_factor(BW, spacing, nbins=NBINS_EF, lado_mm=None, **kw):
    """Resumen del EF y su histograma sobre [-1, 1], ponderado por voxel.

    Con `lado_mm`, sobre el subcubo central de ese lado (ver EF_SUBCUBO_VOX).
    """
    from .elipsoide import factor_elipsoide
    BW = np.asarray(BW, bool)
    if lado_mm is not None:
        BW = subcubo_central(BW, spacing, lado_mm)
    opciones = dict(EF_PRESUPUESTO)
    opciones.update(kw)
    r = factor_elipsoide(BW, _sp3(spacing), campo=True, **opciones)
    campo = r.pop("EF_campo", None)
    out = {k: _num(v) for k, v in r.items()}
    out["n_subcubo"] = int(BW.shape[0])
    out["lado_subcubo_mm"] = float(BW.shape[0] * _sp3(spacing)[0])
    if campo is not None:
        v = np.asarray(campo, float)
        v = v[np.isfinite(v)]
        h, bordes = np.histogram(v, bins=nbins, range=(-1.0, 1.0))
        s = h.sum()
        out["histograma"] = (h / s).tolist() if s else h.tolist()
        out["bordes"] = bordes.tolist()
    return out


def curvaturas(BW, spacing, sigma_vox=SIGMA_VOX, margen=MARGEN_CURVATURA):
    """k1, k2 y areas de la interfaz de la mascara suavizada.

    Se mide sobre la mascara suavizada con una gaussiana y no sobre la
    binaria: la curvatura de una escalera de voxeles es la de las esquinas,
    no la de la trabecula (ver curvatura.py). VOI y candidatos pasan por el
    mismo suavizado, que es lo que hace comparables sus perfiles.
    """
    from scipy import ndimage
    from .curvatura import curvaturas_de_campo
    # `curvaturas_de_campo` orienta la normal hacia donde el campo CRECE, y
    # la convencion del proyecto es que apunte al vacio. La mascara suavizada
    # crece hacia el hueso, asi que se pasa su complemento; si no, una bola de
    # hueso saldria concava (bloque 28).
    campo = 1.0 - ndimage.gaussian_filter(np.asarray(BW, float), sigma_vox)
    return curvaturas_de_campo(campo, _sp3(spacing), 0.5, margen=margen)


def perfiles_curvatura(curvas, nbins=NBINS_CURVATURA):
    """Perfil (k1, k2) de cada estructura con un LIMITE COMUN, y su resumen.

    El limite comun (el mayor de los percentiles 99 de |k|) hace que los
    perfiles se lean sobre los mismos ejes: con un limite por estructura, dos
    manchas identicas parecerian distintas.
    """
    from .curvatura import perfil, resumen
    lim = 0.0
    for c in curvas.values():
        v = np.abs(np.r_[c["k1"], c["k2"]])
        v = v[np.isfinite(v)]
        if v.size:
            lim = max(lim, float(np.percentile(v, 99)))
    lim = lim or 1.0
    out = {}
    for est, c in curvas.items():
        P, _bx, _by = perfil(c["k1"], c["k2"], c["areas"], limite=lim,
                             nbins=nbins)
        out[est] = {"resumen": {k: _num(v) for k, v in
                                resumen(c["k1"], c["k2"], c["areas"]).items()},
                    "perfil": np.round(P, 7).tolist(), "limite": lim,
                    "n_vertices": int(c["n_vertices"]),
                    "n_descartados": int(c["n_descartados"])}
    return out


# ---------------------------------------------------------------------------
# 2. Sensibilidad a la posicion de la superficie (sustituto del umbral)
# ---------------------------------------------------------------------------

def nivel_de_desplazamiento(delta_vox, sigma_vox=SIGMA_VOX):
    """Umbral t sobre la mascara suavizada que desplaza `delta_vox` voxeles
    una interfaz plana: t = Phi(-delta / sigma)."""
    from scipy.special import ndtr
    return float(ndtr(-float(delta_vox) / float(sigma_vox)))


def desplazar_superficie(BW, delta_vox, sigma_vox=SIGMA_VOX):
    """Mascara con la superficie movida `delta_vox` voxeles (> 0 engrosa)."""
    from scipy import ndimage
    campo = ndimage.gaussian_filter(np.asarray(BW, float), sigma_vox)
    return campo > nivel_de_desplazamiento(delta_vox, sigma_vox)


def sensibilidad_superficie(BW, spacing, desplazamientos=DESPLAZAMIENTOS,
                            n_mec=N_MEC_DEF, eje=2, E_s=E_S_DEF,
                            nu_s=NU_S_DEF, apoyo="deslizante",
                            sigma_vox=SIGMA_VOX, progreso=None):
    """Morfometria y E_app de la estructura con la superficie desplazada.

    Cada fila lleva el desplazamiento en voxeles y en micrometros, el umbral
    equivalente t, BV/TV, Tb.Th, Tb.Sp, Tb.N, DA y E_app a `n_mec`^3. El
    ensayo usa el mismo remuestreo y el mismo solver que el de la app.
    """
    from .elastic import remuestrear_bw
    from .morphometry import morfometria
    from .resistencia import ensayo_compresion_eje
    sp = _sp3(spacing)
    filas = []
    for k, d in enumerate(desplazamientos):
        _progreso(progreso, k / len(desplazamientos),
                  f"superficie desplazada {d:+g} voxeles")
        bw = desplazar_superficie(BW, d, sigma_vox)
        m = morfometria(bw, sp, do_mil=True)
        fila = {"desplazamiento_vox": float(d),
                "desplazamiento_um": float(d * sp[0] * 1000.0),
                "t": nivel_de_desplazamiento(d, sigma_vox)}
        fila.update({c: _num(m.get(c)) for c in ("BVTV", "TbTh", "TbSp",
                                                 "TbN", "DA")})
        bwr, spr = remuestrear_bw(bw, sp, n_mec)
        try:
            r = ensayo_compresion_eje(bwr, spr, eje=eje, E_s=E_s, nu_s=nu_s,
                                      apoyo=apoyo)
            fila.update(E_app=_num(r.get("E_app")), ok=bool(r.get("ok")),
                        residuo=_num(r.get("residuo_rel")))
        except Exception as e:                      # se anota, no se oculta
            fila.update(E_app=None, ok=False, error=f"{type(e).__name__}: {e}")
        filas.append(fila)
    return {"filas": filas,
            "parametros": {"sigma_vox": float(sigma_vox), "n_mec": int(n_mec),
                           "eje": int(eje), "E_s_Pa": float(E_s),
                           "nu_s": float(nu_s), "apoyo": apoyo}}


# ---------------------------------------------------------------------------
# 3. Condiciones de contorno sobre la misma rejilla
# ---------------------------------------------------------------------------

def condiciones_contorno(BW, spacing, n_mec=N_MEC_DEF, eje=2, E_s=E_S_DEF,
                         nu_s=NU_S_DEF, C_periodica=None, progreso=None):
    """E/E_s a lo largo de `eje` con cada condicion de contorno de spinpy.

    deslizante  ensayo con caras cargadas libres de deslizar
    empotrado   ensayo con caras cargadas sin desplazamiento lateral
    periodica   homogeneizacion periodica; E del eje del tensor C

    Las tres sobre la mascara remuestreada a `n_mec`^3, para que la unica
    diferencia sea la condicion de contorno. `C_periodica` (6x6, Pa) es el
    tensor que la sesion ya homogeneizo A ESA MISMA REJILLA: la
    homogeneizacion es con mucho lo mas caro (255 s a 32^3 en 4 nucleos) y
    repetirla no anade nada. Sin el, se homogeneiza aqui.
    """
    from .elastic import constantes_ingenieria, homogeneizar, remuestrear_bw
    from .resistencia import ensayo_compresion_eje
    bwr, spr = remuestrear_bw(np.asarray(BW, bool), _sp3(spacing), n_mec)
    out = {"n": int(bwr.shape[0]), "eje": int(eje),
           "BVTV_rejilla": float(bwr.mean())}
    for k, apoyo in enumerate(APOYOS_ENSAYO):
        _progreso(progreso, k / 3, f"ensayo con apoyo {apoyo}")
        try:
            r = ensayo_compresion_eje(bwr, spr, eje=eje, E_s=E_s, nu_s=nu_s,
                                      apoyo=apoyo)
            out[apoyo] = {"E_rel": _num(r["E_app"] / E_s) if r.get("ok")
                          else None, "ok": bool(r.get("ok")),
                          "residuo": _num(r.get("residuo_rel"))}
        except Exception as e:
            out[apoyo] = {"E_rel": None, "ok": False,
                          "error": f"{type(e).__name__}: {e}"}
    if C_periodica is not None:
        Ch, info = np.asarray(C_periodica, float), {"ok": True,
                                                    "origen": "sesion"}
    else:
        _progreso(progreso, 2 / 3, "homogeneizacion periodica")
        Ch, info = homogeneizar(bwr, E_s=E_s, nu_s=nu_s,
                                vox_size=float(spr[0]))
        info["origen"] = "calculada"
    E_eje = None
    if info.get("ok") and np.all(np.isfinite(Ch)):
        ec = constantes_ingenieria(Ch)
        E_eje = _num(ec.get(("Ex", "Ey", "Ez")[int(eje)]))
    out["periodica"] = {"E_rel": None if E_eje is None else E_eje / E_s,
                        "ok": E_eje is not None,
                        "origen": info.get("origen"),
                        "residuo": _num(info.get("residuo_rel"))}
    out["parametros"] = {"E_s_Pa": float(E_s), "nu_s": float(nu_s)}
    return out


# ---------------------------------------------------------------------------
# Todo junto, sobre las estructuras de una sesion
# ---------------------------------------------------------------------------

ESTUDIOS = ("forma", "superficie", "contorno", "perdida", "fallo")


def _tensor_de_sesion(res, codigo):
    """(C en Pa, n, E_s, nu_s) del tensor que la sesion homogeneizo para la
    estructura `codigo`, o (None,)*4 si no hay."""
    for fam, suf in (("spinodoide", ""), ("dual-lattice", "_dual")):
        rec = res.get("elastico" + suf)
        if not isinstance(rec, dict):
            continue
        for k, v in (rec.get("por_estructura") or {}).items():
            est = "voi" if k == "voi" else fam
            C = (v or {}).get("C_Pa")
            if est != codigo or C is None or not rec.get("resolucion"):
                continue
            C = np.asarray(C, float)
            if C.shape == (6, 6) and np.all(np.isfinite(C)):
                return (C, int(rec["resolucion"]),
                        float(rec.get("E_s_Pa") or E_S_DEF),
                        float(rec.get("nu_s") or NU_S_DEF))
    return None, None, None, None


def completar(doc, estructuras, que=("forma", "superficie", "contorno"),
              n_mec=None, rehacer=False, progreso=None):
    """Calcula los estudios que falten y los guarda en `doc["resultados"]`.

    `estructuras`: lista de (codigo, BW, spacing) con codigo "voi",
    "spinodoide" o "dual-lattice", la misma que prepara el informe.
    `n_mec`: rejilla de los ensayos; por omision, la del ensayo de la sesion
    o N_MEC_DEF. Lo que ya existe en el documento no se recalcula salvo con
    `rehacer`. Devuelve la lista de estudios calculados.

    "perdida" y "fallo" (simulaciones in silico) se corren sobre el VOI si
    esta y, si no, sobre la primera estructura: son las mismas funciones que
    usa el visor.
    """
    res = doc.setdefault("resultados", {})
    if n_mec is None:
        n_mec = int((res.get("resistencia") or {}).get("resolucion")
                    or N_MEC_DEF)
    E_s = float((res.get("resistencia") or {}).get("E_s_Pa") or E_S_DEF)
    nu_s = float((res.get("resistencia") or {}).get("nu_s") or NU_S_DEF)
    ests = [(c, np.asarray(b, bool), _sp3(s)) for c, b, s in estructuras]
    hechos = []
    t0 = time.time()

    def hace_falta(clave):
        return rehacer or not isinstance(res.get(clave), dict)

    if "forma" in que and ests and hace_falta("forma"):
        _progreso(progreso, 0.0, "forma: Ellipsoid Factor y curvaturas")
        pe, curvas = {}, {}
        h_max = max(float(sp[0]) for _c, _b, sp in ests)
        lado_ef = EF_SUBCUBO_VOX * h_max
        for c, bw, sp in ests:
            _progreso(progreso, 0.05, f"Ellipsoid Factor: {c}")
            pe[c] = {"EF": ellipsoid_factor(bw, sp, lado_mm=lado_ef)}
            curvas[c] = curvaturas(bw, sp)
        for c, v in perfiles_curvatura(curvas).items():
            pe[c]["curvatura"] = v
        res["forma"] = {"por_estructura": pe,
                        "parametros": {"sigma_vox": SIGMA_VOX,
                                       "margen_vox": MARGEN_CURVATURA,
                                       "nbins": NBINS_CURVATURA,
                                       "EF_lado_subcubo_mm": lado_ef,
                                       "EF_presupuesto": EF_PRESUPUESTO}}
        hechos.append("forma")

    if "superficie" in que and ests and hace_falta("sensibilidad_superficie"):
        _progreso(progreso, 0.3, "sensibilidad a la superficie")
        pe = {c: sensibilidad_superficie(bw, sp, n_mec=n_mec, E_s=E_s,
                                         nu_s=nu_s)
              for c, bw, sp in ests}
        res["sensibilidad_superficie"] = {
            "por_estructura": pe,
            "parametros": next(iter(pe.values()))["parametros"]}
        hechos.append("superficie")

    if "contorno" in que and ests and hace_falta("condiciones_contorno"):
        _progreso(progreso, 0.6, "condiciones de contorno")
        pe = {}
        for c, bw, sp in ests:
            C, n_c, E_c, nu_c = _tensor_de_sesion(res, c)
            if C is None:
                pe[c] = condiciones_contorno(bw, sp, n_mec=n_mec, E_s=E_s,
                                             nu_s=nu_s)
            else:
                # La rejilla y el material del tensor mandan: las tres
                # condiciones tienen que resolver la misma estructura.
                pe[c] = condiciones_contorno(bw, sp, n_mec=n_c, E_s=E_c,
                                             nu_s=nu_c, C_periodica=C)
        res["condiciones_contorno"] = {"por_estructura": pe}
        hechos.append("contorno")

    base = next(((c, bw, sp) for c, bw, sp in ests if c == "voi"),
                ests[0] if ests else None)
    if base is not None and "perdida" in que and hace_falta(
            "simulacion_perdida"):
        from .simulacion import simular_perdida
        _progreso(progreso, 0.8, "simulacion de perdida osea")
        r = simular_perdida(base[1], base[2], guardar_mascaras=False)
        r.pop("mascaras", None)
        r["estructura_codigo"] = base[0]
        res["simulacion_perdida"] = r
        hechos.append("perdida")
    if base is not None and "fallo" in que and hace_falta("fallo_progresivo"):
        from .simulacion import fallo_progresivo
        _progreso(progreso, 0.9, "fallo progresivo")
        r = fallo_progresivo(base[1], base[2], guardar_mascaras=False)
        r.pop("mascaras", None)
        r["estructura_codigo"] = base[0]
        res["fallo_progresivo"] = r
        hechos.append("fallo")

    if hechos:
        res.setdefault("estudios_informe", {})["tiempo_s"] = round(
            time.time() - t0, 1)
    _progreso(progreso, 1.0, "listo")
    return hechos
