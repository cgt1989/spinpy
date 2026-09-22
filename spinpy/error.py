"""
error.py — Funcion de error morfometrico y criterios de fabrica.

Port de:
    computeMorphometricError   (AppFinal_V2.m:3725-3785)
    localFabricCos             (AppFinal_V2.m:6051)
    localFabricAlignWeight     (AppFinal_V2.m:6076)
    localTieneEje              (AppFinal_V2.m:6035)

INVARIANTE C6 — LA RAZON DE SER DE ESTA FUNCION
------------------------------------------------
La version anterior de la app sumaba los terminos disponibles y saltaba EN
SILENCIO los que salieran NaN. Consecuencia: un candidato en el que alguna
metrica no se podia calcular acumulaba MENOS error que uno en el que si se
calculaba, y por tanto parecia mejor. El optimizador quedaba sesgado
justamente hacia las geometrias degeneradas, que son las que mas metricas
hacen fallar.

La correccion es dividir por la suma de pesos EFECTIVAMENTE usados, de modo
que el error es un promedio ponderado y no depende de cuantos terminos
sobrevivan. Si se rompe esto, el sesgo vuelve sin dar ningun sintoma.

TERMINO MECANICO (J1)
---------------------
La app puede anadir Ez_rel y Ez_Ex, que requieren homogeneizacion elastica.
Sus pesos se multiplican por `peso_mecanico`; con 0 la funcion se comporta
como la version sin termino mecanico, bit a bit — que es exactamente lo que
hace MATLAB cuando el ajuste mecanico esta desactivado.

TERMINOS OPCIONALES (correccion O1) — NINGUNO ENTRA POR OMISION
----------------------------------------------------------------
`Estudio_Discriminadores` midio que el objetivo iguala lo que no distingue y
deja fuera lo que si:

    Tb.Sp               0.6 sigma  — es una FORMULA, Tb.Th (1/BV.TV - 1), de
                                     dos terminos que el ajuste ya iguala
    H desv (curvatura)  64 sigma   Po.Dm medio      39 sigma
    Tb.Th CV            23 sigma   frac. placas EF  36 sigma

y la literatura reciente compara distribuciones (Xiao et al. 2025, distancia
de Hellinger entre distribuciones de placas y barras) y el tensor de rigidez
entero (Wang et al. 2025; Otto et al. 2025, perdida logaritmica). El argumento
`pesos` da acceso a todo eso:

  * sustituye el peso de cualquier termino base —`{"TbSp": 0}` lo quita—;
  * activa terminos de `TERMINOS_OPCIONALES` con el peso que se les de.

`pesos=None` (por omision) recorre EXACTAMENTE el mismo camino de siempre, en
el mismo orden de sumas, y el error es bit a bit el validado contra MATLAB
(`validar_error.py`, 1173 pares; el bloque 17 de `tests/` lo repite). Un
termino opcional ausente en alguno de los dos dicts se salta y renormaliza
(C6), igual que los base.

Como se mide cada termino, para que sumen en la misma escala que las
diferencias relativas al cuadrado de los terminos base:

  escalar       ((b - a) / |a|)^2
  distribucion  Wasserstein: (W1(a, b) / |media de a|)^2. Si b es a escalada
                por un factor s, vale (s - 1)^2, lo mismo que un escalar.
                Hellinger: H^2, en [0, 1].
  tensor        d_LE(a, b)^2 / 6, con d_LE la distancia log-euclidea entre
                los tensores en sus ejes materiales. Si b = s a, vale
                (ln s)^2, que para s cerca de 1 es (s - 1)^2.

Los terminos mecanicos opcionales NO se multiplican por `peso_mecanico`: su
peso es el que se les da en `pesos`. `peso_mecanico` sigue gobernando solo
Ez_rel y Ez_Ex, como en MATLAB.
"""

from __future__ import annotations

import numpy as np

# Pesos de AppFinal_V2.m:3757-3758. Los dos ultimos se multiplican por el peso
# del ajuste mecanico.
CAMPOS = ["BVTV", "DA", "BSBV", "TbTh", "TbSp", "TbN", "PoTot", "Ez_rel", "Ez_Ex"]
PESOS_BASE = [3.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.5, 2.0, 1.0]

PENALIZACION = 1e6      # nada medible

# Terminos que se pueden activar con `pesos`. `medida` es la compuerta de
# `morphometry.metricas_forma` que los produce; `coste` decide DONDE los paga
# el ajuste: "forma" en cada evaluacion de la busqueda, "caro" y "mecanico"
# solo entre finalistas (etapa D). El orden de este dict es el orden de suma.
TERMINOS_OPCIONALES = {
    "PoDm":          {"tipo": "escalar", "medida": "poro", "coste": "forma"},
    "TbTh_CV":       {"tipo": "escalar", "medida": "tbth", "coste": "forma"},
    "H_desv":        {"tipo": "escalar", "medida": "curvatura",
                      "coste": "forma"},
    "EF_frac_placa": {"tipo": "escalar", "medida": "ef", "coste": "caro"},
    "dist_TbTh":     {"tipo": "distribucion", "medida": "tbth",
                      "coste": "forma", "cuantiles": "TbTh"},
    "dist_PoDm":     {"tipo": "distribucion", "medida": "poro",
                      "coste": "forma", "cuantiles": "PoDm"},
    "dist_H":        {"tipo": "distribucion", "medida": "curvatura",
                      "coste": "forma", "cuantiles": "H"},
    "C_logE":        {"tipo": "tensor", "medida": None, "coste": "mecanico",
                      "clave": "C_principal_rel"},
    "E1_rel":        {"tipo": "escalar", "medida": None, "coste": "mecanico"},
    "E2_rel":        {"tipo": "escalar", "medida": None, "coste": "mecanico"},
    "E3_rel":        {"tipo": "escalar", "medida": None, "coste": "mecanico"},
    "G23_rel":       {"tipo": "escalar", "medida": None, "coste": "mecanico"},
    "G13_rel":       {"tipo": "escalar", "medida": None, "coste": "mecanico"},
    "G12_rel":       {"tipo": "escalar", "medida": None, "coste": "mecanico"},
}

# Sugerencia, NO valor por omision: fuera el termino que no separa nada y
# dentro los tres baratos que mas separan. El EF queda fuera por coste.
PESOS_DISCRIMINANTES = {"TbSp": 0.0, "PoDm": 1.0, "TbTh_CV": 1.0,
                        "H_desv": 1.0}

DISTANCIAS = ("wasserstein", "hellinger")


def _escalar(v):
    """Devuelve el valor si es un escalar finito utilizable; None si no."""
    if v is None:
        return None
    a = np.asarray(v)
    if a.size != 1:
        return None
    x = float(a.reshape(()))
    return x if np.isfinite(x) else None


def validar_pesos(pesos):
    """Comprueba `pesos` y lo devuelve como dict de floats (o None).

    Un nombre mal escrito se rechaza: con un dict permisivo, `{"PoDM": 1}`
    no activaria nada y el ajuste parecia haber usado el termino.
    """
    if pesos is None:
        return None
    out = {}
    for k, v in dict(pesos).items():
        if k == "DA2":
            raise ValueError("el peso de DA2 se da con `peso_da2`, no en "
                             "`pesos`")
        if k not in CAMPOS and k not in TERMINOS_OPCIONALES:
            raise ValueError(
                f"termino desconocido en pesos: {k!r}. Base: {CAMPOS}; "
                f"opcionales: {list(TERMINOS_OPCIONALES)}")
        w = float(v)
        if not np.isfinite(w) or w < 0:
            raise ValueError(f"peso de {k!r} no valido: {v!r}")
        out[k] = w
    return out


def medidas_necesarias(pesos):
    """Que hay que medir para poder evaluar estos pesos.

    Devuelve {"forma": set, "caro": set, "mecanico": bool}: las compuertas de
    `metricas_forma` que hacen falta en cada evaluacion, las que solo se pagan
    entre finalistas, y si algun termino necesita homogeneizar.
    """
    pesos = validar_pesos(pesos) or {}
    forma, caro, mec = set(), set(), False
    for k, w in pesos.items():
        t = TERMINOS_OPCIONALES.get(k)
        if t is None or w <= 0:
            continue
        if t["coste"] == "forma":
            forma.add(t["medida"])
        elif t["coste"] == "caro":
            caro.add(t["medida"])
        elif t["coste"] == "mecanico":
            mec = True
    return {"forma": forma, "caro": caro, "mecanico": mec}


# ---------------------------------------------------------------------------
# Distancias entre distribuciones dadas por sus cuantiles
# ---------------------------------------------------------------------------

def _niveles(n):
    return (np.arange(int(n)) + 0.5) / int(n)


def _a_niveles(q, n):
    q = np.asarray(q, float).ravel()
    if q.size == n:
        return q
    return np.interp(_niveles(n), _niveles(q.size), q)


def wasserstein_cuantiles(qa, qb):
    """Distancia de Wasserstein-1 entre dos distribuciones 1D.

    W1 = integral_0^1 |F_a^-1(p) - F_b^-1(p)| dp, por la regla del punto medio
    sobre los niveles de Hazen. Si las tablas tienen distinto numero de
    niveles, se lleva la corta a los de la larga.
    """
    qa = np.asarray(qa, float).ravel()
    qb = np.asarray(qb, float).ravel()
    n = max(qa.size, qb.size)
    return float(np.mean(np.abs(_a_niveles(qa, n) - _a_niveles(qb, n))))


def _cdf(x, q):
    q = np.maximum.accumulate(np.asarray(q, float).ravel())
    esc = max(float(np.abs(q).max()), 1e-300)
    # rampa infinitesimal: np.interp necesita abscisas crecientes y una
    # distribucion con valores repetidos (un pico) no las da
    q = q + np.arange(q.size) * 1e-12 * esc
    return np.interp(x, q, _niveles(q.size), left=0.0, right=1.0)


def hellinger_cuantiles(qa, qb, nbins=32):
    """Distancia de Hellinger entre dos distribuciones 1D, en [0, 1].

    Se reconstruyen las funciones de distribucion desde los cuantiles, se
    reparten en `nbins` clases comunes sobre el rango conjunto y
    H = sqrt(1 - sum sqrt(p_a p_b)). 0 si son iguales, 1 si no se solapan.
    """
    qa = np.asarray(qa, float).ravel()
    qb = np.asarray(qb, float).ravel()
    lo = float(min(qa.min(), qb.min()))
    hi = float(max(qa.max(), qb.max()))
    if not hi > lo:
        return 0.0
    bordes = np.linspace(lo, hi, int(nbins) + 1)
    bordes[0] = np.nextafter(lo, -np.inf)
    pa = np.diff(_cdf(bordes, qa))
    pb = np.diff(_cdf(bordes, qb))
    sa, sb = pa.sum(), pb.sum()
    if sa <= 0 or sb <= 0:
        return np.nan
    bc = float(np.sum(np.sqrt((pa / sa) * (pb / sb))))
    return float(np.sqrt(max(0.0, 1.0 - bc)))


def _termino_opcional(nombre, m_a, m_b, distancia):
    """Valor de un termino opcional (sin peso), o None si no es evaluable."""
    t = TERMINOS_OPCIONALES[nombre]
    if t["tipo"] == "escalar":
        v1, v2 = _escalar(m_a.get(nombre)), _escalar(m_b.get(nombre))
        if v1 is None or v2 is None or abs(v1) <= 0:
            return None
        rel = (v2 - v1) / abs(v1)
        return rel * rel

    if t["tipo"] == "distribucion":
        qa = (m_a.get("cuantiles") or {}).get(t["cuantiles"])
        qb = (m_b.get("cuantiles") or {}).get(t["cuantiles"])
        if qa is None or qb is None or len(qa) < 2 or len(qb) < 2:
            return None
        qa, qb = np.asarray(qa, float), np.asarray(qb, float)
        if not (np.all(np.isfinite(qa)) and np.all(np.isfinite(qb))):
            return None
        if distancia == "hellinger":
            h = hellinger_cuantiles(qa, qb)
            return None if not np.isfinite(h) else h * h
        mu = abs(float(np.mean(qa)))
        if mu <= 0:
            return None
        w = wasserstein_cuantiles(qa, qb) / mu
        return w * w

    # tensor
    from .elastic import distancia_log_euclidea
    Ca, Cb = m_a.get(t["clave"]), m_b.get(t["clave"])
    if Ca is None or Cb is None:
        return None
    Ca, Cb = np.asarray(Ca, float), np.asarray(Cb, float)
    if Ca.shape != (6, 6) or Cb.shape != (6, 6) or \
            not (np.all(np.isfinite(Ca)) and np.all(np.isfinite(Cb))):
        return None
    try:
        d = distancia_log_euclidea(Ca, Cb)
    except ValueError:
        return None
    return d * d / 6.0


def error_morfometrico(m_voi, m_spin, peso_mecanico=0.0, peso_da2=0.0,
                       pesos=None, distancia="wasserstein"):
    """Error morfometrico ponderado entre VOI y candidato.

    Devuelve (err, n_usados). `err` es un promedio ponderado de diferencias
    relativas al cuadrado; vale PENALIZACION si no hubo ningun termino
    utilizable.

    peso_mecanico : equivalente de localMechFitWeight(). Escala los pesos de
        Ez_rel y Ez_Ex, que salen de la homogeneizacion periodica
        (`spinpy.fit.metricas_mecanicas`). Si se pone > 0 pero esos campos no
        estan en los dicts, simplemente no se usan y el error se renormaliza
        por el peso realmente aplicado (invariante C6): un termino ausente no
        puede contar como un termino con error cero.

        En el ajuste no se activa dentro de la busqueda sino en una etapa de
        desempate entre finalistas — homogeneizar las ~53 evaluaciones seria
        entre tres y cuatro ordenes de magnitud mas caro. Ver
        `ajustar_spinodoide`.

    peso_da2 : peso del termino DA2 (anisotropia secundaria, autovalor mayor
        entre el intermedio del MIL). NO esta en MATLAB y vale 0 por omision,
        de modo que el error del spinodoide sigue siendo bit a bit el validado.
        Existe para el dual-lattice: el DA escalar no distingue una fabrica
        axisimetrica (DA2 ~ DA, lo que da alargar un solo eje) de una
        triaxial como la del hueso equino (DA 1.45-1.51, DA2 1.12-1.18), y
        sin este termino un preset triaxial con el mismo DA solo empata con
        el axial (Estudio_Familias, 4.2 y recomendacion 5).

    pesos : dict opcional (correccion O1, ver la cabecera del modulo). Sustituye
        pesos base por nombre y activa terminos de `TERMINOS_OPCIONALES`. None
        deja el camino validado contra MATLAB intacto.

    distancia : "wasserstein" (por omision) o "hellinger", para los terminos
        `dist_*`.
    """
    if distancia not in DISTANCIAS:
        raise ValueError(f"distancia debe ser una de {DISTANCIAS}")
    pesos_base = list(PESOS_BASE)
    extra = {}
    if pesos is not None:
        pesos = validar_pesos(pesos)
        for i, c in enumerate(CAMPOS):
            if c in pesos:
                pesos_base[i] = pesos[c]
        extra = {k: w for k, w in pesos.items() if k in TERMINOS_OPCIONALES}

    pesos_l = pesos_base
    pesos_l[7] *= peso_mecanico
    pesos_l[8] *= peso_mecanico
    campos = list(CAMPOS)
    if peso_da2 > 0:
        campos.append("DA2")
        pesos_l.append(float(peso_da2))

    acc = 0.0
    peso_usado = 0.0
    n_usados = 0

    for campo, w in zip(campos, pesos_l):
        if w <= 0:
            continue
        v1 = _escalar(m_voi.get(campo)) if hasattr(m_voi, "get") else None
        v2 = _escalar(m_spin.get(campo)) if hasattr(m_spin, "get") else None
        if v1 is None or v2 is None or abs(v1) <= 0:
            continue
        rel = (v2 - v1) / abs(v1)
        acc += w * rel * rel
        peso_usado += w        # C6: se acumula el peso REALMENTE utilizado
        n_usados += 1

    if extra and hasattr(m_voi, "get") and hasattr(m_spin, "get"):
        for nombre in TERMINOS_OPCIONALES:
            w = extra.get(nombre, 0.0)
            if w <= 0:
                continue
            t = _termino_opcional(nombre, m_voi, m_spin, distancia)
            if t is None:
                continue
            acc += w * t
            peso_usado += w
            n_usados += 1

    if peso_usado <= 0:
        return PENALIZACION, 0
    return acc / peso_usado, n_usados


# ---------------------------------------------------------------------------
# Criterios de fabrica (orientacion)
# ---------------------------------------------------------------------------

def tiene_eje(m):
    """Port de localTieneEje: ¿el struct trae una direccion principal usable?"""
    if not isinstance(m, dict):
        return False
    d = m.get("dir_principal")
    if d is None:
        return False
    d = np.asarray(d, dtype=float).ravel()
    if d.size != 3 or not np.all(np.isfinite(d)) or np.linalg.norm(d) <= 0:
        return False
    return bool(m.get("MIL_valid", True))


def fabric_cos(m_a, m_b):
    """Coseno del angulo entre direcciones principales, en valor absoluto.

    Se compara como EJE y no como vector: el signo de un autovector es
    arbitrario, de ahi el valor absoluto.
    """
    if not tiene_eje(m_a) or not tiene_eje(m_b):
        return np.nan
    a = np.asarray(m_a["dir_principal"], float).ravel()
    b = np.asarray(m_b["dir_principal"], float).ravel()
    c = abs(float((a / np.linalg.norm(a)) @ (b / np.linalg.norm(b))))
    return float(min(1.0, max(0.0, c)))


def peso_alineacion(m_voi, da_min=1.15, da_pleno=1.37, w_max=2.0):
    """Port de localFabricAlignWeight.

    Cuanto vale desempatar por orientacion en ESTE VOI. Devuelve 0 —es decir,
    no desempatar— cuando la direccion principal no significa nada:

      * DA por debajo de 1.15: por debajo del suelo de ruido del MIL (~1.07)
        la direccion principal es ruido, no estructura.
      * MIL acotado: en estructuras laminares el elipsoide degenera y el DA es
        solo una cota inferior.

    Entre 1.15 y 1.37 el peso sube linealmente hasta el maximo.
    """
    if not isinstance(m_voi, dict):
        return 0.0
    da = _escalar(m_voi.get("DA"))
    if da is None:
        return 0.0
    if m_voi.get("MIL_clamped", False):
        return 0.0
    if not tiene_eje(m_voi):
        return 0.0
    if da <= da_min:
        return 0.0
    return float(w_max * min(1.0, (da - da_min) / (da_pleno - da_min)))
