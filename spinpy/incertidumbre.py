"""
incertidumbre.py — El error de un ajuste es una variable aleatoria.

EL PROBLEMA
-----------
El generador es estocastico: los mismos parametros con otra semilla dan otra
estructura y otro error. La busqueda escalonada evalua cada punto con UNA
semilla y se queda con el minimo, y el minimo de muchas tiradas ruidosas es
optimista por construccion (maldicion del ganador). Reportar ese numero solo es
reportar la mejor tirada.

Park et al. (2026, arXiv 2607.11209) lo miden para spinodoides: descriptores
de cono identicos dan rigideces dispersas, parametros con la misma media pueden
tener dispersiones distintas, y un optimo determinista incumple las
restricciones en cuanto se tiene en cuenta esa variabilidad. En este proyecto
ya se vio con la desconexion (`Estudio_Percolacion`: de 0.0 a 86.5 % entre
cinco semillas del mismo punto).

LO QUE HACE ESTE MODULO
-----------------------
  * `informe`: K realizaciones nuevas del ganador (semillas distintas de la de
    busqueda), con el error medio +- sd y la media +- sd de cada metrica.
  * `suelo_autoconsistente`: el error entre realizaciones del MISMO candidato,
    tomando una como objetivo y otra como candidato. Es lo minimo que puede
    dar cualquier ajuste aunque conociera la respuesta: por debajo de el, las
    diferencias entre ajustes no significan nada. Es la misma cifra que midio
    `Validacion_Anexo` a mano, ahora en cada ajuste.
  * `seleccion_robusta`: entre los finalistas, el que minimiza
    media + k sd del error sobre K semillas, en vez del mejor de una sola.

K >= 5 (`K_MIN`). Con menos, la sd de una muestra es tan ruidosa (su error
relativo es ~1/sqrt(2(K-1)), 50 % con K = 3) que media + k sd ordena al azar.

Las funciones no saben de familias: reciben `generar(semilla) -> mascara`,
`medir(mascara) -> dict` y `error(m_a, m_b) -> (err, n)`, de modo que el
spinodoide y el dual-lattice las usan igual.
"""

from __future__ import annotations

import numpy as np

K_MIN = 5
SELECCIONES = ("minimo", "robusta")

NOTA = ("El error de la semilla de busqueda es el minimo entre muchas "
        "evaluaciones ruidosas y por eso optimista (maldicion del ganador); "
        "la media sobre K semillas nuevas es la estimacion honesta. El suelo "
        "autoconsistente es el error entre realizaciones del mismo "
        "candidato: por debajo de el, ningun ajuste puede mejorar y las "
        "diferencias entre ajustes no significan nada.")


def validar_replicas(K):
    """0 desactiva; 1..K_MIN-1 es un error, no un valor que se corrija."""
    K = int(K)
    if K == 0:
        return 0
    if K < K_MIN:
        raise ValueError(
            f"replicas={K}: hacen falta al menos {K_MIN}. Con menos, la sd "
            f"de la muestra tiene un error relativo de ~"
            f"{100 / np.sqrt(2 * max(K - 1, 1)):.0f} % y cualquier criterio "
            f"basado en ella decide al azar. 0 desactiva el informe.")
    return K


def semillas(seed, K):
    """Semillas de las replicas: nunca la de busqueda, que es la afortunada."""
    return [int(seed) + 1 + k for k in range(int(K))]


def resumen(valores):
    v = np.asarray([x for x in valores if x is not None], float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {"media": np.nan, "sd": np.nan, "min": np.nan, "max": np.nan,
                "n": 0}
    return {"media": float(v.mean()),
            "sd": float(v.std(ddof=1)) if v.size > 1 else np.nan,
            "min": float(v.min()), "max": float(v.max()), "n": int(v.size)}


def evaluar(generar, medir, error, m_voi, sems, progreso=None, etiqueta=""):
    """Mide y evalua una realizacion por semilla. Devuelve (metricas, errores)."""
    ms, errs = [], []
    for i, s in enumerate(sems):
        if progreso:
            progreso(i, len(sems), etiqueta)
        m = medir(generar(int(s)))
        e, _ = error(m_voi, m)
        ms.append(m)
        errs.append(float(e))
    return ms, errs


def suelo_autoconsistente(ms, error):
    """Error entre pares ORDENADOS de realizaciones del mismo candidato.

    Ordenados porque el error normaliza por el objetivo y no es simetrico.
    """
    vals = []
    for i, a in enumerate(ms):
        for j, b in enumerate(ms):
            if i != j:
                vals.append(float(error(a, b)[0]))
    out = resumen(vals)
    out["n_pares"] = len(vals)
    return out


def resumen_metricas(ms, claves):
    out = {}
    for k in claves:
        v = []
        for m in ms:
            x = m.get(k)
            if x is not None and np.ndim(x) == 0 and np.isfinite(float(x)):
                v.append(float(x))
        if len(v) >= 2:
            r = resumen(v)
            r["cv_pct"] = (100.0 * r["sd"] / abs(r["media"])
                           if r["media"] else np.nan)
            out[k] = r
    return out


def puntuacion(errs, k):
    """media + k sd; infinito si no hay al menos dos errores finitos."""
    v = np.asarray(errs, float)
    v = v[np.isfinite(v)]
    if v.size < 2:
        return np.inf
    return float(v.mean() + float(k) * v.std(ddof=1))


def informe(generar, medir, error, m_voi, seed, K, claves, progreso=None):
    """Incertidumbre del ganador sobre K semillas nuevas."""
    sems = semillas(seed, K)
    ms, errs = evaluar(generar, medir, error, m_voi, sems, progreso,
                       "incertidumbre")
    return {"K": int(K), "semillas": sems, "error": resumen(errs),
            "errores": errs, "metricas": resumen_metricas(ms, claves),
            "suelo_autoconsistente": suelo_autoconsistente(ms, error),
            "nota": NOTA}


def seleccion_robusta(finalistas, generar_de, medir, error, m_voi, seed, K,
                      k=1.0, progreso=None):
    """Indice del finalista con menor media + k sd del error, y la tabla.

    `generar_de(finalista, semilla) -> mascara`. Todos los finalistas usan las
    MISMAS K semillas: la comparacion es pareada y la diferencia entre ellos no
    se mezcla con la de las semillas.
    """
    sems = semillas(seed, K)
    tabla = []
    for i, t in enumerate(finalistas):
        _, errs = evaluar(lambda s, _t=t: generar_de(_t, s), medir, error,
                          m_voi, sems, progreso,
                          f"robusta {i + 1}/{len(finalistas)}")
        tabla.append({"finalista": i, "error": resumen(errs),
                      "errores": errs, "puntuacion": puntuacion(errs, k)})
    idx = int(np.argmin([f["puntuacion"] for f in tabla])) if tabla else 0
    return idx, tabla
