"""
lote.py — Ajuste y generacion de replicas en lote, con salida tabular.

PARA QUE SIRVE Y PARA QUE NO
-----------------------------
El objetivo del proyecto es aumentar el N en estudios de volumenes oseos
escasos. Conviene ser preciso sobre que hace y que no hace este modulo,
porque la diferencia decide si un resultado es publicable:

  GENERAR K REPLICAS DE UN VOI NO DA N = K. Da N = 1 con K replicas tecnicas.
  La varianza total se descompone en la de ENTRE especimenes —la biologica, que
  fijan los animales— y la de DENTRO, que es el ruido del generador. Las
  replicas dividen la segunda por raiz de K y NO TOCAN la primera. Para una
  afirmacion sobre "el hueso sesamoideo equino" los grados de libertad los
  ponen los caballos, no las realizaciones.

Lo que las replicas SI compran:
  * reducir la incertidumbre de la estimacion DE CADA especimen
  * experimentos numericos controlados que el hueso real no permite (variar
    BV/TV con DA fijo, por ejemplo)
  * rellenar huecos del muestreo
  * alimentar un modelo sustituto estructura -> propiedad

Por eso la tabla que produce este modulo lleva SIEMPRE las columnas `animal`,
`sitio` y `replica` separadas: sin esa estructura no se puede descomponer la
varianza despues, y el error solo se detecta cuando ya se gastaron semanas de
CPU. `spinpy.estadistica` opera directamente sobre ellas.

REPRODUCIBILIDAD
----------------
Cada fila registra la semilla y los parametros exactos que la produjeron, de
modo que cualquier espécimen sintetico se puede regenerar bit a bit.
"""

from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import numpy as np

from . import procedencia
from .dual_lattice import generar_dual_lattice
from .elastic import constantes_ingenieria, homogeneizar, remuestrear_bw
from .fit import ajustar_spinodoide, longitud_caracteristica
from .grf import generar_mascara
from .io import leer_voi
from .morphometry import morfometria

METRICAS = ["BVTV", "PoTot", "BSBV", "BSTV", "BS", "BV", "TV",
            "TbTh", "TbSp", "TbN", "DA", "DA2", "FracPort"]


def _etiquetar(ruta):
    """Deduce animal, sitio y tipo del nombre del archivo.

    Los VOIs equinos siguen `VOI_<sitio>_<tipo>.vtk` dentro de una carpeta `H<n>`;
    los porcinos, `VOI_<tag>.mat`. Si no encaja en ningun patron se devuelve el
    nombre completo como animal, que es preferible a inventarse una estructura
    jerarquica que no existe.
    """
    ruta = Path(ruta)
    animal, sitio, tipo = ruta.stem, "", ""
    m = re.match(r"VOI_(\w+?)_(cubico|PCAaligned)$", ruta.stem, re.IGNORECASE)
    if m:
        sitio, tipo = m.group(1), m.group(2)
        animal = ruta.parent.name
    else:
        m2 = re.match(r"VOI_(\w+)$", ruta.stem)
        if m2:
            animal, sitio, tipo = m2.group(1), "unico", "cubico"
    return animal, sitio, tipo


# ---------------------------------------------------------------------------

def _medir(BW, spacing, mecanica=False, res_homog=16, E_s=20e9, nu_s=0.30):
    """Morfometria y, opcionalmente, constantes elasticas homogeneizadas."""
    fila = {}
    m = morfometria(BW, spacing)
    for c in METRICAS:
        v = m.get(c, np.nan)
        fila[c] = float(v) if np.isscalar(v) else np.nan
    fila["MIL_clamped"] = bool(m.get("MIL_clamped", False))
    d = np.asarray(m.get("dir_principal", [0, 0, 1]), float).ravel()
    fila["dir_z"] = float(abs(d[2])) if d.size == 3 else np.nan

    if mecanica:
        bw, sp = remuestrear_bw(BW, spacing, res_homog)
        C, info = homogeneizar(bw, E_s, nu_s, vox_size=sp)
        if info["ok"]:
            ec = constantes_ingenieria(C)
            fila["Ez_Pa"] = ec["Ez"]
            fila["Ez_rel"] = ec["Ez"] / E_s
            fila["Ez_Ex"] = ec["Ez"] / ec["Ex"]
            fila["E_medio_rel"] = ec["E_medio"] / E_s
        else:
            fila.update({k: np.nan for k in
                         ("Ez_Pa", "Ez_rel", "Ez_Ex", "E_medio_rel")})
        fila["homog_ok"] = bool(info["ok"])
        fila["homog_res"] = int(bw.shape[0])
    return fila


# Parametros que describen un ajuste en la tabla, por familia. Se escriben en
# cada fila para que cualquier espécimen sintetico se pueda regenerar.
PARAMETROS_FAMILIA = {
    "spinodoide": ("densidad", "wave_number_pi", "num_waves"),
    "dual-lattice": ("densidad", "celdas", "irregularidad"),
}


def generar_desde_parametros(par, semilla):
    """Mascara de un ajuste a partir de su dict `parametros`, de cualquier familia.

    La traduccion en si vive en `procedencia.kwargs_generador`, que es la
    UNICA del paquete y comprueba que las dos lecturas del numero de onda
    cuadran: un parametro con el valor del deslizador tomado por radianes
    falla ahi en vez de generar una estructura pi veces mas gruesa.
    """
    kw = procedencia.kwargs_generador(par, semilla)
    if par.get("familia") == "dual-lattice":
        BW, _, _ = generar_dual_lattice(**kw)
        return BW
    BW, _, _ = generar_mascara(**kw)
    return BW


def _trabajo_replica(par, sp_cand, semilla, mecanica, res_homog, E_s, nu_s):
    """Una realizacion sintetica. Aislada para poder paralelizarla."""
    BW = generar_desde_parametros(par, semilla)
    fila = _medir(BW, sp_cand, mecanica, res_homog, E_s, nu_s)
    fila["semilla"] = int(semilla)
    return fila


def dispersion_semillas(parametros, spacing, n_semillas=8, semilla_base=None,
                        m_ref=None, extra=False, progreso=None):
    """Repite la MISMA parametrizacion con K semillas y resume la dispersion.

    PARA QUE SIRVE. El generador es estocastico: dos evaluaciones del mismo
    vector de parametros dan realizaciones distintas. Esa dispersion es el
    SUELO DE RUIDO por debajo del cual ninguna diferencia significa nada — ni
    entre dos ajustes, ni entre Python y MATLAB, ni entre un candidato y el
    VOI. Toda la validacion del port se apoya en ella, y hasta ahora habia que
    cambiar la semilla a mano y anotar los numeros.

    Si se pasa `m_ref` (tipicamente la morfometria del VOI), cada metrica trae
    ademas su **z**: cuantas desviaciones tipicas del propio generador separan
    al VOI de la media de las realizaciones. Es el mismo criterio con el que se
    valido el port —|z| < 2 es indistinguible del ruido— y convierte una
    diferencia porcentual, que no dice nada por si sola, en una afirmacion
    comprobable.

    Un CV pequeno NO significa que la metrica sea fiable: significa que el
    generador la reproduce de forma estable. BS tiene un CV del 0.19% a
    densidad 0.55 y aun asi arrastra el sesgo de marching cubes.
    """
    par = dict(parametros)
    if semilla_base is None:
        semilla_base = int(par.get("seed", 20260720))
    par.pop("seed", None)
    # `parametros` es el dict de llamada al generador, no el de un ajuste: la
    # familia viaja como clave aparte porque las dos firmas no se parecen.
    familia = par.pop("familia", "spinodoide")
    generador = (generar_dual_lattice if familia == "dual-lattice"
                 else generar_mascara)

    filas = []
    for k in range(int(n_semillas)):
        if progreso:
            progreso(k, int(n_semillas), f"Realizacion {k+1}")
        BW, _, _ = generador(seed=int(semilla_base) + k, **par)
        m = morfometria(BW, spacing, extra=extra)
        filas.append({kk: float(m[kk]) for kk in METRICAS
                      if kk in m and np.isscalar(m[kk])})
    if progreso:
        progreso(int(n_semillas), int(n_semillas), "Listo")

    resumen = {}
    for met in filas[0]:
        v = np.array([f[met] for f in filas if np.isfinite(f.get(met, np.nan))])
        if v.size < 2:
            continue
        # ddof=1: se estima la dispersion de la POBLACION de realizaciones a
        # partir de una muestra, no se describe la muestra. Con K=8 la
        # diferencia frente a ddof=0 es del 7%, que no es despreciable cuando
        # el resultado se usa como umbral.
        media, sd = float(v.mean()), float(v.std(ddof=1))
        d = {"media": media, "sd": sd, "min": float(v.min()),
             "max": float(v.max()),
             "cv_pct": 100.0 * sd / abs(media) if media else np.nan,
             "n": int(v.size)}
        if m_ref is not None and np.isscalar(m_ref.get(met, None)):
            ref = float(m_ref[met])
            if np.isfinite(ref):
                d["ref"] = ref
                d["dif_rel_pct"] = 100.0 * (media - ref) / abs(ref) if ref else np.nan
                d["z"] = (media - ref) / sd if sd > 0 else np.nan
        resumen[met] = d

    return {"resumen": resumen, "realizaciones": filas,
            "n_semillas": int(n_semillas), "semilla_base": int(semilla_base),
            "familia": familia}


def cerrar_obreros():
    """Cierra los procesos obreros que `loky` deja vivos para reutilizarlos.

    joblib mantiene sus obreros arrancados unos minutos despues de que
    `Parallel` devuelva, para no pagar el arranque la proxima vez. El precio
    es que HEREDAN LA SALIDA ESTANDAR del proceso padre: si alguien canaliza
    la salida del programa, la tuberia no ve nunca el final aunque el padre
    haya terminado. Compilando se perdieron 27 minutos por esto, con un obrero
    girando al 100 % de una CPU despues de que el guion hubiera acabado.

    Se llama al terminar el lote. Si el backend no es loky no hay nada que
    cerrar y no pasa nada.
    """
    try:
        from joblib.externals.loky import get_reusable_executor
        get_reusable_executor().shutdown(wait=True)
    except Exception:
        pass


def opciones_paralelo(n_jobs=-1):
    """Como paralelizar, segun se corra desde el codigo o desde el ejecutable.

    DESDE EL CODIGO: procesos. Es lo correcto — el trabajo es CPU pura y el
    GIL lo serializaria.

    EMPAQUETADO CON PyInstaller: hilos, y no es una preferencia sino la unica
    opcion que funciona. joblib usa `loky`, que arranca a sus obreros
    reejecutando `sys.executable` con argumentos del estilo `-c "from loky
    ..."`. En un ejecutable congelado `sys.executable` es la propia
    aplicacion, que no acepta `-c`: los obreros mueren nada mas nacer y
    joblib levanta `TerminatedWorkerError`. Medido sobre el .exe construido
    aqui: el obrero muere a los 24 s y el proceso padre se queda colgado.

    Con hilos se pierde parte del paralelismo, pero no todo: generar la
    mascara, la morfometria y la homogeneizacion pasan casi todo su tiempo
    dentro de numpy, scipy y pyamg, que sueltan el GIL. El lote tarda mas en
    el ejecutable que desde el codigo, y conviene saberlo antes de lanzarlo.
    """
    if getattr(sys, "frozen", False):
        return {"n_jobs": n_jobs, "prefer": "threads"}
    return {"n_jobs": n_jobs, "prefer": "processes"}


def correr_lote(vois, n_replicas=10, modo="completo", mecanica=False,
                res_homog=16, E_s=20e9, nu_s=0.30, semilla_base=770000,
                n_jobs=-1, progreso=None, familia="spinodoide",
                opciones_ajuste=None):
    """Ajusta cada VOI y genera `n_replicas` sinteticas de cada ajuste.

    familia : "spinodoide" (por omision, el comportamiento de siempre) o
        "dual-lattice". Las dos se ajustan con la misma busqueda escalonada y
        el mismo modo; la tabla lleva una columna `familia` para que dos lotes
        de familias distintas se puedan concatenar sin mezclarse.

    opciones_ajuste : dict que se pasa tal cual al ajuste (`pesos`,
        `distancia`, `replicas`, `seleccion`, `k_robusto`, `alineacion`...).
        Lo que se uso queda en la tabla de ajustes, columna `objetivo`.

    PROCEDENCIA. Cada fila de las dos tablas lleva las columnas de
    `procedencia.COLUMNAS` (version de spinpy y del formato, familia, esquema,
    numero de onda con sus dos lecturas). En las filas sinteticas `semilla` es
    la de ESA replica y la del ajuste va en `semilla_ajuste`.

    Devuelve un DataFrame ordenado, una fila por espécimen:

        origen    'real' para el VOI medido, 'sintetico' para las replicas
        animal    H1..H4, C1..V5 — la unidad de inferencia biologica
        sitio     proximal / medio / distal
        replica   0 para el real, 1..K para las sinteticas
        semilla   la que genero esa realizacion (NaN en los reales)
        ...       parametros ajustados y metricas

    La fila real de cada VOI se incluye con `replica = 0` para que la
    comparacion real-vs-sintetico salga de la misma tabla y no de dos ficheros
    que puedan desincronizarse.
    """
    import pandas as pd
    from joblib import Parallel, delayed

    filas = []
    ajustes = []
    vois = [Path(v) for v in vois]
    t0 = time.time()

    for iv, ruta in enumerate(vois, start=1):
        animal, sitio, tipo = _etiquetar(ruta)
        if progreso:
            progreso(iv - 1, len(vois), f"Ajustando {ruta.name}")

        VOI, spacing = leer_voi(ruta)
        m_voi = morfometria(VOI, spacing)

        # --- espécimen real (replica 0) ---
        f = _medir(VOI, spacing, mecanica, res_homog, E_s, nu_s)
        f.update(procedencia.columnas(procedencia.bloque(familia)))
        f.update({"origen": "real", "animal": animal, "sitio": sitio,
                  "tipo": tipo, "replica": 0, "semilla": np.nan,
                  "archivo": ruta.name, "familia": familia,
                  # El VOI no se genera: ni esquema ni numero de onda.
                  "esquema": None, "wave_number_pi": None,
                  "wave_number_rad": None})
        filas.append(f)

        # --- ajuste ---
        kw = {"modo": modo, "m_voi": m_voi}
        if familia == "dual-lattice":
            from .fit_dual import ajustar_dual_lattice
            kw.update(opciones_ajuste or {})
            r = ajustar_dual_lattice(VOI, spacing, **kw)
            par = r["parametros"]
            propios = {k: par[k] for k in ("densidad", "celdas",
                                           "irregularidad", "resolucion")}
            propios["estiramiento"] = par["estiramiento"]
        else:
            kw["precision"] = "f32"
            kw.update(opciones_ajuste or {})
            r = ajustar_spinodoide(VOI, spacing, **kw)
            par = r["parametros"]
            propios = {k: par[k] for k in ("densidad", "wave_number_pi",
                                           "num_waves", "resolucion",
                                           "esquema")}
            propios["thetas"] = par["thetas"]
        cols_proc = procedencia.columnas(r["procedencia"])
        inc = r.get("incertidumbre") or {}
        alin = r.get("alineacion") or {}
        ajustes.append({**cols_proc, "familia": familia, "animal": animal,
                        "sitio": sitio, "tipo": tipo, "archivo": ruta.name,
                        "error": r["error"],
                        "error_media_K": (inc.get("error") or {}).get("media"),
                        "error_sd_K": (inc.get("error") or {}).get("sd"),
                        "K": inc.get("K", 0),
                        "suelo_media": (inc.get("suelo_autoconsistente")
                                        or {}).get("media"),
                        "alineada": alin.get("alineada"),
                        "angulo_eje_deg": alin.get("angulo_despues_deg"),
                        "objetivo": repr((r.get("objetivo") or {}).get("pesos")),
                        "n_evaluaciones": r["n_evaluaciones"],
                        "tiempo_s": r["tiempo_s"], **propios})

        lc = longitud_caracteristica(VOI.shape, spacing)
        sp_cand = np.full(3, lc / par["resolucion"])

        # --- replicas sinteticas, en paralelo ---
        if progreso:
            progreso(iv - 1, len(vois),
                     f"Generando {n_replicas} replicas de {ruta.name}")
        semillas = semilla_base + 1000 * iv + np.arange(n_replicas)
        res = Parallel(**opciones_paralelo(n_jobs))(
            delayed(_trabajo_replica)(par, sp_cand, s, mecanica, res_homog,
                                      E_s, nu_s) for s in semillas)
        for k, f in enumerate(res, start=1):
            f.update({c: v for c, v in cols_proc.items() if c != "semilla"})
            f.update({"origen": "sintetico", "animal": animal, "sitio": sitio,
                      "tipo": tipo, "replica": k, "archivo": ruta.name,
                      "error_ajuste": r["error"], "familia": familia,
                      "semilla_ajuste": cols_proc["semilla"]})
            f.update({k2: par[k2] for k2 in PARAMETROS_FAMILIA[familia]})
            filas.append(f)

    if progreso:
        progreso(len(vois), len(vois), "Listo")

    orden = ["origen", "animal", "sitio", "tipo", "replica", "semilla",
             "archivo"] + METRICAS
    df = pd.DataFrame(filas)
    cols = [c for c in orden if c in df.columns]
    df = df[cols + [c for c in df.columns if c not in cols]]
    cerrar_obreros()
    return df, pd.DataFrame(ajustes), {"tiempo_s": time.time() - t0,
                                       "n_vois": len(vois),
                                       "n_replicas": n_replicas}
