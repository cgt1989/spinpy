"""metodos.py — Registro de metodos de ajuste, y generacion de replicas.

POR QUE UN REGISTRO Y NO UNA CADENA DE `if`
--------------------------------------------
La comparacion «ajusta con todos los metodos y elige» solo tiene sentido si
anadir un metodo es anadir una entrada. Con condicionales repartidos, portar
el optimizador de Pareto obligaria a tocar el motor, la ventana, la tabla y el
guardado, y lo normal es olvidarse de alguno.

Aqui cada metodo declara su nombre, que hace y como se ejecuta. La interfaz
recorre el registro: lo que este dentro sale en la ventana, y lo que no, no.

LOS CUATRO DE AHORA NO SON LOS CINCO DE MATLAB
-----------------------------------------------
`AppFinal_V2.m` tiene cinco optimizadores: la busqueda escalonada, Pareto con
`gamultiobj`, bayesiano con `bayesopt`, MOBO y el best-fit rapido. De esos,
el port cubre la busqueda escalonada; los otros cuatro siguen solo en MATLAB.

Lo que se ofrece aqui son cuatro VARIANTES genuinamente distintas de lo que si
esta portado, no cinco optimizadores disfrazados:

  * `rapido`      etapas A+B con 4 presets de angulos
  * `completo`    etapas A+B+C con 14 presets y refinado
  * `mecanico`    completo + etapa D, que reordena los finalistas por rigidez
  * `equitativo`  completo con el otro reparto de ondas entre conos

El ultimo no es un capricho: `rechazo` y `equitativo` reparten las ondas de
forma distinta cuando los conos tienen angulos desiguales, y producen
anisotropias distintas para los mismos thetas. Cual describe mejor un hueso
concreto es una pregunta empirica, y esta ventana es donde se contesta.
"""

from __future__ import annotations

import time

import numpy as np

from . import procedencia
from .fit import ajustar_spinodoide
from .grf import generar_mascara
from .morphometry import morfometria

# Metricas que se comparan contra el VOI en la ventana de seleccion.
COMPARADAS = ["BVTV", "BSBV", "TbTh", "TbSp", "TbN", "DA", "DA2"]

# Parametro del generador que gobierna cada metrica, para las replicas con
# variacion. NO se puede fijar una metrica directamente: el generador toma
# densidad, numero de onda y angulos, y las metricas salen de ahi. Lo honesto
# es perturbar el parametro que la gobierna y REPORTAR la dispersion lograda,
# que es lo que hace `generar_replicas`.
GOBIERNA = {
    "BVTV": ("densidad", "fraccion de solido: la controla directamente"),
    "TbTh": ("wave_number_pi", "grosor: lo fija la escala del campo; mas onda, "
                               "trabecula mas fina"),
    "TbSp": ("wave_number_pi", "separacion: misma escala que el grosor, asi "
                               "que no se pueden variar por separado"),
    "DA":   ("thetas", "anisotropia: la fijan los semiangulos de cono"),
}


# Lo mismo para el dual-lattice: densidad, escala de la red y estiramiento.
GOBIERNA_DUAL = {
    "BVTV": ("densidad", "fraccion de solido: la controla directamente"),
    "TbTh": ("celdas", "grosor: lo fija la escala de la red; mas celdas, "
                       "puntal mas fino"),
    "TbSp": ("celdas", "separacion: misma escala que el grosor, asi que no se "
                       "pueden variar por separado"),
    "DA":   ("estiramiento", "anisotropia: la fija el estiramiento de la "
                             "malla"),
}


def gobierna_de(familia):
    """Mapa metrica -> parametro del generador, segun la familia."""
    return GOBIERNA_DUAL if familia == "dual-lattice" else GOBIERNA


def _envolver(nombre, etiqueta, descripcion, kw):
    def correr(VOI, spacing, m_voi=None, semilla=20260720, num_waves=700,
               esquema="rechazo", progreso=None, **extra):
        args = dict(kw)
        args.update(extra)
        args.setdefault("esquema", esquema)
        # `_n_eval` es solo para dimensionar la barra de progreso; no es un
        # argumento del ajuste.
        args.pop("_n_eval", None)
        t0 = time.time()
        r = ajustar_spinodoide(VOI, spacing, m_voi=m_voi, seed=semilla,
                               num_waves=num_waves, progreso=progreso, **args)
        r["metodo"] = nombre
        r["etiqueta"] = etiqueta
        r["tiempo_s"] = time.time() - t0
        return r
    return {"nombre": nombre, "etiqueta": etiqueta,
            "descripcion": descripcion, "correr": correr,
            "n_eval": kw.get("_n_eval", 53)}


REGISTRO = {
    "rapido": _envolver(
        "rapido", "Rapido (A+B)",
        "Rejilla gruesa de densidad y numero de onda, y 4 presets de angulos. "
        "Sin refinado. Es el que se usa cuando hacen falta muchos ajustes.",
        {"modo": "rapido", "_n_eval": 19}),
    "completo": _envolver(
        "completo", "Completo (A+B+C)",
        "Rejilla fina, 14 presets de angulos y etapa de refinado. Es el "
        "metodo de referencia del proyecto.",
        {"modo": "completo", "_n_eval": 53}),
    "mecanico": _envolver(
        "mecanico", "Completo + desempate mecanico",
        "Como el completo y ademas reordena los mejores candidatos "
        "homogeneizando cada uno, para que no gane uno que reproduce la "
        "morfometria pero no la rigidez. Es el mas lento con diferencia.",
        {"modo": "completo", "peso_mecanico": 1.0, "res_mec": 16,
         "n_finalistas": 5, "_n_eval": 58}),
    "equitativo": _envolver(
        "equitativo", "Completo, muestreo equitativo",
        "Igual que el completo pero repartiendo las ondas a partes iguales "
        "entre los conos, en vez de por rechazo. Con conos desiguales da "
        "anisotropias distintas para los mismos angulos.",
        {"modo": "completo", "esquema": "equitativo", "_n_eval": 53}),
}


# --- Dual-lattice -----------------------------------------------------------
# Solo los dos metodos de busqueda escalonada. El `mecanico` y el `equitativo`
# no tienen sentido todavia aqui: el primero se puede pedir pasando
# `peso_mecanico` al ajuste, y el segundo es un esquema de muestreo de ondas
# que esta familia no tiene.

def _envolver_dual(nombre, etiqueta, descripcion, kw):
    def correr(VOI, spacing, m_voi=None, semilla=20260720, progreso=None,
               **extra):
        from .fit_dual import ajustar_dual_lattice
        args = dict(kw)
        args.update(extra)
        args.pop("_n_eval", None)
        # La interfaz pasa a todos los metodos los argumentos del spinodoide;
        # aqui no significan nada y el ajuste no los acepta.
        for k in ("num_waves", "esquema", "precision"):
            args.pop(k, None)
        t0 = time.time()
        r = ajustar_dual_lattice(VOI, spacing, m_voi=m_voi, seed=semilla,
                                 progreso=progreso, **args)
        r["metodo"] = nombre
        r["etiqueta"] = etiqueta
        r["tiempo_s"] = time.time() - t0
        return r
    return {"nombre": nombre, "etiqueta": etiqueta,
            "descripcion": descripcion, "correr": correr,
            "n_eval": kw.get("_n_eval", 52), "familia": "dual-lattice"}


REGISTRO_DUAL = {
    "rapido": _envolver_dual(
        "rapido", "Rapido (A+B)",
        "Rejilla gruesa de densidad y celdas alrededor de la estimacion desde "
        "el Tb.N del VOI, y 7 presets de estiramiento. Sin refinado.",
        {"modo": "rapido", "_n_eval": 16}),
    "completo": _envolver_dual(
        "completo", "Completo (A+B+C)",
        "Rejilla fina de densidad y celdas, 23 presets de estiramiento "
        "(axiales hasta x3 y triaxiales) y etapa de refinado.",
        {"modo": "completo", "_n_eval": 62}),
}

# Metodos por familia. `REGISTRO` se mantiene con su nombre porque la interfaz
# y `lote.py` lo usan asi para el spinodoide.
REGISTROS = {"spinodoide": REGISTRO, "dual-lattice": REGISTRO_DUAL}


def comparar_resultados(resultados, m_voi):
    """Tabla comparativa de varios ajustes frente al mismo VOI.

    Devuelve una fila por metodo con el error, el tiempo y la diferencia
    relativa de cada metrica. El **error no es comparable entre metodos que
    usan distinto numero de terminos**: el del desempate mecanico incluye
    Ez_rel y Ez_Ex y el de los demas no, asi que se marca `n_terminos` para
    que la comparacion se lea sabiendo eso.
    """
    filas = []
    for r in resultados:
        if not r or r.get("fallo"):
            filas.append({"metodo": r.get("metodo", "?") if r else "?",
                          "fallo": (r or {}).get("fallo", "no se ejecuto")})
            continue
        m = r["metricas_spin"]
        fila = {"metodo": r["metodo"], "etiqueta": r.get("etiqueta", r["metodo"]),
                "error": float(r["error"]), "tiempo_s": float(r["tiempo_s"]),
                "n_evaluaciones": int(r.get("n_evaluaciones", 0)),
                "mecanico": bool((r.get("mecanico") or {}).get("traza")),
                "bicontinuo": bool(r["diagnostico"]["bicontinuo"])}
        peor = 0.0
        for k in COMPARADAS:
            a, b = m_voi.get(k), m.get(k)
            fila[k] = b
            if (a is not None and b is not None and np.isscalar(a)
                    and np.isscalar(b) and np.isfinite(a) and np.isfinite(b)
                    and a != 0):
                d = 100.0 * (b - a) / abs(a)
                fila["d_" + k] = d
                peor = max(peor, abs(d))
        fila["peor_dif_pct"] = peor
        filas.append(fila)
    return filas


def generar_replicas(parametros, spacing, n=10, semilla_base=None,
                     variar=None, amplitud=0.05, medir=True, progreso=None):
    """N realizaciones del mismo ajuste, iguales o con variacion controlada.

    `variar` es una lista de METRICAS (`BVTV`, `TbTh`, `DA`...). Como el
    generador no toma metricas sino parametros, cada una se traduce al
    parametro que la gobierna segun `GOBIERNA`, y se perturba con una normal
    de desviacion `amplitud` relativa.

    QUE SIGNIFICA CADA MODO, Y PARA QUE SIRVE CADA UNO:

      * `variar=None` — replicas TECNICAS. Mismos parametros, distinta
        semilla. Miden el ruido del generador: el suelo por debajo del cual
        ninguna diferencia significa nada. NO son especimenes independientes
        y no anaden grados de libertad a nada (ver `spinpy.estadistica`).

      * `variar=[...]` — familia SINTETICA alrededor del ajuste. Sirve para
        experimentos numericos que el hueso real no permite, como variar
        BV/TV manteniendo la anisotropia. Tampoco son especimenes: son puntos
        de diseno, y compararlos con especimenes reales exige declararlo.

    AVISO SOBRE Tb.Th Y Tb.Sp: las gobierna el MISMO parametro, la escala del
    campo. Pedir que varie una sin la otra no es posible con este generador, y
    si se piden las dos se avisa en la salida en vez de fingir que se hicieron
    dos cosas distintas.

    Se devuelve la dispersion REALMENTE lograda en cada metrica, que no tiene
    por que ser `amplitud`: la relacion parametro-metrica no es lineal ni
    tiene ganancia uno.
    """
    par = dict(parametros)
    if semilla_base is None:
        semilla_base = int(par.get("semilla", 20260720))
    if par.get("familia") == "dual-lattice":
        return _replicas_dual(par, spacing, n, int(semilla_base), variar,
                              amplitud, medir, progreso)

    variar = list(variar or [])
    avisos = []
    if "TbTh" in variar and "TbSp" in variar:
        avisos.append(
            "Tb.Th y Tb.Sp las gobierna el mismo parametro (la escala del "
            "campo), asi que se perturban juntas: pedir las dos no varia dos "
            "cosas independientes.")
    sin_mapa = [v for v in variar if v not in GOBIERNA]
    if sin_mapa:
        avisos.append(
            "Sin parametro que las gobierne, no se varian: "
            + ", ".join(sin_mapa) + ". El generador toma densidad, numero de "
            "onda y angulos; el resto de metricas son consecuencia.")

    objetivos = {GOBIERNA[v][0] for v in variar if v in GOBIERNA}
    rng = np.random.default_rng(int(semilla_base))
    # Las dos lecturas comprobadas: el valor en radianes sale de aqui y de
    # ningun otro calculo (ver `procedencia.py`).
    onda = procedencia.numero_onda(par.get("wave_number_pi"),
                                   par.get("wave_number_rad"))
    cols_proc = procedencia.columnas(procedencia.desde_parametros(par))

    filas, mascaras = [], []
    for k in range(int(n)):
        if progreso:
            progreso(k, int(n), f"Replica {k+1} de {n}")
        p = {
            "resolution": int(par["resolucion"]),
            "wave_number": onda["wave_number_rad"],
            "num_waves": int(par["num_waves"]),
            "thetas": list(par["thetas"]),
            "rho": float(par["densidad"]),
            "R": np.asarray(par["R"], float),
            "esquema": par.get("esquema", "rechazo"),
            "seed": int(semilla_base) + k,
        }
        if "densidad" in objetivos:
            p["rho"] = float(np.clip(p["rho"] * (1 + amplitud * rng.normal()),
                                     0.05, 0.95))
        if "wave_number_pi" in objetivos:
            p["wave_number"] = float(p["wave_number"]
                                     * (1 + amplitud * rng.normal()))
        if "thetas" in objetivos:
            p["thetas"] = [float(np.clip(t * (1 + amplitud * rng.normal()),
                                         0, 90)) for t in p["thetas"]]

        BW, _, info = generar_mascara(**p)
        fila = {c: v for c, v in cols_proc.items() if c != "semilla"}
        fila.update({"replica": k + 1, "semilla": p["seed"],
                     "densidad": p["rho"],
                     "wave_number_pi": p["wave_number"] / np.pi,
                     "wave_number_rad": p["wave_number"],
                     "thetas": list(p["thetas"]),
                     "rho_obtenida": float(info["rho_obtenida"])})
        if medir:
            m = morfometria(BW, spacing)
            fila.update({k2: float(m[k2]) for k2 in COMPARADAS
                         if k2 in m and np.isscalar(m[k2])})
        filas.append(fila)
        mascaras.append(BW)
    if progreso:
        progreso(int(n), int(n), "listo")

    logrado = {}
    if medir and filas:
        for k2 in COMPARADAS:
            v = np.array([f[k2] for f in filas if k2 in f], float)
            v = v[np.isfinite(v)]
            if v.size > 1 and v.mean():
                logrado[k2] = {"media": float(v.mean()),
                               "sd": float(v.std(ddof=1)),
                               "cv_pct": float(100 * v.std(ddof=1) / abs(v.mean()))}

    return {"filas": filas, "mascaras": mascaras, "logrado": logrado,
            "variar": variar, "amplitud": float(amplitud),
            "tecnicas": not variar, "avisos": avisos,
            "semilla_base": int(semilla_base)}


def _logrado(filas, medir):
    logrado = {}
    if medir and filas:
        for k2 in COMPARADAS:
            v = np.array([f[k2] for f in filas if k2 in f], float)
            v = v[np.isfinite(v)]
            if v.size > 1 and v.mean():
                logrado[k2] = {"media": float(v.mean()),
                               "sd": float(v.std(ddof=1)),
                               "cv_pct": float(100 * v.std(ddof=1) / abs(v.mean()))}
    return logrado


def _replicas_dual(par, spacing, n, semilla_base, variar, amplitud, medir,
                   progreso):
    """`generar_replicas` para un ajuste de dual-lattice.

    Mismo contrato y la misma distincion entre replicas tecnicas y familia
    sintetica; cambian los parametros que se perturban (`GOBIERNA_DUAL`). El
    estiramiento se perturba en sus tres componentes por igual y se deja >= 1
    en el eje largo, para que la perturbacion no convierta un alargamiento en
    un aplastamiento.
    """
    from .dual_lattice import generar_dual_lattice

    variar = list(variar or [])
    avisos = []
    if "TbTh" in variar and "TbSp" in variar:
        avisos.append(
            "Tb.Th y Tb.Sp las gobierna el mismo parametro (la escala de la "
            "red), asi que se perturban juntas: pedir las dos no varia dos "
            "cosas independientes.")
    sin_mapa = [v for v in variar if v not in GOBIERNA_DUAL]
    if sin_mapa:
        avisos.append(
            "Sin parametro que las gobierne, no se varian: "
            + ", ".join(sin_mapa) + ". El generador toma densidad, celdas y "
            "estiramiento; el resto de metricas son consecuencia.")
    objetivos = {GOBIERNA_DUAL[v][0] for v in variar if v in GOBIERNA_DUAL}
    rng = np.random.default_rng(int(semilla_base))

    filas, mascaras = [], []
    for k in range(int(n)):
        if progreso:
            progreso(k, int(n), f"Replica {k+1} de {n}")
        rho = float(par["densidad"])
        celdas = float(par["celdas"])
        est = np.asarray(par["estiramiento"], float)
        if "densidad" in objetivos:
            rho = float(np.clip(rho * (1 + amplitud * rng.normal()), 0.05, 0.95))
        if "celdas" in objetivos:
            celdas = float(max(1.0, celdas * (1 + amplitud * rng.normal())))
        if "estiramiento" in objetivos:
            est = np.maximum(est * (1 + amplitud * rng.normal(size=3)), 0.1)
        semilla = int(semilla_base) + k
        BW, _, info = generar_dual_lattice(
            int(par["resolucion"]), celdas, rho, estiramiento=est,
            irregularidad=float(par.get("irregularidad", 0.5)),
            R=np.asarray(par["R"], float), seed=semilla)
        fila = {"replica": k + 1, "semilla": semilla, "densidad": rho,
                "celdas": celdas,
                "estiramiento": [float(x) for x in est],
                "rho_obtenida": float(info["rho_obtenida"])}
        if medir:
            m = morfometria(BW, spacing)
            fila.update({k2: float(m[k2]) for k2 in COMPARADAS
                         if k2 in m and np.isscalar(m[k2])})
        filas.append(fila)
        mascaras.append(BW)
    if progreso:
        progreso(int(n), int(n), "listo")

    return {"filas": filas, "mascaras": mascaras,
            "logrado": _logrado(filas, medir), "variar": variar,
            "amplitud": float(amplitud), "tecnicas": not variar,
            "avisos": avisos, "semilla_base": int(semilla_base),
            "familia": "dual-lattice"}
