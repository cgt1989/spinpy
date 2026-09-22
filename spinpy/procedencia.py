"""
procedencia.py — De donde sale cada numero que se guarda.

POR QUE EXISTE
--------------
Un resultado sin procedencia no se puede reproducir, y lo que no se puede
reproducir no se puede corregir. El caso que lo motiva esta medido.
`Replicacion_AppFinal_H4/resultados.json` guardaba el ajuste del VOI proximal
de H4 como `wave: 15`. Ese campo es el valor del DESLIZADOR, que tanto la app
de MATLAB como el port multiplican por pi antes de generar
(AppFinal_V2.m:1275, `fit.py`). `Estudio_Convergencia/code/run_estudio.py` lo
transcribio como `wave_number=15.0`, sin pi: la estructura estudiada era 3.3
veces mas gruesa que la ajustada, y sobre ella se construyo un problema abierto
—21 veces mas blanda, 13.8 % de islas, sin convergencia— que no describe el
ajuste.

Dos decisiones cierran esa puerta:

  1. El numero de onda se guarda SIEMPRE con sus dos lecturas y con nombres que
     no admiten duda: `wave_number_pi` (multiplo de pi, lo que muestra el
     deslizador) y `wave_number_rad` (lo que recibe el generador). La clave
     suelta `wave` no se vuelve a escribir.
  2. Hay UNA sola traduccion de un dict de parametros a una llamada al
     generador (`kwargs_generador`), y comprueba que las dos lecturas cuadran.
     Si alguien pasa el valor del deslizador como radianes, las dos lecturas
     discrepan en un factor pi y la llamada falla en vez de generar otra cosa.

Cada JSON lleva un bloque `procedencia`; cada CSV, las mismas columnas en CADA
fila. Una tabla se concatena con otras y una cabecera se pierde por el camino;
una columna no.

VERSION DEL FORMATO
-------------------
  1  todo lo escrito antes de este modulo: sin procedencia
  2  bloque `procedencia` y numero de onda con sus dos lecturas

Leer un archivo de version 1 no falla —hay meses de resultados asi— pero AVISA
(`AvisoFormatoAntiguo`) y dice que supuesto se hizo al migrar cada campo. El
aviso no se puede desactivar desde aqui a proposito: quien quiera silenciarlo
tiene que hacerlo con `warnings`, y eso queda escrito en su codigo.
"""

from __future__ import annotations

import copy
import json
import platform
import time
import warnings
from pathlib import Path

import numpy as np

FORMATO = "spinpy/resultado"
VERSION_FORMATO = 2

# Tolerancia relativa al contrastar las dos lecturas del numero de onda. Un
# redondeo de JSON o de CSV mueve la decimoquinta cifra; un pi olvidado mueve
# la primera.
TOL_ONDA = 1e-9

# Columnas de procedencia de una fila de CSV, en este orden.
COLUMNAS = ("spinpy_version", "formato_version", "familia", "esquema",
            "semilla", "wave_number_pi", "wave_number_rad")


class AvisoFormatoAntiguo(UserWarning):
    """Se leyo un archivo anterior al formato con procedencia."""


class ErrorNumeroOnda(ValueError):
    """El numero de onda falta, es ambiguo o sus dos lecturas no cuadran."""


def version_spinpy():
    from . import __version__
    return __version__


# ---------------------------------------------------------------------------
# Numero de onda
# ---------------------------------------------------------------------------

def numero_onda(pi=None, rad=None):
    """Las dos lecturas del numero de onda, comprobadas.

    pi  : multiplo de pi —el valor del deslizador y de `wave_number_pi`—.
    rad : lo que recibe `generar_mascara(wave_number=...)`.

    Con una sola se deduce la otra. Con las dos se exige `rad = pi * pi_`; si
    no cuadran se levanta `ErrorNumeroOnda`, y si el cociente es justamente pi
    el mensaje lo dice, porque es el error que ya ocurrio una vez.
    """
    if pi is None and rad is None:
        raise ErrorNumeroOnda("falta el numero de onda: ni wave_number_pi ni "
                              "wave_number_rad")
    pi = None if pi is None else float(pi)
    rad = None if rad is None else float(rad)
    if pi is None:
        pi = rad / np.pi
    elif rad is None:
        rad = pi * np.pi
    elif not np.isclose(rad, pi * np.pi, rtol=TOL_ONDA, atol=0.0):
        pista = ""
        if pi != 0 and np.isclose(rad, pi, rtol=1e-6):
            pista = (" El valor en radianes es IGUAL al multiplo de pi: falta "
                     "multiplicar por pi (el error de run_estudio.py).")
        raise ErrorNumeroOnda(
            f"las dos lecturas del numero de onda no cuadran: wave_number_pi="
            f"{pi!r} implica {pi * np.pi!r} rad, pero wave_number_rad={rad!r}."
            + pista)
    return {"wave_number_pi": pi, "wave_number_rad": rad}


# ---------------------------------------------------------------------------
# Bloque de procedencia
# ---------------------------------------------------------------------------

def bloque(familia, esquema=None, semilla=None, wave_number_pi=None,
           wave_number_rad=None, **extra):
    """Bloque `procedencia` de un resultado.

    En el dual-lattice no hay esquema de muestreo ni numero de onda: se
    escriben como "no_aplica" y None, no se omiten, para que todas las filas de
    una tabla mezclada tengan las mismas columnas y un hueco no se lea como un
    dato perdido.
    """
    fam = str(familia)
    d = {
        "formato": FORMATO,
        "version_formato": VERSION_FORMATO,
        "spinpy": version_spinpy(),
        "fecha": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "familia": fam,
        "semilla": None if semilla is None else int(semilla),
    }
    if fam == "dual-lattice":
        d.update(esquema="no_aplica", wave_number_pi=None, wave_number_rad=None)
    else:
        d["esquema"] = esquema if esquema is not None else "rechazo"
        if wave_number_pi is None and wave_number_rad is None:
            d.update(wave_number_pi=None, wave_number_rad=None)
        else:
            d.update(numero_onda(wave_number_pi, wave_number_rad))
    d.update(extra)
    return d


def desde_parametros(par):
    """Bloque de procedencia a partir del dict `parametros` de un ajuste."""
    fam = par.get("familia", "spinodoide")
    if fam == "dual-lattice":
        return bloque(fam, semilla=par.get("semilla"))
    return bloque(fam, esquema=par.get("esquema", "rechazo"),
                  semilla=par.get("semilla"),
                  wave_number_pi=par.get("wave_number_pi"),
                  wave_number_rad=par.get("wave_number_rad"))


def columnas(proc):
    """Las columnas de procedencia de una fila de CSV, desde un bloque."""
    return {"spinpy_version": proc.get("spinpy"),
            "formato_version": proc.get("version_formato"),
            "familia": proc.get("familia"),
            "esquema": proc.get("esquema"),
            "semilla": proc.get("semilla"),
            "wave_number_pi": proc.get("wave_number_pi"),
            "wave_number_rad": proc.get("wave_number_rad")}


# ---------------------------------------------------------------------------
# La unica traduccion parametros -> generador
# ---------------------------------------------------------------------------

def kwargs_generador(par, semilla=None):
    """Argumentos del generador a partir del dict `parametros` de un ajuste.

    Es el UNICO sitio que convierte un resultado guardado en una llamada a
    `generar_mascara` o a `generar_dual_lattice`. `semilla` sustituye a la del
    ajuste (replicas); None usa la guardada.

    Un dict que solo trae la clave ambigua `wave` se rechaza: no hay forma de
    saber desde aqui si es el deslizador o ya radianes. `migrar` sabe hacer el
    supuesto y lo declara; esta funcion no supone nada.
    """
    fam = par.get("familia", "spinodoide")
    s = int(par["semilla"] if semilla is None else semilla)
    R = np.asarray(par["R"], float)
    if fam == "dual-lattice":
        return dict(resolution=int(par["resolucion"]),
                    celdas=float(par["celdas"]), rho=float(par["densidad"]),
                    estiramiento=tuple(float(x) for x in par["estiramiento"]),
                    irregularidad=float(par.get("irregularidad", 0.5)),
                    R=R, seed=s)
    if ("wave" in par and "wave_number_pi" not in par
            and "wave_number_rad" not in par):
        raise ErrorNumeroOnda(
            "el dict solo trae la clave ambigua 'wave'. Pasalo antes por "
            "procedencia.migrar, que declara el supuesto (multiplo de pi).")
    onda = numero_onda(par.get("wave_number_pi"), par.get("wave_number_rad"))
    return dict(resolution=int(par["resolucion"]),
                wave_number=onda["wave_number_rad"],
                num_waves=int(par["num_waves"]),
                thetas=[float(t) for t in par["thetas"]],
                rho=float(par["densidad"]), R=R,
                esquema=par.get("esquema", "rechazo"), seed=s)


# ---------------------------------------------------------------------------
# Lectura de archivos, con migracion declarada
# ---------------------------------------------------------------------------

def _version_de(doc):
    if not isinstance(doc, dict):
        return 1
    proc = doc.get("procedencia")
    if isinstance(proc, dict) and "version_formato" in proc:
        return int(proc["version_formato"])
    if "version_formato" in doc:
        return int(doc["version_formato"])
    return 1


def _es_numero(v):
    return isinstance(v, (int, float, np.integer, np.floating)) and \
        not isinstance(v, bool) and np.isfinite(float(v))


def _migrar_nodo(nodo, cuenta, ruta):
    if isinstance(nodo, list):
        return [_migrar_nodo(v, cuenta, f"{ruta}[{i}]")
                for i, v in enumerate(nodo)]
    if not isinstance(nodo, dict):
        return nodo
    out = {k: _migrar_nodo(v, cuenta, f"{ruta}.{k}" if ruta else str(k))
           for k, v in nodo.items()}
    tiene_pi = _es_numero(out.get("wave_number_pi"))
    tiene_rad = _es_numero(out.get("wave_number_rad"))
    if tiene_pi and not tiene_rad:
        out["wave_number_rad"] = float(out["wave_number_pi"]) * np.pi
        cuenta.setdefault("pi_sin_rad", [0, ruta])[0] += 1
    elif _es_numero(out.get("wave")) and not tiene_pi and not tiene_rad:
        out["wave_number_pi"] = float(out["wave"])
        out["wave_number_rad"] = float(out["wave"]) * np.pi
        cuenta.setdefault("wave", [0, ruta])[0] += 1
    elif (_es_numero(out.get("wave_number")) and not tiene_pi
          and not tiene_rad):
        out["wave_number_rad"] = float(out["wave_number"])
        out["wave_number_pi"] = float(out["wave_number"]) / np.pi
        cuenta.setdefault("wave_number", [0, ruta])[0] += 1
    return out


def migrar(doc):
    """Lleva un documento de formato 1 al 2. Devuelve (doc, avisos).

    No modifica la entrada. Cada supuesto se cuenta una vez con el numero de
    sitios donde se aplico y la ruta del primero: una traza con 53 entradas no
    produce 53 avisos, pero tampoco uno que oculte cuantos fueron.
    """
    version = _version_de(doc)
    cuenta = {}
    out = _migrar_nodo(copy.deepcopy(doc), cuenta, "")
    avisos = []
    if version < VERSION_FORMATO:
        avisos.append(
            f"Archivo en formato {version}, anterior a la procedencia: no "
            f"consta con que version de spinpy, esquema de muestreo ni semilla "
            f"se obtuvieron sus numeros.")
    if "wave" in cuenta:
        n, r = cuenta["wave"]
        avisos.append(
            f"Clave ambigua 'wave' en {n} sitio(s) (primero: {r or 'raiz'}). "
            f"Se SUPONE multiplo de pi —la convencion del deslizador y de la "
            f"traza del ajuste— y se anaden wave_number_pi y wave_number_rad. "
            f"Si ese numero ya estaba en radianes, la estructura regenerada "
            f"sera pi veces mas fina: es el error que invalido el problema "
            f"abierto de Estudio_Convergencia.")
    if "pi_sin_rad" in cuenta:
        n, r = cuenta["pi_sin_rad"]
        avisos.append(
            f"wave_number_pi sin wave_number_rad en {n} sitio(s) (primero: "
            f"{r or 'raiz'}): se anade rad = pi x valor. No hay ambiguedad, el "
            f"nombre declara la unidad.")
    if "wave_number" in cuenta:
        n, r = cuenta["wave_number"]
        avisos.append(
            f"'wave_number' sin sus dos lecturas en {n} sitio(s) (primero: "
            f"{r or 'raiz'}): se toma como RADIANES, la convencion de "
            f"generar_mascara, que es lo que de verdad se genero.")
    return out, avisos


def leer_json(ruta):
    """Lee un resultado JSON de spinpy, migrandolo y avisando si es antiguo."""
    doc = json.loads(Path(ruta).read_text(encoding="utf-8"))
    out, avisos = migrar(doc)
    for a in avisos:
        warnings.warn(f"{Path(ruta).name}: {a}", AvisoFormatoAntiguo,
                      stacklevel=2)
    return out


def leer_csv(ruta):
    """Lee una tabla CSV de spinpy (`;` y coma decimal, o `,`), avisando.

    Devuelve un DataFrame con las columnas de procedencia completadas donde se
    pueden deducir sin suponer nada, y con el supuesto declarado donde no.
    """
    import pandas as pd

    ruta = Path(ruta)
    with open(ruta, encoding="utf-8-sig") as f:
        cabecera = f.readline()
    if ";" in cabecera:
        df = pd.read_csv(ruta, sep=";", decimal=",", encoding="utf-8-sig")
    else:
        df = pd.read_csv(ruta, encoding="utf-8-sig")

    avisos = []
    if "formato_version" not in df.columns:
        avisos.append("tabla en formato 1, anterior a la procedencia: no "
                      "consta version de spinpy, esquema ni semilla del ajuste.")
    if "wave_number_pi" in df.columns and "wave_number_rad" not in df.columns:
        df["wave_number_rad"] = df["wave_number_pi"] * np.pi
        avisos.append("columna wave_number_pi sin wave_number_rad: se anade "
                      "rad = pi x valor.")
    elif ("wave" in df.columns and "wave_number_pi" not in df.columns
          and "wave_number_rad" not in df.columns):
        df["wave_number_pi"] = df["wave"]
        df["wave_number_rad"] = df["wave"] * np.pi
        avisos.append("columna ambigua 'wave': se SUPONE multiplo de pi. Si ya "
                      "estaba en radianes, la estructura regenerada sera pi "
                      "veces mas fina.")
    for a in avisos:
        warnings.warn(f"{ruta.name}: {a}", AvisoFormatoAntiguo, stacklevel=2)
    return df


def serializable(obj):
    """Convierte arrays y escalares de numpy en tipos que `json` sabe volcar."""
    if isinstance(obj, dict):
        return {str(k): serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [serializable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return serializable(obj.tolist())
    if isinstance(obj, (np.floating, np.integer, np.bool_)):
        return obj.item()
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj
