"""
spinpy.motores — Motores de elementos finitos que corren DENTRO de spinpy.

Sustituyen a FEBio (un ejecutable externo) en los ensayos de compresion con
malla de ladrillos (hex8) o malla suave (TET10), lineales y no lineales. Todos
reciben el MISMO problema (`problema_de_malla`) y devuelven el desplazamiento
en el orden de nodos de spinpy, de modo que el postproceso (E_app, von Mises,
Pistoia) es comun y la comparacion entre motores mide solo el motor.

  app       el resolvedor propio (`resistencia`); hex8 lineal. Siempre.
  ngsolve   NGSolve (pip: ngsolve, mkl). hex8/TET10, lineal y no lineal.
  fenicsx   FEniCSx (solo conda-forge; compila en tiempo de ejecucion).
  skfem     scikit-fem (pip). hex8/TET10, lineal y no lineal; lento en NL.
  sfepy     SfePy (pip). hex8/TET10 lineal; no lineal solo SVK y con Newton
            de convergencia lineal (tangente inconsistente, medido).

La eleccion y las medidas que la justifican: `comparativa_motores/INFORME.md`.

`disponibles()` dice cuales estan instalados sin importarlos del todo;
`resolver(problema, motor)` corre en el proceso actual y
`resolver_aislado(...)` en un proceso hijo cancelable (lo usa la GUI: un fallo
de memoria o de una biblioteca en C no cierra la aplicacion, y «Detener» mata
el calculo en curso).
"""

from __future__ import annotations

import importlib
import importlib.util
import multiprocessing as mp
import time

from ._comun import (ErrorMotor, NoDisponible, esquinas,           # noqa: F401
                     problema_de_malla, rectificar_aristas)

#: Orden de presentacion. La app va primero: es la linea base validada.
MOTORES = ("app", "ngsolve", "fenicsx", "skfem", "sfepy")

#: Modulo que delata la instalacion de cada motor.
_SONDA = {"app": None, "ngsolve": "ngsolve", "fenicsx": "dolfinx",
          "skfem": "skfem", "sfepy": "sfepy"}

#: Lo que resuelve cada motor. `nl` = materiales no lineales disponibles.
CAPACIDADES = {
    "app": {"mallas": ("hex8",), "nl": (), "plato": False},
    "ngsolve": {"mallas": ("hex8", "tet10"), "nl": ("svk", "neohookeano"),
                "plato": True},
    "fenicsx": {"mallas": ("hex8", "tet10"), "nl": ("svk", "neohookeano"),
                "plato": True},
    "skfem": {"mallas": ("hex8", "tet10"), "nl": ("svk", "neohookeano"),
              "plato": True},
    "sfepy": {"mallas": ("hex8", "tet10"), "nl": ("svk",), "plato": True},
}

ETIQUETAS = {"app": "App (spinpy)", "ngsolve": "NGSolve",
             "fenicsx": "FEniCSx", "skfem": "scikit-fem", "sfepy": "SfePy"}

#: Motor recomendado para lo que la app no resuelve (TET10, no lineal),
#: segun la comparativa: exacto en todas las pruebas, instalable con pip en
#: Windows, el mas rapido con TET10 y el unico con directo multihilo sin
#: compilador. `comparativa_motores/INFORME.md`, seccion de conclusiones.
RECOMENDADO = "ngsolve"


def _modulo(motor):
    return importlib.import_module(f".m_{motor}", __name__)


def instalado(motor):
    s = _SONDA[motor]
    return s is None or importlib.util.find_spec(s) is not None


def disponibles():
    """Motores instalados en este entorno, en el orden de `MOTORES`."""
    return [m for m in MOTORES if instalado(m)]


#: Distribucion de cada motor, para leer su version SIN importarlo.
_DISTRIBUCION = {"ngsolve": "ngsolve", "fenicsx": "fenics-dolfinx",
                 "skfem": "scikit-fem", "sfepy": "sfepy"}


def version(motor):
    """Version del motor. Se lee de los metadatos del paquete: importar la
    biblioteca en el proceso de la GUI carga sus DLL junto a las de VTK y Qt,
    y los motores solo deben cargarse en el proceso hijo que resuelve."""
    dist = _DISTRIBUCION.get(motor)
    if dist is not None:
        try:
            from importlib.metadata import version as _v
            return _v(dist)
        except Exception:
            pass
    try:
        return _modulo(motor).version()
    except Exception:
        return None


def puede(motor, malla="hex8", analisis="lineal", material="svk",
          control="fuerza"):
    """True si `motor` resuelve esa combinacion (sin importarlo)."""
    c = CAPACIDADES[motor]
    if malla not in c["mallas"]:
        return False
    if control == "plato" and not c["plato"]:
        return False
    return analisis == "lineal" or material in c["nl"]


def resolver(problema, motor):
    """Resuelve en este proceso. Devuelve {u, F_reac, meta}."""
    if not instalado(motor):
        raise NoDisponible(f"{ETIQUETAS[motor]} no esta instalado")
    m = problema["meta"]
    if not puede(motor, m["tipo"], m["analisis"], m.get("material", "svk"),
                 m["control"]):
        raise NoDisponible(f"{ETIQUETAS[motor]} no resuelve {m['tipo']} "
                           f"{m['analisis']} ({m.get('material')}, "
                           f"{m['control']})")
    t0 = time.perf_counter()
    out = _modulo(motor).resolver(problema)
    out["meta"]["tiempo_total_s"] = time.perf_counter() - t0
    out["meta"]["version"] = version(motor)
    return out


class Cancelado(RuntimeError):
    """El usuario detuvo el calculo; el proceso hijo ya se termino."""


def _hijo(conexion, problema, motor, hilos):
    import os
    if hilos:
        for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                  "OPENBLAS_NUM_THREADS", "SPINPY_HILOS"):
            os.environ[k] = str(int(hilos))
    try:
        out = resolver(problema, motor)
        out["meta"]["rss_pico_MB"] = rss_pico_MB()
        conexion.send(("ok", out))
    except BaseException as e:                       # noqa: BLE001
        conexion.send(("error", (type(e).__name__, str(e))))
    finally:
        conexion.close()


def rss_pico_MB():
    """Pico de memoria residente de ESTE proceso, en MB (None si no se sabe).

    Linux: VmHWM de /proc/self/status. NO `ru_maxrss`: en Linux ese maximo se
    conserva a traves de exec, y un hijo `spawn` heredaria el del proceso de
    la GUI en el momento del fork (medido: los tres motores de una prueba
    daban el mismo valor, el de la ventana). Windows: PeakWorkingSetSize.
    """
    try:
        with open("/proc/self/status") as fh:
            for linea in fh:
                if linea.startswith("VmHWM:"):
                    return int(linea.split()[1]) / 1024.0
    except OSError:
        pass
    return _rss_windows_MB()


def _rss_windows_MB():
    """Pico de memoria del proceso en Windows (PeakWorkingSetSize)."""
    try:
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD),
                        ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t),
                        ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t),
                        ("PeakPagefileUsage", ctypes.c_size_t)]
        c = PMC()
        c.cb = ctypes.sizeof(PMC)
        ctypes.windll.psapi.GetProcessMemoryInfo(
            ctypes.windll.kernel32.GetCurrentProcess(), ctypes.byref(c),
            c.cb)
        return c.PeakWorkingSetSize / 2 ** 20
    except Exception:
        return None


def resolver_aislado(problema, motor, cancelar=None, hilos=None,
                     intervalo=0.25, tiempo_max=None):
    """`resolver` en un PROCESO HIJO (spawn), cancelable.

    `cancelar`: funcion sin argumentos o `threading.Event`; si se activa, se
    termina el hijo y se lanza `Cancelado`. Un hijo que muere sin responder
    (memoria agotada, violacion de segmento en una biblioteca en C) se
    convierte en `ErrorMotor` en lugar de cerrar la aplicacion.
    """
    if hasattr(cancelar, "is_set"):
        cancelar = cancelar.is_set
    ctx = mp.get_context("spawn")
    rx, tx = ctx.Pipe(duplex=False)
    p = ctx.Process(target=_hijo, args=(tx, problema, motor, hilos),
                    daemon=True)
    t0 = time.perf_counter()
    p.start()
    tx.close()
    try:
        while True:
            if rx.poll(intervalo):
                try:
                    estado, dato = rx.recv()
                except EOFError:
                    estado, dato = "muerto", None
                break
            if not p.is_alive():
                estado, dato = "muerto", None
                break
            if cancelar is not None and cancelar():
                raise Cancelado(f"{ETIQUETAS[motor]}: calculo detenido")
            if tiempo_max and time.perf_counter() - t0 > tiempo_max:
                raise ErrorMotor(f"{ETIQUETAS[motor]}: supero "
                                 f"{tiempo_max:.0f} s")
    finally:
        if p.is_alive():
            p.terminate()
        p.join(5)
    if estado == "ok":
        return dato
    if estado == "muerto":
        raise ErrorMotor(f"{ETIQUETAS[motor]}: el proceso de calculo termino "
                         f"sin responder (codigo {p.exitcode}); suele ser "
                         "memoria agotada")
    tipo, msg = dato
    if tipo == "NoDisponible":
        raise NoDisponible(msg)
    raise ErrorMotor(f"{ETIQUETAS[motor]}: {tipo}: {msg}")
