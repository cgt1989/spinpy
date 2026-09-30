"""
_base.py — Lo comun a los cinco adaptadores de motor.

Cada adaptador (`app`, `skfem`, `sfepy`, `ngsolve`, `fenicsx`) recibe el MISMO
problema discreto en un .npz (nodos, conectividad, caras cargadas, apoyos,
material) y devuelve el desplazamiento nodal en el ORDEN DE NODOS DE SPINPY.
Asi las diferencias entre motores solo pueden venir del motor (montaje,
cuadratura, condiciones de contorno, resolvedor), nunca de la malla.

Este modulo solo usa numpy y la biblioteca estandar: tiene que importarse
igual en el entorno de pip (app, scikit-fem, SfePy, NGSolve) y en el de
conda-forge (FEniCSx), que no tiene spinpy instalado.

Protocolo de una corrida (lo lanza `correr.py` en un PROCESO HIJO, para que
la memoria de pico sea la del motor y no la del lanzador):

    python -m motores.<motor> problema.npz salida.npz

La salida lleva `u` (N, 3) y, en no lineal, la fuerza de reaccion por paso;
el JSON con tiempos y memoria va en `salida.npz` como `meta` (texto).
"""

from __future__ import annotations

import json
import resource
import sys
import time

import numpy as np


class Cronometro:
    """Tiempos por etapa (s, reloj de pared)."""

    def __init__(self):
        self.t = {}
        self._t0 = None
        self._k = None

    def __call__(self, etapa):
        ahora = time.perf_counter()
        if self._k is not None:
            self.t[self._k] = self.t.get(self._k, 0.0) + ahora - self._t0
        self._k, self._t0 = etapa, ahora
        return self

    def fin(self):
        self(None)
        return dict(self.t)


def rss_pico_MB():
    """Pico de memoria residente del proceso (Linux: VmHWM).

    No `ru_maxrss`: en Linux ese maximo sobrevive a exec y el hijo heredaria
    el del lanzador en el momento del fork.
    """
    try:
        with open("/proc/self/status") as fh:
            for linea in fh:
                if linea.startswith("VmHWM:"):
                    return int(linea.split()[1]) / 1024.0
    except OSError:
        pass
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def leer_problema(ruta):
    d = np.load(ruta, allow_pickle=False)
    p = {k: d[k] for k in d.files}
    p["meta"] = json.loads(str(p["meta"]))
    return p


def modos_rigidos(coord):
    """Seis modos de solido rigido (3N, 6), GDL entrelazados x, y, z."""
    n = coord.shape[0]
    x, y, z = coord[:, 0], coord[:, 1], coord[:, 2]
    B = np.zeros((3 * n, 6))
    B[0::3, 0] = 1.0
    B[1::3, 1] = 1.0
    B[2::3, 2] = 1.0
    B[0::3, 3] = -y
    B[1::3, 3] = x
    B[1::3, 4] = -z
    B[2::3, 4] = y
    B[0::3, 5] = z
    B[2::3, 5] = -x
    return B


def modos_rigidos_gdl(xyz, comp):
    """Modos rigidos para GDL con numeracion arbitraria.

    `xyz` (n, 3): posicion del nodo de cada GDL; `comp` (n,): su componente
    (0, 1, 2). Lo necesitan los motores que no entrelazan x, y, z.
    """
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    B = np.zeros((xyz.shape[0], 6))
    for c in range(3):
        B[comp == c, c] = 1.0
    B[comp == 0, 3] = -y[comp == 0]
    B[comp == 1, 3] = x[comp == 1]
    B[comp == 1, 4] = -z[comp == 1]
    B[comp == 2, 4] = y[comp == 2]
    B[comp == 0, 5] = z[comp == 0]
    B[comp == 2, 5] = -x[comp == 2]
    return B


def resolver_scipy(A, b, B=None, modo="iterativo", tol=1e-10):
    """Sistema disperso SDP con scipy/pyamg (para los motores sin resolvedor).

    directo    SuperLU (`scipy.sparse.linalg.splu`) con ordenacion COLAMD.
    iterativo  CG precondicionado con multigrid algebraico de agregacion
               suavizada (pyamg), con los modos rigidos `B` como espacio nulo
               cercano: la misma receta que `resistencia.ensayo_compresion`.
    Devuelve (x, info) con el residuo relativo ||b - A x|| / ||b||.
    """
    from scipy.sparse.linalg import splu
    A = A.tocsr()
    if modo == "directo":
        x = splu(A.tocsc()).solve(b)
        nombre = "SuperLU (scipy)"
        it = None
    else:
        import pyamg
        ml = pyamg.smoothed_aggregation_solver(A, B=B, max_coarse=500)
        res = []
        x = ml.solve(b, tol=tol, maxiter=1000, accel="cg", residuals=res)
        nombre = "CG + AMG agregacion suavizada (pyamg, modos rigidos)"
        it = len(res) - 1
    nb = np.linalg.norm(b)
    r = float(np.linalg.norm(b - A @ x) / (nb if nb > 0 else 1.0))
    return x, {"solver": nombre, "residuo_rel": r, "iteraciones": it}


def gdl_fijos(p):
    """GDL restringidos (indice 3*nodo + componente) del apoyo del problema.

    deslizante: uz = 0 en toda la base; ux = uy = 0 en `ancla_xy`; uy = 0 en
    `ancla_y` (lo minimo para quitar los movimientos de solido rigido).
    empotrado: los tres en toda la base. Con plato, uz prescrito en el techo.
    """
    base = p["base_nodos"].astype(np.int64)
    if p["meta"]["apoyo"] == "empotrado":
        f = np.concatenate([3 * base, 3 * base + 1, 3 * base + 2])
    else:
        a, b = int(p["meta"]["ancla_xy"]), int(p["meta"]["ancla_y"])
        f = np.concatenate([3 * base + 2, [3 * a, 3 * a + 1, 3 * b + 1]])
    return np.unique(f)


def emparejar(coord_motor, coord_spinpy, tol_rel=1e-9):
    """Indice `i` tal que coord_motor[i[k]] == coord_spinpy[k].

    Los motores renumeran nodos y GDL; se vuelve al orden de spinpy por
    coordenadas. Se exige coincidencia a `tol_rel` del tamano de la malla: si
    un nodo no tiene pareja la comparacion no tendria sentido y se aborta.
    """
    from scipy.spatial import cKDTree
    L = float(np.ptp(coord_spinpy, axis=0).max())
    d, i = cKDTree(coord_motor).query(coord_spinpy)
    if d.max() > tol_rel * L:
        raise RuntimeError(f"nodos sin pareja: distancia maxima {d.max():.3e} "
                           f"(L = {L:.3e})")
    return i


def escribir_salida(ruta, u, meta, **extra):
    meta = dict(meta)
    meta["rss_pico_MB"] = rss_pico_MB()
    np.savez(ruta, u=np.asarray(u, float), meta=json.dumps(meta), **extra)


def resolver_paquete(p, motor):
    """Caso del banco -> problema de `spinpy.motores` -> (u, meta, extra)."""
    import os
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from spinpy import motores
    m = p["meta"]
    if m["analisis"] == "lineal":
        m.setdefault("control", "fuerza")
        m.setdefault("cargas", [m["p"]])
    else:
        m.setdefault("control", "plato")
        m.setdefault("cargas", list(m["eps_plato"]))
    os.environ.setdefault("SPINPY_HILOS", os.environ.get("OMP_NUM_THREADS",
                                                         "1"))
    out = motores.resolver(p, motor)
    return out["u"], out["meta"], {"F_reac": out["F_reac"]}


def principal(resolver):
    """Punto de entrada comun: `resolver(problema) -> (u, meta, extra)`."""
    entrada, salida = sys.argv[1], sys.argv[2]
    rss0 = rss_pico_MB()
    p = leer_problema(entrada)
    u, meta, extra = resolver(p)
    meta["rss_importacion_MB"] = rss0
    escribir_salida(salida, u, meta, **(extra or {}))
