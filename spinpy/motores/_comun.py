"""
_comun.py — Lo que comparten los adaptadores de motor.

EL PROBLEMA
-----------
Un motor recibe un dict `problema` (ver `problema_de_malla`) y devuelve
`Solucion`: el desplazamiento nodal (N, 3) en el ORDEN DE NODOS DE SPINPY y la
fuerza de reaccion en el techo por paso de carga. Unidades mm - N - MPa.

  nodos, elems      (N, 3) mm; (M, 8) hex8 C3D8 o (M, 10) TET10 C3D10
  caras_techo       quad4 / tri6 cargadas (normal +z)
  caras_base        quad4 / tri6 de la base (esquinas primero)
  base_nodos        todos los nodos del plano z0 (intermedios incluidos)
  techo_nodos       todos los nodos del plano z1
  meta              tipo, apoyo, E, nu, H, A_bruta, ancla_xy, ancla_y y:
      analisis      'lineal' | 'nl'
      control       'fuerza': traccion p uniforme sobre las caras del techo
                    (la tension aparente por la seccion BRUTA, repartida
                    sobre el area osea, como la app); 'plato': uz = -eps H
                    impuesto en todo el techo, sin friccion.
      cargas        lista de p (MPa, control por fuerza) o de eps (plato),
                    una por paso; el lineal usa la ultima.
      material      'svk' | 'neohookeano' (solo no lineal)
      solver        'directo' | 'iterativo' | 'auto'
  BW, spacing       solo para el motor 'app' (su ensayo parte de la mascara)

La fuerza de reaccion se define igual en todos: F = -sum(f_int,z) sobre los
GDL z del techo, con f_int el vector de fuerzas internas en el estado
resuelto (positiva en compresion). En lineal con control por fuerza es la
carga aplicada; con plato, la que hace el plato.
"""

from __future__ import annotations

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


class NoDisponible(NotImplementedError):
    """El motor no resuelve esta combinacion (elemento, analisis, material)."""


class ErrorMotor(RuntimeError):
    """El motor fallo o no convergio."""


def modos_rigidos(coord):
    """Seis modos de solido rigido (3N, 6), GDL entrelazados x, y, z."""
    n = coord.shape[0]
    return modos_rigidos_gdl(np.repeat(coord, 3, axis=0),
                             np.tile(np.arange(3), n))


def modos_rigidos_gdl(xyz, comp):
    """Modos rigidos para GDL con numeracion arbitraria.

    `xyz` (n, 3): posicion del nodo de cada GDL; `comp` (n,): su componente.
    Tres traslaciones y tres giros infinitesimales: el espacio nulo de la
    rigidez sin apoyos, que el multigrid algebraico necesita representar en
    cada nivel grueso (ver el comentario del AMG en `resistencia.py`).
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
    """Sistema disperso SDP con scipy/pyamg (motores sin resolvedor propio).

    directo    SuperLU (`scipy.sparse.linalg.splu`), ordenacion COLAMD.
    iterativo  CG precondicionado con multigrid algebraico de agregacion
               suavizada (pyamg) y los modos rigidos `B` como espacio nulo
               cercano: la receta de `resistencia.ensayo_compresion`.
    Devuelve (x, info) con el residuo relativo ||b - A x|| / ||b||; si no
    alcanza 100 tol se lanza `ErrorMotor` (un iterativo no avisa solo).
    """
    from scipy.sparse.linalg import splu
    A = A.tocsr()
    if modo == "directo":
        x = splu(A.tocsc()).solve(b)
        nombre, it = "SuperLU (scipy)", None
    else:
        import pyamg
        ml = pyamg.smoothed_aggregation_solver(A, B=B, max_coarse=500)
        res = []
        x = ml.solve(b, tol=tol, maxiter=2000, accel="cg", residuals=res)
        nombre = "CG + AMG agregacion suavizada (pyamg, modos rigidos)"
        it = len(res) - 1
    nb = np.linalg.norm(b)
    r = float(np.linalg.norm(b - A @ x) / (nb if nb > 0 else 1.0))
    if not np.isfinite(r) or (modo != "directo" and r > 100 * tol):
        raise ErrorMotor(f"El resolvedor no convergio (residuo {r:.1e}).")
    return x, {"solver": nombre, "residuo_rel": r, "iteraciones": it}


def emparejar(coord_motor, coord_spinpy, tol_rel=1e-9):
    """Indice `i` tal que coord_motor[i[k]] == coord_spinpy[k].

    Los motores renumeran nodos y GDL; se vuelve al orden de spinpy por
    coordenadas. Si un nodo no tiene pareja se aborta: la solucion no se
    podria atribuir a la malla de spinpy.
    """
    from scipy.spatial import cKDTree
    L = float(np.ptp(coord_spinpy, axis=0).max())
    d, i = cKDTree(coord_motor).query(coord_spinpy)
    if d.max() > tol_rel * L:
        raise ErrorMotor(f"nodos sin pareja: distancia maxima {d.max():.3e} "
                         f"(L = {L:.3e})")
    return i


def esquinas(elems):
    """Nodos de esquina y conectividad renumerada sobre ellos."""
    if elems.shape[1] == 8:
        return np.arange(int(elems.max()) + 1), elems
    esq = np.unique(elems[:, :4])
    return esq, np.searchsorted(esq, elems[:, :4])


def lame(E, nu):
    """(lambda, mu) de Lame a partir de E y nu."""
    return E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))


def rectificar_aristas(nodos, elems):
    """TET10 con cada nodo intermedio en el punto medio de su arista.

    `fem.preparar_tet10` devuelve a su plano exacto los nodos a menos de
    1e-3 h de un plano del cubo; si una esquina se mueve y la otra no, el nodo
    intermedio deja de estar en el punto medio (medido: hasta 8e-4 h). Los
    motores construyen P2 sobre las esquinas, con aristas RECTAS: para que
    todos resuelvan la misma geometria se recoloca cada intermedio. Devuelve
    (nodos, desplazamiento maximo).
    """
    nodos = np.array(nodos, float)
    if elems.shape[1] != 10:
        return nodos, 0.0
    mov = 0.0
    for k, (i, j) in enumerate(((0, 1), (1, 2), (0, 2), (0, 3), (1, 3),
                                (2, 3))):
        med = 0.5 * (nodos[elems[:, i]] + nodos[elems[:, j]])
        mov = max(mov, float(np.abs(nodos[elems[:, 4 + k]] - med).max()))
        nodos[elems[:, 4 + k]] = med
    return nodos, mov


def problema_de_malla(malla, E, nu, apoyo="deslizante", analisis="lineal",
                      control="fuerza", cargas=(1.0,), material="svk",
                      solver="auto", tol=1e-10, BW=None):
    """Problema de motor a partir de una malla de `fem.mallar`.

    E en MPa, cargas en MPa (control por fuerza: tension aparente sobre la
    seccion BRUTA; aqui se convierte en la presion sobre el area osea del
    techo) o adimensionales (plato: deformacion aparente).
    """
    from ..escribe import area_caras
    nodos, elems = malla["nodos"], np.asarray(malla["elems"], np.int64)
    nodos, mov = rectificar_aristas(nodos, elems)
    z = nodos[:, 2]
    z0, z1 = float(z.min()), float(z.max())
    tol_z = 1e-9 * max(z1 - z0, 1.0)
    base = np.nonzero(z <= z0 + tol_z)[0]
    techo = np.nonzero(z >= z1 - tol_z)[0]
    if elems.shape[1] == 8:
        cb = elems[np.all(z[elems[:, 0:4]] <= z0 + tol_z, axis=1), 0:4]
        esq_base = base
    else:
        from ..fem import _caras_planos
        cb = _caras_planos(dict(malla, nodos=nodos))[(2, 0)][0]
        esq_base = np.intersect1d(base, np.unique(elems[:, :4]))
    ct = np.asarray(malla["caras_techo"], np.int64)
    cbx = nodos[esq_base]
    a = int(esq_base[int(np.argmin(cbx[:, 0] + cbx[:, 1]))])
    b = int(esq_base[int(np.argmax(cbx[:, 0] - cbx[:, 1]))])
    A_osea = float(area_caras(nodos, ct).sum())
    A_bruta = float(malla["A_bruta"])
    cargas = [float(c) for c in np.atleast_1d(cargas)]
    if control == "fuerza":
        cargas = [c * A_bruta / A_osea for c in cargas]
    m = {"tipo": "hex8" if elems.shape[1] == 8 else "tet10",
         "apoyo": apoyo, "E": float(E), "nu": float(nu),
         "H": float(z1 - z0), "A_bruta": A_bruta, "A_osea_techo": A_osea,
         "ancla_xy": a, "ancla_y": b, "analisis": analisis,
         "control": control, "cargas": cargas, "material": material,
         "solver": solver, "tol": float(tol), "rectificacion_max_mm": mov,
         "n_nodos": int(nodos.shape[0]), "n_elems": int(elems.shape[0]),
         "n_gdl": int(3 * nodos.shape[0])}
    if BW is not None:
        m["spacing"] = np.asarray(malla["spacing"], float).tolist()
    p = {"nodos": nodos, "elems": elems, "caras_techo": ct, "caras_base": cb,
         "base_nodos": base, "techo_nodos": techo, "meta": m}
    if BW is not None:
        p["BW"] = np.asarray(BW, bool)
    return p


def fijos_basicos(p):
    """[(nodo, componente)] del apoyo (sin el plato)."""
    m = p["meta"]
    base = p["base_nodos"]
    if m["apoyo"] == "empotrado":
        return [(int(n), c) for n in base for c in range(3)]
    return ([(int(n), 2) for n in base]
            + [(m["ancla_xy"], 0), (m["ancla_xy"], 1), (m["ancla_y"], 1)])


def solucion(u, F, meta, **extra):
    """Salida comun de un motor."""
    out = {"u": np.asarray(u, float), "F_reac": np.asarray(F, float),
           "meta": meta}
    out.update(extra)
    return out
