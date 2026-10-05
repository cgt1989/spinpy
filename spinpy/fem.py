"""
fem.py: Ensayos de elementos finitos de spinpy con motores INTERNOS.

Sustituye a la integracion con FEBio (`febio4.exe`, un proceso externo). El
mallado, los protocolos, el postproceso y las tablas son los mismos que se
validaron contra FEBio (`comparativa_febio/`, `comparativa_febio_tet/`); lo que
cambia es quien resuelve: uno o varios de los motores de `spinpy.motores`
(la app, NGSolve, FEniCSx, scikit-fem, SfePy), elegidos por el usuario y
comparados en `comparativa_motores/INFORME.md`.

Dos mallas posibles sobre la MISMA mascara:

  hex8    la malla de voxeles de la app, un elemento por voxel. Con el motor
          'app' es exactamente `resistencia.ensayo_compresion`; con los demas
          motores es el mismo problema discreto (la comparativa lo verifica a
          ~1e-10), y ademas admite plato rigido y no linealidad geometrica.
  tet10   una malla tetraedrica cuadratica SUAVE (`solido.malla_tet10`) que
          quita los escalones. La app no la resuelve; si los demas motores.

Con varios motores sobre la misma malla y el mismo protocolo, la diferencia
entre motores mide la IMPLEMENTACION (debe ser ~0); la diferencia entre hex8
y tet10 con un mismo motor mide solo el efecto de la MALLA.

DIFERENCIAS CON LA RUTA DE FEBIO (declaradas)
---------------------------------------------
* Lineal: se resuelve el problema LINEAL directamente, una sola vez. FEBio
  es no lineal geometricamente y habia que extrapolar a carga nula con dos
  cargas pequenas (u = 2 u(s) - u(2 s)).
* Fuerza de reaccion: sale del vector de fuerzas internas en el techo. FEBio
  4.5 escribia reaccion cero en los `zero displacement` y se reconstruia por
  la integral de volumen de sigma_zz, que aqui queda como COMPROBACION de
  equilibrio (`dF_rel`).
* No lineal con fuerza: la traccion es MUERTA (fija en la configuracion de
  referencia). La presion de FEBio seguia a la cara deformada (seguidora); la
  diferencia es del orden de la rotacion de las caras del techo, pequena a
  las cargas de los protocolos, pero no nula.

LA GUI SOLO LLAMA A ESTE MODULO: todo corre igual desde la CLI y los tests.

Unidades: E y sigma viajan en Pa, como en el resto de spinpy; los motores
trabajan en mm - N - MPa y la conversion se hace aqui.
"""

from __future__ import annotations

import hashlib
import os
import time
from pathlib import Path

import numpy as np

from . import motores
from .motores import Cancelado, ErrorMotor, NoDisponible      # noqa: F401

#: Motor por omision: la app con ladrillos (el validado); para lo que la app
#: no resuelve, el recomendado por la comparativa.
MOTOR_DEF = "app"


# ---------------------------------------------------------------------------
# Magnitudes derivadas del campo de tension
# ---------------------------------------------------------------------------

def von_mises(s):
    """von Mises = sqrt(3 J2) de tensiones Voigt [xx yy zz yz xz xy]."""
    s = np.asarray(s, float)
    return np.sqrt(np.maximum(
        0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2
               + (s[:, 2] - s[:, 0]) ** 2)
        + 3.0 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2), 0.0))


def eps_eff(sig, E, nu):
    """Deformacion efectiva de Pistoia a partir de la tension (Voigt, Pa).

    spinpy la calcula de sigma y eps; aqui eps sale de sigma por la
    flexibilidad isotropa, con distorsiones de ingenieria. Es la misma
    energia: U = sigma . S sigma / 2, eps_eff = sqrt(2 U / E).
    """
    S = np.zeros((6, 6))
    S[:3, :3] = -nu / E
    np.fill_diagonal(S[:3, :3], 1.0 / E)
    S[3, 3] = S[4, 4] = S[5, 5] = 2.0 * (1.0 + nu) / E
    sig = np.asarray(sig, float)
    eps = sig @ S.T
    U = 0.5 * np.einsum("ij,ij->i", sig, eps)
    return np.sqrt(np.maximum(2.0 * U / E, 0.0))



def huella(nodos, elems):
    """SHA-256 de la malla (coordenadas float64 y conectividad int64, en C)."""
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(np.asarray(nodos, np.float64)).tobytes())
    h.update(np.ascontiguousarray(np.asarray(elems, np.int64)).tobytes())
    return h.hexdigest()




# ---------------------------------------------------------------------------
# Mallas
# ---------------------------------------------------------------------------

MALLAS = ("hex8", "tet10")

# Caras de un tetraedro por sus esquinas locales, cada una con la esquina
# opuesta, y aristas C3D10 (ranuras 4..9) para montar las tri6.
_CARAS_TET = ((1, 2, 3, 0), (0, 3, 2, 1), (0, 1, 3, 2), (0, 2, 1, 3))
_ARISTA_C3D10 = {(0, 1): 4, (1, 2): 5, (0, 2): 6, (0, 3): 7, (1, 3): 8,
                 (2, 3): 9}

#: Perdida de volumen de la malla suave, en % del volumen de VOXELES del hueso
#: portante, a partir de la cual el registro va con reservas. DECLARADO ANTES
#: DE MEDIR (2026-09-24), a partir de lo que ya estaba medido en `solido.py`
#: (-5 % frente a voxeles, -1.9 % frente a la superficie cruda, spinodoide a
#: 48^3): la rigidez escala ~ rho^2, asi que un 3 % de volumen son ~6 % de
#: rigidez que no vienen de los escalones. Es NUESTRO criterio.
PERDIDA_VOLUMEN_MAX_PCT = 3.0

#: Hueso descartado al exigir que la malla una base y techo, en % del volumen
#: de la malla, a partir del cual hay reservas. Mismo umbral que el 1 % de
#: `informe.DESCONEXION_RESERVAS`, por coherencia. NUESTRO criterio.
DESCARTE_MAX_PCT = 1.0


def _volumen_tet(nodos, elems):
    P = nodos[elems[:, :4]]
    return np.einsum("ij,ij->i", P[:, 1] - P[:, 0],
                     np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 0])) / 6.0


def caras_borde_tet10(elems):
    """Caras de borde de una malla TET10: (elemento, tri6 saliente).

    Una cara de esquinas que aparece en un solo tetraedro es de borde. Se
    orienta con la normal hacia FUERA del tetraedro (la esquina opuesta queda
    detras), que es lo que necesita la traccion en el techo, y se completa con sus
    nodos intermedios en el orden tri6 (esquinas a, b, c; aristas ab, bc, ca).
    """
    elems = np.asarray(elems, dtype=np.int64)
    n = elems.shape[0]
    tri, dueno = [], []
    for a, b, c, _ in _CARAS_TET:
        tri.append(elems[:, [a, b, c]])
        dueno.append(np.arange(n))
    tri = np.concatenate(tri)
    dueno = np.concatenate(dueno)
    loc = np.concatenate([np.tile(np.array(f[:3]), (n, 1)) for f in _CARAS_TET])
    clave = np.sort(tri, axis=1)
    _, inv, cuenta = np.unique(clave, axis=0, return_inverse=True,
                               return_counts=True)
    borde = cuenta[inv.ravel()] == 1
    tri, dueno, loc = tri[borde], dueno[borde], loc[borde]

    def medio(i, j):
        k = np.array([_ARISTA_C3D10[tuple(sorted((int(p), int(q))))]
                      for p, q in zip(i, j)])
        return elems[dueno, k]

    m_ab = medio(loc[:, 0], loc[:, 1])
    m_bc = medio(loc[:, 1], loc[:, 2])
    m_ca = medio(loc[:, 2], loc[:, 0])
    return dueno, np.column_stack([tri, m_ab, m_bc, m_ca])


def _orientar_saliente(nodos, elems, dueno, caras):
    """Invierte las caras cuya normal apunta hacia dentro de su elemento."""
    P = nodos[caras[:, :3]]
    n = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    c_el = nodos[elems[dueno, :4]].mean(axis=1)
    dentro = np.einsum("ij,ij->i", n, c_el - P[:, 0]) > 0
    caras = caras.copy()
    # (a, b, c, ab, bc, ca) -> (a, c, b, ca, bc, ab)
    caras[dentro] = caras[dentro][:, [0, 2, 1, 5, 4, 3]]
    return caras


def _componentes_portantes(nodos, elems, z0, z1, tol):
    """Elementos de los componentes que tocan base Y techo."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    e4 = elems[:, :4]
    n = nodos.shape[0]
    filas = np.repeat(e4[:, 0], 3)
    cols = e4[:, 1:].ravel()
    G = coo_matrix((np.ones(filas.size), (filas, cols)), shape=(n, n))
    ncomp, lab = connected_components(G, directed=False)
    z = nodos[:, 2]
    en_base = set(np.unique(lab[z <= z0 + tol]))
    en_techo = set(np.unique(lab[z >= z1 - tol]))
    buenos = np.array(sorted(en_base & en_techo), dtype=np.int64)
    usados = np.unique(e4)
    n_comp = int(np.unique(lab[usados]).size)
    return np.isin(lab[e4[:, 0]], buenos), n_comp


def _compactar(nodos, elems):
    usados, inv = np.unique(elems, return_inverse=True)
    return nodos[usados], inv.reshape(elems.shape)


def opciones_malla(**kw):
    """Opciones de la malla suave con sus valores por omision, validadas."""
    from .solido import DECIMADO, SUAVIZADO_BANDA, SUAVIZADO_ITER
    o = {"suavizado": SUAVIZADO_ITER, "banda": SUAVIZADO_BANDA,
         "decimado": DECIMADO, "tam_max_mm": None, "minratio": None,
         "corregir_volumen": True}
    for k, v in kw.items():
        if k not in o:
            raise ValueError(f"opcion de malla desconocida: {k}")
        o[k] = v
    return o


def maxvolume_de_tamano(L):
    """Volumen del tetraedro regular de arista L: L^3 / (6 sqrt 2)."""
    return None if not L else float(L) ** 3 / (6.0 * np.sqrt(2.0))


def preparar_tet10(nodos, elems, forma, spacing):
    """Post-proceso de una malla TET10 dentro del cubo de la mascara.

    El cubo es el de `solido.superficie_cerrada`: planos en -h/2 y
    (n - 1/2) h en cada eje. Se exige que la MALLA una base y techo (se
    descartan los componentes que no toquen los dos planos z, y se reporta el
    volumen perdido), se extraen las caras de borde tri6 salientes, las del
    techo (cargadas) y la capa superficial: elementos con una cara de borde
    que NO este sobre un plano del cubo.

    Devuelve (nodos, elems, vol, caras_techo, superficie, z0, z1, informe).
    Sirve igual para una malla que no salga de voxeles (la validacion con la
    esfera analitica usa esta misma funcion).
    """
    spacing = np.asarray(spacing, float)
    z0, z1 = -0.5 * spacing[2], (forma[2] - 0.5) * spacing[2]
    tol = 1e-6 * float(spacing[2])
    # Los nodos vienen de marching cubes en FLOAT32: a z ~ 5 mm el techo queda
    # a ~4e-7 mm de su plano (medido en el VOI proximal de H4 a 48^3), mas que
    # cualquier tolerancia razonable en unidades de h. Se devuelven a su plano
    # EXACTO, en float64, los nodos a menos de 1e-3 h de el. Una arista con
    # los dos extremos en el plano tiene alli tambien su nodo intermedio.
    nodos = np.array(nodos, dtype=np.float64, copy=True)
    for e in range(3):
        for v in (-0.5 * spacing[e], (forma[e] - 0.5) * spacing[e]):
            cerca = np.abs(nodos[:, e] - v) <= 1e-3 * spacing[e]
            nodos[cerca, e] = v
    # ASTILLAS PLANAS. tetgen deja a veces tetraedros con sus cuatro esquinas
    # sobre una misma cara del cubo (volumen ~1e-10 de la mediana, medido en
    # el espinodoide de referencia a 32^3: 29 elementos). Al devolver los
    # nodos a su plano exacto su volumen pasa a cero y la malla se rechazaba
    # (tambien por la ruta de FEBio). Esas astillas no aportan volumen ni
    # rigidez: sus cuatro caras cubren dos veces el mismo cuadrilatero de la
    # cara del cubo. Se quitan y la cara expuesta pasa a ser la del elemento
    # interior contiguo; se declara cuantas.
    en_mismo_plano = np.zeros(elems.shape[0], bool)
    for e in range(3):
        for v in (-0.5 * spacing[e], (forma[e] - 0.5) * spacing[e]):
            en_mismo_plano |= np.all(nodos[elems[:, :4], e] == v, axis=1)
    n_astillas = int(en_mismo_plano.sum())
    if n_astillas:
        elems = elems[~en_mismo_plano]
        nodos, elems = _compactar(nodos, elems)
    vol = _volumen_tet(nodos, elems)
    if (vol <= 0).any():
        raise RuntimeError(f"{int((vol <= 0).sum())} tetraedros con "
                           "volumen no positivo.")
    V_total = float(vol.sum())
    port, n_comp = _componentes_portantes(nodos, elems, z0, z1, tol)
    elems, vol = elems[port], vol[port]
    nodos, elems = _compactar(nodos, elems)
    V = float(vol.sum())
    dueno, caras_b = caras_borde_tet10(elems)
    caras_b = _orientar_saliente(nodos, elems, dueno, caras_b)
    # Caras sobre los planos del cubo: todas sus esquinas en el plano.
    en_plano = np.zeros(caras_b.shape[0], bool)
    techo = np.zeros(caras_b.shape[0], bool)
    for e in range(3):
        a, b = -0.5 * spacing[e], (forma[e] - 0.5) * spacing[e]
        c = nodos[caras_b[:, :3], e]
        t = 1e-6 * float(spacing[e])
        en_a = np.all(np.abs(c - a) <= t, axis=1)
        en_b = np.all(np.abs(c - b) <= t, axis=1)
        en_plano |= en_a | en_b
        if e == 2:
            techo = en_b
    caras = caras_b[techo]
    sup = np.zeros(elems.shape[0], bool)
    sup[dueno[~en_plano]] = True
    V_caja = float(np.prod(np.asarray(forma) * spacing))
    inf = {"n_componentes": n_comp, "V_malla_total_mm3": V_total,
           "V_malla_mm3": V,
           "descartado_pct": 100.0 * (V_total - V) / V_total,
           "BVTV_malla": V / V_caja,
           "vol_elem_min_mm3": float(vol.min()),
           "astillas_planas_quitadas": n_astillas,
           "vol_elem_mediana_mm3": float(np.median(vol))}
    return nodos, elems, vol, caras, sup, z0, z1, inf


def malla_de_tet10(nodos, elems, forma, spacing, informe=None):
    """Dict de malla (como el de `mallar`) a partir de un TET10 cualquiera."""
    spacing = np.asarray(spacing, float)
    nodos, elems, vol, caras, sup, z0, z1, inf = preparar_tet10(
        np.asarray(nodos, float), np.asarray(elems, np.int64), forma, spacing)
    inf = {**(informe or {}), **inf, "tipo": "tet10",
           "n_nodos": int(nodos.shape[0]), "n_elems": int(elems.shape[0]),
           "n_superficie": int(sup.sum()), "n_caras_techo": int(caras.shape[0])}
    return {"tipo": "tet10", "nodos": nodos, "elems": elems,
            "caras_techo": caras, "vol_elem": vol, "superficie": sup,
            "A_bruta": float(forma[0] * spacing[0] * forma[1] * spacing[1]),
            "H": float(z1 - z0), "z0": z0, "z1": z1, "forma": tuple(forma),
            "spacing": spacing, "eje": 2, "huella": huella(nodos, elems),
            "informe": inf}


def mallar(BW, spacing, tipo="hex8", eje=2, progreso=None, **opciones):
    """Mascara -> malla lista para los motores, con el eje de carga en z.

    Siempre se filtra PRIMERO a lo portante (`resistencia._solo_portante`,
    el mismo filtro del ensayo de la app): un fragmento que no une base y
    techo solo aporta una submatriz casi singular. Para cargar en X o Y se
    PERMUTA la mascara como `ensayo_compresion_eje` y se malla el resultado:
    la carga siempre va en z.

    Devuelve un dict con
      tipo, nodos (N, 3) mm, elems (M, 8|10) base 0,
      caras_techo   quad4 / tri6 cargadas, normal +z,
      vol_elem      volumen de cada elemento (mm^3),
      superficie    elementos de la capa superficial (bool, M): los que
                    tienen una cara en el borde libre del hueso. Las caras
                    sobre los seis planos del cubo NO cuentan (la regla de
                    las isocaps, la misma de `resistencia.capa_superficie`).
      A_bruta, H, z0, z1, forma, spacing, eje, huella (SHA-256), informe.

    Con TET10, `informe` declara la perdida de volumen frente a los voxeles
    del hueso portante y frente a la superficie cruda de marching cubes, y el
    volumen descartado al exigir que la MALLA (no la mascara) una base y
    techo: el suavizado puede cortar un puntal de un voxel.
    """
    from .resistencia import _PERM, _solo_portante, capa_superficie
    from .solido import malla_hex, malla_tet10

    if tipo not in MALLAS:
        raise ValueError(f"malla desconocida: {tipo!r}")
    BW = np.asarray(BW, dtype=bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    p = _PERM[int(eje)]
    BW = np.transpose(BW, p)
    spacing = spacing[list(p)]
    nx, ny, nz = BW.shape

    BW_res = _solo_portante(BW)
    if not BW_res[:, :, 0].any() or not BW_res[:, :, -1].any():
        raise ValueError("No hay hueso que una la base con el techo.")
    frac_port = float(BW_res.sum() / max(BW.sum(), 1))
    V_vox = float(BW_res.sum() * np.prod(spacing))
    A_bruta = float(nx * spacing[0] * ny * spacing[1])
    inf = {"tipo": tipo, "eje": int(eje), "forma": [nx, ny, nz],
           "spacing_mm": spacing.tolist(), "BVTV_voxel": float(BW.mean()),
           "frac_portante_voxel": frac_port, "V_voxel_portante_mm3": V_vox}

    if tipo == "hex8":
        nodos, elems, _ = malla_hex(BW_res, spacing)
        from .escribe import _caras_techo_hex
        caras = _caras_techo_hex(nodos, elems)
        vol = np.full(elems.shape[0], float(np.prod(spacing)))
        sup = capa_superficie(BW_res)[np.nonzero(BW_res)]
        z0, z1 = 0.0, float(nz * spacing[2])
        inf.update(V_malla_mm3=float(vol.sum()), vol_pct_voxel=100.0,
                   descartado_pct=0.0)
    else:
        o = opciones_malla(**opciones)
        objetivo = V_vox if o["corregir_volumen"] else None
        nodos, elems, _s, it = malla_tet10(
            BW_res, spacing, suavizado=o["suavizado"], banda=o["banda"],
            decimado=o["decimado"], progreso=progreso,
            maxvolume=maxvolume_de_tamano(o["tam_max_mm"]),
            minratio=o["minratio"], ejes_planos=(0, 1, 2),
            volumen_objetivo=objetivo)
        if not it["orden_ok"]:
            raise RuntimeError("Orden de nodos TET10 incorrecto: "
                               + it["orden"]["msg"])
        nodos, elems, vol, caras, sup, z0, z1, extra = preparar_tet10(
            nodos, elems, BW.shape, spacing)
        V = extra["V_malla_mm3"]
        inf.update(extra)
        inf.update(
            vol_pct_voxel=100.0 * V / V_vox,
            vol_pct_MC=100.0 * V / it["V_mc"] if it["V_mc"] > 0 else None,
            perdida_pct=100.0 - 100.0 * V / V_vox,
            opciones={k: v for k, v in o.items()},
            maxvolume_mm3=maxvolume_de_tamano(o["tam_max_mm"]),
            tetgen={k: it[k] for k in ("estanca", "bordes_abiertos",
                                       "tri_superficie", "orden",
                                       "correccion_volumen", "tiempo_s")})
        inf["reservas_volumen"] = abs(inf["perdida_pct"]) > PERDIDA_VOLUMEN_MAX_PCT
        inf["reservas_descarte"] = inf["descartado_pct"] > DESCARTE_MAX_PCT

    inf.update(n_nodos=int(nodos.shape[0]), n_elems=int(elems.shape[0]),
               n_gdl=int(3 * nodos.shape[0]),
               n_superficie=int(sup.sum()), n_caras_techo=int(caras.shape[0]))
    return {"tipo": tipo, "nodos": nodos, "elems": elems,
            "caras_techo": caras, "vol_elem": vol, "superficie": sup,
            "A_bruta": A_bruta, "H": float(z1 - z0), "z0": z0, "z1": z1,
            "forma": (nx, ny, nz), "spacing": spacing, "eje": int(eje),
            "huella": huella(nodos, elems), "informe": inf}


def tamano_previsto(BW, spacing, tipo, n=None):
    """Elementos, GDL y memoria previstos del resolvedor directo, sin mallar.

    hex8: EXACTO (voxeles portantes y sus nodos a la resolucion `n`). tet10:
    `tiempos.n_tet10` a partir de los voxeles de superficie, +-4 % en hueso
    trabecular (ver `tiempos`). Cuesta menos de un segundo a 48^3.
    """
    from . import tiempos
    from .resistencia import _solo_portante, capa_superficie
    if n is None:
        n = N_HEX_DEF if tipo == "hex8" else N_TET_DEF
    B0, _sp = _remuestrear(BW, spacing, n)
    B = _solo_portante(B0)
    frac = float(B.sum() / max(B0.sum(), 1))
    if tipo == "hex8":
        i, j, k = np.nonzero(B)
        nn = np.array(B.shape) + 1
        esq = np.array([[a, b, c] for a in (0, 1) for b in (0, 1)
                        for c in (0, 1)])
        ids = np.unique(((i[:, None] + esq[:, 0]) + nn[0] * (
            (j[:, None] + esq[:, 1]) + nn[1] * (k[:, None] + esq[:, 2]))))
        gdl = 3 * int(ids.size)
        n_el = int(i.size)
    else:
        n_el = int(round(tiempos.n_tet10(capa_superficie(B).sum(), B.sum())))
        gdl = int(round(tiempos.gdl_tet10(n_el)))
    return {"tipo": tipo, "n": int(max(B.shape)), "n_elems": n_el,
            "gdl": gdl, "memoria_MB": tiempos.memoria_fem(tipo, gdl),
            "exacto": tipo == "hex8", "frac_portante": frac}


def memoria_equipo_MB():
    """RAM fisica total del equipo en MB (None si no se puede leer)."""
    try:
        import ctypes

        class MEMORYSTATUSEX(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong),
                        ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong),
                        ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong),
                        ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong),
                        ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("sullAvailExtendedVirtual", ctypes.c_ulonglong)]
        m = MEMORYSTATUSEX()
        m.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
        return m.ullTotalPhys / 2 ** 20
    except Exception:
        try:
            return (os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
                    / 2 ** 20)
        except (ValueError, OSError, AttributeError):
            return None


class ErrorMalla(RuntimeError):
    """El mallado TET10 fallo (tetgen abortado o proceso caido)."""


def _mallar_hijo(BW, spacing, tipo, eje, opciones):
    return mallar(BW, spacing, tipo, eje=eje, **opciones)


def mallar_aislado(BW, spacing, tipo="tet10", eje=2, **opciones):
    """`mallar` en un PROCESO HIJO: si tetgen se cae, no se lleva la app.

    tetgen es codigo C sin red: sobre el VOI proximal de H4 remuestreado a
    20^3 (puntales de ~1 voxel, que el suavizado deja casi degenerados) mato
    el interprete con una violacion de segmento. En la GUI eso cerraba la
    ventana entera, con horas de sesion. Aqui el hijo muere solo y se lanza
    `ErrorMalla`, que la etapa registra como fallo y la cadena sigue.
    El hijo cuesta unos segundos de arranque (importar numpy y pyvista);
    con hex8 no hay tetgen y se malla en el mismo proceso.
    """
    if tipo == "hex8":
        return mallar(BW, spacing, tipo, eje=eje, **opciones)
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    from concurrent.futures.process import BrokenProcessPool
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=1, mp_context=ctx) as ex:
        fut = ex.submit(_mallar_hijo, np.asarray(BW, bool),
                        np.asarray(spacing, float), tipo, int(eje),
                        dict(opciones))
        try:
            return fut.result()
        except BrokenProcessPool as e:
            raise ErrorMalla(
                "El mallado TET10 se cayo (tetgen abortado). Suele pasar "
                "cuando la resolucion deja puntales de uno o dos voxeles: "
                "subir la resolucion de la mascara o usar Ladrillos.") from e




# ---------------------------------------------------------------------------
# Protocolos y analisis
# ---------------------------------------------------------------------------

#: Una entrada por protocolo; agregar uno aqui lo agrega a los dialogos.
#: 'app' es el ensayo por omision de la app (`resistencia`). 'tapia2026'
#: reproduce los cuatro valores que definen el ensayo de Tapia, Gonzalez,
#: Vidal & Salinas, Biology 2026;15:722 (los mismos que `visor.PAPER_*`, lo
#: comprueba el bloque 27); alli el mallado fue SOLID187 de 0.05 mm en
#: ANSYS, de ahi `tam_elem_mm`. Si se edita cualquiera de los valores
#: publicados, el registro deja de llamarse Tapia (`nombre_protocolo`).
#: 'homogeneizacion' da el tensor elastico: periodico con hex8 (el mismo
#: problema discreto que `elastic.homogeneizar`) y cotas KUBC/SUBC con tet10.
PROTOCOLOS_FEM = {
    "app": {"tipo": "compresion", "etiqueta": "Ensayo de la app",
            "E_s": 20e9, "nu": 0.30, "sigma_app": 1e6, "carga_N": None,
            "apoyo": "deslizante"},
    "tapia2026": {"tipo": "compresion",
                  "etiqueta": "Protocolo de Tapia et al. (2026)",
                  "E_s": 18e9, "nu": 0.30, "sigma_app": None,
                  "carga_N": 100.0, "apoyo": "empotrado",
                  "tam_elem_mm": 0.05},
    "homogeneizacion": {"tipo": "homogeneizacion",
                        "etiqueta": "Homogeneización (tensor elástico)",
                        "E_s": 20e9, "nu": 0.30, "amplitud": 1e-5,
                        "escala_vacio": 1e-6},
}

#: Claves que definen un protocolo publicado: si alguna cambia, el registro
#: se renombra (p. ej. «Tapia (modificado)»).
CLAVES_PROTOCOLO = ("E_s", "nu", "sigma_app", "carga_N", "apoyo", "amplitud",
                    "escala_vacio")

#: Analisis de un protocolo de compresion. 'lineal_plato', 'nl_plato' y
#: 'nl_pistoia' necesitan el lineal (la deformacion equivalente y la carga de
#: Pistoia salen de el): si no se pidio, se corre igual y se declara.
ANALISIS = ("lineal", "lineal_plato", "nl_fuerza", "nl_plato", "nl_pistoia")



def protocolo(clave, **cambios):
    """Copia de un protocolo con `cambios` aplicados y su nombre de registro."""
    base = PROTOCOLOS_FEM[clave]
    p = dict(base)
    p.update(cambios)
    p["clave"] = clave
    p["modificado"] = any(
        k in base and not _iguales(p.get(k), base.get(k))
        for k in CLAVES_PROTOCOLO)
    p["nombre"] = nombre_protocolo(clave, p["modificado"])
    return p


def _iguales(a, b):
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return float(a) == float(b)
    return a == b


def nombre_protocolo(clave, modificado=False):
    base = {"app": "App", "tapia2026": "Tapia",
            "homogeneizacion": "Homogeneización"}.get(clave, clave)
    return f"{base} (modificado)" if modificado else base


def _sigma_ref(prot, A_bruta):
    """Tension aparente del protocolo en Pa (A_bruta en mm^2).

    Con `carga_N` se deduce de la seccion BRUTA, como `ensayo_compresion`
    con `unidad='mm'`: sigma[Pa] = 1e6 F[N] / A[mm^2].
    """
    if prot.get("carga_N") is not None:
        return 1e6 * float(prot["carga_N"]) / float(A_bruta)
    return float(prot["sigma_app"])


def _uz_medio(malla, u):
    """Desplazamiento vertical medio del techo, ponderado por AREA.

    Es la integral de uz sobre las caras cargadas entre su area. Con quad4
    bilineales cada nodo de una cara pesa A/4; con tri6 las esquinas pesan
    CERO y cada nodo intermedio A/3 (integral de las funciones de forma
    cuadraticas). La media simple de los nodos del techo (la definicion de
    la app con voxeles) pesaria de mas las esquinas de las caras.
    """
    from .escribe import area_caras
    caras = malla["caras_techo"]
    A = area_caras(malla["nodos"], caras)
    uz = u[:, 2]
    if caras.shape[1] == 4:
        m = uz[caras].mean(axis=1)
    else:
        m = uz[caras[:, 3:]].mean(axis=1)
    return float((A * m).sum() / A.sum())


def _uz_medio_nodal(malla, u):
    """Media simple de uz en los nodos del techo: la definicion de la app."""
    z = malla["nodos"][:, 2]
    tolz = 1e-9 * max(float(z.max() - z.min()), 1.0)
    return float(u[z >= z.max() - tolz, 2].mean())


def _fuerza(malla, sigma):
    """F = -sum(sigma_zz V_e) / H (N), con sigma en Pa y mm.

    Integral de volumen de la tension, exacta en el problema discreto lineal
    (desplazamiento virtual v = z e_z). Los motores dan ademas la reaccion
    del vector de fuerzas internas; esta integral queda como COMPROBACION de
    equilibrio independiente del motor (`dF_rel`).
    """
    return float(-(sigma[:, 2] * malla["vol_elem"]).sum() / malla["H"]) * 1e-6


#: Nombre anterior, para scripts y sesiones que lo usan.
PROTOCOLOS_FEBIO = PROTOCOLOS_FEM


# ---------------------------------------------------------------------------
# Correcciones de los artefactos de borde (comparativa_motores/correcciones)
# ---------------------------------------------------------------------------

#: Margen del NUCLEO respecto de las seis caras del VOI (mm). Las trabeculas
#: cortadas por las caras laterales quedan descargadas hasta ~0,6 mm de la
#: cara y la traccion del techo altera ~1,25 mm por debajo; midiendo E y el
#: p99 en el nucleo, frente a la configuracion embebida (el mismo hueso
#: rodeado de hueso), el error de E bajo del -43/-55 % al -5/-10 % con
#: traccion y al +1/+2 % con plato rigido (hex8, 32^3 a 64^3).
NUCLEO_MARGEN_MM = 0.625
#: Espesor de cada franja del extensometro virtual del nucleo (mm).
FRANJA_EXTENSOMETRO_MM = 0.25


def _centroides_uz(nodos, elems, u):
    """Centroide de cada elemento y u_z en el. TET10: las funciones de forma
    cuadraticas valen -1/8 en las esquinas y 1/4 en las aristas."""
    nodos, elems = np.asarray(nodos, float), np.asarray(elems, np.int64)
    if elems.shape[1] == 8:
        return nodos[elems].mean(1), u[elems, 2].mean(1)
    uz = (-0.125 * u[elems[:, :4], 2].sum(1)
          + 0.25 * u[elems[:, 4:], 2].sum(1))
    return nodos[elems[:, :4]].mean(1), uz


def _volumen_recortado(P, lo, hi):
    """Volumen exacto del tetraedro de esquinas P (4, 3) dentro de la caja.

    Se recorta por los seis semiespacios: el poliedro convexo recortado es la
    envolvente de los vertices que quedan dentro y de los cortes de cada
    segmento entre dos vertices con el plano (para un convexo, esos puntos
    estan en el recortado y contienen a sus vertices).
    """
    from scipy.spatial import ConvexHull, QhullError
    V = np.asarray(P, float)
    for e in range(3):
        for lim, signo in ((lo[e], 1.0), (hi[e], -1.0)):
            d = signo * (V[:, e] - lim)
            dentro = d >= 0
            if dentro.all():
                continue
            if not dentro.any():
                return 0.0
            a, b = np.nonzero(dentro)[0], np.nonzero(~dentro)[0]
            t = d[a][:, None] / (d[a][:, None] - d[b][None, :])
            cortes = (V[a][:, None, :] + t[..., None]
                      * (V[b][None, :, :] - V[a][:, None, :])).reshape(-1, 3)
            V = np.vstack([V[dentro], cortes])
    try:
        return float(ConvexHull(V).volume)
    except (QhullError, ValueError):
        return 0.0                       # degenerado: volumen nulo


def _fraccion_en_caja(nodos, elems, lo, hi, dentro):
    """Fraccion del volumen de cada elemento dentro de la caja [lo, hi].

    hex8 con la caja ajustada a la rejilla: 0 o 1 (decide el centroide).
    TET10 (aristas rectas): 1 o 0 si sus cuatro esquinas estan dentro, o
    fuera por un mismo plano; en los que la caja corta, el volumen recortado
    exacto (`_volumen_recortado`) entre el del tetraedro.
    """
    frac = dentro.astype(float)
    if elems.shape[1] == 8:
        return frac
    P = nodos[elems[:, :4]]                                   # (M, 4, 3)
    den = np.all((P >= lo) & (P <= hi), axis=(1, 2))
    fuera = np.zeros(len(P), bool)
    for e in range(3):
        fuera |= np.all(P[:, :, e] <= lo[e], axis=1)
        fuera |= np.all(P[:, :, e] >= hi[e], axis=1)
    frac[den], frac[fuera] = 1.0, 0.0
    vol = np.abs(np.einsum("ij,ij->i", P[:, 1] - P[:, 0],
                           np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]))) / 6
    for m in np.nonzero(~den & ~fuera)[0]:
        frac[m] = min(1.0, _volumen_recortado(P[m], lo, hi) / vol[m])
    return frac


def magnitudes_nucleo(malla, u, sigma, sigma_ref=None,
                      margen=NUCLEO_MARGEN_MM, franja=FRANJA_EXTENSOMETRO_MM):
    """E aparente y p99 de von Mises en el NUCLEO del VOI.

    El nucleo es la caja a `margen` mm o mas de las seis caras del cubo. Se
    mide dentro del ensayo del VOI completo (no se recorta nada):
      sigma_n  = -sum(sigma_zz,e V_e) / V_caja, elementos con centroide en
                 la caja (la tension media de la caja, poros incluidos);
      eps_n    = extensometro virtual: diferencia de u_z medio (ponderado por
                 volumen, en el centroide) entre la franja superior y la
                 inferior de la caja, de `franja` mm cada una, dividida por
                 la distancia entre sus z medios;
      E_app    = sigma_n / eps_n;
      p99      = p99 de von Mises en la capa superficial de la caja,
                 ponderado por volumen. Con `sigma_ref` (Pa, la tension
                 aparente del protocolo) se expresa a esa tension aplicada AL
                 NUCLEO: p99 * sigma_ref / sigma_n, la misma normalizacion
                 con la que se valido.
    `sigma` en Pa (la de `tension_elemental`). Si el VOI es demasiado
    pequeno para el margen, devuelve {"ok": False}.
    """
    from .resistencia import percentil_ponderado
    nodos, elems = malla["nodos"], np.asarray(malla["elems"], np.int64)
    sp = np.asarray(malla["spacing"], float)
    forma = np.asarray(malla["forma"], float)
    # La caja se ajusta a la rejilla de voxeles de la mascara: con hex8 los
    # elementos con centroide dentro la cubren EXACTAMENTE (sin el ajuste, un
    # margen que no es multiplo del voxel sesga la tension media por el
    # cociente entre volumenes); con TET10 la seleccion por centroide es
    # insesgada en promedio.
    k = np.maximum(np.rint(margen / sp), 1.0)
    lo = (-0.5 * sp if malla["tipo"] == "tet10" else np.zeros(3)) + k * sp
    hi = lo + (forma - 2 * k) * sp
    margen = float((k * sp).min())
    if np.any(hi - lo < 4 * franja):
        return {"ok": False, "msg": "VOI demasiado pequeno para el nucleo",
                "margen_mm": margen}
    c, uz = _centroides_uz(nodos, elems, u)
    w = np.asarray(malla["vol_elem"], float)
    dentro = np.all((c >= lo) & (c <= hi), axis=1)
    V = float(np.prod(hi - lo))
    frac = _fraccion_en_caja(nodos, elems, lo, hi, dentro)
    s_n = float(-(sigma[:, 2] * w * frac).sum() / V)
    z = c[:, 2]
    arr = dentro & (z >= hi[2] - franja)
    aba = dentro & (z <= lo[2] + franja)
    if not arr.any() or not aba.any() or s_n <= 0:
        return {"ok": False, "msg": "sin hueso en las franjas del nucleo",
                "margen_mm": margen}
    dz = np.average(z[arr], weights=w[arr]) - np.average(z[aba], weights=w[aba])
    du = np.average(uz[arr], weights=w[arr]) - np.average(uz[aba],
                                                          weights=w[aba])
    eps = abs(du / dz)
    sup = dentro & np.asarray(malla["superficie"], bool)
    vm = von_mises(sigma)
    out = {"ok": True, "E_app": s_n / eps, "eps_app": eps, "sigma_Pa": s_n,
           "margen_mm": margen, "franja_mm": franja,
           "frac_volumen": float(w[dentro].sum() / w.sum()),
           "n_superficie": int(sup.sum())}
    if sup.any():
        p99 = float(percentil_ponderado(vm[sup], w[sup], 99))
        out["vm_p99_superficie"] = (p99 if sigma_ref is None
                                    else p99 * float(sigma_ref) / s_n)
    return out



def estadisticos(malla, sigma, prot, sigma_ref, F_total_Pa_mm2):
    """p99 de superficie y Pistoia de un campo de tensiones, ponderados por volumen.

    Se usan las MISMAS funciones que la app (`criterio_pistoia`,
    `estadisticos_vm`), a las que se pasa `vol_solido`: con hexaedros de un
    voxel los pesos son iguales y el resultado es bit a bit el de siempre.
    Con TET10 la tension de cada elemento es la media de sus puntos de Gauss.
    """
    from .resistencia import criterio_pistoia
    res = {"ok": True, "sigma_app": float(sigma_ref),
           "F_total": float(F_total_Pa_mm2),
           "eps_eff_solido": eps_eff(sigma, prot["E_s"], prot["nu"]),
           "vm_solido": von_mises(sigma),
           "superficie_solido": np.asarray(malla["superficie"], bool),
           "vol_solido": np.asarray(malla["vol_elem"], float)}
    return criterio_pistoia(res)




# ---------------------------------------------------------------------------
# Ensayo con un motor
# ---------------------------------------------------------------------------

def tension_elemental(nodos, elems, u, E, nu):
    """Tension (Voigt, Pa) media por elemento del campo lineal `u` (mm).

    hex8: gradiente en el centro del voxel, la convencion de la app
    (`resistencia._b_centro`); con elementos rectangulares es la media de la
    deformacion sobre el elemento. tet10 de aristas rectas: la deformacion
    es lineal en el elemento y su valor en el centroide es su media de
    volumen (lo mismo que escribia FEBio: la media de sus puntos de Gauss).
    """
    nodos = np.asarray(nodos, float)
    elems = np.asarray(elems, np.int64)
    if elems.shape[1] == 8:
        P = nodos[elems]
        h = P[:, 6] - P[:, 0]
        s = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                      [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
        G = 0.125 * s[None] * (2.0 / h)[:, None, :]
    else:
        L = np.full(4, 0.25)
        dNdL = np.zeros((10, 4))
        for i in range(4):
            dNdL[i, i] = 4 * L[i] - 1
        for k, (i, j) in enumerate(((0, 1), (1, 2), (0, 2), (0, 3), (1, 3),
                                    (2, 3))):
            dNdL[4 + k, i] = 4 * L[j]
            dNdL[4 + k, j] = 4 * L[i]
        dLdxi = np.array([[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         float)
        P = nodos[elems[:, :4]]
        J = np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0],
                      P[:, 3] - P[:, 0]], axis=2)
        G = np.einsum("ak,mkj->maj", dNdL @ dLdxi, np.linalg.inv(J))
    H = np.einsum("mai,maj->mij", np.asarray(u, float)[elems], G)
    eps = np.stack([H[:, 0, 0], H[:, 1, 1], H[:, 2, 2],
                    H[:, 1, 2] + H[:, 2, 1], H[:, 0, 2] + H[:, 2, 0],
                    H[:, 0, 1] + H[:, 1, 0]], axis=1)
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[0, 0] = D[1, 1] = D[2, 2] = lam + 2 * mu
    D[3, 3] = D[4, 4] = D[5, 5] = mu
    return eps @ D.T


def _problema(malla, prot, analisis, control, cargas, material, solver,
              BW=None):
    return motores.problema_de_malla(
        malla, E=prot["E_s"] / 1e6, nu=prot["nu"], apoyo=prot["apoyo"],
        analisis=analisis, control=control, cargas=cargas, material=material,
        solver=solver, BW=BW)


def _corregido(reg):
    """Valores CORREGIDOS de los artefactos de borde, con su metodo.

    Validados frente a la configuracion embebida en
    `comparativa_motores/correcciones/` (hex8, 32^3 a 64^3):
      E_app   nucleo con plato rigido (error +0,8 a +1,5 %); sin el plato,
              nucleo con traccion (-5,5 a -10 %).
      p99     VOI completo con plato (+0,5 a +9,5 %); sin el plato, nucleo
              con traccion (-1,9 a +8,7 %).
    La linea base (traccion, VOI completo) daba -43 a -55 % en E_app y
    +33 a +53 % en el p99.
    """
    lin = reg.get("lineal") or {}
    lp = reg.get("lineal_plato") or {}
    nl, npl = lin.get("nucleo") or {}, lp.get("nucleo") or {}
    out = {}
    if npl.get("ok"):
        out["E_app"], out["metodo_E"] = npl["E_app"], "plato_nucleo"
    elif nl.get("ok"):
        out["E_app"], out["metodo_E"] = nl["E_app"], "traccion_nucleo"
    p99 = (lp.get("pistoia") or {}).get("vm_p99_superficie")
    if p99 is not None:
        out["vm_p99_superficie"], out["metodo_p99"] = p99, "plato_voi"
    elif nl.get("vm_p99_superficie") is not None:
        out["vm_p99_superficie"] = nl["vm_p99_superficie"]
        out["metodo_p99"] = "traccion_nucleo"
    if lin.get("E_app") and "E_app" in out:
        out["cociente_E_lineal"] = out["E_app"] / lin["E_app"]
    return out


def ensayo(malla, prot, carpeta=None, analisis=("lineal",), motor=MOTOR_DEF,
           hilos=None, cancelar=None, material="svk", pasos=1, progreso=None,
           solver="auto", aislado=True, BW=None, prefijo="ensayo", **_viejos):
    """Ensayo de compresion de un protocolo sobre `malla`, con `motor`.

    Corre los `analisis` pedidos (ver `ANALISIS`) y devuelve un registro con
    procedencia. Un analisis que falla se registra y NO detiene los demas;
    solo `Cancelado` sale hacia arriba.

      lineal      el problema lineal a la carga del protocolo: E_app, fuerza
                  (comprobacion de equilibrio), p99 de von Mises en la capa
                  superficial y Pistoia; ademas E_app y p99 en el NUCLEO
                  (`magnitudes_nucleo`), lejos de las caras cortadas.
      lineal_plato el mismo lineal con plato rigido sin friccion en el techo
                  (corrige el artefacto del techo cargado): E_app, p99 y
                  Pistoia a la tension del protocolo, y su nucleo. Con el
                  calcula el registro `corregido` (ver `_corregido`).
      nl_fuerza   no lineal a la carga del protocolo, traccion muerta, en
                  `pasos` incrementos.
      nl_plato    plato rigido sin friccion: lineal del plato y no lineal a la
                  deformacion equivalente a la carga del protocolo (la del
                  lineal con fuerza), como `comparativa_febio/plato.py`.
      nl_pistoia  no lineal a la carga de fallo de Pistoia del lineal.

    E_app se define con uz medio del techo PONDERADO POR AREA (`_uz_medio`);
    con hex8 se da ademas `E_app_nodal` (media simple de nodos, la de la app).
    `BW` (la mascara de la malla hex8) solo lo necesita el motor 'app'.
    `aislado`: cada resolucion en un proceso hijo cancelable (la GUI).
    """
    analisis = [a for a in ANALISIS if a in set(analisis)]
    if not analisis:
        raise ValueError("No se pidio ningun analisis.")
    necesita_lineal = "lineal" in analisis or any(
        a in analisis for a in ("lineal_plato", "nl_plato", "nl_pistoia"))
    sref = _sigma_ref(prot, malla["A_bruta"])
    F_ref_N = sref * malla["A_bruta"] * 1e-6
    H = malla["H"]
    reg = {"protocolo": prot.get("clave"), "nombre": prot.get("nombre"),
           "modificado": bool(prot.get("modificado")),
           "parametros": {k: prot.get(k) for k in CLAVES_PROTOCOLO
                          if k in prot},
           "malla": malla["tipo"], "huella_malla": malla["huella"],
           "informe_malla": malla["informe"], "eje": malla["eje"],
           "sigma_ref_Pa": sref, "F_ref_N": F_ref_N,
           "material_no_lineal": material, "pasos": int(pasos),
           "motor": {"clave": motor, "nombre": motores.ETIQUETAS[motor],
                     "version": motores.version(motor), "solver": None},
           "corridas": [], "fallos": [], "analisis_pedidos": list(analisis),
           "carpeta": str(carpeta) if carpeta else None}
    n_tot = (1 if necesita_lineal else 0) + analisis.count("nl_fuerza") \
        + int(("lineal_plato" in analisis or "nl_plato" in analisis)
              and motores.puede(motor, malla["tipo"], "lineal",
                                control="plato")) \
        + analisis.count("nl_plato") + analisis.count("nl_pistoia")
    hechas = [0]

    def run(nombre, an, control, cargas):
        """Una resolucion; devuelve la salida del motor o None si fallo."""
        if progreso:
            progreso(hechas[0], n_tot, f"{motores.ETIQUETAS[motor]} · "
                     f"{nombre}")
        p = _problema(malla, prot, an, control, cargas, material, solver,
                      BW=BW if motor == "app" else None)
        try:
            out = (motores.resolver_aislado(p, motor, cancelar=cancelar,
                                            hilos=hilos) if aislado
                   else motores.resolver(p, motor))
        except (ErrorMotor, NoDisponible) as e:
            reg["fallos"].append({"corrida": nombre, "msg": str(e)})
            hechas[0] += 1
            return None
        hechas[0] += 1
        mt = out["meta"]
        reg["motor"]["solver"] = reg["motor"]["solver"] or mt.get("solver")
        reg["corridas"].append({
            "corrida": nombre, "tiempo_s": mt.get("tiempo_total_s"),
            "tiempos": mt.get("tiempos"), "memoria_MB": mt.get("rss_pico_MB"),
            "solver": mt.get("solver"), "residuo_rel": mt.get("residuo_rel"),
            "iteraciones": mt.get("iteraciones"),
            "iteraciones_newton": mt.get("iteraciones_newton"),
            "n_gdl": mt.get("n_gdl"), "control": control})
        return out

    lin = None
    if necesita_lineal:
        out = run("lineal", "lineal", "fuerza", [sref / 1e6])
        if out is not None:
            u = out["u"]
            s = tension_elemental(malla["nodos"], malla["elems"], u,
                                  prot["E_s"], prot["nu"])
            uz = _uz_medio(malla, u)
            eps_app = abs(uz) / H
            F = _fuerza(malla, s)
            lin = {"E_app": sref / eps_app, "eps_app": eps_app,
                   "F_N": F, "dF_rel": abs(F - F_ref_N) / F_ref_N,
                   "F_reac_N": float(out["F_reac"][-1]),
                   "u_max_mm": float(np.abs(u).max())}
            if malla["tipo"] == "hex8":
                lin["E_app_nodal"] = sref / (abs(_uz_medio_nodal(malla, u))
                                             / H)
            p = estadisticos(malla, s, prot, sref, sref * malla["A_bruta"])
            lin["pistoia"] = {k: v for k, v in p.items()
                              if isinstance(v, (int, float, bool, str))}
            lin["nucleo"] = magnitudes_nucleo(malla, u, s, sigma_ref=sref)
            reg["_campos_lineal"] = {"u": u, "sigma": s}
        reg["lineal"] = lin

    plato_lin = None
    puede_plato = motores.puede(motor, malla["tipo"], "lineal",
                                control="plato")
    if lin and "lineal_plato" in analisis and not puede_plato:
        reg["lineal_plato"] = {"no_disponible": True, "msg": (
            f"{motores.ETIQUETAS[motor]} no resuelve el plato rigido; el "
            "valor corregido usa el nucleo con traccion")}
    if lin and puede_plato and ("lineal_plato" in analisis
                                or "nl_plato" in analisis):
        eps_p = sref / lin["E_app"]
        plato_lin = run("plato_lineal", "lineal", "plato", [eps_p])
        if plato_lin is not None and "lineal_plato" in analisis:
            u = plato_lin["u"]
            F = float(plato_lin["F_reac"][-1])
            s_p = F * 1e6 / malla["A_bruta"]
            # Lineal: el campo a la tension del protocolo es el del plato
            # escalado por sref / s_p.
            s = tension_elemental(malla["nodos"], malla["elems"], u,
                                  prot["E_s"], prot["nu"]) * (sref / s_p)
            p = estadisticos(malla, s, prot, sref, sref * malla["A_bruta"])
            reg["lineal_plato"] = {
                "E_app": s_p / eps_p, "eps_plato": eps_p, "F_N": F,
                "cociente_plato_fuerza": (s_p / eps_p) / lin["E_app"],
                "pistoia": {k: v for k, v in p.items()
                            if isinstance(v, (int, float, bool, str))},
                "nucleo": magnitudes_nucleo(malla, u * (sref / s_p), s,
                                            sigma_ref=sref)}
            reg["_campos_lineal_plato"] = {"u": u * (sref / s_p), "sigma": s}
    if lin:
        reg["corregido"] = _corregido(reg)

    def escalones(total):
        return [total * (k + 1) / max(1, int(pasos)) for k in
                range(max(1, int(pasos)))]

    if "nl_fuerza" in analisis:
        out = run("nl_fuerza", "nl", "fuerza", escalones(sref / 1e6))
        if out is not None:
            e = abs(_uz_medio(malla, out["u"])) / H
            reg["nl_fuerza"] = {"E_app": sref / e, "eps_app": e,
                                "F_N": float(out["F_reac"][-1])}
            if lin:
                reg["nl_fuerza"]["dE_rel"] = reg["nl_fuerza"]["E_app"] \
                    / lin["E_app"] - 1.0
            reg["_campos_nl_fuerza"] = {"u": out["u"]}

    if "nl_plato" in analisis and lin:
        eps_p = sref / lin["E_app"]
        d = {"eps_plato": eps_p}
        out = plato_lin
        if out is not None:
            d["E_app_lineal"] = out["F_reac"][-1] * 1e6 \
                / malla["A_bruta"] / eps_p
            d["cociente_plato_fuerza_lineal"] = d["E_app_lineal"] \
                / lin["E_app"]
        out = run("nl_plato", "nl", "plato", escalones(eps_p))
        if out is not None:
            d["F_N"] = float(out["F_reac"][-1])
            d["E_app"] = d["F_N"] * 1e6 / malla["A_bruta"] / eps_p
            if "E_app_lineal" in d:
                d["dE_rel"] = d["E_app"] / d["E_app_lineal"] - 1.0
        reg["nl_plato"] = d

    if "nl_pistoia" in analisis and lin and lin["pistoia"].get("ok"):
        k = lin["pistoia"]["factor"]
        d = {"factor": k, "sigma_Pa": k * sref}
        out = run("nl_pistoia", "nl", "fuerza", escalones(k * sref / 1e6))
        if out is not None:
            e = abs(_uz_medio(malla, out["u"])) / H
            d.update(E_app=k * sref / e, eps_app=e,
                     dE_rel=(k * sref / e) / lin["E_app"] - 1.0,
                     F_N=float(out["F_reac"][-1]))
        reg["nl_pistoia"] = d

    if progreso:
        progreso(n_tot, n_tot, "listo")
    reg["ok"] = bool(reg["corridas"]) and not reg["fallos"]
    return reg


# ---------------------------------------------------------------------------
# Un analisis completo: el UNICO camino de calculo de las dos puertas
# ---------------------------------------------------------------------------

#: Resolucion por omision de la malla de ladrillos: la de los ensayos que se
#: reportan (Estudio_Convergencia: 40-48; a 32^3 la rigidez del VOI proximal pierde un
#: ~10 %).
N_HEX_DEF = 40
#: Resolucion por omision de la mascara de la que sale la malla suave. Algo
#: mas fina que la de ladrillos (el suavizado de Taubin se come puntales de
#: uno o dos voxeles, y a 40^3 Tb.Th/h ronda 2), pero no 64: el VOI proximal
#: de H4 a 64^3 son 350 000 TET10 y 1.9 millones de GDL (medido), que con
#: Pardiso no caben en un equipo de 16 GB. A 48^3: 179 000 y 1.0 millones.
N_TET_DEF = 48


def _remuestrear(BW, spacing, n):
    from .elastic import remuestrear_bw
    if n is None or int(n) == max(BW.shape):
        return np.asarray(BW, bool), np.asarray(spacing, float)
    return remuestrear_bw(BW, spacing, int(n))


def ensayo_app(BW, spacing, prot, eje=2):
    """El ensayo LINEAL de la app con los valores del protocolo.

    Es la columna «app (hex8, lineal)» de la tabla. Mismo `resistencia` que
    el panel de Analisis mecanico; nada propio.
    """
    from .resistencia import criterio_pistoia, ensayo_compresion_eje
    kw = {"E_s": prot["E_s"], "nu_s": prot["nu"], "apoyo": prot["apoyo"]}
    if prot.get("carga_N") is not None:
        kw["carga_N"] = prot["carga_N"]
    else:
        kw["sigma0"] = prot["sigma_app"]
    t0 = time.perf_counter()
    r = ensayo_compresion_eje(BW, spacing, eje=eje, **kw)
    out = {"ok": bool(r.get("ok")), "msg": r.get("msg", ""),
           "tiempo_s": time.perf_counter() - t0,
           "solver": r.get("solver"), "residuo_rel": r.get("residuo_rel")}
    if r.get("ok"):
        p = criterio_pistoia(r)
        out.update(E_app=r["E_app"], n_elem=r["n_elem"],
                   frac_portante=r["frac_portante"],
                   pistoia={k: v for k, v in p.items()
                            if isinstance(v, (int, float, bool, str))})
    return out


def analizar(BW, spacing, prot, malla="hex8", analisis=("lineal",), eje=2,
             carpeta=".", n=None, motores_fem=None, hilos=None, cancelar=None,
             material="svk", pasos=1, progreso=None, comparar_app=True,
             opciones_malla=None, etiqueta="estructura", conv_malla=False,
             solver="auto", aislado=True, motor=None, **_viejos):
    """Una estructura, un protocolo de compresion, una malla, un eje, y uno
    o VARIOS motores.

    Es lo que corren la GUI («Informe FEM (Auto)…» y el ensayo del panel) y
    la CLI: remuestrea la mascara a `n` (vecino mas proximo, como la app;
    None = la resolucion de siempre de cada malla), malla UNA vez y resuelve
    la misma malla con cada motor de `motores_fem` (`ensayo`). Si
    `comparar_app`, corre ademas el ensayo lineal de la app sobre la MISMA
    mascara remuestreada a la resolucion de ladrillos (la columna «app»).

    `conv_malla` (solo TET10): repite el lineal con la mitad del tamano de
    elemento y sin decimar, con el primer motor, y guarda
    `convergencia_malla` (dE_rel, dp99_rel).

    Devuelve UNA LISTA de registros, uno por motor, que comparten `app`,
    `estructura`, `n`, `tiempo_malla_s` y la malla (`_malla`). Los que no
    resuelven esa combinacion (p. ej. la app con TET10) se registran con el
    motivo en `fallos` y `no_disponible`.
    """
    if motor is not None and motores_fem is None:
        motores_fem = [motor]
    motores_fem = list(motores_fem or [MOTOR_DEF])
    if n is None:
        n = N_HEX_DEF if malla == "hex8" else N_TET_DEF
    t0 = time.perf_counter()
    B, sp = _remuestrear(BW, spacing, n)
    op = dict(opciones_malla or {})
    if malla == "tet10" and op.get("tam_max_mm") is None and \
            prot.get("tam_elem_mm"):
        op["tam_max_mm"] = prot["tam_elem_mm"]
    m = mallar_aislado(B, sp, malla, eje=eje,
                       **(op if malla == "tet10" else {}))
    t_malla = time.perf_counter() - t0
    BW_hex = None
    if malla == "hex8":
        from .resistencia import _PERM, _solo_portante
        BW_hex = _solo_portante(np.transpose(B, _PERM[int(eje)]))
    app = None
    if comparar_app and "lineal" in analisis:
        Ba, spa = _remuestrear(BW, spacing, N_HEX_DEF if malla == "tet10"
                               else n)
        app = ensayo_app(Ba, spa, prot, eje=eje)
        app["n"] = int(max(Ba.shape))
    regs = []
    for mot in motores_fem:
        nombre = f"{etiqueta}_{prot.get('clave')}_{malla}_{'XYZ'[eje]}_{mot}"
        reg = ensayo(m, prot, Path(carpeta) / nombre, analisis=[
            a for a in analisis if a == "lineal"
            or motores.puede(mot, malla, "nl", material,
                             "plato" if a == "nl_plato" else "fuerza")],
            motor=mot, hilos=hilos, cancelar=cancelar, material=material,
            pasos=pasos, progreso=progreso, solver=solver, aislado=aislado,
            BW=BW_hex) if motores.puede(mot, malla) else None
        if reg is None:
            reg = {"protocolo": prot.get("clave"), "nombre": prot.get("nombre"),
                   "malla": malla, "eje": int(eje), "ok": False,
                   "motor": {"clave": mot, "nombre": motores.ETIQUETAS[mot],
                             "version": motores.version(mot)},
                   "fallos": [{"corrida": "todas", "msg": (
                       f"{motores.ETIQUETAS[mot]} no resuelve {malla}")}],
                   "no_disponible": True, "corridas": [],
                   "analisis_pedidos": list(analisis),
                   "informe_malla": m["informe"], "huella_malla": m["huella"]}
        else:
            omitidos = [a for a in analisis if a not in reg["analisis_pedidos"]]
            for a in omitidos:
                reg["fallos"].append({"corrida": a, "msg": (
                    f"{motores.ETIQUETAS[mot]} no resuelve {a} con "
                    f"{material}")})
            reg["analisis_pedidos"] = list(analisis)
            reg["ok"] = bool(reg["corridas"]) and not reg["fallos"]
        reg.update(estructura=etiqueta, n=int(max(B.shape)),
                   eje_nombre="XYZ"[eje], tiempo_malla_s=t_malla,
                   opciones_malla=op, _malla=m)
        if app is not None:
            reg["app"] = app
        regs.append(reg)

    if conv_malla and malla == "tet10" and regs and regs[0].get("lineal"):
        # Segunda malla: la MITAD del tamano de elemento (de la mediana del
        # volumen de la primera si no se fijo) y sin decimar la superficie,
        # sobre la MISMA mascara. Solo el lineal, con el primer motor.
        V_med = float(np.median(m["vol_elem"]))
        tam1 = op.get("tam_max_mm") or (6.0 * np.sqrt(2.0) * V_med) ** (1 / 3)
        op2 = dict(op, tam_max_mm=0.5 * float(tam1), decimado=0)
        try:
            m2 = mallar_aislado(B, sp, "tet10", eje=eje, **op2)
            r2 = ensayo(m2, prot, analisis=("lineal",), motor=motores_fem[0],
                        hilos=hilos, cancelar=cancelar, solver=solver,
                        aislado=aislado)
            l1, l2 = regs[0]["lineal"], r2.get("lineal")
            d = {"tam_mm": [float(tam1), 0.5 * float(tam1)],
                 "n_elems": [int(m["elems"].shape[0]),
                             int(m2["elems"].shape[0])],
                 "huella_fina": m2["huella"], "fallos": r2["fallos"],
                 "motor": motores_fem[0]}
            if l2:
                d["dE_rel"] = l2["E_app"] / l1["E_app"] - 1.0
                p1 = l1["pistoia"].get("vm_p99_superficie")
                p2 = l2["pistoia"].get("vm_p99_superficie")
                if p1 and p2:
                    d["dp99_rel"] = p2 / p1 - 1.0
            regs[0]["convergencia_malla"] = d
        except ErrorMalla as e:
            regs[0]["convergencia_malla"] = {
                "fallos": [{"corrida": "malla_fina", "msg": str(e)}]}
    return regs


def tabla_motores(registros, referencia=None):
    """Filas de la TABLA COMPARATIVA ENTRE MOTORES.

    Agrupa por (estructura, protocolo, malla, eje, n) y compara cada motor
    con el de `referencia` (por omision, el primero del grupo que resolvio
    el lineal): dE_app, dp99 de superficie y d sigma_fallo relativos. Sobre
    la misma malla las diferencias miden solo la implementacion; lo esperable
    es ~1e-8 o menos (la tolerancia de los resolvedores iterativos).
    """
    grupos = {}
    for r in registros:
        if r.get("tipo") == "homogeneizacion":
            continue
        k = (r.get("estructura"), r.get("nombre"), r.get("malla"),
             r.get("eje_nombre"), r.get("n"))
        grupos.setdefault(k, []).append(r)
    filas = []
    for k, rs in grupos.items():
        ref = None
        for r in rs:
            if r.get("lineal") and (referencia is None or
                                    (r.get("motor") or {}).get("clave")
                                    == referencia):
                ref = r
                break
        for r in rs:
            mt = r.get("motor") or {}
            lin = r.get("lineal") or {}
            p = lin.get("pistoia") or {}
            f = {"estructura": k[0], "protocolo": k[1], "malla": k[2],
                 "eje": k[3], "n": k[4], "motor": mt.get("nombre"),
                 "clave_motor": mt.get("clave"), "version": mt.get("version"),
                 "solver": mt.get("solver"), "ok": r.get("ok"),
                 "no_disponible": bool(r.get("no_disponible")),
                 "fallos": len(r.get("fallos", []))}
            if lin:
                f.update(E_app_MPa=lin["E_app"] / 1e6,
                         vm_p99_sup_MPa=p.get("vm_p99_superficie",
                                              np.nan) / 1e6,
                         sigma_fallo_MPa=p.get("sigma_fallo", np.nan) / 1e6,
                         dF_rel=lin.get("dF_rel"))
                corr = [c for c in r.get("corridas", [])
                        if c["corrida"] == "lineal"]
                if corr:
                    f["tiempo_s"] = corr[0].get("tiempo_s")
                    f["memoria_MB"] = corr[0].get("memoria_MB")
                    f["residuo_rel"] = corr[0].get("residuo_rel")
                if ref is not None and ref is not r:
                    lr = ref["lineal"]
                    pr = lr.get("pistoia") or {}
                    f["ref"] = (ref.get("motor") or {}).get("nombre")
                    f["dE_rel"] = lin["E_app"] / lr["E_app"] - 1.0
                    if p.get("vm_p99_superficie") and \
                            pr.get("vm_p99_superficie"):
                        f["dp99_rel"] = p["vm_p99_superficie"] \
                            / pr["vm_p99_superficie"] - 1.0
                    if p.get("sigma_fallo") and pr.get("sigma_fallo"):
                        f["dfallo_rel"] = p["sigma_fallo"] \
                            / pr["sigma_fallo"] - 1.0
                    cu, cr = r.get("_campos_lineal"), ref.get("_campos_lineal")
                    if cu is not None and cr is not None:
                        f["du_rel"] = float(np.abs(cu["u"] - cr["u"]).max()
                                            / np.abs(cr["u"]).max())
            cor = r.get("corregido") or {}
            if "E_app" in cor:
                f["E_app_corregido_MPa"] = cor["E_app"] / 1e6
            for kk in ("nl_fuerza", "nl_plato", "nl_pistoia"):
                d = r.get(kk) or {}
                if "dE_rel" in d:
                    f[f"dE_{kk}"] = d["dE_rel"]
            filas.append(f)
    return filas


def registro_json(reg):
    """Copia serializable del registro: sin campos ni mallas en memoria."""
    def limpio(v):
        if isinstance(v, dict):
            return {k: limpio(x) for k, x in v.items()
                    if not str(k).startswith("_")}
        if isinstance(v, (list, tuple)):
            return [limpio(x) for x in v]
        if isinstance(v, np.generic):
            return v.item()
        if isinstance(v, np.ndarray):
            return v.tolist()
        return v
    return limpio(reg)


def fila_tabla(reg):
    """Numeros de la tabla comparativa de un registro (MPa y adimensionales).

    `dE_app_implementacion` (solo hex8): el motor frente a la app con la
    DEFINICION de la app (media simple de nodos del techo); debe ser ~0.
    Las demas diferencias frente a la app mezclan malla y definicion y se
    rotulan asi en la GUI.
    """
    f = {"estructura": reg.get("estructura"), "protocolo": reg.get("nombre"),
         "motor": (reg.get("motor") or {}).get("nombre"),
         "malla": reg.get("malla"), "eje": reg.get("eje_nombre"),
         "n": reg.get("n"), "ok": reg.get("ok"),
         "fallos": len(reg.get("fallos", []))}
    lin = reg.get("lineal") or {}
    app = reg.get("app") or {}
    inf = reg.get("informe_malla") or {}
    f["BVTV_malla"] = inf.get("BVTV_malla", inf.get("BVTV_voxel"))
    f["perdida_volumen_pct"] = inf.get("perdida_pct")
    f["descartado_pct"] = inf.get("descartado_pct")
    if lin:
        p = lin.get("pistoia", {})
        f.update(E_app_MPa=lin["E_app"] / 1e6,
                 vm_p99_sup_MPa=(p.get("vm_p99_superficie", np.nan) / 1e6),
                 sigma_fallo_MPa=(p.get("sigma_fallo", np.nan) / 1e6),
                 dF_rel=lin.get("dF_rel"))
        if "E_app_nodal" in lin and app.get("ok"):
            f["dE_app_implementacion"] = lin["E_app_nodal"] / app["E_app"] - 1
    if app.get("ok"):
        pa = app.get("pistoia", {})
        f.update(E_app_app_MPa=app["E_app"] / 1e6,
                 vm_p99_sup_app_MPa=pa.get("vm_p99_superficie", np.nan) / 1e6,
                 sigma_fallo_app_MPa=pa.get("sigma_fallo", np.nan) / 1e6)
    for k in ("nl_fuerza", "nl_plato", "nl_pistoia"):
        d = reg.get(k) or {}
        if "dE_rel" in d:
            f[f"dE_{k}"] = d["dE_rel"]
    if (reg.get("nl_plato") or {}).get("cociente_plato_fuerza_lineal"):
        f["plato_fuerza_lineal"] = reg["nl_plato"]["cociente_plato_fuerza_lineal"]
    if (reg.get("lineal_plato") or {}).get("cociente_plato_fuerza"):
        f["plato_fuerza_lineal"] = reg["lineal_plato"]["cociente_plato_fuerza"]
    cor = reg.get("corregido") or {}
    if "E_app" in cor:
        f["E_app_corregido_MPa"] = cor["E_app"] / 1e6
        f["metodo_E_corregido"] = cor.get("metodo_E")
    if "vm_p99_superficie" in cor:
        f["vm_p99_corregido_MPa"] = cor["vm_p99_superficie"] / 1e6
        f["metodo_p99_corregido"] = cor.get("metodo_p99")
    return f


# ---------------------------------------------------------------------------
# Homogeneizacion
# ---------------------------------------------------------------------------

VOIGT = ("xx", "yy", "zz", "yz", "xz", "xy")


def tensor_de_voigt(e):
    """Voigt [xx yy zz yz xz xy] con distorsiones de INGENIERIA -> tensor."""
    return np.array([[e[0], e[5] / 2, e[4] / 2],
                     [e[5] / 2, e[1], e[3] / 2],
                     [e[4] / 2, e[3] / 2, e[2]]], float)



def _caras_planos(malla):
    """Caras tri6 de borde sobre cada uno de los seis planos del cubo."""
    nodos, elems = malla["nodos"], malla["elems"]
    dueno, cb = caras_borde_tet10(elems)
    cb = _orientar_saliente(nodos, elems, dueno, cb)
    sp, forma = malla["spacing"], malla["forma"]
    out = {}
    for e in range(3):
        a, b = -0.5 * sp[e], (forma[e] - 0.5) * sp[e]
        c = nodos[cb[:, :3], e]
        t = 1e-6 * float(sp[e])
        out[(e, 0)] = (cb[np.all(np.abs(c - a) <= t, axis=1)], a)
        out[(e, 1)] = (cb[np.all(np.abs(c - b) <= t, axis=1)], b)
    return out


def _mayor_componente(malla):
    """La malla TET10 reducida a su mayor componente conexo."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    nodos, elems = malla["nodos"], malla["elems"]
    e4 = elems[:, :4]
    n = nodos.shape[0]
    G = coo_matrix((np.ones(3 * e4.shape[0]),
                    (np.repeat(e4[:, 0], 3), e4[:, 1:].ravel())), shape=(n, n))
    _, lab = connected_components(G, directed=False)
    le = lab[e4[:, 0]]
    vol = malla["vol_elem"]
    tot = np.bincount(le, weights=vol)
    k = int(np.argmax(tot))
    keep = le == k
    m = malla_de_tet10(nodos, elems[keep], malla["forma"], malla["spacing"],
                       informe=malla["informe"])
    m["informe"]["descartado_homog_pct"] = float(
        100.0 * (vol.sum() - vol[keep].sum()) / vol.sum())
    return m


def homogeneizar(BW, spacing, prot, malla="hex8", carpeta=".", n=None,
                 progreso=None, etiqueta="estructura", motores_fem=None,
                 **_viejos):
    """Tensor elastico periodico sobre la malla de ladrillos.

    Con hex8 es `elastic.homogeneizar`, el problema discreto que se valido
    contra FEBio a 1e-6 del maximo de C (`comparativa_febio/porcino/`). Las
    cotas KUBC/SUBC sobre TET10 que daba la ruta de FEBio no estan todavia en
    los motores internos: con tet10 se devuelve el registro con el motivo en
    `fallos` en lugar de un tensor.
    """
    from .elastic import homogeneizar as homog_app
    if n is None:
        n = N_HEX_DEF
    B, sp = _remuestrear(BW, spacing, n)
    reg = {"protocolo": prot.get("clave"), "nombre": prot.get("nombre"),
           "modificado": bool(prot.get("modificado")),
           "tipo": "homogeneizacion", "malla": malla, "estructura": etiqueta,
           "n": int(max(B.shape)),
           "parametros": {k: prot.get(k) for k in CLAVES_PROTOCOLO
                          if k in prot},
           "motor": {"clave": "app", "nombre": motores.ETIQUETAS["app"],
                     "version": motores.version("app")},
           "corridas": [], "fallos": [], "carpeta": str(carpeta)}
    if malla != "hex8":
        reg["fallos"].append({"corrida": "homogeneizacion", "msg": (
            "Las cotas KUBC/SUBC con TET10 no estan disponibles en los "
            "motores internos; usar ladrillos (tensor periodico).")})
        reg["no_disponible"] = True
        reg["ok"] = False
        return reg
    if progreso:
        progreso(0, 1, "homogeneizacion periodica")
    t0 = time.perf_counter()
    C, info = homog_app(B, prot["E_s"], prot["nu"], vox_size=sp,
                        escala_vacio=prot.get("escala_vacio", 1e-6))
    reg["corridas"].append({"corrida": "periodica",
                            "tiempo_s": time.perf_counter() - t0,
                            "solver": info.get("solver"),
                            "residuo_rel": info.get("residuo_rel")})
    if not info.get("ok"):
        reg["fallos"].append({"corrida": "periodica", "msg": info["msg"]})
    reg["C_periodico"] = np.asarray(C)
    reg["C_app"] = np.asarray(C)
    reg["condicion"] = "periodica"
    reg["V_celda_mm3"] = float(np.prod(np.asarray(B.shape) * sp))
    reg["ok"] = not reg["fallos"]
    if progreso:
        progreso(1, 1, "listo")
    return reg


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    """python -m spinpy.fem VOI.vtk [...] --protocolo app --motores ngsolve"""
    import argparse
    import json

    from .io import leer_voi
    ap = argparse.ArgumentParser(
        prog="python -m spinpy.fem",
        description="Ensayos de spinpy con los motores FEM internos.")
    ap.add_argument("vois", nargs="+", help="VOI(s): .vtk, .mat, carpeta o .tif")
    ap.add_argument("--protocolo", nargs="+", default=["app"],
                    choices=[k for k, v in PROTOCOLOS_FEM.items()
                             if v["tipo"] == "compresion"])
    ap.add_argument("--malla", nargs="+", default=["hex8", "tet10"],
                    choices=list(MALLAS))
    ap.add_argument("--motores", nargs="+", default=None,
                    choices=list(motores.MOTORES),
                    help="por omision, todos los instalados")
    ap.add_argument("--analisis", nargs="+", default=["lineal"],
                    choices=list(ANALISIS))
    ap.add_argument("--ejes", default="Z", help="Z, X, Y o XYZ")
    ap.add_argument("--n-hex", type=int, default=N_HEX_DEF)
    ap.add_argument("--n-tet", type=int, default=N_TET_DEF)
    ap.add_argument("--tam", type=float, default=None,
                    help="tamano maximo de tetraedro (mm)")
    ap.add_argument("--sin-correccion", action="store_true",
                    help="no corregir la perdida de volumen del suavizado")
    ap.add_argument("--conv-malla", action="store_true",
                    help="TET10: repetir el lineal con la mitad del tamano")
    ap.add_argument("--material", default="svk", choices=["svk",
                                                          "neohookeano"])
    ap.add_argument("--pasos", type=int, default=1)
    ap.add_argument("--hilos", type=int, default=None)
    ap.add_argument("--salida", default="./fem_resultados")
    a = ap.parse_args(argv)

    mots = a.motores or motores.disponibles()
    salida = Path(a.salida)
    salida.mkdir(parents=True, exist_ok=True)
    print("Motores: " + ", ".join(f"{motores.ETIQUETAS[m]} "
                                  f"{motores.version(m)}" for m in mots))
    op = {"tam_max_mm": a.tam, "corregir_volumen": not a.sin_correccion}
    filas = []
    for ruta in a.vois:
        BW, spc = leer_voi(ruta)
        et = Path(ruta).stem
        for pk in a.protocolo:
            prot = protocolo(pk)
            for mt in a.malla:
                for e in a.ejes.upper():
                    eje = "XYZ".index(e)
                    print(f"== {et} · {prot['nombre']} · {mt} · {e}",
                          flush=True)
                    regs = analizar(
                        BW, spc, prot, malla=mt, analisis=a.analisis, eje=eje,
                        carpeta=salida, n=a.n_hex if mt == "hex8" else a.n_tet,
                        motores_fem=mots, hilos=a.hilos, material=a.material,
                        pasos=a.pasos, opciones_malla=op, etiqueta=et,
                        conv_malla=a.conv_malla, aislado=False)
                    for f in tabla_motores(regs):
                        filas.append(f)
                        print("   " + ", ".join(
                            f"{k}={v:.5g}" if isinstance(v, float)
                            else f"{k}={v}" for k, v in f.items()), flush=True)
                    with open(salida / "resultados_fem.jsonl", "a",
                              encoding="utf-8") as fh:
                        for r in regs:
                            fh.write(json.dumps(registro_json(r),
                                                ensure_ascii=False) + "\n")
    return filas


if __name__ == "__main__":
    main()
