"""
solido.py — Exportacion del spinodoide (o de un VOI) como SOLIDO.

Un STL es solo la piel: describe la frontera, no el interior. Para simular hace
falta una malla VOLUMETRICA. Este modulo produce las dos variantes que tienen
sentido sobre una mascara de voxeles.

HEXAEDRICA (directa desde los voxeles)
--------------------------------------
Cada voxel solido se convierte en un hexaedro de 8 nodos. Es exacta —reproduce
la mascara sin ninguna aproximacion—, instantanea y NO PUEDE FALLAR: no hay
mallador que se atragante. Es el metodo estandar de micro-elementos finitos
para hueso trabecular (Van Rietbergen, Ruegsegger) y es lo que hacen
`localWriteFEBioHex` y `localWriteAbaqusHex` de la app de MATLAB.

Su limitacion: el hexaedro trilineal es un elemento LINEAL, mas rigido a
flexion que uno cuadratico con los mismos grados de libertad, y la superficie
queda escalonada.

TETRAEDRICA TET10 (cuadratica)
------------------------------
Superficie por marching cubes, suavizada, reparada y tetraedralizada. La
superficie sale lisa y el elemento cuadratico representa mejor la flexion, a
cambio de ser mucho mas lento y de introducir tres pasos donde se pierde algo
de material.

TRES COSAS QUE COSTARON DESCUBRIR Y QUE NO SE PUEDEN OMITIR
------------------------------------------------------------
1. SUAVIZAR ES OBLIGATORIO antes de tetraedralizar. Sobre la superficie cruda
   de marching cubes tetgen aborta con `recoversubfaces` o `split_segment`:
   los datos binarios producen miles de triangulos coplanares y aristas
   diminutas que no puede recuperar. No es cosmetica.

2. EL ORDEN DE NODOS DE TETGEN NO ES EL DE C3D10. Medido (distancias
   exactamente cero):

       ranura tetgen   4       5       6       7       8       9
       arista        (2,3)   (0,3)   (0,1)   (1,2)   (1,3)   (0,2)
       C3D10 espera  (0,1)   (1,2)   (0,2)   (0,3)   (1,3)   (2,3)

   De ahi PERM_TETGEN_A_C3D10. Sin la permutacion los elementos quedan
   retorcidos y el analisis corre igualmente devolviendo tensiones sin
   sentido: un error que no da ningun sintoma. Por eso se VERIFICA sobre la
   malla concreta antes de escribir.

3. SISTEMA DE UNIDADES mm-N-MPa. Ni Abaqus ni ANSYS llevan unidades:
   interpretan los numeros en el sistema coherente que elija quien escribe. La
   malla esta en milimetros, asi que el modulo se escribe en MPa. Poner 20e9
   (pascales) declara un material 10^6 veces mas rigido y el analisis corre sin
   protestar.
"""

from __future__ import annotations

import time

import numpy as np
import pyvista as pv

# Esquinas del hexaedro en el orden de C3D8 de Abaqus y SOLID185 de ANSYS:
# cara inferior en sentido antihorario y despues la superior.
ESQUINAS = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                     [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])

# tetgen -> C3D10 / SOLID187
PERM_TETGEN_A_C3D10 = [0, 1, 2, 3, 6, 7, 9, 5, 8, 4]
ARISTAS_C3D10 = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]

# DECIMADO era 0.8. Medido sobre un spinodoide a 48^3 (rho 0.35), descompuesto
# paso a paso frente al volumen del CONTEO de voxeles:
#     superficie cerrada cruda   -3.15 %   (definicion: MC a nivel 0.5, no
#                                           conteo de cubos; no es perdida)
#     + Taubin(20)               -4.3 %
#     + decimate(0.8)            -5.6 %       decimate(0.5): -4.4 %
#     + PyMeshFix                -6.6 %       con 0.5:       -5.0 %
# `volume_preservation=True` de decimate NO cambio ni un decimal, asi que no
# se usa. Con 0.5 la perdida real -la que anaden suavizado, decimado y
# reparacion sobre la superficie cruda- baja del 3.5 % al 1.9 %, a cambio de
# 2.5 veces mas triangulos para tetgen. Por eso la referencia que se reporta y
# que dispara el aviso es la superficie cruda (`vol_pct_MC`), no el conteo.
SUAVIZADO_ITER, SUAVIZADO_BANDA, DECIMADO = 20, 0.05, 0.5


# ---------------------------------------------------------------------------
# Malla hexaedrica
# ---------------------------------------------------------------------------

def malla_hex(BW, spacing):
    """Un hexaedro por voxel solido. Devuelve (nodos, elems, informe).

    Solo se crean los nodos que algun elemento usa: una mascara con BV/TV de
    0.3 generaria de otro modo un 70% de nodos sueltos que ningun elemento
    referencia, y varios preprocesadores los rechazan.
    """
    BW = np.asarray(BW, dtype=bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    nx, ny, nz = BW.shape
    nnx, nny = nx + 1, ny + 1

    i, j, k = np.nonzero(BW)
    n_el = i.size
    if n_el == 0:
        raise ValueError("La mascara no contiene material solido.")

    # Indice global de nodo en la rejilla (nx+1, ny+1, nz+1)
    def nid(a, b, c):
        return a + nnx * b + nnx * nny * c

    conn = np.stack([nid(i + e[0], j + e[1], k + e[2]) for e in ESQUINAS],
                    axis=1)                                   # (n_el, 8)

    usados, conn_compacta = np.unique(conn, return_inverse=True)
    conn_compacta = conn_compacta.reshape(conn.shape)

    ku, resto = np.divmod(usados, nnx * nny)
    ju, iu = np.divmod(resto, nnx)
    nodos = np.stack([iu, ju, ku], axis=1) * spacing[None, :]

    inf = {"tipo": "hex8", "n_nodos": int(nodos.shape[0]),
           "n_elems": int(n_el), "n_gdl": int(3 * nodos.shape[0]),
           "volumen": float(n_el * np.prod(spacing)),
           "BV_voxel": float(BW.sum() * np.prod(spacing))}
    inf["vol_pct_BV"] = 100.0 * inf["volumen"] / inf["BV_voxel"]
    return nodos, conn_compacta, inf


def volumen_hex(nodos, elems):
    """Volumen total, por el jacobiano en el centro de cada hexaedro.

    Sirve de comprobacion independiente: si la conectividad estuviera en un
    orden equivocado los volumenes saldrian negativos.
    """
    P = nodos[elems]                                        # (n, 8, 3)
    dx = np.linalg.norm(P[:, 1] - P[:, 0], axis=1)
    dy = np.linalg.norm(P[:, 3] - P[:, 0], axis=1)
    dz = np.linalg.norm(P[:, 4] - P[:, 0], axis=1)
    return float(np.sum(dx * dy * dz)), int(np.sum(dx * dy * dz <= 0))


# ---------------------------------------------------------------------------
# Malla tetraedrica TET10
# ---------------------------------------------------------------------------

def superficie_cerrada(BW, spacing):
    """Isosuperficie 0.5 cerrada por relleno de vacio.

    El relleno no es solo para cerrar: como en el corte el gradiente es
    puramente normal, las caras de los bordes salen PLANAS en vez de abiertas,
    que es lo que se necesita si despues hay que apoyar la probeta.
    """
    BW = np.asarray(BW, dtype=bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    M = np.pad(BW, 1, constant_values=False)
    g = pv.ImageData(dimensions=M.shape, spacing=tuple(spacing),
                     origin=tuple(-spacing))
    g.point_data["v"] = M.astype(np.float32).flatten(order="F")
    return g.contour([0.5], scalars="v").triangulate().clean()


def _aplanar_tapas(sup, z0, z1, banda, eje=2):
    """Devuelve al plano exacto los nodos que el suavizado saco de las tapas."""
    p = np.array(sup.points, dtype=float, copy=True)
    z = p[:, eje]
    n0 = int((z <= z0 + banda).sum()); n1 = int((z >= z1 - banda).sum())
    p[z <= z0 + banda, eje] = z0
    p[z >= z1 - banda, eje] = z1
    m = sup.copy(); m.points = p
    return m.clean(), {"aplanados_base": n0, "aplanados_techo": n1}


def _planos(BW_shape, spacing, eje):
    """Planos de las dos tapas del eje: los de `superficie_cerrada`."""
    return -0.5 * spacing[eje], (BW_shape[eje] - 0.5) * spacing[eje]


def _area_libre(sup, BW_shape, spacing, ejes):
    """Area de los triangulos que NO estan sobre un plano de tapa."""
    s = sup.triangulate()
    tri = np.asarray(s.faces).reshape(-1, 4)[:, 1:]
    p = np.asarray(s.points, float)
    en_plano = np.zeros(tri.shape[0], bool)
    for e in ejes:
        for v in _planos(BW_shape, spacing, e):
            en_plano |= np.all(np.abs(p[tri, e] - v) <= 1e-9 * spacing[e],
                               axis=1)
    P = p[tri]
    a = 0.5 * np.linalg.norm(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]),
                             axis=1)
    return float(a[~en_plano].sum())


def _desplazar_normal(sup, delta, BW_shape, spacing, ejes):
    """Desplaza la superficie `delta` (mm) por su normal saliente.

    Es la correccion de volumen: el suavizado de Taubin, el decimado y la
    reparacion encogen el solido (medido: -1.9 % frente a la superficie
    cruda). Un desplazamiento uniforme por la normal cambia el volumen en
    delta * area, a primer orden, y delta puede ser NEGATIVO: el suavizado
    encoge lo convexo y agranda lo concavo, y en una estructura dominada por
    poros (una cavidad en un bloque) el solido GANA volumen.

    Los vertices que estaban sobre un plano de tapa vuelven EXACTAMENTE a ese
    plano (solo se mueven dentro de el): recortar al cubo no basta, porque
    con delta < 0 la tapa se meteria hacia dentro, dejaria de ser plana y
    tetgen aborta (medido: violacion de segmento en la cavidad de 16^3).
    """
    p0 = np.asarray(sup.points, float)
    en_plano = []
    for e in ejes:
        a, b = _planos(BW_shape, spacing, e)
        t = 1e-9 * float(spacing[e])
        en_plano.append((e, np.abs(p0[:, e] - a) <= t, a))
        en_plano.append((e, np.abs(p0[:, e] - b) <= t, b))
    # SIN auto_orient_normals: orienta cada cascara hacia fuera de SI MISMA, y
    # la pared de un poro cerrado quedaba apuntando al hueso (la correccion
    # engordaba el solido cuando debia adelgazarlo). La orientacion que deja
    # PyMeshFix es la saliente del solido —su volumen con signo es el del
    # hueso—, y se conserva.
    s = sup.compute_normals(point_normals=True, cell_normals=False,
                            auto_orient_normals=False,
                            consistent_normals=False, split_vertices=False)
    p = p0 + float(delta) * np.asarray(s.point_data["Normals"], float)
    for e in ejes:
        a, b = _planos(BW_shape, spacing, e)
        p[:, e] = np.clip(p[:, e], a, b)
    for e, m_, v in en_plano:
        p[m_, e] = v
    m = sup.copy(); m.points = p
    return m


def malla_tet10(BW, spacing, suavizado=SUAVIZADO_ITER, banda=SUAVIZADO_BANDA,
                decimado=DECIMADO, progreso=None, maxvolume=None,
                minratio=None, ejes_planos=(2,), volumen_objetivo=None):
    """Mascara -> TET10. Devuelve (nodos, elems, superficie, informe).

    Los cuatro ultimos argumentos son opcionales; con su valor por omision el
    resultado es el de siempre.

      maxvolume          volumen maximo de tetraedro (mm^3) para tetgen; None
                         deja que el tamano lo marque la superficie.
      minratio           razon radio-arista maxima de la mejora de calidad
                         de tetgen (su defecto es 2.0).
      ejes_planos        tapas que se devuelven a su plano exacto tras el
                         suavizado. Por omision solo z (base y techo del
                         ensayo). (0, 1, 2) deja las seis caras del cubo
                         planas: lo necesitan las condiciones de borde de la
                         homogeneizacion y la regla que excluye las caras del
                         cubo de la capa superficial.
      volumen_objetivo   si se da (mm^3), la superficie reparada se desplaza
                         por su normal para encerrar ese volumen (dos
                         iteraciones) antes de tetraedralizar.
    """
    import pymeshfix
    import tetgen

    def pr(f, m):
        if progreso:
            progreso(f, m)

    BW = np.asarray(BW, dtype=bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    BV = float(BW.sum() * np.prod(spacing))
    t0 = time.time()

    pr(0.10, "Extrayendo la superficie…")
    sup0 = superficie_cerrada(BW, spacing)
    V_mc = float(sup0.volume)          # referencia honesta de la perdida

    pr(0.30, "Suavizando y decimando…")
    sup = sup0
    if suavizado:
        sup = sup.smooth_taubin(n_iter=int(suavizado), pass_band=float(banda))
    if decimado:
        sup = sup.decimate(float(decimado)).triangulate().clean()

    pr(0.45, "Aplanando las tapas…")
    sup, inf_pl = _aplanar_tapas(sup, -0.5 * spacing[2],
                                 (BW.shape[2] - 0.5) * spacing[2], spacing[2])
    for e in ejes_planos:
        if e != 2:
            sup, _ = _aplanar_tapas(sup, *_planos(BW.shape, spacing, e),
                                    spacing[e], eje=e)

    pr(0.55, "Reparando (estanqueidad)…")
    # remove_smallest_components=False: el defecto de PyMeshFix BORRA todas
    # las cascaras menos la mayor. En una estructura con poros cerrados (una
    # cavidad que no toca las caras del cubo) eso quitaba la pared del poro y
    # tetgen lo rellenaba de tetraedros: el poro salia macizo, sin ningun
    # aviso (medido 2026-09-24 en una cavidad esferica a 16^3: la malla tenia
    # el volumen del cubo entero). Las regiones que quedan dentro de un poro
    # se quitan despues de tetraedralizar (`_quitar_poros`).
    mf = pymeshfix.MeshFix(sup.triangulate().clean())
    mf.repair(remove_smallest_components=False)
    sup = mf.mesh

    correccion = None
    if volumen_objetivo is not None:
        V0 = float(sup.volume)
        Vobj = float(volumen_objetivo)
        antes = sup
        for _ in range(2):
            # Solo se mueve la superficie LIBRE: las caras sobre los planos del
            # cubo vuelven a su plano. Dividir por el area total subestimaba el
            # desplazamiento en una estructura con mucha tapa.
            A = _area_libre(sup, BW.shape, spacing, set(ejes_planos) | {2})
            if A <= 0:
                break
            delta = (Vobj - float(sup.volume)) / A
            sup = _desplazar_normal(sup, delta, BW.shape, spacing,
                                    set(ejes_planos) | {2})
        # El desplazamiento puede CRUZAR entre si puntales finos (delta del
        # orden de 0.1 voxel cuando la perdida es grande): tetgen aborta con
        # una violacion de acceso sobre una superficie que se autointerseca
        # (medido en el VOI proximal de H4 a 32^3). Se repara otra vez.
        mf = pymeshfix.MeshFix(sup.triangulate().clean())
        mf.repair(remove_smallest_components=False)
        sup = mf.mesh
        # Guarda: si el volumen no se acerco al objetivo (orientacion
        # inesperada, superficie rara), se deshace y se declara.
        aplicada = abs(float(sup.volume) - Vobj) < abs(V0 - Vobj)
        if not aplicada:
            sup = antes
        correccion = {"V_antes": V0, "V_despues": float(sup.volume),
                      "V_objetivo": Vobj, "aplicada": bool(aplicada)}

    # Informe de estanqueidad ANTES de tetraedralizar y de escribir: una STL
    # con bordes abiertos no se imprime, y tetgen sobre ella aborta o, peor,
    # rellena lo que no debe. Se comprueba sobre la malla concreta.
    bordes = int(sup.n_open_edges)
    estanca = bool(sup.is_manifold and bordes == 0)
    try:
        n_comp = int(len(np.unique(
            sup.connectivity("all").point_data["RegionId"])))
    except Exception:
        n_comp = -1

    pr(0.70, "Tetraedralizando…")
    tg = tetgen.TetGen(sup)
    opciones = {}
    if maxvolume is not None:
        opciones["maxvolume"] = float(maxvolume)
    if minratio is not None:
        opciones["minratio"] = float(minratio)
    salida = tg.tetrahedralize(order=2, regionattrib=True, **opciones)
    nodos = np.asarray(salida[0], dtype=float)
    elems = np.asarray(salida[1], dtype=np.int64)[:, PERM_TETGEN_A_C3D10]
    nodos, elems, inf_poros = _quitar_poros(nodos, elems, salida[2], BW,
                                            spacing)

    pr(0.95, "Verificando…")
    ok, vinf = verificar_tet10(nodos, elems)
    V_tet = float(np.abs(_vol_tet4(nodos, elems)).sum())
    inf = {"tipo": "tet10", "n_nodos": int(nodos.shape[0]),
           "n_elems": int(elems.shape[0]), "n_gdl": int(3 * nodos.shape[0]),
           "volumen": V_tet, "BV_voxel": BV,
           "vol_pct_BV": 100.0 * V_tet / BV,
           # Frente a la superficie cerrada cruda: es la perdida que de verdad
           # anaden suavizado, decimado y reparacion. La diferencia entre
           # conteo de voxeles y superficie cruda es de definicion.
           "V_mc": V_mc,
           "vol_pct_MC": 100.0 * V_tet / V_mc if V_mc > 0
           else float("nan"),
           "poros": inf_poros,
           "estanca": estanca, "bordes_abiertos": bordes,
           "n_componentes": n_comp,
           "tri_superficie": int(sup.n_cells),
           "orden_ok": bool(ok), "orden": vinf,
           "aplanado": inf_pl, "tiempo_s": time.time() - t0,
           "maxvolume": maxvolume, "minratio": minratio,
           "ejes_planos": [int(e) for e in sorted(set(ejes_planos) | {2})],
           "correccion_volumen": correccion}
    pr(1.0, "Listo.")
    return nodos, elems, sup, inf


def _vol_tet4(nodos, elems):
    P = nodos[elems[:, :4]]
    return np.einsum("ij,ij->i", P[:, 1] - P[:, 0],
                     np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 0])) / 6.0


def _quitar_poros(nodos, elems, atributos, BW, spacing):
    """Quita las regiones de tetgen que caen dentro de un poro cerrado.

    tetgen numera cada region encerrada por la superficie (`regionattrib`).
    Una cavidad cerrada es una region propia; se decide si es hueso o poro
    por VOTO de volumen de sus tetraedros: la fase de la mascara en el voxel
    que contiene cada centroide (voxel i centrado en i*h, la convencion de
    `superficie_cerrada`). Se vota por region y no por elemento porque un
    tetraedro pegado a la superficie suavizada puede caer en un voxel de la
    otra fase sin que su region lo sea.
    """
    at = np.asarray(atributos).ravel()
    if at.size != elems.shape[0]:
        return nodos, elems, {"regiones": None}
    v = np.abs(_vol_tet4(nodos, elems))
    c = nodos[elems[:, :4]].mean(axis=1)
    idx = np.rint(c / spacing[None, :]).astype(np.int64)
    idx = np.clip(idx, 0, np.array(BW.shape) - 1)
    es_hueso = BW[idx[:, 0], idx[:, 1], idx[:, 2]]
    quitar = np.zeros(elems.shape[0], bool)
    regiones = []
    for a in np.unique(at):
        m = at == a
        frac = float(v[m & es_hueso].sum() / max(v[m].sum(), 1e-300))
        poro = frac < 0.5
        regiones.append({"region": float(a), "V_mm3": float(v[m].sum()),
                         "frac_hueso": frac, "poro": bool(poro)})
        quitar |= m & poro
    if quitar.any():
        elems = elems[~quitar]
        usados, inv = np.unique(elems, return_inverse=True)
        nodos, elems = nodos[usados], inv.reshape(elems.shape)
    return nodos, elems, {"regiones": regiones,
                          "n_poros_quitados": int(sum(r["poro"]
                                                      for r in regiones)),
                          "V_poros_quitados_mm3": float(v[quitar].sum())}


def verificar_tet10(nodos, elems, tol_rel=1e-6, n_muestra=2000):
    """Comprueba que cada nodo intermedio esta en el medio de SU arista."""
    nodos = np.asarray(nodos, float)
    elems = np.asarray(elems, dtype=np.int64)
    if elems.shape[1] != 10:
        return False, {"msg": f"{elems.shape[1]} nodos por elemento, no 10."}
    rng = np.random.default_rng(0)
    E = elems[rng.choice(elems.shape[0],
                         size=min(n_muestra, elems.shape[0]), replace=False)]
    esc = float(np.linalg.norm(nodos.max(0) - nodos.min(0)))
    peor = 0.0
    for k, (a, b) in enumerate(ARISTAS_C3D10):
        med = 0.5 * (nodos[E[:, a]] + nodos[E[:, b]])
        peor = max(peor, float(np.linalg.norm(nodos[E[:, 4 + k]] - med,
                                              axis=1).max() / esc))
    ok = peor <= tol_rel
    return ok, {"desviacion_rel": peor,
                "msg": ("Orden de nodos correcto (C3D10)." if ok else
                        "ORDEN DE NODOS INCORRECTO: no exportar.")}
