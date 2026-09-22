"""
elastic.py — Homogeneizacion numerica periodica sobre la rejilla de voxeles.

Port de:
    localHex8Ke              (AppFinal_V2.m:6121-6167)
    localHomogenizeVoxel     (AppFinal_V2.m:6192-6430)
    localEngineeringConstants(AppFinal_V2.m:6452-...)

Resuelve el problema de celda unidad de verdad —elementos hexaedricos sobre la
rejilla de voxeles, condiciones de contorno periodicas, seis casos de
deformacion macroscopica unitaria— y obtiene

    C_ij = (1/|V|) * sum_e (u0_i - chi_i)' ke (u0_j - chi_j)

siendo chi el campo corrector. Referencia: Andreassen & Andreasen,
Comput. Mater. Sci. 83 (2014) 488-495.

SUSTITUYE al tensor heuristico que la app tenia antes (E_eff = E_s*rho^2 con un
reparto lineal segun los angulos y un factor 0.9 sin justificacion), en el que
la rigidez no derivaba de la geometria en ningun momento: dos estructuras con
la misma densidad daban el mismo tensor aunque su topologia fuera distinta.

CONVENCIONES QUE NO SE PUEDEN CAMBIAR
-------------------------------------
  * Voigt [xx yy zz yz xz xy] con deformaciones angulares de INGENIERIA, en
    todo el archivo y en todo el proyecto.
  * El vacio se modela con rigidez 1e-6*E_s en lugar de eliminarse: evita que
    islas de material desconectadas hagan singular la matriz, a costa de un
    sesgo despreciable frente al solido.
  * La rejilla de NODOS tiene el mismo tamano que la de elementos: el nodo
    nx+1 es, por periodicidad, el nodo 1.
  * Los arrays se recorren en orden Fortran (columna primero) para que el
    indice de elemento coincida con el de MATLAB.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, cg, spilu, splu

VOID_SCALE = 1e-6      # rigidez relativa del vacio (valor de la app)
UMBRAL_DIRECTO = 6000  # GDL por debajo de los cuales se resuelve directo
#
# ELECCION DE SOLVER — medida, no heredada
# ----------------------------------------
# MATLAB usa `decomposition(...,'chol')` por debajo de 1.2e4 GDL y `pcg` con
# Cholesky incompleta por encima. Ninguna de las dos piezas existe igual aqui:
#
#   * El LU directo de scipy sufre el relleno intrinseco de las rejillas 3D:
#     el coste crece mucho mas deprisa que el numero de incognitas.
#   * scipy no tiene `ichol`. Su `spilu` es una ILU no simetrica que con el
#     contraste 1e6 entre solido y vacio se degrada tanto que el CG no
#     converge: a 32^3 se rindio tras 2294 s sin alcanzar la tolerancia.
#
# La salida es multigrid algebraico (pyamg), que ataca justamente esta clase de
# problema. Tiempos MEDIDOS sobre el mismo spinodoide, maquina en reposo (una
# medicion anterior con otros procesos compitiendo dio cifras hasta 4 veces
# mayores; estas son las buenas):
#
#     lado    GDL      directo    AMG+CG    factor
#     12^3    5.184       1.8 s     1.8 s     1.0
#     16^3   12.288      14.0 s     5.7 s     2.5
#     20^3   24.000      59.7 s    14.8 s     4.0
#     24^3   41.472     223.9 s    34.3 s     6.5
#     28^3   65.856    1696.6 s    55.6 s    30.5
#     32^3   98.304       ---      99.2 s     ---
#
# Los dos caminos coinciden en C_h dentro de ~1e-15 relativo. La ventaja del
# multigrid crece con el tamano, que es exactamente lo que se necesita: por
# debajo de UMBRAL_DIRECTO empatan y el directo ademas es exacto, asi que se
# conserva ahi.


# ---------------------------------------------------------------------------
# Matriz de rigidez elemental
# ---------------------------------------------------------------------------

def hex8_ke(dx, dy, dz, E, nu):
    """Rigidez de un hexaedro trilineal de 8 nodos, Gauss 2x2x2.

    Se calcula numericamente en lugar de usar una matriz precomputada: mas
    largo, pero verificable y sin riesgo de erratas de transcripcion.
    """
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))

    D = np.array([
        [lam + 2 * mu, lam, lam, 0, 0, 0],
        [lam, lam + 2 * mu, lam, 0, 0, 0],
        [lam, lam, lam + 2 * mu, 0, 0, 0],
        [0, 0, 0, mu, 0, 0],
        [0, 0, 0, 0, mu, 0],
        [0, 0, 0, 0, 0, mu],
    ], dtype=float)

    xn = np.array([-1, 1, 1, -1, -1, 1, 1, -1], dtype=float)
    yn = np.array([-1, -1, 1, 1, -1, -1, 1, 1], dtype=float)
    zn = np.array([-1, -1, -1, -1, 1, 1, 1, 1], dtype=float)

    g = np.array([-1 / np.sqrt(3), 1 / np.sqrt(3)])
    detJ = (dx / 2) * (dy / 2) * (dz / 2)
    invJ = np.diag([2 / dx, 2 / dy, 2 / dz])

    ke = np.zeros((24, 24))
    for xi in g:
        for eta in g:
            for zet in g:
                dNdxi = 0.125 * np.vstack([
                    xn * (1 + eta * yn) * (1 + zet * zn),
                    yn * (1 + xi * xn) * (1 + zet * zn),
                    zn * (1 + xi * xn) * (1 + eta * yn),
                ])
                dNdx = invJ @ dNdxi              # 3x8

                B = np.zeros((6, 24))
                B[0, 0::3] = dNdx[0]                              # eps_xx
                B[1, 1::3] = dNdx[1]                              # eps_yy
                B[2, 2::3] = dNdx[2]                              # eps_zz
                B[3, 1::3] = dNdx[2]; B[3, 2::3] = dNdx[1]        # gamma_yz
                B[4, 0::3] = dNdx[2]; B[4, 2::3] = dNdx[0]        # gamma_xz
                B[5, 0::3] = dNdx[1]; B[5, 1::3] = dNdx[0]        # gamma_xy

                ke += B.T @ D @ B * detJ
    return (ke + ke.T) / 2


# ---------------------------------------------------------------------------

def _edof_periodico(nx, ny, nz):
    """Grados de libertad por elemento, con conectividad periodica.

    Orden de nodos identico al de MATLAB:
        (i,j,k) (i+1,j,k) (i+1,j+1,k) (i,j+1,k)
        (i,j,k+1) (i+1,j,k+1) (i+1,j+1,k+1) (i,j+1,k+1)
    """
    ex, ey, ez = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz),
                             indexing="ij")
    ex = ex.ravel(order="F"); ey = ey.ravel(order="F"); ez = ez.ravel(order="F")

    def nodo(a, b, c):
        return (a % nx) + nx * (b % ny) + nx * ny * (c % nz)

    nodos = np.stack([
        nodo(ex, ey, ez), nodo(ex + 1, ey, ez),
        nodo(ex + 1, ey + 1, ez), nodo(ex, ey + 1, ez),
        nodo(ex, ey, ez + 1), nodo(ex + 1, ey, ez + 1),
        nodo(ex + 1, ey + 1, ez + 1), nodo(ex, ey + 1, ez + 1),
    ], axis=1)                                    # (nel, 8)

    edof = np.empty((nodos.shape[0], 24), dtype=np.int64)
    for a in range(8):
        edof[:, 3 * a + 0] = 3 * nodos[:, a] + 0
        edof[:, 3 * a + 1] = 3 * nodos[:, a] + 1
        edof[:, 3 * a + 2] = 3 * nodos[:, a] + 2
    return edof


def _u0_casos(dx, dy, dz):
    """Desplazamientos nodales de los seis casos de deformacion unitaria."""
    xn = np.array([0, 1, 1, 0, 0, 1, 1, 0], float) * dx
    yn = np.array([0, 0, 1, 1, 0, 0, 1, 1], float) * dy
    zn = np.array([0, 0, 0, 0, 1, 1, 1, 1], float) * dz

    u0 = np.zeros((24, 6))
    for c in range(6):
        e = np.zeros(6); e[c] = 1.0
        # Voigt [xx yy zz yz xz xy] -> tensor (las angulares van a la mitad)
        Eps = np.array([[e[0], e[5] / 2, e[4] / 2],
                        [e[5] / 2, e[1], e[3] / 2],
                        [e[4] / 2, e[3] / 2, e[2]]])
        for a in range(8):
            u0[3 * a:3 * a + 3, c] = Eps @ np.array([xn[a], yn[a], zn[a]])
    return u0


def _ensamblar(keS, edof, escala, ndof, bloques=8):
    """Matriz global dispersa, ensamblada por bloques de elementos.

    Ensamblar de una vez requiere nel*576 tripletas: a 48^3 son 64 millones de
    entradas, del orden de 1 GB. Trocear acota el pico de memoria sin cambiar
    el resultado, porque la suma de matrices dispersas es exacta.
    """
    nel = edof.shape[0]
    ke_plana = keS.ravel()
    K = sparse.csr_matrix((ndof, ndof))
    paso = max(1, nel // bloques)
    for i0 in range(0, nel, paso):
        i1 = min(i0 + paso, nel)
        ed = edof[i0:i1]
        filas = np.repeat(ed, 24, axis=1).ravel()
        cols = np.tile(ed, (1, 24)).ravel()
        vals = (ke_plana[None, :] * escala[i0:i1, None]).ravel()
        K = K + sparse.coo_matrix((vals, (filas, cols)),
                                  shape=(ndof, ndof)).tocsr()
    return K


# ---------------------------------------------------------------------------

def homogeneizar(BW, E_s=1.0, nu_s=0.3, vox_size=1.0, tol=1e-8,
                 escala_vacio=VOID_SCALE, verbose=False):
    """Tensor de rigidez homogeneizado C_h (6x6) en orden de Voigt.

    Devuelve (Ch, info). Si algo impide resolver, Ch es NaN e `info['msg']`
    explica por que — nunca se devuelve un tensor silenciosamente equivocado.

    escala_vacio : rigidez de la fase `False` relativa a la solida. El valor
        por defecto 1e-6 es el de la app y representa vacio. Se expone como
        parametro porque permite VERIFICAR la funcion contra el promedio de
        Backus con dos fases reales: con vacio a 1e-6 casi todas las entradas
        del tensor exacto valen ~0 y la comparacion no distingue una
        implementacion correcta de una equivocada.
    """
    BW = np.asarray(BW, dtype=bool)
    vox = np.atleast_1d(np.asarray(vox_size, dtype=float)).ravel()
    if vox.size == 1:
        vox = np.repeat(vox, 3)
    dx, dy, dz = vox

    nx, ny, nz = BW.shape
    nel = nx * ny * nz
    ndof = 3 * nel

    info = {"ok": False, "msg": "", "n_elem": nel, "n_dof": ndof,
            "solver": "", "rho_solido": float(BW.mean())}
    Ch = np.full((6, 6), np.nan)

    if not BW.any():
        info["msg"] = "La estructura no contiene material solido."
        return Ch, info

    keS = hex8_ke(dx, dy, dz, E_s, nu_s)
    edof = _edof_periodico(nx, ny, nz)

    escala = np.full(nel, float(escala_vacio))
    escala[BW.ravel(order="F")] = 1.0

    K = _ensamblar(keS, edof, escala, ndof)
    K = ((K + K.T) / 2).tocsc()

    u0 = _u0_casos(dx, dy, dz)
    fe = keS @ u0                                     # (24, 6)

    Fmat = np.zeros((ndof, 6))
    for c in range(6):
        contrib = escala[:, None] * fe[:, c][None, :]      # (nel, 24)
        Fmat[:, c] = np.bincount(edof.ravel(), weights=contrib.ravel(),
                                 minlength=ndof)

    # K es singular solo por las tres traslaciones de solido rigido: se fija el
    # nodo 0 y el sistema queda definido positivo.
    libres = np.arange(3, ndof)
    chi = np.zeros((ndof, 6))
    Kff = K[libres][:, libres].tocsc()

    resuelto = False
    if ndof <= UMBRAL_DIRECTO:
        try:
            lu = splu(Kff)
            for c in range(6):
                chi[libres, c] = lu.solve(Fmat[libres, c])
            info["solver"] = "LU directo"
            resuelto = np.all(np.isfinite(chi))
        except Exception:
            resuelto = False

    if not resuelto:
        Kff_csr = Kff.tocsr()
        ml = None
        try:
            import pyamg
            ml = pyamg.smoothed_aggregation_solver(Kff_csr, max_coarse=500)
            nombre = "AMG (agregacion suavizada) + CG"
        except ImportError:
            ml = None
        except Exception:
            ml = None

        peor = 0.0
        if ml is not None:
            for c in range(6):
                b = Fmat[libres, c]
                x = ml.solve(b, tol=tol, maxiter=500, accel="cg")
                chi[libres, c] = x
                nb = np.linalg.norm(b)
                peor = max(peor, np.linalg.norm(b - Kff_csr @ x) /
                           (nb if nb > 0 else 1.0))
        else:
            # Respaldo sin pyamg: ILU y, si falla, Jacobi. Documentado como
            # poco fiable con contraste alto — se conserva para no quedarse
            # sin ninguna via, no porque sea buena.
            M = None
            nombre = ""
            for droptol in (1e-4, 1e-3, 1e-2):
                try:
                    ilu = spilu(Kff.tocsc(), drop_tol=droptol, fill_factor=10)
                    M = LinearOperator(Kff.shape, ilu.solve)
                    nombre = f"CG + ILU (drop_tol {droptol:.0e})"
                    break
                except Exception:
                    M = None
            if M is None:
                d = np.asarray(Kff.diagonal(), dtype=float)
                d[d <= 0] = np.finfo(float).eps
                M = LinearOperator(Kff.shape, lambda v, d=d: v / d)
                nombre = "CG + Jacobi (diagonal)"
            for c in range(6):
                b = Fmat[libres, c]
                x, _ = cg(Kff_csr, b, rtol=tol, maxiter=20000, M=M)
                chi[libres, c] = x
                nb = np.linalg.norm(b)
                peor = max(peor, np.linalg.norm(b - Kff_csr @ x) /
                           (nb if nb > 0 else 1.0))

        info["solver"] = nombre
        info["residuo_rel"] = float(peor)

        # F10: hay que COMPROBAR el residuo. Un solver iterativo no lanza error
        # al no converger: devuelve el mejor iterado, que pasa cualquier
        # comprobacion de isfinite y da por bueno un tensor equivocado.
        if not np.isfinite(peor) or peor > max(1e-6, 100 * tol):
            info["msg"] = (f"El solver iterativo no convergio (residuo relativo "
                           f"peor = {peor:.1e}). El tensor no es fiable; prueba "
                           f"una resolucion menor.")
            return Ch, info
        resuelto = np.all(np.isfinite(chi))

    if not resuelto:
        info["msg"] = "La solucion del sistema no es finita."
        return Ch, info

    # --- Tensor homogeneizado ---------------------------------------------
    Vtot = (nx * dx) * (ny * dy) * (nz * dz)
    difs = [u0[:, c][None, :] - chi[edof, c] for c in range(6)]

    Ch = np.zeros((6, 6))
    for i in range(6):
        Ai = difs[i] @ keS                      # (nel, 24)
        for j in range(i, 6):
            Ch[i, j] = float(np.sum(np.sum(Ai * difs[j], axis=1) * escala) / Vtot)
            Ch[j, i] = Ch[i, j]
    Ch = (Ch + Ch.T) / 2

    info["ok"] = True
    info["msg"] = (f"Homogeneizacion completada. {nel} elementos, {ndof} GDL, "
                   f"fraccion solida {info['rho_solido']:.3f} ({info['solver']}).")
    return Ch, info


# ---------------------------------------------------------------------------

def remuestrear_bw(BW, spacing, n_max):
    """Reduce la mascara a n_max voxeles por eje. Port de localResampleBW.

    Submuestreo por vecino mas proximo, no por promediado: la mascara debe
    seguir siendo binaria. El spacing se corrige por el factor real alcanzado,
    no por el pedido, porque `round` sobre linspace no da exactamente n_max.

    El coste de la homogeneizacion crece con el CUBO del lado, asi que este
    recorte no es cosmetico: es lo que hace viable calcular un tensor sobre un
    VOI de 97^3.
    """
    BW = np.asarray(BW, dtype=bool)
    spacing = np.asarray(spacing, dtype=float).ravel()
    if spacing.size != 3:
        spacing = np.ones(3)
    sz = np.array(BW.shape)
    if sz.max() <= n_max:
        return BW, spacing

    factor = n_max / sz
    idx = [np.round(np.linspace(0, sz[k] - 1,
                                max(2, int(round(sz[k] * factor[k]))))).astype(int)
           for k in range(3)]
    BWr = BW[np.ix_(idx[0], idx[1], idx[2])]
    spc = spacing * (sz / np.array(BWr.shape))
    return BWr, spc


def constantes_ingenieria(C):
    """Constantes de ingenieria a partir de C (6x6), via S = inv(C).

    Orden de Voigt [xx yy zz yz xz xy] con deformaciones angulares de
    ingenieria, de modo que S[3,3] = 1/G_yz directamente.
    """
    C = np.asarray(C, dtype=float)
    S = np.linalg.inv(C)
    ec = {
        "Ex": 1.0 / S[0, 0], "Ey": 1.0 / S[1, 1], "Ez": 1.0 / S[2, 2],
        "Gyz": 1.0 / S[3, 3], "Gxz": 1.0 / S[4, 4], "Gxy": 1.0 / S[5, 5],
        "nu_xy": -S[0, 1] / S[0, 0], "nu_xz": -S[0, 2] / S[0, 0],
        "nu_yz": -S[1, 2] / S[1, 1],
    }
    E = np.array([ec["Ex"], ec["Ey"], ec["Ez"]])
    ec["E_medio"] = float(E.mean())
    ec["anisotropia_E"] = float(E.max() / E.min()) if E.min() > 0 else np.nan
    return ec


# ---------------------------------------------------------------------------
# Comparar tensores: distancia log-euclidea, ejes materiales, clase de simetria
# ---------------------------------------------------------------------------
#
# POR QUE HACE FALTA MAS QUE Ez/Es Y Ez/Ex
# ----------------------------------------
# La etapa D del ajuste compara dos escalares. Wang et al. (2025, Materials &
# Design) igualan las nueve constantes de un hueso y Otto et al. (2025,
# GAMM-Mitteilungen) el tensor entero con una perdida logaritmica, y los dos
# muestran que dos estructuras con el mismo Ez pueden tener cortantes y
# anisotropia muy distintos. Lo que sigue da las piezas para comparar el
# tensor completo sin que la orientacion se cuele en la comparacion.

# Voigt [xx yy zz yz xz xy] -> par de indices del tensor de orden 2.
_VOIGT = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))
_W_MANDEL = np.array([1.0, 1.0, 1.0, np.sqrt(2.0), np.sqrt(2.0), np.sqrt(2.0)])


def voigt_a_tensor4(C):
    """C (6x6, Voigt con deformaciones angulares de ingenieria) -> C_ijkl.

    Con esa convencion C_IJ = C_ijkl sin factores, que es lo que hace la
    conversion trivial y la razon por la que el proyecto la fija.
    """
    C = np.asarray(C, float)
    T = np.zeros((3, 3, 3, 3))
    for I, (i, j) in enumerate(_VOIGT):
        for J, (k, l) in enumerate(_VOIGT):
            v = C[I, J]
            for a, b in ((i, j), (j, i)):
                for c, d in ((k, l), (l, k)):
                    T[a, b, c, d] = v
    return T


def tensor4_a_voigt(T):
    T = np.asarray(T, float)
    C = np.zeros((6, 6))
    for I, (i, j) in enumerate(_VOIGT):
        for J, (k, l) in enumerate(_VOIGT):
            C[I, J] = T[i, j, k, l]
    return C


def en_ejes(C, V):
    """Componentes de C en la base ortonormal cuyas COLUMNAS son V.

    C'_ijkl = V_pi V_qj V_rk V_sl C_pqrs. Un V con determinante -1 da el mismo
    tensor que su rotacion: un tensor de orden par no cambia por inversion.
    """
    V = np.asarray(V, float)
    T = voigt_a_tensor4(C)
    Tp = np.einsum("pi,qj,rk,sl,pqrs->ijkl", V, V, V, V, T, optimize=True)
    return tensor4_a_voigt(Tp)


def a_mandel(C):
    """Voigt de ingenieria -> Mandel: C_M = W C W, W = diag(1,1,1,r2,r2,r2).

    En Mandel el tensor es una matriz 6x6 de verdad —sus autovalores son los de
    la aplicacion lineal— y una rotacion es una matriz ortogonal. En Voigt
    ninguna de las dos cosas es cierta, y por eso el logaritmo se toma aqui.
    """
    return np.asarray(C, float) * np.outer(_W_MANDEL, _W_MANDEL)


def _logm_sdp(M):
    M = (np.asarray(M, float) + np.asarray(M, float).T) / 2.0
    lam, V = np.linalg.eigh(M)
    if not np.all(np.isfinite(lam)) or np.any(lam <= 0):
        raise ValueError("el tensor no es definido positivo: no tiene "
                         "logaritmo real")
    return (V * np.log(lam)) @ V.T


def distancia_log_euclidea(C1, C2):
    """|| log C1 - log C2 ||_F, con los tensores en notacion de Mandel.

    PROPIEDADES QUE LA HACEN LA ADECUADA (y que el bloque 17 comprueba):
      * es simetrica y vale 0 solo si C1 = C2;
      * d(aC, C) = sqrt(6) |ln a|: una estructura el doble de rigida en todo
        esta a la misma distancia que una la mitad, que es lo que pide un
        error RELATIVO;
      * d(C1, C2) = d(S1, S2), porque log(C^-1) = -log C: igualar la rigidez
        y la flexibilidad a la vez, que es el motivo por el que Otto et al.
        (2025) eligen una perdida logaritmica;
      * no cambia si se gira a los DOS con la misma rotacion.

    Lo que NO es: invariante si se gira solo uno. Para comparar sin que cuente
    la orientacion hay que llevar los dos a sus ejes materiales antes
    (`constantes_principales`), que es lo que hace el ajuste.
    """
    return float(np.linalg.norm(_logm_sdp(a_mandel(C1))
                                - _logm_sdp(a_mandel(C2))))


def _autovectores(M, tol):
    """Autovectores y numero de autovalores DISTINTOS (1, 2 o 3)."""
    lam, V = np.linalg.eigh((M + M.T) / 2.0)
    esc = max(float(np.max(np.abs(lam))), 1e-300)
    dif = np.diff(np.sort(lam)) / esc
    return V, 1 + int(np.sum(dif > tol))


def ejes_materiales(C, tol=1e-6):
    """Ejes de simetria material de C (columnas de una matriz 3x3).

    Se toman de los autovectores del tensor DILATACIONAL d_ij = C_ijkk o del
    de VOIGT v_ij = C_ikjk (Cowin 1987): para un material ortotropo los dos son
    diagonales en sus ejes de simetria. Se usa el que tenga mas autovalores
    distintos, y el dilatacional en caso de empate.

    Con DOS autovalores iguales (material transversalmente isotropo) el
    autovector distinto es el eje de simetria y los otros dos son una base
    cualquiera del plano; para ese material cualquier base del plano es
    exacta, asi que no hay nada que decidir. La primera version exigia tres
    autovalores distintos y devolvia los ejes del laboratorio para un laminado
    girado: sus constantes no se recuperaban (bloque 17).

    Si los dos tensores son isotropos —material cubico o isotropo— no fijan
    ejes y se devuelven los del laboratorio. En las estructuras de este
    proyecto, que viven en una rejilla de voxeles, los ejes cubicos coinciden
    con los de la rejilla, asi que es el supuesto razonable; se declara con
    `fuente`.

    Devuelve (V, fuente) con fuente "dilatacional" | "voigt" | "laboratorio".
    """
    T = voigt_a_tensor4(C)
    Vd, nd = _autovectores(np.einsum("ijkk->ij", T), tol)
    Vv, nv = _autovectores(np.einsum("ikjk->ij", T), tol)
    if nd > 1 and nd >= nv:
        return Vd, "dilatacional"
    if nv > 1:
        return Vv, "voigt"
    return np.eye(3), "laboratorio"


def constantes_principales(C):
    """Constantes de ingenieria en los ejes MATERIALES, ordenados E1 >= E2 >= E3.

    Con eso dos tensores se comparan por lo que son y no por como estan
    colocados: el eje mas rigido de uno contra el mas rigido del otro. La
    orientacion se trata aparte, en la alineacion del ajuste.

    `desviacion_ortotropa` es la fraccion de la norma (Mandel) que queda en
    las entradas que un material ortotropo tendria a cero una vez en sus ejes.
    Si no es pequena, las nueve constantes no describen el tensor.
    """
    V, fuente = ejes_materiales(C)
    Cp = en_ejes(C, V)
    ec = constantes_ingenieria(Cp)
    orden = np.argsort([-ec["Ex"], -ec["Ey"], -ec["Ez"]], kind="stable")
    V = V[:, orden]
    Cp = en_ejes(C, V)
    ec = constantes_ingenieria(Cp)
    M = a_mandel(Cp)
    fuera = M.copy()
    fuera[:3, :3] = 0.0
    fuera[3, 3] = fuera[4, 4] = fuera[5, 5] = 0.0
    nM = float(np.linalg.norm(M))
    return {
        "E1": ec["Ex"], "E2": ec["Ey"], "E3": ec["Ez"],
        "G23": ec["Gyz"], "G13": ec["Gxz"], "G12": ec["Gxy"],
        "nu12": ec["nu_xy"], "nu13": ec["nu_xz"], "nu23": ec["nu_yz"],
        "ejes": V, "fuente_ejes": fuente, "C_principal": Cp,
        "desviacion_ortotropa": (float(np.linalg.norm(fuera)) / nM
                                 if nM > 0 else np.nan),
    }


def clase_anisotropia(C, tol=0.05):
    """Clase de simetria elastica de C, con una tolerancia relativa.

    CRITERIO PROPIO, DECLARADO. Otto et al. (2025) proponen un metodo para
    asignar la clase de anisotropia de un tensor; no se ha podido leer su
    texto completo, asi que esto NO lo reproduce. Es un criterio por
    tolerancia sobre las constantes en los ejes materiales:

      sin_simetria_ortotropa  desviacion_ortotropa > tol
      isotropa                E1=E2=E3, G iguales y G = E / (2 (1 + nu))
      cubica                  E1=E2=E3 y G iguales, sin la relacion anterior
      transversal             dos E iguales, con sus G y la relacion en el plano
      ortotropa               el resto

    "Iguales" es |a - b| <= tol * max(|a|, |b|). Con la dispersion entre
    realizaciones de un spinodoide (varios por ciento) una tolerancia menor
    que esa convierte ruido en clase.
    """
    p = constantes_principales(C)
    E = [p["E1"], p["E2"], p["E3"]]
    nu = {(0, 1): p["nu12"], (0, 2): p["nu13"], (1, 2): p["nu23"]}
    # Cortante entre los ejes a y b (base 0): G23 es el de los ejes 1 y 2.
    Gpar = {frozenset((1, 2)): p["G23"], frozenset((0, 2)): p["G13"],
            frozenset((0, 1)): p["G12"]}

    def igual(a, b):
        return abs(a - b) <= tol * max(abs(a), abs(b), 1e-300)

    out = {"clase": "ortotropa", "eje": None, "tol": float(tol),
           "desviacion_ortotropa": p["desviacion_ortotropa"],
           "E": E, "G": [p["G23"], p["G13"], p["G12"]]}
    if not np.isfinite(p["desviacion_ortotropa"]) or \
            p["desviacion_ortotropa"] > tol:
        out["clase"] = "sin_simetria_ortotropa"
        return out

    g = [p["G23"], p["G13"], p["G12"]]
    if igual(E[0], E[1]) and igual(E[1], E[2]) and igual(E[0], E[2]):
        if igual(g[0], g[1]) and igual(g[1], g[2]) and igual(g[0], g[2]):
            nu_m = float(np.mean([p["nu12"], p["nu13"], p["nu23"]]))
            E_m = float(np.mean(E))
            iso = igual(float(np.mean(g)), E_m / (2.0 * (1.0 + nu_m)))
            out["clase"] = "isotropa" if iso else "cubica"
            return out

    # Transversal: el eje de simetria es el que tiene la E distinta. En el
    # plano (a, b) perpendicular: E_a = E_b, G_ac = G_bc y G_ab = E/(2(1+nu_ab)).
    for eje, (a, b) in ((2, (0, 1)), (0, (1, 2)), (1, (0, 2))):
        if not igual(E[a], E[b]):
            continue
        G_ac = Gpar[frozenset((a, eje))]
        G_bc = Gpar[frozenset((b, eje))]
        G_ab = Gpar[frozenset((a, b))]
        nu_ab = nu[tuple(sorted((a, b)))]
        if igual(G_ac, G_bc) and igual(G_ab, E[a] / (2.0 * (1.0 + nu_ab))):
            out["clase"] = "transversal"
            out["eje"] = int(eje)
            return out
    return out


# ---------------------------------------------------------------------------

def backus_laminado(fracciones, E, nu):
    """Tensor EXACTO de un laminado de capas isotropas normales a z.

    Promedio de Backus (1962). Es la solucion analitica contra la que se
    verifica `homogeneizar`: para un apilado de capas planas la solucion de
    elementos finitos debe reproducirla, porque el campo corrector es constante
    dentro de cada capa y el hexaedro trilineal lo representa exactamente.

    fracciones : fracciones de volumen de cada capa (suman 1)
    E, nu      : modulos y coeficientes de Poisson de cada capa
    """
    f = np.asarray(fracciones, float)
    E = np.asarray(E, float)
    nu = np.asarray(nu, float)
    f = f / f.sum()

    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    A = lam + 2 * mu

    prom = lambda v: float(np.sum(f * v))

    inv_A = prom(1.0 / A)
    lam_A = prom(lam / A)
    C33 = 1.0 / inv_A
    C13 = C33 * lam_A
    C11 = prom(4 * mu * (lam + mu) / A) + C33 * lam_A**2
    C12 = prom(2 * mu * lam / A) + C33 * lam_A**2
    C44 = 1.0 / prom(1.0 / mu)      # cortante fuera del plano (yz, xz)
    C66 = prom(mu)                  # cortante en el plano xy

    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = C11
    C[2, 2] = C33
    C[0, 1] = C[1, 0] = C12
    C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = C13
    C[3, 3] = C[4, 4] = C44
    C[5, 5] = C66
    return C
