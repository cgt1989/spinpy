"""
dual_lattice.py — Generador dual-lattice VOXEL-PRIMERO.

Segunda familia de microestructuras del proyecto, junto al spinodoide de
`grf.py`. Es la que Vafaeefar et al. (2022) introducen como modelo de hueso
trabecular y la UNICA de las tres que estudian (gyroide, spinodoide,
dual-lattice) que reproduce la rigidez del hueso al que se calibra:

    E3 (MPa)   hueso 1188   gyroide 887   spinodoide 714   dual-lattice 1313

con BV/TV, Tb.Th y DA igualados en las tres. El spinodoide ajustado a nuestro
VOI equino es ~8 veces mas blando que el hueso (ver Estudio_Discriminadores,
"problema abierto"); el dual-lattice es la familia con la que contrastar si esa brecha es
propia del spinodoide o de cualquier estructura calibrada solo con BV/TV, Tb.Th
y DA.

LA CONSTRUCCION (Vafaeefar et al. 2022, seccion 2.1.3)
------------------------------------------------------
  1. Se reparte una nube de puntos en el dominio y se triangula en tetraedros
     (Delaunay). Un factor de ESTIRAMIENTO por eje alarga las celdas, y con
     ellas las barras: es lo que da la anisotropia.
  2. El ESQUELETO es la red dual: el centroide de cada tetraedro unido al
     centroide de cada una de sus cuatro caras. Dos tetraedros vecinos
     comparten cara, asi que sus barras se encuentran en ella y forman un
     puntal centroide-cara-centroide. Cada nodo tiene exactamente CUATRO
     puntales: es la conectividad nodal 4-N que Vafaeefar senalan como la
     razon de su rigidez, frente a la 3-N del hueso.
  3. El puntal se engruesa con un radio.

Lo que cambia respecto al original, y por que:

  * GIBBON (`dualLattice.m`) construye una MALLA de superficie con barras de
    seccion triangular, y la voxeliza despues. Aqui se va directo a voxeles,
    como en `grf.py`: la mascara es  {x : distancia(x, esqueleto) <= r}.
    Las barras quedan de seccion REDONDA. Es una diferencia de forma, no de
    topologia, y hay que declararla al comparar con cifras de GIBBON.

  * La triangulacion no es de TetGen sino `scipy.spatial.Delaunay` sobre una
    rejilla con perturbacion aleatoria (`irregularidad`). Una rejilla perfecta
    es degenerada para Delaunay (ocho puntos cocirculares por celda); la
    perturbacion lo evita y es ademas la fuente de la aleatoriedad de la
    familia, como las fases en el spinodoide.

  * El radio NO es un parametro de entrada: sale de la DENSIDAD. Se calcula la
    distancia al esqueleto en todos los voxeles y se toma como radio el
    cuantil `rho` de esa distancia. Es el mismo truco que el umbral erf^-1 del
    spinodoide —el nivel se elige para que la fraccion solida sea la pedida—
    y deja la densidad EXACTA en vez de aproximada. Asi el ajuste busca en las
    mismas variables que con el spinodoide: densidad y escala.

LOS PARAMETROS, Y SU EQUIVALENTE EN EL SPINODOIDE
-------------------------------------------------
    rho            densidad relativa               =  rho
    celdas         celdas por lado del cubo unidad ~  numero de onda
                   (espaciado medio = 1/celdas). Fija la escala: mas celdas,
                   trabeculas mas finas y mas juntas.
    estiramiento   (ex, ey, ez), alargamiento      ~  thetas
                   relativo de las celdas por eje. Se normaliza a producto 1,
                   asi que solo cuenta la PROPORCION: (1,1,2) alarga en z sin
                   cambiar el volumen medio de celda.
    irregularidad  amplitud de la perturbacion de  (sin equivalente)
                   la rejilla, en fraccion del espaciado. 0.5 por omision.
    R              rotacion                        =  R
    seed           semilla                         =  seed

EL ESTIRAMIENTO SE APLICA A LA MALLA, NO A LOS PUNTOS
-----------------------------------------------------
Vafaeefar describen "una malla tetraedrica con el factor de estiramiento
deseado en cada direccion": se malla y la malla se estira. Eso fija la
CONECTIVIDAD antes de deformar, y todos los puntales se alargan de forma afin.

La primera version de este modulo hacia otra cosa que parece equivalente y no
lo es: estirar la nube de puntos y DESPUES triangularla. Delaunay trabaja en la
metrica isotropa, asi que ante puntos mas separados en z vuelve a elegir
tetraedros bien formados y la red deshace buena parte del estiramiento. Medido
sobre el VOI proximal de H4 (4.55 celdas, rho 0.28, 96^3, dos semillas):

    estiramiento z    1.0    1.5    2.0    2.5    3.0    4.0    5.0
    DA               1.04   1.12   1.15   1.12   1.12   1.16   1.15

El DA se satura en ~1.15 por mucho que se estire, y el VOI tiene 1.51: con
aquella construccion la familia no podia alcanzar la anisotropia del hueso, el
mismo defecto que la correccion G1 arreglo en el spinodoide. Triangulando la
rejilla sin estirar y estirando los vertices, el esqueleto es exactamente la
imagen afin del isotropo (el centroide de un tetraedro y el de una cara se
conservan por aplicaciones afines), y el DA ya responde. Mismo caso:

    estiramiento z    1.0    1.25   1.5    2.0    2.5    3.0    4.0
    DA               1.04   1.09   1.18   1.38   1.56   1.69   1.85

EL DOMINIO ES UN RECORTE, NO UNA CAJA CERRADA
---------------------------------------------
Los puntos se generan en una region MAYOR que el cubo —la esfera que lo
circunscribe mas dos celdas— y el cubo se recorta de ahi. Asi el borde del
cubo corta puntales como el bisturi corta trabeculas en un VOI, en vez de
fabricar una cascara artificial alrededor. Tambien permite rotar: la red se
genera en su marco, se gira, y el cubo ve siempre la misma red girada con la
misma semilla.

Convencion de ejes y rejilla, identicas a `grf.campo_grf`: dimension 1 = X y
`linspace(0, domain_size, resolution)` en cada eje.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy.spatial import Delaunay, cKDTree

# Muestras por lado de voxel al discretizar los puntales para el arbol k-d. La
# distancia final no es la de la muestra sino la EXACTA al segmento de esa
# muestra (ver `_distancia_esqueleto`); el muestreo solo decide cual es el
# segmento mas proximo, y con tres por voxel se equivoca en casos de empate
# cuya diferencia de distancia es menor que (h/6)^2 / (2 r).
_MUESTRAS_POR_VOXEL = 3

_CARAS = list(combinations(range(4), 3))


def esqueleto(celdas, estiramiento=(1.0, 1.0, 1.0), irregularidad=0.5,
              R=None, seed=None, domain_size=1.0, margen_celdas=2.0):
    """Red dual de una triangulacion de Delaunay que cubre el cubo.

    Devuelve un dict con
        nodos     (m, 3)  centroides de los tetraedros conservados
        caras     (k, 3)  centroides de las caras usadas (compartidas o no)
        seg_a     (4m, 3) extremo nodo de cada semipuntal
        seg_b     (4m, 3) extremo cara de cada semipuntal
        seg_nodo  (4m,)   indice del nodo de cada semipuntal
        seg_cara  (4m,)   indice de la cara de cada semipuntal
        celdas, estiramiento (normalizado), espaciado, irregularidad

    Solo se conservan los tetraedros cuyo centroide cae dentro del cubo
    ampliado en `margen_celdas` espaciados: los de fuera no pueden aportar
    ningun puntal que llegue al cubo con un radio razonable, y los de la
    envolvente de la nube son laminas degeneradas.
    """
    celdas = float(celdas)
    if celdas <= 0:
        raise ValueError(f"celdas debe ser > 0; se recibio {celdas}")
    irregularidad = float(irregularidad)
    if not (0.0 < irregularidad <= 1.0):
        raise ValueError("irregularidad debe estar en (0, 1]: con 0 la rejilla "
                         "es degenerada para Delaunay")
    e = np.asarray(estiramiento, float).ravel()
    if e.size != 3 or np.any(e <= 0):
        raise ValueError("estiramiento debe ser tres valores positivos")
    e = e / np.prod(e) ** (1.0 / 3.0)
    R = np.eye(3) if R is None else np.asarray(R, float)

    L = float(domain_size)
    d = L / celdas
    rng = np.random.default_rng(seed)

    # Radio que hay que cubrir en el marco de la red: la esfera circunscrita
    # al cubo mas el margen. Antes de estirar, cada eje necesita ese radio
    # dividido por su estiramiento.
    radio = np.sqrt(3.0) / 2.0 * L + margen_celdas * d * float(e.max())
    n_eje = np.ceil(radio / (e * d)).astype(int)
    ejes = [np.arange(-k, k + 1, dtype=float) for k in n_eje]
    G = np.stack(np.meshgrid(*ejes, indexing="ij"), axis=-1).reshape(-1, 3)
    G = G + rng.uniform(-0.5, 0.5, G.shape) * irregularidad

    # Se triangula la rejilla ISOTROPA y despues se estiran los vertices: la
    # conectividad se decide antes de deformar. Ver "EL ESTIRAMIENTO SE APLICA
    # A LA MALLA" en la cabecera: triangular la nube ya estirada satura el DA.
    tri = Delaunay(G)
    X = (G * d * e) @ R.T + 0.5 * L
    T = tri.simplices
    V = X[T]                                   # (m, 4, 3)
    C = V.mean(axis=1)

    tope = margen_celdas * d * float(e.max())
    dentro = np.all((C >= -tope) & (C <= L + tope), axis=1)
    # Descarta laminas: volumen casi nulo frente al de una celda media.
    vol = np.abs(np.einsum("ij,ij->i", V[:, 1] - V[:, 0],
                           np.cross(V[:, 2] - V[:, 0], V[:, 3] - V[:, 0]))) / 6.0
    dentro &= vol > 1e-6 * d ** 3
    T, V, C = T[dentro], V[dentro], C[dentro]
    m = len(T)

    # Caras: una clave por terna ordenada de vertices, para que las dos mitades
    # de un puntal compartan EXACTAMENTE el mismo punto de cara.
    ternas = np.sort(np.concatenate([T[:, list(c)] for c in _CARAS]), axis=1)
    claves, inv = np.unique(ternas, axis=0, return_inverse=True)
    inv = inv.ravel()
    F = X[claves].mean(axis=1)

    # `ternas` concatena primero la cara 0 de todos los tetraedros, luego la
    # cara 1, etc.: el nodo de cada semipuntal sigue ese mismo orden.
    seg_nodo = np.tile(np.arange(m), 4)
    seg_cara = inv
    return {
        "nodos": C, "caras": F,
        "seg_a": C[seg_nodo], "seg_b": F[seg_cara],
        "seg_nodo": seg_nodo, "seg_cara": seg_cara,
        "celdas": celdas, "estiramiento": e, "espaciado": d,
        "irregularidad": irregularidad,
    }


def _distancia_esqueleto(pts, A, B, h):
    """Distancia exacta de cada punto al semipuntal mas proximo.

    Se muestrea cada segmento con paso <= h/3, se busca la muestra mas proxima
    con un arbol k-d y se calcula la distancia EXACTA al segmento de esa
    muestra. Devuelve (distancia, indice de segmento).
    """
    AB = B - A
    largo = np.linalg.norm(AB, axis=1)
    n_s = max(2, int(np.ceil(largo.max() / (h / _MUESTRAS_POR_VOXEL))) + 1)
    t = np.linspace(0.0, 1.0, n_s)
    M = (A[:, None, :] + t[None, :, None] * AB[:, None, :]).reshape(-1, 3)
    _, j = cKDTree(M).query(pts, k=1, workers=-1)
    s = j // n_s

    a, ab = A[s], AB[s]
    den = np.maximum((ab * ab).sum(1), 1e-300)
    u = np.clip(((pts - a) * ab).sum(1) / den, 0.0, 1.0)
    dist = np.linalg.norm(pts - (a + u[:, None] * ab), axis=1)
    return dist, s


def campo_dual(resolution, celdas, estiramiento=(1.0, 1.0, 1.0),
               irregularidad=0.5, R=None, seed=None, domain_size=1.0):
    """Distancia de cada voxel al esqueleto, forma (res, res, res), dim 1 = X.

    Es el analogo del GRF: el campo continuo cuyo conjunto de nivel da la
    mascara. Devuelve (D, esqueleto).
    """
    res = int(resolution)
    L = float(domain_size)
    esq = esqueleto(celdas, estiramiento, irregularidad, R=R, seed=seed,
                    domain_size=L)
    eje = np.linspace(0.0, L, res)
    X, Y, Z = np.meshgrid(eje, eje, eje, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    # Solo cuentan los semipuntales que pueden llegar al cubo. El margen de un
    # espaciado basta: un radio mayor que eso ya no es una red de barras.
    tope = esq["espaciado"] * float(esq["estiramiento"].max())
    A, B = esq["seg_a"], esq["seg_b"]
    cerca = np.all((np.minimum(A, B) <= L + tope)
                   & (np.maximum(A, B) >= -tope), axis=1)
    h = L / max(res - 1, 1)
    D, _ = _distancia_esqueleto(pts, A[cerca], B[cerca], h)
    return D.reshape(res, res, res), esq


def generar_dual_lattice(resolution, celdas, rho, estiramiento=(1.0, 1.0, 1.0),
                         irregularidad=0.5, R=None, seed=None, domain_size=1.0,
                         mayor_componente=False):
    """Mascara booleana solida (True = hueso), forma (res, res, res), dim 1 = X.

    Fase solida = {D <= r}, con r el cuantil `rho` de la distancia al
    esqueleto: se toman los k = round(rho * res^3) voxeles mas proximos, de
    modo que la densidad obtenida es la pedida salvo empates exactos de
    distancia.

    Misma firma de retorno que `grf.generar_mascara`: (BW, campo, info).
    """
    if not (0.0 < rho < 1.0):
        raise ValueError(f"rho debe estar en (0,1); se recibio {rho}")
    D, esq = campo_dual(resolution, celdas, estiramiento, irregularidad, R=R,
                        seed=seed, domain_size=domain_size)
    plano = D.ravel()
    k = int(np.clip(round(float(rho) * plano.size), 1, plano.size))
    r = float(np.partition(plano, k - 1)[k - 1])
    BW = D <= r

    if mayor_componente:
        from .morphometry import mayor_componente_6
        BW = mayor_componente_6(BW)

    info = {
        "familia": "dual-lattice",
        "radio_puntal": r,
        "espesor_puntal": 2.0 * r,
        "rho_objetivo": float(rho),
        "rho_obtenida": float(BW.mean()),
        "celdas": float(celdas),
        "estiramiento": [float(x) for x in esq["estiramiento"]],
        "irregularidad": float(irregularidad),
        "n_nodos": int(len(esq["nodos"])),
        "n_caras": int(len(esq["caras"])),
    }
    return BW, D, info
