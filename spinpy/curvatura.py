# -*- coding: utf-8 -*-
"""
curvatura.py — Curvaturas principales de la interfaz hueso-vacio.

QUE MIDE, Y POR QUE ES UNA METRICA APARTE
------------------------------------------
BV/TV, Tb.Th y DA describen CUANTO material hay, como de grueso y hacia donde
mira. Ninguna de las tres dice nada sobre la FORMA LOCAL de la superficie: una
red conectada de trabeculas y un monton de islas sueltas con el mismo volumen,
el mismo espesor y la misma orientacion son indistinguibles para ellas. La
curvatura si las distingue, porque una isla es convexa por todas partes
(k1 > 0 y k2 > 0) mientras que una red es mayoritariamente silla de montar
(k1 > 0 > k2).

Eso es exactamente el problema abierto de este proyecto: el spinodoide ajustado
al VOI acierta BV/TV, Tb.Th y DA y aun asi tiene el 13.8 % de su hueso en islas
aisladas. El perfil de curvaturas es un descriptor barato que si lo ve.

La metrica es la de

    Guo, Y., Sharma, S., & Kumar, S. (2024). Inverse designing surface
        curvatures by deep learning. Advanced Intelligent Systems, 6(6),
        2300789. https://doi.org/10.1002/aisy.202300789

cuya ecuacion (3) define la probabilidad de que la superficie presente el par
(k1, k2) como la fraccion de AREA con esas curvaturas — no la fraccion de
elementos. La diferencia no es cosmetica: una malla de marching cubes tiene
triangulos de areas muy distintas, y contar elementos pesaria mas las zonas
teseladas fino, que son justo las de curvatura alta. `perfil` pondera por area.

CONVENCION DE SIGNO — SE DECLARA AQUI Y NO SE TOCA
--------------------------------------------------
La normal apunta HACIA EL VACIO (fuera del hueso), y las curvaturas son los
valores propios del mapa de Weingarten dn restringido al plano tangente. Con
eso:

    bola de hueso de radio R        k1 = k2 = +1/R      (convexa, positiva)
    poro esferico de radio R        k1 = k2 = -1/R
    trabecula tipica                k1 > 0 > k2         (silla de montar)
    superficie minima               k1 = -k2            (H = 0)

Es la convencion de la literatura de hueso (curvatura positiva = solido
convexo) y la que hace que "H = 0" signifique superficie minima, que es la
afirmacion que el articulo discute a proposito de las TPMS. El signo se
verifica sobre una esfera en el bloque 09 de `tests/`; si alguien lo invierte,
esa prueba falla.

DOS CAMINOS, Y POR QUE HACEN FALTA LOS DOS
-------------------------------------------
1. `curvaturas_implicitas` — EXACTA. Si la superficie es el conjunto de nivel
   de un campo del que se conocen gradiente y hessiano (el GRF de un
   spinodoide, una superficie nodal periodica), la curvatura tiene formula
   cerrada y no hay error de discretizacion ninguno. Es la respuesta contra la
   que se valida lo demas.

2. `curvaturas_malla` — DISCRETA. Un VOI de micro-CT no tiene campo analitico:
   solo hay voxeles. Aqui se ajusta la segunda forma fundamental por triangulo
   (Rusinkiewicz, 2004) y se promedia a los vertices ponderando por area. Es lo
   unico aplicable al hueso real, y por eso hay que validarlo contra el camino
   1 antes de creerselo.

La unidad es 1/mm, coherente con el resto del proyecto (Tb.Th en mm). Una
curvatura medida sobre una mascara con `spacing` en mm sale en 1/mm sin
conversion adicional.

LO QUE NO HACE
--------------
No suaviza. Sobre una mascara BINARIA cruda, la superficie de marching cubes es
una escalera y su curvatura discreta es ruido con estructura: se miden las
esquinas de los voxeles, no la trabecula. Por eso la entrada natural es
`curvaturas_de_campo`, que mide sobre el campo CONTINUO — el GRF en un
spinodoide, o la mascara suavizada en un VOI—, que es lo que hace el articulo
al medir sobre su campo de fase.

Referencia del estimador discreto:
    Rusinkiewicz, S. (2004). Estimating curvatures and their derivatives on
        triangle meshes. En 3DPVT 2004 (pp. 486-493). IEEE.
        https://doi.org/10.1109/TDPVT.2004.1335277
"""
from __future__ import annotations

import numpy as np

__all__ = ["curvaturas_implicitas", "curvaturas_malla", "curvaturas_de_campo",
           "areas_baricentricas", "normales_desde_campo", "perfil", "resumen",
           "mediana_ponderada"]


# ---------------------------------------------------------------------------
# Camino exacto: conjunto de nivel de un campo con derivadas conocidas
# ---------------------------------------------------------------------------

def curvaturas_implicitas(grad, hess):
    """k1 >= k2 de la superficie {f = c} a partir de grad f y hess f.

    Con n = grad/|grad| (que apunta hacia donde el campo CRECE, es decir hacia
    el vacio, porque el solido es {f <= c}), el mapa de Weingarten restringido
    al plano tangente es

        S = P H P / |grad|        con P = I - n n^T

    y sus dos valores propios no nulos son k1 y k2. El tercero es exactamente
    cero por construccion — P anula la direccion normal— y se descarta POR
    MAGNITUD, no por posicion: en una silla perfecta los dos utiles son +a y -a
    y el cero queda EN MEDIO al ordenar. Quedarse con "los dos ultimos" o
    tirar "el de en medio" da curvaturas equivocadas justo en las sillas, que
    son la mayor parte de la superficie de una trabecula.

    grad : (n, 3)      hess : (n, 3, 3)     ->  (k1, k2), cada uno (n,)
    """
    g = np.asarray(grad, float)
    H = np.asarray(hess, float)
    ng = np.linalg.norm(g, axis=1)
    ok = ng > 0

    n = np.zeros_like(g)
    n[ok] = g[ok] / ng[ok, None]

    P = np.eye(3)[None, :, :] - n[:, :, None] * n[:, None, :]
    S = np.einsum("nij,njk,nkl->nil", P, H, P)
    with np.errstate(divide="ignore", invalid="ignore"):
        S = S / np.where(ok, ng, np.nan)[:, None, None]

    w = np.linalg.eigvalsh((S + np.swapaxes(S, 1, 2)) / 2.0)   # (n, 3), asc
    idx = np.argsort(np.abs(w), axis=1)[:, 1:]                 # los dos mayores
    k = np.take_along_axis(w, idx, axis=1)
    return np.max(k, axis=1), np.min(k, axis=1)


# ---------------------------------------------------------------------------
# Camino discreto: malla de triangulos
# ---------------------------------------------------------------------------

def areas_baricentricas(verts, faces):
    """Area asignada a cada vertice: un tercio de la de cada triangulo suyo.

    Es el peso de la ecuacion (3) del articulo. Se usa la reparticion
    baricentrica y no la de Voronoi porque la de Voronoi se vuelve NEGATIVA en
    triangulos obtusos, y marching cubes produce triangulos obtusos a montones;
    un peso negativo en un histograma de probabilidad no significa nada.
    """
    v = np.asarray(verts, float)
    f = np.asarray(faces, int)
    a = np.linalg.norm(np.cross(v[f[:, 1]] - v[f[:, 0]],
                                v[f[:, 2]] - v[f[:, 0]]), axis=1) / 2.0
    w = np.zeros(len(v))
    for c in range(3):
        np.add.at(w, f[:, c], a / 3.0)
    return w


def _normales_vertice(verts, faces):
    """Normales por vertice, promedio de las de cara ponderado por area.

    El producto vectorial sin normalizar YA lleva el area dentro, asi que
    sumarlo tal cual es la media ponderada; normalizar antes de sumar daria el
    mismo peso a un triangulo diminuto de una esquina que a uno grande.
    """
    v = np.asarray(verts, float)
    f = np.asarray(faces, int)
    cr = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    nv = np.zeros_like(v)
    for c in range(3):
        np.add.at(nv, f[:, c], cr)
    ln = np.linalg.norm(nv, axis=1)
    ln[ln == 0] = 1.0
    return nv / ln[:, None]


def normales_desde_campo(verts, campo, spacing, origen=0.0):
    """Normales orientadas HACIA EL VACIO usando el gradiente del campo.

    La orientacion no se puede dejar al orden de los triangulos: marching cubes
    la fija por convenio interno, y una malla con la normal invertida da todas
    las curvaturas con el signo cambiado. Eso convierte una silla en otra silla
    —o sea, pasa desapercibido en lo simetrico— y a la vez convierte una isla
    convexa en un poro y H en -H, que es justo lo que aqui se quiere medir.

    Se interpola el gradiente de `campo` por diferencias centradas en los
    vertices y se toma su direccion: el solido es {campo <= nivel}, luego el
    campo crece hacia el vacio.
    """
    from scipy import ndimage

    v = np.asarray(verts, float)
    campo = np.asarray(campo, float)
    sp = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if sp.size == 1:
        sp = np.repeat(sp, 3)
    org = np.atleast_1d(np.asarray(origen, float)).ravel()
    if org.size == 1:
        org = np.repeat(org, 3)

    idx = ((v - org) / sp).T                       # (3, n), en indices
    g = np.stack([ndimage.map_coordinates(
        np.gradient(campo, sp[e], axis=e), idx, order=1, mode="nearest")
        for e in range(3)], axis=1)
    ln = np.linalg.norm(g, axis=1)
    ln[ln == 0] = 1.0
    return g / ln[:, None]


def curvaturas_malla(verts, faces, normales=None):
    """k1 >= k2 por vertice, ajustando la segunda forma fundamental por cara.

    Metodo de Rusinkiewicz (2004): en cada triangulo se busca la matriz 2x2
    simetrica II que explica como cambia la normal a lo largo de sus tres
    aristas,

        II . (e . u, e . w)^T  =  (dn . u, dn . w)^T       para las 3 aristas,

    seis ecuaciones y tres incognitas que se resuelven por minimos cuadrados.
    El resultado se sube a un tensor 3x3 en el plano de la cara, se acumula en
    los vertices ponderando por area y se diagonaliza en el plano tangente del
    vertice.

    POR QUE ESTE Y NO EL LAPLACIANO DE COTANGENTES. El operador de cotangentes
    da la curvatura MEDIA barata y bien, y la gaussiana por defecto angular,
    pero NO da las curvaturas principales por separado: hay que sacarlas de una
    cuadratica cuyo discriminante se vuelve negativo en cuanto hay ruido, y
    entonces k1 y k2 son complejos. El perfil (k1, k2) del articulo necesita
    las dos curvaturas en TODA la superficie, sin huecos.

    normales : (n, 3) normales por vertice YA orientadas hacia el vacio. Si no
        se dan, se deducen del orden de los triangulos, cuyo signo depende de
        quien construyo la malla. Con campo disponible —siempre, en un
        spinodoide— usar `normales_desde_campo`.

    Devuelve (k1, k2, normales).
    """
    v = np.asarray(verts, float)
    f = np.asarray(faces, int)
    nv = (_normales_vertice(v, f) if normales is None
          else np.asarray(normales, float))

    p0, p1, p2 = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    n0, n1, n2 = nv[f[:, 0]], nv[f[:, 1]], nv[f[:, 2]]

    e = np.stack([p2 - p1, p0 - p2, p1 - p0], axis=1)          # (m, 3, 3)
    dn = np.stack([n2 - n1, n0 - n2, n1 - n0], axis=1)

    cr = np.cross(e[:, 0], e[:, 1])
    area = np.linalg.norm(cr, axis=1) / 2.0
    bueno = area > 0
    nf = np.zeros_like(cr)
    nf[bueno] = cr[bueno] / (2.0 * area[bueno, None])

    u = e[:, 0].copy()
    lu = np.linalg.norm(u, axis=1)
    bueno &= lu > 0
    u[bueno] /= lu[bueno, None]
    w = np.cross(nf, u)                                        # completa (u, w)

    # El sistema se monta a mano y se resuelve en bloque: son seis filas fijas
    # y tres incognitas (a, b, c), y llamar a lstsq una vez por triangulo
    # cuesta minutos en una malla de medio millon de caras.
    eu = np.einsum("mij,mj->mi", e, u)                         # (m, 3)
    ew = np.einsum("mij,mj->mi", e, w)
    du = np.einsum("mij,mj->mi", dn, u)
    dw = np.einsum("mij,mj->mi", dn, w)

    m = len(f)
    A = np.zeros((m, 6, 3))
    b = np.zeros((m, 6))
    for i in range(3):
        A[:, 2 * i, 0] = eu[:, i]
        A[:, 2 * i, 1] = ew[:, i]
        b[:, 2 * i] = du[:, i]
        A[:, 2 * i + 1, 1] = eu[:, i]
        A[:, 2 * i + 1, 2] = ew[:, i]
        b[:, 2 * i + 1] = dw[:, i]

    AtA = np.einsum("mki,mkj->mij", A, A)
    Atb = np.einsum("mki,mk->mi", A, b)
    # Regularizacion minima: en un triangulo degenerado AtA es singular y
    # `solve` levantaria LinAlgError para toda la malla por culpa de uno.
    AtA += 1e-12 * np.eye(3)[None]
    sol = np.linalg.solve(AtA, Atb[..., None])[..., 0]         # (m, 3) a, b, c

    T = (sol[:, 0, None, None] * (u[:, :, None] * u[:, None, :])
         + sol[:, 1, None, None] * (u[:, :, None] * w[:, None, :]
                                    + w[:, :, None] * u[:, None, :])
         + sol[:, 2, None, None] * (w[:, :, None] * w[:, None, :]))

    acc = np.zeros((len(v), 3, 3))
    peso = np.zeros(len(v))
    aw = np.where(bueno, area, 0.0)
    for c in range(3):
        np.add.at(acc, f[:, c], T * aw[:, None, None])
        np.add.at(peso, f[:, c], aw)
    peso[peso == 0] = 1.0
    acc /= peso[:, None, None]

    # Diagonalizacion en el plano tangente del vertice. Mismo criterio que en
    # el camino exacto: se descarta el valor propio nulo por magnitud.
    P = np.eye(3)[None] - nv[:, :, None] * nv[:, None, :]
    S = np.einsum("nij,njk,nkl->nil", P, acc, P)
    ev = np.linalg.eigvalsh((S + np.swapaxes(S, 1, 2)) / 2.0)
    idx = np.argsort(np.abs(ev), axis=1)[:, 1:]
    k = np.take_along_axis(ev, idx, axis=1)
    return np.max(k, axis=1), np.min(k, axis=1), nv


def curvaturas_de_campo(campo, spacing, nivel, exacta=None, origen=0.0,
                        margen=0):
    """Malla del conjunto de nivel y curvaturas sobre ella, en un solo paso.

    `campo` es CONTINUO (el GRF, no la mascara binaria): la curvatura de una
    escalera de voxeles no es la de la trabecula.

    `exacta`, si se da, es una funcion pts -> (grad, hess) con las derivadas
    analiticas en los vertices; entonces se calculan tambien las curvaturas
    exactas y se devuelven las dos series. Asi se valida el estimador discreto
    contra la respuesta cerrada sobre la MISMA malla, que es la unica forma de
    que la comparacion no mezcle el error del estimador con el del mallado.

    `margen` descarta los vertices a menos de esa distancia (en voxeles) de las
    caras del cubo. En el borde la malla se corta y las diferencias de normal
    quedan cojas, asi que la curvatura de esos vertices no es de la superficie
    sino del recorte. Es un descarte, no una correccion: se dice cuantos.

    Devuelve un diccionario con verts, faces, areas, normales, k1, k2, el
    numero de vertices descartados y, si procede, k1_exacta y k2_exacta.
    """
    from skimage import measure

    campo = np.asarray(campo, float)
    sp = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if sp.size == 1:
        sp = np.repeat(sp, 3)
    org = np.atleast_1d(np.asarray(origen, float)).ravel()
    if org.size == 1:
        org = np.repeat(org, 3)

    verts, faces, _, _ = measure.marching_cubes(campo, level=float(nivel),
                                                spacing=tuple(sp))
    nv = normales_desde_campo(verts, campo, sp)
    k1, k2, nv = curvaturas_malla(verts, faces, normales=nv)
    areas = areas_baricentricas(verts, faces)

    fuera = np.zeros(len(verts), dtype=bool)
    if margen > 0:
        lim = (np.array(campo.shape) - 1.0) * sp
        d = np.minimum(verts, lim[None, :] - verts)
        fuera = (d < margen * sp[None, :]).any(axis=1)
        areas = np.where(fuera, 0.0, areas)

    out = {"verts": verts + org[None, :], "faces": faces, "normales": nv,
           "areas": areas, "k1": k1, "k2": k2,
           "n_vertices": int(len(verts)), "n_descartados": int(fuera.sum())}
    if exacta is not None:
        g, H = exacta(verts + org[None, :])
        e1, e2 = curvaturas_implicitas(g, H)
        out["k1_exacta"], out["k2_exacta"] = e1, e2
    return out


# ---------------------------------------------------------------------------
# Perfil de curvaturas — ecuacion (3) del articulo
# ---------------------------------------------------------------------------

def perfil(k1, k2, areas, limite=None, nbins=200):
    """Histograma 2D de (k1, k2) ponderado por AREA y normalizado a 1.

    Devuelve (P, bordes_k1, bordes_k2) con P de forma (nbins, nbins). Es la
    ecuacion (3) del articulo, salvo que alli se serializa la mitad triangular
    para alimentar una red neuronal; aqui se devuelve la matriz entera porque
    lo que hacemos con ella es dibujarla.
    """
    k1 = np.asarray(k1, float)
    k2 = np.asarray(k2, float)
    w = np.asarray(areas, float)
    fin = np.isfinite(k1) & np.isfinite(k2) & np.isfinite(w) & (w > 0)
    if limite is None:
        v = np.abs(np.r_[k1[fin], k2[fin]])
        limite = float(np.percentile(v, 99)) if v.size else 1.0
        limite = limite if limite > 0 else 1.0
    bordes = np.linspace(-limite, limite, nbins + 1)
    P, bx, by = np.histogram2d(k1[fin], k2[fin], bins=[bordes, bordes],
                               weights=w[fin])
    s = P.sum()
    return (P / s if s > 0 else P), bx, by


def resumen(k1, k2, areas):
    """Cifras que resumen un perfil de curvaturas, ponderadas por area.

    H y K son la media y la gaussiana; `silla` es la fraccion de area con
    k1 > 0 > k2, que es la firma de una red conectada, y `convexa` la de area
    con las dos positivas, que es la de las islas sueltas.
    """
    k1 = np.asarray(k1, float)
    k2 = np.asarray(k2, float)
    w = np.asarray(areas, float)
    fin = np.isfinite(k1) & np.isfinite(k2) & np.isfinite(w) & (w > 0)
    k1, k2, w = k1[fin], k2[fin], w[fin]
    A = float(w.sum())
    if A <= 0:
        return {"area": 0.0}
    H = (k1 + k2) / 2.0
    K = k1 * k2
    Hm = float((w * H).sum() / A)
    return {
        "area": A,
        "H_medio": Hm,
        # La DISPERSION importa tanto como el centro: dos superficies con la
        # misma curvatura media pueden ser una lisa y la otra llena de
        # esquinas, y al comparar con un perfil publicado hay que poder decir
        # si lo que difiere es donde esta la mancha o como de ancha es.
        "H_desv": float(np.sqrt((w * (H - Hm) ** 2).sum() / A)),
        "H_mediana": mediana_ponderada(H, w),
        "K_medio": float((w * K).sum() / A),
        "K_integral": float((w * K).sum()),
        "k1_medio": float((w * k1).sum() / A),
        "k2_medio": float((w * k2).sum() / A),
        "silla": float(w[(k1 > 0) & (k2 < 0)].sum() / A),
        "convexa": float(w[(k1 > 0) & (k2 > 0)].sum() / A),
        "concava": float(w[(k1 < 0) & (k2 < 0)].sum() / A),
    }


def mediana_ponderada(x, w):
    """Mediana de x con pesos w. Publica porque el centro de un perfil de
    curvaturas es justo esto y lo piden desde fuera del modulo."""
    o = np.argsort(x)
    x, w = x[o], w[o]
    c = np.cumsum(w)
    return float(x[min(np.searchsorted(c, c[-1] / 2.0), len(x) - 1)])
