"""grises.py: proyeccion de la malla de ladrillos sobre la imagen en GRISES.

La variante «proyI» de INFORME.md necesita un campo continuo cuya
isosuperficie sea la frontera del hueso. En un VOI de micro-CT ese campo es la
propia reconstruccion en grises, interpolada trilinealmente, y la frontera es
la isosuperficie al MISMO umbral con el que se segmento la mascara. Este
modulo reune las tres piezas que faltaban:

  leer_gris          la pila del escaner SIN binarizar (mismo lector y mismo
                     filtro de rebanadas que `voi.leer_pila_tiff`).
  extraer_cubo_gris  el cubo alineado por PCA, con los MISMOS puntos que
                     `voi.extraer_cubo`: por vecino mas proximo reproduce la
                     mascara bit a bit (se comprueba) y por interpolacion
                     trilineal da el campo para proyectar.
  campo_para         (f, grad, T) en coordenadas de la malla remuestreada a n,
                     con la misma correspondencia de indices que
                     `elastic.remuestrear_bw` (vecino mas proximo sobre
                     linspace), y un filtro gaussiano opcional.

`simular_microct` fabrica una reconstruccion sintetica a partir de una
geometria continua (volumen parcial, desenfoque del sistema y ruido) para
validar el flujo contra una superficie conocida antes de usar hueso real.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from spinpy import voi                                    # noqa: E402


# ---------------------------------------------------------------------------
# Lectura sin binarizar y cubo en grises
# ---------------------------------------------------------------------------

def leer_gris(origen, patron="*.tif"):
    """Volumen crudo (nx, ny, nz), dimension 1 = X, y las rebanadas apartadas."""
    p = Path(origen)
    if p.is_file():
        return voi._leer_pila_multipagina(p), []
    rutas, forma, apartados = voi._separar_rebanadas(
        voi.listar_rebanadas(p, patron))
    ny, nx = forma
    im0 = voi._leer_imagen(rutas[0])
    vol = np.zeros((nx, ny, len(rutas)), dtype=im0.dtype)
    for k, r in enumerate(rutas):
        vol[:, :, k] = (voi._leer_imagen(r) if k else im0).T
    return vol, apartados


def _puntos_cubo(spacing, centro_pca, ejes, centro, lado_vox):
    """Los puntos de `voi.extraer_cubo`, en indices CONTINUOS de voxel."""
    h = (lado_vox - 1) // 2
    r = np.arange(lado_vox) - h
    dx, dy, dz = np.meshgrid(r * spacing[0], r * spacing[1], r * spacing[2],
                             indexing="ij")
    delta = np.stack([dx.ravel(), dy.ravel(), dz.ravel()], axis=1)
    pts_mm = (delta + np.asarray(centro_pca, float)) @ ejes.T + centro
    return pts_mm / spacing[None, :]


def extraer_cubo_gris(vol, spacing, centro_pca, ejes, centro, lado_mm=5.0,
                      lado_vox=None, orden=1):
    """Cubo en grises con la geometria de `voi.extraer_cubo`.

    orden 0 da el vecino mas proximo (los mismos voxeles que la mascara);
    orden 1, interpolacion trilineal (fuera del volumen, el valor del borde).
    """
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    if lado_vox is None:
        lado_vox = int(round(lado_mm / spacing[0]))
    if lado_vox % 2 == 0:
        lado_vox += 1
    q = _puntos_cubo(spacing, centro_pca, ejes, centro, lado_vox)
    # Solo se lee la caja que rodea al cubo: la pila entera en float32
    # ocuparia varias veces lo que ocupa en enteros.
    idx0 = np.rint(q).astype(np.int64)
    dentro = np.all((idx0 >= 0) & (idx0 < np.array(vol.shape)), axis=1)
    lo = np.clip(np.floor(q.min(0)).astype(int) - 2, 0, None)
    hi = np.minimum(np.ceil(q.max(0)).astype(int) + 3, np.array(vol.shape))
    vol = vol[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
    q = q - lo[None, :]
    if orden == 0:
        # np.rint, como extraer_cubo (map_coordinates redondea distinto en .5)
        # Fuera del volumen, fondo por debajo de cualquier umbral: lo mismo
        # que hace extraer_cubo con la mascara.
        idx = np.clip(idx0 - lo[None, :], 0, np.array(vol.shape) - 1)
        vals = vol[idx[:, 0], idx[:, 1], idx[:, 2]].astype(np.float32)
        vals[~dentro] = float(vol.min()) - 1.0
    else:
        vals = ndimage.map_coordinates(vol.astype(np.float32), q.T, order=1,
                                       mode="nearest").astype(np.float32)
    cubo = vals.reshape(lado_vox, lado_vox, lado_vox)
    return np.transpose(cubo, (1, 0, 2))


# ---------------------------------------------------------------------------
# Campo para la proyeccion
# ---------------------------------------------------------------------------

def campo_para(G, T, N_malla, spc_malla, sigma=0.0):
    """(f, grad, T) sobre puntos en mm de la malla de lado `N_malla`.

    `elastic.remuestrear_bw` toma el voxel nativo round(j (N-1)/(n-1)) como
    voxel j de la malla, cuyo centro esta en (j + 1/2) spc. La correspondencia
    afin que respeta eso es q = (x / spc - 1/2) (N-1)/(n-1), con q el indice
    continuo en G. Con n = N es q = x / h - 1/2, el centro de cada voxel.

    G: cubo en grises a resolucion nativa (N^3), trilineal. sigma: filtro
    gaussiano previo, en voxeles NATIVOS (0 = sin filtro).
    """
    G = np.asarray(G, np.float32)
    if sigma and sigma > 0:
        G = ndimage.gaussian_filter(G, sigma, mode="nearest")
    N = np.array(G.shape, float)
    n = np.array(N_malla, float) if np.ndim(N_malla) else np.full(3, float(N_malla))
    spc = np.asarray(spc_malla, float) * np.ones(3)
    esc = (N - 1) / np.maximum(n - 1, 1)
    dG = [np.gradient(G, axis=k).astype(np.float32) for k in range(3)]

    def q(P):
        return ((P / spc[None, :]) - 0.5) * esc[None, :]

    def f(P):
        return ndimage.map_coordinates(G, q(P).T, order=1, mode="nearest")

    def gr(P):
        Q = q(P).T
        return np.stack([ndimage.map_coordinates(d, Q, order=1, mode="nearest")
                         for d in dG], axis=1) * (esc / spc)[None, :]

    return f, gr, float(T)


# ---------------------------------------------------------------------------
# micro-CT sintetico
# ---------------------------------------------------------------------------

def simular_microct(indicadora, N, sub=3, psf=0.7, snr=8.0, semilla=0,
                    fondo=0.0, hueso=1.0):
    """Reconstruccion sintetica de lado N de una geometria continua.

    indicadora(P) -> bool sobre puntos en [0, 1)^3. Volumen parcial: media de
    sub^3 muestras por voxel. Desenfoque: gaussiana de sigma `psf` voxeles.
    Ruido: gaussiano con desviacion (hueso - fondo) / snr. Devuelve float32.
    """
    rng = np.random.default_rng(semilla)
    o = (np.arange(sub) + 0.5) / sub
    c = np.arange(N)
    acc = np.zeros((N, N, N), np.float32)
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    for a in o:
        for b in o:
            for d in o:
                P = np.stack([(X + a) / N, (Y + b) / N, (Z + d) / N], -1)
                acc += indicadora(P.reshape(-1, 3)).reshape(N, N, N)
    acc /= sub ** 3
    img = fondo + (hueso - fondo) * acc
    if psf and psf > 0:
        img = ndimage.gaussian_filter(img, psf, mode="nearest")
    if snr and snr > 0:
        img = img + rng.normal(0.0, (hueso - fondo) / snr, img.shape)
    return img.astype(np.float32)
