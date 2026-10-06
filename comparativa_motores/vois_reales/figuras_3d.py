"""
figuras_3d.py: Figuras tridimensionales y de campos del informe de VOIs reales.

    cd comparativa_motores/vois_reales
    python figuras_3d.py CARPETA_DE_VOIS

Superficies por marching cubes (scikit-image) dibujadas con matplotlib, sin
OpenGL. Los mapas y perfiles de tension leen resultados/campos.{json,npz}
(`python estudio.py CARPETA campos`). Solo PNG: con cientos de miles de
triangulos, el SVG no es practicable.
"""

import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI / "informe_vois"))
sys.path.insert(0, str(AQUI.parent.parent))
from estilo_figuras import (AZUL, GRIS, NARANJA, TINTA, TINTA2,   # noqa: E402
                            aplicar, coma)
import matplotlib.pyplot as plt                                    # noqa: E402
from matplotlib.colors import to_rgb                               # noqa: E402
from matplotlib.lines import Line2D                                # noqa: E402
from matplotlib.ticker import FuncFormatter                        # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection            # noqa: E402
from scipy.ndimage import gaussian_filter                          # noqa: E402
from skimage.measure import marching_cubes                         # noqa: E402

from spinpy.elastic import remuestrear_bw                          # noqa: E402
from spinpy.io import leer_voi                                     # noqa: E402

RES = AQUI / "resultados"
FIGS = AQUI / "informe_vois" / "figs"
VOIS = [("VOI_C1.mat", "C1", "porcino"),
        ("VOI_C2.mat", "C2", "porcino"),
        ("VOI_proximal_cubico.vtk", "Prox. cúbico", "equino"),
        ("VOI_proximal_PCAaligned.vtk", "Prox. PCA", "equino"),
        ("VOI_medio_PCAaligned.vtk", "Medio PCA", "equino")]
HUESO, GRIS_CLARO = "#d9d4c7", "#cfccc4"
LUZ = np.array([-0.45, -0.65, 0.62])
LUZ /= np.linalg.norm(LUZ)
aplicar()


def guardar_png(fig, nombre, dpi=260):
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / f"{nombre}.png", dpi=dpi)
    plt.close(fig)


def superficie(M, h, suavizado=0.8):
    """Mascara -> (vertices en mm, caras). Se rellena con un voxel de poro
    para cerrar la superficie en las caras del cubo (trabeculas cortadas)."""
    V = np.pad(M.astype(float), 1)
    if suavizado:
        V = gaussian_filter(V, suavizado)
    if V.max() < 0.5:
        return np.zeros((0, 3)), np.zeros((0, 3), int)
    v, f, _, _ = marching_cubes(V, 0.5, spacing=(h, h, h))
    return v - h + 0.5 * h, f


def sombrear(v, f, color):
    n = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-30
    i = 0.42 + 0.58 * np.abs(n @ LUZ)
    base = np.broadcast_to(np.asarray(color, float), (len(f), 3))
    return np.clip(base * i[:, None], 0, 1)


def dibujar(ax, v, f, colores, alfa=1.0):
    if len(f) == 0:
        return
    pc = Poly3DCollection(v[f], facecolors=colores, linewidths=0,
                          alpha=alfa)
    ax.add_collection3d(pc)


def dibujar_varias(ax, partes, h):
    """Varias mascaras en UNA coleccion: matplotlib ordena por profundidad
    los poligonos de una coleccion, pero las colecciones entre si solo como
    bloques, y una tapaba a la otra."""
    V, F, C, o = [], [], [], 0
    for M, color in partes:
        v, f = superficie(M, h)
        V.append(v)
        F.append(f + o)
        C.append(sombrear(v, f, to_rgb(color)))
        o += len(v)
    v, f = np.vstack(V), np.vstack(F)
    dibujar(ax, v, f, np.vstack(C))


def caja(ax, lo, hi, color, lw=1.0, ls="-"):
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    c = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1])
                  for z in (lo[2], hi[2])])
    for a in range(8):
        for b in range(a + 1, 8):
            if np.sum(c[a] != c[b]) == 1:
                ax.plot(*zip(c[a], c[b]), color=color, lw=lw, ls=ls)


def ejes(ax, L, elev=22, azim=-58):
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_zlim(0, L)
    ax.set_box_aspect((1, 1, 1), zoom=1.1)
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()


# ---------------------------------------------------------------------------
# Galeria de los cinco VOIs, a escala
# ---------------------------------------------------------------------------

def fig_galeria(carpeta):
    mor = json.loads((RES / "morfometria.json").read_text())
    fig = plt.figure(figsize=(6.8, 1.95))
    Lmax = 5.0
    for k, (archivo, nombre, especie) in enumerate(VOIS):
        BW, sp = leer_voi(Path(carpeta) / archivo)
        L = BW.shape[0] * sp[0]
        B, spr = remuestrear_bw(BW, sp, 80)
        v, f = superficie(B, spr[0])
        ax = fig.add_subplot(1, 5, k + 1, projection="3d")
        dibujar(ax, v, f, sombrear(v, f, to_rgb(HUESO)))
        caja(ax, (0, 0, 0), (L, L, L), GRIS, lw=0.6)
        ejes(ax, Lmax)
        ax.set_title(f"{nombre}\n{especie}, {coma(round(L, 1))} mm, "
                     f"BV/TV {coma(round(mor[nombre]['BVTV'], 2))}",
                     fontsize=7.3, loc="center", pad=-2)
        print("galeria", nombre, len(f), "triangulos", flush=True)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=0.86, wspace=0)
    guardar_png(fig, "f_vois_3d")


# ---------------------------------------------------------------------------
# Construccion de la referencia embebida sobre un VOI real
# ---------------------------------------------------------------------------

def _caras_corte(v, f, L, h, ejes_=(0, 1)):
    """Triangulos sobre los planos del cubo (las secciones de trabeculas
    cortadas) en los ejes indicados."""
    c = v[f].mean(axis=1)
    tol = 0.55 * h
    sel = np.zeros(len(f), bool)
    for e in ejes_:
        sel |= (c[:, e] < tol) | (c[:, e] > L - tol)
    return sel


def fig_construccion(carpeta, nombre="Prox. cúbico", N=96, n=64):
    archivo = next(a for a, nm, _ in VOIS if nm == nombre)
    BW, sp = leer_voi(Path(carpeta) / archivo)
    B, sp = remuestrear_bw(BW, sp, N)
    h = float(sp[0])
    m = (N - n) // 2
    LB, Lc = N * h, n * h
    k = int(round(0.625 / h))
    fig = plt.figure(figsize=(6.8, 2.95))

    # (a) bloque con un cuarto frontal retirado: el cubo de ensayo en azul
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    I = np.arange(N)
    cubo = np.zeros_like(B)
    cubo[m:m + n, m:m + n, m:m + n] = True
    corte = (I[:, None, None] >= N // 2) & (I[None, :, None] < N // 2)
    corte = np.broadcast_to(corte, B.shape)
    fuera = B & ~cubo & ~corte
    dentro = B & cubo & ~corte
    dibujar_varias(ax, ((fuera, GRIS_CLARO), (dentro, "#8fbdf0")), h)
    caja(ax, (0, 0, 0), (LB, LB, LB), GRIS, lw=0.6)
    caja(ax, (m * h,) * 3, ((m + n) * h,) * 3, NARANJA, lw=1.0)
    ejes(ax, LB)
    ax.set_title("(a) Bloque: VOI completo a 96³", fontsize=8)

    # (b) cubo de ensayo aislado: secciones cortadas y plato
    ax = fig.add_subplot(1, 3, 2, projection="3d")
    Bi = B[m:m + n, m:m + n, m:m + n]
    v, f = superficie(Bi, h)
    col = sombrear(v, f, to_rgb("#8fbdf0"))
    lat = _caras_corte(v, f, Lc, h, (0, 1))
    techo = _caras_corte(v, f, Lc, h, (2,)) & (v[f].mean(axis=1)[:, 2] > Lc / 2)
    col[lat] = to_rgb(NARANJA)
    col[techo] = to_rgb("#a33f16")
    dibujar(ax, v, f, col)
    caja(ax, (0, 0, 0), (Lc, Lc, Lc), NARANJA, lw=1.0)
    z = Lc * 1.04
    plato = [[(0, 0, z), (Lc, 0, z), (Lc, Lc, z), (0, Lc, z)]]
    ax.add_collection3d(Poly3DCollection(plato, facecolors="#4a4945",
                                         alpha=0.55, linewidths=0))
    for x in (0.25, 0.5, 0.75):
        ax.quiver(x * Lc, 0.5 * Lc, z + 0.32 * Lc, 0, 0, -0.26 * Lc,
                  color=TINTA, lw=1.0, arrow_length_ratio=0.35)
    ejes(ax, Lc * 1.1)
    ax.set_title("(b) Cubo de ensayo aislado (64³)", fontsize=8)

    # (c) nucleo y franjas del extensometro dentro del cubo
    ax = fig.add_subplot(1, 3, 3, projection="3d")
    i = np.arange(n)
    nuc = np.zeros_like(Bi)
    nuc[k:n - k, k:n - k, k:n - k] = True
    s = int(round(0.25 / h))
    franja = np.zeros_like(Bi)
    franja[k:n - k, k:n - k, k:k + s] = True
    franja[k:n - k, k:n - k, n - k - s:n - k] = True
    corte = (i[:, None, None] >= n // 2) & (i[None, :, None] < n // 2)
    corte = np.broadcast_to(corte, Bi.shape)
    dibujar_varias(ax, ((Bi & ~nuc & ~corte, GRIS_CLARO),
                        (Bi & nuc & ~franja & ~corte, "#8fbdf0"),
                        (Bi & franja & ~corte, "#1baf7a")), h)
    caja(ax, (0, 0, 0), (Lc, Lc, Lc), NARANJA, lw=1.0)
    caja(ax, (k * h,) * 3, ((n - k) * h,) * 3, AZUL, lw=1.0, ls="--")
    ejes(ax, Lc)
    ax.set_title("(c) Núcleo a 0,625 mm y extensómetro", fontsize=8)

    fig.legend(handles=[
        Line2D([], [], ls="", marker="s", ms=7, color=GRIS_CLARO,
               label="hueso circundante"),
        Line2D([], [], ls="", marker="s", ms=7, color="#8fbdf0",
               label="hueso del cubo de ensayo"),
        Line2D([], [], ls="", marker="s", ms=7, color=NARANJA,
               label="secciones cortadas (caras laterales)"),
        Line2D([], [], ls="", marker="s", ms=7, color="#a33f16",
               label="secciones cargadas (techo)"),
        Line2D([], [], ls="", marker="s", ms=7, color="#1baf7a",
               label="franjas del extensómetro (0,25 mm)")],
        loc="lower center", ncol=3, fontsize=7, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0, right=1, bottom=0.17, top=0.92, wspace=0.02)
    guardar_png(fig, "f_construccion_3d")


# ---------------------------------------------------------------------------
# Campos de tension: corte y perfiles
# ---------------------------------------------------------------------------

CASOS = (("traccion", "cubo aislado, tracción (hasta la V2.0.2)", NARANJA),
         ("plato", "cubo aislado, plato rígido (F1)", AZUL),
         ("referencia", "mismo cubo dentro del bloque (referencia)", TINTA))


TIT_CORTE = ("(a) Tracción, cubo aislado", "(b) Plato rígido, cubo aislado",
             "(c) Referencia: dentro del bloque")


def fig_campos(nombre="Prox. cúbico"):
    """Corte vertical y = L/2 del cubo de ensayo: tension vertical
    normalizada en los tres ensayos (arriba) y diferencia de cada cubo
    aislado con la referencia (abajo)."""
    cam = json.loads((RES / "campos.json").read_text())
    if nombre not in cam:
        nombre = next(iter(cam))
    cortes = np.load(RES / "campos.npz")
    clave = nombre.replace(" ", "_").replace(".", "")
    h = cam[nombre]["h_mm"]
    Z = {c: cortes[f"{clave}|{c}|zz"].T for c, _, _ in CASOS}
    n = Z["referencia"].shape[0]
    L = n * h
    fig = plt.figure(figsize=(6.8, 4.75))
    cmap = plt.get_cmap("magma_r").copy()
    cmap.set_bad("#ffffff")
    cdif = plt.get_cmap("RdBu_r").copy()
    cdif.set_bad("#ffffff")
    w, x0 = 0.25, 0.005
    for j, (caso, _, _) in enumerate(CASOS):
        ax = fig.add_axes([x0 + j * 0.262, 0.53, w, 0.40])
        im = ax.imshow(np.ma.masked_invalid(Z[caso]), origin="lower",
                       cmap=cmap, vmin=0, vmax=8, extent=(0, L, 0, L),
                       interpolation="nearest")
        _marco(ax)
        ax.set_title(TIT_CORTE[j], fontsize=7.4, loc="left")
    for j, caso in enumerate(("traccion", "plato")):
        ax = fig.add_axes([x0 + j * 0.262, 0.04, w, 0.40])
        D = Z[caso] - Z["referencia"]
        imd = ax.imshow(np.ma.masked_invalid(D), origin="lower", cmap=cdif,
                        vmin=-3, vmax=3, extent=(0, L, 0, L),
                        interpolation="nearest")
        _marco(ax)
        ax.set_title(("(d) Tracción menos referencia",
                      "(e) Plato menos referencia")[j], fontsize=7.4,
                     loc="left")
    for y, im_, et in ((0.53, im, r"$-\sigma_{zz}\,/\,\bar{\sigma}$"),
                       (0.04, imd, r"diferencia de $-\sigma_{zz}\,/\,\bar{\sigma}$")):
        cax = fig.add_axes([0.80, y, 0.014, 0.40])
        cb = fig.colorbar(im_, cax=cax)
        cb.set_label(et, fontsize=7.5)
        cb.ax.yaxis.set_major_formatter(FuncFormatter(coma))
        cb.ax.tick_params(labelsize=7)
    fig.text(0.535, 0.24, f"{nombre}\ncorte vertical en y = L/2\n"
             f"lado del cubo: {coma(round(L, 2))} mm\n"
             "carga vertical (z)\n\nrojo: el cubo aislado carga\n"
             "más que en el bloque\nazul: carga menos", fontsize=7.2,
             color=TINTA2, va="center")
    guardar_png(fig, "f_campos_corte", dpi=300)


def _marco(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for s_ in ax.spines.values():
        s_.set_visible(True)
        s_.set_color(GRIS)


def fig_perfiles():
    """Perfiles (11) de los cubos aislados divididos por los de la
    referencia: la arquitectura, comun a los tres ensayos, se cancela."""
    cam = json.loads((RES / "campos.json").read_text())
    nombres = [nm for _, nm, _ in VOIS if nm in cam]
    fig, axs = plt.subplots(2, len(nombres), figsize=(6.8, 3.9),
                            sharey="row", squeeze=False)
    for j, nm in enumerate(nombres):
        d = cam[nm]
        ref = d["referencia"]
        for caso, et, col in CASOS[:2]:
            q = d[caso]
            lat = np.divide(q["lateral"], ref["lateral"])
            ver = np.divide(q["vertical"], ref["vertical"])
            axs[0, j].plot(q["d_mm"], lat, color=col, lw=1.3, marker="o",
                           ms=2.2)
            axs[1, j].plot(ver, np.array(q["z_mm"]) / (len(q["z_mm"])
                                                      * d["h_mm"]),
                           color=col, lw=1.2)
        axs[0, j].axhline(1, color=TINTA, lw=0.8)
        axs[1, j].axvline(1, color=TINTA, lw=0.8)
        axs[0, j].set_title(nm, fontsize=7.8)
        axs[0, j].set_xlabel("d a la cara lateral (mm)", fontsize=7)
        axs[1, j].set_xlabel("von Mises / referencia", fontsize=7)
        axs[1, j].set_xlim(0.4, 2.6)
        for ax in axs[:, j]:
            ax.tick_params(labelsize=6.5)
    axs[0, 0].set_ylabel(r"$\sigma_{zz}$ de la capa /" "\nreferencia",
                         fontsize=7.5)
    axs[0, 0].set_ylim(0.4, 1.6)
    axs[1, 0].set_ylabel("altura z / H", fontsize=7.5)
    fig.legend(handles=[Line2D([], [], color=c, lw=1.6, label=e)
                        for _, e, c in CASOS[:2]] +
               [Line2D([], [], color=TINTA, lw=0.8,
                       label="igual que el cubo dentro del bloque")],
               loc="upper center", ncol=3, fontsize=7,
               bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.93), h_pad=1.2)
    for ax in axs.ravel():
        for eje in (ax.xaxis, ax.yaxis):
            eje.set_major_formatter(FuncFormatter(coma))
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "f_perfiles.png", dpi=300)
    fig.savefig(FIGS / "f_perfiles.svg")
    plt.close(fig)


if __name__ == "__main__":
    carpeta = sys.argv[1] if len(sys.argv) > 1 else None
    que = sys.argv[2:] or ["galeria", "construccion", "campos"]
    if carpeta and "galeria" in que:
        fig_galeria(carpeta)
    if carpeta and "construccion" in que:
        fig_construccion(carpeta)
    if "campos" in que and (RES / "campos.json").exists():
        fig_campos()
        fig_perfiles()
    print("figuras en", FIGS)
