"""
estilo_figuras.py: Estilo matplotlib de los informes en formato de guía de campo.

    import sys; sys.path.insert(0, ".claude/skills/informe-guia-campo")
    from estilo_figuras import aplicar, guardar, AZUL, NARANJA, AQUA
    aplicar()
    ...
    ax._cat_y = True          # si el eje y lleva etiquetas de texto
    guardar(fig, carpeta_figs, "f_nombre")    # PNG 300 ppp y SVG

Fuentes IBM Plex Sans (con DejaVu Sans de respaldo para letras griegas),
paleta categórica en orden fijo, rejilla tenue y coma decimal en los ejes
lineales. El texto va en tinta, nunca en el color de la serie.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
from matplotlib import font_manager                              # noqa: E402
from matplotlib.ticker import FixedFormatter, FuncFormatter      # noqa: E402

AQUI = Path(__file__).resolve().parent

#: Paleta categórica validada: usar en este orden y no pasar de tres series
#: en dispersión (más series: facetas o «Otros»).
AZUL, NARANJA, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
#: Rampa secuencial azul (claro a oscuro) para resoluciones o magnitudes.
AZULES = ("#86b6ef", "#2a78d6", "#184f95")
GRIS, TINTA, TINTA2 = "#8a8984", "#0b0b0b", "#52514e"


def aplicar():
    for f in (AQUI / "fuentes").glob("PlexSans-*.ttf"):
        font_manager.fontManager.addfont(str(f))
    plt.rcParams.update({
        "font.family": ["IBM Plex Sans", "DejaVu Sans"], "font.size": 8.5,
        "axes.titlesize": 9, "axes.labelsize": 8.5, "axes.titleweight": 600,
        "axes.titlelocation": "left", "axes.spines.top": False,
        "axes.spines.right": False, "axes.edgecolor": "#8a8984",
        "axes.labelcolor": "#2b2a27", "xtick.color": TINTA2,
        "ytick.color": TINTA2, "axes.grid": True, "grid.color": "#e6e5e0",
        "grid.linewidth": 0.6, "legend.frameon": False, "figure.dpi": 150,
        "savefig.dpi": 300, "savefig.bbox": "tight", "lines.linewidth": 1.8,
        "svg.fonttype": "path"})


def coma(v, _pos=None):
    """Número con coma decimal y signo menos tipográfico."""
    return f"{v:g}".replace(".", ",").replace("-", "−")


def guardar(fig, carpeta, nombre):
    """Aplica la coma decimal a los ejes lineales numéricos y guarda PNG y SVG.

    Un eje y con etiquetas de texto (barras horizontales, filas) se marca con
    `ax._cat_y = True` para que no se reescriban sus etiquetas.
    """
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    for ax in fig.axes:
        for eje, escala in ((ax.xaxis, ax.get_xscale()),
                            (ax.yaxis, ax.get_yscale())):
            if eje is ax.yaxis and getattr(ax, "_cat_y", False):
                continue
            if escala == "linear" and not isinstance(
                    eje.get_major_formatter(), FixedFormatter):
                eje.set_major_formatter(FuncFormatter(coma))
    fig.savefig(carpeta / f"{nombre}.png")
    fig.savefig(carpeta / f"{nombre}.svg")
    plt.close(fig)
