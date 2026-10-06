"""
figuras_vois.py: Figuras del informe de las correcciones en VOIs reales.

    cd comparativa_motores/vois_reales
    python figuras_vois.py [CARPETA_DE_VOIS]   # la carpeta, para la Figura 1

Lee resultados/{validacion,practico,morfometria}.json y escribe PNG y SVG en
informe_vois/figs/.
"""

import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI / "informe_vois"))
sys.path.insert(0, str(AQUI.parent.parent))
from estilo_figuras import (AQUA, AZUL, GRIS, NARANJA, TINTA,  # noqa: E402
                            TINTA2, aplicar, guardar)
import matplotlib.pyplot as plt                                  # noqa: E402
from matplotlib.lines import Line2D                              # noqa: E402

RES = AQUI / "resultados"
FIGS = AQUI / "informe_vois" / "figs"
val = json.loads((RES / "validacion.json").read_text())
pra = json.loads((RES / "practico.json").read_text())
mor = json.loads((RES / "morfometria.json").read_text())
VOIS = list(val)
aplicar()


def err(x, r):
    return 100 * (x / r - 1)


def _filas(ax):
    y = np.arange(len(VOIS))[::-1]
    ax.set_yticks(y)
    ax.set_yticklabels([f"{v} ({val[v]['especie']})" for v in VOIS])
    ax._cat_y = True
    ax.grid(axis="y", visible=False)
    ax.axvline(0, color=TINTA2, lw=0.9)
    ax.axvspan(-5, 5, color="#f0efec", lw=0, zorder=0)
    return y


def fig_validacion_E():
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    y = _filas(ax)
    for yi, v in zip(y, VOIS):
        d = val[v]
        ri = d["referencia"]["plato"]["interior"]["E"]
        m = d["metodos"]
        ax.plot(err(m["B0"]["E"], ri), yi + 0.16, "o", ms=6.5, color=NARANJA,
                mec="white", mew=0.8)
        ax.plot(err(m["F1"]["E"], ri), yi + 0.05, "o", ms=6.5, color=AZUL,
                mec="white", mew=0.8)
        if "E" in m["F12"]:
            rn = d["referencia"]["plato"]["nucleo"]["E"]
            ax.plot(err(m["F12"]["E"], rn), yi - 0.06, "o", ms=6.5,
                    color=AQUA, mec="white", mew=0.8)
        rx = d["referencia"]["plato"]["nucleo_exp"]["E"]
        ax.plot(err(m["F12_exp"]["E"], rx), yi - 0.17, "o", ms=6.5,
                mfc="white", mec=AQUA, mew=1.6)
    ax.legend(handles=[
        Line2D([], [], ls="", marker="o", color=NARANJA, label="publicado hasta la V2.0.2"),
        Line2D([], [], ls="", marker="o", color=AZUL, label="F1: plato rígido"),
        Line2D([], [], ls="", marker="o", color=AQUA,
               label="F1 + F2: plato y núcleo (app)"),
        Line2D([], [], ls="", marker="o", mfc="white", mec=AQUA, mew=1.6,
               label="plato y núcleo de 0,25 mm (exploratorio)")],
        loc="lower left", fontsize=7.2, ncol=2, bbox_to_anchor=(0, 1.01))
    ax.set_xlabel("error de E aparente frente al mismo hueso embebido (%)")
    ax.set_xlim(-40, 10)
    fig.tight_layout()
    guardar(fig, FIGS, "f_validacion_E")


def fig_validacion_p99():
    fig, ax = plt.subplots(figsize=(6.6, 2.6))
    y = _filas(ax)
    for yi, v in zip(y, VOIS):
        d = val[v]
        ri = d["referencia"]["plato"]["interior"]["p99"]
        m = d["metodos"]
        ax.plot(err(m["B0"]["p99"], ri), yi + 0.12, "o", ms=6.5,
                color=NARANJA, mec="white", mew=0.8)
        ax.plot(err(m["F1"]["p99"], ri), yi - 0.12, "o", ms=6.5, color=AZUL,
                mec="white", mew=0.8)
    ax.legend(handles=[
        Line2D([], [], ls="", marker="o", color=NARANJA, label="publicado hasta la V2.0.2"),
        Line2D([], [], ls="", marker="o", color=AZUL,
               label="p99 corregido (plato rígido)")],
        loc="lower right", fontsize=7.5)
    ax.set_xlabel("error del p99 de von Mises frente al hueso embebido (%)")
    fig.tight_layout()
    guardar(fig, FIGS, "f_validacion_p99")


def fig_practico():
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True)
    x = np.arange(len(VOIS))
    for ax, malla, tit in ((axs[0], "hex8", "(a) Ladrillos, 40³"),
                           (axs[1], "tet10", "(b) Malla suave, 40³")):
        hoy = [pra.get(f"{v}|{malla}", {}).get("E_app", np.nan) for v in VOIS]
        cor = [pra.get(f"{v}|{malla}", {}).get("E_corregido", np.nan)
               for v in VOIS]
        ax.bar(x - 0.19, hoy, width=0.36, color=NARANJA, label="E publicado hasta la V2.0.2")
        ax.bar(x + 0.19, cor, width=0.36, color=AZUL, label="E corregido")
        ax.set_xticks(x)
        ax.set_xticklabels(VOIS, rotation=25, ha="right", fontsize=7.5)
        ax._cat_x = True
        ax.set_title(tit)
        for xi, v in zip(x, VOIS):
            if "error" in pra.get(f"{v}|{malla}", {}):
                ax.text(xi, 150, "sin malla", rotation=90, ha="center",
                        va="bottom", fontsize=7, color=TINTA2)
    axs[0].set_ylabel("E aparente (MPa)")
    axs[0].legend(loc="upper left", fontsize=7.5)
    fig.tight_layout()
    guardar(fig, FIGS, "f_practico")


def fig_a3():
    fig, ax = plt.subplots(figsize=(6.6, 2.5))
    y = np.arange(len(VOIS))[::-1]
    for yi, v in zip(y, VOIS):
        h = pra.get(f"{v}|hex8", {})
        t = pra.get(f"{v}|tet10", {})
        if not h or not t or "E_app" not in t:
            ax.text(-1, yi, "malla suave no generada (tetgen)", va="center",
                    ha="right", fontsize=7.5, color=TINTA2)
            continue
        b = err(t["E_app"], h["E_app"])
        c = err(t["E_corregido"], h["E_corregido"])
        ax.barh(yi + 0.17, b, height=0.3, color=NARANJA,
                label="hasta la V2.0.2" if yi == y[0] else None)
        ax.barh(yi - 0.17, c, height=0.3, color=AZUL,
                label="corregido" if yi == y[0] else None)
        for val_, dy in ((b, 0.17), (c, -0.17)):
            ax.text(val_ + (1 if val_ >= 0 else -1), yi + dy,
                    f"{val_:+.0f} %".replace("-", "−"), va="center",
                    ha="left" if val_ >= 0 else "right", fontsize=7,
                    color=TINTA2)
    ax.axvline(0, color=TINTA2, lw=0.9)
    ax.axvline(-10, color=GRIS, lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(VOIS)
    ax._cat_y = True
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("E de la malla suave frente a la de ladrillos, ambas a 40³ (%)")
    ax.set_xlim(-55, 4)
    ax.legend(loc="lower left", fontsize=7.5)
    fig.tight_layout()
    guardar(fig, FIGS, "f_a3")


def fig_resolucion():
    fig, ax = plt.subplots(figsize=(6.6, 2.4))
    x = np.arange(len(VOIS))
    r40 = [mor[v]["TbTh"] / (val[v]["lado_voi_mm"] / 40) for v in VOIS]
    r96 = [mor[v]["TbTh"] / val[v]["h_mm"] for v in VOIS]
    ax.bar(x - 0.19, r40, width=0.36, color=NARANJA,
           label="40³ (ladrillos por omisión)")
    ax.bar(x + 0.19, r96, width=0.36, color=AZUL,
           label="96³ (validación embebida)")
    ax.axhline(4, color=TINTA, lw=1.1, ls="--")
    ax.text(2.5, 4.15, "4 elementos por espesor", fontsize=7.5,
            color=TINTA2, ha="center")
    ax.set_xticks(x)
    ax.set_xticklabels(VOIS)
    ax._cat_x = True
    ax.set_ylabel("Tb.Th / tamaño de vóxel")
    ax.legend(loc="upper left", fontsize=7.5)
    fig.tight_layout()
    guardar(fig, FIGS, "f_resolucion")


def fig_perfil_bvtv():
    """BV/TV por capa horizontal frente a la altura relativa."""
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.7), sharey=True)
    grupos = (("porcino", "(a) Porcinos, 3 mm"), ("equino", "(b) Equinos, 5 mm"))
    for ax, (esp, tit) in zip(axs, grupos):
        nombres = [v for v in VOIS if val[v]["especie"] == esp]
        for v, col in zip(nombres, (AZUL, NARANJA, AQUA)):
            p = np.asarray(mor[v]["perfil_z"])
            z = (np.arange(len(p)) + 0.5) / len(p)
            ax.plot(p, z, color=col, lw=1.4, label=v)
        ax.axhspan(16 / 96, 80 / 96, color="#f0efec", lw=0, zorder=0)
        ax.axhline(16 / 96, color=GRIS, lw=0.7, ls="--")
        ax.axhline(80 / 96, color=GRIS, lw=0.7, ls="--")
        ax.set_title(tit)
        ax.set_xlabel("BV/TV de la capa")
        ax.set_xlim(0, 1.02)
        ax.legend(loc="lower right", fontsize=7.3)
    axs[0].set_ylabel("altura relativa z / H")
    axs[0].text(0.02, 0.5, "cubo de\nensayo", fontsize=7, color=TINTA2,
                va="center")
    fig.tight_layout()
    guardar(fig, FIGS, "f_perfil_bvtv")


def fig_esquema_real(carpeta):
    """Corte de un VOI real con el bloque, el cubo de ensayo y el nucleo."""
    from spinpy.elastic import remuestrear_bw
    from spinpy.io import leer_voi
    v = "Prox. cúbico"
    BW, sp = leer_voi(Path(carpeta) / val[v]["archivo"])
    B, sp = remuestrear_bw(BW, sp, 96)
    h = sp[0]
    L = 96 * h
    corte = B[:, 48, :].T                     # plano x-z
    fig, ax = plt.subplots(figsize=(3.6, 3.6))
    ax.imshow(corte, origin="lower", cmap="Greys", extent=(0, L, 0, L),
              vmin=0, vmax=1.6, interpolation="nearest")
    m = 16 * h
    core = m + round(0.625 / h) * h
    for (a, b), col, et, ls in (((0, L), TINTA, "bloque (VOI completo)", "-"),
                                ((m, L - m), NARANJA, "cubo de ensayo", "-"),
                                ((core, L - core), AZUL,
                                 "núcleo (0,625 mm)", "--")):
        ax.plot([a, b, b, a, a], [a, a, b, b, a], color=col, lw=1.8, ls=ls,
                label=et)
    ax.fill_between([core, L - core], L - core - 0.25, L - core,
                    color=AZUL, alpha=0.25, lw=0)
    ax.fill_between([core, L - core], core, core + 0.25, color=AZUL,
                    alpha=0.25, lw=0)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("z (mm), carga vertical")
    ax.grid(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), fontsize=7.2,
              ncol=1)
    fig.tight_layout()
    guardar(fig, FIGS, "f_esquema_real")


if __name__ == "__main__":
    fig_validacion_E()
    fig_validacion_p99()
    fig_practico()
    fig_a3()
    fig_resolucion()
    fig_perfil_bvtv()
    if len(sys.argv) > 1:
        fig_esquema_real(sys.argv[1])
    print("figuras en", FIGS)
