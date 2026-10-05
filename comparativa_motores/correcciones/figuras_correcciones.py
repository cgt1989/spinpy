"""
figuras_correcciones.py: Figuras del informe de correcciones de A1, A2 y A3.

    cd comparativa_motores/correcciones
    python figuras_correcciones.py      # lee resultados/*.json

Escribe PNG y SVG en informe_correcciones/figs/.
"""

import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI / "informe_correcciones"))
from estilo_figuras import (AQUA, AZUL, AZULES, GRIS, NARANJA,  # noqa: E402
                            TINTA2, aplicar, guardar)
import matplotlib.pyplot as plt                                  # noqa: E402

RES = AQUI / "resultados"
FIGS = AQUI / "informe_correcciones" / "figs"
dec = json.loads((RES / "decision.json").read_text())
ref = json.loads((RES / "referencia.json").read_text())
hx = json.loads((RES / "hex8.json").read_text())
val = json.loads((RES / "validacion_app.json").read_text())
aplicar()

METODOS = [("B0", "Línea base: tracción, VOI completo"),
           ("F1", "F1: plato rígido, VOI completo"),
           ("F2", "F2: tracción, núcleo"),
           ("F12", "F1 + F2: plato rígido, núcleo")]
NS = (32, 48, 64)


def puntos(clave, nombre, xlim, xticks, leyenda="lower right"):
    fig, ax = plt.subplots(figsize=(6.6, 2.7))
    y = np.arange(len(METODOS))[::-1]
    ax.axvspan(-5, 5, color="#f0efec", lw=0, zorder=0)
    for n, col, dy in zip(NS, AZULES, (0.18, 0, -0.18)):
        for yi, (m, _) in zip(y, METODOS):
            v = 100 * dec["hex8"][str(n)][m][clave]
            ax.plot([v], [yi + dy], "o", ms=6, color=col, mec="white",
                    mew=0.8, label=f"{n}³" if m == "B0" else None)
    ax.axvline(0, color=TINTA2, lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([t for _, t in METODOS])
    ax._cat_y = True
    ax.grid(axis="y", visible=False)
    ax.set_xlim(*xlim)
    ax.set_xticks(xticks)
    ax.set_xlabel(nombre)
    ax.legend(title="hex8", loc=leyenda, ncol=3, fontsize=7.5,
              title_fontsize=7.5)
    fig.tight_layout()
    return fig


def fig_E():
    guardar(puntos("s_E", "error de E aparente frente a la referencia embebida (%)",
                   (-60, 10), [-60, -50, -40, -30, -20, -10, -5, 0, 5, 10],
                   leyenda="lower left"),
            FIGS, "f_error_E")


def fig_p99():
    guardar(puntos("s_p99", "error del p99 de von Mises frente a la referencia (%)",
                   (-20, 60), [-20, -10, -5, 0, 5, 10, 20, 30, 40, 50, 60]),
            FIGS, "f_error_p99")


def fig_valores():
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.8))
    for ax, reg, tit, met in ((axs[0], "voi", "(a) Región del VOI completo",
                               (("B0", NARANJA, "tracción (app hoy)"),
                                ("F1", AZUL, "plato rígido"))),
                              (axs[1], "nucleo", "(b) Núcleo del VOI",
                               (("F2", NARANJA, "tracción"),
                                ("F12", AZUL, "plato rígido")))):
        r = [ref[f"n{n}"][reg]["E"] for n in NS]
        ax.plot(NS, r, color=TINTA2, ls="--", marker="s", ms=4,
                label="referencia embebida")
        for m, col, et in met:
            ax.plot(NS, [dec["hex8"][str(n)][m]["E"] for n in NS], color=col,
                    marker="o", ms=4.5, label=et)
        ax.set_title(tit)
        ax.set_xticks(NS)
        ax.set_xlim(26, 70)
        ax.set_ylim(0, 760)
        ax.set_xlabel("vóxeles por lado (n)")
        ax.legend(loc="lower right", fontsize=7.5)
    axs[0].set_ylabel("E aparente (MPa)")
    fig.tight_layout()
    guardar(fig, FIGS, "f_valores_E")


def fig_campo():
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.6))
    pares = (("n{}_fuerza", "campo_voi", "tracción, VOI"),
             ("n{}_plato", "campo_voi", "plato, VOI"),
             ("n{}_fuerza", "campo_nucleo", "tracción, núcleo"),
             ("n{}_plato", "campo_nucleo", "plato, núcleo"))
    x = np.arange(len(pares))
    for ax, clave, tit, ylab in ((axs[0], "err_mediana",
                                  "(a) Error local mediano",
                                  "|Δ von Mises| / referencia"),
                                 (axs[1], "cola_recuperada",
                                  "(b) Zonas calientes recuperadas",
                                  "fracción del 1 % superior")):
        for n, col, dx in zip(NS, AZULES, (-0.24, 0, 0.24)):
            v = [hx[k.format(n)][c][clave] for k, c, _ in pares]
            ax.bar(x + dx, v, width=0.22, color=col,
                   label=f"{n}³" if clave == "err_mediana" else None)
        ax.set_xticks(x)
        ax.set_xticklabels([e for _, _, e in pares], rotation=20, ha="right",
                           fontsize=7.5)
        ax._cat_x = True
        ax.set_title(tit)
        ax.set_ylabel(ylab)
        ax.set_ylim(0, 0.65)
    axs[0].legend(title="hex8", fontsize=7.5, title_fontsize=7.5, ncol=3,
                  loc="upper right")
    fig.tight_layout()
    guardar(fig, FIGS, "f_campo")


def fig_a3():
    fig, ax = plt.subplots(figsize=(6.6, 2.6))
    y = np.arange(len(METODOS))[::-1]
    for n, col, dy in ((32, AZULES[0], 0.17), (48, AZULES[1], -0.17)):
        for yi, (m, _) in zip(y, METODOS):
            v = 100 * dec["a3"][str(n)][m]["brecha_E"]
            ax.barh(yi + dy, v, height=0.3, color=col,
                    label=f"{n}³" if m == "B0" else None)
            ax.text(v - 0.8, yi + dy, f"{v:+.1f} %".replace(".", ",")
                    .replace("-", "−"), va="center", ha="right",
                    fontsize=7, color=TINTA2)
    ax.axvline(0, color=TINTA2, lw=0.9)
    ax.axvline(-10, color=GRIS, lw=0.9, ls="--")
    ax.text(-9.4, y[0] + 0.45, "umbral −10 %", fontsize=7.5, color=TINTA2)
    ax.set_yticks(y)
    ax.set_yticklabels([t for _, t in METODOS])
    ax._cat_y = True
    ax.grid(axis="y", visible=False)
    ax.set_xlim(-48, 2)
    ax.set_xlabel("E de TET10 frente a hex8 con el mismo protocolo (%)")
    ax.legend(title="TET10", loc="lower left", fontsize=7.5,
              title_fontsize=7.5)
    fig.tight_layout()
    guardar(fig, FIGS, "f_a3")


def fig_validacion():
    fig, ax = plt.subplots(figsize=(6.6, 2.5))
    casos = [k for k in val]
    x = np.arange(len(casos))
    b = [100 * val[k]["err_E_base"] for k in casos]
    c = [100 * val[k]["err_E_corregido"] for k in casos]
    ax.bar(x - 0.18, b, width=0.34, color=NARANJA, label="E publicado hoy")
    ax.bar(x + 0.18, c, width=0.34, color=AZUL, label="E corregido")
    for xi, v in zip(x + 0.18, c):
        ax.text(xi, v - 2.5 if v < 0 else v + 1, f"{v:+.1f} %".replace(".", ",")
                .replace("-", "−"), ha="center",
                va="top" if v < 0 else "bottom", fontsize=7, color=TINTA2)
    ax.axhline(0, color=TINTA2, lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([k.replace("_n", " ").replace("hex8", "hex8")
                        .replace("tet10", "TET10") + "³" for k in casos])
    ax._cat_x = True
    ax.set_ylabel("error frente a la referencia (%)")
    ax.set_ylim(-80, 8)
    ax.legend(loc="lower right")
    fig.tight_layout()
    guardar(fig, FIGS, "f_validacion")


if __name__ == "__main__":
    fig_E()
    fig_p99()
    fig_valores()
    fig_campo()
    fig_a3()
    fig_validacion()
    print("figuras en", FIGS)
