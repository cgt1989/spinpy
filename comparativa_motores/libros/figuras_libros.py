"""
figuras_libros.py: Figuras del informe de los cambios de la revision de la
biblioteca. Solo lee resultados/*.json.

    python comparativa_motores/libros/figuras_libros.py

Guarda PNG y SVG en informe_libros/figs/.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
RES = AQUI / "resultados"
FIGS = AQUI / "informe_libros" / "figs"
sys.path.insert(0, str(AQUI / "informe_libros"))
from estilo_figuras import (AQUA, AZUL, GRIS, NARANJA, TINTA2,   # noqa: E402
                            aplicar, guardar)

import matplotlib.pyplot as plt                                  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter  # noqa: E402,E501


def ticks_log(ax, eje, valores):
    """Marcas fijas con coma decimal en un eje logaritmico."""
    a = ax.yaxis if eje == "y" else ax.xaxis
    a.set_major_locator(FixedLocator(valores))
    a.set_major_formatter(FixedFormatter(
        [f"{v:g}".replace(".", ",") for v in valores]))
    a.set_minor_formatter(NullFormatter())


def leer(nombre):
    p = RES / nombre
    return json.loads(p.read_text()) if p.exists() else None


def f_flexion():
    d = leer("p6_flexion.json")
    nts = sorted(int(k) for k in d["error_rel"]["hex8"])
    fig, ax = plt.subplots(figsize=(5.2, 2.7))
    for el, c, et in (("hex8", AZUL, "hex8 (trilineal, V2.1.1)"),
                      ("hex8i", NARANJA, "hex8i (modos incompatibles)")):
        e = [100 * d["error_rel"][el][str(n)] for n in nts]
        ax.plot(nts, e, "o-", color=c, label=et)
    ax.axhline(0, color=GRIS, lw=0.8)
    ax.set_xlabel("elementos por canto de la viga")
    ax.set_ylabel("error de la flecha (%)")
    ax.set_xticks(nts)
    ax.legend(loc="lower right")
    guardar(fig, FIGS, "f_flexion")


def f_espinodoide():
    d6 = leer("p6_espinodoide.json")["24"]
    a6 = leer("p6_espinodoide_ajuste.json")
    d7 = leer("p7_hex_orden2.json")
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
    ax = axs[0]
    Einf = a6["Einf_MPa"]
    for el, c, et in (("hex8", AZUL, "hex8"), ("hex8i", NARANJA, "hex8i")):
        g = [d6[f"{el}_m{m}"]["n_dof"] for m in (1, 2, 3, 4)]
        e = [100 * (d6[f"{el}_m{m}"]["E_app"] / 1e6 / Einf - 1)
             for m in (1, 2, 3, 4)]
        ax.loglog(g, e, "o-", color=c, label=et)
    ticks_log(ax, "y", [20, 30, 50, 100])
    ax.set_title("a. Ensayo de la app (tracción)")
    ax.set_xlabel("grados de libertad")
    ax.set_ylabel("error de E aparente (%)")
    ax.legend()
    ax = axs[1]
    E1 = [d7[f"24_orden1_m{m}"] for m in (1, 2, 3)]
    E2 = [d7[f"24_orden2_m{m}"] for m in (1, 2)]
    from scipy.optimize import least_squares
    h = 1 / np.array([1, 2, 3.])
    y = np.array([x["E_app_MPa"] for x in E1])
    s = least_squares(lambda p: y - (p[0] + p[1] * h ** p[2]), [270, 250, 1])
    Einf7 = s.x[0]
    for serie, c, et in ((E1, AZUL, "hex8 (Q1)"), (E2, AQUA, "Q2")):
        ax.loglog([x["n_gdl"] for x in serie],
                  [100 * (x["E_app_MPa"] / Einf7 - 1) for x in serie], "o-",
                  color=c, label=et)
    ticks_log(ax, "y", [10, 20, 30, 50, 100])
    ax.set_title("b. NGSolve con plato rígido")
    ax.set_xlabel("grados de libertad")
    ax.set_ylabel("error de E aparente (%)")
    ax.legend()
    fig.tight_layout()
    guardar(fig, FIGS, "f_espinodoide")
    return Einf, Einf7


def f_material():
    sys.path.insert(0, str(AQUI.parent))
    from comun import uniaxial
    a = np.linspace(0.001, 0.75, 300)
    fig, ax = plt.subplots(figsize=(5.2, 2.7))
    ax.plot(100 * a, [-uniaxial("svk", 1 - x)[0] / 1000 for x in a],
            color=AZUL, label="St. Venant-Kirchhoff")
    ax.plot(100 * a, [-uniaxial("neohookeano", 1 - x)[0] / 1000 for x in a],
            color=NARANJA, label="neo-Hookeano")
    ax.axvline(100 * (1 - 1 / np.sqrt(3)), color=GRIS, lw=0.8, ls="--")
    ax.text(43.5, 16, "máximo de carga\ndel SVK (42 %)", fontsize=7.5,
            color=TINTA2)
    ax.set_xlabel("acortamiento uniaxial (%)")
    ax.set_ylabel("tensión nominal de compresión (GPa)")
    ax.set_ylim(0, 20)
    ax.legend(loc="upper left")
    guardar(fig, FIGS, "f_material")


def f_homogeneizacion():
    d = leer("p1_homogeneizacion_n16.json")
    d24 = leer("p1_homogeneizacion_n24.json")
    fig, axs = plt.subplots(1, 2 if d24 else 1, figsize=(7.2, 2.8),
                            squeeze=False)
    for ax, dd, n in zip(axs[0], (d, d24), (16, 24)):
        if not dd:
            continue
        rhos = sorted(dd, key=float)
        x = np.arange(len(rhos))
        for k, (v, c) in enumerate((("escalar", AZUL),
                                    ("traslaciones", NARANJA),
                                    ("rigidos", AQUA))):
            its = [max(dd[r][v]["iteraciones"]) if dd[r].get(v) else 0
                   for r in rhos]
            ok = [dd[r][v]["ok"] if dd[r].get(v) else False for r in rhos]
            b = ax.bar(x + (k - 1) * 0.27, its, 0.25, color=c,
                       label={"rigidos": "rígidos"}.get(v, v))
            for xi, o, hgt in zip(x + (k - 1) * 0.27, ok, its):
                if not o:
                    ax.text(xi, hgt + 8, "×", ha="center", color=TINTA2,
                            fontsize=9)
        ax.xaxis.set_major_locator(FixedLocator(x))
        ax.xaxis.set_major_formatter(FixedFormatter(
            [f"{float(r):.2f}".replace(".", ",") for r in rhos]))
        ax.set_xlabel("ρ nominal del espinodoide")
        ax.set_ylabel("iteraciones del CG (máx. 500)")
        ax.set_title(f"{'a' if n == 16 else 'b'}. {n}³")
        ax.set_ylim(0, 560)
    axs[0][0].legend(loc="upper right")
    fig.tight_layout()
    guardar(fig, FIGS, "f_homogeneizacion")


def f_pca():
    d = leer("p5_pca.json")
    fig, ax = plt.subplots(figsize=(5.2, 2.7))
    r = [f["cociente_s2_s1"] for f in d]
    ax.semilogy(r, [f["grados_por_pct_eje1"] for f in d], "o-", color=AZUL,
                label="sensibilidad estimada (° por 1 %)")
    ax.semilogy(r, [f["giro_medido_max"] for f in d], "s-", color=NARANJA,
                label="giro medido con ruido (°)")
    ticks_log(ax, "y", [0.3, 1, 3, 10, 30])
    ax.axhline(5, color=GRIS, lw=0.8, ls="--")
    ax.text(0.6, 6, "umbral de aviso", fontsize=7.5, color=TINTA2)
    ax.set_xlabel("cociente s2/s1 de las desviaciones")
    ax.set_ylabel("grados")
    ax.legend(loc="upper left")
    guardar(fig, FIGS, "f_pca")


def f_newton():
    d = leer("p2_newton.json")
    filas = []
    for tipo in ("hex8", "tet10"):
        for mat, et in (("svk", "SVK"), ("neohookeano", "neo")):
            for s in (2, 5, 10, 20):
                k = f"esp_{tipo}_24_{mat}_{s}MPa"
                if f"{k}_v211" in d and f"{k}_glob" in d:
                    filas.append((f"{tipo.upper() if tipo == 'tet10' else tipo}"
                                  f" {et} {s}", d[f"{k}_v211"],
                                  d[f"{k}_glob"]))
    fig, ax = plt.subplots(figsize=(7.2, 2.9))
    x = np.arange(len(filas))
    for j, (c, et) in enumerate(((AZUL, "Newton de la V2.1.1"),
                                 (NARANJA, "Newton globalizado"))):
        its = [sum(f[1 + j]["iter"]) if f[1 + j].get("ok") else 0
               for f in filas]
        ax.bar(x + (j - 0.5) * 0.38, its, 0.36, color=c, label=et)
        for xi, f, h in zip(x + (j - 0.5) * 0.38, filas, its):
            if not f[1 + j].get("ok"):
                ax.text(xi, 1, "×", ha="center", color=TINTA2, fontsize=10)
    ax.xaxis.set_major_locator(FixedLocator(x))
    ax.xaxis.set_major_formatter(FixedFormatter([f[0] for f in filas]))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7.5)
    ax.set_ylabel("iteraciones de Newton")
    ax.set_xlabel("malla, material y tensión aparente (MPa)")
    ax.legend(loc="upper left")
    guardar(fig, FIGS, "f_newton")


def f_cota():
    d = leer("p4_cota_error.json")
    if not d:
        return
    fig, ax = plt.subplots(figsize=(5.2, 2.9))
    for clave, v in d.items():
        c = AZUL if clave.startswith("hex8") else NARANJA
        ax.loglog(v["error_rel_medido"], v["cota_error_rel"], "o", color=c)
        ax.annotate(clave.replace("_", " ").replace("tet10", "TET10") + "³",
                    (v["error_rel_medido"], v["cota_error_rel"]),
                    textcoords="offset points", xytext=(5, 3), fontsize=7.5,
                    color=TINTA2)
    lo = min(min(v["error_rel_medido"], v["cota_error_rel"])
             for v in d.values()) / 3
    hi = max(max(v["error_rel_medido"], v["cota_error_rel"])
             for v in d.values()) * 3
    ax.plot([lo, hi], [lo, hi], color=GRIS, lw=0.8, ls="--")
    ax.text(hi / 2.5, hi / 9, "cota = error", fontsize=7.5, color=TINTA2)
    ax.set_xlabel("error relativo medido frente al directo")
    ax.set_ylabel("cota ε·κ estimada")
    ax.plot([], [], "o", color=AZUL, label="hex8")
    ax.plot([], [], "o", color=NARANJA, label="TET10")
    ax.legend(loc="upper left")
    guardar(fig, FIGS, "f_cota")


def main():
    aplicar()
    f_flexion()
    print("E_inf", f_espinodoide())
    f_material()
    f_homogeneizacion()
    f_pca()
    f_newton()
    f_cota()


if __name__ == "__main__":
    main()
