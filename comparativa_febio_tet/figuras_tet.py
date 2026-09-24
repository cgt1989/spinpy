"""figuras_tet.py — Figuras de la validacion de la malla suave (TET10).

    cd Port_Python
    python comparativa_febio_tet/figuras_tet.py

Solo LEE `resultados/cavidad.jsonl`; no resuelve nada. Deja en `figs/`:
  fig_cavidad_pico.png   pico y p99 de von Mises en la pared frente a la
                         resolucion: hex8, TET10 y la referencia analitica.
  fig_cavidad_suavizado.png  el mismo pico frente a las iteraciones de Taubin.
  fig_cavidad_rigidez.png    E_app de las tres mallas.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AQUI = Path(__file__).resolve().parent
FIGS = AQUI / "figs"
GOODIER = 1.9808
AZUL, NARANJA, GRIS = "#1f4e9c", "#b06000", "#555555"


def leer():
    ult = {}
    for l in open(AQUI / "resultados" / "cavidad.jsonl", encoding="utf-8"):
        f = json.loads(l)
        if f.get("ok"):
            ult[f["etiqueta"]] = f        # la ultima corrida de cada etiqueta
    return ult


def main():
    FIGS.mkdir(exist_ok=True)
    d = leer()
    ref = d.get("ref_esfera96") or d.get("ref_esfera64")
    ref64 = d.get("ref_esfera64")
    ns = sorted({f["n"] for f in d.values() if f.get("tipo_caso") == "voxeles"})

    def serie(pref, clave, ns_):
        x, y = [], []
        for n in ns_:
            f = d.get(f"{pref}_n{n}")
            if f:
                x.append(n * 0.2)
                y.append(f[clave])
        return x, y

    # -- pico y p99 frente a la resolucion --
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=True)
    for ax, clave, tit in ((axs[0], "vm_max_pared", "máximo en la pared"),
                           (axs[1], "vm_p99_pared", "p99 de la capa superficial")):
        for pref, col, et in (("hex8", AZUL, "ladrillos (hex8)"),
                              ("tet10", NARANJA, "malla suave (TET10)")):
            x, y = serie(pref, clave, ns)
            ax.plot(x, y, "o-", color=col, label=et)
        ax.axhline(ref[clave], color=GRIS, lw=1.2,
                   label="esfera analítica (TET10 fino)")
        ax.axhspan(ref[clave] * 0.97, ref[clave] * 1.03, color=GRIS, alpha=0.12,
                   label="±3 % (tolerancia declarada)")
        ax.axhline(GOODIER, color=GRIS, lw=0.8, ls=":",
                   label="Goodier, medio infinito")
        ax.set_title(tit, fontsize=10)
        ax.set_xlabel("r / h  (vóxeles en el radio)")
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("σ von Mises / σ₀")
    h, e = axs[1].get_legend_handles_labels()
    fig.legend(h, e, fontsize=7.5, frameon=False, loc="lower center", ncol=5)
    fig.suptitle("Cavidad esférica (r = 0,2 L), ensayo de la app: el pico no "
                 "converge con ninguna de las dos mallas", fontsize=10)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(FIGS / "fig_cavidad_pico.png", dpi=200)
    plt.close(fig)

    # -- suavizado --
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    for n, mk in ((32, "o"), (40, "s")):
        for dec, ls, et in ((None, "-", "decimado 0,5"),
                            (0, "--", "sin decimar")):
            x, y = [], []
            for it in (20, 60, 120):
                if dec is None:
                    k = f"tet10_n{n}" if it == 20 else f"tet10_n{n}_s{it}"
                else:
                    k = f"tet10_n{n}_s{it}_d0"
                if k in d:
                    x.append(it)
                    y.append(d[k]["vm_max_pared"] / ref["vm_max_pared"] - 1)
            ax.plot(x, [100 * v for v in y], mk + ls,
                    color=AZUL if n == 32 else NARANJA,
                    label=f"n = {n}, {et}")
    ax.axhspan(-3, 3, color=GRIS, alpha=0.12)
    ax.axhline(0, color=GRIS, lw=1)
    ax.set_xlabel("iteraciones de Taubin")
    ax.set_ylabel("pico TET10 − referencia [%]")
    ax.set_title("El pico depende del suavizado, un parámetro de modelado",
                 fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_cavidad_suavizado.png", dpi=200)
    plt.close(fig)

    # -- rigidez --
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    for pref, col, et in (("hex8", AZUL, "ladrillos (hex8)"),
                          ("tet10", NARANJA, "malla suave (TET10)")):
        x, y = serie(pref, "E_app_MPa", ns)
        ax.plot(x, [100 * (v / ref["E_app_MPa"] - 1) for v in y], "o-",
                color=col, label=et)
    if ref64:
        ax.plot([], [], " ", label=f"referencias: 64 → 96 cambia "
                f"{100 * (ref['E_app_MPa'] / ref64['E_app_MPa'] - 1):+.3f} %"
                .replace(".", ","))
    ax.axhline(0, color=GRIS, lw=1)
    ax.set_xlabel("r / h  (vóxeles en el radio)")
    ax.set_ylabel("E_app − referencia [%]")
    ax.set_title("La rigidez sí es estable (< 0,5 %)", fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_cavidad_rigidez.png", dpi=200)
    plt.close(fig)
    print("figuras en", FIGS)


if __name__ == "__main__":
    main()
