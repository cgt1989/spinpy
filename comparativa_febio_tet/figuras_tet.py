"""figuras_tet.py — Figuras de la validacion de la malla suave (TET10).

    cd Port_Python
    python comparativa_febio_tet/figuras_tet.py

Solo LEE `resultados/cavidad.jsonl`; no resuelve nada. Deja en `figs/`:
  fig_cavidad_pico.png   pico y p99 de von Mises en la pared frente a la
                         resolucion: hex8, TET10 y la referencia analitica.
  fig_cavidad_suavizado.png  el mismo pico frente a las iteraciones de Taubin.
  fig_cavidad_rigidez.png    E_app de las tres mallas.
y las mismas tres en ingles en `figs/en/` (las usa el README en ingles). Los
nombres en espanol se quedan en `figs/` porque INFORME.md los cita asi.
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

# Rotulos por idioma. En ingles el separador decimal es el punto.
TEXTOS = {
    "es": {
        "max": "máximo en la pared", "p99": "p99 de la capa superficial",
        "hex8": "ladrillos (hex8)", "tet10": "malla suave (TET10)",
        "ref": "esfera analítica (TET10 fino)",
        "tol": "±3 % (tolerancia declarada)",
        "goodier": "Goodier, medio infinito",
        "rh": "r / h  (vóxeles en el radio)",
        "sup_pico": "Cavidad esférica (r = 0,2 L), ensayo de la app: el pico "
                    "no converge con ninguna de las dos mallas",
        "dec": "decimado 0,5", "sindec": "sin decimar",
        "taubin": "iteraciones de Taubin",
        "y_suav": "pico TET10 − referencia [%]",
        "tit_suav": "El pico depende del suavizado, un parámetro de modelado",
        "refs": "referencias: 64 → 96 cambia {:+.3f} %",
        "y_rig": "E_app − referencia [%]",
        "tit_rig": "La rigidez sí es estable (< 0,5 %)",
        "coma": True,
    },
    "en": {
        "max": "wall maximum", "p99": "surface-layer p99",
        "hex8": "bricks (hex8)", "tet10": "smooth mesh (TET10)",
        "ref": "analytical sphere (fine TET10)",
        "tol": "±3 % (declared tolerance)",
        "goodier": "Goodier, infinite medium",
        "rh": "r / h  (voxels across the radius)",
        "sup_pico": "Spherical cavity (r = 0.2 L), app test: the peak does "
                    "not converge with either mesh",
        "dec": "decimated 0.5", "sindec": "not decimated",
        "taubin": "Taubin iterations",
        "y_suav": "TET10 peak − reference [%]",
        "tit_suav": "The peak depends on smoothing, a modelling parameter",
        "refs": "references: 64 → 96 changes by {:+.3f} %",
        "y_rig": "E_app − reference [%]",
        "tit_rig": "Stiffness is stable (< 0.5 %)",
        "coma": False,
    },
}


def leer():
    ult = {}
    for l in open(AQUI / "resultados" / "cavidad.jsonl", encoding="utf-8"):
        f = json.loads(l)
        if f.get("ok"):
            ult[f["etiqueta"]] = f        # la ultima corrida de cada etiqueta
    return ult


def main():
    d = leer()
    for idioma, destino in (("es", FIGS), ("en", FIGS / "en")):
        destino.mkdir(parents=True, exist_ok=True)
        figuras(d, TEXTOS[idioma], destino)
        print("figuras en", destino)


def figuras(d, T, destino):
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
    for ax, clave, tit in ((axs[0], "vm_max_pared", T["max"]),
                           (axs[1], "vm_p99_pared", T["p99"])):
        for pref, col, et in (("hex8", AZUL, T["hex8"]),
                              ("tet10", NARANJA, T["tet10"])):
            x, y = serie(pref, clave, ns)
            ax.plot(x, y, "o-", color=col, label=et)
        ax.axhline(ref[clave], color=GRIS, lw=1.2,
                   label=T["ref"])
        ax.axhspan(ref[clave] * 0.97, ref[clave] * 1.03, color=GRIS, alpha=0.12,
                   label=T["tol"])
        ax.axhline(GOODIER, color=GRIS, lw=0.8, ls=":",
                   label=T["goodier"])
        ax.set_title(tit, fontsize=10)
        ax.set_xlabel(T["rh"])
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("σ von Mises / σ₀")
    h, e = axs[1].get_legend_handles_labels()
    fig.legend(h, e, fontsize=7.5, frameon=False, loc="lower center", ncol=5)
    fig.suptitle(T["sup_pico"], fontsize=10)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(destino / "fig_cavidad_pico.png", dpi=200)
    plt.close(fig)

    # -- suavizado --
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    for n, mk in ((32, "o"), (40, "s")):
        for dec, ls, et in ((None, "-", T["dec"]),
                            (0, "--", T["sindec"])):
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
    ax.set_xlabel(T["taubin"])
    ax.set_ylabel(T["y_suav"])
    ax.set_title(T["tit_suav"], fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(destino / "fig_cavidad_suavizado.png", dpi=200)
    plt.close(fig)

    # -- rigidez --
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    for pref, col, et in (("hex8", AZUL, T["hex8"]),
                          ("tet10", NARANJA, T["tet10"])):
        x, y = serie(pref, "E_app_MPa", ns)
        ax.plot(x, [100 * (v / ref["E_app_MPa"] - 1) for v in y], "o-",
                color=col, label=et)
    if ref64:
        et = T["refs"].format(
            100 * (ref["E_app_MPa"] / ref64["E_app_MPa"] - 1))
        ax.plot([], [], " ", label=et.replace(".", ",") if T["coma"] else et)
    ax.axhline(0, color=GRIS, lw=1)
    ax.set_xlabel(T["rh"])
    ax.set_ylabel(T["y_rig"])
    ax.set_title(T["tit_rig"], fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(destino / "fig_cavidad_rigidez.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
