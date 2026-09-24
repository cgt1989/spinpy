"""figuras_febio.py — Figuras del informe de la comparacion con FEBio.

    cd Port_Python
    python comparativa_febio/figuras_febio.py

Lee SOLO `resultados/*.jsonl` y los archivos de FEBio de `corridas/`; no
resuelve nada. Escribe PNG a 300 ppp en `figs/`. Colores y tintas de
`spinpy.figuras`, para que el informe tenga el aspecto del resto del proyecto.

  fig1  coincidencia lineal: cada metrica, cada caso, frente a su tolerancia
  fig2  E no lineal / E lineal a 1 MPa frente a BV/TV (fuerza y plato)
  fig3  E no lineal / E lineal frente a la carga, proximal 32 y 48
  fig4  donde se concentra el desvio no lineal (proximal 48, 1 MPa)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                              # noqa: E402
import numpy as np                                           # noqa: E402

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(AQUI.parent))
from spinpy.figuras import COLOR, EJE, REJILLA, TINTA, TINTA_2  # noqa: E402

RES = AQUI / "resultados"
FIGS = AQUI / "figs"
C_FUERZA, C_PLATO = COLOR["voi"], COLOR["spinodoide"]

plt.rcParams.update({"font.size": 9, "axes.edgecolor": EJE,
                     "axes.labelcolor": TINTA, "xtick.color": TINTA_2,
                     "ytick.color": TINTA_2, "axes.spines.top": False,
                     "axes.spines.right": False})


def _leer(nombre):
    ruta = RES / nombre
    if not ruta.exists():
        return []
    return [json.loads(l) for l in ruta.read_text(encoding="utf-8").splitlines()
            if l.strip()]


def _rejilla(ax, eje="y"):
    ax.grid(axis=eje, color=REJILLA, linewidth=0.6)
    ax.set_axisbelow(True)


def _guardar(fig, nombre):
    FIGS.mkdir(exist_ok=True)
    fig.savefig(FIGS / nombre, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  ", FIGS / nombre)


# ---------------------------------------------------------------------------

METRICAS = [("dE_rel", "E_app", 1e-6), ("du_max_rel", "u (máx.)", 1e-5),
            ("dsig_max_rel", "σ (máx.)", 1e-4),
            ("d_vm_p99_superficie_rel", "p99 sup.", 1e-5),
            ("d_factor_rel", "Pistoia", 1e-5), ("dF_rel", "equilibrio", 1e-9)]


def fig1(filas):
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    for i, (k, nombre, tol) in enumerate(METRICAS):
        v = np.array([f[k] for f in filas
                      if f.get(k) is not None and np.isfinite(f[k])])
        v = np.maximum(v, 1e-13)
        x = i + np.linspace(-0.22, 0.22, v.size)
        ax.scatter(x, v, s=22, color=C_FUERZA, edgecolor="white",
                   linewidth=0.6, zorder=3)
        ax.plot([i - 0.35, i + 0.35], [tol, tol], color=TINTA, linewidth=1.4)
    ax.set_yscale("log")
    ax.set_xticks(range(len(METRICAS)))
    ax.set_xticklabels([m[1] for m in METRICAS])
    ax.set_ylabel("diferencia relativa spinpy – FEBio")
    ax.set_xlim(-0.5, len(METRICAS) - 0.5)
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", color=C_FUERZA, label="un caso"),
        Line2D([], [], color=TINTA, linewidth=1.4,
               label="tolerancia predicha")],
        frameon=False, fontsize=7.5, loc="lower left")
    _rejilla(ax)
    _guardar(fig, "fig1_coincidencia_lineal.png")


def fig2(filas, plato):
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    pts = [(f["BVTV"], f["dE_nl_1MPa"], f["caso"]) for f in filas
           if f.get("dE_nl_1MPa") is not None
           and f["caso"] not in ("bloque", "proximal_48_empotrado")]
    for b, d, c in pts:
        m = "o" if c.endswith("_32") else "s"
        ax.scatter(b, 100 * abs(d), s=34, marker=m, color=C_FUERZA,
                   edgecolor="white", linewidth=0.6, zorder=3)
        if c.startswith("espinodoide"):
            ax.annotate("espinodoide ajustado", (b, 100 * abs(d)),
                        xytext=(8, 0), textcoords="offset points",
                        fontsize=7, color=TINTA_2, va="center")
        elif c.startswith("proximal") or c.endswith("_48"):
            nombre = c.split("_")[0]
            if c.endswith("_48"):
                ax.annotate(nombre, (b, 100 * abs(d)), xytext=(8, 0),
                            textcoords="offset points", fontsize=7,
                            color=TINTA_2, va="center")
    for p in plato:
        d = p.get("dE_plato_nl_1MPa")
        if d is not None:
            m = "o" if p["caso"].endswith("_32") else "s"
            ax.scatter(p["BVTV"], 100 * abs(d), s=34, marker=m,
                       color=C_PLATO, edgecolor="white", linewidth=0.6,
                       zorder=3)
    ax.set_yscale("log")
    ax.set_xlabel("BV/TV")
    ax.set_ylabel("|E no lineal / E lineal − 1| a 1 MPa (%)")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=C_FUERZA, label="fuerza impuesta (app)"),
         Line2D([], [], marker="o", ls="", color=C_PLATO, label="plato rígido"),
         Line2D([], [], marker="o", ls="", color=TINTA_2, label="32³"),
         Line2D([], [], marker="s", ls="", color=TINTA_2, label="48³")]
    ax.legend(handles=h, frameon=False, fontsize=7.5, loc="upper right")
    _rejilla(ax)
    _guardar(fig, "fig2_no_lineal_bvtv.png")


def fig3(barrido, plato):
    if not barrido:
        return
    fig, axs = plt.subplots(1, len(barrido), figsize=(3.3 * len(barrido), 3.0),
                            sharey=True, squeeze=False)
    fallo = {p["caso"]: p.get("sigma_fallo_pistoia_MPa") for p in plato}
    for ax, b in zip(axs[0], barrido):
        s = np.array(b["cargas_Pa"]) / 1e6
        for clave, color, nombre in (("fuerza", C_FUERZA, "fuerza impuesta"),
                                     ("plato", C_PLATO, "plato rígido")):
            y = np.array([np.nan if v is None else v for v in b[clave]])
            ax.plot(s, y, color=color, linewidth=2, marker="o", markersize=4,
                    label=nombre)
        ax.axvline(1.0, color=EJE, linewidth=1, ls="--")
        ax.text(0.95, 0.60, "carga de la app", transform=ax.get_xaxis_transform(),
                fontsize=6.5, color=TINTA_2, rotation=90, ha="right")
        for x, v in zip(s, b["fuerza"]):
            if v is None:
                ax.text(x, 0.35, "no converge", transform=ax.get_xaxis_transform(),
                        fontsize=6.5, color=C_FUERZA, ha="center",
                        rotation=90)
        sf = fallo.get(b["caso"])
        if sf:
            ax.axvline(sf, color=TINTA_2, linewidth=1, ls=":")
            ax.text(sf * 0.93, 0.60, "fallo Pistoia",
                    transform=ax.get_xaxis_transform(), fontsize=6.5,
                    color=TINTA_2, rotation=90, ha="right")
        ax.set_xscale("log")
        ax.set_title(b["caso"].replace("_", " ") + "³", fontsize=9)
        ax.set_xlabel("tensión aparente (MPa)")
        _rejilla(ax)
    axs[0][0].set_ylabel("E no lineal / E lineal")
    axs[0][0].legend(frameon=False, fontsize=7.5, loc="lower left")
    _guardar(fig, "fig3_barrido_carga.png")


def fig4(clave="proximal_48"):
    """Desvio no lineal por nodo, en planta del techo y proyectado en x-z."""
    import comparar_febio as cf
    from spinpy.resistencia import _solo_portante
    from spinpy.solido import malla_hex
    d = AQUI / "corridas" / clave
    if not (d / f"{clave}_1000000Pa_u.txt").exists():
        return
    BW, spc, _a, _f = cf.CASOS[clave]()
    nodos, _e, _i = malla_hex(_solo_portante(BW), spc)
    u = {s: cf._ultimo_registro(d / f"{clave}_{s}Pa_u.txt", 3) * (1e6 / s)
         for s in (1000, 2000, 1000000)}
    u_lin = 2 * u[1000] - u[2000]
    du = np.linalg.norm(u[1000000] - u_lin, axis=1) * 1e3     # micras
    ijk = np.round(nodos / spc).astype(int)
    n = BW.shape[0] + 1

    top = np.full((n, n), np.nan)
    m = ijk[:, 2] == ijk[:, 2].max()
    top[ijk[m, 0], ijk[m, 1]] = du[m]
    xz = np.full((n, BW.shape[2] + 1), np.nan)
    for (i, _j, k), v in zip(ijk, du):
        if np.isnan(xz[i, k]) or v > xz[i, k]:
            xz[i, k] = v

    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.2))
    vmax = np.nanmax(du)
    L = n * spc[0]
    for ax, img, tit, yl in ((axs[0], top, "a) techo (planta, z = H)", "y (mm)"),
                             (axs[1], xz, "b) máximo a lo largo de y", "z (mm)")):
        h = ax.imshow(img.T, origin="lower", cmap="Blues", vmin=0, vmax=vmax,
                      extent=(0, L, 0, img.shape[1] * spc[2]))
        ax.set_title(tit, fontsize=9)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel(yl)
        for s in ax.spines.values():
            s.set_visible(True)
    cb = fig.colorbar(h, ax=axs, shrink=0.85)
    cb.set_label("|u no lineal − u lineal| a 1 MPa (µm)")
    _guardar(fig, "fig4_mapa_desvio.png")
    print(f"   fig4: desvio max {vmax:.1f} µm, mediana {np.median(du):.3f} µm")


def main():
    filas = _leer("febio.jsonl")
    plato = _leer("plato.jsonl")
    barrido = _leer("barrido.jsonl")
    fig1(filas)
    fig2(filas, plato)
    fig3(barrido, plato)
    fig4()


if __name__ == "__main__":
    main()
