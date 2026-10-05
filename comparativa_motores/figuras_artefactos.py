"""
figuras_artefactos.py: Analisis complementario y figuras del informe de artefactos.

    cd comparativa_motores
    python artefactos.py            # antes: resuelve y guarda los campos
    python figuras_artefactos.py    # figuras en informe_artefactos/figs/

Lee `resultados/artefactos.json` y los campos por elemento de
`resultados/artefactos_campos/*.npz`, calcula los perfiles de tension junto al
techo y a las caras laterales y el nivel de azar de la persistencia de zonas
calientes, y escribe esos numeros en `resultados/artefactos_extra.json`.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
import numpy as np                                               # noqa: E402
from matplotlib import font_manager                              # noqa: E402
from scipy.spatial import cKDTree                                # noqa: E402

AQUI = Path(__file__).resolve().parent
RES = AQUI / "resultados"
CAMPOS = RES / "artefactos_campos"
FIGS = AQUI / "informe_artefactos" / "figs"
LADO, BANDA = 5.0, 0.30

for f in (AQUI / "informe_artefactos" / "fuentes").glob("PlexSans-*.ttf"):
    font_manager.fontManager.addfont(str(f))
plt.rcParams.update({
    "font.family": ["IBM Plex Sans", "DejaVu Sans"], "font.size": 8.5,
    "axes.titlesize": 9, "axes.labelsize": 8.5, "axes.titleweight": 600,
    "axes.titlelocation": "left", "axes.spines.top": False,
    "axes.spines.right": False, "axes.edgecolor": "#8a8984",
    "axes.labelcolor": "#2b2a27", "xtick.color": "#52514e",
    "ytick.color": "#52514e", "axes.grid": True, "grid.color": "#e6e5e0",
    "grid.linewidth": 0.6, "legend.frameon": False, "figure.dpi": 150,
    "savefig.dpi": 300, "savefig.bbox": "tight", "lines.linewidth": 1.8,
    "svg.fonttype": "path"})

# Paleta categorica (orden fijo, validada) y rampa secuencial azul
AZUL, NARANJA, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GRIS, TINTA, TINTA2 = "#8a8984", "#0b0b0b", "#52514e"


def cargar(k):
    d = np.load(CAMPOS / f"{k}.npz")
    return {x: d[x] for x in d.files}


def pond(v, w, q):
    o = np.argsort(v)
    v, w = v[o], w[o]
    c = np.cumsum(w) - 0.5 * w
    return float(np.interp(q / 100.0, c / w.sum(), v))


# ---------------------------------------------------------------------------
# Analisis complementario
# ---------------------------------------------------------------------------

def perfil(c, coord, bordes):
    """Media ponderada y p95 de von Mises (todo el tejido) por franja."""
    vm, w = c["vm"], c["vol"]
    med, p95 = [], []
    for a, b in zip(bordes[:-1], bordes[1:]):
        m = (coord >= a) & (coord < b)
        med.append(float(np.average(vm[m], weights=w[m])))
        p95.append(pond(vm[m], w[m], 95))
    return np.array(med), np.array(p95)


def dist_lateral(cen):
    return np.minimum(cen[:, :2], LADO - cen[:, :2]).min(1)


def azar_persistencia(ka, kb, n_rep=20, semilla=0):
    """Fraccion de elementos superficiales de `ka` ELEGIDOS AL AZAR (tantos
    como su cola) con un elemento de la cola de `kb` a <= BANDA mm."""
    a, b = cargar(ka), cargar(kb)
    sa, sb = a["sup"], b["sup"]
    ca = a["cen"][sa]
    vb = b["vm"][sb]
    qb = pond(vb, b["vol"][sb], 99)
    arbol = cKDTree(b["cen"][sb][vb >= qb])
    va = a["vm"][sa]
    n = int((va >= pond(va, a["vol"][sa], 99)).sum())
    rng = np.random.default_rng(semilla)
    f = [float((arbol.query(ca[rng.choice(len(ca), n, replace=False)])[0]
                <= BANDA).mean()) for _ in range(n_rep)]
    return float(np.mean(f)), float(np.std(f))


def main():
    FIGS.mkdir(parents=True, exist_ok=True)
    r = json.loads((RES / "artefactos.json").read_text())
    extra = {}

    # Perfiles en z y frente a la distancia a las caras laterales
    bz = np.linspace(0, LADO, 21)
    bl = np.linspace(0, 1.0, 11)
    perf = {}
    for k in ("hex8_n48_fuerza_deslizante", "hex8_n48_plato_deslizante",
              "tet10_n48_fuerza_deslizante", "tet10_n48_plato_deslizante",
              "hex8_n64_fuerza_deslizante"):
        c = cargar(k)
        mz, pz = perfil(c, c["cen"][:, 2], bz)
        ml, pl = perfil(c, dist_lateral(c["cen"]), bl)
        # Fraccion del volumen a <= BANDA de una cara lateral y su tension
        # media relativa a la del nucleo (> 2 BANDA)
        d = dist_lateral(c["cen"])
        w = c["vol"]
        borde, nucleo = d <= BANDA, d > 2 * BANDA
        perf[k] = {"z": mz.tolist(), "z_p95": pz.tolist(), "lat": ml.tolist(),
                   "lat_p95": pl.tolist(),
                   "frac_vol_borde": float(w[borde].sum() / w.sum()),
                   "vm_borde_sobre_nucleo": float(
                       np.average(c["vm"][borde], weights=w[borde])
                       / np.average(c["vm"][nucleo], weights=w[nucleo]))}
    extra["perfiles"] = perf
    extra["bordes_z_mm"] = bz.tolist()
    extra["bordes_lat_mm"] = bl.tolist()

    # Nivel de azar de la persistencia
    az = {}
    for par in r["persistencia"]:
        ka, kb = par.split("|")
        m, s = azar_persistencia(ka, kb)
        az[par] = {"azar_media": m, "azar_de": s}
    extra["persistencia_azar"] = az
    (RES / "artefactos_extra.json").write_text(json.dumps(extra, indent=1))
    for k, v in az.items():
        print(k, f"{r['persistencia'][k]['cola_cerca']:.3f} frente a azar "
              f"{v['azar_media']:.3f} +- {v['azar_de']:.3f}")
    for k, v in perf.items():
        print(k, f"borde {v['frac_vol_borde']:.3f} del volumen; vm borde/"
              f"nucleo {v['vm_borde_sobre_nucleo']:.3f}")

    fig_perfil_z(perf, bz)
    fig_perfil_lateral(perf, bl)
    fig_enriquecimiento(r)
    fig_dp99(r)
    fig_convergencia(r)
    fig_zz(r)
    fig_persistencia(r, az)
    fig_corte()


# ---------------------------------------------------------------------------
# Figuras
# ---------------------------------------------------------------------------

def _coma(v, _pos=None):
    t = f"{v:g}".replace(".", ",").replace("-", "\u2212")
    return t


def _guardar(fig, nombre):
    from matplotlib.ticker import FuncFormatter, FixedFormatter
    for ax in fig.axes:
        for eje, escala in ((ax.xaxis, ax.get_xscale()),
                            (ax.yaxis, ax.get_yscale())):
            if eje is ax.yaxis and getattr(ax, "_cat_y", False):
                continue
            if escala == "linear" and not isinstance(eje.get_major_formatter(),
                                                     FixedFormatter):
                eje.set_major_formatter(FuncFormatter(_coma))
    fig.savefig(FIGS / f"{nombre}.png")
    fig.savefig(FIGS / f"{nombre}.svg")
    plt.close(fig)


def _rotulo(ax, x, y, texto, color):
    ax.annotate(texto, (x[-1], y[-1]), xytext=(4, 0),
                textcoords="offset points", va="center", fontsize=8,
                color=TINTA2)
    ax.plot([x[-1]], [y[-1]], "o", ms=4, color=color)


def fig_perfil_z(perf, bz):
    zc = 0.5 * (bz[:-1] + bz[1:])
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True)
    for ax, malla, tit in ((axs[0], "hex8", "(a) Ladrillos (hex8), 48³"),
                           (axs[1], "tet10", "(b) Malla suave (TET10), 48³")):
        for control, col, et in (("fuerza", NARANJA, "tracción uniforme"),
                                 ("plato", AZUL, "plato rígido")):
            v = perf[f"{malla}_n48_{control}_deslizante"]["z_p95"]
            ax.plot(zc, v, color=col, marker="o", ms=3.5, label=et)
        ax.axvspan(LADO - BANDA, LADO, color="#f0efec", zorder=0, lw=0)
        ax.axvspan(0, BANDA, color="#f0efec", zorder=0, lw=0)
        ax.set_title(tit)
        ax.set_xlabel("altura z (mm): base a la izquierda, techo a la derecha")
        ax.set_xlim(0, LADO)
    axs[0].set_ylabel("p95 de von Mises / σ aparente")
    axs[1].legend(loc="upper left")
    fig.tight_layout()
    _guardar(fig, "f_perfil_z")


def fig_perfil_lateral(perf, bl):
    dc = 0.5 * (bl[:-1] + bl[1:])
    fig, ax = plt.subplots(figsize=(4.4, 2.8))
    for k, col, et in (("hex8_n48_fuerza_deslizante", AZUL, "hex8 48³"),
                       ("hex8_n64_fuerza_deslizante", AQUA, "hex8 64³"),
                       ("tet10_n48_fuerza_deslizante", NARANJA, "TET10 48³")):
        v = np.array(perf[k]["lat"])
        ax.plot(dc, v / v[-3:].mean(), color=col, marker="o", ms=3.5,
                label=et)
    ax.axvspan(0, BANDA, color="#f0efec", zorder=0, lw=0)
    ax.axhline(1, color=GRIS, lw=0.8, ls="--")
    ax.set_xlabel("distancia a la cara lateral más próxima (mm)")
    ax.set_ylabel("von Mises medio / valor del núcleo")
    ax.set_title("Tensión junto a las caras cortadas del VOI")
    ax.legend(loc="lower right")
    fig.tight_layout()
    _guardar(fig, "f_perfil_lateral")


CLASES = [("techo", "Franja del techo cargado"),
          ("union_singular", "Uniones por arista o vértice (hex8)"),
          ("base", "Franja de la base"),
          ("arista_entrante", "Aristas cóncavas de la escalera (hex8)"),
          ("mala_calidad", "Tetraedros de mala calidad (TET10)"),
          ("borde_lateral", "Franja de las caras laterales"),
          ("ancla", "Entorno de los nodos ancla")]


def fig_enriquecimiento(r):
    casos = [("hex8_n32_fuerza_deslizante", "hex8 32³", "#86b6ef"),
             ("hex8_n48_fuerza_deslizante", "hex8 48³", AZUL),
             ("hex8_n64_fuerza_deslizante", "hex8 64³", "#184f95"),
             ("tet10_n48_fuerza_deslizante", "TET10 48³", NARANJA)]
    fig, ax = plt.subplots(figsize=(6.6, 3.1))
    y = np.arange(len(CLASES))[::-1]
    for (k, et, col), dy in zip(casos, (0.24, 0.08, -0.08, -0.24)):
        for yi, (c, _) in zip(y, CLASES):
            d = r[k]["clases"].get(c)
            if d is None or not np.isfinite(d.get("enriquecimiento", np.nan)):
                continue
            e = max(d["enriquecimiento"], 0.02)
            ax.plot([e], [yi + dy], "o", ms=5.5, color=col,
                    mec="white", mew=0.8, label=et if yi == y[0] else None)
    ax.axvline(1, color=TINTA2, lw=0.9)
    ax.text(0.95, y[0] + 0.5, "1 = tanto como su volumen", fontsize=7.5,
            color=TINTA2, ha="right")
    ax.set_xscale("log")
    ax.set_xlim(0.015, 30)
    ax.set_xticks([0.02, 0.1, 0.5, 1, 2, 5, 10, 20])
    ax.set_xticklabels(["≈0", "0,1", "0,5", "1", "2", "5", "10", "20"])
    ax.set_yticks(y)
    ax.set_yticklabels([t for _, t in CLASES])
    ax._cat_y = True
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("enriquecimiento de la cola (> p99) en cada zona, escala log")
    ax.legend(loc="lower right", ncol=4, fontsize=7.5,
              bbox_to_anchor=(1.0, -0.42))
    fig.tight_layout()
    _guardar(fig, "f_enriquecimiento")


def fig_dp99(r):
    casos = [("hex8_n48_fuerza_deslizante", "hex8 48³", AZUL),
             ("tet10_n48_fuerza_deslizante", "TET10 48³", NARANJA)]
    cl = [c for c in CLASES if c[0] not in ("ancla", "arista_entrante")]
    fig, ax = plt.subplots(figsize=(6.6, 2.6))
    y = np.arange(len(cl))[::-1]
    for (k, et, col), dy in zip(casos, (0.17, -0.17)):
        for yi, (c, _) in zip(y, cl):
            d = r[k]["clases"].get(c)
            if d is None:
                continue
            ax.barh(yi + dy, d["dp99_pct"], height=0.3, color=col,
                    label=et if yi == y[0] else None)
            ax.text(d["dp99_pct"] + (0.25 if d["dp99_pct"] >= 0 else -0.25),
                    yi + dy, f"{d['dp99_pct']:+.1f} %".replace(".", ",").replace("-", "\u2212"),
                    va="center", ha="left" if d["dp99_pct"] >= 0 else "right",
                    fontsize=7, color=TINTA2)
    ax.axvline(0, color=TINTA2, lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([t for _, t in cl])
    ax._cat_y = True
    ax.grid(axis="y", visible=False)
    ax.set_xlim(-13, 10)
    ax.set_xlabel("cambio del p99 de superficie al excluir la zona (%)")
    ax.legend(loc="upper right")
    fig.tight_layout()
    _guardar(fig, "f_dp99")


def fig_convergencia(r):
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.8))
    for tipo, ns, col, et in (("hex8", (32, 48, 64), AZUL, "hex8"),
                              ("tet10", (32, 48), NARANJA, "TET10")):
        ks = [f"{tipo}_n{n}_fuerza_deslizante" for n in ns]
        E = [r[k]["E_app_area_MPa"] for k in ks]
        p = [r[k]["p99_sup"] for k in ks]
        axs[0].plot(ns, E, color=col, marker="o", ms=4.5)
        axs[1].plot(ns, p, color=col, marker="o", ms=4.5)
        _rotulo(axs[0], ns, E, et, col)
        _rotulo(axs[1], ns, p, et, col)
    for ax in axs:
        ax.set_xticks((32, 48, 64))
        ax.set_xlim(28, 72)
        ax.set_xlabel("vóxeles por lado (n)")
    axs[0].set_title("(a) Rigidez aparente")
    axs[0].set_ylabel("E aparente (MPa, por área)")
    axs[0].set_ylim(0, 340)
    axs[1].set_title("(b) Tensión citada")
    axs[1].set_ylabel("p99 de von Mises / σ aparente")
    axs[1].set_ylim(0, 95)
    fig.tight_layout()
    _guardar(fig, "f_convergencia")


def fig_zz(r):
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.8))
    for tipo, ns, col, et in (("hex8", (32, 48, 64), AZUL, "hex8"),
                              ("tet10", (32, 48), NARANJA, "TET10")):
        h = [LADO / n for n in ns]
        z = [r[f"{tipo}_n{n}_fuerza_deslizante"]["zz_global"] for n in ns]
        axs[0].plot(h, z, color=col, marker="o", ms=4.5, label=et)
        pend = np.polyfit(np.log(h), np.log(z), 1)[0]
        axs[0].text(0.97, 0.20 if tipo == "hex8" else 0.08,
                    f"{et}: pendiente {pend:.2f}".replace(".", ","),
                    transform=axs[0].transAxes, ha="right", fontsize=7.5,
                    color=TINTA2)
    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_xticks([0.08, 0.1, 0.12, 0.16])
    axs[0].set_xticklabels(["0,08", "0,10", "0,12", "0,16"])
    axs[0].set_yticks([0.25, 0.3, 0.4, 0.5])
    axs[0].set_yticklabels(["0,25", "0,30", "0,40", "0,50"])
    axs[0].minorticks_off()
    axs[0].set_xlabel("tamaño de vóxel h (mm)")
    axs[0].set_ylabel("indicador ZZ global")
    axs[0].set_title("(a) Error de recuperación al refinar")
    axs[0].legend(loc="upper left", fontsize=7.5)
    c = cargar("hex8_n48_fuerza_deslizante")
    s = c["sup"]
    rel = c["eta"][s] / np.maximum(c["vm"][s], 1e-12)
    cola = c["vm"][s] >= pond(c["vm"][s], c["vol"][s], 99)
    bins = np.linspace(0, 1.5, 31)
    axs[1].hist(rel[~cola], bins=bins, density=True, color="#9ec5f4",
                label="resto de la superficie")
    axs[1].hist(rel[cola], bins=bins, density=True, histtype="step",
                color=NARANJA, lw=1.8, label="cola (> p99)")
    axs[1].set_xlabel("indicador ZZ del elemento / su von Mises")
    axs[1].set_ylabel("densidad")
    axs[1].set_title("(b) hex8 48³: incertidumbre local")
    axs[1].legend(loc="upper right", fontsize=7.5)
    fig.tight_layout()
    _guardar(fig, "f_zz")


def fig_persistencia(r, az):
    pares = [("hex8_n32_fuerza_deslizante|hex8_n64_fuerza_deslizante",
              "hex8 32³ → 64³"),
             ("hex8_n48_fuerza_deslizante|hex8_n64_fuerza_deslizante",
              "hex8 48³ → 64³"),
             ("tet10_n32_fuerza_deslizante|tet10_n48_fuerza_deslizante",
              "TET10 32³ → 48³"),
             ("hex8_n48_fuerza_deslizante|tet10_n48_fuerza_deslizante",
              "hex8 48³ → TET10 48³"),
             ("hex8_n48_fuerza_deslizante|hex8_n48_fuerza_empotrado",
              "base deslizante → empotrada"),
             ("hex8_n48_fuerza_deslizante|hex8_n48_plato_deslizante",
              "tracción → plato (hex8)"),
             ("tet10_n48_fuerza_deslizante|tet10_n48_plato_deslizante",
              "tracción → plato (TET10)")]
    fig, ax = plt.subplots(figsize=(6.6, 2.9))
    y = np.arange(len(pares))[::-1]
    for yi, (k, et) in zip(y, pares):
        v = 100 * r["persistencia"][k]["cola_cerca"]
        a = 100 * az[k]["azar_media"]
        ax.barh(yi, v, height=0.55, color=AZUL)
        ax.plot([a, a], [yi - 0.32, yi + 0.32], color=TINTA, lw=1.6)
        ax.text(v + 1, yi, f"{v:.0f} %", va="center", fontsize=7.5,
                color=TINTA2)
    ax.plot([], [], color=TINTA, lw=1.6, label="nivel de azar")
    ax.set_yticks(y)
    ax.set_yticklabels([e for _, e in pares])
    ax._cat_y = True
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, 108)
    ax.set_xlabel("zonas calientes de la primera malla que reaparecen "
                  "a ≤ 0,3 mm en la segunda (%)")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, -0.02))
    fig.tight_layout()
    _guardar(fig, "f_persistencia")


def fig_corte():
    """Rebanada y in [2,2; 2,8] mm del campo hex8 48³: traccion frente a plato."""
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 3.3), sharey=True)
    vmax = None
    for ax, k, tit in ((axs[0], "hex8_n48_fuerza_deslizante",
                        "(a) Tracción uniforme en el techo"),
                       (axs[1], "hex8_n48_plato_deslizante",
                        "(b) Plato rígido en el techo")):
        c = cargar(k)
        m = (c["cen"][:, 1] > 2.2) & (c["cen"][:, 1] < 2.8)
        v = c["vm"][m]
        if vmax is None:
            vmax = float(np.percentile(v, 99.5))
        o = np.argsort(v)
        sc = ax.scatter(c["cen"][m, 0][o], c["cen"][m, 2][o], c=v[o], s=2.2,
                        cmap="Blues", vmin=0, vmax=vmax, linewidths=0,
                        rasterized=True)
        ax.set_aspect("equal")
        ax.set_xlim(0, LADO)
        ax.set_ylim(0, LADO)
        ax.set_title(tit)
        ax.set_xlabel("x (mm)")
        ax.grid(False)
    axs[0].set_ylabel("z (mm)")
    cb = fig.colorbar(sc, ax=axs, shrink=0.85, pad=0.02)
    cb.set_label("von Mises / σ aparente")
    _guardar(fig, "f_corte")


if __name__ == "__main__":
    main()
