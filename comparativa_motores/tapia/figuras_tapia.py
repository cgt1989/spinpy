"""
figuras_tapia.py: Figuras y tablas de la replicacion del protocolo de Tapia
et al. (2026) con cada motor FEM.

    xvfb-run -a python figuras_tapia.py        # los renders 3D necesitan OpenGL

Solo LEE `campos/` y `resultados/tapia.jsonl` (los escribe `replicar_tapia.py`).
Las diferencias de cada motor se miden contra la REFERENCIA: NGSolve con el
resolvedor directo (residuo ~1e-14), por elemento y divididas por el p99 del
campo de referencia en el tejido, como en `comparativa_febio/porcino`.

Figuras (en `figs/`):
  mapa_vm_solvers, mapa_u_solvers      campo de cada motor, una escala por fila
  dif_vm_solvers,  dif_u_solvers       log10(|motor - referencia| / p99)
  resumen_diferencias                  maxima y p99 de la diferencia por motor
  distribucion_diferencias             ECDF de la diferencia por elemento
Tablas en `tablas_tapia.md`.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent.parent))

import matplotlib                                               # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                 # noqa: E402
from matplotlib import cm                                       # noqa: E402
from matplotlib.colors import Normalize                         # noqa: E402
from matplotlib.image import imread                             # noqa: E402

from spinpy import figuras as F                                 # noqa: E402

FIGS = AQUI / "figs"
REND = FIGS / "render"
CAMPOS = AQUI / "campos"
LADO = 700
VOIS = ("C2", "C4")
ESTR = ("voi", "spinodoide", "dual-lattice")
NOM_E = {"voi": "VOI", "spinodoide": "Spinodoide", "dual-lattice": "Dual-lattice"}
MOT = ("app", "ngsolve", "fenicsx", "skfem", "sfepy")
NOM = {"app": "App (spinpy)", "ngsolve": "NGSolve", "fenicsx": "FEniCSx",
       "skfem": "scikit-fem", "sfepy": "SfePy"}
COLOR = {"app": "#52514e", "ngsolve": "#2a78d6", "fenicsx": "#eb6834",
         "skfem": "#1baf7a", "sfepy": "#eda100"}
MARCA = {"app": "o", "ngsolve": "s", "fenicsx": "^", "skfem": "D", "sfepy": "v"}
MAG = {"vm": ("σ von Mises [MPa]", "inferno"),
       "u": ("Deformación total |u| [mm]", "viridis")}
PISO = -16.0                 # log10 de la diferencia por debajo del cual se trunca


def cargar(v, est):
    g = np.load(CAMPOS / f"{v}_{est}.npz")
    d = {"geo": {k: g[k] for k in g.files}, "campos": {}}
    for k in MOT + ("referencia",):
        f = CAMPOS / f"{v}_{est}__{k}.npz"
        if f.exists():
            z = np.load(f)
            d["campos"][k] = {m: z[f"{m}_{k}"] for m in ("vm", "u")}
    return d


def malla3d(geo, vals):
    c = np.full(tuple(int(x) for x in geo["forma"]), np.nan)
    i = geo["ijk"].astype(int)
    c[i[:, 0], i[:, 1], i[:, 2]] = vals
    return c


def render(campo, sp, clim, cmap, ruta):
    from scipy import ndimage
    ruta = Path(ruta)
    if ruta.exists():
        return ruta
    ruta.parent.mkdir(parents=True, exist_ok=True)
    solido = np.isfinite(campo)
    malla = F.malla_figura(solido, sp, suavizar=False)
    _d, (ix, iy, iz) = ndimage.distance_transform_edt(~solido,
                                                      return_indices=True)
    rell = campo[ix, iy, iz]
    idx = np.rint(np.asarray(malla.points, float) / sp).astype(int)
    for k in range(3):
        np.clip(idx[:, k], 0, campo.shape[k] - 1, out=idx[:, k])
    malla["v"] = rell[idx[:, 0], idx[:, 1], idx[:, 2]]
    lo, hi = F._cubo(solido, sp)
    p = F._plotter(LADO)
    try:
        p.add_mesh(malla, scalars="v", cmap=cmap, clim=clim,
                   show_scalar_bar=False, smooth_shading=True, specular=0.15)
        F._camara(p, lo, hi, "iso")
        p.screenshot(str(ruta))
    finally:
        p.close()
    return ruta


def guardar(fig, nombre):
    FIGS.mkdir(exist_ok=True)
    fig.savefig(FIGS / f"{nombre}.png", dpi=220, facecolor="white",
                bbox_inches="tight")
    fig.savefig(FIGS / f"{nombre}.pdf", facecolor="white",
                bbox_inches="tight")
    plt.close(fig)


def filas():
    return [(v, e) for v in VOIS for e in ESTR
            if (CAMPOS / f"{v}_{e}.npz").exists()]


def diferencia(c, k, m):
    """log10(|campo_k - referencia| / p99 referencia) por elemento."""
    ref = c["referencia"][m]
    p99 = float(np.percentile(ref, 99))
    d = np.abs(c[k][m] - ref) / p99
    return np.log10(np.maximum(d, 10.0 ** PISO))


def fig_mapas(m, datos):
    etiq, cmap = MAG[m]
    fr = [(v, e) for v, e in filas() if (v, e) in datos]
    nf, nc = len(fr), len(MOT)
    fig, ejes = plt.subplots(nf, nc, figsize=(2.15 * nc + 0.7, 2.15 * nf + 0.8),
                             facecolor="white", squeeze=False)
    for i, (v, e) in enumerate(fr):
        d = datos[(v, e)]
        ref = d["campos"].get("referencia") or next(iter(d["campos"].values()))
        clim = (float(np.percentile(ref[m], 1)),
                float(np.percentile(ref[m], 99)))
        sp = d["geo"]["spacing"].astype(float)
        for j, k in enumerate(MOT):
            ax = ejes[i, j]
            ax.set_axis_off()
            if k in d["campos"]:
                r = render(malla3d(d["geo"], d["campos"][k][m]), sp, clim,
                           cmap, REND / f"{m}_{v}_{e}_{k}.png")
                ax.imshow(imread(r))
            else:
                ax.text(0.5, 0.5, "sin resultado", ha="center", va="center",
                        transform=ax.transAxes, fontsize=7, color="#888")
            if i == 0:
                ax.set_title(NOM[k], fontsize=8.5, color=F.TINTA)
            if j == 0:
                ax.text(-0.02, 0.5, f"{v} · {NOM_E[e]}",
                        transform=ax.transAxes, rotation=90, ha="right",
                        va="center", fontsize=8, color=F.COLOR[e],
                        fontweight="bold")
            if j == nc - 1:
                cax = ax.inset_axes([1.03, 0.1, 0.05, 0.8])
                b = fig.colorbar(cm.ScalarMappable(norm=Normalize(*clim),
                                                   cmap=cmap), cax=cax)
                b.ax.tick_params(labelsize=6, colors=F.TINTA_2)
                b.outline.set_edgecolor(F.EJE)
    fig.suptitle(f"{etiq}: protocolo de Tapia et al. (40³, empotrado, "
                 f"18 GPa, 100 N). Una escala por fila (p1 a p99 de la "
                 f"referencia).", fontsize=9, color=F.TINTA, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    guardar(fig, f"mapa_{m}_solvers")


def fig_dif(m, datos):
    etiq, _ = MAG[m]
    fr = [(v, e) for v, e in filas() if (v, e) in datos]
    ks = MOT
    nf, nc = len(fr), len(ks)
    fig, ejes = plt.subplots(nf, nc, figsize=(2.15 * nc + 0.9, 2.15 * nf + 1.1),
                             facecolor="white", squeeze=False)
    for i, (v, e) in enumerate(fr):
        d = datos[(v, e)]
        sp = d["geo"]["spacing"].astype(float)
        for j, k in enumerate(ks):
            ax = ejes[i, j]
            ax.set_axis_off()
            if k in d["campos"] and "referencia" in d["campos"]:
                r = render(malla3d(d["geo"], diferencia(d["campos"], k, m)),
                           sp, (PISO, -2.0), "cividis",
                           REND / f"d{m}_{v}_{e}_{k}.png")
                ax.imshow(imread(r))
            if i == 0:
                ax.set_title(NOM[k], fontsize=8.5, color=F.TINTA)
            if j == 0:
                ax.text(-0.02, 0.5, f"{v} · {NOM_E[e]}",
                        transform=ax.transAxes, rotation=90, ha="right",
                        va="center", fontsize=8, color=F.COLOR[e],
                        fontweight="bold")
    cax = fig.add_axes([0.25, 0.035, 0.5, 0.012])
    b = fig.colorbar(cm.ScalarMappable(norm=Normalize(PISO, -2.0),
                                       cmap="cividis"), cax=cax,
                     orientation="horizontal")
    b.set_label("log₁₀ ( |motor − referencia| / p99 de la referencia )",
                fontsize=7.5, color=F.TINTA_2)
    b.ax.tick_params(labelsize=7, colors=F.TINTA_2)
    fig.suptitle(f"{etiq}: diferencia por elemento frente a NGSolve (directo). "
                 f"Por debajo de 10^{PISO:.0f} se trunca.", fontsize=9,
                 color=F.TINTA, y=0.995)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    guardar(fig, f"dif_{m}_solvers")


def estadisticas(datos):
    """{(v,e,k,m): (max, p99, mediana)} de |motor - ref| / p99 ref."""
    out = {}
    for (v, e), d in datos.items():
        if "referencia" not in d["campos"]:
            continue
        for k in MOT:
            if k not in d["campos"]:
                continue
            for m in ("vm", "u"):
                ref = d["campos"]["referencia"][m]
                p99 = float(np.percentile(ref, 99))
                x = np.abs(d["campos"][k][m] - ref) / p99
                out[(v, e, k, m)] = (float(x.max()),
                                     float(np.percentile(x, 99)),
                                     float(np.median(x)))
    return out


def fig_resumen(datos):
    est = estadisticas(datos)
    fr = [f for f in filas() if f in datos]
    fig, ejes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True,
                             facecolor="white")
    for ax, m in zip(ejes, ("vm", "u")):
        x = np.arange(len(fr))
        for j, k in enumerate(MOT):
            y = [max(est.get((v, e, k, m), (np.nan,))[0], 1e-17)
                 if (v, e, k, m) in est else np.nan for v, e in fr]
            ax.plot(x + (j - 2) * 0.09, y, MARCA[k], color=COLOR[k], ms=6,
                    label=NOM[k], mec="white", mew=0.5)
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{v}\n{NOM_E[e]}" for v, e in fr], fontsize=8)
        ax.set_title(MAG[m][0], fontsize=10, color=F.TINTA)
        ax.grid(True, axis="y", color=F.REJILLA, lw=0.7, which="both")
        ax.set_axisbelow(True)
    ejes[0].set_ylabel("máx |motor − referencia| / p99 de la referencia",
                       fontsize=9, color=F.TINTA)
    ejes[0].legend(frameon=False, fontsize=8, ncol=1)
    fig.suptitle("Mayor diferencia por elemento de cada motor frente a "
                 "NGSolve (directo), protocolo de Tapia", fontsize=10,
                 color=F.TINTA)
    fig.tight_layout()
    guardar(fig, "resumen_diferencias")


def fig_ecdf(datos):
    fr = [f for f in filas() if f in datos and "referencia" in datos[f]["campos"]]
    fig, ejes = plt.subplots(2, len(fr), figsize=(2.6 * len(fr) + 0.6, 5.4),
                             sharex=True, sharey=True, facecolor="white",
                             squeeze=False)
    for i, m in enumerate(("vm", "u")):
        for j, (v, e) in enumerate(fr):
            ax = ejes[i, j]
            d = datos[(v, e)]
            for k in MOT:
                if k not in d["campos"]:
                    continue
                x = np.sort(diferencia(d["campos"], k, m))
                ax.plot(x, np.linspace(0, 1, len(x)), color=COLOR[k], lw=1.4,
                        label=NOM[k])
            ax.grid(True, color=F.REJILLA, lw=0.6)
            if i == 0:
                ax.set_title(f"{v} · {NOM_E[e]}", fontsize=8.5,
                             color=F.COLOR[e])
            if j == 0:
                ax.set_ylabel(MAG[m][0] + "\nfracción de elementos",
                              fontsize=8)
            if i == 1:
                ax.set_xlabel("log₁₀ (|Δ| / p99)", fontsize=8)
    ejes[0, 0].legend(frameon=False, fontsize=7)
    fig.suptitle("Distribución acumulada de la diferencia por elemento "
                 "frente a NGSolve (directo)", fontsize=10, color=F.TINTA)
    fig.tight_layout()
    guardar(fig, "distribucion_diferencias")


def tablas(datos):
    filas_j = [json.loads(x) for x in
               (AQUI / "resultados" / "tapia.jsonl").read_text().splitlines()
               if x.strip()]
    ult = {}
    for f in filas_j:
        ult[(f["voi"], f["estructura"], f["clave"])] = f
    est = estadisticas(datos)
    md = ["## Resultados por motor\n"]
    for v, e in [f for f in filas() if f in datos]:
        md.append(f"\n### {v} · {NOM_E[e]}\n")
        md.append("| Motor | Resolvedor | E_app (MPa) | p99 vM (MPa) | máx vM "
                  "(MPa) | |u| máx (mm) | Tiempo motor (s) | Memoria (GB) |"
                  " máx Δσ | máx Δ|u| |")
        md.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for k in MOT + ("referencia",):
            f = ult.get((v, e, k))
            if not f:
                continue
            c = (f.get("corridas") or [{}])[0]
            nombre = NOM.get(k, "NGSolve (directo, ref.)")
            if "E_app_MPa" not in f:
                md.append(f"| {nombre} | FALLO | | | | | | | | |")
                continue
            sv = est.get((v, e, k, "vm"))
            su = est.get((v, e, k, "u"))
            md.append(
                f"| {nombre} | {str(c.get('solver'))[:34]} | "
                f"{f['E_app_MPa']:.3f} | {f['vm_p99_MPa']:.3f} | "
                f"{f['vm_max_MPa']:.2f} | {f['u_max_mm']:.5f} | "
                f"{c.get('tiempo_s') or 0:.1f} | "
                f"{(c.get('memoria_MB') or 0) / 1024:.2f} | "
                f"{'—' if k == 'referencia' or not sv else f'{sv[0]:.1e}'} | "
                f"{'—' if k == 'referencia' or not su else f'{su[0]:.1e}'} |")
    (AQUI / "tablas_tapia.md").write_text("\n".join(md) + "\n",
                                          encoding="utf-8")


def main():
    datos = {f: cargar(*f) for f in filas()}
    tablas(datos)
    for m in ("vm", "u"):
        fig_mapas(m, datos)
        fig_dif(m, datos)
    fig_resumen(datos)
    fig_ecdf(datos)
    print("listo")


if __name__ == "__main__":
    main()
