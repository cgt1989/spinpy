"""
figuras.py: figuras y tablas de la presentacion del protocolo de Tapia con
C1 a C5, sus spinodoides y sus dual-lattice.

    xvfb-run -a python figuras.py

Solo LEE `campos/`, `resultados/` y `../tapia/estructuras/`. Escribe en `figs/`:
  bloque.png         cubo macizo: von Mises/sigma y |u| en un corte, apoyo
                     deslizante frente a empotrado (base fija) con el mismo
                     protocolo
  mapa_vm.png        sigma von Mises [MPa], 3 estructuras x 5 VOIs
  mapa_u.png         deformacion total |u| [mm]
  resumen.json       cifras de la tabla de resultados y de los ajustes
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent.parent))
sys.path.insert(0, str(AQUI.parent / "tapia"))

import matplotlib                                               # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                 # noqa: E402
from matplotlib import cm                                       # noqa: E402
from matplotlib.colors import Normalize                         # noqa: E402
from matplotlib.image import imread                             # noqa: E402

import figuras_tapia as FT                                      # noqa: E402
from spinpy import fem                                          # noqa: E402

FIGS = AQUI / "figs"
REND = FIGS / "render"
CAMPOS = AQUI / "campos"
EST_DIR = AQUI.parent / "tapia" / "estructuras"
VOIS = ("C1", "C2", "C3", "C4", "C5")
ESTR = ("voi", "spinodoide", "dual-lattice")
NOM = {"voi": "VOI", "spinodoide": "Spinodoide", "dual-lattice": "Dual-lattice"}
MAG = {"vm": ("Tensión de von Mises [MPa]", "inferno"),
       "u": ("Deformación total |u| [mm]", "viridis")}
FT.LADO = 600
TINTA, TINTA2 = "#1d2733", "#5b6573"
E_S, NU, F_N, L_MM, N = 18e9, 0.30, 100.0, 3.008, 40


def guardar(fig, nombre):
    FIGS.mkdir(exist_ok=True)
    fig.savefig(FIGS / f"{nombre}.png", dpi=200, facecolor="white",
                bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Cubo macizo
# ---------------------------------------------------------------------------

def fig_bloque():
    sig = F_N / L_MM ** 2                                 # MPa
    h = L_MM / N
    BW = np.ones((N, N, N), bool)
    datos = {}
    for apoyo in ("deslizante", "empotrado"):
        malla = fem.mallar(BW, np.full(3, h), "hex8")
        prot = fem.protocolo("tapia2026", apoyo=apoyo)
        reg = fem.ensayo(malla, prot, analisis=("lineal",), motor="ngsolve",
                         solver="directo", aislado=False)
        u = reg["_campos_lineal"]["u"]
        s = fem.tension_elemental(malla["nodos"], malla["elems"], u, E_S, NU)
        vm = fem.von_mises(s) * 1e-6
        um = np.linalg.norm(u, axis=1)[malla["elems"]].mean(axis=1)
        cen = malla["nodos"][malla["elems"]].mean(axis=1)
        ijk = np.floor(cen / h).astype(int)
        gv = np.full((N, N, N), np.nan)
        gu = np.full((N, N, N), np.nan)
        gv[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = vm
        gu[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = um
        datos[apoyo] = (gv, gu)
    fig, ejes = plt.subplots(2, 2, figsize=(6.6, 6.4), facecolor="white")
    clim_v = (0.95, 1.30)
    umax = max(np.nanmax(datos[a][1]) for a in datos) * 1e3
    clim_u = (0.0, umax)
    for j, apoyo in enumerate(("deslizante", "empotrado")):
        gv, gu = datos[apoyo]
        a = ejes[0, j]
        im = a.imshow((gv[:, N // 2, :] / sig).T, origin="lower",
                      extent=(0, L_MM, 0, L_MM), cmap="inferno",
                      vmin=clim_v[0], vmax=clim_v[1])
        a.set_title({"deslizante": "Apoyo deslizante\n(uz = 0 en la base)",
                     "empotrado": "Base empotrada (Tapia)\n(ux = uy = uz = 0)"}[
            apoyo], fontsize=9.5, color=TINTA)
        b = ejes[1, j]
        im2 = b.imshow(gu[:, N // 2, :].T * 1e3, origin="lower",
                       extent=(0, L_MM, 0, L_MM), cmap="viridis",
                       vmin=clim_u[0], vmax=clim_u[1])
        for ax in (a, b):
            ax.set_xlabel("x [mm]", fontsize=8, color=TINTA2)
            ax.tick_params(labelsize=7, colors=TINTA2)
        a.set_ylabel("z [mm]  (base abajo, carga arriba)", fontsize=8,
                     color=TINTA2)
        b.set_ylabel("z [mm]", fontsize=8, color=TINTA2)
        if j == 1:
            for ax, mm, tit in ((a, im, "von Mises / σ"),
                                (b, im2, "|u| [µm]")):
                cb = fig.colorbar(mm, ax=ax, fraction=0.046, pad=0.03)
                cb.ax.tick_params(labelsize=7)
                cb.set_label(tit, fontsize=8)
    fig.suptitle("Cubo macizo de 3,008 mm, 100 N, E = 18 GPa. Corte en y = L/2",
                 fontsize=10, color=TINTA)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    guardar(fig, "bloque")


# ---------------------------------------------------------------------------
# Mapas de las 15 estructuras
# ---------------------------------------------------------------------------

def cargar(v, est):
    z = np.load(CAMPOS / f"{v}_{est}.npz")
    return {k: z[k] for k in z.files}


def fig_mapas(m):
    etiq, cmap = MAG[m]
    datos = {(v, e): cargar(v, e) for v in VOIS for e in ESTR
             if (CAMPOS / f"{v}_{e}.npz").exists()}
    fig, ejes = plt.subplots(3, 5, figsize=(13.2, 7.4), facecolor="white")
    for j, v in enumerate(VOIS):
        tres = [datos[(v, e)][f"{m}_referencia"] for e in ESTR
                if (v, e) in datos]
        if not tres:
            continue
        lo = 0.0
        hi = float(max(np.percentile(t, 99) for t in tres))
        clim = (lo, hi)
        for i, e in enumerate(ESTR):
            ax = ejes[i, j]
            ax.set_axis_off()
            if (v, e) not in datos:
                continue
            d = datos[(v, e)]
            sp = d["spacing"].astype(float)
            r = FT.render(FT.malla3d(d, d[f"{m}_referencia"]), sp, clim, cmap,
                          REND / f"{m}_{v}_{e}.png")
            ax.imshow(imread(r))
            if j == 0:
                ax.text(-0.02, 0.5, NOM[e], transform=ax.transAxes,
                        rotation=90, ha="right", va="center", fontsize=11,
                        color=TINTA, fontweight="bold")
            if i == 0:
                ax.set_title(v, fontsize=13, color=TINTA, fontweight="bold")
        cax = fig.add_axes([ejes[2, j].get_position().x0 + 0.01,
                            ejes[2, j].get_position().y0 - 0.035,
                            ejes[2, j].get_position().width - 0.02, 0.014])
        cb = fig.colorbar(cm.ScalarMappable(norm=Normalize(*clim), cmap=cmap),
                          cax=cax, orientation="horizontal")
        cb.ax.tick_params(labelsize=8, colors=TINTA2)
        cb.set_label(etiq if j == 2 else "", fontsize=9, color=TINTA)
    fig.subplots_adjust(left=0.05, right=0.995, top=0.95, bottom=0.09,
                        wspace=0.02, hspace=0.04)
    guardar(fig, f"mapa_{m}")


# ---------------------------------------------------------------------------
# Cifras
# ---------------------------------------------------------------------------

def resumen():
    filas = [json.loads(x) for x in
             open(AQUI / "resultados" / "tapia_c1_c5.jsonl", encoding="utf-8")]
    ref = {(f["voi"], f["estructura"]): f for f in filas
           if f["clave"] == "referencia" and "E_app_MPa" in f}
    app = {(f["voi"], f["estructura"]): f for f in filas
           if f["clave"] == "app" and "E_app_MPa" in f}
    tabla, dif = [], []
    for v in VOIS:
        for e in ESTR:
            r = ref.get((v, e))
            if not r:
                continue
            a = app.get((v, e))
            d = {"voi": v, "estructura": e, "BVTV": r["BVTV_malla"],
                 "E_app_MPa": r["E_app_MPa"], "vm_p99_MPa": r["vm_p99_MPa"],
                 "vm_max_MPa": r["vm_max_MPa"],
                 "u_tot_max_mm": r["u_tot_max_mm"],
                 "uz_max_mm": r["uz_max_mm"]}
            if a:
                d["dif_E_app_app_vs_ngsolve"] = abs(
                    a["E_app_MPa"] / r["E_app_MPa"] - 1)
                d["dif_vm_max_app_vs_ngsolve"] = abs(
                    a["vm_max_MPa"] / r["vm_max_MPa"] - 1)
                dif.append(d["dif_E_app_app_vs_ngsolve"])
            tabla.append(d)
    ajustes = []
    for v in VOIS:
        for fam in ("spinodoide", "dual-lattice"):
            p = EST_DIR / f"{v}_{fam}.json"
            if p.exists():
                j = json.load(open(p, encoding="utf-8"))
                ajustes.append({"voi": v, "familia": fam,
                                "error": j["error"],
                                "BVTV_generada": j["BVTV_generada"],
                                "BVTV_voi": j["BVTV_voi"],
                                "tiempo_s": j["tiempo_ajuste_s"]})
    json.dump({"tabla": tabla, "ajustes": ajustes,
               "dif_max_app_vs_ngsolve": max(dif) if dif else None},
              open(AQUI / "figs" / "resumen.json", "w", encoding="utf-8"),
              indent=1, default=float)
    print(f"{len(tabla)} estructuras en la tabla; {len(ajustes)} ajustes")


if __name__ == "__main__":
    FIGS.mkdir(exist_ok=True)
    fig_bloque()
    for m in ("vm", "u"):
        fig_mapas(m)
    resumen()
