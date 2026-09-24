"""figuras_porcino.py — Figuras de la validacion mecanica frente a FEBio.

    cd Port_Python
    python comparativa_febio/porcino/figuras_porcino.py            # es y en
    python comparativa_febio/porcino/figuras_porcino.py --idioma es

Solo LEE: `resultados/*.jsonl`, `campos/*.npz` y, para la respuesta no lineal a
la carga de fallo, el logfile de FEBio que dejo `validar_porcino.py` en
`corridas/`. No resuelve nada.

Los mapas 3D se pintan con las MISMAS piezas que la figura 8 del informe de la
app (`figuras.malla_figura`, `_camara`, `_cubo`, el mapa `inferno`, proyeccion
paralela, cada vertice con el valor del voxel solido mas proximo), para que el
mapa «de la app» sea el que la app dibuja y el de FEBio se pinte con la misma
receta y la MISMA escala de color.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent.parent))
sys.path.insert(0, str(AQUI.parent))

from comparar_febio import _ultimo_registro, _von_mises  # noqa: E402
from spinpy import figuras as F                         # noqa: E402

DIR_RES = AQUI / "resultados"
DIR_CAMPOS = AQUI / "campos"
DIR_CORR = AQUI / "corridas"
DIR_FIGS = AQUI / "figs"
DIR_REND = AQUI / "figs" / "render"
N_CONV = (22, 28, 34, 40)
ORDEN = F.ORDEN
COLOR = F.COLOR
DPI = 300
LADO = 900

T = {
    "es": dict(
        voi="VOI", spinodoide="Spinodoide", dl="Dual-lattice",
        app="spinpy (app)", febio="FEBio", febio_lin="FEBio, lineal",
        febio_nl="FEBio, no lineal", dif_lin="|FEBio − app|",
        dif_nl="FEBio no lineal − app", vm="σ von Mises [MPa]",
        desp="deformación total |u| [mm]",
        dif_rel="log₁₀ |Δ| / p99 del campo", dif_pct="Δ / p99 del campo [%]",
        compresion="Ensayo de la app (40³, deslizante, 20 GPa, 1 MPa)",
        comparado="Protocolo de Tapia et al. (40³, empotrado, 18 GPa, 100 N)",
        homogeneizacion="Tensor periódico (32³)",
        convergencia="Convergencia (22–40³)", fallo="Fallo progresivo (32³)",
        dif_eje="|diferencia relativa| spinpy frente a FEBio",
        tol="tolerancia declarada", paso="paso", E_rel="E_app / E_app,0",
        F_fallo="F de fallo de Pistoia [N]", n="n (vóxeles por lado)",
        E_app="E_app [MPa]", excedencia="1 − F  (fracción de la capa superficial)",
        rig_plato="plato rígido", nl_carga="no lineal, carga del ensayo",
        nl_fallo="no lineal, carga de Pistoia", lineal="lineal",
        const="constante de ingeniería [GPa]",
        paridad_x="σ von Mises spinpy [MPa]",
        paridad_y="σ von Mises FEBio [MPa]",
        metricas={"E": "E_app", "u": "campo u", "sig": "campo σ",
                  "vm": "von Mises por elemento", "desp": "|u| por elemento",
                  "p99": "p99 de superficie", "pist": "factor de Pistoia",
                  "F": "fuerza (equilibrio)", "C": "tensor C (máx.)",
                  "LE": "distancia log-euclídea", "Econv": "E_app (serie)",
                  "Efallo": "E_app por paso", "Ffallo": "F de fallo por paso"},
    ),
    "en": dict(
        voi="VOI", spinodoide="Spinodoid", dl="Dual-lattice",
        app="spinpy (app)", febio="FEBio", febio_lin="FEBio, linear",
        febio_nl="FEBio, nonlinear", dif_lin="|FEBio − app|",
        dif_nl="FEBio nonlinear − app", vm="von Mises σ [MPa]",
        desp="total deformation |u| [mm]",
        dif_rel="log₁₀ |Δ| / field p99", dif_pct="Δ / field p99 [%]",
        compresion="App test (40³, sliding, 20 GPa, 1 MPa)",
        comparado="Tapia et al. protocol (40³, fixed, 18 GPa, 100 N)",
        homogeneizacion="Periodic tensor (32³)",
        convergencia="Convergence (22–40³)", fallo="Progressive failure (32³)",
        dif_eje="|relative difference| spinpy vs FEBio",
        tol="declared tolerance", paso="step", E_rel="E_app / E_app,0",
        F_fallo="Pistoia failure load [N]", n="n (voxels per side)",
        E_app="E_app [MPa]", excedencia="1 − F  (fraction of surface layer)",
        rig_plato="rigid platen", nl_carga="nonlinear, test load",
        nl_fallo="nonlinear, Pistoia load", lineal="linear",
        const="engineering constant [GPa]",
        paridad_x="spinpy von Mises σ [MPa]",
        paridad_y="FEBio von Mises σ [MPa]",
        metricas={"E": "E_app", "u": "u field", "sig": "σ field",
                  "vm": "element von Mises", "desp": "element |u|",
                  "p99": "surface p99", "pist": "Pistoia factor",
                  "F": "force (equilibrium)", "C": "tensor C (max.)",
                  "LE": "log-Euclidean distance", "Econv": "E_app (series)",
                  "Efallo": "E_app per step", "Ffallo": "failure F per step"},
    ),
}


def nombre(est, t):
    return {"voi": t["voi"], "spinodoide": t["spinodoide"],
            "dual-lattice": t["dl"]}[est]


def leer(analisis):
    ruta = DIR_RES / f"{analisis}.jsonl"
    if not ruta.exists():
        return []
    filas = [json.loads(l) for l in open(ruta, encoding="utf-8") if l.strip()]
    # Si un analisis se repitio, vale la ULTIMA corrida de cada clave.
    clave = (lambda f: (f["estructura"], f.get("n"), f.get("paso")))
    ult = {}
    for f in filas:
        ult[clave(f)] = f
    return list(ult.values())


def _leyenda_abajo(fig, handles, labels, ncol, dy=0.0):
    """Leyenda FUERA de los ejes, centrada bajo la figura.

    Dentro de los ejes tapaba datos (tolerancias en la figura de coincidencia,
    la serie del dual-lattice en la de convergencia). `bbox_inches="tight"` al
    guardar amplia el lienzo para que quepa.
    """
    fig.legend(handles, labels, loc="upper center", ncol=ncol,
               bbox_to_anchor=(0.5, dy), fontsize=6.5, frameon=False,
               handlelength=1.8, columnspacing=1.4)


def _guardar(fig, destino):
    destino = Path(destino)
    destino.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destino.with_suffix(".png"), dpi=DPI, facecolor=F.FONDO,
                bbox_inches="tight")
    fig.savefig(destino.with_suffix(".pdf"), facecolor=F.FONDO,
                bbox_inches="tight")
    return destino.with_suffix(".png")


# ---------------------------------------------------------------------------
# Mapas 3D
# ---------------------------------------------------------------------------

def campo3d(z, valores):
    c = np.full(tuple(int(x) for x in z["forma"]), np.nan)
    ijk = z["ijk"].astype(int)
    c[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = valores
    return c


def render(campo, sp, clim, cmap, ruta):
    """Una estructura coloreada por `campo` (NaN fuera del hueso), con la
    receta de `figuras.renders_von_mises` y una escala fija."""
    from scipy import ndimage
    ruta = Path(ruta)
    if ruta.exists():
        return ruta
    ruta.parent.mkdir(parents=True, exist_ok=True)
    solido = np.isfinite(campo)
    malla = F.malla_figura(solido, sp, suavizar=False)
    _d, (ix, iy, iz) = ndimage.distance_transform_edt(
        ~solido, return_indices=True)
    relleno = campo[ix, iy, iz]
    idx = np.rint(np.asarray(malla.points, float) / sp).astype(int)
    for k in range(3):
        np.clip(idx[:, k], 0, campo.shape[k] - 1, out=idx[:, k])
    malla["v"] = relleno[idx[:, 0], idx[:, 1], idx[:, 2]]
    lo, hi = F._cubo(solido, sp)
    p = F._plotter(LADO)
    try:
        p.add_mesh(malla, scalars="v", cmap=cmap, clim=clim,
                   show_scalar_bar=False, smooth_shading=True, specular=0.15)
        try:
            p.enable_anti_aliasing("ssaa")
        except Exception:
            pass
        F._camara(p, lo, hi, "iso")
        p.screenshot(str(ruta))
    finally:
        p.close()
    return ruta


def _nl_fallo_vm(est):
    """von Mises no lineal a la carga de fallo, reescalado a 1 MPa."""
    base = DIR_CORR / "compresion" / est / f"compresion_{est}_nl_fallo_s.txt"
    if not base.exists():
        return None, None
    fila = next((f for f in leer("compresion") if f["estructura"] == est),
                None)
    if fila is None or not fila.get("nl_fallo", {}).get("converge"):
        return None, None
    s = _ultimo_registro(base, 6) * 1e6
    esc = fila["sigma_app_Pa"] / fila["nl_fallo"]["sigma_Pa"]
    return _von_mises(s) * esc, fila["nl_fallo"]["sigma_Pa"]


def panel_mapas(analisis, magnitud, idioma):
    """Filas: estructuras. Columnas: app | FEBio lineal | |Δ| lineal |
    FEBio no lineal | Δ no lineal. Escala comun en las columnas de valor."""
    from matplotlib import cm
    from matplotlib.colors import Normalize
    from matplotlib.image import imread
    t = T[idioma]
    datos = {}
    for est in ORDEN:
        ruta = DIR_CAMPOS / f"{analisis}_{est}.npz"
        if ruta.exists():
            datos[est] = dict(np.load(ruta))
    if not datos:
        return None
    esc = 1e-6 if magnitud == "vm" else 1.0
    col_nl = None
    todos = []
    series = {}
    for est, z in datos.items():
        a = z[f"{magnitud}_spinpy"].astype(float) * esc
        b = z[f"{magnitud}_febio"].astype(float) * esc
        if magnitud == "vm" and analisis == "compresion":
            nl, sg = _nl_fallo_vm(est)
            nl = None if nl is None else nl * esc
            col_nl = ("fallo", sg)
        else:
            nl = z.get(f"{magnitud}_febio_nl")
            nl = None if nl is None else nl.astype(float) * esc
            col_nl = ("ensayo", None)
        series[est] = (a, b, nl)
        todos.append(a)
    todos = np.concatenate(todos)
    clim = (float(np.percentile(todos, 1)), float(np.percentile(todos, 99)))
    cmap = F.CMAP_TENSION if magnitud == "vm" else "viridis"
    lim_nl = 0.0
    for est, (a, b, nl) in series.items():
        if nl is not None:
            ref = np.percentile(a, 99)
            lim_nl = max(lim_nl, float(np.percentile(
                np.abs(nl - a) / ref * 100, 99.5)))
    lim_nl = lim_nl or 1.0
    rutas = {}
    for est, (a, b, nl) in series.items():
        z = datos[est]
        sp = z["spacing"].astype(float)
        ref = float(np.percentile(a, 99))
        dlin = np.log10(np.maximum(np.abs(b - a) / ref, 1e-16))
        pre = DIR_REND / f"{analisis}_{magnitud}_{est}"
        r = [render(campo3d(z, a), sp, clim, cmap, f"{pre}_app.png"),
             render(campo3d(z, b), sp, clim, cmap, f"{pre}_febio.png"),
             render(campo3d(z, dlin), sp, (-12, -4), "cividis",
                    f"{pre}_dlin.png")]
        if nl is not None:
            r.append(render(campo3d(z, nl), sp, clim, cmap,
                            f"{pre}_nl_{col_nl[0]}.png"))
            r.append(render(campo3d(z, (nl - a) / ref * 100), sp,
                            (-lim_nl, lim_nl), "RdBu_r",
                            f"{pre}_dnl_{col_nl[0]}.png"))
        rutas[est] = r
    ncol = max(len(r) for r in rutas.values())
    nfil = len(rutas)
    fig = F._figura(2.05 * ncol, 2.05 * nfil + 0.9)
    tit = [t["app"], t["febio_lin"], t["dif_lin"]]
    if ncol > 3:
        extra = ""
        if col_nl[0] == "fallo":
            extra = (" (σ de fallo)" if idioma == "es"
                     else " (failure σ)")
        tit += [t["febio_nl"] + extra, t["dif_nl"]]
    ancho = 0.9 / ncol
    alto = 0.84 / nfil
    for i, est in enumerate(rutas):
        for j, ruta in enumerate(rutas[est]):
            ax = fig.add_axes([0.02 + j * ancho, 0.14 + (nfil - 1 - i) * alto,
                               ancho * 0.98, alto * 0.98])
            ax.imshow(imread(ruta))
            ax.set_axis_off()
            if i == 0:
                ax.set_title(tit[j], fontsize=7.5, color=F.TINTA)
            if j == 0:
                ax.text(-0.02, 0.5, nombre(est, t), transform=ax.transAxes,
                        rotation=90, ha="right", va="center", fontsize=8.5,
                        color=COLOR[est], fontweight="bold")
    barras = [((0, 2), Normalize(*clim), cmap,
               t["vm"] if magnitud == "vm" else t["desp"]),
              ((2, 3), Normalize(-12, -4), "cividis", t["dif_rel"])]
    if ncol > 3:
        barras.append(((4, 5), Normalize(-lim_nl, lim_nl), "RdBu_r",
                       t["dif_pct"]))
    for (j0, j1), norm, cm_, et in barras:
        cax = fig.add_axes([0.02 + j0 * ancho + 0.1 * ancho, 0.075,
                            (j1 - j0) * ancho - 0.2 * ancho, 0.022])
        b = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm_), cax=cax,
                         orientation="horizontal")
        b.outline.set_edgecolor(F.EJE)
        b.ax.tick_params(labelsize=6.5, colors=F.TINTA_2)
        b.set_label(et, fontsize=7, color=F.TINTA_2)
    if ncol > 3 and col_nl[0] == "fallo":
        # la columna no lineal se muestra reescalada a la carga del ensayo
        pass
    fig.suptitle(t[analisis], fontsize=9, color=F.TINTA, y=1.035)
    return _guardar(fig, DIR_FIGS / idioma
                    / f"mapa_{magnitud}_{analisis}")


# ---------------------------------------------------------------------------
# Graficos de datos
# ---------------------------------------------------------------------------

TOL = {"E": 1e-6, "u": 1e-5, "sig": 1e-4, "vm": 1e-4, "desp": 1e-5,
       "p99": 1e-5, "pist": 1e-5, "F": 1e-9, "C": 1e-6, "LE": 1e-6,
       "Econv": 1e-6, "Efallo": 1e-6, "Ffallo": 1e-5}


def fig_coincidencia(idioma):
    """Todas las diferencias spinpy/FEBio del problema lineal en una figura."""
    t = T[idioma]
    filas = []   # (grupo, metrica, estructura, valor)
    for an in ("compresion", "comparado"):
        for f in leer(an):
            for m, k in (("E", "dE_rel"), ("u", "du_max_rel"),
                         ("sig", "dsig_max_rel"), ("vm", "dvm_max_elem_rel"),
                         ("desp", "ddesp_max_rel"),
                         ("p99", "d_vm_p99_superficie_rel"),
                         ("pist", "d_factor_rel"), ("F", "dF_rel")):
                if f.get(k) is not None:
                    filas.append((an, m, f["estructura"], f[k]))
    for f in leer("homogeneizacion"):
        filas.append(("homogeneizacion", "C", f["estructura"],
                      f["dC_max_rel"]))
        filas.append(("homogeneizacion", "LE", f["estructura"],
                      f["d_log_euclidea"]))
    for an, m, k in (("convergencia", "Econv", "dE_rel"),
                     ("fallo", "Efallo", "dE_rel"),
                     ("fallo", "Ffallo", "dF_rel")):
        por = {}
        for f in leer(an):
            if an == "fallo" and f.get("voxeles_intactos_distintos"):
                continue
            por.setdefault(f["estructura"], []).append(f[k])
        for est, v in por.items():
            filas.append((an, m, est, max(v)))
    if not filas:
        return None
    grupos = [g for g in ("compresion", "comparado", "homogeneizacion",
                          "convergencia", "fallo")
              if any(r[0] == g for r in filas)]
    # Cada grupo lleva su PROPIA fila de titulo: antes el titulo se escribia
    # a media fila de la primera metrica y se montaba sobre la ultima
    # etiqueta del grupo anterior.
    filas_y = []                       # ("titulo", g) o ("metrica", g, m)
    for g in grupos:
        filas_y.append(("titulo", g))
        for m in dict.fromkeys(r[1] for r in filas if r[0] == g):
            filas_y.append(("metrica", g, m))
    n = len(filas_y)
    fig = F._figura(6.6, 0.22 * n + 0.9)
    ax = fig.add_axes([0.36, 0.08, 0.6, 0.86])
    F._estilo(ax, rejilla="x")
    y = {}
    for i, f in enumerate(filas_y):
        yy = n - 1 - i
        if f[0] == "titulo":
            ax.text(-0.55, yy, t[f[1]], transform=ax.get_yaxis_transform(),
                    fontsize=7, color=F.TINTA, fontweight="bold", ha="left",
                    va="center")
            if i:
                ax.axhline(yy + 0.5, color=F.EJE, lw=0.5)
        else:
            y[(f[1], f[2])] = yy
    desp = {"voi": -0.2, "spinodoide": 0.0, "dual-lattice": 0.2}
    for g, m, est, v in filas:
        ax.plot(max(v, 1e-16), y[(g, m)] + desp[est], "o", ms=4,
                color=COLOR[est], mec="white", mew=0.4)
    for (g, m), yy in y.items():
        ax.plot([TOL[m]] * 2, [yy - 0.35, yy + 0.35], color=F.TINTA, lw=1.2)
    ax.set_xscale("log")
    ax.set_xlim(1e-14, 1e-3)
    ax.set_ylim(-0.7, n - 0.4)
    ax.set_yticks(list(y.values()))
    ax.set_yticklabels([t["metricas"][m] for g, m in y], fontsize=7)
    ax.set_xlabel(t["dif_eje"], fontsize=7.5, color=F.TINTA_2)
    from matplotlib.lines import Line2D
    hs = ([Line2D([], [], marker="o", ls="", color=COLOR[e]) for e in ORDEN]
          + [Line2D([], [], color=F.TINTA)])
    _leyenda_abajo(fig, hs, [nombre(e, t) for e in ORDEN] + [t["tol"]], 4,
                   dy=-0.01)
    return _guardar(fig, DIR_FIGS / idioma / "fig_coincidencia")


def fig_paridad(analisis, idioma):
    t = T[idioma]
    ests = [e for e in ORDEN if (DIR_CAMPOS / f"{analisis}_{e}.npz").exists()]
    if not ests:
        return None
    fig = F._figura(2.3 * len(ests), 2.35)
    for j, est in enumerate(ests):
        z = np.load(DIR_CAMPOS / f"{analisis}_{est}.npz")
        a = z["vm_spinpy"].astype(float) / 1e6
        b = z["vm_febio"].astype(float) / 1e6
        ax = fig.add_subplot(1, len(ests), j + 1)
        F._estilo(ax, rejilla="both")
        ax.hexbin(a, b, gridsize=60, bins="log", cmap="Greys", mincnt=1,
                  linewidths=0)
        lim = [0, float(max(a.max(), b.max())) * 1.02]
        ax.plot(lim, lim, color=COLOR[est], lw=0.8)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal")
        d = np.abs(b - a).max() / a.max()
        ax.set_title(f"{nombre(est, t)}  ·  máx |Δ|/máx = {d:.0e}"
                     if idioma == "es" else
                     f"{nombre(est, t)}  ·  max |Δ|/max = {d:.0e}",
                     fontsize=7.5, color=F.TINTA)
        ax.set_xlabel(t["paridad_x"], fontsize=7, color=F.TINTA_2)
        if j == 0:
            ax.set_ylabel(t["paridad_y"], fontsize=7, color=F.TINTA_2)
    fig.tight_layout()
    return _guardar(fig, DIR_FIGS / idioma / f"fig_paridad_{analisis}")


def fig_colas(idioma):
    """Excedencia de von Mises en la capa superficial, los dos programas."""
    t = T[idioma]
    ans = [a for a in ("compresion", "comparado") if leer(a)]
    if not ans:
        return None
    fig = F._figura(3.2 * len(ans), 2.6)
    for j, an in enumerate(ans):
        ax = fig.add_subplot(1, len(ans), j + 1)
        F._estilo(ax, rejilla="both")
        for f in sorted(leer(an), key=lambda f: ORDEN.index(f["estructura"])):
            q = f["cuantiles_vm_superficie"]
            p = np.asarray(q["p"], float) / 100
            ax.plot(q["spinpy_MPa"], 1 - p, color=COLOR[f["estructura"]],
                    lw=1.6, alpha=0.55, label=f"{nombre(f['estructura'], t)}"
                    f" · {t['app']}")
            ax.plot(q["febio_MPa"], 1 - p, ls="none", marker=".", ms=2.5,
                    color=COLOR[f["estructura"]],
                    label=f"{nombre(f['estructura'], t)} · FEBio")
            p99 = f["vm_p99_superficie_spinpy"] * 1e-6
            ax.axvline(p99, color=COLOR[f["estructura"]], lw=0.5, ls=":")
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 1)
        # Eje cortado en p99,9, como la figura 11 de la app: el maximo es el
        # punto que no converge y aplastaba la cola.
        corte = max(float(np.interp(99.9, f["cuantiles_vm_superficie"]["p"],
                                    f["cuantiles_vm_superficie"]["spinpy_MPa"]))
                    for f in leer(an))
        ax.set_xlim(0, 1.05 * corte)
        ax.set_xlabel(t["vm"], fontsize=7, color=F.TINTA_2)
        ax.set_ylabel(t["excedencia"], fontsize=7, color=F.TINTA_2)
        ax.set_title(t[an], fontsize=7.5, color=F.TINTA)
    fig.tight_layout()
    h, l = fig.axes[0].get_legend_handles_labels()
    _leyenda_abajo(fig, h, l, 3)
    return _guardar(fig, DIR_FIGS / idioma / "fig_colas_superficie")


def fig_rigidez(idioma):
    t = T[idioma]
    ans = [a for a in ("compresion", "comparado") if leer(a)]
    if not ans:
        return None
    fig = F._figura(3.4 * len(ans), 2.7)
    for j, an in enumerate(ans):
        ax = fig.add_subplot(1, len(ans), j + 1)
        F._estilo(ax, rejilla="y")
        filas = sorted(leer(an), key=lambda f: ORDEN.index(f["estructura"]))
        series = [(t["app"], lambda f: f["E_app_spinpy_MPa"], 1.0, None),
                  (t["febio_lin"], lambda f: f["E_app_febio_MPa"], 0.75,
                   None),
                  (t["nl_carga"], lambda f: (f.get("nl_ensayo") or {}).get(
                      "E_app_MPa"), 0.5, "//")]
        if an == "compresion":
            series += [(t["nl_fallo"], lambda f: (f.get("nl_fallo") or {}
                                                  ).get("E_app_MPa"), 0.35,
                        "xx"),
                       (t["rig_plato"], lambda f: f.get("E_app_plato_MPa"),
                        0.2, "..")]
        w = 0.8 / len(series)
        for k, (et, fn, alfa, rayado) in enumerate(series):
            for i, f in enumerate(filas):
                v = fn(f)
                if v is None:
                    continue
                ax.bar(i - 0.4 + w * (k + 0.5), v, w * 0.92,
                       color=COLOR[f["estructura"]], alpha=alfa,
                       hatch=rayado, edgecolor="white", lw=0.3,
                       label=et if i == 0 else None)
        ax.set_xticks(range(len(filas)))
        ax.set_xticklabels([nombre(f["estructura"], t) for f in filas],
                           fontsize=7)
        ax.set_ylabel(t["E_app"], fontsize=7, color=F.TINTA_2)
        ax.set_title(t[an], fontsize=7.5, color=F.TINTA)
    fig.tight_layout()
    # Una sola leyenda, con las series del panel mas completo; en gris,
    # porque el color ya es la estructura.
    from matplotlib.patches import Patch
    h, l = fig.axes[0].get_legend_handles_labels()
    gris = [Patch(facecolor="#888888", alpha=x.patches[0].get_alpha(),
                  hatch=x.patches[0].get_hatch(), edgecolor="white")
            for x in h]
    _leyenda_abajo(fig, gris, l, len(l))
    return _guardar(fig, DIR_FIGS / idioma / "fig_rigidez")


def fig_tensor(idioma):
    t = T[idioma]
    filas = sorted(leer("homogeneizacion"),
                   key=lambda f: ORDEN.index(f["estructura"]))
    if not filas:
        return None
    claves = ("Ex", "Ey", "Ez", "Gyz", "Gxz", "Gxy")
    fig = F._figura(6.4, 2.6)
    ax = fig.add_axes([0.08, 0.16, 0.6, 0.74])
    F._estilo(ax, rejilla="y")
    w = 0.8 / len(filas)
    for i, f in enumerate(filas):
        for k, c in enumerate(claves):
            x = k - 0.4 + w * (i + 0.5)
            ax.bar(x, f["constantes_spinpy"][c] / 1e9, w * 0.9,
                   color=COLOR[f["estructura"]], alpha=0.8,
                   label=nombre(f["estructura"], t) if k == 0 else None)
            ax.plot(x, f["constantes_febio"][c] / 1e9, "_", ms=9, mew=1.4,
                    color=F.TINTA, label=("FEBio" if (k == 0 and i == 0)
                                          else None))
    ax.set_xticks(range(len(claves)))
    ax.set_xticklabels(["E_x", "E_y", "E_z", "G_yz", "G_xz", "G_xy"],
                       fontsize=7)
    ax.set_ylabel(t["const"], fontsize=7, color=F.TINTA_2)
    ax.legend(fontsize=6, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.12), ncol=4)
    ax.set_title(t["homogeneizacion"], fontsize=7.5, color=F.TINTA)
    # Mapa de |ΔC| / max|C| del VOI
    f = filas[0]
    Cs, Cf = np.asarray(f["C_spinpy_Pa"]), np.asarray(f["C_febio_Pa"])
    ax2 = fig.add_axes([0.74, 0.2, 0.22, 0.62])
    im = ax2.imshow(np.log10(np.maximum(np.abs(Cf - Cs) / np.abs(Cs).max(),
                                        1e-16)), cmap="cividis", vmin=-12,
                    vmax=-6)
    ax2.set_xticks(range(6))
    ax2.set_yticks(range(6))
    lab = ["xx", "yy", "zz", "yz", "xz", "xy"]
    ax2.set_xticklabels(lab, fontsize=5.5)
    ax2.set_yticklabels(lab, fontsize=5.5)
    ax2.set_title(("log₁₀ |ΔC_ij| / máx|C| · " if idioma == "es" else
                   "log₁₀ |ΔC_ij| / max|C| · ") + nombre(f["estructura"], t),
                  fontsize=6.5, color=F.TINTA)
    b = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
    b.ax.tick_params(labelsize=5.5)
    return _guardar(fig, DIR_FIGS / idioma / "fig_tensor")


def fig_convergencia(idioma):
    t = T[idioma]
    filas = leer("convergencia")
    if not filas:
        return None
    fig = F._figura(3.6, 2.6)
    ax = fig.add_subplot(111)
    F._estilo(ax, rejilla="both")
    for est in ORDEN:
        fs = sorted((f for f in filas if f["estructura"] == est),
                    key=lambda f: f["n"])
        if not fs:
            continue
        n = [f["n"] for f in fs]
        ax.plot(n, [f["E_app_spinpy_MPa"] for f in fs], "-", lw=1.6,
                color=COLOR[est], alpha=0.6,
                label=f"{nombre(est, t)} · {t['app']}")
        ax.plot(n, [f["E_app_febio_MPa"] for f in fs], "o", ms=3.5,
                mfc="white", color=COLOR[est],
                label=f"{nombre(est, t)} · FEBio")
        ses = [(f["n"], f["E_app_sesion_MPa"]) for f in fs
               if "E_app_sesion_MPa" in f]
        if ses:
            ax.plot(*zip(*ses), "x", ms=5, color=F.TINTA,
                    label="sesión de la app" if idioma == "es"
                    else "app session")
    ax.set_xlabel(t["n"], fontsize=7, color=F.TINTA_2)
    ax.set_ylabel(t["E_app"], fontsize=7, color=F.TINTA_2)
    ax.set_title(t["convergencia"], fontsize=7.5, color=F.TINTA)
    ax.set_xticks(list(N_CONV))
    fig.tight_layout()
    h, l = ax.get_legend_handles_labels()
    # la marca de la sesion al final: asi cada columna es una estructura
    ses = [i for i, x in enumerate(l) if x.startswith(("sesión", "app session"))]
    orden = [i for i in range(len(l)) if i not in ses] + ses[:1]
    _leyenda_abajo(fig, [h[i] for i in orden], [l[i] for i in orden], 4)
    return _guardar(fig, DIR_FIGS / idioma / "fig_convergencia")


def fig_fallo(idioma):
    t = T[idioma]
    filas = leer("fallo")
    if not filas:
        return None
    fig = F._figura(6.6, 2.5)
    axs = [fig.add_subplot(1, 2, k + 1) for k in range(2)]
    for ax in axs:
        F._estilo(ax, rejilla="both")
    for est in ORDEN:
        fs = sorted((f for f in filas if f["estructura"] == est),
                    key=lambda f: f["paso"])
        if not fs:
            continue
        p = [f["paso"] for f in fs]
        E0s, E0f = fs[0]["E_app_spinpy_MPa"], fs[0]["E_app_febio_MPa"]
        axs[0].plot(p, [f["E_app_spinpy_MPa"] / E0s for f in fs], "-",
                    lw=1.6, color=COLOR[est], alpha=0.6,
                    label=f"{nombre(est, t)} · {t['app']}")
        axs[0].plot(p, [f["E_app_febio_MPa"] / E0f for f in fs], "o", ms=3.5,
                    mfc="white", color=COLOR[est],
                    label=f"{nombre(est, t)} · FEBio")
        axs[1].plot(p, [f["F_fallo_spinpy_N"] for f in fs], "-", lw=1.6,
                    color=COLOR[est], alpha=0.6)
        axs[1].plot(p, [f["F_fallo_febio_N"] for f in fs], "o", ms=3.5,
                    mfc="white", color=COLOR[est])
    axs[0].set_ylabel(t["E_rel"], fontsize=7, color=F.TINTA_2)
    axs[1].set_ylabel(t["F_fallo"], fontsize=7, color=F.TINTA_2)
    for ax in axs:
        ax.set_xlabel(t["paso"], fontsize=7, color=F.TINTA_2)
    fig.suptitle(t["fallo"], fontsize=7.5, color=F.TINTA)
    fig.tight_layout()
    h, l = axs[0].get_legend_handles_labels()
    _leyenda_abajo(fig, h, l, 3)
    return _guardar(fig, DIR_FIGS / idioma / "fig_fallo")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--idioma", nargs="*", default=["es", "en"])
    ap.add_argument("--sin-3d", action="store_true")
    a = ap.parse_args()
    for idi in a.idioma:
        hechas = [fig_coincidencia(idi), fig_rigidez(idi), fig_colas(idi),
                  fig_tensor(idi), fig_convergencia(idi), fig_fallo(idi),
                  fig_paridad("comparado", idi),
                  fig_paridad("compresion", idi)]
        if not a.sin_3d:
            hechas += [panel_mapas("comparado", "vm", idi),
                       panel_mapas("compresion", "vm", idi),
                       panel_mapas("comparado", "desp", idi)]
        for h in hechas:
            if h is not None:
                print(h)


if __name__ == "__main__":
    main()
