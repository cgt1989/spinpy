"""
figuras.py — Tablas (Markdown) y figuras de la comparativa de motores.

    python comparativa_motores/figuras.py

Lee `resultados/*.jsonl` (los escribe `correr.py`) y deja las figuras en
`figs/` y las tablas en `tablas.md`. Nada se calcula aqui: solo se ordena lo
medido.

Colores: un color FIJO por motor (sigue a la entidad, no al rango); la app, la
linea base, en tinta neutra. La paleta categorica es la de referencia del
proyecto, validada para daltonismo (CVD dE >= 9 entre adyacentes); como tres
tonos quedan por debajo de 3:1 de contraste sobre el fondo, cada motor lleva
ademas su propio marcador y las cifras estan en las tablas.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
RES = AQUI / "resultados"
FIGS = AQUI / "figs"

MOTORES = ("app", "ngsolve", "fenicsx", "skfem", "sfepy")
NOMBRE = {"app": "App (spinpy)", "ngsolve": "NGSolve", "fenicsx": "FEniCSx",
          "skfem": "scikit-fem", "sfepy": "SfePy"}
COLOR = {"app": "#52514e", "ngsolve": "#2a78d6", "fenicsx": "#eb6834",
         "skfem": "#1baf7a", "sfepy": "#eda100"}
MARCA = {"app": "o", "ngsolve": "s", "fenicsx": "^", "skfem": "D",
         "sfepy": "v"}
TINTA, TINTA_2, REJILLA, FONDO = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def leer(grupo):
    f = RES / f"{grupo}.jsonl"
    if not f.exists():
        return []
    return [json.loads(x) for x in f.read_text().splitlines() if x.strip()]


def t_total(d):
    return sum((d.get("tiempos") or {}).values()) if d.get("ok") else None


def t_sin_jit(d):
    t = dict(d.get("tiempos") or {})
    t.pop("jit", None)
    return sum(t.values()) if d.get("ok") else None


def fmt(x, f=".3g"):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return format(x, f)


def sci(x):
    if x is None:
        return "—"
    if x == 0:
        return "0"
    return f"{x:.1e}"


# ---------------------------------------------------------------------------
# Tablas
# ---------------------------------------------------------------------------

def tabla_lineal(filas, titulo):
    """Una fila por (caso, motor, resolvedor)."""
    L = [f"#### {titulo}", "",
         "| Caso | GDL | Motor | Resolvedor | Tiempo (s) | JIT (s) | "
         "Memoria (MB) | Residuo | Iter. | Δu/ref | ΔvM/ref | "
         "E_app (MPa) | p99 sup. (MPa) | σ fallo (MPa) | dF |",
         "|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
         "---:|---:|"]
    for d in filas:
        m = d["meta_caso"]
        if not d.get("ok"):
            L.append(f"| {d['caso']} | {m['n_gdl']} | {NOMBRE[d['motor']]} | "
                     f"{d['solver']} | **falló**: {d.get('motivo', '')[:60]} "
                     "| | | | | | | | | | |")
            continue
        E = d.get("E_app_area_MPa") if m["tipo"] == "tet10" \
            else d.get("E_app_nodal_MPa")
        L.append(
            f"| {d['caso']} | {d.get('n_gdl', m['n_gdl'])} | "
            f"{NOMBRE[d['motor']]} | {d.get('solver') or '—'} | "
            f"{fmt(t_sin_jit(d), '.2f')} | "
            f"{fmt((d.get('tiempos') or {}).get('jit'), '.2f')} | "
            f"{fmt(d.get('rss_pico_MB'), '.0f')} | "
            f"{sci(d.get('residuo_rel'))} | {d.get('iteraciones') or '—'} | "
            f"{sci(d.get('du_rel_ref'))} | {sci(d.get('dvm_rel_ref'))} | "
            f"{fmt(E, '.6g')} | {fmt(d.get('vm_p99_sup_MPa'), '.5g')} | "
            f"{fmt(d.get('sigma_fallo_MPa'), '.5g')} | {sci(d.get('dF_rel'))} |")
    return L + [""]


def tabla_nl(filas):
    L = ["#### No lineal (control por desplazamiento del plato)", "",
         "| Caso | GDL | Motor | Tiempo (s) | Memoria (MB) | "
         "Iteraciones de Newton por paso | F por paso (N) | "
         "Error frente a la solución cerrada | Δu/ref |",
         "|---|---:|---|---:|---:|---|---|---:|---:|"]
    for d in filas:
        m = d["meta_caso"]
        if not d.get("ok"):
            L.append(f"| {d['caso']} | {m['n_gdl']} | {NOMBRE[d['motor']]} | "
                     f"**falló**: {d.get('motivo', '')[:70]} | | | | | |")
            continue
        F = d.get("F_reac_N") or []
        err = "—"
        if d.get("F_exacta_N"):
            e = np.abs(np.array(F) / np.array(d["F_exacta_N"]) - 1).max()
            err = sci(float(e))
        L.append(f"| {d['caso']} | {m['n_gdl']} | {NOMBRE[d['motor']]} | "
                 f"{fmt(t_sin_jit(d), '.2f')} | "
                 f"{fmt(d.get('rss_pico_MB'), '.0f')} | "
                 f"{d.get('iteraciones_newton')} | "
                 f"{', '.join(fmt(x, '.6g') for x in F)} | {err} | "
                 f"{sci(d.get('du_rel_ref'))} |")
    return L + [""]


# ---------------------------------------------------------------------------
# Figuras
# ---------------------------------------------------------------------------

def _estilo(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(TINTA_2)
        ax.spines[s].set_linewidth(0.6)
    ax.tick_params(colors=TINTA_2, labelsize=7, width=0.6)
    ax.grid(True, which="major", color=REJILLA, lw=0.5)
    ax.set_axisbelow(True)


def fig_escalado(filas, destino, clave, etiqueta_y, titulo):
    """Tiempo o memoria frente a GDL, log-log, por motor; (a) iterativo,
    (b) directo. La app en los dos paneles (su resolvedor es el suyo)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for ax, modo, et in zip(axs, ("iterativo", "directo"),
                            ("(a) resolvedor iterativo",
                             "(b) resolvedor directo")):
        _estilo(ax)
        for mot in MOTORES:
            pts = [(d["meta_caso"]["n_gdl"], clave(d)) for d in filas
                   if d["motor"] == mot and d.get("ok") and clave(d)
                   and (d["solver"] == modo or mot == "app")]
            if not pts:
                continue
            pts.sort()
            x, y = zip(*pts)
            ax.plot(x, y, "-", marker=MARCA[mot], color=COLOR[mot], lw=1.6,
                    ms=5, mfc=FONDO, mew=1.3, label=NOMBRE[mot])
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("grados de libertad", fontsize=7.5, color=TINTA_2)
        ax.set_title(et, fontsize=8, color=TINTA, loc="left")
    axs[0].set_ylabel(etiqueta_y, fontsize=7.5, color=TINTA_2)
    axs[0].legend(frameon=False, fontsize=6.5, labelcolor=TINTA)
    fig.suptitle(titulo, fontsize=8.5, color=TINTA, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(destino, dpi=200, facecolor=FONDO)
    plt.close(fig)
    return destino


def fig_exactitud(filas, destino):
    """Diferencia de desplazamientos frente a la referencia directa."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    _estilo(ax)
    for mot in MOTORES:
        pts = [(d["meta_caso"]["n_gdl"], max(d["du_rel_ref"], 1e-16),
                d["solver"]) for d in filas
               if d["motor"] == mot and d.get("ok")
               and d.get("du_rel_ref") is not None
               and list(d.get("referencia") or []) != [mot, d["solver"]]]
        for modo, relleno in (("iterativo", FONDO), ("directo", None),
                              ("propio", FONDO)):
            q = sorted((x, y) for x, y, s in pts if s == modo)
            if not q:
                continue
            x, y = zip(*q)
            ax.plot(x, y, ls="-" if modo != "directo" else ":",
                    marker=MARCA[mot], color=COLOR[mot], lw=1.2, ms=5,
                    mfc=relleno or COLOR[mot], mew=1.2,
                    label=f"{NOMBRE[mot]} ({modo})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(1e-8, color=TINTA_2, lw=0.8, ls="--")
    ax.text(ax.get_xlim()[0], 1.3e-8, " tolerancia de la app (1e-8)",
            fontsize=6, color=TINTA_2, va="bottom")
    ax.set_xlabel("grados de libertad", fontsize=7.5, color=TINTA_2)
    ax.set_ylabel("max|u − u_ref| / max|u_ref|", fontsize=7.5, color=TINTA_2)
    ax.legend(frameon=False, fontsize=5.5, labelcolor=TINTA, ncol=2)
    fig.tight_layout()
    fig.savefig(destino, dpi=200, facecolor=FONDO)
    plt.close(fig)
    return destino


def fig_no_lineal(filas, destino):
    """(a) Curva fuerza-deformacion del espinodoide por motor (se superponen)
    y la lineal de referencia; (b) iteraciones de Newton por paso."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    esp = [d for d in filas if d.get("ok") and not d["meta_caso"].get("bloque")]
    if not esp:
        return None
    mats = sorted({d["meta_caso"]["material"] for d in esp})
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0))
    _estilo(axs[0])
    _estilo(axs[1])
    anchos = 0.8 / max(1, len({d["motor"] for d in esp}))
    for d in esp:
        eps = np.array(d["meta_caso"]["eps_plato"]) * 100
        F = d["F_reac_N"]
        ls = "-" if d["meta_caso"]["material"] == "svk" else "--"
        axs[0].plot(eps, F, ls, marker=MARCA[d["motor"]],
                    color=COLOR[d["motor"]], lw=1.4, ms=5, mfc=FONDO,
                    label=f"{NOMBRE[d['motor']]}, "
                          f"{'SVK' if ls == '-' else 'neo-Hookeano'}")
    axs[0].set_xlabel("deformación del plato (%)", fontsize=7.5,
                      color=TINTA_2)
    axs[0].set_ylabel("fuerza de reacción (N)", fontsize=7.5, color=TINTA_2)
    axs[0].set_title("(a) respuesta del espinodoide", fontsize=8, color=TINTA,
                     loc="left")
    axs[0].legend(frameon=False, fontsize=5.5, labelcolor=TINTA)
    svk = [d for d in esp if d["meta_caso"]["material"] == mats[-1]]
    mots = [m for m in MOTORES if any(d["motor"] == m for d in svk)]
    for j, m in enumerate(mots):
        d = next(d for d in svk if d["motor"] == m)
        it = d.get("iteraciones_newton") or []
        x = np.arange(len(it))
        axs[1].bar(x + (j - (len(mots) - 1) / 2) * anchos, it, anchos * 0.9,
                   color=COLOR[m], label=NOMBRE[m])
    axs[1].set_xlabel("paso de carga", fontsize=7.5, color=TINTA_2)
    axs[1].set_ylabel("iteraciones de Newton", fontsize=7.5, color=TINTA_2)
    axs[1].set_title(f"(b) iteraciones por paso ({mats[-1]})", fontsize=8,
                     color=TINTA, loc="left")
    axs[1].legend(frameon=False, fontsize=6, labelcolor=TINTA)
    fig.tight_layout()
    fig.savefig(destino, dpi=200, facecolor=FONDO)
    plt.close(fig)
    return destino


def main():
    FIGS.mkdir(exist_ok=True)
    ex, nl = leer("exactos"), leer("nl")
    cav, eh, et = leer("cavidad"), leer("espinodoide_hex"), \
        leer("espinodoide_tet")
    L = ["# Tablas de la comparativa de motores (generadas por figuras.py)",
         ""]
    if ex:
        L += tabla_lineal(ex, "Casos con solución exacta (bloque macizo)")
    if cav:
        L += tabla_lineal(cav, "Cavidad esférica")
    if eh:
        L += tabla_lineal(eh, "Espinodoide, malla de ladrillos (hex8)")
    if et:
        L += tabla_lineal(et, "Espinodoide, malla suave (TET10)")
    if nl:
        L += tabla_nl(nl)
    (AQUI / "tablas.md").write_text("\n".join(L), encoding="utf-8")
    if eh:
        fig_escalado(eh, FIGS / "fig1_tiempo_hex8.png", t_sin_jit,
                     "tiempo de pared (s)",
                     "Espinodoide hex8: tiempo de montaje y resolución")
        fig_escalado(eh, FIGS / "fig2_memoria_hex8.png",
                     lambda d: d.get("rss_pico_MB"), "memoria de pico (MB)",
                     "Espinodoide hex8: memoria de pico del proceso")
    if et:
        fig_escalado(et, FIGS / "fig3_tiempo_tet10.png", t_sin_jit,
                     "tiempo de pared (s)",
                     "Espinodoide TET10: tiempo de montaje y resolución")
        fig_escalado(et, FIGS / "fig4_memoria_tet10.png",
                     lambda d: d.get("rss_pico_MB"), "memoria de pico (MB)",
                     "Espinodoide TET10: memoria de pico del proceso")
    todas = [d for d in eh + et + cav if d.get("ok")]
    if todas:
        fig_exactitud(todas, FIGS / "fig5_exactitud.png")
    if nl:
        fig_no_lineal(nl, FIGS / "fig6_no_lineal.png")
    print((AQUI / "tablas.md").read_text()[:3000])


if __name__ == "__main__":
    main()
