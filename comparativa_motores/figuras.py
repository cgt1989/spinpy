"""
figuras.py: tablas (Markdown) y figuras de la comparativa de motores.

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


def modo(d):
    """'directo' o 'iterativo' a partir del nombre del resolvedor.

    En el JSON la clave `solver` guarda el NOMBRE que devuelve el motor (el
    modo pedido lo sobrescribe la meta del motor al guardar la fila)."""
    n = str(d.get("solver") or "")
    if d["motor"] == "app":
        return "propio"
    return "directo" if any(k in n for k in ("SuperLU", "PARDISO", "MUMPS",
                                               "sparsecholesky", "Newton")) \
        else "iterativo"


def t_total(d):
    return sum((d.get("tiempos") or {}).values()) if d.get("ok") else None


def t_sin_jit(d):
    t = dict(d.get("tiempos") or {})
    t.pop("jit", None)
    return sum(t.values()) if d.get("ok") else None


def fmt(x, f=".3g"):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/d"
    return format(x, f)


def sci(x):
    if x is None:
        return "n/d"
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
            f"{NOMBRE[d['motor']]} | {d.get('solver') or 'n/d'} | "
            f"{fmt(t_sin_jit(d), '.2f')} | "
            f"{fmt((d.get('tiempos') or {}).get('jit'), '.2f')} | "
            f"{fmt(d.get('rss_pico_MB'), '.0f')} | "
            f"{sci(d.get('residuo_rel'))} | {d.get('iteraciones') or 'n/d'} | "
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
        err = "n/d"
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

def _ejes_log(ax):
    from matplotlib.ticker import LogLocator, NullFormatter
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_major_locator(LogLocator(base=10))


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
    for ax, modo_p, et in zip(axs, ("iterativo", "directo"),
                            ("(a) resolvedor iterativo",
                             "(b) resolvedor directo")):
        _estilo(ax)
        for mot in MOTORES:
            for var in ("", "sparsecholesky"):
                pts = [(d["meta_caso"]["n_gdl"], clave(d)) for d in filas
                       if d["motor"] == mot and d.get("ok") and clave(d)
                       and (modo(d) == modo_p or mot == "app")
                       and (("sparsecholesky" in str(d.get("solver")))
                            == bool(var))]
                if not pts:
                    continue
                pts = sorted(dict(pts).items())
                x, y = zip(*pts)
                ax.plot(x, y, "--" if var else "-", marker=MARCA[mot],
                        color=COLOR[mot], lw=1.6, ms=5,
                        mfc=COLOR[mot] if var else FONDO, mew=1.3,
                        label=NOMBRE[mot] + (" (Cholesky propio)" if var
                                             else " (resolvedor propio)"
                                             if mot == "app" else ""))
        ax.set_xscale("log")
        ax.set_yscale("log")
        _ejes_log(ax)
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
                modo(d)) for d in filas
               if d["motor"] == mot and d.get("ok")
               and (d.get("du_rel_ref") or 0) > 0]
        for m_, relleno in (("iterativo", FONDO), ("directo", None),
                            ("propio", FONDO)):
            q = sorted(dict((x, y) for x, y, s in pts if s == m_).items())
            if not q:
                continue
            x, y = zip(*q)
            ax.plot(x, y, ls="-" if m_ != "directo" else ":",
                    marker=MARCA[mot], color=COLOR[mot], lw=1.2, ms=5,
                    mfc=relleno or COLOR[mot], mew=1.2,
                    label=f"{NOMBRE[mot]} ({m_})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    _ejes_log(ax)
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


# ---------------------------------------------------------------------------
# Convergencia con solucion exacta (convergencia.py)
# ---------------------------------------------------------------------------

ELEMENTO = {"hex8": "hex8 (Q1)", "tet4": "TET4 (P1)", "tet10": "TET10 (P2)"}
#: Las rectas de cada elemento van en tinta (el color sigue al MOTOR, como en
#: el resto de figuras); el elemento lo distinguen el trazo y la etiqueta.
TRAZO = {"hex8": (TINTA_2, "-"), "tet4": (TINTA_2, "--"),
         "tet10": (TINTA, ":")}


def pendiente(x, y):
    """Pendiente de la recta de minimos cuadrados de log(y) frente a log(x)."""
    return float(np.polyfit(np.log(x), np.log(y), 1)[0])


def serie_h(filas, tipo, clave, motor=None):
    """(h, e) de un elemento: el motor pedido o, si no, la referencia de
    cada malla (los motores coinciden a 1e-9, ver `du_rel_ref`)."""
    pts = {}
    for d in filas:
        if d["estudio"] != "h" or not d.get("ok") or d["elemento"] != tipo:
            continue
        if motor is not None and d["motor"] != motor:
            continue
        if motor is None and (d.get("du_rel_ref") or 0) != 0:
            continue
        pts[d["h"]] = d[clave]
    h = np.array(sorted(pts, reverse=True))
    return h, np.array([pts[x] for x in h])


def tabla_convergencia(filas):
    """Orden observado por elemento y motor (pendiente global y local)."""
    L = ["#### Convergencia con solución exacta (refinamiento h)", "",
         "| Elemento | Motor | n | GDL | error L2 rel. | error H1 rel. | "
         "orden local L2 | orden local H1 | Δu/ref | tiempo (s) |",
         "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for tipo in ELEMENTO:
        for mot in MOTORES:
            fs = sorted((d for d in filas if d["estudio"] == "h"
                         and d.get("ok") and d["elemento"] == tipo
                         and d["motor"] == mot), key=lambda d: d["n"])
            ant = None
            for d in fs:
                oL = oH = "n/d"
                if ant is not None:
                    r = np.log(ant["h"] / d["h"])
                    oL = f"{np.log(ant['e_L2'] / d['e_L2']) / r:.3f}"
                    oH = f"{np.log(ant['e_H1'] / d['e_H1']) / r:.3f}"
                L.append(f"| {ELEMENTO[tipo]} | {NOMBRE[mot]} | {d['n']} | "
                         f"{d['n_gdl']} | {d['e_L2_rel']:.4e} | "
                         f"{d['e_H1_rel']:.4e} | {oL} | {oH} | "
                         f"{sci(d.get('du_rel_ref'))} | "
                         f"{fmt(d.get('tiempo_s'), '.2f')} |")
                ant = d
    L += ["", "#### Convergencia con solución exacta (refinamiento p, "
          "malla fija)", "",
          "| Motor | n | grado p | GDL | error L2 rel. | error H1 rel. | "
          "tiempo (s) |", "|---|---:|---:|---:|---:|---:|---:|"]
    for d in filas:
        if d["estudio"] == "p" and d.get("ok"):
            L.append(f"| {NOMBRE[d['motor']]} | {d['n']} | {d['grado']} | "
                     f"{d['n_gdl']} | {d['e_L2_rel']:.4e} | "
                     f"{d['e_H1_rel']:.4e} | {fmt(d.get('tiempo_s'), '.2f')}"
                     " |")
    return L + [""]


def fig_convergencia(filas, destino):
    """(a, b) Error relativo en L2 y en la seminorma H1 frente al tamano de
    elemento h, con la pendiente de cada elemento; los marcadores de los
    motores se superponen. (c) Refinamiento p con la malla fija: log(e)
    frente al grado es una recta. (d) Error L2 frente a los GDL: h frente a
    p."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.2))
    ax_a, ax_b, ax_c, ax_d = axs.ravel()
    mots_h = [m for m in MOTORES if any(d["estudio"] == "h" and d["motor"] == m
                                        and d.get("ok") for d in filas)]
    for ax, clave, nombre, letra in ((ax_a, "e_L2_rel", "L2", "(a)"),
                                     (ax_b, "e_H1_rel", "H1", "(b)")):
        _estilo(ax)
        ax.grid(True, which="minor", color=REJILLA, lw=0.3)
        for tipo, (col, ls) in TRAZO.items():
            h, e = serie_h(filas, tipo, clave)
            if h.size < 2:
                continue
            m_ = pendiente(h, e)
            ax.plot(h, e, ls, color=col, lw=1.4, zorder=1,
                    label=f"{ELEMENTO[tipo]}: pendiente {m_:.3f}"
                    .replace(".", ","))
            for j, mot in enumerate(mots_h):
                hm, em = serie_h(filas, tipo, clave, mot)
                if hm.size:
                    ax.plot(hm, em, ls="none", marker=MARCA[mot],
                            color=COLOR[mot], ms=7 - 0.9 * j, mfc="none",
                            mew=1.1, zorder=2 + j)
        ax.set_xscale("log")
        ax.set_yscale("log")
        from matplotlib.ticker import FixedLocator, NullFormatter
        hs = sorted({d["h"] for d in filas if d["estudio"] == "h"})
        ax.xaxis.set_major_locator(FixedLocator(hs))
        ax.xaxis.set_minor_locator(FixedLocator([]))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xticklabels([f"1/{round(1 / x)}" for x in hs])
        ax.set_xlabel("tamaño de elemento h (mm)", fontsize=7.5,
                      color=TINTA_2)
        ax.set_ylabel(f"error relativo en {nombre}", fontsize=7.5,
                      color=TINTA_2)
        ax.set_title(f"{letra} refinamiento h, norma {nombre}", fontsize=8,
                     color=TINTA, loc="left")
        ax.legend(frameon=False, fontsize=6.3, labelcolor=TINTA,
                  loc="lower right", handlelength=2.6)
    # motores (comun a a y b), arriba a la izquierda de (a): ahi no hay datos
    ax_a.add_artist(ax_a.legend_)
    ax_a.legend(handles=[Line2D([], [], ls="none", marker=MARCA[m],
                                color=COLOR[m], mfc="none", mew=1.1, ms=6,
                                label=NOMBRE[m]) for m in mots_h],
                frameon=False, fontsize=6.3, labelcolor=TINTA,
                loc="upper left", title="motor", title_fontsize=6.3)
    # (c) refinamiento p
    _estilo(ax_c)
    pp = [d for d in filas if d["estudio"] == "p" and d.get("ok")]
    for mot in [m for m in MOTORES if any(d["motor"] == m for d in pp)]:
        fs = sorted((d for d in pp if d["motor"] == mot),
                    key=lambda d: d["grado"])
        g = np.array([d["grado"] for d in fs])
        for clave, ls, et in (("e_L2_rel", "-", "L2"), ("e_H1_rel", "--",
                                                        "H1")):
            e = np.array([d[clave] for d in fs])
            b = -np.polyfit(g, np.log(e), 1)[0]
            ax_c.plot(g, e, ls, marker=MARCA[mot], color=COLOR[mot], lw=1.4,
                      ms=5, mfc=FONDO if ls == "--" else COLOR[mot],
                      mew=1.2, label=f"{NOMBRE[mot]}, {et}: "
                      f"e ∝ exp(−{b:.2f} p)".replace(".", ","))
    ax_c.set_yscale("log")
    ax_c.set_xlabel("grado del polinomio p (malla fija, h = 0,25 mm)",
                    fontsize=7.5, color=TINTA_2)
    ax_c.set_ylabel("error relativo", fontsize=7.5, color=TINTA_2)
    ax_c.set_title("(c) refinamiento p (TET4 → TET10 → TETN)", fontsize=8,
                   color=TINTA, loc="left")
    ax_c.legend(frameon=False, fontsize=6, labelcolor=TINTA)
    # (d) h frente a p por GDL
    _estilo(ax_d)
    ax_d.grid(True, which="minor", color=REJILLA, lw=0.3)
    for tipo in ("tet4", "tet10"):
        col, ls = TRAZO[tipo]
        fs = sorted((d for d in filas if d["estudio"] == "h" and d.get("ok")
                     and d["elemento"] == tipo and d["motor"] == "ngsolve"),
                    key=lambda d: d["n_gdl"])
        if fs:
            ax_d.plot([d["n_gdl"] for d in fs], [d["e_L2_rel"] for d in fs],
                      ls, color=col, lw=1.4, marker="o", ms=4, mfc=FONDO,
                      label=f"refinamiento h, {ELEMENTO[tipo]}")
    for mot in [m for m in MOTORES if any(d["motor"] == m for d in pp)]:
        fs = sorted((d for d in pp if d["motor"] == mot),
                    key=lambda d: d["grado"])
        ax_d.plot([d["n_gdl"] for d in fs], [d["e_L2_rel"] for d in fs], "-",
                  marker=MARCA[mot], color=COLOR[mot], lw=1.4, ms=5,
                  label=f"refinamiento p, {NOMBRE[mot]}")
    ax_d.set_xscale("log")
    ax_d.set_yscale("log")
    ax_d.set_xlabel("grados de libertad", fontsize=7.5, color=TINTA_2)
    ax_d.set_ylabel("error relativo en L2", fontsize=7.5, color=TINTA_2)
    ax_d.set_title("(d) coste: error frente a GDL", fontsize=8, color=TINTA,
                   loc="left")
    ax_d.legend(frameon=False, fontsize=6, labelcolor=TINTA)
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
    if leer("convergencia"):
        L += tabla_convergencia(leer("convergencia"))
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
    conv = leer("convergencia")
    if conv:
        fig_convergencia(conv, FIGS / "fig5_convergencia.png")
    # Concordancia entre motores sobre una misma malla: ya no es figura del
    # informe (la Figura 5 es la convergencia), queda como complemento.
    todas = [d for d in eh + et + cav if d.get("ok")]
    if todas:
        fig_exactitud(todas, FIGS / "concordancia_motores.png")
    if nl:
        fig_no_lineal(nl, FIGS / "fig6_no_lineal.png")
    print((AQUI / "tablas.md").read_text()[:3000])


if __name__ == "__main__":
    main()
