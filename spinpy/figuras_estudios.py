"""
figuras_estudios.py: figuras 12 a 19 del informe para publicacion.

Todas salen SOLO del documento de la sesion, sin volumenes, asi que se
rehacen igual desde un JSON exportado. Cada una devuelve [png, pdf] o [] si
el documento no tiene su dato.

  fig12_perdida      rigidez y morfometria durante la perdida osea simulada
  fig13_fallo        carga de fallo, rigidez y dano en el fallo progresivo
  fig14_superficie   superficie 3D del modulo de Young direccional E(n)
  fig15_curvaturas   perfil de curvaturas principales (k1, k2) ponderado por
                     area y fracciones de silla, convexa y concava
  fig16_ef           distribucion del Ellipsoid Factor (placas frente a barras)
  fig17_febio        paridad spinpy frente a FEBio
  fig18_superficie   sensibilidad de la morfometria y de E_app a desplazar la
                     superficie (sustituto de un cambio de umbral)
  fig19_contorno     E/E_s con cada condicion de contorno

Estilo, colores y guardado son los de `figuras.py`: el mismo color por
estructura en todo el informe, y el color nunca va solo.
"""

from __future__ import annotations

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .figuras import (COLOR, EJE, FONDO, NOMBRES, ORDEN, REJILLA, TINTA,
                      TINTA_2, _estilo, _f, _figura, _guardar, _i, _n,
                      modulo_direccional)


def _titulo(ax, txt):
    ax.set_title(txt, fontsize=8, color=TINTA, loc="left")


def _etiquetas(ax, x, y):
    ax.set_xlabel(x, fontsize=7.5, color=TINTA_2)
    ax.set_ylabel(y, fontsize=7.5, color=TINTA_2)


def _leyenda(ax, **kw):
    kw.setdefault("loc", "best")
    ax.legend(frameon=False, fontsize=6.3, labelcolor=TINTA, **kw)


def _por_estructura(rec):
    """{codigo: valor} del registro en el orden VOI, spinodoide, dual."""
    pe = (rec or {}).get("por_estructura") or {}
    return {e: pe[e] for e in ORDEN if isinstance(pe.get(e), dict)}


# ---------------------------------------------------------------------------
# Figura 12: perdida osea simulada
# ---------------------------------------------------------------------------

def fig_perdida(doc, destino, idioma="es"):
    from .informe import estructura_de
    i = _i(idioma)
    rec = (doc.get("resultados") or {}).get("simulacion_perdida")
    if not isinstance(rec, dict) or len(rec.get("pasos") or []) < 2:
        return []
    est = estructura_de(rec.get("estructura_codigo") or rec.get("estructura"),
                        doc)
    filas = rec["pasos"]
    col = COLOR[est]
    x = np.array([100.0 * (1.0 - float(f.get("hueso_rel", np.nan)))
                  for f in filas])
    rec_fase = np.array([f.get("fase") == "recuperacion" for f in filas])

    fig = _figura(6.6, 2.7)
    ax = fig.add_subplot(1, 2, 1)
    _estilo(ax, "y")
    for clave, marca, etq in (
            ("E_rel", "o", (r"$E/E_0$", r"$E/E_0$")),
            ("sigma_fallo_rel", "s",
             (r"$\sigma_{fallo}/\sigma_0$", r"$\sigma_{fail}/\sigma_0$"))):
        y = np.array([_f(f, clave) if _f(f, clave) is not None else np.nan
                      for f in filas])
        ls = "-" if clave == "E_rel" else "--"
        ax.plot(x[~rec_fase], y[~rec_fase], ls, marker=marca, color=col,
                ms=3.5, lw=1.2, mfc=col if clave == "E_rel" else FONDO,
                label=etq[i])
        if rec_fase.any():
            ax.plot(x[rec_fase], y[rec_fase], ls, marker=marca, color=TINTA_2,
                    ms=3.2, lw=1.0, mfc=FONDO)
    if rec_fase.any():
        ax.plot([], [], "-", color=TINTA_2, lw=1.0,
                label=("recuperación", "recovery")[i])
    ax.set_ylim(bottom=0)
    _etiquetas(ax, ("hueso retirado [%]", "bone removed [%]")[i],
               ("cociente frente al paso 0", "ratio to step 0")[i])
    _titulo(ax, f"(a) {NOMBRES[est][i]}, "
            + str(rec.get("protocolo", "")).replace("_", " "))
    _leyenda(ax)

    ax = fig.add_subplot(1, 2, 2)
    _estilo(ax, "y")
    f0 = filas[0]
    for clave, etq, marca in (("BVTV", "BV/TV", "o"), ("TbTh", "Tb.Th", "s"),
                              ("TbN", "Tb.N", "^"), ("ConnD", "Conn.D", "D")):
        v0 = _f(f0, clave)
        if not v0:
            continue
        y = [(_f(f, clave) or np.nan) / v0 for f in filas]
        ax.plot([f["paso"] for f in filas], y, "-", marker=marca, ms=3.2,
                lw=1.0, color=TINTA if clave == "BVTV" else TINTA_2,
                mfc=FONDO if clave != "BVTV" else TINTA, label=etq)
    ax.axhline(1.0, color=EJE, lw=0.8)
    # Conn.D puede multiplicarse por decenas al perforarse las placas; en
    # escala lineal aplastaria las demas curvas.
    ax.set_yscale("log")
    _etiquetas(ax, ("paso", "step")[i],
               ("cociente frente al paso 0 (log)", "ratio to step 0 (log)")[i])
    _titulo(ax, ("(b) morfometría en cada paso", "(b) morphometry at each "
                                                   "step")[i])
    _leyenda(ax, ncol=2)
    fig.tight_layout()
    return _guardar(fig, destino)


# ---------------------------------------------------------------------------
# Figura 13: fallo progresivo
# ---------------------------------------------------------------------------

def fig_fallo(doc, destino, idioma="es"):
    from .informe import estructura_de
    i = _i(idioma)
    rec = (doc.get("resultados") or {}).get("fallo_progresivo")
    if not isinstance(rec, dict) or len(rec.get("pasos") or []) < 2:
        return []
    est = estructura_de(rec.get("estructura_codigo") or rec.get("estructura"),
                        doc)
    filas = rec["pasos"]
    col = COLOR[est]
    paso = np.array([int(f["paso"]) for f in filas])
    F = np.array([_f(f, "F_fallo") if _f(f, "F_fallo") is not None
                  else np.nan for f in filas])
    colapso = (rec.get("resumen") or {}).get("colapso_en_paso")

    fig = _figura(6.6, 2.7)
    ax = fig.add_subplot(1, 2, 1)
    _estilo(ax, "y")
    ax.plot(paso, F, "-o", color=col, ms=3.5, lw=1.2,
            label=("carga de fallo de Pistoia", "Pistoia failure load")[i])
    r = rec.get("resumen") or {}
    if r.get("paso_F_max") is not None and _f(r, "F_max") is not None:
        ax.plot([r["paso_F_max"]], [r["F_max"]], "*", color=TINTA, ms=9,
                label=("máximo antes del colapso", "maximum before "
                                                    "collapse")[i])
    if colapso is not None:
        ax.axvspan(colapso - 0.5, paso.max() + 0.5, color=REJILLA, lw=0,
                   alpha=0.7, label=("tras el colapso", "after collapse")[i])
    _etiquetas(ax, ("paso", "step")[i], "F [N]")
    _titulo(ax, f"(a) {NOMBRES[est][i]}")
    _leyenda(ax)

    ax = fig.add_subplot(1, 2, 2)
    _estilo(ax, "y")
    ax.plot(paso, [(_f(f, "E_rel") if _f(f, "E_rel") is not None else np.nan)
                   for f in filas], "-o", color=col, ms=3.5, lw=1.2,
            label=r"$E/E_0$")
    ax.plot(paso, [(_f(f, "dano_acumulado") or 0.0) for f in filas], "--s",
            color=TINTA_2, ms=3.2, lw=1.0, mfc=FONDO,
            label=("tejido roto (fracción)", "broken tissue (fraction)")[i])
    E_col = _f(r, "E_rel_colapso")
    if E_col is not None:
        ax.axhline(E_col, color=TINTA, lw=0.8, ls=":",
                   label=("umbral de colapso", "collapse threshold")[i])
    ax.set_ylim(0, 1.05)
    _etiquetas(ax, ("paso", "step")[i], ("fracción", "fraction")[i])
    _titulo(ax, ("(b) rigidez y daño acumulado", "(b) stiffness and "
                                                   "accumulated damage")[i])
    _leyenda(ax)
    fig.tight_layout()
    return _guardar(fig, destino)


# ---------------------------------------------------------------------------
# Figura 14: superficie del modulo de Young direccional
# ---------------------------------------------------------------------------

def rigideces(doc):
    """{estructura: C en MPa} de los tensores homogeneizados validos."""
    res = doc.get("resultados") or {}
    rig = {}
    for fam, suf in (("spinodoide", ""), ("dual-lattice", "_dual")):
        pe = (res.get("elastico" + suf) or {}).get("por_estructura") or {}
        for k, v in pe.items():
            est = "voi" if k == "voi" else fam
            if est in rig or (v or {}).get("C_Pa") is None:
                continue
            try:
                C = np.asarray(v["C_Pa"], float)
                if C.shape == (6, 6) and np.all(np.isfinite(C)):
                    np.linalg.inv(C)
                    rig[est] = C / 1e6
            except (ValueError, np.linalg.LinAlgError):
                pass
    return {e: rig[e] for e in ORDEN if e in rig}


def fig_superficie_E(doc, destino, idioma="es", n=60):
    """E(n) sobre la esfera de direcciones, como superficie radial en MPa.

    Misma escala de color y mismos ejes en todos los paneles: la superficie de
    una estructura mas blanda sale mas pequena, no igual de grande.
    """
    from matplotlib import cm
    from matplotlib.colors import Normalize
    i = _i(idioma)
    rig = rigideces(doc)
    if not rig:
        return []
    th, ph = np.meshgrid(np.linspace(0, np.pi, n), np.linspace(0, 2 * np.pi,
                                                                2 * n))
    d = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph),
                  np.cos(th)], axis=-1).reshape(-1, 3)
    sup = {e: modulo_direccional(C, d).reshape(th.shape)
           for e, C in rig.items()}
    vmax = max(float(v.max()) for v in sup.values())
    vmin = min(float(v.min()) for v in sup.values())
    norma = Normalize(vmin=vmin, vmax=vmax)
    from matplotlib import colormaps
    mapa = colormaps["viridis"]
    fig = _figura(2.3 * len(sup) + 0.6, 2.7)
    for k, (e, E) in enumerate(sup.items()):
        ax = fig.add_subplot(1, len(sup), k + 1, projection="3d")
        X, Y, Z = (E * np.sin(th) * np.cos(ph), E * np.sin(th) * np.sin(ph),
                   E * np.cos(th))
        ax.plot_surface(X, Y, Z, facecolors=mapa(norma(E)), rstride=1,
                        cstride=1, linewidth=0, antialiased=True,
                        shade=False)
        for s in (ax.set_xlim, ax.set_ylim, ax.set_zlim):
            s(-vmax, vmax)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("x", fontsize=6.5, color=TINTA_2, labelpad=-10)
        ax.set_ylabel("y", fontsize=6.5, color=TINTA_2, labelpad=-10)
        ax.set_zlabel("z", fontsize=6.5, color=TINTA_2, labelpad=-10)
        ax.tick_params(labelsize=0, length=0, pad=-4)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        for eje in (ax.xaxis, ax.yaxis, ax.zaxis):
            eje.pane.set_facecolor(FONDO)
            eje.pane.set_edgecolor(REJILLA)
        ax.view_init(elev=22, azim=-58)
        ax.set_title(f"({'abc'[k]}) {NOMBRES[e][i]}\nEmax/Emin = "
                     f"{_n(E.max() / E.min(), idioma)}", fontsize=7.5,
                     color=TINTA, pad=-2)
    sm = cm.ScalarMappable(norm=norma, cmap=mapa)
    fig.subplots_adjust(left=0.0, right=0.88, wspace=0.05)
    cb = fig.colorbar(sm, cax=fig.add_axes([0.915, 0.18, 0.014, 0.62]))
    cb.set_label(("E(n) [MPa]", "E(n) [MPa]")[i], fontsize=7, color=TINTA_2)
    cb.ax.tick_params(labelsize=6, colors=TINTA_2)
    return _guardar(fig, destino)


# ---------------------------------------------------------------------------
# Figura 15: curvaturas principales
# ---------------------------------------------------------------------------

def fig_curvaturas(doc, destino, idioma="es"):
    from matplotlib.colors import LogNorm
    i = _i(idioma)
    forma = (doc.get("resultados") or {}).get("forma")
    pe = {e: v["curvatura"] for e, v in _por_estructura(forma).items()
          if isinstance(v.get("curvatura"), dict)
          and v["curvatura"].get("perfil")}
    if not pe:
        return []
    fig = _figura(2.1 * len(pe) + 2.6, 2.55)
    # Columnas: un perfil por estructura, la barra de color, un hueco para
    # su etiqueta y el panel de tipos de superficie.
    ancho = [1.0] * len(pe) + [0.05, 0.32, 1.1]
    gs = fig.add_gridspec(1, len(pe) + 3, width_ratios=ancho, wspace=0.3)
    positivos = [np.asarray(v["perfil"], float) for v in pe.values()]
    pmin = min(float(p[p > 0].min()) for p in positivos if (p > 0).any())
    pmax = max(float(p.max()) for p in positivos)
    norma = LogNorm(vmin=max(pmin, pmax * 1e-4), vmax=pmax)
    im = None
    for k, (e, v) in enumerate(pe.items()):
        ax = fig.add_subplot(gs[0, k])
        P = np.asarray(v["perfil"], float)
        L = float(v["limite"])
        P = np.where(P > 0, P, np.nan)
        # P[a, b]: a recorre k1, b recorre k2. imshow pinta filas en y, asi
        # que se traspone para tener k1 en x y k2 en y.
        im = ax.imshow(P.T, origin="lower", extent=(-L, L, -L, L),
                       cmap="magma_r", norm=norma, interpolation="nearest",
                       aspect="equal")
        ax.axhline(0, color=EJE, lw=0.6)
        ax.axvline(0, color=EJE, lw=0.6)
        ax.plot([-L, L], [L, -L], ls="--", color=TINTA_2, lw=0.7)
        ax.tick_params(labelsize=6, colors=TINTA_2)
        ax.set_xlabel(r"$\kappa_1$ [mm$^{-1}$]", fontsize=7, color=TINTA_2)
        if k == 0:
            ax.set_ylabel(r"$\kappa_2$ [mm$^{-1}$]", fontsize=7,
                          color=TINTA_2)
        else:
            ax.tick_params(labelleft=False)
        r = v.get("resumen") or {}
        ax.set_title(f"({'abcd'[k]}) {NOMBRES[e][i]}\nH = "
                     f"{_n(r.get('H_medio', np.nan), idioma)} mm⁻¹",
                     fontsize=7.5, color=TINTA, loc="left")
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[0, len(pe)]))
    cb.ax.set_title(("fracción\nde área", "area\nfraction")[i], fontsize=6,
                    color=TINTA_2, loc="left", pad=4)
    cb.ax.tick_params(labelsize=5.5, colors=TINTA_2)

    ax = fig.add_subplot(gs[0, -1])
    _estilo(ax, "y")
    cats = (("silla", ("silla", "saddle"), TINTA),
            ("convexa", ("convexa", "convex"), "#9a988f"),
            ("concava", ("cóncava", "concave"), FONDO))
    abajo = np.zeros(len(pe))
    xs = np.arange(len(pe))
    for clave, etq, col in cats:
        h = np.array([100.0 * ((v.get("resumen") or {}).get(clave) or 0.0)
                      for v in pe.values()])
        ax.bar(xs, h, 0.6, bottom=abajo, color=col, edgecolor=TINTA,
               linewidth=0.4, label=etq[i])
        abajo += h
    ax.set_xticks(xs)
    ax.set_xticklabels([NOMBRES[e][i] for e in pe], fontsize=6.3,
                       rotation=20, color=TINTA_2)
    for t, e in zip(ax.get_xticklabels(), pe):
        t.set_color(COLOR[e])
    ax.set_ylim(0, 100)
    ax.set_ylabel(("área [%]", "area [%]")[i], fontsize=7, color=TINTA_2)
    _titulo(ax, f"({'abcd'[len(pe)]}) " + ("tipo de superficie",
                                           "surface type")[i])
    ax.legend(frameon=False, fontsize=6, labelcolor=TINTA,
              loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=3)
    return _guardar(fig, destino)


# ---------------------------------------------------------------------------
# Figura 16: Ellipsoid Factor
# ---------------------------------------------------------------------------

# Vafaeefar et al. exigen que los elipsoides cubran mas del 95 % del solido
# para dar el EF por valido (ver elipsoide.factor_elipsoide).
EF_RELLENO_MIN = 0.95


def fig_ellipsoid_factor(doc, destino, idioma="es"):
    i = _i(idioma)
    forma = (doc.get("resultados") or {}).get("forma")
    pe = {e: v["EF"] for e, v in _por_estructura(forma).items()
          if isinstance(v.get("EF"), dict) and v["EF"].get("histograma")}
    if not pe:
        return []
    fig = _figura(6.6, 2.6)
    ax = fig.add_subplot(1, 2, 1)
    _estilo(ax, "y")
    for e, v in pe.items():
        b = np.asarray(v["bordes"], float)
        h = np.asarray(v["histograma"], float) * 100.0
        ax.stairs(h, b, color=COLOR[e], lw=1.3,
                  label=NOMBRES[e][i] + (" *" if (v.get("EF_relleno") or 0)
                                         < EF_RELLENO_MIN else ""))
        if v.get("EF") is not None:
            ax.axvline(v["EF"], color=COLOR[e], lw=0.8, ls="--")
    for x in (-0.2, 0.2):
        ax.axvline(x, color=EJE, lw=0.7, ls=":")
    ax.text(-0.6, ax.get_ylim()[1] * 0.95, ("placa", "plate")[i],
            ha="center", va="top", fontsize=6.5, color=TINTA_2)
    ax.text(0.6, ax.get_ylim()[1] * 0.95, ("barra", "rod")[i], ha="center",
            va="top", fontsize=6.5, color=TINTA_2)
    ax.set_xlim(-1, 1)
    _etiquetas(ax, "EF", ("vóxeles [%]", "voxels [%]")[i])
    _titulo(ax, ("(a) distribución", "(a) distribution")[i])
    _leyenda(ax, loc="upper left", bbox_to_anchor=(0.0, 0.88))

    ax = fig.add_subplot(1, 2, 2)
    _estilo(ax, "y")
    xs = np.arange(len(pe))
    ancho = 0.36
    for j, (clave, etq, rel) in enumerate((
            ("EF_frac_placa", ("placas (EF < −0,2)", "plates (EF < −0.2)"),
             FONDO),
            ("EF_frac_barra", ("barras (EF > +0,2)", "rods (EF > +0.2)"),
             None))):
        for k, (e, v) in enumerate(pe.items()):
            y = 100.0 * (v.get(clave) or 0.0)
            ax.bar(k + (j - 0.5) * ancho, y, ancho,
                   color=rel or COLOR[e], edgecolor=COLOR[e], linewidth=1.0,
                   hatch="////" if rel else None)
            ax.text(k + (j - 0.5) * ancho, y + 1, _n(y, idioma, ".2g"),
                    ha="center", va="bottom", fontsize=6, color=TINTA_2)
    ax.set_xticks(xs)
    ax.set_xticklabels([NOMBRES[e][i] for e in pe], fontsize=7)
    for t, e in zip(ax.get_xticklabels(), pe):
        t.set_color(COLOR[e])
    ax.set_ylabel(("sólido [%]", "solid [%]")[i], fontsize=7.5,
                  color=TINTA_2)
    ax.legend(handles=[Patch(facecolor=FONDO, edgecolor=TINTA_2,
                             hatch="////", label=("placas", "plates")[i]),
                       Patch(facecolor=TINTA_2, edgecolor=TINTA_2,
                             label=("barras", "rods")[i])],
              frameon=False, fontsize=6.3, labelcolor=TINTA)
    _titulo(ax, ("(b) placas y barras", "(b) plates and rods")[i])
    fig.tight_layout()
    return _guardar(fig, destino)


# ---------------------------------------------------------------------------
# Figura 17: paridad spinpy frente a FEBio
# ---------------------------------------------------------------------------

MAGNITUDES_FEBIO = (("E_app_app_MPa", "E_app_MPa", "o", "E_app"),
                    ("vm_p99_sup_app_MPa", "vm_p99_sup_MPa", "s",
                     "p99 von Mises"),
                    ("sigma_fallo_app_MPa", "sigma_fallo_MPa", "^",
                     r"$\sigma_{Pistoia}$"))


def fig_paridad_febio(doc, destino, idioma="es"):
    """Cada magnitud de la app frente a la de FEBio sobre la misma malla.

    Solo registros hex8: con TET10 la malla es otra y la diferencia mezcla
    malla y solver, que no es lo que una grafica de paridad debe mostrar.
    """
    from .febio import fila_tabla
    from .informe import estructura_de
    i = _i(idioma)
    regs = ((doc.get("resultados") or {}).get("febio") or {}).get(
        "registros") or []
    puntos = []
    for reg in regs:
        if str(reg.get("malla", "hex8")).lower() != "hex8":
            continue
        f = fila_tabla(reg)
        est = estructura_de(f.get("estructura"), doc)
        for app, feb, marca, etq in MAGNITUDES_FEBIO:
            a, b = _f(f, app), _f(f, feb)
            if a and b and a > 0 and b > 0:
                puntos.append((est, a, b, marca, etq))
    if not puntos:
        return []
    fig = _figura(6.6, 2.8)
    ax = fig.add_subplot(1, 2, 1)
    _estilo(ax, "both")
    vals = np.array([[p[1], p[2]] for p in puntos])
    lo, hi = vals.min() * 0.7, vals.max() * 1.4
    ax.plot([lo, hi], [lo, hi], color=EJE, lw=0.8, zorder=1)
    for est, a, b, marca, _e in puntos:
        ax.plot(a, b, marca, color=COLOR[est], ms=4.5, mfc=FONDO, mew=1.1,
                zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    _etiquetas(ax, ("spinpy [MPa]", "spinpy [MPa]")[i], "FEBio [MPa]")
    _titulo(ax, ("(a) paridad", "(a) parity")[i])
    mans = [Line2D([], [], ls="", marker=m, color=TINTA_2, mfc=FONDO,
                   label=etq) for _a, _b, m, etq in MAGNITUDES_FEBIO]
    mans += [Line2D([], [], ls="", marker="o", color=COLOR[e], label=
                    NOMBRES[e][i]) for e in ORDEN
             if any(p[0] == e for p in puntos)]
    _leyenda(ax, handles=mans, loc="upper left")

    ax = fig.add_subplot(1, 2, 2)
    _estilo(ax, "y")
    for k, (est, a, b, marca, _e) in enumerate(puntos):
        d = abs(b / a - 1.0)
        ax.plot(k, max(d, 1e-12), marca, color=COLOR[est], ms=4.5, mfc=FONDO,
                mew=1.1)
    ax.set_yscale("log")
    ax.set_xticks([])
    _etiquetas(ax, ("magnitud y estructura", "quantity and structure")[i],
               ("|FEBio/spinpy − 1|", "|FEBio/spinpy − 1|")[i])
    _titulo(ax, ("(b) diferencia relativa", "(b) relative difference")[i])
    fig.tight_layout()
    return _guardar(fig, destino)


# ---------------------------------------------------------------------------
# Figura 18: sensibilidad a la posicion de la superficie
# ---------------------------------------------------------------------------

METRICAS_SENS = (("BVTV", "BV/TV"), ("TbTh", "Tb.Th"), ("DA", "DA"),
                 ("E_app", r"$E_{app}$"))


def fig_sensibilidad_superficie(doc, destino, idioma="es"):
    i = _i(idioma)
    rec = (doc.get("resultados") or {}).get("sensibilidad_superficie")
    pe = {e: v.get("filas") or [] for e, v in _por_estructura(rec).items()}
    pe = {e: f for e, f in pe.items() if len(f) >= 3}
    if not pe:
        return []
    fig = _figura(6.6, 2.5)
    for k, (clave, etq) in enumerate(METRICAS_SENS):
        ax = fig.add_subplot(1, 4, k + 1)
        _estilo(ax, "y")
        for e, filas in pe.items():
            base = next((f for f in filas if f["desplazamiento_vox"] == 0.0),
                        None)
            v0 = _f(base, clave)
            if not v0:
                continue
            x = [f["desplazamiento_um"] for f in filas]
            y = [100.0 * ((_f(f, clave) or np.nan) / v0 - 1.0) for f in filas]
            ax.plot(x, y, "-o", color=COLOR[e], ms=3, lw=1.1,
                    label=NOMBRES[e][i])
        ax.axhline(0, color=EJE, lw=0.8)
        ax.axvline(0, color=EJE, lw=0.8)
        ax.tick_params(labelsize=6)
        ax.set_xlabel(("desplazamiento [µm]", "offset [µm]")[i],
                      fontsize=6.8, color=TINTA_2)
        if k == 0:
            ax.set_ylabel(("cambio frente a la original [%]",
                           "change from the original [%]")[i], fontsize=6.8,
                          color=TINTA_2)
        _titulo(ax, f"({'abcd'[k]}) {etq}")
    h, l = fig.axes[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=6.5, ncol=len(l),
               loc="lower center", bbox_to_anchor=(0.5, 0.98),
               labelcolor=TINTA)
    fig.tight_layout(w_pad=0.8)
    return _guardar(fig, destino)


# ---------------------------------------------------------------------------
# Figura 19: condiciones de contorno
# ---------------------------------------------------------------------------

CONDICIONES = (("deslizante", ("apoyo deslizante", "sliding support"),
                None),
               ("empotrado", ("apoyo empotrado", "fixed support"), "////"),
               ("periodica", ("periódica", "periodic"), "...."))


def fig_condiciones_contorno(doc, destino, idioma="es"):
    i = _i(idioma)
    rec = (doc.get("resultados") or {}).get("condiciones_contorno")
    pe = _por_estructura(rec)
    pe = {e: v for e, v in pe.items()
          if any(_f(v.get(c) or {}, "E_rel") for c, _e, _h in CONDICIONES)}
    if not pe:
        return []
    fig = _figura(4.6, 2.7)
    ax = fig.add_subplot(1, 1, 1)
    _estilo(ax, "y")
    ancho = 0.26
    ymax = 0.0
    for k, (e, v) in enumerate(pe.items()):
        per = _f(v.get("periodica") or {}, "E_rel")
        for j, (c, _etq, patron) in enumerate(CONDICIONES):
            y = _f(v.get(c) or {}, "E_rel")
            if y is None:
                continue
            xpos = k + (j - 1) * ancho
            ax.bar(xpos, y, ancho, color=COLOR[e] if patron is None else
                   FONDO, edgecolor=COLOR[e], linewidth=1.0, hatch=patron)
            ymax = max(ymax, y)
            txt = _n(y, idioma, ".3g")
            if per and c != "periodica":
                txt += f"\n×{_n(y / per, idioma, '.2f')}"
            ax.text(xpos, y, txt, ha="center", va="bottom", fontsize=5.6,
                    color=TINTA_2)
    ax.set_ylim(0, ymax * 1.3)
    ax.set_xticks(range(len(pe)))
    ax.set_xticklabels([f"{NOMBRES[e][i]}\n{v.get('n', '')}³"
                        for e, v in pe.items()], fontsize=7)
    for t, e in zip(ax.get_xticklabels(), pe):
        t.set_color(COLOR[e])
    ax.set_ylabel(r"$E/E_s$", fontsize=7.5, color=TINTA_2)
    ax.legend(handles=[Patch(facecolor=TINTA_2 if p is None else FONDO,
                             edgecolor=TINTA_2, hatch=p, label=etq[i])
                       for _c, etq, p in CONDICIONES],
              frameon=False, fontsize=6.3, labelcolor=TINTA, ncol=3,
              loc="lower center", bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout()
    return _guardar(fig, destino)


FIGURAS = (
    ("perdida", ("fig12_perdida_osea", "fig12_bone_loss"), fig_perdida),
    ("fallo", ("fig13_fallo_progresivo", "fig13_progressive_failure"),
     fig_fallo),
    ("superficie_E", ("fig14_superficie_E", "fig14_E_surface"),
     fig_superficie_E),
    ("curvaturas", ("fig15_curvaturas", "fig15_curvatures"), fig_curvaturas),
    ("ellipsoid", ("fig16_ellipsoid_factor", "fig16_ellipsoid_factor_en"),
     fig_ellipsoid_factor),
    ("febio", ("fig17_febio", "fig17_febio_en"), fig_paridad_febio),
    ("sensibilidad", ("fig18_sensibilidad_superficie",
                      "fig18_surface_sensitivity"),
     fig_sensibilidad_superficie),
    ("contorno", ("fig19_condiciones_contorno", "fig19_boundary_conditions"),
     fig_condiciones_contorno),
)
