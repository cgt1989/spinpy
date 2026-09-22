"""
figuras.py — Figuras del informe para publicacion, en calidad de imprenta.

QUE PRODUCE
-----------
  fig1_morfometria   cociente candidato / VOI de cada metrica (VOI = 1)
  fig2_mecanica      modulos por homogeneizacion (Ex, Ey, Ez) y E_app por eje,
                     VOI frente a candidato; lo NO citable va rayado y lo
                     citable con reservas lleva *
  fig3_estructuras   isosuperficie 3D del VOI y de cada estructura ajustada, en
                     el estilo y la vista principales
  fig4_vistas        rejilla estructura x vista en el estilo principal
  fig5_secciones     cortes centrales XY, XZ e YZ, hueso en tinta, 1 mm de barra
  fig6_anisotropia   cortes polares de la fabrica MIL y de E(n) del tensor C
  fig7_distribucion  espesor trabecular local y tamano de poro (opcional)
  fig8_von_mises     tension de von Mises del analisis comparado en 3D, con
                     una escala comun para todas las estructuras (opcional)
  3d/                galeria: cada estructura en cada estilo y cada vista
                     elegidos, para cambiar la figura sin volver a calcular

Las de datos salen en PNG a 600 ppp —lo que piden la mayoria de revistas para
graficos de linea— y en PDF vectorial. Cada render 3D principal sale en PNG de
2400 px de lado, los de la galeria a 1600 px, y los paneles en PNG a 600 ppp.

Los renders usan proyeccion PARALELA y una camara que solo depende del cubo:
VOI y candidatos, escalados al mismo lado, salen a la misma escala en cualquier
vista. La superficie se puede suavizar con Taubin, que no encoge ni adelgaza;
es solo para la figura y el pie lo declara.

POR QUE COCIENTES EN LA FIGURA 1
--------------------------------
BV/TV es adimensional, Tb.Th va en mm y Conn.D en 1/mm^3. En un solo eje no
caben sin inventar una escala para cada una, y un segundo eje engaña. El
cociente con el VOI las pone todas en la misma escala y responde a la pregunta
del ajuste —cuanto se parece—; el valor absoluto del VOI va en la etiqueta y
todos los numeros estan en la tabla del informe.

COLOR
-----
Identidad por estructura, fija: VOI azul, spinodoide naranja, dual-lattice
aguamarina. Son las tres primeras posiciones de una paleta categorica validada
para daltonismo con todos los pares a la vez (distancia CVD >= 9). El color no
va nunca solo: hay leyenda, etiqueta de valor y, para la citabilidad, rayado.

HILOS
-----
Las figuras de datos se construyen con `matplotlib.figure.Figure` y su propio
lienzo Agg, no con pyplot: pyplot elige un backend global que dentro del visor
es el de Qt. Asi se pueden dibujar en un hilo de trabajo.

El 3D NO: VTK necesita el contexto OpenGL del hilo que lo crea, y el visor ya
tiene uno en el hilo de la interfaz. `renders_3d` se llama desde ese hilo.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

COLOR = {"voi": "#2a78d6", "spinodoide": "#eb6834", "dual-lattice": "#1baf7a"}
TINTA = "#0b0b0b"
TINTA_2 = "#52514e"
REJILLA = "#e1e0d9"
EJE = "#c3c2b7"
FONDO = "#ffffff"

DPI = 600
LADO_3D = 2400
ORDEN = ("voi", "spinodoide", "dual-lattice")

NOMBRES = {"voi": ("VOI", "VOI"), "spinodoide": ("Spinodoide", "Spinodoid"),
           "dual-lattice": ("Dual-lattice", "Dual-lattice")}

METRICAS = (("BVTV", "BV/TV", ""), ("BSBV", "BS/BV", "mm⁻¹"),
            ("TbTh", "Tb.Th", "mm"), ("TbSp", "Tb.Sp", "mm"),
            ("TbN", "Tb.N", "mm⁻¹"), ("DA", "DA", ""),
            ("ConnD", "Conn.D", "mm⁻³"), ("SMI", "SMI", ""))


def _i(idioma):
    return 0 if idioma == "es" else 1


def _n(x, idioma, espec=".3g"):
    x = float(x)
    # Por encima de 100 se redondea a entero en vez de dejar que `g` pase a
    # notacion cientifica: "2950" se lee en una barra, "2,95e+03" no.
    s = format(x, ".0f" if abs(x) >= 100 else espec)
    return s.replace(".", ",") if idioma == "es" else s


def _f(d, clave):
    try:
        x = float(np.asarray((d or {}).get(clave), float).reshape(()))
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def _estilo(ax, rejilla="x"):
    ax.set_facecolor(FONDO)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(EJE)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=TINTA_2, labelsize=7.5, length=3, width=0.6)
    ax.grid(axis=rejilla, color=REJILLA, linewidth=0.6)
    ax.set_axisbelow(True)


def _figura(ancho, alto):
    fig = Figure(figsize=(ancho, alto), facecolor=FONDO)
    FigureCanvasAgg(fig)
    return fig


def _guardar(fig, destino, pdf=True):
    destino = Path(destino)
    destino.parent.mkdir(parents=True, exist_ok=True)
    rutas = [destino.with_suffix(".png")]
    fig.savefig(rutas[0], dpi=DPI, facecolor=FONDO, bbox_inches="tight")
    if pdf:
        rutas.append(destino.with_suffix(".pdf"))
        fig.savefig(rutas[1], facecolor=FONDO, bbox_inches="tight")
    return [str(r) for r in rutas]


# ---------------------------------------------------------------------------
# Figura 1 — morfometria relativa
# ---------------------------------------------------------------------------

def fig_morfometria(doc, destino, idioma="es"):
    """Cociente candidato/VOI por metrica. Devuelve [png, pdf] o [] si no hay."""
    i = _i(idioma)
    m_v = doc.get("morfometria_voi") or {}
    cands = [(f, doc.get(k)) for f, k in (("spinodoide", "morfometria_spin"),
                                          ("dual-lattice", "morfometria_dual"))
             if doc.get(k)]
    filas = []
    for clave, etq, uni in METRICAS:
        v = _f(m_v, clave)
        if v and any(_f(mc, clave) is not None for _f_, mc in cands):
            filas.append((clave, etq, uni, v))
    if not filas or not cands:
        return []

    n, k = len(filas), len(cands)
    fig = _figura(3.5, 0.34 * n * k + 0.35 * n + 0.8)
    ax = fig.add_subplot(111)
    _estilo(ax, "x")
    alto = 0.78 / k
    y0 = np.arange(n)[::-1].astype(float)
    maximo = 1.0
    for j, (fam, mc) in enumerate(cands):
        r = np.array([(_f(mc, c) / v) if _f(mc, c) is not None else np.nan
                      for c, _e, _u, v in filas])
        y = y0 + (k - 1) / 2 * alto - j * alto
        ax.barh(y, r, height=alto * 0.86, color=COLOR[fam],
                label=NOMBRES[fam][i], zorder=2)
        for yy, rr in zip(y, r):
            if np.isfinite(rr):
                # A la derecha de la barra Y de la linea de referencia: una
                # barra algo menor que 1 dejaba el numero pisado por la linea.
                ax.text(max(rr, 1.0), yy, "  " + _n(rr, idioma), va="center",
                        ha="left", fontsize=6.5, color=TINTA_2)
        if np.isfinite(r).any():
            maximo = max(maximo, float(np.nanmax(r)))
    ax.axvline(1.0, color=TINTA, linewidth=1.0, zorder=3)
    etiquetas = [f"{etq}\nVOI {_n(v, idioma)}{(' ' + uni) if uni else ''}"
                 for _c, etq, uni, v in filas]
    ax.set_yticks(y0)
    ax.set_yticklabels(etiquetas, fontsize=7, color=TINTA)
    ax.set_xlim(0, maximo * 1.28)
    ax.set_xlabel(("candidato / VOI", "candidate / VOI")[i], fontsize=7.5,
                  color=TINTA_2)
    manejadores = [Patch(color=COLOR[f], label=NOMBRES[f][i]) for f, _m in cands]
    manejadores.append(Line2D([], [], color=TINTA, linewidth=1.0,
                              label=("VOI (referencia = 1)",
                                     "VOI (reference = 1)")[i]))
    ax.legend(handles=manejadores, frameon=False, fontsize=6.8, ncol=k + 1,
              loc="lower left", bbox_to_anchor=(0.0, 1.0),
              handlelength=1.4, columnspacing=1.0, labelcolor=TINTA)
    return _guardar(fig, destino)


# ---------------------------------------------------------------------------
# Figura 2 — mecanica, con la citabilidad a la vista
# ---------------------------------------------------------------------------

def _estado(items, magnitud, estructura, eje=None):
    for it in items or []:
        if (it["magnitud"] == magnitud and it["estructura"] == estructura
                and it.get("eje") == eje):
            return it["estado"]
    return "citable"


def fig_mecanica(doc, items, destino, idioma="es"):
    """E por homogeneizacion y E_app por eje. Devuelve [png, pdf] o []."""
    i = _i(idioma)
    res = doc.get("resultados") or {}
    paneles = []

    serie = {}
    for fam, suf in (("spinodoide", ""), ("dual-lattice", "_dual")):
        rec = res.get("elastico" + suf)
        if not isinstance(rec, dict):
            continue
        for k, v in (rec.get("por_estructura") or {}).items():
            est = "voi" if k == "voi" else fam
            if est in serie:
                continue
            vals = {e: _f(v, e) for e in ("Ex", "Ey", "Ez")}
            if any(x is not None for x in vals.values()):
                estado = _estado(items, "elastico", est)
                serie[est] = {e: (x / 1e6, estado) for e, x in vals.items()
                              if x is not None}
    if serie:
        paneles.append((("Homogeneización periódica",
                         "Periodic homogenisation")[i],
                        ["Ex", "Ey", "Ez"], serie, None))

    serie = {}
    for fam, suf in (("spinodoide", ""), ("dual-lattice", "_dual")):
        rec = res.get("resistencia" + suf)
        if not isinstance(rec, dict):
            continue
        for k, v in (rec.get("por_estructura") or {}).items():
            est = "voi" if k == "voi" else fam
            for eje, p in (v.get("ejes") or {}).items():
                E = _f(p, "E_app")
                if E is not None:
                    serie.setdefault(est, {})[eje] = (
                        E / 1e6, _estado(items, "E_app", est, eje))
    if serie:
        cats = sorted({e for s in serie.values() for e in s})
        paneles.append((("Ensayo de compresión", "Compression test")[i],
                        cats, serie, ("eje de carga", "loading axis")[i]))
    if not paneles:
        return []

    ETIQ = {"Ex": "$E_x$", "Ey": "$E_y$", "Ez": "$E_z$"}
    fig = _figura(3.3 * len(paneles) + 0.3, 2.8)
    presentes, hay_rayado, hay_reservas = [], False, False
    for n_p, (titulo, cats, serie, xlab) in enumerate(paneles):
        ax = fig.add_subplot(1, len(paneles), n_p + 1)
        _estilo(ax, "y")
        ests = [e for e in ORDEN if e in serie]
        k = len(ests)
        ancho = 0.8 / k
        x = np.arange(len(cats), dtype=float)
        ymax = 0.0
        for j, est in enumerate(ests):
            if est not in presentes:
                presentes.append(est)
            for xi, cat in zip(x - 0.4 + ancho * (j + 0.5), cats):
                if cat not in serie[est]:
                    continue
                val, estado = serie[est][cat]
                ymax = max(ymax, val)
                if estado == "no_citable":
                    hay_rayado = True
                    ax.bar(xi, val, width=ancho * 0.86, facecolor=FONDO,
                           edgecolor=COLOR[est], hatch="////", linewidth=0.8,
                           zorder=2)
                else:
                    ax.bar(xi, val, width=ancho * 0.86, color=COLOR[est],
                           zorder=2)
                marca = "*" if estado == "reservas" else ""
                hay_reservas = hay_reservas or bool(marca)
                ax.text(xi, val, _n(val, idioma) + marca, ha="center",
                        va="bottom", fontsize=6, color=TINTA_2)
        ax.set_xticks(x)
        ax.set_xticklabels([ETIQ.get(c, c) for c in cats], fontsize=7.5,
                           color=TINTA)
        ax.set_ylim(0, ymax * 1.2 if ymax > 0 else 1)
        ax.set_ylabel("E [MPa]", fontsize=7.5, color=TINTA_2)
        if xlab:
            ax.set_xlabel(xlab, fontsize=7.5, color=TINTA_2)
        ax.set_title(titulo, fontsize=8, color=TINTA, loc="left")
    manejadores = [Patch(color=COLOR[e], label=NOMBRES[e][i])
                   for e in presentes]
    if hay_rayado:
        manejadores.append(Patch(facecolor=FONDO, edgecolor=TINTA_2,
                                 hatch="////", label=("no citable",
                                                      "not citable")[i]))
    if hay_reservas:
        manejadores.append(Line2D([], [], linestyle="none", marker="$*$",
                                  color=TINTA_2, label=("con reservas",
                                                        "with caveats")[i]))
    fig.legend(handles=manejadores, frameon=False, fontsize=6.8,
               ncol=len(manejadores), loc="lower center",
               bbox_to_anchor=(0.5, 0.98), labelcolor=TINTA)
    fig.tight_layout()
    return _guardar(fig, destino)


# ---------------------------------------------------------------------------
# Figura 3 — estructuras en 3D
# ---------------------------------------------------------------------------

def malla_cerrada(BW, spacing):
    """Isosuperficie del solido, CERRADA contra las caras del cubo.

    Contornear la mascara tal cual deja abiertas las caras donde el hueso toca
    el borde: es lo correcto para MEDIR (esas tapas no son interfase osea), pero
    en una figura hace que la pieza parezca mucho mas porosa de lo que es. Se
    rodea con una capa de fondo antes de contornear. Cambia la figura, no la
    medida.
    """
    import pyvista as pv
    BW = np.asarray(BW, dtype=bool)
    lleno = np.zeros(np.array(BW.shape) + 2, dtype=np.float32)
    lleno[1:-1, 1:-1, 1:-1] = BW
    grid = pv.ImageData(dimensions=lleno.shape,
                        spacing=tuple(float(s) for s in spacing),
                        origin=tuple(-float(s) for s in spacing))
    grid.point_data["v"] = lleno.flatten(order="F")
    malla = grid.contour([0.5], scalars="v")
    if malla.n_points:
        malla = malla.compute_normals(auto_orient_normals=True,
                                      consistent_normals=False)
    return malla



# Estilos de render. `color=None` toma el color de identidad de la estructura;
# `escalar="z"` colorea por altura. Gris y hueso son los de las revistas de
# hueso y biomateriales: sin color que distraiga de la forma, y legibles en
# escala de grises impresa. La altura usa cividis, secuencial y uniforme en
# luminancia, que se lee igual con cualquier tipo de daltonismo.
ESTILOS_3D = {
    "gris": {"color": "#c9c8c3", "nombre": ("gris", "grey")},
    "hueso": {"color": "#e4d8bf", "nombre": ("hueso (marfil)", "bone (ivory)")},
    "identidad": {"color": None,
                  "nombre": ("color por estructura", "colour by structure")},
    "altura": {"escalar": "z", "cmap": "cividis",
               "nombre": ("color por altura (Z)", "colour by height (Z)")},
}

# Vistas: direccion de la camara DESDE el centro del cubo y vector arriba. Las
# isometricas giran de 90 en 90 grados alrededor de Z, que es el eje axial del
# VOI; las ortogonales miran a lo largo de un eje, con Z arriba salvo la
# superior, que mira desde +Z con Y arriba.
VISTAS_3D = {
    "iso": ((1, 1, 1), (0, 0, 1),
            ("isométrica", "isometric")),
    "iso_90": ((-1, 1, 1), (0, 0, 1),
               ("isométrica girada 90°", "isometric rotated 90°")),
    "iso_180": ((-1, -1, 1), (0, 0, 1),
                ("isométrica girada 180°", "isometric rotated 180°")),
    "iso_270": ((1, -1, 1), (0, 0, 1),
                ("isométrica girada 270°", "isometric rotated 270°")),
    "frontal": ((0, -1, 0), (0, 0, 1), ("frontal (XZ)", "front (XZ)")),
    "lateral": ((1, 0, 0), (0, 0, 1), ("lateral (YZ)", "side (YZ)")),
    "superior": ((0, 0, 1), (0, 1, 0), ("superior (XY)", "top (XY)")),
}

LADO_GALERIA = 1600
# Suavizado de Taubin de las figuras (nunca de las medidas). La seccion de
# modelos del informe lee estos dos valores.
TAUBIN_ITER = 40
TAUBIN_BANDA = 0.03
# Por omision: el estilo y la vista de la figura 3, y lo que se renderiza para
# la galeria y el panel de vistas.
OPCIONES_3D = {"estilos": ["gris", "identidad"],
               "vistas": ["iso", "iso_180", "frontal", "lateral", "superior"],
               "estilo_principal": "gris", "vista_principal": "iso",
               "suavizar": True, "distribuciones": True, "von_mises": True}


def opciones_3d(op=None):
    """Opciones completas y validas: lo desconocido se descarta, lo que
    falta toma el valor por omision, y el estilo y la vista principales entran
    siempre en la galeria."""
    o = dict(OPCIONES_3D)
    o.update({k: v for k, v in (op or {}).items() if k in OPCIONES_3D})
    o["estilos"] = [e for e in o["estilos"] if e in ESTILOS_3D]
    o["vistas"] = [v for v in o["vistas"] if v in VISTAS_3D]
    if o["estilo_principal"] not in ESTILOS_3D:
        o["estilo_principal"] = OPCIONES_3D["estilo_principal"]
    if o["vista_principal"] not in VISTAS_3D:
        o["vista_principal"] = OPCIONES_3D["vista_principal"]
    if o["estilo_principal"] not in o["estilos"]:
        o["estilos"].insert(0, o["estilo_principal"])
    if o["vista_principal"] not in o["vistas"]:
        o["vistas"].insert(0, o["vista_principal"])
    for k in ("suavizar", "distribuciones", "von_mises"):
        o[k] = bool(o[k])
    return o


def malla_figura(BW, spacing, suavizar=True):
    """`malla_cerrada` y, si se pide, suavizada con Taubin.

    El contorno de una mascara binaria sale en escalera, y en un render a
    2400 px la escalera es lo primero que se ve. Taubin alterna un paso que
    suaviza y otro que infla, asi que NO encoge la pieza ni adelgaza las
    trabeculas finas, que es lo que haria un Laplaciano o un desenfoque del
    volumen antes de contornear. Es solo para la figura: toda medida usa la
    malla sin tocar, y el pie de la figura lo declara.
    """
    malla = malla_cerrada(BW, spacing)
    if suavizar and malla.n_points:
        try:
            malla = malla.smooth_taubin(n_iter=TAUBIN_ITER,
                                        pass_band=TAUBIN_BANDA)
            malla = malla.compute_normals(auto_orient_normals=True,
                                          consistent_normals=False)
        except Exception:
            pass
    return malla


def _cubo(BW, spacing):
    """Caja nominal del cubo de voxeles, en mm. La misma para todas las
    estructuras escaladas al VOI, y por eso la misma camara."""
    n = np.asarray(np.shape(BW), float)
    s = np.asarray(spacing, float)
    return -0.5 * s, (n - 0.5) * s


def _camara(p, lo, hi, vista):
    """Proyeccion PARALELA con una escala que solo depende del cubo.

    En perspectiva el tamano aparente depende de la distancia y cada vista
    enfoca distinto; en paralela, dos estructuras del mismo lado salen a la
    misma escala en cualquier vista, que es lo que permite ponerlas lado a
    lado en un panel y compararlas a ojo.
    """
    d, arriba, _nom = VISTAS_3D[vista]
    d = np.asarray(d, float)
    d /= np.linalg.norm(d)
    c = 0.5 * (lo + hi)
    diag = float(np.linalg.norm(hi - lo))
    arriba = np.asarray(arriba, float)
    u = arriba - np.dot(arriba, d) * d
    u /= np.linalg.norm(u)
    r = np.cross(d, u)
    esquinas = np.array([[x, y, z] for x in (lo[0], hi[0])
                         for y in (lo[1], hi[1]) for z in (lo[2], hi[2])]) - c
    alto = np.ptp(esquinas @ u)
    ancho = np.ptp(esquinas @ r)
    cam = p.camera
    cam.focal_point = tuple(c)
    cam.position = tuple(c + 3.0 * diag * d)
    cam.up = tuple(u)
    cam.parallel_projection = True
    cam.parallel_scale = 0.5 * max(alto, ancho) * 1.04
    p.renderer.ResetCameraClippingRange()
    # Sin render explicito, `screenshot` devuelve el ultimo fotograma: todas
    # las vistas de un mismo plotter salian iguales a la primera.
    p.render()


def _anadir(p, malla, estilo, estructura, lo, hi):
    e = ESTILOS_3D[estilo]
    if e.get("escalar") == "z":
        m = malla.copy()
        m["z"] = m.points[:, 2]
        p.add_mesh(m, scalars="z", cmap=e["cmap"], clim=(lo[2], hi[2]),
                   show_scalar_bar=False, smooth_shading=True,
                   specular=0.2, specular_power=15)
    else:
        p.add_mesh(malla, color=e["color"] or COLOR[estructura],
                   smooth_shading=True, specular=0.2, specular_power=15)
    diag = float(np.linalg.norm(hi - lo))
    # Oclusion ambiental: sin sombra de contacto no se sabe que trabecula pasa
    # por delante de cual. Radio a escala de la escena, como en el visor.
    try:
        p.enable_ssao(radius=0.045 * diag, bias=0.0009 * diag,
                      kernel_size=128, blur=True)
    except Exception:
        pass
    try:
        p.enable_anti_aliasing("ssaa")
    except Exception:
        pass


def _plotter(lado):
    import pyvista as pv
    p = pv.Plotter(off_screen=True, window_size=(lado, lado))
    p.set_background(FONDO)
    return p


def render_3d(BW, spacing, color, destino, lado=LADO_3D):
    """Un render isometrico de un color. Se conserva por compatibilidad."""
    destino = Path(destino)
    destino.parent.mkdir(parents=True, exist_ok=True)
    lo, hi = _cubo(BW, spacing)
    p = _plotter(lado)
    try:
        p.add_mesh(malla_cerrada(BW, spacing), color=color,
                   smooth_shading=True, specular=0.25, specular_power=15)
        try:
            p.enable_anti_aliasing("ssaa")
        except Exception:
            pass
        _camara(p, lo, hi, "iso")
        p.screenshot(str(destino))
    finally:
        p.close()
    return str(destino)


def renders_3d(estructuras, carpeta, opciones=None):
    """Renders de cada estructura. Hilo de la interfaz (OpenGL).

    Devuelve {"principal": {estructura: png}, "galeria": {(est, estilo,
    vista): png}}. La principal —estilo y vista elegidos, a 2400 px— es la
    de la figura 3 y conserva el nombre `fig3_<estructura>.png`; la galeria,
    a 1600 px, va a `3d/` para que el autor elija otra sin volver a calcular.
    Cada estructura se contornea UNA vez y cada estilo usa un solo plotter
    para todas sus vistas: solo cambia la camara.
    """
    o = opciones_3d(opciones)
    carpeta = Path(carpeta)
    (carpeta / "3d").mkdir(parents=True, exist_ok=True)
    out = {"principal": {}, "galeria": {}}
    for est, BW, sp in estructuras:
        malla = malla_figura(BW, sp, o["suavizar"])
        lo, hi = _cubo(BW, sp)
        p = _plotter(LADO_3D)
        try:
            _anadir(p, malla, o["estilo_principal"], est, lo, hi)
            _camara(p, lo, hi, o["vista_principal"])
            ruta = carpeta / f"fig3_{est}.png"
            p.screenshot(str(ruta))
            out["principal"][est] = str(ruta)
        finally:
            p.close()
        for estilo in o["estilos"]:
            p = _plotter(LADO_GALERIA)
            try:
                _anadir(p, malla, estilo, est, lo, hi)
                for vista in o["vistas"]:
                    _camara(p, lo, hi, vista)
                    ruta = carpeta / "3d" / f"{est}_{estilo}_{vista}.png"
                    p.screenshot(str(ruta))
                    out["galeria"][(est, estilo, vista)] = str(ruta)
            finally:
                p.close()
    return out


def panel_3d(rutas, estructuras, destino, idioma="es"):
    """Reune los renders en un panel rotulado (a), (b), (c)."""
    from matplotlib.image import imread
    i = _i(idioma)
    k = len(rutas)
    fig = _figura(2.4 * k, 2.65)
    for j, (ruta, est) in enumerate(zip(rutas, estructuras)):
        ax = fig.add_subplot(1, k, j + 1)
        ax.imshow(imread(ruta))
        ax.set_axis_off()
        ax.set_title(f"({'abc'[j]}) {NOMBRES[est][i]}", fontsize=8.5,
                     color=TINTA, loc="left")
    fig.tight_layout()
    return _guardar(fig, destino, pdf=False)[0]


def panel_vistas(galeria, estructuras, estilo, vistas, destino, idioma="es"):
    """Rejilla estructura x vista en un solo estilo (figura 4).

    Filas = estructuras y columnas = vistas, para que cada columna compare
    las estructuras desde el mismo sitio. Con proyeccion paralela y el mismo
    cubo, todas las celdas estan a la misma escala.
    """
    from matplotlib.image import imread
    i = _i(idioma)
    vistas = [v for v in vistas
              if all((e, estilo, v) in galeria for e in estructuras)]
    if not vistas or not estructuras:
        return None
    filas, cols = len(estructuras), len(vistas)
    fig = _figura(1.55 * cols + 0.2, 1.62 * filas + 0.25)
    letras = "abcdefghijklmnopqrstuvwxyz"
    for a, est in enumerate(estructuras):
        for b, vista in enumerate(vistas):
            ax = fig.add_subplot(filas, cols, a * cols + b + 1)
            ax.imshow(imread(galeria[(est, estilo, vista)]))
            ax.set_axis_off()
            if a == 0:
                ax.set_title(VISTAS_3D[vista][2][i], fontsize=6.5,
                             color=TINTA_2)
            if b == 0:
                ax.text(-0.04, 0.5, NOMBRES[est][i], transform=ax.transAxes,
                        rotation=90, ha="right", va="center", fontsize=7.5,
                        color=TINTA)
            ax.text(0.02, 0.98, f"({letras[a * cols + b]})",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=6.5, color=TINTA,
                    bbox=dict(boxstyle="square,pad=0.15", fc=FONDO,
                              ec="none", alpha=0.85))
    fig.subplots_adjust(left=0.04, right=0.995, top=0.93, bottom=0.01,
                        wspace=0.02, hspace=0.04)
    return _guardar(fig, destino, pdf=False)[0]


# ---------------------------------------------------------------------------
# Figura 5 — secciones ortogonales
# ---------------------------------------------------------------------------

def fig_secciones(estructuras, destino, idioma="es"):
    """Cortes centrales XY, XZ e YZ de cada estructura, a la misma escala.

    Es la figura clasica de micro-CT y la que un revisor pide para ver la
    segmentacion: el 3D ensena la forma, el corte ensena el grosor y el poro
    sin la ambiguedad de la perspectiva. Hueso en tinta sobre blanco, que se
    imprime bien en escala de grises. Sin matplotlib 3D ni VTK: corre en un
    hilo de trabajo. La mascara tiene dimension 1 = X, 2 = Y, 3 = Z.
    """
    if not estructuras:
        return []
    from matplotlib.colors import ListedColormap
    i = _i(idioma)
    cmap = ListedColormap([FONDO, TINTA])
    planos = (("XY", ("corte XY (Z central)", "XY section (mid Z)")),
              ("XZ", ("corte XZ (Y central)", "XZ section (mid Y)")),
              ("YZ", ("corte YZ (X central)", "YZ section (mid X)")))
    filas = len(estructuras)
    fig = _figura(6.3, 2.2 * filas + 0.2)
    letras = "abcdefghijklmnopqrstuvwxyz"
    for a, (est, BW, sp) in enumerate(estructuras):
        BW = np.asarray(BW, bool)
        sp = np.asarray(sp, float)
        L = np.asarray(BW.shape) * sp
        cx, cy, cz = (s // 2 for s in BW.shape)
        cortes = {"XY": (BW[:, :, cz].T, (0, L[0], 0, L[1]), ("X", "Y")),
                  "XZ": (BW[:, cy, :].T, (0, L[0], 0, L[2]), ("X", "Z")),
                  "YZ": (BW[cx, :, :].T, (0, L[1], 0, L[2]), ("Y", "Z"))}
        for b, (clave, titulos) in enumerate(planos):
            img, ext, (eh, ev) = cortes[clave]
            ax = fig.add_subplot(filas, 3, a * 3 + b + 1)
            ax.imshow(img, cmap=cmap, origin="lower", extent=ext,
                      interpolation="nearest", vmin=0, vmax=1)
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ax.spines.values():
                s.set_color(EJE)
                s.set_linewidth(0.6)
            ax.set_xlabel(eh, fontsize=6.5, color=TINTA_2, labelpad=1)
            ax.set_ylabel(ev, fontsize=6.5, color=TINTA_2, labelpad=1)
            if a == 0:
                ax.set_title(titulos[i], fontsize=7, color=TINTA_2)
            if b == 0:
                ax.text(-0.2, 0.5, NOMBRES[est][i], transform=ax.transAxes,
                        rotation=90, ha="right", va="center", fontsize=8,
                        color=TINTA)
            ax.text(0.02, 0.98, f"({letras[a * 3 + b]})",
                    transform=ax.transAxes, ha="left", va="top", fontsize=6.5,
                    color=TINTA, bbox=dict(boxstyle="square,pad=0.15",
                                           fc=FONDO, ec="none", alpha=0.85))
            if b == 0:
                # Barra de escala de 1 mm (o de un quinto del lado si el cubo
                # es menor que 2 mm), abajo a la izquierda.
                barra = 1.0 if L[0] >= 2.0 else float(f"{L[0] / 5:.1g}")
                x0, y0 = 0.06 * ext[1], 0.06 * ext[3]
                ax.plot([x0, x0 + barra], [y0, y0], color="#d1453b",
                        linewidth=2.2, solid_capstyle="butt")
                ax.text(x0 + barra / 2, y0 + 0.035 * ext[3],
                        (f"{barra:g} mm").replace(".", ",")
                        if idioma == "es" else f"{barra:g} mm",
                        ha="center", va="bottom", fontsize=6,
                        color="#d1453b", fontweight="bold",
                        bbox=dict(boxstyle="square,pad=0.12", fc=FONDO,
                                  ec="none", alpha=0.85))
    fig.tight_layout(pad=0.4)
    return _guardar(fig, destino, pdf=False)


# ---------------------------------------------------------------------------
# Figura 6 — anisotropia: fabrica (MIL) y modulo de Young direccional
# ---------------------------------------------------------------------------

PLANOS = (("XY", (0, 1)), ("XZ", (0, 2)), ("YZ", (1, 2)))


def _mil_plano(m, ejes, ang):
    """MIL(n) = 1/sqrt(n^T M n) del tensor de fabrica, en un plano.

    M se reconstruye de sus autovalores y autovectores, que es lo que guarda la
    morfometria: con la convencion 1/MIL(n)^2 = n^T M n del tensor MIL. Devuelve
    None si no hay tensor.
    """
    lam = m.get("eigenvalues") if m else None
    vec = m.get("eigenvectors") if m else None
    if lam is None or vec is None:
        return None
    lam = np.asarray(lam, float)
    vec = np.asarray(vec, float)
    if lam.shape != (3,) or vec.shape != (3, 3) or not np.all(np.isfinite(lam)):
        return None
    M = vec @ np.diag(lam) @ vec.T
    n = np.zeros((ang.size, 3))
    n[:, ejes[0]] = np.cos(ang)
    n[:, ejes[1]] = np.sin(ang)
    q = np.einsum("ij,jk,ik->i", n, M, n)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(q > 0, 1.0 / np.sqrt(q), np.nan)


def modulo_direccional(C, n):
    """E(n) de un tensor de rigidez en Voigt [xx yy zz yz xz xy], en las
    unidades de C.

    Con deformaciones de cizalla INGENIERILES (la convencion de todo el
    proyecto), una traccion uniaxial s n(x)n tiene sigma_voigt = s a con
    a = [n1², n2², n3², n2 n3, n1 n3, n1 n2], y la deformacion normal a lo
    largo de n es eps_voigt . a. Por tanto 1/E(n) = a^T S a con S = C^-1.
    """
    S = np.linalg.inv(np.asarray(C, float))
    n = np.atleast_2d(np.asarray(n, float))
    a = np.column_stack([n[:, 0] ** 2, n[:, 1] ** 2, n[:, 2] ** 2,
                         n[:, 1] * n[:, 2], n[:, 0] * n[:, 2],
                         n[:, 0] * n[:, 1]])
    return 1.0 / np.einsum("ij,jk,ik->i", a, S, a)


def fig_anisotropia(doc, destino, idioma="es"):
    """Cortes polares de la fabrica y de E(n) en XY, XZ e YZ.

    Fila 1: MIL(n) en mm, del tensor de fabrica de cada estructura. Fila 2:
    modulo de Young direccional E(n) en MPa, del tensor C homogeneizado. Las
    dos filas preguntan lo mismo —en que direccion es «fuerte» la estructura—
    una por la geometria y otra por la mecanica, y ponerlas juntas deja ver
    si coinciden. Devuelve [png, pdf] o [] si no hay nada que dibujar.
    """
    i = _i(idioma)
    res = doc.get("resultados") or {}
    morf = {"voi": doc.get("morfometria_voi"),
            "spinodoide": doc.get("morfometria_spin"),
            "dual-lattice": doc.get("morfometria_dual")}
    rig = {}
    for fam, suf in (("spinodoide", ""), ("dual-lattice", "_dual")):
        pe = (res.get("elastico" + suf) or {}).get("por_estructura") or {}
        for k, v in pe.items():
            est = "voi" if k == "voi" else fam
            if est not in rig and v.get("C_Pa") is not None:
                try:
                    C = np.asarray(v["C_Pa"], float)
                    if C.shape == (6, 6) and np.all(np.isfinite(C)):
                        np.linalg.inv(C)
                        rig[est] = C / 1e6
                except (ValueError, np.linalg.LinAlgError):
                    pass
    ang = np.linspace(0.0, 2.0 * np.pi, 361)
    filas = []
    mil = {e: [_mil_plano(m, ej, ang) for _p, ej in PLANOS]
           for e, m in morf.items() if m}
    mil = {e: c for e, c in mil.items() if all(x is not None for x in c)}
    if mil:
        filas.append((("Fábrica MIL [mm]", "MIL fabric [mm]")[i], mil))
    if rig:
        curvas = {}
        for e, C in rig.items():
            cs = []
            for _p, ej in PLANOS:
                n = np.zeros((ang.size, 3))
                n[:, ej[0]] = np.cos(ang)
                n[:, ej[1]] = np.sin(ang)
                cs.append(modulo_direccional(C, n))
            curvas[e] = cs
        filas.append((("Módulo de Young E(n) [MPa]",
                       "Young's modulus E(n) [MPa]")[i], curvas))
    if not filas:
        return []
    fig = _figura(6.6, 2.35 * len(filas) + 0.35)
    presentes = []
    for a, (titulo, curvas) in enumerate(filas):
        rmax = max(float(np.nanmax(c)) for cs in curvas.values() for c in cs)
        for b, (plano, _ej) in enumerate(PLANOS):
            ax = fig.add_subplot(len(filas), 3, a * 3 + b + 1, polar=True)
            for est in ORDEN:
                if est not in curvas:
                    continue
                if est not in presentes:
                    presentes.append(est)
                ax.plot(ang, curvas[est][b], color=COLOR[est], linewidth=1.3)
            ax.set_ylim(0, rmax * 1.08)
            ax.set_facecolor(FONDO)
            ax.grid(color=REJILLA, linewidth=0.6)
            ax.spines["polar"].set_color(EJE)
            ax.tick_params(labelsize=5.5, colors=TINTA_2, pad=1)
            # Dos anillos rotulados, en la diagonal: con los de matplotlib
            # por omision los numeros se pisaban entre si.
            paso = float(f"{rmax / 2:.1g}")
            ax.set_yticks([paso, 2 * paso] if 2 * paso <= rmax * 1.08
                          else [paso])
            ax.yaxis.set_major_formatter(
                lambda x, _p: _n(x, idioma, ".2g"))
            ax.set_rlabel_position(45)
            h, v = plano[0], plano[1]
            ax.set_xticks(np.radians([0, 90, 180, 270]))
            ax.set_xticklabels([f"+{h}", f"+{v}", f"−{h}", f"−{v}"],
                               fontsize=6.5, color=TINTA_2)
            if b == 0:
                ax.text(-0.42, 0.5, titulo, transform=ax.transAxes,
                        rotation=90, ha="center", va="center", fontsize=7.5,
                        color=TINTA)
            if a == 0:
                ax.set_title(("plano ", "plane ")[i] + plano, fontsize=7.5,
                             color=TINTA, pad=18)
    manejadores = [Line2D([], [], color=COLOR[e], linewidth=1.5,
                          label=NOMBRES[e][i]) for e in presentes]
    fig.legend(handles=manejadores, frameon=False, fontsize=7,
               ncol=len(manejadores), loc="lower center",
               bbox_to_anchor=(0.5, 0.985), labelcolor=TINTA)
    fig.tight_layout(pad=0.6, w_pad=1.2, h_pad=1.4)
    return _guardar(fig, destino)


# ---------------------------------------------------------------------------
# Figura 7 — distribuciones de espesor trabecular y de tamano de poro
# ---------------------------------------------------------------------------

def distribuciones(estructuras):
    """{estructura: {"TbTh": valores, "PoDm": valores}} en mm.

    Espesor local por esferas inscritas (Hildebrand y Ruegsegger) sobre el
    hueso y sobre el poro, el mismo algoritmo que Po.Dm en la morfometria. Es
    la parte cara de la figura: ~70 s para un VOI de 188^3. Hilo de trabajo.
    """
    from .espesor import espesor_local
    out = {}
    for est, BW, sp in estructuras:
        BW = np.asarray(BW, bool)
        sp = np.asarray(sp, float)
        e = espesor_local(BW, sp)
        p = espesor_local(~BW, sp)
        out[est] = {"TbTh": np.asarray(e[BW], np.float32),
                    "PoDm": np.asarray(p[~BW], np.float32),
                    "h": float(np.max(sp))}
    return out


def fig_distribuciones(dist, destino, idioma="es"):
    """Histogramas (densidad) de Tb.Th local y Po.Dm, con la media marcada.

    El valor medio de Tb.Th de la tabla sale del modelo de placas
    (2 BV/BS); aqui es el espesor local, que ve la FORMA de la distribucion:
    dos estructuras con la misma media pueden tener una muy ancha y otra muy
    estrecha, y eso la tabla no lo dice.
    """
    if not dist:
        return []
    i = _i(idioma)
    fig = _figura(6.6, 2.45)
    presentes = []
    for b, (clave, titulo) in enumerate((
            ("TbTh", ("Espesor trabecular local", "Local trabecular "
                                                  "thickness")[i]),
            ("PoDm", ("Tamaño de poro local", "Local pore size")[i]))):
        ax = fig.add_subplot(1, 2, b + 1)
        _estilo(ax, "y")
        todos = np.concatenate([d[clave] for d in dist.values()
                                if d[clave].size])
        if not todos.size:
            continue
        # El espesor local toma valores discretos, multiplos del voxel. Con
        # clases mas estrechas que el voxel mas grueso el histograma sale en
        # peine; se usa ese voxel como ancho de clase, centrada en sus
        # multiplos.
        w = max(float(d.get("h", 0.0)) for d in dist.values()) or (
            float(np.percentile(todos, 99.5)) / 45)
        tope = float(np.percentile(todos, 99.5))
        bordes = (np.arange(0, int(np.ceil(tope / w)) + 2) - 0.5) * w
        bordes[0] = 0.0
        for est in ORDEN:
            if est not in dist or not dist[est][clave].size:
                continue
            if est not in presentes:
                presentes.append(est)
            v = dist[est][clave]
            ax.hist(v, bins=bordes, density=True, histtype="stepfilled",
                    color=COLOR[est], alpha=0.18, linewidth=0)
            ax.hist(v, bins=bordes, density=True, histtype="step",
                    color=COLOR[est], linewidth=1.2)
            media = float(np.mean(v))
            ax.axvline(media, color=COLOR[est], linewidth=1.0,
                       linestyle=(0, (3, 2)))
        ax.set_xlim(0.0, bordes[-1])
        ax.set_xlabel(titulo + " [mm]", fontsize=7.5, color=TINTA_2)
        ax.set_ylabel(("densidad [1/mm]", "density [1/mm]")[i], fontsize=7.5,
                      color=TINTA_2)
    manejadores = [Patch(facecolor=COLOR[e], alpha=0.5, edgecolor=COLOR[e],
                         label=NOMBRES[e][i]) for e in presentes]
    manejadores.append(Line2D([], [], color=TINTA_2, linewidth=1.0,
                              linestyle=(0, (3, 2)),
                              label=("media", "mean")[i]))
    fig.legend(handles=manejadores, frameon=False, fontsize=7,
               ncol=len(manejadores), loc="lower center",
               bbox_to_anchor=(0.5, 0.97), labelcolor=TINTA)
    fig.tight_layout()
    return _guardar(fig, destino)


# ---------------------------------------------------------------------------
# Figura 8 — tension de von Mises sobre la superficie, escala comun
# ---------------------------------------------------------------------------

CMAP_TENSION = "inferno"


def renders_von_mises(campos, carpeta, vista="iso", lado=LADO_GALERIA):
    """Render 3D de cada estructura coloreada por von Mises. Hilo de la
    interfaz (OpenGL).

    `campos`: {estructura: (campo en Pa con NaN fuera del hueso, spacing)} del
    analisis comparado. La escala de color es COMUN —percentiles 1 y 99 de
    todo el tejido de todas las estructuras— porque dos mapas con escalas
    propias se ven igual de «calientes» aunque uno soporte el doble de
    tension. Devuelve ({estructura: png}, (min, max) en MPa) o ({}, None).
    """
    from scipy import ndimage
    ests = [e for e in ORDEN if e in campos]
    if not ests:
        return {}, None
    vals = []
    for e in ests:
        c = np.asarray(campos[e][0], float)
        vals.append(c[np.isfinite(c) & (c > 0)] / 1e6)
    todos = np.concatenate([v for v in vals if v.size]) if vals else []
    if not len(todos):
        return {}, None
    clim = (float(np.percentile(todos, 1)), float(np.percentile(todos, 99)))
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    out = {}
    for e in ests:
        campo, sp = campos[e]
        campo = np.asarray(campo, float)
        sp = np.asarray(sp, float)
        solido = np.isfinite(campo)
        malla = malla_figura(solido, sp, suavizar=False)
        if not malla.n_points:
            continue
        # Cada vertice toma el valor del voxel SOLIDO mas proximo. No se usa
        # `espesor.muestrear_en_puntos`: su maximo 3x3x3 es inocuo para un
        # espesor, pero en una tension inflaria justo los picos.
        _d, (ix, iy, iz) = ndimage.distance_transform_edt(
            ~solido, return_indices=True)
        relleno = campo[ix, iy, iz] / 1e6
        idx = np.rint(np.asarray(malla.points, float) / sp).astype(int)
        for k in range(3):
            np.clip(idx[:, k], 0, campo.shape[k] - 1, out=idx[:, k])
        malla["vm"] = relleno[idx[:, 0], idx[:, 1], idx[:, 2]]
        lo, hi = _cubo(solido, sp)
        p = _plotter(lado)
        try:
            p.add_mesh(malla, scalars="vm", cmap=CMAP_TENSION, clim=clim,
                       show_scalar_bar=False, smooth_shading=True,
                       specular=0.15)
            try:
                p.enable_anti_aliasing("ssaa")
            except Exception:
                pass
            _camara(p, lo, hi, vista)
            ruta = carpeta / f"fig8_{e}.png"
            p.screenshot(str(ruta))
            out[e] = str(ruta)
        finally:
            p.close()
    return out, clim


def panel_von_mises(rutas, clim, destino, idioma="es", carga_N=100.0):
    """Renders de von Mises lado a lado con UNA barra de color comun."""
    from matplotlib import cm
    from matplotlib.colors import Normalize
    from matplotlib.image import imread
    i = _i(idioma)
    ests = [e for e in ORDEN if e in rutas]
    if not ests:
        return None
    k = len(ests)
    fig = _figura(2.3 * k + 0.9, 2.6)
    for j, est in enumerate(ests):
        ax = fig.add_axes([0.01 + j * (0.86 / k), 0.02, 0.86 / k - 0.01,
                           0.86])
        ax.imshow(imread(rutas[est]))
        ax.set_axis_off()
        ax.set_title(f"({'abc'[j]}) {NOMBRES[est][i]}", fontsize=8.5,
                     color=TINTA, loc="left")
    cax = fig.add_axes([0.9, 0.14, 0.022, 0.64])
    barra = fig.colorbar(cm.ScalarMappable(norm=Normalize(*clim),
                                           cmap=CMAP_TENSION), cax=cax)
    barra.outline.set_edgecolor(EJE)
    barra.ax.tick_params(labelsize=6.5, colors=TINTA_2)
    barra.set_label(("σ von Mises [MPa]", "von Mises σ [MPa]")[i]
                    + f", {carga_N:g} N", fontsize=7, color=TINTA_2)
    return _guardar(fig, destino, pdf=False)[0]
