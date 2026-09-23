"""
figura_metodo.py — Figura 0 del informe: el metodo de principio a fin.

QUE DIBUJA
----------
Una fila por etapa del trabajo, en formato 16:9 para que sirva igual en un
articulo (doble columna) que en una diapositiva:

  VOI           rebanada de micro-CT -> la misma rebanada segmentada -> pila
                3D con el cubo recortado en su marco PCA -> VOI -> morfometria
                que persigue el ajuste
  spinodoide    direcciones de onda en sus conos -> GRF en las caras del cubo
                -> histograma del campo y umbral phi0 -> solido -> candidato
                frente al VOI
  dual-lattice  red dual 4-N -> distancia al esqueleto -> umbral r como
                cuantil -> solido -> candidato frente al VOI

Solo aparecen las filas de lo que hay: una familia, dos, o ninguna fila de
VOI si la sesion no tiene VOI.

NADA ES UN DIBUJO
-----------------
Cada panel sale del dato que describe:

  * la rebanada es el archivo del escaner que pasa por el CENTRO del cubo, y
    la segmentada es esa misma rebanada de la mascara que entro en la pila,
    asi que las dos coinciden pixel a pixel; el contorno es la interseccion
    exacta del cubo con ese plano;
  * el cubo se dibuja con las esquinas que salen de la MISMA transformacion
    que `voi.extraer_cubo` (marco PCA -> mm), no con una caja aproximada;
  * las direcciones de onda se recalculan con la semilla del ajuste y el mismo
    orden de consumo del generador que `grf.campo_grf`: son las que genero la
    estructura (el bloque 24 lo comprueba bit a bit);
  * los campos y los umbrales son los que devuelve el generador, y el solido
    es exactamente {campo <= umbral} (tambien comprobado);
  * las tablas son la morfometria de la sesion, la misma de la seccion 2.

EL CONTEXTO DEL VOI
-------------------
Un VOI recortado de una pila ya no lleva la pila: el visor la descarta para
liberar memoria. Por eso el contexto (pila reducida, rebanada, esquinas) se
toma EN EL MOMENTO DEL RECORTE con `contexto_pila`, y el documento de la
sesion guarda `registro_recorte`, con lo necesario para rehacerlo desde el
disco (`rehacer_contexto`) cuando el informe se reconstruye con la CLI. Si el
VOI se cargo de un .vtk, no hay pila: la fila del VOI muestra su corte central
y lo dice.

HILOS
-----
`renders_metodo` usa VTK y va en el hilo de la interfaz, como `renders_3d`.
Todo lo demas (contexto, datos de familia, composicion) es numpy/matplotlib.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import figuras as F

# Lado maximo, en celdas, de la pila reducida que se dibuja en la fila del
# VOI. Es contexto (donde esta el cubo dentro del hueso), no medida: 160 da
# una silueta nitida y cabe en unos MB.
LADO_REDUCIDO_MAX = 160
# Nivel de la isosuperficie de la pila reducida (media por bloques): por
# debajo de 0,5 para que las trabeculas finas no desaparezcan al promediar.
NIVEL_REDUCIDO = 0.3
# Subcubo central de la red dual que se dibuja (fraccion del lado): la red
# entera a 96^3 es una maraña que no se lee.
SUBCUBO_ESQUELETO = (0.2, 0.8)
LADO_RENDER = 1200
COLOR_VOI_3D = F.ESTILOS_3D["hueso"]["color"]


def _t(es, en, idioma):
    return es if idioma == "es" else en


# ---------------------------------------------------------------------------
# Contexto del VOI (pila, rebanada, cubo)
# ---------------------------------------------------------------------------

def _lado_vox(lado_mm, h):
    """El mismo lado en voxeles que `voi.extraer_cubo` (impar)."""
    n = int(round(float(lado_mm) / float(h)))
    return n + 1 if n % 2 == 0 else n


def esquinas_cubo(ejes, centro, centro_pca, lado_mm, spacing):
    """Esquinas (8, 3) en mm, en el sistema de la pila, del cubo que recorta
    `voi.extraer_cubo` con esos argumentos.

    `extraer_cubo` muestrea centros de voxel a (i - h_v) * h del centro, con
    i = 0..lado_vox-1; la caja que ocupan esos voxeles llega medio voxel mas
    alla del ultimo centro. Mismo mapeo marco PCA -> mm: (delta + c) @ ejes.T
    + centro.
    """
    sp = np.atleast_1d(np.asarray(spacing, float)).ravel()
    sp = np.repeat(sp, 3) if sp.size == 1 else sp[:3]
    n = _lado_vox(lado_mm, sp[0])
    semi = ((n - 1) / 2.0 + 0.5) * sp
    s = np.array([[a, b, c] for a in (-1, 1) for b in (-1, 1) for c in (-1, 1)],
                 float) * semi
    return (s + np.asarray(centro_pca, float)) @ np.asarray(ejes, float).T \
        + np.asarray(centro, float)


def corte_cubo_plano_z(esquinas, z0):
    """Poligono (m, 2) de la interseccion del cubo con el plano z = z0, en mm
    (x, y), ordenado por angulo. Vacio si el plano no corta el cubo."""
    E = np.asarray(esquinas, float)
    pts = []
    for a in range(8):
        for b in range(a + 1, 8):
            if bin(a ^ b).count("1") != 1:          # solo aristas
                continue
            za, zb = E[a, 2] - z0, E[b, 2] - z0
            if za == 0.0:
                pts.append(E[a, :2])
            if za * zb < 0.0:
                t = za / (za - zb)
                pts.append(E[a, :2] + t * (E[b, :2] - E[a, :2]))
    if len(pts) < 3:
        return np.zeros((0, 2))
    P = np.unique(np.round(np.array(pts), 12), axis=0)
    c = P.mean(axis=0)
    return P[np.argsort(np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0]))]


def reducir(mask, lado_max=LADO_REDUCIDO_MAX):
    """Media por bloques p^3 (p entero) para que el lado mayor quede <= lado_max.
    Devuelve (volumen float32, p)."""
    mask = np.asarray(mask)
    p = max(1, int(np.ceil(max(mask.shape) / float(lado_max))))
    n = (np.array(mask.shape) // p) * p
    m = mask[:n[0], :n[1], :n[2]].reshape(n[0] // p, p, n[1] // p, p,
                                          n[2] // p, p)
    return m.mean(axis=(1, 3, 5), dtype=np.float32), p


def _rebanada_gris(origen, patron, k, n_rebanadas):
    """Rebanada k del archivo original, (x, y), o None si no se puede leer o
    la serie del disco no es la que se cargo."""
    from . import voi as V
    try:
        p = Path(origen)
        if p.is_file():
            from PIL import Image
            im = Image.open(str(p))
            im.seek(int(k))
            a = np.asarray(im)
        else:
            rutas, _f, _ap = V._separar_rebanadas(V.listar_rebanadas(p, patron))
            if len(rutas) != int(n_rebanadas):
                return None
            a = V._leer_imagen(rutas[int(k)])
        a = np.asarray(a, float)
        if a.ndim == 3:
            a = a[..., 0]
        return a.T
    except Exception:
        return None


def contexto_pila(mask, spacing, ejes, centro, centro_pca, lado_mm,
                  origen=None, patron="*.tif", reducido=None):
    """Lo que la fila del VOI necesita de la pila, tomado antes de soltarla.

    `reducido`: (volumen, p) de `reducir` si ya se calculo (el visor lo hace
    en el hilo de trabajo antes de abrir el dialogo de recorte).
    """
    mask = np.asarray(mask, bool)
    sp = np.atleast_1d(np.asarray(spacing, float)).ravel()
    sp = np.repeat(sp, 3) if sp.size == 1 else sp[:3]
    esq = esquinas_cubo(ejes, centro, centro_pca, lado_mm, sp)
    c_mm = np.asarray(centro_pca, float) @ np.asarray(ejes, float).T \
        + np.asarray(centro, float)
    k = int(np.clip(np.rint(c_mm[2] / sp[2]), 0, mask.shape[2] - 1))
    vol, p = reducido if reducido is not None else reducir(mask)
    gris = (_rebanada_gris(origen, patron, k, mask.shape[2])
            if origen else None)
    if gris is not None and gris.shape != mask.shape[:2]:
        gris = None
    return {"hueso": vol, "paso": int(p), "spacing": sp, "esquinas": esq,
            "k": k, "n_rebanadas": int(mask.shape[2]),
            "binaria": mask[:, :, k].copy(), "gris": gris,
            "poligono": corte_cubo_plano_z(esq, k * sp[2]),
            "lado_mm": float(lado_mm),
            "lado_vox": _lado_vox(lado_mm, sp[0]),
            "origen": str(origen) if origen else None}


def registro_recorte(pila, ejes, centro, centro_pca, lado_mm):
    """Lo que el documento guarda para rehacer el contexto desde el disco.

    `pila`: {"origen", "patron", "umbral", "tam_voxel_mm"} de la carga.
    """
    return {"origen": str(pila.get("origen")), "patron": pila.get("patron"),
            "umbral": pila.get("umbral"),
            "tam_voxel_mm": pila.get("tam_voxel_mm"),
            "ejes": np.asarray(ejes, float).tolist(),
            "centro": np.asarray(centro, float).tolist(),
            "centro_pca": np.asarray(centro_pca, float).tolist(),
            "lado_mm": float(lado_mm)}


def rehacer_contexto(rec):
    """Contexto desde `registro_recorte`, releyendo la pila (lento: la pila
    entera). None si el origen ya no esta o no se deja leer."""
    if not isinstance(rec, dict) or not rec.get("origen"):
        return None
    try:
        from .voi import leer_pila_tiff
        BW, sp, _info = leer_pila_tiff(rec["origen"],
                                       patron=rec.get("patron") or "*.tif",
                                       umbral=rec.get("umbral"),
                                       tam_voxel=rec.get("tam_voxel_mm"),
                                       unidad="mm")
        return contexto_pila(BW, sp, rec["ejes"], rec["centro"],
                             rec["centro_pca"], rec["lado_mm"],
                             origen=rec["origen"],
                             patron=rec.get("patron") or "*.tif")
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Datos de cada familia
# ---------------------------------------------------------------------------

def datos_familia(par, BW=None, campo=None, info=None):
    """Mascara, campo, umbral y la aleatoriedad visible de un ajuste.

    Con `BW/campo/info` ya regenerados (el informe los tiene de `_paquete`)
    no se vuelve a llamar al generador.
    """
    from . import procedencia
    kw = procedencia.kwargs_generador(par)
    fam = par.get("familia", "spinodoide")
    if fam == "dual-lattice":
        from . import dual_lattice as dl
        if BW is None:
            BW, campo, info = dl.generar_dual_lattice(**kw)
        esq = dl.esqueleto(kw["celdas"], kw["estiramiento"],
                           kw["irregularidad"], R=kw["R"], seed=kw["seed"])
        return {"familia": fam, "BW": BW, "campo": campo,
                "umbral": float(info["radio_puntal"]), "esqueleto": esq,
                "par": par, "info": info}
    from . import grf
    if BW is None:
        BW, campo, info = grf.generar_mascara(**kw)
    # Mismo generador y mismo orden de consumo que `campo_grf`: primero las
    # direcciones. Son las de la estructura, no unas parecidas.
    rng = np.random.default_rng(kw["seed"])
    dirs = grf.wave_directions(kw["num_waves"], kw["thetas"], R=kw["R"],
                               esquema=kw["esquema"], rng=rng)
    return {"familia": fam, "BW": BW, "campo": campo,
            "umbral": float(info["levelset"]), "dirs": dirs, "par": par,
            "info": info}


# ---------------------------------------------------------------------------
# Renders (VTK: hilo de la interfaz)
# ---------------------------------------------------------------------------

def _recortar_blanco(a, margen=8):
    a = np.asarray(a)
    m = np.where((a[..., :3] < 245).any(axis=2))
    if not m[0].size:
        return a
    y0, y1, x0, x1 = m[0].min(), m[0].max(), m[1].min(), m[1].max()
    return a[max(0, y0 - margen):y1 + margen + 1,
             max(0, x0 - margen):x1 + margen + 1]


def _superficie(m):
    """`extract_surface` sin el aviso de pyvista reciente, y sin romper las
    versiones que aun no aceptan `algorithm`."""
    try:
        return m.extract_surface(algorithm="dataset_surface")
    except TypeError:
        return m.extract_surface()


def _foto(p):
    p.render()
    a = p.screenshot(return_img=True)
    p.close()
    return _recortar_blanco(a)


def _caja(p, lo, hi, color="#9a9994", ancho=2):
    import pyvista as pv
    p.add_mesh(pv.Box(bounds=(lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]))
               .extract_all_edges(), color=color, line_width=ancho)


def _aa(p):
    try:
        p.enable_anti_aliasing("ssaa")
    except Exception:
        pass


def render_solido(BW, spacing, color, suavizar=True, lado=LADO_RENDER):
    sp = np.atleast_1d(np.asarray(spacing, float)).ravel()
    sp = np.repeat(sp, 3) if sp.size == 1 else sp[:3]
    malla = F.malla_figura(BW, sp, suavizar=suavizar)
    lo, hi = F._cubo(BW, sp)
    p = F._plotter(lado)
    p.add_mesh(malla, color=color, smooth_shading=True, specular=0.25,
               specular_power=15)
    try:
        d = float(np.linalg.norm(hi - lo))
        p.enable_ssao(radius=0.045 * d, bias=0.0009 * d, kernel_size=128,
                      blur=True)
    except Exception:
        pass
    _aa(p)
    _caja(p, lo, hi, color="#c3c2b7", ancho=1)
    F._camara(p, lo, hi, "iso")
    return _foto(p)


def render_hueso(ctx, lado=LADO_RENDER):
    import pyvista as pv
    vol, paso, sp = ctx["hueso"], ctx["paso"], ctx["spacing"]
    h = sp * paso
    lleno = np.zeros(np.array(vol.shape) + 2, np.float32)
    lleno[1:-1, 1:-1, 1:-1] = vol
    # Centro del bloque p^3 = (p - 1)/2 voxeles mas alla de su esquina.
    org = -h + 0.5 * (paso - 1) * sp
    g = pv.ImageData(dimensions=lleno.shape, spacing=tuple(h), origin=tuple(org))
    g.point_data["v"] = lleno.ravel(order="F")
    malla = g.contour([NIVEL_REDUCIDO], scalars="v")
    try:
        malla = malla.smooth_taubin(n_iter=30, pass_band=0.05)
    except Exception:
        pass
    malla = malla.compute_normals(auto_orient_normals=True,
                                  consistent_normals=False)
    p = F._plotter(lado)
    p.add_mesh(malla, color=COLOR_VOI_3D, smooth_shading=True, opacity=0.35,
               specular=0.2)
    E = ctx["esquinas"]
    for a in range(8):
        for b in range(a + 1, 8):
            if bin(a ^ b).count("1") == 1:
                p.add_mesh(pv.Line(E[a], E[b]), color=F.COLOR["voi"],
                           line_width=5)
    p.add_mesh(_superficie(pv.PolyData(E).delaunay_3d()),
               color=F.COLOR["voi"], opacity=0.35)
    _aa(p)
    lo, hi = np.array(malla.bounds[::2]), np.array(malla.bounds[1::2])
    c = 0.5 * (lo + hi)
    diag = float(np.linalg.norm(hi - lo))
    d = np.array([1.0, -0.6, 0.55])
    d /= np.linalg.norm(d)
    p.camera.focal_point = tuple(c)
    p.camera.position = tuple(c + 3 * diag * d)
    p.camera.up = (0, 0, 1)
    p.camera.parallel_projection = True
    p.camera.parallel_scale = 0.44 * diag
    p.renderer.ResetCameraClippingRange()
    return _foto(p)


def render_campo(campo, lado_mm, cmap, clim, umbral, color_iso, lado=LADO_RENDER):
    """Campo en las tres caras visibles del cubo desde la vista isometrica, con
    la isolinea del umbral: es la frontera del solido en esas caras."""
    import pyvista as pv
    n = campo.shape[0]
    h = float(lado_mm) / n
    g = pv.ImageData(dimensions=campo.shape, spacing=(h, h, h))
    g.point_data["c"] = np.asarray(campo, np.float32).ravel(order="F")
    sup = _superficie(g)
    p = F._plotter(lado)
    p.add_mesh(sup, scalars="c", cmap=cmap, clim=clim, show_scalar_bar=False,
               lighting=False)
    iso = sup.contour([float(umbral)], scalars="c")
    if iso.n_points:
        p.add_mesh(iso, color=color_iso, line_width=2.5)
    _aa(p)
    lo, hi = np.zeros(3), np.full(3, (n - 1) * h)
    _caja(p, lo, hi, color="#c3c2b7", ancho=1)
    F._camara(p, lo, hi, "iso")
    return _foto(p)


def render_conos(thetas, R, dirs, color, lado=LADO_RENDER):
    """Esfera unidad, conos de semiangulo theta_i sobre los ejes de R (las
    columnas), y las direcciones de onda (+n y -n: cos es par)."""
    import pyvista as pv
    R = np.asarray(R, float)
    p = F._plotter(lado)
    p.add_mesh(pv.Sphere(radius=1.0, theta_resolution=80, phi_resolution=80),
               color="#f1f0ec", opacity=0.25, smooth_shading=True)
    tonos = ("#c9542a", "#e0874f", "#f0b88c")
    for i, t in enumerate(np.asarray(thetas, float)):
        eje = R[:, i]
        p.add_mesh(pv.Line(-1.3 * eje, 1.3 * eje), color="#777777",
                   line_width=2)
        if t <= 0:
            continue
        alto, radio = np.cos(np.radians(t)), np.sin(np.radians(t))
        for s in (1, -1):
            p.add_mesh(pv.Cone(center=s * eje * alto / 2, direction=-s * eje,
                               height=alto, radius=radio, resolution=80,
                               capping=False),
                       color=tonos[i], opacity=0.22, smooth_shading=True)
    pts = np.vstack([dirs, -dirs])
    p.add_mesh(pv.PolyData(pts).glyph(geom=pv.Sphere(radius=0.022),
                                      scale=False, orient=False),
               color=color, smooth_shading=True)
    _aa(p)
    p.camera.focal_point = (0, 0, 0)
    p.camera.position = (3.2, -4.2, 2.6)
    p.camera.up = (0, 0, 1)
    p.camera.parallel_projection = True
    p.camera.parallel_scale = 1.45
    return _foto(p)


def render_esqueleto(esq, lado_mm, color, lado=LADO_RENDER):
    import pyvista as pv
    lo, hi = SUBCUBO_ESQUELETO
    A, B = esq["seg_a"], esq["seg_b"]
    ok = np.all((np.minimum(A, B) <= hi + 0.3) & (np.maximum(A, B) >= lo - 0.3),
                axis=1)
    A, B = A[ok] * lado_mm, B[ok] * lado_mm
    n = len(A)
    lineas = np.hstack([np.full((n, 1), 2), np.arange(n)[:, None],
                        n + np.arange(n)[:, None]]).ravel()
    red = pv.PolyData(np.vstack([A, B]), lines=lineas)
    L0, L1 = lo * lado_mm, hi * lado_mm
    caja = (L0, L1, L0, L1, L0, L1)
    tubo = red.tube(radius=0.005 * lado_mm).clip_box(caja, invert=False)
    p = F._plotter(lado)
    p.add_mesh(_superficie(tubo), color=color, smooth_shading=True)
    nod = esq["nodos"] * lado_mm
    nod = nod[np.all((nod >= L0) & (nod <= L1), axis=1)]
    if len(nod):
        p.add_mesh(pv.PolyData(nod).glyph(
            geom=pv.Sphere(radius=0.011 * lado_mm), scale=False, orient=False),
            color="#222222", smooth_shading=True)
    _aa(p)
    lo3, hi3 = np.full(3, L0), np.full(3, L1)
    _caja(p, lo3, hi3, color="#c3c2b7", ancho=1)
    F._camara(p, lo3, hi3, "iso")
    return _foto(p)


def renders_metodo(voi, familias, lado_mm, ctx=None, suavizar=True,
                   lado=LADO_RENDER):
    """Todas las imagenes 3D de la figura. Hilo de la interfaz (OpenGL).

    voi: (BW, spacing) o None. familias: {familia: datos_familia(...)}.
    Devuelve {clave: imagen RGB}; lo que falla se omite y se anota en
    "_fallos" (la figura se compone con lo que haya).
    """
    out, fallos = {}, []

    def hacer(clave, fn, *a, **k):
        try:
            out[clave] = fn(*a, **k)
        except Exception as e:                       # noqa: BLE001
            fallos.append(f"{clave}: {type(e).__name__}: {e}")

    if ctx is not None:
        hacer("hueso", render_hueso, ctx, lado)
    if voi is not None:
        hacer("voi", render_solido, voi[0], voi[1], COLOR_VOI_3D, suavizar,
              lado)
    for fam, d in familias.items():
        col = F.COLOR[fam]
        n = d["BW"].shape[0]
        sp = np.full(3, float(lado_mm) / n)
        if fam == "dual-lattice":
            hacer(fam + ":esqueleto", render_esqueleto, d["esqueleto"],
                  lado_mm, col, lado)
            D = np.asarray(d["campo"], float) * lado_mm      # a mm
            hacer(fam + ":campo", render_campo, D, lado_mm, "viridis_r",
                  (0.0, float(np.percentile(D, 99))),
                  d["umbral"] * lado_mm, "#ffffff", lado)
        else:
            par = d["par"]
            hacer(fam + ":conos", render_conos, par["thetas"], par["R"],
                  d["dirs"], col, lado)
            G = np.asarray(d["campo"], float)
            v = float(np.percentile(np.abs(G), 99.5))
            hacer(fam + ":campo", render_campo, G, lado_mm, "RdBu_r", (-v, v),
                  d["umbral"], "#222222", lado)
        hacer(fam + ":solido", render_solido, d["BW"], sp, col, suavizar, lado)
    out["_fallos"] = fallos
    return out


# ---------------------------------------------------------------------------
# Composicion (matplotlib)
# ---------------------------------------------------------------------------

METRICAS = (("BVTV", "BV/TV", "", 3), ("TbTh", "Tb.Th", "mm", 3),
            ("TbSp", "Tb.Sp", "mm", 3), ("BSBV", "BS/BV", "mm⁻¹", 1),
            ("DA", "DA", "", 2))
FONDO_TARJETA = "#f4f3ef"
ANCHO, ALTO_FILA = 13.33, 2.45


def _num(x, d, idioma):
    s = f"{x:.{d}f}"
    return s.replace(".", ",") if idioma == "es" else s


def _v(morf, k):
    try:
        x = float(np.asarray((morf or {}).get(k), float).reshape(()))
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def _mostrar(ax, a):
    ax.imshow(a, interpolation="lanczos")
    ax.set_aspect("equal")
    ax.axis("off")


def _titulo(ax, letra, texto, fs):
    ax.text(0.0, 1.02, letra, transform=ax.transAxes, fontsize=fs + 1,
            weight="bold", color=F.TINTA, va="bottom", ha="left")
    ax.annotate(texto, xy=(0.0, 1.02), xycoords="axes fraction",
                xytext=(fs * 1.25, 0), textcoords="offset points",
                fontsize=fs, color=F.TINTA, va="bottom", ha="left")


def _flecha(fig, p0, p1, color=F.EJE, lw=1.4, ms=12):
    from matplotlib.patches import FancyArrowPatch
    fig.patches.append(FancyArrowPatch(p0, p1, transform=fig.transFigure,
                                       arrowstyle="-|>", mutation_scale=ms,
                                       color=color, lw=lw))


def _flecha_ejes(fig, a1, a2):
    b1, b2 = a1.get_position(), a2.get_position()
    y = 0.5 * (b1.y0 + b1.y1)
    _flecha(fig, (b1.x1 + 0.004, y), (b2.x0 - 0.004, y))


def _banda(fig, y0, y1, texto, color):
    from matplotlib.patches import FancyBboxPatch
    fig.patches.append(FancyBboxPatch(
        (0.008, y0), 0.022, y1 - y0,
        boxstyle="round,pad=0,rounding_size=0.006",
        transform=fig.transFigure, fc=color, ec="none"))
    fig.text(0.019, 0.5 * (y0 + y1), texto, rotation=90, ha="center",
             va="center", color="white", fontsize=11, weight="bold")


def _fondo_tarjeta(ax):
    from matplotlib.patches import FancyBboxPatch
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.02, 0.04), 0.96, 0.92,
                                boxstyle="round,pad=0,rounding_size=0.05",
                                transform=ax.transAxes, fc=FONDO_TARJETA,
                                ec="none"))


def _tarjeta_voi(ax, morf, idioma, fs):
    _fondo_tarjeta(ax)
    ax.text(0.08, 0.86, _t("Objetivo del ajuste", "Fitting target", idioma),
            fontsize=fs, weight="bold", color=F.COLOR["voi"],
            transform=ax.transAxes)
    y = 0.72
    for k, nom, u, d in METRICAS:
        x = _v(morf, k)
        ax.text(0.08, y, nom, fontsize=fs - 0.5, color=F.TINTA_2,
                transform=ax.transAxes)
        ax.text(0.92, y, "—" if x is None else f"{_num(x, d, idioma)} {u}".strip(),
                fontsize=fs - 0.5, color=F.TINTA, ha="right",
                transform=ax.transAxes)
        y -= 0.13
    ax.text(0.08, 0.06, _t("misma ruta de medida\npara VOI y candidatos",
                           "same measurement path\nfor VOI and candidates",
                           idioma),
            fontsize=fs - 2.5, color=F.TINTA_2, transform=ax.transAxes,
            va="bottom")


def _tabla_candidato(ax, morf, morf_voi, color, idioma, fs):
    _fondo_tarjeta(ax)
    ax.text(0.08, 0.86, _t("Candidato vs VOI", "Candidate vs VOI", idioma)
            if morf_voi else _t("Candidato", "Candidate", idioma),
            fontsize=fs, weight="bold", color=color, transform=ax.transAxes)
    y = 0.72
    for k, nom, _u, d in METRICAS:
        x, t = _v(morf, k), _v(morf_voi, k)
        ax.text(0.08, y, nom, fontsize=fs - 0.5, color=F.TINTA_2,
                transform=ax.transAxes)
        ax.text(0.62, y, "—" if x is None else _num(x, d, idioma),
                fontsize=fs - 0.5, color=F.TINTA, ha="right",
                transform=ax.transAxes)
        if x is not None and t:
            dd = f"{100.0 * (x - t) / t:+.1f} %"
            ax.text(0.94, y, dd.replace(".", ",") if idioma == "es" else dd,
                    fontsize=fs - 1.5, color=F.TINTA_2, ha="right",
                    transform=ax.transAxes)
        y -= 0.13


def _histograma(ax, valores, umbral, color, xlabel, texto, fs, gauss=False):
    F._estilo(ax, rejilla="y")
    h, b = np.histogram(np.asarray(valores, float).ravel(), bins=70,
                        density=True)
    c = 0.5 * (b[1:] + b[:-1])
    w = b[1] - b[0]
    iz = c <= umbral
    ax.bar(c[iz], h[iz], width=w, color=color, lw=0)
    ax.bar(c[~iz], h[~iz], width=w, color=F.REJILLA, lw=0)
    if gauss:
        x = np.linspace(b[0], b[-1], 300)
        ax.plot(x, np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi), color=F.TINTA,
                lw=0.9, ls=(0, (3, 2)))
    ax.axvline(umbral, color=F.TINTA, lw=1.0)
    ax.text(0.98, 0.96, texto, transform=ax.transAxes, fontsize=fs - 1,
            ha="right", va="top", color=F.TINTA, linespacing=1.35)
    ax.set_ylim(0, h.max() * 1.45)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.tick_params(labelsize=fs - 1.5)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=fs - 0.5, color=F.TINTA_2)


def _recuadro(img2d, poligono_px, margen=0.06):
    """Ventana (y0, y1, x0, x1) que contiene el hueso de la rebanada y el
    poligono del cubo, con margen."""
    ys, xs = np.nonzero(img2d)
    if len(poligono_px):
        xs = np.concatenate([xs, poligono_px[:, 0]])
        ys = np.concatenate([ys, poligono_px[:, 1]])
    if not len(xs):
        return 0, img2d.shape[0], 0, img2d.shape[1]
    lado = max(np.ptp(xs), np.ptp(ys)) * (1 + 2 * margen)
    cx, cy = 0.5 * (xs.min() + xs.max()), 0.5 * (ys.min() + ys.max())
    y0 = int(max(0, cy - lado / 2))
    x0 = int(max(0, cx - lado / 2))
    return (y0, int(min(img2d.shape[0], cy + lado / 2)), x0,
            int(min(img2d.shape[1], cx + lado / 2)))


def _rebanada(ax, ctx, gris, fs, barra=True):
    """Rebanada (gris o segmentada) recortada al hueso, con el corte del cubo
    en azul y 1 mm de barra. Pixel = (fila y, columna x), como el archivo."""
    from matplotlib.patches import Polygon, Rectangle
    sp = ctx["spacing"]
    B = ctx["binaria"].T                         # (y, x)
    pol = ctx["poligono"] / sp[:2] if len(ctx["poligono"]) else np.zeros((0, 2))
    y0, y1, x0, x1 = _recuadro(B, pol)
    if gris:
        g = ctx["gris"].T
        lo, hi = np.percentile(g, (0.5, 99.7))
        ax.imshow(np.clip((g - lo) / max(hi - lo, 1e-12), 0, 1), cmap="gray",
                  interpolation="lanczos")
    else:
        ax.imshow(B, cmap="gray_r", interpolation="nearest", vmin=0, vmax=1)
    if len(pol):
        ax.add_patch(Polygon(pol, closed=True, fill=False,
                             ec=F.COLOR["voi"], lw=1.4))
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    if barra:
        px = 1.0 / sp[0]
        ancho = x1 - x0
        xb, yb = x1 - 0.14 * ancho - px, y1 - 0.06 * ancho
        col = "white" if gris else F.TINTA
        ax.add_patch(Rectangle((xb, yb), px, 0.012 * ancho, color=col))
        ax.text(xb + px / 2, yb - 0.01 * ancho, "1 mm", color=col,
                ha="center", va="bottom", fontsize=fs - 2)
    ax.axis("off")


def componer_metodo(imagenes, familias, doc, destino, idioma="es", ctx=None,
                    voi=None, lado_mm=None):
    """Compone y guarda la figura 0. Devuelve [png, pdf] o [] si no hay
    nada que dibujar.

    imagenes: `renders_metodo`. familias: {familia: datos_familia}. doc: el
    documento de la sesion (morfometria). voi: (BW, spacing) o None.
    """
    fams = [f for f in ("spinodoide", "dual-lattice") if f in familias]
    hay_voi = voi is not None and "voi" in imagenes
    if not fams and not hay_voi:
        return []
    fs = 10
    filas = ([("voi", None)] if hay_voi else []) + [(f, None) for f in fams]
    nf = len(filas)
    alto = ALTO_FILA * nf + 0.15
    fig = F._figura(ANCHO, alto)
    xs = [0.045, 0.235, 0.425, 0.615, 0.805]
    w = 0.17
    # Posiciones en pulgadas, pasadas a fraccion: cada fila mide lo mismo sea
    # cual sea el numero de filas.
    arriba_in, hueco_in, alto_in = 0.30, 0.52, ALTO_FILA - 0.52
    letras = iter("abcdefghijklmnopqrstuvwxyz")
    mv = doc.get("morfometria_voi") if doc else None
    nombres = {"voi": ("VOI", "VOI"), "spinodoide": ("Spinodoide", "Spinodoid"),
               "dual-lattice": ("Dual-lattice", "Dual-lattice")}
    li = 0 if idioma == "es" else 1
    tarjetas = {}
    for i, (fila, _x) in enumerate(filas):
        y1 = 1.0 - (arriba_in + i * ALTO_FILA) / alto
        y0 = y1 - alto_in / alto
        ax = [fig.add_axes([x, y0, w, y1 - y0]) for x in xs]
        color = F.COLOR[fila]
        _banda(fig, y0, y1, nombres[fila][li], color)
        cadena = []
        if fila == "voi":
            if ctx is not None:
                hay_gris = ctx.get("gris") is not None
                _rebanada(ax[0], ctx, hay_gris, fs)
                h_um = ctx["spacing"][0] * 1000.0
                _titulo(ax[0], next(letras), _t(
                    f"micro-CT, {_num(h_um, 1, idioma)} µm/vóxel",
                    f"micro-CT, {_num(h_um, 1, idioma)} µm/voxel", idioma)
                    if hay_gris else _t("rebanada de la pila", "stack slice",
                                        idioma), fs)
                _rebanada(ax[1], ctx, False, fs, barra=False)
                _titulo(ax[1], next(letras), _t("segmentación", "segmentation",
                                                idioma), fs)
                if "hueso" in imagenes:
                    _mostrar(ax[2], imagenes["hueso"])
                else:
                    ax[2].axis("off")
                _titulo(ax[2], next(letras), _t("pila 3D, marco PCA, VOI",
                                                "3D stack, PCA frame, VOI",
                                                idioma), fs)
                cadena = [0, 1, 2, 3]
            else:
                for a in ax[:2]:
                    a.axis("off")
                ax[0].text(0.0, 0.5, _t(
                    "VOI leído de un archivo ya\nrecortado: la sesión no tiene\n"
                    "la pila de micro-CT de la\nque salió.",
                    "VOI read from an already\ncropped file: the session does\n"
                    "not hold the micro-CT stack\nit came from.", idioma),
                    transform=ax[0].transAxes, fontsize=fs - 1,
                    color=F.TINTA_2, va="center", linespacing=1.4)
                BWv = np.asarray(voi[0], bool)
                ax[2].imshow(BWv[:, BWv.shape[1] // 2, :].T, cmap="gray_r",
                             origin="lower", interpolation="nearest",
                             vmin=0, vmax=1)
                ax[2].set_xticks([])
                ax[2].set_yticks([])
                for s in ax[2].spines.values():
                    s.set_color(F.EJE)
                _titulo(ax[2], next(letras), _t("corte central XZ del VOI",
                                                "central XZ section of the VOI",
                                                idioma), fs)
                cadena = [2, 3]
            _mostrar(ax[3], imagenes["voi"])
            lado_txt = _num(float(lado_mm), 2, idioma) if lado_mm else "?"
            _titulo(ax[3], next(letras), _t(f"VOI cúbico {lado_txt} mm",
                                            f"cubic VOI {lado_txt} mm",
                                            idioma), fs)
            _tarjeta_voi(ax[4], mv, idioma, fs)
            tarjetas["voi"] = ax[4]
        else:
            d = familias[fila]
            if fila == "spinodoide":
                k1, k2 = fila + ":conos", fila + ":campo"
                t1 = _t(r"ondas en conos $\theta_i$",
                        r"waves in cones $\theta_i$", idioma)
                t2 = r"GRF $\phi(\mathbf{x})$"
                _histograma(ax[2], d["campo"], d["umbral"], F.COLOR[fila], "",
                            r"$\phi_0=\sqrt{2}\,\mathrm{erf}^{-1}(2\rho-1)$",
                            fs, gauss=True)
                t4 = r"$\{\phi\leq\phi_0\}$"
            else:
                k1, k2 = fila + ":esqueleto", fila + ":campo"
                t1 = _t("red dual 4-N", "4-N dual lattice", idioma)
                t2 = _t(r"distancia $d(\mathbf{x})$",
                        r"distance $d(\mathbf{x})$", idioma)
                L = float(lado_mm) if lado_mm else 1.0
                _histograma(ax[2], np.asarray(d["campo"], float) * L,
                            d["umbral"] * L, F.COLOR[fila],
                            r"$d$ (mm)" if lado_mm else r"$d$",
                            _t(r"$r$ = cuantil $\rho$ de $d$",
                               r"$r$ = $\rho$-quantile of $d$", idioma), fs)
                t4 = r"$\{d\leq r\}$"
            for j, clave in ((0, k1), (1, k2), (3, fila + ":solido")):
                if clave in imagenes:
                    _mostrar(ax[j], imagenes[clave])
                else:
                    ax[j].axis("off")
            for j, t in ((0, t1), (1, t2), (2, _t("umbral por densidad",
                                                   "density threshold",
                                                   idioma)), (3, t4)):
                _titulo(ax[j], next(letras), t, fs)
            morf = doc.get({"spinodoide": "morfometria_spin",
                            "dual-lattice": "morfometria_dual"}[fila]) \
                if doc else None
            _tabla_candidato(ax[4], morf, mv, F.COLOR[fila], idioma, fs)
            tarjetas[fila] = ax[4]
            cadena = [0, 1, 2, 3]
        for a, b in zip(cadena[:-1], cadena[1:]):
            _flecha_ejes(fig, ax[a], ax[b])
    if "voi" in tarjetas and fams:
        b0 = tarjetas["voi"].get_position()
        b1 = tarjetas[fams[0]].get_position()
        x = 0.5 * (b0.x0 + b0.x1)
        _flecha(fig, (x, b0.y0 - 0.005), (x, b1.y1 + 0.005), color=F.TINTA,
                ms=10)
        fig.text(x + 0.006, 0.5 * (b0.y0 + b1.y1),
                 _t("ajuste", "fit", idioma), fontsize=fs - 1.5,
                 color=F.TINTA_2, va="center")
    return F._guardar(fig, destino)


def figura_metodo(voi, familias, doc, carpeta, lado_mm, ctx=None,
                  suavizar=True, lado=LADO_RENDER):
    """Renders + las dos versiones (ES/EN). Hilo de la interfaz.

    Devuelve ({"es": [png, pdf], "en": [...]}, fallos de render).
    """
    img = renders_metodo(voi, familias, lado_mm, ctx, suavizar, lado)
    fallos = img.pop("_fallos", [])
    out = {}
    for idioma, nombre in (("es", "fig0_metodo"), ("en", "fig0_method")):
        out[idioma] = componer_metodo(img, familias, doc,
                                      Path(carpeta) / nombre, idioma, ctx,
                                      voi, lado_mm)
    return out, fallos
