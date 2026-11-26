"""escena_cargas.py - La escena 3D del analisis con cargas: estructura, cubo, cargas y campos.

POR QUE ESTA EN UN MODULO APARTE
--------------------------------
Es solo pyvista, sin Qt, como `escena_recorte`: la ventana (`DialogoCargas`)
le pasa su `QtInteractor`, y las pruebas y las capturas le pasan un
`pv.Plotter` fuera de pantalla. Asi se puede ver lo que dibuja, y probar lo que
hace el raton, sin abrir la interfaz.

EL MARCO ES EL DE LA ESTRUCTURA
-------------------------------
x, y, z son los ejes de la mascara, en mm, con el voxel (i, j, k) ocupando
[i, i+1] x [j, j+1] x [k, k+1] por el espaciado. Es el mismo marco en el que
`spinpy.cargas` define regiones y nodos, de modo que lo que se pincha y lo que se
resuelve no pueden diferir. (El recorte, en cambio, dibuja en el marco PCA del
hueso: aqui el VOI ya es un cubo recortado y sus ejes son los del ensayo.)

QUE SE VE
---------
* la estructura (hueso, espinodoide o dual-lattice), traslucida si se quiere;
* el cubo como caja, con el nombre de cada cara (+x ... -z) para pinchar por
  cara;
* cada CARGA como una flecha de su color, con su region resaltada: la flecha
  llega a la superficie si empuja hacia dentro y sale de ella si tira; un
  momento es un arco con punta; un desplazamiento impuesto, una flecha verde;
* cada APOYO como un parche azul sobre las caras de hueso que sujeta;
* el RESULTADO: un campo coloreado sobre la superficie, con la deformada
  amplificada y un corte opcional.

DOS FORMAS DE PINTAR UN CAMPO
-----------------------------
`suave` (por omision): la isosuperficie del hueso, coloreada con el valor del
campo en el voxel mas proximo promediado con sus vecinos (3x3x3), y deformada
con la interpolacion trilineal del desplazamiento. Es la que se ve bien.
`voxeles`: las caras de voxel expuestas, cada una con el valor EXACTO de su
elemento y movida por los desplazamientos de sus nodos: es lo que el solver
resolvio, sin suavizar. Para citar un valor se mira esta.

EL RATON
--------
`activar_pinchazo` coloca un observador sobre el boton izquierdo: un clic que no
se arrastra pincha (un arrastre sigue orbitando la camara). Solo se pincha la
estructura y la caja; las flechas y las regiones resaltadas no son pinchables.
"""

from __future__ import annotations

import numpy as np
import pyvista as pv

from . import cargas as C

COLOR_HUESO = "#d9c7a0"
COLOR_CUBO = "#ff8a1f"
FONDO = ("#141a22", "#2c3745")           # (abajo, arriba)
COLORES_CARGA = ("#e0382a", "#f08a24", "#d4b017", "#d6336c", "#9c5de0")
COLORES_APOYO = ("#2f80d1", "#19b394", "#4f5fd6", "#2aa7c9")
COLOR_DESP = "#27ae60"
COLOR_SEL = "#ffe066"

_MATERIAL = dict(smooth_shading=True, ambient=0.2, diffuse=0.8,
                 specular=0.25, specular_power=16)


# ---------------------------------------------------------------------------
# Mallas de la estructura
# ---------------------------------------------------------------------------

def malla_estructura(BW, spacing):
    """Isosuperficie 0.5 de la mascara, CERRADA en los bordes y en el marco de cargas.

    Se rodea de ceros para que cierre en las caras del cubo (se ve un bloque
    macizo, como lo ve el solver) y se lleva al marco donde el voxel i ocupa
    [i, i+1]: los centros de voxel quedan en (i + 0.5) * spacing.
    """
    BW = np.asarray(BW, bool)
    if not BW.any():
        return pv.PolyData()
    sp = C._spacing3(spacing)
    pad = np.pad(BW, 1).astype(np.float32)
    g = pv.ImageData(dimensions=pad.shape, spacing=(1.0, 1.0, 1.0))
    g.point_data["v"] = pad.flatten(order="F")
    m = g.contour([0.5], scalars="v")
    if not m.n_points:
        return pv.PolyData()
    # indice de relleno -> indice de voxel (-1) -> marco de nodos (+0.5)
    m.points = (m.points - 1.0 + 0.5) * sp
    return m.compute_normals(auto_orient_normals=True,
                             consistent_normals=False)


def malla_caras(caras, shape, spacing, q=None, escala_normal=0.0):
    """Las caras de voxel `q` (todas si None) como PolyData de cuadrilateros.

    `escala_normal`: fraccion del espaciado que se separa cada cara a lo largo
    de su normal, para que un resalte no quede enterrado en la superficie.
    """
    if q is None:
        q = np.arange(caras["cent"].shape[0])
    q = np.asarray(q, int)
    if q.size == 0:
        return pv.PolyData()
    sp = C._spacing3(spacing)
    ids, inv = np.unique(caras["nod"][q].ravel(), return_inverse=True)
    pts = C._coord_nodos(ids, shape, spacing)
    f = np.empty((q.size, 5), np.int64)
    f[:, 0] = 4
    f[:, 1:] = inv.reshape(-1, 4)
    m = pv.PolyData(pts, f.ravel())
    if escala_normal:
        # desplazar cada nodo la media de las normales de sus caras
        nn = np.zeros_like(pts)
        np.add.at(nn, inv, np.repeat(caras["nor"][q], 4, axis=0))
        ln = np.linalg.norm(nn, axis=1, keepdims=True)
        m.points = pts + escala_normal * float(sp.min()) * nn / np.where(
            ln > 0, ln, 1.0)
    return m


# ---------------------------------------------------------------------------
# Escena
# ---------------------------------------------------------------------------

class EscenaCargas:
    """Escena de cargas sobre un plotter de pyvista (Qt o fuera de pantalla)."""

    def __init__(self, plotter, BW, spacing, fondo=FONDO, caras=None):
        self.p = plotter
        self.BW = np.asarray(BW, bool)
        self.shape = self.BW.shape
        self.sp = C._spacing3(spacing)
        self.ext = C.extension(self.shape, self.sp)
        self.caras = (caras if caras is not None
                      else C.caras_expuestas(self.BW, self.sp))
        self.malla = malla_estructura(self.BW, self.sp)
        self._nodos_hueso = None
        self.traslucido = True
        self.modo_campo = "suave"
        self._resultado = None            # (res, campo, opciones)
        self._nombres = set()             # actores de cargas y apoyos
        self._mango = None
        self._pinchazo = None
        self._escala = float(np.linalg.norm(self.ext))
        self.flechas = {}                 # indice de carga -> (origen, direccion)
        try:
            self.p.set_background(fondo[0], top=fondo[1])
        except Exception:
            self.p.set_background(fondo[0])
        try:
            self.p.enable_lightkit()
        except Exception:
            pass
        try:
            self.p.add_axes(xlabel="x", ylabel="y", zlabel="z",
                            line_width=2, color="white")
        except Exception:
            pass
        self.poner_cubo()
        self.pintar_estructura()

    # -- datos auxiliares ---------------------------------------------------
    @property
    def nodos_hueso(self):
        if self._nodos_hueso is None:
            nn = np.array(self.shape) + 1
            i, j, k = np.nonzero(self.BW)
            from .resistencia import ESQUINAS
            ids = [(i + c[0]) + nn[0] * (j + c[1]) + nn[0] * nn[1] * (k + c[2])
                   for c in ESQUINAS]
            self._nodos_hueso = np.unique(np.stack(ids, axis=1))
        return self._nodos_hueso

    def seleccion(self, region):
        return C.seleccionar(region, self.caras, self.shape, self.sp,
                             self.nodos_hueso)

    def _no_pinchable(self, actor):
        try:
            actor.SetPickable(False)
        except Exception:
            pass

    # -- estructura y cubo ----------------------------------------------------
    def pintar_estructura(self, render=False):
        if self._resultado is not None or not self.malla.n_points:
            return
        op = 0.35 if self.traslucido else 1.0
        a = self.p.add_mesh(self.malla, name="estructura", color=COLOR_HUESO,
                            opacity=op, render=False, **_MATERIAL)
        a.SetPickable(True)
        if render:
            self.p.render()

    def fijar_traslucido(self, activo):
        self.traslucido = bool(activo)
        self.pintar_estructura(render=True)

    def poner_cubo(self, etiquetas=True):
        b = (0, self.ext[0], 0, self.ext[1], 0, self.ext[2])
        caja = pv.Box(bounds=b)
        a = self.p.add_mesh(caja, name="cubo", color=COLOR_CUBO, opacity=0.04,
                            lighting=False, render=False)
        a.SetPickable(True)
        ar = self.p.add_mesh(caja.extract_all_edges(), name="cubo_aristas",
                             color=COLOR_CUBO, line_width=2, lighting=False,
                             render=False)
        self._no_pinchable(ar)
        if etiquetas:
            c = self.ext / 2.0
            pts, txt = [], []
            for nombre, (ax, s) in C.CARAS.items():
                q = c.copy()
                q[ax] = self.ext[ax] if s > 0 else 0.0
                pts.append(q)
                txt.append(nombre)
            try:
                self.p.add_point_labels(
                    np.array(pts), txt, name="cubo_caras", font_size=13,
                    text_color="#ffd9b0", shape=None, show_points=False,
                    always_visible=True, render=False)
            except Exception:
                pass

    # -- cargas y apoyos --------------------------------------------------------
    def _borrar_caso(self):
        for n in self._nombres:
            try:
                self.p.remove_actor(n, render=False)
            except Exception:
                pass
        self._nombres = set()
        self.flechas = {}

    def _add(self, malla, nombre, **kw):
        a = self.p.add_mesh(malla, name=nombre, render=False, **kw)
        self._no_pinchable(a)
        self._nombres.add(nombre)
        return a

    def direccion_efectiva(self, carga, sel):
        """Vector unitario con el que se dibuja la flecha de una carga.

        Para una direccion `normal_*` es la normal media de las caras de la
        region (la fuerza de cada cara sigue la suya; la flecha resume).
        """
        t = carga["tipo"]
        if t == "momento":
            return C._unitario(carga["eje"], "eje")
        if t == "desplazamiento":
            u = carga.get("u")
            if u is not None:
                v = np.array([0.0 if x is None else x for x in u])
            else:
                dv = C.direccion_vector(carga.get("direccion", "+z"))
                v = np.zeros(3) if dv is None else dv * float(
                    carga.get("magnitud", 1.0))
            n = np.linalg.norm(v)
            return v / n if n > 0 else np.array([0.0, 0.0, -1.0])
        d = carga.get("direccion", "normal_entrante" if t == "presion"
                      else "normal_entrante")
        v = C.direccion_vector(d)
        if v is not None:
            m = float(carga.get("magnitud", 1.0))
            return v if (t != "presion" or m >= 0) else -v
        q = sel["quads"]
        nm = (-self.caras["nor"][q].mean(0) if q is not None and len(q)
              else np.array([0.0, 0.0, -1.0]))
        if str(d).strip().lower() == "normal_saliente":
            nm = -nm
        n = np.linalg.norm(nm)
        if n < 1e-9 and q is not None and len(q):
            nm = -self.caras["nor"][q[0]]
            n = 1.0
        return nm / n if n > 0 else np.array([0.0, 0.0, -1.0])

    def _centro_region(self, sel):
        X = C._coord_nodos(sel["nodos"], self.shape, self.sp)
        w = sel["peso"]
        return (X * w[:, None]).sum(0) / w.sum()

    def _flecha(self, origen, direccion, largo, color, nombre):
        v = np.asarray(direccion, float)
        a = pv.Arrow(start=tuple(origen), direction=tuple(v), scale=largo,
                     tip_length=0.28, tip_radius=0.09, shaft_radius=0.035)
        self._add(a, nombre, color=color, lighting=False)

    def _arco(self, centro, eje, radio, color, nombre):
        e = np.asarray(eje, float)
        a = np.array([1.0, 0, 0]) if abs(e[0]) < 0.9 else np.array([0, 1.0, 0])
        u = np.cross(e, a)
        u /= np.linalg.norm(u)
        v = np.cross(e, u)
        t = np.linspace(0.15 * np.pi, 1.7 * np.pi, 40)
        pts = centro + radio * (np.outer(np.cos(t), u) + np.outer(np.sin(t), v))
        self._add(pv.lines_from_points(pts), nombre, color=color,
                  line_width=5, lighting=False)
        tang = -np.sin(t[-1]) * u + np.cos(t[-1]) * v
        cono = pv.Cone(center=tuple(pts[-1] + 0.3 * radio * tang),
                       direction=tuple(tang), height=0.6 * radio,
                       radius=0.22 * radio, resolution=24)
        self._add(cono, nombre + "_punta", color=color, lighting=False)

    def poner_caso(self, caso, activa=None, activo_apoyo=None, render=True):
        """Dibuja las cargas y los apoyos del caso; `activa` resalta una carga.

        Una region vacia (no toca hueso) se omite en silencio aqui: la ventana
        la senala en su lista, y el solver la rechaza con un mensaje.
        """
        self._borrar_caso()
        L = 0.22 * self._escala
        for i, c in enumerate(caso.get("cargas", [])):
            try:
                sel = self.seleccion(c["region"])
            except C.ErrorCaso:
                continue
            if sel["nodos"].size == 0:
                continue
            col = (COLOR_DESP if c["tipo"] == "desplazamiento"
                   else COLORES_CARGA[i % len(COLORES_CARGA)])
            if sel["quads"] is not None and sel["quads"].size:
                self._add(malla_caras(self.caras, self.shape, self.sp,
                                      sel["quads"], 0.05), "carga_reg_%d" % i,
                          color=COLOR_SEL if i == activa else col,
                          opacity=0.85, lighting=False)
            else:
                X = C._coord_nodos(sel["nodos"], self.shape, self.sp)
                self._add(pv.PolyData(X), "carga_reg_%d" % i, color=col,
                          render_points_as_spheres=True, point_size=9,
                          lighting=False)
            o = self._centro_region(sel)
            d = self.direccion_efectiva(c, sel)
            esc = L * (1.35 if i == activa else 1.0)
            if c["tipo"] == "momento":
                self._arco(o, d, 0.35 * esc, col, "carga_flecha_%d" % i)
                pos_et = o + d * 0.45 * esc
            else:
                # empuja hacia dentro: la punta toca la superficie; tira: sale
                dentro = dentro_flecha(d, sel, self)
                ori = o - d * esc if dentro else o
                self._flecha(ori, d, esc, col, "carga_flecha_%d" % i)
                pos_et = ori - d * 0.08 * esc if dentro else o + d * 1.08 * esc
            self.flechas[i] = (o, d)
            try:
                self.p.add_point_labels(
                    [pos_et], ["C%d" % (i + 1)], name="carga_et_%d" % i,
                    font_size=15, text_color=col, shape=None,
                    show_points=False, always_visible=True, render=False)
                self._nombres.add("carga_et_%d" % i)
            except Exception:
                pass
        for i, a in enumerate(caso.get("apoyos", [])):
            try:
                sel = self.seleccion(a["region"])
            except C.ErrorCaso:
                continue
            if sel["nodos"].size == 0:
                continue
            col = COLORES_APOYO[i % len(COLORES_APOYO)]
            if sel["quads"] is not None and sel["quads"].size:
                self._add(malla_caras(self.caras, self.shape, self.sp,
                                      sel["quads"], 0.05), "apoyo_reg_%d" % i,
                          color=COLOR_SEL if i == activo_apoyo else col,
                          opacity=0.8, lighting=False)
            else:
                X = C._coord_nodos(sel["nodos"], self.shape, self.sp)
                self._add(pv.PolyData(X), "apoyo_reg_%d" % i, color=col,
                          render_points_as_spheres=True, point_size=9,
                          lighting=False)
        if render:
            self.p.render()

    # -- resultado -----------------------------------------------------------------
    def _valores_suaves(self, campo):
        from scipy import ndimage
        ok = np.isfinite(campo)
        num = ndimage.uniform_filter(np.where(ok, campo, 0.0), 3, mode="nearest")
        den = ndimage.uniform_filter(ok.astype(float), 3, mode="nearest")
        with np.errstate(invalid="ignore", divide="ignore"):
            suave = np.where(den > 0.05, num / den, np.nan)
        P = self.malla.points / self.sp
        idx = np.clip(np.floor(P).astype(int), 0, np.array(self.shape) - 1)
        return suave[idx[:, 0], idx[:, 1], idx[:, 2]]

    def mostrar_resultado(self, res, campo, clim=None, deformada=False,
                          escala=1.0, corte=None, barra=True, titulo=None):
        """Colorea la superficie con `res['campos'][campo]`.

        `corte` = (eje 0-2, posicion en mm, lado +1/-1) para ver el interior.
        Devuelve el rango (min, max) usado.
        """
        self._resultado = (res, campo, dict(
            clim=clim, deformada=deformada, escala=escala, corte=corte,
            barra=barra, titulo=titulo))
        vals = np.asarray(res["campos"][campo], float)
        if self.modo_campo == "voxeles":
            m, s = self._malla_voxeles(res, vals, deformada, escala)
        else:
            m = self.malla.copy()
            s = self._valores_suaves(vals)
            if deformada:
                m.points = C.desplazar_puntos(res, m.points, escala)
        if clim is None:
            fin = s[np.isfinite(s)]
            clim = ((float(fin.min()), float(fin.max())) if fin.size
                    else (0.0, 1.0))
            if clim[0] == clim[1]:
                clim = (clim[0] - 1e-12, clim[1] + 1e-12)
        m["valor"] = s
        if corte is not None and m.n_points:
            ax, pos, lado = corte
            n = np.zeros(3)
            n[ax] = 1.0
            o = np.zeros(3)
            o[ax] = pos
            c = m.clip(normal=n, origin=o, invert=(lado < 0))
            if c.n_points:
                m = c
        self.p.remove_actor("estructura", render=False)
        # El titulo va como texto aparte: dentro de la barra horizontal se
        # solapa con sus cifras.
        titulo = titulo or C.ETIQUETA_CAMPO.get(campo, campo)
        try:
            from .idioma import _
            titulo = _(titulo)
        except Exception:
            pass
        self.p.add_text(titulo, position="upper_right", font_size=11,
                        color="white", name="titulo_campo", render=False)
        self.p.add_mesh(
            m, name="campo", scalars="valor", cmap="turbo", clim=clim,
            nan_color="#8d8d8d", show_scalar_bar=barra, render=False,
            scalar_bar_args={"title": " ", "n_labels": 5, "color": "white",
                             "vertical": False, "height": 0.06,
                             "position_y": 0.02, "fmt": "%.3g"},
            **_MATERIAL).SetPickable(True)
        self.p.render()
        return clim

    def _malla_voxeles(self, res, vals, deformada, escala):
        caras = self.caras
        v = caras["vox"]
        s_cara = vals[v[:, 0], v[:, 1], v[:, 2]]
        m = malla_caras(caras, self.shape, self.sp)
        if deformada:
            ids = np.unique(caras["nod"].ravel())
            nn = np.array(self.shape) + 1
            c, r = np.divmod(ids, nn[0] * nn[1])
            b, a = np.divmod(r, nn[0])
            du = res["u_malla"][a, b, c].astype(float)
            m.points = m.points + escala * du
        sc = np.zeros(m.n_cells)
        sc[:] = s_cara
        m.cell_data["valor"] = sc
        m = m.cell_data_to_point_data()
        return m, np.asarray(m.point_data["valor"], float)

    def quitar_resultado(self):
        self._resultado = None
        for n in ("campo", "titulo_campo"):
            self.p.remove_actor(n, render=False)
        try:
            self.p.remove_scalar_bar(render=False)
        except Exception:
            pass
        self.pintar_estructura()
        self.p.render()

    def fijar_modo_campo(self, modo):
        self.modo_campo = "voxeles" if modo == "voxeles" else "suave"
        if self._resultado is not None:
            res, campo, o = self._resultado
            self.mostrar_resultado(res, campo, **o)

    # -- raton ------------------------------------------------------------------------
    def activar_pinchazo(self, al_pinchar):
        """`al_pinchar(punto, shift)` se llama con el clic que no se arrastra."""
        self.desactivar_pinchazo()
        try:
            from vtkmodules.vtkRenderingCore import vtkCellPicker
        except Exception:
            return False
        picker = vtkCellPicker()
        picker.SetTolerance(0.002)
        estado = {"xy": None}
        iren = self.p.iren

        def pulsar(_o, _e):
            estado["xy"] = iren.interactor.GetEventPosition()

        def soltar(_o, _e):
            if estado["xy"] is None:
                return
            x0, y0 = estado["xy"]
            x, y = iren.interactor.GetEventPosition()
            estado["xy"] = None
            if abs(x - x0) > 4 or abs(y - y0) > 4:
                return                    # fue un arrastre: la camara orbito
            self._pinchar_pantalla(x, y, al_pinchar,
                                   bool(iren.interactor.GetShiftKey()), picker)

        try:
            a = iren.add_observer("LeftButtonPressEvent", pulsar)
            b = iren.add_observer("LeftButtonReleaseEvent", soltar)
        except Exception:
            return False
        self._pinchazo = (a, b, picker, al_pinchar)
        return True

    def _pinchar_pantalla(self, x, y, al_pinchar, shift=False, picker=None):
        if picker is None:
            from vtkmodules.vtkRenderingCore import vtkCellPicker
            picker = vtkCellPicker()
            picker.SetTolerance(0.002)
        if not picker.Pick(x, y, 0, self.p.renderer):
            return None
        pt = np.array(picker.GetPickPosition(), float)
        al_pinchar(pt, shift)
        return pt

    def desactivar_pinchazo(self):
        if self._pinchazo is None:
            return
        for o in self._pinchazo[:2]:
            try:
                self.p.iren.remove_observer(o)
            except Exception:
                pass
        self._pinchazo = None

    def activar_mango(self, origen, direccion, largo, al_mover):
        """Esfera en la punta de la flecha: arrastrarla gira y alarga la fuerza."""
        self.quitar_mango()
        origen = np.asarray(origen, float)
        punta = origen + np.asarray(direccion, float) * largo

        def mover(c):
            v = np.asarray(c, float) - origen
            n = float(np.linalg.norm(v))
            if n > 1e-9:
                al_mover(v / n, n)

        try:
            self._mango = self.p.add_sphere_widget(
                mover, center=tuple(punta), radius=0.03 * self._escala,
                color=COLOR_SEL, selected_color="#ffffff",
                interaction_event="always", test_callback=False)
            return True
        except Exception:
            self._mango = None
            return False

    def quitar_mango(self):
        if self._mango is None:
            return
        self._mango = None
        try:
            self.p.clear_sphere_widgets()
        except Exception:
            pass

    # -- camara ----------------------------------------------------------------------
    def vista(self, nombre):
        f = {"xy": self.p.view_xy, "xz": self.p.view_xz,
             "yz": self.p.view_yz}.get(nombre)
        if f is not None:
            f()
            if nombre == "xz":
                self.p.camera.up = (0, 0, 1)
        else:
            self.p.view_isometric()
        self.p.reset_camera()
        self.p.render()

    def captura(self, ruta):
        """Guarda la vista actual en PNG; devuelve la ruta."""
        self.p.screenshot(str(ruta))
        return ruta


def dentro_flecha(d, sel, escena):
    """True si la flecha dibujada empuja hacia dentro de la superficie."""
    q = sel["quads"]
    if q is None or not len(q):
        return True
    return float(d @ escena.caras["nor"][q].mean(0)) <= 0
