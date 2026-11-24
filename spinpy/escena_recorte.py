"""escena_recorte.py - La escena 3D del recorte del VOI: hueso, cubo y tramos.

POR QUE ESTA EN UN MODULO APARTE
--------------------------------
Es solo pyvista, sin Qt: la ventana (`visor.DialogoRecorte`) le pasa su
`QtInteractor`, y las pruebas y las capturas le pasan un `pv.Plotter` fuera de
pantalla. Asi se puede ver lo que dibuja sin abrir la interfaz.

TODO SE DIBUJA EN EL MARCO PCA
------------------------------
Eje x = PC1 (eje mayor del hueso), y = PC2, z = PC3. En ese marco el cubo es
una caja alineada con los ejes, que es lo que permite recortar el hueso con una
caja y dibujarlo sin matrices; y el hueso aparece siempre tumbado a lo largo
de x, oriente como oriente se escaneo, que es lo que no hace una vista de
cortes en los ejes del escaner.

QUE SE VE
---------
* el hueso entero, traslucido, a resolucion reducida (la de la figura del
  metodo, `figura_metodo.reducir`): es contexto, no medida;
* el contenido del cubo, opaco y de otro color, a resolucion COMPLETA y tal y
  como lo recorta `extraer_cubo`: lo que se ve ahi es lo que se va a medir;
* el cubo como caja con aristas y caras tenues, con su lado en la propia
  escena, o el cilindro (radio y altura) con las aristas de sus tapas;
* los tres tramos sugeridos, numerados.

La malla del hueso sale de un volumen reducido por bloques, asi que es una
APROXIMACION de la superficie (trabeculas finas pueden perderse). No se mide
nada sobre ella.
"""

from __future__ import annotations

import numpy as np
import pyvista as pv

from .voi import _spacing3, cilindro_a_pca, cubo_a_pca, mundo_a_pca

COLOR_HUESO = "#d9c7a0"
COLOR_CUBO = "#ff8a1f"
COLOR_TRAMO = "#4fb3ff"
FONDO = ("#141a22", "#2c3745")           # (abajo, arriba)

# Material con un punto de especularidad, igual que los paneles principales:
# con material puramente difuso las trabeculas se funden unas con otras.
_MATERIAL = dict(smooth_shading=True, ambient=0.2, diffuse=0.8,
                 specular=0.25, specular_power=16)


def malla_hueso_pca(reducido, spacing, ejes, centro):
    """Superficie del hueso en el marco PCA, desde el volumen reducido.

    `reducido` es `(volumen float32, p)` de `figura_metodo.reducir`: media por
    bloques de p voxeles. El centro del bloque i esta en (i*p + (p-1)/2) voxeles,
    y se contornea a 0.5, la misma convencion que `malla_de_mascara` del visor.
    """
    vol, p = reducido
    sp = _spacing3(spacing)
    grid = pv.ImageData(dimensions=vol.shape, spacing=tuple(p * sp),
                        origin=tuple(0.5 * (p - 1) * sp))
    grid.point_data["v"] = np.asarray(vol, np.float32).flatten(order="F")
    m = grid.contour([0.5], scalars="v")
    if not m.n_points:
        return pv.PolyData()
    m.points = mundo_a_pca(m.points, ejes, centro)
    return m.compute_normals(auto_orient_normals=True,
                             consistent_normals=False)


def _isosuperficie_cerrada(volumen):
    """Isosuperficie 0.5 de un bool 3D, en indices de voxel, CERRADA en los bordes.

    Se rodea de una capa de ceros para que la superficie cierre en las caras del
    volumen: se ve un bloque macizo, como la probeta que se va a medir, en vez
    de una piel abierta por los seis lados. El -1 deshace el relleno; el borde
    queda medio voxel por fuera del ultimo centro, que es lo que ocupa el voxel.
    Devuelve None si no hay superficie.
    """
    pad = np.pad(np.asarray(volumen, bool), 1).astype(np.float32)
    grid = pv.ImageData(dimensions=pad.shape, spacing=(1, 1, 1))
    grid.point_data["v"] = pad.flatten(order="F")
    m = grid.contour([0.5], scalars="v")
    if not m.n_points:
        return None
    m.points = m.points - 1.0
    return m


def _con_normales(m):
    return m.compute_normals(auto_orient_normals=True,
                             consistent_normals=False)


def malla_cubo_pca(cubo, spacing, centro_pca):
    """Superficie del contenido del cubo, en el marco PCA y a resolucion plena."""
    cubo = np.asarray(cubo, bool)
    if not cubo.any():
        return pv.PolyData()
    m = _isosuperficie_cerrada(cubo)
    if m is None:
        return pv.PolyData()
    m.points = cubo_a_pca(m.points, cubo.shape[0], spacing, centro_pca)
    return _con_normales(m)


def malla_cilindro_pca(volumen, dominio, spacing, centro_pca, eje):
    """Superficie del contenido del cilindro, en el marco PCA y a resolucion plena.

    `volumen` y `dominio` son los de `voi.extraer_cilindro`. Se cierra contra la
    pared curva, no contra la caja: lo que se ve es el cilindro de hueso que se
    va a medir.
    """
    volumen = np.asarray(volumen, bool) & np.asarray(dominio, bool)
    if not volumen.any():
        return pv.PolyData()
    m = _isosuperficie_cerrada(volumen)
    if m is None:
        return pv.PolyData()
    m.points = cilindro_a_pca(m.points, volumen.shape, spacing, centro_pca,
                              eje)
    return _con_normales(m)


def malla_region_pca(recorte, spacing):
    """Superficie del contenido de un recorte de `voi.recortar` (cubo o cilindro)."""
    reg = recorte["region"]
    if reg["forma"] == "cilindro":
        return malla_cilindro_pca(recorte["cubo"], recorte["dominio"], spacing,
                                  recorte["centro_pca"], reg.get("eje", 0))
    return malla_cubo_pca(recorte["cubo"], spacing, recorte["centro_pca"])


def _limites(centro_pca, lado_mm):
    c = np.asarray(centro_pca, float)
    h = 0.5 * float(lado_mm)
    return (c[0] - h, c[0] + h, c[1] - h, c[1] + h, c[2] - h, c[2] + h)


class EscenaRecorte:
    """Escena del recorte sobre un plotter de pyvista (Qt o fuera de pantalla)."""

    def __init__(self, plotter, malla_hueso, lo, hi, fondo=FONDO):
        self.p = plotter
        self.malla = malla_hueso
        self.lo, self.hi = np.asarray(lo, float), np.asarray(hi, float)
        self.traslucido = True
        self.seccionado = False
        self._centro = None
        self._lado = None
        self._contenido = None
        self._mango = None
        self._escala = float(np.linalg.norm(self.hi - self.lo))

        try:
            self.p.set_background(fondo[0], top=fondo[1])
        except Exception:
            self.p.set_background(fondo[0])
        try:
            self.p.enable_lightkit()
        except Exception:
            pass
        try:
            self.p.add_axes(xlabel="PC1", ylabel="PC2", zlabel="PC3",
                            line_width=2, color="white")
        except Exception:
            pass
        self._pintar_hueso()

    # -- hueso ---------------------------------------------------------------
    def _malla_visible(self):
        m = self.malla
        if self.seccionado and self._centro is not None and m.n_points:
            # Se quita la mitad que da a la camara (z > z del cubo) para ver
            # el cubo DENTRO del hueso.
            cort = m.clip(normal=(0, 0, 1), origin=(0, 0, self._centro[2]),
                          invert=True)
            if cort.n_points:
                return cort
        return m

    def _pintar_hueso(self, render=False):
        m = self._malla_visible()
        if not m.n_points:
            return
        # Opaco al seccionar: con el corte ya se ve el interior, y la
        # transparencia encima de un corte solo ensucia.
        op = 0.28 if (self.traslucido and not self.seccionado) else 1.0
        self.p.add_mesh(m, name="hueso", color=COLOR_HUESO, opacity=op,
                        render=render, **_MATERIAL)

    def fijar_traslucido(self, activo):
        self.traslucido = bool(activo)
        self._pintar_hueso(render=True)

    def fijar_seccion(self, activo):
        self.seccionado = bool(activo)
        self._pintar_hueso()
        self._pintar_contenido()
        self.p.render()

    # -- cubo y cilindro ---------------------------------------------------------
    def poner_cubo(self, centro_pca, lado_mm, texto=None, render=True):
        """Mueve la caja del cubo. Barato: se puede llamar mientras se arrastra."""
        self._centro = np.asarray(centro_pca, float)
        self._lado = float(lado_mm)
        caja = pv.Box(bounds=_limites(self._centro, self._lado))
        self._dibujar_forma(caja, _limites(self._centro, self._lado), texto,
                            render)

    def poner_cilindro(self, centro_pca, radio_mm, altura_mm, eje, texto=None,
                       render=True):
        """Mueve el cilindro (su eje es el eje PCA `eje`: 0 = PC1, 1 = PC2...)."""
        self._centro = np.asarray(centro_pca, float)
        self._lado = 2.0 * float(radio_mm)
        direccion = [0.0, 0.0, 0.0]
        direccion[int(eje)] = 1.0
        cil = pv.Cylinder(center=tuple(self._centro), direction=direccion,
                          radius=float(radio_mm), height=float(altura_mm),
                          resolution=64, capping=True)
        semi = np.full(3, float(radio_mm))
        semi[int(eje)] = 0.5 * float(altura_mm)
        b = (self._centro[0] - semi[0], self._centro[0] + semi[0],
             self._centro[1] - semi[1], self._centro[1] + semi[1],
             self._centro[2] - semi[2], self._centro[2] + semi[2])
        self._dibujar_forma(cil, b, texto, render)

    def _dibujar_forma(self, malla, limites, texto, render):
        """Caja o cilindro, con aristas y etiqueta. Los actores se reemplazan
        por nombre, asi que cambiar de forma no deja restos de la anterior."""
        self.p.add_mesh(malla, name="caja", color=COLOR_CUBO, opacity=0.05,
                        lighting=False, render=False)
        # Aristas: las del cubo son todas; del cilindro, solo las de las tapas
        # (angulo de 30 grados: la pared curva es una superficie suave).
        aristas = (malla.extract_feature_edges(feature_angle=30.0)
                   if malla.n_cells > 6 else malla.extract_all_edges())
        self.p.add_mesh(aristas, name="aristas", color=COLOR_CUBO,
                        line_width=3, lighting=False, render=False)
        if texto:
            # La etiqueta cuelga de una esquina superior para no tapar el
            # contenido.
            b = limites
            self.p.add_point_labels(
                [[b[0], b[3], b[5]]], [texto], name="etiqueta_cubo",
                font_size=16, text_color=COLOR_CUBO, shape=None,
                show_points=False, always_visible=True, render=False)
        if self.seccionado:
            self._pintar_hueso()
            self._pintar_contenido()
        if self._mango is not None:
            self.mover_mango(self._centro)
        if render:
            self.p.render()

    def poner_contenido(self, malla_cubo, render=True):
        """Pinta el contenido ya mallado del cubo (ver `malla_cubo_pca`)."""
        self._contenido = malla_cubo
        self._pintar_contenido()
        if render:
            self.p.render()

    def _pintar_contenido(self):
        m = self._contenido
        if m is not None and self.seccionado and self._centro is not None \
                and m.n_points:
            # El mismo plano que al hueso: el cubo aparece cortado y se ve su
            # interior en naranja contra el hueso crema.
            cort = m.clip(normal=(0, 0, 1), origin=(0, 0, self._centro[2]),
                          invert=True)
            m = cort if cort.n_points else m
        if m is None or not m.n_points:
            self.p.remove_actor("contenido", render=False)
        else:
            self.p.add_mesh(m, name="contenido", color=COLOR_CUBO,
                            render=False, **_MATERIAL)

    # -- tramos --------------------------------------------------------------
    def poner_tramos(self, centros_pca, etiquetas, render=True):
        c = np.asarray(centros_pca, float)
        self.p.add_mesh(pv.PolyData(c), name="tramos", color=COLOR_TRAMO,
                        render_points_as_spheres=True, point_size=14,
                        lighting=False, render=False)
        self.p.add_point_labels(
            c, list(etiquetas), name="etiquetas_tramos", font_size=18,
            text_color=COLOR_TRAMO, shape=None, show_points=False,
            always_visible=True, point_size=1, render=False)
        if render:
            self.p.render()

    # -- mango para arrastrar el cubo -----------------------------------------
    def activar_mango(self, al_mover, al_soltar=None):
        """Esfera que se arrastra con el raton; mueve el cubo.

        `al_mover(centro)` se llama durante el arrastre y `al_soltar(centro)`
        al soltar. Falla en silencio si el widget no esta disponible (depende
        del backend de VTK): los mandos siguen sirviendo.
        """
        if self._centro is None:
            return False
        radio = 0.035 * self._escala

        def mover(c):
            al_mover(np.asarray(c, float))

        try:
            self._mango = self.p.add_sphere_widget(
                mover, center=tuple(self._centro), radius=radio,
                color=COLOR_CUBO, selected_color="#ffd166",
                interaction_event="always", test_callback=False)
            if al_soltar is not None:
                self._mango.AddObserver(
                    "EndInteractionEvent",
                    lambda *_a: al_soltar(
                        np.asarray(self._mango.GetCenter(), float)))
            return True
        except Exception:
            self._mango = None
            return False

    def quitar_mango(self):
        self._mango = None
        try:
            self.p.clear_sphere_widgets()
            self.p.render()
        except Exception:
            pass

    def mover_mango(self, centro_pca):
        """Lleva el mango al centro del cubo cuando este se mueve por otro lado."""
        if self._mango is None:
            return
        try:
            self._mango.SetCenter(*[float(v) for v in centro_pca])
        except Exception:
            pass

    # -- camara ----------------------------------------------------------------
    def vista(self, nombre):
        """'iso', 'long12' (PC1-PC2), 'long13' (PC1-PC3) o 'transversal' (PC2-PC3)."""
        f = {"long12": self.p.view_xy, "long13": self.p.view_xz,
             "transversal": self.p.view_yz}.get(nombre)
        if f is not None:
            f()
            if nombre == "long13":
                # view_xz mira desde -y; se endereza para que PC3 quede arriba.
                self.p.camera.up = (0, 0, 1)
        else:
            self.p.view_isometric()
        self.p.reset_camera()
        self.p.render()
