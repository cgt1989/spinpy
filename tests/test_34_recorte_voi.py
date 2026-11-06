"""
Bloque 34 - Recorte del VOI: tamano absoluto, posicion del cubo y su dibujo.

QUE SE VERIFICA Y POR QUE
-------------------------
El dialogo de recorte dibuja el cubo sobre el hueso en 3D y en cortes planos, y
deja pedirlo por su volumen absoluto. Es geometria que no da sintomas cuando
falla: un cubo dibujado girado o desplazado respecto al que se recorta tiene
el mismo aspecto que uno correcto, y quien lo coloca mirando la pantalla recorta
otra region sin enterarse. Por eso se comprueba contra `extraer_cubo`, que es
lo que de verdad se mide:

  34.1  `dimensionar_cubo`: el volumen pedido da el mismo cubo que el lado
        equivalente, el lado real es un numero IMPAR de voxeles y la
        desviacion respecto a lo pedido se declara.
  34.2  `cubo_a_pca`: cada voxel del array que devuelve `extraer_cubo`, llevado
        al marco PCA y de vuelta a la pila, cae sobre un voxel de la pila con
        el mismo valor. Sin la correccion de la transposicion (1, 0, 2) esto
        falla en cuanto el hueso no es simetrico.
  34.3  `GuiaPCA` es la inversa de si misma y reproduce `centro_en` hasta la
        resolucion de su histograma.
  34.4  `perfil_bvtv` calcula el mismo cubo que `extraer_cubo` y se puede
        cancelar.
  34.5  `seccion_pca` pinta el corte por el centro del cubo: sobre un bloque
        macizo de lado conocido, el area del corte es la del bloque.
  34.6  `escena_recorte.malla_cubo_pca` dibuja el contenido del cubo donde
        esta: sobre un bloque asimetrico, la caja de la malla en el marco PCA
        es la interseccion bloque-cubo. Y la escena se compone sin errores.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  dimensionar_cubo   exacto (aritmetica de enteros) y 1e-12 en los reales
  cubo_a_pca         0 voxeles de discrepancia
  GuiaPCA            ida y vuelta 1e-9 mm; frente a centro_en, 0.25 mm
                     (cuatro cajas del histograma: la ventana de centro_en
                     tiene bordes exactos y la de GuiaPCA los redondea a caja)
  perfil_bvtv        exacto (mismo calculo)
  seccion_pca        exacto sobre un bloque alineado con los ejes
  malla_cubo_pca     1e-9 mm: el isovalor 0.5 cae medio voxel fuera del
                     ultimo voxel macizo, y ese valor se conoce de antemano
"""
import numpy as np
import pytest

from spinpy import grf
from spinpy import voi as V

BLOQUE = "34 · recorte del VOI"
REF = "Definicion de `voi.extraer_cubo` (port de VOIsTIFF.m)."

H = 0.0515


def _hueso(n=(96, 80, 64), h=H, semilla=3):
    """Elipsoide girado y poroso, SIN simetria: detecta ejes cruzados."""
    rng = np.random.default_rng(semilla)
    idx = np.indices(n).reshape(3, -1).T * h
    c = np.array(n) * h / 2
    R = grf.euler_R(12, 25, 40)
    q = (idx - c) @ R
    dentro = ((q / np.array([1.6, 1.0, 0.8])) ** 2).sum(1) <= 1.0
    poros = rng.random(len(idx)) < 0.30
    return (dentro & ~poros).reshape(n)


@pytest.fixture(scope="module")
def marco():
    mask = _hueso()
    ejes, centro, proy = V.marco_pca(mask, H)
    return mask, ejes, centro, proy


# -- 34.1 --------------------------------------------------------------------

def test_dimensionar_por_lado_y_por_volumen(registro):
    d_lado = V.dimensionar_cubo(0.051489, lado_mm=5.0)
    d_vol = V.dimensionar_cubo(0.051489, volumen_mm3=125.0)      # 5^3 mm3
    ok = (d_lado["lado_vox"] == d_vol["lado_vox"] == 97
          and d_lado["lado_vox"] % 2 == 1
          and np.isclose(d_lado["lado_real_mm"], 97 * 0.051489, rtol=1e-12)
          and np.isclose(d_lado["volumen_real_mm3"],
                         (97 * 0.051489) ** 3, rtol=1e-12)
          and np.isclose(d_vol["lado_pedido_mm"], 5.0, rtol=1e-12))
    registro.anotar(BLOQUE, "volumen pedido = lado equivalente", REF,
                    97, d_vol["lado_vox"], "exacto", "lado_vox impar", ok)
    assert ok
    # la desviacion se declara: 97 voxeles no son exactamente 5000 um
    assert d_lado["desviacion_rel"] == pytest.approx(97 * 0.051489 / 5.0 - 1)


def test_dimensionar_coincide_con_lo_que_recorta(marco):
    """El lado anunciado es el del array que devuelve `extraer_cubo`."""
    mask, ejes, centro, proy = marco
    c = V.centros_por_tercios(proy)[1]
    for lado in (0.9, 1.4, 2.01):
        cubo, _f = V.extraer_cubo(mask, H, c, ejes, centro, lado)
        assert cubo.shape[0] == V.dimensionar_cubo(H, lado_mm=lado)["lado_vox"]


def test_dimensionar_rechaza_entradas_incoherentes():
    with pytest.raises(ValueError):
        V.dimensionar_cubo(H)
    with pytest.raises(ValueError):
        V.dimensionar_cubo(H, lado_mm=1.0, volumen_mm3=1.0)
    for malo in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            V.dimensionar_cubo(H, volumen_mm3=malo)
        with pytest.raises(ValueError):
            V.dimensionar_cubo(H, lado_mm=malo)


# -- 34.2 --------------------------------------------------------------------

def test_cubo_a_pca_cae_sobre_el_mismo_voxel(registro, marco):
    mask, ejes, centro, proy = marco
    c = V.centros_por_tercios(proy)[0]
    cubo, _f = V.extraer_cubo(mask, H, c, ejes, centro, 1.0)
    n = cubo.shape[0]
    idx = np.indices((n, n, n)).reshape(3, -1).T
    mundo = V.pca_a_mundo(V.cubo_a_pca(idx, n, H, c), ejes, centro)
    vox = np.rint(mundo / H).astype(int)
    dentro = np.all((vox >= 0) & (vox < np.array(mask.shape)), axis=1)
    esperado = cubo[idx[dentro, 0], idx[dentro, 1], idx[dentro, 2]]
    leido = mask[vox[dentro, 0], vox[dentro, 1], vox[dentro, 2]]
    fallos = int((esperado != leido).sum())
    registro.anotar(BLOQUE, "cubo_a_pca reproduce extraer_cubo", REF, 0,
                    fallos, "0 voxeles", f"{int(dentro.sum())} voxeles "
                    "comparados", fallos == 0)
    assert dentro.sum() > 1000 and cubo.any()
    assert fallos == 0


def test_cubo_a_pca_detecta_la_transposicion(marco):
    """Sin la permutacion (1, 0, 2) el mapeo NO coincide: la prueba anterior
    discrimina de verdad y no pasaria con la identidad."""
    mask, ejes, centro, proy = marco
    c = V.centros_por_tercios(proy)[0]
    cubo, _f = V.extraer_cubo(mask, H, c, ejes, centro, 1.0)
    n = cubo.shape[0]
    idx = np.indices((n, n, n)).reshape(3, -1).T
    ingenuo = (idx - 0.5 * (n - 1)) * H + c            # sin transponer
    vox = np.rint(V.pca_a_mundo(ingenuo, ejes, centro) / H).astype(int)
    dentro = np.all((vox >= 0) & (vox < np.array(mask.shape)), axis=1)
    leido = mask[vox[dentro, 0], vox[dentro, 1], vox[dentro, 2]]
    assert (cubo[idx[dentro, 0], idx[dentro, 1], idx[dentro, 2]]
            != leido).sum() > 0


def test_mundo_y_pca_son_inversas(marco):
    _m, ejes, centro, _p = marco
    p = np.random.default_rng(0).normal(size=(20, 3))
    np.testing.assert_allclose(
        V.mundo_a_pca(V.pca_a_mundo(p, ejes, centro), ejes, centro), p,
        atol=1e-12)


# -- 34.3 --------------------------------------------------------------------

def test_guia_ida_y_vuelta_y_frente_a_centro_en(registro, marco):
    _m, _e, _c, proy = marco
    g = V.GuiaPCA(proy)
    peor_rt, peor_ce = 0.0, 0.0
    for t in (0.0, 0.07, 0.25, 0.5, 0.81, 1.0):
        for o2, o3 in ((0, 0), (0.3, -0.2), (-0.5, 0.4)):
            c = g.centro(t, o2, o3)
            t2, a2, a3 = g.parametros(c)
            peor_rt = max(peor_rt, abs(t2 - t) * g.L, abs(a2 - o2),
                          abs(a3 - o3))
            ref = V.centro_en(proy, t, o2, o3)
            peor_ce = max(peor_ce, float(np.abs(c - ref).max()))
    registro.anotar(BLOQUE, "GuiaPCA ida y vuelta", REF, 0, peor_rt,
                    "1e-9 mm", "centro -> parametros -> centro", peor_rt < 1e-9)
    registro.anotar(BLOQUE, "GuiaPCA frente a centro_en", REF, 0, peor_ce,
                    "0.25 mm", "max |diferencia| en PC1, PC2, PC3",
                    peor_ce < 0.25)
    assert peor_rt < 1e-9
    assert peor_ce < 0.25


def test_guia_acota_t_y_rechaza_hueso_sin_extension():
    g = V.GuiaPCA(np.random.default_rng(1).normal(size=(500, 3)))
    assert g.parametros(g.centro(5.0))[0] == 1.0
    assert g.centro(-3.0)[0] == pytest.approx(g.lo[0])
    with pytest.raises(ValueError):
        V.GuiaPCA(np.zeros((10, 3)))


# -- 34.4 --------------------------------------------------------------------

def test_perfil_es_el_cubo_que_se_recorta(registro, marco):
    mask, ejes, centro, proy = marco
    g = V.GuiaPCA(proy)
    pf = V.perfil_bvtv(mask, H, ejes, centro, g, 0.8, n=5)
    cubo, fuera = V.extraer_cubo(mask, H, g.centro(0.5), ejes, centro, 0.8)
    ok = bool(np.isclose(pf["bvtv"][2], cubo.mean(), rtol=0, atol=0)
              and np.isclose(pf["fuera"][2], fuera, rtol=0, atol=0))
    registro.anotar(BLOQUE, "perfil_bvtv = extraer_cubo", REF,
                    float(cubo.mean()), float(pf["bvtv"][2]), "exacto",
                    "t = 0.5", ok)
    assert ok and len(pf["t"]) == 5 and pf["x_mm"][0] == pytest.approx(g.lo[0])


def test_perfil_se_puede_cancelar(marco):
    mask, ejes, centro, proy = marco
    llamadas = []

    def cancelar():
        llamadas.append(1)
        return len(llamadas) > 2

    assert V.perfil_bvtv(mask, H, ejes, centro, V.GuiaPCA(proy), 0.8, n=9,
                         cancelar=cancelar) is None


# -- 34.5 --------------------------------------------------------------------

def test_seccion_de_un_bloque_macizo(registro):
    h = 0.1
    mask = np.zeros((60, 50, 40), bool)
    mask[10:40, 5:35, 8:28] = True                  # 30 x 30 x 20 voxeles
    ejes, centro = np.eye(3), np.zeros(3)
    lo = np.zeros(3)
    hi = (np.array(mask.shape) - 1) * h
    c = np.array([2.5, 2.0, 1.8])
    img, (a0, a1, b0, b1) = V.seccion_pca(mask, h, ejes, centro, c, 0, 1, lo,
                                          hi)
    # plano z = 1.8 mm (voxel 18, dentro del bloque): 30 x 30 voxeles
    esperado = 30 * 30
    registro.anotar(BLOQUE, "seccion_pca de un bloque macizo", REF, esperado,
                    int(img.sum()), "exacto",
                    "area del corte en voxeles", img.sum() == esperado)
    assert img.sum() == esperado
    assert img.shape == (50, 60) and (a0, b0) == (0.0, 0.0)
    # el mismo plano fuera del bloque esta vacio
    vacio, _ = V.seccion_pca(mask, h, ejes, centro, [2.5, 2.0, 0.2], 0, 1,
                             lo, hi)
    assert not vacio.any()


# -- 34.6 --------------------------------------------------------------------

def test_malla_del_cubo_esta_donde_el_cubo(registro):
    """Bloque con tres extensiones distintas, ejes del escaner = ejes PCA: si
    el dibujo del cubo cruzara dos ejes, la caja saldria permutada."""
    from spinpy.escena_recorte import malla_cubo_pca
    h = 0.1
    mask = np.zeros((70, 60, 50), bool)
    mask[20:41, 10:31, 5:26] = True        # voxeles 20..40, 10..30, 5..25
    ejes, centro = np.eye(3), np.zeros(3)
    c = np.array([30.0, 20.0, 15.0]) * h   # cubo de 31 voxeles: 15..45, 5..35, 0..30
    cubo, fuera = V.extraer_cubo(mask, h, c, ejes, centro, 31 * h)
    m = malla_cubo_pca(cubo, h, c)
    esperado = np.array([[19.5, 40.5], [9.5, 30.5], [4.5, 25.5]]) * h
    obtenido = np.array(m.bounds).reshape(3, 2)
    err = float(np.abs(obtenido - esperado).max())
    registro.anotar(BLOQUE, "malla_cubo_pca: caja = bloque dentro del cubo",
                    REF, 0.0, err, "1e-9 mm",
                    "caja de la malla en el marco PCA", err < 1e-9)
    assert fuera == 0.0 and err < 1e-9


def test_malla_del_cubo_vacio_no_falla():
    from spinpy.escena_recorte import malla_cubo_pca
    assert malla_cubo_pca(np.zeros((9, 9, 9), bool), 0.1,
                          np.zeros(3)).n_points == 0


def test_escena_se_compone(render_disponible, marco):
    if not render_disponible:
        pytest.skip("VTK no puede dibujar fuera de pantalla en este equipo")
    import pyvista as pv
    from spinpy import escena_recorte as ER
    from spinpy import figura_metodo as FM
    mask, ejes, centro, proy = marco
    g = V.GuiaPCA(proy)
    lo, hi = V.rango_pca(proy)
    malla = ER.malla_hueso_pca(FM.reducir(mask), H, ejes, centro)
    assert malla.n_points > 0
    # el hueso mallado cabe en la caja de las proyecciones PCA, con la holgura
    # de un bloque de la reduccion
    b = np.array(malla.bounds).reshape(3, 2)
    holgura = FM.reducir(mask)[1] * H
    assert (b[:, 0] >= lo - holgura).all() and (b[:, 1] <= hi + holgura).all()

    p = pv.Plotter(off_screen=True, window_size=(160, 120))
    e = ER.EscenaRecorte(p, malla, lo, hi)
    c = g.centro(0.5)
    cubo, _f = V.extraer_cubo(mask, H, c, ejes, centro, 1.0)
    e.poner_cubo(c, cubo.shape[0] * H, texto="1 mm", render=False)
    e.poner_contenido(ER.malla_cubo_pca(cubo, H, c), render=False)
    e.poner_tramos(V.centros_por_tercios(proy), ["1", "2", "3"], render=False)
    for vista in ("iso", "long12", "long13", "transversal"):
        e.vista(vista)
    e.fijar_seccion(True)
    e.fijar_traslucido(False)
    for nombre in ("hueso", "caja", "aristas", "contenido", "tramos"):
        assert nombre in p.renderer.actors, nombre
    img = p.screenshot(return_img=True)
    p.close()
    assert img.std() > 1.0          # algo se dibujo, no un lienzo vacio
