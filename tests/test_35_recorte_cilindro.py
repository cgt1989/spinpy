"""
Bloque 35 - Recorte cilindrico del VOI: geometria, dominio y medida.

QUE SE VERIFICA Y POR QUE
-------------------------
El cilindro se recorta con la misma regla que el cubo, pero su resultado es una
caja envolvente MAS un dominio. Esa pareja es la fuente del error mas facil de
cometer sin darse cuenta: medir sobre la caja entera cuenta el exterior del
cilindro como poro, y el BV/TV sale ~ pi/4 de lo real sin ningun aviso. Por eso
se comprueba contra geometrias de respuesta conocida:

  35.1  `dimensionar_cilindro`: la altura es un numero IMPAR de voxeles, el
        dominio tiene el numero de voxeles que cuenta una rejilla independiente
        y el volumen discreto se acerca al analitico pi r^2 h.
  35.2  `extraer_cilindro` sobre un bloque macizo: el dominio es el cilindro, el
        BV/TV es 1 (no pi/4) y la caja tiene el orden de ejes documentado, para
        los tres ejes posibles. Un bloque sin simetria detecta ejes cruzados.
  35.3  `cilindro_a_pca`: cada voxel de la caja, llevado al marco PCA y de
        vuelta a la pila, cae sobre un voxel con el mismo valor.
  35.4  `morfometria_dominio`: sobre laminas de espesor conocido, el BV/TV y el
        Tb.Th restringidos al dominio coinciden con la solucion continua, y la
        medida ingenua sobre la caja entera NO (la prueba discrimina).
  35.5  El dibujo: la malla del contenido del cilindro esta donde el cilindro y
        la escena se compone con un cilindro.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  dimensionar_cilindro   exacto (enteros); volumen discreto vs analitico 1 %
  extraer_cilindro       BV/TV del bloque macizo: exacto (1.0)
  cilindro_a_pca         0 voxeles de discrepancia
  morfometria_dominio    BV/TV 1 % y Tb.Th 6 % frente a la solucion continua.
                         El 6 % no es holgura gratuita: la pared curva corta las
                         laminas paralelas al eje y el area se subestima ~ 4 %
                         (medido), que es del orden del 8.5 % que ya acepta
                         marching cubes en este proyecto (correccion C4).
  malla del cilindro     1e-9 mm
"""
import numpy as np
import pytest

from spinpy import morphometry as M
from spinpy import voi as V

BLOQUE = "35 · recorte cilíndrico"
REF = "Definicion de `voi.extraer_cilindro` y solucion continua de laminas."

H = 0.02


def _lam(n=120, periodo=12, grosor=4):
    """Laminas normales a x: espesor `grosor` voxeles cada `periodo`."""
    x = np.arange(n)
    return ((x % periodo) < grosor)[:, None, None] * np.ones((n, n, n), bool)


# -- 35.1 --------------------------------------------------------------------

def test_dimensionar_cilindro(registro):
    d = V.dimensionar_cilindro(H, radio_mm=0.8, altura_mm=1.0)
    # radio 0.8 mm = 40 voxeles -> caja de 81; altura 1.0 mm = 50 -> 51 (impar)
    r = 0.8 / H
    n = 81
    u = np.arange(n) - 40
    esperado_plano = int(((u[:, None] ** 2 + u[None, :] ** 2) <= r * r).sum())
    ok = (d["lado_caja"] == 81 and d["alto_vox"] == 51
          and d["n_vox_dominio"] == esperado_plano * 51)
    rel = abs(d["volumen_real_mm3"] / d["volumen_analitico_mm3"] - 1.0)
    registro.anotar(BLOQUE, "dimensionar_cilindro: dominio y altura impar",
                    REF, esperado_plano * 51, d["n_vox_dominio"], "exacto",
                    "voxeles del dominio", ok)
    registro.anotar(BLOQUE, "volumen discreto vs analitico", REF,
                    d["volumen_analitico_mm3"], d["volumen_real_mm3"], "1 %",
                    f"desviacion {100 * rel:.3f} %", rel < 0.01)
    assert ok and rel < 0.01
    assert d["altura_real_mm"] == pytest.approx(51 * H)
    assert d["desviacion_rel_altura"] == pytest.approx(51 * H / 1.0 - 1.0)


def test_dimensionar_cilindro_rechaza_entradas_incoherentes():
    for malo in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            V.dimensionar_cilindro(H, malo, 1.0)
        with pytest.raises(ValueError):
            V.dimensionar_cilindro(H, 1.0, malo)
    with pytest.raises(ValueError):
        V.extraer_cilindro(np.ones((9, 9, 9), bool), H, np.zeros(3),
                           np.eye(3), np.zeros(3), 0.05, 0.05, eje=3)


# -- 35.2 --------------------------------------------------------------------

@pytest.mark.parametrize("eje", [0, 1, 2])
def test_bloque_macizo_da_bvtv_uno(registro, eje):
    """Un bloque macizo lleno el cilindro: BV/TV = 1 dentro del dominio, aunque
    la caja envolvente tenga pi/4 de su volumen ocupado."""
    n = 90
    mask = np.ones((n, n, n), bool)
    c = np.full(3, n * H / 2)
    r = V.extraer_cilindro(mask, H, c, np.eye(3), np.zeros(3), 0.6, 0.8,
                           eje=eje)
    vol, dom = r["volumen"], r["dominio"]
    ocupacion_caja = vol.mean()
    registro.anotar(BLOQUE, f"BV/TV de un bloque macizo (eje {eje})", REF, 1.0,
                    r["bvtv"], "exacto",
                    f"la caja entera daria {ocupacion_caja:.3f}",
                    r["bvtv"] == 1.0)
    assert r["bvtv"] == 1.0 and r["fuera"] == 0.0
    assert ocupacion_caja < 0.85                       # ~ pi/4, no 1
    assert vol.shape == (dom.shape[0], dom.shape[0], r["dim"]["alto_vox"])
    assert vol.sum() == dom.sum()                      # nada fuera del dominio


@pytest.mark.parametrize("eje", [0, 1, 2])
def test_ejes_de_la_caja(eje):
    """El eje del cilindro es la ULTIMA dimension de la caja. Se comprueba con
    un bloque que solo cubre un tramo del eje: el hueso ocupa la mitad de la
    altura del cilindro y solo si la dimension 2 es el eje sale 0.5."""
    n = 90
    mask = np.zeros((n, n, n), bool)
    sl = [slice(None)] * 3
    sl[eje] = slice(0, n // 2)
    mask[tuple(sl)] = True
    c = np.full(3, n * H / 2)
    r = V.extraer_cilindro(mask, H, c, np.eye(3), np.zeros(3), 0.5, 0.8,
                           eje=eje)
    perfil = r["volumen"].sum(axis=(0, 1)) / r["dominio"].sum(axis=(0, 1))
    assert perfil[0] == 1.0 or perfil[-1] == 1.0
    assert 0.0 in (perfil[0], perfil[-1])
    assert 0.4 < r["bvtv"] < 0.6


def test_cilindro_fuera_del_volumen_se_declara():
    mask = np.ones((60, 60, 60), bool)
    c = np.array([0.0, 0.6, 0.6])              # el eje del cilindro en el borde
    r = V.extraer_cilindro(mask, H, c, np.eye(3), np.zeros(3), 0.4, 0.6, eje=0)
    assert 0.3 < r["fuera"] < 0.7
    assert r["bvtv"] == pytest.approx(1.0 - r["fuera"], abs=1e-12)


# -- 35.3 --------------------------------------------------------------------

@pytest.mark.parametrize("eje", [0, 1, 2])
def test_cilindro_a_pca_cae_sobre_el_mismo_voxel(registro, eje):
    rng = np.random.default_rng(5)
    mask = rng.random((70, 64, 58)) < 0.4        # sin simetria ninguna
    c = np.array([0.6, 0.55, 0.5])
    r = V.extraer_cilindro(mask, H, c, np.eye(3), np.zeros(3), 0.3, 0.5,
                           eje=eje)
    vol, dom = r["volumen"], r["dominio"]
    idx = np.indices(vol.shape).reshape(3, -1).T
    mundo = V.cilindro_a_pca(idx, vol.shape, H, c, eje)
    vox = np.rint(mundo / H).astype(int)
    dentro = np.all((vox >= 0) & (vox < np.array(mask.shape)), axis=1)
    en_dom = dom[idx[:, 0], idx[:, 1], idx[:, 2]]
    sel = dentro & en_dom
    esperado = vol[idx[sel, 0], idx[sel, 1], idx[sel, 2]]
    leido = mask[vox[sel, 0], vox[sel, 1], vox[sel, 2]]
    fallos = int((esperado != leido).sum())
    registro.anotar(BLOQUE, f"cilindro_a_pca reproduce extraer_cilindro "
                    f"(eje {eje})", REF, 0, fallos, "0 voxeles",
                    f"{int(sel.sum())} voxeles comparados", fallos == 0)
    assert sel.sum() > 1000 and fallos == 0


# -- 35.4 --------------------------------------------------------------------

def _continuo(R, alto, c, periodo=12, grosor=4, n=120, h=H):
    """BV, BS y Tb.Th de las laminas en un cilindro continuo de eje paralelo a
    las laminas (la pared las corta)."""
    def cuerda(x):
        return 2 * np.sqrt(np.clip(R * R - (x - c) ** 2, 0, None))
    BV = BS = 0.0
    for k in range(n // periodo):
        x1, x2 = (periodo * k - 0.5) * h, (periodo * k + grosor - 0.5) * h
        xs = np.linspace(x1, x2, 20001)
        BV += np.trapezoid(cuerda(xs), xs) * alto
        BS += sum(cuerda(x) * alto for x in (x1, x2) if abs(x - c) < R)
    return BV, BS, 2 * BV / BS, BV / (np.pi * R * R * alto)


def test_morfometria_dominio_frente_a_la_solucion_continua(registro):
    mask = _lam()
    c = np.full(3, 120 * H / 2)
    R, alto = 0.8, 51 * H
    BV, BS, tbth, bvtv = _continuo(R, alto, c[0])
    r = V.extraer_cilindro(mask, H, c, np.eye(3), np.zeros(3), R, 51 * H,
                           eje=1)                  # eje paralelo a las laminas
    m = M.morfometria_dominio(r["volumen"], r["dominio"], H)
    e_bv = abs(m["BVTV"] / bvtv - 1.0)
    e_th = abs(m["TbTh"] / tbth - 1.0)
    ingenuo = M.morfometria(r["volumen"], H, do_mil=False)
    registro.anotar(BLOQUE, "BV/TV del dominio frente a la solucion continua",
                    REF, bvtv, m["BVTV"], "1 %", f"error {100 * e_bv:.2f} %",
                    e_bv < 0.01)
    registro.anotar(BLOQUE, "Tb.Th del dominio frente a la solucion continua",
                    REF, tbth, m["TbTh"], "6 %", f"error {100 * e_th:.2f} %",
                    e_th < 0.06)
    assert e_bv < 0.01 and e_th < 0.06
    # la medida sobre la caja entera esta MUY fuera: la prueba discrimina
    assert abs(ingenuo["BVTV"] / bvtv - 1.0) > 0.15
    assert m["dominio"] is True and not m["MIL_valid"]
    assert np.isnan(m["DA"]) and np.isnan(m["SMI"])


def test_morfometria_dominio_eje_perpendicular_a_las_laminas():
    """Con el eje normal a las laminas la pared no las corta de canto: el
    area es la de laminas circulares y el Tb.Th sale exacto."""
    mask = _lam()
    c = np.full(3, 120 * H / 2)
    r = V.extraer_cilindro(mask, H, c, np.eye(3), np.zeros(3), 0.8, 1.0,
                           eje=0)
    m = M.morfometria_dominio(r["volumen"], r["dominio"], H)
    assert m["TbTh"] == pytest.approx(4 * H, rel=0.01)


def test_morfometria_dominio_con_dominio_completo_es_la_de_la_caja():
    """Si el dominio es toda la caja, la medida es la de `morfometria`."""
    mask = _lam(n=60)
    dom = np.ones(mask.shape, bool)
    m = M.morfometria_dominio(mask, dom, H)
    ref = M.morfometria(mask, H, do_mil=False)
    for k in ("BVTV", "BS", "BSBV", "TbTh", "TbSp", "TbN"):
        assert m[k] == pytest.approx(ref[k], rel=1e-12), k


def test_morfometria_dominio_casos_limite():
    dom = np.zeros((9, 9, 9), bool)
    dom[2:7, 2:7, 2:7] = True
    vacio = np.zeros(dom.shape, bool)
    m = M.morfometria_dominio(vacio, dom, H)
    assert m["BVTV"] == 0.0 and np.isnan(m["TbTh"])
    with pytest.raises(ValueError):
        M.morfometria_dominio(vacio, np.zeros_like(dom), H)
    with pytest.raises(ValueError):
        M.morfometria_dominio(vacio, np.ones((3, 3, 3), bool), H)


# -- 35.5 --------------------------------------------------------------------

@pytest.mark.parametrize("eje", [0, 1, 2])
def test_malla_del_cilindro_esta_donde_el_cilindro(registro, eje):
    """Bloque macizo y asimetrico dentro del cilindro: la caja de la malla, en
    el marco PCA, es la del cilindro recortado por el bloque."""
    from spinpy.escena_recorte import malla_cilindro_pca
    mask = np.ones((80, 80, 80), bool)
    c = np.array([0.8, 0.8, 0.8])
    r = V.extraer_cilindro(mask, H, c, np.eye(3), np.zeros(3), 0.3, 0.5,
                           eje=eje)
    m = malla_cilindro_pca(r["volumen"], r["dominio"], H, c, eje)
    semi = np.full(3, r["dim"]["lado_caja"] * H / 2)
    semi[eje] = r["dim"]["alto_vox"] * H / 2
    esperado = np.stack([c - semi, c + semi], axis=1)
    obtenido = np.array(m.bounds).reshape(3, 2)
    # el isovalor 0.5 cae medio voxel fuera del ultimo voxel macizo; en el
    # plano la pared curva toca la caja solo en cuatro puntos, y en la
    # direccion del eje la tapa es plana
    err_eje = abs(obtenido[eje] - esperado[eje]).max()
    registro.anotar(BLOQUE, f"malla_cilindro_pca: tapas (eje {eje})", REF, 0.0,
                    float(err_eje), "1e-9 mm", "extension de la malla en el "
                    "eje del cilindro", err_eje < 1e-9)
    assert err_eje < 1e-9
    for k in range(3):
        if k != eje:
            assert abs(obtenido[k] - esperado[k]).max() < 1e-9


def test_malla_del_cilindro_vacio_no_falla():
    from spinpy.escena_recorte import malla_cilindro_pca
    dom = V.dominio_cilindro(9, 9, 4.0)
    assert malla_cilindro_pca(np.zeros_like(dom), dom, 0.1, np.zeros(3),
                              0).n_points == 0


def test_escena_con_cilindro(render_disponible):
    if not render_disponible:
        pytest.skip("VTK no puede dibujar fuera de pantalla en este equipo")
    import pyvista as pv
    from spinpy import escena_recorte as ER
    mask = _lam(n=100)
    c = np.full(3, 100 * H / 2)
    ejes, centro = np.eye(3), np.zeros(3)
    malla = ER.malla_hueso_pca((mask.astype(np.float32), 1), H, ejes, centro)
    lo, hi = np.zeros(3), np.full(3, 100 * H)
    p = pv.Plotter(off_screen=True, window_size=(160, 120))
    e = ER.EscenaRecorte(p, malla, lo, hi)
    r = V.extraer_cilindro(mask, H, c, ejes, centro, 0.5, 0.8, eje=1)
    e.poner_cilindro(c, 0.5, r["dim"]["altura_real_mm"], 1, texto="r 0.5",
                     render=False)
    e.poner_contenido(ER.malla_cilindro_pca(r["volumen"], r["dominio"], H, c,
                                            1), render=False)
    for nombre in ("caja", "aristas", "contenido"):
        assert nombre in p.renderer.actors, nombre
    # cambiar de forma reemplaza la caja, no la acumula
    e.poner_cubo(c, 0.6, render=False)
    e.poner_cilindro(c, 0.4, 0.6, 2, render=False)
    e.fijar_seccion(True)
    img = p.screenshot(return_img=True)
    p.close()
    assert img.std() > 1.0


# -- 35.6 --------------------------------------------------------------------
# Lo que la ventana usa para despachar entre cubo y cilindro. Es logica pura:
# se prueba aqui para no depender de Qt.

def test_recortar_despacha_por_forma():
    rng = np.random.default_rng(2)
    mask = rng.random((70, 70, 70)) < 0.5
    c = np.full(3, 70 * H / 2)
    cu = V.recortar(mask, H, np.eye(3), np.zeros(3),
                    {"forma": "cubo", "centro_pca": c, "lado_mm": 0.6})
    ref, fuera = V.extraer_cubo(mask, H, c, np.eye(3), np.zeros(3), 0.6)
    assert cu["dominio"] is None and np.array_equal(cu["cubo"], ref)
    assert cu["bvtv"] == float(ref.mean()) and cu["fuera"] == fuera
    cil = V.recortar(mask, H, np.eye(3), np.zeros(3),
                     {"forma": "cilindro", "centro_pca": c, "radio_mm": 0.3,
                      "altura_mm": 0.6, "eje": 1})
    ref = V.extraer_cilindro(mask, H, c, np.eye(3), np.zeros(3), 0.3, 0.6,
                             eje=1)
    assert np.array_equal(cil["cubo"], ref["volumen"])
    assert np.array_equal(cil["dominio"], ref["dominio"])
    assert cil["bvtv"] == ref["bvtv"]
    with pytest.raises(ValueError):
        V.recortar(mask, H, np.eye(3), np.zeros(3),
                   {"forma": "esfera", "centro_pca": c})


def test_semiejes_de_los_cortes():
    cil = {"forma": "cilindro", "radio_mm": 0.3, "altura_mm": 1.0, "eje": 2}
    # planos que contienen el eje: rectangulo R x H/2; perpendicular: circulo
    assert V.semiejes_corte(cil, 0, 2) == ((0.3, 0.5), False)
    assert V.semiejes_corte(cil, 2, 1) == ((0.5, 0.3), False)
    assert V.semiejes_corte(cil, 0, 1) == ((0.3, 0.3), True)
    cub = {"forma": "cubo", "lado_mm": 0.8}
    assert V.semiejes_corte(cub, 0, 1) == ((0.4, 0.4), False)


def test_misma_region_distingue_forma_y_medidas():
    a = {"forma": "cilindro", "radio_mm": 0.3, "altura_mm": 1.0, "eje": 0}
    assert V.misma_region(a, dict(a))
    assert not V.misma_region(a, {**a, "eje": 1})
    assert not V.misma_region(a, {**a, "radio_mm": 0.31})
    assert not V.misma_region(a, {"forma": "cubo", "lado_mm": 1.0})


def test_registro_del_recorte_declara_el_cilindro():
    from spinpy import figura_metodo as FM
    pila = {"origen": "/x", "patron": "*.tif", "umbral": None,
            "tam_voxel_mm": 0.02}
    reg = {"forma": "cilindro", "radio_mm": 0.3, "altura_mm": 1.0, "eje": 2}
    rec = FM.registro_recorte(pila, np.eye(3), np.zeros(3), np.zeros(3), None,
                              region=reg)
    assert rec["forma"] == "cilindro" and rec["lado_mm"] is None
    assert (rec["radio_mm"], rec["altura_mm"], rec["eje"]) == (0.3, 1.0, 2)
    assert FM.rehacer_contexto(rec) is None      # la figura 0 es del cubo
    # el cubo conserva su registro de siempre, sin claves nuevas
    cu = FM.registro_recorte(pila, np.eye(3), np.zeros(3), np.zeros(3), 1.0)
    assert "forma" not in cu and cu["lado_mm"] == 1.0
