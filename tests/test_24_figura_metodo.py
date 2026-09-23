"""
Bloque 24 — Figura 0 del informe: el metodo (`spinpy.figura_metodo`).

QUE SE VERIFICA
  La figura promete que ningun panel es un dibujo. Se comprueba cada promesa
  contra el dato del que sale:

  direcciones  las direcciones de onda que dibuja son BIT A BIT las que uso
               `grf.campo_grf` para generar la estructura (misma semilla,
               mismo orden de consumo del generador).
  umbral       el solido dibujado es exactamente {campo <= umbral}, en las
               dos familias: el umbral del panel es el del generador.
  cubo         las esquinas del cubo, llevadas al marco PCA, son ±(lado_vox/2)
               voxeles, con lado_vox el del cubo que devuelve
               `voi.extraer_cubo`: la caja que se ve es la que se recorto.
  corte        el poligono del cubo en una rebanada perpendicular a un eje del
               cubo es un cuadrado de ese lado, y fuera del cubo no hay
               poligono; con el cubo girado, cada vertice cae en la superficie.
  rebanada     sobre una pila TIFF sintetica: la rebanada gris es el archivo
               que pasa por el centro del cubo, la segmentada es esa rebanada
               de la mascara, y `rehacer_contexto` (la CLI) da lo mismo que el
               contexto tomado al recortar.
  pie          el pie describe solo lo que la figura lleva (con o sin pila,
               una o dos familias) y no deja marcadores sin rellenar.
  figura       se compone con dos familias y pila, con una y sin pila, en ES y
               EN; con mas filas la figura es mas alta (necesita VTK).
  tiempo       el termino de la figura 0 en `tiempos` reproduce las pasadas
               medidas (2026-09-23) dentro de ±25 %, que es el ruido medido
               entre pasadas (±30 %); marcarla nunca abarata el informe.

Tolerancias declaradas antes de medir: todo es exacto salvo el tiempo.
"""
from pathlib import Path

import numpy as np
import pytest

from spinpy import figura_metodo as FM
from spinpy import figuras as F
from spinpy import grf, informe, procedencia, tiempos as T, voi as V

BLOQUE = "24 Figura del metodo"
REF = "figura_metodo.py (figura 0 del informe)"


def _par_spin():
    return {"familia": "spinodoide", "resolucion": 24, "densidad": 0.3,
            "wave_number_pi": 15.0, "wave_number_rad": 15.0 * np.pi,
            "num_waves": 200, "thetas": [15.0, 15.0, 60.0],
            "R": grf.euler_R(0, 10, 20).tolist(), "esquema": "rechazo",
            "semilla": 11}


def _par_dual():
    return {"familia": "dual-lattice", "resolucion": 24, "densidad": 0.3,
            "celdas": 3.0, "estiramiento": [1.0, 1.0, 1.6],
            "irregularidad": 0.5, "R": np.eye(3).tolist(), "semilla": 5}


def test_direcciones_del_generador(registro):
    par = _par_spin()
    kw = procedencia.kwargs_generador(par)
    _G, dirs_gen, _f = grf.campo_grf(kw["resolution"], kw["wave_number"],
                                     kw["num_waves"], kw["thetas"], R=kw["R"],
                                     esquema=kw["esquema"], seed=kw["seed"])
    d = FM.datos_familia(par)
    ok = np.array_equal(d["dirs"], dirs_gen)
    registro.anotar(BLOQUE, "direcciones de onda = las del generador", REF,
                    None, None, "bit a bit", "misma semilla y orden", ok)
    assert ok


def test_umbral_reproduce_mascara(registro):
    s, dl = FM.datos_familia(_par_spin()), FM.datos_familia(_par_dual())
    ok_s = np.array_equal(s["campo"] <= s["umbral"], s["BW"])
    ok_d = np.array_equal(dl["campo"] <= dl["umbral"], dl["BW"])
    ok = ok_s and ok_d
    registro.anotar(BLOQUE, "solido = {campo <= umbral} en las dos familias",
                    REF, None, None, "exacto", "umbral del generador", ok,
                    nota=f"spinodoide {ok_s}, dual {ok_d}")
    assert ok


def _pila_elipsoide(n=(60, 56, 48), h=0.05, semilla=3):
    """Hueso sintetico: elipsoide girado, poroso, en una pila (x, y, z)."""
    rng = np.random.default_rng(semilla)
    idx = np.indices(n).reshape(3, -1).T * h
    c = np.array(n) * h / 2
    R = grf.euler_R(0, 20, 35)
    q = (idx - c) @ R
    dentro = ((q / np.array([1.3, 0.9, 0.7])) ** 2).sum(1) <= 1.0
    poros = rng.random(len(idx)) < 0.35
    return (dentro & ~poros).reshape(n)


def test_esquinas_son_el_recorte(registro):
    h = 0.05
    mask = _pila_elipsoide(h=h)
    ejes, centro, proy = V.marco_pca(mask, h)
    c = V.centros_por_tercios(proy)[1]
    lado = 0.8
    cubo, _f = V.extraer_cubo(mask, h, c, ejes, centro, lado)
    E = FM.esquinas_cubo(ejes, centro, c, lado, h)
    en_pca = (E - centro) @ ejes - c
    semi = cubo.shape[0] / 2.0 * h
    ok = bool(np.allclose(np.abs(en_pca), semi, atol=1e-9, rtol=0))
    registro.anotar(BLOQUE, "esquinas del cubo = caja de extraer_cubo", REF,
                    semi, float(np.abs(en_pca).max()), "1e-9 mm",
                    f"lado_vox {cubo.shape[0]}", ok)
    assert ok


def test_corte_del_cubo(registro):
    semi = 0.5
    E = FM.esquinas_cubo(np.eye(3), np.zeros(3), np.zeros(3),
                         2 * semi - 0.1, 0.1)       # lado_vox = 9 -> semi .45
    semi = 4.5 * 0.1
    P = FM.corte_cubo_plano_z(E, 0.1)
    area = 0.5 * abs(np.dot(P[:, 0], np.roll(P[:, 1], 1))
                     - np.dot(P[:, 1], np.roll(P[:, 0], 1)))
    fuera = FM.corte_cubo_plano_z(E, 5.0)
    # Cubo girado: cada vertice del corte esta en la superficie del cubo.
    R = grf.euler_R(15, 30, 40)
    Eg = FM.esquinas_cubo(R, np.zeros(3), np.zeros(3), 0.9, 0.1)
    Pg = FM.corte_cubo_plano_z(Eg, 0.05)
    loc = np.column_stack([Pg, np.full(len(Pg), 0.05)]) @ R
    en_sup = np.allclose(np.abs(loc).max(axis=1), semi, atol=1e-9)
    ok = (len(P) == 4 and abs(area - (2 * semi) ** 2) < 1e-12
          and len(fuera) == 0 and len(Pg) >= 3 and en_sup)
    registro.anotar(BLOQUE, "corte del cubo con la rebanada", REF,
                    (2 * semi) ** 2, float(area), "1e-12 mm²",
                    "cuadrado dentro, nada fuera, vertices en la superficie",
                    ok, nota=f"{len(Pg)} vertices con el cubo girado")
    assert ok


def _escribir_pila(carpeta, mask, rng):
    """Rebanadas de 16 bits: hueso ~30000, fondo ~5000, con ruido."""
    from PIL import Image
    carpeta.mkdir(parents=True, exist_ok=True)
    gris = np.where(mask, 30000, 5000) + rng.integers(0, 2000, mask.shape)
    gris = gris.astype(np.uint16)
    for k in range(mask.shape[2]):
        Image.fromarray(gris[:, :, k].T).save(carpeta / f"s_rec{k:05d}.tif")
    return gris


def test_rebanada_y_cli(registro, tmp_path):
    h = 0.05
    mask = _pila_elipsoide(h=h)
    rng = np.random.default_rng(1)
    gris = _escribir_pila(tmp_path / "pila", mask, rng)
    BW, sp, info = V.leer_pila_tiff(tmp_path / "pila", patron="s_rec*.tif",
                                    umbral=17500, tam_voxel=h, unidad="mm")
    ejes, centro, proy = V.marco_pca(BW, sp)
    c = V.centros_por_tercios(proy)[1]
    ctx = FM.contexto_pila(BW, sp, ejes, centro, c, 0.8,
                           origen=tmp_path / "pila", patron="s_rec*.tif")
    c_mm = c @ ejes.T + centro
    k = int(np.rint(c_mm[2] / h))
    rec = FM.registro_recorte({"origen": str(tmp_path / "pila"),
                               "patron": "s_rec*.tif", "umbral": 17500.0,
                               "tam_voxel_mm": h}, ejes, centro, c, 0.8)
    ctx2 = FM.rehacer_contexto(rec)
    ok = (ctx["k"] == k and ctx["gris"] is not None
          and np.array_equal(ctx["gris"], gris[:, :, k].astype(float))
          and np.array_equal(ctx["binaria"], BW[:, :, k])
          and np.array_equal(BW, mask)
          and ctx2 is not None and ctx2["k"] == ctx["k"]
          and np.array_equal(ctx2["gris"], ctx["gris"])
          and np.array_equal(ctx2["binaria"], ctx["binaria"])
          and np.allclose(ctx2["esquinas"], ctx["esquinas"], atol=1e-12))
    registro.anotar(BLOQUE, "rebanada del centro del cubo y rehacer desde "
                    "la CLI", REF, None, None, "exacto",
                    "gris = archivo, binaria = mascara, CLI = recorte", ok,
                    nota=f"k = {ctx['k']} de {ctx['n_rebanadas']}")
    assert ok


def _prep(fams=("spinodoide",), ctx=None, voi=True):
    ests = ([("voi", np.zeros((8, 8, 8), bool), np.full(3, 0.05))]
            if voi else [])
    return {"estructuras": ests, "generados": {f: None for f in fams},
            "contexto_voi": ctx, "lado_mm": 5.0, "opciones_3d": None}


def test_pie_describe_lo_que_hay(registro):
    ctx = {"k": 9, "n_rebanadas": 40, "paso": 3, "gris": np.zeros((2, 2)),
           "binaria": np.zeros((2, 2), bool)}
    casos = {
        "pila, dos familias": (_prep(("spinodoide", "dual-lattice"), ctx),
                               ["10 de 40", "3³", "φ₀", "4-N", "{d ≤ r}"],
                               ["ya recortado"]),
        "sin pila, spinodoide": (_prep(("spinodoide",)),
                                 ["ya recortado", "φ₀"], ["4-N", "rebanada"]),
        "sin VOI, dual": (_prep(("dual-lattice",), voi=False),
                          ["4-N"], ["Fila superior", "φ₀"]),
    }
    fallos = []
    for nombre, (prep, si, no) in casos.items():
        for idioma in ("es", "en"):
            t = informe.pie_metodo(prep, idioma)
            if idioma == "es":
                fallos += [f"{nombre}: falta {x}" for x in si if x not in t]
                fallos += [f"{nombre}: sobra {x}" for x in no if x in t]
            import re
            resto = re.findall(r"\{[a-z_]+\}", t)
            if resto:
                fallos.append(f"{nombre}/{idioma}: sin rellenar {resto}")
    ok = not fallos
    registro.anotar(BLOQUE, "el pie describe solo lo que hay", REF, None,
                    None, "exacto", "tres variantes, ES y EN", ok,
                    nota="; ".join(fallos) or "ok")
    assert ok


def test_opciones_y_tiempo(registro):
    o = F.opciones_3d(None)
    off = F.opciones_3d({"metodo": 0})
    medidos = {(1, False): (29.8, 37.8), (1, True): (44.2, 50.0),
               (2, True): (48.2, 65.7), (2, False): (46.3, 56.6)}
    malos = []
    for (n, pila), (a, b) in medidos.items():
        m = T.metodo(n, pila)
        if not (0.75 * a <= m <= 1.25 * b):
            malos.append(f"{n} fam, pila {pila}: {m:.1f} fuera de {a}-{b}")
    plan = dict(familias=["spinodoide"], vox_voi=97 ** 3, n_fit=96,
                n_estilos=1, n_vistas=1)
    base = T.total(T.estimar(plan))
    con = T.total(T.estimar(dict(plan, metodo=True, pila=True)))
    ok = (o["metodo"] is True and off["metodo"] is False and not malos
          and con > base)
    registro.anotar(BLOQUE, "opcion de figura y su tiempo", REF, None, None,
                    "±25 %", "cuatro casos medidos; marcarla encarece", ok,
                    nota="; ".join(malos) or f"+{con - base:.0f} s")
    assert ok


def test_figura_se_compone(registro, tmp_path, monkeypatch):
    pytest.importorskip("pyvista")
    monkeypatch.setattr(F, "DPI", 80)            # la composicion, no la nitidez
    h = 0.05
    mask = _pila_elipsoide(h=h)
    ejes, centro, proy = V.marco_pca(mask, h)
    c = V.centros_por_tercios(proy)[1]
    cubo, _f = V.extraer_cubo(mask, h, c, ejes, centro, 0.8)
    ctx = FM.contexto_pila(mask, h, ejes, centro, c, 0.8)
    fams = {"spinodoide": FM.datos_familia(_par_spin()),
            "dual-lattice": FM.datos_familia(_par_dual())}
    lado = cubo.shape[0] * h
    doc = {"morfometria_voi": {"BVTV": 0.3, "TbTh": 0.2, "TbSp": 0.4,
                               "BSBV": 10.0, "DA": 1.4},
           "morfometria_spin": {"BVTV": 0.31, "TbTh": 0.21, "TbSp": 0.41,
                                "BSBV": 9.5, "DA": 1.3},
           "morfometria_dual": None}
    try:
        img = FM.renders_metodo((cubo, np.full(3, h)), fams, lado, ctx,
                                lado=200)
    except Exception as e:                          # sin OpenGL
        pytest.skip(f"sin render: {e}")
    fallos = img.pop("_fallos")
    from PIL import Image
    tres = FM.componer_metodo(img, fams, doc, tmp_path / "tres", "es", ctx,
                              (cubo, np.full(3, h)), lado)
    dos = FM.componer_metodo(img, {"spinodoide": fams["spinodoide"]}, doc,
                             tmp_path / "dos", "en", None,
                             (cubo, np.full(3, h)), lado)
    alto3 = Image.open(tres[0]).size[1]
    alto2 = Image.open(dos[0]).size[1]
    ok = (not fallos and len(tres) == 2 and len(dos) == 2
          and all(Path(r).exists() for r in tres + dos) and alto3 > alto2)
    registro.anotar(BLOQUE, "la figura se compone (3 y 2 filas)", REF, None,
                    None, "exacto", "PNG y PDF; mas filas, mas alta", ok,
                    nota=f"alto {alto3} / {alto2} px; fallos {fallos}")
    assert ok
