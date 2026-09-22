"""
Bloque 8 — Entrada por pila TIFF: escala fisica, convencion de ejes y umbral.

QUE SE VERIFICA Y POR QUE ESTO MERECE PRUEBAS PROPIAS
------------------------------------------------------
Un .vtk trae el tamano de voxel dentro. Una pila de imagenes de micro-CT no
trae NADA: la escala hay que sacarla del `*_rec.log` del escaner o pedirsela a
quien carga. Es la clase de error que no da sintomas — con un tamano de voxel
equivocado el volumen se ve idéntico, la segmentacion es la misma y BV/TV
tambien (es adimensional), pero Tb.Th, Tb.Sp, TV, Conn.D y la rigidez aparente
salen escalados por un factor constante y nada avisa. Por eso se comprueba:

  8.1  el tamano de voxel se lee del log real de H4 y vale 51.489 um, que es
       el valor documentado del proyecto. Se comprueba ademas que se llega a el
       desde el .log, desde la carpeta de rebanadas y desde la carpeta padre,
       porque NRecon deja el log en cualquiera de los dos sitios.
  8.2  la conversion de unidades es exacta: 51.489 um, 0.051489 mm y
       0.0051489 cm tienen que dar EXACTAMENTE el mismo spacing interno.
  8.3  la convencion de ejes se conserva: dimension 1 = X, igual que
       `readVTKVOI` y que el resto de spinpy. Mezclarlas es el bug B1/B2 del
       historial. Se verifica con un volumen asimetrico en las tres
       dimensiones, que es el unico que lo detecta.
  8.4  el umbral automatico recupera la particion real cuando las dos
       poblaciones estan separadas (Otsu), y reconoce una pila ya binaria.
  8.5  sin escala y sin log se LANZA excepcion. No hay valor por defecto: es
       preferible parar a suponer 1 mm/voxel.
  8.6  los archivos que el escaner deja junto a la reconstruccion y que NO son
       rebanadas se apartan, y se apartan CON NOMBRE. Regresion: apuntar a
       `H4/` con `*.tif` fallaba con «H 1-4_rec_spr.tif: forma (1944, 936)
       distinta de (200, 936)», porque el lector tomaba como referencia la
       primera imagen del orden natural —una proyeccion `_pp1`— y abortaba al
       llegar a la siguiente. Filtrar por nombre no vale: `_prev_00000235`
       tambien lleva numero de rebanada.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  tamano de voxel   exacto a 1e-12 (es una lectura de texto, no una medida)
  unidades          exacto a 1e-15 relativo
  ejes              exacto (comparacion de formas, sin tolerancia)
  Otsu              BV/TV recuperado a 0.005 absoluto con dos gaussianas
                    separadas 5 sigma; con menos separacion Otsu no tiene por
                    que acertar y no se le exige.
"""

from pathlib import Path

import numpy as np
import pytest

from spinpy.voi import UNIDADES_MM, leer_pila_tiff, tam_voxel_desde_log

BLOQUE = "8 · entrada por pila TIFF"
REF = ("SkyScan/NRecon `*_rec.log`, campo «Image Pixel Size (um)»; "
       "Otsu N. IEEE Trans Syst Man Cybern 1979;9:62.")

# Valor documentado del proyecto para H1-H4 (log SkyScan *_rec.log y cabecera de la app)
TAM_H4_UM = 51.489
RAIZ = Path(__file__).resolve().parents[2]


def _escribir_pila(carpeta, vol):
    from skimage import io as skio
    carpeta.mkdir(parents=True, exist_ok=True)
    for k in range(vol.shape[2]):
        # .T deshace la trasposicion que aplica el lector: se guarda como
        # (fila, columna) = (y, x), que es lo que escribe cualquier escaner.
        skio.imsave(carpeta / f"sl_{k:05d}.tif", vol[:, :, k].T,
                    check_contrast=False)


@pytest.mark.parametrize("desde", ["log", "carpeta_rebanadas", "carpeta_padre"])
def test_81_tam_voxel_desde_log_real(registro, desde):
    origen = {"log": RAIZ / "H4" / "H 1-4_rec.log",
              "carpeta_rebanadas": RAIZ / "H4" / "Segmentadas",
              "carpeta_padre": RAIZ / "H4"}[desde]
    if not Path(origen).exists():
        pytest.skip(f"no esta el dato real: {origen}")

    mm, info = tam_voxel_desde_log(origen)
    esperado = TAM_H4_UM * 1e-3
    ok = mm is not None and abs(mm - esperado) < 1e-12
    registro.anotar(BLOQUE, f"tamano de voxel de H4 desde {desde}", REF,
                    esperado, mm, "1e-12 mm", "lectura exacta del log", ok,
                    nota=(info.get("clave") or ""))
    assert ok, f"desde {desde} se leyo {mm} mm en vez de {esperado}"


def test_82_unidades_equivalentes(tmp_path, registro):
    vol = (np.arange(6 * 5 * 3).reshape(6, 5, 3) % 2).astype(np.uint8)
    _escribir_pila(tmp_path / "u", vol)

    spacings = []
    for val, uni in ((TAM_H4_UM, "um"), (TAM_H4_UM * 1e-3, "mm"),
                     (TAM_H4_UM * 1e-4, "cm"), (TAM_H4_UM * 1e3, "nm")):
        _bw, sp, _i = leer_pila_tiff(tmp_path / "u", tam_voxel=val, unidad=uni)
        spacings.append(float(sp[0]))

    esperado = TAM_H4_UM * 1e-3
    peor = max(abs(s - esperado) / esperado for s in spacings)
    ok = peor < 1e-15
    registro.anotar(BLOQUE, "um/mm/cm/nm dan el mismo spacing interno", REF,
                    0.0, peor, "1e-15 relativo",
                    "la unidad de entrada no puede cambiar el resultado", ok)
    assert ok, f"las unidades no coinciden: {spacings}"
    # y el diccionario de unidades no puede perder factores por un retoque
    assert UNIDADES_MM["um"] == 1e-3 and UNIDADES_MM["cm"] == 10.0


def test_83_convencion_de_ejes(tmp_path, registro):
    # Asimetrico en las tres dimensiones: es el unico volumen que distingue
    # (x,y,z) de (y,x,z) y de (z,y,x).
    nx, ny, nz = 7, 11, 5
    vol = np.zeros((nx, ny, nz), dtype=np.uint8)
    vol[0, :, :] = 1                      # una cara entera en x = 0
    vol[:, 0, :] = 1
    _escribir_pila(tmp_path / "e", vol)

    BW, _sp, _i = leer_pila_tiff(tmp_path / "e", tam_voxel=1.0, unidad="mm")
    ok = BW.shape == (nx, ny, nz) and np.array_equal(BW, vol > 0)
    registro.anotar(BLOQUE, "dimension 1 = X, volumen recuperado sin permutar",
                    REF, 0.0, 0.0 if ok else 1.0, "exacto",
                    "forma y contenido identicos al original", ok,
                    nota=f"{BW.shape} vs {(nx, ny, nz)}")
    assert BW.shape == (nx, ny, nz), f"forma {BW.shape}, esperada {(nx, ny, nz)}"
    assert np.array_equal(BW, vol > 0), "el volumen sale permutado o reflejado"


def test_84_umbral_automatico(tmp_path, registro):
    rng = np.random.default_rng(1)
    mask = rng.random((30, 30, 10)) < 0.40
    # dos gaussianas separadas 5 sigma: es el regimen en el que Otsu tiene que
    # acertar. Con menos separacion no se le exige nada.
    gris = np.where(mask, rng.normal(42000, 2500, mask.shape),
                    rng.normal(9500, 2500, mask.shape))
    _escribir_pila(tmp_path / "g", gris.clip(0, 65535).astype(np.uint16))
    # DOS poblaciones -> Otsu clasico. El defecto del lector es 3 clases
    # porque el caso real es una reconstruccion con aire, medula y hueso; aqui
    # se pide 2 explicitamente, que es lo correcto para este dato.
    _bw, _sp, info = leer_pila_tiff(tmp_path / "g", tam_voxel=1.0, unidad="um",
                                    clases_otsu=2)

    err = abs(info["BVTV"] - mask.mean())
    ok = err < 0.005 and info["umbral_metodo"].startswith("Otsu automatico, 2")
    registro.anotar(BLOQUE, "Otsu recupera BV/TV con poblaciones separadas",
                    REF, float(mask.mean()), info["BVTV"], "0.005 absoluto",
                    "dos gaussianas a 5 sigma", ok,
                    nota=f"umbral {info['umbral']:.0f}")
    assert ok, f"BV/TV {info['BVTV']:.4f} vs {mask.mean():.4f} ({info})"

    # y una pila ya binaria se reconoce como tal, sin aplicar Otsu
    _escribir_pila(tmp_path / "b", mask.astype(np.uint8))
    _bw, _sp, info_b = leer_pila_tiff(tmp_path / "b", tam_voxel=1.0, unidad="um")
    assert info_b["umbral_metodo"].startswith("ya binaria")
    assert abs(info_b["BVTV"] - mask.mean()) < 1e-12


def test_85_sin_escala_falla(tmp_path, registro):
    # Ni cuadrada ni de 3 columnas: `imsave` interpretaria una rebanada (3,3)
    # como RGB y fallaria al ESCRIBIR el dato de prueba, no al leerlo.
    vol = np.ones((6, 5, 2), dtype=np.uint8)
    _escribir_pila(tmp_path / "s", vol)
    with pytest.raises(ValueError) as exc:
        leer_pila_tiff(tmp_path / "s")           # sin tam_voxel y sin log
    ok = "tamano de voxel" in str(exc.value)
    registro.anotar(BLOQUE, "sin escala y sin log se para en vez de suponer",
                    REF, 1.0, 1.0 if ok else 0.0, "exacto",
                    "debe lanzar ValueError explicando que falta", ok)
    assert ok, f"mensaje poco claro: {exc.value}"

    with pytest.raises(ValueError):
        leer_pila_tiff(tmp_path / "s", tam_voxel=10, unidad="pulgadas")
    with pytest.raises(ValueError):
        leer_pila_tiff(tmp_path / "s", tam_voxel=0, unidad="um")


def test_86_aparta_lo_que_no_es_rebanada(tmp_path, registro):
    """Regresion del fallo real de H4: proyecciones y previsualizaciones."""
    d = tmp_path / "mezcla"
    serie = np.zeros((9, 7, 6), dtype=np.uint16)
    serie[2:5, 1:4, :] = 1000
    _escribir_pila(d, serie)

    # Lo que deja NRecon junto a la reconstruccion, con OTRAS formas. Se les
    # pone numero en el nombre a proposito: si el filtro fuese por nombre,
    # estos pasarian igual.
    from skimage import io as skio
    intrusos = {"pp1.tif": np.zeros((5, 12), np.uint16),
                "prev_00000235.tif": np.zeros((12, 12), np.uint16),
                "rec_spr.tif": np.zeros((20, 12), np.uint16)}
    for nombre, arr in intrusos.items():
        skio.imsave(d / nombre, arr, check_contrast=False)

    BW, _sp, info = leer_pila_tiff(d, tam_voxel=TAM_H4_UM, unidad="um")

    ok = (BW.shape == serie.shape
          and info["n_rebanadas"] == serie.shape[2]
          and len(info["apartados"]) == len(intrusos)
          and {a["nombre"] for a in info["apartados"]} == set(intrusos))
    registro.anotar(BLOQUE, "aparta proyecciones y previsualizaciones, con nombre",
                    REF, float(serie.shape[2]), float(info["n_rebanadas"]),
                    "exacto", "la serie dominante por forma es la buena", ok,
                    nota=f"apartados: {[a['nombre'] for a in info['apartados']]}")
    assert BW.shape == serie.shape, f"forma {BW.shape}, esperada {serie.shape}"
    assert len(info["apartados"]) == len(intrusos), info["apartados"]
    assert np.array_equal(BW, serie > 0), "el volumen no es el de la serie"

    # Y si NO hay serie dominante, se para en vez de adivinar.
    d2 = tmp_path / "sin_mayoria"
    d2.mkdir()
    for i, forma in enumerate([(9, 7), (12, 12), (20, 12), (5, 12)]):
        skio.imsave(d2 / f"a_{i:03d}.tif", np.zeros(forma, np.uint16),
                    check_contrast=False)
    with pytest.raises(ValueError) as exc:
        leer_pila_tiff(d2, tam_voxel=1.0, unidad="mm")
    assert "dominante" in str(exc.value), exc.value


def test_87_otsu_tres_clases(tmp_path, registro):
    """Con aire, medula y hueso, 2 clases cuenta la medula como hueso.

    Es el caso real de una reconstruccion de micro-CT y el motivo de que el
    defecto del lector sean 3 clases. Medido sobre H4 frente a la segmentacion
    manual (tercios 0.28 / 0.55 / 0.77): 2 clases da 0.61 / 0.79 / 0.94 y
    3 clases 0.35 / 0.59 / 0.81.
    """
    rng = np.random.default_rng(7)
    n = (40, 30, 12)
    # tres poblaciones bien separadas, con el aire dominando como en un escaner
    r = rng.random(n)
    clase = np.where(r < 0.70, 0, np.where(r < 0.88, 1, 2))   # aire/medula/hueso
    centros = {0: 2000.0, 1: 20000.0, 2: 45000.0}
    gris = np.zeros(n)
    for c, mu in centros.items():
        gris[clase == c] = rng.normal(mu, 1500, int((clase == c).sum()))
    _escribir_pila(tmp_path / "t3", gris.clip(0, 65535).astype(np.uint16))

    hueso = float((clase == 2).mean())
    _b2, _s2, i2 = leer_pila_tiff(tmp_path / "t3", tam_voxel=1.0, unidad="um",
                                  clases_otsu=2)
    _b3, _s3, i3 = leer_pila_tiff(tmp_path / "t3", tam_voxel=1.0, unidad="um",
                                  clases_otsu=3)

    e2, e3 = abs(i2["BVTV"] - hueso), abs(i3["BVTV"] - hueso)
    ok = e3 < 0.01 and e3 < e2
    registro.anotar(BLOQUE, "Otsu 3 clases separa hueso de medula", REF,
                    hueso, i3["BVTV"], "0.01 absoluto",
                    "y ademas mejor que 2 clases", ok,
                    nota=f"2 clases: {i2['BVTV']:.4f} (u={i2['umbral']:.0f}); "
                         f"3 clases: {i3['BVTV']:.4f} (u={i3['umbral']:.0f})")
    assert e3 < 0.01, f"3 clases dio BV/TV {i3['BVTV']:.4f} vs {hueso:.4f}"
    assert e3 < e2, "3 clases no mejora a 2 clases en un volumen de 3 poblaciones"
    assert i3["umbral"] > i2["umbral"], "el umbral alto debe ser mayor"
