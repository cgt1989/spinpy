"""
Bloque 19 — Informe para publicacion (`spinpy.informe`).

QUE SE VERIFICA
  Que el informe no deja pasar como citable ninguna de las situaciones que los
  estudios del proyecto midieron como engañosas, que no marca como problema lo
  que no lo es, y que el paquete de reproduccion detecta de verdad una
  estructura distinta.

  limpio        una sesion sin problemas: E_app y tensor citables; la carga de
                Pistoia siempre con reservas; el maximo de von Mises nunca.
  von Mises     los tres sabores del pico se gradan por separado: maximo no
                citable, percentil 99 sobre todo el tejido con reservas,
                percentil 99 de la capa superficial citable, y este ultimo baja
                a reservas si la capa tiene menos de `N_SUPERFICIE_MIN`
                elementos.
  umbrales      cada comprobacion salta en su umbral declarado (n < 40,
                residuo > 1e-6, desconexion 1 % / 5 %, K < 5, thetas
                degenerados, DA < 1.10, BV/TV >= 0.83, numero de onda sin pi).
  textos        cada motivo tiene texto en los dos idiomas y se rellena; el
                metodo lleva las dos lecturas del numero de onda, los DOI y
                los tres ORCID.
  autoria       ORCID y version coinciden con CITATION.cff.
  reproduccion  la huella coincide con `generar_mascara` llamado a mano, la
                verificacion pasa, y falla al cambiar la semilla o quitar el pi.
  carpeta       todo en una carpeta: Markdown y PDF en los dos idiomas, las
                figuras de datos en PNG de al menos 1500 px y en PDF, los
                cortes ortogonales y las distribuciones de espesor y poro de
                cada estructura, y cada figura enlazada desde el Markdown
                existe. El render 3D no se
                prueba aqui (necesita OpenGL); lo cubre la prueba con un VOI.
  chequeo       la lista de chequeo de reporte (Bouxsein 2010, Erdemir 2012)
                cambia de estado cuando cambia el dato, da los mismos estados
                en ES y EN, cita sus guias y va al final del informe.
  README        la insignia y el texto dicen tantos bloques como archivos
                tests/test_NN_*.py hay.
"""
import copy
import json
import re
from pathlib import Path

import numpy as np
import pytest

import spinpy
from spinpy import informe, lista_chequeo, procedencia
from spinpy.grf import generar_mascara

BLOQUE = "19 Informe para publicacion"
REF = "informe.py (Estudio_Convergencia, F10, Estudio_MIL, Estudio_Familias)"
RAIZ = Path(__file__).resolve().parents[1]


def _par():
    return {"familia": "spinodoide", "resolucion": 24, "densidad": 0.3,
            "wave_number_pi": 15.0, "wave_number_rad": 15.0 * np.pi,
            "num_waves": 200, "thetas": [15.0, 15.0, 60.0],
            "thetas_degenerados": False, "R": np.eye(3).tolist(),
            "esquema": "rechazo", "semilla": 11}


def _doc():
    par = _par()
    m = {"BVTV": 0.30, "TbTh": 0.195, "DA": 1.50, "BSBV": 10.3}
    return {
        "procedencia": procedencia.bloque("spinodoide", semilla=11,
                                          wave_number_pi=15.0),
        "voi": {"nombre": "sintetico", "ruta": None, "sha256": None,
                "forma": [97, 97, 97], "spacing_mm": [0.051489] * 3},
        "morfometria_voi": dict(m),
        "morfometria_spin": dict(m, BVTV=0.31),
        "morfometria_dual": None,
        "resultados": {
            "morfometria": {"resolucion": 97, "extra": False},
            "ajuste": {"procedencia": procedencia.desde_parametros(par),
                       "parametros": par, "error": 0.05,
                       "diagnostico": {"bicontinuo": True,
                                       "poro_conexo_spin": 0.99},
                       "incertidumbre": {"K": 5,
                                         "error": {"media": 0.06, "sd": 0.01},
                                         "suelo_autoconsistente": {"media": 0.03}},
                       "objetivo": {"pesos": None, "peso_mecanico": 0.0}},
            "elastico": {"resolucion": 48, "E_s_Pa": 2e10, "nu_s": 0.3,
                         "por_estructura": {
                             "voi": {"Ex": 1e9, "Ey": 1e9, "Ez": 2e9,
                                     "solver": "AMG (agregacion suavizada) + CG",
                                     "residuo_rel": 1e-8},
                             "spin": {"Ex": 1e8, "Ey": 1e8, "Ez": 2e8,
                                      "solver": "AMG (agregacion suavizada) + CG",
                                      "residuo_rel": 1e-8}}},
            "resistencia": {"resolucion": 48, "apoyo": "libre", "E_s_Pa": 2e10,
                            "nu_s": 0.3, "por_estructura": {
                                "voi": {"ejes": {"Z": {
                                    "E_app": 1e9, "sigma_fallo": 5e6,
                                    "residuo_rel": 1e-8, "frac_portante": 0.999,
                                    "vm_max_fallo": 1e8}}}}},
        },
    }


def _busca(items, magnitud, estructura="voi"):
    for it in items:
        if it["magnitud"] == magnitud and it["estructura"] == estructura:
            return it
    raise AssertionError(f"sin item {magnitud}/{estructura}")


def test_sesion_limpia(registro):
    items = informe.comprobar(_doc())
    E = _busca(items, "E_app")
    ela = _busca(items, "elastico", "spinodoide")
    fallo = _busca(items, "fallo")
    vm = _busca(items, "vm_max")
    aj = _busca(items, "error_ajuste", "spinodoide")
    ok = (E["estado"] == informe.CITABLE and ela["estado"] == informe.CITABLE
          and aj["estado"] == informe.CITABLE
          and fallo["motivos"] == ["pistoia_calibracion"]
          and vm["estado"] == informe.NO_CITABLE
          and "vm_maximo" in vm["motivos"])
    registro.anotar(BLOQUE, "sesion limpia: sin falsos avisos", REF, None,
                    None, "exacto", "citable salvo Pistoia (reservas) y "
                    "maximo de von Mises (no citable)", ok,
                    nota=json.dumps(informe.resumen(items)))
    assert ok


def _con(ruta, valor):
    d = _doc()
    nodo = d
    for k in ruta[:-1]:
        nodo = nodo[k]
    nodo[ruta[-1]] = valor
    return d


CASOS = [
    ("n = 32 < 40", ("resultados", "elastico", "resolucion"), 32,
     "elastico", "voi", "resolucion_mecanica", informe.RESERVAS),
    ("n = 40 no avisa", ("resultados", "elastico", "resolucion"), 40,
     "elastico", "voi", None, informe.CITABLE),
    ("residuo 1e-5", ("resultados", "elastico", "por_estructura", "voi",
                      "residuo_rel"), 1e-5,
     "elastico", "voi", "residuo", informe.NO_CITABLE),
    ("desconexion 3 %", ("resultados", "resistencia", "por_estructura", "voi",
                         "ejes", "Z", "frac_portante"), 0.97,
     "E_app", "voi", "desconexion", informe.RESERVAS),
    ("desconexion 10 %", ("resultados", "resistencia", "por_estructura", "voi",
                          "ejes", "Z", "frac_portante"), 0.90,
     "E_app", "voi", "desconexion_alta", informe.NO_CITABLE),
    ("K = 1", ("resultados", "ajuste", "incertidumbre", "K"), 1,
     "error_ajuste", "spinodoide", "sin_replicas", informe.RESERVAS),
    ("thetas degenerados", ("resultados", "ajuste", "parametros",
                            "thetas_degenerados"), True,
     "parametros", "spinodoide", "thetas_degenerados", informe.RESERVAS),
    ("numero de onda sin pi", ("resultados", "ajuste", "parametros",
                               "wave_number_rad"), 15.0,
     "parametros", "spinodoide", "onda_inconsistente", informe.NO_CITABLE),
    ("VOI no trabecular", ("morfometria_voi", "BVTV"), 0.90,
     "morfometria", "voi", "voi_no_trabecular", informe.RESERVAS),
]


@pytest.mark.parametrize("nombre,ruta,valor,mag,est,motivo,estado", CASOS,
                         ids=[c[0] for c in CASOS])
def test_umbrales(registro, nombre, ruta, valor, mag, est, motivo, estado):
    it = _busca(informe.comprobar(_con(ruta, valor)), mag, est)
    ok = (it["estado"] == estado and
          (motivo in it["motivos"] if motivo else not it["motivos"]))
    registro.anotar(BLOQUE, f"umbral: {nombre}", REF, None, None, "exacto",
                    f"{motivo or 'sin motivo'} -> {estado}", ok,
                    nota=f"{it['estado']} {it['motivos']}")
    assert ok


def test_da_bajo_suelo(registro):
    items = informe.comprobar(_con(("morfometria_voi", "DA"), 1.05))
    it = _busca(items, "da")
    ok = it["motivos"] == ["da_sin_eje"] and not any(
        x["magnitud"] == "da" for x in informe.comprobar(_doc()))
    registro.anotar(BLOQUE, "DA 1.05 sin eje; DA 1.50 no avisa", REF, 1.10,
                    1.05, "umbral avisos.DA_MIN_EJE", "da_sin_eje", ok)
    assert ok


def test_textos_bilingues(registro):
    faltan = sorted(set(informe.GRAVEDAD) ^ set(informe.MOTIVOS))
    for c in informe.MOTIVOS:
        for idioma in ("es", "en"):
            informe.texto_motivo(c, {}, idioma)        # no revienta sin datos
    doc = _doc()
    es, _r = informe.parrafos_metodos(doc, "es")
    en, refs = informe.parrafos_metodos(doc, "en")
    es, en = " ".join(es), " ".join(en)
    items = informe.comprobar(doc)
    paq = informe.paquete_reproduccion(doc, regenerar=False)
    md = informe.informe_markdown(doc, items, paq, "en")
    ok = (not faltan
          and "15π (47.1239 rad)" in en and "15π (47,1239 rad)" in es
          and "10.1038/s41524-020-0341-6" in md
          and all(a["orcid"] in md for a in informe.AUTORES)
          and "[PENDIENTE]" in md
          and "wave:" not in md and "'wave'" not in md
          and "kumar2020" in refs and "pistoia2002" in refs)
    registro.anotar(BLOQUE, "textos ES/EN, dos lecturas de beta, DOI, ORCID",
                    REF, None, None, "exacto", "presentes", ok,
                    nota=f"motivos sin texto: {faltan}")
    assert ok


def test_autoria_coincide_con_citation(registro):
    cff = (RAIZ / "CITATION.cff").read_text(encoding="utf-8")
    orcid_cff = set(re.findall(r"^\s*orcid:\s*'https://orcid.org/([\dX-]+)'",
                               cff, re.M))
    version_cff = re.search(r"^version:\s*(\S+)", cff, re.M).group(1)
    orcid_mod = {a["orcid"] for a in informe.AUTORES}
    ok = orcid_cff == orcid_mod and version_cff == spinpy.__version__
    registro.anotar(BLOQUE, "ORCID y version = CITATION.cff", "CITATION.cff",
                    None, None, "igual", "una sola autoria", ok,
                    nota=f"{sorted(orcid_cff)} v{version_cff}")
    assert ok


def test_reproduccion_bit_a_bit(registro, tmp_path):
    doc = _doc()
    r = informe.escribir_informe(doc, tmp_path, figuras=False,
                                 render_3d=False, pdf=False)
    paq = json.loads((tmp_path / "reproduccion.json").read_text(
        encoding="utf-8"))
    ref, _c, _i = generar_mascara(resolution=24, wave_number=15 * np.pi,
                                  num_waves=200, thetas=[15, 15, 60], rho=0.3,
                                  R=np.eye(3), esquema="rechazo", seed=11)
    huella_directa = informe.huella_mascara(ref)
    v_ok = informe.verificar_reproduccion(paq)

    otra = copy.deepcopy(paq)
    otra["ajustes"]["spinodoide"]["parametros"]["semilla"] = 12
    v_semilla = informe.verificar_reproduccion(otra)

    sin_pi = copy.deepcopy(paq)
    sin_pi["ajustes"]["spinodoide"]["parametros"]["wave_number_rad"] = 15.0
    v_pi = informe.verificar_reproduccion(sin_pi)

    vacio = informe.verificar_reproduccion({"ajustes": {}})

    ok = (all(Path(f).exists() for f in r["archivos"])
          and paq["ajustes"]["spinodoide"]["sha256_mascara"] == huella_directa
          and v_ok["ok"] and not v_semilla["ok"] and not v_pi["ok"]
          and "onda_inconsistente" in v_pi["ajustes"]["spinodoide"]["motivo"]
          and not vacio["ok"])
    registro.anotar(BLOQUE, "paquete de reproduccion", REF, None, None,
                    "bit a bit", "pasa igual; falla con otra semilla, sin pi "
                    "o sin nada que comprobar", ok,
                    nota=huella_directa[:16])
    assert ok


def _ancho_png(ruta):
    with open(ruta, "rb") as f:
        return int.from_bytes(f.read(24)[16:20], "big")


def test_carpeta_completa(registro, tmp_path):
    r = informe.escribir_informe(_doc(), tmp_path, render_3d=False)
    A = informe.ARCHIVOS
    hay = {k: (tmp_path / v).exists() for k, v in A.items()}
    try:
        import reportlab  # noqa: F401
        con_pdf = True
    except ImportError:
        con_pdf = False
    pngs = sorted((tmp_path / "figuras").glob("*.png"))
    pdfs_fig = sorted((tmp_path / "figuras").glob("*.pdf"))
    enlaces = {}
    for idioma in ("es", "en"):
        md = (tmp_path / A["md_" + idioma]).read_text(encoding="utf-8")
        enlaces[idioma] = re.findall(r"^!\[[^\]]*\]\(([^)]+)\)$", md, re.M)
    ok = (hay["md_es"] and hay["md_en"] and hay["paquete"] and hay["sesion"]
          and ((hay["pdf_es"] and hay["pdf_en"]
                and (tmp_path / A["pdf_es"]).read_bytes()[:4] == b"%PDF")
               if con_pdf else "pdf_sin_reportlab" in r["avisos"])
          # fig1, fig2 y fig7 (distribuciones) en PNG y PDF por idioma, y
          # fig5 (cortes) en PNG. fig6 no sale: este documento no trae tensor
          # de fabrica ni tensor C.
          and len(pngs) == 8 and len(pdfs_fig) == 6
          and all(_ancho_png(p) >= 1500 for p in pngs)
          and all(len(v) == 4 and all((tmp_path / e).exists() for e in v)
                  for v in enlaces.values())
          and enlaces["es"] != enlaces["en"])
    registro.anotar(BLOQUE, "carpeta del informe completa", REF, None, None,
                    "exacto", "MD+PDF ES/EN, figuras PNG>=1500 px y PDF, "
                    "enlaces validos", ok,
                    nota=f"{len(pngs)} png, anchos "
                         f"{[_ancho_png(p) for p in pngs]}, pdf={con_pdf}")
    assert ok


def test_sabores_del_pico_de_von_mises(registro):
    """Maximo, p99 global y p99 de la capa NO valen lo mismo.

    Lo que el estudio de convergencia valida es el p99 de la capa superficial
    (-0.3 % frente a la solucion cerrada). El maximo no converge y el p99
    global esta dominado por el material a granel. El informe tiene que
    separarlos, porque las tres claves se parecen y salen del mismo calculo.
    """
    doc = _doc()
    z = doc["resultados"]["resistencia"]["por_estructura"]["voi"]["ejes"]["Z"]
    z.update({"vm_p99": 4e7, "vm_p99_fallo": 8e7,
              "vm_p99_superficie": 5e7, "vm_p99_superficie_fallo": 9e7,
              "vm_n_superficie": 3000})
    items = informe.comprobar(doc)
    mx = _busca(items, "vm_max")
    glob = _busca(items, "vm_p99")
    sup = _busca(items, "vm_p99_superficie")
    ok = (mx["estado"] == informe.NO_CITABLE
          and glob["estado"] == informe.RESERVAS
          and "vm_p99_global" in glob["motivos"]
          and sup["estado"] == informe.CITABLE)
    registro.anotar(BLOQUE, "von Mises: maximo / p99 global / p99 de la capa",
                    REF, None, None, "exacto",
                    "no citable / reservas / citable", ok,
                    nota=f"max={mx['estado']} global={glob['estado']} "
                         f"capa={sup['estado']}")
    assert ok, (mx["estado"], glob["estado"], sup["estado"])

    # Muestra pequena: el mismo numero deja de ser citable. El umbral es
    # NUESTRO criterio, no una medida, y por eso se comprueba en su valor
    # declarado y justo por debajo.
    z["vm_n_superficie"] = informe.N_SUPERFICIE_MIN - 1
    sup = _busca(informe.comprobar(doc), "vm_p99_superficie")
    ok = (sup["estado"] == informe.RESERVAS
          and "vm_superficie_pocos" in sup["motivos"])
    registro.anotar(BLOQUE, "capa con pocos elementos: el p99 baja a reservas",
                    REF, float(informe.N_SUPERFICIE_MIN),
                    float(informe.N_SUPERFICIE_MIN - 1), 0.0,
                    "umbral declarado N_SUPERFICIE_MIN", ok,
                    nota=str(sup["motivos"]))
    assert ok, sup

    z["vm_n_superficie"] = informe.N_SUPERFICIE_MIN
    sup = _busca(informe.comprobar(doc), "vm_p99_superficie")
    assert sup["estado"] == informe.CITABLE, sup


# ---------------------------------------------------------------------------
# Lista de chequeo de reporte (Bouxsein 2010, Erdemir 2012)
# ---------------------------------------------------------------------------

def _fila_chequeo(filas, guia, item_es):
    for f in filas:
        if f["guia"] == guia and f["item"].startswith(item_es):
            return f
    raise AssertionError(f"sin fila {guia}/{item_es}")


def test_lista_chequeo_sigue_a_los_datos(registro):
    """Cada estado sale de la sesion: cambia cuando cambia el dato."""
    B, E = "bouxsein2010", "erdemir2012"
    filas = lista_chequeo.evaluar(_doc(), "es")
    base = {
        "segmentacion": _fila_chequeo(filas, B, "Segmentación")["estado"],
        "convergencia": _fila_chequeo(filas, E, "Convergencia")["estado"],
        "validacion": _fila_chequeo(filas, E, "Validación")["estado"],
        "vóxel": _fila_chequeo(filas, B, "Tamaño de vóxel")["estado"],
    }
    doc = _doc()
    doc["voi"]["recorte"] = {"umbral": 18983.0}
    seg = _fila_chequeo(lista_chequeo.evaluar(doc, "es"), B, "Segmentación")
    doc = _doc()
    doc["resultados"]["convergencia"] = {"puntos": [
        {"n": n, "ok": True, "E_app": 1e9, "residuo": 1e-9}
        for n in (24, 32, 40)]}
    conv = _fila_chequeo(lista_chequeo.evaluar(doc, "es"), E, "Convergencia")
    doc = _doc()
    for k in ("elastico", "resistencia"):
        doc["resultados"].pop(k)
    sin_fe = [f for f in lista_chequeo.evaluar(doc, "es") if f["guia"] == E]
    ok = (base == {"segmentacion": lista_chequeo.NO_CUMPLE,
                   "convergencia": lista_chequeo.NO_CUMPLE,
                   "validacion": lista_chequeo.NO_CUMPLE,
                   "vóxel": lista_chequeo.CUMPLE}
          and seg["estado"] == lista_chequeo.PARCIAL
          and "18983" in seg["donde"]
          and conv["estado"] == lista_chequeo.CUMPLE
          and len(sin_fe) == 1 and sin_fe[0]["estado"] ==
          lista_chequeo.NO_APLICA)
    registro.anotar(BLOQUE, "lista de chequeo: el estado sigue a los datos",
                    "Bouxsein 2010; Erdemir 2012", None, None, "exacto",
                    "umbral, convergencia y ausencia de FE cambian el estado",
                    ok, nota=f"{base} seg={seg['estado']} "
                             f"conv={conv['estado']} sin_fe={len(sin_fe)}")
    assert ok


def test_lista_chequeo_en_el_informe(registro):
    """ES y EN dan los mismos estados; la seccion y sus guias estan citadas."""
    doc = _doc()
    es = lista_chequeo.evaluar(doc, "es")
    en = lista_chequeo.evaluar(doc, "en")
    items = informe.comprobar(doc)
    paq = informe.paquete_reproduccion(doc, regenerar=False)
    md_es = informe.informe_markdown(doc, items, paq, "es")
    md_en = informe.informe_markdown(doc, items, paq, "en")
    c = lista_chequeo.resumen(es)
    ok = ([f["estado"] for f in es] == [f["estado"] for f in en]
          and all(f["estado"] in lista_chequeo.ORDEN for f in es)
          and all(f["accion"] for f in es
                  if f["estado"] in (lista_chequeo.PARCIAL,
                                     lista_chequeo.NO_CUMPLE))
          and "Lista de chequeo de reporte" in md_es
          and "Reporting checklist" in md_en
          and all(doi in md_en for doi in ("10.1002/jbmr.141",
                                           "10.1016/j.jbiomech.2011.11.038",
                                           "10.1016/j.medengphy.2021.03.011"))
          and md_en.rstrip().splitlines()[-1].startswith("| Erdemir 2012"))
    registro.anotar(BLOQUE, "lista de chequeo: ES = EN, citada, al final",
                    "Bouxsein 2010; Erdemir 2012", None, float(len(es)),
                    "exacto", "mismos estados en los dos idiomas", ok,
                    nota=json.dumps(c))
    assert ok


def test_insignia_readme_cuenta_los_bloques(registro):
    """La insignia y el texto del README dicen cuantos bloques hay."""
    n = len(list((RAIZ / "tests").glob("test_[0-9][0-9]_*.py")))
    es = (RAIZ / "README.md").read_text(encoding="utf-8")
    en = (RAIZ / "README.en.md").read_text(encoding="utf-8")
    ok = (f"-{n}%20bloques-" in es and f"\n{n} bloques con tolerancias" in es
          and f"-{n}%20blocks-" in en and f"\n{n} blocks with tolerances" in en)
    registro.anotar(BLOQUE, "README: numero de bloques de verificacion",
                    "README.md, README.en.md", float(n), None, "exacto",
                    "insignia y texto = archivos tests/test_NN_*.py", ok)
    assert ok, f"el README no dice {n} bloques"
