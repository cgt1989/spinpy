"""
Bloque 23 — Citabilidad de convergencia, dispersion y simulaciones, y las
figuras 9-11 del informe.

QUE SE VERIFICA
  convergencia  cada veredicto de `resistencia.estudio_convergencia` tiene su
                estado: serie plana o en meseta citable, sin meseta / sin
                extrapolar / deriva de densidad con reservas, no monotona o
                con menos de tres mallas no citable. Los umbrales saltan en
                el valor declarado (`CONV_ERROR_EXTRAPOLADO_MAX`,
                `CONV_DERIVA_RHO_MAX`) y el residuo que cuenta es el PEOR de
                la serie. Y sobre una estructura real: el veredicto que el
                dialogo escribe y el estado del informe coinciden.
  dispersion    K < 5 con reservas; la semilla de la busqueda como primera
                realizacion (maldicion del ganador) con reservas; replicas
                generadas siempre con reservas (no son especimenes).
  simulaciones  paso 0 sin resolver no citable; pasos fallidos con reservas,
                pero una estructura SIN CAMINO PORTANTE no es un fallo (su
                rigidez es cero); la pendiente log-log siempre con reservas;
                el fallo progresivo siempre con reservas, y los pasos tras el
                colapso no cuentan como fallidos.
  cuantiles     el cuantil 99 que se guarda para la figura 11 es EXACTAMENTE
                el `vm_p99_superficie` citable, sobre un ensayo real.
  figuras       9, 10 y 11 salen en PNG y PDF cuando hay datos y no salen
                cuando no; la separacion z es (media - VOI) / sd; el informe
                las ordena por numero (fig10 detras de fig9, no de fig1) y
                cada enlace existe.
"""
import copy
import json
from pathlib import Path

import numpy as np
import pytest

from spinpy import informe, resistencia
from spinpy import figuras as F

from test_19_informe import _doc

BLOQUE = "23 Estudios en el informe"
REF = ("informe.py (Estudio_Convergencia, Estudio_Discriminadores, "
       "Estudio_Simulaciones)")


def _malla_barras(n=24, paso=8, grosor=3):
    """Barras en x, y, z que atraviesan el cubo entero: tocan base y techo,
    asi que el ensayo de compresion tiene camino portante."""
    BW = np.zeros((n, n, n), dtype=bool)
    for a in range(paso // 2, n, paso):
        for b in range(paso // 2, n, paso):
            s = slice(a - grosor // 2, a + grosor // 2 + 1)
            t = slice(b - grosor // 2, b + grosor // 2 + 1)
            BW[:, s, t] = BW[s, :, t] = BW[s, t, :] = True
    return BW


def _busca(items, magnitud, estructura):
    for it in items:
        if it["magnitud"] == magnitud and it["estructura"] == estructura:
            return it
    raise AssertionError(f"sin item {magnitud}/{estructura}")


# ---------------------------------------------------------------------------
# Convergencia
# ---------------------------------------------------------------------------

def _conv(E=(800e6, 900e6, 950e6, 970e6), rho=(0.30,) * 4, **extra):
    ns = (22, 28, 34, 40)[:len(E)]
    puntos = [{"n": n, "h": 5.0 / n, "rho": r, "ok": True, "E_app": e,
               "frac_portante": 0.999, "residuo": 1e-8, "solver": "AMG"}
              for n, e, r in zip(ns, E, rho)]
    rec = {"puntos": puntos, "eje": "z", "estructura_codigo": "voi"}
    rec.update(extra)
    return rec


def _items_conv(rec):
    d = _doc()
    d["resultados"]["convergencia"] = rec
    return informe.comprobar(d)


CASOS_CONV = [
    ("meseta: error 4.9 %", dict(monotona=True, deriva_rho_rel=0.0,
                                 dispersion_rel=0.18,
                                 error_estimado_rel=0.049,
                                 E_extrapolado=1.02e9),
     None, informe.CITABLE),
    ("sin meseta: error 5 % justo", dict(monotona=True, deriva_rho_rel=0.0,
                                         dispersion_rel=0.18,
                                         error_estimado_rel=0.05,
                                         E_extrapolado=1.02e9),
     "conv_sin_meseta", informe.RESERVAS),
    ("deriva de densidad 2 % justo", dict(monotona=True, deriva_rho_rel=0.02,
                                          dispersion_rel=0.18),
     "conv_deriva_rho", informe.RESERVAS),
    ("deriva 1.9 % y sin orden", dict(monotona=True, deriva_rho_rel=0.019,
                                      dispersion_rel=0.18,
                                      salto_final_rel=0.02),
     "conv_sin_extrapolar", informe.RESERVAS),
    ("no monotona", dict(monotona=False, deriva_rho_rel=0.0,
                         dispersion_rel=0.4),
     "conv_no_monotona", informe.NO_CITABLE),
    ("oscila 3 % justo", dict(monotona=False, deriva_rho_rel=0.0,
                              dispersion_rel=0.03),
     "conv_oscila", informe.RESERVAS),
    ("oscila 3.01 %", dict(monotona=False, deriva_rho_rel=0.0,
                           dispersion_rel=0.0301),
     "conv_no_monotona", informe.NO_CITABLE),
    ("serie plana", dict(monotona=True, deriva_rho_rel=0.0,
                         dispersion_rel=1e-9),
     None, informe.CITABLE),
]


@pytest.mark.parametrize("nombre,extra,motivo,estado", CASOS_CONV,
                         ids=[c[0] for c in CASOS_CONV])
def test_convergencia_umbrales(registro, nombre, extra, motivo, estado):
    it = _busca(_items_conv(_conv(**extra)), "convergencia", "voi")
    ok = (it["estado"] == estado and
          (motivo in it["motivos"] if motivo else not it["motivos"]))
    registro.anotar(BLOQUE, f"convergencia: {nombre}", REF, None, None,
                    "exacto", f"{motivo or 'sin motivo'} -> {estado}", ok,
                    nota=f"{it['estado']} {it['motivos']}")
    assert ok


def test_convergencia_pocos_y_residuo(registro):
    pocos = _busca(_items_conv(_conv(E=(8e8, 9e8))), "convergencia", "voi")
    rec = _conv(monotona=True, deriva_rho_rel=0.0, dispersion_rel=0.18,
                error_estimado_rel=0.01, E_extrapolado=9.8e8)
    # El punto MAS GRUESO es el que falla el residuo: igualmente no citable.
    rec["puntos"][0]["residuo"] = 3e-6
    peor = _busca(_items_conv(rec), "convergencia", "voi")
    # n = 40 en la malla mas fina: no hay aviso de resolucion.
    ok = (pocos["motivos"][:1] == ["conv_pocos_puntos"]
          and pocos["estado"] == informe.NO_CITABLE
          and "residuo" in peor["motivos"]
          and peor["estado"] == informe.NO_CITABLE
          and "resolucion_mecanica" not in peor["motivos"])
    registro.anotar(BLOQUE, "convergencia: 2 mallas y peor residuo", REF,
                    None, None, "exacto", "no citable en ambos", ok,
                    nota=f"{pocos['motivos']} / {peor['motivos']}")
    assert ok


def test_convergencia_estructura_de(registro):
    d = _doc()
    casos = {"VOI": "voi", "Dual-lattice": "dual-lattice",
             "Spinodoid": "spinodoide", None: "voi",
             "dual-lattice": "dual-lattice"}
    sin_voi = dict(d, voi=None)
    ok = (all(informe.estructura_de(k, d) == v for k, v in casos.items())
          and informe.estructura_de(None, sin_voi) == "spinodoide")
    registro.anotar(BLOQUE, "estructura de un registro viejo por su etiqueta",
                    REF, None, None, "exacto", "codigo correcto", ok)
    assert ok


def test_convergencia_real_coincide_con_veredicto(registro):
    """Sobre una estructura de verdad: lo que el dialogo dice y lo que el
    informe grada no pueden discrepar."""
    BW = _malla_barras(30, paso=10, grosor=3)
    rec = resistencia.estudio_convergencia(BW, 0.05, resoluciones=(16, 20, 24),
                                           eje=2)
    rec["estructura_codigo"] = "voi"
    it = _busca(_items_conv(rec), "convergencia", "voi")
    v = rec.get("veredicto", "")
    if "NO es monotona" in v:
        esperado = "conv_no_monotona"
    elif "oscila sin tendencia" in v:
        esperado = "conv_oscila"
    elif "no varia" in v or "citable declarando" in v:
        esperado = None
    elif "Hacen falta" in v:
        esperado = "conv_pocos_puntos"
    else:
        esperado = {"conv_sin_meseta", "conv_sin_extrapolar",
                    "conv_deriva_rho"}
    motivos = [m for m in it["motivos"] if m.startswith("conv_")]
    if isinstance(esperado, set):
        ok = len(motivos) == 1 and motivos[0] in esperado
    else:
        ok = motivos == ([esperado] if esperado else [])
    registro.anotar(BLOQUE, "convergencia real: veredicto == estado", REF,
                    None, None, "exacto", str(esperado), ok,
                    nota=f"{v[:60]} | {it['motivos']}")
    assert ok


# ---------------------------------------------------------------------------
# Dispersion y replicas
# ---------------------------------------------------------------------------

def test_dispersion_y_replicas(registro):
    d = _doc()
    d["resultados"]["dispersion"] = {"resumen": {"BVTV": {"media": 0.3,
                                                          "sd": 0.01}},
                                     "n_semillas": 8, "semilla_base": 12}
    limpio = _busca(informe.comprobar(d), "dispersion", "spinodoide")
    d["resultados"]["dispersion"]["n_semillas"] = 4
    d["resultados"]["dispersion"]["semilla_base"] = 11   # la de la busqueda
    malo = _busca(informe.comprobar(d), "dispersion", "spinodoide")
    d["resultados"]["replicas"] = {"n": 20, "logrado": {}}
    rep = _busca(informe.comprobar(d), "replicas", "spinodoide")
    ok = (limpio["estado"] == informe.CITABLE and not limpio["motivos"]
          and set(malo["motivos"]) == {"sin_replicas", "semilla_busqueda"}
          and rep["motivos"] == ["replicas_no_especimenes"]
          and "20" in informe.texto_motivos(rep, "es"))
    registro.anotar(BLOQUE, "dispersion: K y semilla de busqueda; replicas",
                    REF, None, None, "exacto",
                    "limpio citable; K=4+semilla 11 reservas", ok,
                    nota=f"{malo['motivos']} / {rep['motivos']}")
    assert ok


# ---------------------------------------------------------------------------
# Simulaciones
# ---------------------------------------------------------------------------

def _paso(i, E, ok=True, residuo=1e-8, **kw):
    f = {"paso": i, "fase": "inicial" if i == 0 else "perdida", "ok": ok,
         "E_app": E, "residuo": residuo, "frac_portante": 0.999,
         "n_mec": 48, "BVTV": 0.3 - 0.02 * i}
    f.update(kw)
    return f


def _sim(filas, tbth_h=2.5, pendiente=None):
    r = {"E_rel_final": 0.5, "BVTV_rel_final": 0.8, "TbTh_h_mec": tbth_h}
    if pendiente is not None:
        r["pendiente_loglog"] = pendiente
    return {"protocolo": "adelgazamiento", "pasos": filas, "resumen": r,
            "parametros": {"n_mec": 48, "eje": 2}, "estructura": "VOI",
            "estructura_codigo": "voi"}


def test_simulacion_perdida(registro):
    d = _doc()
    base = [_paso(0, 1e9), _paso(1, 8e8), _paso(2, 6e8)]
    d["resultados"]["simulacion_perdida"] = _sim(base)
    limpio = _busca(informe.comprobar(d), "sim_rigidez", "voi")

    # Un paso sin camino portante es un resultado (E = 0), no un fallo.
    sin_camino = base + [_paso(3, float("nan"), ok=False, sin_camino=True)]
    d["resultados"]["simulacion_perdida"] = _sim(sin_camino)
    sc = _busca(informe.comprobar(d), "sim_rigidez", "voi")

    fallido = base + [_paso(3, 5e8, residuo=2e-6)]
    d["resultados"]["simulacion_perdida"] = _sim(fallido, tbth_h=1.6,
                                                 pendiente=3.1)
    items = informe.comprobar(d)
    fa = _busca(items, "sim_rigidez", "voi")
    pe = _busca(items, "sim_pendiente", "voi")

    d["resultados"]["simulacion_perdida"] = _sim(
        [_paso(0, float("nan"), ok=False)] + base[1:])
    p0 = _busca(informe.comprobar(d), "sim_rigidez", "voi")
    ok = (limpio["estado"] == informe.CITABLE and not limpio["motivos"]
          and not sc["motivos"]
          and set(fa["motivos"]) == {"sim_pasos_fallidos", "sim_tbth_h"}
          and fa["estado"] == informe.RESERVAS
          and "sim_pendiente" in pe["motivos"]
          and p0["estado"] == informe.NO_CITABLE
          and "sim_paso0" in p0["motivos"])
    registro.anotar(BLOQUE, "simulacion de perdida: paso 0, fallidos, "
                    "sin camino, Tb.Th/h 1.6 < 1.7, pendiente", REF, None,
                    None, "exacto", "segun GRAVEDAD", ok,
                    nota=f"{fa['motivos']} / {p0['motivos']}")
    assert ok


def test_fallo_progresivo(registro):
    d = _doc()
    filas = [dict(_paso(0, 1e9), tras_colapso=False),
             dict(_paso(1, 7e8), tras_colapso=False),
             dict(_paso(2, 3e8), tras_colapso=True),
             dict(_paso(3, float("nan"), ok=False), tras_colapso=True)]
    d["resultados"]["fallo_progresivo"] = {
        "tipo": "fallo_progresivo", "pasos": filas,
        "resumen": {"F_max": 120.0, "E_rel_final": 0.3, "TbTh_h_mec": 3.0},
        "parametros": {"n_mec": 48, "eje": 2, "rigidez_danada": 0.05},
        "estructura": "Spinodoide", "estructura_codigo": "spinodoide"}
    it = _busca(informe.comprobar(d), "fallo_progresivo", "spinodoide")
    ok = (set(it["motivos"]) == {"pistoia_calibracion", "fallo_modelo_dano"}
          and it["estado"] == informe.RESERVAS and it["eje"] == "z"
          and "0,05" in informe.texto_motivos(it, "es"))
    registro.anotar(BLOQUE, "fallo progresivo: siempre con reservas; tras el "
                    "colapso no cuenta", REF, None, None, "exacto",
                    "pistoia + modelo de dano", ok, nota=str(it["motivos"]))
    assert ok


def test_textos_nuevos(registro):
    nuevos = [c for c in informe.GRAVEDAD if c.startswith(
        ("conv_", "sim_", "semilla_", "replicas_", "fallo_modelo"))]
    malos = []
    for c in nuevos:
        for idioma in ("es", "en"):
            t = informe.texto_motivo(c, {}, idioma)
            if not t or t != t.strip():
                malos.append((c, idioma))
    ok = len(nuevos) == 13 and not malos
    registro.anotar(BLOQUE, "textos ES/EN de los 13 motivos nuevos", REF, 13,
                    len(nuevos), "exacto", "todos con texto", ok,
                    nota=str(malos))
    assert ok


# ---------------------------------------------------------------------------
# Cuantiles y figuras
# ---------------------------------------------------------------------------

def test_cuantil_99_es_el_citable(registro):
    BW = _malla_barras()
    r = resistencia.ensayo_compresion(BW, [0.05] * 3, unidad="mm",
                                      carga_N=100.0)
    sup = resistencia.estadisticos_vm(r)
    c = resistencia.cuantiles_vm_superficie(r)
    k = c["p"].index(99.0)
    ok = (c["valor"][k] == sup["vm_p99_superficie"]
          and c["n"] == sup["vm_n_superficie"]
          and np.all(np.diff(c["valor"]) >= 0)
          and len(c["p"]) == len(resistencia.PROB_CUANTILES_VM))
    registro.anotar(BLOQUE, "cuantil 99 guardado == vm_p99_superficie", REF,
                    sup["vm_p99_superficie"], c["valor"][k], "exacto",
                    "igualdad bit a bit", ok)
    assert ok


def _doc_completo():
    d = _doc()
    res = d["resultados"]
    res["convergencia"] = _conv(monotona=True, deriva_rho_rel=0.0,
                                dispersion_rel=0.18,
                                error_estimado_rel=0.03,
                                E_extrapolado=9.9e8, orden=1.9)
    res["ajuste"]["incertidumbre"].update(
        errores=[0.05, 0.07, 0.06, 0.055, 0.065],
        metricas={"BVTV": {"media": 0.31, "sd": 0.005},
                  "TbTh": {"media": 0.15, "sd": 0.002}})
    res["dispersion"] = {"resumen": {"SMI": {"media": 1.8, "sd": 0.07,
                                             "z": 14.0}},
                         "n_semillas": 8, "semilla_base": 40}
    p = list(resistencia.PROB_CUANTILES_VM)
    res["analisis_comparado"] = {"carga_N": 100.0, "resolucion": 48,
                                 "por_estructura": {
        k: {"E_app": 1e3, "vm_p99_superficie": 2.0 * s,
            "vm_n_superficie": 5000, "residuo": 1e-8, "frac_portante": 1.0,
            "cuantiles_vm_superficie": {
                "p": p, "valor": list(s * np.linspace(0.1, 2.5, len(p))),
                "n": 5000}}
        for k, s in (("voi", 1.0), ("spin", 1.6))}}
    return d


def test_separacion_z(registro):
    d = _doc_completo()
    z = F.separaciones(d)["spinodoide"]
    # VOI: BVTV 0.30, TbTh 0.195 (de _doc)
    ok = (abs(z["BVTV"][0] - (0.31 - 0.30) / 0.005) < 1e-12
          and abs(z["TbTh"][0] - (0.15 - 0.195) / 0.002) < 1e-12
          and z["SMI"] == (14.0, "dispersion")
          and z["BVTV"][1] == "ajuste")
    registro.anotar(BLOQUE, "separacion z = (media - VOI) / sd", REF, 2.0,
                    z["BVTV"][0], "1e-12", "a mano", ok)
    assert ok


def test_figuras_y_orden(registro, tmp_path):
    d = _doc_completo()
    items = informe.comprobar(d)
    r9 = F.fig_convergencia(d, items, tmp_path / "f9")
    r10 = F.fig_incertidumbre(d, items, tmp_path / "f10")
    r11 = F.fig_von_mises_superficie(d, items, tmp_path / "f11", "en")
    vacio = _doc()
    vacias = (F.fig_convergencia(vacio, [], tmp_path / "v9")
              + F.fig_von_mises_superficie(vacio, [], tmp_path / "v11"))
    sin_inc = copy.deepcopy(vacio)
    sin_inc["resultados"]["ajuste"]["incertidumbre"] = {}
    vacias += F.fig_incertidumbre(sin_inc, [], tmp_path / "v10")

    orden = sorted(["figuras/fig10_x.png", "figuras/fig1_x.png",
                    "figuras/fig9_x.png", "figuras/fig11_x.png",
                    "figuras/fig2_x.png"], key=informe._numero_figura)

    out = informe.escribir_informe(d, tmp_path / "inf", regenerar=False,
                                   render_3d=False, pdf=False)
    md = (tmp_path / "inf" / informe.ARCHIVOS["md_es"]).read_text("utf-8")
    enlaces = [l.split("](")[1].rstrip(")") for l in md.splitlines()
               if l.startswith("![")]
    nums = [informe._numero_figura(e)[0] for e in enlaces]
    ok = (all(Path(x).exists() for x in r9 + r10 + r11)
          and [Path(x).suffix for x in r9] == [".png", ".pdf"]
          and len(r10) == 2 and len(r11) == 2 and not vacias
          and orden == ["figuras/fig1_x.png", "figuras/fig2_x.png",
                        "figuras/fig9_x.png", "figuras/fig10_x.png",
                        "figuras/fig11_x.png"]
          and {9, 10, 11} <= set(nums) and nums == sorted(nums)
          and all((tmp_path / "inf" / e).exists() for e in enlaces))
    registro.anotar(BLOQUE, "figuras 9-11: con datos salen, sin datos no; "
                    "orden numerico en el informe", REF, None, None,
                    "exacto", "PNG+PDF, fig9 < fig10 < fig11", ok,
                    nota=json.dumps({"enlaces": nums,
                                     "resumen": out["resumen"]}))
    assert ok
