"""
Bloque 16 — Procedencia de los resultados (`spinpy.procedencia`).

QUE SE VERIFICA
  Que el error que invalido el problema abierto de Estudio_Convergencia —el
  valor del deslizador (15) transcrito como radianes, sin pi— ya no puede
  pasar en silencio por el paquete, y que los archivos de antes se siguen
  leyendo pero con aviso.

  version      spinpy.__version__ coincide con pyproject.toml.
  onda         una lectura da la otra; las dos se contrastan; el caso
               "rad = valor del deslizador" se rechaza y el mensaje lo nombra.
  traduccion   `kwargs_generador` es la unica: rechaza la clave suelta `wave`
               y regenera bit a bit lo que genera `generar_mascara` con 15 pi.
  migracion    un JSON y un CSV de formato 1 se leen con AvisoFormatoAntiguo,
               con las dos lecturas anadidas y el supuesto declarado.
  ajuste       el resultado de un ajuste trae su bloque, sus dos lecturas y
               una traza sin la clave `wave`.
"""
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pytest

import spinpy
from spinpy import procedencia
from spinpy.grf import generar_mascara
from spinpy.lote import generar_desde_parametros

BLOQUE = "16 Procedencia"
REF = "procedencia.py (el pi de Estudio_Convergencia y Estudio_Percolacion)"


def test_version_coincide_con_pyproject(registro):
    txt = (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text(
        encoding="utf-8")
    v = re.search(r'^version\s*=\s*"([^"]+)"', txt, re.M).group(1)
    ok = v == spinpy.__version__
    registro.anotar(BLOQUE, "version del paquete", "pyproject.toml", None,
                    None, "igual", "una sola fuente", ok,
                    nota=f"{spinpy.__version__} / {v}")
    assert ok


def test_numero_onda_dos_lecturas(registro):
    a = procedencia.numero_onda(pi=15)
    b = procedencia.numero_onda(rad=15 * np.pi)
    c = procedencia.numero_onda(pi=15, rad=15 * np.pi)
    # 15 pi / pi no es 15 bit a bit: se compara con la tolerancia del modulo.
    ok = (a == c and a["wave_number_rad"] == 15 * np.pi
          and np.isclose(b["wave_number_pi"], 15.0, rtol=1e-12, atol=0)
          and b["wave_number_rad"] == 15 * np.pi)
    with pytest.raises(procedencia.ErrorNumeroOnda, match="falta multiplicar"):
        procedencia.numero_onda(pi=15, rad=15.0)
    with pytest.raises(procedencia.ErrorNumeroOnda):
        procedencia.numero_onda()
    registro.anotar(BLOQUE, "lecturas del numero de onda", REF,
                    15 * np.pi, a["wave_number_rad"], "exacto",
                    "rad = pi x deslizador; 15 como rad se rechaza", ok)
    assert ok


def _par(**extra):
    p = {"familia": "spinodoide", "resolucion": 24, "densidad": 0.3,
         "wave_number_pi": 15.0, "num_waves": 200, "thetas": [15, 15, 60],
         "R": np.eye(3).tolist(), "esquema": "rechazo", "semilla": 11}
    p.update(extra)
    return p


def test_una_sola_traduccion(registro):
    kw = procedencia.kwargs_generador(_par())
    ref, _, _ = generar_mascara(resolution=24, wave_number=15 * np.pi,
                                num_waves=200, thetas=[15, 15, 60], rho=0.3,
                                R=np.eye(3), esquema="rechazo", seed=11)
    BW = generar_desde_parametros(_par(), 11)
    ok = (kw["wave_number"] == 15 * np.pi and np.array_equal(BW, ref))
    solo_wave = _par()
    del solo_wave["wave_number_pi"]
    solo_wave["wave"] = 15
    with pytest.raises(procedencia.ErrorNumeroOnda, match="ambigua"):
        procedencia.kwargs_generador(solo_wave)
    with pytest.raises(procedencia.ErrorNumeroOnda):
        procedencia.kwargs_generador(_par(wave_number_rad=15.0))
    registro.anotar(BLOQUE, "kwargs_generador regenera 15 pi", REF, None,
                    None, "mascara identica", "unica traduccion", ok)
    assert ok


def test_migra_json_antiguo_con_aviso(tmp_path, registro):
    viejo = {"params": {"wave": 15, "dens": 0.3},
             "traza": [{"wave_number_pi": 10.0}, {"wave_number_pi": 12.0},
                       {"wave_number_pi": 14.0}]}
    ruta = tmp_path / "resultados.json"
    ruta.write_text(json.dumps(viejo), encoding="utf-8")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        doc = procedencia.leer_json(ruta)
    textos = [str(x.message) for x in w
              if issubclass(x.category, procedencia.AvisoFormatoAntiguo)]
    ok = (doc["params"]["wave_number_rad"] == 15 * np.pi
          and doc["params"]["wave_number_pi"] == 15.0
          and doc["traza"][2]["wave_number_rad"] == 14.0 * np.pi
          and any("formato 1" in t for t in textos)
          and any("ambigua 'wave'" in t for t in textos)
          and any("3 sitio(s)" in t for t in textos))
    # Un documento que ya es de formato 2 no avisa.
    nuevo = {"procedencia": procedencia.bloque("spinodoide", semilla=1,
                                                wave_number_pi=15)}
    ruta2 = tmp_path / "nuevo.json"
    ruta2.write_text(json.dumps(nuevo), encoding="utf-8")
    with warnings.catch_warnings(record=True) as w2:
        warnings.simplefilter("always")
        procedencia.leer_json(ruta2)
    ok &= not [x for x in w2
               if issubclass(x.category, procedencia.AvisoFormatoAntiguo)]
    registro.anotar(BLOQUE, "JSON de formato 1", REF, None, None,
                    "migra y avisa; el de formato 2 no avisa",
                    "AvisoFormatoAntiguo", ok, nota=" | ".join(textos)[:300])
    assert ok, textos


def test_migra_csv_antiguo_con_aviso(tmp_path, registro):
    ruta = tmp_path / "lote_ajustes.csv"
    ruta.write_text("archivo;densidad;wave\nVOI_a.vtk;0,3;15\n",
                    encoding="utf-8-sig")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        df = procedencia.leer_csv(ruta)
    ok = (abs(float(df["wave_number_rad"][0]) - 15 * np.pi) < 1e-12
          and len([x for x in w if issubclass(
              x.category, procedencia.AvisoFormatoAntiguo)]) == 2)
    registro.anotar(BLOQUE, "CSV de formato 1", REF, 15 * np.pi,
                    float(df["wave_number_rad"][0]), "1e-12",
                    "migra y avisa", ok)
    assert ok


@pytest.mark.lento
def test_ajuste_lleva_procedencia(registro):
    VOI, _, _ = generar_mascara(resolution=32, wave_number=12 * np.pi,
                                num_waves=400, thetas=[15, 15, 45], rho=0.32,
                                seed=5)
    r = spinpy.ajustar_spinodoide(VOI, np.full(3, 0.1), modo="rapido",
                                  num_waves=400, replicas=0)
    p, pr = r["parametros"], r["procedencia"]
    ok = (pr["version_formato"] == procedencia.VERSION_FORMATO
          and pr["spinpy"] == spinpy.__version__
          and pr["esquema"] == "rechazo" and pr["semilla"] == p["semilla"]
          and p["wave_number_rad"] == p["wave_number_pi"] * np.pi
          and pr["wave_number_rad"] == p["wave_number_rad"]
          and all("wave" not in t for t in r["traza"])
          and all(t["wave_number_rad"] == t["wave_number_pi"] * np.pi
                  for t in r["traza"]))
    BW = generar_desde_parametros(p, p["semilla"])
    registro.anotar(BLOQUE, "resultado de un ajuste", REF, None, None,
                    "bloque, dos lecturas, traza sin 'wave'", "ajuste rapido",
                    ok, nota=f"regenera {BW.shape}")
    assert ok
