"""
Bloque 26 — El apoyo del ensayo de compresion viaja como CLAVE, no como texto.

EL FALLO QUE CIERRA
  El visor pasaba `cmb_apoyo.currentText()` a `ensayo_compresion`. Con la
  interfaz en ingles, `idioma_textos` traduce "empotrado" -> "fixed", y el
  ensayo decidia con `startswith("empotr")`: elegir "fixed" ejecutaba EN
  SILENCIO el apoyo DESLIZANTE, y el registro exportado decia "fixed". El
  apoyo cambia E_app y la carga de fallo sin dar ningun sintoma.

LO QUE SE VERIFICA
  (1) `normalizar_apoyo` acepta solo las claves canonicas (con mayusculas o
      espacios alrededor) y LANZA ValueError con cualquier otra cosa, incluidas
      las traducciones del propio diccionario de idioma.
  (2) `ensayo_compresion`, `ensayo_compresion_eje`, `escribir_febio` y las dos
      simulaciones rechazan "fixed" antes de calcular nada.
  (3) Las dos claves producen ensayos DISTINTOS sobre un cubo macizo (si el
      apoyo no llegara al ensayo, darian lo mismo) y el registro guarda la
      clave canonica.
  (4) El visor ya no lee el texto del combo, y el orden de sus elementos es
      el de `resistencia.APOYOS` (el helper lee el INDICE).
  (5) El informe nombra el apoyo EJECUTADO: un registro antiguo que dice
      "fixed" se describe como deslizante, que es lo que se calculo.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  Todo exacto: son cadenas, excepciones e igualdades de flotantes producidos
  por el mismo codigo. (3) solo exige E_app(empotrado) > E_app(deslizante):
  la coaccion lateral solo puede rigidizar (bloque 05 da las cotas).
"""
import ast
from pathlib import Path

import numpy as np
import pytest

from spinpy import informe
from spinpy.escribe import escribir_febio
from spinpy.idioma_textos import EN
from spinpy.resistencia import (APOYOS, ensayo_compresion,
                                ensayo_compresion_eje, normalizar_apoyo)
from spinpy.simulacion import fallo_progresivo, simular_perdida
from spinpy.solido import malla_hex

BLOQUE = "26 Apoyo canonico"
REF = "correccion del apoyo traducido (visor en ingles)"
VISOR = Path(__file__).resolve().parent.parent / "visor.py"


def _anotar_bool(registro, prueba, ok, criterio="exacto", nota=""):
    registro.anotar(BLOQUE, prueba, REF, 1, int(bool(ok)), 0.0, criterio,
                    bool(ok), nota=nota)


def test_normalizar_acepta_solo_claves_canonicas(registro):
    buenos = {"deslizante": "deslizante", "empotrado": "empotrado",
              "Empotrado": "empotrado", " DESLIZANTE ": "deslizante"}
    for entrada, esperado in buenos.items():
        ok = normalizar_apoyo(entrada) == esperado
        _anotar_bool(registro, f"acepta {entrada!r}", ok)
        assert ok, entrada
    # Las traducciones del diccionario van PRIMERO: son el caso real.
    malos = [EN["empotrado"], EN["deslizante"], "empotr", "libre", "", None]
    for entrada in malos:
        with pytest.raises(ValueError):
            normalizar_apoyo(entrada)
        _anotar_bool(registro, f"rechaza {entrada!r}", True)


def test_ensayos_y_exportacion_rechazan_fixed(registro, tmp_path):
    BW = np.ones((4, 4, 4), bool)
    sp = [0.1] * 3
    casos = {
        "ensayo_compresion": lambda: ensayo_compresion(BW, sp, apoyo="fixed"),
        "ensayo_compresion_eje": lambda: ensayo_compresion_eje(
            BW, sp, eje=0, apoyo="fixed"),
        "simular_perdida": lambda: simular_perdida(
            BW, sp, pasos=1, n_mec=4, apoyo="fixed"),
        "fallo_progresivo": lambda: fallo_progresivo(
            BW, sp, pasos=1, n_mec=4, apoyo="fixed"),
    }
    nodos, elems, _ = malla_hex(BW, np.asarray(sp))
    casos["escribir_febio"] = lambda: escribir_febio(
        nodos, elems, tmp_path / "c.feb", apoyo="fixed")
    for nombre, f in casos.items():
        with pytest.raises(ValueError):
            f()
        _anotar_bool(registro, f"{nombre}(apoyo='fixed') lanza ValueError",
                     True)
    assert not (tmp_path / "c.feb").exists()


def test_las_dos_claves_dan_ensayos_distintos(registro):
    BW = np.ones((8, 8, 8), bool)
    kw = dict(E_s=10.0, nu_s=0.3, sigma0=1e-3)
    r_d = ensayo_compresion(BW, [1.0] * 3, apoyo="deslizante", **kw)
    r_e = ensayo_compresion(BW, [1.0] * 3, apoyo="empotrado", **kw)
    r_E = ensayo_compresion(BW, [1.0] * 3, apoyo=" Empotrado", **kw)
    assert r_d["ok"] and r_e["ok"] and r_E["ok"]
    ok = r_e["E_app"] > r_d["E_app"]
    registro.anotar(BLOQUE, "E_app empotrado > deslizante (cubo macizo)",
                    "bloque 05: la coaccion lateral rigidiza",
                    r_d["E_app"], r_e["E_app"], "desigualdad",
                    "E_emp > E_desl", ok)
    assert ok
    ok = (r_E["E_app"] == r_e["E_app"] and r_E["apoyo"] == "empotrado"
          and r_d["apoyo"] == "deslizante")
    _anotar_bool(registro, "el registro guarda la clave canonica", ok)
    assert ok


def test_visor_lee_el_indice_del_combo(registro):
    fuente = VISOR.read_text(encoding="utf-8")
    ok = "cmb_apoyo.currentText()" not in fuente
    _anotar_bool(registro, "visor.py no lee el texto de cmb_apoyo", ok)
    assert ok

    # El combo principal se llena con una lista literal: su orden es el que
    # `apoyo_de_combo` traduce por indice con APOYOS.
    listas = []
    for nodo in ast.walk(ast.parse(fuente)):
        if (isinstance(nodo, ast.Call)
                and isinstance(nodo.func, ast.Attribute)
                and nodo.func.attr == "addItems"
                and isinstance(nodo.func.value, ast.Attribute)
                and nodo.func.value.attr == "cmb_apoyo"
                and nodo.args and isinstance(nodo.args[0], ast.List)):
            listas.append(tuple(ast.literal_eval(nodo.args[0])))
    ok = listas == [APOYOS]
    _anotar_bool(registro, "orden del combo == resistencia.APOYOS", ok,
                 nota=str(listas))
    assert ok, listas


def test_informe_nombra_el_apoyo_ejecutado(registro):
    casos = [("empotrado", "es", "empotrado"), ("empotrado", "en", "fixed"),
             ("deslizante", "en", "sliding"),
             # Registro anterior a la correccion: "fixed" ejecutaba deslizante.
             ("fixed", "es", "deslizante"), ("fixed", "en", "sliding"),
             (None, "en", "—")]
    for apoyo, idioma, esperado in casos:
        obt = informe.nombre_apoyo(apoyo, idioma)
        ok = obt == esperado
        _anotar_bool(registro, f"nombre_apoyo({apoyo!r}, {idioma!r})", ok,
                     nota=f"{obt!r} (esperado {esperado!r})")
        assert ok, (apoyo, idioma, obt)
