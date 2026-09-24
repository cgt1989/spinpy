"""
Bloque 25 — Exportacion a FEBio: el archivo plantea el MISMO ensayo que spinpy.

REFERENCIA
  `comparativa_febio/PREDICCIONES.md` y `RESULTADOS.md`: el ensayo de
  compresion de spinpy frente a FEBio 4.5 sobre los VOIs de H4 y el
  espinodoide ajustado. Esta prueba es la version pequena y repetible.

LO QUE SE VERIFICA
  (1) Sin FEBio: el .feb que escribe `escribe.escribir_febio` lleva la malla
      entera, la presion cuya integral sobre el area osea del techo es
      sigma_app * A_bruta, y ancla las MISMAS esquinas que el apoyo deslizante
      de `resistencia.ensayo_compresion`.
  (2) Con FEBio instalado (si no, se salta): sobre una giroide pequena, la
      solucion de FEBio extrapolada a carga nula coincide con la de spinpy.
      Se extrapola porque FEBio es no lineal geometricamente: con una sola
      carga pequena, los voladizos del techo dejan un desvio proporcional a
      la carga (3e-4 a 1 kPa en el VOI proximal de H4).

TOLERANCIAS DECLARADAS ANTES DE MEDIR (las de PREDICCIONES.md)
  (1) conteos exactos; fuerza 1e-12 relativa (suma en coma flotante).
  (2) E_app 1e-6, desplazamientos 1e-5 relativo al maximo.
"""
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

from spinpy.escribe import escribir_febio
from spinpy.resistencia import _solo_portante, ensayo_compresion
from spinpy.solido import malla_hex

BLOQUE = "25 Exportacion FEBio"
REF = "comparativa_febio/PREDICCIONES.md"
FEBIO = Path(r"C:\Program Files\FEBioStudio2\bin\febio4.exe")


def _giroide(n, nivel=-0.6):
    x = np.linspace(0.0, 1.0, n)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    f = (np.sin(4 * np.pi * X) * np.cos(4 * np.pi * Y)
         + np.sin(4 * np.pi * Y) * np.cos(4 * np.pi * Z)
         + np.sin(4 * np.pi * Z) * np.cos(4 * np.pi * X))
    return f > nivel


def _caso(n=12, h=0.1):
    BW = _solo_portante(_giroide(n))
    spc = np.full(3, h)
    nodos, elems, _ = malla_hex(BW, spc)
    return BW, spc, nodos, elems, float(n * h * n * h)


def test_archivo_plantea_el_mismo_ensayo(registro, tmp_path):
    BW, spc, nodos, elems, A = _caso()
    sigma = 1e6
    inf = escribir_febio(nodos, elems, tmp_path / "g.feb", sigma_app=sigma,
                         A_bruta=A)
    raiz = ET.parse(tmp_path / "g.feb").getroot()
    malla = raiz.find("Mesh")
    n_nod = len(malla.find("Nodes"))
    n_el = len(malla.find("Elements"))
    caras = malla.find("Surface")
    p = float(raiz.find("Loads/surface_load/pressure").text)

    # Area osea del techo, contada a mano: voxeles en la ultima capa.
    A_osea = float(BW[:, :, -1].sum() * spc[0] * spc[1])
    F = p * A_osea                                   # N
    F_esp = sigma / 1e6 * A

    # Esquinas del apoyo deslizante, con la regla de ensayo_compresion.
    z = nodos[:, 2]
    base = np.nonzero(z <= z.min() + 1e-9)[0]
    cb = nodos[base]
    a = base[int(np.argmin(cb[:, 0] + cb[:, 1]))] + 1
    b = base[int(np.argmax(cb[:, 0] - cb[:, 1]))] + 1
    ns = {e.get("name"): e.text.replace(",", " ").split()
          for e in malla.findall("NodeSet")}

    casos = [("nodos", nodos.shape[0], n_nod),
             ("elementos", elems.shape[0], n_el),
             ("caras cargadas = voxeles del techo", int(BW[:, :, -1].sum()),
              len(caras)),
             ("ancla (x,y) = esquina de spinpy", a, int(ns["ancla_xy"][0])),
             ("ancla y = esquina de spinpy", b, int(ns["ancla_y"][0]))]
    for nombre, esp, obt in casos:
        registro.anotar(BLOQUE, nombre, REF, esp, obt, 0.0, "exacto",
                        esp == obt)
        assert esp == obt, (nombre, esp, obt)

    err = abs(F - F_esp) / F_esp
    registro.anotar(BLOQUE, "presion x area osea = sigma_app x A_bruta", REF,
                    F_esp, F, err, "1e-12 relativo", err < 1e-12)
    assert err < 1e-12
    assert inf["F_N"] == pytest.approx(F_esp, rel=1e-12)


def _correr(ruta):
    subprocess.run([str(FEBIO), "-i", ruta.name, "-silent"], cwd=ruta.parent,
                   stdin=subprocess.DEVNULL, capture_output=True, check=True)
    txt = ruta.with_name(ruta.stem + "_u.txt").read_text()
    ultimo = txt.split("*Step")[-1]
    filas = [l.split() for l in ultimo.splitlines()
             if l.strip() and not l.startswith("*") and "=" not in l]
    a = np.array(filas, float)
    return a[np.argsort(a[:, 0]), 1:].ravel()


@pytest.mark.skipif(not FEBIO.exists(), reason="FEBio 4 no esta instalado")
def test_febio_resuelve_lo_mismo(registro, tmp_path):
    BW, spc, nodos, elems, A = _caso()
    res = ensayo_compresion(BW, spc)
    u = {}
    for s in (1e3, 2e3):
        ruta = tmp_path / f"g{int(s)}.feb"
        escribir_febio(nodos, elems, ruta, sigma_app=s, A_bruta=A)
        u[s] = _correr(ruta) * (1e6 / s)
    u_f = 2 * u[1e3] - u[2e3]                       # extrapolado a carga nula

    z = nodos[:, 2]
    techo = z >= z.max() - 1e-9
    E_f = 1e6 / (abs(u_f[2::3][techo].mean()) / (z.max() - z.min()))
    dE = abs(E_f - res["E_app"]) / res["E_app"]
    du = np.abs(u_f - res["u"]).max() / np.abs(res["u"]).max()

    registro.anotar(BLOQUE, "E_app FEBio (extrapolado) = spinpy", REF,
                    res["E_app"], E_f, dE, "1e-6 relativo", dE < 1e-6)
    registro.anotar(BLOQUE, "desplazamientos FEBio = spinpy", REF, 0.0, du,
                    du, "1e-5 del maximo", du < 1e-5)
    assert dE < 1e-6
    assert du < 1e-5
