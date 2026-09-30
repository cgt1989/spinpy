"""estudio_grises.py: ladrillos proyectados sobre los grises del micro-CT.

    python comparativa_febio_tet/hex_mejorado/estudio_grises.py sintetico
    python comparativa_febio_tet/hex_mejorado/estudio_grises.py h4
    python comparativa_febio_tet/hex_mejorado/estudio_grises.py h4 --pila RUTA --n 24 32 48 97

SINTETICO. El giroide de INFORME.md §3 «escaneado» a N = 75 voxeles por lado
(Tb.Th/h nativo ~ 6, el del VOI proximal de H4), con volumen parcial,
desenfoque del sistema y ruido, y segmentado por Otsu. Se malla la mascara
remuestreada a n y se proyecta sobre la isosuperficie de los grises al mismo
umbral. La referencia es la geometria continua (hex8I a 128^3 del giroide,
688,55 MPa): mide cuanto se pierde por pasar de la superficie exacta a la de
una imagen con ruido.

H4. La pila del escaner (por omision ../H4 junto al repositorio, como en
`probar_voi_h4.py`), segmentada con Otsu de 3 clases como la aplicacion, el
cubo del tercio indicado alineado por PCA y, para cada n, las variantes hex8,
hex8I y proyI con y sin filtro gaussiano previo de los grises. Sin geometria
exacta, la referencia es la propia variante a resolucion nativa.

PREDICCIONES ESCRITAS ANTES DE CORRER EL CASO REAL (2026-09-30)
----------------------------------------------------------------
G1  Sintetico sin ruido: proyI sobre los grises queda a <= 3 % de la
    referencia continua desde n = 24, como sobre la superficie exacta.
G2  Sintetico con ruido (SNR 8): sin filtro la proyeccion sigue al ruido y
    empeora; con un filtro de sigma 1 voxel nativo vuelve a <= 3 %.
G3  H4: la dispersion de E_app de proyI entre n = 24 y la nativa es menor
    que la de hex8 (hoy +-3 % desde n = 40 y peor por debajo).
G4  H4: hex8I y proyI a resolucion nativa difieren <= 3 % (a Tb.Th/h ~ 6 la
    superficie ya pesa poco en la rigidez).

Una linea JSON por corrida en la salida estandar y en `resultados/`.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[1]
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(RAIZ))
import giroide                                       # noqa: E402
import grises                                        # noqa: E402
import hex_mejorado as hm                            # noqa: E402
from spinpy import voi                               # noqa: E402
from spinpy.elastic import remuestrear_bw            # noqa: E402

E_REF_GIROIDE = 688.55436      # hex8I a 128^3 sobre el giroide continuo


def _otsu(v, clases=2):
    from skimage.filters import threshold_multiotsu, threshold_otsu
    v = np.asarray(v, np.float64).ravel()
    if clases == 2:
        return float(threshold_otsu(v))
    return float(threshold_multiotsu(v, classes=clases)[-1])


def _emitir(d, fh=None):
    lin = json.dumps({k: (round(v, 6) if isinstance(v, float) else v)
                      for k, v in d.items()})
    print(lin, flush=True)
    if fh:
        fh.write(lin + "\n")
        fh.flush()


def barrido(G, T, M, h0, ns, sigmas, variantes, extra, fh=None, ref=None):
    """Para cada n: remuestrea la mascara, malla y resuelve cada variante."""
    N = M.shape[0]
    for n in ns:
        B, spc = remuestrear_bw(M, np.full(3, h0), n)
        for v in variantes:
            sig = sigmas if v.startswith("proy") else [None]
            for s in sig:
                kw = {"hex8": dict(suave=False),
                      "hex8I": dict(suave=False, incompatible=True),
                      "proy": dict(), "proyI": dict(incompatible=True)}[v]
                if v.startswith("proy"):
                    kw["campo"] = grises.campo_para(G, T, B.shape, spc, s)
                t0 = time.perf_counter()
                o = hm.correr(B, spc[0], **kw)
                d = {**extra, "n": int(B.shape[0]), "N_nativo": int(N),
                     "n_rel": float(B.shape[0]) / N, "variante": v,
                     "sigma_filtro": s, "BVTV_mascara": float(B.mean()),
                     **o, "t_total": time.perf_counter() - t0}
                if ref:
                    d["err_pct"] = 100 * (o["E_app"] / ref - 1)
                _emitir(d, fh)


def sintetico(a):
    cel = 3
    ind = lambda P: giroide.g(*(P * 2 * np.pi * cel + 0.3).T) > giroide.T
    salida = AQUI / "resultados" / "grises_sintetico.jsonl"
    with open(salida, "a", encoding="utf-8") as fh:
        for snr in a.snr:
            G = grises.simular_microct(ind, a.N, psf=a.psf,
                                       snr=(snr if snr > 0 else None))
            T = _otsu(G, 2)
            M = G > T
            extra = {"caso": "sintetico", "snr": snr, "psf": a.psf,
                     "umbral": T, "BVTV_nativo": float(M.mean())}
            barrido(G, T, M, 1.0 / a.N, a.n or [24, 32, 48, a.N],
                    a.sigma, a.variantes, extra, fh, ref=E_REF_GIROIDE)


def h4(a):
    pila = Path(a.pila) if a.pila else RAIZ.parent / "H4"
    if not pila.exists():
        sys.exit(f"No encuentro la pila del escaner en {pila}. "
                 f"Indicala con --pila.")
    vol, apartados = grises.leer_gris(pila, a.patron)
    if a.tam_voxel:
        h0 = float(a.tam_voxel) * voi.UNIDADES_MM[a.unidad]
        origen_h = f"indicado ({a.tam_voxel} {a.unidad})"
    else:
        h0, _log = voi.tam_voxel_desde_log(pila)
        if h0 is None:
            sys.exit("Sin tamano de voxel: usa --tam-voxel.")
        origen_h = "log del escaner"
    sp = np.full(3, h0)
    T = float(a.umbral) if a.umbral is not None else voi._umbral_otsu(vol, 3)
    M = vol > T
    ejes, centro, proy = voi.marco_pca(M, sp)
    centros = voi.centros_por_tercios(proy)
    cubo, fuera = voi.extraer_cubo(M, sp, centros[a.tercio], ejes, centro,
                                   a.lado_mm)
    G0 = grises.extraer_cubo_gris(vol, sp, centros[a.tercio], ejes, centro,
                                  a.lado_mm, orden=0)
    identica = bool(np.array_equal(G0 > T, cubo))
    G = grises.extraer_cubo_gris(vol, sp, centros[a.tercio], ejes, centro,
                                 a.lado_mm, orden=1)
    del vol, M
    extra = {"caso": "h4", "pila": str(pila), "tercio": a.tercio,
             "umbral": T, "h_mm": h0, "h_origen": origen_h,
             "fuera": fuera, "mascara_identica": identica,
             "BVTV_nativo": float(cubo.mean()),
             "BVTV_trilineal": float((G > T).mean()),
             "apartados": len(apartados)}
    if a.vtk_referencia and Path(a.vtk_referencia).exists():
        from spinpy import leer_voi
        R, _ = leer_voi(a.vtk_referencia)
        if R.shape == cubo.shape:
            extra["coincide_vtk"] = float((R == cubo).mean())
    print(json.dumps({k: v for k, v in extra.items()}), flush=True)
    if not identica:
        sys.exit("La mascara en grises no reproduce la de voi.extraer_cubo.")
    N = cubo.shape[0]
    salida = AQUI / "resultados" / "grises_h4.jsonl"
    with open(salida, "a", encoding="utf-8") as fh:
        barrido(G, T, cubo, h0, a.n or [24, 32, 40, 48, 64, N], a.sigma,
                a.variantes, extra, fh)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("caso", choices=["sintetico", "h4"])
    ap.add_argument("--n", nargs="*", type=int)
    ap.add_argument("--sigma", nargs="*", type=float, default=[0.0, 1.0])
    ap.add_argument("--variantes", nargs="*",
                    default=["hex8", "hex8I", "proyI"])
    # sintetico
    ap.add_argument("--N", type=int, default=75)
    ap.add_argument("--psf", type=float, default=0.7)
    ap.add_argument("--snr", nargs="*", type=float, default=[0.0, 8.0])
    # h4
    ap.add_argument("--pila")
    ap.add_argument("--patron", default="*.tif")
    ap.add_argument("--tam-voxel", type=float)
    ap.add_argument("--unidad", default="um")
    ap.add_argument("--umbral", type=float)
    ap.add_argument("--tercio", type=int, default=0)
    ap.add_argument("--lado-mm", type=float, default=5.0)
    ap.add_argument("--vtk-referencia",
                    default=str(RAIZ.parent / "H4" / "Segmentadas"
                                / "VOI_proximal_cubico.vtk"))
    a = ap.parse_args()
    (sintetico if a.caso == "sintetico" else h4)(a)


if __name__ == "__main__":
    main()
