"""diagnostico_mascara.py: a 32^3, que parte del error viene de la mascara gruesa?

Micro-CT sintetico del giroide (N = 75, Otsu), malla a n = 24, 32, 48 con
cinco combinaciones de mascara (que puntos se muestrean) y superficie de
proyeccion (grises o exacta). Salida: resultados/diagnostico_mascara.log.
"""
from pathlib import Path
import sys, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import giroide, grises, hex_mejorado as hm
from spinpy.elastic import remuestrear_bw
REF = 688.55436
N = 75
ind = lambda P: giroide.g(*(P * 6 * np.pi + 0.3).T) > giroide.T
G = grises.simular_microct(ind, N, psf=0.7, snr=None)
from skimage.filters import threshold_otsu
T = float(threshold_otsu(G.ravel())); M = G > T
for n in (24, 32, 48):
    B, spc = remuestrear_bw(M, np.full(3, 1 / N), n)
    j = np.arange(n); P = np.stack(np.meshgrid(j, j, j, indexing="ij"), -1).reshape(-1, 3)
    # posiciones fisicas de los voxeles nativos elegidos por remuestrear_bw
    idx = np.round(np.linspace(0, N - 1, n))
    Bc = ind(np.stack(np.meshgrid((idx + .5) / N, (idx + .5) / N, (idx + .5) / N, indexing="ij"), -1).reshape(-1, 3)).reshape(n, n, n)
    Bd = giroide.mascara(n, 3)
    # grises muestreados (trilineal) en los centros de la malla y umbral
    f, gr, _ = grises.campo_para(G, T, n, spc)
    Bg = (f((P + .5) * spc) > T).reshape(n, n, n)
    # campo exacto en coordenadas de la malla con el MISMO estiramiento que el remuestreo
    s = (N - 1) / (n - 1)
    fe, ge, Te = giroide.campo(n, 3)
    fx = lambda X: fe(((X / spc - .5) * s + .5) / N)
    gx = lambda X: ge(((X / spc - .5) * s + .5) / N) * (s / spc / N)
    for et, BW, campo in [("masc_remuestreo+grises", B, (f, gr, T)),
                          ("masc_remuestreo+exacta", B, (fx, gx, Te)),
                          ("masc_continuo_mismos_pts+exacta", Bc, (fx, gx, Te)),
                          ("masc_grises_trilineal+grises", Bg, (f, gr, T)),
                          ("masc_directa(estudio1)+exacta", Bd, giroide.campo(n, 3))]:
        o = hm.correr(BW, 1 / n if "directa" in et else spc[0], campo=campo, incompatible=True)
        print(n, et, "BV/TV", round(BW.mean(), 4), "err", round(100 * (o["E_app"] / REF - 1), 1), "dV", round(o["dV_pct"], 2), flush=True)
    print(n, "difiere remuestreo vs continuo mismos puntos:", int((B != Bc).sum()), "vs directa:", int((B != Bd).sum()))
