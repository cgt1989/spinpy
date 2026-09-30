"""diagnostico_mayoria.py: mascara gruesa por fraccion de volumen (mayoria) frente a vecino mas proximo.

Umbral ideal 0,5, sin y con ruido (SNR 8, filtro sigma 1). Salida:
resultados/diagnostico_mayoria.log.
"""
from pathlib import Path
import sys, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import giroide, grises, hex_mejorado as hm
from spinpy.elastic import remuestrear_bw
REF = 688.55436
N = 75
ind = lambda P: giroide.g(*(P * 6 * np.pi + 0.3).T) > giroide.T
for snr in (None, 8.0):
    G = grises.simular_microct(ind, N, psf=0.7, snr=snr)
    T = 0.5; M = G > T
    for n in (20, 24, 32, 40, 48):
        B, spc = remuestrear_bw(M, np.full(3, 1 / N), n)
        f, gr, _ = grises.campo_para(G, T, n, spc, sigma=1.0 if snr else 0.0)
        # mascara por mayoria: fraccion de 4^3 submuestras con f > T
        j = np.arange(n); P = np.stack(np.meshgrid(j, j, j, indexing="ij"), -1).reshape(-1, 3).astype(float)
        o = (np.arange(4) + .5) / 4; frac = np.zeros(len(P))
        for a in o:
            for b in o:
                for c in o:
                    frac += f((P + [a, b, c]) * spc) > T
        Bm = (frac / 64 > 0.5).reshape(n, n, n)
        for et, BW, kw in [("hex8", B, dict(suave=False)), ("hex8I", B, dict(suave=False, incompatible=True)),
                           ("proyI_nearest", B, dict(campo=(f, gr, T), incompatible=True)),
                           ("hex8I_mayoria", Bm, dict(suave=False, incompatible=True)),
                           ("proyI_mayoria", Bm, dict(campo=(f, gr, T), incompatible=True))]:
            r = hm.correr(BW, spc[0], **kw)
            print(snr, n, et, "BV/TV", round(BW.mean(), 4), "err", round(100 * (r["E_app"] / REF - 1), 1), "dV", round(r["dV_pct"], 2), flush=True)
