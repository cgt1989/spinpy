"""calidad_tet10.py: Calidad de las mallas TET10 del espinodoide (Tabla 3 del informe).

Relacion radio-arista, diedro minimo y elementos diminutos, sobre las cuatro
esquinas de cada tetraedro. Escribe resultados/calidad_tet10.json.
"""
import sys, json, numpy as np
from pathlib import Path
AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI)); sys.path.insert(0, str(AQUI.parent))
import casos

def metricas(nodos, elems):
    P = nodos[elems[:, :4]]                      # (M,4,3)
    a, b, c, d = P[:, 0], P[:, 1], P[:, 2], P[:, 3]
    aristas = np.stack([b-a, c-a, d-a, c-b, d-b, d-c], 1)
    L = np.linalg.norm(aristas, axis=2)
    V = np.abs(np.einsum('ij,ij->i', b-a, np.cross(c-a, d-a))) / 6
    # circunradio
    u, v, w = b-a, c-a, d-a
    uu, vv, ww = (u*u).sum(1), (v*v).sum(1), (w*w).sum(1)
    cr = np.linalg.norm(uu[:, None]*np.cross(v, w) + vv[:, None]*np.cross(w, u)
                        + ww[:, None]*np.cross(u, v), axis=1) / (12*V)
    re = cr / L.min(1)                           # relacion radio-arista (>= 0.612)
    # angulo diedro minimo
    caras = [(0,1,2,3),(0,1,3,2),(0,2,3,1),(1,2,3,0)]
    n = []
    for i, j, k, o in caras:
        nn = np.cross(P[:, j]-P[:, i], P[:, k]-P[:, i])
        s = np.sign(np.einsum('ij,ij->i', nn, P[:, o]-P[:, i]))
        n.append(-nn * s[:, None] / np.linalg.norm(nn, axis=1)[:, None])
    dmin = np.full(len(elems), 180.0)
    for x in range(4):
        for y in range(x+1, 4):
            ang = 180 - np.degrees(np.arccos(np.clip((n[x]*n[y]).sum(1), -1, 1)))
            dmin = np.minimum(dmin, ang)
    q = lambda z, p: float(np.percentile(z, p))
    return {"n": int(len(elems)),
            "radio_arista": {"mediana": q(re, 50), "p95": q(re, 95), "p99": q(re, 99), "max": float(re.max()),
                             "frac_gt_2": float((re > 2).mean()), "frac_gt_10": float((re > 10).mean())},
            "diedro_min_grados": {"min": float(dmin.min()), "p1": q(dmin, 1), "p5": q(dmin, 5), "mediana": q(dmin, 50),
                                  "frac_lt_5": float((dmin < 5).mean()), "frac_lt_1": float((dmin < 1).mean())},
            "vol_rel_mediana": {"min": float(V.min()/np.median(V)),
                                "n_lt_1e-6": int((V < 1e-6*np.median(V)).sum())}}

out = {}
for n in (32, 48):
    c = casos.espinodoide_tet(n)
    out[f"espinodoide_tet_n{n}"] = metricas(c["nodos"], c["elems"])
    print(n, json.dumps(out[f"espinodoide_tet_n{n}"]), flush=True)
json.dump(out, open(AQUI / 'resultados' / 'calidad_tet10.json', 'w'), indent=1)
