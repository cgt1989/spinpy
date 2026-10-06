"""
p6_flexion.py: Hexaedro trilineal (hex8) frente a hexaedro con modos
incompatibles (hex8i) en una viga en voladizo de voxeles.

    python comparativa_motores/libros/p6_flexion.py

Viga de L = 10 mm, seccion cuadrada t x t con t = 1 mm, empotrada en x = 0 y
con una fuerza transversal P = 1 N repartida en la cara x = L. Se resuelve con
nt = 1, 2, 3, 4, 6 y 8 elementos por canto (cubos de lado t/nt) y se compara
la flecha media de la cara libre con la de la malla mas fina de hex8i
extrapolada por Richardson. La viga de un voxel o dos de canto es el caso de
una trabecula en flexion con la resolucion por omision de la app.

Escribe resultados/p6_flexion.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parents[1]))
from spinpy.elastic import hex8_ke                               # noqa: E402

E, NU, L, T, P = 1000.0, 0.3, 10.0, 1.0, 1.0


def voladizo(nt, elemento):
    h = T / nt
    nx, ny, nz = int(round(L / h)), nt, nt
    ke = hex8_ke(h, h, h, E, NU, elemento)
    idx = np.arange((nx + 1) * (ny + 1) * (nz + 1)).reshape(
        nx + 1, ny + 1, nz + 1)
    i, j, k = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz),
                          indexing="ij")
    i, j, k = i.ravel(), j.ravel(), k.ravel()
    nod = np.stack([idx[i, j, k], idx[i + 1, j, k], idx[i + 1, j + 1, k],
                    idx[i, j + 1, k], idx[i, j, k + 1], idx[i + 1, j, k + 1],
                    idx[i + 1, j + 1, k + 1], idx[i, j + 1, k + 1]], 1)
    edof = (3 * nod[:, :, None] + np.arange(3)).reshape(len(nod), 24)
    ndof = 3 * idx.size
    K = sparse.coo_matrix((np.tile(ke.ravel(), len(edof)),
                           (np.repeat(edof, 24, 1).ravel(),
                            np.tile(edof, (1, 24)).ravel())),
                          shape=(ndof, ndof)).tocsr()
    f = np.zeros(ndof)
    libre = idx[-1].ravel()
    # fuerza en z repartida por area tributaria en la cara x = L
    w = np.ones((ny + 1, nz + 1))
    w[0, :] *= 0.5; w[-1, :] *= 0.5; w[:, 0] *= 0.5; w[:, -1] *= 0.5
    f[3 * libre + 2] = -P * (w / w.sum()).ravel()
    fijo = np.concatenate([3 * idx[0].ravel() + c for c in range(3)])
    lib = np.setdiff1d(np.arange(ndof), fijo)
    u = np.zeros(ndof)
    u[lib] = spsolve(K[lib][:, lib].tocsc(), f[lib])
    return float(-(u[3 * libre + 2] * (w / w.sum()).ravel()).sum())


def main():
    nts = (1, 2, 3, 4, 6, 8)
    res = {el: {nt: voladizo(nt, el) for nt in nts} for el in ("hex8", "hex8i")}
    # referencia: Richardson con las dos mallas mas finas de hex8i, orden 2
    a, b = res["hex8i"][6], res["hex8i"][8]
    ref = b + (b - a) / ((8 / 6) ** 2 - 1)
    euler = P * L ** 3 / (3 * E * T ** 4 / 12)
    out = {"referencia_richardson": ref, "euler_bernoulli": euler,
           "flecha": res,
           "error_rel": {el: {nt: v / ref - 1 for nt, v in d.items()}
                         for el, d in res.items()}}
    (AQUI / "resultados").mkdir(exist_ok=True)
    (AQUI / "resultados" / "p6_flexion.json").write_text(
        json.dumps(out, indent=1))
    print(f"referencia {ref:.5f} mm (Euler-Bernoulli {euler:.5f})")
    for nt in nts:
        print(f"nt={nt}:  hex8 {out['error_rel']['hex8'][nt]:+.3%}   "
              f"hex8i {out['error_rel']['hex8i'][nt]:+.3%}")


if __name__ == "__main__":
    main()
