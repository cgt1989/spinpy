"""estudio.py: Reproduce las tablas de INFORME.md.

    python comparativa_febio_tet/hex_mejorado/estudio.py cavidad
    python comparativa_febio_tet/hex_mejorado/estudio.py giroide --n 16 20 24 32 40 48 64
    python comparativa_febio_tet/hex_mejorado/estudio.py giroide --n 96 128 --solo-voxel
    python comparativa_febio_tet/hex_mejorado/estudio.py viga

Una linea JSON por corrida en la salida estandar.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(AQUI.parent))
sys.path.insert(0, str(AQUI.parents[1]))
import giroide                                   # noqa: E402
import hex_mejorado as hm                        # noqa: E402

VARIANTES = {
    "hex8": dict(suave=False),
    "hex8I": dict(suave=False, incompatible=True),
    "taubin10": dict(iters=10),
    "taubin20": dict(iters=20),
    "proy": dict(),
    "proyI": dict(incompatible=True),
}


def emitir(d):
    print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v)
                      for k, v in d.items()}), flush=True)


def cavidad(ns, variantes):
    from cavidad import geometria
    for n in ns:
        BW, spc, c, r = geometria(n)
        cc = c + 0.5 * spc[0]            # centro en coordenadas de la malla hex
        f = lambda P: np.linalg.norm(P - cc, axis=1)
        g = lambda P: (P - cc) / np.maximum(
            np.linalg.norm(P - cc, axis=1), 1e-30)[:, None]
        for v in variantes:
            kw = dict(VARIANTES[v])
            if v.startswith("proy"):
                kw["campo"] = (f, g, r)
            emitir({**hm.correr(BW, spc[0], **kw), "caso": "cavidad",
                    "n": n, "variante": v})


def giro(ns, variantes, celdas=3):
    for n in ns:
        BW = giroide.mascara(n, celdas)
        for v in variantes:
            kw = dict(VARIANTES[v])
            if v.startswith("proy"):
                kw.update(campo=giroide.campo(n, celdas), tope=0.5)
            emitir({**hm.correr(BW, 1.0 / n, **kw), "caso": "giroide",
                    "n": n, "variante": v, "BVTV_voxel": float(BW.mean())})


def viga():
    """Voladizo L/t = 10, carga en la punta, frente a Timoshenko."""
    from scipy import sparse
    from scipy.sparse.linalg import spsolve
    from spinpy.resistencia import _matriz_D
    from spinpy.solido import malla_hex
    E_, nu = 1.0, 0.3
    D = _matriz_D(E_, nu)
    for nt in (1, 2, 4):
        L, t = 10.0, 1.0
        h = t / nt
        X, el, _ = malla_hex(np.ones((int(L / h), nt, nt), bool), h)
        for inc in (False, True):
            X0 = np.array(hm.ESQUINAS, float)[None] * h
            ke = (hm.ke_lote_I if inc else hm.ke_lote)(X0, D)[0]
            ne, n = len(el), 3 * len(X)
            ed = (3 * el[:, :, None] + np.arange(3)).reshape(ne, 24)
            K = sparse.coo_matrix((np.tile(ke.ravel(), ne),
                                   (np.repeat(ed, 24, 1).ravel(),
                                    np.tile(ed, (1, 24)).ravel())),
                                  shape=(n, n)).tocsr()
            F = np.zeros(n)
            punta = np.nonzero(np.abs(X[:, 0] - L) < 1e-9)[0]
            F[3 * punta + 2] = -1.0 / punta.size
            emp = np.nonzero(X[:, 0] < 1e-9)[0]
            lib = np.setdiff1d(np.arange(n),
                               np.r_[3 * emp, 3 * emp + 1, 3 * emp + 2])
            u = np.zeros(n)
            u[lib] = spsolve(K[lib][:, lib].tocsc(), F[lib])
            ref = L ** 3 / (3 * E_ * t ** 4 / 12) + \
                L / (5 / 6 * E_ / (2 * (1 + nu)) * t * t)
            emitir({"caso": "viga", "elem_espesor": nt,
                    "variante": "hex8I" if inc else "hex8",
                    "flecha_rel": float(-u[3 * punta + 2].mean() / ref)})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("caso", choices=["cavidad", "giroide", "viga"])
    ap.add_argument("--n", nargs="*", type=int)
    ap.add_argument("--variantes", nargs="*", default=None)
    ap.add_argument("--solo-voxel", action="store_true")
    a = ap.parse_args()
    var = a.variantes or (["hex8", "hex8I"] if a.solo_voxel else
                          ["hex8", "hex8I", "proy", "proyI"])
    if a.caso == "cavidad":
        cavidad(a.n or [16, 24, 32, 40, 48], var)
    elif a.caso == "giroide":
        giro(a.n or [16, 20, 24, 32, 40, 48, 64], var)
    else:
        viga()


if __name__ == "__main__":
    main()
