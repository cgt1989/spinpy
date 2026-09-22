"""
validar_elastico.py — Verifica la homogeneizacion portada.

Dos verificaciones independientes, en este orden:

  1. CONTRA LA SOLUCION ANALITICA (no necesita MATLAB).
     Para un laminado de capas planas el promedio de Backus (1962) da el tensor
     EXACTO, y la solucion de elementos finitos debe reproducirlo: el campo
     corrector es constante dentro de cada capa y el hexaedro trilineal lo
     representa sin error de discretizacion.

     Se usan contrastes moderados entre las dos fases ademas del vacio a 1e-6.
     Con vacio puro casi todas las entradas del tensor exacto valen ~0 y la
     comparacion no distinguiria una implementacion correcta de una
     equivocada; con contraste moderado todas las entradas son informativas.

  2. CONTRA MATLAB sobre estructuras reales.
     La homogeneizacion es DETERMINISTA —no interviene el generador aleatorio—
     asi que aqui si cabe exigir coincidencia numerica, como en la morfometria
     sobre VOIs.

Uso:
    python validar_elastico.py --exportar   genera casos_elastico.mat
    matlab -batch "cd('...'); validar_elastico_matlab"
    python validar_elastico.py              compara
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

sys.path.insert(0, str(Path(__file__).parent))
from spinpy import generar_mascara                                  # noqa: E402
from spinpy.elastic import (backus_laminado, constantes_ingenieria,  # noqa: E402
                            homogeneizar)

RAIZ = Path(__file__).parent
SALIDA = RAIZ / "resultados"
CASOS = SALIDA / "casos_elastico.mat"
REF_MAT = SALIDA / "elastico_matlab.json"

E_S, NU_S = 1.0, 0.3


def construir_casos():
    """Estructuras deterministas de prueba, guardadas para que MATLAB las lea."""
    casos = {}

    N = 16
    lam = np.zeros((N, N, N), bool); lam[:, :, :N // 2] = True
    casos["laminado"] = lam

    N = 16
    c = (np.arange(N) + 0.5) / N - 0.5
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    casos["esfera"] = (X**2 + Y**2 + Z**2) <= 0.32**2

    barras = np.zeros((N, N, N), bool)
    paso = N // 2
    for i in range(0, N, paso):
        barras[i:i + 3, :, :] = True
        barras[:, i:i + 3, :] = True
    casos["barras"] = barras

    BW, _, _ = generar_mascara(16, 10 * np.pi, 300, [15, 15, 45], 0.45, seed=7)
    casos["spinodoide"] = BW

    d = {k: v.astype(np.uint8) for k, v in casos.items()}
    d["nombres"] = np.array(list(casos.keys()), dtype=object)
    d["E_s"] = E_S
    d["nu_s"] = NU_S
    d["voxSize"] = 1.0 / 16
    savemat(str(CASOS), d)
    print(f"{len(casos)} casos -> {CASOS}")
    for k, v in casos.items():
        print(f"   {k:12s} {v.shape}  fraccion solida = {v.mean():.4f}")
    print(f"\nAhora:  matlab -batch \"cd('{RAIZ}'); validar_elastico_matlab\"")
    return casos


def parte_analitica():
    print("=" * 78)
    print(" 1. CONTRA LA SOLUCION ANALITICA (promedio de Backus)")
    print("=" * 78)
    print(f"   {'caso':30s} {'solver':13s} {'dif.rel.max':>12s} {'|ceros|':>10s}")
    peor = 0.0
    filas = []
    for N, esp, ev in [(12, 6, 0.3), (16, 8, 0.3), (16, 4, 0.1), (12, 6, 1e-6)]:
        BW = np.zeros((N, N, N), bool); BW[:, :, :esp] = True
        Ch, info = homogeneizar(BW, E_S, NU_S, vox_size=1.0 / N, escala_vacio=ev)
        f = esp / N
        Cex = backus_laminado([f, 1 - f], [E_S, ev * E_S], [NU_S, NU_S])
        m = np.abs(Cex) > 1e-9
        rel = float((np.abs(Ch - Cex)[m] / np.abs(Cex)[m]).max())
        cero = float(np.abs(Ch[~m]).max()) if (~m).any() else 0.0
        # El contraste 1e6 empeora el condicionamiento en seis ordenes, asi que
        # perder seis cifras respecto al caso moderado es lo esperado.
        tol = 1e-7 if ev < 1e-3 else 1e-12
        ok = rel <= tol
        peor = max(peor, rel / tol)
        etq = f"{N}^3  f={f:.2f}  fase2={ev:g}"
        print(f"   [{'OK  ' if ok else 'FALLA'}] {etq:24s} {info['solver']:13s} "
              f"{rel:12.3e} {cero:10.2e}")
        filas.append({"caso": etq, "dif_rel_max": rel, "ceros": cero,
                      "tol": tol, "ok": bool(ok)})
    return filas


def parte_matlab():
    print()
    print("=" * 78)
    print(" 2. CONTRA MATLAB sobre estructuras reales (determinista)")
    print("=" * 78)
    if not REF_MAT.exists():
        print(f"   Falta {REF_MAT.name}: ejecuta validar_elastico_matlab en MATLAB.")
        return None

    mat = json.loads(REF_MAT.read_text(encoding="utf-8"))
    d = loadmat(str(CASOS))
    filas = []
    print(f"   {'caso':14s} {'rho':>7s} {'solver':13s} {'dif.rel.max C':>14s} "
          f"{'dif Ez':>10s}")
    for r in mat["resultados"]:
        nombre = r["nombre"] if isinstance(r["nombre"], str) else r["nombre"][0]
        BW = np.asarray(d[nombre]).astype(bool)
        Cm = np.asarray(r["C"], dtype=float).reshape(6, 6)
        Cp, info = homogeneizar(BW, E_S, NU_S, vox_size=float(d["voxSize"]))
        esc = np.abs(Cm).max()
        rel = float(np.abs(Cp - Cm).max() / esc)
        ecm = constantes_ingenieria(Cm)
        ecp = constantes_ingenieria(Cp)
        dez = abs(ecp["Ez"] - ecm["Ez"]) / abs(ecm["Ez"])
        print(f"   {nombre:14s} {BW.mean():7.4f} {info['solver']:13s} "
              f"{rel:14.3e} {dez:10.3e}")
        filas.append({"caso": nombre, "rho": float(BW.mean()),
                      "dif_rel_max": rel, "dif_Ez": float(dez),
                      "Ez_matlab": ecm["Ez"], "Ez_python": ecp["Ez"]})
    return filas


def main():
    if "--exportar" in sys.argv:
        construir_casos()
        return 0
    a = parte_analitica()
    b = parte_matlab()
    res = {"analitica": a, "matlab": b}
    (SALIDA / "validacion_elastica.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n   Resultados -> {SALIDA / 'validacion_elastica.json'}")
    fallos = sum(1 for f in a if not f["ok"])
    return 1 if fallos else 0


if __name__ == "__main__":
    sys.exit(main())
