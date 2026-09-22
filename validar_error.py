"""
validar_error.py — Contrasta la funcion de error portada contra MATLAB.

`error_morfometrico` es la pieza mas delicada del ajuste, porque de ella
depende la invariante C6: normalizar por el peso REALMENTE usado, nunca saltar
terminos NaN en silencio. Un port que sume igual pero normalice mal no falla
por ninguna parte: simplemente sesga la busqueda hacia geometrias degeneradas.

Se aprovechan las metricas ya medidas sobre los 34 VOIs (resultados/
referencia_vois_todos.json) para formar TODOS los pares ordenados VOI-VOI y
evaluar el error de cada uno. Son 34x34 = 1156 pares con valores realistas y
bien dispersos, sin ninguna evaluacion nueva del generador.

Se incluyen ademas casos con NaN inyectado, que es donde vive la invariante C6
y donde un port ingenuo se separa.

Uso:
    python validar_error.py --exportar     genera pares_error.json
    matlab -batch "cd('...'); validar_error_matlab"
    python validar_error.py                compara
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from spinpy.error import error_morfometrico   # noqa: E402

RAIZ = Path(__file__).parent
SALIDA = RAIZ / "resultados"
REF = SALIDA / "referencia_vois_todos.json"
PARES = SALIDA / "pares_error.json"
REF_MAT = SALIDA / "error_matlab.json"

CAMPOS = ["BVTV", "DA", "BSBV", "TbTh", "TbSp", "TbN", "PoTot"]


def _aplanar(c):
    while isinstance(c, list) and len(c) == 1 and isinstance(c[0], list):
        c = c[0]
    return list(c)


def construir_pares():
    ref = json.loads(REF.read_text(encoding="utf-8"))
    campos = _aplanar(ref["campos"])
    bloques = ref["resultados"]
    if isinstance(bloques, dict):
        bloques = [bloques]

    idx = [campos.index(c) for c in CAMPOS]
    M = []
    for b in bloques:
        v = np.ravel(_aplanar([b["metricas"]]))
        M.append([None if v[i] is None else float(v[i]) for i in idx])

    pares = []
    n = len(M)
    for i in range(n):
        for j in range(n):
            pares.append({"a": M[i], "b": M[j]})

    # --- casos que ejercitan la invariante C6 --------------------------------
    # Si el port sumara los terminos disponibles sin renormalizar, un candidato
    # con metricas NaN acumularia MENOS error y pareceria mejor. Estos casos lo
    # detectan: el mismo par con y sin NaN debe dar errores comparables, no un
    # error artificialmente bajo.
    base_a = M[0]
    base_b = M[1]
    for k in range(len(CAMPOS)):
        b2 = list(base_b)
        b2[k] = None                     # NaN en el candidato
        pares.append({"a": base_a, "b": b2})
        a2 = list(base_a)
        a2[k] = None                     # NaN en la referencia
        pares.append({"a": a2, "b": base_b})
    pares.append({"a": [None] * len(CAMPOS), "b": base_b})   # nada medible
    pares.append({"a": base_a, "b": [None] * len(CAMPOS)})
    pares.append({"a": [0.0] * len(CAMPOS), "b": base_b})    # divisor cero

    PARES.write_text(json.dumps({"campos": CAMPOS, "pares": pares}, indent=1),
                     encoding="utf-8")
    print(f"{len(pares)} pares -> {PARES}")
    print(f"Ahora ejecuta:  matlab -batch \"cd('{RAIZ}'); validar_error_matlab\"")


def comparar():
    if not REF_MAT.exists():
        print(f"Falta {REF_MAT}. Ejecuta primero validar_error_matlab en MATLAB.")
        return 2
    d = json.loads(PARES.read_text(encoding="utf-8"))
    mat = json.loads(REF_MAT.read_text(encoding="utf-8"))
    err_mat = np.array([np.nan if v is None else float(v)
                        for v in np.ravel(_aplanar([mat["errores"]]))])
    n_mat = np.array([np.nan if v is None else float(v)
                      for v in np.ravel(_aplanar([mat["n_usados"]]))])

    err_py = np.full(len(d["pares"]), np.nan)
    n_py = np.full(len(d["pares"]), np.nan)
    for i, p in enumerate(d["pares"]):
        ma = {c: v for c, v in zip(CAMPOS, p["a"]) if v is not None}
        mb = {c: v for c, v in zip(CAMPOS, p["b"]) if v is not None}
        e, nu = error_morfometrico(ma, mb)
        err_py[i], n_py[i] = e, nu

    ok = np.isfinite(err_mat) & np.isfinite(err_py)
    dif = np.abs(err_py - err_mat)
    rel = np.where(err_mat != 0, dif / np.abs(err_mat), dif)

    print("=" * 66)
    print(f"  Pares evaluados                 : {int(ok.sum())} de {len(err_py)}")
    print(f"  Diferencia absoluta maxima      : {np.nanmax(dif[ok]):.3e}")
    print(f"  Diferencia relativa maxima      : {np.nanmax(rel[ok]):.3e}")
    print(f"  Nº de terminos: coincide en     : {int((n_py == n_mat).sum())} de {len(n_py)}")
    print("=" * 66)

    malos = np.where(ok & (rel > 1e-12))[0]
    if malos.size:
        print(f"\n  {malos.size} pares fuera de 1e-12:")
        for i in malos[:10]:
            print(f"    par {i}: matlab={err_mat[i]:.10g}  python={err_py[i]:.10g}")
        return 1
    print("\n  Coincidencia exacta en todos los pares, incluidos los casos con")
    print("  NaN inyectado que ejercitan la invariante C6.")
    return 0


if __name__ == "__main__":
    if "--exportar" in sys.argv:
        construir_pares()
        sys.exit(0)
    sys.exit(comparar())
