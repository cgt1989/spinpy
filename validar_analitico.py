"""
validar_analitico.py — Verificacion del port contra soluciones manufacturadas.

Mismo espiritu que los experimentos M1-M5 de Validacion_Anexo: geometrias cuyo
volumen, area y anisotropia se conocen en forma cerrada, para comprobar que la
morfometria portada mide lo que dice medir SIN necesidad de MATLAB.

Se incluye ademas la comparacion de los dos esquemas de muestreo de ondas
(rechazo vs equitativo), que es la discrepancia de fondo entre AppFinal_V2 y el
repositorio TPMS-Scaffolds-generator.

Uso:  python validar_analitico.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from spinpy import morfometria, tensor_mil, wave_directions   # noqa: E402

RAIZ = Path(__file__).parent
SALIDA = RAIZ / "resultados"
SALIDA.mkdir(exist_ok=True)

registros = []
_fallos = 0


def check(nombre, obtenido, esperado, tol_rel, nota=""):
    """Compara y registra. tol_rel es tolerancia RELATIVA."""
    global _fallos
    if esperado == 0:
        err = abs(obtenido)
    else:
        err = abs(obtenido - esperado) / abs(esperado)
    ok = err <= tol_rel
    if not ok:
        _fallos += 1
    registros.append({"prueba": nombre, "obtenido": float(obtenido),
                      "esperado": float(esperado), "error_rel": float(err),
                      "tol_rel": float(tol_rel), "ok": bool(ok), "nota": nota})
    marca = "OK  " if ok else "FALLA"
    print(f"  [{marca}] {nombre:38s} obt={obtenido:>12.6g} "
          f"esp={esperado:>12.6g} err={100*err:>7.3f}%  (tol {100*tol_rel:.1f}%)")
    return ok


def rejilla(n, h):
    """Centros de voxel de una rejilla n^3 con paso h, centrada en el origen."""
    c = (np.arange(n) + 0.5) * h - n * h / 2.0
    return np.meshgrid(c, c, c, indexing="ij")


# ===========================================================================
print("=" * 78)
print(" A. ESFERA  —  BV exacto, BS exacto, isotropia")
print("=" * 78)

N, H = 120, 1.0 / 120
X, Y, Z = rejilla(N, H)
r = 0.30
BW = (X**2 + Y**2 + Z**2) <= r**2
m = morfometria(BW, [H, H, H])

TV = (N * H) ** 3
BV_ex = 4.0 / 3.0 * np.pi * r**3
BS_ex = 4.0 * np.pi * r**2

check("esfera BV/TV", m["BVTV"], BV_ex / TV, 0.01)
check("esfera BV", m["BV"], BV_ex, 0.01)
# Marching cubes sobre datos BINARIOS sobreestima el area (escalonado). Es un
# sesgo conocido y sistematico: por eso la app mide VOI y candidato con el
# MISMO metodo, para que se cancele al restar (correccion C4).
check("esfera BS", m["BS"], BS_ex, 0.10,
      "MC binario sobreestima area; sesgo comun a ambos lados de la comparacion")
check("esfera Tb.Th = 2BV/BS", m["TbTh"], 2 * m["BV"] / m["BS"], 1e-12,
      "identidad de Parfitt, no 4BV/BS")
check("esfera DA ~ 1 (isotropa)", m["DA"], 1.0, 0.10,
      "suelo de ruido del MIL: DA ~ 1.07")

# ===========================================================================
print()
print("=" * 78)
print(" B. LAMINAS PERPENDICULARES A Z  —  BV/TV y BS exactos, DA degenerado")
print("=" * 78)

N, H = 96, 1.0 / 96
BW = np.zeros((N, N, N), dtype=bool)
periodo, espesor, desfase = 12, 4, 2   # 4 de cada 12 planos son solidos
for k in range(N):
    if ((k - desfase) % periodo) < espesor and k >= desfase:
        BW[:, :, k] = True
m = morfometria(BW, [H, H, H])

bvtv_ex = espesor / periodo
n_interfaces = 2 * (N // periodo)     # dos caras por lamina

# Valor esperado CONSISTENTE CON MARCHING CUBES, no el area geometrica ideal.
# Dos efectos de borde, ambos deliberados y compartidos con MATLAB:
#   (1) la malla no se cierra en la frontera del array, asi que una lamina que
#       tocase k=0 perderia esa cara. El desfase=2 evita el caso, y las 16
#       caras quedan en el interior.
#   (2) la superficie se extiende entre los centros de voxel extremos, de modo
#       que su extension lateral es (N-1)*H y no N*H.
# Sin (2) el area "ideal" seria 16.0 y se leeria un error del 2% que en
# realidad no existe.
lado_mc = (N - 1) * H
BS_ex = n_interfaces * lado_mc**2

check("laminas BV/TV", m["BVTV"], bvtv_ex, 1e-12)
check("laminas BS (consistente con MC)", m["BS"], BS_ex, 1e-6)
check("laminas Tb.Th = 2BV/BS", m["TbTh"], 2 * m["BV"] / BS_ex, 1e-9)

# Sesgo fisico residual: Tb.Th verdadero de un apilado de placas es el espesor.
# La diferencia es el efecto (2) de arriba, y es el mismo en el VOI y en el
# candidato, de modo que se cancela en gran medida al restar (correccion C4).
sesgo = (m["TbTh"] - espesor * H) / (espesor * H)
print(f"       Tb.Th medido={m['TbTh']:.6g}  espesor real={espesor*H:.6g}  "
      f"sesgo de borde={100*sesgo:+.2f}%  (esperado +{100*((N/(N-1))**2-1):.2f}%)")
print(f"       DA = {m['DA']:.4g}   clamped = {m['MIL_clamped']}   "
      f"dir_principal = {np.round(m['dir_principal'], 3)}")
if m["DA"] > 1.5:
    print("  [OK  ] laminas DA >> 1 y direccion principal en el plano de placas")
else:
    print("  [FALLA] laminas: DA no refleja la anisotropia laminar")
    _fallos += 1

# ===========================================================================
print()
print("=" * 78)
print(" C. BARRAS PARALELAS A Z  —  columnar, direccion principal = Z")
print("=" * 78)

N, H = 96, 1.0 / 96
X, Y, Z = rejilla(N, H)
paso, rad = 1.0 / 4, 0.055
xm = np.mod(X + 0.5, paso) - paso / 2
ym = np.mod(Y + 0.5, paso) - paso / 2
BW = (xm**2 + ym**2) <= rad**2
m = morfometria(BW, [H, H, H])

n_barras = int(round(1.0 / paso)) ** 2
bvtv_ex = n_barras * np.pi * rad**2 * 1.0 / 1.0
check("barras BV/TV", m["BVTV"], bvtv_ex, 0.03)
check("barras BS", m["BS"], n_barras * 2 * np.pi * rad * 1.0, 0.10,
      "MC binario sobreestima; se excluyen las tapas por ser el borde abierto")

dz = abs(float(m["dir_principal"][2]))
print(f"       DA = {m['DA']:.4g}   |dir_principal . z| = {dz:.4f}")
if dz > 0.95:
    print("  [OK  ] barras: la direccion de MIL maximo es el eje de las barras")
else:
    print("  [FALLA] barras: direccion principal mal orientada")
    _fallos += 1

# ===========================================================================
print()
print("=" * 78)
print(" D. FRACCION PORTANTE  —  cilindro pasante + isla suelta")
print("=" * 78)

N, H = 48, 1.0 / 48
X, Y, Z = rejilla(N, H)
pasante = (X**2 + Y**2) <= 0.12**2                 # cilindro que atraviesa Z
isla = ((X - 0.3) ** 2 + (Y - 0.3) ** 2 + Z**2) <= 0.06**2
BW = pasante | isla
m = morfometria(BW, [H, H, H], do_mil=False)

n_pas = int(pasante.sum())
n_tot = int(BW.sum())
check("FracPort con isla suelta", m["FracPort"], n_pas / n_tot, 0.02,
      "solo cuenta el material en componentes que unen z=0 con z=fin")

m2 = morfometria(pasante, [H, H, H], do_mil=False)
check("FracPort sin isla", m2["FracPort"], 1.0, 1e-12)

# ===========================================================================
print()
print("=" * 78)
print(" E. MIL — invariancia del DA frente al remuestreo y factor F5")
print("=" * 78)

N, H = 64, 1.0 / 64
X, Y, Z = rejilla(N, H)
BW = (np.sin(2 * np.pi * X / 0.25) * np.sin(2 * np.pi * Y / 0.25)
      + np.sin(2 * np.pi * Z / 0.5)) > 0.2
fab = tensor_mil(BW, [H, H, H])

# MIL en unidades de longitud: debe ser del orden del tamano de trabecula, no
# de la mitad. Con el bug F5 (dividir por totalInt en vez de 0.5*totalInt) los
# valores salen exactamente a la mitad.
mil_med = float(np.nanmean(fab["MIL"]))
bvtv = float(BW.mean())
print(f"       BV/TV = {bvtv:.4f}   MIL medio = {mil_med:.5f} mm   "
      f"DA = {fab['DA']:.4g}")
if 0.5 * H < mil_med < 20 * H:
    print("  [OK  ] MIL medio en escala de longitud plausible (factor F5 aplicado)")
else:
    print("  [FALLA] MIL medio fuera de escala: revisar el factor 0.5*totalInt")
    _fallos += 1

# El DA es un cociente de autovalores: escalar el spacing no debe cambiarlo.
fab2 = tensor_mil(BW, [2 * H, 2 * H, 2 * H])
check("DA invariante al escalado isotropo", fab2["DA"], fab["DA"], 1e-9)
check("MIL escala linealmente con spacing",
      float(np.nanmean(fab2["MIL"])), 2 * mil_med, 1e-9)

# ===========================================================================
print()
print("=" * 78)
print(" F. ESQUEMAS DE MUESTREO — GIBBON (rechazo) vs profesor (equitativo)")
print("=" * 78)
print("   Fraccion de ondas cuyo angulo minimo cae bajo cada theta activo.")
print()

comparacion = []
for thetas in [(15, 45, 0), (15, 15, 0), (0, 0, 60), (30, 60, 90), (15, 90, 0)]:
    fila = {"thetas": list(thetas)}
    for esquema in ("rechazo", "equitativo"):
        rng = np.random.default_rng(20260720)
        W = wave_directions(40000, thetas, esquema=esquema, rng=rng)
        ang = np.degrees(np.arccos(np.clip(np.abs(W @ np.eye(3)), -1, 1)))
        rep = [float(np.mean(ang[:, i] < t)) if t > 0 else 0.0
               for i, t in enumerate(thetas)]
        fila[esquema] = rep
    dif = max(abs(a - b) for a, b in zip(fila["rechazo"], fila["equitativo"]))
    fila["dif_max"] = float(dif)
    comparacion.append(fila)
    print(f"   thetas={str(thetas):14s}")
    print(f"       rechazo (GIBBON)    = {[round(v,3) for v in fila['rechazo']]}")
    print(f"       equitativo (prof.)  = {[round(v,3) for v in fila['equitativo']]}"
          f"     dif_max = {dif:.3f}")

print()
print("   Los dos esquemas COINCIDEN solo cuando los conos activos tienen el")
print("   mismo angulo. Con angulos distintos reparten las ondas de forma")
print("   diferente y producen anisotropias distintas para los mismos thetas.")

# ===========================================================================
print()
print("=" * 78)
res = {"pruebas": registros, "muestreo": comparacion,
       "n_fallos": _fallos, "n_pruebas": len(registros)}
(SALIDA / "validacion_analitica.json").write_text(
    json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
print(f" {len(registros)} comprobaciones numericas, {_fallos} fallo(s).")
print(f" Resultados -> {SALIDA / 'validacion_analitica.json'}")
print("=" * 78)
sys.exit(1 if _fallos else 0)
