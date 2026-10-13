"""
sanidad_cubo.py: chequeo de sanidad de TODOS los analisis mecanicos con un cubo
macizo, contra soluciones exactas.

    python sanidad_cubo.py            # ~5 min

Un cubo macizo de 3,008 mm de lado (el tamano de los VOIs C1 a C5), con el
protocolo de Tapia: E_s = 18 GPa, nu = 0,3, F = 100 N sobre la seccion bruta.
Tiene solucion exacta, asi que cualquier desviacion es del metodo, no de la
estructura.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  lineal, resolvedor directo       1e-8  relativa
  no lineal (Newton, SVK)          1e-6  relativa frente a la solucion cerrada
  homogeneizacion periodica        1e-8  relativa
  convergencia entre resoluciones  1e-6  relativa
  prisma de lados distintos        1e-8 (un cubo no detecta que se confundan
                                   ejes, areas o alturas)
  rigidizacion por base empotrada  se REPORTA, no se juzga: no hay solucion
                                   cerrada; solo se exigen equilibrio y el
                                   acuerdo entre motores (1e-8)

SOLUCION EXACTA (apoyo deslizante: uz = 0 en la base, resto libre)
  sigma = F/A,  eps_z = -sigma/E,  eps_x = eps_y = nu sigma/E
  uz(z) = -eps z;  ux = eps_x (x - x_min);  uy = eps_y (y - y_min)
  von Mises = sigma en TODO el bloque.
  No lineal (Saint Venant-Kirchhoff, tension uniaxial): E lam (lam^2 - 1) / 2
  = -sigma con lam = 1 - eps_app.

EJES. Base: 'deslizante' fija uz y lo minimo contra el movimiento rigido;
'empotrado' (Tapia) fija ux, uy, uz. Cara superior: x e y libres, fuerza
nodal igual en -z (no es un plato rigido). Caras laterales: libres.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent.parent))

from spinpy import fem, motores                                 # noqa: E402
from spinpy.elastic import constantes_ingenieria, homogeneizar  # noqa: E402
from spinpy.resistencia import (EPS_CRITICA, ensayo_triaxial,   # noqa: E402
                                estudio_convergencia)

L_MM, E_S, NU, F_N = 3.008, 18e9, 0.30, 100.0
N_HEX, N_TET = 16, 8
TOL_LIN, TOL_NL, TOL_HOM, TOL_CONV = 1e-8, 1e-6, 1e-8, 1e-6
FILAS: list[dict] = []


def anotar(grupo, prueba, esperado, obtenido, tol, nota=""):
    """Error relativo |obtenido/esperado - 1| (o absoluto si esperado = 0)."""
    esperado, obtenido = float(esperado), float(obtenido)
    err = (abs(obtenido - esperado) if esperado == 0
           else abs(obtenido / esperado - 1.0))
    ok = bool(err <= tol) if tol is not None else None
    FILAS.append({"grupo": grupo, "prueba": prueba, "esperado": esperado,
                  "obtenido": obtenido, "error": err, "tol": tol, "ok": ok,
                  "nota": nota})
    marca = {True: "OK ", False: "FALLA", None: "info"}[ok]
    print(f"  [{marca}] {prueba:58s} esp {esperado:+.6e}  obt {obtenido:+.6e}"
          f"  err {err:.2e}", flush=True)


def svk_uniaxial(sigma_pa):
    """eps_app exacto del SVK con traccion nominal sigma (compresion) libre lateral."""
    from scipy.optimize import brentq
    lam = brentq(lambda l: E_S * l * (l * l - 1) / 2 + sigma_pa, 0.5, 1.0,
                 xtol=1e-16, rtol=1e-15)
    return 1.0 - lam


def correr(tipo, apoyo, motor, analisis, n):
    h = L_MM / n
    BW = np.ones((n, n, n), bool)
    malla = (fem.mallar_aislado(BW, np.full(3, h), tipo) if tipo != "hex8"
             else fem.mallar(BW, np.full(3, h), tipo))
    prot = fem.protocolo("tapia2026", apoyo=apoyo)
    reg = fem.ensayo(malla, prot, analisis=analisis, motor=motor,
                     solver="auto" if motor == "app" else "directo",
                     aislado=False, BW=BW if motor == "app" else None)
    return malla, prot, reg


def campos(malla, reg):
    u = reg["_campos_lineal"]["u"]
    s = fem.tension_elemental(malla["nodos"], malla["elems"], u, E_S, NU)
    return u, s, fem.von_mises(s)


def chequeos_lineales(tipo, apoyo, motor, n):
    g = f"{tipo} {apoyo} {motor}"
    print(f"\n== {g} (n = {n})", flush=True)
    t0 = time.perf_counter()
    nl = motor != "skfem" and not (motor == "app")
    analisis = (("lineal", "lineal_plato", "nl_fuerza", "nl_pistoia") if nl
                else ("lineal", "lineal_plato"))
    malla, prot, reg = correr(tipo, apoyo, motor, analisis, n)
    A = L_MM ** 2
    sig = 1e6 * F_N / A                                 # Pa
    eps = sig / E_S
    u, s, vm = campos(malla, reg)
    nod = malla["nodos"]
    zt, zb = nod[:, 2].max(), nod[:, 2].min()
    base, techo = nod[:, 2] <= zb + 1e-9, nod[:, 2] >= zt - 1e-9
    lin = reg["lineal"]
    anotar(g, "equilibrio: reaccion = fuerza aplicada", F_N, lin["F_reac_N"],
           TOL_LIN)
    anotar(g, "equilibrio: integral de volumen de szz = fuerza", F_N,
           lin["F_N"], TOL_LIN)
    # equilibrio de SECCION: la tension media de cada capa horizontal es
    # -sigma. Solo con ladrillos: una capa de hex8 ES una seccion horizontal;
    # los tetraedros agrupados por la altura de su centroide no lo son (la
    # primera version de esta prueba lo hizo y fallaba igual con los dos
    # motores, un defecto de la prueba, no del solver). En TET10 rige el
    # equilibrio global de arriba.
    cen = nod[malla["elems"]].mean(axis=1)
    if tipo == "hex8":
        capas = np.unique(np.round(cen[:, 2], 9))
        vol = malla["vol_elem"]
        dev = max(abs((s[cen[:, 2].round(9) == c, 2]
                       * vol[cen[:, 2].round(9) == c]).sum()
                      / vol[cen[:, 2].round(9) == c].sum() / -sig - 1.0)
                  for c in capas)
        anotar(g, "equilibrio de seccion: szz medio por capa = -sigma", 1.0,
               1.0 + dev, TOL_LIN)
    if apoyo == "deslizante":
        anotar(g, "E_app = E_s", E_S, lin["E_app"], TOL_LIN)
        anotar(g, "von Mises uniforme: media = sigma", sig * 1e-6,
               vm.mean() * 1e-6, TOL_LIN)
        anotar(g, "von Mises uniforme: dispersion (max-min)/sigma", 0.0,
               (vm.max() - vm.min()) / sig, TOL_LIN)
        ex = np.sqrt(2 * (NU * eps * L_MM) ** 2 + (eps * L_MM) ** 2)
        anotar(g, "deformacion total |u| maxima (esquina superior)", ex,
               np.linalg.norm(u, axis=1).max(), TOL_LIN)
        anotar(g, "uz del techo = -eps L", -eps * L_MM,
               u[techo, 2].mean(), TOL_LIN)
        esq = techo & (nod[:, 0] >= nod[:, 0].max() - 1e-9)
        anotar(g, "ux del borde libre del techo = nu eps L (eje libre)",
               NU * eps * L_MM, u[esq, 0].max(), TOL_LIN)
        pl = reg.get("lineal_plato") or {}
        if pl and not pl.get("no_disponible"):
            anotar(g, "plato rigido: E_app = E_s", E_S, pl["E_app"], TOL_LIN)
            anotar(g, "plato rigido / fuerza = 1", 1.0,
                   pl["cociente_plato_fuerza"], TOL_LIN)
            c = reg.get("corregido") or {}
            if c.get("E_app"):
                anotar(g, "E corregido (nucleo con plato) = E_s", E_S,
                       c["E_app"], TOL_LIN)
        if nl:
            nf = reg["nl_fuerza"]
            anotar(g, "no lineal a 100 N: E_app frente a SVK cerrado",
                   sig / svk_uniaxial(sig), nf["E_app"], TOL_NL)
            npst = reg["nl_pistoia"]
            sfall = EPS_CRITICA * E_S
            anotar(g, "Pistoia: sigma de fallo = eps_crit E_s", sfall,
                   npst["sigma_Pa"], TOL_LIN)
            anotar(g, "no lineal a la carga de fallo: E_app frente a SVK",
                   sfall / svk_uniaxial(sfall), npst["E_app"], TOL_NL)
    else:
        anotar(g, "base: |ux|, |uy|, |uz| maximos = 0 (ejes fijos)", 0.0,
               np.abs(u[base]).max(), 1e-12)
        e0 = techo & (nod[:, 0] >= nod[:, 0].max() - 1e-9)
        anotar(g, "techo: ux del borde NO es nulo (eje libre), mm",
               NU * eps * L_MM, u[e0, 0].max(), None,
               "se reporta frente al valor del caso deslizante")
        anotar(g, "rigidizacion E_app / E_s - 1 (base empotrada)", 0.0,
               lin["E_app"] / E_S - 1.0, None,
               "no hay solucion cerrada: se reporta")
        anotar(g, "von Mises maximo / sigma (esquinas de la base)", 1.0,
               vm.max() / sig, None, "se reporta")
        c = np.abs(cen - cen.mean(axis=0)).max(axis=1) < 1.5 * (
            L_MM / n)
        anotar(g, "von Mises en el centro / sigma", 1.0,
               vm[c].mean() / sig, None, "se reporta")
    print(f"   ({time.perf_counter() - t0:.0f} s)", flush=True)
    return lin["E_app"], vm.max() * 1e-6, np.abs(u).max()


def homogeneizacion():
    g = "homogeneizacion periodica"
    print(f"\n== {g}", flush=True)
    BW = np.ones((12, 12, 12), bool)
    Ch, info = homogeneizar(BW, E_s=E_S, nu_s=NU, vox_size=L_MM / 12)
    ec = constantes_ingenieria(Ch)
    G = E_S / (2 * (1 + NU))
    for k in ("Ex", "Ey", "Ez"):
        anotar(g, f"{k} = E_s", E_S, ec[k], TOL_HOM)
    for k in ("Gyz", "Gxz", "Gxy"):
        anotar(g, f"{k} = E/(2(1+nu))", G, ec[k], TOL_HOM)
    for k in ("nu_xy", "nu_xz", "nu_yz"):
        anotar(g, f"{k} = nu", NU, ec[k], TOL_HOM)
    lam = E_S * NU / ((1 + NU) * (1 - 2 * NU))
    anotar(g, "C11 = lambda + 2 mu", lam + 2 * G, Ch[0, 0], TOL_HOM)
    anotar(g, "C12 = lambda", lam, Ch[0, 1], TOL_HOM)


def prisma():
    """Prisma de lados distintos: un cubo es simetrico y NO detecta que se
    confundan ejes, areas o alturas (A = Lx Ly, H = Lz). Exacto con apoyo
    deslizante; con base empotrada solo equilibrio."""
    g = "prisma 12x16x24 (lados distintos)"
    print(f"\n== {g}", flush=True)
    h = 0.1
    nx, ny, nz = 12, 16, 24
    Lx, Ly, Lz = nx * h, ny * h, nz * h
    BW = np.ones((nx, ny, nz), bool)
    sig = 1e6 * F_N / (Lx * Ly)                          # Pa
    eps = sig / E_S
    for apoyo in ("deslizante", "empotrado"):
        malla = fem.mallar(BW, np.full(3, h), "hex8")
        prot = fem.protocolo("tapia2026", apoyo=apoyo)
        reg = fem.ensayo(malla, prot, analisis=("lineal",), motor="ngsolve",
                         solver="directo", aislado=False)
        u, s, vm = campos(malla, reg)
        nod = malla["nodos"]
        techo = nod[:, 2] >= nod[:, 2].max() - 1e-9
        lin = reg["lineal"]
        gg = f"{g} {apoyo}"
        anotar(gg, "equilibrio: reaccion = fuerza aplicada", F_N,
               lin["F_reac_N"], TOL_LIN)
        if apoyo == "deslizante":
            anotar(gg, "E_app = E_s (area Lx Ly, altura Lz)", E_S,
                   lin["E_app"], TOL_LIN)
            anotar(gg, "uz del techo = -sigma Lz / E", -eps * Lz,
                   u[techo, 2].mean(), TOL_LIN)
            anotar(gg, "von Mises = sigma = F/(Lx Ly)", sig * 1e-6,
                   vm.mean() * 1e-6, TOL_LIN)
            ex = np.sqrt((NU * eps * Lx) ** 2 + (NU * eps * Ly) ** 2
                         + (eps * Lz) ** 2)
            anotar(gg, "|u| maxima en la esquina superior lejana", ex,
                   np.linalg.norm(u, axis=1).max(), TOL_LIN)
        else:
            anotar(gg, "base: ejes fijos (|u| = 0)", 0.0,
                   np.abs(u[nod[:, 2] <= nod[:, 2].min() + 1e-9]).max(), 1e-12)
            anotar(gg, "rigidizacion E_app / E_s - 1", 0.0,
                   lin["E_app"] / E_S - 1.0, None, "se reporta")
    # compresion en los tres ejes sobre el prisma (area y altura propias)
    r = ensayo_triaxial(BW, np.full(3, h), E_s=E_S, nu_s=NU,
                        apoyo="deslizante", carga_N=100.0)
    for nom in "XYZ":
        anotar(g, f"compresion en {nom}: E_app = E_s", E_S,
               r[nom]["E_app"], TOL_LIN)


def triaxial_y_convergencia():
    g = "app: triaxial y convergencia"
    print(f"\n== {g}", flush=True)
    BW = np.ones((16, 16, 16), bool)
    sp = np.full(3, L_MM / 16)
    r = ensayo_triaxial(BW, sp, E_s=E_S, nu_s=NU, apoyo="deslizante",
                        carga_N=100.0)
    for ax, nom in enumerate("XYZ"):
        anotar(g, f"compresion en {nom}: E_app = E_s", E_S,
               r[nom]["E_app"], TOL_LIN)
    c = estudio_convergencia(BW, sp, resoluciones=(8, 10, 12, 14, 16), E_s=E_S,
                             nu_s=NU, apoyo="deslizante", carga_N=100.0)
    for f in c["puntos"]:
        E = f.get("E_app")
        if E:
            anotar(g, f"convergencia a {f.get('n', '?')}^3: E_app = E_s",
                   E_S, E, TOL_CONV)


def main():
    t0 = time.perf_counter()
    combos = []
    for apoyo in ("deslizante", "empotrado"):
        for motor in ("app", "ngsolve", "skfem"):
            combos.append(("hex8", apoyo, motor, N_HEX))
        for motor in ("ngsolve", "skfem"):
            combos.append(("tet10", apoyo, motor, N_TET))
    cruz = {}
    for tipo, apoyo, motor, n in combos:
        if not motores.puede(motor, tipo):
            print(f"\n== {tipo} {apoyo} {motor}: el motor no resuelve {tipo}")
            continue
        cruz[(tipo, apoyo, motor)] = chequeos_lineales(tipo, apoyo, motor, n)
    # acuerdo entre motores en el caso empotrado (sin solucion cerrada)
    g = "acuerdo entre motores (hex8, empotrado)"
    base = cruz.get(("hex8", "empotrado", "ngsolve"))
    for m in ("app", "skfem"):
        o = cruz.get(("hex8", "empotrado", m))
        if base and o:
            anotar(g, f"E_app {m} frente a NGSolve", base[0], o[0], TOL_LIN)
            anotar(g, f"von Mises maximo {m} frente a NGSolve", base[1],
                   o[1], TOL_LIN)
            anotar(g, f"|u| maximo {m} frente a NGSolve", base[2], o[2],
                   TOL_LIN)
    homogeneizacion()
    triaxial_y_convergencia()
    prisma()
    ok = [f for f in FILAS if f["ok"] is True]
    mal = [f for f in FILAS if f["ok"] is False]
    print(f"\nRESUMEN: {len(ok)} pruebas OK, {len(mal)} fallan, "
          f"{sum(f['ok'] is None for f in FILAS)} informativas "
          f"({time.perf_counter() - t0:.0f} s)")
    for f in mal:
        print("  FALLA:", f["grupo"], "|", f["prueba"], f"err {f['error']:.2e}")
    (AQUI / "resultados").mkdir(exist_ok=True)
    json.dump({"E_s": E_S, "nu": NU, "F_N": F_N, "L_mm": L_MM,
               "filas": FILAS}, open(AQUI / "resultados" / "sanidad_cubo.json",
                                     "w", encoding="utf-8"),
              indent=1, default=float, ensure_ascii=False)
    return 1 if mal else 0


if __name__ == "__main__":
    sys.exit(main())
