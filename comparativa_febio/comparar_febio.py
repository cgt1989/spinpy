"""comparar_febio.py — El ensayo de compresion de spinpy frente a FEBio 4.

    cd Port_Python
    python comparativa_febio/comparar_febio.py              # todos los casos
    python comparativa_febio/comparar_febio.py --casos bloque proximal_32

Para cada caso resuelve `resistencia.ensayo_compresion` en spinpy, escribe el
MISMO problema con `escribe.escribir_febio`, lo corre en `febio4`, lee sus
desplazamientos, reacciones y tensiones, y compara. Las tolerancias estan en
`PREDICCIONES.md`, escritas antes de la primera corrida.

DOS DECISIONES QUE HAY QUE SABER PARA LEER LOS NUMEROS
------------------------------------------------------
* FEBio es NO LINEAL geometricamente (su `isotropic elastic` es St.
  Venant-Kirchhoff en grandes desplazamientos) y spinpy es lineal. Para
  contrastar el problema lineal se corre FEBio a 1 kPa y a 2 kPa y se
  extrapola a carga nula, u_lin = 2 u(1 kPa) - u(2 kPa) por unidad de carga:
  el desvio no lineal es proporcional a la carga y asi se cancela a primer
  orden. PREDICCIONES.md suponia que 1 kPa bastaba (el bloque macizo daba
  3e-8); en el VOI proximal a 32^3 dio 3e-4, por voladizos del techo.
* Se corre ademas a 1 MPa, la carga de referencia de la app, y se reporta
  cuanto se aparta la respuesta no lineal de la lineal: eso no valida ni
  invalida a spinpy, mide si la hipotesis lineal se sostiene.
* Los estadisticos (percentil de superficie, Pistoia) se calculan con las
  MISMAS funciones de spinpy sobre los dos campos. Lo que se contrasta es la
  solucion de elementos finitos, no la implementacion de un percentil.

Resultados: una linea por caso en `resultados/febio.jsonl`, con el bloque de
procedencia. Los .feb y las salidas de FEBio quedan en `corridas/<caso>/`
(no se versionan: se regeneran).
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))
from spinpy import __version__ as VERSION_SPINPY          # noqa: E402
from spinpy import generar_mascara, leer_voi, procedencia  # noqa: E402
from spinpy.elastic import remuestrear_bw                  # noqa: E402
from spinpy.escribe import escribir_febio                  # noqa: E402
from spinpy.resistencia import (_solo_portante, criterio_pistoia,  # noqa: E402
                                ensayo_compresion)
from spinpy.solido import malla_hex                        # noqa: E402

FEBIO = Path(r"C:\Program Files\FEBioStudio2\bin\febio4.exe")
VOIS = AQUI.parent.parent / "H4" / "Segmentadas"
E_S, NU_S = 20e9, 0.30
SIGMA_SPINPY = 1e6         # el valor por omision de la app
SIGMA_LIN = (1e3, 2e3)     # dos cargas pequenas para extrapolar a cero

# El espinodoide ajustado al VOI proximal cubico de H4. VERBATIM de
# Estudio_Percolacion/code/percolacion.py (AJUSTADO, SEMILLA, LC_VOI); el
# numero de onda es el del deslizador y se multiplica por pi, como fit.py.
AJUSTADO = dict(rho=0.30847554, wave_pi=15.0, num_waves=800,
                thetas=(15.0, 15.0, 60.0), seed=20260720, Lc=4.994433)


# ---------------------------------------------------------------------------
# Casos
# ---------------------------------------------------------------------------

def _voi(nombre, n):
    BW, spc = leer_voi(VOIS / f"VOI_{nombre}_cubico.vtk")
    BW, spc = remuestrear_bw(BW, spc, n)
    return BW, spc


def _espinodoide(n):
    a = AJUSTADO
    BW, _, _ = generar_mascara(resolution=n, wave_number=a["wave_pi"] * np.pi,
                               num_waves=a["num_waves"],
                               thetas=list(a["thetas"]), rho=a["rho"],
                               seed=a["seed"])
    return np.asarray(BW, bool), np.full(3, a["Lc"] / n)


CASOS = {
    "bloque": lambda: (np.ones((8, 8, 12), bool), np.full(3, 0.05),
                       "deslizante", "solido"),
    **{f"{v}_{n}": (lambda v=v, n=n: (*_voi(v, n), "deslizante", "VOI"))
       for n in (32, 48) for v in ("proximal", "medio", "distal")},
    "espinodoide_48": lambda: (*_espinodoide(48), "deslizante", "spinodoide"),
    "proximal_48_empotrado": lambda: (*_voi("proximal", 48), "empotrado",
                                      "VOI"),
}


# ---------------------------------------------------------------------------
# FEBio
# ---------------------------------------------------------------------------

def _ultimo_registro(ruta, ncol):
    """Ultimo bloque de datos de un archivo de logfile de FEBio.

    FEBio escribe un bloque por paso convergido, INCLUIDO el paso 0 (todo
    ceros). Se toma el ultimo y se ordena por el identificador, que es la
    primera columna.
    """
    bloques = re.split(r"^\*Step.*$", Path(ruta).read_text(), flags=re.M)
    filas = [l.split() for l in bloques[-1].splitlines()
             if l.strip() and not l.startswith("*")]
    a = np.array(filas, float)
    if a.shape[1] != ncol + 1:
        raise RuntimeError(f"{ruta}: {a.shape[1]} columnas, se esperaban "
                           f"{ncol + 1}")
    return a[np.argsort(a[:, 0]), 1:]


def correr_febio(feb):
    t0 = time.perf_counter()
    r = subprocess.run([str(FEBIO), "-i", feb.name, "-silent"], cwd=feb.parent,
                       stdin=subprocess.DEVNULL, capture_output=True,
                       text=True, errors="replace")
    t = time.perf_counter() - t0
    log = feb.with_suffix(".log").read_text(errors="replace")
    if r.returncode != 0 or "N O R M A L   T E R M I N A T I O N" not in log:
        raise RuntimeError(f"FEBio no termino bien en {feb} "
                           f"(codigo {r.returncode}); ver {feb.with_suffix('.log')}")
    m = re.search(r"version\s+(\d+\.\d+\.\d+)", log)
    return t, (m.group(1) if m else "?")


def _eps_eff(sig, E, nu):
    """Deformacion efectiva de Pistoia a partir de la tension (Voigt, Pa).

    spinpy la calcula de sigma y eps; aqui eps sale de sigma por la
    flexibilidad isotropa, con distorsiones de ingenieria. Es la misma
    energia: U = sigma . S sigma / 2.
    """
    S = np.zeros((6, 6))
    S[:3, :3] = -nu / E
    np.fill_diagonal(S[:3, :3], 1.0 / E)
    S[3, 3] = S[4, 4] = S[5, 5] = 2.0 * (1.0 + nu) / E
    eps = sig @ S.T
    U = 0.5 * np.einsum("ij,ij->i", sig, eps)
    return np.sqrt(np.maximum(2.0 * U / E, 0.0))


def _von_mises(s):
    return np.sqrt(np.maximum(
        0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2
               + (s[:, 2] - s[:, 0]) ** 2)
        + 3.0 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2), 0.0))


def _rel(a, b):
    return float(abs(a - b) / abs(b)) if b else float("nan")


# ---------------------------------------------------------------------------
# Un caso
# ---------------------------------------------------------------------------

def comparar(clave, dir_corridas):
    BW, spc, apoyo, familia = CASOS[clave]()
    nx, ny, nz = BW.shape
    fila = {"caso": clave, "familia": familia, "n": [nx, ny, nz],
            "spacing_mm": spc.tolist(), "apoyo": apoyo,
            "BVTV": float(BW.mean())}
    print(f"\n== {clave}: {BW.shape}, BV/TV {BW.mean():.4f}, {apoyo}")

    t0 = time.perf_counter()
    res = ensayo_compresion(BW, spc, E_s=E_S, nu_s=NU_S, sigma0=SIGMA_SPINPY,
                            apoyo=apoyo)
    fila["t_spinpy_s"] = time.perf_counter() - t0
    if not res["ok"]:
        raise RuntimeError(f"spinpy no resolvio {clave}: {res['msg']}")
    fila.update(solver_spinpy=res.get("solver"),
                residuo_spinpy=res.get("residuo_rel"),
                n_elem=res["n_elem"], n_dof=res["n_dof"],
                frac_portante=res["frac_portante"])

    # La malla que resolvio spinpy: los mismos voxeles portantes, recorridos
    # en el mismo orden. Si el numero de nodos no casa, el campo u no se
    # puede comparar nodo a nodo y no se sigue.
    nodos, elems, _ = malla_hex(_solo_portante(BW), spc)
    if 3 * nodos.shape[0] != res["n_dof"] or elems.shape[0] != res["n_elem"]:
        raise RuntimeError("La malla exportada no es la que resolvio spinpy")

    d = dir_corridas / clave
    d.mkdir(parents=True, exist_ok=True)
    A_bruta = float(nx * spc[0] * ny * spc[1])

    def febio(sigma):
        """Corre FEBio a `sigma` (Pa) y devuelve u (mm) y sigma (Pa)
        reescalados a la carga de spinpy, o None si no converge."""
        feb = d / f"{clave}_{sigma:.0f}Pa.feb"
        info = escribir_febio(nodos, elems, feb, E_s=E_S, nu_s=NU_S,
                              sigma_app=sigma, A_bruta=A_bruta, apoyo=apoyo)
        try:
            t, ver = correr_febio(feb)
        except RuntimeError as e:
            print(f"   FEBio a {sigma:g} Pa: {e}")
            return None, info
        esc = SIGMA_SPINPY / sigma
        u = _ultimo_registro(feb.with_name(feb.stem + "_u.txt"), 3).ravel()
        s = _ultimo_registro(feb.with_name(feb.stem + "_s.txt"), 6)
        fila["version_febio"] = ver
        fila.setdefault("t_febio_s", {})[f"{sigma:g}"] = t
        return (u * esc, s * esc * 1e6), info

    r1, info = febio(SIGMA_LIN[0])
    r2, _ = febio(SIGMA_LIN[1])
    if r1 is None or r2 is None:
        raise RuntimeError(f"FEBio no convergio a carga pequena en {clave}")
    (u1, s1), (u2, s2) = r1, r2
    # Extrapolacion a carga nula: el termino no lineal de primer orden se
    # cancela. La diferencia entre las dos cargas mide ese termino.
    u_f, s_f = 2 * u1 - u2, 2 * s1 - s2
    fila["desvio_nl_1kPa"] = float(np.abs(u1 - u_f).max() / np.abs(u_f).max())

    # --- equilibrio: la fuerza que aplico FEBio es la de spinpy -----------
    # FEBio 4.5 escribe Rz = 0 en los apoyos `zero displacement`, asi que no
    # se suman reacciones. En el problema discreto, la integral de volumen de
    # la tension es exacta: sum_e sigma_zz V_e = sum_nodos z f_z = -F H (la
    # base esta en z = 0 y no contribuye). La media de Gauss por V_e ES la
    # integral del elemento.
    F_spinpy_N = res["F_total"] * 1e-6               # Pa*mm^2 -> N
    H = float(nodos[:, 2].max() - nodos[:, 2].min())
    fila["F_spinpy_N"] = F_spinpy_N
    fila["F_febio_N"] = float(-s_f[:, 2].sum() * np.prod(spc) / H) * 1e-6
    fila["dF_rel"] = _rel(fila["F_febio_N"], F_spinpy_N)

    # --- E_app, con la misma definicion que spinpy ------------------------
    z = nodos[:, 2]
    techo = z >= z.max() - 1e-9 * max(z.max() - z.min(), 1.0)
    eps_app = abs(float(np.mean(u_f[2::3][techo]))) / float(z.max() - z.min())
    E_f = SIGMA_SPINPY / eps_app
    u_s_ref = res["u"]
    fila.update(E_app_spinpy_MPa=res["E_app"] / 1e6, E_app_febio_MPa=E_f / 1e6,
                dE_rel=_rel(E_f, res["E_app"]))

    # --- respuesta NO lineal a la carga de referencia de la app -----------
    nl, _ = febio(SIGMA_SPINPY)
    if nl is None:
        fila["no_converge_1MPa"] = True
    else:
        u_nl = nl[0]
        e_nl = abs(float(np.mean(u_nl[2::3][techo]))) / float(z.max() - z.min())
        fila["E_app_febio_nl_1MPa"] = SIGMA_SPINPY / e_nl / 1e6
        fila["dE_nl_1MPa"] = _rel(SIGMA_SPINPY / e_nl, res["E_app"])
        fila["du_nl_1MPa"] = float(np.abs(u_nl - u_s_ref).max()
                                   / np.abs(u_s_ref).max())

    # --- campos -----------------------------------------------------------
    u_s = res["u"]
    fila["du_max_rel"] = float(np.abs(u_f - u_s).max() / np.abs(u_s).max())
    sig_s = res["sig"]
    vm_s = res["vm"]
    fila["dsig_max_rel"] = float(np.abs(s_f - sig_s).max() / vm_s.max())
    vm_f = _von_mises(s_f)
    fila["dvm_max_elem_rel"] = float(np.abs(vm_f - vm_s).max() / vm_s.max())

    # --- estadisticos citables, con las MISMAS funciones -------------------
    res_f = {"ok": True, "sigma_app": res["sigma_app"],
             "F_total": res["F_total"], "E_app": E_f,
             "eps_eff_solido": _eps_eff(s_f, E_S, NU_S), "vm_solido": vm_f,
             "superficie_solido": res["superficie_solido"]}
    p_s, p_f = criterio_pistoia(res), criterio_pistoia(res_f)
    for k in ("vm_p99_superficie", "factor", "sigma_fallo", "vm_max",
              "vm_p99"):
        fila[f"{k}_spinpy"] = p_s.get(k)
        fila[f"{k}_febio"] = p_f.get(k)
        fila[f"d_{k}_rel"] = _rel(p_f[k], p_s[k])
    fila["vm_n_superficie"] = p_s.get("vm_n_superficie")
    fila["F_fallo_spinpy_N"] = p_s["F_fallo"] * 1e-6
    fila["F_fallo_febio_N"] = p_f["F_fallo"] * 1e-6

    if familia == "spinodoide":
        a = AJUSTADO
        fila["procedencia"] = procedencia.bloque(
            "spinodoide", esquema="rechazo", semilla=a["seed"],
            wave_number_pi=a["wave_pi"], num_waves=a["num_waves"],
            thetas=list(a["thetas"]), rho=a["rho"], Lc_mm=a["Lc"])
    else:
        fila["procedencia"] = procedencia.bloque(familia)
    fila["procedencia"]["febio"] = fila["version_febio"]
    fila["feb"] = {k: v for k, v in info.items() if k != "ruta"}

    print(f"   E_app  spinpy {fila['E_app_spinpy_MPa']:.6f}  FEBio "
          f"{fila['E_app_febio_MPa']:.6f} MPa   dE {fila['dE_rel']:.2e}")
    print(f"   F      dF {fila['dF_rel']:.2e}   du {fila['du_max_rel']:.2e}   "
          f"dsig {fila['dsig_max_rel']:.2e}")
    print(f"   p99sup d {fila['d_vm_p99_superficie_rel']:.2e}   Pistoia d "
          f"{fila['d_factor_rel']:.2e}   vm_max d {fila['d_vm_max_rel']:.2e}")
    print(f"   no lineal: desvio a 1 kPa {fila['desvio_nl_1kPa']:.2e};  a 1 MPa "
          f"dE {fila.get('dE_nl_1MPa', float('nan')):.2e}  du "
          f"{fila.get('du_nl_1MPa', float('nan')):.2e}")
    print(f"   tiempos spinpy {fila['t_spinpy_s']:.1f} s  FEBio "
          f"{sum(fila['t_febio_s'].values()):.1f} s   ({res.get('solver')}, residuo "
          f"{res.get('residuo_rel', float('nan')):.1e})")
    return fila


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--casos", nargs="*", default=list(CASOS))
    ap.add_argument("--salida", default=str(AQUI / "resultados" / "febio.jsonl"))
    ap.add_argument("--corridas", default=str(AQUI / "corridas"))
    a = ap.parse_args()
    if not FEBIO.exists():
        sys.exit(f"No se encuentra FEBio en {FEBIO}")
    salida = Path(a.salida)
    salida.parent.mkdir(parents=True, exist_ok=True)
    print(f"spinpy {VERSION_SPINPY}  ->  {salida}")
    for clave in a.casos:
        fila = comparar(clave, Path(a.corridas))
        with open(salida, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(fila, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
