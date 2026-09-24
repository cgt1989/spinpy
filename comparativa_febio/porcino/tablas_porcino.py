"""tablas_porcino.py — Tablas Markdown del informe, leidas de `resultados/`.

    python comparativa_febio/porcino/tablas_porcino.py > tablas.md

Las cifras del INFORME y del README salen de aqui y no se teclean: asi un
numero del texto no puede discrepar del JSONL que lo produjo. Decimales con
coma (el informe es en castellano).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
from figuras_porcino import ORDEN, leer  # noqa: E402

NOM = {"voi": "VOI", "spinodoide": "Spinodoide", "dual-lattice": "Dual-lattice"}


def c(x, esp=".4g"):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return format(float(x), esp).replace(".", ",")


def e(x):
    """Diferencia relativa en notacion 1e-9 -> 1·10⁻⁹ compacta."""
    if x is None:
        return "—"
    x = float(x)
    if x == 0:
        return "0"
    m, p = f"{x:.0e}".split("e")
    sup = str(int(p)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{m}·10{sup}"


def pct(x, d=2):
    return "—" if x is None else (f"{100 * float(x):+.{d}f} %").replace(".", ",")


def ordenar(filas):
    return sorted(filas, key=lambda f: (ORDEN.index(f["estructura"]),
                                        f.get("n") or 0, f.get("paso") or 0))


def tabla(cab, filas):
    s = "| " + " | ".join(cab) + " |\n|" + "---|" * len(cab) + "\n"
    return s + "".join("| " + " | ".join(r) + " |\n" for r in filas)


def lineal(an):
    fs = ordenar(leer(an))
    return tabla(
        ["estructura", "BV/TV", "E_app spinpy (MPa)", "ΔE_app", "Δu máx.",
         "Δσ máx.", "Δ von Mises", "Δ p99 sup.", "Δ Pistoia", "ΔF",
         "Δ sesión"],
        [[NOM[f["estructura"]], c(f["BVTV"], ".3f"),
          c(f["E_app_spinpy_MPa"], ".1f"), e(f["dE_rel"]),
          e(f["du_max_rel"]), e(f["dsig_max_rel"]),
          e(f["dvm_max_elem_rel"]), e(f["d_vm_p99_superficie_rel"]),
          e(f["d_factor_rel"]), e(f["dF_rel"]),
          e(f.get("d_sesion_E_rel"))] for f in fs])


def magnitudes(an):
    fs = ordenar(leer(an))
    cab = ["estructura", "E_app (MPa)", "p99 sup. (MPa)", "máx. (MPa)",
           "n capa sup.", "desplaz. máx. (µm)"]
    if an == "compresion":
        cab += ["σ fallo (MPa)", "F fallo (N)"]
    filas = []
    for f in fs:
        r = [NOM[f["estructura"]], c(f["E_app_spinpy_MPa"], ".1f"),
             c(f["vm_p99_superficie_spinpy"] / 1e6, ".2f"),
             c(f["vm_max_spinpy"] / 1e6, ".1f"),
             str(f["vm_n_superficie"]),
             c(f["desp_max_spinpy_mm"] * 1e3, ".2f")]
        if an == "compresion":
            r += [c(f["sigma_fallo_spinpy"] / 1e6, ".2f"),
                  c(f["F_fallo_spinpy_N"], ".1f")]
        filas.append(r)
    return tabla(cab, filas)


def no_lineal():
    fs = ordenar(leer("compresion"))
    filas = []
    for f in fs:
        a, b = f.get("nl_ensayo") or {}, f.get("nl_fallo") or {}
        pn = f.get("plato_nl_fallo") or {}
        filas.append([
            NOM[f["estructura"]], pct(a.get("dE_rel"), 3),
            c(b.get("sigma_Pa", 0) / 1e6, ".2f"), pct(b.get("dE_rel")),
            pct(b.get("d_p99_sup_rel")),
            c(100 * b["frac_sobre_eps_crit"], ".2f") + " %"
            if b.get("frac_sobre_eps_crit") is not None else "—",
            c(f.get("E_app_plato_MPa"), ".1f"),
            c(f.get("plato_sobre_fuerza"), ".3f"), pct(pn.get("dE_rel"))])
    return tabla(["estructura", "ΔE no lineal a 1 MPa", "σ fallo (MPa)",
                  "ΔE no lineal a σ fallo", "Δ p99 sup. a σ fallo",
                  "tejido > 0,7 % a σ fallo", "E_app plato (MPa)",
                  "plato / fuerza", "ΔE plato no lineal a ε fallo"], filas)


def no_lineal_tapia():
    fs = ordenar(leer("comparado"))
    return tabla(["estructura", "σ aparente (MPa)", "ΔE no lineal a 100 N",
                  "Δ p99 sup.", "Δu máx."],
                 [[NOM[f["estructura"]], c(f["sigma_app_Pa"] / 1e6, ".3f"),
                   pct((f.get("nl_ensayo") or {}).get("dE_rel"), 3),
                   pct((f.get("nl_ensayo") or {}).get("d_p99_sup_rel"), 3),
                   pct((f.get("nl_ensayo") or {}).get("du_max_rel"), 3)]
                  for f in fs])


def homog():
    fs = ordenar(leer("homogeneizacion"))
    filas = []
    for f in fs:
        s, fb = f["constantes_spinpy"], f["constantes_febio"]
        filas.append([NOM[f["estructura"]], c(f["BVTV"], ".3f")]
                     + [c(s[k] / 1e6, ".1f") for k in ("Ex", "Ey", "Ez")]
                     + [c(s["Gyz"] / 1e6, ".1f"), e(f["dC_max_rel"]),
                        e(f["d_log_euclidea"]), e(f["asimetria_febio"]),
                        e(f.get("d_sesion_C_rel"))])
    return tabla(["estructura", "BV/TV", "E_x (MPa)", "E_y (MPa)",
                  "E_z (MPa)", "G_yz (MPa)", "ΔC máx. (rel.)",
                  "d log-euclídea", "asimetría FEBio", "Δ sesión"], filas)


def conv():
    fs = ordenar(leer("convergencia"))
    return tabla(["estructura", "n", "h (µm)", "BV/TV", "E_app spinpy (MPa)",
                  "ΔE FEBio", "E_app sesión (MPa)"],
                 [[NOM[f["estructura"]], str(f["n"]),
                   c(f["h_mm"] * 1e3, ".1f"), c(f["BVTV"], ".4f"),
                   c(f["E_app_spinpy_MPa"], ".1f"), e(f["dE_rel"]),
                   c(f.get("E_app_sesion_MPa"), ".1f")] for f in fs])


def fallo():
    fs = ordenar(leer("fallo"))
    return tabla(["estructura", "paso", "daño acum.", "E/E₀ spinpy",
                  "ΔE FEBio", "F fallo spinpy (N)", "ΔF FEBio",
                  "rotos spinpy / FEBio", "intactos distintos"],
                 [[NOM[f["estructura"]], str(f["paso"]),
                   c(100 * f["dano_acumulado"], ".1f") + " %",
                   c(f["E_rel_spinpy"], ".3f"), e(f["dE_rel"]),
                   c(f["F_fallo_spinpy_N"], ".1f"), e(f["dF_rel"]),
                   f"{f['rotos_spinpy']} / {f['rotos_febio']}",
                   str(f["voxeles_intactos_distintos"])] for f in fs])


def tiempos():
    out = []
    for an in ("compresion", "comparado", "homogeneizacion", "convergencia",
               "fallo"):
        fs = leer(an)
        if fs:
            t = sum(sum((f.get("t_febio_s") or {}).values()) for f in fs)
            n = sum(len(f.get("t_febio_s") or {}) for f in fs)
            out.append([an, str(n), c(t / 60, ".1f")])
    return tabla(["análisis", "corridas de FEBio", "minutos de FEBio"], out)


if __name__ == "__main__":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    for tit, fn in (("Ensayo de la app: lineal", lambda: lineal("compresion")),
                    ("Ensayo de la app: magnitudes",
                     lambda: magnitudes("compresion")),
                    ("Ensayo de la app: no lineal y plato", no_lineal),
                    ("Tapia: lineal", lambda: lineal("comparado")),
                    ("Tapia: magnitudes", lambda: magnitudes("comparado")),
                    ("Tapia: no lineal", no_lineal_tapia),
                    ("Homogeneizacion", homog), ("Convergencia", conv),
                    ("Fallo progresivo", fallo), ("Tiempos", tiempos)):
        print(f"### {tit}\n\n{fn()}")
