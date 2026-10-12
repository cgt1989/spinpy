"""
informe_tiempos.py: Tablas y figuras del informe de TIEMPOS de los motores FEM.

    python comparativa_motores/informe_tiempos.py [carpeta_jsonl] [carpeta_ref]

Lee los `.jsonl` que escribe `correr.py --salida <carpeta>` (por omision
`tiempos_sesion/`) y, si existe, la campana publicada (`resultados/`) para
contrastar. Escribe en `tiempos_sesion/`: `tablas_tiempos.md` y las figuras
`t_*.png`. Nada se calcula aqui salvo agregados (minimos, razones,
exponentes de escala por regresion log-log): solo se ordena lo medido.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
ENT = Path(sys.argv[1]) if len(sys.argv) > 1 else AQUI / "tiempos_sesion"
REF = Path(sys.argv[2]) if len(sys.argv) > 2 else AQUI / "resultados"

MOTORES = ("app", "ngsolve", "fenicsx", "skfem", "sfepy")
NOMBRE = {"app": "App (spinpy)", "ngsolve": "NGSolve", "fenicsx": "FEniCSx",
          "skfem": "scikit-fem", "sfepy": "SfePy"}
COLOR = {"app": "#52514e", "ngsolve": "#2a78d6", "fenicsx": "#eb6834",
         "skfem": "#1baf7a", "sfepy": "#eda100"}
MARCA = {"app": "o", "ngsolve": "s", "fenicsx": "^", "skfem": "D",
         "sfepy": "v"}
TINTA, REJILLA, FONDO = "#0b0b0b", "#e4e3df", "#fcfcfb"


def leer(carpeta, grupo):
    f = Path(carpeta) / f"{grupo}.jsonl"
    if not f.exists():
        return []
    return [json.loads(x) for x in f.read_text(encoding="utf-8").splitlines()
            if x.strip()]


def clase(fila):
    """'propio' (app), 'iterativo' (CG) o 'directo' (factorizacion)."""
    if fila["motor"] == "app":
        return "propio"
    return "iterativo" if "CG" in fila["solver"] else "directo"


def t_total(f):
    return f.get("tiempo_total_s") or sum(f.get("tiempos", {}).values())


def tabla(cabecera, filas):
    s = "| " + " | ".join(cabecera) + " |\n"
    s += "|" + "|".join("---" if i == 0 else "---:"
                        for i in range(len(cabecera))) + "|\n"
    for r in filas:
        s += "| " + " | ".join(r) + " |\n"
    return s


def fmt_t(x):
    if x is None:
        return "n/d"
    return f"{x:.2f}" if x < 10 else (f"{x:.1f}" if x < 1000 else f"{x:.0f}")


def casos_de(filas):
    orden = []
    for f in filas:
        if f["caso"] not in orden:
            orden.append(f["caso"])
    return orden


def matriz(filas, cas, columnas):
    """Filas por caso; una columna (motor, clase) con el tiempo total."""
    out = []
    for c in cas:
        d = {}
        for f in filas:
            k = (f["motor"], clase(f))
            if f["caso"] == c and (k not in d or (f.get("ok") and (
                    not d[k].get("ok") or t_total(f) < t_total(d[k])))):
                d[k] = f
        g = next(f for f in filas if f["caso"] == c)
        r = [c, f"{g['meta_caso']['n_gdl']:,}".replace(",", " ")]
        for k in columnas:
            f = d.get(k)
            r.append("n/d" if f is None else
                     (fmt_t(t_total(f)) if f.get("ok") else "FALLO"))
        out.append(r)
    return out


COLS = [("app", "propio"), ("ngsolve", "directo"), ("ngsolve", "iterativo"),
        ("fenicsx", "directo"), ("fenicsx", "iterativo"),
        ("skfem", "directo"), ("skfem", "iterativo"),
        ("sfepy", "directo"), ("sfepy", "iterativo")]
CAB = (["Caso", "GDL"] + [f"{NOMBRE[m].split(' ')[0]} {c[:3]}."
                          for m, c in COLS])


def exponente(x, y):
    x, y = np.log(np.asarray(x, float)), np.log(np.asarray(y, float))
    if len(x) < 2:
        return None
    return float(np.polyfit(x, y, 1)[0])


def mejor_por_motor(filas, caso):
    out = {}
    for f in filas:
        if f["caso"] == caso and f.get("ok"):
            t = t_total(f)
            if f["motor"] not in out or t < out[f["motor"]][0]:
                out[f["motor"]] = (t, clase(f), f.get("rss_pico_MB"))
    return out


def figura_escala(filas, caso_prefijo, titulo, nombre):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ejes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True,
                             facecolor=FONDO)
    for ax, cl in zip(ejes, ("directo", "iterativo")):
        for m in MOTORES:
            pts = sorted((f["meta_caso"]["n_gdl"], t_total(f))
                         for f in filas if f["motor"] == m and f.get("ok")
                         and f["caso"].startswith(caso_prefijo)
                         and (clase(f) == cl or (m == "app"
                                                 and cl == "iterativo")))
            if pts:
                x, y = zip(*pts)
                ax.plot(x, y, marker=MARCA[m], color=COLOR[m], lw=1.6,
                        ms=6, label=NOMBRE[m])
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title("Resolvedor " + cl + (" (con la app)" if cl ==
                                           "iterativo" else ""),
                     color=TINTA, fontsize=10)
        ax.set_xlabel("Grados de libertad", color=TINTA)
        ax.grid(True, color=REJILLA, lw=0.7, which="both")
        ax.set_facecolor(FONDO)
    ejes[0].set_ylabel("Tiempo total del motor (s)", color=TINTA)
    ejes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(titulo, color=TINTA, fontsize=11)
    fig.tight_layout()
    fig.savefig(ENT / nombre, dpi=160, facecolor=FONDO)
    plt.close(fig)


def figura_etapas(filas, caso, nombre):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sel = []
    for m in MOTORES:
        c = mejor_por_motor(filas, caso).get(m)
        if c:
            f = next(x for x in filas if x["caso"] == caso
                     and x["motor"] == m and x.get("ok")
                     and clase(x) == c[1])
            sel.append((f"{NOMBRE[m]}\n({c[1]})", f["tiempos"]))
    if not sel:
        return
    etapas = []
    for _, t in sel:
        for k in t:
            if k not in etapas:
                etapas.append(k)
    cols = ["#6b7b8c", "#2a78d6", "#eda100", "#1baf7a", "#eb6834",
            "#9b59b6"]
    fig, ax = plt.subplots(figsize=(8.5, 4), facecolor=FONDO)
    base = np.zeros(len(sel))
    for i, e in enumerate(etapas):
        v = np.array([t.get(e, 0.0) for _, t in sel])
        ax.barh(range(len(sel)), v, left=base, color=cols[i % len(cols)],
                label=e, height=0.6)
        base += v
    ax.set_yticks(range(len(sel)))
    ax.set_yticklabels([s for s, _ in sel], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Tiempo (s)", color=TINTA)
    ax.set_title(f"Etapas de la configuracion mas rapida de cada motor: "
                 f"{caso}", color=TINTA, fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.set_facecolor(FONDO)
    ax.grid(True, axis="x", color=REJILLA, lw=0.7)
    fig.tight_layout()
    fig.savefig(ENT / nombre, dpi=160, facecolor=FONDO)
    plt.close(fig)


def main():
    grupos = ["exactos", "cavidad", "espinodoide_hex", "espinodoide_tet",
              "nl"]
    D = {g: leer(ENT, g) for g in grupos}
    P = {g: leer(REF, g) for g in grupos} if REF.exists() else {}
    md = []
    md.append("## Tablas de tiempos\n")
    md.append("Tiempo total del motor en segundos (malla propia del motor, "
              "montaje, resolucion y extraccion; sin importar la "
              "biblioteca). `dir.` = factorizacion directa, `ite.` = CG con "
              "multigrid algebraico; la app usa su propio CG + AMG.\n")
    for g in ("exactos", "cavidad", "espinodoide_hex", "espinodoide_tet"):
        if not D[g]:
            continue
        md.append(f"\n### {g}\n")
        md.append(tabla(CAB, matriz(D[g], casos_de(D[g]), COLS)))

    if D["nl"]:
        md.append("\n### nl (Newton, resolvedor directo)\n")
        filas = D["nl"]
        cas = casos_de(filas)
        cols = [("ngsolve", "directo"), ("fenicsx", "directo"),
                ("skfem", "directo"), ("sfepy", "directo")]
        md.append(tabla(["Caso", "GDL"] + [NOMBRE[m] for m, _ in cols],
                        matriz(filas, cas, cols)))

    # mas rapido por motor y caso
    md.append("\n## Configuracion mas rapida de cada motor\n")
    for g in ("espinodoide_hex", "espinodoide_tet"):
        if not D[g]:
            continue
        md.append(f"\n### {g}\n")
        filas_t = []
        for c in casos_de(D[g]):
            mj = mejor_por_motor(D[g], c)
            gdl = next(f for f in D[g] if f["caso"] == c)["meta_caso"][
                "n_gdl"]
            r = [c, f"{gdl:,}".replace(",", " ")]
            for m in MOTORES:
                r.append("n/d" if m not in mj else
                         f"{fmt_t(mj[m][0])} ({mj[m][1][:3]}.; "
                         f"{mj[m][2] / 1024:.1f} GB)")
            filas_t.append(r)
        md.append(tabla(["Caso", "GDL"] + [NOMBRE[m] for m in MOTORES],
                        filas_t))

    # exponentes de escala
    md.append("\n## Exponente de escala\n")
    md.append("Pendiente de la regresion log-log del tiempo total frente a "
              "los GDL (hexaedros; todos los tamanos del espinodoide con "
              "resultado). Un exponente de 1 es escala lineal.\n\n")
    filas_e = []
    for m in MOTORES:
        for cl in ("propio", "directo", "iterativo"):
            pts = sorted((f["meta_caso"]["n_gdl"], t_total(f))
                         for f in D["espinodoide_hex"]
                         if f["motor"] == m and clase(f) == cl
                         and f.get("ok"))
            if len(pts) >= 3:
                x, y = zip(*pts)
                filas_e.append([NOMBRE[m], cl, str(len(pts)),
                                f"{exponente(x, y):.2f}"])
    md.append(tabla(["Motor", "Resolvedor", "Puntos", "Exponente"], filas_e))

    # contraste con la campana publicada
    if P.get("espinodoide_hex") and D["espinodoide_hex"]:
        md.append("\n## Contraste con la campana publicada\n")
        md.append("Cociente entre el tiempo medido ahora y el publicado "
                  "(`resultados/`) para la misma combinacion de caso, motor "
                  "y resolvedor (hexaedros del espinodoide). Mide la "
                  "reproducibilidad entre campanas, no el error del motor.\n")
        filas_c = []
        for g in ("espinodoide_hex", "espinodoide_tet"):
            for c in casos_de(D[g]):
                r = [c]
                for m in MOTORES:
                    for cl in ("propio", "directo", "iterativo"):
                        a = [f for f in D[g] if f["caso"] == c
                             and f["motor"] == m and clase(f) == cl
                             and f.get("ok")]
                        b = [f for f in P.get(g, []) if f["caso"] == c
                             and f["motor"] == m and clase(f) == cl
                             and f.get("ok")]
                        if a and b:
                            r.append(f"{t_total(a[0]) / t_total(b[0]):.2f}")
                        elif (m, cl) in COLS:
                            r.append("n/d")
                filas_c.append(r)
        cab = ["Caso"] + [f"{NOMBRE[m].split(' ')[0]} {c[:3]}."
                          for m, c in COLS]
        md.append(tabla(cab, filas_c))

    # fallos
    fall = [f for g in grupos for f in D[g] if not f.get("ok")]
    md.append("\n## Corridas fallidas o omitidas\n")
    if fall:
        md.append(tabla(["Caso", "Motor", "Resolvedor", "Motivo"],
                        [[f["caso"], f["motor"], f["solver"],
                          str(f.get("motivo"))] for f in fall]))
    else:
        md.append("Ninguna corrida fallo.\n")

    (ENT / "tablas_tiempos.md").write_text("\n".join(md), encoding="utf-8")
    if D["espinodoide_hex"]:
        figura_escala(D["espinodoide_hex"], "espinodoide_hex",
                      "Espinodoide con hexaedros: tiempo frente a tamano",
                      "t_escala_hex.png")
        figura_etapas(D["espinodoide_hex"], "espinodoide_hex_n64",
                      "t_etapas_hex64.png")
    if D["espinodoide_tet"]:
        figura_escala(D["espinodoide_tet"], "espinodoide_tet",
                      "Espinodoide con TET10: tiempo frente a tamano",
                      "t_escala_tet.png")
        figura_etapas(D["espinodoide_tet"], "espinodoide_tet_n32",
                      "t_etapas_tet32.png")
    print("listo:", ENT / "tablas_tiempos.md")


if __name__ == "__main__":
    main()
