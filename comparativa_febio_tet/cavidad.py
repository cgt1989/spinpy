"""cavidad.py — ¿Converge la malla suave (TET10) al valor exacto en una cavidad?

    cd Port_Python
    python comparativa_febio_tet/cavidad.py                  # n = 16 24 32
    python comparativa_febio_tet/cavidad.py --n 16 24 32 48

EL PROBLEMA
-----------
`Estudio_Convergencia` mostro, en una celda PERIODICA con una cavidad
esferica, que el pico de von Mises con voxeles no converge: se pasa un +5.6 %
del valor de Goodier y salta ~10 puntos al refinar. La pregunta que motiva la
malla suave es la objecion de regularidad: ¿quitar los escalones hace
converger el pico?

Con tetraedros no hay nodos emparejados en caras opuestas, asi que aqui el
problema es FINITO: un cubo de lado L = 1 mm con una cavidad centrada de radio
r = 0.2 L (la r/L del estudio C1), sometido al ensayo de la app (1 MPa sobre
la seccion bruta, apoyo deslizante, E 20 GPa, nu 0.30). La capa superficial
es exactamente la pared de la cavidad: las caras del cubo se excluyen por la
regla de las isocaps.

LA REFERENCIA
-------------
Goodier (medio infinito) NO es la solucion de este problema finito: la
seccion neta a la altura de la cavidad es un 12.6 % menor que la bruta. La
referencia "exacta" es el MISMO problema resuelto en FEBio con TET10 sobre la
ESFERA ANALITICA (no desde voxeles), a dos refinamientos. Goodier se reporta
solo como control de orden de magnitud.

PREDICCIONES Y TOLERANCIAS, ESCRITAS ANTES DE LA PRIMERA CORRIDA (2026-09-24)
-----------------------------------------------------------------------------
P1  La referencia analitica esta convergida: entre sus dos refinamientos el
    maximo de von Mises en la pared de la cavidad cambia <= 1 %, y E_app
    <= 0.1 %. Si no, la referencia no sirve y el resto no se juzga.
P2  (la pregunta) Con la malla suave desde voxeles, el MAXIMO de von Mises en
    la pared se acerca a la referencia al refinar la mascara, y a la mayor n
    corrida queda a <= 3 % de ella. Con ladrillos (hex8) el maximo no
    converge (es lo ya medido); se reporta para contrastar.
P3  El p99 de la capa superficial (hazen, ponderado por volumen) de la malla
    suave queda a <= 3 % del de la referencia a la mayor n.
P4  E_app de la malla suave queda a <= 1 % de la referencia a la mayor n
    (la corrige la perdida de volumen, `corregir_volumen=True`).
Si P2 falla, la malla suave NO responde a la objecion de regularidad en este
caso y hay que decirlo asi.

Resultados: una linea por corrida en `resultados/cavidad.jsonl`.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))
from spinpy import __version__ as VERSION_SPINPY           # noqa: E402
from spinpy import febio, procedencia                       # noqa: E402
from spinpy.solido import PERM_TETGEN_A_C3D10, verificar_tet10  # noqa: E402

R_REL = 0.20
L = 1.0
NU = 0.30


def goodier(nu=NU):
    """von Mises en el ecuador de una cavidad esferica, medio infinito."""
    kt = (27.0 - 15.0 * nu) / (2.0 * (7.0 - 5.0 * nu))
    sff = (15.0 * nu - 3.0) / (2.0 * (7.0 - 5.0 * nu))
    return float(np.sqrt(0.5 * ((kt - sff) ** 2 + sff ** 2 + kt ** 2)))


def geometria(n):
    """Mascara (solido = fuera de la esfera) y la esfera en coordenadas de malla.

    Convencion de `superficie_cerrada`: el voxel i tiene su centro en i*h y
    el cubo va de -h/2 a (n - 1/2) h.
    """
    h = L / n
    c = np.arange(n) * h
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    centro = np.full(3, (n / 2.0 - 0.5) * h)
    r = R_REL * L
    d = np.sqrt((X - centro[0]) ** 2 + (Y - centro[1]) ** 2
                + (Z - centro[2]) ** 2)
    return d > r, np.full(3, h), centro, r


def malla_analitica(n, res_esfera, tam):
    """TET10 sobre el cubo con la esfera ANALITICA como hueco."""
    import pyvista as pv
    import tetgen
    _, spc, centro, r = geometria(n)
    a, b = -0.5 * spc[0], (n - 0.5) * spc[0]
    cubo = pv.Cube(bounds=(a, b, a, b, a, b)).triangulate().subdivide(3)
    esfera = pv.Sphere(radius=r, center=centro, theta_resolution=res_esfera,
                       phi_resolution=res_esfera)
    tg = tetgen.TetGen(pv.merge([esfera, cubo]).clean())
    tg.add_hole(list(centro))
    nod, el = tg.tetrahedralize(order=2, minratio=1.5,
                                maxvolume=febio.maxvolume_de_tamano(tam))[:2]
    nod = np.asarray(nod, float)
    el = np.asarray(el, np.int64)[:, PERM_TETGEN_A_C3D10]
    ok, v = verificar_tet10(nod, el)
    if not ok:
        raise RuntimeError(v["msg"])
    return febio.malla_de_tet10(nod, el, (n, n, n), spc,
                                informe={"origen": "esfera analitica",
                                         "res_esfera": res_esfera,
                                         "tam_max_mm": tam})


def medir(malla, carpeta, etiqueta):
    prot = febio.protocolo("app")
    t0 = time.perf_counter()
    reg = febio.ensayo(malla, prot, carpeta / etiqueta, analisis=("lineal",),
                       prefijo=etiqueta)
    lin = reg["lineal"]
    if lin is None:
        return {"etiqueta": etiqueta, "ok": False, "fallos": reg["fallos"]}
    s = reg["_campos_lineal"]["sigma"]
    vm = febio.von_mises(s) / prot["sigma_app"]
    sup = malla["superficie"]
    p = lin["pistoia"]
    return {"etiqueta": etiqueta, "ok": True, "malla": malla["tipo"],
            "n_elems": int(malla["elems"].shape[0]),
            "n_nodos": int(malla["nodos"].shape[0]),
            "n_superficie": int(sup.sum()),
            "BVTV_malla": float(malla["vol_elem"].sum() / L ** 3),
            "E_app_MPa": lin["E_app"] / 1e6, "dF_rel": lin["dF_rel"],
            "vm_max_pared": float(vm[sup].max()),
            "vm_p99_pared": float(p["vm_p99_superficie"]) / prot["sigma_app"],
            "tiempo_s": time.perf_counter() - t0,
            "memoria_MB": max((c["memoria_MB"] or 0)
                              for c in reg["corridas"]),
            "febio": reg["febio"]["version"],
            "informe_malla": {k: v for k, v in malla["informe"].items()
                              if isinstance(v, (int, float, str, bool))}}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--n", nargs="*", type=int, default=[16, 24, 32])
    ap.add_argument("--referencia", nargs="*", type=int, default=[64, 96],
                    help="resoluciones de la esfera analitica")
    ap.add_argument("--suavizado", nargs="*", type=int, default=[None],
                    help="iteraciones de Taubin de la malla suave (barrido)")
    ap.add_argument("--decimado", type=float, default=None)
    ap.add_argument("--solo-tet10", action="store_true")
    ap.add_argument("--salida", default=str(AQUI / "resultados"
                                            / "cavidad.jsonl"))
    ap.add_argument("--corridas", default=str(AQUI / "corridas" / "cavidad"))
    a = ap.parse_args()
    salida = Path(a.salida)
    salida.parent.mkdir(parents=True, exist_ok=True)
    carpeta = Path(a.corridas)
    print(f"spinpy {VERSION_SPINPY}; Goodier (medio infinito) vm/sigma0 = "
          f"{goodier():.4f}")

    def anotar(fila):
        fila["procedencia"] = procedencia.bloque("cavidad_esferica")
        fila["r_rel"] = R_REL
        with open(salida, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(fila, ensure_ascii=False) + "\n")
        if fila.get("ok"):
            print(f"  {fila['etiqueta']:22s} elems {fila['n_elems']:7d}  "
                  f"E {fila['E_app_MPa']:9.2f} MPa  vm_max {fila['vm_max_pared']:.4f}"
                  f"  p99 {fila['vm_p99_pared']:.4f}  BV/TV {fila['BVTV_malla']:.4f}"
                  f"  {fila['tiempo_s']:.0f} s  {fila['memoria_MB']:.0f} MB",
                  flush=True)
        else:
            print(f"  {fila['etiqueta']}: FALLO {fila.get('fallos')}", flush=True)

    # Referencia: n fija (32) para el cubo; se refina la esfera y el tamano.
    for k, res in enumerate(a.referencia):
        m = malla_analitica(32, res, tam=0.06 / (1 + k))
        anotar({**medir(m, carpeta, f"ref_esfera{res}"), "n": 32,
                "tipo_caso": "referencia"})
    for n in a.n:
        BW, spc, _, _ = geometria(n)
        if not a.solo_tet10:
            m = febio.mallar(BW, spc, "hex8")
            anotar({**medir(m, carpeta, f"hex8_n{n}"), "n": n,
                    "r_h": R_REL * n, "tipo_caso": "voxeles"})
        for it in a.suavizado:
            op = {}
            if it is not None:
                op["suavizado"] = it
            if a.decimado is not None:
                op["decimado"] = a.decimado
            et = f"tet10_n{n}" + (f"_s{it}" if it is not None else "") + \
                (f"_d{a.decimado:g}" if a.decimado is not None else "")
            m = febio.mallar(BW, spc, "tet10", **op)
            anotar({**medir(m, carpeta, et), "n": n, "r_h": R_REL * n,
                    "tipo_caso": "voxeles", "opciones_malla": op})


if __name__ == "__main__":
    main()
