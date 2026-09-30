"""
correr.py — Campana de la comparativa de motores FEM.

    cd comparativa_motores
    python correr.py exactos cavidad espinodoide_hex espinodoide_tet nl

Cada corrida (caso x motor x resolvedor) va en un PROCESO HIJO: la memoria de
pico es la del motor y un fallo (memoria agotada, violacion de segmento) no
arrastra la campana. El lanzador vigila la memoria residente del hijo y lo
detiene por encima de `--rss-max` GB, anotando el fallo.

FEniCSx corre con el Python de su entorno de conda-forge (`--python-fenicsx`
o la variable FENICSX_PYTHON); los demas, con el del lanzador.

Salida: una linea JSON por corrida en `resultados/<grupo>.jsonl` con tiempos
por etapa, memoria, resolvedor, residuo, iteraciones, las magnitudes del
postproceso comun y la diferencia de desplazamientos frente a la
REFERENCIA del caso (ver `_referencia`).
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
import casos                                                    # noqa: E402
import comun                                                    # noqa: E402

MOTORES = ("app", "skfem", "sfepy", "ngsolve", "fenicsx")
#: Orden de preferencia de la referencia: una solucion por resolvedor
#: DIRECTO (residuo ~1e-14), de cualquiera de los motores; si ninguno corrio
#: en directo, la iterativa de menor residuo.
PREF_REF = ("ngsolve", "fenicsx", "sfepy", "skfem")
GDL_MAX_SUPERLU = 1.5e5


def _rss_hijo_MB(pid):
    try:
        with open(f"/proc/{pid}/status") as fh:
            for linea in fh:
                if linea.startswith("VmRSS:"):
                    return int(linea.split()[1]) / 1024.0
    except OSError:
        pass
    return 0.0


def correr_uno(caso_npz, motor, solver, py, tmp, tiempo_max, rss_max_MB,
               entorno):
    salida = Path(tmp) / f"{motor}_{solver}.npz"
    if salida.exists():
        salida.unlink()
    cmd = [py, "-m", f"motores.mot_{motor}", str(caso_npz), str(salida)]
    t0 = time.perf_counter()
    p = subprocess.Popen(cmd, cwd=AQUI, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True, env=entorno)
    pico, motivo = 0.0, None
    while p.poll() is None:
        time.sleep(0.2)
        pico = max(pico, _rss_hijo_MB(p.pid))
        if pico > rss_max_MB:
            p.kill()
            motivo = f"memoria > {rss_max_MB / 1024:.0f} GB"
        elif time.perf_counter() - t0 > tiempo_max:
            p.kill()
            motivo = f"tiempo > {tiempo_max:.0f} s"
    out, err = p.communicate()
    t = time.perf_counter() - t0
    if motivo or p.returncode != 0 or not salida.exists():
        ultima = [x for x in (err or "").strip().splitlines()
                  if x.strip()][-1:] or [""]
        return {"ok": False, "motivo": motivo or ultima[0][:300],
                "codigo": p.returncode, "tiempo_pared_s": t}, None, None
    d = np.load(salida, allow_pickle=False)
    meta = json.loads(str(d["meta"]))
    extra = {k: d[k] for k in d.files if k not in ("u", "meta")}
    meta.update(ok=True, tiempo_pared_s=t, rss_vigilado_MB=pico)
    return meta, d["u"], extra


def _referencia(corridas):
    """(motor, solver) de la referencia del caso."""
    for mot in PREF_REF:
        k = (mot, "directo")
        if k in corridas and corridas[k][0].get("ok"):
            return k
    ok = [(k, v[0].get("residuo_rel") or 1.0) for k, v in corridas.items()
          if v[0].get("ok") and k[0] != "app"]
    if ok:
        return min(ok, key=lambda x: x[1])[0]
    return ("app", "propio") if ("app", "propio") in corridas else None


def campana(nombre_grupo, lista, motores, solvers, py_fx, tiempo_max,
            rss_max_MB, hilos):
    res = AQUI / "resultados" / f"{nombre_grupo}.jsonl"
    res.parent.mkdir(exist_ok=True)
    entorno = dict(os.environ, OMP_NUM_THREADS=str(hilos),
                   MKL_NUM_THREADS=str(hilos),
                   OPENBLAS_NUM_THREADS=str(hilos), PYTHONWARNINGS="ignore")
    for fabrica in lista:
        t0 = time.perf_counter()
        caso = fabrica()
        t_caso = time.perf_counter() - t0
        m = caso["meta"]
        print(f"== {m['nombre']}: {m['n_elems']} elementos, {m['n_gdl']} GDL"
              f" ({t_caso:.1f} s)", flush=True)
        tmp = tempfile.mkdtemp(prefix="cm_", dir=os.environ.get("CM_TMP"))
        corridas = {}
        try:
            for motor in motores:
                if motor == "app" and not (m["tipo"] == "hex8"
                                           and m["analisis"] == "lineal"):
                    continue
                if motor == "sfepy" and m["analisis"] == "nl" and \
                        m["material"] != "svk":
                    continue
                for solver in (["propio"] if motor == "app" else solvers):
                    if m["analisis"] == "nl" and solver != "directo":
                        continue
                    if solver == "directo" and \
                            m["n_gdl"] > m.get("gdl_max_directo", 2e6):
                        continue
                    # SuperLU (scipy) no tiene paralelismo ni ordenacion de
                    # diseccion anidada: por encima de ~1.5e5 GDL en 3D su
                    # relleno lo deja fuera de la hora por corrida.
                    if solver == "directo" and motor in ("skfem", "sfepy") \
                            and m["n_gdl"] > GDL_MAX_SUPERLU:
                        continue
                    caso["meta"]["solver"] = solver
                    npz = comun.guardar(caso, Path(tmp) / "caso.npz")
                    py = py_fx if motor == "fenicsx" else sys.executable
                    meta, u, extra = correr_uno(npz, motor, solver, py, tmp,
                                                tiempo_max, rss_max_MB,
                                                entorno)
                    corridas[(motor, solver)] = (meta, u, extra)
                    tt = meta.get("tiempos", {})
                    print(f"   {motor:8s} {solver:9s} "
                          + (f"{sum(tt.values()):8.2f} s  "
                             f"{meta.get('rss_pico_MB', 0):7.0f} MB  "
                             f"res {meta.get('residuo_rel') or 0:.1e}"
                             if meta.get("ok") else
                             f"FALLO: {meta.get('motivo')}"), flush=True)
            ref = _referencia(corridas)
            u_ref = corridas[ref][1] if ref else None
            for (motor, solver), (meta, u, extra) in corridas.items():
                fila = {"grupo": nombre_grupo, "caso": m["nombre"],
                        "meta_caso": {k: v for k, v in m.items()
                                      if k != "solver"},
                        "motor": motor, "solver": solver,
                        "referencia": list(ref) if ref else None,
                        "hilos": hilos, "maquina": platform.processor()
                        or platform.machine(), **meta}
                if meta.get("ok"):
                    if m["analisis"] == "lineal":
                        pp = comun.postproceso(caso, u)
                        fila.update({k: v for k, v in pp.items()
                                     if not k.startswith("_")})
                        if u_ref is not None:
                            fila["du_rel_ref"] = comun.dif_rel(u, u_ref)
                            pr = comun.postproceso(caso, u_ref)
                            fila["dvm_rel_ref"] = comun.dif_rel(pp["_vm"],
                                                                pr["_vm"])
                        if ("app", "propio") in corridas and \
                                corridas[("app", "propio")][1] is not None:
                            fila["du_rel_app"] = comun.dif_rel(
                                u, corridas[("app", "propio")][1])
                    else:
                        fila["F_reac_N"] = extra["F_reac"].tolist()
                        if m["tipo"] == "hex8" and m.get("bloque"):
                            fila["F_exacta_N"] = [
                                -comun.uniaxial(m["material"], 1 - e)[0]
                                * m["A_bruta"] for e in m["eps_plato"]]
                        if u_ref is not None:
                            fila["du_rel_ref"] = comun.dif_rel(u, u_ref)
                with open(res, "a", encoding="utf-8") as fh:
                    fh.write(json.dumps(fila, ensure_ascii=False,
                                        default=float) + "\n")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


def grupos(a):
    n_hex = a.n_hex
    n_tet = a.n_tet
    eps_bl = [0.05, 0.10, 0.15, 0.20]
    eps_sp = [0.0025, 0.005, 0.01, 0.02]
    return {
        "exactos": [casos.bloque_hex, casos.bloque_tet],
        "cavidad": [lambda: casos.cavidad_hex(32),
                    lambda: casos.cavidad_tet(64, 0.06)],
        "espinodoide_hex": [lambda n=n: _dir_max(casos.espinodoide_hex(n),
                                                 a.gdl_directo_hex)
                            for n in n_hex],
        "espinodoide_tet": [lambda n=n: _dir_max(casos.espinodoide_tet(n),
                                                 a.gdl_directo_tet)
                            for n in n_tet],
        "nl": ([lambda mat=mat: casos.no_lineal(
                    _marca(casos.comun.caso_hex(np.ones((4, 4, 6), bool),
                                                0.05, "bloque_hex")),
                    mat, eps_bl) for mat in ("svk", "neohookeano")]
               + [lambda mat=mat: casos.no_lineal(
                    casos.espinodoide_hex(a.n_nl), mat, eps_sp)
                  for mat in ("svk", "neohookeano")]),
    }


def _marca(c):
    c["meta"]["bloque"] = True
    return c


def _dir_max(c, g):
    c["meta"]["gdl_max_directo"] = g
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("grupos", nargs="+")
    ap.add_argument("--motores", nargs="+", default=list(MOTORES))
    ap.add_argument("--solvers", nargs="+", default=["iterativo", "directo"])
    ap.add_argument("--n-hex", nargs="+", type=int,
                    default=[24, 32, 48, 64, 80])
    ap.add_argument("--n-tet", nargs="+", type=int, default=[32, 48])
    ap.add_argument("--n-nl", type=int, default=24)
    ap.add_argument("--gdl-directo-hex", type=float, default=4e5)
    ap.add_argument("--gdl-directo-tet", type=float, default=1.2e6)
    ap.add_argument("--hilos", type=int, default=os.cpu_count())
    ap.add_argument("--tiempo-max", type=float, default=3600.0)
    ap.add_argument("--rss-max", type=float, default=13.0, help="GB")
    ap.add_argument("--python-fenicsx", default=os.environ.get(
        "FENICSX_PYTHON", "/tmp/claude-0/mm/fx/bin/python"))
    a = ap.parse_args()
    g = grupos(a)
    for nombre in a.grupos:
        campana(nombre, g[nombre], a.motores, a.solvers, a.python_fenicsx,
                a.tiempo_max, a.rss_max * 1024, a.hilos)


if __name__ == "__main__":
    main()
