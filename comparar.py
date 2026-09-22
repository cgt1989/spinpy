"""
comparar.py — Contraste Python vs MATLAB contra la banda de ruido estocastico.

QUE SE COMPARA Y POR QUE ASI
----------------------------
El generador de spinodoides es ESTOCASTICO: `spinodoid` remuestrea las
direcciones de onda en cada llamada, de modo que dos evaluaciones del mismo
vector de parametros dan realizaciones distintas y metricas distintas. Existe
por tanto una varianza irreducible por debajo de la cual dos numeros NO son
distinguibles. Esa varianza esta medida en
`Validacion_Anexo/resultados/ruido_estocastico.csv`.

En consecuencia, la pregunta correcta NO es "?dan Python y MATLAB el mismo
numero?" —no pueden, ni siquiera MATLAB consigo mismo— sino:

    ?la diferencia entre las medias de las dos implementaciones cae dentro
     de la dispersion que la propia implementacion de MATLAB exhibe consigo
     misma?

El criterio operativo es la diferencia de medias expresada en unidades de la
desviacion tipica de MATLAB (una z). |z| < 2 significa indistinguible al nivel
de ruido del propio generador.

Se contrastan tres cosas:
  1. Python 'rechazo' vs MATLAB voxel-primero  -> fidelidad del port.
  2. Python 'equitativo' vs MATLAB             -> magnitud del desacuerdo por
     la convencion de muestreo del repositorio del profesor.
  3. MATLAB voxel-primero vs MATLAB ruta-malla -> cuanto cuesta el rodeo por
     la malla, que es lo que el port se ahorra.

Uso:  python comparar.py         (requiere resultados/referencia_matlab.json)
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from spinpy import generar_mascara, morfometria   # noqa: E402

RAIZ = Path(__file__).parent
SALIDA = RAIZ / "resultados"
REF = SALIDA / "referencia_matlab.json"
RUIDO = RAIZ.parent / "Validacion_Anexo" / "resultados" / "ruido_estocastico.csv"

CAMPOS = ["BVTV", "DA", "DA2", "BSBV", "BS", "TbTh", "TbSp", "TbN",
          "PoTot", "FracPort"]


def cargar_banda_ruido():
    """CV% por metrica medido en MATLAB (run_ruido_estocastico.m)."""
    if not RUIDO.exists():
        return {}
    banda = {}
    with open(RUIDO, encoding="utf-8") as fh:
        for fila in csv.DictReader(fh):
            banda.setdefault(fila["metrica"], []).append(float(fila["cv_pct"]))
    return {k: float(np.mean(v)) for k, v in banda.items()}


def _aplanar_campos(c):
    """jsonencode de MATLAB envuelve el cell una vez de mas: [[...]] -> [...]."""
    while isinstance(c, list) and len(c) == 1 and isinstance(c[0], list):
        c = c[0]
    return list(c)


def _vec(v):
    """Vector float aplanado; los null de MATLAB (NaN) pasan a np.nan."""
    v = _aplanar_campos(v) if isinstance(v, list) else [v]
    return np.array([np.nan if x is None else float(x) for x in np.ravel(v)])


def estadisticas(valores):
    """media/std ignorando NaN, columna a columna."""
    A = np.asarray(valores, dtype=float)
    if A.ndim == 1:
        A = A[None, :]
    med = np.full(A.shape[1], np.nan)
    sd = np.full(A.shape[1], np.nan)
    for j in range(A.shape[1]):
        v = A[:, j][np.isfinite(A[:, j])]
        if v.size:
            med[j] = v.mean()
        if v.size > 1:
            sd[j] = v.std(ddof=1)
    return med, sd


CACHE = SALIDA / "python_cache.json"


def _clave(par, esquema, n_rep, resolucion, semilla0, ip):
    return (f"{par['nombre']}|{esquema}|{n_rep}|{resolucion}|{semilla0}|{ip}"
            f"|{par['dens']}|{par['wave']}|{int(par['nw'])}"
            f"|{[int(t) for t in np.ravel(par['th'])]}")


def _cache_leer():
    if CACHE.exists():
        try:
            return json.loads(CACHE.read_text(encoding="utf-8"))
        except Exception:
            pass
    return {}


def _cache_escribir(d):
    CACHE.write_text(json.dumps(d, indent=1), encoding="utf-8")


def corre_python(par, esquema, n_rep, resolucion, spacing, semilla0, ip,
                 usar_cache=True):
    """n_rep realizaciones con el port, mismas semillas nominales que MATLAB.

    El resultado se cachea en resultados/python_cache.json: el lado Python se
    puede precalcular mientras MATLAB genera su referencia (--precalcular),
    ya que los parametros no dependen de la salida de MATLAB.
    """
    cache = _cache_leer() if usar_cache else {}
    k = _clave(par, esquema, n_rep, resolucion, semilla0, ip)
    if k in cache:
        return np.asarray(cache[k]["M"], dtype=float), float(cache[k]["t"])

    M = np.full((n_rep, len(CAMPOS)), np.nan)
    t0 = time.time()
    for r in range(n_rep):
        BW, _, _ = generar_mascara(
            resolution=resolucion,
            wave_number=par["wave"] * np.pi,
            num_waves=int(par["nw"]),
            thetas=par["th"],
            rho=par["dens"],
            esquema=esquema,
            seed=semilla0 + 1000 * ip + r + 1,
        )
        m = morfometria(BW, spacing)
        for q, c in enumerate(CAMPOS):
            v = m.get(c, np.nan)
            if np.isscalar(v) or (isinstance(v, float)):
                M[r, q] = float(v)
    t = (time.time() - t0) / n_rep

    cache = _cache_leer()
    cache[k] = {"M": M.tolist(), "t": t}
    _cache_escribir(cache)
    return M, t


def precalcular(n_rep):
    """Calcula y cachea el lado Python sin esperar a MATLAB.

    Los parametros de los puntos no dependen de la salida de MATLAB, asi que
    esta fase puede correr EN PARALELO con gen_referencia.m. Solo necesita la
    lista de puntos, que se lee de cualquier referencia_matlab.json previa.
    """
    if not REF.exists():
        print(f"Falta {REF}: hace falta una corrida previa (aunque sea SPIN_NREP=1)"
              " para conocer la lista de puntos.")
        return 2
    ref = json.loads(REF.read_text(encoding="utf-8"))
    resolucion = int(ref["resolucion"])
    spacing = [float(ref["spacing_mm"])] * 3
    semilla0 = int(ref["semilla0"])
    resultados = ref["resultados"]
    if isinstance(resultados, dict):
        resultados = [resultados]

    print(f"Precalculando lado Python: {len(resultados)} puntos x {n_rep} "
          f"realizaciones x 2 esquemas, resolucion={resolucion}")
    for ip, blk in enumerate(resultados, start=1):
        par = blk["parametros"]
        for esquema in ("rechazo", "equitativo"):
            _, t = corre_python(par, esquema, n_rep, resolucion, spacing,
                                semilla0, ip)
            print(f"  {blk['punto']:16s} {esquema:11s} {t:5.1f} s/realizacion")
    print(f"Cache -> {CACHE}")
    return 0


def main():
    if not REF.exists():
        print(f"Falta {REF}.\nEjecuta primero:  matlab -batch \"cd('{RAIZ}'); gen_referencia\"")
        return 2

    ref = json.loads(REF.read_text(encoding="utf-8"))
    resolucion = int(ref["resolucion"])
    spacing = [float(ref["spacing_mm"])] * 3
    semilla0 = int(ref["semilla0"])
    banda = cargar_banda_ruido()

    resultados = ref["resultados"]
    if isinstance(resultados, dict):
        resultados = [resultados]

    print("=" * 92)
    print(" CONTRASTE PORT PYTHON  vs  REFERENCIA MATLAB")
    print(f" resolucion={resolucion}  spacing={spacing[0]:.6g} mm  "
          f"nrep_matlab={ref['nrep']}")
    print(" z = (media_py - media_mat) / std_mat   ->   |z| < 2 = indistinguible")
    print("=" * 92)

    salida = {"meta": {"resolucion": resolucion, "spacing_mm": spacing[0],
                       "nrep_matlab": ref["nrep"], "banda_ruido_cv_pct": banda},
              "puntos": []}
    z_rechazo, z_equit = [], []

    for ip, blk in enumerate(resultados, start=1):
        par = blk["parametros"]
        nombre = blk["punto"]
        campos_mat = _aplanar_campos(blk["voxel"]["campos"])
        med_mat = _vec(blk["voxel"]["media"])
        std_mat = _vec(blk["voxel"]["std"])
        med_mal = _vec(blk["malla"]["media"])

        # Reordenar por si MATLAB emitiera otro orden de campos. Las metricas
        # que la copia de applib no calcula (DA2, FracPort: correcciones K1/K3
        # posteriores a la extraccion) se marcan NaN y quedan fuera del
        # contraste; ver nota sobre applib desincronizado en el README.
        def _tomar(v, c):
            return v[campos_mat.index(c)] if c in campos_mat else np.nan
        med_mat = np.array([_tomar(med_mat, c) for c in CAMPOS])
        std_mat = np.array([_tomar(std_mat, c) for c in CAMPOS])
        med_mal = np.array([_tomar(med_mal, c) for c in CAMPOS])

        n_rep = int(ref["nrep"])
        print(f"\n--- {nombre} : dens={par['dens']} wave={par['wave']}pi "
              f"nw={int(par['nw'])} thetas={[int(t) for t in np.ravel(par['th'])]} ---")

        M_rec, t_rec = corre_python(par, "rechazo", n_rep, resolucion,
                                    spacing, semilla0, ip)
        M_equ, t_equ = corre_python(par, "equitativo", n_rep, resolucion,
                                    spacing, semilla0, ip)
        med_rec, std_rec = estadisticas(M_rec)
        med_equ, _ = estadisticas(M_equ)

        print(f"    {'metrica':9s} {'MATLAB':>11s} {'+-std':>9s} "
              f"{'PY rechazo':>11s} {'z':>7s} | {'PY equitat':>11s} {'z':>7s} "
              f"| {'MAT malla':>10s}")
        for q, c in enumerate(CAMPOS):
            s = std_mat[q]
            zr = (med_rec[q] - med_mat[q]) / s if np.isfinite(s) and s > 0 else np.nan
            ze = (med_equ[q] - med_mat[q]) / s if np.isfinite(s) and s > 0 else np.nan
            if np.isfinite(zr):
                z_rechazo.append(abs(zr))
            if np.isfinite(ze):
                z_equit.append(abs(ze))
            marca = " " if (not np.isfinite(zr) or abs(zr) < 2) else "*"
            print(f"  {marca} {c:9s} {med_mat[q]:>11.5g} {s:>9.3g} "
                  f"{med_rec[q]:>11.5g} {zr:>7.2f} | {med_equ[q]:>11.5g} "
                  f"{ze:>7.2f} | {med_mal[q]:>10.5g}")

        print(f"    tiempo/realizacion:  python={t_rec:.1f} s   "
              f"matlab={float(np.ravel(blk['voxel']['tiempo_s'])[0]):.1f} s")

        salida["puntos"].append({
            "punto": nombre, "parametros": par, "campos": CAMPOS,
            "matlab_voxel_media": med_mat.tolist(),
            "matlab_voxel_std": std_mat.tolist(),
            "matlab_malla_media": med_mal.tolist(),
            "python_rechazo_media": med_rec.tolist(),
            "python_rechazo_std": std_rec.tolist(),
            "python_equitativo_media": med_equ.tolist(),
            "z_rechazo": [float((med_rec[q] - med_mat[q]) / std_mat[q])
                          if np.isfinite(std_mat[q]) and std_mat[q] > 0 else None
                          for q in range(len(CAMPOS))],
            "z_equitativo": [float((med_equ[q] - med_mat[q]) / std_mat[q])
                             if np.isfinite(std_mat[q]) and std_mat[q] > 0 else None
                             for q in range(len(CAMPOS))],
            "tiempo_py_s": t_rec,
        })

    print("\n" + "=" * 92)
    print(" RESUMEN")
    print("=" * 92)
    zr = np.asarray(z_rechazo)
    ze = np.asarray(z_equit)
    print(f"  Python 'rechazo'    vs MATLAB : |z| mediano={np.median(zr):.2f}  "
          f"max={zr.max():.2f}   fuera de banda (|z|>=2): "
          f"{int((zr >= 2).sum())}/{zr.size}")
    print(f"  Python 'equitativo' vs MATLAB : |z| mediano={np.median(ze):.2f}  "
          f"max={ze.max():.2f}   fuera de banda (|z|>=2): "
          f"{int((ze >= 2).sum())}/{ze.size}")
    salida["resumen"] = {
        "rechazo": {"z_mediano": float(np.median(zr)), "z_max": float(zr.max()),
                    "fuera_banda": int((zr >= 2).sum()), "n": int(zr.size)},
        "equitativo": {"z_mediano": float(np.median(ze)), "z_max": float(ze.max()),
                       "fuera_banda": int((ze >= 2).sum()), "n": int(ze.size)},
    }

    (SALIDA / "comparacion.json").write_text(
        json.dumps(salida, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n  Resultados -> {SALIDA / 'comparacion.json'}")
    return 0


if __name__ == "__main__":
    if "--precalcular" in sys.argv:
        n = 8
        for a in sys.argv[1:]:
            if a.isdigit():
                n = int(a)
        sys.exit(precalcular(n))
    sys.exit(main())
