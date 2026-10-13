"""
validacion_exp.py: contraste de lo predicho con lo medido en un ensayo fisico.

POR QUE EXISTE
--------------
La comparacion con FEBio y entre motores prueba que spinpy resuelve bien el
problema que plantea (malla de voxeles, tejido homogeneo e isotropo, condicion
de carga ideal, criterio de Pistoia). No prueba que ese problema represente el
hueso. Eso solo lo contrasta un ensayo fisico del mismo espécimen, y este modulo
es la herramienta para hacer ese contraste y para reportarlo sin adornos.

ESTADO ACTUAL (V2.2.0)
----------------------
El repositorio NO incluye ensayos mecanicos del hueso al que corresponden los
VOIs de `datos/`. Mientras no existan, los valores mecanicos absolutos de
spinpy (E_app, carga de fallo) son resultados de simulacion verificados
numericamente, no predicciones validadas experimentalmente. Este modulo deja
listo el procedimiento: cuando haya datos medidos se entregan en un CSV
(plantilla en `datos/porcino_koria2020/ensayos_experimentales_PLANTILLA.csv`)
y se obtiene el contraste con `spinpy-validar-exp`.

QUE CALCULA
-----------
Dadas parejas (predicho, medido) de una misma magnitud:

- Error: RMSE, MAE, error relativo medio absoluto y sesgo relativo.
- Concordancia de Lin (CCC), que penaliza tanto la dispersion como el sesgo
  sistematico; la correlacion de Pearson sola no detecta que todo este
  desplazado por un factor.
- Regresion medido ~ predicho con intervalos de confianza de la pendiente y el
  intercepto (t de Student, n-2 g.l.). Pendiente 1 e intercepto 0 es el ideal.
- Bland-Altman: sesgo y limites de acuerdo al 95 %, con el IC de cada limite
  (Bland y Altman 1986, 1999). Se hace en escala logaritmica (razon medido /
  predicho), porque el error de la mecanica trabecular es multiplicativo.
- Spearman, como medida de ordenacion que no depende de la escala.
- Calibracion del modulo del tejido E_s: la rigidez aparente es lineal en E_s,
  de modo que, conocido el E_s con que se calculo E_app, se estima el E_s que
  mejor ajusta lo medido (minimos cuadrados por el origen), con su error en
  validacion cruzada dejando uno fuera. Distingue un error de escala (E_s mal
  supuesto, que se corrige) de un error de forma (la geometria no explica las
  diferencias entre especimenes, que no se corrige).

LO QUE NO HACE
--------------
- No decide por usted que el acuerdo sea "suficiente". Solo emite veredicto si
  se le pasan criterios declarados ANTES de mirar los datos (`criterios=`); sin
  ellos informa "no evaluado". Fijar el umbral despues de ver el resultado es
  el error que este modulo existe para evitar.
- No corrige por medidas repetidas en un mismo animal. Con varios VOI por
  animal los pares no son independientes y los intervalos salen demasiado
  estrechos: use una sola pareja por animal o un modelo mixto.
- Con n < 8 los intervalos son muy anchos; el modulo lo avisa.

Referencias metodologicas
-------------------------
Lin LI (1989) A concordance correlation coefficient to evaluate
reproducibility. Biometrics 45:255-268.
Bland JM, Altman DG (1986) Statistical methods for assessing agreement between
two methods of clinical measurement. Lancet 327:307-310.
Bland JM, Altman DG (1999) Measuring agreement in method comparison studies.
Stat Methods Med Res 8:135-160.
Erdemir A, Guess TM, Halloran J, Tadepalli SC, Morrison TM (2012) Considerations
for reporting finite element analysis studies in biomechanics. J Biomech
45:625-633.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

N_MIN = 3          # menos que esto no admite regresion con IC
N_AVISO = 8        # por debajo, los IC son demasiado anchos para concluir

CLAVES_CRITERIOS = ("ccc_min", "sesgo_rel_max", "mape_max", "pendiente_ic_incluye_1")


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------

def _pares(pred, med):
    p = np.asarray(pred, float).ravel()
    m = np.asarray(med, float).ravel()
    if p.shape != m.shape:
        raise ValueError(f"predicho y medido deben tener el mismo largo "
                         f"({p.size} frente a {m.size})")
    if not (np.all(np.isfinite(p)) and np.all(np.isfinite(m))):
        raise ValueError("hay valores no finitos en predicho o medido")
    if p.size < N_MIN:
        raise ValueError(f"se necesitan al menos {N_MIN} parejas, hay {p.size}")
    return p, m


def ccc_lin(p, m):
    """Coeficiente de concordancia de Lin, con varianzas poblacionales (1/n)."""
    mp, mm = p.mean(), m.mean()
    vp, vm = p.var(), m.var()
    cov = float(np.mean((p - mp) * (m - mm)))
    den = vp + vm + (mp - mm) ** 2
    return float(2.0 * cov / den) if den > 0 else float("nan")


def _regresion(p, m, alfa):
    n = p.size
    r = stats.linregress(p, m)
    t = stats.t.ppf(1 - alfa / 2, n - 2)
    # intercept_stderr existe desde scipy 1.6
    s_int = float(getattr(r, "intercept_stderr", np.nan))
    return {
        "pendiente": float(r.slope),
        "pendiente_ic": [float(r.slope - t * r.stderr),
                         float(r.slope + t * r.stderr)],
        "intercepto": float(r.intercept),
        "intercepto_ic": [float(r.intercept - t * s_int),
                          float(r.intercept + t * s_int)],
        "r": float(r.rvalue),
        "r2": float(r.rvalue ** 2),
        "pearson_p": float(r.pvalue),
    }


def _bland_altman_log(p, m, alfa):
    """Bland-Altman sobre log(medido/predicho); devuelve factores multiplicativos."""
    if np.any(p <= 0) or np.any(m <= 0):
        return None
    d = np.log(m / p)
    n = d.size
    sesgo, sd = float(d.mean()), float(d.std(ddof=1))
    t = stats.t.ppf(1 - alfa / 2, n - 1)
    z = stats.norm.ppf(1 - alfa / 2)
    # Bland y Altman 1999: ee de los limites = sd * sqrt(1/n + z^2 / (2(n-1)))
    ee_lim = sd * math.sqrt(1.0 / n + z ** 2 / (2.0 * (n - 1)))
    ee_sesgo = sd / math.sqrt(n)
    inf, sup = sesgo - z * sd, sesgo + z * sd
    return {
        "factor_sesgo": math.exp(sesgo),
        "factor_sesgo_ic": [math.exp(sesgo - t * ee_sesgo),
                            math.exp(sesgo + t * ee_sesgo)],
        "limite_inf": math.exp(inf),
        "limite_sup": math.exp(sup),
        "limite_inf_ic": [math.exp(inf - t * ee_lim), math.exp(inf + t * ee_lim)],
        "limite_sup_ic": [math.exp(sup - t * ee_lim), math.exp(sup + t * ee_lim)],
        "sd_log": sd,
    }


def calibrar_Es(pred, med, Es_usado=None):
    """Factor de escala que mejor ajusta lo medido, y el E_s que implica.

    La rigidez aparente de un modelo lineal es proporcional al modulo del
    tejido E_s. Si `pred` se calculo con un E_s conocido (`Es_usado`), el
    factor k de `medido ~ k * pred` (minimos cuadrados por el origen) es el
    cociente entre el E_s que mejor ajusta y el usado: Es_calibrado = k *
    Es_usado.

    Se devuelve tambien el error relativo en validacion cruzada dejando uno
    fuera (cada espécimen predicho con el k de los demas) y el R2 por el
    origen. Un error de escala (E_s mal supuesto) se corrige y deja un error
    de validacion cruzada pequeño; un error de forma (la geometria no explica
    las diferencias entre especimenes) no se corrige y el error de validacion
    cruzada sigue alto.

    Solo vale para magnitudes lineales en E_s (E_app); NO para la carga de
    fallo de Pistoia, que depende de la deformacion critica y no de E_s.
    """
    x, y = _pares(pred, med)
    k = float(np.sum(x * y) / np.sum(x * x))
    n = x.size
    err = np.empty(n)
    for i in range(n):
        s = np.arange(n) != i
        ki = float(np.sum(x[s] * y[s]) / np.sum(x[s] * x[s]))
        err[i] = (ki * x[i] - y[i]) / y[i]
    ss_res = float(np.sum((y - k * x) ** 2))
    ss_tot = float(np.sum(y ** 2))
    out = {
        "k": k,
        "r2_origen": 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
        "error_rel_loo_medio_abs": float(np.mean(np.abs(err))),
        "error_rel_loo": [float(e) for e in err],
    }
    if Es_usado:
        out["Es_usado"] = float(Es_usado)
        out["Es_calibrado"] = k * float(Es_usado)
    return out


# ---------------------------------------------------------------------------
# Contraste
# ---------------------------------------------------------------------------

def acuerdo(pred, med, ids=None, alfa=0.05, criterios=None, magnitud=""):
    """Contraste completo de una magnitud. Devuelve un diccionario serializable.

    criterios (opcional, declarados antes de medir), claves admitidas:
      ccc_min                 CCC de Lin minimo
      sesgo_rel_max           |sesgo relativo medio| maximo
      mape_max                error relativo medio absoluto maximo
      pendiente_ic_incluye_1  True: el IC de la pendiente debe contener a 1
    Sin criterios, "veredicto" es "no evaluado".
    """
    p, m = _pares(pred, med)
    n = p.size
    rel = (p - m) / m if np.all(m != 0) else np.full(n, np.nan)
    rho = stats.spearmanr(p, m)
    res = {
        "magnitud": magnitud,
        "n": int(n),
        "alfa": alfa,
        "ids": list(ids) if ids is not None else None,
        "rmse": float(np.sqrt(np.mean((p - m) ** 2))),
        "mae": float(np.mean(np.abs(p - m))),
        "mape": float(np.nanmean(np.abs(rel))),
        "sesgo_rel": float(np.nanmean(rel)),
        "ccc": ccc_lin(p, m),
        "spearman": float(rho.statistic),
        "spearman_p": float(rho.pvalue),
        "regresion": _regresion(p, m, alfa),
        "bland_altman_log": _bland_altman_log(p, m, alfa),
        "avisos": [],
    }
    if n < N_AVISO:
        res["avisos"].append(
            f"n = {n}: con menos de {N_AVISO} parejas los intervalos son muy "
            f"anchos y el contraste no permite concluir equivalencia.")
    if res["bland_altman_log"] is None:
        res["avisos"].append("Hay valores no positivos: no se calcula el "
                             "Bland-Altman logaritmico.")
    res["criterios"] = dict(criterios) if criterios else None
    res["veredicto"], res["detalle_veredicto"] = _veredicto(res, criterios)
    return res


def _veredicto(res, criterios):
    if not criterios:
        return "no evaluado", ["No se declararon criterios de aceptacion "
                               "antes de medir."]
    desconocidas = set(criterios) - set(CLAVES_CRITERIOS)
    if desconocidas:
        raise ValueError(f"criterios desconocidos: {sorted(desconocidas)}; "
                         f"admitidos: {CLAVES_CRITERIOS}")
    det, ok = [], True

    def _eval(nombre, cumple, texto):
        nonlocal ok
        ok = ok and cumple
        det.append(f"{'cumple' if cumple else 'NO cumple'}: {texto}")

    if "ccc_min" in criterios:
        _eval("ccc", res["ccc"] >= criterios["ccc_min"],
              f"CCC = {res['ccc']:.3f} (minimo {criterios['ccc_min']})")
    if "sesgo_rel_max" in criterios:
        _eval("sesgo", abs(res["sesgo_rel"]) <= criterios["sesgo_rel_max"],
              f"|sesgo relativo| = {abs(res['sesgo_rel']):.3f} "
              f"(maximo {criterios['sesgo_rel_max']})")
    if "mape_max" in criterios:
        _eval("mape", res["mape"] <= criterios["mape_max"],
              f"error relativo medio absoluto = {res['mape']:.3f} "
              f"(maximo {criterios['mape_max']})")
    if criterios.get("pendiente_ic_incluye_1"):
        lo, hi = res["regresion"]["pendiente_ic"]
        _eval("pendiente", lo <= 1.0 <= hi,
              f"IC de la pendiente [{lo:.3f}, {hi:.3f}] contiene a 1")
    return ("cumple" if ok else "no cumple"), det


# ---------------------------------------------------------------------------
# Entrada / salida
# ---------------------------------------------------------------------------

COLUMNAS = ("id", "magnitud", "predicho", "medido")


def leer_csv(ruta):
    """Lee el CSV de contraste y agrupa por magnitud.

    Columnas obligatorias: id, magnitud, predicho, medido. Las filas con
    `medido` vacio se ignoran (es lo que trae la plantilla sin rellenar).
    Devuelve {magnitud: (ids, predicho, medido)}.
    """
    grupos = {}
    with open(ruta, newline="", encoding="utf-8") as f:
        lector = csv.DictReader(f)
        falta = [c for c in COLUMNAS if c not in (lector.fieldnames or [])]
        if falta:
            raise ValueError(f"faltan columnas en el CSV: {falta}")
        for fila in lector:
            if not (fila["medido"] or "").strip():
                continue
            g = grupos.setdefault(fila["magnitud"].strip(), ([], [], []))
            g[0].append(fila["id"].strip())
            g[1].append(float(fila["predicho"]))
            g[2].append(float(fila["medido"]))
    return grupos


def informe_md(resultados):
    """Resumen en Markdown de uno o varios contrastes."""
    L = ["# Contraste con ensayos fisicos", ""]
    for r in resultados:
        L.append(f"## {r['magnitud'] or 'magnitud'} (n = {r['n']})")
        L.append("")
        reg, ba = r["regresion"], r["bland_altman_log"]
        L.append(f"- RMSE {r['rmse']:.4g}; MAE {r['mae']:.4g}; error relativo "
                 f"medio absoluto {100 * r['mape']:.1f} %; sesgo relativo "
                 f"{100 * r['sesgo_rel']:+.1f} %.")
        L.append(f"- CCC de Lin {r['ccc']:.3f}; Spearman {r['spearman']:.3f}; "
                 f"R2 {reg['r2']:.3f}.")
        L.append(f"- medido = {reg['pendiente']:.3f} x predicho "
                 f"{reg['intercepto']:+.4g} "
                 f"(IC {100 * (1 - r['alfa']):.0f} %, pendiente "
                 f"[{reg['pendiente_ic'][0]:.3f}, {reg['pendiente_ic'][1]:.3f}]).")
        if ba:
            L.append(f"- Razon medido/predicho: factor medio "
                     f"{ba['factor_sesgo']:.3f}; limites de acuerdo "
                     f"[{ba['limite_inf']:.3f}, {ba['limite_sup']:.3f}].")
        L.append(f"- Veredicto: **{r['veredicto']}**.")
        for d in r["detalle_veredicto"]:
            L.append(f"  - {d}")
        for a in r["avisos"]:
            L.append(f"- Aviso: {a}")
        L.append("")
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="spinpy-validar-exp",
        description="Contrasta lo predicho por spinpy con ensayos fisicos.")
    ap.add_argument("csv", help="CSV con id, magnitud, predicho, medido")
    ap.add_argument("--json", help="guarda los resultados en este JSON")
    ap.add_argument("--ccc-min", type=float)
    ap.add_argument("--sesgo-rel-max", type=float)
    ap.add_argument("--mape-max", type=float)
    ap.add_argument("--pendiente-incluye-1", action="store_true")
    a = ap.parse_args(argv)

    crit = {}
    if a.ccc_min is not None:
        crit["ccc_min"] = a.ccc_min
    if a.sesgo_rel_max is not None:
        crit["sesgo_rel_max"] = a.sesgo_rel_max
    if a.mape_max is not None:
        crit["mape_max"] = a.mape_max
    if a.pendiente_incluye_1:
        crit["pendiente_ic_incluye_1"] = True

    grupos = leer_csv(a.csv)
    if not grupos:
        print("El CSV no tiene ninguna fila con 'medido' rellenado: no hay "
              "nada que contrastar todavia.", file=sys.stderr)
        return 2
    res = [acuerdo(p, m, ids=i, criterios=crit or None, magnitud=mag)
           for mag, (i, p, m) in grupos.items()]
    print(informe_md(res))
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=2, ensure_ascii=False),
                                encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
