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
  desplazado por un factor. Se da con su IC por la transformacion z de
  Fisher (Lin 1989); los criterios comparan el LIMITE INFERIOR, no el
  estimador puntual.
- Regresion medido ~ predicho con intervalos de confianza de la pendiente y el
  intercepto (t de Student, n-2 g.l.). Pendiente 1 e intercepto 0 es el ideal.
- Bland-Altman: sesgo y limites de acuerdo al 95 %, con el IC de cada limite
  (Bland y Altman 1986, 1999). Se hace en escala logaritmica (razon medido /
  predicho), porque el error de la mecanica trabecular es multiplicativo. El
  IC de los limites usa la aproximacion de Bland y Altman (1999), que es
  anticonservadora con n pequeño (cobertura simulada ~0,92 con n = 5-8 y ~0,94
  con n = 20 para un nominal de 0,95); el modulo lo avisa si n < 20.
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
  el error que este modulo existe para evitar. Un criterio que no se puede
  calcular (NaN) da "no evaluable", nunca "no cumple" en silencio.
- Los mismos criterios se aplican a todas las magnitudes del CSV. Si E_app y la
  carga de fallo necesitan umbrales distintos, use un CSV por magnitud.
- Convenciones de signo: `sesgo_rel` = (predicho - medido) / medido (negativo
  si el modelo subestima); `factor_sesgo` = medido / predicho (mayor que 1 si
  el modelo subestima). Son las dos caras de lo mismo, con referencias opuestas.
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
N_BA_AVISO = 20    # por debajo, el IC de los limites de Bland-Altman es aproximado

CLAVES_CRITERIOS = ("ccc_min", "sesgo_rel_max", "mape_max", "pendiente_equiv_delta")


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


def ccc_ic(p, m, alfa=0.05):
    """IC del CCC de Lin por la transformacion z de Fisher (Lin 1989).

    Devuelve (inferior, superior) o None si no se puede calcular (r = 0,
    CCC = +-1 o varianza nula). Es una aproximacion asintotica: con n < 20 es
    orientativa.
    """
    n = p.size
    if n < 4 or p.std() == 0 or m.std() == 0:
        return None
    r = float(np.corrcoef(p, m)[0, 1])
    rc = ccc_lin(p, m)
    if not np.isfinite(rc) or abs(r) < 1e-12 or abs(rc) >= 1 or abs(rc) < 1e-12:
        return None
    u = (p.mean() - m.mean()) / math.sqrt(p.std() * m.std())
    d = (1 - rc ** 2)
    var_z = (1.0 / (n - 2)) * (
        (1 - r ** 2) * rc ** 2 * d / (r ** 2 * d ** 2)
        + 4 * rc ** 3 * (1 - rc) * u ** 2 / (r * d ** 2)
        - 2 * rc ** 4 * u ** 4 / (r ** 2 * d ** 2))
    if not np.isfinite(var_z) or var_z <= 0:
        return None
    z = math.atanh(rc)
    h = stats.norm.ppf(1 - alfa / 2) * math.sqrt(var_z)
    return (math.tanh(z - h), math.tanh(z + h))


def _regresion(p, m, alfa):
    n = p.size
    if np.ptp(p) == 0:
        # Predicho constante: la pendiente no esta definida. Segun la version
        # de SciPy, linregress devuelve NaN o lanza ValueError; aqui se trata
        # igual en todas.
        nan = float("nan")
        return {"pendiente": nan, "pendiente_ic": [nan, nan],
                "intercepto": nan, "intercepto_ic": [nan, nan],
                "r": nan, "r2": nan, "pearson_p": nan}
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
    if np.any(x <= 0) or np.any(y <= 0):
        raise ValueError("la calibracion de E_s exige predicho y medido "
                         "estrictamente positivos")
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
        # Factor geometrico: el error de la mecanica trabecular es
        # multiplicativo y los minimos cuadrados absolutos los dominan los
        # especimenes mas rigidos. Si k y k_geometrico difieren mucho, mandan
        # pocos especimenes.
        "k_geometrico": float(np.exp(np.mean(np.log(y / x)))),
        # R2 NO centrado (por el origen): sale alto aunque no haya relacion
        # entre especimenes. No lo use como evidencia de acuerdo; mire el
        # error de validacion cruzada.
        "r2_origen": 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
        "error_rel_loo_medio_abs": float(np.mean(np.abs(err))),
        "error_rel_loo": [float(e) for e in err],
    }
    if Es_usado:
        out["Es_usado"] = float(Es_usado)
        out["Es_calibrado"] = k * float(Es_usado)
        out["Es_calibrado_geometrico"] = out["k_geometrico"] * float(Es_usado)
    return out


# ---------------------------------------------------------------------------
# Contraste
# ---------------------------------------------------------------------------

def acuerdo(pred, med, ids=None, alfa=0.05, criterios=None, magnitud=""):
    """Contraste completo de una magnitud. Devuelve un diccionario serializable.

    criterios (opcional, declarados antes de medir), claves admitidas:
      ccc_min                 minimo del LIMITE INFERIOR del IC del CCC de Lin
      sesgo_rel_max           |sesgo relativo medio| maximo
      mape_max                error relativo medio absoluto maximo
      pendiente_equiv_delta   delta: el IC de la pendiente debe caer ENTERO
                              dentro de [1 - delta, 1 + delta] (equivalencia)
    Sin criterios, "veredicto" es "no evaluado"; si algun criterio no se puede
    calcular (NaN), es "no evaluable".
    """
    p, m = _pares(pred, med)
    n = p.size
    avisos = []
    if np.all(m > 0):
        rel = (p - m) / m
    else:
        rel = np.full(n, np.nan)
        avisos.append("Hay valores medidos <= 0: el error relativo (MAPE, "
                      "sesgo relativo) no esta definido y queda en NaN.")
    if np.ptp(p) == 0 or np.ptp(m) == 0:
        rho = None
        avisos.append("Una de las series es constante: la regresion y la "
                      "correlacion de Spearman no estan definidas.")
    else:
        rho = stats.spearmanr(p, m)
    ic_ccc = ccc_ic(p, m, alfa)
    res = {
        "magnitud": magnitud,
        "n": int(n),
        "alfa": alfa,
        "ids": list(ids) if ids is not None else None,
        "rmse": float(np.sqrt(np.mean((p - m) ** 2))),
        "mae": float(np.mean(np.abs(p - m))),
        "mape": float(np.mean(np.abs(rel))),
        "sesgo_rel": float(np.mean(rel)),
        "ccc": ccc_lin(p, m),
        "ccc_ic": list(ic_ccc) if ic_ccc else None,
        "spearman": float(rho.statistic) if rho else float("nan"),
        "spearman_p": float(rho.pvalue) if rho else float("nan"),
        "regresion": _regresion(p, m, alfa),
        "bland_altman_log": _bland_altman_log(p, m, alfa),
        "avisos": avisos,
    }
    if not np.isfinite(res["ccc"]):
        res["avisos"].append("CCC no definido (ambas series constantes): no "
                             "hay variacion que concordar.")
    elif ic_ccc is None:
        res["avisos"].append("No se pudo calcular el IC del CCC (acuerdo casi "
                             "perfecto, correlacion nula o n < 4).")
    if n < N_BA_AVISO and res["bland_altman_log"] is not None:
        res["avisos"].append(
            f"n = {n}: el IC de los limites de Bland-Altman es aproximado y "
            f"anticonservador con n < {N_BA_AVISO}.")
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
    """'cumple', 'no cumple', 'no evaluable' o 'no evaluado' (sin criterios)."""
    if not criterios:
        return "no evaluado", ["No se declararon criterios de aceptacion "
                               "antes de medir."]
    desconocidas = set(criterios) - set(CLAVES_CRITERIOS)
    if desconocidas:
        raise ValueError(f"criterios desconocidos: {sorted(desconocidas)}; "
                         f"admitidos: {CLAVES_CRITERIOS}")
    det, ok, indef = [], True, False

    def _eval(valor, cumple, texto):
        nonlocal ok, indef
        if valor is None or not np.isfinite(valor):
            indef = True
            det.append(f"no evaluable: {texto}")
            return
        ok = ok and bool(cumple)
        det.append(f"{'cumple' if cumple else 'NO cumple'}: {texto}")

    if "ccc_min" in criterios:
        # Se compara el limite inferior del IC, no el estimador puntual.
        ic = res.get("ccc_ic")
        v = ic[0] if ic else res["ccc"]
        que = ("limite inferior del IC del CCC" if ic
               else "CCC puntual (no hubo IC)")
        _eval(v, v >= criterios["ccc_min"] if v == v else False,
              f"{que} = {v:.3f} (minimo {criterios['ccc_min']})")
    if "sesgo_rel_max" in criterios:
        v = abs(res["sesgo_rel"])
        _eval(v, v <= criterios["sesgo_rel_max"],
              f"|sesgo relativo| = {v:.3f} (maximo {criterios['sesgo_rel_max']})")
    if "mape_max" in criterios:
        v = res["mape"]
        _eval(v, v <= criterios["mape_max"],
              f"error relativo medio absoluto = {v:.3f} "
              f"(maximo {criterios['mape_max']})")
    if "pendiente_equiv_delta" in criterios:
        # Equivalencia, no ausencia de diferencia: el IC ENTERO de la
        # pendiente debe caer dentro de [1 - delta, 1 + delta]. Un IC ancho
        # (pocos datos, mucho ruido) no puede cumplirlo.
        d = float(criterios["pendiente_equiv_delta"])
        lo, hi = res["regresion"]["pendiente_ic"]
        _eval(lo, lo >= 1 - d and hi <= 1 + d,
              f"IC de la pendiente [{lo:.3f}, {hi:.3f}] dentro de "
              f"[{1 - d:.3f}, {1 + d:.3f}]")
    if indef:
        return "no evaluable", det
    return ("cumple" if ok else "no cumple"), det


# ---------------------------------------------------------------------------
# Entrada / salida
# ---------------------------------------------------------------------------

COLUMNAS = ("id", "magnitud", "predicho", "medido")


def _num(texto, fila, id_, campo):
    t = (texto or "").strip().replace(",", ".")
    try:
        return float(t)
    except ValueError:
        raise ValueError(f"fila {fila} (id '{id_}'): '{campo}' no es un numero "
                         f"({texto!r})") from None


def leer_csv(ruta):
    """Lee el CSV de contraste y agrupa por magnitud.

    Columnas obligatorias: id, magnitud, predicho, medido. Acepta separador
    coma o punto y coma (Excel en español), coma decimal y BOM de UTF-8. Las
    filas con `medido` vacio se ignoran (es lo que trae la plantilla sin
    rellenar); si `medido` esta y `predicho` no, es un error con el numero de
    fila. Devuelve {magnitud: (ids, predicho, medido)}.
    """
    texto = Path(ruta).read_text(encoding="utf-8-sig")
    cabecera = texto.splitlines()[0] if texto else ""
    sep = ";" if cabecera.count(";") > cabecera.count(",") else ","
    grupos = {}
    lector = csv.DictReader(texto.splitlines(), delimiter=sep)
    falta = [c for c in COLUMNAS if c not in (lector.fieldnames or [])]
    if falta:
        raise ValueError(f"faltan columnas en el CSV: {falta}")
    for n, fila in enumerate(lector, start=2):
        if not (fila["medido"] or "").strip():
            continue
        id_ = fila["id"].strip()
        g = grupos.setdefault(fila["magnitud"].strip(), ([], [], []))
        g[0].append(id_)
        g[1].append(_num(fila["predicho"], n, id_, "predicho"))
        g[2].append(_num(fila["medido"], n, id_, "medido"))
    return grupos


def _a_json(x):
    """NaN e inf no son JSON valido: se escriben como null."""
    if isinstance(x, dict):
        return {k: _a_json(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_a_json(v) for v in x]
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def informe_md(resultados):
    """Resumen en Markdown de uno o varios contrastes."""
    L = ["# Contraste con ensayos fisicos", ""]
    for r in resultados:
        L.append(f"## {r['magnitud'] or 'magnitud'} (n = {r['n']})")
        L.append("")
        reg, ba = r["regresion"], r["bland_altman_log"]
        L.append(f"- RMSE {r['rmse']:.4g}; MAE {r['mae']:.4g}; error relativo "
                 f"medio absoluto {100 * r['mape']:.1f} %; sesgo relativo "
                 f"(predicho - medido) / medido {100 * r['sesgo_rel']:+.1f} %.")
        ic = r.get("ccc_ic")
        L.append(f"- CCC de Lin {r['ccc']:.3f}"
                 + (f" (IC {100 * (1 - r['alfa']):.0f} %: {ic[0]:.3f} a "
                    f"{ic[1]:.3f})" if ic else "")
                 + f"; Spearman {r['spearman']:.3f}; R2 {reg['r2']:.3f}.")
        L.append(f"- medido = {reg['pendiente']:.3f} x predicho "
                 f"{reg['intercepto']:+.4g} "
                 f"(IC {100 * (1 - r['alfa']):.0f} %, pendiente "
                 f"[{reg['pendiente_ic'][0]:.3f}, {reg['pendiente_ic'][1]:.3f}]).")
        if ba:
            L.append(f"- Razon medido / predicho (mayor que 1 si el modelo "
                     f"subestima): factor medio {ba['factor_sesgo']:.3f}; "
                     f"limites de acuerdo [{ba['limite_inf']:.3f}, "
                     f"{ba['limite_sup']:.3f}].")
        cal = r.get("calibracion")
        if cal:
            L.append(f"- Calibracion de E_s: k = {cal['k']:.3f} "
                     f"(geometrico {cal['k_geometrico']:.3f}); error relativo "
                     f"medio en validacion cruzada "
                     f"{100 * cal['error_rel_loo_medio_abs']:.1f} %"
                     + (f"; E_s calibrado {cal['Es_calibrado']:.4g}"
                        if "Es_calibrado" in cal else "") + ".")
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
    ap.add_argument("--ccc-min", type=float,
                    help="minimo del limite inferior del IC del CCC")
    ap.add_argument("--sesgo-rel-max", type=float)
    ap.add_argument("--mape-max", type=float)
    ap.add_argument("--pendiente-delta", type=float,
                    help="el IC de la pendiente debe caer en [1-d, 1+d]")
    ap.add_argument("--Es-usado", type=float, dest="es_usado",
                    help="E_s (MPa) con que se calculo el predicho; activa la "
                         "calibracion de E_s en las magnitudes E_app*")
    a = ap.parse_args(argv)

    crit = {}
    if a.ccc_min is not None:
        crit["ccc_min"] = a.ccc_min
    if a.sesgo_rel_max is not None:
        crit["sesgo_rel_max"] = a.sesgo_rel_max
    if a.mape_max is not None:
        crit["mape_max"] = a.mape_max
    if a.pendiente_delta is not None:
        crit["pendiente_equiv_delta"] = a.pendiente_delta

    grupos = leer_csv(a.csv)
    if not grupos:
        print("El CSV no tiene ninguna fila con 'medido' rellenado: no hay "
              "nada que contrastar todavia.", file=sys.stderr)
        return 2
    if crit and len(grupos) > 1:
        print("Aviso: los mismos criterios se aplican a todas las magnitudes "
              f"({', '.join(grupos)}).", file=sys.stderr)
    res = []
    for mag, (i, p, m) in grupos.items():
        r = acuerdo(p, m, ids=i, criterios=crit or None, magnitud=mag)
        if a.es_usado and mag.lower().startswith("e_app"):
            try:
                r["calibracion"] = calibrar_Es(p, m, Es_usado=a.es_usado)
            except ValueError as e:
                r["avisos"].append(f"Calibracion de E_s omitida: {e}")
        res.append(r)
    print(informe_md(res))
    if a.json:
        Path(a.json).write_text(
            json.dumps(_a_json(res), indent=2, ensure_ascii=False,
                       allow_nan=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
