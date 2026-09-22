"""
estadistica.py — Descomposicion de varianza y equivalencia real vs sintetico.

EL ERROR QUE ESTE MODULO EXISTE PARA IMPEDIR
---------------------------------------------
Generar K replicas de cada VOI y despues tratarlas como K especimenes
independientes es pseudorreplicacion (Hurlbert 1984). Infla los grados de
libertad, estrecha artificialmente los intervalos de confianza y convierte
cualquier diferencia en "significativa".

La varianza total se descompone en

    sigma^2_total = sigma^2_entre + sigma^2_dentro

donde `entre` es la variabilidad biologica entre especimenes y `dentro` la del
generador. Promediar K replicas divide sigma^2_dentro por K y deja
sigma^2_entre INTACTA. El N efectivo para una inferencia biologica lo fijan los
especimenes, no las realizaciones.

`n_efectivo` lo calcula explicitamente para que aparezca en cualquier informe y
sea dificil confundirlo con el numero de filas de la tabla.

POR QUE TOST Y NO UNA t DE STUDENT
-----------------------------------
Para sostener que los sinteticos SUSTITUYEN a los reales hace falta demostrar
EQUIVALENCIA, y "no se rechazo la hipotesis nula" no es eso: con N pequeno no
se rechaza nunca y con N grande se rechaza siempre, asi que la prueba clasica
responde a la pregunta contraria a la que interesa.

TOST invierte la carga: se declara de antemano un margen y se demuestra que la
diferencia cae dentro. Y ese margen no hay que inventarlo — el suelo de ruido
del propio generador, ya medido en Validacion_Anexo, es el limite fisico por
debajo del cual dos realizaciones de los MISMOS parametros tampoco se
distinguen. Pedir mas seria pedir que el sintetico se parezca al real mas de lo
que el real se parece a si mismo.
"""

from __future__ import annotations

import numpy as np


def descomponer_varianza(df, metrica, col_grupo="animal", col_replica="replica",
                         solo=None):
    """Componentes de varianza entre y dentro de espécimen, ICC y N efectivo.

    Modelo de efectos aleatorios de una via sobre los especimenes, con el
    estimador clasico de ANOVA (Searle):

        sigma^2_dentro = CM_dentro
        sigma^2_entre  = (CM_entre - CM_dentro) / n0

    con n0 el tamano medio armonico de grupo. `sigma^2_entre` puede salir
    NEGATIVA cuando la variabilidad entre especimenes es menor que el ruido; se
    acota a cero y se marca, porque un valor negativo no es interpretable pero
    tampoco es un fallo del calculo.
    """
    d = df if solo is None else df[solo]
    d = d[[col_grupo, col_replica, metrica]].dropna(subset=[metrica])
    grupos = [g[metrica].to_numpy(dtype=float) for _, g in d.groupby(col_grupo)]
    grupos = [g for g in grupos if g.size > 0]

    a = len(grupos)
    N = int(sum(g.size for g in grupos))
    if a < 2 or N <= a:
        return {"metrica": metrica, "n_grupos": a, "n_total": N,
                "msg": "Hacen falta al menos dos grupos con replicas."}

    medias = np.array([g.mean() for g in grupos])
    n_i = np.array([g.size for g in grupos], dtype=float)
    gran = float(np.concatenate(grupos).mean())

    SC_entre = float(np.sum(n_i * (medias - gran) ** 2))
    SC_dentro = float(sum(((g - g.mean()) ** 2).sum() for g in grupos))
    CM_entre = SC_entre / (a - 1)
    CM_dentro = SC_dentro / (N - a)

    n0 = (N - float(np.sum(n_i ** 2)) / N) / (a - 1)
    var_entre = (CM_entre - CM_dentro) / n0 if n0 > 0 else np.nan
    acotada = bool(np.isfinite(var_entre) and var_entre < 0)
    if acotada:
        var_entre = 0.0

    total = var_entre + CM_dentro
    icc = var_entre / total if total > 0 else np.nan

    return {
        "metrica": metrica,
        "n_grupos": a,                     # especimenes: la unidad biologica
        "n_total": N,                      # filas de la tabla
        "n_por_grupo_medio": float(np.mean(n_i)),
        "media": gran,
        "var_entre": float(var_entre),
        "var_dentro": float(CM_dentro),
        "sd_entre": float(np.sqrt(max(var_entre, 0))),
        "sd_dentro": float(np.sqrt(max(CM_dentro, 0))),
        "cv_entre_pct": float(100 * np.sqrt(max(var_entre, 0)) / abs(gran))
                        if gran else np.nan,
        "cv_dentro_pct": float(100 * np.sqrt(max(CM_dentro, 0)) / abs(gran))
                         if gran else np.nan,
        "ICC": float(icc),
        "n_efectivo": a,     # NO n_total: las replicas no anaden especimenes
        "var_entre_acotada_a_cero": acotada,
    }


def descomponer_anidada(df, metrica, col_animal="animal", col_sitio="sitio",
                        col_replica="replica", solo=None):
    """Descomposicion ANIDADA: animal / sitio(animal) / replica.

    POR QUE NO BASTA UN SOLO NIVEL. Agrupar por `animal` mete en el mismo saco
    la variacion entre SITIOS del mismo animal y el ruido del generador, que
    son cosas distintas: la primera es biologia y la segunda es numerica.
    Sobre los VOIs equinos —4 animales x 3 sitios x K replicas— esa mezcla
    inflaba el CV "dentro" hasta duplicar el "entre" y hacia el ICC
    ininterpretable.

    El diseno real tiene tres niveles y el modelo tiene que tener tres:

        sigma^2_total = sigma^2_animal + sigma^2_sitio(animal) + sigma^2_replica

    ANOVA anidada clasica (Sokal & Rohlf) sobre un diseno balanceado. Los
    componentes pueden salir negativos cuando un nivel aporta menos que el
    ruido del siguiente; se acotan a cero y se marca, porque un valor negativo
    no es interpretable pero tampoco es un fallo de calculo.

    El N para inferir entre animales sigue siendo el numero de ANIMALES.
    """
    d = df if solo is None else df[solo]
    d = d[[col_animal, col_sitio, metrica]].dropna(subset=[metrica])
    if d.empty:
        return {"metrica": metrica, "msg": "Sin datos."}

    animales = list(dict.fromkeys(d[col_animal]))
    a = len(animales)
    grupos = {}
    for an in animales:
        sub = d[d[col_animal] == an]
        grupos[an] = [g[metrica].to_numpy(float)
                      for _, g in sub.groupby(col_sitio)]

    s = int(np.mean([len(v) for v in grupos.values()]))
    r = float(np.mean([len(x) for v in grupos.values() for x in v]))
    N = int(sum(x.size for v in grupos.values() for x in v))
    if a < 2 or s < 1 or r < 2:
        return {"metrica": metrica, "n_animales": a, "n_total": N,
                "msg": "Diseno insuficiente para anidar."}

    todos = np.concatenate([x for v in grupos.values() for x in v])
    gran = float(todos.mean())

    med_an = {an: float(np.concatenate(v).mean()) for an, v in grupos.items()}
    SC_an = sum(np.concatenate(v).size * (med_an[an] - gran) ** 2
                for an, v in grupos.items())
    SC_si = sum(x.size * (x.mean() - med_an[an]) ** 2
                for an, v in grupos.items() for x in v)
    SC_re = sum(((x - x.mean()) ** 2).sum()
                for v in grupos.values() for x in v)

    gl_an = a - 1
    gl_si = sum(len(v) for v in grupos.values()) - a
    gl_re = N - sum(len(v) for v in grupos.values())
    if min(gl_an, gl_si, gl_re) < 1:
        return {"metrica": metrica, "n_animales": a, "n_total": N,
                "msg": "Grados de libertad insuficientes."}

    CM_an, CM_si, CM_re = SC_an / gl_an, SC_si / gl_si, SC_re / gl_re
    v_re = CM_re
    v_si = (CM_si - CM_re) / r
    v_an = (CM_an - CM_si) / (s * r)
    acot = [bool(v_si < 0), bool(v_an < 0)]
    v_si, v_an = max(v_si, 0.0), max(v_an, 0.0)
    tot = v_an + v_si + v_re

    cv = lambda v: 100 * np.sqrt(max(v, 0)) / abs(gran) if gran else np.nan
    return {
        "metrica": metrica, "media": gran,
        "n_animales": a, "n_sitios_por_animal": s,
        "n_replicas": int(round(r)), "n_total": N,
        "cv_animal_pct": cv(v_an), "cv_sitio_pct": cv(v_si),
        "cv_replica_pct": cv(v_re),
        "pct_var_animal": 100 * v_an / tot if tot else np.nan,
        "pct_var_sitio": 100 * v_si / tot if tot else np.nan,
        "pct_var_replica": 100 * v_re / tot if tot else np.nan,
        "n_efectivo_animal": a,       # para inferir entre animales
        "n_efectivo_especimen": a * s,  # para inferir entre especimenes
        "acotado_a_cero": acot,
    }


def tost(x, y, margen_rel, alfa=0.05):
    """Prueba de equivalencia por dos t unilaterales (Schuirmann).

    Concluye equivalencia si el intervalo de confianza (1-2*alfa) de la
    diferencia de medias cae ENTERO dentro de +-margen. `margen_rel` es
    relativo a la media de `x` (el grupo de referencia), que es como se declara
    de forma natural a partir de un coeficiente de variacion.
    """
    from scipy import stats

    x = np.asarray(x, float); x = x[np.isfinite(x)]
    y = np.asarray(y, float); y = y[np.isfinite(y)]
    if x.size < 2 or y.size < 2:
        return {"msg": "Hacen falta al menos dos observaciones por grupo."}

    mx, my = float(x.mean()), float(y.mean())
    margen = abs(margen_rel * mx)
    dif = my - mx
    se = float(np.sqrt(x.var(ddof=1) / x.size + y.var(ddof=1) / y.size))
    if se == 0:
        return {"equivalente": bool(abs(dif) <= margen), "dif": dif,
                "margen": margen, "msg": "Varianza nula en ambos grupos."}

    gl = (se ** 4) / ((x.var(ddof=1) / x.size) ** 2 / (x.size - 1)
                      + (y.var(ddof=1) / y.size) ** 2 / (y.size - 1))
    t1 = (dif + margen) / se
    t2 = (dif - margen) / se
    p1 = 1 - stats.t.cdf(t1, gl)      # H0: dif <= -margen
    p2 = stats.t.cdf(t2, gl)          # H0: dif >= +margen
    p = max(p1, p2)
    tc = stats.t.ppf(1 - alfa, gl)
    ic = (dif - tc * se, dif + tc * se)

    return {"n_x": int(x.size), "n_y": int(y.size),
            "media_x": mx, "media_y": my, "dif": dif,
            "dif_rel_pct": 100 * dif / mx if mx else np.nan,
            "margen": margen, "margen_rel_pct": 100 * margen_rel,
            "ic_inf": ic[0], "ic_sup": ic[1], "p_tost": float(p),
            "equivalente": bool(p < alfa),
            "gl": float(gl)}


def informe(df, metricas, margen_rel=None, col_grupo="animal"):
    """Descompone la varianza y contrasta equivalencia real vs sintetico.

    `margen_rel` por metrica; si falta, se usa el CV DENTRO observado, que es
    el suelo de ruido del propio generador y por tanto el margen mas exigente
    que tiene sentido pedir.
    """
    import pandas as pd

    reales = df["origen"] == "real"
    sint = df["origen"] == "sintetico"
    filas = []
    for met in metricas:
        if met not in df.columns:
            continue
        var = descomponer_varianza(df, met, col_grupo=col_grupo, solo=sint)
        mr = (margen_rel or {}).get(met)
        if mr is None:
            mr = max(var.get("cv_dentro_pct", np.nan) / 100
                     if np.isfinite(var.get("cv_dentro_pct", np.nan)) else 0.05,
                     0.02)
        eq = tost(df.loc[reales, met], df.loc[sint, met], mr)
        filas.append({**{k: var.get(k) for k in
                         ("metrica", "n_grupos", "n_total", "n_efectivo",
                          "cv_entre_pct", "cv_dentro_pct", "ICC")},
                      **{f"eq_{k}": eq.get(k) for k in
                         ("media_x", "media_y", "dif_rel_pct",
                          "margen_rel_pct", "p_tost", "equivalente")}})
    return pd.DataFrame(filas)
