"""
simulacion.py — Simulaciones in silico sobre el gemelo digital de un VOI.

Dos familias de experimento que no se pueden hacer sobre el animal sin
sacrificar uno por punto de la curva:

  1. PERDIDA OSEA en pasos controlados, con cuatro protocolos que se
     distinguen por QUE hueso se pierde primero.
  2. FALLO PROGRESIVO: se ablanda el tejido que rompe, se vuelve a cargar y
     se repite, hasta que la estructura no aguanta.

En cada paso se mide la morfometria completa y se resuelve un ensayo de
compresion (`resistencia.ensayo_compresion_eje`), de modo que cada curva es
rigidez y carga de fallo frente a masa osea y arquitectura.

LA REGLA QUE HACE COMPARABLES LOS PROTOCOLOS
--------------------------------------------
Todos retiran, en cada paso, la MISMA cantidad de hueso: una fraccion fija del
volumen oseo INICIAL. Lo unico que cambia es el orden en que se eligen los
voxeles. Asi, a igual paso, dos protocolos tienen exactamente el mismo BV/TV,
y cualquier diferencia de rigidez entre ellos es de ARQUITECTURA, no de masa.
Es la comparacion que interesa: la literatura clinica repite que la densidad
sola no predice la fractura, y esto permite medir cuanto pesa lo demas.

    protocolo            se retira primero                   se parece a
    -------------------  ----------------------------------  ----------------------
    adelgazamiento       lo mas cercano a la medula          envejecimiento, perdida
                         (distancia al poro)                 lenta y uniforme
    trabeculas_finas     lo de menor espesor local           OVX, remodelado acelerado:
                         (Hildebrand-Ruegsegger)             las trabeculas finas se
                                                             perforan y desaparecen
    desuso               la superficie menos deformada en    inmovilizacion,
                         el ensayo de ESE paso               microgravedad: el tejido
                         (idea del mecanostato)              que no trabaja se reabsorbe
    recuperacion         perdida por trabeculas finas y      farmaco anabolico tras la
                         despues engrosamiento desde la      perdida
                         superficie hasta la masa inicial

La RECUPERACION es la que da el resultado con mas lectura clinica: se vuelve a
la MISMA masa osea y se compara la rigidez con la de partida. Si no se
recupera, la diferencia es la arquitectura que se perdio para siempre: una
trabecula que desaparece no vuelve por engrosar las que quedan.

LO QUE ESTO NO ES
-----------------
No es un modelo de remodelado calibrado. No hay tasas, ni celulas, ni tiempo:
un "paso" es una cantidad de hueso, no un mes. Los protocolos son REGLAS DE
ORDEN que capturan la firma geometrica de cada situacion, y lo que se puede
afirmar con ellos es comparativo —este protocolo cuesta mas rigidez que aquel a
igual masa—, no una prediccion de cuanto hueso pierde un paciente en un ano.

DOS RESOLUCIONES, Y POR QUE
--------------------------
  * La GEOMETRIA se modifica a la resolucion de trabajo (`n_trabajo`, 96 por
    omision o la del VOI si es menor). Por debajo no hay voxeles suficientes
    para retirar un 5 % con orden: a 32^3 una trabecula tiene uno o dos
    voxeles de grosor y "adelgazarla" es borrarla.
  * La MECANICA se resuelve a `n_mec` (32 por omision), remuestreando la
    geometria de cada paso por vecino mas proximo, como todo el proyecto.
    Los valores absolutos de E_app a esa resolucion no son citables (ver
    `resistencia.estudio_convergencia`); los COCIENTES respecto al paso 0, si,
    porque todos los pasos comparten el sesgo.

Para el desuso, el campo de deformacion de la malla mecanica se lleva a la de
trabajo voxel a voxel y se dilata con un maximo 3x3x3, la misma idea que
`espesor.muestrear_en_puntos`: sin eso, el hueso que el remuestreo dejo fuera
de la malla mecanica tendria "deformacion desconocida" y se retiraria el
primero por un artefacto, no por no trabajar.

DETERMINISMO
------------
Los empates de puntuacion (todos los voxeles de superficie tienen la misma
distancia al poro) se rompen con un ruido minimo de semilla fija. Misma
entrada y misma semilla dan exactamente la misma secuencia de mascaras.
"""

from __future__ import annotations

import time

import numpy as np
from scipy import ndimage

from .elastic import remuestrear_bw
from .espesor import espesor_local
from .morphometry import morfometria
from .resistencia import (E_S_DEF, EPS_CRITICA, FRAC_CRITICA, NU_S_DEF,
                          criterio_pistoia, ensayo_compresion_eje)

PROTOCOLOS = ("adelgazamiento", "trabeculas_finas", "desuso", "recuperacion")
N_TRABAJO_DEF = 96
N_MEC_DEF = 32
SEMILLA_DEF = 20260720

# Por debajo de esta fraccion de hueso no se sigue: no queda estructura que
# ensayar, y la morfometria de un puñado de voxeles no significa nada.
BVTV_MIN = 0.03

# `ensayo_compresion` da F_total = sigma0 * A_bruta con sigma0 en Pa y el area
# en unidades de `spacing` al cuadrado. Con el spacing en MM (unidad interna
# del proyecto) eso son Pa*mm2 = 1e-6 N. Sin este factor la carga de fallo
# sale un millon de veces mayor: medido en la primera prueba del dialogo,
# "61 MN" para un cubo de 4 mm.
PA_MM2_A_N = 1e-6

# Rigidez que conserva un elemento roto en el fallo progresivo. Se ABLANDA en
# vez de borrarse: borrar puntales de uno o dos elementos dejaba fragmentos
# casi como mecanismos y la rigidez del VOI proximal de H4 subia y bajaba entre
# pasos (1 -> 0.15 -> 0.07 -> 0.001 -> 0.11), cosa que un dano no puede hacer.
# El 5 % es una convencion, no una medida del tejido: se declara y viaja con el
# resultado.
RIGIDEZ_DANADA = 0.05

# Rigidez relativa por debajo de la cual la estructura se da por COLAPSADA y la
# carga de fallo deja de leerse. Hace falta porque, una vez rota una banda que
# cruza la seccion, la carga que calcula el criterio SIGUE SUBIENDO paso a paso
# —la sostiene el 5 % residual del tejido ablandado, no resistencia real—.
# Medido en el VOI proximal de H4: la rigidez cae a 0.15 con un 2 % de dano y
# la "carga de fallo" pasa de 165 N a 215 N con la rigidez en 0.07. El 50 % es
# una convencion declarada, no una medida.
E_REL_COLAPSO = 0.5

# Metricas escalares que se guardan por paso. Conn.D y SMI necesitan
# `extra=True`, que a 96^3 cuesta unos segundos: aqui se pagan porque la
# perdida de CONECTIVIDAD es justo lo que distingue a los protocolos.
METRICAS_PASO = ("BVTV", "TbTh", "TbSp", "TbN", "BSBV", "DA", "ConnD", "SMI")


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------

def _progreso(cb, i, n, msg):
    if cb:
        cb(int(i), int(n), str(msg))


def _desempate(forma, semilla):
    """Ruido minimo y reproducible para romper empates de puntuacion."""
    rng = np.random.default_rng(int(semilla))
    return rng.random(forma).astype(np.float64) * 1e-9


def _retirar(BW, puntuacion, k):
    """Quita de BW los k voxeles solidos de MENOR puntuacion."""
    BW = BW.copy()
    idx = np.flatnonzero(BW)
    k = int(min(max(k, 0), idx.size))
    if k == 0:
        return BW, 0
    p = puntuacion.ravel()[idx]
    sel = idx[np.argpartition(p, k - 1)[:k]]
    BW.ravel()[sel] = False
    return BW, k


def _anadir(BW, puntuacion, k):
    """Anade a BW los k voxeles VACIOS de menor puntuacion."""
    BW = BW.copy()
    idx = np.flatnonzero(~BW)
    k = int(min(max(k, 0), idx.size))
    if k == 0:
        return BW, 0
    p = puntuacion.ravel()[idx]
    sel = idx[np.argpartition(p, k - 1)[:k]]
    BW.ravel()[sel] = True
    return BW, k


def _indices_remuestreo(n_fino, n_grueso):
    """Para cada indice de la rejilla fina, el de la gruesa mas proximo.

    `remuestrear_bw` toma en la gruesa los indices round(linspace(0, n-1, m))
    de la fina. Aqui se hace el camino inverso: cada voxel fino se asigna a la
    muestra gruesa cuyo indice fino queda mas cerca.
    """
    muestras = np.round(np.linspace(0, n_fino - 1, n_grueso)).astype(int)
    fino = np.arange(n_fino)
    j = np.searchsorted(muestras, fino)
    j = np.clip(j, 1, n_grueso - 1)
    izq, der = muestras[j - 1], muestras[j]
    return np.where(np.abs(fino - izq) <= np.abs(der - fino), j - 1, j)


def campo_a_rejilla_fina(campo_grueso, forma_fina):
    """Lleva un campo por elemento de la malla mecanica a la de trabajo.

    Los NaN (material no portante o fuera de la malla) valen -1, y se dilata
    con un maximo 3x3x3 ANTES de llevarlo a la rejilla fina: asi un voxel fino
    que el remuestreo dejo fuera recoge la deformacion de su vecino portante.
    Solo queda en -1 lo que no tiene ningun vecino que trabaje, que es
    precisamente lo que no transmite carga.
    """
    c = np.nan_to_num(np.asarray(campo_grueso, float), nan=-1.0)
    c = ndimage.maximum_filter(c, size=3, mode="nearest")
    ix = [_indices_remuestreo(forma_fina[a], c.shape[a]) for a in range(3)]
    return c[np.ix_(ix[0], ix[1], ix[2])]


def _medir(BW, sp):
    m = morfometria(BW, sp, extra=True)
    return {k: (float(m[k]) if np.ndim(m.get(k)) == 0 and m.get(k) is not None
                and np.isfinite(m[k]) else float("nan"))
            for k in METRICAS_PASO} | {
        "dir_principal": [float(x) for x in m.get("dir_principal", [0, 0, 1])]}


def _ensayar(BW, sp, n_mec, eje, E_s, nu_s, apoyo, frac, eps_crit):
    """Ensayo de compresion y criterio de Pistoia a la resolucion mecanica."""
    bw, spr = remuestrear_bw(BW, sp, int(n_mec))
    r = ensayo_compresion_eje(bw, spr, eje=int(eje), E_s=E_s, nu_s=nu_s,
                              apoyo=apoyo)
    out = {"ok": bool(r.get("ok")), "msg": r.get("msg", ""),
           "E_app": float("nan"), "sigma_fallo": float("nan"),
           "F_fallo": float("nan"), "frac_portante": r.get("frac_portante"),
           "residuo": r.get("residuo_rel"), "n_mec": int(bw.shape[0])}
    if not r.get("ok"):
        # Sin camino portante la estructura no aguanta carga: su rigidez
        # aparente es CERO, no desconocida. Se marca aparte para no mezclarlo
        # con un fallo del solver.
        out["sin_camino"] = any(t in out["msg"] for t in
                                ("base y techo", "capa superior",
                                 "no contiene material"))
        return out, None, bw
    p = criterio_pistoia(r, frac=frac, eps_crit=eps_crit)
    out.update({"E_app": float(r["E_app"]),
                "sigma_fallo": float(p.get("sigma_fallo", np.nan)),
                "F_fallo": float(p.get("F_fallo", np.nan)) * PA_MM2_A_N,
                "factor_pistoia": float(p.get("factor", np.nan))})
    return out, r, bw


# ---------------------------------------------------------------------------
# Puntuaciones: QUE se pierde primero en cada protocolo
# ---------------------------------------------------------------------------

def puntuacion_adelgazamiento(BW, sp):
    """Distancia al poro: la superficie sale primero, capa a capa."""
    return ndimage.distance_transform_edt(BW, sampling=sp)


def puntuacion_trabeculas_finas(BW, sp):
    """Espesor local de la trabecula que contiene el voxel.

    Todos los voxeles de un puntal fino comparten un espesor pequeno, asi que
    salen JUNTOS: el puntal desaparece entero en vez de adelgazar. Dentro del
    mismo espesor se desempata por distancia al poro, para que un puntal a
    medio retirar pierda antes la piel que el eje.
    """
    esp = espesor_local(BW, sp)
    d = ndimage.distance_transform_edt(BW, sampling=sp)
    return esp + 1e-3 * d / max(float(d.max()), 1e-12) * float(np.min(sp))


def puntuacion_desuso(BW, sp, campo_eps_fino):
    """Deformacion efectiva, con prioridad absoluta para la SUPERFICIE.

    El remodelado ocurre en las superficies oseas, no dentro del tejido. Se
    ordena primero por capa —superficie (distancia al poro <= 1.5 voxeles)
    antes que interior— y dentro de cada capa por deformacion creciente: lo
    menos cargado de la superficie se va primero. Solo cuando la superficie
    poco cargada se agota se entra en el interior.
    """
    d = ndimage.distance_transform_edt(BW, sampling=sp)
    superficie = d <= 1.5 * float(np.min(sp))
    e = np.asarray(campo_eps_fino, float)
    validos = e[BW & (e > 0)]
    escala = float(np.percentile(validos, 99)) if validos.size else 1.0
    rango = np.clip(e / max(escala, 1e-30), -1.0, 1.0)      # [-1, 1]
    return rango + 4.0 * (~superficie)


# ---------------------------------------------------------------------------
# Perdida osea
# ---------------------------------------------------------------------------

def simular_perdida(BW, spacing, protocolo="adelgazamiento", pasos=8,
                    perdida_paso=0.05, n_trabajo=N_TRABAJO_DEF,
                    n_mec=N_MEC_DEF, eje=2, E_s=E_S_DEF, nu_s=NU_S_DEF,
                    apoyo="deslizante", semilla=SEMILLA_DEF,
                    frac=FRAC_CRITICA, eps_crit=EPS_CRITICA,
                    guardar_mascaras=True, progreso=None):
    """Perdida osea en `pasos` pasos de `perdida_paso` del hueso INICIAL.

    protocolo : 'adelgazamiento' | 'trabeculas_finas' | 'desuso' |
                'recuperacion' (perdida por trabeculas finas y engrosamiento
                de vuelta a la masa inicial en el mismo numero de pasos)

    Devuelve un dict con
        pasos      lista de dicts, uno por paso (el 0 es la estructura intacta)
        mascaras   lista de mascaras a la resolucion de trabajo (si se pide)
        spacing    spacing de esas mascaras, en mm
        resumen    cocientes finales y pendiente log-log de E frente a BV/TV
    """
    if protocolo not in PROTOCOLOS:
        raise ValueError(f"protocolo desconocido: {protocolo!r}; "
                         f"use uno de {PROTOCOLOS}")
    t0 = time.time()
    BW = np.asarray(BW, dtype=bool)
    sp = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if sp.size == 1:
        sp = np.repeat(sp, 3)

    BW, sp = remuestrear_bw(BW, sp, int(min(n_trabajo, max(BW.shape))))
    N0 = int(BW.sum())
    if N0 == 0:
        raise ValueError("La estructura no contiene hueso.")
    k_paso = max(1, int(round(float(perdida_paso) * N0)))
    ruido = _desempate(BW.shape, semilla)

    fases = [("perdida", int(pasos))]
    if protocolo == "recuperacion":
        fases.append(("recuperacion", int(pasos)))
    total = 1 + sum(n for _, n in fases)

    filas, mascaras = [], []
    actual = BW
    ultimo_ensayo = None

    def registrar(i, fase, retirados):
        nonlocal ultimo_ensayo
        _progreso(progreso, len(filas), total,
                  f"paso {i} ({fase}): midiendo y ensayando")
        m = _medir(actual, sp)
        mec, r, bw_mec = _ensayar(actual, sp, n_mec, eje, E_s, nu_s, apoyo,
                                  frac, eps_crit)
        ultimo_ensayo = (r, bw_mec)
        fila = {"paso": i, "fase": fase, "voxeles": int(actual.sum()),
                "hueso_rel": float(actual.sum() / N0),
                "cambio_voxeles": int(retirados), **m, **mec}
        filas.append(fila)
        if guardar_mascaras:
            mascaras.append(actual.copy())

    registrar(0, "inicial", 0)
    i = 0
    for fase, n in fases:
        for _ in range(n):
            i += 1
            if fase == "perdida":
                if actual.mean() <= BVTV_MIN:
                    break
                if protocolo == "adelgazamiento":
                    punt = puntuacion_adelgazamiento(actual, sp)
                elif protocolo in ("trabeculas_finas", "recuperacion"):
                    punt = puntuacion_trabeculas_finas(actual, sp)
                else:
                    r, bw_mec = ultimo_ensayo
                    if r is None:
                        # Ya no hay camino portante: nada "trabaja", y la
                        # regla del mecanostato no distingue. Se sigue
                        # adelgazando para completar la curva.
                        punt = puntuacion_adelgazamiento(actual, sp)
                    else:
                        fino = campo_a_rejilla_fina(r["campo_eps_eff"],
                                                    actual.shape)
                        punt = puntuacion_desuso(actual, sp, fino)
                actual, hecho = _retirar(actual, punt + ruido, k_paso)
                registrar(i, "perdida", -hecho)
            else:
                faltan = N0 - int(actual.sum())
                if faltan <= 0:
                    break
                # Engrosamiento: se anade el poro mas cercano al hueso, capa
                # a capa. Nunca se crean trabeculas nuevas en medio de la
                # medula, que es lo que hace un anabolico real sobre las
                # superficies existentes.
                d_poro = ndimage.distance_transform_edt(~actual, sampling=sp)
                actual, hecho = _anadir(actual, d_poro + ruido,
                                        min(k_paso, faltan))
                registrar(i, "recuperacion", hecho)

    out = {"protocolo": protocolo, "pasos": filas,
           "mascaras": mascaras if guardar_mascaras else [],
           "spacing": [float(x) for x in sp],
           "parametros": {"pasos": int(pasos), "perdida_paso": float(perdida_paso),
                          "n_trabajo": int(BW.shape[0]), "n_mec": int(n_mec),
                          "eje": int(eje), "apoyo": apoyo, "E_s": float(E_s),
                          "nu_s": float(nu_s), "semilla": int(semilla),
                          "frac": float(frac), "eps_crit": float(eps_crit)},
           "tiempo_s": round(time.time() - t0, 1)}
    out["resumen"] = resumir(filas)
    out["resumen"]["TbTh_h_mec"] = tbth_por_elemento(
        filas[0].get("TbTh"), BW.shape[0] * float(sp[0]), filas[0].get("n_mec"))
    _progreso(progreso, total, total, "listo")
    return out


# Elementos de la malla mecanica por trabecula por debajo de los cuales los
# cocientes de rigidez pueden ser artefacto. El estudio de convergencia midio
# que la resolucion efectiva es Tb.Th/h, no n, y que el VOI proximal de H4
# converge a +-3 % desde Tb.Th/h ~1.7 (n = 40-48).
TBTH_H_MIN = 1.7


def tbth_por_elemento(tbth_mm, lado_mm, n_mec):
    """Tb.Th / h de la malla mecanica; NaN si falta algun dato."""
    try:
        h = float(lado_mm) / float(n_mec)
        v = float(tbth_mm) / h
    except (TypeError, ValueError, ZeroDivisionError):
        return float("nan")
    return v if np.isfinite(v) else float("nan")


def resumir(filas):
    """Cocientes respecto al paso 0 y pendiente log-log de E frente a BV/TV.

    La pendiente es DESCRIPTIVA. Estudio_Percolacion explica por que no debe leerse como
    el exponente de Gibson-Ashby: cerca del umbral de rigidez el exponente
    aparente se infla. Sirve para comparar protocolos entre si sobre la misma
    estructura, que es lo que aqui se hace.
    """
    if not filas:
        return {}
    f0 = filas[0]
    E0, S0 = f0.get("E_app"), f0.get("sigma_fallo")
    for f in filas:
        f["E_rel"] = (f["E_app"] / E0 if E0 and np.isfinite(E0) and E0 > 0
                      and np.isfinite(f["E_app"]) else
                      (0.0 if f.get("sin_camino") else float("nan")))
        f["sigma_fallo_rel"] = (f["sigma_fallo"] / S0 if S0 and np.isfinite(S0)
                                and S0 > 0 and np.isfinite(f["sigma_fallo"])
                                else (0.0 if f.get("sin_camino") else float("nan")))
    per = [f for f in filas if f["fase"] in ("inicial", "perdida")
           and f["E_rel"] > 0 and np.isfinite(f["E_rel"])]
    res = {"E_rel_final": filas[-1]["E_rel"],
           "BVTV_rel_final": filas[-1]["BVTV"] / f0["BVTV"] if f0["BVTV"] else float("nan"),
           "ConnD_rel_final": (filas[-1]["ConnD"] / f0["ConnD"]
                               if f0.get("ConnD") not in (None, 0) and np.isfinite(f0["ConnD"])
                               else float("nan")),
           "sin_camino_desde": next((f["paso"] for f in filas if f.get("sin_camino")), None)}
    if len(per) >= 3:
        x = np.log([f["BVTV"] for f in per])
        y = np.log([f["E_rel"] for f in per])
        if np.ptp(x) > 0:
            res["pendiente_loglog"] = float(np.polyfit(x, y, 1)[0])
    rec = [f for f in filas if f["fase"] == "recuperacion"]
    if rec:
        # El minimo entre los pasos de perdida VALIDOS. Tomar el ultimo paso
        # tal cual daba NaN en cuanto ese paso no convergia, y el dialogo
        # decia "tocaba fondo en nan" (primera prueba del dialogo).
        fondo = [f["E_rel"] for f in filas if f["fase"] == "perdida"
                 and np.isfinite(f["E_rel"])]
        res["E_rel_minimo"] = float(min(fondo)) if fondo else float("nan")
        res["E_rel_recuperado"] = rec[-1]["E_rel"]
        res["hueso_rel_recuperado"] = rec[-1]["hueso_rel"]
    return res


# ---------------------------------------------------------------------------
# Fallo progresivo
# ---------------------------------------------------------------------------

def fallo_progresivo(BW, spacing, pasos=10, n_mec=N_MEC_DEF, eje=2,
                     E_s=E_S_DEF, nu_s=NU_S_DEF, apoyo="deslizante",
                     frac=FRAC_CRITICA, eps_crit=EPS_CRITICA,
                     rigidez_danada=RIGIDEZ_DANADA,
                     guardar_mascaras=True, progreso=None):
    """Ablanda en cada paso el tejido que el criterio de Pistoia da por roto.

    COMO
      1. Se resuelve el ensayo lineal. Sobre el tejido todavia INTACTO se busca
         el factor de carga k al que la fraccion `frac` alcanza `eps_crit`
         (Pistoia). El tejido ya roto no cuenta: ya fallo.
      2. Los elementos intactos con k * eps_eff >= eps_crit son los que rompen
         a esa carga: su rigidez pasa a `rigidez_danada` * E_s. No se borran,
         para que la estructura siga conectada y la rigidez solo pueda bajar.
      3. Se vuelve a resolver la estructura danada. Repetir.

    Las mascaras guardadas son el tejido INTACTO de cada paso, que es lo que
    tiene sentido mirar en 3D: el roto sigue en la malla, pero ya no trabaja.

    La carga de fallo de cada paso es la de la estructura YA danada. Si el
    dano se reparte, la carga puede mantenerse o subir unos pasos y la maxima
    estima la resistencia ultima. Si se concentra en una banda que cruza la
    seccion, la rigidez cae de golpe: es la firma de una estructura fragil, y
    la estructura se da por colapsada (`E_REL_COLAPSO`). Despues del colapso la
    carga calculada vuelve a subir por el tejido ablandado, asi que `F_max` se
    toma solo entre los pasos anteriores y los demas se marcan `tras_colapso`.

    Se hace entero a la resolucion mecanica: el dano es un conjunto de
    ELEMENTOS de esa malla, y llevarlo a otra rejilla solo anadiria
    interpolacion sin informacion nueva.

    LIMITES DECLARADOS. Es elastico lineal por tramos con ablandamiento de
    elementos en un solo escalon: no hay plasticidad, ni dano gradual, ni
    contacto. La rigidez residual del 5 % es una convencion. Los dos parametros de Pistoia se calibraron en radio distal
    humano. La serie sirve para comparar estructuras y protocolos, no para dar
    una carga de rotura absoluta.
    """
    t0 = time.time()
    BW = np.asarray(BW, dtype=bool)
    sp = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if sp.size == 1:
        sp = np.repeat(sp, 3)
    actual, spr = remuestrear_bw(BW, sp, int(n_mec))
    N0 = int(actual.sum())
    if N0 == 0:
        raise ValueError("La estructura no contiene hueso.")
    tbth0 = float(morfometria(actual, spr, do_mil=False).get("TbTh", np.nan))

    modulo = np.ones(actual.shape, dtype=np.float64)
    danado = np.zeros(actual.shape, dtype=bool)
    filas, mascaras = [], []
    for i in range(int(pasos) + 1):
        _progreso(progreso, i, pasos + 1, f"paso {i}: ensayando")
        r = ensayo_compresion_eje(actual, spr, eje=int(eje), E_s=E_s,
                                  nu_s=nu_s, apoyo=apoyo, modulo_rel=modulo)
        intactos = int((actual & ~danado).sum())
        fila = {"paso": i, "voxeles": intactos,
                "dano_acumulado": float(1.0 - intactos / N0),
                "BVTV": float(actual.mean()), "ok": bool(r.get("ok")),
                "msg": r.get("msg", ""), "frac_portante": r.get("frac_portante"),
                "residuo": r.get("residuo_rel"), "E_app": float("nan"),
                "sigma_fallo": float("nan"), "F_fallo": float("nan"),
                "rotos": 0}
        if guardar_mascaras:
            mascaras.append(actual & ~danado)
        if not r.get("ok"):
            filas.append(fila)
            break
        campo = np.asarray(r["campo_eps_eff"], float)
        sano = actual & ~danado & np.isfinite(campo)
        e = campo[sano]
        if e.size == 0:
            fila["ok"] = False
            fila["msg"] = "No queda tejido intacto que cargar."
            filas.append(fila)
            break
        # Pistoia sobre el tejido intacto: el percentil (1 - frac) y el factor
        # que lo lleva a eps_crit. Mismo calculo que `criterio_pistoia`, pero
        # sin contar el tejido que ya rompio.
        umbral = float(np.percentile(e, 100.0 * (1.0 - frac)))
        k = eps_crit / umbral if umbral > 0 else float("nan")
        fila.update({"E_app": float(r["E_app"]),
                     "sigma_fallo": float(k * r["sigma_app"]),
                     "F_fallo": float(k * r["F_total"]) * PA_MM2_A_N})
        rotos = sano & (k * np.nan_to_num(campo, nan=0.0)
                        >= eps_crit * (1.0 - 1e-12))
        fila["rotos"] = int(rotos.sum())
        filas.append(fila)
        if i < int(pasos):
            danado |= rotos
            modulo[rotos] = float(rigidez_danada)

    E0 = filas[0]["E_app"]
    F0 = filas[0]["F_fallo"]
    for f in filas:
        f["E_rel"] = f["E_app"] / E0 if np.isfinite(f["E_app"]) and E0 else (
            0.0 if not f["ok"] else float("nan"))
        f["F_rel"] = f["F_fallo"] / F0 if np.isfinite(f["F_fallo"]) and F0 else (
            0.0 if not f["ok"] else float("nan"))
    # Colapso: primer paso sin solucion o con la rigidez por debajo del
    # umbral. La carga de fallo solo se lee ANTES de el (ver E_REL_COLAPSO).
    colapso = next((f["paso"] for f in filas if not f["ok"]
                    or (np.isfinite(f["E_rel"]) and f["E_rel"] < E_REL_COLAPSO)),
                   None)
    for f in filas:
        f["tras_colapso"] = colapso is not None and f["paso"] >= colapso
    validas = [f for f in filas if f["ok"] and np.isfinite(f["F_fallo"])
               and not f["tras_colapso"]]
    pico = max(validas, key=lambda f: f["F_fallo"]) if validas else None
    resumen = {
        "F_max": pico["F_fallo"] if pico else float("nan"),
        "paso_F_max": pico["paso"] if pico else None,
        "F_max_rel": pico["F_rel"] if pico else float("nan"),
        "dano_en_F_max": pico["dano_acumulado"] if pico else float("nan"),
        "E_rel_final": filas[-1]["E_rel"],
        "colapso_en_paso": colapso,
        "E_rel_colapso": E_REL_COLAPSO,
        "TbTh_h_mec": tbth0 / float(spr[0]) if np.isfinite(tbth0) else float("nan"),
    }
    return {"tipo": "fallo_progresivo", "pasos": filas,
            "mascaras": mascaras if guardar_mascaras else [],
            "spacing": [float(x) for x in spr],
            "parametros": {"pasos": int(pasos), "n_mec": int(actual.shape[0]),
                           "eje": int(eje), "apoyo": apoyo, "E_s": float(E_s),
                           "nu_s": float(nu_s), "frac": float(frac),
                           "eps_crit": float(eps_crit),
                           "rigidez_danada": float(rigidez_danada)},
            "resumen": resumen, "tiempo_s": round(time.time() - t0, 1)}
