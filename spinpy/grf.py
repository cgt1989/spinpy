"""
grf.py — Generador de spinodoides VOXEL-PRIMERO.

Port a Python del campo aleatorio gaussiano (GRF) de Kumar et al. (2020), tal
como lo implementa `spinodoid.m` de GIBBON, con las dos convenciones de
muestreo de direcciones de onda que coexisten en la literatura y en el codigo
disponible en este proyecto.

POR QUE VOXEL-PRIMERO
---------------------
`AppFinal_V2.m` genera la malla triangular con `spinodoid()`, la escala al
tamano fisico del VOI y despues la voxeliza (`localVoxelizeMesh`) para medirla
con `localMorphometryFromVoxels`. El repositorio del profesor
(TPMS-Scaffolds-generator) va aun mas lejos: malla tetraedrica con CGAL y
metricas sobre la malla.

Para un bucle de ajuste ese camino es innecesariamente caro: mallar en cada
evaluacion domina el coste. El campo GRF ya esta definido sobre una rejilla, y
la mascara booleana sale de un umbral. Aqui se corta ahi: rejilla -> umbral ->
mascara. La morfometria (morphometry.py) opera sobre la mascara, que es
exactamente la representacion que la app usa para comparar VOI y candidato.

LAS DOS CONVENCIONES DE MUESTREO
--------------------------------
Los `thetas` son semiangulos de cono alrededor de los ejes x, y, z. Un vector
de onda es admisible si cae dentro de alguno de los conos activos. Como
repartir `numWaves` vectores entre los conos NO esta fijado por la teoria, y
las dos implementaciones disponibles resuelven distinto:

  'rechazo'    (GIBBON / AppFinal_V2, spinodoid.m:178-206)
      Se muestrea un candidato isotropo en la esfera y se acepta si cae en
      CUALQUIER cono. El reparto entre conos resulta proporcional al angulo
      solido de cada uno: un cono ancho recibe muchas mas ondas que uno
      estrecho.

  'equitativo' (TPMS-Scaffolds-generator, src/TPMS.py:223)
      `numWaves // n_conos_activos` ondas por cono, cada una muestreada
      uniformemente dentro de su propio cono. El reparto es 50/50 aunque los
      conos tengan angulos muy distintos.

Coinciden solo si todos los conos activos tienen el mismo angulo (o hay uno
solo, o el caso es isotropo). Con thetas=(15,45,0) el rechazo manda ~90% de
las ondas al cono de 45 grados y el equitativo un 50%. Producen anisotropias
DISTINTAS para los mismos thetas nominales; ver validar_analitico.py.

Ninguna es "la correcta": son dos definiciones de la misma familia. Lo que no
se puede es cruzar resultados entre ambas sin declarar cual se uso.

Convencion de ejes: la mascara devuelta tiene dimension 1 = X, igual que
`readVTKVOI` en AppFinal_V2.m y que `voxelize_to_mask` en el repositorio del
profesor.
"""

from __future__ import annotations

import numpy as np
from scipy.special import erfinv as _erfinv_exacta


# ---------------------------------------------------------------------------
# Umbral del level-set
# ---------------------------------------------------------------------------

def _erfinv_winitzki(y: float) -> float:
    """Aproximacion de Winitzki a erf^-1 (a = 0.147).

    Es la que usa src/TPMS.py del repositorio del profesor, donde aparece con
    el nombre `erf` (esta bien empleada, pero el nombre confunde: es la
    inversa). Error relativo del orden de 2e-3. Se incluye para poder
    reproducir sus numeros exactamente; por defecto se usa la exacta de scipy.
    """
    a = 0.147
    s = np.sign(y)
    ln = np.log(1.0 - y * y)
    first = 2.0 / (np.pi * a) + ln / 2.0
    return float(s * np.sqrt(np.sqrt(first * first - ln / a) - first))


def level_set(rho: float, erfinv: str = "exacta") -> float:
    """Umbral phi0 tal que P(GRF <= phi0) = rho para un GRF N(0,1).

    phi0 = sqrt(2) * erf^-1(2*rho - 1)     [Kumar et al. 2020]

    Es identico en GIBBON (`spinodoid.m:264`) y en el repositorio del profesor
    (`TPMS.py:272`, salvo por la aproximacion de erfinv).
    """
    if not (0.0 < rho < 1.0):
        raise ValueError(f"rho debe estar en (0,1); se recibio {rho}")
    if erfinv == "exacta":
        return float(np.sqrt(2.0) * _erfinv_exacta(2.0 * rho - 1.0))
    if erfinv == "winitzki":
        return float(np.sqrt(2.0) * _erfinv_winitzki(2.0 * rho - 1.0))
    raise ValueError("erfinv debe ser 'exacta' o 'winitzki'")


# ---------------------------------------------------------------------------
# Muestreo de direcciones de onda
# ---------------------------------------------------------------------------

def _waves_rechazo(num_waves, thetas, R, rng, lote=4096):
    """Muestreo por rechazo, convencion GIBBON.

    Replica `spinodoid.m:178-206`: candidato isotropo `randn(1,3)` normalizado,
    aceptado si el angulo minimo respecto de ALGUN eje rotado es menor que el
    theta de ese eje. Se vectoriza por lotes; el bucle `while` original acepta
    de uno en uno, lo que es equivalente en distribucion.
    """
    thetas = np.asarray(thetas, dtype=float)
    ejes = R @ np.eye(3)          # columnas = ejes rotados
    aceptados = []
    n = 0
    while n < num_waves:
        c = rng.standard_normal((lote, 3))
        c /= np.linalg.norm(c, axis=1, keepdims=True)
        # angulo al eje i considerando +eje y -eje (de ahi el valor absoluto)
        cosang = np.abs(c @ ejes)                       # (lote, 3)
        ang = np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))
        ok = np.any(ang < thetas[None, :], axis=1)
        v = c[ok]
        if v.size:
            aceptados.append(v)
            n += len(v)
    return np.vstack(aceptados)[:num_waves]


def _muestrear_cono(eje: int, theta_deg: float, n: int, rng):
    """Muestreo uniforme dentro de un cono, convencion del profesor.

    Replica `src/TPMS.py:233-254` (`_sample_cone`), incluida la rotacion del
    cono del eje z al eje pedido y el signo aleatorio que cubre las dos ramas.
    """
    ct = rng.uniform(np.cos(np.radians(theta_deg)), 1.0, n)
    ph = rng.uniform(0.0, 2.0 * np.pi, n)
    st = np.sqrt(np.maximum(0.0, 1.0 - ct**2))
    v = np.column_stack([st * np.cos(ph), st * np.sin(ph), ct])
    if eje == 0:      # z -> x
        v = np.column_stack([v[:, 2], v[:, 1], -v[:, 0]])
    elif eje == 1:    # z -> y
        v = np.column_stack([v[:, 0], v[:, 2], -v[:, 1]])
    v = v * rng.choice([-1.0, 1.0], size=n)[:, None]
    return v


def _waves_equitativo(num_waves, thetas, R, rng):
    """Reparto equitativo entre conos activos, convencion del profesor.

    Replica `src/TPMS.py:210-231` (`_sample_wave_vectors`). La rotacion R no
    existe en el original: aqui se aplica a posteriori sobre las direcciones,
    que es la extension natural y coincide con rotar los conos.
    """
    thetas = np.asarray(thetas, dtype=float)
    if np.all(thetas <= 0):
        # Isotropo: toda la esfera (TPMS.py:213-217)
        u = rng.uniform(-1.0, 1.0, num_waves)
        ph = rng.uniform(0.0, 2.0 * np.pi, num_waves)
        r = np.sqrt(np.maximum(0.0, 1.0 - u**2))
        v = np.column_stack([r * np.cos(ph), r * np.sin(ph), u])
        return v @ R.T

    activos = [(i, t) for i, t in enumerate(thetas) if t > 0]
    por_cono, resto = divmod(num_waves, len(activos))
    trozos = []
    for k, (eje, t) in enumerate(activos):
        m = por_cono + (1 if k < resto else 0)
        trozos.append(_muestrear_cono(eje, t, m, rng))
    return np.vstack(trozos) @ R.T


def wave_directions(num_waves, thetas, R=None, esquema="rechazo", rng=None):
    """Direcciones de onda unitarias, (num_waves, 3).

    esquema : 'rechazo' (GIBBON/AppFinal_V2) | 'equitativo' (profesor)
    """
    if rng is None:
        rng = np.random.default_rng()
    R = np.eye(3) if R is None else np.asarray(R, dtype=float)
    if R.shape != (3, 3):
        raise ValueError("R debe ser 3x3")
    if abs(np.linalg.det(R) - 1.0) > 1e-8 or np.linalg.norm(R.T @ R - np.eye(3)) > 1e-8:
        raise ValueError("R debe pertenecer a SO(3)")

    thetas = np.asarray(thetas, dtype=float)
    if np.any(thetas < 0) or np.any(thetas > 90):
        raise ValueError("thetas debe estar entre 0 y 90 grados")

    if esquema == "rechazo":
        if np.all(thetas <= 0):
            # GIBBON no admite thetas todos nulos (nunca aceptaria un
            # candidato). El caso isotropo se representa con thetas=90.
            thetas = np.full(3, 90.0)
        return _waves_rechazo(num_waves, thetas, R, rng)
    if esquema == "equitativo":
        return _waves_equitativo(num_waves, thetas, R, rng)
    raise ValueError("esquema debe ser 'rechazo' o 'equitativo'")


# ---------------------------------------------------------------------------
# Region degenerada del espacio de angulos
# ---------------------------------------------------------------------------

_N_FIBONACCI = 20000
_ESFERA = None


def _esfera_fibonacci():
    """Puntos casi uniformes en la esfera, fijos: sin azar, sin semilla."""
    global _ESFERA
    if _ESFERA is None:
        k = np.arange(_N_FIBONACCI) + 0.5
        z = 1.0 - 2.0 * k / _N_FIBONACCI
        r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
        phi = np.pi * (1.0 + np.sqrt(5.0)) * k
        _ESFERA = np.abs(np.column_stack([r * np.cos(phi), r * np.sin(phi), z]))
    return _ESFERA


def region_degenerada(thetas, esquema="rechazo", num_waves=None):
    """¿Dan estos thetas EXACTAMENTE la misma estructura que (90, 90, 90)?

    POR QUE PASA
    ------------
    Con muestreo por rechazo se aceptan, uniformes, los candidatos que caen en
    la UNION de los conos. Si esa union cubre toda la esfera, se aceptan todos
    los candidatos: el generador aleatorio se consume igual que en el caso
    isotropo y la realizacion es la misma bit a bit. Una region entera del
    espacio nominal es un solo punto, con gradiente nulo para cualquier
    optimizador.

    CUANDO LA UNION CUBRE LA ESFERA
    -------------------------------
    Una direccion unitaria d escapa del cono i si |d_i| <= cos(theta_i). Hay
    un conjunto de medida positiva que escapa de los tres si y solo si existe
    d con d_i^2 < cos^2(theta_i) para todo i y sum d_i^2 = 1, es decir si
    todos los cosenos son positivos y sum cos^2(theta_i) > 1. Por tanto:

        degenerado  <=>  algun theta >= 90   o   sum cos^2(theta_i) <= 1

    La condicion que se habia escrito —los tres por encima de 54.74 grados,
    el angulo de (1,1,1) a un eje— es solo el caso simetrico (3 cos^2 54.74 =
    1). La general es mas ancha: (70, 70, 30) tambien es isotropo, con
    sum cos^2 = 0.98. El bloque 17 de `tests/` lo comprueba generando.

    Con muestreo EQUITATIVO no ocurre: cada cono recibe su tercio de ondas
    concentradas junto a su eje, cubra o no la union la esfera. Tampoco con
    todos los thetas a 0 en equitativo, que ya es el caso isotropo por
    definicion y no una degeneracion.

    CERCA DEL BORDE TAMBIEN
    -----------------------
    Fuera de la region, pero cerca, la parte de la esfera que la union NO
    cubre es diminuta, y si ninguno de los candidatos cae en ella la
    realizacion vuelve a ser exactamente la isotropa. Medido: (60, 60, 40),
    con sum cos^2 = 1.087, dio a 32^3 y 300 ondas la misma mascara que
    (90, 90, 90). Por eso se da tambien `fraccion_sin_cubrir` (sobre 20 000
    puntos fijos de Fibonacci, resolucion 5e-5) y, si se pasa `num_waves`,
    `prob_realizacion_isotropa` = (1 - fraccion)^num_waves: la probabilidad de
    que ninguna onda aprovechable caiga fuera de la union. Por encima de 0.5
    el gradiente respecto a estos angulos es, en la practica, nulo.

    Devuelve datos, no frases: {"degenerado", "motivo", "suma_cos2",
    "holgura", "fraccion_sin_cubrir", "prob_realizacion_isotropa"}. `motivo`
    es "" | "theta_90" | "suma_cos2" | "todos_cero" | "equitativo". `holgura`
    = sum cos^2 - 1: negativa o cero dentro de la region.
    """
    th = np.asarray(thetas, dtype=float).ravel()
    c = np.cos(np.radians(np.clip(th, 0.0, 90.0)))
    s = float((c ** 2).sum())
    fuera = np.all(_esfera_fibonacci() <= np.where(th <= 0, 1.0, c)[None, :],
                   axis=1)
    frac = 0.0 if np.any(th >= 90.0) else float(fuera.mean())
    out = {"degenerado": False, "motivo": "", "suma_cos2": s,
           "holgura": s - 1.0, "fraccion_sin_cubrir": frac,
           "prob_realizacion_isotropa": None}
    if esquema != "rechazo":
        out["motivo"] = "equitativo"
        return out
    if num_waves is not None:
        out["prob_realizacion_isotropa"] = float((1.0 - frac) ** int(num_waves))
    if np.all(th <= 0):
        out.update(degenerado=True, motivo="todos_cero")
    elif np.any(th >= 90.0):
        out.update(degenerado=True, motivo="theta_90")
    elif s <= 1.0 + 1e-12:
        out.update(degenerado=True, motivo="suma_cos2")
    if out["degenerado"]:
        out["fraccion_sin_cubrir"] = 0.0
        if num_waves is not None:
            out["prob_realizacion_isotropa"] = 1.0
    return out


def canonicalizar_thetas(thetas, esquema="rechazo"):
    """Los thetas que representan la estructura: (90, 90, 90) si son
    degenerados, ellos mismos si no. Sirve de clave para no evaluar dos veces
    un mismo punto con nombres distintos."""
    if region_degenerada(thetas, esquema)["degenerado"]:
        return [90.0, 90.0, 90.0]
    return [float(t) for t in np.asarray(thetas, float).ravel()]


# ---------------------------------------------------------------------------
# Campo y mascara
# ---------------------------------------------------------------------------

def euler_R(rx_deg, ry_deg, rz_deg):
    """Matriz de rotacion R = Rz*Ry*Rx, la convencion de `generateSpinodoid`.

    Replica AppFinal_V2.m:1283-1289. Los deslizadores de rotacion de la app
    solo cubren +-90 grados y no pueden representar todas las orientaciones de
    un VOI; por eso existe el mecanismo `pendingR` (correccion F4), que deja a
    los optimizadores forzar una R arbitraria. Aqui R se acepta directamente.
    """
    rx, ry, rz = (np.radians(float(a)) for a in (rx_deg, ry_deg, rz_deg))
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(rx), -np.sin(rx)],
                   [0, np.sin(rx), np.cos(rx)]])
    Ry = np.array([[np.cos(ry), 0, np.sin(ry)],
                   [0, 1, 0],
                   [-np.sin(ry), 0, np.cos(ry)]])
    Rz = np.array([[np.cos(rz), -np.sin(rz), 0],
                   [np.sin(rz), np.cos(rz), 0],
                   [0, 0, 1]])
    return Rz @ Ry @ Rx


def _evaluar_campo(pts, dirs, fases, wave_number, motor, precision, bloque):
    """Suma de cosenos sobre la nube de puntos.

    La expresion  sum_i cos(beta*<n_i,x> + gamma_i)  es un PRODUCTO MATRICIAL
    disfrazado de bucle: <n_i, x_j> es la matriz (puntos x ondas). Formularlo
    como matmul lo deja en manos de BLAS en lugar de la difusion de numpy, que
    es lo que hacia la version anterior.

    Medido a 64^3 con 700 ondas (8 nucleos):
        difusion numpy f64  5.87 s   (version original)
        matmul   numpy f64  3.66 s   x1.6   dif. maxima 2e-14
        matmul   numpy f32  1.10 s   x5.3   dif. maxima 1e-5

    Sobre la precision reducida: el campo se usa para decidir GRF <= phi0, de
    modo que solo importan los voxeles a menos de ~1e-5 del umbral. En los
    casos medidos NO cambio ninguno de los 262144. Por la densidad de valores
    del GRF cerca del umbral cabe esperar del orden de 1-3 voxeles en otras
    realizaciones, es decir 1e-5 relativo en BV/TV: tres ordenes de magnitud
    por debajo del ruido estocastico del propio generador.

    Aun asi el valor por defecto es f64: la precision reducida se pide
    explicitamente, porque es un compromiso y no debe entrar en silencio.

    POR QUE NO SE USA TORCH
    -----------------------
    Una version con torch multihilo medida aparte daba 0.52 s (x10). Se
    descarto: torch y el numpy de MKL enlazan runtimes de OpenMP distintos y
    en este entorno chocan con `OMP: Error #15`. El unico rodeo disponible
    (KMP_DUPLICATE_LIB_OK) esta documentado por Intel como no soportado y
    capaz de "producir resultados incorrectos en silencio". En una herramienta
    cuyo valor es la fidelidad numerica verificada, un factor 2 adicional no
    compensa ese riesgo.
    """
    if motor not in ("numpy",):
        raise ValueError(
            f"motor '{motor}' no disponible; solo 'numpy'. Ver la nota sobre "
            "torch y el conflicto de OpenMP en la documentacion de esta funcion.")

    n_pts = pts.shape[0]
    nw = dirs.shape[0]
    dt = np.float32 if precision == "f32" else np.float64

    if bloque is None:
        # ~64 MB por bloque intermedio
        bloque = max(1024, int(8e6 / max(nw, 1)))

    D = np.ascontiguousarray((wave_number * dirs.T).astype(dt))   # (3, nw)
    F = fases.astype(dt)
    P = np.ascontiguousarray(pts.astype(dt))

    out = np.empty(n_pts, dtype=dt)
    for i in range(0, n_pts, bloque):
        A = P[i:i + bloque] @ D
        A += F
        np.cos(A, out=A)
        out[i:i + bloque] = A.sum(1)
    return out


def campo_grf(resolution, wave_number, num_waves, thetas, R=None,
              esquema="rechazo", seed=None, domain_size=1.0,
              motor="numpy", precision="f64", bloque=None):
    """Evalua el GRF sobre la rejilla del dominio.

    Replica `spinodoid.m:212-232`:
        rejilla  = linspace(0, domainSize, resolution) en los tres ejes,
                   ndgrid -> dimension 1 = X
        GRF(x)   = sum_i sqrt(2/N) * cos(waveNumber * <n_i, x> + gamma_i)

    Las fases se muestrean uniformes en [0, 2*pi). GIBBON usa `rand_angle`, que
    pese a su construccion rebuscada es uniforme en [0, 2*pi) — se comprobo
    analiticamente — de modo que `rng.uniform` es equivalente en distribucion.

    El campo se evalua por bloques de PUNTOS con un producto matricial (ver
    `_evaluar_campo`), no por bloques de ondas con difusion. Materializar
    (res^3, num_waves) de una vez ocuparia ~1.5 GB a 64^3 x 700.

    precision : 'f64' (por defecto, exacta) | 'f32' (~3x, ver _evaluar_campo)

    El orden de consumo del generador aleatorio —primero las direcciones,
    despues las fases— se mantiene identico al de la version anterior, de modo
    que una misma semilla sigue dando la misma realizacion.

    Devuelve (GRF, dirs, fases) con GRF de forma (res, res, res), dim 1 = X.
    """
    rng = np.random.default_rng(seed)
    res = int(resolution)
    dirs = wave_directions(num_waves, thetas, R=R, esquema=esquema, rng=rng)
    fases = rng.uniform(0.0, 2.0 * np.pi, num_waves)

    eje = np.linspace(0.0, float(domain_size), res)
    X, Y, Z = np.meshgrid(eje, eje, eje, indexing="ij")   # dim 1 = X (ndgrid)
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    suma = _evaluar_campo(pts, dirs, fases, wave_number, motor, precision, bloque)
    GRF = (np.sqrt(2.0 / num_waves) * suma).astype(np.float64).reshape(res, res, res)
    return GRF, dirs, fases


def generar_mascara(resolution, wave_number, num_waves, thetas, rho,
                    R=None, esquema="rechazo", seed=None, domain_size=1.0,
                    erfinv="exacta", mayor_componente=False,
                    motor="numpy", precision="f64"):
    """Mascara booleana solida (True = hueso), forma (res, res, res), dim 1 = X.

    Fase solida = {GRF <= phi0}, la misma convencion que `spinodoid.m` usa al
    llamar a `isosurface(..., levelset)` con `isocaps 'enclose','below'`.

    AVISO sobre el codigo del profesor: alli `implicit_func` devuelve
    GRF + phi0 y la fase solida sale de `pgal.Difference(cuboid, tpms)`, es
    decir de DOS inversiones de signo encadenadas. Si se toma su
    `implicit_func` y se umbraliza directamente en phi < 0 se obtiene la fase
    PORO, no la osea. Aqui la convencion es explicita y directa.

    mayor_componente : si True, conserva solo la mayor componente conexa por
        caras (conectividad 6), analogo a `localLargestComponent6`. La app
        filtra la mayor componente de SUPERFICIE via `tesgroup`, que no es
        exactamente lo mismo; por eso el valor por defecto es False y el filtro
        se declara explicitamente cuando se quiere.
    """
    GRF, dirs, fases = campo_grf(resolution, wave_number, num_waves, thetas,
                                 R=R, esquema=esquema, seed=seed,
                                 domain_size=domain_size,
                                 motor=motor, precision=precision)
    phi0 = level_set(rho, erfinv=erfinv)
    BW = GRF <= phi0

    if mayor_componente:
        from .morphometry import mayor_componente_6
        BW = mayor_componente_6(BW)

    info = {
        "levelset": phi0,
        "rho_objetivo": float(rho),
        "rho_obtenida": float(BW.mean()),
        "esquema": esquema,
        "n_ondas": int(num_waves),
        # Las dos lecturas del numero de onda, con nombre que declara la
        # unidad: `wave_number` es radianes (lo que recibio el generador) y
        # `wave_number_pi` su multiplo de pi (lo que muestra el deslizador).
        # Ver `procedencia.py` para el error que motiva guardar las dos.
        "wave_number": float(wave_number),
        "wave_number_rad": float(wave_number),
        "wave_number_pi": float(wave_number) / np.pi,
        "thetas": [float(t) for t in np.asarray(thetas).ravel()],
    }
    reg = region_degenerada(thetas, esquema, num_waves=num_waves)
    info["thetas_degenerados"] = bool(reg["degenerado"])
    info["prob_realizacion_isotropa"] = reg["prob_realizacion_isotropa"]
    return BW, GRF, info


def derivadas_grf(pts, dirs, fases, wave_number, num_waves=None, bloque=None):
    """Gradiente y hessiano ANALITICOS del GRF en una nube de puntos.

    El campo es una suma de cosenos, asi que sus derivadas son exactas y no
    hace falta diferenciarlo numericamente:

        phi(x)      = sqrt(2/N) sum_q cos(beta <n_q, x> + gamma_q)
        grad phi    = -sqrt(2/N) beta   sum_q n_q          sin(...)
        hess phi    = -sqrt(2/N) beta^2 sum_q n_q (x) n_q  cos(...)

    PARA QUE SIRVE. La curvatura de la isosuperficie sale de estas dos
    cantidades (`curvatura.curvaturas_implicitas`), de modo que un spinodoide
    tiene curvatura EXACTA, sin error de mallado ni de diferencias finitas. Es
    la respuesta cerrada contra la que se valida el estimador discreto que
    luego hay que aplicar al hueso real, donde no hay campo analitico.

    `dirs` y `fases` son los que devuelve `campo_grf`; la misma semilla da las
    mismas y por tanto el mismo campo. `num_waves` se toma de `dirs` si no se
    da: el factor sqrt(2/N) depende de N y ponerlo a mano invita a que no
    coincida con el campo del que salieron las direcciones.

    Devuelve (grad, hess) con formas (n, 3) y (n, 3, 3).
    """
    P = np.ascontiguousarray(np.asarray(pts, dtype=np.float64))
    D = np.asarray(dirs, dtype=np.float64)
    F = np.asarray(fases, dtype=np.float64)
    N = int(D.shape[0] if num_waves is None else num_waves)
    b = float(wave_number)
    amp = np.sqrt(2.0 / N)

    n_pts = P.shape[0]
    if bloque is None:
        bloque = max(256, int(2e6 / max(D.shape[0], 1)))

    grad = np.empty((n_pts, 3))
    hess = np.empty((n_pts, 3, 3))
    Dt = np.ascontiguousarray((b * D).T)                        # (3, nw)
    for i in range(0, n_pts, bloque):
        A = P[i:i + bloque] @ Dt
        A += F
        s = np.sin(A)
        c = np.cos(A)
        grad[i:i + bloque] = -amp * b * (s @ D)
        # sum_q n_q (x) n_q cos(...): se contrae en un solo einsum por bloque.
        hess[i:i + bloque] = -amp * b * b * np.einsum(
            "pq,qi,qj->pij", c, D, D, optimize=True)
    return grad, hess
