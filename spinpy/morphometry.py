"""
morphometry.py — Port de la morfometria de AppFinal_V2.m.

Traduccion linea a linea de:
    localMorphometryFromVoxels   (AppFinal_V2.m:5890-5996)
    localMILTensor               (AppFinal_V2.m:5733-5866)
    localLoadBearingFraction     (AppFinal_V2.m:6012-6030)
    localLargestComponent6       (AppFinal_V2.m:6502-6572)

Se respetan las invariantes documentadas del proyecto, que son correcciones ya
pagadas en el codigo MATLAB y que un port ingenuo reintroduce en silencio:

  * Tb.Th = 2*BV/BS  (modelo de placas de Parfitt, como CTAn/Scanco).
    NO 4*BV/BS, que es el modelo de barras y duplica el valor.
  * DA sale del TENSOR MIL (Harrigan & Mann 1984), nunca de la covarianza de
    la nube de puntos. El MIL tiene un suelo de ruido de DA ~ 1.07: por debajo
    de esa diferencia no hay significado.
  * Las estructuras laminares degeneran el ajuste del elipsoide; se acota el
    autovalor y se marca el resultado como COTA INFERIOR (`MIL_clamped`).
  * La superficie se mide con marching cubes sobre la mascara binaria, igual
    que el VOI real, para que el sesgo del metodo se cancele al restar.

CONVENCION DE EJES — el punto donde mas facil es equivocarse
------------------------------------------------------------
MATLAB llama `isosurface(permute(BW,[2 1 3]), 0.5)` porque `isosurface`
devuelve los vertices en convencion meshgrid (columna, fila, pagina); la
permutacion previa deja V(:,1) = X. Es el bug B1/B2 del historial del proyecto.

En Python NO hace falta permutar: `skimage.measure.marching_cubes` devuelve los
vertices en orden de indice del array, asi que con BW de forma (nx, ny, nz) los
vertices ya salen como (X, Y, Z). Se usa `method='lorensen'` (marching cubes
clasico) por ser el mas proximo al de MATLAB; el 'lewiner' por defecto resuelve
las ambiguedades de otra manera y da areas ligeramente distintas.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from skimage import measure

# Conectividad 6 (solo caras), la misma que usa el mallado hexaedrico de la
# homogeneizacion: dos voxeles que se tocan por arista o vertice no comparten
# nodos y no se transmiten carga.
_CONN6 = ndimage.generate_binary_structure(3, 1)


# ---------------------------------------------------------------------------
# Utilidades de componentes conexas
# ---------------------------------------------------------------------------

def mayor_componente_6(BW):
    """Mayor componente conexa por caras. Port de localLargestComponent6."""
    BW = np.asarray(BW, dtype=bool)
    lab, n = ndimage.label(BW, structure=_CONN6)
    if n <= 1:
        return BW
    tam = np.bincount(lab.ravel())
    tam[0] = 0
    return lab == int(tam.argmax())


def fraccion_portante(BW):
    """Fraccion de material en componentes que unen las dos caras en Z.

    Port de localLoadBearingFraction (K3). Vale 1 cuando todo el material forma
    parte de caminos que atraviesan la probeta y baja al aparecer islas y
    munones ciegos. Se mide, no se optimiza (peso 0 en el error).
    """
    BW = np.asarray(BW, dtype=bool)
    n_tot = int(BW.sum())
    if n_tot == 0:
        return 0.0
    lab, n = ndimage.label(BW, structure=_CONN6)
    if n == 0:
        return 0.0
    ini = set(np.unique(lab[:, :, 0])) - {0}
    fin = set(np.unique(lab[:, :, -1])) - {0}
    pasantes = ini & fin
    if not pasantes:
        return 0.0
    return float(np.isin(lab, list(pasantes)).sum() / n_tot)


# ---------------------------------------------------------------------------
# Superficie
# ---------------------------------------------------------------------------

def area_superficie(BW, spacing, metodo="lorensen"):
    """Area de la isosuperficie 0.5 de la mascara, en unidades fisicas.

    Equivale a las lineas 5942-5950 de AppFinal_V2.m. Devuelve 0.0 si no hay
    interfaz (mascara vacia o llena).
    """
    BW = np.asarray(BW, dtype=bool)
    if not BW.any() or BW.all():
        return 0.0
    verts, faces, _, _ = measure.marching_cubes(
        BW.astype(np.float32), level=0.5,
        spacing=tuple(float(s) for s in spacing), method=metodo)
    a = verts[faces[:, 0]]
    b = verts[faces[:, 1]]
    c = verts[faces[:, 2]]
    return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum())


def area_superficie_partes(BW, spacing, metodo="lorensen"):
    """Separa el area en INTERNA (interfaz hueso-vacio) y EXTERNA (las tapas).

    POR QUE HACEN FALTA LAS DOS, Y POR QUE NO SON LA MISMA COSA
    -----------------------------------------------------------
    La superficie que se mide sobre un VOI tiene dos partes de naturaleza
    distinta y no se pueden sumar sin pensar:

      INTERNA   la interfaz real hueso-vacio, que es donde vive el hueso: donde
                se remodela, donde se adhieren las celulas, donde ocurre todo.
                Es la que entra en BS/BV y en Tb.Th, y la unica que tiene
                sentido biologico. Es lo que este modulo llamaba `BS` a secas,
                y lo sigue llamando: no se cambia el nombre de una metrica que
                ya esta en informes y en el error del ajuste.

      EXTERNA   las seis caras por donde el cubo CORTA las trabeculas. No es
                interfaz de nada: es la huella del bisturi. Si se sumara a la
                interna, BS/BV subiria al reducir el VOI —mas corte por unidad
                de volumen— y Tb.Th bajaria, sin que la estructura hubiera
                cambiado en absoluto.

    Esa es exactamente la invariante del proyecto ("el area excluye las seis
    tapas"). Lo que anade esta funcion es poder MEDIR la que se excluye, en vez
    de solo descartarla: su tamano dice cuanto del VOI es borde, que es un
    diagnostico de si la ventana es bastante grande para lo que se mide.

    COMO SE SEPARAN, CON UNA SOLA DEFINICION DE AREA
    ------------------------------------------------
    Tentador y equivocado: contar voxeles de hueso en las seis caras y
    multiplicar por el area del voxel. Eso mezcla dos definiciones de area —una
    de marching cubes y otra de conteo de cuadraditos— que difieren en torno al
    8.5 % (la correccion C4 de este proyecto), y la resta de dos definiciones
    distintas no significa nada.

    Aqui se usa la MISMA marching cubes dos veces: sobre la mascara tal cual,
    que deja la superficie abierta en el borde y da la interna; y sobre la
    mascara rodeada de una capa de vacio, que obliga al algoritmo a cerrar las
    tapas y da el total. La externa es la diferencia, y las dos salen de la
    misma regla.

    Devuelve (interna, externa, total), en unidades de `spacing` al cuadrado.
    """
    BW = np.asarray(BW, dtype=bool)
    interna = area_superficie(BW, spacing, metodo=metodo)
    # Una sola capa basta: marching cubes solo necesita ver vacio al otro lado.
    cerrada = np.pad(BW, 1, mode="constant", constant_values=False)
    total = area_superficie(cerrada, spacing, metodo=metodo)
    # El max() no es cosmetico: con una mascara vacia o llena las dos llamadas
    # devuelven 0 por caminos distintos, y un -1e-15 se propagaria a un BS/BV
    # negativo.
    return interna, max(total - interna, 0.0), total


def tamano_poro(BW, spacing, n_radios=None, campo=False):
    """Tamano de poro Po.Dm: espesor local de la fase PORO, no de la osea.

    Es el metodo de Hildebrand y Ruegsegger aplicado al complemento de la
    mascara, que es lo que usan CTAn y Scanco para Tb.Sp y lo que la literatura
    de andamios llama distribucion de tamano de poro.

    SOBRE "HACERLO RAPIDO CON LA TRANSFORMADA DE DISTANCIA"
    -------------------------------------------------------
    La transformada de distancia ya esta dentro: es el primer paso del metodo.
    Lo que NO vale es quedarse ahi. La distancia al hueso mas proximo, d(x), no
    es el tamano del poro en x: vale el radio del poro solo en el eje medial y
    cae a cero al acercarse a la pared. Un poro cilindrico de diametro D daria
    valores entre 0 y D/2 segun donde se mire, y su media saldria del orden de
    D/3 en vez de D. El paso que falta -y que cuesta las iteraciones de
    `espesor.espesor_local`- es quedarse con la mayor esfera que CONTIENE el
    punto, no con la centrada en el. Ver la cabecera de `espesor.py`, donde
    esto esta medido.

    QUE SE DEVUELVE, Y POR QUE UNA DISTRIBUCION Y NO UN NUMERO
    -----------------------------------------------------------
    Tb.Sp, tal como lo calcula `morfometria`, sale del modelo de placas de
    Parfitt: Tb.Th * (1/BV.TV - 1). Es una formula, no una medida: supone que
    la estructura son placas paralelas. Po.Dm se MIDE sobre la geometria, y en
    hueso trabecular los dos no tienen por que coincidir; la diferencia entre
    ambos es en si misma un indicador de cuanto se aparta la estructura del
    modelo de placas.

    Devuelve un diccionario con Po.Dm media, mediana, desviacion, percentiles y
    el numero de voxeles de poro medidos. Con `campo=True` en la llamada a
    `espesor.espesor_local` se puede recuperar el campo entero si hace falta
    pintarlo; aqui solo se resumen las cifras.

    EL AVISO DE BORDE, Y POR QUE NO ES EL OBVIO
    --------------------------------------------
    La transformada de distancia mide la distancia al HUESO, y fuera del array
    no hay hueso: por tanto trata el exterior como si el poro siguiera. Para un
    poro que de verdad continua fuera del VOI eso es lo correcto; para uno que
    se cerraba justo despues, Po.Dm sale de mas.

    La primera version de este aviso contaba los voxeles de poro cuya COMPONENTE
    toca una cara del cubo. No sirve: en cualquier estructura de celda abierta
    -o sea, en todo hueso trabecular- el poro es una sola componente que toca
    las seis caras, y el aviso salia 100 % siempre. Medido sobre el VOI
    proximal de H4: 100 %, sin decir nada de nadie.

    Lo que se cuenta ahora es la fraccion de poro cuya ESFERA INSCRITA se sale
    de la ventana, es decir los puntos cuya distancia al borde es menor que su
    radio: son exactamente aquellos cuyo espesor se apoya en la extrapolacion.
    Y ademas se dan las cifras restringidas a los que NO se salen
    (`PoDm_interior`): si las dos versiones se parecen, la ventana no manda; si
    difieren, el VOI es pequeno para los poros que contiene y hay que decirlo.

    Sigue siendo un AVISO y no una correccion: no hay forma honesta de saber
    cuanto seguia el poro fuera del VOI.
    """
    from .espesor import espesor_local, estadisticas

    BW = np.asarray(BW, dtype=bool)
    poro = ~BW
    if not poro.any():
        return {"PoDm": np.nan, "PoDm_mediana": np.nan, "PoDm_sd": np.nan,
                "PoDm_p05": np.nan, "PoDm_p95": np.nan, "PoDm_n": 0,
                "PoDm_frac_ventana": np.nan, "PoDm_interior": np.nan,
                "PoDm_interior_mediana": np.nan, "PoDm_interior_n": 0}

    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)

    esp = espesor_local(poro, spacing, n_radios=n_radios)
    st = estadisticas(esp, poro)

    # Radio de la esfera inscrita en cada punto y distancia a la cara mas
    # proxima del cubo. Donde el radio es mayor, la esfera se sale y el
    # espesor medido se apoya en suponer que el poro continua fuera.
    r = ndimage.distance_transform_edt(poro, sampling=spacing)
    # Se pliega eje a eje y no con `np.stack` + `minimum.reduce`: los tres
    # vectores tienen longitudes distintas en cuanto el VOI no es cubico, y
    # apilarlos revienta.
    ejes = [np.minimum(np.arange(n), n - 1 - np.arange(n)) * sp
            for n, sp in zip(poro.shape, spacing)]
    d_borde = np.minimum(np.minimum(ejes[0][:, None, None],
                                    ejes[1][None, :, None]),
                         ejes[2][None, None, :])
    sale = poro & (d_borde < r)
    dentro = poro & ~sale

    n_poro = int(poro.sum())
    st_in = estadisticas(esp, dentro) if dentro.any() else {}

    out = {
        "PoDm": st.get("media", np.nan),
        "PoDm_mediana": st.get("mediana", np.nan),
        "PoDm_sd": st.get("sd", np.nan),
        "PoDm_p05": st.get("p05", np.nan),
        "PoDm_p95": st.get("p95", np.nan),
        "PoDm_n": n_poro,
        # Aviso: fraccion de poro cuyo espesor se apoya en la extrapolacion
        # fuera de la ventana.
        "PoDm_frac_ventana": float(int(sale.sum()) / max(n_poro, 1)),
        # Las mismas cifras sobre el poro que NO se sale. Si se parecen a las
        # de arriba, la ventana no manda.
        "PoDm_interior": st_in.get("media", np.nan),
        "PoDm_interior_mediana": st_in.get("mediana", np.nan),
        "PoDm_interior_n": int(dentro.sum()),
    }
    if campo:
        # El campo entero, para quien necesite la DISTRIBUCION y no solo sus
        # cifras (`metricas_forma`). Cero en el hueso.
        out["PoDm_campo"] = esp
    return out


# ---------------------------------------------------------------------------
# Metricas de forma para la funcion objetivo: distribuciones, no solo medias
# ---------------------------------------------------------------------------
#
# QUE SON Y POR QUE VIVEN AQUI
# ----------------------------
# `Estudio_Discriminadores` midio que las metricas que el ajuste iguala (BV/TV,
# Tb.Th, DA) no separan el VOI del spinodoide ajustado, y que las que si lo
# hacen son de FORMA: dispersion de la curvatura media (64 sigma), tamano de
# poro (39), coeficiente de variacion del espesor (23), fraccion de placas por
# Ellipsoid Factor (36). Para que puedan entrar en el objetivo tienen que
# medirse por el MISMO camino en el VOI y en el candidato (correccion C4), y el
# sitio de eso es este modulo.
#
# Ademas de las cifras se devuelven los CUANTILES de cada distribucion
# (`cuantiles`), que es lo que necesita una distancia entre distribuciones
# (Wasserstein, Hellinger) al estilo de Xiao et al. (2025). Se guardan los
# cuantiles y no las muestras: 101 numeros por distribucion, que caben en un
# JSON, y con ellos la distancia de Wasserstein-1 es exacta salvo la regla del
# punto medio.

# Suavizado para la curvatura de una mascara binaria, en voxeles. Es el de
# `Estudio_Discriminadores` —el minimo que quita la escalera de marching cubes
# sin comerse una trabecula de 3.3 voxeles— y se aplica igual a los dos lados.
SIGMA_CURVATURA = 1.0

# Niveles de probabilidad guardados por distribucion.
N_CUANTILES = 101


def niveles_cuantiles(n=N_CUANTILES):
    """Niveles (k + 0.5) / n: posiciones de Hazen, la convencion de `prctile`
    de MATLAB que usa todo el proyecto."""
    return (np.arange(int(n)) + 0.5) / int(n)


def cuantiles(valores, pesos=None, n=N_CUANTILES):
    """Cuantiles de una distribucion en los niveles de `niveles_cuantiles`.

    Con `pesos` (por ejemplo el area de cada vertice en una curvatura) cada
    valor cuenta en proporcion a su peso, con la misma posicion de Hazen: el
    centro de su tramo de peso acumulado. Devuelve un array, o None si no hay
    ningun valor utilizable.
    """
    v = np.asarray(valores, float).ravel()
    ok = np.isfinite(v)
    if pesos is not None:
        w = np.asarray(pesos, float).ravel()
        ok &= np.isfinite(w) & (w > 0)
    if not ok.any():
        return None
    p = niveles_cuantiles(n)
    if pesos is None:
        return np.percentile(v[ok], 100.0 * p, method="hazen")
    v, w = v[ok], w[ok]
    o = np.argsort(v, kind="stable")
    v, w = v[o], w[o]
    c = (np.cumsum(w) - 0.5 * w) / w.sum()
    return np.interp(p, c, v)


def curvatura_de_mascara(BW, spacing, sigma=SIGMA_CURVATURA, margen=2):
    """Curvaturas principales de una mascara BINARIA, por el camino del VOI.

    Binarizar -> suavizar con una gaussiana de `sigma` voxeles -> mallar el
    nivel 0.5 -> Rusinkiewicz con normales tomadas del campo suavizado. Es el
    unico camino disponible para un VOI de micro-CT y por eso se usa tambien
    en el candidato, aunque del spinodoide se conozca el campo exacto: medir
    uno por el campo y otro por la mascara compararia dos definiciones.

    Los vertices a menos de `margen` voxeles de las caras se descartan (area
    cero): alli la malla esta cortada y la curvatura es del recorte.

    Devuelve (k1, k2, areas), en 1/mm si `spacing` va en mm.
    """
    from skimage import measure
    from .curvatura import (areas_baricentricas, curvaturas_malla,
                            normales_desde_campo)

    BW = np.asarray(BW, dtype=bool)
    sp = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if sp.size == 1:
        sp = np.repeat(sp, 3)
    campo = ndimage.gaussian_filter(BW.astype(np.float32), float(sigma))
    verts, faces, _, _ = measure.marching_cubes(campo, level=0.5,
                                                spacing=tuple(sp))
    # El campo suavizado CRECE hacia el hueso; la normal tiene que apuntar al
    # VACIO (convencion de `curvatura.py`), de ahi el signo.
    nv = normales_desde_campo(verts, -campo, sp)
    k1, k2, _ = curvaturas_malla(verts, faces, normales=nv)
    areas = areas_baricentricas(verts, faces)
    if margen > 0:
        lim = (np.array(campo.shape) - 1.0) * sp
        d = np.minimum(verts, lim[None, :] - verts)
        areas = np.where((d < margen * sp[None, :]).any(axis=1), 0.0, areas)
    return k1, k2, areas


def metricas_forma(BW, spacing, tbth=False, poro=False, curvatura=False,
                   ef=False, n_cuantiles=N_CUANTILES):
    """Metricas de forma que puede usar la funcion objetivo, a demanda.

    Cada una tiene su compuerta porque su coste es muy distinto (medido):

      tbth       espesor local del solido: TbTh_local, TbTh_CV y cuantiles.
                 Una transformada de distancia por radio.
      poro       tamano de poro Po.Dm (`tamano_poro`) y sus cuantiles.
                 0.7 s a 64^3, 12 s a 128^3.
      curvatura  H_desv, H_medio, fracciones concava/silla y cuantiles de H,
                 ponderados por area. Una gaussiana y una marching cubes.
      ef         Ellipsoid Factor (`elipsoide.factor_elipsoide`). 373 s a
                 48^3: NUNCA dentro de un bucle de ajuste; el ajuste lo deja
                 para los finalistas.

    Las cifras de espesor local NO son citables como Tb.Th (sesgo de ~30 % a
    3-4 voxeles, ver `espesor.py`); sirven para comparar VOI y candidato
    medidos con el mismo voxel, que es lo que hace el ajuste.
    """
    BW = np.asarray(BW, dtype=bool)
    sp = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if sp.size == 1:
        sp = np.repeat(sp, 3)
    out = {"cuantiles": {}}
    trivial = (not BW.any()) or BW.all()

    if tbth:
        out.update(TbTh_local=np.nan, TbTh_local_sd=np.nan, TbTh_CV=np.nan)
        if not trivial:
            from .espesor import espesor_local
            v = espesor_local(BW, sp)[BW].astype(float)
            v = v[np.isfinite(v) & (v > 0)]
            if v.size > 1:
                media, sd = float(v.mean()), float(v.std(ddof=1))
                out.update(TbTh_local=media, TbTh_local_sd=sd,
                           TbTh_CV=sd / media if media > 0 else np.nan)
                out["cuantiles"]["TbTh"] = cuantiles(v, n=n_cuantiles).tolist()

    if poro:
        tp = tamano_poro(BW, sp, campo=not trivial)
        campo_poro = tp.pop("PoDm_campo", None)
        out.update(tp)
        if campo_poro is not None:
            v = np.asarray(campo_poro, float)[~BW]
            v = v[np.isfinite(v) & (v > 0)]
            if v.size:
                out["cuantiles"]["PoDm"] = cuantiles(v, n=n_cuantiles).tolist()

    if curvatura:
        out.update(H_desv=np.nan, H_medio=np.nan, curv_concava=np.nan,
                   curv_silla=np.nan)
        if not trivial:
            from .curvatura import resumen
            k1, k2, a = curvatura_de_mascara(BW, sp)
            r = resumen(k1, k2, a)
            if r.get("area", 0.0) > 0:
                out.update(H_desv=r["H_desv"], H_medio=r["H_medio"],
                           curv_concava=r["concava"], curv_silla=r["silla"])
                q = cuantiles((np.asarray(k1) + np.asarray(k2)) / 2.0, a,
                              n=n_cuantiles)
                if q is not None:
                    out["cuantiles"]["H"] = q.tolist()

    if ef:
        from .elipsoide import factor_elipsoide
        out.update(factor_elipsoide(BW, sp))

    return out


# ---------------------------------------------------------------------------
# Tensor MIL
# ---------------------------------------------------------------------------

def tensor_mil(BW, spacing, n_dirs=61, n_lines=32, min_chord_frac=0.35,
               muestreo="voxel", sigma_vox=0.7):
    """Mean Intercept Length y tensor de fabrica. Port de localMILTensor.

    MUESTREO (correccion M1, OPCIONAL; por omision no cambia nada)
      "voxel"      lo de siempre: cada muestra de la recta lee el voxel que la
                   contiene. Identico bit a bit a localMILTensor de MATLAB.
      "suavizado"  la mascara se suaviza con una gaussiana de `sigma_vox`
                   voxeles y cada muestra lee ese campo por interpolacion
                   trilineal, con umbral 0.5.

    Por que existe (Estudio_MIL): sobre voxeles, una recta que corre OBLICUA a
    lo largo de un puntal cruza una y otra vez la escalera de su superficie y
    cuenta intersecciones que no existen. El MIL cae sobre todo a lo largo del
    eje de la estructura, el DA baja y la direccion principal se tuerce hacia
    el eje de rejilla mas proximo. Medido en un dual-lattice girado de forma
    exacta (48^3, estiramiento 2.5, cinco semillas): con "voxel" el error de
    direccion llega a 14.5 +- 1.9 grados y el DA va de 1.48 (alineado) a 1.13
    segun el giro; con "suavizado" el error no pasa de 1.6 +- 0.9 grados y el
    DA se queda en 1.55-1.64. El DA de la estructura continua es ~1.93: NINGUN
    muestreo de voxeles lo recupera a 48^3, asi que el DA absoluto sigue sin
    ser comparable entre resoluciones; lo que arregla M1 es la DEPENDENCIA DE
    LA ORIENTACION. VOI y candidato deben medirse con el mismo muestreo (C4).

    Direcciones cuasi-uniformes en el hemisferio por espiral de Fibonacci; para
    cada direccion se lanza una rejilla de n_lines^2 rectas paralelas, se
    acumula longitud osea y numero de transicciones solido<->poro, y

        MIL(d) = longitud_osea_total / (0.5 * n_transiciones)

    El factor 0.5 es la correccion F5: `totalInt` cuenta la entrada Y la salida
    de cada cuerda, es decir 2 por cuerda, y el MIL estandar (CTAn/Scanco)
    divide por el numero de cuerdas. Sin el factor se devuelve la mitad del
    valor. El DA no se ve afectado por ser un cociente, pero el MIL expuesto si.

    Despues se ajusta por minimos cuadrados el elipsoide 1/MIL^2 = n' M n y se
    diagonaliza. El autovalor MINIMO corresponde al MIL MAXIMO, es decir a la
    direccion principal de la estructura.

    Devuelve un dict con las mismas claves que el struct `fab` de MATLAB.
    """
    BW = np.asarray(BW, dtype=bool)
    spacing = np.asarray(spacing, dtype=float).ravel()
    sz = np.array(BW.shape, dtype=int)

    fab = {"DA": np.nan, "MIL": None, "dirs": None,
           "eigenvalues": np.array([1.0, 1.0, 1.0]),
           "eigenvectors": np.eye(3),
           "dir_principal": np.array([0.0, 0.0, 1.0]),
           "valid": False, "clamped": False}

    if muestreo not in ("voxel", "suavizado"):
        raise ValueError(f"muestreo debe ser 'voxel' o 'suavizado'; se recibio {muestreo!r}")
    fab["muestreo"] = muestreo

    n_bone = int(BW.sum())
    if BW.ndim < 3 or n_bone == 0 or n_bone == BW.size:
        return fab   # solido homogeneo o vacio: la anisotropia no esta definida

    campo_suave = None
    if muestreo == "suavizado":
        # sigma en VOXELES de cada eje: el artefacto es de la rejilla, no del mm
        campo_suave = ndimage.gaussian_filter(BW.astype(np.float64), float(sigma_vox))

    L = sz * spacing
    ctr = L / 2.0
    step = 0.5 * spacing.min()          # 2 muestras por voxel
    R = 0.5 * float(np.linalg.norm(L))
    min_chord = min_chord_frac * float(L.min())

    # --- Direcciones: espiral de Fibonacci sobre el hemisferio --------------
    k = np.arange(n_dirs) + 0.5
    cosT = k / n_dirs
    sinT = np.sqrt(np.maximum(0.0, 1.0 - cosT**2))
    phi = np.pi * (1.0 + np.sqrt(5.0)) * k
    dirs = np.column_stack([sinT * np.cos(phi), sinT * np.sin(phi), cosT])
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)

    # tt = -R:step:R  (colon de MATLAB: valores <= R)
    n_t = int(np.floor((2.0 * R) / step + 1e-9)) + 1
    tt = -R + step * np.arange(n_t)

    off = np.linspace(-R, R, n_lines)
    OU, OV = np.meshgrid(off, off, indexing="ij")
    OU = OU.ravel()
    OV = OV.ravel()

    MIL = np.full(n_dirs, np.nan)

    for q in range(n_dirs):
        d = dirs[q]
        tmp = np.array([0.0, 0.0, 1.0]) if abs(d[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        u = np.cross(d, tmp); u /= np.linalg.norm(u)
        v = np.cross(d, u);   v /= np.linalg.norm(v)

        origins = ctr[None, :] + OU[:, None] * u[None, :] + OV[:, None] * v[None, :]

        # P: (n_lineas, n_t, 3)
        P = origins[:, None, :] + tt[None, :, None] * d[None, None, :]
        idx = np.floor(P / spacing[None, None, :]).astype(np.int64)

        dentro = np.ones(idx.shape[:2], dtype=bool)
        for ax in range(3):
            dentro &= (idx[:, :, ax] >= 0) & (idx[:, :, ax] < sz[ax])

        # Cuerdas demasiado cortas: muy sesgadas por el borde, se descartan
        largo_ok = dentro.sum(axis=1) * step >= min_chord
        if not largo_ok.any():
            continue

        vals = np.zeros(dentro.shape, dtype=bool)
        sel = dentro & largo_ok[:, None]
        if sel.any():
            if campo_suave is None:
                ii = idx[:, :, 0][sel]
                jj = idx[:, :, 1][sel]
                kk = idx[:, :, 2][sel]
                vals[sel] = BW[ii, jj, kk]
            else:
                # centro del voxel i en (i + 0.5) * spacing
                c = (P[sel] / spacing[None, :] - 0.5).T
                vals[sel] = ndimage.map_coordinates(
                    campo_suave, c, order=1, mode="nearest") >= 0.5

        # Primer y ultimo indice dentro del volumen, por linea
        first = np.argmax(dentro, axis=1)
        last = dentro.shape[1] - 1 - np.argmax(dentro[:, ::-1], axis=1)
        util = largo_ok & (last - first + 1 >= 2) & dentro.any(axis=1)
        if not util.any():
            continue

        cols = np.arange(dentro.shape[1])[None, :]
        seg = (cols >= first[:, None]) & (cols <= last[:, None]) & util[:, None]

        total_bone = float((vals & seg).sum()) * step
        cambio = (vals[:, :-1] != vals[:, 1:]) & seg[:, :-1] & seg[:, 1:]
        total_int = int(cambio.sum())

        if total_int > 0:
            MIL[q] = total_bone / (0.5 * total_int)

    good = np.isfinite(MIL) & (MIL > 0)
    if int(good.sum()) < 9:
        return fab   # datos insuficientes para ajustar el elipsoide

    n = dirs[good]
    y = 1.0 / (MIL[good] ** 2)
    A = np.column_stack([n[:, 0]**2, n[:, 1]**2, n[:, 2]**2,
                         2 * n[:, 1] * n[:, 2],
                         2 * n[:, 0] * n[:, 2],
                         2 * n[:, 0] * n[:, 1]])
    m, *_ = np.linalg.lstsq(A, y, rcond=None)

    M = np.array([[m[0], m[5], m[4]],
                  [m[5], m[1], m[3]],
                  [m[4], m[3], m[2]]])
    M = (M + M.T) / 2.0

    lam, Vm = np.linalg.eigh(M)
    if not np.all(np.isfinite(lam)):
        return fab

    # Laminares: MIL -> inf en el plano de las placas, 1/MIL^2 -> 0, y el
    # ajuste puede dar un autovalor negativo. Se acota y se marca: el DA
    # resultante es una COTA INFERIOR.
    if np.any(lam <= 0):
        y_pos = y[y > 0]
        if y_pos.size == 0:
            return fab
        lam = np.maximum(lam, 0.5 * float(y_pos.min()))
        fab["clamped"] = True

    orden = np.argsort(lam)
    lam_s = lam[orden]
    Vs = Vm[:, orden]

    fab["eigenvalues"] = lam_s
    fab["eigenvectors"] = Vs
    fab["dir_principal"] = Vs[:, 0] / np.linalg.norm(Vs[:, 0])
    fab["DA"] = float(np.sqrt(lam_s[2] / lam_s[0]))
    fab["MIL"] = MIL
    fab["dirs"] = dirs
    fab["valid"] = True
    return fab


# ---------------------------------------------------------------------------
# Motor unico de morfometria
# ---------------------------------------------------------------------------

def morfometria_malla(BW, spacing, suavizado=20, banda=0.05, do_mil=True):
    """Morfometria sobre la superficie SUAVIZADA Y REPARADA, no sobre voxeles.

    Segundo modo de medida, que convive con `morfometria` y no la sustituye.
    Construye la isosuperficie cerrada (`solido.superficie_cerrada`), la suaviza
    con Taubin, la repara con PyMeshFix y la recorta a la caja del VOI. Sobre
    esa malla mide BV (volumen encerrado) y BS (area SIN las seis tapas del
    cubo), y de ahi las metricas de Parfitt. DA, DA2 y la fraccion portante
    siguen saliendo de los voxeles: son topologicas o del tensor MIL, y no
    tienen version de malla aqui (esa seria otra funcion).

    POR QUE EXISTE — MEDIDO SOBRE UNA ESFERA DE AREA CONOCIDA
      marching cubes crudo (lo que mide `morfometria`)   area +8.5 %
      cerrada + PyMeshFix sin suavizar                    area +8.5 %
      cerrada + Taubin(20) + PyMeshFix                    area +0.6..0.9 %
                                                          volumen +0.1..0.4 %
    PyMeshFix por si solo no quita ningun sesgo: solo cierra y repara. Lo que
    elimina la escalera de los voxeles es el suavizado de Taubin, que a
    diferencia del laplaciano conserva el volumen. La pareja es lo que da una
    superficie "de hueso" y no "de cubos".

    LO QUE CUESTA, Y SE DEVUELVE PARA QUE SE VEA
      * En trabeculas de uno o dos voxeles PyMeshFix se come material al
        reparar: un spinodoide a 48^3 pierde ~1 % frente a la superficie
        cerrada cruda. `perdida_vs_mc_pct` lo mide en cada llamada.
      * El volumen de la superficie cerrada cruda ya es ~3 % menor que el
        conteo de voxeles en estructuras finas: marching cubes coloca la
        superficie a mitad de camino entre centros solidos y vacios, y el
        conteo cuenta cubos enteros. Eso NO es una perdida, es una diferencia
        de definicion, y sobre una esfera el volumen de marching cubes es el
        mas proximo al exacto. `dif_BV_vs_voxel_pct` la reporta.
    LA CAJA DE MEDIDA ES LA DE CENTROS DE VOXEL, MEDIO VOXEL POR DENTRO
      El suavizado REDONDEA las aristas donde una trabecula cortada toca la
      cara del cubo: la tapa retrocede hacia dentro, deja de ser plana y la
      punta redondeada se cuenta como interfaz. Medido en una losa, eso
      inflaba BS un 27 %. Devolver los vertices al plano lo arregla en la
      losa pero en un spinodoide colapsa triangulos en las puntas y deja la
      malla no-manifold. La solucion robusta es recortar la malla suavizada
      MEDIO VOXEL POR DENTRO, en los planos de los centros de voxel: el
      recortador de superficies cerradas elimina la zona redondeada y fabrica
      el mismo las tapas planas, sin tocar vertices. Medido: losa +0.9 % de
      area y +0.1 % de volumen; spinodoide BS_malla/BS_voxel = 0.91-0.92, que
      es justo el sesgo de voxel que el suavizado quita.
      Consecuencia: `TV` aqui es ((n-1)*sp)^3, la caja de centros, que es
      tambien la de la superficie abierta de skimage que mide `morfometria`.
      Se devuelve `TV_voxel` para la referencia; BV/TV es comparable entre
      modos porque los dos son densidades sobre su propia caja.

    LA REGLA DE ORO SE MANTIENE (correccion C4): VOI y candidato se miden con
    el MISMO modo. Mezclar una tabla de voxeles con una de malla es mezclar
    dos definiciones de BS que difieren un 9 %.

    Devuelve las mismas claves que `morfometria` mas: modo="malla", BS_tapas,
    V_mc, V_rep, TV_voxel, BV_voxel, BS_voxel, dif_BV_vs_voxel_pct,
    perdida_vs_mc_pct, n_tri, estanca, bordes_abiertos (de la malla reparada),
    bordes_abiertos_recorte, reparado_tras_recorte, suavizado, banda.
    """
    import pyvista as pv                       # noqa: F401  (usado por solido)
    try:
        import pymeshfix
    except ImportError as e:
        raise ImportError("El modo de medida sobre malla necesita pymeshfix: "
                          "pip install pymeshfix  (extra 'malla')") from e
    from .solido import superficie_cerrada

    BW = np.asarray(BW, dtype=bool)
    spacing = np.atleast_1d(np.asarray(spacing, dtype=float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    n = np.array(BW.shape, dtype=float)
    TV_voxel = float(np.prod(n * spacing))
    TV = float(np.prod((n - 1.0) * spacing))      # caja de centros de voxel

    # La parte topologica y el MIL, de los voxeles, y de paso los valores de
    # voxel de BV y BS para el diagnostico.
    mv = morfometria(BW, spacing, do_mil=do_mil)
    m = {"modo": "malla", "TV": TV, "TV_voxel": TV_voxel,
         "suavizado": int(suavizado), "banda": float(banda),
         "DA": mv["DA"], "DA2": mv["DA2"], "MIL_valid": mv["MIL_valid"],
         "MIL_clamped": mv["MIL_clamped"], "FracPort": mv["FracPort"],
         "BV_voxel": mv["BV"], "BS_voxel": mv["BS"]}
    nan = {"BV": np.nan, "BS": np.nan, "BVTV": np.nan, "PoTot": np.nan,
           "BSBV": np.nan, "BSTV": np.nan, "BSPV": np.nan,
           "BS_interna": np.nan, "BS_externa": np.nan, "BS_total": np.nan,
           "BS_frac_externa": np.nan,
           "TbTh": np.nan, "TbSp": np.nan,
           "TbN": np.nan, "BS_tapas": np.nan, "V_mc": np.nan,
           "V_rep": np.nan,
           "dif_BV_vs_voxel_pct": np.nan, "perdida_vs_mc_pct": np.nan,
           "n_tri": 0, "estanca": False, "bordes_abiertos": -1,
           "bordes_abiertos_recorte": -1, "reparado_tras_recorte": False}
    if not BW.any() or BW.all():
        m.update(nan)
        return m

    # Planos de la caja de centros de voxel. La normal indica el lado que se
    # queda: +e en el plano bajo, -e en el alto.
    lo = np.zeros(3)
    hi = (n - 1.0) * spacing

    def recortar(malla):
        """Las seis caras de una vez con vtkClipClosedSurface.

        Recortar plano a plano con `clip_closed_surface` de PyVista falla en
        mallas grandes: sobre el VOI real de 97^3 (470 000 triangulos) el
        tapado de un recorte deja aristas abiertas que el siguiente rechaza
        como "non-manifold". Con una coleccion de seis planos VTK tapa una
        sola vez y de forma consistente. Aun asi puede dejar unas pocas
        aristas abiertas (34 sobre 461 000 triangulos en ese VOI). NO se
        rematan por defecto: rematarlas con PyMeshFix costo 112 s y se comio
        un 0.7 % del volumen para cerrar un 0.007 % de las aristas. Se
        reportan en `bordes_abiertos_recorte` y solo se repara si superan el
        0.1 % de los triangulos, que es cuando el volumen dejaria de ser
        fiable; entonces `reparado_tras_recorte` lo dice.
        """
        import vtk
        planos = vtk.vtkPlaneCollection()
        for ax in range(3):
            for origen, signo in ((lo, 1.0), (hi, -1.0)):
                p = vtk.vtkPlane()
                o = [0.0, 0.0, 0.0]; nrm = [0.0, 0.0, 0.0]
                o[ax] = float(origen[ax]); nrm[ax] = signo
                p.SetOrigin(o); p.SetNormal(nrm)
                planos.AddItem(p)
        f = vtk.vtkClipClosedSurface()
        f.SetInputData(malla)
        f.SetClippingPlanes(planos)
        f.SetGenerateFaces(True)
        f.SetGenerateOutline(False)
        f.SetTolerance(1e-9)
        f.Update()
        r = pv.wrap(f.GetOutput()).triangulate().clean()
        abiertos = int(r.n_open_edges) if r.n_cells else 0
        remate = False
        if r.n_cells and abiertos > 0.001 * r.n_cells:
            mfr = pymeshfix.MeshFix(r)
            mfr.repair()
            r = mfr.mesh.triangulate().clean()
            remate = True
        return r, abiertos, remate

    sup = superficie_cerrada(BW, spacing)
    V_mc = float(sup.volume)                 # cruda, cerrada, caja completa

    s = sup
    if suavizado:
        s = s.smooth_taubin(n_iter=int(suavizado), pass_band=float(banda))
    mf = pymeshfix.MeshFix(s.triangulate().clean())
    mf.repair()
    rep = mf.mesh                            # suavizada + reparada, cerrada
    V_rep = float(abs(rep.volume))           # misma caja completa que V_mc
    # La estanqueidad que importa es la de ESTA malla: es la superficie de
    # la estructura. El recorte a la caja es un paso de medida.
    m["estanca"] = bool(rep.is_manifold and rep.n_open_edges == 0)
    m["bordes_abiertos"] = int(rep.n_open_edges)

    s, abiertos, remate = recortar(rep)
    m["bordes_abiertos_recorte"] = int(abiertos)
    m["reparado_tras_recorte"] = bool(remate)
    if s.n_cells == 0:
        m.update(nan)
        return m

    # Tapas: celdas cuyo centro cae en uno de los seis planos de la caja
    c = np.asarray(s.cell_centers().points, dtype=float)
    tol = 1e-6 * float(np.min(spacing))
    en_tapa = np.zeros(s.n_cells, dtype=bool)
    for ax in range(3):
        en_tapa |= np.abs(c[:, ax] - lo[ax]) < tol
        en_tapa |= np.abs(c[:, ax] - hi[ax]) < tol
    areas = np.asarray(s.compute_cell_sizes(length=False, area=True,
                                            volume=False).cell_data["Area"],
                       dtype=float)
    BS_tapas = float(np.abs(areas[en_tapa]).sum())
    BS = float(np.abs(areas).sum() - BS_tapas)
    BV = float(abs(s.volume))

    m.update({"BV": BV, "BS": BS, "BS_tapas": BS_tapas, "V_mc": V_mc,
              "V_rep": V_rep, "n_tri": int(s.n_cells)})
    m["BVTV"] = BV / TV
    m["PoTot"] = (1.0 - m["BVTV"]) * 100.0
    # Dos numeros distintos, los dos en la caja COMPLETA y sin recortar nada:
    #   dif_BV_vs_voxel: DEFINICION, MC a nivel 0.5 frente a conteo de cubos.
    #   perdida_vs_mc:   lo que suavizar + reparar quitan a la superficie
    #                    cruda. Medido en el VOI proximal de H4: -0.4 %.
    # (Una version anterior recortaba tambien la cruda para comparar en la
    # caja de centros: el recorte de una superficie escalonada deja cientos
    # de aristas abiertas y la referencia salia rota. No hace falta recortar
    # para medir una perdida: basta con comparar las dos mallas cerradas.)
    m["dif_BV_vs_voxel_pct"] = (100.0 * (V_mc / mv["BV"] - 1.0)
                                if mv["BV"] > 0 else np.nan)
    m["perdida_vs_mc_pct"] = (100.0 * (V_rep / V_mc - 1.0) if V_mc > 0
                              else np.nan)
    if BS > 0 and BV > 0:
        m["BSBV"] = BS / BV
        m["BSTV"] = BS / TV
        # Las mismas tres metricas que en el modo voxel, pero calculadas con
        # LA SUPERFICIE DE ESTE MODO. Copiarlas del otro mezclaria dos
        # definiciones de BS que difieren un 9 % (correccion C4): aqui la
        # interna es la malla suavizada sin tapas y la externa son las tapas
        # que el recortador fabrica, que ya estaban medidas por separado.
        if TV > BV:
            m["BSPV"] = BS / (TV - BV)
        m["BS_interna"] = BS
        m["BS_externa"] = BS_tapas
        m["BS_total"] = BS + BS_tapas
        if (BS + BS_tapas) > 0:
            m["BS_frac_externa"] = BS_tapas / (BS + BS_tapas)
        m["TbTh"] = 2.0 * BV / BS               # Parfitt, modelo de placas
        m["TbSp"] = m["TbTh"] * ((1.0 / m["BVTV"]) - 1.0)
        m["TbN"] = m["BVTV"] / m["TbTh"]
    else:
        for k in ("BSBV", "BSTV", "TbTh", "TbSp", "TbN"):
            m[k] = np.nan
    return m


def conectividad(BW, spacing):
    """Densidad de conectividad (Conn.D), Odgaard & Gundersen 1993.

        Conn = 1 - chi        Conn.D = Conn / TV     [1/mm^3]

    siendo chi la caracteristica de Euler de la fase solida. Cuenta las
    conexiones REDUNDANTES de la red: cuantas trabeculas se podrian cortar
    antes de partirla en dos.

    POR QUE IMPORTA Y POR QUE NO LA SUSTITUYE BV/TV. En perdida osea las
    trabeculas no solo adelgazan: se desconectan. Dos estructuras con el mismo
    BV/TV y el mismo Tb.Th pueden tener redes completamente distintas, y la
    resistencia depende de la red, no solo de cuanto material hay. Es la
    metrica que separa "hueso mas fino" de "hueso roto".

    Se usa conectividad 26 para el solido, que es la convencion de Odgaard
    (emparejada con 6 para el fondo) y la que aplican CTAn y Scanco.
    """
    BW = np.asarray(BW, dtype=bool)
    spacing = np.asarray(spacing, dtype=float).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    TV = float(np.prod(np.array(BW.shape, float) * spacing))
    if not BW.any() or TV <= 0:
        return {"euler": np.nan, "Conn": np.nan, "ConnD": np.nan}
    chi = float(measure.euler_number(BW, connectivity=3))
    conn = 1.0 - chi
    return {"euler": chi, "Conn": conn, "ConnD": conn / TV}


def indice_smi(BW, spacing, dr_rel=0.05, metodo="lorensen"):
    """Structure Model Index (Hildebrand & Ruegsegger 1997).

        SMI = 6 * BV * (dBS/dr) / BS^2

    Mide si la estructura se parece mas a placas o a barras dilatandola un
    poco y viendo como cambia la superficie: en una placa el area apenas
    cambia (SMI ~ 0), en una barra crece linealmente (SMI ~ 3) y en una esfera
    crece mas deprisa (SMI ~ 4).

    LA CRITICA, QUE HAY QUE REPORTAR JUNTO AL NUMERO. Salmon et al. (2015)
    mostraron que el SMI esta confundido por la CONCAVIDAD: las superficies
    concavas aportan area negativa al dilatar, asi que en hueso denso el indice
    sale negativo y deja de tener la interpretacion placa-barra. No es un
    artefacto de calculo, es una limitacion de la definicion. Se incluye porque
    se sigue pidiendo en la literatura, no porque sea fiable a BV/TV alto.

    La dilatacion se hace desplazando los vertices de la isosuperficie a lo
    largo de sus normales, no con morfologia sobre voxeles: un paso de un voxel
    entero seria demasiado grosero para una derivada.

    SESGO DE DISCRETIZACION, MEDIDO. Sobre geometrias de respuesta conocida:

        esfera   +3.45  (exacto 4)   -13.8%
        cilindro +2.68  (exacto 3)   -10.8%
        losa      0.00  (exacto 0)     exacto

    El deficit no es un fallo del calculo: marching cubes sobreestima BS un
    8.9% sobre datos binarios y el SMI va como 1/BS^2, asi que cabe esperar
    -15.7%. La losa sale exacta porque su dBS/dr es cero y el sesgo de BS no
    tiene por donde entrar. Toda implementacion sobre voxeles arrastra esto.
    """
    BW = np.asarray(BW, dtype=bool)
    spacing = np.asarray(spacing, dtype=float).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    if not BW.any() or BW.all():
        return {"SMI": np.nan, "dBSdr": np.nan}

    verts, faces, normales, _ = measure.marching_cubes(
        BW.astype(np.float32), level=0.5,
        spacing=tuple(float(s) for s in spacing), method=metodo)

    def area(v):
        a, b, c = v[faces[:, 0]], v[faces[:, 1]], v[faces[:, 2]]
        return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum())

    # Las normales de marching_cubes apuntan hacia FUERA del solido en la
    # convencion de skimage, asi que dilatar es sumar en su direccion. Con el
    # signo al reves el SMI sale negativo en todo —una esfera daba -3.36 en vez
    # de +4—, que es indistinguible del negativo legitimo de las superficies
    # concavas y por tanto un error que no salta a la vista.
    n = normales / np.maximum(np.linalg.norm(normales, axis=1, keepdims=True),
                              1e-30)
    dr = float(dr_rel * np.min(spacing))
    A0 = area(verts)
    A1 = area(verts + dr * n)
    dBSdr = (A1 - A0) / dr

    BV = float(BW.sum() * np.prod(spacing))
    smi = 6.0 * BV * dBSdr / (A0 ** 2) if A0 > 0 else np.nan
    return {"SMI": float(smi), "dBSdr": float(dBSdr), "BS_smi": A0}


def morfometria(BW, spacing, do_mil=True, modelo_tbth="plate", do_poro=False,
                metodo_mc="lorensen", extra=False, do_ef=False,
                muestreo_mil="voxel"):
    """Port de localMorphometryFromVoxels.

    Es el UNICO sitio donde se calculan BV/TV, BS/BV, Tb.Th, Tb.Sp, Tb.N y DA,
    igual que en la app: cualquier metrica nueva debe pasar por aqui, porque el
    sesgo del metodo solo se cancela si VOI y candidato recorren el mismo
    camino (correccion C4).

    spacing en mm/voxel. Devuelve un dict con las claves del struct de MATLAB.
    """
    BW = np.asarray(BW, dtype=bool)
    spacing = np.asarray(spacing, dtype=float).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    sz = np.array(BW.shape, dtype=float)

    m = {}
    m["BVTV"] = float(BW.sum() / BW.size)
    m["PoTot"] = (1.0 - m["BVTV"]) * 100.0
    m["TV"] = float(np.prod(sz * spacing))
    m["BV"] = m["BVTV"] * m["TV"]

    # Valores por defecto si algo falla (misma forma de struct en toda salida)
    m.update({"BS": np.nan, "BSBV": np.nan, "BSTV": np.nan, "BSPV": np.nan,
              "BS_interna": np.nan, "BS_externa": np.nan, "BS_total": np.nan,
              "BS_frac_externa": np.nan,
              "TbTh": np.nan, "TbSp": np.nan, "TbN": np.nan,
              "DA": np.nan, "DA2": np.nan, "FracPort": np.nan,
              "dir_principal": np.array([0.0, 0.0, 1.0]),
              "eigenvalues": np.array([1.0, 1.0, 1.0]),
              "eigenvectors": np.eye(3),
              "MIL_valid": False, "MIL_clamped": False})

    if m["BVTV"] <= 0.0 or m["BVTV"] >= 1.0:
        return m

    m["FracPort"] = fraccion_portante(BW)

    BS = area_superficie(BW, spacing, metodo=metodo_mc)
    if BS > 0:
        m["BS"] = BS
        m["BSBV"] = BS / m["BV"]
        m["BSTV"] = BS / m["TV"]
        # Parfitt (placas), como CTAn/Scanco. 'rod' duplicaria Tb.Th.
        m["TbTh"] = (4.0 if modelo_tbth == "rod" else 2.0) * m["BV"] / BS
        m["TbSp"] = m["TbTh"] * ((1.0 / m["BVTV"]) - 1.0)
        m["TbN"] = m["BVTV"] / m["TbTh"]

    # S/V de la fase PORO. En hueso se usa BS/BV -superficie especifica del
    # solido, Parfitt- pero la literatura de materiales porosos para implantes
    # pide la misma superficie dividida por el volumen de POROS, que es la que
    # gobierna transporte y adhesion celular. Es la misma interfaz dividida por
    # el otro volumen, asi que sale gratis y no se toca BS/BV.
    if np.isfinite(m.get("BS", np.nan)) and m["BVTV"] < 1.0:
        m["BSPV"] = m["BS"] / (m["TV"] * (1.0 - m["BVTV"]))

    if extra:
        # Conn.D y SMI se calculan solo si se piden: el SMI necesita otra
        # marching cubes completa y en el bucle de ajuste eso son 53
        # evaluaciones de mas sin que ninguna entre en la funcion de error.
        m.update(conectividad(BW, spacing))
        m.update(indice_smi(BW, spacing, metodo=metodo_mc))

        # La separacion interna/externa cuesta una marching cubes mas, y el
        # tamano de poro una transformada de distancia por radio. Ninguna de
        # las dos entra en la funcion de error, asi que van aqui: en el bucle
        # de ajuste serian cientos de evaluaciones regaladas.
        interna, externa, total = area_superficie_partes(
            BW, spacing, metodo=metodo_mc)
        m["BS_interna"] = interna
        m["BS_externa"] = externa
        m["BS_total"] = total
        # Cuanto del area medida es huella del corte. Si esto es grande, el VOI
        # es pequeno para la estructura que contiene.
        m["BS_frac_externa"] = (externa / total) if total > 0 else np.nan
        if do_poro:
            m.update(tamano_poro(BW, spacing))
        if do_ef:
            # Ellipsoid Factor (Doube 2015): placa/barra sin la confusion por
            # concavidad del SMI. Es, con diferencia, la medida mas cara del
            # modulo -una optimizacion por semilla-, por eso tiene su propia
            # compuerta y nunca entra en el bucle de ajuste.
            from .elipsoide import factor_elipsoide
            m.update(factor_elipsoide(BW, spacing))

    if do_mil:
        # muestreo_mil="suavizado" es la correccion M1 (opcional, ver
        # tensor_mil); VOI y candidato deben usar el mismo valor (C4).
        fab = tensor_mil(BW, spacing, muestreo=muestreo_mil)
        if fab["valid"]:
            m["DA"] = fab["DA"]
            m["dir_principal"] = fab["dir_principal"]
            m["eigenvalues"] = fab["eigenvalues"]
            m["eigenvectors"] = fab["eigenvectors"]
            m["MIL_valid"] = True
            m["MIL_clamped"] = fab["clamped"]
            m["MIL"] = fab["MIL"]
            m["MIL_dirs"] = fab["dirs"]
            # K1: DA solo usa el autovalor mayor y el menor, y no distingue
            # una estructura columnar de una laminar. DA2 = sqrt(l2/l1) es el
            # dato que falta y sale gratis del mismo tensor.
            m["DA2"] = float(np.sqrt(fab["eigenvalues"][1] / fab["eigenvalues"][0]))

    return m
