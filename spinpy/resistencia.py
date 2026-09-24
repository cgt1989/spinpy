"""
resistencia.py — Ensayo de compresion axial y carga de fallo aparente.

Cubre lo que a la rigidez le falta: para biomecanica, saber que un hueso es
rigido es media historia; la pregunta que suele importar es a que carga rompe.

CRITERIO DE PISTOIA
-------------------
Pistoia et al. (Bone 30:842, 2002). Sobre una UNICA resolucion LINEAL se busca
el factor por el que hay que escalar la carga para que un porcentaje del tejido
supere una deformacion critica:

    el hueso falla cuando el 2% del volumen oseo supera el 0.7% de
    deformacion efectiva

Como el problema es lineal, la deformacion escala con la carga y el factor sale
de una division: no hace falta analisis no lineal ni iterar. El criterio se
valido originalmente contra ensayos de radio distal humano y predice la carga
de rotura mucho mejor que la densidad sola.

La deformacion efectiva es la definicion energetica que usa esa literatura:

    eps_eff = sqrt(2 * U / E_s)      con U la densidad de energia de deformacion

Para un estado uniaxial se reduce a |eps_zz|, lo que permite comprobarla contra
una solucion exacta.

LOS DOS PARAMETROS SON CONVENCIONES, NO CONSTANTES FISICAS
-----------------------------------------------------------
El 2% y el 0.7% se calibraron en radio distal humano a una resolucion concreta.
No son universales: dependen de la especie, del sitio y —esto importa— de la
RESOLUCION de la malla, porque el percentil de una cola depende de cuantos
elementos haya. Se exponen como argumentos y se reportan siempre junto al
resultado, en lugar de esconderlos como constantes.

CONDICIONES DE CONTORNO
-----------------------
  'empotrado'  base totalmente fija: reproduce un ensayo con friccion entre
               probeta y platos, introduce coaccion lateral y rigidiza.
  'deslizante' uz = 0 en la base y solo lo minimo para eliminar los movimientos
               de solido rigido: compresion sin friccion, estado uniaxial.

La eleccion cambia E_app y la carga de fallo, asi que hay que declararla.

TRES CAMPOS DE SALIDA, NO UNO
------------------------------
Del mismo campo de tensiones salen `campo_eps_eff` (deformacion efectiva de
Pistoia, energetica) y `campo_vm` (von Mises, solo la parte desviadora). No son
el mismo mapa reescalado: coinciden salvo un factor en estado uniaxial y se
separan en los nudos trabeculares, donde el estado es triaxial. Uno predice la
carga aparente de fallo y el otro localiza donde plastificaria el tejido.

El tercero es `campo_desp`, la DEFORMACION TOTAL: el modulo del vector
desplazamiento, |u| = sqrt(ux^2 + uy^2 + uz^2), en las unidades de `spacing`.
Es la variable "Total Deformation" que reporta ANSYS y la que usan los
estudios comparativos microCT-FEA junto a la von Mises. No es una tension ni
una deformacion unitaria: es un desplazamiento, y su maximo esta siempre en la
cara cargada porque acumula todo lo que se ha movido la probeta. Mide RIGIDEZ
GLOBAL, mientras la von Mises localiza la concentracion de tension.

TRES EJES
---------
`ensayo_compresion_eje` corre el mismo ensayo en X, Y o Z permutando la
estructura, y `ensayo_triaxial` los tres de una vez. El cociente E_max/E_min
que sale de ahi es una anisotropia MECANICA, contrastable con Ex/Ey/Ez de la
homogeneizacion periodica y distinta del DA del MIL, que es geometrico.
"""

from __future__ import annotations

import time

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, cg, splu

from .elastic import UMBRAL_DIRECTO, VOID_SCALE, hex8_ke

E_S_DEF, NU_S_DEF = 20e9, 0.30
SIGMA0_DEF = 1e6          # 1 MPa de tension aparente de referencia

# Claves CANONICAS del apoyo del ensayo de compresion. Son las unicas que
# entienden `ensayo_compresion` y `escribe.escribir_febio`, y las que se
# guardan en los registros exportados. El visor traduce el texto del combo
# ("empotrado" -> "fixed" con la interfaz en ingles), asi que NUNCA debe
# pasar `currentText()`: antes se decidia con startswith("empotr") y "fixed"
# ejecutaba en silencio el apoyo DESLIZANTE mientras el registro decia
# "fixed". Ahora un texto desconocido lanza ValueError.
APOYOS = ("deslizante", "empotrado")


def normalizar_apoyo(apoyo):
    """Devuelve la clave canonica de `apoyo` o lanza ValueError.

    Admite mayusculas y espacios alrededor, nada mas: ni traducciones
    ("fixed", "sliding") ni prefijos. Un apoyo mal nombrado cambia E_app y la
    carga de fallo sin dar ningun sintoma, asi que vale mas un error.
    """
    clave = str(apoyo).strip().lower()
    if clave not in APOYOS:
        raise ValueError(
            "apoyo desconocido: %r (se admite %s)"
            % (apoyo, " o ".join(repr(a) for a in APOYOS)))
    return clave
FRAC_CRITICA = 0.02       # 2% del volumen oseo
EPS_CRITICA = 0.007       # 0.7% de deformacion efectiva

# Umbrales del veredicto de `estudio_convergencia`. Los lee tambien
# `informe.comprobar`, para que el informe y el dialogo no puedan discrepar.
# Son CRITERIO NUESTRO, no una medida: por encima de un 2 % de deriva de
# densidad el remuestreo ya cambia la estructura (no se extrapola), y un error
# estimado frente a Richardson por debajo del 5 % es lo que se da por meseta.
CONV_DERIVA_RHO_MAX = 0.02
CONV_ERROR_EXTRAPOLADO_MAX = 0.05
# Una serie NO monotona cuya dispersion total (max - min) / media no pasa de
# este valor oscila sin tendencia dentro de la banda en la que el VOI proximal
# de H4 se dio por convergido (+-3 % desde n = 40, Estudio_Convergencia). No
# se extrapola —el orden no tiene sentido—, pero tampoco es la serie que «no
# converge en absoluto»: medido en el VOI porcino V1, 4041/4072/4054/4021 MPa
# de 22^3 a 40^3, 1.3 % de dispersion. Usar la banda ENTERA como dispersion
# total es la lectura estricta del +-3 %; es CRITERIO NUESTRO.
CONV_OSCILACION_MAX = 0.03

ESQUINAS = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                     [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])


def _b_centro(dx, dy, dz):
    """Matriz B en el centro del elemento (xi = eta = zeta = 0)."""
    xn = np.array([-1, 1, 1, -1, -1, 1, 1, -1], float)
    yn = np.array([-1, -1, 1, 1, -1, -1, 1, 1], float)
    zn = np.array([-1, -1, -1, -1, 1, 1, 1, 1], float)
    dNdx = np.diag([2 / dx, 2 / dy, 2 / dz]) @ (0.125 * np.vstack([xn, yn, zn]))
    B = np.zeros((6, 24))
    B[0, 0::3] = dNdx[0]
    B[1, 1::3] = dNdx[1]
    B[2, 2::3] = dNdx[2]
    B[3, 1::3] = dNdx[2]; B[3, 2::3] = dNdx[1]
    B[4, 0::3] = dNdx[2]; B[4, 2::3] = dNdx[0]
    B[5, 0::3] = dNdx[1]; B[5, 1::3] = dNdx[0]
    return B


def _matriz_D(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[0, 0] = D[1, 1] = D[2, 2] = lam + 2 * mu
    D[3, 3] = D[4, 4] = D[5, 5] = mu
    return D


def _modos_rigidos(coord):
    """Los seis modos de solido rigido, como matriz (ndof, 6).

    Tres traslaciones y tres rotaciones infinitesimales. Son el espacio nulo de
    la matriz de rigidez sin apoyos, y el multigrid algebraico los necesita
    para construir bien los niveles gruesos en problemas de elasticidad.
    """
    n = coord.shape[0]
    x, y, z = coord[:, 0], coord[:, 1], coord[:, 2]
    B = np.zeros((3 * n, 6))
    B[0::3, 0] = 1.0                      # traslacion x
    B[1::3, 1] = 1.0                      # traslacion y
    B[2::3, 2] = 1.0                      # traslacion z
    B[0::3, 3] = -y; B[1::3, 3] = x       # giro alrededor de z
    B[1::3, 4] = -z; B[2::3, 4] = y       # giro alrededor de x
    B[0::3, 5] = z;  B[2::3, 5] = -x      # giro alrededor de y
    return B


def _solo_portante(BW):
    """Conserva el material en componentes que unen la base con el techo.

    Un fragmento oseo que no llega a ninguno de los dos platos no transmite
    carga: en el modelo solo aporta una submatriz casi singular sujeta por la
    rigidez del vacio.
    """
    from scipy import ndimage
    c6 = ndimage.generate_binary_structure(3, 1)
    lab, n = ndimage.label(BW, structure=c6)
    if n == 0:
        return BW
    ini = set(np.unique(lab[:, :, 0])) - {0}
    fin = set(np.unique(lab[:, :, -1])) - {0}
    p = ini & fin
    return np.isin(lab, list(p)) if p else BW


def _pct(v, q):
    """Percentil con la convencion de `prctile` de MATLAB.

    Posiciones de trazado (i-0.5)/n, que numpy llama 'hazen'. El defecto de
    numpy, (i-1)/(n-1), NO es lo mismo: sobre los VOI de H4 desplaza el p99
    entre un 0.4 y un 0.8 %.
    """
    v = np.asarray(v, float)
    return float(np.percentile(v, q, method="hazen")) if v.size else float("nan")


def _pesos_iguales(w):
    return w is None or (np.size(w) > 0 and bool(np.all(w == np.ravel(w)[0])))


def percentil_ponderado(v, w, q, metodo="hazen"):
    """Percentil ponderado por `w` (volumen de cada elemento).

    Con voxeles todos los elementos pesan lo mismo y el percentil por conteo
    YA esta ponderado por volumen. Con tetraedros no: «el 2 % del volumen
    oseo» (Pistoia) o el p99 de la capa superficial tienen que contar volumen,
    no elementos.

    DEFINICION. Se ordenan los valores; cada elemento ocupa un tramo de peso
    w_i y se situa en su punto MEDIO, c_i = S_i - w_i / 2 (S_i, peso
    acumulado). Entre esos puntos se interpola linealmente y fuera se toma el
    extremo:
      'hazen'   posicion c_i / W: con pesos iguales es (i - 0.5)/n, la de
                `prctile` de MATLAB (Q5 de Hyndman y Fan).
      'linear'  posicion (c_i - c_1)/(c_n - c_1): con pesos iguales es
                (i - 1)/(n - 1), el defecto de numpy (Q7), que es el que usan
                el umbral de Pistoia y el vm_p99 de todo el tejido.

    Con pesos iguales (o `w` None) se llama a `np.percentile` con el mismo
    metodo, de modo que el resultado es BIT A BIT el de siempre: el bloque 27
    lo comprueba, y comprueba tambien que la formula ponderada con pesos
    iguales coincide con numpy a redondeo.
    """
    v = np.asarray(v, float).ravel()
    if v.size == 0:
        return float("nan") if np.ndim(q) == 0 else np.full(np.shape(q), np.nan)
    if _pesos_iguales(w):
        r = np.percentile(v, q, method=metodo)
        return float(r) if np.ndim(r) == 0 else r
    w = np.asarray(w, float).ravel()
    if w.size != v.size or (w < 0).any():
        raise ValueError("pesos incompatibles con los valores")
    return _percentil_formula(v, w, q, metodo)


def _percentil_formula(v, w, q, metodo):
    orden = np.argsort(v, kind="stable")
    v, w = v[orden], w[orden]
    S = np.cumsum(w)
    c = S - 0.5 * w
    if metodo == "hazen":
        pos = c / S[-1]
    elif metodo == "linear":
        pos = (c - c[0]) / (c[-1] - c[0]) if v.size > 1 else np.zeros(1)
    else:
        raise ValueError(f"metodo desconocido: {metodo!r}")
    r = np.interp(np.asarray(q, float) / 100.0, pos, v)
    return float(r) if np.ndim(r) == 0 else r


def _media(v, w):
    if _pesos_iguales(w):
        return float(v.mean())
    return float(np.average(v, weights=w))


def capa_superficie(BW):
    """Elementos solidos con al menos un vecino POR CARA en el vacio.

    Es la capa donde vive la escalera de la voxelizacion y donde se concentra
    la tension. Importa porque el maximo de von Mises NO CONVERGE con la
    resolucion —ni siquiera contra solucion cerrada: en la cavidad esferica
    periodica el pico binario se pasa un +5.6 % del valor de Goodier a la
    resolucion mas fina y va dando saltos de ~10 puntos por el camino
    (`Estudio_Convergencia`)—, mientras que un percentil 99 TOMADO SOBRE ESTA
    CAPA se queda a -0.3 % del exacto. Un p99 sobre TODO el tejido tampoco
    sirve en estructuras densas: lo domina el material a granel, donde la
    tension es la nominal, y no dice nada del pico.

    OJO — ESTA DEFINICION NO ES PERIODICA, y la de
    `Estudio_Convergencia/code/conv_nucleo.py` SI LO ES.
    Alli la celda se homogeneiza con condiciones periodicas, asi que un vecino
    que sale por una cara entra por la opuesta y `np.roll` es lo correcto.
    Aqui el ensayo tiene base y techo: envolver compararia la cara cargada con
    la apoyada y marcaria como superficie libre justo los elementos donde la
    condicion de contorno concentra tension por si misma. Fuera del dominio se
    supone material, de modo que las seis caras del VOI no cuentan como
    interfaz —la misma regla que excluye las isocaps del area superficial en
    `morphometry`—. Las dos definiciones coinciden en el interior y solo se
    separan en el borde; no son intercambiables, y por eso no se comparte el
    codigo.
    """
    B = np.asarray(BW, dtype=bool)
    if B.ndim != 3:
        raise ValueError("capa_superficie espera una mascara 3D.")
    P = np.pad(B, 1, mode="constant", constant_values=True)
    vecino_vacio = np.zeros(B.shape, dtype=bool)
    for eje in (0, 1, 2):
        for desp in (1, -1):
            vecino_vacio |= ~np.roll(P, desp, axis=eje)[1:-1, 1:-1, 1:-1]
    return B & vecino_vacio


def estadisticos_vm(res):
    """Estadisticos de von Mises en la capa superficial de un ensayo resuelto.

    Toma el diccionario que devuelve `ensayo_compresion`. Usa `vm_solido` y la
    mascara `superficie_solido` que ese mismo ensayo deja preparada, las dos en
    el mismo orden y en doble precision. Para resultados ANTIGUOS, que no
    llevan esa mascara, la deduce de `campo_vm`; sale lo mismo salvo que el
    campo se guarda en float32, asi que los valores pueden diferir en la
    septima cifra. Se avisa aqui porque si no, un `vm_max_superficie` que no
    coincide exactamente con `vm_max` parece un error y no lo es.

    Devuelve `vm_p99_superficie` —el estadistico que el estudio de
    convergencia valida— y, para poder juzgarlo, `vm_n_superficie`: un
    percentil de la cola de una muestra pequena es ruido. Un dominio
    completamente macizo no tiene superficie libre interior y devuelve NaN con
    n = 0; es correcto, no hay concentracion que medir.

    Los percentiles de aqui usan la convencion de `prctile` (ver `_pct`).

    Si `res` trae `vol_solido` —el volumen de cada elemento, alineado con
    `vm_solido`; lo ponen los ensayos en FEBio con tetraedros— el percentil,
    la media y la fraccion se PONDERAN por volumen (`percentil_ponderado`).
    Con volumenes iguales, o sin la clave, el resultado es bit a bit el de
    siempre.
    """
    vm = np.asarray(res.get("vm_solido", []), float)
    sup = res.get("superficie_solido")
    w = res.get("vol_solido")
    if sup is None:
        campo = res.get("campo_vm")
        if campo is None:
            return {}
        campo = np.asarray(campo, float)
        dentro = np.isfinite(campo)
        vm = campo[dentro]
        sup = capa_superficie(dentro)[dentro]
        w = None
    sup = np.asarray(sup, bool).ravel()
    if vm.size == 0 or sup.size != vm.size:
        return {}
    if w is not None:
        w = np.asarray(w, float).ravel()
        if w.size != vm.size:
            return {}
    ok = np.isfinite(vm)
    vs = vm[sup & ok]
    ws = None if w is None else w[sup & ok]
    if _pesos_iguales(ws) and _pesos_iguales(w):
        frac = (float(vs.size) / float(ok.sum()) if ok.any()
                else float("nan"))
    else:
        frac = float(ws.sum() / w[ok].sum()) if ok.any() else float("nan")
    return {"vm_p99_superficie": (percentil_ponderado(vs, ws, 99)
                                  if vs.size else float("nan")),
            "vm_max_superficie": float(vs.max()) if vs.size else float("nan"),
            "vm_media_superficie": _media(vs, ws) if vs.size
            else float("nan"),
            "vm_n_superficie": int(vs.size),
            "vm_frac_superficie": frac}


# Probabilidades (en %) a las que se guardan los cuantiles de von Mises de la
# capa superficial: cada punto entero y, en la cola, cada decima. La cola es
# lo que se cita (p99), asi que es donde hace falta resolucion; 111 numeros
# por estructura caben en el JSON sin guardar el campo.
PROB_CUANTILES_VM = tuple([float(q) for q in range(0, 99)]
                          + [round(99.0 + 0.1 * k, 1) for k in range(11)])


def cuantiles_vm_superficie(res, probs=PROB_CUANTILES_VM, escala=1.0):
    """Distribucion de von Mises en la capa superficial, como cuantiles.

    Misma muestra que `estadisticos_vm` —el tejido de `superficie_solido`— y
    misma convencion de percentil (hazen, la de `prctile`), de modo que el
    cuantil 99 de la lista coincide exactamente con `vm_p99_superficie`. Es lo
    que permite dibujar la cola de la que sale el valor citable a partir del
    JSON exportado, sin guardar el campo de cientos de MB.

    `escala` multiplica los valores (1e-6 para pasar de Pa a MPa). Devuelve
    {"p": [...], "valor": [...], "n": int} o {} si no hay capa superficial.
    """
    vm = np.asarray(res.get("vm_solido", []), float)
    sup = res.get("superficie_solido")
    w = res.get("vol_solido")
    if sup is None:
        campo = res.get("campo_vm")
        if campo is None:
            return {}
        campo = np.asarray(campo, float)
        dentro = np.isfinite(campo)
        vm = campo[dentro]
        sup = capa_superficie(dentro)[dentro]
        w = None
    sup = np.asarray(sup, bool).ravel()
    if vm.size == 0 or sup.size != vm.size:
        return {}
    m = sup & np.isfinite(vm)
    vs = vm[m]
    if vs.size == 0:
        return {}
    ws = None if w is None else np.asarray(w, float).ravel()[m]
    q = np.asarray(percentil_ponderado(vs, ws, list(probs))) * float(escala)
    return {"p": [float(x) for x in probs], "valor": [float(x) for x in q],
            "n": int(vs.size)}


def ensayo_compresion(BW, spacing, E_s=E_S_DEF, nu_s=NU_S_DEF,
                      sigma0=SIGMA0_DEF, apoyo="deslizante",
                      escala_vacio=None, tol=1e-8, solo_portante=True,
                      carga_N=None, unidad="mm", progreso=None,
                      modulo_rel=None):
    """Compresion axial en z. Devuelve campos y magnitudes aparentes.

    `modulo_rel` (opcional, misma forma que BW): rigidez relativa de cada
    voxel de hueso, en (0, 1]. Con None —el valor por omision— todo el hueso
    tiene E_s y el resultado es identico al de siempre. Existe para el fallo
    progresivo (`spinpy.simulacion.fallo_progresivo`), que ABLANDA el tejido
    roto en lugar de borrarlo: borrar puntales de uno o dos elementos dejaba
    fragmentos casi como mecanismos y la rigidez subia y bajaba entre pasos
    (medido sobre el VOI proximal de H4 a 32^3). La tension se escala con el
    mismo factor, asi que la energia y eps_eff de un elemento ablandado bajan
    en proporcion.

    La tension aparente `sigma0` se define sobre la seccion BRUTA del VOI
    (lado x lado), no sobre el area osea. Es la convencion de la app y la que
    hace comparables los resultados con la morfometria: aplicarla sobre el area
    real daria una tension mucho mayor y dependiente del BV/TV.

    `carga_N` es la alternativa: fuerza TOTAL en newtons, de la que se deduce
    sigma0 = carga_N / A_bruta. Existe porque los estudios comparativos
    microCT-FEA aplican una fuerza fija (100 N es lo habitual) y no una tension.
    OJO AL COMPARAR: a fuerza fija, dos VOIs de distinta seccion reciben
    tensiones aparentes distintas; a tension fija reciben fuerzas distintas.
    Si los dos volumenes tienen el mismo lado -que es el caso al comparar un
    VOI con su spinodoide ajustado, escalado al tamano fisico del VOI- las dos
    convenciones coinciden. Si no, hay que declarar cual se uso.

    UNIDADES — `unidad` SOLO importa si se da `carga_N`
    ---------------------------------------------------
    El resto de la funcion es indiferente a las unidades: la deformacion es
    adimensional, asi que las tensiones salen en las mismas unidades que
    `E_s` y los desplazamientos en las mismas que `spacing`. Una fuerza, en
    cambio, mezcla las dos: sigma = F / A obliga a que F, E_s y spacing
    pertenezcan al mismo sistema.

    En esta aplicacion `spacing` esta en MILIMETROS y `E_s` en PASCALES, que
    no es un sistema coherente. Para que sigma salga en Pa a partir de un area
    en mm^2 hay que multiplicar la fuerza en newtons por 1e6:

        sigma[Pa] = F[N] / A[m^2] = 1e6 * F[N] / A[mm^2]

    `unidad='mm'` (por defecto) aplica ese factor; `unidad='m'` no lo aplica y
    corresponde a un sistema SI puro. Sin esto, pedir 100 N con el spacing en
    mm aplica en realidad 1e-4 N y todo sale un millon de veces mas pequeno:
    von Mises de 0.000 MPa y desplazamientos de 0.00000 mm.

    EL VACIO NO SE MALLA — y por que importa
    -----------------------------------------
    `homogeneizar` rellena los poros con un material de rigidez 1e-6*E_s para
    que la celda periodica nunca sea singular. Aqui eso NO funciona: medido
    sobre el VOI proximal de H4 a 32^3 (BV/TV 0.28), el solver se quedaba en un
    residuo de 3.3e-1 —sin converger en absoluto— porque el contraste de 1e6
    entre hueso y relleno arruina el condicionamiento en una estructura poco
    densa.

    La salida es la de la practica habitual en micro-elementos finitos (Van
    Rietbergen): mallar SOLO el hueso y renumerar los nodos. El sistema queda
    bien condicionado mientras el hueso este conectado, que es justo lo que
    garantiza `solo_portante`. Ademas es mas rapido: a BV/TV 0.28 hay un 72%
    menos de elementos.

    `escala_vacio` se conserva por compatibilidad; si se le da un valor, se
    vuelve al modelo antiguo con relleno blando.
    """
    def pr(f, m):
        if progreso:
            progreso(f, m)

    apoyo = normalizar_apoyo(apoyo)
    BW = np.asarray(BW, dtype=bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    dx, dy, dz = spacing
    nx, ny, nz = BW.shape
    nel = nx * ny * nz
    v_elem = float(dx * dy * dz)
    nnx, nny, nnz = nx + 1, ny + 1, nz + 1
    ndof = 3 * nnx * nny * nnz

    out = {"ok": False, "msg": "", "apoyo": apoyo, "E_s": E_s,
           "sigma_app": float(sigma0), "rho": float(BW.mean())}
    if not BW.any():
        out["msg"] = "La estructura no contiene material solido."
        return out
    if not BW[:, :, -1].any():
        out["msg"] = ("La capa superior no contiene material: no hay "
                      "superficie por la que introducir la carga.")
        return out

    pr(0.05, "Ensamblando…")
    keS = hex8_ke(dx, dy, dz, E_s, nu_s)

    BW_res = _solo_portante(BW) if solo_portante else BW
    out["frac_portante"] = float(BW_res.sum() / max(BW.sum(), 1))
    if not BW_res.any() or not BW_res[:, :, -1].any():
        out["msg"] = "Tras filtrar, no queda material que una base y techo."
        return out

    nod = lambda a, b, c: a + nnx * b + nnx * nny * c

    if escala_vacio is None:
        # Solo los voxeles de hueso son elementos, y solo sus nodos existen.
        i_, j_, k_ = np.nonzero(BW_res)
        conn = np.stack([nod(i_ + c[0], j_ + c[1], k_ + c[2])
                         for c in ESQUINAS], axis=1)
        usados, conn = np.unique(conn, return_inverse=True)
        conn = conn.reshape(-1, 8)
        n_nod = usados.size
        ndof = 3 * n_nod
        escala = np.ones(conn.shape[0])
        if modulo_rel is not None:
            escala = np.clip(np.asarray(modulo_rel, float)[i_, j_, k_],
                             1e-12, None)
        # coordenadas de los nodos que sobreviven
        kk, resto = np.divmod(usados, nnx * nny)
        jj, ii = np.divmod(resto, nnx)
        coord = np.stack([ii, jj, kk], axis=1) * spacing[None, :]
    else:
        ei, ej, ek = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz),
                                 indexing="ij")
        ei = ei.ravel(order="F"); ej = ej.ravel(order="F"); ek = ek.ravel(order="F")
        conn = np.stack([nod(ei + c[0], ej + c[1], ek + c[2])
                         for c in ESQUINAS], axis=1)
        escala = np.full(conn.shape[0], float(escala_vacio))
        solido_F = BW_res.ravel(order="F")
        escala[solido_F] = (1.0 if modulo_rel is None else np.clip(
            np.asarray(modulo_rel, float).ravel(order="F")[solido_F],
            1e-12, None))
        n_nod = nnx * nny * nnz
        ndof = 3 * n_nod
        kk, resto = np.divmod(np.arange(n_nod), nnx * nny)
        jj, ii = np.divmod(resto, nnx)
        coord = np.stack([ii, jj, kk], axis=1) * spacing[None, :]

    nelem = conn.shape[0]
    out.update({"n_elem": int(nelem), "n_dof": int(ndof)})

    edof = np.empty((nelem, 24), dtype=np.int64)
    for a in range(8):
        for k in range(3):
            edof[:, 3 * a + k] = 3 * conn[:, a] + k

    ke_plana = keS.ravel()
    K = sparse.csr_matrix((ndof, ndof))
    paso = max(1, nelem // 8)
    for i0 in range(0, nelem, paso):
        i1 = min(i0 + paso, nelem)
        ed = edof[i0:i1]
        K = K + sparse.coo_matrix(
            ((ke_plana[None, :] * escala[i0:i1, None]).ravel(),
             (np.repeat(ed, 24, axis=1).ravel(), np.tile(ed, (1, 24)).ravel())),
            shape=(ndof, ndof)).tocsr()
    K = (K + K.T) / 2

    # --- carga sobre los nodos de la cara superior -------------------------
    A_bruta = float((nx * dx) * (ny * dy))
    if carga_N is not None:
        # Fuerza total impuesta: la tension aparente se deduce de la seccion
        # bruta, no al reves. El factor de unidad convierte newtons al sistema
        # implicito de esta llamada (ver el docstring).
        factor = {"mm": 1e6, "m": 1.0}.get(str(unidad).lower())
        if factor is None:
            raise ValueError("unidad debe ser 'mm' o 'm'; se recibio "
                             f"{unidad!r}")
        F_total = float(carga_N) * factor
        sigma0 = F_total / A_bruta if A_bruta > 0 else 0.0
        out["sigma_app"] = float(sigma0)
        out["carga_N"] = float(carga_N)
        out["unidad"] = str(unidad).lower()
    else:
        F_total = float(sigma0) * A_bruta
    z = coord[:, 2]
    tolz = 1e-9 * max(z.max() - z.min(), 1.0)
    techo = np.nonzero(z >= z.max() - tolz)[0]
    base = np.nonzero(z <= z.min() + tolz)[0]
    if techo.size == 0 or base.size == 0:
        out["msg"] = "No hay nodos en la base o en el techo."
        return out

    # Reparto por AREA TRIBUTARIA, no por numero de nodos. Un nodo de esquina
    # sirve a un solo elemento y uno interior a cuatro; darles la misma fuerza
    # sobrecarga los bordes y estropea el campo. Se noto porque un bloque
    # macizo dejaba de tener deformacion uniforme: eps_eff variaba un 125%
    # cuando deberia ser constante, y E_app caia un 2.7%.
    #
    # Se recorre cada elemento cuya cara superior esta en el techo y se le da
    # a cada uno de sus cuatro nodos superiores un cuarto de su parte. Ese es
    # el equivalente nodal consistente de una presion uniforme.
    z_nodo_sup = coord[conn[:, 4], 2]
    elem_techo = np.nonzero(z_nodo_sup >= z.max() - tolz)[0]
    if elem_techo.size == 0:
        out["msg"] = "Ningun elemento llega al techo."
        return out
    F = np.zeros(ndof)
    f_por_nodo = F_total / (elem_techo.size * 4)
    for a in range(4, 8):
        np.add.at(F, 3 * conn[elem_techo, a] + 2, -f_por_nodo)

    # --- contorno ----------------------------------------------------------
    if apoyo == "empotrado":
        fijos = np.concatenate([3 * base, 3 * base + 1, 3 * base + 2])
    else:
        # Apoyo deslizante: uz = 0 en toda la base y SOLO lo minimo para
        # eliminar los movimientos de solido rigido que quedan (traslacion en
        # x e y, y giro alrededor de z).
        #
        # Los dos nodos tienen que ser ESQUINAS GEOMETRICAS, no dos indices
        # cualesquiera de la lista. Al renumerar los nodos del hueso el orden
        # deja de tener relacion con la posicion, y coger `base[0]` y `base[1]`
        # puede dar dos nodos alineados en y: entonces fijar uy en ambos no
        # elimina el giro y ademas sobrerestringe. Se noto porque un bloque
        # macizo dejaba de dar deformacion uniforme (E_app 19.47 en vez de 20).
        cb = coord[base]
        a = base[int(np.argmin(cb[:, 0] + cb[:, 1]))]            # esquina (min x, min y)
        b = base[int(np.argmax(cb[:, 0] - cb[:, 1]))]            # esquina (max x, min y)
        fijos = np.concatenate([3 * base + 2,
                                [3 * a, 3 * a + 1, 3 * b + 1]])
    fijos = np.unique(fijos)
    libres = np.setdiff1d(np.arange(ndof), fijos)

    # --- resolucion --------------------------------------------------------
    pr(0.35, f"Resolviendo {libres.size} grados de libertad…")
    u = np.zeros(ndof)
    Kff = K[libres][:, libres].tocsr()
    Ff = F[libres]
    resuelto = False
    if libres.size <= UMBRAL_DIRECTO:
        try:
            u[libres] = splu(Kff.tocsc()).solve(Ff)
            out["solver"] = "LU directo"
            resuelto = bool(np.all(np.isfinite(u)))
        except Exception:
            resuelto = False
    if not resuelto:
        try:
            import pyamg
            # ESPACIO NULO CERCANO: los seis modos de solido rigido.
            #
            # Sin esto el multigrid agrega como si el problema fuera de tipo
            # Poisson, y en elasticidad eso funciona mal: sobre el VOI proximal
            # de H4 (BV/TV 0.28) el residuo se quedaba en 7e-4 a 32^3 y EMPEORABA
            # al refinar —1.8e-2 a 48^3, 6.3e-2 a 97^3—, que es la firma de una
            # jerarquia mal construida, no de un sistema singular. El solver
            # directo resolvia el mismo caso sin problema.
            #
            # Las traslaciones y rotaciones deben poder representarse EXACTAMENTE
            # en cada nivel grueso; si no, el suavizador no puede eliminar los
            # modos de baja energia y la convergencia se degrada con el tamano.
            B = _modos_rigidos(coord)[libres, :]
            ml = pyamg.smoothed_aggregation_solver(Kff, B=B, max_coarse=500)
            u[libres] = ml.solve(Ff, tol=tol, maxiter=500, accel="cg")
            out["solver"] = "AMG (con modos rigidos) + CG"
        except ImportError:
            d = np.asarray(Kff.diagonal(), float)
            d[d <= 0] = np.finfo(float).eps
            M = LinearOperator(Kff.shape, lambda v, d=d: v / d)
            u[libres], _ = cg(Kff, Ff, rtol=tol, maxiter=20000, M=M)
            out["solver"] = "CG + Jacobi"
        nb = np.linalg.norm(Ff)
        res = np.linalg.norm(Ff - Kff @ u[libres]) / (nb if nb > 0 else 1.0)
        out["residuo_rel"] = float(res)
        # Un solver iterativo no lanza error al no converger: devuelve el mejor
        # iterado, que pasa cualquier isfinite y daria un campo equivocado.
        if not np.isfinite(res) or res > max(1e-6, 100 * tol):
            out["msg"] = (f"El solver no convergio (residuo {res:.1e}). "
                          f"El campo no es fiable.")
            return out
        resuelto = True

    # --- deformaciones y energia -------------------------------------------
    pr(0.85, "Recuperando deformaciones…")
    B0 = _b_centro(dx, dy, dz)
    D = _matriz_D(E_s, nu_s)
    eps = u[edof] @ B0.T                      # (nelem, 6)
    sig = (eps @ D.T) * escala[:, None]

    # Densidad de energia de deformacion, y de ahi la deformacion efectiva de
    # Pistoia. En estado uniaxial esta se reduce a |eps_zz|, lo que da una
    # comprobacion exacta contra un bloque macizo.
    U = 0.5 * np.einsum("ij,ij->i", sig, eps)
    eps_eff = np.sqrt(np.maximum(2.0 * U / E_s, 0.0))

    # Tension equivalente de von Mises, del MISMO campo de tensiones. Sale
    # gratis: `sig` ya esta calculada para la energia.
    #
    # NO es una reetiquetacion de eps_eff. La deformacion efectiva de Pistoia
    # viene de la energia TOTAL, que incluye la parte hidrostatica; von Mises
    # descarta esa parte y solo mide la desviadora. Los dos campos coinciden en
    # un estado uniaxial puro salvo un factor, y se separan justo donde importa:
    # en los nudos trabeculares, donde el estado es triaxial. Por eso el mapa de
    # von Mises localiza la plastificacion y el de Pistoia predice la carga
    # aparente de fallo; se reportan los dos, no uno en lugar del otro.
    s = sig
    vm = np.sqrt(np.maximum(
        0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2
               + (s[:, 2] - s[:, 0]) ** 2)
        + 3.0 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2), 0.0))

    solido = (escala >= 1.0)

    # Campo 3D de deformacion efectiva, para poder pintarlo sobre la geometria.
    # Solo tiene valor donde hay hueso PORTANTE: el material filtrado no se
    # resolvio, y rellenarlo con ceros lo pintaria como si estuviera descargado
    # cuando en realidad no se sabe. Se deja NaN, que no se colorea.
    def a_campo(v):
        c = np.full(BW.shape, np.nan, dtype=np.float32)
        if escala_vacio is None:
            c[np.nonzero(BW_res)] = v
        else:
            c.ravel(order="F")[:] = np.where(solido, v, np.nan)
        return c

    out["campo_eps_eff"] = a_campo(eps_eff)
    out["campo_vm"] = a_campo(vm)

    # Capa superficial del MODELO que se acaba de resolver, en el mismo orden
    # que `vm_solido` / `eps_eff_solido`. Se calcula aqui, y no mas tarde a
    # partir del campo, porque el campo se guarda en float32 y porque aqui se
    # sabe con certeza que geometria se mallo: si `solo_portante` retiro un
    # fragmento, el hueco que deja ES superficie libre en el problema resuelto.
    # El `[solido]` final no sobra: con `modulo_rel` el tejido ablandado tiene
    # escala < 1 y queda fuera de `vm_solido`, asi que la mascara tiene que
    # perder esos mismos elementos o deja de estar alineada. Ablandado NO es
    # vacio —sigue transmitiendo carga—, de modo que cuenta como hueso al
    # decidir quien tiene superficie libre y solo se cae de la MUESTRA.
    if escala_vacio is None:
        out["superficie_solido"] = capa_superficie(
            BW_res)[np.nonzero(BW_res)][solido]
    else:
        out["superficie_solido"] = capa_superficie(BW).ravel(order="F")[solido]

    # DEFORMACION TOTAL: modulo del desplazamiento nodal, promediado a los
    # ocho nodos de cada elemento para poder pintarlo sobre la misma rejilla
    # que los otros dos campos. El promedio por elemento es lo que hace que
    # las tres salidas sean comparables entre si; el maximo NODAL, que es el
    # que reporta ANSYS como "Total Deformation", se devuelve aparte sin
    # promediar para no perderlo.
    desp_nodal = np.linalg.norm(u.reshape(-1, 3), axis=1)
    desp_elem = desp_nodal[conn].mean(axis=1)
    out["campo_desp"] = a_campo(desp_elem)
    out["desp_solido"] = desp_elem
    out["desp_max"] = float(desp_nodal.max())
    out["desp_media"] = float(desp_elem.mean())
    out["desp_max_elem"] = float(desp_elem.max())

    out.update({
        "ok": True, "u": u, "eps": eps, "sig": sig,
        "eps_eff": eps_eff, "SED": U, "solido": solido, "vm": vm,
        "F_total": F_total, "A_bruta": A_bruta, "v_elem": v_elem,
        "eps_eff_solido": eps_eff[solido], "vm_solido": vm[solido],
    })

    # rigidez aparente
    uz = u[2::3]
    H = float(z.max() - z.min())
    eps_app = abs(float(np.mean(uz[techo]))) / max(H, 1e-300)
    out["eps_app"] = eps_app
    out["E_app"] = float(sigma0) / max(eps_app, 1e-300)
    pr(1.0, "Listo.")
    return out


def criterio_pistoia(res, frac=FRAC_CRITICA, eps_crit=EPS_CRITICA):
    """Carga de fallo aparente a partir de un ensayo lineal ya resuelto.

    Devuelve el factor por el que hay que escalar la carga para que la fraccion
    `frac` del volumen oseo alcance `eps_crit` de deformacion efectiva, y la
    tension y carga de fallo que resultan.

    Como el problema es LINEAL, eps escala con la carga y el factor es una
    division. Toda la maquinaria no lineal que uno esperaria aqui sobra: esa es
    justamente la gracia del criterio.
    """
    if not res.get("ok"):
        return {"ok": False, "msg": res.get("msg", "El ensayo no se resolvio.")}

    e = np.asarray(res["eps_eff_solido"], float)
    # `vol_solido` (opcional, alineado con eps_eff_solido/vm_solido): volumen
    # de cada elemento. Con voxeles no hace falta; con tetraedros de tamanos
    # distintos «el 2 % del volumen» ya no es «el 2 % de los elementos».
    w = res.get("vol_solido")
    w = None if w is None else np.asarray(w, float).ravel()
    fin = np.isfinite(e)
    e = e[fin]
    we = None if w is None else w[fin]
    if e.size == 0:
        return {"ok": False, "msg": "No hay tejido oseo con deformacion."}

    # Percentil (1-frac): el valor que supera exactamente esa fraccion del
    # tejido. Con elementos iguales el percentil por conteo YA esta ponderado
    # por volumen; si no, lo pondera `percentil_ponderado` (mismo tipo 7).
    umbral = float(percentil_ponderado(e, we, 100.0 * (1.0 - frac),
                                       metodo="linear"))
    if umbral <= 0:
        return {"ok": False, "msg": "Deformacion nula en el tejido."}

    k = eps_crit / umbral
    # von Mises del tejido, escalada a la carga de fallo: es la tension que
    # alcanzaria el material en el instante que el criterio declara la rotura.
    # Contrastarla con el limite elastico del tejido (~150-200 MPa en hueso
    # mineralizado) es la comprobacion de coherencia que el criterio no hace.
    vms = np.asarray(res.get("vm_solido", []), float)
    fv = np.isfinite(vms)
    vms = vms[fv]
    wv = None if (w is None or w.size != fv.size) else w[fv]
    vm_out = {}
    if vms.size:
        # NOTA DE CONVENCION: estas cuatro claves conservan `np.percentile` con
        # el defecto de numpy, que NO es el `prctile` de MATLAB que declara el
        # resto del proyecto. Se dejan como estaban a proposito, para no mover
        # numeros que ya estan en JSON exportados y en figuras; la diferencia
        # medida en los VOI de H4 es del 0.4-0.8 %. Las claves de superficie
        # que se anaden debajo SI usan la convencion declarada, porque son
        # nuevas y porque asi coinciden con `Estudio_Convergencia`, que es de
        # donde sale su validacion.
        p99 = percentil_ponderado(vms, wv, 99.0, metodo="linear")
        vm_out = {"vm_max": float(vms.max()), "vm_media": _media(vms, wv),
                  "vm_p99": float(p99),
                  "vm_max_fallo": float(k * vms.max()),
                  "vm_p99_fallo": float(k * p99)}

        # Capa superficial: el unico sabor del pico que converge. El maximo se
        # sigue reportando —esta en los resultados de siempre y hay quien lo
        # pide—, pero el numero citable es este.
        sup = estadisticos_vm(res)
        vm_out.update(sup)
        p99s = sup.get("vm_p99_superficie")
        if p99s is not None and np.isfinite(p99s):
            vm_out["vm_p99_superficie_fallo"] = float(k * p99s)

    return {
        "ok": True,
        "factor": float(k),
        **vm_out,
        "sigma_fallo": float(k * res["sigma_app"]),
        "F_fallo": float(k * res["F_total"]),
        "eps_eff_p": umbral,
        "eps_eff_max": float(e.max()),
        "eps_eff_media": _media(e, we),
        "ponderado_volumen": not _pesos_iguales(we),
        "frac": float(frac), "eps_crit": float(eps_crit),
        "n_elem_solido": int(e.size),
        "E_app": res.get("E_app"),
        # El criterio se calibro en radio distal humano: los dos parametros
        # viajan con el resultado para que nadie los de por universales.
        "nota": (f"Criterio de Pistoia con frac={frac:.1%} y "
                 f"eps_crit={eps_crit:.2%}; ambos son convenciones calibradas "
                 f"en radio distal humano, no constantes fisicas."),
    }


# ---------------------------------------------------------------------------
# Ensayo en un eje cualquiera, y en los tres
# ---------------------------------------------------------------------------

# Permutaciones CICLICAS que llevan el eje pedido a la posicion z. Se eligen
# ciclicas y no cualquier reordenacion porque det(P) = +1: la terna sigue siendo
# dextrogira, asi que las componentes de `eps` y `sig` del marco del ensayo
# conservan el signo de las tensiones tangenciales. Con el material isotropo de
# aqui daria igual, pero dejarlo bien cuesta lo mismo y no obliga a revisarlo si
# alguna vez se mete un tejido anisotropo.
_PERM = {0: (1, 2, 0), 1: (2, 0, 1), 2: (0, 1, 2)}
EJES = ("X", "Y", "Z")


def ensayo_compresion_eje(BW, spacing, eje=2, **kw):
    """`ensayo_compresion` a lo largo de X, Y o Z (`eje` = 0, 1 o 2).

    El ensayo esta escrito para comprimir en z —el filtro portante, la cara
    cargada, el apoyo y `E_app` miran todos a la tercera dimension—. En vez de
    duplicar esa logica con indices parametrizados, que es donde se cuelan los
    errores que no dan error, se PERMUTA la estructura para que el eje pedido
    ocupe la posicion z, se resuelve el mismo codigo de siempre y se deshace la
    permutacion en los campos 3D.

    QUE MARCO USA CADA SALIDA — importa al leerla:

      * `E_app`, `eps_app`, `sigma_fallo`, `frac_portante` y todo lo escalar son
        invariantes o ya estan referidos al eje ensayado: se leen tal cual.
      * `campo_eps_eff`, `campo_vm` y `campo_desp` se devuelven en el marco
        ORIGINAL, listos para pintar sobre la geometria sin tocar nada.
      * `eps`, `sig` y `u` se quedan en el marco DEL ENSAYO (donde la carga va
        en z). Deshacer esa rotacion componente a componente seria facil de
        equivocar y ninguna parte del visor los usa, asi que se marcan con
        `marco_permutado` en vez de convertirlos a medias.
    """
    eje = int(eje)
    if eje not in _PERM:
        raise ValueError("eje debe ser 0 (X), 1 (Y) o 2 (Z)")
    p = _PERM[eje]

    BW = np.asarray(BW, dtype=bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)

    if kw.get("modulo_rel") is not None:
        kw = dict(kw)
        kw["modulo_rel"] = np.transpose(np.asarray(kw["modulo_rel"]), p)
    out = ensayo_compresion(np.transpose(BW, p), spacing[list(p)], **kw)
    out["eje"] = eje
    out["eje_nombre"] = EJES[eje]

    inv = tuple(int(i) for i in np.argsort(p))
    for clave in ("campo_eps_eff", "campo_vm", "campo_desp"):
        if clave in out and out[clave] is not None:
            out[clave] = np.transpose(out[clave], inv)
    out["marco_permutado"] = (eje != 2)
    return out


def ensayo_triaxial(BW, spacing, ejes=(0, 1, 2), progreso=None, **kw):
    """Ensayo de compresion en los tres ejes. Devuelve un dict por eje.

    PARA QUE SIRVE: `E_app` de un solo eje no dice si la estructura es
    anisotropa; con los tres sale un cociente E_max/E_min que es una medida
    MECANICA de anisotropia, independiente del DA del tensor MIL —que es
    puramente geometrico— y comparable con Ex, Ey, Ez de la homogeneizacion
    periodica de `spinpy.elastic`.

    Que las dos rutas coincidan no es trivial: la homogeneizacion impone
    contorno PERIODICO sobre la celda y este ensayo impone platos reales sobre
    una probeta finita. Coinciden solo si la estructura es grande frente a la
    trabecula. Discrepar es informacion —efecto de tamano o de borde—, no un
    fallo; por eso conviene mirar las dos y no sustituir una por la otra.

    La clave `resumen` trae E_app por eje y la anisotropia, ya calculados.
    """
    res = {}
    ejes = tuple(int(e) for e in ejes)
    for i, e in enumerate(ejes):
        if progreso:
            progreso(i / max(len(ejes), 1), f"Comprimiendo en {EJES[e]}…")
        res[EJES[e]] = ensayo_compresion_eje(BW, spacing, eje=e, **kw)
    if progreso:
        progreso(1.0, "Listo.")

    E = {EJES[e]: res[EJES[e]].get("E_app") for e in ejes
         if res[EJES[e]].get("ok")}
    resumen = {"E_app": E, "ejes": [EJES[e] for e in ejes]}
    v = [x for x in E.values() if x and np.isfinite(x) and x > 0]
    if len(v) >= 2:
        resumen["anisotropia_E_app"] = float(max(v) / min(v))
        resumen["eje_mas_rigido"] = max(E, key=lambda k: E[k])
    res["resumen"] = resumen
    return res


# ---------------------------------------------------------------------------
# Convergencia de malla
# ---------------------------------------------------------------------------

def estudio_convergencia(BW, spacing, resoluciones=(20, 24, 28, 32, 40),
                         eje=2, remuestrear=None, progreso=None, **kw):
    """`E_app` a varias resoluciones, para decidir si el valor es citable.

    POR QUE HACE FALTA. En micro-elementos finitos sobre estructuras poco
    densas `E_app` NO converge sin más al refinar: medido sobre el VOI proximal
    de H4 (BV/TV 0.28) salieron 700, 797 y 365 MPa a 32^3, 48^3 y 64^3. Un
    valor tomado de una sola resolución no significa nada, y esa advertencia ya
    estaba escrita en la interfaz sin que hubiera forma de resolverla.

    QUE SIGNIFICA UNA SERIE NO MONOTONA. No es ruido numérico —el problema es
    determinista y el residuo del solver se comprueba en cada punto—: es que al
    refinar cambia la GEOMETRIA. Cada remuestreo redibuja qué trabéculas
    quedan conectadas, y una trabécula de uno o dos vóxeles puede aparecer,
    desaparecer o partirse. El modelo no está convergiendo hacia una estructura
    fija; está resolviendo estructuras distintas. Por eso se devuelve también
    `rho` y `frac_portante` de cada punto: si varían, la comparación de E_app
    entre resoluciones no es una comparación de discretización.

    Se devuelve además el orden de convergencia observado entre los tres
    últimos puntos y la extrapolación de Richardson, PERO solo cuando la serie
    es monótona y la densidad estable; en caso contrario se marcan como no
    aplicables en vez de dar un número que invitaría a citarlo.
    """
    if remuestrear is None:
        from .elastic import remuestrear_bw as remuestrear

    BW = np.asarray(BW, dtype=bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)

    puntos = []
    resoluciones = [int(n) for n in resoluciones]
    for i, n in enumerate(resoluciones):
        if progreso:
            progreso(i / len(resoluciones), f"Resolviendo a {n}³…")
        bw, spr = remuestrear(BW, spacing, n)
        r = ensayo_compresion_eje(bw, spr, eje=eje, **kw)
        p = {"n": n, "h": float(spr[0]), "rho": float(bw.mean()),
             "ok": bool(r.get("ok"))}
        if r.get("ok"):
            p.update({"E_app": float(r["E_app"]),
                      "n_elem": int(r["n_elem"]),
                      "frac_portante": float(r.get("frac_portante", 1.0)),
                      "residuo": float(r.get("residuo_rel", 0.0)),
                      "solver": r.get("solver", "")})
        else:
            p["msg"] = r.get("msg", "")
        puntos.append(p)
    if progreso:
        progreso(1.0, "Listo.")

    out = {"puntos": puntos, "eje": EJES[int(eje)]}
    val = [p for p in puntos if p["ok"]]
    if len(val) < 3:
        out["veredicto"] = ("Hacen falta al menos tres resoluciones resueltas "
                            "para hablar de convergencia.")
        return out

    E = np.array([p["E_app"] for p in val], float)
    h = np.array([p["h"] for p in val], float)
    rho = np.array([p["rho"] for p in val], float)

    out["dispersion_rel"] = float((E.max() - E.min()) / E.mean())
    out["deriva_rho_rel"] = float((rho.max() - rho.min()) / rho.mean())
    d = np.diff(E)
    monotona = bool(np.all(d > 0) or np.all(d < 0))
    out["monotona"] = monotona

    # SERIE PLANA. Un bloque macizo da el mismo E_app a toda resolucion salvo
    # el error de redondeo, y esas diferencias de 1e-11 pasan el test de
    # monotonia por casualidad: basta con que el ruido tenga el mismo signo.
    # Extrapolar de ahi da un "orden de convergencia" que es puro artefacto
    # del redondeo. Se corta antes: si no hay variacion, ya se convergio y no
    # hay nada que extrapolar.
    if out["dispersion_rel"] < 1e-6:
        out["salto_final_rel"] = float(abs(E[-1] - E[-2]) / abs(E[-1]))
        out["veredicto"] = (
            f"E_app no varia con la malla (dispersion "
            f"{out['dispersion_rel']:.1e}): ya esta convergido. Es lo que se "
            f"espera de una estructura que el elemento representa exactamente, "
            f"como un bloque macizo.")
        return out
    # El ultimo salto relativo es la señal practica: si refinar ya casi no
    # mueve el valor, el valor esta cerca de su limite.
    out["salto_final_rel"] = float(abs(E[-1] - E[-2]) / abs(E[-1]))

    if monotona and out["deriva_rho_rel"] < CONV_DERIVA_RHO_MAX:
        # Richardson con los tres ultimos puntos. Solo tiene sentido si la
        # serie es monotona: con una serie que sube y baja, el "orden" que sale
        # de la formula es un artefacto sin interpretacion.
        E1, E2, E3 = E[-3], E[-2], E[-1]
        h1, h2, h3 = h[-3], h[-2], h[-1]
        den = E2 - E1
        r21 = h1 / h2
        if abs(den) > 0 and abs((E3 - E2) / den) not in (0.0, 1.0):
            try:
                pord = np.log(abs(den / (E3 - E2))) / np.log(r21)
                out["orden"] = float(pord)
                out["E_extrapolado"] = float(
                    E3 + (E3 - E2) / ((h2 / h3) ** pord - 1.0))
                out["error_estimado_rel"] = float(
                    abs(out["E_extrapolado"] - E3) / abs(E3))
            except (ValueError, ZeroDivisionError, FloatingPointError):
                pass

    if not monotona and out["dispersion_rel"] <= CONV_OSCILACION_MAX:
        out["veredicto"] = (
            f"La serie no es monotona pero oscila sin tendencia dentro de un "
            f"{100*out['dispersion_rel']:.1f}% (banda {100*CONV_OSCILACION_MAX:.0f}"
            f"%). No se extrapola; cita el valor de la malla mas fina "
            f"declarando esa banda como incertidumbre de discretizacion.")
    elif not monotona:
        out["veredicto"] = (
            f"La serie NO es monotona (dispersion {100*out['dispersion_rel']:.0f}%). "
            f"No hay convergencia que extrapolar: al refinar cambia la "
            f"geometria conectada, no solo la discretizacion. NO cites un "
            f"E_app de una sola resolucion.")
    elif out.get("error_estimado_rel", 1.0) < CONV_ERROR_EXTRAPOLADO_MAX:
        out["veredicto"] = (
            f"Serie monotona y el ultimo refinado mueve el valor un "
            f"{100*out['salto_final_rel']:.1f}%. El extrapolado es citable "
            f"declarando la malla y el error estimado.")
    else:
        out["veredicto"] = (
            f"Serie monotona pero aun lejos del limite (el ultimo refinado "
            f"mueve el valor un {100*out['salto_final_rel']:.1f}%). Hace falta "
            f"refinar mas antes de citar nada.")
    return out
