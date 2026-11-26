"""
cargas.py - Analisis mecanico con cargas y apoyos ELEGIDOS por el usuario.

POR QUE EXISTE
--------------
`resistencia.ensayo_compresion` resuelve UN caso: compresion axial entre dos
platos. Es el ensayo de referencia, pero no el unico que interesa: una carga
oblicua, una cizalla, una flexion o una fuerza aplicada en un punto de la
superficie (un contacto articular, un tornillo) son preguntas distintas y la
respuesta cambia. Aqui el caso de carga es un DATO: una lista de cargas y una
lista de apoyos, cada uno con su region, su direccion y su magnitud. La
ventana `DialogoCargas` (`dialogo_cargas.py`) los construye con el raton; este
modulo no sabe nada de Qt y se puede usar desde un script o un cuaderno.

SE APLICA A CUALQUIER ESTRUCTURA
--------------------------------
El solver recibe una mascara booleana y un `spacing`. Da igual que sea el VOI
del hueso, un espinodoide o un dual-lattice: es la misma malla de hexaedros
trilineales (solo los voxeles de hueso son elementos) y el mismo material, asi
que los resultados de las tres estructuras se pueden comparar.

NO CAMBIA NADA DE LO EXISTENTE
------------------------------
Reutiliza `hex8_ke`, el ensamblado por bloques y la matriz B de
`resistencia`/`elastic`, pero no las modifica. Con una presion sobre la cara +z,
un apoyo deslizante en -z y la estabilizacion automatica, este modulo reproduce
`ensayo_compresion` (el bloque 37 lo comprueba a precision de maquina).

UNIDADES
--------
Longitudes en mm, fuerzas en N, momentos en N*mm, presiones y tensiones en MPa
y modulo de Young en PASCALES (`E_s`, como el resto de la app). Internamente el
sistema es Pa-mm-(Pa*mm^2), el mismo que explica `ensayo_compresion`: una fuerza
en newtons se multiplica por 1e6. Las claves de salida que son tensiones
terminan en `_MPa`.

EL CASO DE CARGA
----------------
Es un diccionario serializable (JSON), que es lo que guarda la sesion:

    caso = {
      "material": {"E_s": 20e9, "nu": 0.3},
      "cargas":  [carga, ...],
      "apoyos":  [apoyo, ...],
      "estabilizar": True,        # ver mas abajo
      "solo_portante": True,
    }

REGION. Donde actua una carga o un apoyo. Siempre en el marco de la estructura
(x, y, z = ejes de la mascara, en mm desde la esquina del volumen):

    {"forma": "cara", "cara": "+z"}                  una cara del cubo
    {"forma": "esfera", "centro": [x,y,z], "radio": r}      un parche
    {"forma": "caja", "lo": [x,y,z], "hi": [x,y,z]}
    {"forma": "cilindro", "centro": [...], "radio": r, "eje": 2, "altura": h}
    {"forma": "nodos", "ids": [id, ...]}     nodos sueltos (carga o apoyo puntual)

y, opcionalmente, `"ambito"`: "superficie" (por omision) o "volumen", y
`"filtro_normal": {"dir": [x,y,z], "angulo": 60}`.

  * superficie: la region son las CARAS DE VOXEL EXPUESTAS (hueso a un lado,
    vacio o el exterior del cubo al otro) cuyo centro cae dentro. Una carga se
    reparte por area, como una presion: es el equivalente nodal consistente. Es
    lo que pide pinchar con el raton sobre la superficie del hueso.
  * volumen: todos los nodos de hueso dentro de la region. Pensado para
    apoyos y desplazamientos impuestos; una fuerza se reparte por nodo.
  * filtro_normal: solo las caras cuya normal saliente forma menos de `angulo`
    grados con `dir`. Sobre una trabecula, el parche de una esfera incluye
    caras de los dos lados; el filtro deja las que dan al punto pinchado.

CARGAS (`"tipo"`):
  "fuerza"         `magnitud` (N) y `direccion`. Es la fuerza TOTAL de la region.
  "presion"        `magnitud` (MPa), positiva = comprime (empuja hacia dentro,
                   normal entrante). `area_ref`: "bruta" (la seccion del cubo
                   que cubre la cara, como la tension aparente de la app) o
                   "expuesta" (el area de las caras de hueso seleccionadas).
  "momento"        `magnitud` (N*mm) alrededor de `eje` (vector), por `centro`
                   (por omision, el baricentro ponderado de la region). Es un
                   PAR: fuerza neta nula.
  "desplazamiento" `u` = [ux, uy, uz] en mm con None en el componente libre, o
                   `magnitud` y `direccion`. Es un dato cinematico, no una
                   fuerza; la reaccion sale del calculo.

`direccion` es un vector, uno de "+x" "-x" "+y" "-y" "+z" "-z", o
"normal_entrante" / "normal_saliente" (la normal de CADA cara de la region).

APOYOS (`"tipo"`):
  "empotrado"      los tres desplazamientos nulos.
  "deslizante"     solo el normal a la cara (o el de `gdl`): la region desliza
                   sin friccion sobre un plato rigido.
  "personalizado"  los grados de libertad de `gdl`, por ejemplo ["x", "z"].

ESTABILIZAR. Si los apoyos no impiden los seis movimientos de solido rigido
(caso tipico: compresion con apoyo deslizante, o cargas autoequilibradas sin
apoyo alguno), la matriz es singular. Con `estabilizar` se anaden las
restricciones MINIMAS (seis grados de libertad elegidos para que sean
independientes) y despues se comprueba que sus reacciones son nulas, que es lo
que ocurre si y solo si las cargas estan equilibradas. Si no lo estan, el
resultado dice que falta un apoyo en lugar de entregar un campo falso.

LA CARGA DEBE LLEGAR A UN APOYO
-------------------------------
Un fragmento de hueso que lleva carga pero no toca ningun apoyo es un
mecanismo, y uno que toca un apoyo pero no lleva carga no transmite nada. Se
trabaja solo con los fragmentos (conexion por caras) que tienen carga y
apoyo; si una carga cae sobre un fragmento sin apoyo se rechaza el caso con un
mensaje, no se descarta en silencio.

PISTOIA Y CARGAS LOCALES
------------------------
El factor de fallo que se devuelve usa el criterio de Pistoia (2 % del volumen
oseo por encima de 0.7 % de deformacion efectiva). Esos dos numeros se
calibraron con compresion de radio distal; con una fuerza concentrada en un
punto la cola de deformacion la dicta el contacto y el criterio no esta
validado. El resultado lo declara en `nota`.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, cg, splu

from .elastic import UMBRAL_DIRECTO, _ensamblar, hex8_ke
from .resistencia import (EPS_CRITICA, ESQUINAS, FRAC_CRITICA, _b_centro,
                          _matriz_D, _modos_rigidos, percentil_ponderado)

E_S_DEF, NU_S_DEF = 20e9, 0.30
N_A_SISTEMA = 1e6          # newtons -> Pa*mm^2 (ver cabecera)

CARAS = {"+x": (0, 1), "-x": (0, -1), "+y": (1, 1), "-y": (1, -1),
         "+z": (2, 1), "-z": (2, -1)}
EJES = ("x", "y", "z")

TIPOS_CARGA = ("fuerza", "presion", "momento", "desplazamiento")
TIPOS_APOYO = ("empotrado", "deslizante", "personalizado")
FORMAS = ("cara", "esfera", "caja", "cilindro", "nodos")
AMBITOS = ("superficie", "volumen")

CAMPOS = ("vm", "eps_eff", "desp", "ux", "uy", "uz", "s1", "s3", "tau_max",
          "sxx", "syy", "szz")
ETIQUETA_CAMPO = {
    "vm": "von Mises (MPa)", "eps_eff": "Deformacion efectiva",
    "desp": "Desplazamiento total (mm)", "ux": "Desplazamiento x (mm)",
    "uy": "Desplazamiento y (mm)", "uz": "Desplazamiento z (mm)",
    "s1": "Tension principal maxima (MPa)",
    "s3": "Tension principal minima (MPa)",
    "tau_max": "Cortante maximo (MPa)", "sxx": "Tension xx (MPa)",
    "syy": "Tension yy (MPa)", "szz": "Tension zz (MPa)"}


class ErrorCaso(ValueError):
    """El caso de carga esta mal formado o no se puede resolver."""


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------

def _spacing3(spacing):
    sp = np.atleast_1d(np.asarray(spacing, float)).ravel()
    return np.repeat(sp, 3) if sp.size == 1 else sp[:3]


def extension(shape, spacing):
    """Tamano del volumen en mm (x, y, z)."""
    return np.asarray(shape, float) * _spacing3(spacing)


def _vec(v, nombre="vector"):
    a = np.asarray(v, float).ravel()
    if a.size != 3 or not np.all(np.isfinite(a)):
        raise ErrorCaso(f"{nombre}: se esperaban tres numeros finitos")
    return a


def _unitario(v, nombre="direccion"):
    a = _vec(v, nombre)
    n = float(np.linalg.norm(a))
    if n < 1e-12:
        raise ErrorCaso(f"{nombre}: el vector es nulo")
    return a / n


def direccion_vector(d):
    """Vector unitario de una direccion dada como vector o como '+x'... '-z'.

    Las direcciones `normal_entrante` y `normal_saliente` dependen de cada
    cara y las resuelve el reparto, no esta funcion: devuelve None.
    """
    if isinstance(d, str):
        k = d.strip().lower()
        if k in ("normal_entrante", "normal_saliente"):
            return None
        if k in CARAS:
            ax, s = CARAS[k]
            v = np.zeros(3)
            v[ax] = s
            return v
        raise ErrorCaso(f"direccion desconocida: {d!r}")
    return _unitario(d)


# ---------------------------------------------------------------------------
# Caras de voxel expuestas
# ---------------------------------------------------------------------------

def caras_expuestas(S, spacing):
    """Caras de voxel con hueso a un lado y vacio (o el exterior) al otro.

    Devuelve un diccionario de arrays alineados, uno por cara:
      cent (m,3)  centro, en mm
      nor  (m,3)  normal SALIENTE del hueso
      area (m,)   mm^2
      nod  (m,4)  ids de los cuatro nodos en la rejilla completa
      vox  (m,3)  indices del voxel de hueso al que pertenece
    El id de nodo es `a + nnx*b + nnx*nny*c`, el mismo que usa el ensayo.
    """
    S = np.asarray(S, bool)
    sp = _spacing3(spacing)
    nn = np.array(S.shape) + 1
    cent, nor, area, nod, vox = [], [], [], [], []
    for ax in range(3):
        a1, a2 = (ax + 1) % 3, (ax + 2) % 3
        P = np.pad(S, [(1, 1) if a == ax else (0, 0) for a in range(3)])
        lo = np.take(P, range(0, S.shape[ax] + 1), axis=ax)
        hi = np.take(P, range(1, S.shape[ax] + 2), axis=ax)
        for signo, mask in ((1, lo & ~hi), (-1, hi & ~lo)):
            idx = np.nonzero(mask)
            if idx[0].size == 0:
                continue
            p = np.stack(idx, axis=1).astype(np.int64)     # (m,3) con p[ax]=plano
            v = p.copy()
            if signo == 1:
                v[:, ax] -= 1
            e1 = np.zeros(3, np.int64); e1[a1] = 1
            e2 = np.zeros(3, np.int64); e2[a2] = 1
            esq = [p, p + e1, p + e1 + e2, p + e2]
            nod.append(np.stack(
                [q[:, 0] + nn[0] * q[:, 1] + nn[0] * nn[1] * q[:, 2]
                 for q in esq], axis=1))
            cent.append((p + 0.5 * (e1 + e2)) * sp)
            n = np.zeros((p.shape[0], 3)); n[:, ax] = signo
            nor.append(n)
            area.append(np.full(p.shape[0], sp[a1] * sp[a2]))
            vox.append(v)
    if not cent:
        z3 = np.zeros((0, 3))
        return {"cent": z3, "nor": z3, "area": np.zeros(0),
                "nod": np.zeros((0, 4), np.int64),
                "vox": np.zeros((0, 3), np.int64)}
    return {"cent": np.concatenate(cent), "nor": np.concatenate(nor),
            "area": np.concatenate(area), "nod": np.concatenate(nod),
            "vox": np.concatenate(vox)}


def _coord_nodos(ids, shape, spacing):
    nn = np.array(shape) + 1
    c, r = np.divmod(ids, nn[0] * nn[1])
    b, a = np.divmod(r, nn[0])
    return np.stack([a, b, c], axis=1) * _spacing3(spacing)


def _dentro(region, X, ext, sp, tol):
    """Mascara de los puntos X (n,3) que caen en la geometria de `region`."""
    forma = region.get("forma")
    if forma == "cara":
        cara = str(region.get("cara", "")).strip()
        if cara not in CARAS:
            raise ErrorCaso(f"cara desconocida: {cara!r}")
        ax, s = CARAS[cara]
        plano = ext[ax] if s > 0 else 0.0
        capa = float(region.get("capa", 0))
        return np.abs(X[:, ax] - plano) <= capa * sp[ax] + tol
    if forma == "esfera":
        c = _vec(region["centro"], "centro")
        return np.linalg.norm(X - c, axis=1) <= float(region["radio"]) + tol
    if forma == "caja":
        lo, hi = _vec(region["lo"], "lo"), _vec(region["hi"], "hi")
        return np.all((X >= np.minimum(lo, hi) - tol)
                      & (X <= np.maximum(lo, hi) + tol), axis=1)
    if forma == "cilindro":
        c = _vec(region["centro"], "centro")
        ax = int(region.get("eje", 2))
        if ax not in (0, 1, 2):
            raise ErrorCaso("eje del cilindro: 0, 1 o 2")
        d = X - c
        axial = np.abs(d[:, ax])
        rad = np.linalg.norm(np.delete(d, ax, axis=1), axis=1)
        return ((rad <= float(region["radio"]) + tol)
                & (axial <= 0.5 * float(region["altura"]) + tol))
    raise ErrorCaso(f"forma de region desconocida: {forma!r}")


def seleccionar(region, caras, shape, spacing, nodos_hueso=None):
    """Caras (o nodos) de la region. Ver la cabecera para la semantica.

    Devuelve {"quads": indices o None, "nodos": ids unicos, "peso": peso de
    cada nodo (area tributaria en mm^2, o 1 por nodo en el ambito volumen),
    "area": area total seleccionada}. `nodos_hueso` son los ids de los nodos
    que pertenecen a algun voxel de hueso; solo se usa en el ambito volumen.
    """
    sp = _spacing3(spacing)
    ext = extension(shape, spacing)
    tol = 1e-9 * float(ext.max())
    ambito = str(region.get("ambito", "superficie")).lower()
    if ambito not in AMBITOS:
        raise ErrorCaso(f"ambito desconocido: {ambito!r}")

    if region.get("forma") == "nodos":
        # Nodos elegidos uno a uno (carga o apoyo puntual exacto). Solo
        # cuentan los que pertenecen a hueso; el resto se descarta.
        if nodos_hueso is None:
            raise ErrorCaso("region de nodos sin nodos de hueso")
        ids = np.intersect1d(np.asarray(region.get("ids", []), np.int64),
                             np.asarray(nodos_hueso, np.int64))
        return {"quads": None, "nodos": ids, "peso": np.ones(ids.size),
                "area": float(ids.size)}

    if ambito == "volumen":
        if nodos_hueso is None:
            raise ErrorCaso("ambito volumen sin nodos de hueso")
        X = _coord_nodos(nodos_hueso, shape, spacing)
        m = _dentro(region, X, ext, sp, tol)
        ids = np.asarray(nodos_hueso)[m]
        return {"quads": None, "nodos": ids, "peso": np.ones(ids.size),
                "area": float(ids.size)}

    if region.get("forma") == "cara":
        ax, s = CARAS[str(region.get("cara", "")).strip()]
        # En una cara solo cuentan las caras de voxel de ESA cara del cubo,
        # no las del hueso que miran en el mismo sentido en otro plano.
        plano = ext[ax] if s > 0 else 0.0
        m = ((caras["nor"][:, ax] * s > 0.5)
             & (np.abs(caras["cent"][:, ax] - plano) <= tol))
    else:
        m = _dentro(region, caras["cent"], ext, sp, tol)
    f = region.get("filtro_normal")
    if f:
        d = _unitario(f["dir"], "filtro_normal.dir")
        ang = float(f.get("angulo", 60.0))
        m = m & (caras["nor"] @ d >= np.cos(np.radians(ang)) - 1e-12)
    q = np.nonzero(m)[0]
    nods = caras["nod"][q]
    ids, inv = np.unique(nods.ravel(), return_inverse=True)
    peso = np.zeros(ids.size)
    np.add.at(peso, inv, np.repeat(caras["area"][q] / 4.0, 4))
    return {"quads": q, "nodos": ids, "peso": peso,
            "area": float(caras["area"][q].sum())}


def area_bruta_cara(region, shape, spacing):
    """Seccion del cubo cubierta por una cara; None si la region no es una cara."""
    if region.get("forma") != "cara":
        return None
    ax, _ = CARAS[str(region.get("cara", "")).strip()]
    ext = extension(shape, spacing)
    return float(np.prod(np.delete(ext, ax)))


# ---------------------------------------------------------------------------
# Validacion y plantillas
# ---------------------------------------------------------------------------

def _gdl(g):
    if g is None:
        return None
    out = []
    for c in g:
        k = str(c).strip().lower()
        if k not in EJES:
            raise ErrorCaso(f"gdl desconocido: {c!r} (x, y o z)")
        if EJES.index(k) not in out:
            out.append(EJES.index(k))
    return out


def validar_region(region, etiqueta="region"):
    """Lanza `ErrorCaso` si la region esta mal formada (forma, cara, medidas)."""
    if not isinstance(region, dict):
        raise ErrorCaso(f"{etiqueta}: debe ser un diccionario")
    f = region.get("forma")
    if f not in FORMAS:
        raise ErrorCaso(f"{etiqueta}: forma desconocida {f!r}")
    if str(region.get("ambito", "superficie")).lower() not in AMBITOS:
        raise ErrorCaso(f"{etiqueta}: ambito desconocido")
    try:
        if f == "cara":
            if str(region.get("cara", "")).strip() not in CARAS:
                raise ErrorCaso(f"{etiqueta}: cara desconocida "
                                f"{region.get('cara')!r} (+x ... -z)")
            if int(region.get("capa", 0)) < 0:
                raise ErrorCaso(f"{etiqueta}: capas negativas")
        elif f == "esfera":
            _vec(region["centro"], "centro")
            if not float(region["radio"]) > 0:
                raise ErrorCaso(f"{etiqueta}: el radio debe ser positivo")
        elif f == "caja":
            _vec(region["lo"], "lo")
            _vec(region["hi"], "hi")
        elif f == "cilindro":
            _vec(region["centro"], "centro")
            if int(region.get("eje", 2)) not in (0, 1, 2):
                raise ErrorCaso(f"{etiqueta}: el eje es 0, 1 o 2")
            if not (float(region["radio"]) > 0 and float(region["altura"]) > 0):
                raise ErrorCaso(f"{etiqueta}: radio y altura positivos")
        else:
            np.asarray(region.get("ids", []), np.int64)
        fn = region.get("filtro_normal")
        if fn:
            _unitario(fn["dir"], "filtro_normal.dir")
    except KeyError as e:
        raise ErrorCaso(f"{etiqueta}: falta {e}") from e
    except (TypeError, ValueError) as e:
        if isinstance(e, ErrorCaso):
            raise
        raise ErrorCaso(f"{etiqueta}: {e}") from e


def validar_caso(caso):
    """Lanza `ErrorCaso` si el caso esta mal formado. Devuelve el caso."""
    if not isinstance(caso, dict):
        raise ErrorCaso("el caso de carga debe ser un diccionario")
    mat = caso.get("material", {})
    if not float(mat.get("E_s", E_S_DEF)) > 0:
        raise ErrorCaso("E_s debe ser positivo")
    nu = float(mat.get("nu", NU_S_DEF))
    if not -1 < nu < 0.5:
        raise ErrorCaso("nu debe estar en (-1, 0.5)")
    cargas, apoyos = caso.get("cargas", []), caso.get("apoyos", [])
    if not cargas:
        raise ErrorCaso("no hay ninguna carga: anade al menos una")
    for i, c in enumerate(cargas):
        t = c.get("tipo")
        if t not in TIPOS_CARGA:
            raise ErrorCaso(f"carga {i + 1}: tipo desconocido {t!r}")
        if "region" not in c:
            raise ErrorCaso(f"carga {i + 1}: falta la region")
        validar_region(c["region"], f"carga {i + 1}")
        if t in ("fuerza", "presion", "momento"):
            m = float(c.get("magnitud", 0.0))
            if not np.isfinite(m):
                raise ErrorCaso(f"carga {i + 1}: magnitud no finita")
        if t == "fuerza":
            d = c.get("direccion", "normal_entrante")
            if direccion_vector(d) is None and not str(d).startswith("normal"):
                raise ErrorCaso(f"carga {i + 1}: direccion invalida")
        if t == "momento" and "eje" not in c:
            raise ErrorCaso(f"carga {i + 1}: el momento necesita `eje`")
        if t == "desplazamiento" and "u" not in c and "magnitud" not in c:
            raise ErrorCaso(f"carga {i + 1}: desplazamiento sin `u`")
    for i, a in enumerate(apoyos):
        t = a.get("tipo")
        if t not in TIPOS_APOYO:
            raise ErrorCaso(f"apoyo {i + 1}: tipo desconocido {t!r}")
        if "region" not in a:
            raise ErrorCaso(f"apoyo {i + 1}: falta la region")
        validar_region(a["region"], f"apoyo {i + 1}")
        _gdl(a.get("gdl"))
        if t == "personalizado" and not a.get("gdl"):
            raise ErrorCaso(f"apoyo {i + 1}: `personalizado` necesita `gdl`")
    return caso


def caso_vacio():
    return {"material": {"E_s": E_S_DEF, "nu": NU_S_DEF}, "cargas": [],
            "apoyos": [], "estabilizar": True, "solo_portante": True}


PLANTILLAS = ("compresion", "traccion", "cizalla", "flexion", "torsion",
              "puntual")


def plantilla(nombre, shape, spacing, eje=2, sigma_MPa=1.0, F_N=100.0,
              desp_mm=None):
    """Caso de carga de partida para un ensayo clasico, sobre el cubo entero.

    `eje` es el eje del ensayo (0, 1, 2). Las plantillas son el PUNTO DE
    PARTIDA que la ventana ofrece: se editan o se sustituyen. Compresion y
    traccion reproducen el ensayo de la app (presion sobre la seccion bruta,
    apoyo deslizante); cizalla, flexion y torsion actuan sobre la cara +eje con
    la cara -eje empotrada.
    """
    ext = extension(shape, spacing)
    eje = int(eje)
    mas, menos = "+" + EJES[eje], "-" + EJES[eje]
    caso = caso_vacio()
    a = (eje + 1) % 3
    if nombre in ("compresion", "traccion"):
        s = abs(float(sigma_MPa)) * (1 if nombre == "compresion" else -1)
        caso["cargas"] = [{"tipo": "presion", "magnitud": s,
                           "region": {"forma": "cara", "cara": mas},
                           "area_ref": "bruta"}]
        caso["apoyos"] = [{"tipo": "deslizante",
                           "region": {"forma": "cara", "cara": menos}}]
    elif nombre == "cizalla":
        d = np.zeros(3); d[a] = 1
        caso["cargas"] = [{"tipo": "fuerza", "magnitud": float(F_N),
                           "direccion": d.tolist(),
                           "region": {"forma": "cara", "cara": mas}}]
        caso["apoyos"] = [{"tipo": "empotrado",
                           "region": {"forma": "cara", "cara": menos}}]
    elif nombre == "flexion":
        eje_m = np.zeros(3); eje_m[a] = 1
        L = float(np.min(np.delete(ext, eje)))
        caso["cargas"] = [{"tipo": "momento", "magnitud": float(F_N) * L,
                           "eje": eje_m.tolist(),
                           "region": {"forma": "cara", "cara": mas}}]
        caso["apoyos"] = [{"tipo": "empotrado",
                           "region": {"forma": "cara", "cara": menos}}]
    elif nombre == "torsion":
        eje_m = np.zeros(3); eje_m[eje] = 1
        L = float(np.min(np.delete(ext, eje)))
        caso["cargas"] = [{"tipo": "momento", "magnitud": float(F_N) * L,
                           "eje": eje_m.tolist(),
                           "region": {"forma": "cara", "cara": mas}}]
        caso["apoyos"] = [{"tipo": "empotrado",
                           "region": {"forma": "cara", "cara": menos}}]
    elif nombre == "puntual":
        c = ext / 2.0
        c[eje] = ext[eje]
        caso["cargas"] = [{
            "tipo": "fuerza", "magnitud": float(F_N),
            "direccion": "normal_entrante",
            "region": {"forma": "esfera", "centro": c.tolist(),
                       "radio": 0.15 * float(np.min(ext)),
                       "filtro_normal": {"dir": _unit_eje(eje, 1),
                                         "angulo": 60.0}}}]
        caso["apoyos"] = [{"tipo": "empotrado",
                           "region": {"forma": "cara", "cara": menos}}]
    else:
        raise ErrorCaso(f"plantilla desconocida: {nombre!r}")
    if desp_mm is not None and nombre in ("compresion", "traccion", "cizalla"):
        u = [None, None, None]
        u[eje if nombre != "cizalla" else a] = (
            (-abs(desp_mm) if nombre == "compresion" else abs(desp_mm)))
        caso["cargas"] = [{"tipo": "desplazamiento", "u": u,
                           "region": {"forma": "cara", "cara": mas}}]
    return caso


def _unit_eje(eje, s):
    v = [0.0, 0.0, 0.0]
    v[eje] = float(s)
    return v


# ---------------------------------------------------------------------------
# Ensamblado de cargas y restricciones
# ---------------------------------------------------------------------------

def _fuerza_nodal(carga, sel, caras, shape, spacing):
    """Fuerzas nodales (dict id -> vector, en N) de una carga de fuerza o presion."""
    t = carga["tipo"]
    d = carga.get("direccion", "normal_entrante" if t == "presion" else None)
    if d is None:
        raise ErrorCaso("la fuerza necesita direccion")
    if carga["region"].get("forma") == "nodos" and (
            t == "presion" or direccion_vector(d) is None):
        raise ErrorCaso("una presion o una direccion normal necesitan caras "
                        "de hueso; sobre nodos sueltos usa una fuerza con "
                        "un vector")
    if t == "presion":
        M = float(carga["magnitud"])
        ref = str(carga.get("area_ref", "bruta")).lower()
        A = area_bruta_cara(carga["region"], shape, spacing) \
            if ref == "bruta" else None
        if A is None:
            A = sel["area"]
        total = M * A                       # MPa * mm^2 = N
    else:
        total = float(carga["magnitud"])
    if sel["area"] <= 0 or sel["nodos"].size == 0:
        raise ErrorCaso("la region de la carga no toca ninguna cara de "
                        "hueso: muevela o agrandala")
    v = direccion_vector(d)
    ids, w = sel["nodos"], sel["peso"]
    if v is not None:
        F = np.outer(w / w.sum(), v * total)
        return ids, F
    # normal de cada cara: la fuerza de cada cara sigue SU normal
    q = sel["quads"]
    if q is None:
        raise ErrorCaso("la direccion normal solo existe en el ambito "
                        "superficie")
    signo = -1.0 if str(d).strip().lower() == "normal_entrante" else 1.0
    nor = signo * caras["nor"][q]
    fq = (caras["area"][q] / sel["area"] * total)[:, None] * nor / 4.0
    pos = np.searchsorted(ids, caras["nod"][q].ravel())
    F = np.zeros((ids.size, 3))
    np.add.at(F, pos, np.repeat(fq, 4, axis=0))
    return ids, F


def _par_nodal(carga, sel, shape, spacing):
    """Fuerzas nodales de un par M alrededor de un eje: fuerza neta nula."""
    if sel["nodos"].size == 0:
        raise ErrorCaso("la region del momento no toca ningun nodo de hueso")
    e = _unitario(carga["eje"], "eje")
    ids, w = sel["nodos"], sel["peso"]
    X = _coord_nodos(ids, shape, spacing)
    c = (_vec(carga["centro"], "centro") if carga.get("centro") is not None
         else (X * w[:, None]).sum(0) / w.sum())
    r = X - c
    rp = r - np.outer(r @ e, e)
    den = float((w * (rp * rp).sum(1)).sum())
    if den <= 0:
        raise ErrorCaso("el momento actua sobre un solo punto del eje: "
                        "no hay brazo de palanca")
    F = (np.cross(e, r) * w[:, None]) * (float(carga["magnitud"]) / den)
    return ids, F


def _u_impuesto(carga):
    if carga.get("u") is not None:
        u = list(carga["u"])
        if len(u) != 3:
            raise ErrorCaso("u debe tener tres componentes (None = libre)")
        return [None if x is None else float(x) for x in u]
    d = carga.get("direccion", "+z")
    v = direccion_vector(d)
    if v is None:
        raise ErrorCaso("un desplazamiento impuesto no puede ser normal")
    m = float(carga["magnitud"])
    if np.count_nonzero(np.abs(v) > 1e-12) == 1:
        return [m * v[k] if abs(v[k]) > 1e-12 else None for k in range(3)]
    return [m * v[k] for k in range(3)]


def _gdl_apoyo(apoyo, sel, caras):
    g = _gdl(apoyo.get("gdl"))
    if g:
        return g
    t = apoyo["tipo"]
    if t == "empotrado":
        return [0, 1, 2]
    # deslizante: el eje normal a la cara; en otra region, el de la normal media
    r = apoyo["region"]
    if r.get("forma") == "cara":
        return [CARAS[str(r["cara"]).strip()][0]]
    q = sel["quads"]
    if q is not None and len(q):
        return [int(np.argmax(np.abs(caras["nor"][q]).mean(0)))]
    return [2]


# ---------------------------------------------------------------------------
# Estabilizacion
# ---------------------------------------------------------------------------

def _dofs_estabilizadores(coord, fijos_loc):
    """Grados de libertad MINIMOS que completan los seis modos de solido rigido.

    `coord` (n,3) son los nodos de UN cuerpo; `fijos_loc` los grados de libertad
    (indices locales 3*nodo+k) que ya estan restringidos. Devuelve los que
    faltan, elegidos entre nodos repartidos por todo el cuerpo (extremos y
    muestra) con un QR con pivote, para que los seis queden bien
    condicionados. Vacio si ya estaba determinado.
    """
    n = coord.shape[0]
    L = float(np.ptp(coord, axis=0).max()) or 1.0
    R = _modos_rigidos(coord - coord.mean(0))
    R[:, 3:] /= L
    if len(fijos_loc):
        U, s, Vt = np.linalg.svd(R[np.asarray(fijos_loc)], full_matrices=True)
        rango = int((s > 1e-8 * s[0]).sum()) if s.size else 0
        N = Vt[rango:].T
    else:
        N = np.eye(6)
    if N.shape[1] == 0:
        return []
    cand = set(int(i) for i in np.argmin(coord, axis=0)) \
        | set(int(i) for i in np.argmax(coord, axis=0))
    rng = np.random.default_rng(0)
    cand |= set(int(i) for i in rng.choice(n, size=min(n, 80), replace=False))
    cand = np.array(sorted(cand))
    dofs = (3 * cand[:, None] + np.arange(3)[None, :]).ravel()
    M = R[dofs] @ N
    from scipy.linalg import qr
    _, _, piv = qr(M.T, pivoting=True)
    elegidos = [int(dofs[p]) for p in piv[:N.shape[1]]]
    chk = np.linalg.svd(R[np.concatenate([np.asarray(fijos_loc, int),
                                          elegidos]).astype(int)],
                        compute_uv=False)
    if chk.size < 6 or chk[-1] < 1e-8 * chk[0]:
        raise ErrorCaso("no se pudo estabilizar la pieza (geometria "
                        "degenerada)")
    return elegidos


# ---------------------------------------------------------------------------
# Resolucion
# ---------------------------------------------------------------------------

def _componentes_caras(BW):
    from scipy import ndimage
    return ndimage.label(BW, structure=ndimage.generate_binary_structure(3, 1))


def resolver_caso(BW, spacing, caso, tol=1e-8, progreso=None,
                  campos=CAMPOS):
    """Resuelve el caso de carga sobre la estructura `BW` (elasticidad lineal).

    Devuelve un diccionario con `ok` y `msg`; con `ok` verdadero trae:

      campos      {nombre: array 3D float32, NaN donde no hay hueso portante}
      stats       por campo: max, min, media y p99 (tensiones en MPa)
      reacciones  una por apoyo o desplazamiento impuesto: fuerza (N) y su
                  modulo, el area y el numero de nodos
      equilibrio  residuo relativo de fuerzas y de momentos (deberia ser ~1e-9
                  con solver directo y ~`tol` con el iterativo)
      cargas      la fuerza total y el par total que aplica cada carga
      rigidez     desplazamiento medio de las regiones cargadas y N/mm
      pistoia     factor de fallo, con su advertencia sobre cargas locales
      avisos      lista de textos; si no esta vacia, hay algo que leer
      u_malla     (nx+1, ny+1, nz+1, 3) float32, para dibujar la deformada

    Un caso mal formado lanza `ErrorCaso`; uno bien formado que no se puede
    resolver (carga sin apoyo, fuerzas no equilibradas, solver sin converger)
    devuelve `ok=False` y un `msg` que dice que cambiar.
    """
    def pr(f, m):
        if progreso:
            progreso(f, m)

    BW = np.asarray(BW, bool)
    if BW.ndim != 3:
        raise ErrorCaso("la estructura debe ser una mascara 3D")
    caso = validar_caso(caso)
    sp = _spacing3(spacing)
    dx, dy, dz = sp
    shape = BW.shape
    nx, ny, nz = shape
    nn = np.array(shape) + 1
    E_s = float(caso.get("material", {}).get("E_s", E_S_DEF))
    nu_s = float(caso.get("material", {}).get("nu", NU_S_DEF))
    out = {"ok": False, "msg": "", "avisos": [], "rho": float(BW.mean()),
           "E_s": E_s, "nu_s": nu_s}
    if not BW.any():
        out["msg"] = "La estructura no contiene material solido."
        return out

    pr(0.03, "Seleccionando las regiones…")
    caras = caras_expuestas(BW, sp)
    nodo_id = lambda a, b, c: a + nn[0] * b + nn[0] * nn[1] * c
    i_, j_, k_ = np.nonzero(BW)
    conn_g = np.stack([nodo_id(i_ + c[0], j_ + c[1], k_ + c[2])
                       for c in ESQUINAS], axis=1)
    nodos_hueso = np.unique(conn_g)
    lab, ncomp = _componentes_caras(BW)
    lab_vox = lab[i_, j_, k_]

    def comps_de(sel):
        if sel["quads"] is not None:
            v = caras["vox"][sel["quads"]]
            return set(np.unique(lab[v[:, 0], v[:, 1], v[:, 2]]).tolist())
        m = np.isin(conn_g, sel["nodos"]).any(axis=1)
        return set(np.unique(lab_vox[m]).tolist())

    # --- cargas -> fuerzas nodales y restricciones ---------------------------
    fuerzas = []        # (ids, F(n,3) en N, indice de carga)
    impuestos = []      # (ids, u3, indice)
    sel_cargas = []
    for n, c in enumerate(caso["cargas"]):
        sel = seleccionar(c["region"], caras, shape, sp, nodos_hueso)
        if sel["nodos"].size == 0:
            out["msg"] = (f"La region de la carga {n + 1} no toca ninguna "
                          f"cara de hueso.")
            return out
        sel_cargas.append(sel)
        if c["tipo"] in ("fuerza", "presion"):
            ids, F = _fuerza_nodal(c, sel, caras, shape, sp)
            fuerzas.append((ids, F, n))
        elif c["tipo"] == "momento":
            ids, F = _par_nodal(c, sel, shape, sp)
            fuerzas.append((ids, F, n))
        else:
            impuestos.append((sel["nodos"], _u_impuesto(c), n))

    sel_apoyos, gdl_apoyos = [], []
    for n, a in enumerate(caso["apoyos"]):
        sel = seleccionar(a["region"], caras, shape, sp, nodos_hueso)
        if sel["nodos"].size == 0:
            out["msg"] = (f"La region del apoyo {n + 1} no toca ninguna cara "
                          f"de hueso.")
            return out
        sel_apoyos.append(sel)
        gdl_apoyos.append(_gdl_apoyo(a, sel, caras))

    # --- que fragmentos entran -----------------------------------------------
    con_carga, con_apoyo = set(), set()
    for sel in sel_cargas:
        con_carga |= comps_de(sel)
    for sel in sel_apoyos:
        con_apoyo |= comps_de(sel)
    for ids, u3, n in impuestos:
        con_apoyo |= comps_de(sel_cargas[n])
    estabilizar = bool(caso.get("estabilizar", True))
    solo_portante = bool(caso.get("solo_portante", True))
    if solo_portante:
        if estabilizar:
            entran = set(con_carga) | set(con_apoyo & con_carga)
        else:
            entran = con_carga & con_apoyo
        sin_apoyo = con_carga - con_apoyo
        if sin_apoyo and not estabilizar:
            out["msg"] = ("Hay carga sobre un fragmento de hueso que no toca "
                          "ningun apoyo (mecanismo). Mueve el apoyo o la "
                          "carga, o activa la estabilizacion si las cargas "
                          "estan equilibradas.")
            return out
    else:
        entran = set(range(1, ncomp + 1))
    if not entran:
        out["msg"] = ("Ninguna carga llega a un apoyo: no hay fragmento de "
                      "hueso que una las regiones.")
        return out
    elem_ok = np.isin(lab_vox, list(entran))
    out["frac_portante"] = float(elem_ok.sum() / BW.sum())
    if out["frac_portante"] < 0.999:
        out["avisos"].append(
            "El %.1f %% del hueso no une cargas y apoyos y no se resuelve."
            % (100 * (1 - out["frac_portante"])))

    pr(0.08, "Ensamblando…")
    conn_ok = conn_g[elem_ok]
    usados, conn = np.unique(conn_ok, return_inverse=True)
    conn = conn.reshape(-1, 8)
    n_nod, nelem = usados.size, conn.shape[0]
    ndof = 3 * n_nod
    coord = _coord_nodos(usados, shape, sp)
    out.update({"n_elem": int(nelem), "n_dof": int(ndof)})
    ijk_elem = np.stack([i_[elem_ok], j_[elem_ok], k_[elem_ok]], axis=1)

    def locales(ids):
        """Posicion compacta de los ids de rejilla que sobreviven."""
        pos = np.searchsorted(usados, ids)
        pos[pos >= usados.size] = 0
        return pos, usados[pos] == ids

    # fuerzas
    F = np.zeros(ndof)
    perdida = []
    desglose = []
    for ids, Fn, n in fuerzas:
        pos, ok = locales(ids)
        if not ok.all():
            perdida.append(n)
        for k in range(3):
            np.add.at(F, 3 * pos[ok] + k, Fn[ok, k] * N_A_SISTEMA)
        sel = sel_cargas[n]
        X = _coord_nodos(ids, shape, sp)
        desglose.append({
            "indice": n, "tipo": caso["cargas"][n]["tipo"],
            "F_N": Fn.sum(0).tolist(),
            "M_Nmm": np.cross(X, Fn).sum(0).tolist(),
            "area_mm2": sel["area"], "n_nodos": int(ids.size)})
    if perdida:
        out["msg"] = ("La carga %s cae en parte sobre hueso que no se "
                      "resuelve; acerca el apoyo." % ", ".join(
                          str(p + 1) for p in perdida))
        return out

    # restricciones: apoyos primero, desplazamientos impuestos encima
    fijo = {}                      # dof -> (valor, etiqueta)
    for n, (sel, g) in enumerate(zip(sel_apoyos, gdl_apoyos)):
        pos, ok = locales(sel["nodos"])
        for k in g:
            for d in 3 * pos[ok] + k:
                fijo[int(d)] = (0.0, ("apoyo", n))
    for ids, u3, n in impuestos:
        pos, ok = locales(ids)
        for k in range(3):
            if u3[k] is None:
                continue
            for d in 3 * pos[ok] + k:
                fijo[int(d)] = (float(u3[k]), ("desp", n))

    # cuerpos (por nodos compartidos) y estabilizadores
    estab = []
    cuerpo = np.zeros(n_nod, int)
    cuerpos_libres = set()         # cuerpos sin ninguna restriccion propia
    if nelem:
        from scipy.sparse.csgraph import connected_components
        filas = np.repeat(conn, 8, axis=1).ravel()
        cols = np.tile(conn, (1, 8)).ravel()
        G = sparse.coo_matrix((np.ones(filas.size, np.int8), (filas, cols)),
                              shape=(n_nod, n_nod))
        nb, cuerpo = connected_components(G, directed=False)
        for b in range(nb):
            nodos_b = np.nonzero(cuerpo == b)[0]
            loc = {int(g): i for i, g in enumerate(nodos_b)}
            fl = [3 * loc[d // 3] + d % 3 for d in fijo if d // 3 in loc]
            falta = _dofs_estabilizadores(coord[nodos_b], fl)
            if not fl:
                cuerpos_libres.add(b)
            for d in falta:
                estab.append(int(3 * nodos_b[d // 3] + d % 3))
    if estab and not estabilizar:
        out["msg"] = ("Los apoyos no impiden el movimiento de solido rigido "
                      "de la pieza. Anade un apoyo, o activa la "
                      "estabilizacion automatica.")
        return out
    for d in estab:
        fijo.setdefault(d, (0.0, ("estab", 0)))
    out["n_estabilizadores"] = len(estab)

    fijos = np.array(sorted(fijo), dtype=np.int64)
    uf = np.array([fijo[int(d)][0] for d in fijos])
    libres = np.setdiff1d(np.arange(ndof), fijos)

    edof = np.empty((nelem, 24), dtype=np.int64)
    for a in range(8):
        for k in range(3):
            edof[:, 3 * a + k] = 3 * conn[:, a] + k
    keS = hex8_ke(dx, dy, dz, E_s, nu_s)
    K = _ensamblar(keS, edof, np.ones(nelem), ndof)
    K = (K + K.T) / 2

    pr(0.35, "Resolviendo %d grados de libertad…" % libres.size)
    u = np.zeros(ndof)
    u[fijos] = uf
    Kff = K[libres][:, libres].tocsr()
    rhs = F[libres] - K[libres][:, fijos] @ uf
    resuelto = False
    if libres.size <= UMBRAL_DIRECTO:
        try:
            u[libres] = splu(Kff.tocsc()).solve(rhs)
            out["solver"] = "LU directo"
            resuelto = bool(np.all(np.isfinite(u)))
        except Exception:
            resuelto = False
    if not resuelto:
        try:
            import pyamg
            B = _modos_rigidos(coord)[libres, :]
            ml = pyamg.smoothed_aggregation_solver(Kff, B=B, max_coarse=500)
            u[libres] = ml.solve(rhs, tol=tol, maxiter=500, accel="cg")
            out["solver"] = "AMG (con modos rigidos) + CG"
        except ImportError:
            d = np.asarray(Kff.diagonal(), float)
            d[d <= 0] = np.finfo(float).eps
            M = LinearOperator(Kff.shape, lambda v, d=d: v / d)
            u[libres], _ = cg(Kff, rhs, rtol=tol, maxiter=20000, M=M)
            out["solver"] = "CG + Jacobi"
        nb_ = np.linalg.norm(rhs)
        res = np.linalg.norm(rhs - Kff @ u[libres]) / (nb_ if nb_ > 0 else 1.0)
        out["residuo_rel"] = float(res)
        if not np.isfinite(res) or res > max(1e-6, 100 * tol):
            out["msg"] = ("El solver no convergio (residuo %.1e). El campo "
                          "no es fiable." % res)
            return out

    # --- reacciones y equilibrio -----------------------------------------------
    pr(0.8, "Reacciones y equilibrio…")
    Rint = K @ u - F                       # reaccion en los grados fijos
    R = np.zeros(ndof)
    R[fijos] = Rint[fijos]
    Ftot = np.abs(F).sum() / N_A_SISTEMA
    if estab:
        re = np.abs(Rint[np.array(estab)]).sum() / N_A_SISTEMA
        ref = max(Ftot, np.abs(R).sum() / N_A_SISTEMA, 1e-300)
        out["reaccion_estabilizadores_rel"] = float(re / ref)
        if re / ref > 1e-3:
            sueltas = []
            for ids, _Fn, n in fuerzas:
                pos, ok = locales(ids)
                if set(cuerpo[pos[ok]].tolist()) & cuerpos_libres:
                    sueltas.append("C%d" % (n + 1))
            if sueltas:
                out["msg"] = (
                    "La carga %s actua sobre un fragmento de hueso que no "
                    "toca ningun apoyo: no se puede equilibrar. Mueve el "
                    "apoyo o la carga." % ", ".join(sueltas))
            else:
                out["msg"] = (
                    "Las cargas no estan equilibradas con los apoyos "
                    "definidos (los restos de reaccion son el %.1f %% de la "
                    "carga): falta un apoyo, o un par que compense."
                    % (100 * re / ref))
            return out
    reac = []
    for etq, lista in (("apoyo", caso["apoyos"]), ("desp", caso["cargas"])):
        for n in range(len(lista)):
            if etq == "desp" and lista[n]["tipo"] != "desplazamiento":
                continue
            ds = [d for d, (_, e) in fijo.items() if e == (etq, n)]
            if not ds:
                continue
            Fv = np.zeros(3)
            Mv = np.zeros(3)
            ds = np.array(ds)
            nodo, k = ds // 3, ds % 3
            np.add.at(Fv, k, R[ds] / N_A_SISTEMA)
            r = np.cross(coord[nodo], np.eye(3)[k]) * (
                R[ds] / N_A_SISTEMA)[:, None]
            Mv = r.sum(0)
            sel = (sel_apoyos[n] if etq == "apoyo" else sel_cargas[n])
            reac.append({"origen": etq, "indice": n, "F_N": Fv.tolist(),
                         "modulo_N": float(np.linalg.norm(Fv)),
                         "M_Nmm": Mv.tolist(), "area_mm2": sel["area"],
                         "n_nodos": int(sel["nodos"].size)})
    out["reacciones"] = reac
    out["cargas"] = desglose
    Fap = np.zeros(3)
    Map = np.zeros(3)
    for d in desglose:
        Fap += d["F_N"]
        Map += d["M_Nmm"]
    nodo_all, k_all = np.divmod(fijos, 3)
    Rt = np.zeros(3)
    np.add.at(Rt, k_all, R[fijos] / N_A_SISTEMA)
    Mt = (np.cross(coord[nodo_all], np.eye(3)[k_all])
          * (R[fijos] / N_A_SISTEMA)[:, None]).sum(0)
    esc_F = max(np.linalg.norm(Fap), np.linalg.norm(Rt), 1e-300)
    esc_M = max(np.linalg.norm(Map), np.linalg.norm(Mt),
                esc_F * float(np.ptp(coord, axis=0).max()) * 1e-3, 1e-300)
    out["equilibrio"] = {
        "F_aplicada_N": Fap.tolist(), "F_reaccion_N": Rt.tolist(),
        "residuo_F": float(np.linalg.norm(Fap + Rt) / esc_F),
        "residuo_M": float(np.linalg.norm(Map + Mt) / esc_M)}

    # --- campos -----------------------------------------------------------------
    pr(0.88, "Recuperando tensiones…")
    B0 = _b_centro(dx, dy, dz)
    D = _matriz_D(E_s, nu_s)
    eps = u[edof] @ B0.T
    sig = (eps @ D.T)                      # Pa
    U = 0.5 * np.einsum("ij,ij->i", sig, eps)
    eps_eff = np.sqrt(np.maximum(2.0 * U / E_s, 0.0))
    s = sig / 1e6                          # MPa
    vm = np.sqrt(np.maximum(
        0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2
               + (s[:, 2] - s[:, 0]) ** 2)
        + 3.0 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2), 0.0))
    T = np.empty((nelem, 3, 3))
    T[:, 0, 0], T[:, 1, 1], T[:, 2, 2] = s[:, 0], s[:, 1], s[:, 2]
    T[:, 1, 2] = T[:, 2, 1] = s[:, 3]
    T[:, 0, 2] = T[:, 2, 0] = s[:, 4]
    T[:, 0, 1] = T[:, 1, 0] = s[:, 5]
    if any(c in campos for c in ("s1", "s3", "tau_max")):
        ev = np.linalg.eigvalsh(T)         # ascendente
    else:
        ev = np.zeros((nelem, 3))
    un = u.reshape(-1, 3)
    por_elem = lambda c: un[conn, c].mean(axis=1)
    valores = {
        "vm": vm, "eps_eff": eps_eff,
        "desp": np.linalg.norm(un, axis=1)[conn].mean(axis=1),
        "ux": por_elem(0), "uy": por_elem(1), "uz": por_elem(2),
        "s1": ev[:, 2], "s3": ev[:, 0], "tau_max": 0.5 * (ev[:, 2] - ev[:, 0]),
        "sxx": s[:, 0], "syy": s[:, 1], "szz": s[:, 2]}
    out["campos"], out["stats"] = {}, {}
    for nombre in campos:
        v = valores[nombre]
        c = np.full(shape, np.nan, np.float32)
        c[ijk_elem[:, 0], ijk_elem[:, 1], ijk_elem[:, 2]] = v
        out["campos"][nombre] = c
        out["stats"][nombre] = {
            "max": float(v.max()), "min": float(v.min()),
            "media": float(v.mean()),
            "p99": float(percentil_ponderado(v, None, 99.0, metodo="hazen"))}
    out["desp_max_nodal"] = float(np.linalg.norm(un, axis=1).max())

    ug = np.zeros((int(nn[0]), int(nn[1]), int(nn[2]), 3), np.float32)
    a_, b_, c_ = np.unravel_index(usados, tuple(nn), order="F")
    ug[a_, b_, c_] = un
    out["u_malla"] = ug
    out["elem_ijk"] = ijk_elem
    out["spacing"] = sp

    # rigidez de las regiones cargadas
    rig = []
    for ids, Fn, n in fuerzas:
        pos, ok = locales(ids)
        um = un[pos[ok]].mean(0)
        Ft = Fn.sum(0)
        nf = float(np.linalg.norm(Ft))
        if nf > 0:
            dirf = Ft / nf
            ud = float(um @ dirf)
            rig.append({"indice": n, "u_medio_mm": um.tolist(),
                        "u_en_direccion_mm": ud,
                        "rigidez_N_mm": (nf / abs(ud)) if abs(ud) > 0
                        else float("inf")})
        else:
            rig.append({"indice": n, "u_medio_mm": um.tolist()})
    for ids, u3, n in impuestos:
        sel_r = [r for r in reac if r["origen"] == "desp"
                 and r["indice"] == n]
        if not sel_r:
            continue
        dsp = max([abs(x) for x in u3 if x is not None] or [0.0])
        rig.append({"indice": n, "desp_impuesto_mm": dsp,
                    "reaccion_N": sel_r[0]["modulo_N"],
                    "rigidez_N_mm": (sel_r[0]["modulo_N"] / dsp)
                    if dsp > 0 else float("inf")})
    out["rigidez"] = rig

    # Pistoia, con su advertencia
    out["pistoia"] = _pistoia(eps_eff, caso)
    out["interno"] = {
        "sig_MPa": s, "eps_eff": eps_eff,
        "u_medio": {str(n): un[locales(ids)[0][locales(ids)[1]]].mean(0)
                    for ids, _F, n in fuerzas}}
    out["ok"] = True
    pr(1.0, "Listo.")
    return out


def _pistoia(eps_eff, caso):
    """Factor de escala de todas las cargas hasta el criterio de Pistoia."""
    e = np.asarray(eps_eff, float)
    e = e[np.isfinite(e)]
    umbral = float(percentil_ponderado(e, None, 100 * (1 - FRAC_CRITICA),
                                       metodo="linear")) if e.size else 0.0
    if umbral <= 0:
        return None
    k = EPS_CRITICA / umbral
    esc = []
    for i, c in enumerate(caso["cargas"]):
        if c["tipo"] == "desplazamiento":
            u = c.get("u")
            m = (max(abs(x) for x in u if x is not None)
                 if u is not None and any(x is not None for x in u)
                 else abs(float(c.get("magnitud", 0.0))))
            esc.append({"indice": i, "valor": m, "fallo": k * m, "unidad": "mm"})
        else:
            m = abs(float(c.get("magnitud", 0.0)))
            uni = {"fuerza": "N", "presion": "MPa", "momento": "N*mm"}[c["tipo"]]
            esc.append({"indice": i, "valor": m, "fallo": k * m, "unidad": uni})
    return {
        "factor": float(k), "eps_eff_p": umbral, "cargas_fallo": esc,
        "frac": FRAC_CRITICA, "eps_crit": EPS_CRITICA,
        "nota": ("Criterio de Pistoia (2 % del volumen oseo sobre 0.7 % de "
                 "deformacion efectiva), calibrado con compresion de radio "
                 "distal. Con cargas concentradas o desplazamientos "
                 "impuestos no esta validado: tomalo como un indice "
                 "relativo, no como una carga de rotura.")}


def desplazar_puntos(res, puntos, escala=1.0):
    """Puntos (n,3) en mm desplazados por el campo de la solucion (interp. trilineal).

    Sirve para dibujar la deformada de la superficie del hueso. `escala`
    amplifica el desplazamiento; la deformada real es invisible a escala 1.
    """
    ug = res["u_malla"]
    sp = res["spacing"]
    P = np.asarray(puntos, float)
    g = P / sp
    hi = np.array(ug.shape[:3]) - 1
    i0 = np.clip(np.floor(g).astype(int), 0, hi - 1)
    f = np.clip(g - i0, 0.0, 1.0)
    out = np.zeros_like(P)
    for dx_ in (0, 1):
        for dy_ in (0, 1):
            for dz_ in (0, 1):
                w = (np.where(dx_, f[:, 0], 1 - f[:, 0])
                     * np.where(dy_, f[:, 1], 1 - f[:, 1])
                     * np.where(dz_, f[:, 2], 1 - f[:, 2]))
                out += w[:, None] * ug[i0[:, 0] + dx_, i0[:, 1] + dy_,
                                       i0[:, 2] + dz_]
    return P + escala * out


# ---------------------------------------------------------------------------
# Ayudas para elegir con el raton
# ---------------------------------------------------------------------------
#
# La ventana pincha un punto de la escena 3D; estas funciones lo convierten en
# una region. Viven aqui, y no en la ventana, para poder comprobarlas sin
# interfaz y porque son las que deciden QUE se carga.

def cara_del_cubo(punto, shape, spacing):
    """La cara del cubo (+x ... -z) mas cercana al punto, y su distancia en mm."""
    ext = extension(shape, spacing)
    p = np.asarray(punto, float)
    mejor, dmin = None, np.inf
    for nombre, (ax, s) in CARAS.items():
        d = abs(p[ax] - (ext[ax] if s > 0 else 0.0))
        if d < dmin:
            mejor, dmin = nombre, float(d)
    return mejor, dmin


def cara_expuesta_cercana(caras, punto):
    """Indice de la cara de voxel expuesta mas cercana al punto, o None."""
    if caras["cent"].shape[0] == 0:
        return None
    d = np.linalg.norm(caras["cent"] - np.asarray(punto, float), axis=1)
    return int(np.argmin(d))


def normal_en(caras, punto):
    """Normal saliente de la cara de hueso mas cercana al punto pinchado."""
    q = cara_expuesta_cercana(caras, punto)
    return None if q is None else caras["nor"][q].copy()


def nodo_cercano(caras, punto, shape, spacing):
    """Id (en la rejilla) del nodo de superficie mas cercano al punto."""
    q = cara_expuesta_cercana(caras, punto)
    if q is None:
        return None
    ids = caras["nod"][q]
    X = _coord_nodos(ids, shape, spacing)
    return int(ids[int(np.argmin(np.linalg.norm(X - np.asarray(punto, float),
                                                axis=1)))])


def region_parche(punto, radio, caras, angulo=60.0):
    """Region esferica sobre la superficie del hueso, centrada en el pinchazo.

    El centro se lleva al de la cara de hueso mas cercana (el punto pinchado
    cae sobre la isosuperficie suavizada, que no coincide exactamente con las
    caras de voxel), y se filtra por la normal de esa cara para que el parche
    no recoja las caras del otro lado de una trabecula.
    """
    q = cara_expuesta_cercana(caras, punto)
    if q is None:
        raise ErrorCaso("la estructura no tiene superficie expuesta")
    reg = {"forma": "esfera", "centro": caras["cent"][q].tolist(),
           "radio": float(radio)}
    if angulo is not None and angulo < 180:
        reg["filtro_normal"] = {"dir": caras["nor"][q].tolist(),
                                "angulo": float(angulo)}
    return reg


# ---------------------------------------------------------------------------
# Barrido de direccion
# ---------------------------------------------------------------------------

def barrido_direccion(BW, spacing, caso, indice=0, plano=(0, 2), n=24,
                      progreso=None):
    """Respuesta de la estructura al girar UNA carga dentro de un plano.

    Gira la carga `indice` (fuerza, presion o el eje de un momento) en el plano
    de los ejes `plano`, `n` pasos de 360/n grados, y devuelve la rigidez, el
    von Mises y el factor de Pistoia de cada direccion: la «rosa» mecanica de
    la estructura.

    POR QUE CUESTA DOS SOLUCIONES Y NO n. El problema es lineal y el reparto
    de la carga sobre la region no depende de su direccion, asi que la
    respuesta a una carga de direccion cos(t) e1 + sin(t) e2 es la combinacion
    cos(t) R1 + sin(t) R2 de las respuestas a e1 y e2. Se resuelven esas dos
    y los n angulos salen de combinar sus tensiones (el bloque 37 lo comprueba
    contra una solucion directa por angulo). La no linealidad de la
    deformacion efectiva y de von Mises se evalua sobre la tension combinada,
    no sobre las respuestas sueltas.

    Los apoyos NO giran: si el apoyo es deslizante, una componente tangencial
    de la carga no tiene reaccion y esa direccion falla. Para barrer en un
    plano completo hace falta un apoyo empotrado.

    Devuelve un diccionario con `ok`, los angulos en grados y, por angulo,
    `rigidez_N_mm` (fuerza / desplazamiento medio de la region en la direccion
    de la fuerza; solo para fuerza y presion), `vm_max`, `vm_p99`, `factor`
    (Pistoia) y `u_dir_mm`.
    """
    caso = validar_caso(caso)
    cargas = caso["cargas"]
    if not 0 <= indice < len(cargas):
        raise ErrorCaso("indice de carga fuera de rango")
    c0 = cargas[indice]
    if c0["tipo"] not in ("fuerza", "presion", "momento"):
        raise ErrorCaso("solo se puede girar una fuerza, una presion o un "
                        "momento")
    a, b = int(plano[0]), int(plano[1])
    if a == b or a not in (0, 1, 2) or b not in (0, 1, 2):
        raise ErrorCaso("el plano del barrido necesita dos ejes distintos")
    e = np.eye(3)
    clave = "eje" if c0["tipo"] == "momento" else "direccion"

    def con_dir(v):
        cc = {**caso, "cargas": [dict(c) for c in cargas]}
        cc["cargas"][indice][clave] = [float(x) for x in v]
        return cc

    base = []
    for k, v in enumerate((e[a], e[b])):
        if progreso:
            progreso(0.5 * k, "Resolviendo la direccion %d de 2…" % (k + 1))
        r = resolver_caso(BW, spacing, con_dir(v), campos=())
        if not r["ok"]:
            return {"ok": False, "msg": r["msg"]}
        base.append(r)
    E_s, nu_s = base[0]["E_s"], base[0]["nu_s"]
    Dinv = np.linalg.inv(_matriz_D(E_s, nu_s))
    s1 = base[0]["interno"]["sig_MPa"]
    s2 = base[1]["interno"]["sig_MPa"]
    u1 = base[0]["interno"]["u_medio"].get(str(indice))
    u2 = base[1]["interno"]["u_medio"].get(str(indice))
    mag = abs(float(c0.get("magnitud", 0.0)))
    th = np.linspace(0.0, 360.0, int(n), endpoint=False)
    fila = {k: np.zeros(th.size) for k in
            ("rigidez_N_mm", "u_dir_mm", "vm_max", "vm_p99", "factor")}
    for i, t in enumerate(np.radians(th)):
        c, s_ = np.cos(t), np.sin(t)
        s = c * s1 + s_ * s2                                  # MPa
        Pa = s * 1e6
        U = 0.5 * np.einsum("ij,ij->i", Pa, Pa @ Dinv.T)
        ee = np.sqrt(np.maximum(2.0 * U / E_s, 0.0))
        vm = np.sqrt(np.maximum(
            0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2
                   + (s[:, 2] - s[:, 0]) ** 2)
            + 3.0 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2), 0.0))
        fila["vm_max"][i] = vm.max()
        fila["vm_p99"][i] = percentil_ponderado(vm, None, 99.0,
                                                metodo="hazen")
        umb = float(percentil_ponderado(ee, None, 100 * (1 - FRAC_CRITICA),
                                        metodo="linear"))
        fila["factor"][i] = EPS_CRITICA / umb if umb > 0 else np.nan
        if c0["tipo"] != "momento" and u1 is not None:
            um = c * u1 + s_ * u2
            d = c * e[a] + s_ * e[b]
            ud = float(um @ d)
            fila["u_dir_mm"][i] = ud
            fila["rigidez_N_mm"][i] = (
                mag / abs(ud) if abs(ud) > 1e-300 else np.inf)
    if progreso:
        progreso(1.0, "Listo.")
    return {"ok": True, "msg": "", "theta_deg": th, "plano": [a, b],
            "indice": int(indice), "tipo": c0["tipo"], "magnitud": mag,
            **fila, "n_elem": base[0]["n_elem"],
            "solver": base[0].get("solver", "")}


# ---------------------------------------------------------------------------
# Texto, resumen y archivos
# ---------------------------------------------------------------------------

def _t(texto):
    """Traduce un texto de las listas de la interfaz (import tardio: `idioma`
    solo importa Qt dentro de una funcion, pero asi el nucleo no depende de el
    al importarse)."""
    from .idioma import _
    return _(texto)


def _num(x, esp=".4g"):
    return format(float(x), esp)


def describir_region(region):
    """Frase corta de una region, para las listas de la ventana."""
    f = region.get("forma")
    amb = ("" if str(region.get("ambito", "superficie")) == "superficie"
           else _t(", volumen"))
    if f == "cara":
        capa = int(region.get("capa", 0))
        t = _t("cara %s") % region.get("cara")
        if capa:
            t += _t(", %d capas") % capa
        return t + amb
    if f == "esfera":
        c = region["centro"]
        t = _t("parche de radio %s mm en (%s, %s, %s)") % (
            _num(region["radio"], ".3g"), *(_num(x, ".3g") for x in c))
        if region.get("filtro_normal"):
            t += _t(", solo caras a %g° de la normal") % float(
                region["filtro_normal"].get("angulo", 60))
        return t + amb
    if f == "caja":
        return _t("caja") + amb
    if f == "cilindro":
        return _t("cilindro de radio %s mm, eje %s") % (
            _num(region["radio"], ".3g"), EJES[int(region.get("eje", 2))]) + amb
    if f == "nodos":
        return _t("%d nodo(s)") % len(region.get("ids", []))
    return str(f)


def describir_direccion(d):
    if isinstance(d, str):
        return {"normal_entrante": _t("normal entrante"),
                "normal_saliente": _t("normal saliente")}.get(d, d)
    v = np.asarray(d, float)
    n = np.linalg.norm(v)
    return "(%s, %s, %s)" % tuple(_num(x, ".3g") for x in (v / n if n else v))


def describir_carga(c, n=None):
    t = c["tipo"]
    pre = ("C%d: " % (n + 1)) if n is not None else ""
    r = describir_region(c["region"])
    if t == "fuerza":
        return pre + _t("fuerza %s N hacia %s sobre %s") % (
            _num(c["magnitud"]), describir_direccion(
                c.get("direccion", "normal_entrante")), r)
    if t == "presion":
        return pre + _t("presion %s MPa (%s) sobre %s") % (
            _num(c["magnitud"]), c.get("area_ref", "bruta"), r)
    if t == "momento":
        return pre + _t("momento %s N·mm alrededor de %s sobre %s") % (
            _num(c["magnitud"]), describir_direccion(c["eje"]), r)
    u = c.get("u")
    if u is None:
        d = c.get("direccion", "+z")
        return pre + _t("desplazamiento %s mm hacia %s sobre %s") % (
            _num(c["magnitud"]), describir_direccion(d), r)
    return pre + _t("desplazamiento (%s) mm sobre %s") % (
        ", ".join(_t("libre") if x is None else _num(x) for x in u), r)


def describir_apoyo(a, n=None):
    pre = ("A%d: " % (n + 1)) if n is not None else ""
    g = _gdl(a.get("gdl"))
    t = a["tipo"]
    extra = (" (%s)" % ",".join(EJES[k] for k in (g or []))
             if (t == "personalizado" or g) else "")
    nombre = {"empotrado": _t("empotrado"), "deslizante": _t("deslizante"),
              "personalizado": _t("personalizado")}[t]
    return pre + _t("%s%s sobre %s") % (nombre, extra,
                                        describir_region(a["region"]))


def _json_seguro(o):
    if isinstance(o, dict):
        return {str(k): _json_seguro(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_json_seguro(v) for v in o]
    if isinstance(o, np.ndarray):
        return _json_seguro(o.tolist())
    if isinstance(o, (np.floating,)):
        v = float(o)
        return v if np.isfinite(v) else None
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


FORMATO_CASO = "spinpy.cargas/caso"


def guardar_caso(caso, ruta, estructura=None):
    """Escribe el caso en JSON. `estructura` (opcional) anota sobre que se uso."""
    import json
    validar_caso(caso)
    doc = {"formato": FORMATO_CASO, "version": 1, "caso": _json_seguro(caso)}
    if estructura:
        doc["estructura"] = _json_seguro(estructura)
    with open(ruta, "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=2, ensure_ascii=False)
    return ruta


def cargar_caso(ruta):
    """Lee un caso de `guardar_caso`, validandolo. Lanza `ErrorCaso` si no sirve."""
    import json
    try:
        with open(ruta, encoding="utf-8") as f:
            doc = json.load(f)
    except (OSError, ValueError) as e:
        raise ErrorCaso("no se pudo leer el caso: %s" % e) from e
    if not isinstance(doc, dict) or doc.get("formato") != FORMATO_CASO:
        raise ErrorCaso("el archivo no es un caso de carga de spinpy")
    return validar_caso(doc["caso"])


def resumen(res, caso, estructura=None):
    """Lo que se conserva de una solucion: todo menos los arrays grandes.

    Es lo que va a la sesion, al JSON exportado y al informe. Un numero
    sin su caso de carga no se puede reproducir, asi que viajan juntos.
    """
    from . import __version__
    if not res.get("ok"):
        return {"ok": False, "msg": res.get("msg", ""), "caso": caso}
    claves = ("n_elem", "n_dof", "solver", "residuo_rel", "frac_portante",
              "n_estabilizadores", "avisos", "stats", "reacciones", "cargas",
              "equilibrio", "rigidez", "pistoia", "rho", "E_s", "nu_s",
              "desp_max_nodal")
    out = {k: res[k] for k in claves if k in res}
    out.update({"ok": True, "caso": caso, "spinpy": __version__})
    if estructura:
        out["estructura"] = estructura
    return _json_seguro(out)


def exportar_json(res, caso, ruta, estructura=None):
    import json
    with open(ruta, "w", encoding="utf-8") as f:
        json.dump(resumen(res, caso, estructura), f, indent=2,
                  ensure_ascii=False)
    return ruta


def exportar_csv(res, caso, ruta):
    """Tabla de cargas, reacciones y estadisticos en un CSV legible."""
    import csv
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["seccion", "nombre", "Fx_N", "Fy_N", "Fz_N", "modulo_N",
                    "Mx_Nmm", "My_Nmm", "Mz_Nmm", "area_mm2", "n_nodos"])
        for d in res.get("cargas", []):
            c = caso["cargas"][d["indice"]]
            w.writerow(["carga", describir_carga(c, d["indice"]),
                        *d["F_N"], float(np.linalg.norm(d["F_N"])),
                        *d["M_Nmm"], d["area_mm2"], d["n_nodos"]])
        for r in res.get("reacciones", []):
            lista = caso["apoyos"] if r["origen"] == "apoyo" else caso["cargas"]
            txt = (describir_apoyo(lista[r["indice"]], r["indice"])
                   if r["origen"] == "apoyo"
                   else describir_carga(lista[r["indice"]], r["indice"]))
            w.writerow(["reaccion", txt, *r["F_N"], r["modulo_N"],
                        *r["M_Nmm"], r["area_mm2"], r["n_nodos"]])
        w.writerow([])
        w.writerow(["campo", "max", "min", "media", "p99"])
        for k, v in res.get("stats", {}).items():
            w.writerow([ETIQUETA_CAMPO.get(k, k), v["max"], v["min"],
                        v["media"], v["p99"]])
        eq = res.get("equilibrio") or {}
        w.writerow([])
        w.writerow(["residuo_fuerzas", eq.get("residuo_F")])
        w.writerow(["residuo_momentos", eq.get("residuo_M")])
    return ruta


def exportar_vtk(res, ruta):
    """Campos y desplazamientos en un `.vti` (ImageData) que abre ParaView.

    Los campos van como datos de CELDA (un valor por voxel, NaN donde no hay
    hueso resuelto) y el desplazamiento nodal como dato de PUNTO.
    """
    import pyvista as pv
    ug = res["u_malla"]
    n = np.array(ug.shape[:3]) - 1
    g = pv.ImageData(dimensions=tuple(int(x) + 1 for x in n),
                     spacing=tuple(float(x) for x in res["spacing"]),
                     origin=(0.0, 0.0, 0.0))
    for k, c in res["campos"].items():
        g.cell_data[k] = np.asarray(c, np.float32).flatten(order="F")
    g.point_data["u_mm"] = np.stack(
        [ug[..., k].flatten(order="F") for k in range(3)], axis=1)
    g.save(str(ruta))
    return ruta


def escala_deformada(res, fraccion=0.08):
    """Factor que amplifica la deformada hasta `fraccion` del lado menor del cubo.

    La deformada real de un hueso a cargas fisiologicas es de micras: invisible.
    El factor se muestra siempre junto a la figura, porque una deformada
    amplificada no es la deformacion real.
    """
    sp = res.get("spacing")
    ug = res.get("u_malla")
    d = float(res.get("desp_max_nodal", 0.0))
    if sp is None or ug is None or d <= 0:
        return 1.0
    lado = float(np.min((np.array(ug.shape[:3]) - 1) * np.asarray(sp)))
    return float(fraccion * lado / d)
