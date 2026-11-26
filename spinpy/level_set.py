"""
level_set.py — Generador level-set-foam VOXEL-PRIMERO (tercera familia).

Espuma de celda abierta como conjunto de nivel de campos de distancia de
elipsoides empaquetados (Sonon, Francois y Massart 2015; Kilingar et al.
2019, ambos comprobados; la documentacion de morphoGen cita ademas a Sonon et
al. 2012 y a Degdeg y Dirrenberger 2026, que no se pudieron comprobar en
ninguna base bibliografica). La implementacion es la de **morphoGen**
(https://github.com/jdirr/morphoGen, Dirrenberger, GPL-3.0-or-later): este
modulo no copia nada de ella, la llama. Los puntales siguen las aristas de la
teselacion de las semillas, con seccion de curvatura controlada (`curvatura`)
y engrosamiento hacia los nodos (`c1`, `c2`).

POR QUE ESTA FAMILIA
--------------------
Se midio sobre cuatro VOI porcinos (`docs/familias/`): con la densidad fija al
BV/TV medido y un objetivo que incluye Conn.D, espesor local, curvatura y
tamano de poro, esta familia reproduce mejor que el spinodoide y el
dual-lattice lo que el ajuste NO optimiza (SMI, numero de poros, rigidez
elastica periodica).

MISMO CONTRATO QUE LAS OTRAS DOS FAMILIAS
-----------------------------------------
    generar_level_set(resolution, n_celdas, rho, ...) -> (BW, campo, info)

BW es booleana (res, res, res) con dim 1 = X, True = hueso; `campo` es la
funcion implicita (negativa dentro del solido) y `info` lleva `rho_objetivo` y
`rho_obtenida`. La densidad es EXACTA por construccion: morphoGen resuelve el
espesor del puntal para la fraccion volumetrica pedida, igual que el
dual-lattice toma el cuantil de la distancia. Por eso aqui, a diferencia del
spinodoide, la densidad no es un parametro libre de la busqueda.

DIFERENCIAS CON LAS OTRAS DOS (declaradas)
------------------------------------------
* **No admite rotacion `R`.** La teselacion es periodica y sus ejes han de
  coincidir con los de la caja. La anisotropia y su orientacion las fija el
  `estiramiento` por eje (que eje se alarga): es el equivalente de los angulos
  de Euler de las otras familias y es lo que el ajuste elige.
* **Es opcional.** morphoGen exige Python >= 3.12 y es GPL-3.0-or-later,
  mientras que el nucleo de spinpy es MIT y admite 3.10. Por eso se importa
  aqui, perezosamente, y el resto del paquete funciona sin ella. Ver
  `disponible()` y `LEEME` del extra `spinpy[level-set]`.
* **Hay combinaciones inviables.** El empaquetado denso exige una caja de tres
  diametros de semilla por eje; con pocas semillas muy alargadas no cabe.
  `es_viable` aplica la misma regla que morphoGen, y el ajuste trata esas
  combinaciones como no medibles, no como error.
"""

from __future__ import annotations

import math

import numpy as np

FAMILIA = "level-set-foam"

# Fraccion volumetrica de las semillas (empaquetado denso de Lubachevsky-
# Stillinger). Es el valor de las plantillas de morphoGen y el que se uso en
# toda la validacion.
FRACCION_EMPAQUETADO = 0.6


class MorphoGenNoDisponible(ImportError):
    """morphoGen no esta instalado o no es utilizable en este Python."""


def _morphogen():
    """Importa morphoGen o explica como instalarlo. Es el UNICO import."""
    try:
        import morphogen as mg
    except ImportError as exc:                       # pragma: no cover
        raise MorphoGenNoDisponible(
            "La familia level-set-foam necesita morphoGen (Python >= 3.12, "
            "licencia GPL-3.0-or-later): pip install \"spinpy[level-set]\" "
            "o pip install morphogen. Ver README, seccion «Familias de "
            "microestructura».") from exc
    return mg


def disponible():
    """True si morphoGen se puede importar; sin lanzar nada si no."""
    try:
        _morphogen()
        return True
    except MorphoGenNoDisponible:
        return False


def version_morphogen():
    """Version de morphoGen instalada, o None."""
    try:
        return str(_morphogen().__version__)
    except MorphoGenNoDisponible:
        return None


def _forma_semilla(n_celdas, estiramiento, domain_size=1.0):
    """Semieje corto `a` y factores (b/a, c/a) de las semillas elipsoidales.

    El estiramiento (ex, ey, ez) son las proporciones de los semiejes en x, y,
    z; se normaliza respecto a x porque morphoGen toma `a` en x y los factores
    de forma en y y z. El numero de semillas fija el radio GEOMETRICO
    (a b c)^(1/3), asi que cambiar el estiramiento no cambia cuantas celdas hay.
    """
    ex, ey, ez = (float(v) for v in estiramiento)
    if min(ex, ey, ez) <= 0.0:
        raise ValueError("estiramiento: tres numeros positivos")
    b, c = ey / ex, ez / ex
    vol = float(domain_size) ** 3
    rg = (FRACCION_EMPAQUETADO * vol
          / (4.0 / 3.0 * math.pi * float(n_celdas))) ** (1.0 / 3.0)
    a = rg / (b * c) ** (1.0 / 3.0)
    return a, (b, c)


def es_viable(n_celdas, estiramiento, domain_size=1.0):
    """True si el empaquetado denso cabe: la regla de morphoGen.

    Cada eje de la caja, dividido por el factor de forma de ese eje, ha de
    medir al menos tres diametros de la semilla (6 a). Se comprueba ANTES de
    generar para que la interfaz y el ajuste no gasten segundos en un caso que
    va a fallar.
    """
    try:
        a, (b, c) = _forma_semilla(n_celdas, estiramiento, domain_size)
    except ValueError:
        return False
    lados = np.array([1.0, 1.0 / b, 1.0 / c]) * float(domain_size)
    return bool(np.all(lados >= 6.0 * a))


def _es_identidad(R):
    return R is None or np.allclose(np.asarray(R, float), np.eye(3), atol=1e-9)


def generar_level_set(resolution, n_celdas, rho, estiramiento=(1.0, 1.0, 1.0),
                      curvatura=0.0, c1=0.0, c2=0.0, R=None, seed=None,
                      domain_size=1.0, mayor_componente=False):
    """Mascara booleana solida (True = hueso), forma (res, res, res), dim 1 = X.

    Misma firma de retorno que `grf.generar_mascara` y
    `dual_lattice.generar_dual_lattice`: (BW, campo, info).

    n_celdas      numero de semillas del empaquetado (~ numero de celdas)
    rho           fraccion de solido, EXACTA (morphoGen resuelve el espesor)
    estiramiento  (ex, ey, ez): proporciones de los semiejes de las semillas;
                  el eje mayor es el de las trabeculas alineadas
    curvatura     K_c de las secciones del puntal: K_c < 0 convexas, 0 < K_c < 1
                  concavas (bordes de Plateau)
    c1, c2        engrosamiento hacia los nodos, A(xi) = c1 xi^4 + c2 xi^2 + 1
    R             solo None o la identidad (ver cabecera)
    """
    if not (0.0 < rho < 1.0):
        raise ValueError(f"rho debe estar en (0,1); se recibio {rho}")
    if not _es_identidad(R):
        raise ValueError(
            "level-set-foam no admite rotacion R: la teselacion es periodica "
            "y sus ejes son los de la caja. La orientacion se da con el "
            "estiramiento por eje.")
    if not float(curvatura) < 1.0:
        raise ValueError("curvatura: K_c < 1")
    if n_celdas < 8:
        raise ValueError("n_celdas: al menos 8 semillas")
    if not es_viable(n_celdas, estiramiento, domain_size):
        raise ValueError(
            "level-set-foam: el empaquetado denso no cabe con "
            f"n_celdas={n_celdas} y estiramiento={tuple(estiramiento)} "
            "(hacen falta tres diametros de semilla por eje); use mas celdas "
            "o un estiramiento menor.")

    mg = _morphogen()
    a, aspecto = _forma_semilla(n_celdas, estiramiento, domain_size)
    dominio = mg.Domain.cube(float(domain_size), periodic=True)
    modelo = mg.models.LevelSetFoam(
        seed_radius=a, aspect=aspecto, packing="dense",
        packing_fraction=FRACCION_EMPAQUETADO, volume_fraction=float(rho),
        curvature=float(curvatura), c1=float(c1), c2=float(c2))
    espuma = modelo.generate(dominio, seed=None if seed is None else int(seed))
    campo = np.asarray(espuma.field(shape=int(resolution)).data)
    BW = campo <= 0.0

    if mayor_componente:
        from .morphometry import mayor_componente_6
        BW = mayor_componente_6(BW)

    info = {
        "familia": FAMILIA,
        "esquema": "no_aplica",
        "rho_objetivo": float(rho),
        "rho_obtenida": float(BW.mean()),
        "n_celdas": float(n_celdas),
        "n_semillas": int(len(getattr(espuma, "centres", ()))) or None,
        "estiramiento": [float(v) for v in estiramiento],
        "curvatura": float(curvatura),
        "c1": float(c1),
        "c2": float(c2),
        "espesor_adimensional": float(getattr(espuma, "thickness", float("nan"))),
        "morphogen": version_morphogen(),
    }
    return BW, campo, info
