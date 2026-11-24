"""
_plantilla.py: Punto de partida para un motor FEM nuevo.

COMO USARLA (la guia completa es `docs/AGREGAR_UN_MOTOR_FEM.md`)
-----------------------------------------------------------------
1. Copia este archivo a `m_<clave>.py` (por ejemplo `m_fenics_legacy.py`).
2. Rellena `NOMBRE`, `LICENCIA`, `version()` y `resolver()`.
3. Anade una ficha a `REGISTRO` en `spinpy/motores/__init__.py`.
4. Anade el nombre del motor a `spinpy/idioma_textos.py` y a `DINAMICAS` en
   `idioma_revisar.py`.
5. Ejecuta `python -m pytest tests/test_29_motores.py tests/test_36_contrato_
   motores.py`: el bloque 29 recorre `motores.MOTORES` y pone a prueba el motor
   nuevo sin que escribas ninguna prueba.

Este archivo, tal cual, no es un motor: `resolver` lanza `NoDisponible`. No
esta en `REGISTRO`, de modo que ni la interfaz ni las pruebas lo cargan.

REGLAS DE UN MOTOR
------------------
* NO importes la biblioteca del motor al nivel del modulo. `disponibles()` la
  sondea con `importlib.util.find_spec`, y el proceso de la interfaz no debe
  cargar sus DLL (chocan con las de VTK y Qt). Importa dentro de `resolver`:
  corre en un proceso hijo (`resolver_aislado`).
* Unidades: mm, N, MPa.
* Devuelve el desplazamiento en el ORDEN DE NODOS DE SPINPY, no en el del
  motor. Si el motor renumera, usa `_comun.emparejar` (por coordenadas).
* Si el motor no resuelve una combinacion (elemento, analisis, material),
  lanza `NoDisponible`. Nunca devuelvas el resultado de otro modelo.
* Si no converge, lanza `ErrorMotor` con la causa.

EL PROBLEMA QUE RECIBE `resolver(p)`
------------------------------------
Es un dict; la descripcion completa esta en la cabecera de `_comun.py`:

    p["nodos"]       (N, 3) mm
    p["elems"]       (M, 8) hex8 en orden C3D8, o (M, 10) TET10 en orden C3D10
    p["caras_techo"] caras cargadas (quad4 / tri6), normal +z
    p["caras_base"]  caras de la base
    p["base_nodos"], p["techo_nodos"]   todos los nodos de los planos z0 y z1
    p["meta"]        tipo ('hex8' | 'tet10'), apoyo, E, nu, H, A_bruta,
                     analisis ('lineal' | 'nl'), control ('fuerza' | 'plato'),
                     cargas, material ('svk' | 'neohookeano'), solver, ...

LO QUE DEVUELVE
---------------
`_comun.solucion(u, F, meta)`: `u` (N, 3) en mm; `F` la fuerza de reaccion en
el techo por paso de carga (N, positiva en compresion); `meta` un dict libre
con al menos `{"motor": "<clave>"}`. Ver `m_skfem.py` para un ejemplo corto
en Python puro, y `m_ngsolve.py` para uno con solver directo multihilo.
"""

from __future__ import annotations

from ._comun import NoDisponible

NOMBRE = "Mi motor"            # como se muestra en la interfaz
LICENCIA = "MIT"               # licencia de la BIBLIOTECA del motor


def version():
    """Version de la biblioteca (se llama en el proceso hijo, con ella ya
    importable)."""
    raise NoDisponible("plantilla: sin biblioteca")


def resolver(p):
    """Resuelve el problema `p` y devuelve `solucion(u, F, meta)`."""
    # import mi_biblioteca            # <- aqui, nunca arriba
    # 1. Malla: p["nodos"], p["elems"] -> la del motor.
    # 2. Apoyo: `fijos_basicos(p)` da [(nodo, componente)] fijos.
    # 3. Carga: p["meta"]["control"] == "fuerza" -> presion sobre p["caras_techo"];
    #           == "plato" -> uz impuesto en p["techo_nodos"].
    # 4. Resolver: lineal (K u = f) o Newton (p["meta"]["analisis"] == "nl").
    # 5. u en el orden de nodos de spinpy; F = -suma de f_int,z en el techo.
    raise NoDisponible("plantilla: el motor no esta implementado")
