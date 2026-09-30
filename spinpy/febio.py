"""
febio.py — Nombre anterior de `spinpy.fem`, conservado por compatibilidad.

Desde V1.1.0 spinpy ya no ejecuta FEBio: los ensayos con malla de ladrillos o
suave, lineales y no lineales, los resuelven los motores internos
(`spinpy.motores`, elegidos en `spinpy.fem`). La EXPORTACION del ensayo a un
archivo .feb para abrirlo en FEBio Studio sigue en `escribe.escribir_febio` y
`escribe.escribir_febio_ensayo`: escribir un archivo no es depender del
programa.

Las validaciones historicas contra FEBio (`comparativa_febio/`,
`comparativa_febio_tet/`) se reproducen con spinpy V1.0.2 y FEBio 4.5.
"""

from .fem import *                                              # noqa: F401,F403
from .fem import (_caras_planos, _componentes_portantes,         # noqa: F401
                  _mayor_componente, _orientar_saliente, _remuestrear,
                  _uz_medio, _uz_medio_nodal, _fuerza, _sigma_ref)
