"""
mot_ngsolve.py: Envoltorio de linea de comandos del motor 'ngsolve' de spinpy.

El calculo es `spinpy.motores.m_ngsolve`: la campana mide el MISMO codigo que
ejecuta la aplicacion. Este archivo solo traduce el caso del banco (.npz) y
escribe la salida (ver `_base.principal`).
"""

from ._base import principal, resolver_paquete

if __name__ == "__main__":
    principal(lambda p: resolver_paquete(p, "ngsolve"))
