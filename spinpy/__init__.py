"""
spinpy — nucleo voxel-primero del pipeline de spinodoides, portado a Python.

Cubre la parte de AppFinal_V2.m que va del generador a la morfometria:

    generar_mascara(...)  ->  BW (nx, ny, nz)  ->  morfometria(BW, spacing)

Lee VOIs en VTK legacy (.vtk) y en .mat de MATLAB.

Ajusta un spinodoide a un VOI con la busqueda escalonada de fitSpinodoidToVOI
(etapas A/B/C, correcciones G1/G2/K2) mediante ajustar_spinodoide.

Homogeneiza elasticamente sobre la rejilla de voxeles (homogeneizar).

Mide las curvaturas principales de la interfaz hueso-vacio
(curvaturas_de_campo, perfil_curvaturas), que es lo unico que distingue
una red conectada de un monton de islas con el mismo BV/TV.

Ensaya a compresion en X, Y o Z y estima la carga de fallo con el criterio de
Pistoia, devolviendo los campos de deformacion efectiva y de von Mises
(ensayo_compresion, ensayo_compresion_eje, ensayo_triaxial, criterio_pistoia).

Lo que NO cubre (sigue solo en MATLAB): exportacion FEBio y los otros cuatro
optimizadores (Pareto, bayesiano, MOBO y el best-fit rapido con termino
mecanico).
"""

# Una sola fuente de la version dentro del paquete; `tests/test_16` comprueba
# que coincide con `pyproject.toml`. Viaja en el bloque de procedencia de cada
# resultado (`spinpy.procedencia`).
__version__ = "0.2.0"

from .grf import (campo_grf, canonicalizar_thetas, derivadas_grf, euler_R,
                  generar_mascara, level_set, region_degenerada,
                  wave_directions)
# `perfil` y `resumen` se renombran AL SALIR del modulo: dentro de
# `curvatura.py` se leen bien, pero sueltos en el espacio de nombres del
# paquete no dirian de que son el perfil ni el resumen.
from .curvatura import (curvaturas_de_campo, curvaturas_implicitas,
                        curvaturas_malla)
from .curvatura import perfil as perfil_curvaturas
from .curvatura import resumen as resumen_curvaturas
from .avisos import desalineacion, voi_no_trabecular
from .dual_lattice import campo_dual, generar_dual_lattice
from .elipsoide import factor_elipsoide
from .elastic import (backus_laminado, clase_anisotropia,
                      constantes_ingenieria, constantes_principales,
                      distancia_log_euclidea, hex8_ke, homogeneizar)
from .error import (PESOS_DISCRIMINANTES, TERMINOS_OPCIONALES,
                    error_morfometrico, fabric_cos, peso_alineacion)
from . import incertidumbre, procedencia
from .fit import (ajustar_spinodoide, dir_a_euler, longitud_caracteristica,
                  metricas_mecanicas)
from .fit_dual import ajustar_dual_lattice
from .io import leer_mat_voi, leer_voi, leer_vtk_voi
from .lote import METRICAS, correr_lote, dispersion_semillas
from .voi import (apilar, escribir_vtk_voi, extraer_cubo, marco_pca,
                  vois_por_tercios)
from .morphometry import metricas_forma, morfometria_malla
from .morphometry import (area_superficie, conectividad, fraccion_portante,
                          indice_smi, mayor_componente_6, morfometria,
                          tensor_mil)
from .resistencia import (capa_superficie, criterio_pistoia,
                          ensayo_compresion, ensayo_compresion_eje,
                          ensayo_triaxial, estadisticos_vm,
                          estudio_convergencia)
from .simulacion import fallo_progresivo, simular_perdida

__all__ = [
    "__version__", "procedencia", "incertidumbre",
    "campo_grf", "derivadas_grf", "euler_R", "generar_mascara",
    "level_set", "wave_directions", "region_degenerada",
    "canonicalizar_thetas",
    "clase_anisotropia", "constantes_principales", "distancia_log_euclidea",
    "PESOS_DISCRIMINANTES", "TERMINOS_OPCIONALES", "metricas_forma",
    "curvaturas_de_campo", "curvaturas_implicitas", "curvaturas_malla",
    "perfil_curvaturas", "resumen_curvaturas",
    "desalineacion", "voi_no_trabecular",
    "campo_dual", "generar_dual_lattice", "ajustar_dual_lattice",
    "factor_elipsoide",
    "backus_laminado", "constantes_ingenieria", "hex8_ke", "homogeneizar",
    "error_morfometrico", "fabric_cos", "peso_alineacion",
    "ajustar_spinodoide", "dir_a_euler", "longitud_caracteristica",
    "metricas_mecanicas",
    "leer_mat_voi", "leer_voi", "leer_vtk_voi",
    "METRICAS", "correr_lote", "dispersion_semillas",
    "apilar", "escribir_vtk_voi", "extraer_cubo", "marco_pca",
    "vois_por_tercios",
    "area_superficie", "conectividad", "fraccion_portante", "indice_smi",
    "mayor_componente_6", "morfometria", "tensor_mil",
    "capa_superficie", "estadisticos_vm",
    "criterio_pistoia", "ensayo_compresion", "ensayo_compresion_eje",
    "ensayo_triaxial", "estudio_convergencia",
    "simular_perdida", "fallo_progresivo",
]
