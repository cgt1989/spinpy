"""
informe.py — Informe para publicacion: metodos, citabilidad y reproduccion.

POR QUE EXISTE
--------------
Entre "la aplicacion calculo un numero" y "ese numero entra en un articulo"
hay tres pasos que hasta ahora dependian de la memoria de quien escribe:

  1. DESCRIBIR lo que se hizo con los valores que de verdad se usaron. Un
     parrafo de metodos redactado a mano desde los deslizadores es justo como
     se transcribio `wave: 15` sin pi (ver `procedencia.py`).
  2. SABER si el numero se puede citar. Los estudios del proyecto midieron
     varias situaciones que producen numeros finitos, de buen aspecto y sin
     ningun error, que no significan lo que parecen. Estan escritas en los
     informes; quien use la aplicacion sin haberlos leido no las conoce.
  3. PERMITIR que otro lo reproduzca, y comprobar que lo consigue.

Este modulo hace las tres cosas sobre el MISMO documento que escribe
"Exportar resultados" en el visor, asi que funciona igual desde la interfaz
que desde un JSON ya exportado:

    spinpy-informe resultados_sesion.json --carpeta informe
    spinpy-informe --verificar informe/reproduccion.json --voi VOI.vtk

Todo va a UNA carpeta: el informe en Markdown y PDF (ES y EN), las figuras en
`figuras/` (PNG a 600 ppp y PDF vectorial; 3D en PNG de 2400 px), el paquete
de reproduccion y el documento de la sesion, que permite rehacer el informe
sin abrir el visor.

No recalcula nada de lo que describe: el informe tiene que hablar de los
numeros que se vieron en pantalla. Lo unico que ejecuta es el generador, para
anotar la huella SHA-256 de cada estructura ajustada, que es lo que el paquete
de reproduccion compara despues.

COMPROBACIONES DE CITABILIDAD
-----------------------------
Devuelven datos —un estado y codigos de motivo—, no frases, igual que
`avisos.py`. Los textos, en los dos idiomas, estan en `MOTIVOS`. Tres estados:

  citable      ninguna comprobacion lo objeta
  reservas     se puede citar declarando la limitacion en el texto
  no_citable   no describe lo que parece; no debe aparecer como resultado

Cada umbral lleva al lado el estudio que lo midio. Dos NO son una medida sino
un criterio nuestro, y se declaran como tal: `DESCONEXION_RESERVAS` y
`DESCONEXION_MAX`. Lo medido es que el ajuste de H4 proximal tiene un 0.5 % de
hueso desconectado y la estructura patologica de beta = 15 un 13.8 %; donde
poner la frontera entre ambos es una decision, no un resultado.

EL TEXTO DE METODOS
-------------------
Se redacta con los valores del documento, frase a frase segun lo que exista:
sin ajuste no hay parrafo de ajuste, y sin homogeneizacion no se nombra el
solver. El numero de onda aparece SIEMPRE con sus dos lecturas (multiplo de pi
y radianes). En espanol los decimales van con coma; en ingles, con punto.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import string
import sys
import time
from pathlib import Path

import numpy as np

from . import procedencia
from .avisos import DA_MIN_EJE, desalineacion, voi_no_trabecular
from .incertidumbre import K_MIN
from .resistencia import (CONV_DERIVA_RHO_MAX, CONV_ERROR_EXTRAPOLADO_MAX,
                          CONV_OSCILACION_MAX, EPS_CRITICA, FRAC_CRITICA)
from .simulacion import TBTH_H_MIN as SIM_TBTH_H_MIN

# ---------------------------------------------------------------------------
# Autoria y archivo del software
# ---------------------------------------------------------------------------

# Nombres como constan en los registros publicos de ORCID. `tests/test_19`
# comprueba que los ORCID coinciden con `CITATION.cff`.
AUTORES = (
    {"nombre": "Carlos", "apellidos": "González-Torres",
     "orcid": "0000-0002-7765-6388"},
    {"nombre": "David", "apellidos": "Ortiz-Puerta",
     "orcid": "0000-0001-6285-3066"},
    {"nombre": "Mauricio A.", "apellidos": "Sarabia-Vallejos",
     "orcid": "0000-0001-5128-796X"},
)

PENDIENTE = "[PENDIENTE]"
# El DOI se rellena al depositar la version en Zenodo.
REPOSITORIO = "https://github.com/cgt1989/spinpy"
DOI_ARCHIVO = PENDIENTE

# ---------------------------------------------------------------------------
# Umbrales
# ---------------------------------------------------------------------------

# Estudio_Convergencia: en el VOI proximal de H4 la rigidez a 32^3 queda ~10 %
# por debajo; desde n = 40 se mueve dentro de +-3 %.
N_MECANICA_MIN = 40
# Correccion F10: un solver iterativo no falla, devuelve su mejor iterado.
RESIDUO_MAX = 1e-6
# Estudio_MIL: por debajo de Tb.Th/h ~ 3 la gaussiana del muestreo suavizado
# (M1) erosiona puntales de dos voxeles.
TBTH_H_M1_MIN = 3.0
# Un percentil 99 es un punto de la cola: por debajo de este numero de
# elementos en la capa superficial lo que se lee es la muestra, no la
# estructura. CRITERIO NUESTRO, no una medida.
N_SUPERFICIE_MIN = 200
# CRITERIO NUESTRO, no una medida (ver la cabecera).
DESCONEXION_RESERVAS = 0.01
DESCONEXION_MAX = 0.05

CITABLE, RESERVAS, NO_CITABLE = "citable", "reservas", "no_citable"
_ORDEN = {CITABLE: 0, RESERVAS: 1, NO_CITABLE: 2}

GRAVEDAD = {
    "no_resuelto": NO_CITABLE,
    "residuo": NO_CITABLE,
    "desconexion_alta": NO_CITABLE,
    "vm_maximo": NO_CITABLE,
    "vm_p99_global": RESERVAS,
    "vm_superficie_pocos": RESERVAS,
    "onda_inconsistente": NO_CITABLE,
    "residuo_no_registrado": RESERVAS,
    "resolucion_mecanica": RESERVAS,
    "desconexion": RESERVAS,
    "pistoia_calibracion": RESERVAS,
    "voi_no_trabecular": RESERVAS,
    "da_sin_eje": RESERVAS,
    "m1_resolucion": RESERVAS,
    "thetas_degenerados": RESERVAS,
    "sin_replicas": RESERVAS,
    "no_bicontinuo": RESERVAS,
    "desalineado": RESERVAS,
    "sin_procedencia": RESERVAS,
    # Estudio de convergencia en malla
    "conv_pocos_puntos": NO_CITABLE,
    "conv_no_monotona": NO_CITABLE,
    "conv_oscila": RESERVAS,
    "conv_deriva_rho": RESERVAS,
    "conv_sin_meseta": RESERVAS,
    "conv_sin_extrapolar": RESERVAS,
    # Dispersion entre semillas y replicas generadas
    "semilla_busqueda": RESERVAS,
    "replicas_no_especimenes": RESERVAS,
    # Simulaciones in silico
    "sim_paso0": NO_CITABLE,
    "sim_pasos_fallidos": RESERVAS,
    "sim_tbth_h": RESERVAS,
    "sim_pendiente": RESERVAS,
    "fallo_modelo_dano": RESERVAS,
}

# (espanol, ingles). Los marcadores se rellenan con `datos` del item; uno que
# falte sale como nan en vez de romper el informe.
MOTIVOS = {
    "no_resuelto": (
        "el cálculo no se resolvió: no hay valor que citar",
        "the computation did not solve: there is no value to cite"),
    "residuo": (
        "residuo relativo del solver {residuo:.1e} > 1e-6: el solver devolvió "
        "su mejor iterado, no una solución (corrección F10)",
        "solver relative residual {residuo:.1e} > 1e-6: the solver returned "
        "its best iterate, not a solution (correction F10)"),
    "residuo_no_registrado": (
        "el residuo del solver iterativo no quedó registrado (resultado "
        "anterior a este informe): repetir el cálculo para poder comprobarlo",
        "the iterative solver residual was not recorded (result older than "
        "this report): rerun the computation to check it"),
    "resolucion_mecanica": (
        "resuelto a {n}³, por debajo de 40³: en el VOI proximal de H4, 32³ "
        "subestima la rigidez ~10 % y desde 40³ se mueve ±3 % "
        "(Estudio_Convergencia)",
        "solved at {n}³, below 40³: on the H4 proximal VOI, 32³ "
        "underestimates stiffness by ~10 % and from 40³ on it moves ±3 % "
        "(Estudio_Convergencia)"),
    "desconexion": (
        "el {desconexion_pct:.1f} % del hueso no une base y techo y se retiró "
        "antes de resolver: declararlo junto al valor",
        "{desconexion_pct:.1f} % of the bone does not connect bottom and top "
        "and was removed before solving: state it next to the value"),
    "desconexion_alta": (
        "el {desconexion_pct:.1f} % del hueso está desconectado (> 5 %): la "
        "rigidez describe una estructura fragmentada",
        "{desconexion_pct:.1f} % of the bone is disconnected (> 5 %): the "
        "stiffness describes a fragmented structure"),
    "vm_maximo": (
        "el máximo de von Mises no converge con la resolución, ni siquiera con "
        "solución cerrada (+5,6 % sobre Goodier en la cavidad esférica): "
        "informar el percentil 99 de la capa superficial (-0,3 %)",
        "the von Mises maximum does not converge with resolution, not even "
        "with a closed-form answer (+5.6 % over Goodier on the spherical "
        "cavity): report the 99th percentile of the surface layer (-0.3 %)"),
    "vm_p99_global": (
        "percentil 99 sobre TODO el tejido: lo que el estudio de convergencia "
        "valida (-0,3 % frente a la solución cerrada) es el percentil 99 de la "
        "capa superficial, y en una estructura densa el global está dominado "
        "por el material a granel, donde la tensión es la nominal",
        "99th percentile over ALL the tissue: what the convergence study "
        "validates (-0.3 % against the closed-form answer) is the 99th "
        "percentile of the surface layer, and in a dense structure the global "
        "one is dominated by bulk material, where the stress is the nominal "
        "one"),
    "vm_superficie_pocos": (
        "la capa superficial tiene {n_superficie:.0f} elementos (< {n_min:.0f}"
        "): un percentil 99 sobre tan pocos describe la muestra, no la "
        "estructura",
        "the surface layer has {n_superficie:.0f} elements (< {n_min:.0f}): a "
        "99th percentile over so few describes the sample, not the structure"),
    "pistoia_calibracion": (
        "criterio de Pistoia con parámetros (2 % del tejido, 0,7 % de "
        "deformación) calibrados en radio distal humano, no en este hueso: "
        "la carga es una estimación comparativa, no una predicción absoluta",
        "Pistoia criterion with parameters (2 % of tissue, 0.7 % strain) "
        "calibrated on the human distal radius, not on this bone: the load "
        "is a comparative estimate, not an absolute prediction"),
    "voi_no_trabecular": (
        "BV/TV del VOI {BVTV:.2f} >= 0,83: hueso compacto con poros aislados; "
        "ninguna familia resulta bicontinua ahí (Estudio_Familias)",
        "VOI BV/TV {BVTV:.2f} >= 0.83: compact bone with isolated pores; no "
        "family is bicontinuous there (Estudio_Familias)"),
    "da_sin_eje": (
        "DA {DA:.3f} < 1,10: la dirección principal es ruido (suelo del MIL "
        "~1,07); no afirmar orientación ni diferencias de DA menores que ~0,07",
        "DA {DA:.3f} < 1.10: the principal direction is noise (MIL floor "
        "~1.07); do not claim an orientation or DA differences below ~0.07"),
    "m1_resolucion": (
        "muestreo MIL suavizado con Tb.Th/h = {tbth_h:.1f} < 3: la gaussiana "
        "erosiona puntales finos (Estudio_MIL)",
        "smoothed MIL sampling with Tb.Th/h = {tbth_h:.1f} < 3: the Gaussian "
        "erodes thin struts (Estudio_MIL)"),
    "thetas_degenerados": (
        "thetas en la región degenerada: la estructura es la isótropa y los "
        "ángulos ajustados no son identificables; no informarlos como "
        "resultado",
        "thetas in the degenerate region: the structure is the isotropic one "
        "and the fitted angles are not identifiable; do not report them as a "
        "result"),
    "sin_replicas": (
        "error de una sola realización (K = {K} < 5): el generador es "
        "estocástico y la búsqueda se queda con la mejor tirada; repetir con "
        "réplicas >= 5",
        "error from a single realisation (K = {K} < 5): the generator is "
        "stochastic and the search keeps the luckiest draw; rerun with "
        "replicas >= 5"),
    "no_bicontinuo": (
        "el candidato no es bicontinuo (poro conexo {poro:.2f} < 0,9): es un "
        "sólido con poros aislados, no una trabécula",
        "the candidate is not bicontinuous (connected pore {poro:.2f} < 0.9): "
        "it is a solid with isolated pores, not a trabecular structure"),
    "desalineado": (
        "eje del candidato a {angulo:.0f}° del VOI (> 30°): el ensayo lo carga "
        "en otra dirección que al hueso",
        "candidate axis at {angulo:.0f}° from the VOI (> 30°): the test loads "
        "it along a different direction than the bone"),
    "onda_inconsistente": (
        "el número de onda no trae sus dos lecturas coherentes (o solo la "
        "clave ambigua 'wave'): la estructura descrita puede ser pi veces más "
        "gruesa que la ajustada",
        "the wave number lacks two consistent readings (or carries only the "
        "ambiguous 'wave' key): the structure described may be pi times "
        "coarser than the fitted one"),
    "sin_procedencia": (
        "documento sin bloque de procedencia (formato 1): no consta la versión "
        "de spinpy, el esquema ni la semilla",
        "document without a provenance block (format 1): spinpy version, "
        "scheme and seed are not recorded"),
    "conv_pocos_puntos": (
        "solo {n_ok:.0f} malla(s) resuelta(s): con menos de tres no hay "
        "estudio de convergencia",
        "only {n_ok:.0f} mesh(es) solved: with fewer than three there is no "
        "convergence study"),
    "conv_no_monotona": (
        "E_app no es monótono al refinar (dispersión {disp_pct:.0f} %): cada "
        "remuestreo cambia qué trabéculas quedan conectadas; ningún E_app de "
        "una sola malla es citable",
        "E_app is not monotonic under refinement (spread {disp_pct:.0f} %): "
        "each resampling changes which trabeculae stay connected; no "
        "single-mesh E_app is citable"),
    "conv_oscila": (
        "E_app oscila sin tendencia al refinar, dentro de un {disp_pct:.1f} % "
        "(banda {max_pct:.0f} %, criterio nuestro): no se extrapola; citar el "
        "valor de la malla más fina con esa banda como incertidumbre",
        "E_app oscillates without a trend under refinement, within "
        "{disp_pct:.1f} % (band {max_pct:.0f} %, our criterion): no "
        "extrapolation; cite the finest-mesh value with that band as its "
        "uncertainty"),
    "conv_deriva_rho": (
        "la densidad cambia un {deriva_pct:.1f} % entre mallas (> {max_pct:.0f}"
        " %, criterio nuestro): el remuestreo cambia la estructura, no solo "
        "la discretización, y no se extrapola",
        "density changes by {deriva_pct:.1f} % between meshes (> "
        "{max_pct:.0f} %, our criterion): resampling changes the structure, "
        "not only the discretisation, and no extrapolation is made"),
    "conv_sin_meseta": (
        "la malla más fina está a un {err_pct:.1f} % del valor extrapolado "
        "(Richardson; umbral {max_pct:.0f} %, criterio nuestro): aún no está "
        "en la meseta; citar el extrapolado con su error o refinar",
        "the finest mesh is {err_pct:.1f} % from the extrapolated value "
        "(Richardson; threshold {max_pct:.0f} %, our criterion): not yet on "
        "the plateau; cite the extrapolated value with its error or refine"),
    "conv_sin_extrapolar": (
        "serie monótona pero sin orden de convergencia calculable: el último "
        "refinado mueve E_app un {salto_pct:.1f} % y no hay estimación del "
        "error de discretización",
        "monotonic series but no computable order of convergence: the last "
        "refinement moves E_app by {salto_pct:.1f} % and there is no estimate "
        "of the discretisation error"),
    "semilla_busqueda": (
        "la primera realización usa la semilla de la búsqueda ({semilla:.0f}),"
        " la que ganó el ajuste: sesga la media hacia el VOI (maldición del "
        "ganador); las réplicas del ajuste usan semilla+1…semilla+K",
        "the first realisation uses the search seed ({semilla:.0f}), the one "
        "that won the fit: it biases the mean towards the VOI (winner's "
        "curse); the fit replicas use seed+1…seed+K"),
    "replicas_no_especimenes": (
        "{n:.0f} realizaciones de UN ajuste, no {n:.0f} especímenes: no "
        "añaden grados de libertad a una comparación biológica",
        "{n:.0f} realisations of ONE fit, not {n:.0f} specimens: they add no "
        "degrees of freedom to a biological comparison"),
    "sim_paso0": (
        "el ensayo del paso 0 no se resolvió: todos los cocientes se refieren "
        "a él",
        "the step-0 test did not solve: every ratio refers to it"),
    "sim_pasos_fallidos": (
        "{n_fallidos:.0f} paso(s) sin solución válida (solver sin resolver o "
        "residuo > 1e-6) quedan fuera de la curva",
        "{n_fallidos:.0f} step(s) without a valid solution (unsolved or "
        "residual > 1e-6) are left out of the curve"),
    "sim_tbth_h": (
        "Tb.Th/h de la malla mecánica = {tbth_h:.2f} < {min:.1f}: con menos "
        "elementos por trabécula los cocientes de rigidez pueden ser "
        "artefacto (Estudio_Convergencia)",
        "Tb.Th/h of the mechanical mesh = {tbth_h:.2f} < {min:.1f}: with "
        "fewer elements per trabecula the stiffness ratios may be an artefact "
        "(Estudio_Convergencia)"),
    "sim_pendiente": (
        "pendiente log-log DESCRIPTIVA: cerca del umbral de rigidez el "
        "exponente aparente se infla (Estudio_Percolacion); sirve para "
        "comparar protocolos, no como exponente de Gibson-Ashby",
        "DESCRIPTIVE log-log slope: near the rigidity threshold the apparent "
        "exponent is inflated (Estudio_Percolacion); it compares protocols, "
        "it is not a Gibson-Ashby exponent"),
    "fallo_modelo_dano": (
        "el tejido roto se ablanda a {rigidez:.2f} E_s en lugar de retirarse: "
        "es un modelo de daño elegido por nosotros (borrarlo desconectaba "
        "puntales), no una ley medida",
        "broken tissue is softened to {rigidez:.2f} E_s instead of removed: "
        "a damage model of our choosing (deleting it disconnected struts), "
        "not a measured law"),
}

ESTADOS = {CITABLE: ("citable", "citable"),
           RESERVAS: ("con reservas", "with caveats"),
           NO_CITABLE: ("no citable", "not citable")}

MAGNITUDES = {
    "procedencia": ("procedencia del documento", "document provenance"),
    "morfometria": ("morfometría (BV/TV, Tb.Th, DA…)",
                    "morphometry (BV/TV, Tb.Th, DA…)"),
    "da": ("DA y dirección principal", "DA and principal direction"),
    "parametros": ("parámetros del ajuste", "fitted parameters"),
    "error_ajuste": ("error del ajuste", "fit error"),
    "elastico": ("tensor C y constantes elásticas",
                 "stiffness tensor C and elastic constants"),
    "E_app": ("módulo aparente E_app", "apparent modulus E_app"),
    "fallo": ("tensión y carga de fallo (Pistoia)",
              "failure stress and load (Pistoia)"),
    "vm_max": ("tensión de von Mises máxima", "maximum von Mises stress"),
    "vm_p99": ("percentil 99 de von Mises en todo el tejido",
               "99th percentile of von Mises over all the tissue"),
    "vm_p99_superficie": ("percentil 99 de von Mises en la capa superficial",
                          "99th percentile of von Mises on the surface layer"),
    "convergencia": ("E_app frente a la malla (convergencia)",
                     "E_app against mesh size (convergence)"),
    "dispersion": ("dispersión entre semillas", "seed-to-seed scatter"),
    "replicas": ("réplicas generadas", "generated replicas"),
    "sim_rigidez": ("rigidez relativa en la pérdida ósea simulada",
                    "relative stiffness under simulated bone loss"),
    "sim_pendiente": ("pendiente log-log E/E0 frente a BV/TV",
                      "log-log slope of E/E0 against BV/TV"),
    "sim_fallo": ("tensión de fallo relativa en la pérdida ósea simulada",
                  "relative failure stress under simulated bone loss"),
    "fallo_progresivo": ("carga máxima del fallo progresivo",
                         "peak load of progressive failure"),
}

ESTRUCTURAS = {"documento": ("documento", "document"), "voi": ("VOI", "VOI"),
               "spinodoide": ("spinodoide", "spinodoid"),
               "dual-lattice": ("dual-lattice", "dual-lattice")}

# Referencias con DOI: copiadas de `paper/paper.bib`, donde se contrastaron con
# los registros de las editoriales. Las que no tienen DOI comprobado lo dicen.
REFERENCIAS = {
    "kumar2020": ("Kumar S, Tan S, Zheng L, Kochmann DM (2020). Inverse-"
                  "designed spinodoid metamaterials. npj Computational "
                  "Materials 6:73.", "10.1038/s41524-020-0341-6"),
    "soyarslan2018": ("Soyarslan C, Bargmann S, Pradas M, Weissmüller J "
                      "(2018). 3D stochastic bicontinuous microstructures: "
                      "Generation, topology and elasticity. Acta Materialia "
                      "149:326–340.", "10.1016/j.actamat.2018.01.005"),
    # Crossref la fecha en 2023 (vol. 138, febrero); la clave conserva 2022,
    # el año de aceptación, con el que la cita el resto del proyecto.
    "vafaeefar2022": ("Vafaeefar M, Moerman KM, Kavousi M, Vaughan TJ (2023). "
                      "A morphological, topological and mechanical "
                      "investigation of gyroid, spinodoid and dual-lattice "
                      "algorithms as structural models of trabecular bone. "
                      "Journal of the Mechanical Behavior of Biomedical "
                      "Materials 138:105584.", "10.1016/j.jmbbm.2022.105584"),
    # Las siete siguientes, para la seccion de modelos. Las tres primeras
    # estan en paper/paper.bib; las otras cuatro se comprobaron en Crossref.
    "hildebrand1997thickness": ("Hildebrand T, Rüegsegger P (1997). A new "
                                "method for the model-independent assessment "
                                "of thickness in three-dimensional images. "
                                "Journal of Microscopy 185(1):67–75.",
                                "10.1046/j.1365-2818.1997.1340694.x"),
    "salmon2015": ("Salmon PL, Ohlsson C, Shefelbine SJ, Doube M (2015). "
                   "Structure model index does not measure rods and plates in "
                   "trabecular bone. Frontiers in Endocrinology 6:162.",
                   "10.3389/fendo.2015.00162"),
    "virtanen2020": ("Virtanen P, Gommers R, Oliphant TE, et al. (2020). SciPy "
                     "1.0: fundamental algorithms for scientific computing in "
                     "Python. Nature Methods 17:261–272.",
                     "10.1038/s41592-019-0686-2"),
    "doube2015": ("Doube M (2015). The Ellipsoid Factor for quantification of "
                  "rods, plates, and intermediate forms in 3D geometries. "
                  "Frontiers in Endocrinology 6:15.",
                  "10.3389/fendo.2015.00015"),
    "taubin1995": ("Taubin G (1995). Curve and surface smoothing without "
                   "shrinkage. Proceedings of IEEE International Conference "
                   "on Computer Vision, 852–857.",
                   "10.1109/ICCV.1995.466848"),
    "hyndman1996": ("Hyndman RJ, Fan Y (1996). Sample quantiles in "
                    "statistical packages. The American Statistician "
                    "50(4):361–365.", "10.1080/00031305.1996.10473566"),
    "vanrietbergen1995": ("van Rietbergen B, Weinans H, Huiskes R, Odgaard A "
                          "(1995). A new method to determine trabecular bone "
                          "elastic properties and loading using "
                          "micromechanical finite-element models. Journal of "
                          "Biomechanics 28(1):69–81.",
                          "10.1016/0021-9290(95)80008-5"),
    "parfitt1987": ("Parfitt AM, Drezner MK, Glorieux FH, Kanis JA, "
                    "Malluche H, Meunier PJ, Ott SM, Recker RR (1987). Bone "
                    "histomorphometry: Standardization of nomenclature, "
                    "symbols, and units. Journal of Bone and Mineral Research "
                    "2(6):595–610.", "10.1002/jbmr.5650020617"),
    "bouxsein2010": ("Bouxsein ML, Boyd SK, Christiansen BA, Guldberg RE, "
                     "Jepsen KJ, Müller R (2010). Guidelines for assessment of "
                     "bone microstructure in rodents using micro-computed "
                     "tomography. Journal of Bone and Mineral Research "
                     "25(7):1468–1486.", "10.1002/jbmr.141"),
    "lorensen1987": ("Lorensen WE, Cline HE (1987). Marching cubes: A high "
                     "resolution 3D surface construction algorithm. ACM "
                     "SIGGRAPH Computer Graphics 21(4):163–169.",
                     "10.1145/37402.37422"),
    "lewiner2003": ("Lewiner T, Lopes H, Vieira AW, Tavares G (2003). "
                    "Efficient implementation of marching cubes' cases with "
                    "topological guarantees. Journal of Graphics Tools "
                    "8(2):1–15.", "10.1080/10867651.2003.10487582"),
    "vanderwalt2014": ("van der Walt S, Schönberger JL, Nunez-Iglesias J, "
                       "Boulogne F, Warner JD, Yager N, Gouillart E, Yu T "
                       "(2014). scikit-image: image processing in Python. "
                       "PeerJ 2:e453.", "10.7717/peerj.453"),
    "harrigan1984": ("Harrigan TP, Mann RW (1984). Characterization of "
                     "microstructural anisotropy in orthotropic materials "
                     "using a second rank tensor. Journal of Materials "
                     "Science 19(3):761–767.", "10.1007/BF00540446"),
    "odgaard1993": ("Odgaard A, Gundersen HJG (1993). Quantification of "
                    "connectivity in cancellous bone, with special emphasis "
                    "on 3-D reconstructions. Bone 14(2):173–182.",
                    "10.1016/8756-3282(93)90245-6"),
    "hildebrand1997smi": ("Hildebrand T, Rüegsegger P (1997). Quantification "
                          "of bone microarchitecture with the Structure Model "
                          "Index. Computer Methods in Biomechanics and "
                          "Biomedical Engineering 1(1):15–23.",
                          "10.1080/01495739708936692"),
    "andreassen2014": ("Andreassen E, Andreasen CS (2014). How to determine "
                       "composite material properties using numerical "
                       "homogenization. Computational Materials Science "
                       "83:488–495.", "10.1016/j.commatsci.2013.09.006"),
    "bell2022": ("Bell N, Olson LN, Schroder J (2022). PyAMG: Algebraic "
                 "Multigrid Solvers in Python. Journal of Open Source Software "
                 "7(72):4142.", "10.21105/joss.04142"),
    "pistoia2002": ("Pistoia W, van Rietbergen B, Lochmüller E-M, Lill CA, "
                    "Eckstein F, Rüegsegger P (2002). Estimation of distal "
                    "radius failure load with micro-finite element analysis "
                    "models based on three-dimensional peripheral quantitative "
                    "computed tomography images. Bone 30(6):842–848.",
                    "10.1016/S8756-3282(02)00736-6"),
    "harris2020": ("Harris CR, Millman KJ, van der Walt SJ, et al. (2020). "
                   "Array programming with NumPy. Nature 585:357–362.",
                   "10.1038/s41586-020-2649-2"),
}


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------

class _Formateador(string.Formatter):
    """`str.format` con coma decimal en espanol y nan para lo que falte."""

    def __init__(self, idioma):
        super().__init__()
        self.idioma = idioma

    def get_value(self, key, args, kwargs):
        if isinstance(key, str):
            v = kwargs.get(key)
            return float("nan") if v is None else v
        return super().get_value(key, args, kwargs)

    def format_field(self, valor, espec):
        es_num = isinstance(valor, (int, float, np.integer, np.floating)) \
            and not isinstance(valor, bool)
        s = super().format_field(valor, espec)
        if es_num and self.idioma == "es":
            s = s.replace(".", ",")
        return s


def _f(d, clave):
    """Escalar finito de un dict, o None."""
    if not isinstance(d, dict):
        return None
    try:
        x = float(np.asarray(d.get(clave), float).reshape(()))
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def _i(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _idx(idioma):
    return 0 if idioma == "es" else 1


def sha256_archivo(ruta):
    """SHA-256 de un archivo, o None si la ruta no es un archivo legible.

    Una carpeta de rebanadas no tiene una huella unica y devuelve None: el
    informe lo declara en vez de inventar una.
    """
    try:
        p = Path(ruta)
        if not p.is_file():
            return None
        h = hashlib.sha256()
        with open(p, "rb") as f:
            for bloque in iter(lambda: f.read(1 << 20), b""):
                h.update(bloque)
        return h.hexdigest()
    except (TypeError, OSError):
        return None


def huella_mascara(BW):
    """SHA-256 de una mascara binaria: forma + bits en orden Fortran.

    Orden Fortran porque es el de los indices de elemento del proyecto; la
    forma entra en la huella para que dos mascaras con los mismos bits y otra
    forma no coincidan.
    """
    BW = np.asarray(BW, dtype=bool)
    h = hashlib.sha256()
    h.update(repr(tuple(int(s) for s in BW.shape)).encode("ascii"))
    h.update(np.packbits(BW.ravel(order="F")).tobytes())
    return h.hexdigest()


def _regenerar(par, completo=False):
    """Mascara de un ajuste desde sus parametros, por la unica traduccion.

    completo=True devuelve (BW, campo, info): la figura 0 dibuja el campo y
    el umbral, y asi no se llama dos veces al generador.
    """
    kw = procedencia.kwargs_generador(par)
    if par.get("familia") == "dual-lattice":
        from .dual_lattice import generar_dual_lattice
        BW, c, inf = generar_dual_lattice(**kw)
    else:
        from .grf import generar_mascara
        BW, c, inf = generar_mascara(**kw)
    return (BW, c, inf) if completo else BW


def _lado_voi(doc):
    voi = doc.get("voi") or {}
    forma, sp = voi.get("forma"), voi.get("spacing_mm")
    try:
        return float(forma[0]) * float(np.ravel(sp)[0])
    except (TypeError, IndexError, ValueError):
        return None


def _morf(doc, estructura):
    return doc.get({"voi": "morfometria_voi", "spinodoide": "morfometria_spin",
                    "dual-lattice": "morfometria_dual"}[estructura])


# ---------------------------------------------------------------------------
# Citabilidad
# ---------------------------------------------------------------------------

def _item(bloque, magnitud, estructura, valor, motivos, datos=None, eje=None):
    motivos = list(dict.fromkeys(motivos))
    estado = CITABLE
    for m in motivos:
        if _ORDEN[GRAVEDAD[m]] > _ORDEN[estado]:
            estado = GRAVEDAD[m]
    return {"bloque": bloque, "magnitud": magnitud, "eje": eje,
            "estructura": estructura,
            "valor": {k: v for k, v in (valor or {}).items() if v is not None},
            "estado": estado, "motivos": motivos, "datos": dict(datos or {})}


def _items_vm(bloque, d, est, mot, dat, eje=None):
    """Items de citabilidad para las claves de von Mises de un registro.

    Un solo sitio decide como se juzga cada sabor del pico, porque hay tres y
    se parecen demasiado:

      vm_max*            no converge con la resolucion -> no citable
      vm_p99*            percentil sobre TODO el tejido -> con reservas
      vm_p99_superficie* percentil sobre la capa libre  -> citable

    El orden de las comprobaciones importa: `vm_p99_superficie` empieza por
    `vm_p99`, asi que se mira antes el caso largo.
    """
    n_sup = _f(d, "vm_n_superficie")
    out = []
    for kk in d:
        k = str(kk)
        if not k.startswith("vm_"):
            continue
        extra = ["pistoia_calibracion"] if k.endswith("_fallo") else []
        if k.startswith("vm_max"):
            out.append(_item(bloque, "vm_max", est, {kk: _f(d, kk)},
                             mot + ["vm_maximo"] + extra, dat, eje))
        elif k.startswith("vm_p99_superficie"):
            m = list(extra)
            da = dict(dat)
            if n_sup is not None and n_sup < N_SUPERFICIE_MIN:
                m.append("vm_superficie_pocos")
                da["n_superficie"] = n_sup
                da["n_min"] = float(N_SUPERFICIE_MIN)
            out.append(_item(bloque, "vm_p99_superficie", est,
                             {kk: _f(d, kk)}, mot + m, da, eje))
        elif k.startswith("vm_p99"):
            out.append(_item(bloque, "vm_p99", est, {kk: _f(d, kk)},
                             mot + ["vm_p99_global"] + extra, dat, eje))
    return out


def _motivos_mecanicos(n, residuo, frac_portante, ok=True):
    mot, dat = [], {}
    if not ok:
        mot.append("no_resuelto")
    if residuo is not None:
        dat["residuo"] = residuo
        if residuo > RESIDUO_MAX:
            mot.append("residuo")
    if n is not None:
        dat["n"] = n
        if n < N_MECANICA_MIN:
            mot.append("resolucion_mecanica")
    if frac_portante is not None:
        des = 1.0 - frac_portante
        dat["desconexion_pct"] = 100.0 * des
        if des > DESCONEXION_MAX:
            mot.append("desconexion_alta")
        elif des > DESCONEXION_RESERVAS:
            mot.append("desconexion")
    return mot, dat


def comprobar(doc):
    """Lista de items de citabilidad de un documento de resultados.

    Cada item: {"bloque", "magnitud", "eje", "estructura", "valor", "estado",
    "motivos", "datos"}. `magnitud`, `estructura` y `motivos` son codigos;
    `etiqueta_item` y `texto_motivos` los convierten en texto.
    """
    items = []
    res = doc.get("resultados") or {}
    m_voi = doc.get("morfometria_voi") or None

    proc = doc.get("procedencia")
    if (not isinstance(proc, dict) or
            (_i(proc.get("version_formato")) or 1) < procedencia.VERSION_FORMATO):
        items.append(_item("procedencia", "procedencia", "documento", {},
                           ["sin_procedencia"]))

    voi_nt = voi_no_trabecular(m_voi) if m_voi else {"aviso": False}
    lado = _lado_voi(doc)
    muestreo = (doc.get("muestreo_mil")
                or (res.get("morfometria") or {}).get("muestreo_mil")
                or "voxel")

    # --- morfometria ---
    for est in ("voi", "spinodoide", "dual-lattice"):
        m = _morf(doc, est)
        if not m:
            continue
        mot, dat = [], {}
        if est == "voi" and voi_nt["aviso"]:
            mot.append("voi_no_trabecular")
            dat["BVTV"] = voi_nt["BVTV"]
        if est == "voi":
            h = (lado / float((doc.get("voi") or {}).get("forma", [0])[0])
                 if lado else None)
        else:
            n_m = _i((res.get("morfometria") or {}).get("resolucion"))
            h = lado / n_m if (lado and n_m) else None
        tbth = _f(m, "TbTh")
        if muestreo == "suavizado" and h and tbth:
            dat["tbth_h"] = tbth / h
            if dat["tbth_h"] < TBTH_H_M1_MIN:
                mot.append("m1_resolucion")
        items.append(_item("morfometria", "morfometria", est,
                           {"BV/TV": _f(m, "BVTV"), "Tb.Th [mm]": tbth,
                            "DA": _f(m, "DA")}, mot, dat))
        da = _f(m, "DA")
        if da is not None and da < DA_MIN_EJE:
            items.append(_item("morfometria", "da", est, {"DA": da},
                               ["da_sin_eje"], {"DA": da}))

    # --- ajustes ---
    for clave, fam in (("ajuste", "spinodoide"), ("ajuste_dual", "dual-lattice")):
        rec = res.get(clave)
        if not isinstance(rec, dict) or not rec.get("parametros"):
            continue
        par = rec["parametros"]
        mot, dat = [], {}
        if fam == "spinodoide":
            try:
                if ("wave" in par and "wave_number_pi" not in par
                        and "wave_number_rad" not in par):
                    raise procedencia.ErrorNumeroOnda("clave 'wave'")
                procedencia.numero_onda(par.get("wave_number_pi"),
                                        par.get("wave_number_rad"))
            except procedencia.ErrorNumeroOnda:
                mot.append("onda_inconsistente")
            if par.get("thetas_degenerados"):
                mot.append("thetas_degenerados")
            valor = {"rho": _f(par, "densidad"),
                     "beta/pi": _f(par, "wave_number_pi"),
                     "N": _f(par, "num_waves")}
        else:
            valor = {"rho": _f(par, "densidad"), "celdas": _f(par, "celdas")}
        diag = rec.get("diagnostico") or {}
        if diag.get("bicontinuo") is False:
            mot.append("no_bicontinuo")
            dat["poro"] = _f(diag, "poro_conexo_spin")
        if voi_nt["aviso"]:
            mot.append("voi_no_trabecular")
            dat["BVTV"] = voi_nt["BVTV"]
        items.append(_item("ajuste", "parametros", fam, valor, mot, dat))

        inc = rec.get("incertidumbre") or {}
        K = _i(inc.get("K")) or 0
        err = inc.get("error") or {}
        valor = {"error": _f(rec, "error")}
        if K >= K_MIN:
            valor.update(media=_f(err, "media"), sd=_f(err, "sd"))
        items.append(_item("ajuste", "error_ajuste", fam, valor,
                           [] if K >= K_MIN else ["sin_replicas"], {"K": K}))

    # --- mecanica ---
    def desal(fam):
        m_c = _morf(doc, fam)
        if not (m_voi and m_c):
            return [], {}
        d = desalineacion(m_voi, m_c)
        if d["aviso"]:
            return ["desalineado"], {"angulo": d["angulo_deg"]}
        return [], {}

    for fam, suf in (("spinodoide", ""), ("dual-lattice", "_dual")):
        rec = res.get("elastico" + suf)
        if isinstance(rec, dict):
            n = _i(rec.get("resolucion"))
            for k, v in (rec.get("por_estructura") or {}).items():
                est = "voi" if k == "voi" else fam
                ok = _f(v, "Ex") is not None
                residuo = _f(v, "residuo_rel")
                mot, dat = _motivos_mecanicos(n, residuo, None, ok)
                if (ok and residuo is None
                        and "LU" not in str(v.get("solver", ""))):
                    mot.append("residuo_no_registrado")
                if est != "voi":
                    m2, d2 = desal(fam)
                    mot += m2
                    dat.update(d2)
                items.append(_item(
                    "elastico", "elastico", est,
                    {e + " [MPa]": (None if _f(v, e) is None else _f(v, e) / 1e6)
                     for e in ("Ex", "Ey", "Ez")}, mot, dat))

        rec = res.get("resistencia" + suf)
        if isinstance(rec, dict):
            n = _i(rec.get("resolucion"))
            for k, v in (rec.get("por_estructura") or {}).items():
                est = "voi" if k == "voi" else fam
                for eje, p in sorted((v.get("ejes") or {}).items()):
                    mot, dat = _motivos_mecanicos(
                        n, _f(p, "residuo_rel"), _f(p, "frac_portante"),
                        _f(p, "E_app") is not None)
                    if est != "voi":
                        m2, d2 = desal(fam)
                        mot += m2
                        dat.update(d2)
                    E = _f(p, "E_app")
                    items.append(_item("resistencia", "E_app", est,
                                       {"E_app [MPa]": None if E is None
                                        else E / 1e6}, mot, dat, eje))
                    if "sigma_fallo" in p:
                        s = _f(p, "sigma_fallo")
                        items.append(_item(
                            "resistencia", "fallo", est,
                            {"sigma_fallo [MPa]": None if s is None
                             else s / 1e6},
                            mot + ["pistoia_calibracion"], dat, eje))
                    items += _items_vm("resistencia", p, est, mot, dat, eje)

        rec = res.get("analisis_comparado" + suf)
        if isinstance(rec, dict):
            n = _i(rec.get("resolucion"))
            for k, v in (rec.get("por_estructura") or {}).items():
                est = "voi" if k == "voi" else fam
                mot, dat = _motivos_mecanicos(
                    n, _f(v, "residuo"), _f(v, "frac_portante"),
                    _f(v, "E_app") is not None)
                if est != "voi":
                    m2, d2 = desal(fam)
                    mot += m2
                    dat.update(d2)
                items.append(_item("analisis_comparado", "E_app", est,
                                   {"E_app": _f(v, "E_app")}, mot, dat))
                items += _items_vm("analisis_comparado", v, est, mot, dat)
                # Las de von Mises ya se han juzgado arriba con su propio
                # criterio; sin excluirlas, `vm_p99_superficie_fallo` entraria
                # ademas aqui por contener "fallo" y saldria dos veces.
                claves_fallo = [kk for kk in v if "fallo" in str(kk)
                                and not str(kk).startswith("vm_")]
                if claves_fallo:
                    items.append(_item(
                        "analisis_comparado", "fallo", est,
                        {kk: _f(v, kk) for kk in claves_fallo},
                        mot + ["pistoia_calibracion"], dat))

    items += _items_convergencia(doc, res)
    items += _items_dispersion(doc, res)
    items += _items_simulaciones(doc, res)
    return items


def estructura_de(codigo, doc):
    """Codigo de estructura ("voi", "spinodoide", "dual-lattice") de un
    registro de convergencia o de simulacion.

    Los registros nuevos traen `estructura_codigo`. Los anteriores solo
    guardaban la etiqueta de pantalla, que esta traducida; entonces se deduce:
    "VOI" es el VOI y cualquier otra cosa el candidato, de la familia cuyo
    nombre aparezca en la etiqueta, y si ninguna, el spinodoide. Sin etiqueta,
    el VOI si el documento lo tiene, porque es lo que el visor elige primero.
    """
    if codigo in ESTRUCTURAS and codigo != "documento":
        return codigo
    txt = str(codigo or "").lower()
    if not txt:
        return "voi" if doc.get("voi") else "spinodoide"
    if txt.strip() == "voi":
        return "voi"
    return "dual-lattice" if "dual" in txt else "spinodoide"


def _items_convergencia(doc, res):
    rec = res.get("convergencia")
    if not isinstance(rec, dict) or not rec.get("puntos"):
        return []
    est = estructura_de(rec.get("estructura_codigo")
                        or rec.get("estructura"), doc)
    val = [p for p in rec["puntos"] if p.get("ok")]
    mot, dat = [], {"n_ok": float(len(val))}
    valor = {}
    if len(val) < 3:
        mot.append("conv_pocos_puntos")
    if val:
        fino = max(val, key=lambda p: int(p.get("n", 0)))
        E = _f(fino, "E_app")
        valor["E_app [MPa]"] = None if E is None else E / 1e6
        # El residuo que manda es el PEOR de la serie: un punto mal resuelto
        # tuerce la curva entera, no solo su propio valor.
        residuos = [_f(p, "residuo") for p in val
                    if _f(p, "residuo") is not None]
        m2, d2 = _motivos_mecanicos(_i(fino.get("n")),
                                    max(residuos) if residuos else None,
                                    _f(fino, "frac_portante"))
        mot += m2
        dat.update(d2)
    if len(val) >= 3:
        disp = _f(rec, "dispersion_rel")
        if disp is not None and disp < 1e-6:
            pass                          # serie plana: ya convergida
        elif rec.get("monotona") is False:
            osc = disp is not None and disp <= CONV_OSCILACION_MAX
            mot.append("conv_oscila" if osc else "conv_no_monotona")
            dat["disp_pct"] = 100.0 * (disp if disp is not None
                                       else float("nan"))
            dat["max_pct"] = 100.0 * CONV_OSCILACION_MAX
        else:
            deriva = _f(rec, "deriva_rho_rel")
            err = _f(rec, "error_estimado_rel")
            if deriva is not None and deriva >= CONV_DERIVA_RHO_MAX:
                mot.append("conv_deriva_rho")
                dat["deriva_pct"] = 100.0 * deriva
                dat["max_pct"] = 100.0 * CONV_DERIVA_RHO_MAX
            elif err is None:
                mot.append("conv_sin_extrapolar")
                salto = _f(rec, "salto_final_rel")
                dat["salto_pct"] = (100.0 * salto if salto is not None
                                    else float("nan"))
            elif err >= CONV_ERROR_EXTRAPOLADO_MAX:
                mot.append("conv_sin_meseta")
                dat["err_pct"] = 100.0 * err
                dat["max_pct"] = 100.0 * CONV_ERROR_EXTRAPOLADO_MAX
            Ex = _f(rec, "E_extrapolado")
            if Ex is not None:
                valor["E_extrap [MPa]"] = Ex / 1e6
    return [_item("convergencia", "convergencia", est, valor, mot, dat,
                  rec.get("eje"))]


def _semilla_ajuste(rec):
    if not isinstance(rec, dict):
        return None
    for d in (rec.get("procedencia") or {}, rec.get("parametros") or {}):
        for k in ("semilla", "seed"):
            if d.get(k) is not None:
                return _i(d.get(k))
    return None


def _items_dispersion(doc, res):
    items = []
    for fam, suf in (("spinodoide", ""), ("dual-lattice", "_dual")):
        rec = res.get("dispersion" + suf)
        if isinstance(rec, dict) and rec.get("resumen"):
            K = _i(rec.get("n_semillas")) or 0
            mot, dat = [], {"K": K}
            if K < K_MIN:
                mot.append("sin_replicas")
            base = _i(rec.get("semilla_base"))
            if base is not None and base == _semilla_ajuste(
                    res.get("ajuste" + suf)):
                mot.append("semilla_busqueda")
                dat["semilla"] = float(base)
            items.append(_item("dispersion", "dispersion", fam, {"K": K},
                               mot, dat))
        rec = res.get("replicas" + suf)
        if isinstance(rec, dict) and rec.get("n"):
            n = _i(rec.get("n"))
            items.append(_item("replicas", "replicas", fam, {"n": n},
                               ["replicas_no_especimenes"], {"n": float(n)}))
    return items


def _pasos_fallidos(filas, sin_camino_vale=True):
    """Pasos sin solucion valida. `sin_camino` NO es un fallo: una estructura
    sin camino portante tiene rigidez cero, y eso es un resultado."""
    n = 0
    for f in filas:
        if sin_camino_vale and f.get("sin_camino"):
            continue
        r = _f(f, "residuo")
        if not f.get("ok") or (r is not None and r > RESIDUO_MAX):
            n += 1
    return n


def _eje_txt(par):
    return {0: "x", 1: "y", 2: "z"}.get(_i((par or {}).get("eje")))


def _items_simulaciones(doc, res):
    items = []
    rec = res.get("simulacion_perdida")
    if isinstance(rec, dict) and rec.get("pasos"):
        est = estructura_de(rec.get("estructura_codigo")
                            or rec.get("estructura"), doc)
        filas, r = rec["pasos"], rec.get("resumen") or {}
        par = rec.get("parametros") or {}
        f0 = filas[0]
        mot, dat = _motivos_mecanicos(_i(f0.get("n_mec") or par.get("n_mec")),
                                      _f(f0, "residuo"),
                                      _f(f0, "frac_portante"))
        if not f0.get("ok") or _f(f0, "E_app") is None:
            mot.append("sim_paso0")
        n_f = _pasos_fallidos(filas[1:])
        if n_f:
            mot.append("sim_pasos_fallidos")
            dat["n_fallidos"] = float(n_f)
        tbth_h = _f(r, "TbTh_h_mec")
        if tbth_h is not None and tbth_h < SIM_TBTH_H_MIN:
            mot.append("sim_tbth_h")
            dat.update(tbth_h=tbth_h, min=SIM_TBTH_H_MIN)
        eje = _eje_txt(par)
        items.append(_item("simulacion", "sim_rigidez", est,
                           {"E/E0": _f(r, "E_rel_final"),
                            "BV/TV / BV/TV0": _f(r, "BVTV_rel_final")},
                           mot, dat, eje))
        if _f(r, "pendiente_loglog") is not None:
            items.append(_item("simulacion", "sim_pendiente", est,
                               {"pendiente": _f(r, "pendiente_loglog")},
                               mot + ["sim_pendiente"], dat, eje))
        s_fin = _f(filas[-1], "sigma_fallo_rel")
        if s_fin is not None:
            items.append(_item("simulacion", "sim_fallo", est,
                               {"sigma/sigma0": s_fin},
                               mot + ["pistoia_calibracion"], dat, eje))

    rec = res.get("fallo_progresivo")
    if isinstance(rec, dict) and rec.get("pasos"):
        est = estructura_de(rec.get("estructura_codigo")
                            or rec.get("estructura"), doc)
        filas, r = rec["pasos"], rec.get("resumen") or {}
        par = rec.get("parametros") or {}
        f0 = filas[0]
        mot, dat = _motivos_mecanicos(_i(par.get("n_mec")), _f(f0, "residuo"),
                                      _f(f0, "frac_portante"))
        if not f0.get("ok") or _f(f0, "E_app") is None:
            mot.append("sim_paso0")
        # Solo cuentan los pasos ANTES del colapso: despues, no resolver es
        # justamente lo que se espera de una estructura rota.
        antes = [f for f in filas[1:] if not f.get("tras_colapso")]
        n_f = _pasos_fallidos(antes, sin_camino_vale=False)
        if n_f:
            mot.append("sim_pasos_fallidos")
            dat["n_fallidos"] = float(n_f)
        tbth_h = _f(r, "TbTh_h_mec")
        if tbth_h is not None and tbth_h < SIM_TBTH_H_MIN:
            mot.append("sim_tbth_h")
            dat.update(tbth_h=tbth_h, min=SIM_TBTH_H_MIN)
        mot += ["pistoia_calibracion", "fallo_modelo_dano"]
        dat["rigidez"] = _f(par, "rigidez_danada")
        items.append(_item("simulacion", "fallo_progresivo", est,
                           {"F_max [N]": _f(r, "F_max"),
                            "E/E0": _f(r, "E_rel_final")}, mot, dat,
                           _eje_txt(par)))
    return items


def resumen(items):
    out = {CITABLE: 0, RESERVAS: 0, NO_CITABLE: 0}
    for it in items:
        out[it["estado"]] += 1
    return out


def texto_motivo(codigo, datos, idioma="es"):
    return _Formateador(idioma).format(MOTIVOS[codigo][_idx(idioma)],
                                       **(datos or {}))


def texto_motivos(item, idioma="es"):
    if not item["motivos"]:
        return "—"
    return "; ".join(texto_motivo(c, item["datos"], idioma)
                     for c in item["motivos"])


def etiqueta_item(item, idioma="es"):
    i = _idx(idioma)
    txt = (ESTRUCTURAS[item["estructura"]][i] + " — "
           + MAGNITUDES[item["magnitud"]][i])
    return txt + (f" [{item['eje']}]" if item.get("eje") else "")


def _valor_txt(valor, idioma):
    fm = _Formateador(idioma)
    if not valor:
        return "—"
    return "; ".join(f"{k} = {fm.format_field(float(v), '.4g')}"
                     for k, v in valor.items())


# ---------------------------------------------------------------------------
# Texto de metodos
# ---------------------------------------------------------------------------

class _Redactor:
    def __init__(self, idioma):
        self.idioma = idioma
        self.orden = []
        self._fm = _Formateador(idioma)

    def t(self, es, en, **v):
        return self._fm.format(es if self.idioma == "es" else en, **v)

    def n(self, x, espec=".4g"):
        try:
            return self._fm.format_field(float(x), espec)
        except (TypeError, ValueError):
            return "—"

    def c(self, *claves):
        nums = []
        for k in claves:
            if k not in self.orden:
                self.orden.append(k)
            nums.append(self.orden.index(k) + 1)
        return "[" + ", ".join(str(x) for x in nums) + "]"


def _solver_txt(s, idioma):
    s = str(s or "")
    if "AMG" in s:
        return ("gradiente conjugado precondicionado con multimalla algebraica"
                if idioma == "es" else
                "conjugate gradients preconditioned with algebraic multigrid")
    if "LU" in s:
        return ("factorización LU directa" if idioma == "es"
                else "direct LU factorisation")
    if "Jacobi" in s:
        return ("gradiente conjugado con precondicionador de Jacobi"
                if idioma == "es" else
                "conjugate gradients with a Jacobi preconditioner")
    return s or "—"


def _autores_cortos(idioma):
    ap = [a["apellidos"] for a in AUTORES]
    y = " y " if idioma == "es" else " and "
    return ", ".join(ap[:-1]) + y + ap[-1]


def parrafos_metodos(doc, idioma="es", r=None):
    """Devuelve (parrafos, claves_de_referencia_en_orden_de_cita).

    `r` permite compartir la numeracion de citas con otra seccion (la de
    modelos): la lista devuelta es la misma que va creciendo en `r`.
    """
    r = r or _Redactor(idioma)
    P = []
    proc = doc.get("procedencia") or {}
    res = doc.get("resultados") or {}
    voi = doc.get("voi") or {}
    from . import __version__

    P.append(r.t(
        "Los análisis se realizaron con spinpy {v} ({aut}; repositorio {rep}; "
        "DOI {doi}), en Python {py} con NumPy {npv} {c}.",
        "Analyses were performed with spinpy {v} ({aut}; repository {rep}; "
        "DOI {doi}), on Python {py} with NumPy {npv} {c}.",
        v=proc.get("spinpy") or __version__, aut=_autores_cortos(idioma),
        rep=REPOSITORIO, doi=DOI_ARCHIVO,
        py=proc.get("python") or platform.python_version(),
        npv=proc.get("numpy") or np.__version__, c=r.c("harris2020")))

    m_voi = doc.get("morfometria_voi")
    hay_cand = bool(doc.get("morfometria_spin") or doc.get("morfometria_dual"))
    forma, sp = voi.get("forma"), voi.get("spacing_mm")
    if forma and sp:
        P.append(r.t(
            "El volumen de interés ({nom}) se analizó como imagen binaria de "
            "{dims} vóxeles, con un tamaño de vóxel de {um} µm ({lado} mm de "
            "lado).",
            "The volume of interest ({nom}) was analysed as a binary image of "
            "{dims} voxels with a voxel size of {um} µm ({lado} mm side).",
            nom=voi.get("nombre") or "VOI",
            dims="×".join(str(int(x)) for x in forma),
            um=r.n(float(np.ravel(sp)[0]) * 1000.0, ".4g"),
            lado=r.n(_lado_voi(doc), ".4g")))

    if m_voi or hay_cand:
        muestreo = (doc.get("muestreo_mil")
                    or (res.get("morfometria") or {}).get("muestreo_mil")
                    or "voxel")
        txt = r.t(
            "La morfometría siguió la nomenclatura de Parfitt et al. y las "
            "recomendaciones de Bouxsein et al. {c1}. BV/TV se obtuvo por "
            "conteo de vóxeles y la superficie ósea (BS) con marching cubes "
            "{c2}, excluyendo las caras de corte del volumen; el espesor "
            "trabecular se calculó con el modelo de placas, Tb.Th = 2·BV/BS, y "
            "de él Tb.Sp y Tb.N. El grado de anisotropía (DA) se obtuvo del "
            "tensor de longitud media de intercepción (MIL) {c3} con muestreo "
            "{mu}; diferencias de DA menores que ~0,07 quedan por debajo del "
            "suelo de ruido del estimador.",
            "Morphometry followed the nomenclature of Parfitt et al. and the "
            "guidelines of Bouxsein et al. {c1}. BV/TV was obtained by voxel "
            "counting and bone surface (BS) by marching cubes {c2}, excluding "
            "the cut faces of the volume; trabecular thickness used the plate "
            "model, Tb.Th = 2·BV/BS, from which Tb.Sp and Tb.N were derived. "
            "The degree of anisotropy (DA) was obtained from the mean "
            "intercept length (MIL) tensor {c3} with {mu} sampling; DA "
            "differences below ~0.07 lie under the estimator's noise floor.",
            c1=r.c("parfitt1987", "bouxsein2010"),
            c2=r.c("lorensen1987", "lewiner2003", "vanderwalt2014"),
            c3=r.c("harrigan1984"),
            mu=(r.t("suavizado (gaussiana σ = 0,7 vóxeles e interpolación "
                    "trilineal)", "smoothed (Gaussian σ = 0.7 voxels and "
                    "trilinear interpolation)")
                if muestreo == "suavizado" else r.t("por vóxeles", "voxel")))
        if _f(m_voi or {}, "ConnD") is not None or \
                (res.get("morfometria") or {}).get("extra"):
            txt += " " + r.t(
                "La densidad de conectividad (Conn.D) se obtuvo de la "
                "característica de Euler {c1} y el índice de modelo "
                "estructural (SMI) según {c2}.",
                "Connectivity density (Conn.D) was obtained from the Euler "
                "characteristic {c1} and the structure model index (SMI) "
                "following {c2}.",
                c1=r.c("odgaard1993"), c2=r.c("hildebrand1997smi"))
        if hay_cand:
            txt += " " + r.t(
                "Cada candidato sintético se escaló al tamaño físico del VOI y "
                "se midió por el mismo camino que el VOI, de modo que el sesgo "
                "de marching cubes sobre datos binarios afecta por igual a "
                "ambos.",
                "Each synthetic candidate was scaled to the physical size of "
                "the VOI and measured through the same code path as the VOI, "
                "so the bias of marching cubes on binary data affects both "
                "equally.")
        P.append(txt)

    for clave, fam in (("ajuste", "spinodoide"), ("ajuste_dual", "dual-lattice")):
        rec = res.get(clave)
        if not isinstance(rec, dict) or not rec.get("parametros"):
            continue
        par = rec["parametros"]
        if fam == "spinodoide":
            try:
                onda = procedencia.numero_onda(par.get("wave_number_pi"),
                                               par.get("wave_number_rad"))
            except procedencia.ErrorNumeroOnda:
                onda = {"wave_number_pi": None, "wave_number_rad": None}
            th = par.get("thetas") or [None] * 3
            sep = "; " if idioma == "es" else ", "
            txt = r.t(
                "Se ajustó al VOI una microestructura espinodal {c}: el "
                "conjunto de nivel de un campo aleatorio gaussiano de N = {nw} "
                "ondas planas cuyas direcciones se muestrean {esq} conos de "
                "semiángulos θ = ({th})°, con número de onda β = {bpi}π "
                "({brad} rad) sobre el cubo unidad, umbralizado a la densidad "
                "relativa ρ = {rho}. La estructura se generó a {res}³ vóxeles "
                "con semilla {sem} y se escaló al lado del VOI.",
                "A spinodoid microstructure {c} was fitted to the VOI: the "
                "level set of a Gaussian random field of N = {nw} plane waves "
                "whose directions are sampled {esq} cones of half-angles "
                "θ = ({th})°, with wave number β = {bpi}π ({brad} rad) on the "
                "unit cube, thresholded at relative density ρ = {rho}. The "
                "structure was generated at {res}³ voxels with seed {sem} and "
                "scaled to the side of the VOI.",
                c=r.c("kumar2020", "soyarslan2018"),
                nw=_i(par.get("num_waves")),
                esq=(r.t("por rechazo dentro de la unión de",
                         "by rejection within the union of")
                     if par.get("esquema", "rechazo") == "rechazo"
                     else r.t("de forma equitativa entre",
                              "evenly across")),
                th=sep.join(r.n(t, ".4g") for t in th),
                bpi=r.n(onda["wave_number_pi"], ".6g"),
                brad=r.n(onda["wave_number_rad"], ".6g"),
                rho=r.n(par.get("densidad"), ".4g"),
                res=_i(par.get("resolucion")), sem=_i(par.get("semilla")))
        else:
            est = par.get("estiramiento") or [None] * 3
            sep = "; " if idioma == "es" else ", "
            txt = r.t(
                "Se ajustó al VOI una red dual de una teselación de Delaunay "
                "(dual-lattice) {c} con {cel} celdas por lado, estiramiento "
                "({est}), irregularidad {irr} y densidad relativa ρ = {rho}, "
                "generada a {res}³ vóxeles con semilla {sem} y escalada al "
                "lado del VOI.",
                "A dual lattice of a Delaunay tessellation {c} was fitted to "
                "the VOI, with {cel} cells per side, stretch ({est}), "
                "irregularity {irr} and relative density ρ = {rho}, generated "
                "at {res}³ voxels with seed {sem} and scaled to the side of "
                "the VOI.",
                c=r.c("vafaeefar2022"), cel=r.n(par.get("celdas"), ".4g"),
                est=sep.join(r.n(e, ".3g") for e in est),
                irr=r.n(par.get("irregularidad"), ".3g"),
                rho=r.n(par.get("densidad"), ".4g"),
                res=_i(par.get("resolucion")), sem=_i(par.get("semilla")))

        obj = rec.get("objetivo") or {}
        pesos = obj.get("pesos")
        txt += " " + r.t(
            "Los parámetros se eligieron con la búsqueda escalonada de spinpy "
            "minimizando un error morfométrico ponderado{pes}; el error del "
            "candidato devuelto fue {err}.",
            "Parameters were chosen by spinpy's staged search minimising a "
            "weighted morphometric error{pes}; the error of the returned "
            "candidate was {err}.",
            pes=(r.t(" con pesos modificados ({p})", " with modified weights "
                     "({p})", p=json.dumps(pesos, ensure_ascii=False))
                 if pesos else
                 r.t(" con los pesos por omisión", " with the default weights")),
            err=r.n(rec.get("error"), ".4g"))
        inc = rec.get("incertidumbre") or {}
        K = _i(inc.get("K")) or 0
        if K >= K_MIN:
            e = inc.get("error") or {}
            s = inc.get("suelo_autoconsistente") or {}
            txt += " " + r.t(
                "Su incertidumbre se estimó con K = {K} realizaciones de "
                "semillas nuevas, distintas de la de búsqueda: error {m} ± {sd} "
                "(media ± DE), frente a un suelo autoconsistente de {f}.",
                "Its uncertainty was estimated from K = {K} realisations with "
                "fresh seeds, different from the search seed: error {m} ± {sd} "
                "(mean ± SD), against a self-consistent floor of {f}.",
                K=K, m=r.n(e.get("media"), ".4g"), sd=r.n(e.get("sd"), ".3g"),
                f=r.n(s.get("media"), ".4g"))
        else:
            txt += " " + r.t(
                "El error corresponde a una única realización del generador "
                "estocástico.",
                "The error corresponds to a single realisation of the "
                "stochastic generator.")
        if rec.get("alineacion"):
            txt += " " + r.t(
                "La orientación del candidato se impuso a partir del marco del "
                "tensor MIL del VOI; la matriz de rotación usada consta en el "
                "paquete de reproducción.",
                "The candidate's orientation was imposed from the frame of the "
                "VOI's MIL tensor; the rotation matrix used is stored in the "
                "reproduction package.")
        P.append(txt)

    tbth_v = _f(m_voi or {}, "TbTh")
    lado = _lado_voi(doc)
    for suf in ("", "_dual"):
        rec = res.get("elastico" + suf)
        if not isinstance(rec, dict):
            continue
        n = _i(rec.get("resolucion"))
        pe = rec.get("por_estructura") or {}
        solver = next((v.get("solver") for v in pe.values()
                       if _f(v, "Ex") is not None), "")
        amg = "AMG" in str(solver)
        txt = r.t(
            "El tensor de rigidez efectivo se obtuvo por homogeneización "
            "periódica {c} sobre una rejilla de {n}³ elementos hexaédricos "
            "lineales, con tejido isótropo (E_s = {E} GPa, ν_s = {nu}) y los "
            "poros modelados con rigidez 10⁻⁶·E_s. El sistema se resolvió con "
            "{sol}{c2} y solo se aceptaron soluciones con residuo relativo "
            "≤ 10⁻⁶.",
            "The effective stiffness tensor was obtained by periodic "
            "homogenisation {c} on a grid of {n}³ linear hexahedral elements, "
            "with isotropic tissue (E_s = {E} GPa, ν_s = {nu}) and pores "
            "modelled with stiffness 10⁻⁶·E_s. The system was solved with "
            "{sol}{c2}, and only solutions with a relative residual ≤ 10⁻⁶ "
            "were accepted.",
            c=r.c("andreassen2014"), n=n,
            E=r.n((_f(rec, "E_s_Pa") or np.nan) / 1e9, ".4g"),
            nu=r.n(rec.get("nu_s"), ".3g"),
            sol=_solver_txt(solver, idioma),
            c2=(" " + r.c("bell2022")) if amg else "")
        if tbth_v and lado and n:
            txt += " " + r.t(
                "A esa resolución el VOI tiene Tb.Th/h = {x}.",
                "At that resolution the VOI has Tb.Th/h = {x}.",
                x=r.n(tbth_v / (lado / n), ".3g"))
        P.append(txt)

    # Con las dos familias cada ensayo se describe UNA vez: si las dos se
    # ensayaron igual la frase sale identica y se deja una sola, y el rango de
    # la capa superficial se junta en uno.
    textos, n_sup = [], []
    for suf in ("", "_dual"):
        rec = res.get("resistencia" + suf)
        if not isinstance(rec, dict):
            continue
        ejes = sorted({e for v in (rec.get("por_estructura") or {}).values()
                       for e in (v.get("ejes") or {})})
        textos.append(r.t(
            "Se simuló un ensayo de compresión uniaxial en {ejes} sobre {n}³ "
            "elementos hexaédricos (E_s = {E} GPa, ν_s = {nu}; apoyo: «{ap}»), "
            "mallando solo el hueso que une las caras de carga. La carga de "
            "fallo se estimó con el criterio de Pistoia et al. {c}: fallo "
            "cuando el {fr} % del tejido supera una deformación efectiva del "
            "{eps} %, parámetros calibrados en radio distal humano.",
            "A uniaxial compression test was simulated along {ejes} on {n}³ "
            "hexahedral elements (E_s = {E} GPa, ν_s = {nu}; support: "
            "\"{ap}\"), meshing only the bone that connects the loaded faces. "
            "Failure load was estimated with the criterion of Pistoia et al. "
            "{c}: failure when {fr} % of the tissue exceeds an effective "
            "strain of {eps} %, parameters calibrated on the human distal "
            "radius.",
            ejes=", ".join(ejes) or "—", n=_i(rec.get("resolucion")),
            E=r.n((_f(rec, "E_s_Pa") or np.nan) / 1e9, ".4g"),
            nu=r.n(rec.get("nu_s"), ".3g"), ap=rec.get("apoyo") or "—",
            c=r.c("pistoia2002"), fr=r.n(FRAC_CRITICA * 100, ".3g"),
            eps=r.n(EPS_CRITICA * 100, ".3g")))

        # El pico de tension SOLO se declara si de verdad se midio sobre la
        # capa superficial. Es una frase de metodos, no un adorno: dice que
        # numero se esta dando y por que no es el maximo.
        n_sup += [x for x in (_f(pp, "vm_n_superficie")
                              for v in (rec.get("por_estructura") or {})
                              .values()
                              for pp in (v.get("ejes") or {}).values()) if x]
    for t in dict.fromkeys(textos):
        P.append(t)
    if n_sup:
        P.append(r.t(
            "El pico de tensión se informa como percentil 99 de la tensión "
            "de von Mises sobre la capa superficial —los elementos de "
            "hueso con al menos una cara en el vacío, entre {lo} y {hi} "
            "elementos según la estructura—, y no como máximo: sobre esta "
            "malla de vóxeles el máximo no converge con la resolución, "
            "mientras que ese percentil sí.",
            "The stress peak is reported as the 99th percentile of the von "
            "Mises stress over the surface layer —the bone elements with "
            "at least one face in the void, between {lo} and {hi} elements "
            "depending on the structure—, not as a maximum: on this voxel "
            "mesh the maximum does not converge with resolution, whereas "
            "that percentile does.",
            lo=r.n(min(n_sup), ".0f"), hi=r.n(max(n_sup), ".0f")))

    recs = [res.get("analisis_comparado" + suf) for suf in ("", "_dual")]
    recs = [x for x in recs if isinstance(x, dict)]
    textos = []
    for rec in recs:
        varios = len(recs) > 1
        textos.append(r.t(
            "Se realizó además un análisis mecánico comparado del VOI y "
            + ("de cada candidato" if varios else "el candidato")
            + " según el protocolo de {prot}, a {n}³ elementos, con "
            "E_s = {E}, ν_s = {nu}, una carga de {F} N y apoyo «{ap}».",
            "A comparative mechanical analysis of the VOI and "
            + ("each candidate" if varios else "the candidate")
            + " was also performed following the protocol of {prot}, at {n}³ "
            "elements, with E_s = {E}, ν_s = {nu}, a load of {F} N and support "
            "\"{ap}\".",
            prot=rec.get("protocolo") or "—", n=_i(rec.get("resolucion")),
            E=r.n(rec.get("E_s_Pa"), ".4g"), nu=r.n(rec.get("nu_s"), ".3g"),
            F=r.n(rec.get("carga_N"), ".4g"), ap=rec.get("apoyo") or "—"))
    for t in dict.fromkeys(textos):
        P.append(t)

    return P, list(r.orden)


# ---------------------------------------------------------------------------
# Reproduccion
# ---------------------------------------------------------------------------

def _paquete(doc, regenerar=True, campos=None):
    """(paquete, {familia: mascara regenerada}). Las mascaras sirven a las figuras.

    `campos`: dict que, si se pasa, recibe {familia: (campo, info)} del
    generador (figura 0).
    """
    from . import __version__
    res = doc.get("resultados") or {}
    voi = doc.get("voi") or {}
    paq = {
        "formato": "spinpy/reproduccion", "version": 1,
        "procedencia": {"spinpy": __version__,
                        "python": platform.python_version(),
                        "numpy": np.__version__,
                        "fecha": time.strftime("%Y-%m-%dT%H:%M:%S")},
        "voi": {"nombre": voi.get("nombre"), "ruta": voi.get("ruta"),
                "sha256": voi.get("sha256"), "forma": voi.get("forma"),
                "spacing_mm": voi.get("spacing_mm")},
        "ajustes": {},
    }
    mascaras = {}
    for clave, fam in (("ajuste", "spinodoide"), ("ajuste_dual", "dual-lattice")):
        rec = res.get(clave)
        if not isinstance(rec, dict) or not rec.get("parametros"):
            continue
        par = dict(rec["parametros"])
        par.setdefault("familia", fam)
        a = {"parametros": par, "procedencia": rec.get("procedencia"),
             "sha256_mascara": None, "forma": None, "BVTV": None}
        if regenerar:
            try:
                BW, campo, info = _regenerar(par, completo=True)
                if campos is not None:
                    campos[fam] = (campo, info)
                a.update(sha256_mascara=huella_mascara(BW),
                         forma=[int(s) for s in BW.shape],
                         BVTV=float(BW.mean()))
                mascaras[fam] = BW
            except procedencia.ErrorNumeroOnda as e:
                a["error"] = f"onda_inconsistente: {e}"
        paq["ajustes"][fam] = a
    return procedencia.serializable(paq), mascaras


def paquete_reproduccion(doc, regenerar=True):
    """Parametros exactos, huella del VOI y huella de cada estructura ajustada."""
    return _paquete(doc, regenerar)[0]


def verificar_reproduccion(paq, ruta_voi=None):
    """Regenera cada estructura y compara huellas bit a bit.

    Devuelve {"ok", "ajustes": {familia: {"ok", "esperado", "obtenido",
    "motivo"}}, "voi": {...} | None, "avisos": [codigos]}. `ok` es False si
    no habia nada que comprobar: no se da por reproducido lo que no se probo.
    """
    from . import __version__
    out = {"ok": True, "ajustes": {}, "voi": None, "avisos": []}
    guardada = (paq.get("procedencia") or {}).get("spinpy")
    if guardada != __version__:
        out["avisos"].append({"codigo": "version_distinta",
                              "guardada": guardada, "instalada": __version__})
    comprobados = 0
    for fam, a in (paq.get("ajustes") or {}).items():
        esperado = a.get("sha256_mascara")
        if not esperado:
            out["ajustes"][fam] = {"ok": None, "motivo": "sin_huella"}
            continue
        try:
            obtenido = huella_mascara(_regenerar(a["parametros"]))
            ok = obtenido == esperado
            out["ajustes"][fam] = {"ok": ok, "esperado": esperado,
                                   "obtenido": obtenido, "motivo": ""}
        except procedencia.ErrorNumeroOnda as e:
            ok = False
            out["ajustes"][fam] = {"ok": False, "esperado": esperado,
                                   "obtenido": None,
                                   "motivo": f"onda_inconsistente: {e}"}
        comprobados += 1
        out["ok"] = out["ok"] and ok
    if ruta_voi is not None:
        esperado = (paq.get("voi") or {}).get("sha256")
        obtenido = sha256_archivo(ruta_voi)
        ok = bool(esperado) and obtenido == esperado
        out["voi"] = {"ok": ok, "esperado": esperado, "obtenido": obtenido}
        comprobados += 1
        out["ok"] = out["ok"] and ok
    if comprobados == 0:
        out["ok"] = False
        out["avisos"].append({"codigo": "nada_que_verificar"})
    return out


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Comparacion VOI / candidatos
# ---------------------------------------------------------------------------

# Metricas de la tabla comparativa: (clave, etiqueta, unidad, en el ajuste).
# Las seis primeras son las que el objetivo del ajuste intenta igualar; el
# resto se miden despues y dicen lo que el ajuste NO controla.
COMPARA_MORF = (
    ("BVTV", "BV/TV", "", True), ("BSBV", "BS/BV", "mm⁻¹", True),
    ("TbTh", "Tb.Th", "mm", True), ("TbSp", "Tb.Sp", "mm", True),
    ("TbN", "Tb.N", "mm⁻¹", True), ("DA", "DA", "", True),
    ("ConnD", "Conn.D", "mm⁻³", False), ("SMI", "SMI", "", False),
    ("PoDm", "Po.Dm", "mm", False), ("EF_frac_placa", "EF placas", "",
                                     False),
)
# Metricas con signo: su diferencia relativa no significa nada (un SMI de
# -0,97 frente a +1,01 daria «+204 %»). Se da la diferencia absoluta.
DIF_ABSOLUTA = {"SMI"}
FAMILIAS_INF = (("spinodoide", "", "morfometria_spin"),
                ("dual-lattice", "_dual", "morfometria_dual"))
MARCA = {CITABLE: "", RESERVAS: "*", NO_CITABLE: "†"}


def _estado_de(items, bloque, magnitud, est, eje=None, clave=None):
    """El peor estado de los items que casan, o None si no hay ninguno."""
    peor = None
    for it in items:
        if (it["bloque"] != bloque or it["magnitud"] != magnitud
                or it["estructura"] != est or it.get("eje") != eje):
            continue
        if clave is not None and clave not in (it.get("valor") or {}):
            continue
        if peor is None or _ORDEN[it["estado"]] > _ORDEN[peor]:
            peor = it["estado"]
    return peor


def _dif(v, c):
    if v is None or c is None or v == 0:
        return None
    return 100.0 * (c - v) / abs(v)


def seccion_comparacion(doc, items, idioma="es"):
    """Tabla VOI frente a cada familia ajustada: morfometria y mecanica.

    Una columna por estructura presente y, junto a cada candidato, su
    diferencia relativa con el VOI. Deliberadamente NO compara el error del
    ajuste entre familias: el objetivo del dual-lattice suma un termino DA2
    que el del spinodoide no tiene, asi que sus cifras no estan en la misma
    escala. Lo que se resume es la media de |diferencia| en las seis metricas
    que las dos familias intentan igualar, calculada aqui y declarada como tal.

    Devuelve lineas de Markdown, o [] si no hay ningun candidato.
    """
    i = _idx(idioma)
    T = lambda es, en: es if i == 0 else en      # noqa: E731
    fmt = _Formateador(idioma)
    res = doc.get("resultados") or {}
    m_voi = doc.get("morfometria_voi") or {}
    fams = [(f, suf, doc.get(k)) for f, suf, k in FAMILIAS_INF
            if doc.get(k) or isinstance(res.get("ajuste" + suf), dict)]
    if not fams or not m_voi:
        return []
    nombres = {"spinodoide": T("Spinodoide", "Spinodoid"),
               "dual-lattice": "Dual-lattice"}

    def num(x, espec=".4g"):
        return "—" if x is None else fmt.format_field(float(x), espec)

    def pct(d):
        return "—" if d is None else fmt.format_field(d, "+.1f") + " %"

    cab = ["", "VOI"]
    for f, _s, _m in fams:
        cab += [nombres[f], "Δ"]
    L = []

    # -- morfometria --
    L.append("### " + T("Morfometría", "Morphometry"))
    L.append("")
    L.append("| " + " | ".join(cab) + " |")
    L.append("|" + "---|" * len(cab))
    medias = {f: [] for f, _s, _m in fams}
    mas_cerca = {f: 0 for f, _s, _m in fams}
    n_comparadas = 0
    for clave, etq, uni, en_ajuste in COMPARA_MORF:
        v = _f(m_voi, clave)
        vals = [(f, _f(m or {}, clave)) for f, _s, m in fams]
        if v is None and all(c is None for _f_, c in vals):
            continue
        fila = [etq + (f" [{uni}]" if uni else ""), num(v)]
        difs = {}
        for f, c in vals:
            if clave in DIF_ABSOLUTA:
                d = None if (v is None or c is None) else c - v
                difs[f] = d
                fila += [num(c), "—" if d is None
                         else fmt.format_field(d, "+.3g")]
                continue
            d = _dif(v, c)
            difs[f] = d
            fila += [num(c), pct(d)]
            if en_ajuste and d is not None:
                medias[f].append(abs(d))
        L.append("| " + " | ".join(fila) + " |")
        validas = {f: abs(d) for f, d in difs.items() if d is not None}
        if len(fams) > 1 and len(validas) == len(fams):
            n_comparadas += 1
            mejor = min(validas.values())
            for f, d in validas.items():
                if d == mejor:
                    mas_cerca[f] += 1
    L.append("")
    resumen_morf = []
    for f, _s, _m in fams:
        if medias[f]:
            resumen_morf.append(T(
                f"{nombres[f]}: media de |Δ| {num(np.mean(medias[f]), '.1f')} % "
                f"en {len(medias[f])} métricas del ajuste",
                f"{nombres[f]}: mean |Δ| {num(np.mean(medias[f]), '.1f')} % "
                f"over {len(medias[f])} fitted metrics"))
    if resumen_morf:
        L.append("; ".join(resumen_morf) + ".")
        L.append("")

    # -- mecanica --
    filas = []

    def voi_y_cands(bloque_res, extraer):
        """[(valor VOI, estado VOI), (valor, estado, resolucion) por familia]"""
        voi_val, voi_est, cands, resol = None, None, {}, {}
        for f, suf, _m in fams:
            rec = res.get(bloque_res + suf)
            if not isinstance(rec, dict):
                cands[f] = (None, None)
                continue
            resol[f] = _i(rec.get("resolucion"))
            pe = rec.get("por_estructura") or {}
            if voi_val is None and "voi" in pe:
                voi_val, voi_est = extraer(pe["voi"], "voi")
            cands[f] = (extraer(pe["spin"], f) if "spin" in pe
                        else (None, None))
        return voi_val, voi_est, cands, resol

    resoluciones = {}

    def anadir(etq, bloque_res, extraer):
        v, ev, cands, resol = voi_y_cands(bloque_res, extraer)
        if v is None and all(c[0] is None for c in cands.values()):
            return
        for f, n in resol.items():
            resoluciones.setdefault(bloque_res, set()).add(n)
        fila = [etq, num(v) + MARCA.get(ev, "")]
        for f, _s, _m in fams:
            c, ec = cands.get(f, (None, None))
            fila += [num(c) + MARCA.get(ec, ""), pct(_dif(v, c))]
        filas.append(fila)

    for comp in ("Ex", "Ey", "Ez"):
        anadir(f"E{comp[1]} " + T("(homogeneización) [MPa]",
                                  "(homogenisation) [MPa]"),
               "elastico",
               lambda d, est, comp=comp: (
                   None if _f(d, comp) is None else _f(d, comp) / 1e6,
                   _estado_de(items, "elastico", "elastico", est)))
    ejes = sorted({e for f, suf, _m in fams
                   for v in ((res.get("resistencia" + suf) or {})
                             .get("por_estructura") or {}).values()
                   for e in (v.get("ejes") or {})})
    for eje in ejes:
        anadir(T(f"E_app ensayo, eje {eje} [MPa]",
                 f"E_app test, {eje} axis [MPa]"), "resistencia",
               lambda d, est, eje=eje: (
                   None if _f((d.get("ejes") or {}).get(eje), "E_app") is None
                   else _f(d["ejes"][eje], "E_app") / 1e6,
                   _estado_de(items, "resistencia", "E_app", est, eje)))
        anadir(T(f"σ fallo Pistoia, eje {eje} [MPa]",
                 f"Pistoia failure σ, {eje} axis [MPa]"), "resistencia",
               lambda d, est, eje=eje: (
                   None if _f((d.get("ejes") or {}).get(eje),
                              "sigma_fallo") is None
                   else _f(d["ejes"][eje], "sigma_fallo") / 1e6,
                   _estado_de(items, "resistencia", "fallo", est, eje)))
    anadir(T("E_app protocolo comparado [MPa]",
             "E_app compared protocol [MPa]"), "analisis_comparado",
           lambda d, est: (None if _f(d, "E_app") is None
                           else _f(d, "E_app") / 1e6,
                           _estado_de(items, "analisis_comparado", "E_app",
                                      est)))
    anadir(T("σ von Mises p99 capa superficial, 100 N [MPa]",
             "von Mises σ p99 surface layer, 100 N [MPa]"),
           "analisis_comparado",
           lambda d, est: (_f(d, "vm_p99_superficie"),
                           _estado_de(items, "analisis_comparado",
                                      "vm_p99_superficie", est,
                                      clave="vm_p99_superficie")))
    if filas:
        L.append("### " + T("Mecánica", "Mechanics"))
        L.append("")
        L.append("| " + " | ".join(cab) + " |")
        L.append("|" + "---|" * len(cab))
        for fila in filas:
            L.append("| " + " | ".join(fila) + " |")
        L.append("")
        L.append(T("Marcas: * con reservas, † no citable (motivo en la "
                   "sección 3). Δ es la diferencia relativa con el VOI; en "
                   "el SMI, que tiene signo, la diferencia absoluta.",
                   "Marks: * with caveats, † not citable (reason in section "
                   "3). Δ is the relative difference from the VOI; for the "
                   "SMI, which is signed, the absolute difference."))
        L.append("")
        distintas = [b for b, ns in resoluciones.items() if len(ns) > 1]
        if distintas:
            L.append(T(
                "**Atención:** las familias se ensayaron a resoluciones "
                "distintas en " + ", ".join(distintas) + "; esas Δ no son "
                "comparables entre familias.",
                "**Caution:** the families were tested at different "
                "resolutions in " + ", ".join(distintas) + "; those Δ are "
                "not comparable across families."))
            L.append("")

    # -- orientacion --
    ori = []
    for f, suf, _m in fams:
        o = (res.get("analisis_comparado" + suf) or {}).get("orientacion")
        if isinstance(o, dict) and _f(o, "angulo_deg") is not None:
            ori.append(T(f"{nombres[f]} a {num(o['angulo_deg'], '.0f')}°",
                         f"{nombres[f]} at {num(o['angulo_deg'], '.0f')}°"))
    if ori:
        L.append(T("Ángulo entre el eje principal (MIL) del VOI y el del "
                   "candidato, medido sobre la malla ensayada: ",
                   "Angle between the principal (MIL) axis of the VOI and "
                   "that of the candidate, measured on the tested mesh: ")
                 + "; ".join(ori) + T(". Por encima de 30° la comparación "
                                      "mecánica mide sobre todo la "
                                      "desalineación.",
                                      ". Above 30° the mechanical comparison "
                                      "mostly measures the misalignment."))
        L.append("")

    # -- en palabras sencillas --
    if len(fams) > 1 and n_comparadas:
        partes = [T(f"{nombres[f]} está más cerca del VOI en {mas_cerca[f]}",
                    f"{nombres[f]} is closer to the VOI in {mas_cerca[f]}")
                  for f, _s, _m in fams]
        L.append("> " + T(
            f"**En palabras sencillas.** De las {n_comparadas} métricas "
            "morfométricas medidas en las dos familias, "
            + " y ".join(partes) + ". Parecerse en la morfometría no "
            "garantiza parecerse en la rigidez: compárese la tabla de "
            "mecánica. El error de cada ajuste no se compara entre familias "
            "porque el del dual-lattice incluye un término de anisotropía "
            "(DA2) que el del spinodoide no tiene.",
            f"**In plain words.** Of the {n_comparadas} morphometric metrics "
            "measured in both families, " + " and ".join(partes) + ". "
            "Resembling the VOI in morphometry does not guarantee resembling "
            "it in stiffness: see the mechanics table. The fit errors are "
            "not compared across families because the dual-lattice one "
            "includes an anisotropy term (DA2) that the spinodoid one lacks."))
        L.append("")
    return L


def informe_markdown(doc, items, paq, idioma="es", nombre_paquete=None,
                     figuras=None, galeria=None):
    """Markdown del informe. `figuras`: [(ruta relativa, pie)] ya en su idioma;
    `galeria`: {"n", "estilos", "vistas"} de los renders de `figuras/3d/`."""
    i = _idx(idioma)
    T = lambda es, en: es if i == 0 else en      # noqa: E731
    red = _Redactor(idioma)
    P, refs = parrafos_metodos(doc, idioma, red)
    # Las citas de la seccion de modelos continuan la numeracion de los
    # metodos, y `refs` es la misma lista, asi que las recoge sola. Se redacta
    # aqui, antes que nada, para saber si existe y en que numero cae.
    sec_modelos = 5 if figuras else 4
    from .informe_modelos import seccion_modelos
    modelos = seccion_modelos(doc, idioma, red, sec_modelos)
    voi = doc.get("voi") or {}
    from . import __version__
    L = []
    L.append("# " + T("Informe para publicación", "Publication report")
             + (f" — {voi['nombre']}" if voi.get("nombre") else ""))
    L.append("")
    L.append(T(f"Generado con spinpy {__version__} el "
               f"{time.strftime('%Y-%m-%d %H:%M')}.",
               f"Generated with spinpy {__version__} on "
               f"{time.strftime('%Y-%m-%d %H:%M')}."))
    L.append("")

    L.append("## 1. " + T("Métodos", "Methods"))
    L.append("")
    L.append(T("Texto listo para adaptar a la sección de métodos. Cada frase "
               "sale de los valores que se usaron en la sesión; los números "
               "entre corchetes remiten a la lista de referencias.",
               "Text ready to adapt into the methods section. Every sentence "
               "comes from the values used in the session; bracketed numbers "
               "refer to the reference list.")
             + (T(f" Las ecuaciones de cada cálculo, sus condiciones de "
                  f"contorno y tolerancias están en la sección {sec_modelos}.",
                  f" The equations of every computation, their boundary "
                  f"conditions and tolerances are in section {sec_modelos}.")
                if modelos else ""))
    L.append("")
    for p in P:
        L.append(p)
        L.append("")

    # Siempre presente, aunque vacia: la numeracion de las secciones es fija
    # porque los pies de figura remiten a «la sección 3».
    comp = seccion_comparacion(doc, items, idioma)
    dos = bool(doc.get("morfometria_spin")) and bool(
        doc.get("morfometria_dual"))
    L.append("## 2. " + (T("Comparación VOI / Spinodoide / Dual-lattice",
                           "Comparison VOI / Spinodoid / Dual-lattice")
                         if dos else T("Comparación con el VOI",
                                       "Comparison with the VOI")))
    L.append("")
    L += comp or [T("No hay estructuras ajustadas que comparar con el VOI.",
                    "There are no fitted structures to compare with the "
                    "VOI."), ""]
    L.append("## 3. " + T("Citabilidad de los resultados",
                          "Citability of the results"))
    L.append("")
    c = resumen(items)
    L.append(T(f"**{c[CITABLE]}** citables, **{c[RESERVAS]}** con reservas, "
               f"**{c[NO_CITABLE]}** no citables.",
               f"**{c[CITABLE]}** citable, **{c[RESERVAS]}** with caveats, "
               f"**{c[NO_CITABLE]}** not citable."))
    L.append("")
    L.append("| " + " | ".join(T(*h) for h in (
        ("Resultado", "Result"), ("Valor", "Value"), ("Estado", "Status"),
        ("Motivo", "Reason"))) + " |")
    L.append("|---|---|---|---|")
    for it in sorted(items, key=lambda x: -_ORDEN[x["estado"]]):
        L.append("| " + " | ".join((
            etiqueta_item(it, idioma), _valor_txt(it["valor"], idioma),
            ESTADOS[it["estado"]][i], texto_motivos(it, idioma))) + " |")
    L.append("")
    L.append("> " + T(
        "**En palabras sencillas.** Un número puede salir del programa sin ningún "
        "error y aun así no significar lo que parece. «Con reservas» quiere "
        "decir que se puede usar si se explica su limitación en el artículo; "
        "«no citable», que no debe presentarse como resultado.",
        "**In plain words.** A number can come out of the program without any "
        "error and still not mean what it seems. \"With caveats\" means it "
        "can be used if its limitation is explained in the paper; \"not "
        "citable\" means it must not be presented as a result."))
    L.append("")

    sec = 4
    if figuras:
        L.append(f"## {sec}. " + T("Figuras", "Figures"))
        L.append("")
        L.append(T("Los archivos están en `figuras/`: gráficos en PNG a 600 ppp "
                   "y en PDF vectorial; estructuras 3D en PNG de 2400 px de "
                   "lado.",
                   "The files are in `figuras/`: charts as 600 dpi PNG and "
                   "vector PDF; 3D structures as 2400 px PNG."))
        L.append("")
        if galeria and galeria.get("n"):
            from . import figuras as F
            L.append(T(
                f"En `figuras/3d/` hay {galeria['n']} renders más, a 1600 px, "
                "para elegir otra presentación sin volver a calcular: estilos "
                + ", ".join(F.ESTILOS_3D[e]["nombre"][0]
                            for e in galeria["estilos"])
                + "; vistas "
                + ", ".join(F.VISTAS_3D[v][2][0] for v in galeria["vistas"])
                + ". Cada archivo se llama `estructura_estilo_vista.png`.",
                f"`figuras/3d/` holds {galeria['n']} more renders, at 1600 px, "
                "to choose another presentation without recomputing: styles "
                + ", ".join(F.ESTILOS_3D[e]["nombre"][1]
                            for e in galeria["estilos"])
                + "; views "
                + ", ".join(F.VISTAS_3D[v][2][1] for v in galeria["vistas"])
                + ". Each file is named `structure_style_view.png`."))
            L.append("")
        for ruta, pie in figuras:
            L.append(f"![{pie}]({ruta})")
            L.append("")
        sec += 1

    if modelos:
        L += modelos
        sec += 1

    L.append(f"## {sec}. " + T("Referencias", "References"))
    L.append("")
    for k, clave in enumerate(refs, 1):
        txt, doi = REFERENCIAS[clave]
        L.append(f"{k}. {txt}" + (f" https://doi.org/{doi}" if doi else ""))
    L.append("")

    L.append(f"## {sec + 1}. " + T("Cómo citar el software",
                                   "How to cite the software"))
    L.append("")
    autores = ", ".join(
        a["apellidos"] + " " + "".join(p[0] for p in
                                       a["nombre"].replace(".", " ").split())
        for a in AUTORES)
    L.append(T(f"{autores}. spinpy, versión {__version__}. Repositorio: "
               f"{REPOSITORIO}. DOI: {DOI_ARCHIVO}.",
               f"{autores}. spinpy, version {__version__}. Repository: "
               f"{REPOSITORIO}. DOI: {DOI_ARCHIVO}."))
    L.append("")
    for a in AUTORES:
        L.append(f"- {a['nombre']} {a['apellidos']} — "
                 f"https://orcid.org/{a['orcid']}")
    L.append("")

    L.append(f"## {sec + 2}. " + T("Reproducción", "Reproduction"))
    L.append("")
    nombre_paquete = nombre_paquete or "reproduccion.json"
    sha_voi = (paq.get("voi") or {}).get("sha256")
    L.append(T(
        f"`{nombre_paquete}` contiene los parámetros exactos de cada ajuste "
        f"—incluidas la matriz de rotación y las dos lecturas del número de "
        f"onda—, la huella SHA-256 del VOI "
        f"({'`' + sha_voi + '`' if sha_voi else 'no disponible: el VOI no se cargó desde un archivo único'}) "
        f"y la de cada estructura regenerada:",
        f"`{nombre_paquete}` holds the exact parameters of each fit "
        f"—including the rotation matrix and both readings of the wave "
        f"number—, the SHA-256 fingerprint of the VOI "
        f"({'`' + sha_voi + '`' if sha_voi else 'not available: the VOI was not loaded from a single file'}) "
        f"and that of every regenerated structure:"))
    L.append("")
    L.append("| " + T("Familia", "Family") + " | " + T("Rejilla", "Grid")
             + " | SHA-256 |")
    L.append("|---|---|---|")
    for fam, a in (paq.get("ajustes") or {}).items():
        forma = "×".join(str(s) for s in (a.get("forma") or [])) or "—"
        L.append(f"| {ESTRUCTURAS[fam][i]} | {forma} | "
                 f"`{a.get('sha256_mascara') or '—'}` |")
    L.append("")
    L.append(T("Para comprobar que otra instalación obtiene exactamente las "
               "mismas estructuras:",
               "To check that another installation obtains exactly the same "
               "structures:"))
    L.append("")
    L.append("```bash")
    L.append(f"spinpy-informe --verificar {nombre_paquete} --voi RUTA_DEL_VOI")
    L.append("```")
    L.append("")
    L.append("> " + T(
        "**En palabras sencillas.** El generador usa números aleatorios, pero con "
        "la misma semilla y la misma versión produce siempre la misma "
        "estructura. La huella SHA-256 es como una firma de esa estructura: "
        "si otra persona obtiene la misma firma, obtuvo exactamente la misma "
        "estructura, vóxel a vóxel.",
        "**In plain words.** The generator uses random numbers, but with the same "
        "seed and version it always produces the same structure. The SHA-256 "
        "fingerprint is like a signature of that structure: if someone else "
        "gets the same signature, they got exactly the same structure, voxel "
        "by voxel."))
    L.append("")
    return "\n".join(L)


# ---------------------------------------------------------------------------
# Carpeta del informe
# ---------------------------------------------------------------------------

# Nombres FIJOS dentro de la carpeta que elige el usuario: el informe cita el
# paquete por su nombre y el visor los usa para avisar antes de sobrescribir.
ARCHIVOS = {
    "md_es": "informe_publicacion_es.md",
    "md_en": "publication_report_en.md",
    "pdf_es": "informe_publicacion_es.pdf",
    "pdf_en": "publication_report_en.pdf",
    "paquete": "reproduccion.json",
    "sesion": "resultados_sesion.json",
}
CARPETA_FIGURAS = "figuras"

AVISOS = {
    "pdf_sin_reportlab": (
        "no se generaron los PDF: falta reportlab en este entorno",
        "the PDFs were not generated: reportlab is missing in this environment"),
    "figuras_sin_matplotlib": (
        "no se generaron las figuras de datos: falta matplotlib",
        "data figures were not generated: matplotlib is missing"),
    "figuras_3d": (
        "no se generaron las figuras 3D ({detalle})",
        "3D figures were not generated ({detalle})"),
    "voi_no_disponible": (
        "sin render 3D del VOI: no estaba cargado y su ruta no se pudo leer",
        "no 3D render of the VOI: it was not loaded and its path could not be "
        "read"),
    "contexto_voi_no_disponible": (
        "figura 0 sin la pila de micro-CT: el VOI se recortó de una pila que "
        "ya no está en la ruta guardada o no se deja leer",
        "figure 0 without the micro-CT stack: the VOI was cropped from a "
        "stack that is no longer at the saved path or cannot be read"),
}

PIES = {
    "morfometria": (
        "Morfometría de cada estructura ajustada relativa al VOI (línea negra "
        "= 1). Cada barra es el cociente candidato/VOI; bajo cada métrica, su "
        "valor en el VOI.",
        "Morphometry of each fitted structure relative to the VOI (black line "
        "= 1). Each bar is the candidate/VOI ratio; below each metric, its "
        "value in the VOI."),
    "mecanica": (
        "Módulos elásticos efectivos por homogeneización periódica y módulo "
        "aparente en compresión. Las barras rayadas no son citables y las "
        "marcadas con * se citan con reservas (sección 3).",
        "Effective elastic moduli from periodic homogenisation and apparent "
        "compressive modulus. Hatched bars are not citable and those marked * "
        "carry caveats (section 3)."),
    "estructuras": (
        "Isosuperficies del VOI y de las estructuras ajustadas, a la misma "
        "escala física ({lado} mm de lado); vista {vista}, estilo {estilo}, "
        "proyección paralela.",
        "Isosurfaces of the VOI and of the fitted structures at the same "
        "physical scale ({lado} mm side); {vista} view, {estilo} style, "
        "parallel projection."),
    "vistas": (
        "Cada estructura desde varias direcciones (filas: estructuras; "
        "columnas: vistas), estilo {estilo}, proyección paralela y la misma "
        "escala en todas las celdas ({lado} mm de lado). Z es el eje axial.",
        "Each structure from several directions (rows: structures; columns: "
        "views), {estilo} style, parallel projection and the same scale in "
        "every cell ({lado} mm side). Z is the axial direction."),
    "secciones": (
        "Cortes centrales de la imagen binaria de cada estructura, hueso en "
        "negro, a la misma escala ({lado} mm de lado); barra de escala en "
        "rojo. Son los vóxeles tal cual se midieron, sin suavizar.",
        "Central sections of the binary image of each structure, bone in "
        "black, at the same scale ({lado} mm side); scale bar in red. These "
        "are the voxels exactly as measured, without smoothing."),
    "anisotropia": (
        "Anisotropía en tres planos. Arriba, longitud media de intercepción "
        "MIL(n) del tensor de fábrica (geometría); abajo, módulo de Young "
        "direccional E(n) del tensor de rigidez homogeneizado (mecánica). Si "
        "los lóbulos de las dos filas apuntan igual, la dirección rígida es "
        "la de la fábrica. Z es el eje axial.",
        "Anisotropy in three planes. Top, mean intercept length MIL(n) from "
        "the fabric tensor (geometry); bottom, directional Young's modulus "
        "E(n) from the homogenised stiffness tensor (mechanics). If the lobes "
        "of both rows point the same way, the stiff direction is the fabric "
        "one. Z is the axial direction."),
    "distribuciones": (
        "Distribución del espesor trabecular local y del tamaño de poro local "
        "(esferas inscritas) de cada estructura; línea discontinua, la media. "
        "Dos estructuras con el mismo Tb.Th medio pueden tener distribuciones "
        "muy distintas.",
        "Distribution of local trabecular thickness and local pore size "
        "(inscribed spheres) of each structure; dashed line, the mean. Two "
        "structures with the same mean Tb.Th can have very different "
        "distributions."),
    "von_mises": (
        "Tensión de von Mises sobre la superficie ósea bajo el protocolo del "
        "análisis comparado ({carga} N), vista {vista}, con la misma escala "
        "de color para todas las estructuras (percentiles 1 y 99 del tejido). "
        "Los picos locales no son citables: el valor que converge es el "
        "percentil 99 de la capa superficial (sección 3).",
        "Von Mises stress on the bone surface under the compared-analysis "
        "protocol ({carga} N), {vista} view, with the same colour scale for "
        "every structure (tissue 1st and 99th percentiles). Local peaks are "
        "not citable: the converging value is the surface-layer 99th "
        "percentile (section 3)."),
    "convergencia": (
        "Convergencia en malla del módulo aparente en compresión ({est}, eje "
        "{eje}). (a) E_app frente a la resolución efectiva Tb.Th/h (elementos "
        "por espesor trabecular); banda, ±3 % de la malla más fina; línea "
        "discontinua, extrapolación de Richardson, que solo se calcula si la "
        "serie es monótona y la densidad estable. (b) Cambio de la densidad y "
        "de la fracción portante frente a la malla más fina: por encima de "
        "±{deriva} % el remuestreo cambia la estructura y la serie ya no mide "
        "solo la discretización. Estado de citabilidad en la sección 3.",
        "Mesh convergence of the apparent compressive modulus ({est}, {eje} "
        "axis). (a) E_app against the effective resolution Tb.Th/h (elements "
        "per trabecular thickness); band, ±3 % of the finest mesh; dashed "
        "line, Richardson extrapolation, computed only when the series is "
        "monotonic and density stable. (b) Change of density and load-bearing "
        "fraction against the finest mesh: beyond ±{deriva} % resampling "
        "changes the structure and the series no longer measures "
        "discretisation alone. Citability status in section 3."),
    "incertidumbre": (
        "(a) Error morfométrico de K realizaciones nuevas del candidato "
        "ganador (semillas distintas de la de búsqueda), su media ± sd, el "
        "error de la búsqueda (×: la realización ganadora, optimista por "
        "construcción) y el suelo autoconsistente (error entre pares de "
        "realizaciones del mismo candidato; un error dentro del suelo no se "
        "distingue del ruido del generador). (b) Separación de cada métrica "
        "respecto del VOI en desviaciones típicas de las réplicas, escala "
        "lineal dentro de ±2 y logarítmica fuera; la banda gris (|z| < 2) es "
        "indistinguible. Con K = {K} la sd tiene un error relativo de "
        "~{err_sd} %: una |z| cercana a 2 no decide nada.",
        "(a) Morphometric error of K fresh realisations of the winning "
        "candidate (seeds other than the search seed), their mean ± sd, the "
        "search error (×: the winning realisation, optimistic by "
        "construction) and the self-consistent floor (error between pairs of "
        "realisations of the same candidate; an error inside the floor cannot "
        "be told from generator noise). (b) Separation of each metric from "
        "the VOI in replicate standard deviations, linear within ±2 and "
        "logarithmic beyond; the grey band (|z| < 2) is indistinguishable. "
        "With K = {K} the sd has a relative error of ~{err_sd} %: a |z| near "
        "2 decides nothing."),
    "vm_superficie": (
        "Fracción de la capa superficial (elementos de hueso con una cara "
        "hacia el poro; las seis caras del VOI no cuentan) cuya tensión de "
        "von Mises supera σ, bajo el protocolo del análisis comparado "
        "({carga} N). El punto marca el percentil 99 (posiciones (i − 0,5)/n, "
        "la convención de prctile), que es el valor citable. El eje llega "
        "al percentil 99,9: el máximo, que no converge con la malla, queda "
        "fuera a propósito. * con reservas, † no citable (sección 3).",
        "Fraction of the surface layer (bone elements with a face towards the "
        "pore; the six VOI faces do not count) whose von Mises stress exceeds "
        "σ, under the compared-analysis protocol ({carga} N). The dot marks "
        "the 99th percentile (plotting positions (i − 0.5)/n, the prctile "
        "convention), which is the citable value. The axis stops at the "
        "99.9th percentile: the maximum, which does not converge with the "
        "mesh, is left out on purpose. * with caveats, † not citable "
        "(section 3)."),
    "suavizado": (
        " Superficie suavizada con Taubin solo para la figura; todas las "
        "medidas usan la malla sin suavizar.",
        " Surface smoothed with Taubin for the figure only; every measurement "
        "uses the unsmoothed mesh."),
}


# Figura 0: el pie se arma por partes (`pie_metodo`) segun lo que la figura
# lleve de verdad, para no describir paneles que no estan.
PIES_METODO = {
    "cabecera": (
        "Método, de la imagen al candidato.",
        "Method, from image to candidate."),
    "voi_pila": (
        "Fila superior: rebanada {k} de {n} de la pila de micro-CT, la que "
        "pasa por el centro del VOI, {gris}, con el corte exacto del cubo en "
        "azul; la misma rebanada en la máscara segmentada que entró en la "
        "pila; la pila reducida por medias de bloques de {p}³ vóxeles (solo "
        "como contexto, no se mide sobre ella) con el cubo recortado en el "
        "marco de sus ejes principales; el VOI, de {lado} mm de lado, y a la "
        "derecha la morfometría que persigue el ajuste.",
        "Top row: slice {k} of {n} of the micro-CT stack, the one through the "
        "centre of the VOI, {gris}, with the exact section of the cube in "
        "blue; the same slice in the segmented mask that entered the stack; "
        "the stack reduced by {p}³-voxel block means (context only, nothing "
        "is measured on it) with the cube cropped in its principal-axis "
        "frame; the VOI, {lado} mm side, and on the right the morphometry "
        "the fit pursues."),
    "gris": ("en niveles de gris del archivo original",
             "in grey levels from the original file"),
    "sin_gris": ("ya binaria (el archivo original no se pudo releer)",
                 "already binary (the original file could not be re-read)"),
    "voi_archivo": (
        "Fila superior: corte central XZ del VOI, leído de un archivo ya "
        "recortado (la sesión no tiene la pila de micro-CT de la que salió), "
        "el VOI, de {lado} mm de lado, y a la derecha la morfometría que "
        "persigue el ajuste.",
        "Top row: central XZ section of the VOI, read from an already "
        "cropped file (the session does not hold the micro-CT stack it came "
        "from), the VOI, {lado} mm side, and on the right the morphometry the "
        "fit pursues."),
    "familias": (
        "Una fila por familia ajustada, con los parámetros y la semilla de "
        "esta sesión:",
        "One row per fitted family, with the parameters and seed of this "
        "session:"),
    "spinodoide": (
        "las direcciones de onda del generador (y sus opuestas) en los conos "
        "de semiángulo θ sobre los ejes de R, el campo gaussiano φ en las "
        "caras del cubo con la isolínea φ₀, su histograma con la N(0,1) "
        "discontinua y el umbral φ₀ que deja la fracción ρ, y el sólido "
        "{φ ≤ φ₀};",
        "the generator's wave directions (and their opposites) in the cones "
        "of half-angle θ about the axes of R, the Gaussian field φ on the "
        "faces of the cube with the φ₀ isoline, its histogram with the dashed "
        "N(0,1) and the threshold φ₀ that leaves the fraction ρ, and the "
        "solid {φ ≤ φ₀};"),
    "dual": (
        "la red dual 4-N de la triangulación de Delaunay (subcubo central, "
        "{sub} % del lado), la distancia d al esqueleto en las caras del cubo "
        "con la isolínea r, el umbral r como cuantil ρ de d, y el sólido "
        "{{d ≤ r}};",
        "the 4-N dual lattice of the Delaunay triangulation (central "
        "sub-cube, {sub} % of the side), the distance d to the skeleton on "
        "the faces of the cube with the r isoline, the threshold r as the "
        "ρ-quantile of d, and the solid {{d ≤ r}};"),
    "tablas": (
        "las tablas dan la morfometría de la sesión (sección 2) y su "
        "diferencia con el VOI.",
        "the tables give the session morphometry (section 2) and its "
        "difference from the VOI."),
}


def texto_aviso(aviso, idioma="es"):
    cod, _sep, detalle = str(aviso).partition(":")
    if cod not in AVISOS:
        return str(aviso)
    return AVISOS[cod][_idx(idioma)].format(detalle=detalle.strip())


def _espaciado(sp):
    sp = np.atleast_1d(np.asarray(sp, float)).ravel()
    return np.repeat(sp, 3) if sp.size == 1 else sp[:3]


def preparar(doc, carpeta, VOI=None, spacing=None, regenerar=True,
             cargar_voi=True, opciones_3d=None, campos_vm=None,
             contexto_voi=None):
    """Etapa 1 (hilo de trabajo): citabilidad, paquete y estructuras a dibujar.

    `opciones_3d`: estilos, vistas, principal y suavizado de las figuras 3D
    (`figuras.OPCIONES_3D`); None toma los de por omision. `campos_vm`:
    {estructura: (campo de von Mises en Pa, spacing)} del analisis comparado,
    para la figura 8; solo existe dentro del visor, no en un JSON.
    `contexto_voi`: `figura_metodo.contexto_pila` tomado al recortar el VOI
    (figura 0). Sin el, si el documento trae `voi.recorte` y la pila sigue en
    el disco, se rehace releyendola.
    """
    from . import figuras as F
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    avisos = []
    o = F.opciones_3d(opciones_3d)
    generados = {} if o["metodo"] else None
    paq, mascaras = _paquete(doc, regenerar=regenerar, campos=generados)
    rec = (doc.get("voi") or {}).get("recorte")
    if o["metodo"] and contexto_voi is None and rec and cargar_voi:
        from .figura_metodo import rehacer_contexto
        contexto_voi = rehacer_contexto(rec)
        if contexto_voi is None:
            avisos.append("contexto_voi_no_disponible")
    if VOI is None and cargar_voi:
        ruta = (doc.get("voi") or {}).get("ruta")
        try:
            from .io import leer_voi
            VOI, spacing = leer_voi(ruta)
        except Exception:
            VOI = None
            if doc.get("voi"):
                avisos.append("voi_no_disponible")
    lado = _lado_voi(doc)
    estructuras = []
    if VOI is not None:
        sp = _espaciado(spacing)
        estructuras.append(("voi", np.asarray(VOI, dtype=bool), sp))
        lado = lado or float(VOI.shape[0] * sp[0])
    for fam, BW in mascaras.items():
        estructuras.append((fam, BW, np.full(3, lado / BW.shape[0]
                                             if lado else 1.0)))
    return {"carpeta": carpeta, "doc": doc, "items": comprobar(doc),
            "paquete": paq, "estructuras": estructuras, "lado_mm": lado,
            "figuras": {"es": [], "en": []}, "archivos_figuras": [],
            "opciones_3d": opciones_3d, "galeria": None,
            "campos_vm": campos_vm, "generados": generados,
            "contexto_voi": contexto_voi,
            "avisos": avisos}


def figuras_datos(prep):
    """Etapa 2 (hilo de trabajo): figuras de datos en ES y EN."""
    try:
        from . import figuras as F
    except ImportError:
        prep["avisos"].append("figuras_sin_matplotlib")
        return prep
    d = prep["carpeta"] / CARPETA_FIGURAS
    o = F.opciones_3d(prep.get("opciones_3d"))
    dist = None           # se calcula una vez y sirve a los dos idiomas
    for idioma in ("es", "en"):
        i = _idx(idioma)
        for clave, nombre, hacer in (
                ("morfometria", ("fig1_morfometria", "fig1_morphometry")[i],
                 lambda dest: F.fig_morfometria(prep["doc"], dest, idioma)),
                ("mecanica", ("fig2_mecanica", "fig2_mechanics")[i],
                 lambda dest: F.fig_mecanica(prep["doc"], prep["items"], dest,
                                             idioma))):
            rutas = hacer(d / nombre)
            if rutas:
                prep["archivos_figuras"] += rutas
                prep["figuras"][idioma].append(
                    (f"{CARPETA_FIGURAS}/{Path(rutas[0]).name}",
                     PIES[clave][i]))
        # Cortes: sin VTK, asi que van aqui y no con los renders 3D.
        rutas = F.fig_secciones(prep["estructuras"],
                                d / ("fig5_secciones", "fig5_sections")[i],
                                idioma)
        if rutas:
            prep["archivos_figuras"] += rutas
            prep["figuras"][idioma].append(
                (f"{CARPETA_FIGURAS}/{Path(rutas[0]).name}",
                 PIES["secciones"][i].format(lado=_lado_txt(prep, idioma))))
        rutas = F.fig_anisotropia(prep["doc"], d / ("fig6_anisotropia",
                                                    "fig6_anisotropy")[i],
                                  idioma)
        if rutas:
            prep["archivos_figuras"] += rutas
            prep["figuras"][idioma].append(
                (f"{CARPETA_FIGURAS}/{Path(rutas[0]).name}",
                 PIES["anisotropia"][i]))
        _figuras_estudios(prep, d, idioma)
        if o["distribuciones"] and prep["estructuras"]:
            if dist is None:
                dist = F.distribuciones(prep["estructuras"])
            rutas = F.fig_distribuciones(
                dist, d / ("fig7_distribuciones", "fig7_distributions")[i],
                idioma)
            if rutas:
                prep["archivos_figuras"] += rutas
                prep["figuras"][idioma].append(
                    (f"{CARPETA_FIGURAS}/{Path(rutas[0]).name}",
                     PIES["distribuciones"][i]))
    return prep


def _figuras_estudios(prep, d, idioma):
    """Figuras 9-11: convergencia, incertidumbre y cola de von Mises.

    Salen solo del documento —sin volumenes—, asi que tambien se rehacen
    desde un JSON exportado. Cada una aparece solo si su etapa se corrio.
    """
    from . import figuras as F
    i = _idx(idioma)
    doc, items = prep["doc"], prep["items"]
    fm = _Formateador(idioma)
    res = doc.get("resultados") or {}

    def poner(rutas, pie):
        if rutas:
            prep["archivos_figuras"].extend(rutas)
            prep["figuras"][idioma].append(
                (f"{CARPETA_FIGURAS}/{Path(rutas[0]).name}", pie))

    rec = res.get("convergencia") or {}
    est = estructura_de(rec.get("estructura_codigo")
                        or rec.get("estructura"), doc)
    poner(F.fig_convergencia(doc, items, d / ("fig9_convergencia",
                                               "fig9_convergence")[i],
                             idioma),
          PIES["convergencia"][i].format(
              est=ESTRUCTURAS[est][i], eje=rec.get("eje", "z"),
              deriva=fm.format_field(100.0 * CONV_DERIVA_RHO_MAX, "g")))

    Ks = [_i(((res.get(c) or {}).get("incertidumbre") or {}).get("K"))
          for c in ("ajuste", "ajuste_dual")]
    Ks = [k for k in Ks if k]
    K = min(Ks) if Ks else K_MIN
    poner(F.fig_incertidumbre(doc, items, d / ("fig10_incertidumbre",
                                                "fig10_uncertainty")[i],
                              idioma),
          PIES["incertidumbre"][i].format(
              K=K, err_sd=fm.format_field(
                  100.0 / np.sqrt(2.0 * max(K - 1, 1)), ".0f")))

    carga = _f(res.get("analisis_comparado")
               or res.get("analisis_comparado_dual") or {}, "carga_N") or 100.0
    poner(F.fig_von_mises_superficie(
              doc, items, d / ("fig11_von_mises_superficie",
                               "fig11_von_mises_surface")[i], idioma),
          PIES["vm_superficie"][i].format(
              carga=fm.format_field(float(carga), ".0f")))


def _numero_figura(ruta):
    """fig10 va detras de fig9, no detras de fig1: orden numerico, no de
    texto. Lo que no empieza por figN va al final, en orden alfabetico."""
    nombre = Path(ruta).name
    digitos = ""
    for c in nombre[3:] if nombre.startswith("fig") else "":
        if not c.isdigit():
            break
        digitos += c
    return (int(digitos) if digitos else 10 ** 6, nombre)


def _lado_txt(prep, idioma):
    lado = prep["lado_mm"]
    return _Formateador(idioma).format_field(
        float(lado) if lado else float("nan"), ".3g")


def figuras_3d(prep):
    """Etapa 3: renders 3D. Llamar desde el hilo que puede crear OpenGL."""
    if not prep["estructuras"]:
        return prep
    try:
        from . import figuras as F
        d = prep["carpeta"] / CARPETA_FIGURAS
        o = F.opciones_3d(prep.get("opciones_3d"))
        r = F.renders_3d(prep["estructuras"], d, o)
        ests = [e[0] for e in prep["estructuras"]]
        sueltas = [r["principal"][e] for e in ests]
        prep["archivos_figuras"] += sueltas + list(r["galeria"].values())
        prep["galeria"] = {"n": len(r["galeria"]), "estilos": o["estilos"],
                           "vistas": o["vistas"]}
        for idioma, nombre, nombre4 in (
                ("es", "fig3_estructuras", "fig4_vistas"),
                ("en", "fig3_structures", "fig4_views")):
            i = _idx(idioma)
            ruta = F.panel_3d(sueltas, ests, d / nombre, idioma)
            prep["archivos_figuras"].append(ruta)
            suave = PIES["suavizado"][i] if o["suavizar"] else ""
            prep["figuras"][idioma].append((
                f"{CARPETA_FIGURAS}/{Path(ruta).name}",
                PIES["estructuras"][i].format(
                    lado=_lado_txt(prep, idioma),
                    estilo=F.ESTILOS_3D[o["estilo_principal"]]["nombre"][i],
                    vista=F.VISTAS_3D[o["vista_principal"]][2][i]) + suave))
            if len(o["vistas"]) > 1:
                ruta = F.panel_vistas(r["galeria"], ests,
                                      o["estilo_principal"], o["vistas"],
                                      d / nombre4, idioma)
                if ruta:
                    prep["archivos_figuras"].append(ruta)
                    prep["figuras"][idioma].append((
                        f"{CARPETA_FIGURAS}/{Path(ruta).name}",
                        PIES["vistas"][i].format(
                            lado=_lado_txt(prep, idioma),
                            estilo=F.ESTILOS_3D[o["estilo_principal"]]
                            ["nombre"][i]) + suave))
        campos = prep.get("campos_vm")
        if o["von_mises"] and campos:
            rutas_vm, clim = F.renders_von_mises(campos, d,
                                                 o["vista_principal"])
            if rutas_vm:
                prep["archivos_figuras"] += list(rutas_vm.values())
                carga = _f((prep["doc"].get("resultados") or {}).get(
                    "analisis_comparado") or (prep["doc"].get("resultados")
                                              or {}).get(
                    "analisis_comparado_dual") or {}, "carga_N") or 100.0
                for idioma, nombre in (("es", "fig8_von_mises"),
                                       ("en", "fig8_von_mises_en")):
                    i = _idx(idioma)
                    ruta = F.panel_von_mises(rutas_vm, clim, d / nombre,
                                             idioma, carga)
                    if ruta:
                        prep["archivos_figuras"].append(ruta)
                        prep["figuras"][idioma].append((
                            f"{CARPETA_FIGURAS}/{Path(ruta).name}",
                            PIES["von_mises"][i].format(
                                carga=_Formateador(idioma).format_field(
                                    float(carga), ".0f"),
                                vista=F.VISTAS_3D[o["vista_principal"]][2][i])))
    except Exception as e:
        prep["avisos"].append(f"figuras_3d: {type(e).__name__}: {e}")
    _figura_metodo(prep)
    return prep


def pie_metodo(prep, idioma="es"):
    """Pie de la figura 0, con lo que de verdad lleva: con o sin pila, con o
    sin VOI, una familia o dos."""
    from . import figura_metodo as FM
    i = _idx(idioma)
    fm = _Formateador(idioma)
    ctx = prep.get("contexto_voi")
    fams = [f for f in ("spinodoide", "dual-lattice")
            if f in (prep.get("generados") or {})]
    hay_voi = any(e[0] == "voi" for e in prep["estructuras"])
    partes = [PIES_METODO["cabecera"][i]]
    if hay_voi and ctx is not None:
        partes.append(PIES_METODO["voi_pila"][i].format(
            k=int(ctx["k"]) + 1, n=int(ctx["n_rebanadas"]),
            p=int(ctx["paso"]), lado=_lado_txt(prep, idioma),
            gris=PIES_METODO["gris" if ctx.get("gris") is not None
                             else "sin_gris"][i]))
    elif hay_voi:
        partes.append(PIES_METODO["voi_archivo"][i].format(
            lado=_lado_txt(prep, idioma)))
    if fams:
        partes.append(PIES_METODO["familias"][i])
    if "spinodoide" in fams:
        partes.append(PIES_METODO["spinodoide"][i])
    if "dual-lattice" in fams:
        lo, hi = FM.SUBCUBO_ESQUELETO
        partes.append(PIES_METODO["dual"][i].format(
            sub=fm.format_field(100.0 * (hi - lo), ".0f")))
    if fams:
        partes.append(PIES_METODO["tablas"][i])
    o = _opciones_fig(prep)
    texto = " ".join(partes)
    return texto + (PIES["suavizado"][i] if o["suavizar"] else "")


def _opciones_fig(prep):
    from . import figuras as F
    return F.opciones_3d(prep.get("opciones_3d"))


def _figura_metodo(prep):
    """Figura 0 (hilo de la interfaz: renders VTK). Si falla, se avisa y el
    resto del informe sigue."""
    o = _opciones_fig(prep)
    gen = prep.get("generados")
    if not o["metodo"] or gen is None:
        return prep
    try:
        from . import figura_metodo as FM
        doc = prep["doc"]
        res = doc.get("resultados") or {}
        familias = {}
        for (fam, BW, _sp) in prep["estructuras"]:
            if fam == "voi" or fam not in gen:
                continue
            clave = "ajuste" if fam == "spinodoide" else "ajuste_dual"
            par = dict((res.get(clave) or {}).get("parametros") or {})
            par.setdefault("familia", fam)
            campo, info = gen[fam]
            familias[fam] = FM.datos_familia(par, BW=BW, campo=campo,
                                             info=info)
        voi = next(((BW, sp) for (e, BW, sp) in prep["estructuras"]
                    if e == "voi"), None)
        if not familias and voi is None:
            return prep
        d = prep["carpeta"] / CARPETA_FIGURAS
        rutas, fallos = FM.figura_metodo(voi, familias, doc, d,
                                         prep["lado_mm"],
                                         ctx=prep.get("contexto_voi"),
                                         suavizar=o["suavizar"])
        for f in fallos:
            prep["avisos"].append(f"figuras_3d: {f}")
        for idioma in ("es", "en"):
            if rutas.get(idioma):
                prep["archivos_figuras"] += rutas[idioma]
                prep["figuras"][idioma].append(
                    (f"{CARPETA_FIGURAS}/{Path(rutas[idioma][0]).name}",
                     pie_metodo(prep, idioma)))
    except Exception as e:                                # noqa: BLE001
        prep["avisos"].append(f"figuras_3d: fig0: {type(e).__name__}: {e}")
    finally:
        # Los campos del generador ocupan lo suyo y ya no hacen falta.
        prep["generados"] = {k: None for k in (gen or {})} if gen else gen
    return prep


def componer(prep, pdf=True):
    """Etapa 4 (hilo de trabajo): JSON, Markdown y PDF en la carpeta."""
    c = prep["carpeta"]
    doc, items, paq = prep["doc"], prep["items"], prep["paquete"]
    escritos = []

    def escribe(nombre, texto):
        ruta = c / nombre
        ruta.write_text(texto, encoding="utf-8")
        escritos.append(str(ruta))
        return ruta

    escribe(ARCHIVOS["sesion"],
            json.dumps(procedencia.serializable(doc), indent=1,
                       ensure_ascii=False, default=str))
    escribe(ARCHIVOS["paquete"], json.dumps(paq, indent=2, ensure_ascii=False))
    # Por NUMERO de figura: los cortes (fig5) se dibujan antes que los
    # renders (fig3, fig4), pero en el informe van despues; y fig10 va detras
    # de fig9, cosa que el orden de texto no respeta.
    mds = {idioma: escribe(ARCHIVOS["md_" + idioma], informe_markdown(
               doc, items, paq, idioma, ARCHIVOS["paquete"],
               sorted(prep["figuras"][idioma],
                      key=lambda f: _numero_figura(f[0])),
               galeria=prep.get("galeria")))
           for idioma in ("es", "en")}
    if pdf:
        try:
            from . import md2pdf
        except ImportError:
            prep["avisos"].append("pdf_sin_reportlab")
        else:
            for idioma, etiqueta, cabecera, llano in (
                    ("es", "Figura", "spinpy — informe para publicación",
                     "EN PALABRAS SENCILLAS"),
                    ("en", "Figure", "spinpy — publication report",
                     "IN PLAIN WORDS")):
                ruta = c / ARCHIVOS["pdf_" + idioma]
                md2pdf.componer(str(mds[idioma]), str(ruta), cabecera,
                                etiqueta, llano)
                escritos.append(str(ruta))
    return {"carpeta": str(c), "archivos": escritos + prep["archivos_figuras"],
            "items": items, "resumen": resumen(items), "paquete": paq,
            "avisos": list(prep["avisos"])}


def escribir_informe(doc, carpeta, regenerar=True, VOI=None, spacing=None,
                     figuras=True, render_3d=True, pdf=True, opciones_3d=None,
                     campos_vm=None, contexto_voi=None):
    """Todas las etapas seguidas, en el hilo que llama. Para scripts y la CLI.

    Devuelve {"carpeta", "archivos", "items", "resumen", "paquete", "avisos"}.
    """
    prep = preparar(doc, carpeta, VOI, spacing, regenerar,
                    cargar_voi=render_3d, opciones_3d=opciones_3d,
                    campos_vm=campos_vm, contexto_voi=contexto_voi)
    if figuras:
        figuras_datos(prep)
    if render_3d:
        figuras_3d(prep)
    return componer(prep, pdf)


# ---------------------------------------------------------------------------
# Linea de comandos
# ---------------------------------------------------------------------------

def main(argv=None):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    p = argparse.ArgumentParser(
        prog="spinpy-informe",
        description="Informe para publicacion desde un JSON de 'Exportar "
                    "resultados' del visor, o verificacion de un paquete de "
                    "reproduccion.")
    p.add_argument("entrada", type=Path,
                   help="resultados .json del visor; con --verificar, el "
                        "archivo _reproduccion.json")
    p.add_argument("--carpeta", type=Path, default=None,
                   help="carpeta de salida (por omision, informe_<entrada> "
                        "junto a la entrada)")
    p.add_argument("--sin-pdf", action="store_true")
    p.add_argument("--sin-figuras", action="store_true")
    p.add_argument("--sin-3d", action="store_true")
    p.add_argument("--verificar", action="store_true",
                   help="regenera las estructuras y compara sus huellas")
    p.add_argument("--voi", type=Path, default=None,
                   help="con --verificar: comprueba tambien la huella del VOI")
    p.add_argument("--sin-regenerar", action="store_true",
                   help="no ejecuta el generador (el paquete queda sin huellas)")
    a = p.parse_args(argv)

    if a.verificar:
        paq = json.loads(a.entrada.read_text(encoding="utf-8"))
        r = verificar_reproduccion(paq, a.voi)
        for fam, d in r["ajustes"].items():
            print(f"  {fam:13s} {'OK' if d['ok'] else 'DISTINTO' if d['ok'] is False else 'sin huella'}"
                  + (f"  ({d['motivo']})" if d.get("motivo") else ""))
        if r["voi"] is not None:
            print(f"  {'VOI':13s} {'OK' if r['voi']['ok'] else 'DISTINTO'}")
        for av in r["avisos"]:
            print(f"  aviso: {av}")
        print("REPRODUCIDO" if r["ok"] else "NO REPRODUCIDO")
        return 0 if r["ok"] else 1

    doc = procedencia.leer_json(a.entrada)
    carpeta = a.carpeta or a.entrada.with_name("informe_" + a.entrada.stem)
    r = escribir_informe(doc, carpeta, regenerar=not a.sin_regenerar,
                         figuras=not a.sin_figuras, render_3d=not a.sin_3d,
                         pdf=not a.sin_pdf)
    for av in r["avisos"]:
        print(f"  aviso: {texto_aviso(av)}")
    c = r["resumen"]
    print(f"{c[CITABLE]} citables, {c[RESERVAS]} con reservas, "
          f"{c[NO_CITABLE]} no citables")
    for it in r["items"]:
        if it["estado"] == NO_CITABLE:
            print(f"  NO CITABLE  {etiqueta_item(it)}: {texto_motivos(it)}")
    for f in r["archivos"]:
        print(f"  -> {f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
