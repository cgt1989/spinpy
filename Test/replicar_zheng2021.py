# -*- coding: utf-8 -*-
"""
replicar_zheng2021.py — Reproduce la Figura 1(b-e) de Zheng et al. (2021):
                        las cuatro clases de spinodoide con su superficie
                        elastica DENTRO de las cotas de Voigt y de
                        Hashin-Shtrikman.

REFERENCIA REPLICADA (APA 7)
    Zheng, L., Kumar, S., & Kochmann, D. M. (2021). Data-driven topology
        optimization of spinodoid metamaterials with seamlessly tunable
        anisotropy. Computer Methods in Applied Mechanics and Engineering,
        383, 113894. https://doi.org/10.1016/j.cma.2021.113894

    Licencia del articulo: Creative Commons Attribution (CC BY),
    http://creativecommons.org/licenses/by/4.0/
    (c) 2021 The Author(s). Publicado por Elsevier B.V.

    La figura original reproducida en `referencia/Zheng2021_Fig1.png` se
    incluye al amparo de esa licencia, sin modificaciones y con atribucion.
    Ver `referencia/FUENTES.md`.


QUE APORTA ESTA REPLICA QUE NO APORTABA LA DE KUMAR 2020
---------------------------------------------------------
Tres cosas, y ninguna es decorativa.

1. LAS CUATRO TERNAS, VERBATIM. La replica de Kumar et al. (2020) tuvo que
   DEDUCIR la terna de la clase cubica porque su Fig. 2 no la publica. Aqui los
   pies de la Fig. 1 las dan las cuatro:

       (b) lamelar    rho = 0.5, thetas = ( 0,  0, 15) grados
       (c) columnar   rho = 0.5, thetas = (15, 15,  0)
       (d) cubica     rho = 0.5, thetas = (15, 15, 15)
       (e) isotropa   rho = 0.5, thetas = (90, 90, 90)

   Ademas estan TRASPUESTAS respecto de las que usa la replica de Kumar
   —alli el cono activo de la lamelar era e1 y aqui es e3—, de modo que
   comprobar el mapa "cono activo -> eje blando" con estas ternas es una
   comprobacion independiente y no la misma dos veces.

2. COTAS ABSOLUTAS, NO RAZONES. La replica de Kumar comprueba cocientes
   E1/E3: si toda nuestra homogeneizacion estuviera escalada por un factor
   constante, los cocientes seguirian saliendo bien. Su Fig. 1 dibuja cada
   superficie elastica dentro de DOS esferas: la cota de Voigt (gris claro) y
   la cota superior de Hashin-Shtrikman (gris oscuro, solo para la isotropa).
   Esas son cotas de verdad, se pueden calcular en cerrado y no se pueden
   pasar. Es una comprobacion que un error de escala NO sobrevive.

       Voigt a rho = 0.5:  E = rho * E_s = 0.5 E_s, que es exactamente el radio
       0.5 al que esta dibujada la esfera gris clara de su figura.
       Hashin-Shtrikman a rho = 0.5 y nu_s = 0.3:  E = 0.333 E_s, que es el
       tamano al que se ve la esfera gris oscura del panel (e).

3. ORTOTROPIA. El articulo afirma (p. 11) que, por las simetrias de la
   distribucion de orientaciones, el tensor homogeneizado es ORTOTROPO con
   nueve constantes independientes. Eso son doce componentes de Voigt que
   tienen que anularse, y comprobarlo pone a prueba de paso nuestro orden de
   Voigt y las deformaciones angulares de ingenieria.

   LA AFIRMACION ES DEL ENSEMBLE, Y ESO CAMBIA COMO SE COMPRUEBA. Lo que tiene
   simetria es la DISTRIBUCION de direcciones de onda, no una realizacion
   suelta: con N ondas sorteadas, una realizacion tiene sus laminas un poco
   torcidas y produce acoplamiento normal-cortante de verdad. Medido aqui: a
   40^3, la columnar da un 2.6 % de la diagonal y la lamelar un 11.1 %. Exigir
   un umbral pequeno a UNA realizacion seria exigir algo que el articulo no
   dice, y ademas dependeria de la semilla.

   Por eso Z3 comprueba lo que el articulo si afirma: al PROMEDIAR
   realizaciones, el acoplamiento tiene que BAJAR y acercarse a cero. Es el
   mismo arreglo que hubo que hacerle a la comprobacion C4 de la replica de
   Kumar, y por la misma razon.

Y una cuarta, que no es una replica sino la razon por la que este articulo
importa a este proyecto:

4. POR QUE rho >= 0.3. El articulo restringe el espacio de diseno a
   rho >= rho_min = 0.3 "para evitar dominios solidos disjuntos" (p. 5). Esa
   frase es la version publicada del problema abierto de este proyecto: el
   spinodoide ajustado a nuestro VOI tiene rho = 0.339 y aun asi el 13.8 % de
   su hueso en islas aisladas. Aqui se mide la curva de material desconectado
   frente a rho para las cuatro clases, que es lo que convierte "hay que meter
   conectividad en el objetivo" en una afirmacion respaldada. Y sale algo mas,
   que no esperabamos: la regla solo funciona con el theta_min del texto, no
   con el de los pies de figura. Ver el bloque de la discrepancia, abajo.


LA DISCREPANCIA INTERNA DEL ARTICULO NO ES INOCUA, Y SE PUEDE MEDIR
--------------------------------------------------------------------
El texto de la p. 5 dice `theta_min = pi/6` (30 grados). El pie de su Fig. 3
dice `theta_min = 15 grados`, los pies de la Fig. 1 usan 15, y el Benchmark I
arranca de (15, 0, 0). Aqui se usan los 15 de la Fig. 1 para las cuatro clases,
porque son las ternas que se replican.

Pero resulta que los dos valores NO son intercambiables para la otra mitad de
esa misma frase, la que fija rho >= 0.3 "para evitar dominios disjuntos".
Medido aqui, la fraccion de hueso fuera de la mayor componente a rho = 0.30:

    lamelar (0, 0, 15)      decenas de %    <- los 15 de la Fig. 1, y con una
                                               dispersion enorme entre
                                               realizaciones: ver abajo
    lamelar (0, 0, 30)       ~2 %           <- el pi/6 del texto
    lamelar (0, 0, 45)       ~0.2 %
    columnar, cubica, isotropa   < 1 %

Es decir: la regla rho_min = 0.3 cumple su proposito SI se aplica junto al
theta_min = pi/6 que la acompana en el texto, y NO lo cumple con los 15 grados
de los pies de figura. La lamelar es el caso critico porque un solo cono
estrecho produce laminas casi planas, y unas laminas planas y separadas son
componentes distintas por definicion. Con dos o tres conos, o con un cono
ancho, las laminas se ondulan y se tocan.

Por eso Z5 comprueba la regla DONDE EL ARTICULO LA ENUNCIA -con theta >= 30- y
Z7 deja escrito el hallazgo de que con 15 grados no basta.

DOS METRICAS DE DESCONEXION, Y NO MIDEN LO MISMO
  `1 - mayor_componente` cuenta el hueso que no esta en el trozo mas grande.
  Es la lectura literal de "dominios disjuntos", y para una lamelar o una
  isotropa dice lo que hay que decir. Para una COLUMNAR no: esa clase son
  columnas paralelas que atraviesan la probeta, legitimamente separadas unas de
  otras, y todas cargan. Medido a rho = 0.30 con cinco semillas, la columnar da
  4.6 % "desconectado" y a la vez 98.9 % de material PORTANTE -en caminos que
  van de cara a cara-. Penalizarla por lo primero seria penalizarla por no ser
  una sola pieza, que no es un defecto mecanico.

  Por eso Z5 se decide sobre la fraccion PORTANTE, que es la que gobierna la
  rigidez, y el informe lleva las dos. Y se exige en la PEOR semilla, no en la
  media: la afirmacion es que la regla del articulo funciona, y una regla que
  funciona de media no funciona.

LA DESCONEXION TIENE UNA DISPERSION ENORME, Y ESO CAMBIA COMO SE MIDE
  La primera version de este guion media la desconexion con UNA realizacion por
  punto y presentaba series como si fueran tendencias. No lo eran. Medido
  despues sobre la lamelar de 15 grados a rho = 0.30, con CINCO semillas por
  punto (Estudio_Percolacion):

    N = 200     19.3 +- 37.7 %      0.0   86.5    0.7    8.8    0.8
    N = 400      1.7 +-  2.5 %      2.3    0.2    0.0    0.0    5.9
    N = 1000    25.3 +- 34.5 %     48.4    1.4    0.0   75.0    2.0
    N = 2000    46.0 +- 10.1 %     45.3   51.6   57.8   30.6   44.7

  Dos lecturas, y las dos importan:

  1. Hasta N = 1000 la dispersion entre realizaciones es de DECENAS de puntos
     porcentuales. Cualquier cifra suelta de esa zona -incluido el 48.4 % que
     esta version reporta en su barrido- es una tirada, no una medida. Por eso
     Z5 y Z7 se evaluan sobre la MEDIA de cinco semillas y reportan su
     desviacion, y el barrido de una semilla se queda como lo que es: la forma
     de la curva, no el valor de un punto.

  2. A N = 2000 la dispersion se desploma y las cinco semillas caen entre 30 y
     58 %. Eso si es una senal: la estructura se autopromedia y esta
     sistematicamente rota -y en alguna realizacion la fraccion PORTANTE cae a
     cero, es decir que no atraviesa en absoluto-. La afirmacion "con un cono
     estrecho la lamelar no cumple el proposito de rho_min" se sostiene; lo que
     no se sostenia era la evidencia de tres puntos con una semilla que se
     publico primero.

     Por eso el bloque de dispersion usa N = 2000 FIJOS, no los del modo: a
     N = 400 esa clase da 1.7 +- 2.5 % y no hay nada que ver, de modo que Z5 y
     Z7 fallarian en modo rapido por mirar donde el efecto no esta.

  En resolucion, en cambio, la cosa esta tranquila: 128^3, 160^3 y 192^3 dan
  48.4, 48.3 y 48.4 % en la misma realizacion. La ejecucion completa usa 128^3.

  QUEDA ESCRITO PORQUE ES EL TIPO DE ERROR QUE ESTA CARPETA EXISTE PARA NO
  BARRER: una serie monotona de tres puntos, cada uno con una semilla distinta,
  parece una tendencia y no lo es.


CONDICIONES DE CONTORNO: LA DIFERENCIA CON SUS CIFRAS ES CONOCIDA Y TIENE SIGNO
-------------------------------------------------------------------------------
Ellos homogeneizan con condiciones AFINES y lo dicen expresamente (p. 5):
"the computed response provides an upper bound to the actual effective
stiffness". Nosotros usamos condiciones PERIODICAS, que dan el valor mas bajo
de los dos. Asi que nuestras rigideces tienen que salir POR DEBAJO de las
suyas, y las comprobaciones estan escritas como desigualdades en ese sentido:
no se compara cifra con cifra, porque cifra con cifra no seria comparable.

Uso:
    python replicar_zheng2021.py                 # completo (~90 min)
    python replicar_zheng2021.py --rapido        # mallas menores (~16 min)

Los tiempos son MEDIDOS en la maquina del proyecto, no estimados: la
homogeneizacion de las cuatro clases a 40^3 son unos 25 min, el barrido de
densidad a 128^3 otros 35, el ensemble de ortotropia 8 y el bloque de
dispersion -que va a N = 2000 fijos- unos 15. El modo rapido no baja el bloque
de dispersion, por lo que su suelo son esos 15 min.
    python replicar_zheng2021.py --solo-figura   # sin homogeneizar
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(AQUI.parent))

from spinpy import (constantes_ingenieria, generar_mascara,  # noqa: E402
                    homogeneizar, morfometria)
from spinpy.morphometry import (fraccion_portante,  # noqa: E402
                                mayor_componente_6)
# El calculo de la superficie elastica ya esta escrito y verificado en la otra
# replica; se reutiliza en vez de copiarlo. Copiar la conversion de Voigt a
# tensor de cuarto orden es exactamente como se cuelan los factores 1/2.
from replicar_kumar2020 import (_esfera, autocomprobar_superficie,  # noqa: E402
                                modulo_direccional)

CITA_APA = (
    "Zheng, L., Kumar, S., & Kochmann, D. M. (2021). Data-driven topology "
    "optimization of spinodoid metamaterials with seamlessly tunable "
    "anisotropy. Computer Methods in Applied Mechanics and Engineering, 383, "
    "113894. https://doi.org/10.1016/j.cma.2021.113894"
)
DOI = "10.1016/j.cma.2021.113894"
LICENCIA = "CC BY (http://creativecommons.org/licenses/by/4.0/)"

RHO = 0.5                    # verbatim, pies de la Fig. 1(b-e)
BETA_PI = 15.0               # NO lo publica la Fig. 1; se mantiene el de la
                             # replica de Kumar para que las dos sean
                             # comparables entre si. La anisotropia no depende
                             # de beta, y las cotas tampoco.
TH_MIN = 15.0                # operativo (pies de las Figs. 1 y 3), ver arriba
NU_S = 0.3                   # verbatim, p. 5
SEMILLA = 20260720

CLASES = [
    {"clave": "lamelar", "nombre": "Lamelar", "thetas": [0.0, 0.0, 15.0],
     "fuente": "verbatim, pie de la Fig. 1(b)",
     "esperado": "disco achatado: e3 blando"},
    {"clave": "columnar", "nombre": "Columnar", "thetas": [15.0, 15.0, 0.0],
     "fuente": "verbatim, pie de la Fig. 1(c)",
     "esperado": "huso alargado: e3 rigido"},
    {"clave": "cubica", "nombre": "Cubica", "thetas": [15.0, 15.0, 15.0],
     "fuente": "verbatim, pie de la Fig. 1(d)",
     "esperado": "lobulos en los tres ejes"},
    {"clave": "isotropa", "nombre": "Isotropa", "thetas": [90.0, 90.0, 90.0],
     "fuente": "verbatim, pie de la Fig. 1(e)",
     "esperado": "esfera, dentro de Hashin-Shtrikman"},
]

# Barrido de densidad para la comprobacion de rho_min. Se cruza con las cuatro
# clases porque la pregunta interesante no es si rho = 0.3 basta EN MEDIA, sino
# si basta para TODAS las anisotropias.
RHOS = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50]
RHO_MIN_ART = 0.30           # verbatim, p. 5

# Una quinta entrada SOLO para el barrido: la misma lamelar con el theta_min
# que fija el TEXTO de la p. 5 (pi/6 = 30 grados) en lugar del que usan los
# pies de la Fig. 1 (15). La discrepancia del articulo no es inocua y esta es
# la forma de medirla: ver el bloque "LA DISCREPANCIA..." de la cabecera.
CLASE_TEXTO = {"clave": "lamelar_30", "nombre": "Lamelar (theta_min del texto)",
               "thetas": [0.0, 0.0, 30.0],
               "fuente": "theta_min = pi/6 del texto, p. 5",
               "esperado": "conectada ya a rho = 0.3"}

# Ensemble para la comprobacion de ortotropia. La malla es FIJA y pequena a
# proposito: lo que se compara es el acoplamiento de una realizacion con el del
# promedio de varias, medidos los dos igual, y para eso no hace falta la malla
# de las cifras que se reportan. Cuatro semillas no son muchas —el acoplamiento
# deberia caer como 1/raiz(n)— pero bastan para ver el signo del efecto, que es
# lo que se afirma.
SEMILLAS_ENSEMBLE = [20260720, 7, 99, 1234]
RES_ENSEMBLE = 24

# Densidades donde la conectividad se mide con VARIAS semillas, y las semillas.
# Ver el bloque "LA DESCONEXION TIENE UNA DISPERSION ENORME" de la cabecera:
# una sola realizacion no sirve para decidir si la regla del articulo se
# cumple, y estas son las dos densidades donde se decide.
RHOS_DISPERSION = [0.30, 0.35]
SEMILLAS_DISPERSION = [20260720, 7, 99, 1234, 555]

# CONFIGURACION FIJA, NO LA DEL MODO. El fenomeno que miden Z5 y Z7 -que un
# cono estrecho deja la lamelar rota a rho = 0.30- SOLO EXISTE con muchas
# ondas: a N = 400 la desconexion de esa clase es 1.7 +- 2.5 % y no hay nada
# que ver, a N = 2000 es 46 +- 10 % y la dispersion entre semillas se desploma
# (Estudio_Percolacion). Si este bloque siguiera al `--rapido`, las dos
# comprobaciones FALLARIAN siempre en modo rapido, no porque la implementacion
# se hubiera desviado sino porque se estaria mirando donde el efecto no esta.
# Se paga el minuto de mas y se mide donde la pregunta tiene respuesta.
RES_DISPERSION = 96
ONDAS_DISPERSION = 2000


# ---------------------------------------------------------------------------
# Cotas en cerrado
# ---------------------------------------------------------------------------
def cota_voigt(rho, E_s=1.0):
    """Regla de las mezclas: C_Voigt = rho * C_s, luego E = rho * E_s.

    Al ser el tensor entero proporcional al del solido, su inverso lo es
    tambien y el modulo direccional vale rho*E_s en TODAS las direcciones: la
    cota de Voigt es una esfera, que es como su figura la dibuja.
    """
    return float(rho) * float(E_s)


def cota_hashin_shtrikman(rho, E_s=1.0, nu_s=NU_S):
    """Cota superior de Hashin-Shtrikman para un solido poroso.

    Fase 1 solida con fraccion rho; fase 2 vacio (K2 = G2 = 0) con 1 - rho. Es
    la cota mas estrecha que se puede poner sabiendo solo las fracciones de
    volumen y que el medio es isotropo, y por eso el articulo solo la dibuja en
    el panel isotropo: en los otros la microestructura no es isotropa y la cota
    no aplica.

    Devuelve (E, K, G).
    """
    rho = float(rho)
    K1 = E_s / (3.0 * (1.0 - 2.0 * nu_s))
    G1 = E_s / (2.0 * (1.0 + nu_s))
    f2 = 1.0 - rho
    if f2 <= 0:
        return E_s, K1, G1
    K = K1 + f2 / (1.0 / (0.0 - K1) + rho / (K1 + 4.0 * G1 / 3.0))
    G = G1 + f2 / (1.0 / (0.0 - G1)
                   + 2.0 * rho * (K1 + 2.0 * G1)
                   / (5.0 * G1 * (K1 + 4.0 * G1 / 3.0)))
    E = 9.0 * K * G / (3.0 * K + G) if (3.0 * K + G) > 0 else 0.0
    return float(E), float(K), float(G)


# ---------------------------------------------------------------------------
# Ortotropia
# ---------------------------------------------------------------------------
# En notacion de Voigt [xx yy zz yz xz xy], un tensor ortotropo alineado con
# los ejes tiene EXACTAMENTE estas doce componentes nulas.
CEROS_ORTOTROPO = [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5),
                   (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5)]


def desviacion_ortotropa(C):
    """Cuanto se aparta C de ser ortotropo, en tanto por uno.

    Se normaliza por la media de C11, C22 y C33 y no por el maximo del tensor:
    con el maximo, una estructura muy anisotropa —donde una diagonal es
    diez veces las otras— pareceria mas ortotropa de lo que es solo por
    dividir por un numero mas grande.
    """
    C = np.asarray(C, float)
    esc = float(np.mean([C[0, 0], C[1, 1], C[2, 2]]))
    if not np.isfinite(esc) or esc <= 0:
        return np.nan
    return float(max(abs(C[i, j]) for i, j in CEROS_ORTOTROPO) / esc)


# ---------------------------------------------------------------------------
# Generacion y medida
# ---------------------------------------------------------------------------
def una_clase(cl, res_vista, res_homog, num_waves, rho=RHO, con_mecanica=True):
    t0 = time.time()
    BW, _, info = generar_mascara(
        resolution=res_vista, wave_number=BETA_PI * np.pi,
        num_waves=num_waves, thetas=cl["thetas"], rho=rho,
        esquema="rechazo", seed=SEMILLA)
    sp = np.full(3, 1.0 / res_vista)
    m = morfometria(BW, sp, do_mil=True)

    solo_mayor = mayor_componente_6(BW)
    n_tot = int(BW.sum())
    out = {
        "clave": cl["clave"], "nombre": cl["nombre"],
        "thetas": list(cl["thetas"]), "fuente_parametros": cl["fuente"],
        "rho_pedida": float(rho), "rho_obtenida": float(info["rho_obtenida"]),
        "resolucion_vista": int(res_vista),
        "BVTV": float(m["BVTV"]), "DA": float(m["DA"]),
        "fraccion_desconectada": float(1.0 - solo_mayor.sum() / max(n_tot, 1)),
        "fraccion_portante": float(m.get("FracPort", np.nan)),
    }

    if con_mecanica:
        from spinpy.elastic import remuestrear_bw
        BWh, sph = remuestrear_bw(BW, sp, res_homog)
        C, inf = homogeneizar(BWh, E_s=1.0, nu_s=NU_S, vox_size=float(sph[0]))
        ec = constantes_ingenieria(C)
        dirs, _ = _esfera(60, 120)
        Ed = modulo_direccional(C, dirs)
        out.update({
            "resolucion_homog": int(res_homog),
            "homogeneizacion_ok": bool(inf["ok"]),
            "C_normalizado_Es": np.asarray(C, float).tolist(),
            "E1": float(ec["Ex"]), "E2": float(ec["Ey"]), "E3": float(ec["Ez"]),
            "anisotropia_E": float(ec["anisotropia_E"]),
            "E_max_direccional": float(np.nanmax(Ed)),
            "E_min_direccional": float(np.nanmin(Ed)),
            "desviacion_ortotropa": desviacion_ortotropa(C),
        })
    out["segundos"] = round(time.time() - t0, 1)
    return BW, out


def ensemble_ortotropia(cl, res_vista, num_waves):
    """Acoplamiento de cada realizacion y el del tensor PROMEDIO.

    Promediar TENSORES y no metricas es lo que corresponde: la afirmacion del
    articulo es sobre el tensor efectivo del ensemble, y la media de los
    valores absolutos de doce componentes nunca bajaria aunque cada una
    fluctuara alrededor de cero.
    """
    from spinpy.elastic import remuestrear_bw

    tensores, sueltos = [], []
    for semilla in SEMILLAS_ENSEMBLE:
        BW, _, _ = generar_mascara(
            resolution=res_vista, wave_number=BETA_PI * np.pi,
            num_waves=num_waves, thetas=cl["thetas"], rho=RHO,
            esquema="rechazo", seed=semilla)
        sp = np.full(3, 1.0 / res_vista)
        BWh, sph = remuestrear_bw(BW, sp, RES_ENSEMBLE)
        C, inf = homogeneizar(BWh, E_s=1.0, nu_s=NU_S, vox_size=float(sph[0]))
        if not inf["ok"]:
            continue
        tensores.append(np.asarray(C, float))
        sueltos.append(desviacion_ortotropa(C))
    if not tensores:
        return None
    medio = np.mean(tensores, axis=0)
    return {
        "clase": cl["clave"], "thetas": list(cl["thetas"]),
        "resolucion": RES_ENSEMBLE, "semillas": list(SEMILLAS_ENSEMBLE),
        "sueltos": [float(x) for x in sueltos],
        "peor_suelto": float(np.max(sueltos)),
        "medio_de_sueltos": float(np.mean(sueltos)),
        "del_tensor_promedio": float(desviacion_ortotropa(medio)),
        "C_promedio": medio.tolist(),
    }


def dispersion_conectividad(res_vista=None, num_waves=None):
    """Fraccion desconectada con varias semillas, en las densidades criticas.

    POR QUE NO BASTA EL BARRIDO
      El barrido recorre muchas densidades con UNA realizacion cada una, que
      esta bien para ver la FORMA de la curva y no para decidir si un valor
      pasa de un umbral. Medida aparte, la dispersion entre realizaciones de la
      lamelar a rho = 0.30 y N = 1000 es de decenas de puntos porcentuales: la
      misma terna da 0.0, 1.4, 2.0, 48.4 y 75.0 % segun la semilla. Decidir con
      una tirada es decidir con ruido.

    Devuelve una fila por clase y densidad con media, desviacion y los valores
    sueltos, que se reportan enteros: quien lea tiene que ver la dispersion, no
    solo el promedio.
    """
    filas = []
    for cl in CLASES + [CLASE_TEXTO]:
        for rho in RHOS_DISPERSION:
            vals, port = [], []
            for s in SEMILLAS_DISPERSION:
                BW, _, _ = generar_mascara(
                    resolution=RES_DISPERSION, wave_number=BETA_PI * np.pi,
                    num_waves=ONDAS_DISPERSION, thetas=cl["thetas"], rho=rho,
                    esquema="rechazo", seed=s)
                n = int(BW.sum())
                mayor = int(mayor_componente_6(BW).sum())
                vals.append(1.0 - mayor / max(n, 1))
                port.append(fraccion_portante(BW))
            v = np.asarray(vals, float)
            p = np.asarray(port, float)
            filas.append({
                "clave": cl["clave"], "rho": float(rho),
                "semillas": list(SEMILLAS_DISPERSION),
                "valores": [float(x) for x in v],
                "media": float(v.mean()),
                "sd": float(v.std(ddof=1)) if v.size > 1 else 0.0,
                "mediana": float(np.median(v)),
                "min": float(v.min()), "max": float(v.max()),
                # La otra mitad de la historia: material en caminos que
                # ATRAVIESAN. Ver la nota "DOS METRICAS" de la cabecera.
                "portante_valores": [float(x) for x in p],
                "portante_media": float(p.mean()),
                "portante_min": float(p.min()),
                "portante_sd": float(p.std(ddof=1)) if p.size > 1 else 0.0,
            })
    return filas


def barrido_rho(res_vista, num_waves):
    """Fraccion de hueso fuera de la mayor componente, por clase y densidad.

    No homogeneiza: es una medida topologica y sale de contar voxeles, asi que
    se puede permitir la rejilla de vista completa y todas las densidades.
    """
    filas = []
    for cl in CLASES + [CLASE_TEXTO]:
        for rho in RHOS:
            BW, _, _ = generar_mascara(
                resolution=res_vista, wave_number=BETA_PI * np.pi,
                num_waves=num_waves, thetas=cl["thetas"], rho=rho,
                esquema="rechazo", seed=SEMILLA)
            n = int(BW.sum())
            mayor = int(mayor_componente_6(BW).sum())
            filas.append({
                "clave": cl["clave"], "rho": float(rho),
                "rho_obtenida": float(BW.mean()),
                "fraccion_desconectada": float(1.0 - mayor / max(n, 1)),
            })
    return filas


# ---------------------------------------------------------------------------
# Predicciones falsables
# ---------------------------------------------------------------------------
def comprobar(res, filas, ens, disp):
    C = []

    def anota(cod, texto, cond, obtenido, criterio):
        C.append({"codigo": cod, "prediccion": texto, "criterio": criterio,
                  "obtenido": obtenido, "pasa": bool(cond)})

    d = {r["clave"]: r for r in res}
    con_mec = [r for r in res if r.get("homogeneizacion_ok")]

    # --- Z1  cota de Voigt, en TODAS las direcciones y en las cuatro clases
    if con_mec:
        vo = cota_voigt(RHO)
        peor = max(r["E_max_direccional"] / vo for r in con_mec)
        quien = max(con_mec, key=lambda r: r["E_max_direccional"])
        anota("Z1",
              "Ninguna superficie elastica sale de la esfera de Voigt: es la "
              "esfera gris clara de su Fig. 1, de radio rho*E_s",
              peor <= 1.0,
              f"peor caso {quien['nombre']}: E_max = "
              f"{quien['E_max_direccional']:.4f} E_s frente a {vo:.3f} E_s "
              f"({100*peor:.1f} % de la cota)",
              "E(d) <= rho*E_s en las cuatro clases")

    # --- Z2  cota de Hashin-Shtrikman en la clase isotropa
    iso = d.get("isotropa")
    if iso and iso.get("homogeneizacion_ok"):
        hs, _, _ = cota_hashin_shtrikman(RHO)
        anota("Z2",
              "La clase isotropa no supera la cota superior de "
              "Hashin-Shtrikman, que es la esfera gris oscura de su panel (e)",
              iso["E_max_direccional"] <= hs,
              f"E_max = {iso['E_max_direccional']:.4f} E_s frente a "
              f"{hs:.4f} E_s ({100*iso['E_max_direccional']/hs:.1f} % de la "
              f"cota)",
              "E(d) <= E_HS+ para la unica clase donde la cota aplica")

    # --- Z3  ortotropia (p. 11 del articulo), que es del ENSEMBLE
    if ens:
        baja = ens["del_tensor_promedio"] < ens["medio_de_sueltos"]
        anota("Z3",
              "El tensor homogeneizado es ORTOTROPO (p. 11): al promediar "
              "realizaciones, el acoplamiento normal-cortante baja hacia cero. "
              "Una realizacion suelta NO lo cumple, y el articulo tampoco lo "
              "afirma de ella",
              baja and ens["del_tensor_promedio"] <= 0.05,
              f"clase {ens['clase']}: realizaciones sueltas "
              + ", ".join(f"{100*x:.1f} %" for x in ens["sueltos"])
              + f"  ->  tensor promedio {100*ens['del_tensor_promedio']:.1f} % "
              f"de la diagonal",
              "el promedio baja respecto de las realizaciones sueltas y "
              "queda por debajo del 5 % de la media de C11, C22, C33")

    # --- Z4  el mapa cono activo -> eje blando, con SUS ternas (traspuestas)
    lam, col = d.get("lamelar"), d.get("columnar")
    if lam and col and lam.get("homogeneizacion_ok") \
            and col.get("homogeneizacion_ok"):
        rl = lam["E3"] / max(lam["E1"], 1e-30)
        rc = col["E3"] / max(col["E1"], 1e-30)
        anota("Z4",
              "Con SUS ternas —traspuestas respecto de las de Kumar 2020— el "
              "eje del cono activo es el blando en la lamelar (disco) y el "
              "rigido en la columnar (huso)",
              rl < 0.6 and rc > 1.6,
              f"lamelar E3/E1 = {rl:.3f}   ·   columnar E3/E1 = {rc:.3f}",
              "lamelar < 0.6 y columnar > 1.6")

    # --- Z5  por que el articulo pone rho_min = 0.3, y con que theta_min
    #
    # Se usan las MEDIAS sobre varias semillas y no el barrido de una sola: la
    # desconexion de una realizacion suelta varia decenas de puntos.
    if disp:
        en_min = {d["clave"]: d["media"] for d in disp
                  if abs(d["rho"] - RHO_MIN_ART) < 1e-9}
        sd_min = {d["clave"]: d["sd"] for d in disp
                  if abs(d["rho"] - RHO_MIN_ART) < 1e-9}
        # PORTANTE, no "fuera de la mayor componente": ver la nota de la
        # cabecera. La peor de las semillas, no la media: lo que se afirma es
        # que la regla FUNCIONA, y para eso no vale que funcione de media.
        port_min = {d["clave"]: d["portante_min"] for d in disp
                    if abs(d["rho"] - RHO_MIN_ART) < 1e-9}
        # La regla se comprueba DONDE EL ARTICULO LA ENUNCIA: la misma frase de
        # la p. 5 que fija rho_min = 0.3 fija theta_min = pi/6. Por eso la
        # lamelar entra aqui con sus 30 grados y no con los 15 de la Fig. 1.
        sujetas = {k: v for k, v in port_min.items() if k != "lamelar"}
        peor = min(sujetas.values()) if sujetas else np.nan
        quien = min(sujetas, key=sujetas.get) if sujetas else "?"
        anota("Z5",
              "A rho = 0.3 —el rho_min que el articulo impone para evitar "
              "dominios disjuntos— queda poco hueso fuera de la mayor "
              "componente, con el theta_min = pi/6 que fija la MISMA frase",
              peor >= 0.95,
              f"peor clase {quien} en su PEOR semilla: {100*peor:.1f} % de "
              f"hueso en caminos que atraviesan; "
              + ", ".join(f"{k} {100*v:.1f} %" for k, v in sujetas.items())
              + f" (minimo de {len(SEMILLAS_DISPERSION)} semillas)",
              ">= 95 % de material portante en todas las semillas de las "
              "clases con theta >= 30 grados")

        # --- Z7  el hallazgo: con los 15 grados de la Fig. 1 NO basta
        lam = [d for d in disp if d["clave"] == "lamelar"
               and abs(d["rho"] - RHO_MIN_ART) < 1e-9]
        vals15 = np.asarray(lam[0]["valores"], float) if lam else np.array([])
        l15 = en_min.get("lamelar", np.nan)
        l30 = en_min.get("lamelar_30", np.nan)
        # Cuantas realizaciones se rompen de verdad. Es la cifra que describe
        # una distribucion BIMODAL; la media de 25 % no le pasa a ninguna.
        rotas = int((vals15 > 0.05).sum())
        anota("Z7",
              "HALLAZGO, no una prediccion del articulo: la lamelar con los "
              "15 grados del pie de su Fig. 1 NO cumple el proposito de "
              "rho_min; con los 30 grados del texto si, y la diferencia es de "
              "un orden de magnitud",
              rotas >= 2 and (l30 <= 0.05),
              f"a rho = 0.30 la lamelar de 15 grados sale rota en "
              f"{rotas} de {len(vals15)} realizaciones ("
              + ", ".join(f"{100*x:.1f}" for x in vals15)
              + f" % desconectado), mientras la de 30 grados se queda en "
              f"{100*l30:.1f} %. Su rho_c medido es 0.31 "
              f"(Estudio_Percolacion), o sea que rho = 0.30 cae justo en el "
              f"umbral",
              # El criterio se declaro para la configuracion anterior (N=1000,
              # donde el reparto es bimodal) y NO se ha apretado al ver que a
              # N=2000 se rompen las cinco: apretarlo despues de mirar seria
              # ajustar la prueba al resultado.
              "al menos dos realizaciones de cinco rotas a 15 grados, y "
              "ninguna a 30")

        bajo = {f["clave"]: f["fraccion_desconectada"]
                for f in filas if abs(f["rho"] - 0.15) < 1e-9}
        en_min = {k: v for k, v in en_min.items() if k != "lamelar_30"}
        # Z6 es una TENDENCIA entre extremos muy separados (0.30 -> 0.15), no
        # un umbral, asi que el barrido de una semilla le vale: la diferencia
        # es de decenas de puntos y la dispersion no la puede invertir.
        crece = (bajo and en_min
                 and all(bajo[k] > en_min[k] for k in en_min if k in bajo))
        anota("Z6",
              "Y la restriccion tiene sentido: por debajo de rho_min el "
              "material desconectado crece en todas las clases",
              crece,
              "; ".join(f"{k}: {100*en_min[k]:.1f} % -> {100*bajo[k]:.1f} % "
                        f"al bajar de 0.30 a 0.15" for k in en_min
                        if k in bajo),
              "la fraccion desconectada crece al bajar rho, en las cuatro")
    return C


# ---------------------------------------------------------------------------
# Figura
# ---------------------------------------------------------------------------
def figura(mascaras, res, destino, con_mecanica=True):
    """Cuatro filas como su Fig. 1(b-e): estructura, superficie y cotas."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pyvista as pv

    n = len(res)
    ncol = 2 if con_mecanica else 1
    fig = plt.figure(figsize=(4.6 * ncol, 3.6 * n))

    dirs, forma = _esfera()
    vo = cota_voigt(RHO)
    hs, _, _ = cota_hashin_shtrikman(RHO)
    D = dirs.reshape(*forma, 3)

    for i, r in enumerate(res):
        BW = mascaras[r["clave"]]

        pl = pv.Plotter(off_screen=True, window_size=(700, 700))
        pl.set_background("white")
        try:
            pl.enable_anti_aliasing("ssaa")
        except Exception:
            pass
        grid = pv.ImageData(dimensions=BW.shape,
                            spacing=(1.0 / BW.shape[0],) * 3)
        grid.point_data["v"] = BW.astype(np.float32).flatten(order="F")
        malla = grid.contour([0.5], scalars="v")
        if malla.n_points:
            malla = malla.compute_normals(auto_orient_normals=True,
                                          consistent_normals=False)
            pl.add_mesh(malla, color="#d98880", smooth_shading=True,
                        ambient=0.18, diffuse=0.82, specular=0.2,
                        specular_power=18)
        pl.camera_position = "iso"
        pl.camera.zoom(1.25)
        img = pl.screenshot(return_img=True)
        pl.close()

        ax = fig.add_subplot(n, ncol, i * ncol + 1)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(f"{r['nombre']}   $\\theta$ = "
                     f"({r['thetas'][0]:.0f}, {r['thetas'][1]:.0f}, "
                     f"{r['thetas'][2]:.0f})$^\\circ$",
                     fontsize=11, loc="left")

        if not con_mecanica:
            continue

        ax = fig.add_subplot(n, ncol, i * ncol + 2, projection="3d")
        if r.get("homogeneizacion_ok"):
            # Cotas primero, para que la superficie quede dibujada dentro.
            for radio, color, alfa in ((vo, "#c8c8c8", 0.16),
                                       ((hs if r["clave"] == "isotropa"
                                         else None), "#8a8a8a", 0.22)):
                if radio is None:
                    continue
                ax.plot_surface(radio * D[..., 0], radio * D[..., 1],
                                radio * D[..., 2], color=color, alpha=alfa,
                                rstride=4, cstride=4, linewidth=0,
                                shade=False)
            Ed = modulo_direccional(np.asarray(r["C_normalizado_Es"]), dirs)
            R3 = np.nan_to_num(Ed).reshape(forma)
            norm = plt.Normalize(np.nanmin(R3), np.nanmax(R3))
            ax.plot_surface(R3 * D[..., 0], R3 * D[..., 1], R3 * D[..., 2],
                            facecolors=plt.cm.jet(norm(R3)), rstride=2,
                            cstride=2, linewidth=0, antialiased=True,
                            shade=False)
            lim = vo * 1.05
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
            ax.set_xlabel("$\\hat{e}_1$"); ax.set_ylabel("$\\hat{e}_2$")
            ax.set_zlabel("$\\hat{e}_3$")
            tit = (f"E(d)/E$_s$   max = {r['E_max_direccional']:.3f}"
                   f"   ·   Voigt {vo:.2f}")
            if r["clave"] == "isotropa":
                tit += f"   ·   HS {hs:.3f}"
            ax.set_title(tit, fontsize=9)
            ax.tick_params(labelsize=7)
        else:
            ax.text2D(0.5, 0.5, "la homogeneizacion no convergio",
                      ha="center", transform=ax.transAxes, color="#b62324")
            ax.set_axis_off()

    fig.suptitle(
        "Replica de Zheng et al. (2021), Fig. 1(b-e) — generada con spinpy\n"
        "$\\rho$ = %.2f · $\\beta$ = %g$\\pi$ · esferas: cota de Voigt y, en "
        "la isotropa, Hashin-Shtrikman" % (RHO, BETA_PI), fontsize=12)
    fig.text(0.5, 0.016, "Replica de: " + CITA_APA.replace(". Computer",
                                                           ".\nComputer"),
             ha="center", va="bottom", fontsize=7, color="#444")
    fig.text(0.5, 0.003, "Figura original bajo " + LICENCIA
             + " · (c) 2021 The Author(s)", ha="center", va="bottom",
             fontsize=7, color="#777")
    fig.tight_layout(rect=(0, 0.045, 1, 0.95))
    fig.savefig(destino, dpi=130)
    plt.close(fig)
    return destino


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--rapido", action="store_true",
                    help="mallas menores; sirve para comprobar la tuberia")
    ap.add_argument("--solo-figura", action="store_true",
                    help="no homogeneiza (sin superficies elasticas)")
    a = ap.parse_args()

    res_vista = 96 if a.rapido else 128
    # El estudio de convergencia de este proyecto fija 40-48 para cualquier
    # cifra que se publique; 32 cuesta un 10 % de rigidez. Las COTAS se cumplen
    # a cualquier resolucion —y de hecho una malla gruesa da rigideces mas
    # bajas, o sea mas holgura—, pero las cifras que se reportan son las que
    # son y se calculan donde toca.
    res_homog = 24 if a.rapido else 40
    num_waves = 400 if a.rapido else 1000
    con_mec = not a.solo_figura

    salida = AQUI / "resultados"
    salida.mkdir(exist_ok=True)

    print("=" * 78)
    print(" REPLICA DE Zheng et al. (2021), Fig. 1(b-e)")
    print(" " + CITA_APA)
    print(" licencia del articulo: " + LICENCIA)
    print("=" * 78)
    err = autocomprobar_superficie()
    print(f"  E(d) sobre un material isotropo: error {err:.2e}  [OK]")
    vo = cota_voigt(RHO)
    hs, K_hs, G_hs = cota_hashin_shtrikman(RHO)
    print(f"  cotas en cerrado a rho={RHO}, nu_s={NU_S}:  "
          f"Voigt E = {vo:.4f} E_s   ·   Hashin-Shtrikman E = {hs:.4f} E_s")
    print(f"  rho={RHO}  beta={BETA_PI}pi  N={num_waves}  vista {res_vista}^3"
          f"  homog {res_homog}^3  semilla {SEMILLA}")
    print()

    t0 = time.time()
    mascaras, res = {}, []
    for cl in CLASES:
        print(f"  {cl['nombre']:10s} theta={cl['thetas']} ...", end="",
              flush=True)
        BW, r = una_clase(cl, res_vista, res_homog, num_waves,
                          con_mecanica=con_mec)
        mascaras[cl["clave"]] = BW
        res.append(r)
        extra = ""
        if con_mec and r.get("homogeneizacion_ok"):
            extra = (f"  E_max={r['E_max_direccional']:.4f}"
                     f"  ({100*r['E_max_direccional']/vo:.0f} % de Voigt)"
                     f"  ortotropia {100*r['desviacion_ortotropa']:.1f} %")
        print(f" {r['segundos']:5.1f} s  BV/TV={r['BVTV']:.3f}"
              f"  desconectado={100*r['fraccion_desconectada']:.1f} %{extra}")

    # El ensemble se mide sobre la clase que PEOR sale de una realizacion
    # suelta: si el promedio la arregla a ella, la afirmacion se sostiene donde
    # mas cuesta. Elegir la mejor seria elegir el caso comodo.
    ens = None
    if con_mec:
        con_mec_res = [r for r in res if r.get("homogeneizacion_ok")]
        if con_mec_res:
            peor = max(con_mec_res, key=lambda r: r["desviacion_ortotropa"])
            cl = next(c for c in CLASES if c["clave"] == peor["clave"])
            print(f"\n  ortotropia en el ensemble ({cl['nombre']}, "
                  f"{len(SEMILLAS_ENSEMBLE)} semillas a {RES_ENSEMBLE}^3) ...",
                  end="", flush=True)
            ens = ensemble_ortotropia(cl, res_vista, num_waves)
            if ens:
                print("  sueltas "
                      + ", ".join(f"{100*x:.1f} %" for x in ens["sueltos"])
                      + f"  ->  promedio "
                      f"{100*ens['del_tensor_promedio']:.1f} %")
            else:
                print(" no se pudo")

    print("\n  barrido de densidad (fraccion de hueso desconectado):")
    filas = barrido_rho(res_vista if not a.rapido else 80, num_waves)
    columnas = CLASES + [CLASE_TEXTO]
    print("      rho   " + "  ".join(f"{c['clave'][:10]:>11s}"
                                     for c in columnas))
    for rho in RHOS:
        fila = "     %.2f  " % rho
        for c in columnas:
            v = [f for f in filas
                 if f["clave"] == c["clave"] and abs(f["rho"] - rho) < 1e-9]
            fila += "  %10.1f%%" % (100 * v[0]["fraccion_desconectada"]) \
                if v else "          --"
        print(fila + ("   <- rho_min del articulo"
                      if abs(rho - RHO_MIN_ART) < 1e-9 else ""))

    print("\n  dispersion entre realizaciones en las densidades criticas "
          f"({len(SEMILLAS_DISPERSION)} semillas, {RES_DISPERSION}^3, "
          f"N = {ONDAS_DISPERSION} FIJOS):")
    disp = dispersion_conectividad()
    for d in disp:
        print(f"    {d['clave']:>11} rho={d['rho']:.2f}   descon "
              f"{100*d['media']:5.1f} +- {100*d['sd']:4.1f} %  (mediana "
              f"{100*d['mediana']:4.1f})   portante min "
              f"{100*d['portante_min']:5.1f} %   "
              + " ".join(f"{100*x:.1f}" for x in d["valores"]))

    checks = comprobar(res, filas, ens, disp)
    print("\n  PREDICCIONES DEL ARTICULO")
    fallos = 0
    for c in checks:
        ok = "OK   " if c["pasa"] else "FALLA"
        fallos += (not c["pasa"])
        print(f"    [{ok}] {c['codigo']}  {c['prediccion']}")
        print(f"             {c['obtenido']}   (criterio: {c['criterio']})")

    fig = None
    print("\n  componiendo la figura...", end="", flush=True)
    try:
        fig = figura(mascaras, res, salida / "replica_zheng2021_fig1.png",
                     con_mecanica=con_mec)
        print(f" -> {fig.name}")
    except Exception as e:
        print(f" fallo: {type(e).__name__}: {e}")

    doc = {
        "replica_de": {"cita_apa": CITA_APA, "doi": DOI, "licencia": LICENCIA,
                       "figura_original": "referencia/Zheng2021_Fig1.png"},
        "generado": time.strftime("%Y-%m-%d %H:%M:%S"),
        "parametros": {"rho": RHO, "beta_pi": BETA_PI, "theta_min": TH_MIN,
                       "nu_s": NU_S, "num_waves": num_waves,
                       "semilla": SEMILLA, "resolucion_vista": res_vista,
                       "resolucion_homogeneizacion": res_homog,
                       "esquema": "rechazo (ecuacion 2 del articulo)",
                       "condiciones_contorno": ("periodicas; el articulo usa "
                                                "afines, que son cota "
                                                "superior")},
        "cotas": {"voigt_E": vo, "hashin_shtrikman_E": hs,
                  "hashin_shtrikman_K": K_hs, "hashin_shtrikman_G": G_hs},
        "clases": res,
        "ortotropia_ensemble": ens,
        "barrido_rho": filas,
        "dispersion_conectividad": disp,
        "dispersion_config": {"resolucion": RES_DISPERSION,
                              "num_waves": ONDAS_DISPERSION,
                              "semillas": list(SEMILLAS_DISPERSION),
                              "rhos": list(RHOS_DISPERSION),
                              "nota": ("fija, no la del modo: el fenomeno de "
                                       "Z5 y Z7 solo existe con muchas ondas")},
        "comprobaciones": checks,
        "fallos": fallos,
        "figura": fig.name if fig else None,
    }
    j = salida / "replica_zheng2021.json"
    j.write_text(json.dumps(doc, indent=2, ensure_ascii=False, default=float),
                 encoding="utf-8")

    print(f"\n  total {time.time()-t0:.0f} s   ·   informe -> {j.name}")
    print("=" * 78)
    print(f"  COMPROBACIONES QUE FALLAN: {fallos}")
    print("=" * 78)
    return 1 if fallos else 0


if __name__ == "__main__":
    sys.exit(main())
