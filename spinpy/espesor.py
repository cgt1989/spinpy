"""
espesor.py — Campo de espesor local (Tb.Th punto a punto).

Metodo de Hildebrand y Ruegsegger (1997), que es el que usan CTAn y Scanco:
el espesor local en un punto es el DIAMETRO DE LA MAYOR ESFERA que contiene
ese punto y cabe entera dentro de la estructura.

    Th(y) = 2 * max { r(x) : |y - x| <= r(x) }

siendo r(x) la distancia de x al fondo, es decir el radio de la mayor esfera
CENTRADA en x que cabe.

POR QUE NO BASTA LA TRANSFORMADA DE DISTANCIA
----------------------------------------------
Es tentador colorear directamente por la distancia al fondo y llamarlo espesor.
No lo es: r(x) solo da el espesor en el eje medial, y cae a cero al acercarse a
la superficie. Como la superficie es justo lo que se ve al renderizar, colorear
por r(x) mostraria una estructura uniformemente delgada por fuera, con el valor
bueno escondido en el interior. El maximo sobre las esferas que CONTIENEN el
punto —no las centradas en el— es lo que reparte el espesor del eje medial a
todo el grosor de la trabecula.

SESGO DE DISCRETIZACION — MEDIDO, Y GRANDE A ESTA ESCALA
---------------------------------------------------------
El metodo mide distancias entre CENTROS de voxel, asi que la frontera efectiva
del objeto digitalizado cae medio voxel por dentro y el espesor sale corto. El
sesgo depende de cuantos voxeles cruzan la estructura (cilindro, medido):

    diametro (voxeles)     2      3      4      6     10     20
    error                0.0%  -33%   -29%  -5.7%  -11%   -5.7%

Una losa plana sale EXACTA (0.00%), porque no hay curvatura que discretizar.

Esto importa mucho aqui: los VOIs equinos tienen Tb.Th ~0.17 mm con voxel de
0.0515 mm, es decir unos 3.3 voxeles de diametro, justo en la zona donde el
sesgo ronda el 30%. LOS VALORES ABSOLUTOS DE ESTE CAMPO NO SON CITABLES COMO
Tb.Th.

Para lo que si sirve:
  * ver DONDE la estructura es gruesa o fina, que es informacion relativa y el
    sesgo es sistematico;
  * comparar VOI y spinodoide medidos con el MISMO voxel, porque entonces el
    sesgo afecta por igual a los dos y se cancela en gran medida al restar
    —el mismo argumento de la correccion C4.

Para un numero citable esta el Tb.Th de Parfitt que devuelve `morfometria`.

DIFERENCIA CON EL Tb.Th DE LA MORFOMETRIA
------------------------------------------
`morfometria` devuelve el Tb.Th de Parfitt, 2*BV/BS, que es un MODELO DE PLACAS
paralelas: un unico numero por VOI, obtenido de un cociente global. Este modulo
da un campo, y su media NO tiene por que coincidir con aquel valor — son dos
definiciones distintas, y las dos son estandar. La de Parfitt es la que reportan
CTAn y Scanco por defecto y la que usa el ajuste; esta sirve para VER donde la
estructura es gruesa o fina.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage


def espesor_local(BW, spacing=1.0, n_radios=None):
    """Campo de espesor local, en las mismas unidades que `spacing`.

    Devuelve un array del tamano de BW con el espesor en cada voxel solido y
    cero en el fondo.

    El algoritmo recorre los radios de mayor a menor. Para cada radio r toma
    los centros que admiten una esfera de al menos ese radio y marca todo lo
    que queda a distancia <= r de alguno; lo que se marca por primera vez tiene
    espesor 2r. Al ir de mayor a menor, cada voxel se queda con el mayor radio
    que lo cubre, que es la definicion.

    `n_radios` limita cuantos niveles se prueban (uno por nivel cuesta una
    transformada de distancia). Con None se usan todos los radios distintos que
    aparecen, redondeados a la rejilla.
    """
    BW = np.asarray(BW, dtype=bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    if not BW.any():
        return np.zeros(BW.shape, dtype=np.float32)

    # r(x): radio de la mayor esfera centrada en x que cabe en el solido
    r = ndimage.distance_transform_edt(BW, sampling=spacing)

    rmax = float(r.max())
    if rmax <= 0:
        return np.zeros(BW.shape, dtype=np.float32)

    paso = float(np.min(spacing))
    if n_radios is None:
        radios = np.arange(rmax, 0, -paso)
    else:
        radios = np.linspace(rmax, paso, int(n_radios))

    esp = np.zeros(BW.shape, dtype=np.float32)
    pendientes = BW.copy()

    for rad in radios:
        centros = r >= rad
        if not centros.any():
            continue
        # Distancia a los centros admisibles: lo que queda a <= rad de alguno
        # esta cubierto por una esfera de radio rad.
        d = ndimage.distance_transform_edt(~centros, sampling=spacing)
        cubre = (d <= rad) & pendientes
        if cubre.any():
            esp[cubre] = 2.0 * rad
            pendientes &= ~cubre
        if not pendientes.any():
            break

    # Lo que no cubrio ningun radio (voxeles sueltos) se queda con 2*r(x)
    if pendientes.any():
        esp[pendientes] = (2.0 * r[pendientes]).astype(np.float32)
    return esp


def estadisticas(esp, BW=None):
    """Resumen del campo: media, mediana, desviacion y percentiles."""
    v = np.asarray(esp, float)
    v = v[v > 0] if BW is None else v[np.asarray(BW, bool)]
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {}
    return {"media": float(v.mean()), "mediana": float(np.median(v)),
            "sd": float(v.std(ddof=1)) if v.size > 1 else 0.0,
            "p05": float(np.percentile(v, 5)),
            "p95": float(np.percentile(v, 95)),
            "min": float(v.min()), "max": float(v.max())}


def muestrear_en_puntos(esp, spacing, puntos, origen=None):
    """Valor del campo en coordenadas arbitrarias, por vecino mas proximo.

    Los vertices de la isosuperficie caen ENTRE voxeles, y con frecuencia del
    lado del fondo, donde el espesor vale cero. Colorear con eso pintaria toda
    la piel de azul. Se dilata antes el campo con un maximo 3x3x3 para que cada
    vertice recoja el valor del solido contiguo.
    """
    esp = np.asarray(esp, np.float32)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    origen = np.zeros(3) if origen is None else np.asarray(origen, float)

    campo = ndimage.maximum_filter(esp, size=3, mode="nearest")
    idx = np.rint((np.asarray(puntos, float) - origen) / spacing).astype(int)
    for k in range(3):
        np.clip(idx[:, k], 0, esp.shape[k] - 1, out=idx[:, k])
    return campo[idx[:, 0], idx[:, 1], idx[:, 2]]
