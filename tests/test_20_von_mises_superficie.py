"""
Bloque 20 — Capa superficial de von Mises: el unico sabor del pico que se cita.

REFERENCIA
  `Estudio_Convergencia/INFORME.md` (141 puntos FE sobre cuatro familias de
  experimento, cruzados contra `localVonMisesField` a 1e-12 en C_h y E_app).
  Lo que ese estudio establece:
    · el MAXIMO de von Mises no converge con la resolucion en ninguna
      geometria medida, ni siquiera en una con solucion cerrada: en la cavidad
      esferica periodica el pico binario se pasa un +5.6 % del valor de Goodier
      a la resolucion mas fina y va dando saltos de ~10 puntos por el camino;
    · el percentil 99 TOMADO SOBRE LA CAPA SUPERFICIAL se queda a -0.3 % del
      exacto.
  De ahi `resistencia.capa_superficie` y `resistencia.estadisticos_vm`.

LO QUE SE VERIFICA Y LO QUE NO
  NO se reproduce aqui la cavidad esferica de Goodier: eso es homogeneizacion
  periodica y vive en `Estudio_Convergencia`, no en el ensayo de compresion.
  Lo que si se verifica, y de forma exacta, es la DEFINICION de la capa y la
  propiedad que motiva usarla:

    (1) geometria de la capa sobre casos construidos a mano, donde el numero
        de elementos de superficie se cuenta con los dedos;
    (2) la regla de las isocaps: las seis caras del VOI NO son interfaz, la
        misma convencion con la que `morphometry` excluye las isocaps del area
        superficial;
    (3) NO PERIODICIDAD. `Estudio_Convergencia/code/conv_nucleo.py` define su
        capa con `np.roll`, que envuelve, porque alli la celda es periodica.
        Aqui el ensayo tiene base y techo y envolver marcaria como superficie
        libre justo los elementos que la condicion de contorno ya carga de
        mas. Esta prueba falla si alguien copia la version periodica;
    (4) alineacion de `superficie_solido` con `vm_solido`, tambien con
        `modulo_rel` (fallo progresivo), donde el tejido ablandado sale de la
        muestra pero NO cuenta como vacio;
    (5) la propiedad que justifica todo: al refinar la malla sobre la MISMA
        geometria continua, el p99 de la capa se mueve menos que el maximo.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  (1)-(4): exactas, son conteos enteros y longitudes de vector.
  (5)  criterio, no tolerancia numerica: |Dp99_sup| < |Dmax| en variacion
       relativa entre las dos resoluciones. No se fija un numero porque el
       estudio no da uno para esta geometria; lo que se afirma es el ORDEN,
       que es lo que decide cual de los dos se cita.
"""
import numpy as np

from spinpy.resistencia import (capa_superficie, criterio_pistoia,
                                ensayo_compresion, estadisticos_vm)

BLOQUE = "20 Capa superficial de von Mises"
REF = "Estudio_Convergencia/INFORME.md"


def _giroide(n, nivel=-0.6):
    """Giroide voxelizada: superficie libre por todas partes."""
    x = np.linspace(0.0, 1.0, n)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    f = (np.sin(4 * np.pi * X) * np.cos(4 * np.pi * Y)
         + np.sin(4 * np.pi * Y) * np.cos(4 * np.pi * Z)
         + np.sin(4 * np.pi * Z) * np.cos(4 * np.pi * X))
    return f > nivel


def test_geometria_de_la_capa(registro):
    casos = []

    # Macizo: no hay superficie libre INTERIOR. Las seis caras del dominio son
    # isocaps y no cuentan (regla 2).
    B = np.ones((6, 6, 6), bool)
    casos.append(("cubo macizo: sin superficie libre interior", 0,
                  int(capa_superficie(B).sum())))

    # Un voxel aislado: el es toda la capa.
    B = np.zeros((5, 5, 5), bool)
    B[2, 2, 2] = True
    casos.append(("voxel aislado: la capa es el mismo voxel", 1,
                  int(capa_superficie(B).sum())))

    # Losa que cruza el dominio entero en x,y con vacio arriba y abajo: solo
    # sus dos caras z son interfaz. Los cuatro cantos tocan el borde del
    # dominio y NO cuentan -> exactamente 2*nx*ny.
    nx = ny = 6
    B = np.zeros((nx, ny, 7), bool)
    B[:, :, 2:5] = True
    casos.append(("losa: solo las dos caras z, los cantos son isocaps",
                  2 * nx * ny, int(capa_superficie(B).sum())))

    for nombre, esp, obt in casos:
        registro.anotar(BLOQUE, nombre, REF, esp, obt, 0.0,
                        "conteo exacto de elementos", obt == esp)
        assert obt == esp, (nombre, obt, esp)


def test_la_capa_no_es_periodica(registro):
    """Con `np.roll` la capa z=0 se marcaria; con base y techo reales, no."""
    n = 5
    B = np.ones((n, n, n), bool)
    B[:, :, -1] = False          # vacio solo arriba
    sup = capa_superficie(B)

    # La capa de abajo tiene material encima y el borde del dominio debajo.
    obt = int(sup[:, :, 0].sum())
    registro.anotar(BLOQUE, "no periodica: la base no es superficie libre",
                    REF, 0, obt, 0.0,
                    "conteo exacto; con np.roll saldria n*n", obt == 0)
    assert obt == 0, "la capa se esta calculando con envolvente periodica"

    # La ultima capa con material si lo es: tiene vacio justo encima.
    obt = int(sup[:, :, n - 2].sum())
    registro.anotar(BLOQUE, "no periodica: la capa bajo el vacio si lo es",
                    REF, n * n, obt, 0.0, "conteo exacto", obt == n * n)
    assert obt == n * n


def test_alineacion_con_vm_solido(registro):
    """La mascara tiene que indexar `vm_solido`, tambien con `modulo_rel`."""
    n = 24
    BW = _giroide(n)
    x = np.linspace(0.0, 1.0, n)
    X = np.meshgrid(x, x, x, indexing="ij")[0]

    for etiqueta, kw in (("sin modulo_rel", {}),
                         ("con modulo_rel (fallo progresivo)",
                          {"modulo_rel": np.where(X > 0.5, 0.05, 1.0)})):
        r = ensayo_compresion(BW, [0.05] * 3, apoyo="deslizante", **kw)
        assert r["ok"], r["msg"]
        n_vm = int(np.size(r["vm_solido"]))
        n_sup = int(np.size(r["superficie_solido"]))
        registro.anotar(BLOQUE, f"alineacion mascara/muestra — {etiqueta}",
                        REF, n_vm, n_sup, 0.0,
                        "misma longitud exacta", n_sup == n_vm)
        assert n_sup == n_vm, (etiqueta, n_sup, n_vm)

        # Con la mascara alineada los estadisticos salen de `vm_solido`, que
        # es float64; si se dedujeran del campo (float32) el maximo no
        # coincidiria en la ultima cifra.
        est = estadisticos_vm(r)
        vm = np.asarray(r["vm_solido"], float)
        sup = np.asarray(r["superficie_solido"], bool)
        esp = float(vm[sup].max())
        obt = float(est["vm_max_superficie"])
        registro.anotar(BLOQUE,
                        f"maximo de la capa en doble precision — {etiqueta}",
                        REF, esp, obt, 0.0,
                        "igualdad exacta en float64", obt == esp)
        assert obt == esp


def test_orden_de_los_sabores_del_pico(registro):
    """p99 global <= p99 de la capa <= maximo: la concentracion esta fuera."""
    BW = _giroide(24)
    r = ensayo_compresion(BW, [0.05] * 3, apoyo="deslizante")
    assert r["ok"], r["msg"]
    p = criterio_pistoia(r)
    ok = p["vm_p99"] <= p["vm_p99_superficie"] <= p["vm_max"]
    registro.anotar(BLOQUE, "p99 global <= p99 capa <= maximo", REF, None,
                    p["vm_p99_superficie"], None, "desigualdad", ok,
                    nota="p99 %.4g  capa %.4g  max %.4g  n_capa %d" % (
                        p["vm_p99"], p["vm_p99_superficie"], p["vm_max"],
                        p["vm_n_superficie"]))
    assert ok, (p["vm_p99"], p["vm_p99_superficie"], p["vm_max"])


def test_el_p99_de_la_capa_es_mas_estable_que_el_maximo(registro):
    """La propiedad que decide cual de los dos se cita.

    Misma geometria continua a dos resoluciones. El criterio se declara en la
    cabecera: en variacion RELATIVA, el p99 de la capa tiene que moverse menos
    que el maximo. Si esto falla es un hallazgo sobre el estadistico, no un
    error de la prueba.
    """
    vals = {}
    for n in (24, 36):
        BW = _giroide(n)
        r = ensayo_compresion(BW, [1.2 / n] * 3, apoyo="deslizante")
        assert r["ok"], (n, r["msg"])
        p = criterio_pistoia(r)
        vals[n] = (p["vm_max"], p["vm_p99_superficie"], p["vm_n_superficie"])

    d_max = abs(vals[36][0] / vals[24][0] - 1.0)
    d_sup = abs(vals[36][1] / vals[24][1] - 1.0)
    ok = d_sup < d_max
    registro.anotar(
        BLOQUE, "al refinar, el p99 de la capa se mueve menos que el maximo",
        REF, d_max, d_sup, None, "|D p99_capa| < |D max|, relativo", ok,
        nota="24^3 max %.4g capa %.4g (n %d) | 36^3 max %.4g capa %.4g (n %d)"
             " | D_max %.1f %%  D_capa %.1f %%" % (
                 vals[24][0], vals[24][1], vals[24][2],
                 vals[36][0], vals[36][1], vals[36][2],
                 100 * d_max, 100 * d_sup))
    assert ok, (d_max, d_sup)
