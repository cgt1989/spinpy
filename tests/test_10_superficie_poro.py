"""
Bloque 10 — Superficie interna/externa, S/V del poro y tamano de poro.

QUE SE VERIFICA
  Tres metricas que se anadieron para el lenguaje de los materiales porosos:
  la separacion del area en interfaz real y huella del corte, la razon
  superficie/volumen referida al poro, y el tamano de poro medido por esferas
  inscritas. Las tres se comprueban sobre figuras cuya respuesta se conoce en
  cerrado.

REFERENCIAS
  Bouxsein ML et al. Guidelines for assessment of bone microstructure in
    rodents using micro-CT. J Bone Miner Res 2010;25:1468-86. BS es la interfaz
    osea; las caras de corte del VOI no lo son.
  Hildebrand T, Ruegsegger P. A new method for the model-independent
    assessment of thickness in three-dimensional images. J Microsc
    1997;185:67-75. El tamano de poro es el espesor local de la fase poro.
  Lorensen WE, Cline HE. Marching cubes. SIGGRAPH Comput Graph 1987;21:163-9.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  No se inventan aqui: salen de sesgos que este proyecto ya tenia MEDIDOS Y
  ESCRITOS antes de escribir esta prueba, y lo que se comprueba es que las
  metricas nuevas los heredan y no anaden ninguno propio.

    superficie de una figura CURVA      marching cubes sobre binario
                                        sobreestima ~8.5 % (invariante del
                                        proyecto, medido sobre una esfera)
                                        -> se exige [0, +12 %], y el signo
                                           importa tanto como el margen
    superficie de una cara PLANA        sin curvatura que discretizar, pero la
                                        triangulacion bisela las aristas del
                                        cubo -> se exige [-3 %, +3 %]
    tamano de poro de una losa          `espesor.py` documenta error 0.00 % en
                                        losas -> se exige 1 %
    tamano de poro de una esfera        `espesor.py` documenta -5.7 % a 20
                                        voxeles de diametro -> se exige
                                        [-10 %, 0 %]
    interna + externa = total           identidad, 1e-9

  La suma exacta es la comprobacion mas importante de las cinco: si alguna vez
  alguien calcula la externa contando voxeles en vez de con la misma marching
  cubes, esa identidad deja de cumplirse al 8.5 %, que es justo el error que la
  funcion existe para evitar.
"""
import numpy as np

from conftest import bola, cilindro, losa
from spinpy.morphometry import area_superficie_partes, morfometria, tamano_poro

BLOQUE = "10 Superficie interna/externa y tamano de poro"
REF = "Bouxsein et al. 2010; Hildebrand & Ruegsegger 1997"

N = 40
SP = 0.1                      # mm, del orden del voxel de un micro-CT
ESP = [SP, SP, SP]


def test_cubo_macizo_todo_externo(registro):
    """Un cubo lleno no tiene interfaz: toda su area son las seis tapas."""
    BW = np.ones((N, N, N), dtype=bool)
    interna, externa, total = area_superficie_partes(BW, ESP)

    ok = interna == 0.0
    registro.anotar(BLOQUE, "cubo macizo: area interna", REF, 0.0, interna,
                    "= 0 exactamente", "sin interfaz hueso-vacio", ok)
    assert ok, interna

    exacto = 6.0 * (N * SP) ** 2
    err = (externa - exacto) / exacto
    ok2 = abs(err) <= 0.03
    registro.anotar(BLOQUE, "cubo macizo: area externa", REF, exacto, externa,
                    "error relativo <= 3 %",
                    "seis caras planas; la triangulacion bisela las aristas",
                    ok2)
    assert ok2, (externa, exacto)


def test_esfera_interior_todo_interno(registro):
    """Una esfera que no toca el borde no aporta tapas."""
    R = 12.0
    BW = bola(n=N, radio=R)
    interna, externa, total = area_superficie_partes(BW, ESP)

    ok = externa == 0.0
    registro.anotar(BLOQUE, "esfera interior: area externa", REF, 0.0, externa,
                    "= 0 exactamente", "no toca ninguna cara del cubo", ok)
    assert ok, externa

    exacto = 4.0 * np.pi * (R * SP) ** 2
    err = (interna - exacto) / exacto
    ok2 = 0.0 <= err <= 0.12
    registro.anotar(BLOQUE, "esfera interior: area interna", REF, exacto,
                    interna, "entre 0 y +12 %",
                    "hereda el sesgo +8.5 % de marching cubes binario", ok2)
    assert ok2, (interna, exacto, err)


def test_cilindro_separa_lateral_y_tapas(registro):
    """El cilindro tiene las dos partes a la vez y se conocen las dos."""
    rad, n_xy = 10.0, N
    BW = cilindro(n_xy=n_xy, largo=N, diametro=2 * rad)
    interna, externa, total = area_superficie_partes(BW, ESP)

    exacto_i = 2.0 * np.pi * (rad * SP) * (N * SP)
    exacto_e = 2.0 * np.pi * (rad * SP) ** 2
    for nom, obt, exa in (("lateral (interna)", interna, exacto_i),
                          ("tapas (externa)", externa, exacto_e)):
        err = (obt - exa) / exa
        ok = -0.03 <= err <= 0.12
        registro.anotar(BLOQUE, f"cilindro: {nom}", REF, exa, obt,
                        "entre -3 % y +12 %",
                        "superficie curva con el sesgo de marching cubes", ok)
        assert ok, (nom, obt, exa)


def test_identidad_interna_mas_externa(registro):
    """interna + externa = total, por construccion y sin holgura.

    Es la prueba que impide que alguien calcule la externa contando voxeles:
    esa via daria un total distinto del que da marching cubes.
    """
    BW = cilindro(n_xy=N, largo=N, diametro=16.0)
    interna, externa, total = area_superficie_partes(BW, ESP)
    err = abs((interna + externa) - total)
    ok = err <= 1e-9 * max(total, 1.0)
    registro.anotar(BLOQUE, "interna + externa = total", REF, total,
                    interna + externa, "identidad, 1e-9",
                    "una sola definicion de area para las dos partes", ok)
    assert ok, (interna, externa, total)


def test_tamano_poro_losa_exacto(registro):
    """Un hueco plano de espesor conocido: Po.Dm tiene que dar ese espesor.

    `espesor.py` documenta error 0.00 % en losas —no hay curvatura que
    discretizar— asi que aqui no hay excusa para fallar.
    """
    for t in (4, 8, 12):
        BW = np.ones((N, N, N), dtype=bool)
        k0 = (N - t) // 2
        BW[:, :, k0:k0 + t] = False
        r = tamano_poro(BW, ESP)
        exacto = t * SP
        err = (r["PoDm_mediana"] - exacto) / exacto
        ok = abs(err) <= 0.01
        registro.anotar(BLOQUE, f"hueco plano de {t} voxeles: Po.Dm mediana",
                        REF, exacto, r["PoDm_mediana"],
                        "error relativo <= 1 %",
                        "las losas no tienen sesgo de discretizacion", ok)
        assert ok, (t, r["PoDm_mediana"], exacto)

        # El aviso de ventana se dispara -la esfera se sale por los lados- y
        # sin embargo la medida es exacta, porque el poro CONTINUA de verdad
        # fuera. Es la prueba de que el aviso es un aviso y no un sesgo: si
        # alguien lo convierte en correccion, esta prueba se lo dice.
        ok2 = abs(r["PoDm_interior"] - r["PoDm"]) <= 0.01 * exacto
        registro.anotar(BLOQUE,
                        f"hueco plano de {t} voxeles: interior == total",
                        REF, r["PoDm"], r["PoDm_interior"],
                        "difieren menos del 1 %",
                        "el aviso de ventana no implica sesgo", ok2)
        assert ok2, (r["PoDm"], r["PoDm_interior"])


def test_tamano_poro_esferico(registro):
    """Poro esferico cerrado dentro de hueso macizo: Po.Dm = su diametro."""
    R = 12.0
    BW = ~bola(n=N, radio=R)
    r = tamano_poro(BW, ESP)
    exacto = 2.0 * R * SP
    err = (r["PoDm_mediana"] - exacto) / exacto
    ok = -0.10 <= err <= 0.0
    registro.anotar(BLOQUE, "poro esferico D=24 voxeles: Po.Dm mediana", REF,
                    exacto, r["PoDm_mediana"], "entre -10 % y 0 %",
                    "sesgo de esferas inscritas documentado en espesor.py", ok)
    assert ok, (r["PoDm_mediana"], exacto, err)

    ok2 = r["PoDm_frac_ventana"] == 0.0
    registro.anotar(BLOQUE, "poro esferico: aviso de ventana", REF,
                    0.0, r["PoDm_frac_ventana"], "= 0",
                    "la esfera inscrita cabe entera dentro del cubo", ok2)
    assert ok2, r["PoDm_frac_ventana"]


def test_sv_del_poro_contra_el_del_solido(registro):
    """BS/PV y BS/BV son la MISMA superficie sobre volumenes distintos.

    De ahi sale una identidad que no depende de la geometria:

        (BS/PV) / (BS/BV) = BV/PV = BV.TV / (1 - BV.TV)

    Si alguna vez alguien calcula BS/PV con otra superficie -por ejemplo
    incluyendo las tapas- esta razon deja de cumplirse.
    """
    BW = bola(n=N, radio=14.0)
    m = morfometria(BW, ESP, do_mil=False)
    razon = m["BSPV"] / m["BSBV"]
    exacto = m["BVTV"] / (1.0 - m["BVTV"])
    err = (razon - exacto) / exacto
    ok = abs(err) <= 1e-9
    registro.anotar(BLOQUE, "(BS/PV)/(BS/BV) = BV/PV", REF, exacto, razon,
                    "identidad, 1e-9",
                    "misma interfaz, distinto volumen de referencia", ok)
    assert ok, (razon, exacto)


def test_podm_frente_a_tbsp_en_placas(registro):
    """En una pila de placas, el modelo de Parfitt y la medida coinciden.

    Tb.Sp sale de Tb.Th*(1/BV.TV - 1), que SUPONE placas paralelas. Sobre una
    pila de placas de verdad esa suposicion es exacta, asi que Po.Dm medido y
    Tb.Sp calculado tienen que dar lo mismo. Es la unica geometria donde se
    puede exigir: en hueso trabecular la diferencia entre los dos es
    informacion, no error.
    """
    # La pila EMPIEZA Y TERMINA en hueso, a proposito: si el ultimo hueco
    # quedara abierto contra la cara del cubo, su esfera inscrita no estaria
    # acotada por arriba y Po.Dm mediria el borde del array en vez del poro.
    # Medido al escribir esta prueba con la pila abierta: 1.40 en vez de 1.00.
    # No es un fallo de la metrica -es lo que avisa `PoDm_frac_borde`- sino una
    # geometria de prueba mal puesta.
    t_hueso, t_poro = 6, 10
    n_poros = 2
    nz = t_hueso * (n_poros + 1) + t_poro * n_poros
    BW = np.zeros((N, N, nz), dtype=bool)
    for k in range(n_poros + 1):
        z0 = k * (t_hueso + t_poro)
        BW[:, :, z0:z0 + t_hueso] = True
    r = tamano_poro(BW, ESP)
    exacto = t_poro * SP
    err = (r["PoDm_mediana"] - exacto) / exacto
    ok = abs(err) <= 0.05
    registro.anotar(BLOQUE, "pila de placas: Po.Dm = separacion real", REF,
                    exacto, r["PoDm_mediana"], "error relativo <= 5 %",
                    "donde el modelo de placas es exacto, los dos coinciden",
                    ok)
    assert ok, (r["PoDm_mediana"], exacto)
