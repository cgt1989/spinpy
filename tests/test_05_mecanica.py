"""
Bloque 5 — Elasticidad: soluciones cerradas y cotas universales.

REFERENCIAS
  Hill R. The elastic behaviour of a crystalline aggregate. Proc Phys Soc A
  1952;65:349. Cotas de Voigt y Reuss: para CUALQUIER microestructura,
  C_Reuss <= C_eff <= C_Voigt en el sentido de las formas cuadraticas.
  Hashin Z, Shtrikman S. A variational approach to the theory of the elastic
  behaviour of multiphase materials. J Mech Phys Solids 1963;11:127-140.
  Cotas mas estrechas para materiales ISOTROPOS.
  Gibson LJ, Ashby MF. Cellular Solids, 2a ed., Cambridge 1997. E/E_s ~ C rho^n.
  Kumar S, Tan S, Zheng L, Kochmann DM. npj Comput Mater 2020;6. Escalado de
  la rigidez de los spinodoides con la densidad.
  Andreassen E, Andreasen CS. Comput Mater Sci 2014;83:488. El homogeneizador
  que se verifica.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  cubo macizo, homogeneizacion:  C_h = C_isotropo a 1e-9 relativo.
      El hexaedro trilineal reproduce EXACTAMENTE un campo de deformacion
      uniforme, asi que la unica fuente de error es la aritmetica.
  cubo macizo, ensayo deslizante:  E_app = E_s a 1e-9; eps_eff uniforme
      = sigma0/E_s y von Mises uniforme = sigma0 a 1e-7.
      El ensayo resuelve con un solver ITERATIVO (tol = 1e-8), asi que la
      uniformidad elemento a elemento solo puede exigirse al orden de esa
      tolerancia; 1e-7 es diez veces la del solver. (La primera version de
      esta prueba pedia 1e-9 y salio 1.0013e-9: era la tolerancia la que
      estaba mal declarada, no el codigo; se deja constancia.)
  cubo macizo, ensayo empotrado:  E_s <= E_app <= M_oed, con
      M_oed = E(1-nu)/((1+nu)(1-2nu)) el modulo edometrico (deformacion
      lateral impedida). Desigualdad dura, holgura 1e-9.
  spinodoide cualquiera:  C_h simetrico; C_h semidefinido positivo;
      rho*C_s - C_h semidefinido positivo (Voigt). Holgura 1e-9 * ||C||.
  spinodoide ISOTROPO (conos de 90 grados):  K_h <= 1.05 K_HS+ y
      G_h <= 1.05 G_HS+, con K_h, G_h el promedio de Hill del tensor. El 5 %
      es la holgura por la anisotropia residual de una celda finita: la cota
      es para el material isotropo, no para una realizacion.
  Gibson-Ashby:  el exponente n ajustado sobre rho en {0.25, 0.35, 0.50}
      cae en [1, 3]. Es una cota de plausibilidad, no un valor de
      literatura: Kumar 2020 reporta escalados cercanos a lineales a densidad
      media, y las espumas abiertas clasicas dan n ~ 2. Se REPORTA n.

RESULTADO DE LA PRIMERA EJECUCION, DEJADO AQUI A PROPOSITO
  Gibson-Ashby FALLA: n = 3.88 con rho = [0.26, 0.35, 0.51] y
  E/E_s = [0.017, 0.066, 0.224]. El exponente local es mayor a baja densidad
  (4.7 entre 0.26 y 0.35; 3.3 entre 0.35 y 0.51).
  Para descartar la discretizacion se repitio con la MISMA geometria (numero
  de onda 8*pi por celda, conos de 90 grados) a 20^3, 28^3 y 36^3 y rho en
  {0.31, 0.39, 0.51}: n = 3.74, 4.08, 4.11, y E/E_s a rho=0.31 converge
  0.035 -> 0.027 -> 0.025 (resultados/gibson_ashby_vs_celda.json). El
  exponente es robusto a la malla; NO es un artefacto del homogeneizador.
  Lo que ese estudio no separa es el tamano de RVE: con 8*pi caben ~4
  longitudes de onda por celda. El siguiente paso es repetirlo con mas
  longitudes de onda por celda (16*pi a 40^3 o mas) antes de atribuir
  n ~ 4 a la familia. La cota [1, 3] se declaro desde la literatura de
  espumas, no de spinodoides, y por el criterio de esta suite no se mueve:
  la prueba queda roja como hallazgo.

EXPLICADO EN 2026-09-11: ES PERCOLACION DE RIGIDEZ
  El paso que esta nota pedia -mas longitudes de onda por celda- se dio, y de
  paso se midio lo que faltaba. Ver `Estudio_Percolacion/INFORME.md`.

  El exponente se reproduce solo: un barrido independiente de 14 puntos a 32^3
  da n = 4.18, junto al 3.88 de aqui. No se va con la malla ni con el tamano de
  celda. Lo que lo explica es que la densidad hay que contarla desde un UMBRAL:

      E = C rho^n                      n = 4.18    residuo log 0.1496
      E = C (rho - 0.16)^t             t = 2.42    residuo log 0.0978
      E = C (rho - rho_c)^t libre      t = 1.66    residuo log 0.0693,
                                                   con rho_c = 0.222

  y sobre los TRES PUNTOS DE ESTA PRUEBA, la forma con el umbral medido ajusta
  a 0.0167 frente a 0.1029, con t = 2.06 -el exponente de espuma de celda
  abierta de manual-.

  El 0.16 no se ajusto a estos datos: es el umbral de CONECTIVIDAD, medido
  aparte contando voxeles portantes a 128^3 con tres semillas, sin tocar
  ninguna ecuacion de elasticidad. Que los datos elasticos prefieran 0.222 no
  es una contradiccion: el umbral de RIGIDEZ esta por encima del de
  conectividad, porque el conglomerado que percola primero esta lleno de
  extremos colgantes que no cargan. Corroboracion independiente: a rho = 0.25
  la estructura tiene el 95 % de su material en caminos que atraviesan y la
  homogeneizacion NO alcanza la puerta de residuo.

  LA PRUEBA SE QUEDA ROJA IGUAL, y a proposito: lo que mide -que un spinodoide
  no sigue E ~ rho^2- es CIERTO. Lo que ha cambiado no es el veredicto sino que
  ahora se sabe por que.
"""
import numpy as np
import pytest

from spinpy import generar_mascara
from spinpy.elastic import constantes_ingenieria, homogeneizar
from spinpy.resistencia import ensayo_compresion

BLOQUE = "05 Elasticidad (cerradas y cotas)"
E_S, NU = 1.0, 0.30
LAMBDA = E_S * NU / ((1 + NU) * (1 - 2 * NU))
MU = E_S / (2 * (1 + NU))
K_S = E_S / (3 * (1 - 2 * NU))
M_OED = E_S * (1 - NU) / ((1 + NU) * (1 - 2 * NU))


def C_isotropo(lam, mu):
    C = np.zeros((6, 6))
    C[:3, :3] = lam
    C[np.arange(3), np.arange(3)] += 2 * mu
    C[3:, 3:] = np.eye(3) * mu
    return C


def hill_KG(C):
    """Promedio de Hill (Voigt+Reuss)/2 de K y G para un C (6x6) de Voigt."""
    S = np.linalg.inv(C)
    Kv = (C[0, 0] + C[1, 1] + C[2, 2] + 2 * (C[0, 1] + C[0, 2] + C[1, 2])) / 9
    Gv = ((C[0, 0] + C[1, 1] + C[2, 2]) - (C[0, 1] + C[0, 2] + C[1, 2])
          + 3 * (C[3, 3] + C[4, 4] + C[5, 5])) / 15
    Kr = 1 / (S[0, 0] + S[1, 1] + S[2, 2] + 2 * (S[0, 1] + S[0, 2] + S[1, 2]))
    Gr = 15 / (4 * (S[0, 0] + S[1, 1] + S[2, 2])
               - 4 * (S[0, 1] + S[0, 2] + S[1, 2])
               + 3 * (S[3, 3] + S[4, 4] + S[5, 5]))
    return 0.5 * (Kv + Kr), 0.5 * (Gv + Gr)


def hashin_shtrikman_sup(f, K1, G1):
    """Cotas superiores de HS para solido (K1,G1) con fraccion f y poros."""
    f2 = 1 - f
    K = K1 + f2 / (-1 / K1 + f / (K1 + 4 * G1 / 3))
    G = G1 + f2 / (-1 / G1 + 2 * f * (K1 + 2 * G1) / (5 * G1 * (K1 + 4 * G1 / 3)))
    return K, G


def test_homogeneizacion_cubo_macizo(registro):
    Ch, info = homogeneizar(np.ones((12, 12, 12), bool), E_s=E_S, nu_s=NU,
                            vox_size=1.0)
    assert info["ok"], info["msg"]
    C = C_isotropo(LAMBDA, MU)
    err = float(np.max(np.abs(Ch - C)) / np.max(np.abs(C)))
    registro.anotar(BLOQUE, "homogeneizacion cubo macizo: max|Ch - C_iso|",
                    "Andreassen & Andreasen 2014", 0.0, err, 1e-9,
                    "exacto para deformacion uniforme", err <= 1e-9)
    assert err <= 1e-9
    ec = constantes_ingenieria(Ch)
    for k, esp in (("Ez", E_S), ("Gxy", MU), ("nu_xy", NU)):
        e = abs(ec[k] - esp) / esp
        registro.anotar(BLOQUE, f"constantes de ingenieria: {k}", "Hooke",
                        esp, ec[k], 1e-9, "relativo", e <= 1e-9)
        assert e <= 1e-9


def test_ensayo_cubo_macizo_deslizante(registro):
    s0 = 1e-3
    r = ensayo_compresion(np.ones((12, 12, 12), bool), [1.0] * 3, E_s=E_S,
                          nu_s=NU, sigma0=s0, apoyo="deslizante")
    assert r["ok"], r["msg"]
    e = abs(r["E_app"] - E_S) / E_S
    registro.anotar(BLOQUE, "ensayo deslizante cubo macizo: E_app", "Hooke",
                    E_S, r["E_app"], 1e-9, "relativo", e <= 1e-9)
    assert e <= 1e-9
    TOL_SOLVER = 1e-7          # 10 x la tolerancia del solver iterativo (1e-8)
    eps = np.asarray(r["eps_eff_solido"])
    d = float(np.max(np.abs(eps - s0 / E_S)) / (s0 / E_S))
    registro.anotar(BLOQUE, "eps_eff uniforme = sigma0/E_s (Pistoia uniaxial)",
                    "Pistoia 2002, def. energetica", s0 / E_S, float(eps.mean()),
                    TOL_SOLVER, "max desviacion relativa; solver tol 1e-8",
                    d <= TOL_SOLVER,
                    nota="en estado uniaxial eps_eff = |eps_zz|")
    assert d <= TOL_SOLVER
    vm = np.asarray(r["vm_solido"])
    d = float(np.max(np.abs(vm - s0)) / s0)
    registro.anotar(BLOQUE, "von Mises uniforme = sigma0", "def. von Mises",
                    s0, float(vm.mean()), TOL_SOLVER,
                    "max desviacion relativa; solver tol 1e-8", d <= TOL_SOLVER)
    assert d <= TOL_SOLVER


def test_ensayo_cubo_macizo_empotrado_entre_cotas(registro):
    r = ensayo_compresion(np.ones((12, 12, 12), bool), [1.0] * 3, E_s=E_S,
                          nu_s=NU, sigma0=1e-3, apoyo="empotrado")
    assert r["ok"], r["msg"]
    E = r["E_app"]
    ok = (E_S - 1e-9) <= E <= (M_OED + 1e-9)
    registro.anotar(BLOQUE, "ensayo empotrado: E_s <= E_app <= M_oed",
                    "modulo edometrico E(1-nu)/((1+nu)(1-2nu))", M_OED, E,
                    "desigualdad", f"[{E_S:.4f}, {M_OED:.4f}]", ok,
                    nota="la coaccion lateral solo puede rigidizar")
    assert ok, E


@pytest.mark.lento
def test_cotas_de_voigt_en_spinodoides(registro):
    """rho*C_s - C_h y C_h deben ser semidefinidos positivos (Hill 1952)."""
    Cs = C_isotropo(LAMBDA, MU)
    for rho in (0.30, 0.50):
        BW, _, _ = generar_mascara(rho=rho, wave_number=8 * np.pi,
                                   num_waves=300, thetas=[15, 15, 45],
                                   resolution=20, seed=1)
        Ch, info = homogeneizar(BW, E_s=E_S, nu_s=NU, vox_size=1.0 / 20)
        assert info["ok"], info["msg"]
        f = float(BW.mean())
        holg = 1e-9 * np.max(np.abs(Cs))
        sim = float(np.max(np.abs(Ch - Ch.T)))
        l_h = float(np.linalg.eigvalsh(Ch).min())
        l_v = float(np.linalg.eigvalsh(f * Cs - Ch).min())
        registro.anotar(BLOQUE, f"spinodoide rho={f:.2f}: Ch simetrico",
                        "Hill 1952", 0.0, sim, holg, "max|Ch-Ch^T|",
                        sim <= holg)
        registro.anotar(BLOQUE, f"spinodoide rho={f:.2f}: Ch >= 0 (Reuss)",
                        "Hill 1952", None, l_h, f">= -{holg:.1e}",
                        "autovalor minimo", l_h >= -holg)
        registro.anotar(BLOQUE, f"spinodoide rho={f:.2f}: rho*Cs - Ch >= 0 (Voigt)",
                        "Hill 1952", None, l_v, f">= -{holg:.1e}",
                        "autovalor minimo", l_v >= -holg)
        assert sim <= holg and l_h >= -holg and l_v >= -holg


@pytest.mark.lento
def test_hashin_shtrikman_isotropo(registro):
    """Conos de 90 grados = muestreo isotropo. K y G de Hill bajo HS+."""
    BW, _, _ = generar_mascara(rho=0.40, wave_number=8 * np.pi, num_waves=400,
                               thetas=[90, 90, 90], resolution=20, seed=2)
    Ch, info = homogeneizar(BW, E_s=E_S, nu_s=NU, vox_size=1.0 / 20)
    assert info["ok"], info["msg"]
    f = float(BW.mean())
    K, G = hill_KG(Ch)
    Ksup, Gsup = hashin_shtrikman_sup(f, K_S, MU)
    okK = 0 <= K <= 1.05 * Ksup
    okG = 0 <= G <= 1.05 * Gsup
    registro.anotar(BLOQUE, f"HS+: K_hill <= 1.05 K_HS+ (rho={f:.2f})",
                    "Hashin & Shtrikman 1963", Ksup, K, "5 %",
                    "cota superior isotropa", okK)
    registro.anotar(BLOQUE, f"HS+: G_hill <= 1.05 G_HS+ (rho={f:.2f})",
                    "Hashin & Shtrikman 1963", Gsup, G, "5 %",
                    "cota superior isotropa", okG)
    assert okK and okG, (K, Ksup, G, Gsup)


@pytest.mark.lento
def test_escalado_gibson_ashby(registro):
    rhos, Es = [], []
    for rho, seed in ((0.25, 11), (0.35, 12), (0.50, 13)):
        BW, _, _ = generar_mascara(rho=rho, wave_number=8 * np.pi,
                                   num_waves=400, thetas=[90, 90, 90],
                                   resolution=20, seed=seed)
        Ch, info = homogeneizar(BW, E_s=E_S, nu_s=NU, vox_size=1.0 / 20)
        assert info["ok"], info["msg"]
        rhos.append(float(BW.mean()))
        Es.append(constantes_ingenieria(Ch)["E_medio"])
    n, logC = np.polyfit(np.log(rhos), np.log(Es), 1)
    ok = 1.0 <= n <= 3.0
    registro.anotar(BLOQUE, "Gibson-Ashby: exponente n de E/E_s ~ rho^n",
                    "Gibson & Ashby 1997; Kumar 2020", None, float(n),
                    "1 <= n <= 3", "plausibilidad", ok,
                    nota="rho=%s  E/E_s=%s" % (
                        [round(r, 3) for r in rhos], [round(e, 4) for e in Es]))
    assert ok, n
