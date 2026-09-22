"""
fit.py — Ajuste de un spinodoide a un VOI por busqueda escalonada.

Port de `fitSpinodoidToVOI` (AppFinal_V2.m:3407-3720), incluidas las
correcciones G1, G2 y K2, que son la parte mas reciente del trabajo sobre la
app y la que motivo la version V2.

CORRECCION L1 — LA ORIENTACION SE MIDE, NO SE SUPONE
-----------------------------------------------------
Es la PRIMERA correccion que no viene de MATLAB: `AppFinal_V2.m` no la lleva.
Quien contraste contra MATLAB tiene que pasar `alinear_fabrica=False`, y asi lo
hace `validar_ajuste.py`.

  L1. `R_fit` salia de `dir_a_euler(m_voi['dir_principal'])`, que lleva el eje
      z DEL CANDIDATO sobre la direccion principal del VOI. Eso presupone que
      el eje z del candidato es su eje rigido, y no lo es: con muestreo por
      rechazo el eje rigido es el del cono ESTRECHO, y con muestreo equitativo
      la regla se INVIERTE. Ademas, con dos conos estrechos iguales
      —[15,15,X], que son 11 de los 14 presets— no hay eje rigido sino un
      PLANO, y su direccion principal es degenerada.

      Consecuencia medida sobre el VOI proximal de H4: el ajuste de menor error
      dejaba el eje rigido del candidato a 86 grados del hueso y el modulo
      aparente en la direccion de carga salia 57 veces menor. El DA coincidia
      (1.456 frente a 1.480) porque el DA es un ESCALAR y no distingue "rigido
      en Z" de "rigido en Y", y ningun termino del error mira la orientacion.

      -> `alinear_por_fabrica`: se genera el ganador sin rotar, se MIDE su
         tensor de fabrica, se construye la rotacion minima que lleva esa
         direccion sobre la del VOI, se regenera y se vuelve a medir para
         comprobar donde quedo.

CORRECCION L2 — LA ALINEACION SE IMPONE, NO SE INTENTA
-------------------------------------------------------
L1 se RENDIA en tres casos y devolvia entonces la rotacion clasica, que es
precisamente la que dejaba el eje perpendicular: si el candidato tenia la
fabrica degenerada en un plano (los once presets [15,15,X]), si tras rotar el
eje quedaba peor que antes (el generador remuestrea las ondas), o si el
candidato era casi isotropo. `avisos.py` lo midio despues en el banco: en 10 de
12 VOIs equinos el eje quedaba a mas de 30 grados.

  -> `alinear_marco_fabrica` (por omision, `alineacion="impuesta"`):
       * alinea el MARCO MIL entero —cada autovector del candidato con el
         correspondiente del VOI—, no un solo eje. Un candidato plano tiene
         normal bien definida aunque no tenga eje principal, y alinear el
         marco la coloca;
       * compara el eje que SIGNIFICA algo: la direccion principal, o la
         normal del plano si el candidato es plano;
       * itera —regenerar, medir, corregir— hasta quedar dentro de `tol_deg`
         o agotar `iter_max`;
       * y se queda con la MEJOR orientacion medida entre la clasica y las
         iteradas. Nunca devuelve algo peor que la clasica, y dice si llego a
         la tolerancia (`alineada`).
     `alineacion="eje"` conserva L1 tal cual.

LO QUE ESTAS CORRECCIONES NO ARREGLAN
--------------------------------------
El hueco de RIGIDEZ entre hueso y spinodoide es de dos causas, y la
orientacion solo toca una. La otra es de la familia: el exponente de
Gibson-Ashby medido para los spinodoides es n ~ 4.1, frente a n ~ 2 del hueso
trabecular. A rho ~ 0.26 eso es un factor de 3 a 6 que ninguna rotacion cierra.

LAS TRES ETAPAS
---------------
  A. densidad x numero de onda   (angulos fijos en el valor de partida)
  B. presets de angulos conicos  (con los mejores parametros de A)
  C. refinado de densidad y numero de ondas   (solo en modo completo)
  D. desempate entre finalistas  (solo si hay terminos mecanicos o caros)

Las correcciones que NO se pueden perder al portar:

  G2 — el numero de onda se recorre en TODO su rango util [8, 25] en el modo
       completo, con independencia de donde estuviera el deslizador. La version
       anterior usaba una ventana [base-5, base, base+5]; con el valor por
       defecto (15) eso daba [10 15 20] y NUNCA alcanzaba 25, de modo que en
       trabeculas muy finas el ajuste sobreestimaba Tb.Th en torno al 43 %.

  G1 — los angulos conicos se BARREN en la etapa B. Antes quedaban fijos en el
       valor de los deslizadores y el spinodoide salia casi isotropo
       (DA ~1.1-1.2), incapaz de alcanzar el DA 1.4-1.5 del hueso proximal
       equino.

  K2 — la lista de presets cubre las PERMUTACIONES, no solo un octante. thetas
       no es simetrico: el angulo grande marca el eje en torno al cual se
       concentran los vectores de onda, y por tanto el eje BLANDO. Sin las
       permutaciones solo se podia orientar la estructura de una manera.
       Ademas se desempata por orientacion entre los presets cuyo error ya
       esta dentro del 5 % del mejor (ver `_desempatar_theta`).

FUNCION OBJETIVO E INCERTIDUMBRE (correcciones O1 y U1)
-------------------------------------------------------
`pesos` y `distancia` pasan tal cual a `error.error_morfometrico` (ver su
cabecera). Donde se paga cada termino lo decide su coste: los de forma en cada
evaluacion, el Ellipsoid Factor y los mecanicos solo entre finalistas (etapa
D). Con `pesos=None` la busqueda es la de siempre, bit a bit.

`replicas` (5 por omision) mide el ganador con K semillas nuevas y reporta el
error medio +- sd y el suelo autoconsistente (`spinpy.incertidumbre`). No
cambia que candidato gana. `seleccion="robusta"` si lo cambia: elige entre los
finalistas por media + `k_robusto` sd sobre esas K semillas.

POR QUE LA ORIENTACION NO ES UN TERMINO DEL ERROR
--------------------------------------------------
Se probo y se descarto: sobre los tres VOI de H4 con 5 replicas mejoraba la
rigidez pero subia el error de las siete metricas de 0.00726 a 0.00990 (+36 %).
Un sumando nuevo compra su objetivo con el de los demas y no hay forma de
acotar el intercambio desde dentro de la suma. Aqui la orientacion es criterio
SECUNDARIO entre candidatos ya buenos: la morfometria queda acotada por
construccion y la formula del error no se toca.
"""

from __future__ import annotations

import time

import numpy as np

from . import incertidumbre, procedencia
from .error import (DISTANCIAS, TERMINOS_OPCIONALES, error_morfometrico,
                    fabric_cos, medidas_necesarias, peso_alineacion,
                    validar_pesos)
from .grf import (canonicalizar_thetas, euler_R, generar_mascara,
                  region_degenerada)
from .morphometry import metricas_forma, morfometria

SEMILLA_POR_DEFECTO = 20260720      # F12: resultados reproducibles

# Material base de la app (localBaseMaterial): matriz mineralizada.
E_S_DEF, NU_S_DEF = 20e9, 0.30

# Replicas del ganador para el informe de incertidumbre (U1).
REPLICAS_DEF = incertidumbre.K_MIN

# Claves de `metricas_mecanicas` que se copian a los dicts que ve el error.
CLAVES_MECANICAS = ("Ez_rel", "Ez_Ex", "E1_rel", "E2_rel", "E3_rel",
                    "G23_rel", "G13_rel", "G12_rel", "C_principal_rel")

# Metricas de las que el informe de incertidumbre da media +- sd.
CLAVES_INCERTIDUMBRE = ["BVTV", "BSBV", "TbTh", "TbSp", "TbN", "DA", "DA2"]


def metricas_mecanicas(BW, spacing, res_mec=16, E_s=E_S_DEF, nu_s=NU_S_DEF):
    """Constantes elasticas de una estructura, por homogeneizacion periodica.

    Ez_rel y Ez_Ex son los dos campos mecanicos historicos que
    `error_morfometrico` usa con `peso_mecanico > 0`: la rigidez axial
    normalizada y la anisotropia elastica. Ademas se devuelven, en los EJES
    MATERIALES y ordenadas E1 >= E2 >= E3, las constantes de ingenieria
    (E*_rel, G*_rel, nu*), el tensor normalizado en esos ejes
    (`C_principal_rel`, para la distancia log-euclidea), la clase de
    anisotropia y el residuo del solver.

    La estructura se REMUESTREA a `res_mec` antes de homogeneizar porque el
    coste crece con el cubo del lado: a 16^3 son ~6 s y a 32^3 ~99 s. En un
    desempate entre finalistas eso es la diferencia entre medio minuto y media
    hora.

    Devuelve {} si la homogeneizacion falla, y no ceros: `error_morfometrico`
    salta los campos ausentes y renormaliza por el peso realmente usado
    (invariante C6), mientras que un cero fingiria una estructura sin rigidez y
    arrastraria el ajuste hacia geometrias degeneradas.

    Ez_rel DEPENDE FUERTEMENTE DE `res_mec`, y hay que contar con ello. Medido
    sobre el mismo VOI sintetico (48^3, BV/TV 0.345): Ez_rel = 0.108 a 12^3 y
    0.064 a 16^3. No es un error de calculo sino el mismo fenomeno que hace no
    convergente a E_app en el ensayo de compresion —al remuestrear cambia que
    trabeculas quedan conectadas— y por eso el numero **no es citable como
    rigidez del espécimen**.

    Para lo que si sirve es para comparar: en el ajuste, VOI y candidatos se
    homogeneizan todos al MISMO `res_mec`, de modo que el sesgo actua sobre los
    dos lados y se cancela en buena parte al restar. Es exactamente el
    argumento de la correccion C4, y es la unica razon por la que un `res_mec`
    barato es aceptable aqui y no lo seria para reportar un modulo.
    """
    from .elastic import (clase_anisotropia, constantes_ingenieria,
                          constantes_principales, homogeneizar, remuestrear_bw)

    bw, sp = remuestrear_bw(np.asarray(BW, bool), spacing, int(res_mec))
    C, info = homogeneizar(bw, E_s, nu_s, vox_size=sp)
    if not info.get("ok"):
        return {}
    ec = constantes_ingenieria(C)
    if not np.isfinite(ec["Ez"]) or ec["Ex"] <= 0:
        return {}
    out = {"Ez_rel": float(ec["Ez"] / E_s),
           "Ez_Ex": float(ec["Ez"] / ec["Ex"]),
           "Ez_Pa": float(ec["Ez"]), "Ex_Pa": float(ec["Ex"])}
    try:
        p = constantes_principales(C)
        out.update({
            "E1_rel": float(p["E1"] / E_s), "E2_rel": float(p["E2"] / E_s),
            "E3_rel": float(p["E3"] / E_s),
            "G23_rel": float(p["G23"] / E_s), "G13_rel": float(p["G13"] / E_s),
            "G12_rel": float(p["G12"] / E_s),
            "nu12": float(p["nu12"]), "nu13": float(p["nu13"]),
            "nu23": float(p["nu23"]),
            "C_principal_rel": (np.asarray(p["C_principal"]) / E_s).tolist(),
            "desviacion_ortotropa": float(p["desviacion_ortotropa"]),
            "clase_anisotropia": clase_anisotropia(C)["clase"],
        })
    except (np.linalg.LinAlgError, ValueError):
        pass
    if "residuo_rel" in info:
        out["residuo_homog"] = float(info["residuo_rel"])
    return out

# Tope de la rejilla de densidades. DESVIACION DELIBERADA RESPECTO A MATLAB,
# que usa 0.80 (`min(0.80, BVTV*1.3)` en fitSpinodoidToVOI).
#
# POR QUE SE SUBE. Con 0.80, cinco de los doce VOIs equinos del lote quedaban
# pegados al tope y no podian alcanzarse: los distales tienen BV/TV de 0.82 a
# 0.92 y H1 distal llega a 0.9164. El error de ese ajuste se disparaba a 0.162,
# veinte veces el de los proximales, simplemente porque la busqueda no podia
# llegar donde estaba la respuesta.
#
# POR QUE 0.95 Y NO MAS. Medido sobre el generador (64^3, wave 12pi, 700 ondas),
# la fase PORO se fragmenta al subir la densidad:
#
#     rho pedida   0.70    0.80    0.85    0.90    0.95    0.98
#     poro conexo  99.9%   83.4%   23.1%    7.7%    7.5%    6.8%
#     Tb.Th (mm)   0.141   0.204   0.265   0.384   0.724   1.695
#
# El hueso trabecular es BICONTINUO —la medula forma una red conectada—, asi
# que por encima de ~0.85 la estructura deja de ser trabecula y pasa a ser
# cortical con poros aislados. A 0.98 el Tb.Th de 1.7 mm ya es un bloque macizo.
#
# Aun asi el tope se pone en 0.95 y no en 0.85, porque LOS VOIs REALES SE
# COMPORTAN IGUAL: medida su fase poro, H1 medio (BV/TV 0.893) conserva un
# 36.5% conexo y H1 distal (0.916) un 41.9%, con mas de 145 poros aislados cada
# uno. El modelo no esta fallando al representarlos; a esas densidades el
# espécimen real tampoco es bicontinuo. Acotar en 0.85 excluiria especimenes
# que existen.
#
# El ajuste AVISA (`bicontinuo` en el resultado) cuando el candidato ganador
# tiene la fase poro fragmentada, para que la decision de si ese ajuste es
# interpretable la tome una persona.
RHO_MAX = 0.95
RHO_MIN = 0.20


# ---------------------------------------------------------------------------

def dir_a_euler(dir_vec):
    """Port de dirVectorToEulerApprox (AppFinal_V2.m:5414).

    Angulos de Euler que llevan el eje z sobre `dir_vec`. Se fuerza la
    componente z positiva porque la direccion principal es un EJE: su signo es
    arbitrario y quedarse con el hemisferio superior evita rotaciones de 180
    grados que no cambian la geometria.
    """
    d = np.asarray(dir_vec, dtype=float).ravel()
    if d.size != 3 or np.linalg.norm(d) == 0:
        return 0.0, 0.0, 0.0
    d = d / np.linalg.norm(d)
    if d[2] < 0:
        d = -d
    theta = np.arccos(np.clip(d[2], -1.0, 1.0))
    phi = np.arctan2(d[1], d[0])
    rx = 0.0
    ry = float(np.clip(np.degrees(theta), -90.0, 90.0))
    rz = float(np.clip(np.degrees(phi), -180.0, 180.0))
    return rx, ry, rz


def _rot_entre_ejes(a, b):
    """Rotacion minima que lleva el EJE a sobre el eje b (Rodrigues).

    Son ejes, no vectores: el signo de un autovector es arbitrario, asi que se
    elige el sentido de `a` que da el giro mas corto. Sin eso se puede obtener
    una rotacion de 180 grados que es geometricamente equivalente pero mete un
    giro innecesario en las demas componentes.
    """
    a = np.asarray(a, float).ravel()
    b = np.asarray(b, float).ravel()
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return np.eye(3)
    a, b = a / na, b / nb
    if float(a @ b) < 0:
        a = -a
    v = np.cross(a, b)
    s = float(np.linalg.norm(v))
    c = float(a @ b)
    if s < 1e-12:
        return np.eye(3)          # ya alineados (el caso opuesto se descarto)
    vx = np.array([[0.0, -v[2], v[1]],
                   [v[2], 0.0, -v[0]],
                   [-v[1], v[0], 0.0]])
    return np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))


def _angulo_ejes(a, b):
    """Angulo en grados entre dos EJES (0 a 90)."""
    a = np.asarray(a, float).ravel()
    b = np.asarray(b, float).ravel()
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return float("nan")
    c = abs(float((a / na) @ (b / nb)))
    return float(np.degrees(np.arccos(min(1.0, c))))


def alinear_por_fabrica(m_voi, res, sp_cand, dens, wave, nw, thetas,
                        esquema="rechazo", seed=SEMILLA_POR_DEFECTO,
                        precision="f32", da_min=1.10, da2_max=1.06):
    """Rotacion que alinea la fabrica MEDIDA del candidato con la del VOI.

    Es la correccion L1, que se conserva como `alineacion="eje"`. Por omision
    el ajuste usa `alinear_marco_fabrica` (L2), que no se rinde en los casos
    en que esta si (ver la cabecera del modulo).

    CORRECCION L1 — POR QUE NO VALE SUPONER QUE EL EJE RIGIDO ES Z
    ---------------------------------------------------------------
    Hasta aqui, `R_fit` salia de `dir_a_euler(m_voi['dir_principal'])`, que
    construye la rotacion que lleva **el eje z del spinodoide** sobre la
    direccion principal del VOI. Eso presupone que el eje z del candidato es su
    eje rigido, y NO LO ES en general:

      * Con muestreo por RECHAZO, el eje rigido es el del cono ESTRECHO: el
        cono ancho se lleva casi todas las ondas por angulo solido, el campo
        oscila rapido en esa direccion y esa direccion queda blanda. Medido:
        thetas [45,45,15] -> rigido en Z; [45,15,45] -> rigido en Y;
        [15,45,45] -> rigido en X. Tres de tres.
      * Con muestreo EQUITATIVO la regla SE INVIERTE: el cono estrecho recibe
        su tercio de ondas concentrado junto a su eje, y es ese eje el que se
        ablanda. Medido: [45,45,15] -> rigido en Y; [90,90,15] -> rigido en X.
      * Y con dos conos estrechos iguales ([15,15,X]) no hay eje rigido sino un
        PLANO rigido: la direccion principal es degenerada y alinearla no
        significa nada.

    Consecuencia medida sobre el VOI proximal de H4: el ajuste de menor error
    dejaba el eje rigido del candidato a 86 grados del que tiene el hueso, y el
    modulo aparente en la direccion de carga salia 57 veces menor que el del
    VOI. El DA coincidia (1.456 frente a 1.480) porque el DA es un ESCALAR y no
    distingue "rigido en Z" de "rigido en Y".

    Como ningun termino de `error_morfometrico` depende de la orientacion, la
    correccion no cambia que candidato gana: solo lo coloca bien. Por eso se
    aplica DESPUES de elegir los parametros y cuesta dos generaciones extra.

    COMO SE HACE
      1. Se genera el candidato SIN rotar y se mide su tensor de fabrica.
      2. Se construye la rotacion minima que lleva esa direccion medida sobre
         la del VOI.
      3. Se regenera con ella y se VUELVE A MEDIR, porque el generador
         remuestrea las ondas y el resultado no es una rotacion exacta del
         anterior. El angulo que queda se devuelve para que conste.

    CUANDO NO SE APLICA
      Si el VOI o el candidato son casi isotropos (DA < `da_min`) la direccion
      principal es ruido, y si el candidato tiene la fabrica degenerada en un
      plano (DA2 < `da2_max`) no hay un eje que alinear. En los dos casos se
      devuelve la rotacion clasica y se dice por que.

    Devuelve (R, info).
    """
    u_voi = np.asarray(m_voi.get("dir_principal", [0.0, 0.0, 1.0]), float)
    rx, ry, rz = dir_a_euler(u_voi)
    R_clasica = euler_R(rx, ry, rz)
    info = {"metodo": "clasica (eje z del candidato sobre la del VOI)",
            "aplicada": False, "motivo": "", "angulo_antes_deg": float("nan"),
            "angulo_despues_deg": float("nan"),
            "DA_candidato": float("nan"), "DA2_candidato": float("nan")}

    da_voi = float(m_voi.get("DA", np.nan))
    if not np.isfinite(da_voi) or da_voi < da_min:
        info["motivo"] = (f"el VOI es casi isotropo (DA {da_voi:.3f} < "
                          f"{da_min}): su direccion principal es ruido")
        return R_clasica, info

    # 1. candidato SIN rotar
    BW0, _, _ = generar_mascara(
        resolution=res, wave_number=float(wave) * np.pi, num_waves=int(nw),
        thetas=np.clip(np.round(thetas), 0, 90), rho=float(dens),
        R=np.eye(3), esquema=esquema, seed=seed, precision=precision)
    m0 = morfometria(BW0, sp_cand)
    u0 = np.asarray(m0.get("dir_principal", [0.0, 0.0, 1.0]), float)
    da0 = float(m0.get("DA", np.nan))
    da2_0 = float(m0.get("DA2", np.nan))
    info["DA_candidato"] = da0
    info["DA2_candidato"] = da2_0

    if not np.isfinite(da0) or da0 < da_min:
        info["motivo"] = (f"el candidato es casi isotropo (DA {da0:.3f}): no "
                          f"hay eje que alinear")
        return R_clasica, info
    if np.isfinite(da2_0) and da2_0 < da2_max:
        info["motivo"] = (
            f"la fabrica del candidato es degenerada (DA2 {da2_0:.3f} < "
            f"{da2_max}): es rigido en un PLANO, no en un eje, y ninguna "
            f"rotacion lo convierte en axial")
        return R_clasica, info

    # 2. rotacion que lleva la fabrica medida sobre la del VOI
    R = _rot_entre_ejes(u0, u_voi)

    # 3. regenerar y volver a medir: el generador remuestrea, asi que hay que
    #    comprobar donde queda de verdad el eje, no darlo por hecho.
    BW1, _, _ = generar_mascara(
        resolution=res, wave_number=float(wave) * np.pi, num_waves=int(nw),
        thetas=np.clip(np.round(thetas), 0, 90), rho=float(dens),
        R=R, esquema=esquema, seed=seed, precision=precision)
    m1 = morfometria(BW1, sp_cand)
    u1 = np.asarray(m1.get("dir_principal", [0.0, 0.0, 1.0]), float)

    info.update({
        "metodo": "fabrica medida (correccion L1)",
        "aplicada": True,
        "angulo_antes_deg": _angulo_ejes(u0, u_voi),
        "angulo_despues_deg": _angulo_ejes(u1, u_voi),
        "dir_candidato_sin_rotar": [float(x) for x in u0],
        "dir_candidato_rotado": [float(x) for x in u1],
        "dir_voi": [float(x) for x in u_voi],
    })
    # Si tras rotar el eje esta MAS lejos que con la rotacion clasica, no se
    # fuerza: el remuestreo puede haber cambiado la fabrica del candidato.
    if info["angulo_despues_deg"] > info["angulo_antes_deg"]:
        info["aplicada"] = False
        info["motivo"] = ("tras rotar, el eje quedo mas desalineado que antes "
                          "(el generador remuestrea las ondas); se conserva la "
                          "rotacion clasica")
        return R_clasica, info
    return R, info


def _rot_marcos(Va, Vb):
    """Rotacion propia que lleva cada columna de Va sobre la de Vb (como ejes).

    El signo de los autovectores es libre; se invierte el tercero si hace
    falta para que el determinante sea +1 y la rotacion no refleje.
    """
    Va = np.asarray(Va, float)
    Vb = np.asarray(Vb, float)
    s = np.sign(np.linalg.det(Va) * np.linalg.det(Vb)) or 1.0
    return Vb @ np.diag([1.0, 1.0, s]) @ Va.T


def alinear_marco_fabrica(m_voi, generar, medir, R_clasica, m_clasica=None,
                          da_min=1.10, da2_plano=1.06, tol_deg=10.0,
                          iter_max=3):
    """Correccion L2: impone la orientacion del marco MIL del VOI al candidato.

    generar(R) -> mascara del ganador con la rotacion R
    medir(mascara) -> dict de `morfometria` (con `eigenvectors`, `DA`, `DA2`)
    m_clasica : la medida del ganador con `R_clasica`, si ya se tiene (el
        ajuste la trae de la busqueda); si no, se genera.

    Sirve a las dos familias: en el dual-lattice la rotacion es exacta y suele
    bastar una iteracion; en el spinodoide el generador remuestrea las ondas al
    rotar y la iteracion corrige lo que quede.

    El eje que se compara es la direccion principal (autovector del menor
    autovalor del MIL, el de MIL maximo) salvo que el candidato sea PLANO
    (DA2 < `da2_plano`: los dos primeros autovalores iguales), en cuyo caso se
    compara la NORMAL del plano, que es lo que ese candidato si tiene bien
    definido.

    No se aplica —y se dice— si el VOI o el candidato son casi isotropos: sin
    eje no hay orientacion que imponer ni que perder.

    Devuelve (R, info). `info` conserva las claves de L1 (`aplicada`,
    `motivo`, `angulo_antes_deg`, `angulo_despues_deg`, `metodo`) y anade
    `alineada`, `iteraciones`, `eje_comparado` e `historial_deg`.
    """
    info = {"metodo": "marco MIL impuesto (correccion L2)", "aplicada": False,
            "alineada": False, "motivo": "", "tol_deg": float(tol_deg),
            "angulo_antes_deg": float("nan"),
            "angulo_despues_deg": float("nan"), "iteraciones": 0,
            "eje_comparado": "", "historial_deg": [],
            "DA_candidato": float("nan"), "DA2_candidato": float("nan")}

    da_voi = float(m_voi.get("DA", np.nan))
    V_voi = np.asarray(m_voi.get("eigenvectors", np.eye(3)), float)
    if (not np.isfinite(da_voi) or da_voi < da_min
            or not m_voi.get("MIL_valid", True) or V_voi.shape != (3, 3)):
        info["motivo"] = (f"el VOI es casi isotropo o sin fabrica valida (DA "
                          f"{da_voi:.3f} < {da_min}): no hay orientacion que "
                          f"imponer")
        return R_clasica, info

    m0 = medir(generar(np.eye(3)))
    da0 = float(m0.get("DA", np.nan))
    da2_0 = float(m0.get("DA2", np.nan))
    info["DA_candidato"], info["DA2_candidato"] = da0, da2_0
    if not np.isfinite(da0) or da0 < da_min:
        info["motivo"] = (f"el candidato es casi isotropo (DA {da0:.3f}): su "
                          f"orientacion no cambia nada")
        return R_clasica, info

    plano = bool(np.isfinite(da2_0) and da2_0 < da2_plano)
    col = 2 if plano else 0
    info["eje_comparado"] = "normal del plano" if plano else "direccion principal"

    def angulo(m):
        V = np.asarray(m.get("eigenvectors", np.eye(3)), float)
        return _angulo_ejes(V[:, col], V_voi[:, col])

    if m_clasica is None:
        m_clasica = medir(generar(R_clasica))
    candidatos = [(angulo(m_clasica), np.asarray(R_clasica, float), "clasica")]
    info["angulo_antes_deg"] = candidatos[0][0]

    V0 = np.asarray(m0.get("eigenvectors", np.eye(3)), float)
    R = _rot_marcos(V0, V_voi)
    for it in range(int(iter_max)):
        m_it = medir(generar(R))
        a = angulo(m_it)
        candidatos.append((a, R, f"iteracion {it + 1}"))
        info["iteraciones"] = it + 1
        if a <= tol_deg:
            break
        # Correccion desde donde quedo de verdad: el marco medido con R se
        # lleva al del VOI y se compone con R.
        V_it = np.asarray(m_it.get("eigenvectors", np.eye(3)), float)
        R = _rot_marcos(V_it, V_voi) @ R

    info["historial_deg"] = [float(c[0]) for c in candidatos]
    validos = [c for c in candidatos if np.isfinite(c[0])]
    mejor = min(validos, key=lambda c: c[0]) if validos else candidatos[0]
    info["angulo_despues_deg"] = float(mejor[0])
    info["aplicada"] = mejor[2] != "clasica"
    info["alineada"] = bool(np.isfinite(mejor[0]) and mejor[0] <= tol_deg)
    if not info["alineada"]:
        info["motivo"] = (f"tras {info['iteraciones']} iteracion(es) el eje "
                          f"quedo a {mejor[0]:.1f} grados (> {tol_deg}); se "
                          f"usa la mejor orientacion medida ({mejor[2]})")
    return mejor[1], info


def longitud_caracteristica(voi_shape, spacing):
    """Port de localCharacteristicLength: lado fisico medio del VOI, en mm.

    Correccion F1: el spinodoide se genera en un cubo unidad adimensional. Si
    no se escala al tamano fisico del VOI antes de medirlo, Tb.Th, Tb.Sp y
    BS/BV salen en unidades distintas y restarlas no significa nada.
    """
    lc = float(np.mean(np.asarray(voi_shape, float) *
                       np.asarray(spacing, float).ravel()))
    return lc if np.isfinite(lc) and lc > 0 else 1.0


def resolucion_comparacion(voi_shape):
    """Port de localVoxelCountForComparison: 24 <= n <= 96, tomado del VOI."""
    s = np.asarray(voi_shape, dtype=int)
    if s.size != 3:
        return 64
    return int(np.clip(s.min(), 24, 96))


# ---------------------------------------------------------------------------

def _desempatar_theta(presets, errores, metricas, err_in, theta_in, m_in,
                      m_voi, tol=1.05):
    """Port de localDesempatarTheta (K2).

    Primero aplica el criterio de siempre —minimo error—. Despues, SOLO si
    desempatar por orientacion tiene sentido en este VOI, busca entre los
    candidatos cuyo error esta dentro del `tol` del mejor y se queda con el que
    mejor alinea su eje principal con el del VOI.

    No cuesta ninguna evaluacion extra: los presets ya estan todos evaluados, y
    UNA VEZ CADA UNO. Ese detalle importa. Comparar contra un error que ya fue
    elegido como minimo —y por tanto es una realizacion afortunada— hace
    inservible cualquier cota de este tipo: una tirada fresca de los mismos
    parametros la excede casi siempre (maldicion del ganador). Aqui no se
    reevalua nada.
    """
    info = {"n_elegibles": 0, "cos_ganador": np.nan, "cos_min_err": np.nan,
            "aplicado": False, "activo": False}
    t_out, e_out, m_out = theta_in, err_in, m_in

    for k, e in enumerate(errores):
        if e < e_out:
            e_out, t_out, m_out = e, presets[k], metricas[k]

    if peso_alineacion(m_voi) <= 0:
        return t_out, e_out, m_out, info
    info["activo"] = True
    info["cos_min_err"] = fabric_cos(m_voi, m_out)

    validos = [e for e in ([err_in] + list(errores)) if np.isfinite(e)]
    if not validos:
        return t_out, e_out, m_out, info
    e_best = min(validos)
    if not np.isfinite(e_best) or e_best <= 0:
        return t_out, e_out, m_out, info

    mejor_cos = fabric_cos(m_voi, m_out)
    if np.isnan(mejor_cos):
        mejor_cos = -np.inf
    n_eleg = 0

    for k, e in enumerate(errores):
        if not np.isfinite(e) or e > tol * e_best:
            continue
        n_eleg += 1
        c = fabric_cos(m_voi, metricas[k])
        if np.isnan(c):
            continue
        if c > mejor_cos:
            mejor_cos = c
            t_out, e_out, m_out = presets[k], e, metricas[k]
            info["aplicado"] = True

    info["n_elegibles"] = n_eleg
    info["cos_ganador"] = mejor_cos
    return t_out, e_out, m_out, info


# ---------------------------------------------------------------------------

def fraccion_poro_conexa(BW):
    """Fraccion del poro que pertenece a su mayor componente conexa (caras).

    Mide si la estructura sigue siendo BICONTINUA. En hueso trabecular la
    medula forma una red conectada; cuando este valor cae, la estructura ha
    pasado a ser un solido con poros aislados y deja de ser trabecula, por bien
    que coincida el BV/TV.
    """
    from scipy import ndimage
    poro = ~np.asarray(BW, dtype=bool)
    n_poro = int(poro.sum())
    if n_poro == 0:
        return 0.0, 0
    c6 = ndimage.generate_binary_structure(3, 1)
    lab, n = ndimage.label(poro, structure=c6)
    if n == 0:
        return 0.0, 0
    t = np.bincount(lab.ravel()); t[0] = 0
    return float(t.max() / n_poro), int(n)


def validar_opciones(pesos, distancia, replicas, seleccion, alineacion):
    """Comprueba las opciones de objetivo, incertidumbre y alineacion.

    Compartida con `fit_dual`: una opcion mal escrita tiene que fallar al
    empezar, no despues de media hora de busqueda.
    """
    pesos = validar_pesos(pesos)
    if distancia not in DISTANCIAS:
        raise ValueError(f"distancia debe ser una de {DISTANCIAS}")
    K = incertidumbre.validar_replicas(replicas)
    if seleccion not in incertidumbre.SELECCIONES:
        raise ValueError(f"seleccion debe ser una de "
                         f"{incertidumbre.SELECCIONES}")
    if seleccion == "robusta" and K == 0:
        raise ValueError(f"la seleccion robusta necesita replicas >= "
                         f"{incertidumbre.K_MIN}")
    if alineacion not in ("impuesta", "eje"):
        raise ValueError("alineacion debe ser 'impuesta' o 'eje'")
    return pesos, K


def claves_incertidumbre(pesos, usa_mec):
    claves = list(CLAVES_INCERTIDUMBRE)
    for k in (pesos or {}):
        if TERMINOS_OPCIONALES.get(k, {}).get("tipo") == "escalar":
            claves.append(k)
    if usa_mec:
        claves += ["Ez_rel", "Ez_Ex"]
    return list(dict.fromkeys(claves))


def ajustar_spinodoide(VOI, spacing, modo="completo", m_voi=None,
                       esquema="rechazo", precision="f32", seed=SEMILLA_POR_DEFECTO,
                       num_waves=700, resolucion=None, progreso=None,
                       rho_max=RHO_MAX, peso_mecanico=0.0, res_mec=16,
                       n_finalistas=5, E_s=E_S_DEF, nu_s=NU_S_DEF,
                       alinear_fabrica=True, pesos=None,
                       distancia="wasserstein", replicas=REPLICAS_DEF,
                       seleccion="minimo", k_robusto=1.0,
                       alineacion="impuesta"):
    """Busca los parametros de spinodoide que mejor reproducen el VOI.

    modo : 'completo' (etapas A+B+C, 14 presets) | 'rapido' (A+B, 4 presets)

    alinear_fabrica : al terminar la busqueda, la orientacion del ganador se
        corrige usando su fabrica MEDIDA en vez de suponer que su eje rigido es
        z. No cambia que candidato gana —ningun termino del error depende de la
        orientacion— solo lo coloca bien. Ponerlo a False reproduce el
        comportamiento anterior, que es lo que hay que usar para contrastar
        contra `AppFinal_V2.m`, porque MATLAB no lleva esta correccion.

    alineacion : 'impuesta' (L2, por omision: marco MIL entero, iterado, nunca
        peor que la clasica) | 'eje' (L1, la de antes). Ver la cabecera.

    peso_mecanico : si es > 0 se anade una ETAPA D que reordena los
        `n_finalistas` mejores candidatos incluyendo Ez_rel y Ez_Ex, obtenidos
        por homogeneizacion a `res_mec`.

        POR QUE UNA ETAPA APARTE Y NO DENTRO DE LA BUSQUEDA. Homogeneizar es
        entre tres y cuatro ordenes de magnitud mas caro que medir morfometria:
        las ~53 evaluaciones del modo completo pasarian de segundos a horas. La
        etapa D paga la homogeneizacion solo donde puede cambiar la decision
        —entre candidatos que ya son buenos morfometricamente— y son
        `n_finalistas`+1 homogeneizaciones en total, contando la del VOI.

        LO QUE ESTO NO ES. No es una busqueda con el termino mecanico dentro:
        si el optimo mecanico esta en una region que la busqueda morfometrica
        descarto pronto, la etapa D no lo encuentra. Es un desempate informado,
        y `traza_mecanica` deja ver cuanto reordeno para que se pueda juzgar si
        el desempate estaba decidiendo algo o solo confirmando.

    pesos, distancia : terminos opcionales del objetivo (correccion O1, ver
        `error.error_morfometrico`). Los de forma se miden en cada evaluacion;
        el Ellipsoid Factor y los mecanicos (C_logE, E*_rel, G*_rel) activan
        la etapa D y solo se pagan entre finalistas. Ojo con el EF: son varios
        minutos por estructura, tambien sobre el VOI.

    replicas : K semillas nuevas con las que se mide el ganador para el
        informe `incertidumbre` (0 lo desactiva; 1 a 4 es un error).

    seleccion : 'minimo' (el de siempre) | 'robusta': entre los finalistas,
        el de menor media + `k_robusto` sd del error sobre las K semillas.

    precision : 'f32' por defecto en la busqueda. Son ~50 evaluaciones y la
        diferencia frente a f64 no cambia ningun voxel en los casos medidos
        (ver spinpy.grf._evaluar_campo). El candidato ganador se REMIDE al
        final en f64, de modo que las metricas que se reportan son exactas.

    Devuelve un dict con los parametros ganadores, sus metricas y la traza.
    """
    pesos, K = validar_opciones(pesos, distancia, replicas, seleccion,
                                alineacion)
    nec = medidas_necesarias(pesos)
    forma = {k: True for k in nec["forma"]}
    caro = {k: True for k in nec["caro"]}
    usa_mec = bool(peso_mecanico > 0 or nec["mecanico"])
    etapa_d = bool(usa_mec or caro)

    VOI = np.asarray(VOI, dtype=bool)
    spacing = np.asarray(spacing, dtype=float).ravel()
    if m_voi is None:
        m_voi = morfometria(VOI, spacing)
    if forma:
        m_voi = dict(m_voi)
        m_voi.update(metricas_forma(VOI, spacing, **forma))

    lc = longitud_caracteristica(VOI.shape, spacing)
    res = int(resolucion) if resolucion else resolucion_comparacion(VOI.shape)
    sp_cand = np.full(3, lc / res)

    bvtv = float(m_voi["BVTV"])
    rx, ry, rz = dir_a_euler(m_voi["dir_principal"])
    R_fit = euler_R(rx, ry, rz)

    thetas0 = np.array([15.0, 15.0, 15.0])

    if modo == "rapido":
        dens = np.linspace(max(RHO_MIN, bvtv * 0.85), min(rho_max, bvtv * 1.15), 3)
        waves = np.unique(np.clip([10.0, 15.0, 20.0], 8, 25))
        presets = [thetas0, [15, 15, 45], [15, 45, 15], [45, 15, 15]]
        refinar = False
    else:
        # G2: rango util completo, no una ventana alrededor del deslizador
        dens = np.linspace(max(RHO_MIN, bvtv * 0.7), min(rho_max, bvtv * 1.3), 5)
        waves = np.unique(np.round(np.linspace(8, 25, 6)))
        # G1 + K2: familia [15 15 X] y sus permutaciones
        presets = [thetas0, [20, 20, 20],
                   [15, 15, 30], [15, 15, 45], [15, 15, 60],
                   [15, 30, 15], [15, 45, 15], [15, 60, 15],
                   [30, 15, 15], [45, 15, 15], [60, 15, 15],
                   [15, 30, 50], [15, 45, 60], [45, 60, 15]]
        refinar = True

    n_total = len(dens) * len(waves) + len(presets) + (9 if refinar else 0)
    traza = []
    t0 = time.time()
    contador = {"i": 0}

    def generar(d, w, nw, th, R, prec, s=None):
        """w es el MULTIPLO DE PI (deslizador); aqui y solo aqui se pasa a rad."""
        BW, _, _ = generar_mascara(
            resolution=res, wave_number=float(w) * np.pi, num_waves=int(nw),
            thetas=np.clip(np.round(th), 0, 90), rho=float(d), R=R,
            esquema=esquema, seed=seed if s is None else int(s),
            precision=prec)
        return BW

    def medir(BW, completo=False):
        m = morfometria(BW, sp_cand)
        if forma:
            m.update(metricas_forma(BW, sp_cand, **forma))
        if completo:
            if caro:
                m.update(metricas_forma(BW, sp_cand, **caro))
            if usa_mec:
                mec = metricas_mecanicas(BW, sp_cand, res_mec, E_s, nu_s)
                m.update(mec)
        return m

    def err(ma, mb, con_mec=False):
        return error_morfometrico(ma, mb,
                                  peso_mecanico=peso_mecanico if con_mec else 0.0,
                                  pesos=pesos, distancia=distancia)

    def pr_extra(i, n, etapa):
        if progreso:
            progreso(contador["i"], n_total, f"{etapa} ({i + 1}/{n})")

    def evaluar(d, w, nw, th, etapa):
        contador["i"] += 1
        if progreso:
            progreso(contador["i"], n_total, etapa)
        BW = generar(d, w, nw, th, R_fit, precision)
        m = medir(BW)
        e, n_us = err(m_voi, m)
        th_g = np.clip(np.round(th), 0, 90)
        traza.append({"etapa": etapa, "dens": float(d),
                      "wave_number_pi": float(w),
                      "wave_number_rad": float(w) * np.pi,
                      "nw": int(nw), "thetas": [float(t) for t in np.ravel(th)],
                      "degenerado": bool(
                          region_degenerada(th_g, esquema)["degenerado"]),
                      "err": float(e), "n_terminos": int(n_us)})
        return m, e

    def finalistas_de_traza(n):
        # Los finalistas salen de la traza, no de una lista aparte: asi el
        # criterio es exactamente "los mejores por morfometria" y se ve en la
        # propia traza cuales fueron. Se descartan duplicados de parametros,
        # que abundan porque las etapas comparten puntos; los thetas se
        # comparan en forma CANONICA, porque dos nombres de la region
        # degenerada son la misma estructura.
        vistos, fin = set(), []
        for t in sorted(traza, key=lambda t: t["err"]):
            clave = (round(t["dens"], 9), round(t["wave_number_pi"], 9),
                     t["nw"], tuple(canonicalizar_thetas(
                         np.clip(np.round(t["thetas"]), 0, 90), esquema)))
            if clave in vistos:
                continue
            vistos.add(clave)
            fin.append(t)
            if len(fin) >= int(n):
                break
        return fin

    # ---------- ETAPA A ----------
    best = {"err": np.inf, "d": dens[0], "w": waves[0], "nw": num_waves,
            "th": thetas0, "m": None}
    for d in dens:
        for w in waves:
            m, e = evaluar(d, w, num_waves, thetas0, "A")
            if e < best["err"]:
                best.update(err=e, d=d, w=w, nw=num_waves, m=m)

    # ---------- ETAPA B ----------
    thB, errB, mB = [], [], []
    vistos_B = {}
    for th in presets:
        th_try = np.clip(np.round(np.asarray(th, float)), 15, 60)
        # Un preset de la region degenerada es la estructura isotropa con
        # otro nombre: se reutiliza su evaluacion en vez de gastar otra y de
        # que dos "candidatos" identicos compitan en el desempate.
        clave = tuple(canonicalizar_thetas(th_try, esquema))
        if clave in vistos_B:
            m, e = vistos_B[clave]
        else:
            m, e = evaluar(best["d"], best["w"], best["nw"], th_try, "B")
            vistos_B[clave] = (m, e)
        thB.append(th_try); errB.append(e); mB.append(m)

    th_best, err_best, m_best, info_k2 = _desempatar_theta(
        thB, errB, mB, best["err"], best["th"], best["m"], m_voi)
    best.update(err=err_best, th=th_best, m=m_best)

    # ---------- ETAPA C ----------
    if refinar:
        dens_ref = np.unique(np.clip(best["d"] * np.array([0.92, 1.0, 1.08]),
                                     RHO_MIN, rho_max))
        nw_ref = np.unique(np.clip(
            np.round([num_waves - 200, num_waves, num_waves + 200]), 400, 1500).astype(int))
        for d in dens_ref:
            for nw in nw_ref:
                m, e = evaluar(d, best["w"], nw, best["th"], "C")
                if e < best["err"]:
                    best.update(err=e, d=d, nw=nw, m=m)

    # ---------- ETAPA D: desempate entre finalistas ----------
    traza_mec = []
    m_voi_mec = {}
    finalistas = []
    if etapa_d:
        if progreso:
            progreso(contador["i"], n_total, "D")
        m_voi = dict(m_voi)
        if usa_mec:
            m_voi_mec = metricas_mecanicas(VOI, spacing, res_mec, E_s, nu_s)
            m_voi.update({k: v for k, v in m_voi_mec.items()
                          if k in CLAVES_MECANICAS})
        if caro:
            m_voi.update(metricas_forma(VOI, spacing, **caro))

        finalistas = finalistas_de_traza(n_finalistas)
        mejor_mec = None
        for k, t in enumerate(finalistas):
            if progreso:
                progreso(contador["i"], n_total,
                         f"D ({k+1}/{len(finalistas)})")
            # Regenerar es EXACTO, no aproximado: la semilla es la misma en
            # toda la busqueda, asi que los mismos parametros dan la misma
            # realizacion. Por eso no hace falta haber guardado las mascaras,
            # que a 64^3 y 53 evaluaciones no cabrian en memoria.
            BWc = generar(t["dens"], t["wave_number_pi"], t["nw"],
                          t["thetas"], R_fit, precision)
            mc = medir(BWc, completo=True)
            err_c, n_us = err(m_voi, mc, con_mec=True)
            traza_mec.append({
                "orden_morfo": k, "err_morfo": t["err"], "err_con_mec": float(err_c),
                "n_terminos": int(n_us), "dens": t["dens"],
                "wave_number_pi": t["wave_number_pi"],
                "wave_number_rad": t["wave_number_rad"],
                "nw": t["nw"], "thetas": t["thetas"],
                "Ez_rel": mc.get("Ez_rel"), "Ez_Ex": mc.get("Ez_Ex"),
                "homogeneizo": bool(usa_mec and "Ez_rel" in mc)})
            if mejor_mec is None or err_c < mejor_mec[0]:
                mejor_mec = (err_c, t, mc)

        if mejor_mec is not None:
            _, t, mc = mejor_mec
            best.update(err=mejor_mec[0], d=t["dens"], w=t["wave_number_pi"],
                        nw=t["nw"], th=np.asarray(t["thetas"], float), m=mc)

    # ---------- seleccion robusta (U1) ----------
    info_sel = {"criterio": seleccion, "k": float(k_robusto), "K": int(K)}
    if seleccion == "robusta":
        base_fin = finalistas or finalistas_de_traza(n_finalistas)
        idx, tabla = incertidumbre.seleccion_robusta(
            base_fin,
            lambda t, s: generar(t["dens"], t["wave_number_pi"], t["nw"],
                                 t["thetas"], R_fit, precision, s),
            lambda BW: medir(BW, completo=etapa_d),
            lambda a, b: err(a, b, con_mec=etapa_d),
            m_voi, seed, K, k_robusto, progreso=pr_extra)
        for f, t in zip(tabla, base_fin):
            f.update(dens=t["dens"], wave_number_pi=t["wave_number_pi"],
                     nw=t["nw"], thetas=t["thetas"])
        t = base_fin[idx]
        antes = (float(best["d"]), float(best["w"]), int(best["nw"]),
                 tuple(float(x) for x in np.ravel(best["th"])))
        despues = (float(t["dens"]), float(t["wave_number_pi"]), int(t["nw"]),
                   tuple(float(x) for x in t["thetas"]))
        info_sel.update(tabla=tabla, elegido=idx, cambio=antes != despues)
        best.update(err=tabla[idx]["error"]["media"], d=t["dens"],
                    w=t["wave_number_pi"], nw=t["nw"],
                    th=np.asarray(t["thetas"], float),
                    m=None if despues != antes else best["m"])

    # ---------- orientacion: L2 (impuesta) o L1 (eje) ----------
    # Va aqui, con los parametros ya elegidos, porque ningun termino del error
    # depende de la orientacion: corregirla antes no cambiaria la seleccion y
    # costaria varias generaciones por preset en vez de unas pocas en total.
    R_final = R_fit
    info_alin = {"aplicada": False, "metodo": "clasica (desactivada)",
                 "motivo": "alinear_fabrica=False"}
    if alinear_fabrica:
        if alineacion == "eje":
            R_final, info_alin = alinear_por_fabrica(
                m_voi, res, sp_cand, best["d"], best["w"], best["nw"],
                best["th"], esquema=esquema, seed=seed, precision=precision)
        else:
            R_final, info_alin = alinear_marco_fabrica(
                m_voi,
                lambda R: generar(best["d"], best["w"], best["nw"], best["th"],
                                  R, precision),
                lambda BW: morfometria(BW, sp_cand), R_fit,
                m_clasica=best["m"])

    # ---------- remedida final en doble precision ----------
    BW = generar(best["d"], best["w"], best["nw"], best["th"], R_final, "f64")
    m_final = medir(BW, completo=etapa_d)
    err_final, _ = err(m_voi, m_final, con_mec=True)

    # ---------- incertidumbre del ganador (U1) ----------
    info_inc = {"K": 0}
    if K > 0:
        info_inc = incertidumbre.informe(
            lambda s: generar(best["d"], best["w"], best["nw"], best["th"],
                              R_final, precision, s),
            lambda BWr: medir(BWr, completo=etapa_d),
            lambda a, b: err(a, b, con_mec=True),
            m_voi, seed, K, claves_incertidumbre(pesos, usa_mec),
            progreso=pr_extra)
        # Dos referencias frente a la media de las replicas, que NO son lo
        # mismo: el error con el que la busqueda eligio (optimista, maldicion
        # del ganador) y el de la realizacion que se devuelve, con la R final
        # y en f64 (una tirada mas; en H4 proximal dio 0.034 frente a una
        # media de 0.018 +- 0.007).
        info_inc["error_busqueda"] = float(best["err"])
        info_inc["error_realizacion_devuelta"] = float(err_final)
        info_inc["precision"] = precision

    # Diagnostico: ¿el ganador sigue siendo bicontinuo, y toco el tope?
    fp_spin, n_poros = fraccion_poro_conexa(BW)
    fp_voi, n_poros_voi = fraccion_poro_conexa(VOI)
    en_tope = bool(abs(best["d"] - min(rho_max, bvtv * 1.3)) < 1e-9)
    th_final = [float(t) for t in np.ravel(best["th"])]

    parametros = {
        "familia": "spinodoide",
        "densidad": float(best["d"]),
        # Las dos lecturas, siempre juntas (ver `procedencia.py`).
        "wave_number_pi": float(best["w"]),
        "wave_number_rad": float(best["w"]) * np.pi,
        "num_waves": int(best["nw"]),
        "thetas": th_final,
        "thetas_degenerados": bool(region_degenerada(
            np.clip(np.round(th_final), 0, 90), esquema)["degenerado"]),
        # R es la que SE USO para generar la mascara devuelta. Con la
        # correccion de orientacion activa no coincide con los angulos de
        # Euler clasicos, que se conservan aparte para poder reproducir el
        # comportamiento anterior y para colocar los deslizadores.
        "R": R_final.tolist(),
        "R_clasica": R_fit.tolist(),
        "euler_deg": [rx, ry, rz],
        "resolucion": res,
        "spacing_mm": float(sp_cand[0]),
        "esquema": esquema,
        "semilla": int(seed),
        "rho_max": float(rho_max),
    }

    return {
        "familia": "spinodoide",
        "procedencia": procedencia.desde_parametros(parametros),
        "parametros": parametros,
        "objetivo": {
            "pesos": pesos, "distancia": distancia,
            "peso_mecanico": float(peso_mecanico),
            "medidas_por_evaluacion": sorted(forma),
            "medidas_finalistas": sorted(caro) + (["mecanica"] if usa_mec
                                                  else []),
            "etapa_d": etapa_d,
        },
        "diagnostico": {
            "en_tope_densidad": en_tope,
            "poro_conexo_spin": fp_spin,
            "n_poros_spin": n_poros,
            "poro_conexo_voi": fp_voi,
            "n_poros_voi": n_poros_voi,
            "bicontinuo": bool(fp_spin >= 0.9),
            "aviso": ("El candidato NO es bicontinuo: la fase poro esta "
                      "fragmentada, asi que la estructura es un solido con "
                      "poros aislados y no una trabecula."
                      if fp_spin < 0.9 else ""),
        },
        "alineacion": info_alin,
        "error": float(err_final),
        "error_busqueda": float(best["err"]),
        "incertidumbre": info_inc,
        "seleccion": info_sel,
        "metricas_voi": m_voi,
        "metricas_spin": m_final,
        "mascara": BW,
        # Las rejillas son DETERMINISTAS: salen de formulas sobre BV/TV y sobre
        # la direccion principal del VOI, sin intervencion del azar. Son por
        # tanto lo unico que puede contrastarse EXACTAMENTE contra MATLAB, ya
        # que las dos plataformas evaluan realizaciones distintas en cada punto.
        "rejillas": {
            "densidad": [float(x) for x in dens],
            "onda": [float(x) for x in waves],
            "presets_theta": [[float(t) for t in np.ravel(p)] for p in presets],
            "n_eval_previstas": int(n_total),
        },
        "k2": info_k2,
        # Vacio si no hubo etapa D. Cuando la hubo, `orden_morfo` frente al
        # orden por `err_con_mec` dice de un vistazo si el desempate cambio
        # algo: si el ganador sigue siendo el 0, la etapa D confirmo la
        # morfometria y no decidio nada.
        "mecanico": {
            "peso": float(peso_mecanico),
            "res_mec": int(res_mec),
            "voi": m_voi_mec,
            "traza": traza_mec,
            "reordeno": bool(traza_mec and min(
                range(len(traza_mec)),
                key=lambda i: traza_mec[i]["err_con_mec"]) != 0),
        },
        "traza": traza,
        "n_evaluaciones": contador["i"],
        "tiempo_s": time.time() - t0,
    }
