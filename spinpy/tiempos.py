"""
tiempos.py — Cuanto va a tardar el informe automatico, etapa por etapa.

POR QUE UN MODELO Y NO UNA TABLA
--------------------------------
El coste de cada etapa depende del tamano de lo que procesa, y no de forma
lineal: generar un spinodoide crece como n^2.8, homogeneizar como n^4.1. Una
tabla fija de minutos por etapa se equivocaria por un orden de magnitud en
cuanto cambiara la resolucion. Cada etapa se modela con la ley de potencias de
su calculo dominante, ajustada a tiempos MEDIDOS.

DE DONDE SALEN LOS COEFICIENTES
-------------------------------
Medidos el 2026-09-16 en el equipo de desarrollo (VOI_V1 de talo porcino,
188^3 a 16 um, y su spinodoide ajustado), a dos o tres tamanos por calculo:

  generar spinodoide     48^3 1.25 s    96^3 8.75 s                 ~ n^2.81
  generar dual-lattice   48^3 0.60 s    96^3 1.30 s                 ~ n^1.11
  morfometria basica     64^3 1.17 s    96^3 1.90 s   188^3 4.46 s  ~ vox^0.43
  + Conn.D, SMI          96^3 2.34 s   188^3 8.45 s                  ~ vox^0.64
  + Po.Dm (con los dos)  96^3 ~18 s    188^3 58.6 s                  ~ vox^0.6
  EF                     32^3 184 s  (y 373 s a 48^3, ~30 min a 97^3) ~ vox^0.7
  homogeneizacion        16^3 6.4 s     24^3 31.6 s   32^3 ~112 s    ~ n^4.13
  ensayo, un eje         24^3 2.6 s     32^3 8.0 s    48^3 26.6 s    ~ n^3.34
  analisis comparado     0.87 x ensayo de un eje
  perdida osea           2 pasos a 24^3: 18.3 s (5.2 s de geometria por paso)
  fallo progresivo       2 pasos a 24^3: 9.0 s
  ajuste (3 metodos)     96^3 1359 s, 32^3 131 s  = 127.6 evaluaciones
                         equivalentes de generar + medir
  ajuste dual (2)        32^3 71 s  = ~80 evaluaciones equivalentes
  figuras 3D             contorno + Taubin 10.3 s a 188^3; ~2 s por render
                         de galeria y ~4 s por principal
  distribuciones (fig 7) espesor local de hueso + poro: 188^3 73 s,
                         96^3 4.5 s                          ~ vox^1.38

Figura 0 (metodo), medida el 2026-09-23 sobre H4 proximal (VOI 97^3,
candidatos 96^3), dos pasadas por caso, con ruido de +-30 % entre pasadas:

  renders     1 familia 8.8-11.5 s, 2 familias 14.2-18.2 s
              (VOI ~3 s, pila reducida ~1 s, ~6 s por familia)
  componer    ES + EN a 600 ppp: 2 filas 21.0-23.6 s, 3 filas 29.3-49.4 s
              (~12 s por fila); 1 familia + pila 44.2-50.0 s en total
Los campos y las mascaras ya los regenera el paquete de reproduccion: la
figura no vuelve a llamar al generador.

El ajuste del spinodoide y el tensor se calibraron con la ejecucion completa
del informe automatico sobre ese VOI, asi que coincidir ahi no prueba nada.
Lo que NO se uso para ajustar: el ensayo de esa ejecucion (modelo 35 s, real
28 s), el ajuste dual de un VOI sintetico de 32^3 (68 s frente a 71 s; su
coeficiente sale solo de contar evaluaciones) y el total de la ejecucion del
cerdo (1773 s frente a 1755 s). En etapas de pocos segundos el coste fijo
domina y el modelo se queda corto en unos segundos, que no importan.

Y EN OTRO EQUIPO
----------------
Los coeficientes son de UN equipo. `factores` corrige cada etapa con lo que
tardo de verdad en ejecuciones anteriores del mismo equipo (el visor guarda
el cociente real/estimado en QSettings al terminar cada etapa). La primera
estimacion en un equipo nuevo puede errar por un factor 2; se declara asi.
"""

from __future__ import annotations

VOX_VOI = 6644672          # 188^3, el VOI de la calibracion


def gen(familia, n):
    n = float(n)
    if familia == "dual-lattice":
        return 1.30 * (n / 96.0) ** 1.11
    return 8.75 * (n / 96.0) ** 2.81


def morfo(vox, extra=False, poro=False, ef=False):
    r = float(vox) / VOX_VOI
    if poro:
        t = 58.6 * r ** 0.6
    elif extra:
        t = 8.45 * r ** 0.64
    else:
        t = 4.46 * r ** 0.43
    if ef:
        t += 184.0 * (float(vox) / 32768.0) ** 0.7
    return t


def homog(n):
    return 112.0 * (float(n) / 32.0) ** 4.13


def ensayo(n):
    return 8.0 * (float(n) / 32.0) ** 3.34


def comparado(n):
    return 0.87 * ensayo(n) + 0.5


# Evaluaciones equivalentes de un ajuste completo por familia: los metodos
# comparables mas las replicas de incertidumbre.
EVAL_AJUSTE = {"spinodoide": 127.6, "dual-lattice": 80.0}


def ajuste(familia, n_fit, vox_voi):
    return (EVAL_AJUSTE.get(familia, 127.6)
            * (gen(familia, n_fit) + morfo(int(n_fit) ** 3))
            + morfo(vox_voi))


def resoluciones_convergencia(n_fe):
    """Las mismas que `Visor.convergencia_fe`: menos de tres, no se lanza."""
    return sorted({n for n in (int(n_fe * f) for f in (0.55, 0.7, 0.85, 1.0))
                   if n >= 12})


def distribuciones(vox):
    return 73.0 * (float(vox) / VOX_VOI) ** 1.38


def metodo(n_familias, con_pila=False):
    """Figura 0: renders (VOI, pila y cuatro por familia) y dos composiciones."""
    filas = 1 + int(n_familias)
    # La fila de la pila cuesta mas al componer (rebanada en gris a 600 ppp):
    # 1 s de render y ~8 s de composicion por encima de una fila normal.
    return (4.0 + 6.0 * n_familias + (9.0 if con_pila else 0.0)
            + 12.0 * filas)


def informe(familias, n_fit, vox_voi, n_estilos, n_vistas, suavizar=True,
            con_distribuciones=False, con_von_mises=False, con_metodo=False,
            con_pila=False):
    t = 0.6 + sum(gen(f, n_fit) for f in familias)         # regenerar
    t += 12.0                                                # figuras de datos
    estructuras = [vox_voi] + [int(n_fit) ** 3] * len(familias)
    for vox in estructuras:
        malla = 10.3 * (float(vox) / VOX_VOI) ** 0.9
        t += (malla if suavizar else 0.4 * malla)
        t += 4.0 + 2.0 * n_estilos * n_vistas
        if con_distribuciones:
            t += distribuciones(vox)
    if con_von_mises:
        t += 3.0 * len(estructuras) + 2.0
    if con_metodo:
        t += metodo(len(familias), con_pila)
    return t + 12.0                                          # PDF


def estimar(plan, factores=None):
    """Segundos por etapa: {clave: {familia o None: s}}.

    `plan` describe lo que se va a correr (ver `DialogoInformeAuto.plan`);
    `factores` es {clave: real/estimado} de ejecuciones anteriores.
    """
    f = dict(factores or {})
    fams = list(plan["familias"])
    vox_voi = int(plan["vox_voi"])
    n_fit = int(plan["n_fit"])
    out = {}

    def poner(clave, fam, s):
        out.setdefault(clave, {})[fam] = float(s) * float(f.get(clave, 1.0))

    if plan.get("ajuste"):
        for fam in fams:
            poner("ajuste", fam, ajuste(fam, n_fit, vox_voi))
    if plan.get("morfometria"):
        kw = dict(extra=plan.get("extra"), poro=plan.get("poro"),
                  ef=plan.get("ef"))
        s = sum(morfo(int(plan["n_cand"]) ** 3, **kw) for _x in fams)
        if plan.get("medir_voi", True):
            s += morfo(vox_voi, **kw)
        poner("morfometria", None, s)
    if plan.get("dispersion"):
        n = int(plan["n_resm"])
        for fam in fams:
            poner("dispersion", fam, int(plan["K"]) * (
                gen(fam, n) + morfo(n ** 3, extra=plan.get("extra"))))
    if plan.get("elastico"):
        for fam in fams:
            poner("elastico", fam, 2 * homog(plan["n_homog"]))
    if plan.get("ensayo"):
        for fam in fams:
            poner("ensayo", fam,
                  2 * int(plan["n_ejes"]) * ensayo(plan["n_fe"]) + 1.0)
    if plan.get("convergencia"):
        rs = resoluciones_convergencia(int(plan["n_fe"]))
        poner("convergencia", None,
              sum(ensayo(r) for r in rs) if len(rs) >= 3 else 0.0)
    if plan.get("comparado"):
        for fam in fams:
            poner("comparado", fam, 2 * comparado(plan["n_fe"]))
    n_mec = int(plan.get("n_fe", 32))
    n_trab = min(96, int(plan.get("n_sim", 96)))
    if plan.get("perdida"):
        pasos = int(plan["pasos"]) * (2 if plan.get("protocolo")
                                      == "recuperacion" else 1)
        poner("perdida", fams[0], (pasos + 1) * ensayo(n_mec)
              + pasos * 5.2 * (n_trab / 96.0) ** 3)
    if plan.get("fallo"):
        poner("fallo", fams[0], (int(plan["pasos"]) + 1) * ensayo(n_mec) * 1.15)
    poner("informe", None, informe(fams, n_fit, vox_voi,
                                   int(plan.get("n_estilos", 2)),
                                   int(plan.get("n_vistas", 5)),
                                   plan.get("suavizar", True),
                                   plan.get("distribuciones", False),
                                   plan.get("von_mises", False),
                                   plan.get("metodo", False),
                                   plan.get("pila", False)))
    return out


def total(est):
    return sum(sum(v.values()) for v in est.values())


def texto(s, idioma="es"):
    """Duracion legible y redondeada a lo que el modelo puede prometer."""
    s = float(s)
    if s < 45:
        return "< 1 min"
    m = s / 60.0
    if m < 10:
        return f"~{m:.0f} min"
    if m < 60:
        return f"~{5 * round(m / 5):.0f} min"
    h = int(m // 60)
    r = 10 * round((m - 60 * h) / 10)
    if r == 60:
        h, r = h + 1, 0
    return f"~{h} h {r:02d} min" if r else f"~{h} h"


def actualizar_factor(anterior, real, estimado, peso=0.5):
    """Media movil del cociente real/estimado, acotada a [0.2, 5]: una etapa
    rara (el equipo ocupado con otra cosa) no descalibra todo."""
    if not estimado or estimado <= 0 or real is None or real <= 0:
        return anterior
    c = min(5.0, max(0.2, float(real) / float(estimado)))
    if anterior is None:
        return c
    return (1 - peso) * float(anterior) + peso * c
