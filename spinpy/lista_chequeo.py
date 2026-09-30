"""
lista_chequeo.py: seccion «Lista de chequeo de reporte» del informe para
publicacion.

POR QUE EXISTE
--------------
El informe ya dice que se hizo y si cada numero se puede citar. Lo que un
revisor contrasta ademas es si el metodo esta descrito con el detalle que
piden las guias de reporte del campo. Aqui se hace ese contraste item por
item, con los datos de la sesion, y se dice que falta y donde ponerlo.

LAS GUIAS
---------
  Bouxsein et al. 2010 (J Bone Miner Res 25:1468; PMID 20533309). Guia de
      morfometria osea por micro-TC. Esta escrita para hueso de ROEDOR ex vivo;
      aqui se aplica por extension a otras especies y el informe lo declara.
  Erdemir et al. 2012 (J Biomech 45:625; PMID 22236526). Parametros de
      reporte de estudios de elementos finitos en biomecanica, organizados en
      identificacion del modelo, estructura del modelo, estructura de la
      simulacion, verificacion, validacion, disponibilidad y multiescala. Los
      items de abajo siguen esos dominios, contrastados con el texto completo
      (PMC3278509).
  Oefner et al. 2021 (Med Eng Phys 92:25; PMID 34167708) publicaron despues
      una lista de verificacion y validacion para biomecanica ortopedica y
      traumatologica. Se CITA como lectura recomendada, pero sus items no se
      evaluan: no se contrastaron con su texto completo.

QUE SIGNIFICA CADA ESTADO
-------------------------
El estado se refiere a lo que ESTE INFORME ya documenta, no al manuscrito. Lo
que spinpy no puede saber (el protocolo de adquisicion, la especie) sale como
«no cumple» con la accion concreta que corresponde: no se presume cumplido
algo que el programa no ha visto.

  cumple      el informe lo documenta; se dice en que seccion
  parcial     el informe documenta una parte; la accion dice el resto
  no_cumple   la sesion no lo contiene; la accion dice que aportar
  no_aplica   el item no corresponde a esta sesion

Como `informe.comprobar`, la evaluacion devuelve datos (estado y textos ya en
su idioma); `seccion` los compone en Markdown.
"""

from __future__ import annotations

import numpy as np

CUMPLE, PARCIAL, NO_CUMPLE, NO_APLICA = ("cumple", "parcial", "no_cumple",
                                         "no_aplica")
ORDEN = (CUMPLE, PARCIAL, NO_CUMPLE, NO_APLICA)
ESTADOS = {CUMPLE: ("Cumple", "Met"),
           PARCIAL: ("Cumple parcialmente", "Partially met"),
           NO_CUMPLE: ("No cumple", "Not met"),
           NO_APLICA: ("No aplica", "Not applicable")}

TITULOS = {"bouxsein2010": ("Bouxsein et al. 2010: morfometría por micro-TC",
                            "Bouxsein et al. 2010: micro-CT morphometry"),
           "erdemir2012": ("Erdemir et al. 2012: elementos finitos",
                           "Erdemir et al. 2012: finite elements")}

# Dominios de Erdemir et al. 2012, con sus nombres en el articulo.
DOMINIOS = {
    "identificacion": ("identificación del modelo", "model identification"),
    "estructura": ("estructura del modelo", "model structure"),
    "simulacion": ("estructura de la simulación", "simulation structure"),
    "verificacion": ("verificación", "verification"),
    "validacion": ("validación", "validation"),
    "disponibilidad": ("disponibilidad", "availability"),
    "multiescala": ("multiescala", "multiscale"),
}
DOMINIO_BOUXSEIN = {
    "adquisicion": ("adquisición", "acquisition"),
    "procesado": ("procesado de imagen", "image processing"),
    "analisis": ("análisis", "analysis"),
    "resultados": ("resultados", "outcomes"),
}

# Conjunto minimo de variables trabeculares de Bouxsein et al. 2010.
MINIMO_TRABECULAR = ("BVTV", "TbN", "TbTh", "TbSp")


def _fila(guia, dominio, item, estado, donde, accion=""):
    return {"guia": guia, "dominio": dominio, "item": item, "estado": estado,
            "donde": donde, "accion": accion}


def _entero(x, idioma):
    """Entero con separador de miles: espacio en espanol, coma en ingles."""
    s = f"{int(x):,}"
    return s.replace(",", " ") if idioma == "es" else s


def _recs(res, base):
    """Registros de las dos familias (`base` y `base_dual`) que existan."""
    return [x for x in (res.get(base), res.get(base + "_dual"))
            if isinstance(x, dict)]


def _ejes_voi(rec):
    """Resultados por eje del VOI en un registro de resistencia."""
    pe = (rec.get("por_estructura") or {}).get("voi") or {}
    return [v for v in (pe.get("ejes") or {}).values() if isinstance(v, dict)]


def _residuos(res):
    """Todos los residuos relativos que la sesion guardo."""
    from .informe import _f
    out = []
    for rec in _recs(res, "resistencia") + _recs(res, "analisis_comparado"):
        for v in (rec.get("por_estructura") or {}).values():
            for e in ((v or {}).get("ejes") or {}).values():
                x = _f(e, "residuo_rel")
                if x is not None:
                    out.append(x)
    for rec in _recs(res, "elastico"):
        for v in (rec.get("por_estructura") or {}).values():
            x = _f(v, "residuo_rel")
            if x is not None:
                out.append(x)
    conv = res.get("convergencia")
    if isinstance(conv, dict):
        out += [x for x in (_f(p, "residuo") for p in conv.get("puntos") or []
                            if p.get("ok")) if x is not None]
    return out


def evaluar(doc, idioma="es", r=None, secciones=None):
    """Lista de filas de la lista de chequeo, con los textos en `idioma`.

    `secciones`: {"metodos", "comparacion", "citabilidad", "modelos",
    "cita", "reproduccion"} -> numero de seccion del informe, para decir
    donde esta cada cosa. "modelos" puede ser None si no hay esa seccion.
    """
    from . import informe as I          # import tardio: informe importa esto
    r = r or I._Redactor(idioma)
    t = r.t
    S = {"metodos": 1, "comparacion": 2, "citabilidad": 3, "modelos": None,
         "cita": None, "reproduccion": None}
    S.update(secciones or {})

    def sec(*claves):
        nums = [str(S[k]) for k in claves if S.get(k)]
        return "§" + ", §".join(nums) if nums else ""

    res = doc.get("resultados") or {}
    voi = doc.get("voi") or {}
    m_voi = doc.get("morfometria_voi") or {}
    hay_cand = bool(doc.get("morfometria_spin") or doc.get("morfometria_dual"))
    forma, sp = voi.get("forma"), voi.get("spacing_mm")
    um = None
    try:
        um = float(np.ravel(sp)[0]) * 1000.0
    except (TypeError, IndexError, ValueError):
        pass
    lado = I._lado_voi(doc)
    sens = (res.get("sensibilidad_superficie") or {}).get(
        "por_estructura") or {}
    sens_voi = "voi" in sens
    dims = "×".join(str(int(x)) for x in forma) if forma else None

    F = []
    B, E = "bouxsein2010", "erdemir2012"
    D = DOMINIO_BOUXSEIN

    # ------------------------------------------------------------------
    # Bouxsein et al. 2010: morfometria por micro-TC
    # ------------------------------------------------------------------
    if not (m_voi or forma):
        F.append(_fila(B, D["resultados"], t("Toda la guía", "Whole guideline"),
                       NO_APLICA, t("La sesión no contiene un VOI de imagen.",
                                    "The session holds no image VOI.")))
    else:
        F.append(_fila(
            B, D["adquisicion"],
            t("Equipo, medio de escaneo, potencial y corriente del tubo, "
              "tiempo de integración, filtro y reconstrucción",
              "Scanner, scanning medium, tube potential and current, "
              "integration time, filter and reconstruction"),
            NO_CUMPLE,
            t("No consta en la sesión: spinpy recibe la imagen ya "
              "reconstruida.",
              "Not recorded in the session: spinpy receives the image "
              "already reconstructed."),
            t("Declararlos en el manuscrito. En equipos SkyScan figuran en el "
              "archivo .log de la reconstrucción.",
              "State them in the manuscript. On SkyScan systems they are in "
              "the reconstruction .log file.")))
        if um is not None:
            F.append(_fila(B, D["adquisicion"], t("Tamaño de vóxel",
                                                  "Voxel size"), CUMPLE,
                           t("{s}: {u} µm.", "{s}: {u} µm.",
                             s=sec("metodos"), u=r.n(um, ".4g"))))
        else:
            F.append(_fila(B, D["adquisicion"], t("Tamaño de vóxel",
                                                  "Voxel size"), NO_CUMPLE,
                           t("La sesión no registra el tamaño de vóxel.",
                             "The session does not record the voxel size."),
                           t("Declarar el tamaño de vóxel en µm.",
                             "State the voxel size in µm.")))
        rec_pila = voi.get("recorte") if isinstance(voi.get("recorte"),
                                                    dict) else None
        if dims:
            donde = t("{s}: {d} vóxeles", "{s}: {d} voxels", s=sec("metodos"),
                      d=dims)
            if lado:
                donde += t(", {l} mm de lado", ", {l} mm side",
                           l=r.n(lado, ".4g"))
            if rec_pila:
                donde += t("; recortado de la pila en el marco de sus ejes "
                           "principales", "; cropped from the stack in the "
                           "frame of its principal axes")
            F.append(_fila(
                B, D["procesado"], t("Tamaño y localización del VOI",
                                     "Size and location of the VOI"),
                PARCIAL, donde + ".",
                t("Describir la localización anatómica del VOI y el criterio "
                  "con que se delimitó.",
                  "Describe the anatomical location of the VOI and the "
                  "criterion used to delimit it.")))
        F.append(_fila(
            B, D["procesado"], t("Filtrado de la imagen",
                                 "Image filtering"), NO_CUMPLE,
            t("No consta en la sesión.", "Not recorded in the session."),
            t("Declarar el filtrado aplicado antes de segmentar (tipo y "
              "parámetros) o que no se aplicó ninguno.",
              "State the filtering applied before segmentation (type and "
              "parameters) or that none was applied.")))
        umbral = I._f(rec_pila, "umbral") if rec_pila else None
        accion_seg = t(
            "Declarar el método de segmentación y justificar el umbral. La "
            "segmentación es la mayor fuente de incertidumbre del proceso: "
            "conviene informar la sensibilidad de las métricas al umbral.",
            "State the segmentation method and justify the threshold. "
            "Segmentation is the largest source of uncertainty in the "
            "pipeline: report the sensitivity of the metrics to the "
            "threshold.")
        if umbral is not None:
            F.append(_fila(B, D["procesado"], t("Segmentación",
                                                "Segmentation"), PARCIAL,
                           t("Umbral global de {u} (unidades de la imagen), "
                             "registrado al cargar la pila.",
                             "Global threshold of {u} (image units), "
                             "recorded when the stack was loaded.",
                             u=r.n(umbral, ".6g")), accion_seg))
        elif sens_voi:
            F.append(_fila(B, D["procesado"], t("Segmentación",
                                                "Segmentation"), PARCIAL,
                           t("El VOI llegó ya binarizado; {s} cuantifica la "
                             "sensibilidad de las métricas a desplazar su "
                             "superficie, como sustituto del umbral.",
                             "The VOI arrived already binarised; {s} "
                             "quantifies the sensitivity of the metrics to "
                             "moving its surface, as a threshold surrogate.",
                             s=sec("metodos")),
                           t("Declarar el método de segmentación y el umbral "
                             "usados.",
                             "State the segmentation method and threshold "
                             "used.")))
        else:
            F.append(_fila(B, D["procesado"], t("Segmentación",
                                                "Segmentation"), NO_CUMPLE,
                           t("El VOI llegó ya binarizado y la sesión no "
                             "registra cómo.",
                             "The VOI arrived already binarised and the "
                             "session does not record how."), accion_seg))
        if not m_voi:
            F.append(_fila(B, D["analisis"],
                           t("Algoritmos de morfometría",
                             "Morphometry algorithms"), NO_CUMPLE,
                           t("No se midió la morfometría del VOI.",
                             "The VOI morphometry was not measured."),
                           t("Medir la morfometría del VOI.",
                             "Measure the VOI morphometry.")))
        elif I._f(m_voi, "TbTh_local") is not None:
            F.append(_fila(B, D["analisis"],
                           t("Algoritmos 3D directos", "Direct 3D algorithms"),
                           CUMPLE,
                           t("{s}: BV/TV por conteo de vóxeles, BS por "
                             "marching cubes y espesor local por esferas "
                             "inscritas.",
                             "{s}: BV/TV by voxel counting, BS by marching "
                             "cubes and local thickness by inscribed "
                             "spheres.", s=sec("metodos"))))
        else:
            F.append(_fila(
                B, D["analisis"], t("Algoritmos 3D directos",
                                    "Direct 3D algorithms"), PARCIAL,
                t("{s}: Tb.Th, Tb.Sp y Tb.N se derivan del modelo de placas "
                  "(Tb.Th = 2·BV/BS).",
                  "{s}: Tb.Th, Tb.Sp and Tb.N are derived from the plate "
                  "model (Tb.Th = 2·BV/BS).", s=sec("metodos")),
                t("La guía recomienda algoritmos 3D que no supongan la forma "
                  "de la estructura: informar también el espesor local por "
                  "esferas inscritas, que spinpy calcula, o justificar el "
                  "modelo de placas.",
                  "The guideline recommends 3D algorithms that assume no "
                  "structure shape: also report the local thickness by "
                  "inscribed spheres, which spinpy computes, or justify the "
                  "plate model.")))
        if m_voi:
            faltan = [k for k in MINIMO_TRABECULAR if I._f(m_voi, k) is None]
            item = t("Conjunto mínimo: BV/TV, Tb.N, Tb.Th, Tb.Sp",
                     "Minimal set: BV/TV, Tb.N, Tb.Th, Tb.Sp")
            if faltan:
                F.append(_fila(B, D["resultados"], item, NO_CUMPLE,
                               t("Faltan: {k}.", "Missing: {k}.",
                                 k=", ".join(faltan)),
                               t("Medir y tabular las cuatro variables.",
                                 "Measure and tabulate the four variables.")))
            elif hay_cand:
                F.append(_fila(B, D["resultados"], item, CUMPLE,
                               t("{s}: tabla del VOI con unidades.",
                                 "{s}: VOI table with units.",
                                 s=sec("comparacion"))))
            else:
                F.append(_fila(
                    B, D["resultados"], item, PARCIAL,
                    t("Medidas, pero el informe solo las tabula al "
                      "compararlas con un candidato ajustado.",
                      "Measured, but the report only tabulates them when "
                      "comparing with a fitted candidate."),
                    t("Tabular en el manuscrito las cuatro variables del VOI "
                      "con sus unidades.",
                      "Tabulate the four VOI variables with their units in "
                      "the manuscript.")))
            F.append(_fila(B, D["resultados"], t("Nomenclatura y unidades",
                                                 "Nomenclature and units"),
                           CUMPLE,
                           t("{s}: nomenclatura de Parfitt et al. y Bouxsein "
                             "et al., longitudes en mm.",
                             "{s}: nomenclature of Parfitt et al. and "
                             "Bouxsein et al., lengths in mm.",
                             s=sec("metodos", "comparacion"))))

    # ------------------------------------------------------------------
    # Erdemir et al. 2012: elementos finitos
    # ------------------------------------------------------------------
    D = DOMINIOS
    ela = _recs(res, "elastico")
    fe = _recs(res, "resistencia")
    comp = _recs(res, "analisis_comparado")
    conv = res.get("convergencia") if isinstance(res.get("convergencia"),
                                                 dict) else None
    febio = (res.get("febio") or {}).get("registros") \
        if isinstance(res.get("febio"), dict) else None
    motores_fem = ""
    if febio:
        from .informe import nombre_motor
        nombres = []
        for reg in febio:
            n = nombre_motor(reg)
            if n not in nombres and n != "?":
                nombres.append(n)
        motores_fem = ", ".join(nombres)
    if not (ela or fe or comp or conv):
        F.append(_fila(E, D["identificacion"], t("Toda la guía",
                                                 "Whole guideline"),
                       NO_APLICA,
                       t("La sesión no contiene análisis por elementos "
                         "finitos.",
                         "The session holds no finite element analysis.")))
        return F

    analisis = []
    if fe:
        analisis.append(t("ensayo de compresión uniaxial",
                          "uniaxial compression test"))
    if comp:
        analisis.append(t("análisis comparado", "compared analysis"))
    if ela:
        analisis.append(t("homogeneización periódica",
                          "periodic homogenisation"))
    if conv and not (fe or comp):
        analisis.append(t("estudio de convergencia",
                          "convergence study"))
    F.append(_fila(
        E, D["identificacion"],
        t("Descripción del modelo y formulación matemática",
          "Model description and mathematical formulation"), CUMPLE,
        t("{s}: elasticidad lineal estática, tejido isótropo; {a}.",
          "{s}: static linear elasticity, isotropic tissue; {a}.",
          s=sec("metodos", "modelos"), a=", ".join(analisis))))
    F.append(_fila(
        E, D["identificacion"],
        t("Estructura anatómica, especie, estado del espécimen y patología",
          "Anatomical structure, species, specimen state and pathology"),
        NO_CUMPLE,
        t("La sesión solo registra el nombre del VOI ({n}).",
          "The session only records the VOI name ({n}).",
          n=voi.get("nombre") or "VOI"),
        t("Declarar especie, hueso y región, condición del espécimen (ex "
          "vivo, conservación) y cualquier patología.",
          "State species, bone and region, specimen condition (ex vivo, "
          "preservation) and any pathology.")))
    if lado and um is not None:
        F.append(_fila(E, D["identificacion"], t("Escala espacial",
                                                 "Spatial scale"), CUMPLE,
                       t("{s}: cubo de {l} mm de lado, vóxel de {u} µm.",
                         "{s}: cube of {l} mm side, {u} µm voxel.",
                         s=sec("metodos"), l=r.n(lado, ".4g"),
                         u=r.n(um, ".4g"))))
    else:
        F.append(_fila(E, D["identificacion"], t("Escala espacial",
                                                 "Spatial scale"), NO_CUMPLE,
                       t("La sesión no registra el tamaño físico del VOI.",
                         "The session does not record the physical size of "
                         "the VOI."),
                       t("Declarar el lado del VOI y el tamaño de vóxel.",
                         "State the VOI side and the voxel size.")))
    F.append(_fila(
        E, D["identificacion"], t("Limitaciones", "Limitations"), PARCIAL,
        t("{s}: cada resultado va marcado como citable, con reservas o no "
          "citable, con su motivo.",
          "{s}: every result is marked citable, with caveats or not "
          "citable, with its reason.", s=sec("citabilidad")),
        t("Discutir en el manuscrito la limitación principal y su efecto "
          "sobre las conclusiones.",
          "Discuss the main limitation in the manuscript and its effect on "
          "the conclusions.")))

    cc = []
    for rec in fe:
        cc.append(t("apoyo «{a}» en el ensayo de compresión",
                    "support \"{a}\" in the compression test",
                    a=I.nombre_apoyo(rec.get("apoyo"), idioma)))
    for rec in comp:
        cc.append(t("{F} N con apoyo «{a}» en el análisis comparado",
                    "{F} N with support \"{a}\" in the compared analysis",
                    F=r.n(rec.get("carga_N"), ".4g"),
                    a=I.nombre_apoyo(rec.get("apoyo"), idioma)))
    if ela:
        cc.append(t("condiciones periódicas en la homogeneización",
                    "periodic conditions in the homogenisation"))
    if isinstance(res.get("condiciones_contorno"), dict):
        cc.append(t("comparación del módulo con las tres condiciones sobre "
                    "la misma rejilla",
                    "comparison of the modulus under the three conditions "
                    "on the same grid"))
    F.append(_fila(E, D["estructura"],
                   t("Condiciones de contorno y carga",
                     "Boundary and loading conditions"), CUMPLE,
                   t("{s}: {c}.", "{s}: {c}.", s=sec("metodos", "modelos"),
                     c="; ".join(dict.fromkeys(cc)))))
    salidas = []
    if fe or comp or conv:
        salidas.append(t("módulo aparente E_app", "apparent modulus E_app"))
    if fe or comp:
        salidas.append(t("percentil 99 de von Mises en la capa superficial",
                         "99th percentile of von Mises on the surface layer"))
    if fe:
        salidas.append(t("carga de fallo de Pistoia", "Pistoia failure load"))
    if ela:
        salidas.append(t("tensor de rigidez y constantes de ingeniería",
                         "stiffness tensor and engineering constants"))
    F.append(_fila(E, D["estructura"], t("Variables de salida",
                                         "Output variables"), CUMPLE,
                   t("{s}: {v}.", "{s}: {v}.",
                     s=sec("comparacion", "citabilidad"),
                     v=", ".join(salidas))))
    F.append(_fila(
        E, D["estructura"], t("Origen de la geometría",
                              "Source of the geometry"), PARCIAL,
        t("{s}: imagen binaria{d} remuestreada a la rejilla de elementos.",
          "{s}: binary image{d} resampled onto the element grid.",
          s=sec("metodos"),
          d=t(" de {x} vóxeles", " of {x} voxels", x=dims) if dims else ""),
        t("Declarar modalidad de imagen, protocolo de adquisición y "
          "segmentación (ítems de Bouxsein et al.).",
          "State the imaging modality, acquisition protocol and "
          "segmentation (items of Bouxsein et al.).")))

    rejillas = sorted({n for n in (I._i(x.get("resolucion"))
                                   for x in ela + fe + comp) if n})
    n_elem = sorted({n for n in (I._i(e.get("n_elem"))
                                 for rec in fe for e in _ejes_voi(rec)) if n})
    donde = t("{s}: hexaedros lineales de 8 nodos, uno por celda de la "
              "rejilla", "{s}: 8-node linear hexahedra, one per grid cell",
              s=sec("metodos", "modelos"))
    if rejillas:
        donde += t(", rejillas de {g}", ", grids of {g}",
                   g=", ".join(f"{n}³" for n in rejillas))
        if lado:
            donde += t(" (arista de {h} mm)", " (edge of {h} mm)",
                       h=" / ".join(r.n(lado / n, ".3g") for n in rejillas))
    if n_elem:
        extremos = dict.fromkeys((n_elem[0], n_elem[-1]))
        donde += t("; {e} elementos en el ensayo del VOI",
                   "; {e} elements in the VOI test",
                   e=" – ".join(_entero(x, idioma) for x in extremos))
    F.append(_fila(E, D["estructura"], t("Malla", "Mesh"), CUMPLE,
                   donde + "."))

    mat = ela + fe
    E_s = next((I._f(x, "E_s_Pa") for x in mat
                if I._f(x, "E_s_Pa") is not None), None)
    nu_s = next((I._f(x, "nu_s") for x in mat + comp
                 if I._f(x, "nu_s") is not None), None)
    if E_s is None and comp:
        E_s = I._f(comp[0], "E_s_Pa")
    F.append(_fila(
        E, D["estructura"], t("Propiedades del material",
                              "Material properties"), PARCIAL,
        t("{s}: elástico lineal isótropo, E_s = {E} GPa, ν_s = {nu}.",
          "{s}: isotropic linear elastic, E_s = {E} GPa, ν_s = {nu}.",
          s=sec("metodos"),
          E=r.n(E_s / 1e9 if E_s is not None else None, ".4g"),
          nu=r.n(nu_s, ".3g")),
        t("Justificar E_s y ν_s con una referencia o una medida del tejido. "
          "Con tejido homogéneo el módulo aparente es proporcional a E_s.",
          "Justify E_s and ν_s with a reference or a tissue measurement. "
          "With homogeneous tissue the apparent modulus is proportional to "
          "E_s.")))
    F.append(_fila(E, D["estructura"], t("Interacciones y contacto",
                                         "Interactions and contact"),
                   NO_APLICA,
                   t("Un solo material, sin contacto entre componentes.",
                     "A single material, no contact between components.")))

    proc = doc.get("procedencia") or {}
    from . import __version__
    sw = t("{s}: spinpy V{v}, Python {py}, NumPy {npv}",
           "{s}: spinpy V{v}, Python {py}, NumPy {npv}", s=sec("metodos"),
           v=proc.get("spinpy") or __version__,
           py=proc.get("python") or "?", npv=proc.get("numpy") or "?")
    if febio:
        sw += t("; motores FEM ({m}) en {s}", "; FEM engines ({m}) in {s}",
                m=motores_fem or "?", s=sec("modelos"))
    F.append(_fila(E, D["simulacion"], t("Software y versión",
                                         "Software and version"), CUMPLE,
                   sw + "."))
    solver = ""
    for rec in fe + comp + ela:
        for v in (rec.get("por_estructura") or {}).values():
            v = v or {}
            solver = solver or v.get("solver") or next(
                (e.get("solver") for e in (v.get("ejes") or {}).values()
                 if e.get("solver")), "")
    F.append(_fila(E, D["simulacion"],
                   t("Estrategia de solución y algoritmos",
                     "Solution strategy and algorithms"), CUMPLE,
                   t("{s}: análisis estático lineal resuelto con {sol}.",
                     "{s}: static linear analysis solved with {sol}.",
                     s=sec("metodos"),
                     sol=I._solver_txt(solver, idioma))))
    resid = _residuos(res)
    if resid:
        F.append(_fila(E, D["simulacion"],
                       t("Tolerancias de convergencia",
                         "Convergence tolerances"), CUMPLE,
                       t("{s}: criterio de aceptación del residuo; residuo "
                         "relativo máximo en esta sesión: {x}.",
                         "{s}: residual acceptance criterion; largest "
                         "relative residual in this session: {x}.",
                         s=sec("metodos", "modelos"),
                         x=r.n(max(resid), ".1e"))))
    else:
        F.append(_fila(E, D["simulacion"],
                       t("Tolerancias de convergencia",
                         "Convergence tolerances"), PARCIAL,
                       t("{s}: criterio de aceptación del residuo; la "
                         "sesión no guardó los residuos obtenidos.",
                         "{s}: residual acceptance criterion; the session "
                         "did not store the residuals obtained.",
                         s=sec("modelos") or sec("metodos")),
                       t("Informar el residuo relativo de cada solución.",
                         "Report the relative residual of each solution.")))
    F.append(_fila(
        E, D["simulacion"], t("Posprocesado", "Post-processing"), CUMPLE,
        (t("{s}: percentil 99 de von Mises en la capa superficial en lugar "
           "del máximo, y criterio de Pistoia.",
           "{s}: 99th percentile of von Mises on the surface layer instead "
           "of the maximum, and Pistoia criterion.", s=sec("modelos"))
         if (fe or comp) else
         t("{s}: constantes de ingeniería a partir del tensor.",
           "{s}: engineering constants from the tensor.",
           s=sec("modelos") or sec("metodos")))))

    if febio:
        F.append(_fila(E, D["verificacion"],
                       t("Verificación del código", "Code verification"),
                       CUMPLE,
                       t("{s}: la sesión se resolvió además con otros "
                         "motores FEM ({m}) y se tabulan sus diferencias.",
                         "{s}: the session was also solved with other FEM "
                         "engines ({m}) and their differences are tabulated.",
                         m=motores_fem or "?", s=sec("modelos"))))
    else:
        F.append(_fila(
            E, D["verificacion"], t("Verificación del código",
                                    "Code verification"), PARCIAL,
            t("spinpy está verificado frente a soluciones analíticas, frente "
              "a FEBio 4.5 y entre sus motores FEM internos, pero esa "
              "verificación está en el repositorio, no en este informe.",
              "spinpy is verified against analytical solutions, against "
              "FEBio 4.5 and across its internal FEM engines, but that "
              "verification is in the repository, not in this report."),
            t("Citar en el manuscrito la verificación del software.",
              "Cite the software verification in the manuscript.")))
    puntos = [p for p in (conv or {}).get("puntos") or [] if p.get("ok")]
    item_conv = t("Convergencia de malla", "Mesh convergence")
    accion_conv = t("Ejecutar el estudio de convergencia de malla de la "
                    "aplicación e informarlo con el resultado principal.",
                    "Run the application's mesh convergence study and "
                    "report it with the main result.")
    if len(puntos) >= 3:
        ns = sorted({I._i(p.get("n")) for p in puntos if I._i(p.get("n"))})
        F.append(_fila(E, D["verificacion"], item_conv, CUMPLE,
                       t("E_app en {k} resoluciones ({n}); veredicto en {s}.",
                         "E_app at {k} resolutions ({n}); verdict in {s}.",
                         k=len(puntos), n=", ".join(f"{x}³" for x in ns),
                         s=sec("citabilidad"))))
    elif puntos:
        F.append(_fila(E, D["verificacion"], item_conv, PARCIAL,
                       t("Solo {k} resoluciones válidas.",
                         "Only {k} valid resolutions.", k=len(puntos)),
                       accion_conv))
    else:
        F.append(_fila(E, D["verificacion"], item_conv, NO_CUMPLE,
                       t("La sesión no contiene un estudio de convergencia.",
                         "The session holds no convergence study."),
                       accion_conv))
    F.append(_fila(
        E, D["verificacion"], t("Reproducción en otra plataforma",
                                "Reproduction on another platform"), PARCIAL,
        t("{s}: paquete de reproducción con huellas SHA-256.",
          "{s}: reproduction package with SHA-256 fingerprints.",
          s=sec("reproduccion")),
        t("Verificar el paquete en otra instalación (spinpy-informe "
          "--verificar) y declararlo.",
          "Verify the package on another installation (spinpy-informe "
          "--verificar) and state it.")))

    F.append(_fila(
        E, D["validacion"], t("Validación experimental",
                              "Experimental validation"), NO_CUMPLE,
        t("spinpy no contrasta sus predicciones con ensayos físicos; la "
          "coincidencia entre motores FEM es verificación, no validación.",
          "spinpy does not compare its predictions with physical tests; "
          "agreement between FEM engines is verification, not validation."),
        t("Declarar la ausencia de validación experimental, o comparar con un "
          "ensayo mecánico del mismo espécimen o con datos publicados "
          "comparables.",
          "State the absence of experimental validation, or compare with a "
          "mechanical test of the same specimen or with comparable "
          "published data.")))
    F.append(_fila(
        E, D["validacion"], t("Supuestos del modelo", "Model assumptions"),
        PARCIAL,
        t("{s}: tejido homogéneo, isótropo y elástico lineal; poros con "
          "rigidez residual; malla de vóxeles.",
          "{s}: homogeneous, isotropic, linear elastic tissue; pores with "
          "residual stiffness; voxel mesh.", s=sec("modelos") or
          sec("metodos")),
        t("Discutir cómo afectan estos supuestos al resultado principal.",
          "Discuss how these assumptions affect the main result.")))
    K = max([I._i(((x.get("incertidumbre") or {}).get("K"))) or 0
             for x in _recs(res, "ajuste")] or [0])
    n_disp = max([I._i(x.get("n_semillas")) or 0
                  for x in _recs(res, "dispersion")] or [0])
    accion_sens = t("Añadir un análisis de sensibilidad a las entradas que "
                    "no controla el generador: umbral de segmentación, E_s y "
                    "ν_s.",
                    "Add a sensitivity analysis to the inputs the generator "
                    "does not control: segmentation threshold, E_s and "
                    "ν_s.")
    from .incertidumbre import K_MIN
    if sens:
        accion_sens = t("Añadir la sensibilidad a E_s y ν_s, las entradas "
                        "del modelo que no se han variado.",
                        "Add the sensitivity to E_s and ν_s, the model "
                        "inputs that have not been varied.")
    if K >= K_MIN or n_disp or sens:
        partes = []
        if K >= K_MIN or n_disp:
            partes.append(t("variabilidad del generador estocástico con {k} "
                            "semillas", "variability of the stochastic "
                            "generator with {k} seeds", k=max(K, n_disp)))
        if sens:
            partes.append(t("sensibilidad a la posición de la superficie",
                            "sensitivity to the position of the surface"))
        F.append(_fila(E, D["validacion"], t("Incertidumbre y sensibilidad",
                                             "Uncertainty and sensitivity"),
                       PARCIAL, t("{s}: {p}.", "{s}: {p}.", s=sec("metodos"),
                                  p="; ".join(partes)),
                       accion_sens))
    else:
        F.append(_fila(E, D["validacion"], t("Incertidumbre y sensibilidad",
                                             "Uncertainty and sensitivity"),
                       NO_CUMPLE,
                       t("La sesión no contiene análisis de incertidumbre.",
                         "The session holds no uncertainty analysis."),
                       accion_sens))

    donde = t("{s}: repositorio {rep}, código con licencia MIT, autores con "
              "ORCID.",
              "{s}: repository {rep}, MIT-licensed code, authors with ORCID.",
              s=sec("cita"), rep=I.REPOSITORIO)
    if I.DOI_ARCHIVO == I.PENDIENTE:
        F.append(_fila(E, D["disponibilidad"],
                       t("Código, licencia y contacto",
                         "Code, licence and contact"), PARCIAL, donde,
                       t("Depositar la versión usada en un archivo con DOI "
                         "(por ejemplo, Zenodo) y citarlo.",
                         "Deposit the version used in an archive with a DOI "
                         "(e.g. Zenodo) and cite it.")))
    else:
        F.append(_fila(E, D["disponibilidad"],
                       t("Código, licencia y contacto",
                         "Code, licence and contact"), CUMPLE, donde))

    if ela:
        tbn = I._f(m_voi, "TbN")
        donde = t("{s}: homogeneización periódica sobre el VOI completo",
                  "{s}: periodic homogenisation over the whole VOI",
                  s=sec("metodos"))
        if lado and tbn:
            donde += t("; el lado equivale a {x} espaciados trabeculares "
                       "(lado·Tb.N)",
                       "; the side spans {x} trabecular spacings (side·Tb.N)",
                       x=r.n(lado * tbn, ".3g"))
        F.append(_fila(
            E, D["multiescala"], t("Representatividad del elemento de "
                                   "volumen", "Representativeness of the "
                                   "volume element"), PARCIAL, donde + ".",
            t("Justificar que el VOI es representativo (tamaño frente a la "
              "escala trabecular, sensibilidad al tamaño del VOI) y el uso de "
              "condiciones periódicas sobre una muestra de hueso.",
              "Justify that the VOI is representative (size relative to the "
              "trabecular scale, sensitivity to VOI size) and the use of "
              "periodic conditions on a bone sample.")))
    return F


def resumen(filas):
    c = {e: 0 for e in ORDEN}
    for f in filas:
        c[f["estado"]] += 1
    return c


def seccion(doc, idioma, r, numero, secciones=None):
    """Lineas de Markdown de la seccion, citando con el redactor `r`."""
    es = idioma == "es"
    i = 0 if es else 1

    def T(a, b):
        return a if es else b

    filas = evaluar(doc, idioma, r, secciones)
    L = [f"## {numero}. " + T("Lista de chequeo de reporte",
                              "Reporting checklist"), ""]
    L.append(r.t(
        "Esta sección contrasta el informe con dos guías de reporte. Para la "
        "morfometría por micro-TC se usan las recomendaciones de Bouxsein et "
        "al. {c1}, formuladas para hueso de roedor ex vivo y aplicadas aquí "
        "por extensión a otras especies. Para los modelos de elementos "
        "finitos se usan los parámetros de reporte de Erdemir et al. {c2}, "
        "agrupados en sus mismos dominios. Oefner et al. {c3} publicaron "
        "después una lista de verificación y validación para biomecánica "
        "ortopédica y traumatológica; sus ítems no se evalúan aquí, pero "
        "conviene revisarla antes del envío.",
        "This section checks the report against two reporting guidelines. "
        "Micro-CT morphometry is checked against the recommendations of "
        "Bouxsein et al. {c1}, written for ex vivo rodent bone and applied "
        "here by extension to other species. Finite element models are "
        "checked against the reporting parameters of Erdemir et al. {c2}, "
        "grouped in their own domains. Oefner et al. {c3} later published a "
        "verification and validation checklist for orthopaedic and trauma "
        "biomechanics; its items are not evaluated here, but it is worth "
        "reviewing before submission.",
        c1=r.c("bouxsein2010"), c2=r.c("erdemir2012"), c3=r.c("oefner2021")))
    L.append("")
    L.append(T(
        "El estado de cada ítem se refiere a lo que este informe ya "
        "documenta, no al manuscrito. Lo que spinpy no puede conocer, como el "
        "protocolo de adquisición o la especie, figura como «No cumple» junto "
        "con lo que hay que aportar.",
        "The status of each item refers to what this report already "
        "documents, not to the manuscript. What spinpy cannot know, such as "
        "the acquisition protocol or the species, appears as \"Not met\" "
        "together with what must be supplied."))
    L.append("")
    c = resumen(filas)
    pl = {e: "n" if c[e] != 1 else "" for e in ORDEN}
    L.append(T(f"**{c[CUMPLE]}** cumple{pl[CUMPLE]}, **{c[PARCIAL]}** "
               f"cumple{pl[PARCIAL]} parcialmente, **{c[NO_CUMPLE]}** no "
               f"cumple{pl[NO_CUMPLE]}, **{c[NO_APLICA]}** no "
               f"aplica{pl[NO_APLICA]}.",
               f"**{c[CUMPLE]}** met, **{c[PARCIAL]}** partially met, "
               f"**{c[NO_CUMPLE]}** not met, **{c[NO_APLICA]}** not "
               f"applicable."))
    L.append("")
    # Una tabla por guia: la columna del dominio queda corta y el nombre de
    # la guia no se repite en cada fila.
    for k, (guia, titulo) in enumerate(TITULOS.items(), 1):
        propias = [f for f in filas if f["guia"] == guia]
        L.append(f"### {numero}.{k} " + titulo[i])
        L.append("")
        L.append("| " + " | ".join(T(*h) for h in (
            ("Dominio", "Domain"), ("Ítem", "Item"), ("Estado", "Status"),
            ("En este informe", "In this report"),
            ("Qué falta", "What is missing"))) + " |")
        L.append("|---|---|---|---|---|")
        for f in propias:
            L.append("| " + " | ".join((
                f["dominio"][i], f["item"], ESTADOS[f["estado"]][i],
                f["donde"], f["accion"] or "")) + " |")
        L.append("")
    return L
