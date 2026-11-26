"""
informe_cargas.py - Seccion «Casos de carga definidos por el usuario» del informe.

POR QUE UNA SECCION APARTE
--------------------------
El ensayo de compresion de siempre tiene su protocolo declarado en la seccion de
modelos (`informe_modelos`). Un caso de «Cargas con el raton…» no: las cargas, los
apoyos y su direccion los eligio quien uso la ventana, y un numero de von Mises o
de rigidez sin ese caso no se puede reproducir ni interpretar. Esta seccion lo
declara completo: la estructura, cada carga y cada apoyo con su region, el
solver, el equilibrio y las reacciones, y el caso entero en JSON para volver a
resolverlo con `spinpy.cargas.resolver_caso`.

QUE NO AFIRMA
-------------
No pasa por las puertas de citabilidad de `informe.comprobar`: esas puertas
validan el ensayo axial de la app contra una malla y un protocolo fijos. Aqui el
control es el equilibrio de fuerzas y momentos de cada solucion (la suma de
reacciones iguala a la de cargas) y el bloque 37 de la suite, que comprueba el
solver contra soluciones cerradas (compresion, traccion lateral, torsion y
flexion de una barra) y contra `ensayo_compresion`. El factor de Pistoia se
declara con su advertencia: se calibro con compresion de radio distal.

Sale solo si la sesion tiene casos (`resultados["cargas_casos"]`), como el resto
de subsecciones: no se describe un calculo que no se hizo.
"""

from __future__ import annotations

import json

from . import cargas as C


def _dir_txt(d, es):
    if isinstance(d, str):
        return {"normal_entrante": ("normal entrante", "inward normal"),
                "normal_saliente": ("normal saliente", "outward normal")
                }.get(d, (d, d))[0 if es else 1]
    v = [float(x) for x in d]
    return "(%s)" % ", ".join("%.3g" % x for x in v)


def _region_txt(r, es):
    f = r.get("forma")
    amb = "" if r.get("ambito", "superficie") == "superficie" else (
        ", volumen" if es else ", volume")
    if f == "cara":
        return (("cara %s del cubo" if es else "%s face of the cube")
                % r.get("cara")) + amb
    if f == "esfera":
        c = r["centro"]
        t = (("parche de radio %.3g mm centrado en (%.3g, %.3g, %.3g) mm"
              if es else
              "patch of radius %.3g mm centred at (%.3g, %.3g, %.3g) mm")
             % (r["radio"], *c))
        fn = r.get("filtro_normal")
        if fn:
            t += ((", caras a menos de %g° de la normal" if es else
                   ", faces within %g° of the normal") % fn.get("angulo", 60))
        return t + amb
    if f == "caja":
        return ("caja" if es else "box") + amb
    if f == "cilindro":
        return (("cilindro de radio %.3g mm" if es
                 else "cylinder of radius %.3g mm") % r["radio"]) + amb
    if f == "nodos":
        return ("%d nodo(s) elegidos" if es else "%d chosen node(s)"
                ) % len(r.get("ids", []))
    return str(f)


def _carga_txt(c, i, es):
    t = c["tipo"]
    reg = _region_txt(c["region"], es)
    if t == "fuerza":
        return (("C%d: fuerza de %.4g N hacia %s sobre %s" if es
                 else "C%d: force of %.4g N towards %s on %s")
                % (i + 1, c["magnitud"], _dir_txt(
                    c.get("direccion", "normal_entrante"), es), reg))
    if t == "presion":
        return (("C%d: presión de %.4g MPa (%s) sobre %s" if es
                 else "C%d: pressure of %.4g MPa (%s) on %s")
                % (i + 1, c["magnitud"], c.get("area_ref", "bruta"), reg))
    if t == "momento":
        return (("C%d: momento de %.4g N·mm alrededor de %s sobre %s" if es
                 else "C%d: moment of %.4g N·mm about %s on %s")
                % (i + 1, c["magnitud"], _dir_txt(c["eje"], es), reg))
    u = c.get("u")
    if u is not None:
        comp = ", ".join(("libre" if es else "free") if x is None
                         else "%.4g" % x for x in u)
        return (("C%d: desplazamiento impuesto (%s) mm sobre %s" if es
                 else "C%d: prescribed displacement (%s) mm on %s")
                % (i + 1, comp, reg))
    return (("C%d: desplazamiento de %.4g mm hacia %s sobre %s" if es
             else "C%d: displacement of %.4g mm towards %s on %s")
            % (i + 1, c["magnitud"], _dir_txt(c.get("direccion", "+z"), es),
               reg))


def _apoyo_txt(a, i, es):
    g = C._gdl(a.get("gdl"))
    gtxt = (" (%s)" % ",".join(C.EJES[k] for k in g)) if g else ""
    nombre = {"empotrado": ("empotrado", "fixed"),
              "deslizante": ("deslizante", "sliding"),
              "personalizado": ("personalizado", "custom")}[a["tipo"]]
    return ("A%d: %s%s %s %s" % (i + 1, nombre[0 if es else 1], gtxt,
                                 "sobre" if es else "on",
                                 _region_txt(a["region"], es)))


def seccion_cargas(doc, idioma, r, numero):
    """Lineas de Markdown de la seccion; [] si la sesion no tiene casos."""
    casos = [x for x in ((doc.get("resultados") or {}).get("cargas_casos")
                         or []) if isinstance(x, dict) and x.get("ok")]
    if not casos:
        return []
    es = idioma == "es"

    def T(a, b):
        return a if es else b

    L = [f"## {numero}. " + T("Casos de carga definidos por el usuario",
                              "User-defined load cases"), ""]
    L.append(r.t(
        "Además del ensayo de compresión axial, la sesión resolvió {n} "
        "caso(s) de carga elegidos con el ratón. Cada caso es una lista de "
        "cargas y de apoyos, con su región, dirección y magnitud, aplicada "
        "sobre una estructura de hexaedros trilineales de ocho nodos (un "
        "elemento por vóxel de hueso, integración de Gauss 2×2×2) con "
        "material elástico lineal e isótropo. Una fuerza o una presión se "
        "reparte entre los nodos de las caras de vóxel expuestas de su "
        "región según el área tributaria de cada nodo, el equivalente nodal "
        "consistente de una presión; un momento se aplica como un par de "
        "fuerza neta nula; un desplazamiento impuesto es una condición "
        "cinemática. Si los apoyos no impiden el movimiento de sólido rígido "
        "se añaden las restricciones mínimas y se exige que sus reacciones "
        "sean nulas, de modo que un caso con las cargas desequilibradas se "
        "rechaza. Solo se resuelven los fragmentos de hueso (conexión por "
        "caras) que unen cargas y apoyos.",
        "In addition to the axial compression test, the session solved {n} "
        "load case(s) chosen with the mouse. Each case is a list of loads "
        "and supports, with their region, direction and magnitude, applied "
        "to a mesh of eight-node trilinear hexahedra (one element per bone "
        "voxel, 2×2×2 Gauss integration) with a linear isotropic elastic "
        "material. A force or pressure is spread over the nodes of the "
        "exposed voxel faces of its region according to each node's "
        "tributary area, the consistent nodal equivalent of a pressure; a "
        "moment is applied as a couple with zero net force; a prescribed "
        "displacement is a kinematic condition. If the supports do not "
        "prevent rigid-body motion, minimal constraints are added and their "
        "reactions are required to vanish, so a case with unbalanced loads "
        "is rejected. Only the bone fragments (face connectivity) that "
        "link loads to supports are solved.", n=len(casos)))
    L.append("")

    for k, c in enumerate(casos, 1):
        caso = c["caso"]
        est = c.get("estructura") or {}
        mat = caso.get("material", {})
        L.append(f"### {numero}.{k} " + T("Caso %d" % k, "Case %d" % k)
                 + (f" · {est['etiqueta']}" if est.get("etiqueta") else ""))
        L.append("")
        if est:
            L.append(r.t(
                "Estructura: {forma} vóxeles, espaciado {sp} mm, BV/TV "
                "{rho}. Material: E = {E} GPa, ν = {nu}.",
                "Structure: {forma} voxels, spacing {sp} mm, BV/TV {rho}. "
                "Material: E = {E} GPa, ν = {nu}.",
                forma=" × ".join(str(x) for x in est["forma"]),
                sp=" × ".join(r.n(x, ".3g") for x in est["spacing_mm"]),
                rho=r.n(est.get("BV_TV"), ".3f"),
                E=r.n(float(mat.get("E_s", C.E_S_DEF)) / 1e9, ".4g"),
                nu=r.n(mat.get("nu", C.NU_S_DEF), ".3g")))
            L.append("")
        for i, cc in enumerate(caso["cargas"]):
            L.append("- " + _carga_txt(cc, i, es))
        for i, aa in enumerate(caso["apoyos"]):
            L.append("- " + _apoyo_txt(aa, i, es))
        L.append("")
        eq = c.get("equilibrio") or {}
        L.append(r.t(
            "El sistema ({ne} elementos, {nd} grados de libertad) se resolvió "
            "con {solver}. Residuo de equilibrio: fuerzas {rf}, momentos "
            "{rm} (suma de reacciones frente a la de cargas, relativa).",
            "The system ({ne} elements, {nd} degrees of freedom) was solved "
            "with {solver}. Equilibrium residual: forces {rf}, moments {rm} "
            "(sum of reactions against that of the loads, relative).",
            ne=c.get("n_elem", "n/d"), nd=c.get("n_dof", "n/d"),
            solver=c.get("solver", "n/d"),
            rf=r.n(eq.get("residuo_F"), ".1e"),
            rm=r.n(eq.get("residuo_M"), ".1e")))
        L.append("")
        if c.get("reacciones"):
            L.append("| " + T("Origen", "Origin") + " | Fx (N) | Fy (N) | "
                     "Fz (N) | |F| (N) |")
            L.append("|---|---:|---:|---:|---:|")
            for x in c["reacciones"]:
                etq = ("A%d" if x["origen"] == "apoyo" else "C%d") % (
                    x["indice"] + 1)
                L.append("| %s | %s | %s | %s | %s |" % (
                    etq, *(r.n(v, ".4g") for v in x["F_N"]),
                    r.n(x["modulo_N"], ".4g")))
            L.append("")
        st = c.get("stats") or {}
        if "vm" in st:
            L.append(r.t(
                "Tensión de von Mises en el hueso resuelto: máximo {mx} MPa, "
                "percentil 99 {p99} MPa. El máximo de un campo con cargas "
                "concentradas depende de la malla y no es un valor "
                "citable; el percentil 99 es más estable.",
                "von Mises stress in the solved bone: maximum {mx} MPa, "
                "99th percentile {p99} MPa. The maximum of a field with "
                "concentrated loads depends on the mesh and is not a "
                "citable value; the 99th percentile is more stable.",
                mx=r.n(st["vm"]["max"], ".4g"),
                p99=r.n(st["vm"]["p99"], ".4g")))
            L.append("")
        for g in c.get("rigidez") or []:
            if g.get("rigidez_N_mm") not in (None, float("inf")):
                L.append(r.t(
                    "Rigidez de C{n}: {k} N/mm.",
                    "Stiffness of C{n}: {k} N/mm.", n=g["indice"] + 1,
                    k=r.n(g["rigidez_N_mm"], ".4g")))
                L.append("")
        p = c.get("pistoia")
        if p:
            L.append(r.t(
                "Factor de escala hasta el criterio de Pistoia: {f}. {nota}",
                "Scale factor to the Pistoia criterion: {f}. {nota}",
                f=r.n(p["factor"], ".3g"),
                nota=T(
                    "Se calibró con compresión de radio distal; con cargas "
                    "concentradas o desplazamientos impuestos no está "
                    "validado y se da como índice relativo.",
                    "It was calibrated with distal radius compression; with "
                    "concentrated loads or prescribed displacements it is "
                    "not validated and is given as a relative index.")))
            L.append("")
        for a in c.get("avisos") or []:
            L.append(f"> {a}")
            L.append("")
        L.append(T("Caso de carga completo (JSON), para volver a resolverlo "
                   "con `spinpy.cargas.resolver_caso`:",
                   "Full load case (JSON), to solve it again with "
                   "`spinpy.cargas.resolver_caso`:"))
        L.append("")
        L.append("```json")
        L.append(json.dumps(caso, indent=1, ensure_ascii=False))
        L.append("```")
        L.append("")
    return L
