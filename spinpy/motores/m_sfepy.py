"""
m_sfepy.py — Motor SfePy (Python con extensiones C/Cython; BSD-3).

Usa los terminos de SfePy (`dw_lin_elastic`, `dw_surface_ltr`,
`dw_tl_he_svk`), su Newton y su capa de resolvedores:
  directo    `ls.scipy_direct` (SuperLU; UMFPACK no esta en pip);
  iterativo  `ls.scipy_iterative` (CG) con el precondicionador pyamg con modos
             rigidos inyectado por `setup_precond`, el mecanismo de SfePy para
             eso (`ls.pyamg` no admite el espacio nulo cercano).

No lineal: solo St. Venant-Kirchhoff. MEDIDO (comparativa_motores/): la
tangente de `dw_tl_he_svk` en SfePy 2026.3 difiere un 37 % de la derivada
numerica de su propio residuo, asi que Newton converge LINEALMENTE (22-27
iteraciones por paso en un bloque; no converge en 30 al 2 % en un
espinodoide). El residuo es correcto: cuando converge, la fuerza coincide con
la solucion cerrada. El neo-Hookeano de SfePy es la variante desacoplada
(mu/2 (J^-2/3 I1 - 3) + K/2 (J - 1)^2), que no es el de FEBio: no disponible.
"""

from __future__ import annotations

import logging

import numpy as np

from ._comun import (Cronometro, ErrorMotor, NoDisponible, emparejar,
                     modos_rigidos, solucion)

NOMBRE = "SfePy"
LICENCIA = "BSD-3-Clause"


def version():
    import sfepy
    return sfepy.__version__


def _sfepy():
    logging.disable(logging.WARNING)
    from sfepy.base.base import output
    output.set_output(quiet=True)
    import sfepy.discrete as sd
    return sd


def _problema(p):
    _sfepy()
    from sfepy.discrete import FieldVariable
    from sfepy.discrete.conditions import EssentialBC
    from sfepy.discrete.fem import FEDomain, Field, Mesh
    m = p["meta"]
    nodos, elems = p["nodos"], p["elems"]
    if elems.shape[1] == 8:
        esq, conn, desc, orden = np.arange(nodos.shape[0]), elems, "3_8", 1
    else:
        esq = np.unique(elems[:, :4])
        conn, desc, orden = np.searchsorted(esq, elems[:, :4]), "3_4", 2
    mesh = Mesh.from_data("m", nodos[esq], None, [conn.astype(np.int32)],
                          [np.ones(conn.shape[0], np.int32)], [desc])
    dom = FEDomain("d", mesh)
    z0, z1 = float(nodos[:, 2].min()), float(nodos[:, 2].max())
    tol = 1e-9 * max(m["H"], 1.0)
    omega = dom.create_region("Omega", "all")
    top = dom.create_region("Top", f"vertices in (z > {z1 - tol!r})", "facet")
    base = dom.create_region("Base", f"vertices in (z < {z0 + tol!r})",
                             "facet")
    ra = dom.create_region(
        "A", f"vertex {int(np.searchsorted(esq, m['ancla_xy']))}", "vertex")
    rb = dom.create_region(
        "B", f"vertex {int(np.searchsorted(esq, m['ancla_y']))}", "vertex")
    field = Field.from_args("u", np.float64, 3, omega, approx_order=orden)
    u = FieldVariable("u", "unknown", field)
    v = FieldVariable("v", "test", field, primary_var_name="u")
    if m["apoyo"] == "empotrado":
        ebcs = [EssentialBC("base", base, {"u.all": 0.0})]
    else:
        ebcs = [EssentialBC("base", base, {"u.2": 0.0}),
                EssentialBC("a", ra, {"u.[0,1]": 0.0}),
                EssentialBC("b", rb, {"u.1": 0.0})]
    return omega, top, field, u, v, ebcs


def _solver_lineal(modo, field, u, tol):
    from sfepy.solvers.ls import ScipyDirect, ScipyIterative
    if modo == "directo":
        return ScipyDirect({"method": "superlu"}), "SuperLU (ls.scipy_direct)"
    B_todo = modos_rigidos(field.get_coor())

    def precond(mtx, _ctx):
        import pyamg
        eq = u.eq_map.eq
        B = np.zeros((mtx.shape[0], 6))
        B[eq[eq >= 0]] = B_todo[eq >= 0]
        return pyamg.smoothed_aggregation_solver(
            mtx.tocsr(), B=B, max_coarse=500).aspreconditioner(cycle="V")

    return (ScipyIterative({"method": "cg", "i_max": 3000, "eps_r": tol,
                            "eps_a": 1e-300, "setup_precond": precond}),
            "CG (ls.scipy_iterative) + pyamg con modos rigidos")


def _reaccion(pb, expr, st, techo, **mats):
    """F = -sum(f_int,z) en el techo. SfePy anula las filas con condicion
    esencial al evaluar en modo 'weak': se evalua sin condiciones."""
    from sfepy.discrete.conditions import Conditions
    uvec = st()
    pb.time_update(ebcs=Conditions([]))
    pb.get_variables().init_state(uvec)
    r = pb.evaluate(expr, mode="weak", dw_mode="vector",
                    copy_materials=False, **mats)
    return float(-np.asarray(r).reshape(-1, 3)[techo, 2].sum()), uvec


def resolver(p):
    sd = _sfepy()
    from sfepy.discrete.conditions import Conditions, EssentialBC
    from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
    from sfepy.solvers.nls import Newton
    from sfepy.terms import Term
    m = p["meta"]
    if m["analisis"] == "nl" and m["material"] != "svk":
        raise NoDisponible("SfePy: solo St. Venant-Kirchhoff equivalente al "
                           "de FEBio")
    t = Cronometro()
    t("malla")
    omega, top, field, u, v, ebcs = _problema(p)
    plato = m["control"] == "plato"
    coor = field.get_coor()
    techo = np.nonzero(coor[:, 2] > coor[:, 2].max()
                       - 1e-9 * max(m["H"], 1.0))[0]
    mat = sd.Material("m", D=stiffness_from_youngpoisson(3, m["E"], m["nu"]))
    iv = sd.Integral("iv", order=2)
    modo = m.get("solver", "auto")
    if modo == "auto":
        modo = "directo" if (m["analisis"] == "nl" or field.n_nod * 3 < 6000) \
            else "iterativo"
    ls, nombre = _solver_lineal(modo, field, u, m.get("tol", 1e-10))
    lineal = m["analisis"] == "lineal"
    termino = "dw_lin_elastic" if lineal else "dw_tl_he_svk"
    t("montaje")
    tint = Term.new(f"{termino}(m.D, v, u)", iv, omega, m=mat, v=v, u=u)
    F, iters, st, uvec = [], [], None, None
    cargas = m["cargas"][-1:] if lineal else m["cargas"]
    status = {}
    t("solucion")
    for c in cargas:
        eq = tint
        bc = list(ebcs)
        if plato:
            bc.append(EssentialBC("plato", top, {"u.2": -c * m["H"]}))
        else:
            carga = sd.Material("c", val=np.array([[0.0], [0.0], [-c]]))
            eq = tint - Term.new("dw_surface_ltr(c.val, v)",
                                 sd.Integral("is", order=4), top, c=carga,
                                 v=v)
        pb = sd.Problem("p", equations=sd.Equations([sd.Equation("eq", eq)]))
        pb.set_bcs(ebcs=Conditions(bc))
        F_ref = m["E"] * m["A_bruta"]
        conf = ({"i_max": 1, "eps_a": 1e-300, "is_linear": True} if lineal
                else {"i_max": m.get("max_iter", 30), "eps_a": 1e-10 * F_ref,
                      "eps_r": 1e-10, "eps_mode": "or", "macheps": 1e-16,
                      "lin_red": None, "ls_red": 0.5, "ls_min": 1e-5,
                      "check": 0})
        pb.set_solver(Newton(conf, lin_solver=ls, status=status))
        st = pb.solve(state0=uvec, save_results=False)
        if not lineal and status.get("condition", 1) != 0:
            raise ErrorMotor(f"Newton de SfePy no convergio en la carga "
                             f"{c:g} ({status.get('n_iter')} iteraciones)")
        iters.append(int(status.get("n_iter", 1)))
        f, uvec = _reaccion(pb, f"{termino}.2.Omega(m.D, v, u)", st, techo,
                            m=mat)
        F.append(f)
    tiempos = t.fin()
    uu = uvec.reshape(-1, 3)[emparejar(coor, p["nodos"])]
    return solucion(uu, F, {"motor": "sfepy", "tiempos": tiempos,
                            "solver": ("" if lineal else "Newton (SfePy) + ")
                            + nombre,
                            "iteraciones_newton": None if lineal else iters,
                            "n_gdl": int(field.n_nod * 3)})
