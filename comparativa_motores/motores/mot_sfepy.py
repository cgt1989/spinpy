"""
mot_sfepy.py — Adaptador de SfePy (Python + extensiones en C/Cython).

SfePy trae terminos de elasticidad ya escritos (`dw_lin_elastic`,
`dw_surface_ltr`, `dw_tl_he_svk`), su propio Newton y una capa de
resolvedores lineales. Se usa todo de SfePy salvo el precondicionador:
  directo    `ls.scipy_direct` (SuperLU; UMFPACK no esta en pip);
  iterativo  `ls.scipy_iterative` con CG y el precondicionador pyamg con
             modos rigidos inyectado por `setup_precond`, el mecanismo que
             SfePy ofrece para eso. `ls.pyamg` de SfePy no admite el espacio
             nulo cercano, y sin el la agregacion trata el problema como uno
             escalar (ver `resistencia.py`, comentario del AMG).

No lineal: solo St. Venant-Kirchhoff. El neo-Hookeano de SfePy
(`dw_tl_he_neohook` + `dw_tl_bulk_penalty`) es la variante desacoplada,
W = mu/2 (J^-2/3 I1 - 3) + K/2 (J - 1)^2, que NO es la de FEBio
(W = mu/2 (I1 - 3) - mu ln J + lam/2 (ln J)^2): compararla con los demas
mezclaria motor y modelo, asi que se declara no disponible.
"""

from __future__ import annotations

import logging

import numpy as np

logging.disable(logging.WARNING)
from sfepy.base.base import output                           # noqa: E402

output.set_output(quiet=True)
from sfepy.discrete import (Equation, Equations, FieldVariable,  # noqa: E402
                            Integral, Material, Problem)
from sfepy.discrete.conditions import Conditions, EssentialBC  # noqa: E402
from sfepy.discrete.fem import FEDomain, Field, Mesh          # noqa: E402
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson  # noqa: E402
from sfepy.solvers.ls import ScipyDirect, ScipyIterative       # noqa: E402
from sfepy.solvers.nls import Newton                           # noqa: E402
from sfepy.terms import Term                                   # noqa: E402

from ._base import (Cronometro, emparejar, modos_rigidos,      # noqa: E402
                    principal)


def _problema(p, t):
    m = p["meta"]
    nodos, elems = p["nodos"], p["elems"]
    t("malla")
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
    ia = int(np.searchsorted(esq, m["ancla_xy"]))
    ib = int(np.searchsorted(esq, m["ancla_y"]))
    ra = dom.create_region("A", f"vertex {ia}", "vertex")
    rb = dom.create_region("B", f"vertex {ib}", "vertex")
    field = Field.from_args("u", np.float64, 3, omega, approx_order=orden)
    u = FieldVariable("u", "unknown", field)
    v = FieldVariable("v", "test", field, primary_var_name="u")
    if m["apoyo"] == "empotrado":
        ebcs = [EssentialBC("base", base, {"u.all": 0.0})]
    else:
        ebcs = [EssentialBC("base", base, {"u.2": 0.0}),
                EssentialBC("a", ra, {"u.[0,1]": 0.0}),
                EssentialBC("b", rb, {"u.1": 0.0})]
    return dom, omega, top, field, u, v, ebcs


def _solver_lineal(m, field, u):
    modo = m.get("solver", "iterativo")
    if modo == "directo":
        return ScipyDirect({"method": "superlu"}), "SuperLU (ls.scipy_direct)"
    coor = field.get_coor()
    B_todo = modos_rigidos(coor)

    def precond(mtx, _ctx):
        import pyamg
        eq = u.eq_map.eq
        B = np.zeros((mtx.shape[0], 6))
        B[eq[eq >= 0]] = B_todo[eq >= 0]
        ml = pyamg.smoothed_aggregation_solver(mtx.tocsr(), B=B,
                                               max_coarse=500)
        return ml.aspreconditioner(cycle="V")

    ls = ScipyIterative({"method": "cg", "i_max": 1000,
                         "eps_r": m.get("tol", 1e-10), "eps_a": 1e-300,
                         "setup_precond": precond})
    return ls, "CG (ls.scipy_iterative) + pyamg con modos rigidos"


def _u_spinpy(field, u_vec, nodos):
    i = emparejar(field.get_coor(), nodos)
    return u_vec.reshape(-1, 3)[i]


def resolver(p):
    m = p["meta"]
    t = Cronometro()
    dom, omega, top, field, u, v, ebcs = _problema(p, t)
    if m["analisis"] == "nl":
        return _no_lineal(p, t, dom, omega, top, field, u, v, ebcs)
    if m["analisis"] != "lineal":
        raise NotImplementedError(m["analisis"])
    t("montaje")
    mat = Material("m", D=stiffness_from_youngpoisson(3, m["E"], m["nu"]))
    carga = Material("c", val=np.array([[0.0], [0.0], [-m["p"]]]))
    iv = Integral("iv", order=2)
    ist = Integral("is", order=2)
    t1 = Term.new("dw_lin_elastic(m.D, v, u)", iv, omega, m=mat, v=v, u=u)
    t2 = Term.new("dw_surface_ltr(c.val, v)", ist, top, c=carga, v=v)
    pb = Problem("lineal", equations=Equations([Equation("eq", t1 - t2)]))
    pb.set_bcs(ebcs=Conditions(ebcs))
    ls, nombre = _solver_lineal(m, field, u)
    status = {}
    pb.set_solver(Newton({"i_max": 1, "eps_a": 1e-300, "is_linear": True},
                         lin_solver=ls, status=status))
    t("solucion")
    st = pb.solve(save_results=False)
    tiempos = t.fin()
    uu = _u_spinpy(field, st.get_state_parts()["u"], p["nodos"])
    # Residuo del sistema reducido, recalculado aqui (SfePy informa la norma
    # absoluta del residuo no lineal, no la relativa del sistema lineal).
    return uu, {"motor": "sfepy", "solver": nombre,
                "tiempos": tiempos, "n_gdl": int(field.n_nod * 3),
                "estado_newton": {k: (float(v) if isinstance(v, float) else
                                      None)
                                  for k, v in status.items()
                                  if k in ("err", "err0", "time")}}, None


def _no_lineal(p, t, dom, omega, top, field, u, v, ebcs):
    m = p["meta"]
    if m["material"] != "svk":
        raise NotImplementedError("SfePy: solo St. Venant-Kirchhoff "
                                  "equivalente al de FEBio")
    t("montaje")
    mat = Material("m", D=stiffness_from_youngpoisson(3, m["E"], m["nu"]))
    iv = Integral("iv", order=2)
    term = Term.new("dw_tl_he_svk(m.D, v, u)", iv, omega, m=mat, v=v, u=u)
    pb = Problem("nl", equations=Equations([Equation("eq", term)]))
    ls, nombre = _solver_lineal(dict(m, solver=m.get("solver", "directo")),
                                field, u)
    F_ref = m["E"] * m["A_bruta"]
    status = {}
    # eps_mode 'or': SfePy exige por omision las DOS tolerancias (absoluta
    # y relativa); con la relativa a 1e-12 el residuo se estanca en el piso
    # de redondeo y agota las iteraciones (medido en el bloque 4x4x6).
    nls = Newton({"i_max": m.get("max_iter", 30), "eps_a": 1e-10 * F_ref,
                  "eps_r": 1e-10, "eps_mode": "or", "macheps": 1e-16,
                  "lin_red": None,
                  "ls_red": 0.5, "ls_min": 1e-5, "check": 0},
                 lin_solver=ls, status=status)
    pb.set_solver(nls)
    coor = field.get_coor()
    z1 = float(coor[:, 2].max())
    techo = np.nonzero(coor[:, 2] > z1 - 1e-9 * max(m["H"], 1.0))[0]
    fuerzas, iters = [], []
    st = uvec = None
    t("solucion")
    for eps in m["eps_plato"]:
        bc = ebcs + [EssentialBC("plato", top, {"u.2": -eps * m["H"]})]
        pb.set_bcs(ebcs=Conditions(bc))
        st = pb.solve(state0=None if st is None else uvec,
                      save_results=False)
        if status.get("condition", 1) != 0:
            raise RuntimeError(f"Newton de SfePy no convergio a eps={eps}")
        iters.append(int(status["n_iter"]))
        # Reaccion: SfePy anula las filas de los GDL con condicion esencial
        # al evaluar en modo 'weak', asi que se evalua el residuo del estado
        # convergido con las condiciones de contorno retiradas.
        uvec = st()
        pb.time_update(ebcs=Conditions([]))
        pb.get_variables().init_state(uvec)
        r = pb.evaluate("dw_tl_he_svk.2.Omega(m.D, v, u)", mode="weak",
                        dw_mode="vector", m=mat, copy_materials=False)
        fuerzas.append(float(-np.asarray(r).reshape(-1, 3)[techo, 2].sum()))
    tiempos = t.fin()
    uu = _u_spinpy(field, st.get_state_parts()["u"], p["nodos"])
    return uu, {"motor": "sfepy", "solver": "Newton (SfePy) + " + nombre,
                "tiempos": tiempos, "iteraciones_newton": iters,
                "n_gdl": int(field.n_nod * 3)}, {"F_reac": np.array(fuerzas)}


if __name__ == "__main__":
    principal(resolver)
