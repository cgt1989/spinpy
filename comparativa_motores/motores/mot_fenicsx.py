"""
mot_fenicsx.py — Adaptador de FEniCSx (DOLFINx + UFL + FFCx + PETSc).

Se ejecuta con el Python del entorno de conda-forge: FEniCSx no se distribuye
en PyPI. Las formas se escriben en UFL; FFCx las traduce a C y CFFI las
COMPILA la primera vez que se usan (hace falta un compilador de C en tiempo de
ejecucion; queda en cache en ~/.cache/fenics). El tiempo de esa compilacion se
mide aparte ("jit").

  lineal     a(u, v) = sigma(u) : eps(v) dx,  L(v) = -p v_z ds(techo)
  no lineal  Pi(u) = W(F) dx, residuo = derivative(Pi), tangente =
             derivative(residuo); Newton de PETSc (SNES newtonls).

Resolvedores de PETSc: directo = MUMPS; iterativo = CG + GAMG con los seis
modos rigidos como espacio nulo cercano (la receta de los ejemplos de
elasticidad de DOLFINx).
"""

from __future__ import annotations

import os
import time

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

import basix.ufl
import dolfinx
import ufl
from dolfinx import fem, mesh as dmesh
from dolfinx.fem import petsc as fpetsc

from ._base import Cronometro, emparejar, principal

# C3D8 (spinpy) -> orden de vertices de DOLFINx (producto tensorial)
PERM_HEX = [0, 1, 3, 2, 4, 5, 7, 6]


def _malla(p):
    nodos, elems = p["nodos"], p["elems"]
    if elems.shape[1] == 8:
        celda, cells, orden = "hexahedron", elems[:, PERM_HEX], 1
        esq = np.arange(nodos.shape[0])
    else:
        esq = np.unique(elems[:, :4])
        celda, cells, orden = ("tetrahedron",
                               np.searchsorted(esq, elems[:, :4]), 2)
    dom = ufl.Mesh(basix.ufl.element("Lagrange", celda, 1, shape=(3,)))
    msh = dmesh.create_mesh(MPI.COMM_SELF, cells.astype(np.int64),
                            dom, nodos[esq])
    return msh, orden


def _contorno(p, msh, V, extra_techo=None):
    m = p["meta"]
    z = p["nodos"][:, 2]
    z0, z1 = float(z.min()), float(z.max())
    tol = 1e-9 * max(m["H"], 1.0)
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    f_base = dmesh.locate_entities_boundary(
        msh, tdim - 1, lambda x: np.abs(x[2] - z0) < tol)
    f_techo = dmesh.locate_entities_boundary(
        msh, tdim - 1, lambda x: np.abs(x[2] - z1) < tol)
    bcs = []
    cero = PETSc.ScalarType(0.0)
    if m["apoyo"] == "empotrado":
        for c in range(3):
            d = fem.locate_dofs_topological(V.sub(c), tdim - 1, f_base)
            bcs.append(fem.dirichletbc(cero, d, V.sub(c)))
    else:
        d = fem.locate_dofs_topological(V.sub(2), tdim - 1, f_base)
        bcs.append(fem.dirichletbc(cero, d, V.sub(2)))
        for nodo, comps in ((m["ancla_xy"], (0, 1)), (m["ancla_y"], (1,))):
            X = p["nodos"][nodo]
            v = dmesh.locate_entities(
                msh, 0, lambda x, X=X: np.linalg.norm(x.T - X, axis=1) < tol)
            for c in comps:
                d = fem.locate_dofs_topological(V.sub(c), 0, v)
                bcs.append(fem.dirichletbc(cero, d, V.sub(c)))
    mt = dmesh.meshtags(msh, tdim - 1, np.sort(f_techo),
                        np.full(f_techo.size, 1, np.int32))
    if extra_techo is not None:
        d = fem.locate_dofs_topological(V.sub(2), tdim - 1, f_techo)
        bcs.append(fem.dirichletbc(extra_techo, d, V.sub(2)))
    return bcs, mt, f_techo


def _nulo_cercano(V, A):
    """Seis modos rigidos como MatNullSpace (vectores ortonormalizados)."""
    x = V.tabulate_dof_coordinates()
    bs = V.dofmap.index_map_bs
    n = V.dofmap.index_map.size_local
    vecs = [A.createVecLeft() for _ in range(6)]
    arrs = [np.zeros((n, bs)) for _ in range(6)]
    for c in range(3):
        arrs[c][:, c] = 1.0
    X = x[:n]
    arrs[3][:, 0], arrs[3][:, 1] = -X[:, 1], X[:, 0]
    arrs[4][:, 1], arrs[4][:, 2] = -X[:, 2], X[:, 1]
    arrs[5][:, 0], arrs[5][:, 2] = X[:, 2], -X[:, 0]
    for v, a in zip(vecs, arrs):
        v.setArray(a.ravel())
    # Gram-Schmidt: MatNullSpace exige base ortonormal.
    for i, v in enumerate(vecs):
        for w in vecs[:i]:
            v.axpy(-v.dot(w), w)
        v.normalize()
    return PETSc.NullSpace().create(vectors=vecs)


def _ksp(A, V, modo, tol):
    ksp = PETSc.KSP().create(MPI.COMM_SELF)
    ksp.setOperators(A)
    if modo == "directo":
        ksp.setType("preonly")
        ksp.getPC().setType("cholesky")
        ksp.getPC().setFactorSolverType("mumps")
        nombre = "Cholesky MUMPS (PETSc)"
    else:
        A.setNearNullSpace(_nulo_cercano(V, A))
        ksp.setType("cg")
        # Residuo NO precondicionado: la misma medida ||b - A x|| / ||b||
        # que la tolerancia de los demas motores.
        ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
        ksp.setTolerances(rtol=tol, atol=0.0, max_it=2000)
        pc = ksp.getPC()
        pc.setType("gamg")
        opts = PETSc.Options()
        opts["mg_levels_ksp_type"] = "chebyshev"
        opts["mg_levels_pc_type"] = "jacobi"
        opts["pc_gamg_threshold"] = 0.01
        ksp.setFromOptions()
        nombre = "CG + GAMG (PETSc, modos rigidos)"
    return ksp, nombre


def _a_spinpy(V, uh, nodos):
    i = emparejar(V.tabulate_dof_coordinates(), nodos)
    return uh.x.array.reshape(-1, 3)[i]


def _material(m):
    E, nu = m["E"], m["nu"]
    return E / 2 / (1 + nu), E * nu / (1 + nu) / (1 - 2 * nu)


def resolver(p):
    m = p["meta"]
    t = Cronometro()
    t("malla")
    msh, orden = _malla(p)
    V = fem.functionspace(msh, ("Lagrange", orden, (3,)))
    if m["analisis"] == "nl":
        return _no_lineal(p, t, msh, V, orden)
    if m["analisis"] != "lineal":
        raise NotImplementedError(m["analisis"])
    bcs, mt, _ = _contorno(p, msh, V)
    mu, lam = _material(m)
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    def eps(w):
        return ufl.sym(ufl.grad(w))

    def sig(w):
        return 2 * mu * eps(w) + lam * ufl.tr(eps(w)) * ufl.Identity(3)

    ds = ufl.Measure("ds", domain=msh, subdomain_data=mt)
    t("jit")
    a = fem.form(ufl.inner(sig(u), eps(v)) * ufl.dx)
    L = fem.form(-m["p"] * v[2] * ds(1))
    t("montaje")
    A = fpetsc.assemble_matrix(a, bcs=bcs)
    A.assemble()
    b = fpetsc.assemble_vector(L)
    fpetsc.apply_lifting(b, [a], bcs=[bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    fpetsc.set_bc(b, bcs)
    t("solucion")
    ksp, nombre = _ksp(A, V, m.get("solver", "iterativo"), m.get("tol", 1e-10))
    uh = fem.Function(V)
    ksp.solve(b, uh.x.petsc_vec)
    uh.x.scatter_forward()
    r = b.duplicate()
    A.mult(uh.x.petsc_vec, r)
    r.axpy(-1.0, b)
    info = {"residuo_rel": float(r.norm() / b.norm()),
            "iteraciones": int(ksp.getIterationNumber()),
            "razon_ksp": int(ksp.getConvergedReason())}
    if info["razon_ksp"] < 0:
        raise RuntimeError(f"PETSc KSP no convergio ({info['razon_ksp']})")
    t("post")
    uu = _a_spinpy(V, uh, p["nodos"])
    tiempos = t.fin()
    return uu, {"motor": "fenicsx", "solver": nombre, **info,
                "tiempos": tiempos, "n_gdl": int(V.dofmap.index_map.size_local
                                                 * 3),
                "version": dolfinx.__version__}, None


def _no_lineal(p, t, msh, V, orden):
    m = p["meta"]
    mu, lam = _material(m)
    uz = fem.Constant(msh, PETSc.ScalarType(0.0))
    bcs, _, f_techo = _contorno(p, msh, V, extra_techo=None)
    tdim = msh.topology.dim
    d_techo = fem.locate_dofs_topological(V.sub(2), tdim - 1, f_techo)
    bcs.append(fem.dirichletbc(uz, d_techo, V.sub(2)))
    uh = fem.Function(V)
    v = ufl.TestFunction(V)
    Id = ufl.Identity(3)
    F = Id + ufl.grad(uh)
    C = F.T * F
    if m["material"] == "svk":
        Eg = 0.5 * (C - Id)
        W = lam / 2 * ufl.tr(Eg) ** 2 + mu * ufl.inner(Eg, Eg)
    else:
        J = ufl.det(F)
        W = mu / 2 * (ufl.tr(C) - 3) - mu * ufl.ln(J) + lam / 2 * ufl.ln(J) ** 2
    md = {"quadrature_degree": 2} if orden == 1 else {}
    dx = ufl.dx(metadata=md)
    Pi = W * dx
    R = ufl.derivative(Pi, uh, v)
    t("jit")
    R_form = fem.form(R)
    prob = fpetsc.NonlinearProblem(
        R, uh, bcs=bcs, petsc_options_prefix="nl_",
        petsc_options={"snes_type": "newtonls", "snes_linesearch_type": "none",
                       "snes_rtol": 1e-10, "snes_atol": 1e-12 * m["E"]
                       * m["A_bruta"], "snes_stol": 0.0, "snes_max_it": 30,
                       "ksp_type": "preonly", "pc_type": "lu",
                       "pc_factor_mat_solver_type": "mumps"})
    t("solucion")
    techo_z = d_techo
    fuerzas, iters = [], []
    for eps in m["eps_plato"]:
        uz.value = -eps * m["H"]
        prob.solve()
        snes = prob.solver
        if snes.getConvergedReason() <= 0:
            raise RuntimeError(f"SNES no convergio ({snes.getConvergedReason()})"
                               f" a eps={eps}")
        iters.append(int(snes.getIterationNumber()))
        rv = fpetsc.assemble_vector(R_form)
        rv.ghostUpdate(addv=PETSc.InsertMode.ADD,
                       mode=PETSc.ScatterMode.REVERSE)
        arr = rv.getArray()
        # d_techo son indices de GDL «desbloqueados» de V.sub(2)
        fuerzas.append(float(-arr[techo_z].sum()))
    t("post")
    uu = _a_spinpy(V, uh, p["nodos"])
    tiempos = t.fin()
    return uu, {"motor": "fenicsx", "solver": "Newton SNES + MUMPS (PETSc)",
                "tiempos": tiempos, "iteraciones_newton": iters,
                "n_gdl": int(V.dofmap.index_map.size_local * 3),
                "version": dolfinx.__version__}, {"F_reac": np.array(fuerzas)}


if __name__ == "__main__":
    principal(resolver)
