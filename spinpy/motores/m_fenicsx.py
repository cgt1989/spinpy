"""
m_fenicsx.py: Motor FEniCSx (DOLFINx + UFL + FFCx + PETSc; LGPL-3.0).

No esta en PyPI: se instala con conda-forge (`conda install -c conda-forge
fenics-dolfinx`). Las formas se escriben en UFL; FFCx las traduce a C y CFFI
las COMPILA la primera vez (hace falta un compilador de C en tiempo de
ejecucion; quedan en cache en ~/.cache/fenics). Ese tiempo se mide aparte
("jit").

  lineal     a(u, v) = sigma(u) : eps(v) dx,  L(v) = -p v_z ds(techo)
  no lineal  Pi(u) = W(F) dx + p u_z ds(techo); residuo = derivative(Pi),
             tangente derivada por UFL; Newton de PETSc (SNES newtonls).

Resolvedores de PETSc: directo = Cholesky de MUMPS; iterativo = CG + GAMG
con los seis modos rigidos como espacio nulo cercano.
"""

from __future__ import annotations

import numpy as np

from ._comun import Cronometro, ErrorMotor, emparejar, lame, solucion

NOMBRE = "FEniCSx"
LICENCIA = "LGPL-3.0"
# C3D8 (spinpy) -> orden de vertices de DOLFINx (producto tensorial)
PERM_HEX = [0, 1, 3, 2, 4, 5, 7, 6]


def version():
    import dolfinx
    return dolfinx.__version__


def _malla(p):
    import basix.ufl
    import ufl
    from dolfinx import mesh as dmesh
    from mpi4py import MPI
    nodos, elems = p["nodos"], p["elems"]
    if elems.shape[1] == 8:
        celda, cells, orden = "hexahedron", elems[:, PERM_HEX], 1
        esq = np.arange(nodos.shape[0])
    else:
        esq = np.unique(elems[:, :4])
        celda, cells, orden = ("tetrahedron",
                               np.searchsorted(esq, elems[:, :4]), 2)
    dom = ufl.Mesh(basix.ufl.element("Lagrange", celda, 1, shape=(3,)))
    return dmesh.create_mesh(MPI.COMM_SELF, cells.astype(np.int64), dom,
                             nodos[esq]), orden


def _contorno(p, msh, V):
    from dolfinx import fem, mesh as dmesh
    from petsc4py import PETSc
    m = p["meta"]
    z = p["nodos"][:, 2]
    z0, z1 = float(z.min()), float(z.max())
    tol = 1e-9 * max(m["H"], 1.0)
    fd = msh.topology.dim - 1
    msh.topology.create_connectivity(fd, fd + 1)
    f_base = dmesh.locate_entities_boundary(
        msh, fd, lambda x: np.abs(x[2] - z0) < tol)
    f_techo = dmesh.locate_entities_boundary(
        msh, fd, lambda x: np.abs(x[2] - z1) < tol)
    cero = PETSc.ScalarType(0.0)
    bcs = []
    if m["apoyo"] == "empotrado":
        for c in range(3):
            d = fem.locate_dofs_topological(V.sub(c), fd, f_base)
            bcs.append(fem.dirichletbc(cero, d, V.sub(c)))
    else:
        d = fem.locate_dofs_topological(V.sub(2), fd, f_base)
        bcs.append(fem.dirichletbc(cero, d, V.sub(2)))
        for nodo, comps in ((m["ancla_xy"], (0, 1)), (m["ancla_y"], (1,))):
            X = p["nodos"][nodo]
            v = dmesh.locate_entities(
                msh, 0, lambda x, X=X: np.linalg.norm(x.T - X, axis=1) < tol)
            for c in comps:
                d = fem.locate_dofs_topological(V.sub(c), 0, v)
                bcs.append(fem.dirichletbc(cero, d, V.sub(c)))
    mt = dmesh.meshtags(msh, fd, np.sort(f_techo),
                        np.full(f_techo.size, 1, np.int32))
    d_techo = fem.locate_dofs_topological(V.sub(2), fd, f_techo)
    return bcs, mt, d_techo


def _nulo_cercano(V, A):
    from petsc4py import PETSc
    x = V.tabulate_dof_coordinates()
    n = V.dofmap.index_map.size_local
    X = x[:n]
    arrs = [np.zeros((n, 3)) for _ in range(6)]
    for c in range(3):
        arrs[c][:, c] = 1.0
    arrs[3][:, 0], arrs[3][:, 1] = -X[:, 1], X[:, 0]
    arrs[4][:, 1], arrs[4][:, 2] = -X[:, 2], X[:, 1]
    arrs[5][:, 0], arrs[5][:, 2] = X[:, 2], -X[:, 0]
    vecs = []
    for a in arrs:
        v = A.createVecLeft()
        v.setArray(a.ravel())
        for w in vecs:
            v.axpy(-v.dot(w), w)
        v.normalize()
        vecs.append(v)
    return PETSc.NullSpace().create(vectors=vecs)


def _ksp(A, V, modo, tol):
    from mpi4py import MPI
    from petsc4py import PETSc
    ksp = PETSc.KSP().create(MPI.COMM_SELF)
    ksp.setOperators(A)
    if modo == "directo":
        ksp.setType("preonly")
        ksp.getPC().setType("cholesky")
        ksp.getPC().setFactorSolverType("mumps")
        return ksp, "Cholesky MUMPS (PETSc)"
    A.setNearNullSpace(_nulo_cercano(V, A))
    ksp.setType("cg")
    ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
    ksp.setTolerances(rtol=tol, atol=0.0, max_it=3000)
    ksp.getPC().setType("gamg")
    opts = PETSc.Options()
    opts["mg_levels_ksp_type"] = "chebyshev"
    opts["mg_levels_pc_type"] = "jacobi"
    opts["pc_gamg_threshold"] = 0.01
    ksp.setFromOptions()
    return ksp, "CG + GAMG (PETSc, modos rigidos)"


def resolver(p):
    import ufl
    from dolfinx import fem
    from dolfinx.fem import petsc as fpetsc
    from petsc4py import PETSc
    m = p["meta"]
    t = Cronometro()
    t("malla")
    msh, orden = _malla(p)
    V = fem.functionspace(msh, ("Lagrange", orden, (3,)))
    bcs, mt, d_techo = _contorno(p, msh, V)
    plato = m["control"] == "plato"
    uz = fem.Constant(msh, PETSc.ScalarType(0.0))
    carga = fem.Constant(msh, PETSc.ScalarType(0.0))
    if plato:
        bcs.append(fem.dirichletbc(uz, d_techo, V.sub(2)))
    lam, mu = lame(m["E"], m["nu"])
    ds = ufl.Measure("ds", domain=msh, subdomain_data=mt)
    v = ufl.TestFunction(V)
    uh = fem.Function(V)
    modo = m.get("solver", "auto")
    info = {}
    if m["analisis"] == "lineal":
        if modo == "auto":
            modo = "iterativo" if 3 * V.dofmap.index_map.size_local > 3e5 \
                else "directo"
        u = ufl.TrialFunction(V)

        def eps(w):
            return ufl.sym(ufl.grad(w))

        t("jit")
        a = fem.form(ufl.inner(2 * mu * eps(u) + lam * ufl.tr(eps(u))
                               * ufl.Identity(3), eps(v)) * ufl.dx)
        L = fem.form(-carga * v[2] * ds(1))
        if plato:
            uz.value = -m["cargas"][-1] * m["H"]
        else:
            carga.value = m["cargas"][-1]
        t("montaje")
        A = fpetsc.assemble_matrix(a, bcs=bcs)
        A.assemble()
        b = fpetsc.assemble_vector(L)
        fpetsc.apply_lifting(b, [a], bcs=[bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD,
                      mode=PETSc.ScatterMode.REVERSE)
        fpetsc.set_bc(b, bcs)
        t("solucion")
        ksp, info["solver"] = _ksp(A, V, modo, m.get("tol", 1e-10))
        ksp.solve(b, uh.x.petsc_vec)
        uh.x.scatter_forward()
        if ksp.getConvergedReason() < 0 and m.get("solver", "auto") == "auto":
            # MEDIDO: con la malla TET10 de un espinodoide real, CG + GAMG
            # diverge con precondicionador indefinido (codigo -8); en modo
            # automatico se repite con el directo en lugar de fallar.
            info["fallo_iterativo"] = int(ksp.getConvergedReason())
            ksp, info["solver"] = _ksp(A, V, "directo", m.get("tol", 1e-10))
            ksp.solve(b, uh.x.petsc_vec)
            uh.x.scatter_forward()
        if ksp.getConvergedReason() < 0:
            raise ErrorMotor(f"PETSc KSP no convergio "
                             f"({ksp.getConvergedReason()})")
        r = b.duplicate()
        A.mult(uh.x.petsc_vec, r)
        r.axpy(-1.0, b)
        info.update(residuo_rel=float(r.norm() / b.norm()),
                    iteraciones=int(ksp.getIterationNumber()))
        fint = fpetsc.assemble_vector(fem.form(ufl.action(
            ufl.inner(2 * mu * eps(u) + lam * ufl.tr(eps(u)) * ufl.Identity(3),
                      eps(v)) * ufl.dx, uh)))
        F = [float(-fint.getArray()[d_techo].sum())]
        iters = None
    else:
        Id = ufl.Identity(3)
        Fd = Id + ufl.grad(uh)
        C = Fd.T * Fd
        if m["material"] == "svk":
            Eg = 0.5 * (C - Id)
            W = lam / 2 * ufl.tr(Eg) ** 2 + mu * ufl.inner(Eg, Eg)
        else:
            J = ufl.det(Fd)
            W = (mu / 2 * (ufl.tr(C) - 3) - mu * ufl.ln(J)
                 + lam / 2 * ufl.ln(J) ** 2)
        dx = ufl.dx(metadata={"quadrature_degree": 2} if orden == 1 else {})
        Rint = ufl.derivative(W * dx, uh, v)
        R = Rint + carga * v[2] * ds(1)
        t("jit")
        Rint_f = fem.form(Rint)
        prob = fpetsc.NonlinearProblem(
            R, uh, bcs=bcs, petsc_options_prefix="spinpy_nl_",
            petsc_options={
                "snes_type": "newtonls", "snes_linesearch_type": "none",
                "snes_rtol": 1e-10, "snes_atol": 1e-12 * m["E"]
                * m["A_bruta"], "snes_stol": 0.0,
                "snes_max_it": m.get("max_iter", 30),
                "ksp_type": "preonly", "pc_type": "lu",
                "pc_factor_mat_solver_type": "mumps"})
        F, iters = [], []
        t("solucion")
        for c in m["cargas"]:
            if plato:
                uz.value = -c * m["H"]
            else:
                carga.value = c
            prob.solve()
            snes = prob.solver
            if snes.getConvergedReason() <= 0:
                raise ErrorMotor(f"SNES no convergio "
                                 f"({snes.getConvergedReason()}) en la "
                                 f"carga {c:g}")
            iters.append(int(snes.getIterationNumber()))
            rv = fpetsc.assemble_vector(Rint_f)
            rv.ghostUpdate(addv=PETSc.InsertMode.ADD,
                           mode=PETSc.ScatterMode.REVERSE)
            F.append(float(-rv.getArray()[d_techo].sum()))
        info["solver"] = "Newton SNES + MUMPS (PETSc)"
    t("post")
    uu = uh.x.array.reshape(-1, 3)[emparejar(V.tabulate_dof_coordinates(),
                                             p["nodos"])]
    tiempos = t.fin()
    return solucion(uu, F, {"motor": "fenicsx", "tiempos": tiempos,
                            "iteraciones_newton": iters,
                            "n_gdl": int(3 * V.dofmap.index_map.size_local),
                            **info})
