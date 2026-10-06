"""
p4_cota_error.py: Cota del error del iterativo con kappa estimado en el CG.

    python comparativa_motores/libros/p4_cota_error.py [n ...]

RESULTADO NEGATIVO: la cota no detecta el error (ver el informe). Por eso
el CG con Lanczos NO se integro en la aplicacion: vive en este guion, que
sustituye temporalmente `resolver_scipy` de los motores al importarse.

El CG + pyamg de scikit-fem se detenia en el residuo pedido con un error de
desplazamiento de hasta 0,9 % en la malla TET10 del espinodoide a 32^3
(comparativa_motores/INFORME.md, sec. 6). Aqui se resuelve el mismo caso
(ensayo lineal de la app, fuerza en el techo) con scikit-fem iterativo, que
ahora estima kappa con los coeficientes del CG, y con NGSolve directo como
referencia, y se compara el error medido con la cota eps_P * kappa_est.
Tambien en hex8, donde el iterativo es fiable, para ver que la cota no da
falsas alarmas.

Escribe resultados/p4_cota_error.json.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parents[1]))
sys.path.insert(0, str(AQUI.parent))
from spinpy import fem, motores                                 # noqa: E402
from spinpy.motores import _comun, m_skfem                      # noqa: E402
from spinpy.motores._comun import ErrorMotor, problema_de_malla  # noqa: E402
import casos                                                    # noqa: E402


def gc_lanczos(A, b, M, tol=1e-10, maxiter=2000):
    """CG precondicionado que ademas estima el numero de condicion.

    Los coeficientes alfa y beta del CG definen la matriz tridiagonal de
    Lanczos del operador precondicionado M A; sus autovalores extremos
    aproximan los de M A desde dentro (Quarteroni, Sacco y Saleri 2007, secs.
    4.4 y 5.11), asi que kappa_est no supera el kappa verdadero. Con el
    residuo relativo eps, el error relativo de la solucion queda acotado por
    eps kappa (Quarteroni et al., sec. 4.6.2): parar por residuo no basta si
    el sistema esta mal condicionado, que es lo que paso con el CG + pyamg de
    scikit-fem en las mallas TET10 (comparativa_motores/INFORME.md). Para un
    metodo precondicionado, la cota usa el residuo precondicionado
    ||M r|| / ||M b|| y el kappa de M A.
    Devuelve (x, residuos, kappa_est, residuo_precondicionado_rel).
    """
    x = np.zeros_like(b)
    r = b.copy()
    z = M @ r
    p = z.copy()
    rz = float(r @ z)
    nb = float(np.linalg.norm(b)) or 1.0
    res = [float(np.linalg.norm(r)) / nb]
    nz0 = float(np.linalg.norm(z)) or 1.0
    alfas, betas = [], []
    for _ in range(maxiter):
        if res[-1] < tol:
            break
        Ap = A @ p
        alfa = rz / float(p @ Ap)
        x += alfa * p
        r -= alfa * Ap
        z = M @ r
        rz_n = float(r @ z)
        beta = rz_n / rz
        rz = rz_n
        p = z + beta * p
        alfas.append(alfa)
        betas.append(beta)
        res.append(float(np.linalg.norm(r)) / nb)
    zrel = float(np.linalg.norm(z)) / nz0
    k = len(alfas)
    kappa = None
    if k >= 2:
        d = np.empty(k)
        d[0] = 1 / alfas[0]
        for j in range(1, k):
            d[j] = 1 / alfas[j] + betas[j - 1] / alfas[j - 1]
        e = np.sqrt(np.maximum(betas[:k - 1], 0)) / np.asarray(alfas[:k - 1])
        from scipy.linalg import eigvalsh_tridiagonal
        w = eigvalsh_tridiagonal(d, e)
        if w[0] > 0:
            kappa = float(w[-1] / w[0])
    return x, res, kappa, zrel


def resolver_scipy_lanczos(A, b, B=None, modo="iterativo", tol=1e-10):
    """Sistema disperso SDP con scipy/pyamg (motores sin resolvedor propio).

    directo    SuperLU (`scipy.sparse.linalg.splu`), ordenacion COLAMD.
    iterativo  CG precondicionado con multigrid algebraico de agregacion
               suavizada (pyamg) y los modos rigidos `B` como espacio nulo
               cercano: la receta de `resistencia.ensayo_compresion`.
    Devuelve (x, info) con el residuo relativo ||b - A x|| / ||b||; si no
    alcanza 100 tol se lanza `ErrorMotor` (un iterativo no avisa solo).
    """
    from scipy.sparse.linalg import splu
    A = A.tocsr()
    if modo == "directo":
        x = splu(A.tocsc()).solve(b)
        nombre, it = "SuperLU (scipy)", None
    else:
        import pyamg
        ml = pyamg.smoothed_aggregation_solver(A, B=B, max_coarse=500)
        x, res, kappa, zrel = gc_lanczos(A, b,
                                         ml.aspreconditioner(cycle="V"),
                                         tol=tol, maxiter=2000)
        nombre = "CG + AMG agregacion suavizada (pyamg, modos rigidos)"
        it = len(res) - 1
    nb = np.linalg.norm(b)
    r = float(np.linalg.norm(b - A @ x) / (nb if nb > 0 else 1.0))
    if not np.isfinite(r) or (modo != "directo" and r > 100 * tol):
        raise ErrorMotor(f"El resolvedor no convergio (residuo {r:.1e}).")
    info = {"solver": nombre, "residuo_rel": r, "iteraciones": it}
    if modo != "directo":
        # cota del error relativo con el kappa estimado (optimista: kappa_est
        # aproxima el verdadero desde abajo)
        info["kappa_est"] = kappa
        info["residuo_precond_rel"] = zrel
        info["cota_error_rel"] = None if kappa is None else float(zrel * kappa)
    return x, info


# los motores importan resolver_scipy por nombre: se sustituye en los dos
_comun.resolver_scipy = resolver_scipy_lanczos
m_skfem.resolver_scipy = resolver_scipy_lanczos

ARCHIVO = AQUI / "resultados" / "p4_cota_error.json"


def caso(tipo, n, out):
    BW, sp = casos.espinodoide(n)
    malla = fem.mallar(BW, sp, tipo)
    sol = {}
    for motor, solver in (("ngsolve", "directo"), ("skfem", "iterativo")):
        p = problema_de_malla(malla, 20000.0, 0.3, cargas=[1.0],
                              solver=solver, tol=1e-10)
        t0 = time.perf_counter()
        r = motores.resolver(p, motor)
        sol[motor] = (r["u"], r["meta"], time.perf_counter() - t0)
    u_ref, u = sol["ngsolve"][0], sol["skfem"][0]
    m = sol["skfem"][1]
    d = {"n_gdl": int(u.size), "error_rel_medido":
         float(np.linalg.norm(u - u_ref) / np.linalg.norm(u_ref)),
         "error_max_rel": float(np.abs(u - u_ref).max() / np.abs(u_ref).max()),
         "residuo_rel": m.get("residuo_rel"),
         "residuo_precond_rel": m.get("residuo_precond_rel"),
         "kappa_est": m.get("kappa_est"),
         "cota_error_rel": m.get("cota_error_rel"),
         "iteraciones": m.get("iteraciones"),
         "t_directo_s": sol["ngsolve"][2], "t_iterativo_s": sol["skfem"][2]}
    out[f"{tipo}_{n}"] = d
    print(tipo, n, {k: v for k, v in d.items()}, flush=True)
    ARCHIVO.write_text(json.dumps(out, indent=1))


def main(ns):
    ARCHIVO.parent.mkdir(exist_ok=True)
    out = json.loads(ARCHIVO.read_text()) if ARCHIVO.exists() else {}
    for n in ns:
        caso("hex8", n, out)
        caso("tet10", n, out)


if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]] or [24, 32])
