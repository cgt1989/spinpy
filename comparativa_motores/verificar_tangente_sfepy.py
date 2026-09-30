"""verificar_tangente_sfepy.py — Tangente de dw_tl_he_svk frente a diferencias finitas de su residuo.

Ejecutar desde comparativa_motores/ despues de crear /tmp/claude-0/cm/c.npz con un caso
no lineal (bloque 4x4x6). Resultado medido: error relativo 0.37 con dw_tl_he_svk y 2e-11
con dw_lin_elastic (control).
"""
import sys, numpy as np
sys.path.insert(0,'.'); sys.argv=['x']
from motores import mot_sfepy as ms
from motores._base import leer_problema, Cronometro
p = leer_problema('/tmp/claude-0/cm/c.npz'); m = p['meta']
dom, omega, top, field, u, v, ebcs = ms._problema(p, Cronometro())
mat = ms.Material("m", D=ms.stiffness_from_youngpoisson(3, m["E"], m["nu"]))
iv = ms.Integral("iv", order=2)
for nombre in ["dw_tl_he_svk(m.D, v, u)", "dw_lin_elastic(m.D, v, u)"]:
    term = ms.Term.new(nombre, iv, omega, m=mat, v=v, u=u)
    pb = ms.Problem("nl", equations=ms.Equations([ms.Equation("eq", term)]))
    pb.time_update(ebcs=ms.Conditions([]))
    pb.update_materials(); pb.set_default_state()
    rng = np.random.default_rng(0)
    x = 0.02 * m['H'] * rng.standard_normal(field.n_nod*3)
    ev = pb.get_evaluator()
    K = ev.eval_tangent_matrix(x, mtx=pb.equations.create_matrix_graph()).toarray()
    h = 1e-7 * m['H']; dxv = rng.standard_normal(x.size)
    r1 = ev.eval_residual(x + h*dxv); r0 = ev.eval_residual(x - h*dxv)
    fd = (r1 - r0)/(2*h)
    print(nombre, 'error relativo tangente vs DF:', np.linalg.norm(K@dxv - fd)/np.linalg.norm(fd))
