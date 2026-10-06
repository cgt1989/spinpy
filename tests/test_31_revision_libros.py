"""
Bloque 31 — Cambios propuestos por la revision de la biblioteca de elementos
finitos y analisis numerico (docs/revision_libros/).

LO QUE SE VERIFICA
  (1) R6. Hexaedro con modos incompatibles ('hex8i'): seis modos rigidos
      exactos, prueba de la parcela con deformacion constante exacta, y en
      un voladizo de voxeles con 1 y 2 elementos de canto la flecha queda
      cerca de la de referencia, donde el trilineal es demasiado rigido.
  (2) R1. Homogeneizacion de un espinodoide isotropo de rho 0,25 a 16^3 con
      el iterativo: con los modos rigidos como espacio casi nulo converge y
      coincide con el LU; con el espacio escalar de la V2.1.1, no converge.
  (3) R1. Respaldo LU: a rho 0,15 el iterativo no converge con ningun
      espacio casi nulo y la funcion devuelve el tensor del LU.
  (4) R2. Newton globalizado en NGSolve: el bloque macizo con plato da la
      fuerza de la solucion cerrada; en el espinodoide a 24^3 con
      neo-Hookeano y 5 MPa en un solo incremento converge, y el Newton de la
      V2.1.1 (globalizar=False) no.
  (5) R3. El material no lineal por omision es el neo-Hookeano.
  (7) R5. `voi.estabilidad_pca` avisa con un elipsoide casi de revolucion y
      no con uno de ejes bien separados.
  (8) R7. Hexaedros Q2 de NGSolve ('orden_hex' = 2): en el bloque macizo con
      plato, E_app = E_s.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  (1) modos rigidos: |autovalor| < 1e-10; parcela 1e-12 relativa; voladizo:
      hex8i a menos del 2 % de la referencia con 1 y 2 elementos de canto
      (medido -0,74 % y -0,63 %); hex8 a mas del 10 % (medido -35 % y -12 %).
  (2) residuo <= 1e-6 (criterio de la app) y tensor a 1e-6 del LU.
  (4) bloque: 1e-9 de la solucion cerrada.
  (8) 1e-8 relativa.
"""
from __future__ import annotations

import inspect

import numpy as np
import pytest
from scipy import sparse
from scipy.sparse.linalg import spsolve

from spinpy import elastic, fem, generar_mascara, motores
from spinpy.elastic import hex8_ke

BLOQUE = "31 Revision de la biblioteca"
REF = "docs/revision_libros/"
NG = "ngsolve" in motores.disponibles()


def _anotar(registro, prueba, esp, obt, err, crit, ok, nota=""):
    registro.anotar(BLOQUE, prueba, REF, esp, obt, err, crit, bool(ok),
                    nota=nota)


def _voladizo(nt, elemento, L=10.0, T=1.0, E=1000.0, nu=0.3, P=1.0):
    h = T / nt
    nx = int(round(L / h))
    ke = hex8_ke(h, h, h, E, nu, elemento)
    idx = np.arange((nx + 1) * (nt + 1) ** 2).reshape(nx + 1, nt + 1, nt + 1)
    i, j, k = [a.ravel() for a in np.meshgrid(np.arange(nx), np.arange(nt),
                                              np.arange(nt), indexing="ij")]
    nod = np.stack([idx[i, j, k], idx[i + 1, j, k], idx[i + 1, j + 1, k],
                    idx[i, j + 1, k], idx[i, j, k + 1], idx[i + 1, j, k + 1],
                    idx[i + 1, j + 1, k + 1], idx[i, j + 1, k + 1]], 1)
    edof = (3 * nod[:, :, None] + np.arange(3)).reshape(len(nod), 24)
    nd = 3 * idx.size
    K = sparse.coo_matrix((np.tile(ke.ravel(), len(edof)),
                           (np.repeat(edof, 24, 1).ravel(),
                            np.tile(edof, (1, 24)).ravel())),
                          shape=(nd, nd)).tocsr()
    w = np.ones((nt + 1, nt + 1))
    w[0] *= 0.5; w[-1] *= 0.5; w[:, 0] *= 0.5; w[:, -1] *= 0.5
    w = (w / w.sum()).ravel()
    f = np.zeros(nd)
    f[3 * idx[-1].ravel() + 2] = -P * w
    fijo = np.concatenate([3 * idx[0].ravel() + c for c in range(3)])
    lib = np.setdiff1d(np.arange(nd), fijo)
    u = np.zeros(nd)
    u[lib] = spsolve(K[lib][:, lib].tocsc(), f[lib])
    return float(-(u[3 * idx[-1].ravel() + 2] * w).sum())


def test_hex8i_modos_rigidos_y_parcela(registro):
    k = hex8_ke(1.0, 1.0, 1.0, 1.0, 0.3, "hex8i")
    w = np.linalg.eigvalsh(k)
    n0 = int(np.sum(np.abs(w) < 1e-10))
    X = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    k8 = hex8_ke(1.0, 1.0, 1.0, 1.0, 0.3, "hex8")
    peor = 0.0
    for c in range(6):
        e = np.zeros(6); e[c] = 1.0
        Eps = np.array([[e[0], e[5] / 2, e[4] / 2], [e[5] / 2, e[1], e[3] / 2],
                        [e[4] / 2, e[3] / 2, e[2]]])
        u = (X @ Eps.T).ravel()
        peor = max(peor, abs(u @ k @ u - u @ k8 @ u) / (u @ k8 @ u))
    ok = n0 == 6 and w[6] > 1e-6 and peor < 1e-12
    _anotar(registro, "hex8i: 6 modos rigidos y parcela exacta", 6.0,
            float(n0), peor, "1e-12", ok)
    assert ok, (n0, peor)


def test_hex8i_voladizo(registro):
    ref = 4.00227          # Richardson con hex8i a 6 y 8 elementos de canto
    filas = []
    for nt in (1, 2):
        e8 = _voladizo(nt, "hex8") / ref - 1
        ei = _voladizo(nt, "hex8i") / ref - 1
        filas.append((nt, e8, ei))
        ok = abs(ei) < 0.02 and abs(e8) > 0.10
        _anotar(registro, f"voladizo, {nt} elemento(s) de canto: hex8i",
                ref, ref * (1 + ei), abs(ei), "2 %", ok,
                nota=f"hex8 {e8:+.1%}")
    assert all(abs(ei) < 0.02 and abs(e8) > 0.10 for _, e8, ei in filas), filas


def _espinodoide_iso(rho, n=16):
    BW, _, _ = generar_mascara(resolution=n, rho=rho, thetas=(90, 90, 90),
                               wave_number=10 * np.pi, num_waves=700, seed=1)
    return BW


@pytest.mark.lento
def test_homogeneizacion_espacio_nulo(registro, monkeypatch):
    BW = _espinodoide_iso(0.25)
    monkeypatch.setattr(elastic, "UMBRAL_DIRECTO", 10 ** 9)
    C_lu, _ = elastic.homogeneizar(BW, 1.0, 0.3)
    monkeypatch.setattr(elastic, "UMBRAL_DIRECTO", 0)
    monkeypatch.setattr(elastic, "UMBRAL_RESPALDO_LU", 0)
    C, inf = elastic.homogeneizar(BW, 1.0, 0.3, espacio_nulo="rigidos")
    _, inf_esc = elastic.homogeneizar(BW, 1.0, 0.3, espacio_nulo="escalar")
    dif = float(np.abs(C - C_lu).max() / np.abs(C_lu).max())
    ok = inf["ok"] and dif < 1e-6 and not inf_esc["ok"]
    _anotar(registro, "rho 0,25 a 16^3: CG con modos rigidos frente a LU",
            0.0, dif, dif, "1e-6", ok,
            nota=f"escalar: {inf_esc.get('msg', '')[:60]}")
    assert ok, (inf, dif, inf_esc.get("msg"))


@pytest.mark.lento
def test_homogeneizacion_respaldo_lu(registro, monkeypatch):
    BW = _espinodoide_iso(0.15)
    monkeypatch.setattr(elastic, "UMBRAL_DIRECTO", 10 ** 9)
    C_lu, _ = elastic.homogeneizar(BW, 1.0, 0.3)
    monkeypatch.setattr(elastic, "UMBRAL_DIRECTO", 0)
    C, inf = elastic.homogeneizar(BW, 1.0, 0.3)
    dif = float(np.abs(C - C_lu).max() / np.abs(C_lu).max())
    ok = inf["ok"] and "respaldo" in inf["solver"] and dif < 1e-10
    _anotar(registro, "rho 0,15: respaldo LU cuando el CG no converge",
            0.0, dif, dif, "1e-10", ok, nota=inf["solver"][:60])
    assert ok, inf


def _nl(malla, control, cargas, material, globalizar):
    p = motores.problema_de_malla(malla, E=20000.0, nu=0.3, analisis="nl",
                                  control=control, cargas=cargas,
                                  material=material, solver="directo")
    p["meta"]["globalizar"] = globalizar
    return p, motores.resolver(p, "ngsolve")


@pytest.mark.skipif(not NG, reason="NGSolve no instalado")
def test_newton_globalizado_bloque(registro):
    malla = fem.mallar(np.ones((4, 4, 6), bool), np.full(3, 0.05), "hex8")
    lmb, mu = 20000 * 0.3 / (1.3 * 0.4), 20000 / 2.6
    peor = 0.0
    for e in (0.05, 0.10, 0.20):
        p, out = _nl(malla, "plato", [e], "neohookeano", True)
        lz = 1 - e
        from scipy.optimize import brentq
        a = brentq(lambda a: mu * (a * a - 1) + lmb * np.log(a * a * lz),
                   0.5, 2.0, xtol=1e-15)
        P = lz * (mu * (1 - lz ** -2) + lmb * np.log(a * a * lz) / lz ** 2)
        F = out["F_reac"][-1] / p["meta"]["A_bruta"]
        peor = max(peor, abs(F + P) / abs(P))
    ok = peor < 1e-9
    _anotar(registro, "bloque neo-Hookeano con plato, Newton globalizado",
            0.0, peor, peor, "1e-9", ok)
    assert ok, peor


@pytest.mark.skipif(not NG, reason="NGSolve no instalado")
def test_newton_globalizado_espinodoide(registro):
    BW, _, _ = generar_mascara(resolution=24, wave_number=12 * np.pi,
                               num_waves=700, thetas=(30, 30, 90), rho=0.30,
                               seed=1)
    malla = fem.mallar(BW, np.full(3, 5.0 / 24), "hex8")
    _, out = _nl(malla, "fuerza", [5.0], "neohookeano", True)
    conv = out["meta"]["newton"]
    try:
        _nl(malla, "fuerza", [5.0], "neohookeano", False)
        v211 = True
    except motores._comun.ErrorMotor:
        v211 = False
    ok = not v211
    _anotar(registro, "espinodoide 24^3, neo-Hookeano, 5 MPa, un incremento",
            1.0, 1.0, None, "converge; V2.1.1 no", ok, nota=str(conv))
    assert ok


def test_material_por_omision(registro):
    d = inspect.signature(fem.ensayo).parameters["material"].default
    d2 = inspect.signature(fem.analizar).parameters["material"].default
    ok = d == d2 == "neohookeano"
    _anotar(registro, "material no lineal por omision", None, None, None,
            "neohookeano", ok, nota=f"{d}, {d2}")
    assert ok


def test_estabilidad_pca(registro):
    from spinpy.voi import estabilidad_pca, marco_pca

    def elipsoide(b, n=48):
        x, y, z = np.mgrid[-1:1:n * 1j, -1:1:n * 1j, -1:1:n * 1j]
        return x ** 2 + (y / b) ** 2 + (z / 0.5) ** 2 < 0.64

    a = estabilidad_pca(marco_pca(elipsoide(0.6), 1.0)[2])
    b = estabilidad_pca(marco_pca(elipsoide(0.99), 1.0)[2])
    ok = (not a["aviso"]) and b["aviso"] and 1 in b["eje_inestable"]
    _anotar(registro, "aviso de ejes PCA degenerados", None,
            b["grados_por_pct"][0], None, "aviso solo con b = 0,99", ok)
    assert ok, (a, b)


@pytest.mark.skipif(not NG, reason="NGSolve no instalado")
def test_hex_orden2_bloque(registro):
    malla = fem.mallar(np.ones((4, 4, 6), bool), np.full(3, 0.25), "hex8")
    p = motores.problema_de_malla(malla, E=20000.0, nu=0.3, control="plato",
                                  cargas=[1e-3], solver="directo")
    p["meta"]["orden_hex"] = 2
    out = motores.resolver(p, "ngsolve")
    E = abs(out["F_reac"][-1]) / (p["meta"]["A_bruta"] * 1e-3)
    err = abs(E / 20000.0 - 1)
    ok = err < 1e-8
    _anotar(registro, "hexaedros Q2: bloque macizo con plato, E_app = E_s",
            20000.0, E, err, "1e-8", ok)
    assert ok, E
