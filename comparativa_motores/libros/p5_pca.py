"""
p5_pca.py: Sensibilidad de los ejes de `voi.marco_pca` cuando dos varianzas
casi coinciden.

    python comparativa_motores/libros/p5_pca.py

Elipsoides de semiejes (1, b, 0,5) con b de 0,6 a 0,99, a 80^3. Para cada
uno se compara la sensibilidad estimada por `voi.estabilidad_pca` (grados por
1 % de perturbacion de la covarianza) con el giro del eje principal medido al
perturbar la segmentacion con ruido gaussiano (cinco semillas). Escribe
resultados/p5_pca.json.
"""

import json
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parents[1]))
from spinpy.voi import estabilidad_pca, marco_pca               # noqa: E402


def elipsoide(b, ruido=0.0, seed=0, n=80):
    x, y, z = np.mgrid[-1:1:n * 1j, -1:1:n * 1j, -1:1:n * 1j]
    f = x ** 2 + (y / b) ** 2 + (z / 0.5) ** 2
    if ruido:
        f = f + ruido * np.random.default_rng(seed).standard_normal(f.shape)
    return f < 0.64


def main():
    filas = []
    for b in (0.6, 0.8, 0.9, 0.95, 0.97, 0.99):
        e0, _, p0 = marco_pca(elipsoide(b), 1.0)
        giros = []
        for s in range(1, 6):
            e, _, _ = marco_pca(elipsoide(b, 0.15, s), 1.0)
            giros.append(float(np.degrees(np.arccos(
                min(1.0, abs(e[:, 0] @ e0[:, 0]))))))
        est = estabilidad_pca(p0)
        filas.append({"b": b, "cociente_s2_s1": est["cocientes"][0],
                      "grados_por_pct_eje1": est["grados_por_pct"][0],
                      "giro_medido_max": max(giros),
                      "aviso_eje1": 1 in est["eje_inestable"]})
        print(filas[-1])
    (AQUI / "resultados").mkdir(exist_ok=True)
    (AQUI / "resultados" / "p5_pca.json").write_text(json.dumps(filas,
                                                                indent=1))


if __name__ == "__main__":
    main()
