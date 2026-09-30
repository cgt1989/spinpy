"""giroide.py: Giroide esqueletico (BV/TV 0,28) con superficie analitica, para el estudio de hex_mejorado."""
import numpy as np

def g(x, y, z):
    return np.sin(x)*np.cos(y) + np.sin(y)*np.cos(z) + np.sin(z)*np.cos(x)
# nivel para BV/TV 0.28 con muestreo fino
r = np.random.default_rng(0).uniform(0, 2*np.pi, (3, 2_000_000))
T = float(np.quantile(g(*r), 1 - 0.28))
def campo(n, celdas):
    """(f, grad, T) en coordenadas de malla (x = 0 .. 1, h = 1/n)."""
    w = 2*np.pi*celdas
    def f(P):
        a = P*w + 0.3
        return g(a[:, 0], a[:, 1], a[:, 2])
    def gr(P):
        x, y, z = (P*w + 0.3).T
        return w*np.stack([np.cos(x)*np.cos(y) - np.sin(z)*np.sin(x),
                           -np.sin(x)*np.sin(y) + np.cos(y)*np.cos(z),
                           -np.sin(y)*np.sin(z) + np.cos(z)*np.cos(x)], 1)
    return f, gr, T
def mascara(n, celdas):
    c = (np.arange(n) + 0.5) / n * 2*np.pi*celdas + 0.3   # desfase: sin simetria con el borde
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    return g(X, Y, Z) > T
if __name__ == "__main__":
    from spinpy.espesor import espesor_local, estadisticas
    print("T", T)
    for cel in (2, 3, 4):
        BW = mascara(32, cel)
        e = estadisticas(espesor_local(BW, 1.0), BW)
        print(cel, BW.mean(), round(e["media"], 2))
