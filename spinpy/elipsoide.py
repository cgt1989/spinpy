"""
elipsoide.py — Ellipsoid Factor (EF): placas frente a barras, punto a punto.

Metodo de Doube (2015), el que implementa BoneJ y el que usan Vafaeefar et al.
(2022) para comparar hueso trabecular con gyroides, espinodoides y
dual-lattices:

    en cada punto x de la estructura se busca el MAYOR ELIPSOIDE (en volumen)
    que CONTIENE a x y cabe entero dentro del solido. Con sus semiejes
    ordenados a <= b <= c,

        EF(x) = a/b - b/c                      EF en [-1, +1]

    placa  (disco: a << b = c)       EF -> -1
    esfera (a = b = c)               EF =  0
    barra  (aguja: a = b << c)       EF -> +1

POR QUE HACE FALTA SI YA ESTA EL SMI
-----------------------------------
El SMI (`morphometry.indice_smi`) responde a la misma pregunta —placa o
barra— pero por un camino indirecto: dilata la superficie y mira como cambia
el area. Salmon et al. (2015) mostraron que ese camino esta CONFUNDIDO POR LA
CONCAVIDAD: las zonas concavas restan area al dilatar y el indice baja sin que
la estructura sea mas laminar. En este proyecto eso no es teorico: el VOI
proximal de H4 tiene un 11.7 % de area concava frente al 1.1 % del espinodoide
ajustado (`Estudio_Discriminadores`), asi que su SMI de 0.77 puede estar
empujado hacia "placa" por la concavidad y no por las placas.

El EF no mira la superficie: mide la FORMA de los huecos que el solido deja
llenar. No tiene termino de concavidad y es LOCAL —da un valor por voxel—, de
modo que ademas de la media da la distribucion: cuanto del hueso es placa y
cuanto es barra, no solo hacia donde tiende el promedio.

COMO SE BUSCA EL ELIPSOIDE MAXIMO
---------------------------------
El problema no es convexo (maximizar un volumen evitando un conjunto de
puntos), asi que toda implementacion es heuristica. BoneJ dilata el elipsoide
en pasos pequenos, detecta contactos muestreando su superficie y lo gira y lo
"sacude" al azar. Aqui se hace lo mismo con dos diferencias que lo hacen mas
determinista:

  1. SEMILLAS: los maximos locales de la transformada de distancia (la cresta
     del eje medial), uno por bloque de lado ~2 radios, que es donde viven los
     elipsoides grandes. Tras la primera ronda se siembran tambien los voxeles
     que ningun elipsoide cubre, con la condicion de que el elipsoide que
     crezca desde ellos los siga conteniendo.

  2. ESCALA EXACTA: fijados el centro c, la orientacion R y la FORMA (las
     proporciones a0 : a1 : a2), el mayor elipsoide admisible tiene formula
     cerrada. Un punto p, en ejes del elipsoide q = R^T (p - c), queda dentro
     del elipsoide de escala s si  sum (q_j / (s a_j))^2 < 1, asi que

         s^2 = min_p  sum_j (q_j / a_j)^2

     es la escala que toca el primer punto sin tragarselo. El volumen,
     (4/3) pi s^3 a0 a1 a2, queda como funcion SIN RESTRICCIONES de ocho
     numeros —centro (3), giro (3) y dos proporciones— y se maximiza con
     Nelder-Mead.

     La primera version hacia otra cosa, y conviene dejar escrito por que no
     funciona: crecer los semiejes de uno en uno desde la esfera inscrita.
     Una esfera que ya toca un punto en una direccion generica se lo traga en
     cuanto crece CUALQUIER eje, porque el elipsoide nuevo contiene al viejo.
     Medido: el cilindro oblicuo se quedaba en su esfera de partida (3.74 x
     3.74 x 3.74) y la esfera de radio 12 en la de 11.18, con el centro
     clavado medio voxel fuera de sitio. Cambiar la forma a escala libre es lo
     que BoneJ consigue contrayendo antes de dilatar; aqui sale de la formula.

  3. ORIENTACION: la inicial sale del analisis de componentes principales
     del solido alrededor de la semilla (en una placa el eje de menor
     varianza es la normal; en una barra el de mayor es su eje), y el giro es
     variable del optimizador. Sin esto un cilindro oblicuo a los ejes daria
     EF proximo a 0: es exactamente lo que comprueba el bloque 11.

QUE SIGNIFICA "CABE DENTRO": LA DISCRETIZACION
----------------------------------------------
Como en `espesor.py`, se trabaja con CENTROS de voxel: un elipsoide es
admisible si no contiene el centro de ningun voxel de fondo. Durante la
optimizacion solo se miran los voxeles de fondo contiguos al solido mas la capa
exterior entera (en un arbol k-d); al final cada elipsoide se valida contra
TODOS los voxeles de su caja y, si se ha colado alguno, se encoge lo justo. Lo
que se reporta nunca contiene fondo.

EL BORDE DEL VOI
----------------
Fuera del array no se sabe que hay. Se trata como FONDO (una capa de vacio
alrededor), igual que BoneJ trata el borde de la imagen: un elipsoide no puede
salirse del cubo. Consecuencia, que conviene tener presente al leer el numero:
una barra que atraviesa todo el VOI queda acotada por su longitud dentro del
cubo, y su EF es  1 - radio/semilongitud  en vez de 1. Como VOI y candidato se
miden en cubos del mismo lado fisico, el sesgo es comun y se cancela al
comparar (la misma regla de la correccion C4).

COSTE, Y POR QUE EL TOPE DE EVALUACIONES NO SE BAJA A LA LIGERA
---------------------------------------------------------------
Es la medida mas cara del proyecto. Medido sobre un espinodoide columnar de
48^3 (rho 0.35, 41 774 voxeles de hueso, 1317 elipsoides):

    evaluaciones  tiempo   EF mediana  frac. placa  relleno
    300 (sin cache) 409 s   0.401       0.169        0.987
    300             373 s   0.401       0.169        0.987
    150             175 s   0.412       0.147        0.986
     80              99 s   0.404       0.133        0.994

Las dos filas de 300 son identicas: la cache de vecinos no cambia nada, solo
ahorra consultas. Bajar el tope SI cambia: la mediana apenas se mueve, pero la
fraccion de placa cae de 0.169 a 0.133. Las placas necesitan mas iteraciones
para aplanarse que las barras para alargarse, asi que un tope bajo sesga HACIA
BARRA, que es justamente la direccion de la pregunta que esta medida sirve para
contestar. Por eso el valor por omision se queda en 300.

Con MENOS SEMILLAS pasa lo mismo. Mismo espinodoide, 300 evaluaciones,
cambiando el lado del bloque que aporta una semilla (`bloque_radios`):

    bloque_radios  tiempo   elipsoides  frac. placa  frac. barra  relleno
    2 (omision)     373 s     1317        0.169        0.664       0.987
    3               362 s     1100        0.171        0.676       0.984
    4               252 s      801        0.152        0.703       0.971

A 3 no se ahorra casi nada; a 4 se ahorra un tercio y la fraccion de placa
vuelve a caer. Las placas grandes necesitan semillas en su plano medio para que
algun elipsoide las encuentre enteras. No hay atajo barato sin sesgo.

REPRODUCIBILIDAD
----------------
Nelder-Mead es determinista y el unico azar —el desempate entre semillas de
igual radio y el submuestreo si hay demasiadas— usa la semilla del proyecto,
20260720: la misma mascara da siempre el mismo EF.

Referencias:
  Doube M. The ellipsoid factor for quantification of rods, plates, and
    intermediate forms in 3D geometries. Front Endocrinol 2015;6:15.
  Salmon PL et al. Structure Model Index does not measure rods and plates in
    trabecular bone. Front Endocrinol 2015;6:162.
  Vafaeefar M, Moerman KM, Kavousi M, Vaughan TJ. A morphological, topological
    and mechanical investigation of gyroid, spinodoid and dual-lattice
    algorithms as structural models of trabecular bone. JMBBM 2022.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

SEMILLA = 20260720

# Margen para que un punto justo en la superficie del elipsoide no cuente como
# dentro por redondeo.
_HOLGURA = 1.0 - 1e-9


class _Buscador:
    """Estado compartido por todas las semillas de una mascara."""

    def __init__(self, BW, spacing, max_evaluaciones, rng):
        self.BW = BW
        self.sp = spacing
        self.n = np.array(BW.shape)
        self.h = float(spacing.min())
        self.rng = rng
        self.max_eval = int(max_evaluaciones)

        # Capa de vacio alrededor: el exterior del VOI cuenta como fondo. Se
        # mete ENTERA en el arbol, no solo donde toca al solido: si no, un
        # elipsoide alargado podria salirse del cubo por una zona de la capa
        # que no linda con hueso.
        self.P = np.pad(BW, 1, mode="constant", constant_values=False)
        frontera = (ndimage.binary_dilation(
            self.P, structure=np.ones((3, 3, 3), bool)) & ~self.P)
        frontera[[0, -1], :, :] = True
        frontera[:, [0, -1], :] = True
        frontera[:, :, [0, -1]] = True
        self.pts = (np.argwhere(frontera) - 1).astype(float) * spacing
        self.arbol = cKDTree(self.pts)
        self.cota = float(np.linalg.norm((self.n + 2) * spacing))

    # -- utilidades --------------------------------------------------------

    def vecinos(self, c, radio):
        ids = self.arbol.query_ball_point(c, radio)
        return self.pts[ids] if len(ids) else np.empty((0, 3))

    def escala(self, c, R, forma, cache):
        """Escala s del mayor elipsoide c, R, s*forma sin puntos dentro.

        Solo cuentan los puntos a menos de `rq` del centro; si el elipsoide
        resultante llega mas lejos, un punto de fuera podria estar dentro, y
        se repite con la bola ampliada.

        LA CACHE, Y POR QUE NO CAMBIA EL RESULTADO. La primera version
        consultaba el arbol k-d en CADA evaluacion de Nelder-Mead. Medido sobre
        un espinodoide de 48^3: 398 000 consultas, 214 de los 409 s. Aqui se
        consulta una vez por semilla, alrededor de su centro inicial y con un
        radio `Rc` holgado, y en cada evaluacion se usan esos puntos. Mientras
        la bola (c, rq) quepa dentro de la cacheada —desplazamiento + rq <= Rc—
        los puntos cacheados CONTIENEN a los de la bola, y tomar el minimo
        sobre un superconjunto de puntos reales del fondo solo puede anadir
        restricciones que tambien existen: la escala es la misma. Si la bola
        se sale, se reconsulta con un radio mayor.
        """
        amax = float(forma.max())
        rq = cache["rq"]
        while True:
            desp = float(np.linalg.norm(c - cache["c"]))
            if desp + rq > cache["Rc"] and cache["Rc"] < self.cota:
                cache["Rc"] = min(max(2.0 * cache["Rc"], 1.25 * (desp + rq)),
                                  self.cota)
                cache["G"] = self.vecinos(cache["c"], cache["Rc"])
            G = cache["G"]
            if len(G):
                v = ((((G - c) @ R) / forma) ** 2).sum(axis=1)
                s = float(np.sqrt(v.min()))
            else:
                s = np.inf
            if s * amax <= rq or rq >= self.cota:
                cache["rq"] = rq
                return min(s, self.cota / amax)
            rq = min(max(2.0 * rq, 1.01 * s * amax), self.cota)

    def orientacion_inicial(self, idx, r0):
        """Ejes principales del solido alrededor de la semilla."""
        w = int(np.ceil(2.0 * r0 / self.h)) + 1
        lo = np.maximum(idx - w, 0)
        hi = np.minimum(idx + w + 1, self.n)
        sub = self.BW[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
        X = np.argwhere(sub).astype(float) * self.sp
        if len(X) < 5:
            return np.eye(3)
        _, V = np.linalg.eigh(np.cov(X.T))
        if np.linalg.det(V) < 0:
            V[:, 0] = -V[:, 0]
        return V

    def desde_semilla(self, idx, r0, anclada):
        """Maximiza el volumen sobre centro, giro y forma (8 variables).

        x = [desplazamiento del centro / h (3), vector de giro (3),
             log(a1/a0), log(a2/a0)]
        """
        c0 = idx.astype(float) * self.sp
        R0 = self.orientacion_inicial(idx, r0)
        rq0 = max(2.0 * r0, 3.0 * self.h)
        cache = {"c": c0, "rq": rq0, "Rc": 2.0 * rq0}
        cache["G"] = self.vecinos(c0, cache["Rc"])

        def desplegar(x):
            c = c0 + x[:3] * self.h
            R = Rotation.from_rotvec(x[3:6]).as_matrix() @ R0
            forma = np.exp(np.array([0.0, x[6], x[7]]))
            return c, R, forma

        def objetivo(x):
            c, R, forma = desplegar(x)
            s = self.escala(c, R, forma, cache)
            if not np.isfinite(s) or s <= 0.0:
                return 1e3
            if anclada:
                q = (c0 - c) @ R
                if float(((q / (s * forma)) ** 2).sum()) > 1.0:
                    return 1e3            # ya no contiene su semilla
            return -(3.0 * np.log(s) + x[6] + x[7])

        x0 = np.zeros(8)
        # Simplex inicial con pasos de la escala de cada variable: medio
        # voxel en el centro, ~17 grados en el giro, un factor 2 en la forma.
        simplex = np.vstack([x0] + [x0 + d * np.eye(8)[i] for i, d in
                                    enumerate([0.5] * 3 + [0.3] * 3
                                              + [0.7] * 2)])
        r = minimize(objetivo, x0, method="Nelder-Mead",
                     options={"initial_simplex": simplex,
                              "maxfev": self.max_eval, "xatol": 1e-3,
                              "fatol": 1e-4})
        c, R, forma = desplegar(r.x)
        s = self.escala(c, R, forma, cache)
        if not np.isfinite(s) or r.fun >= 1e3:
            return c0, np.eye(3), np.full(3, r0 * _HOLGURA)
        return c, R, forma * (s * _HOLGURA)

    def validar_y_cubrir(self, c, R, ax):
        """Valida contra TODOS los voxeles de su caja y devuelve los cubiertos.

        Encoge el elipsoide lo justo si contiene algun centro de fondo —de la
        mascara o de la capa exterior— o si se sale de la caja con capa.
        """
        lo_lim = -self.sp
        hi_lim = self.n * self.sp
        ext = np.sqrt(((R * ax[None, :]) ** 2).sum(axis=1))
        s = float(np.min(np.minimum((c - lo_lim) / ext, (hi_lim - c) / ext)))
        if s < 1.0:
            ax = ax * (s * _HOLGURA)
            ext = ext * (s * _HOLGURA)

        i0 = np.maximum(np.floor((c - ext) / self.sp).astype(int), -1)
        i1 = np.minimum(np.ceil((c + ext) / self.sp).astype(int), self.n)
        rng_ = [np.arange(i0[k], i1[k] + 1) for k in range(3)]
        I, J, K = np.meshgrid(*rng_, indexing="ij")
        I, J, K = I.ravel(), J.ravel(), K.ravel()
        X = np.column_stack([I, J, K]).astype(float) * self.sp
        Q = (X - c) @ R
        v = ((Q / ax) ** 2).sum(axis=1)
        solido = self.P[I + 1, J + 1, K + 1]

        colados = (v < 1.0) & ~solido
        if colados.any():
            f = np.sqrt(float(v[colados].min())) * _HOLGURA
            ax = ax * f
            v = v / (f * f)
        cubre = (v <= 1.0) & solido
        return ax, I[cubre], J[cubre], K[cubre]


def _elegir_semillas(candidatos, r, bloque, max_n, rng):
    """Un voxel por bloque: el de mayor radio inscrito (empates al azar)."""
    idx = np.argwhere(candidatos)
    if len(idx) == 0:
        return idx, np.empty(0)
    rv = r[idx[:, 0], idx[:, 1], idx[:, 2]]
    dims = tuple(-(-np.array(candidatos.shape) // bloque))
    bid = np.ravel_multi_index(tuple((idx // bloque).T), dims)
    orden = np.lexsort((rng.random(len(rv)), -rv))
    _, prim = np.unique(bid[orden], return_index=True)
    sel = orden[prim]
    if len(sel) > max_n:
        sel = rng.choice(sel, max_n, replace=False)
    return idx[sel], rv[sel]


def factor_elipsoide(BW, spacing=1.0, max_evaluaciones=300, max_semillas=3000,
                     relleno_objetivo=0.98, rondas=4, semilla=SEMILLA,
                     campo=False, bloque_radios=2.0):
    """Ellipsoid Factor de una mascara binaria.

    Devuelve un diccionario con el resumen sobre los voxeles cubiertos:

      EF, EF_mediana, EF_sd, EF_p05, EF_p95
      EF_frac_placa    fraccion del solido con EF < -0.2
      EF_frac_barra    fraccion del solido con EF > +0.2
      EF_relleno       fraccion del solido cubierta por algun elipsoide.
                       Vafaeefar et al. exigen > 95 % para dar el EF por
                       valido; por debajo, el numero describe solo una parte.
      EF_n_elipsoides

    Con `campo=True` anade `EF_campo`: array del tamano de BW con el EF de cada
    voxel y NaN en el fondo y en lo no cubierto.

    Los umbrales de placa y barra (+-0.2) no son de Doube; se eligen simetricos
    y fuera de la banda en que una esfera digitalizada ya oscila, para que las
    dos fracciones digan algo mas que el signo del ruido.

    Percentiles con la convencion de `prctile` de MATLAB (method="hazen"),
    como en el resto del proyecto.
    """
    BW = np.asarray(BW, dtype=bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)

    vacio = {"EF": np.nan, "EF_mediana": np.nan, "EF_sd": np.nan,
             "EF_p05": np.nan, "EF_p95": np.nan, "EF_frac_placa": np.nan,
             "EF_frac_barra": np.nan, "EF_relleno": np.nan,
             "EF_n_elipsoides": 0}
    if not BW.any() or BW.all():
        if campo:
            vacio["EF_campo"] = np.full(BW.shape, np.nan, np.float32)
        return vacio

    rng = np.random.default_rng(semilla)
    B = _Buscador(BW, spacing, max_evaluaciones, rng)

    # Radio inscrito con la capa exterior, para que cerca de las caras no se
    # suponga que el solido sigue.
    r = ndimage.distance_transform_edt(B.P, sampling=spacing)[1:-1, 1:-1, 1:-1]
    cresta = BW & (r >= ndimage.maximum_filter(r, size=3) - 1e-12)
    r_med = float(np.median(r[cresta]))
    bloque = max(2, int(np.ceil(float(bloque_radios) * r_med / B.h)))

    vol_asig = np.zeros(BW.shape, np.float64)
    ef_asig = np.full(BW.shape, np.nan, np.float32)
    n_tot = int(BW.sum())
    n_elip = 0

    candidatos = cresta
    for ronda in range(int(rondas)):
        semillas, radios = _elegir_semillas(candidatos, r, bloque,
                                            max_semillas, rng)
        if len(semillas) == 0:
            break
        for idx, r0 in zip(semillas, radios):
            c, R, ax = B.desde_semilla(idx, float(r0), anclada=ronda > 0)
            ax, I, J, K = B.validar_y_cubrir(c, R, ax)
            if len(I) == 0:
                continue
            a, b, cc = np.sort(ax)
            ef = a / b - b / cc
            vol = float(np.prod(ax))
            # El mayor elipsoide que contiene cada voxel gana, sin importar el
            # orden en que se encontraron: se compara volumen, no turno.
            gana = vol > vol_asig[I, J, K]
            vol_asig[I[gana], J[gana], K[gana]] = vol
            ef_asig[I[gana], J[gana], K[gana]] = ef
            n_elip += 1
        cubierto = vol_asig > 0
        if cubierto.sum() / n_tot >= relleno_objetivo:
            break
        candidatos = BW & ~cubierto

    cubierto = vol_asig > 0
    v = ef_asig[cubierto].astype(float)
    out = dict(vacio)
    out["EF_n_elipsoides"] = n_elip
    out["EF_relleno"] = float(cubierto.sum() / n_tot)
    if v.size:
        out.update({
            "EF": float(v.mean()),
            "EF_mediana": float(np.median(v)),
            "EF_sd": float(v.std(ddof=1)) if v.size > 1 else 0.0,
            "EF_p05": float(np.percentile(v, 5, method="hazen")),
            "EF_p95": float(np.percentile(v, 95, method="hazen")),
            "EF_frac_placa": float((v < -0.2).mean()),
            "EF_frac_barra": float((v > 0.2).mean()),
        })
    if campo:
        out["EF_campo"] = ef_asig
    return out
