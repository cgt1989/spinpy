"""
febio.py — FEBio 4 como segundo resolvedor de los ensayos de spinpy.

La app resuelve sus ensayos con su propio codigo (`resistencia`, `elastic`),
sobre hexaedros de un voxel. Este modulo manda el MISMO problema a FEBio
(`febio4.exe`, un proceso externo) con dos mallas posibles:

  hex8    la malla de voxeles que resuelve la app, un elemento por voxel. Es
          el CONTROL: `comparativa_febio/` la valido contra la app a <= 1e-6
          en E_app, u, sigma, p99 de superficie y Pistoia. Lo que aporta es lo
          que la app no tiene (no linealidad geometrica, plato rigido), pero
          conserva los escalones.
  tet10   una malla tetraedrica cuadratica SUAVE (`solido.malla_tet10`): quita
          los escalones. Es lo mas parecido al metodo de Tapia et al. (2026),
          que mallaron SOLID187 de 0.05 mm en ANSYS sobre el STL.

Con las dos sobre la misma estructura y el mismo protocolo se separan los
efectos: app vs FEBio-hex8 mide la IMPLEMENTACION (debe dar ~0); FEBio-hex8 vs
FEBio-tet10 mide solo el efecto de la MALLA, con el mismo programa, material y
carga.

LA GUI SOLO LLAMA A ESTE MODULO, como a `resistencia`: todo corre igual desde
la CLI y desde los tests.

HECHOS MEDIDOS QUE ESTE CODIGO RESPETA (no son supuestos)
---------------------------------------------------------
* FEBio 4.5.0 termina con violacion de acceso (0xC0000005) si `Output` solo
  tiene logfile: se escribe SIEMPRE un plotfile.
* FEBio 4.5 escribe reaccion CERO en los GDL de un `zero displacement`: la
  fuerza se obtiene de la integral de volumen sum(sigma_zz V_e) = -F H, exacta
  en el problema discreto (la base contribuye con z = z_base).
* FEBio es NO LINEAL geometricamente ('isotropic elastic' es St.
  Venant-Kirchhoff): el problema lineal se obtiene extrapolando a carga nula,
  u = 2 u(s) - u(2 s). Con una sola carga pequena queda un desvio de 4e-4 en
  el VOI proximal de H4, por los voladizos del techo.

Unidades de los .feb: mm - N - MPa. Aqui, como en el resto de spinpy, E y
sigma viajan en Pa y se convierten al escribir y al leer.
"""

from __future__ import annotations

import glob
import hashlib
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

# Nombre del ejecutable y carpeta de la copia empaquetada con el instalador.
EJECUTABLE = "febio4.exe" if os.name == "nt" else "febio4"
CARPETA_EMPAQUETADA = "febio"

# Instalaciones de FEBio Studio: la 2.x deja el solver en su carpeta bin.
PATRONES_ESTUDIO = (r"C:\Program Files\FEBioStudio*\bin",
                    r"C:\Program Files (x86)\FEBioStudio*\bin")

# Texto con el que FEBio cierra una corrida correcta.
TERMINACION_NORMAL = "N O R M A L   T E R M I N A T I O N"


class ErrorFEBio(RuntimeError):
    """FEBio no termino bien. `log` lleva la cola del .log para mostrarla."""

    def __init__(self, msg, log=""):
        super().__init__(msg)
        self.log = log


class Cancelado(RuntimeError):
    """El usuario detuvo la corrida; el proceso de FEBio ya se mato."""


# ---------------------------------------------------------------------------
# Localizar el ejecutable
# ---------------------------------------------------------------------------

def _candidatas_empaquetadas():
    """Donde puede estar la copia que lleva el instalador.

    Con PyInstaller los datos van a `sys._MEIPASS` (modo onefile) o junto al
    ejecutable (modo carpeta); desde las fuentes, junto a `Port_Python/`.
    """
    bases = []
    if getattr(sys, "_MEIPASS", None):
        bases.append(Path(sys._MEIPASS))
    if getattr(sys, "frozen", False):
        bases.append(Path(sys.executable).resolve().parent)
    bases.append(Path(__file__).resolve().parents[1])
    return [b / CARPETA_EMPAQUETADA / EJECUTABLE for b in bases]


def localizar(ruta_usuario=None):
    """Ruta a `febio4.exe`, o None si no hay ninguno.

    Orden: (1) la copia EMPAQUETADA con spinpy —es la que se valido con esta
    version—; (2) la ruta que eligio el usuario (la GUI la guarda en
    QSettings y la pasa aqui; este modulo no depende de Qt); (3) la variable
    de entorno SPINPY_FEBIO; (4) una instalacion de FEBio Studio.
    """
    candidatas = list(_candidatas_empaquetadas())
    if ruta_usuario:
        candidatas.append(Path(ruta_usuario))
    if os.environ.get("SPINPY_FEBIO"):
        candidatas.append(Path(os.environ["SPINPY_FEBIO"]))
    for patron in PATRONES_ESTUDIO:
        for d in sorted(glob.glob(patron), reverse=True):
            candidatas.append(Path(d) / EJECUTABLE)
    for c in candidatas:
        if c.is_file():
            return c
    return None


def origen(exe):
    """'empaquetado', 'usuario/entorno' o 'FEBio Studio', para el registro."""
    if exe is None:
        return None
    exe = Path(exe).resolve()
    if any(exe == c.resolve() for c in _candidatas_empaquetadas()
           if c.is_file()):
        return "empaquetado"
    if "febiostudio" in str(exe).lower():
        return "FEBio Studio"
    return "usuario"


_VERSIONES = {}


def version(exe=None):
    """Version de FEBio ('4.5.0'), leida del .log de un modelo minimo.

    `febio4.exe` no tiene opcion de linea de comandos que la imprima (-v y -h
    dan 'Invalid command line option'), asi que se corre un hexaedro solo, en
    una carpeta temporal: ~0.3 s. Se guarda en cache por ruta.
    """
    exe = Path(exe) if exe else localizar()
    if exe is None:
        return None
    clave = str(exe)
    if clave not in _VERSIONES:
        from .escribe import escribir_febio
        nodos = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                          [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
        with tempfile.TemporaryDirectory() as d:
            feb = Path(d) / "version.feb"
            escribir_febio(nodos, np.arange(8)[None, :], feb, sigma_app=1e3)
            try:
                _VERSIONES[clave] = correr(feb, exe=exe)["version"]
            except (ErrorFEBio, OSError):
                _VERSIONES[clave] = None
    return _VERSIONES[clave]


# ---------------------------------------------------------------------------
# Correr y leer
# ---------------------------------------------------------------------------

def _cola(texto, n=40):
    return "\n".join(texto.splitlines()[-n:])


def _duracion_s(txt):
    m = re.search(r"Total elapsed time[ .]*:\s*[\d:]+\s*\(([\d.eE+-]+)\s*sec",
                  txt)
    return float(m.group(1)) if m else None


def correr(feb, exe=None, hilos=None, cancelar=None, intervalo=0.25,
           tiempo_max=None):
    """Corre `febio4 -i <feb> -silent` en la carpeta del .feb.

    `hilos`: numero de hilos de OpenMP/MKL (Pardiso). None deja el defecto de
    FEBio (todos los nucleos). La GUI pasa nucleos - 1 para seguir
    respondiendo.

    `cancelar`: funcion sin argumentos (o `threading.Event`) que se consulta
    cada `intervalo` s; si devuelve True se MATA el proceso y se lanza
    `Cancelado`. Es lo que hace «Detener tras la etapa actual»: FEBio no tiene
    forma de interrumpirse limpiamente desde fuera.

    Devuelve {tiempo_s, version, memoria_MB (pico que reporta FEBio),
    iteraciones, log}. Lanza `ErrorFEBio` con la cola del .log si FEBio no
    llega a la terminacion normal.
    """
    feb = Path(feb)
    exe = Path(exe) if exe else localizar()
    if exe is None:
        raise ErrorFEBio("No se encontro FEBio (febio4.exe).")
    env = dict(os.environ)
    if hilos:
        env["OMP_NUM_THREADS"] = str(int(hilos))
        env["MKL_NUM_THREADS"] = str(int(hilos))
    if hasattr(cancelar, "is_set"):
        cancelar = cancelar.is_set
    log = feb.with_suffix(".log")
    if log.exists():
        log.unlink()

    flags = getattr(subprocess, "CREATE_NO_WINDOW", 0)
    t0 = time.perf_counter()
    p = subprocess.Popen([str(exe), "-i", feb.name, "-silent"], cwd=feb.parent,
                         stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL, env=env,
                         creationflags=flags)
    try:
        while True:
            try:
                p.wait(timeout=intervalo)
                break
            except subprocess.TimeoutExpired:
                pass
            if cancelar is not None and cancelar():
                p.kill()
                p.wait()
                raise Cancelado(f"Corrida de FEBio detenida: {feb.name}")
            if tiempo_max and time.perf_counter() - t0 > tiempo_max:
                p.kill()
                p.wait()
                raise ErrorFEBio(f"FEBio supero {tiempo_max:.0f} s en {feb.name}")
    finally:
        if p.poll() is None:
            p.kill()
    t = time.perf_counter() - t0

    txt = log.read_text(errors="replace") if log.exists() else ""
    if p.returncode != 0 or TERMINACION_NORMAL not in txt:
        raise ErrorFEBio(f"FEBio no termino bien en {feb.name} (codigo "
                         f"{p.returncode}); ver {log}", _cola(txt))
    m = re.search(r"version\s+(\d+\.\d+\.\d+)", txt)
    mem = re.search(r"Peak memory\s*:\s*([\d.]+)\s*MB", txt)
    it = re.search(r"Total number of equilibrium iterations[ .]*:\s*(\d+)", txt)
    return {"tiempo_s": t, "tiempo_febio_s": _duracion_s(txt),
            "version": m.group(1) if m else "?",
            "memoria_MB": float(mem.group(1)) if mem else None,
            "iteraciones": int(it.group(1)) if it else None,
            "log": str(log)}


def ultimo_registro(ruta, ncol):
    """Ultimo bloque de datos de un archivo de logfile de FEBio.

    FEBio escribe un bloque por paso convergido, INCLUIDO el paso 0 (todo
    ceros), con cabeceras '*Step', '*Time', '*Data'. Se toma el ultimo y se
    ordena por el identificador, que es la primera columna. Movido de
    `comparativa_febio/comparar_febio.py` (`_ultimo_registro`), que lo sigue
    teniendo igual.
    """
    bloques = re.split(r"^\*Step.*$", Path(ruta).read_text(), flags=re.M)
    filas = [l.split() for l in bloques[-1].splitlines()
             if l.strip() and not l.startswith("*") and "=" not in l]
    a = np.array(filas, float)
    if a.ndim != 2 or a.shape[1] != ncol + 1:
        raise ErrorFEBio(f"{ruta}: {a.shape} columnas, se esperaban "
                         f"{ncol + 1}")
    return a[np.argsort(a[:, 0]), 1:]


def leer(feb):
    """Desplazamientos y tensiones del ultimo paso de un .feb ya corrido.

    Devuelve {"u": (N, 3) en mm, "sigma": (M, 6) en Pa}, con la tension de
    Cauchy en el orden Voigt del proyecto [xx yy zz yz xz xy] y por
    elemento: la MEDIA de sus puntos de integracion, que es lo que da el
    logfile de FEBio. En un hexaedro rectangular es el valor del centro; en un
    TET10 de lados rectos con la regla de 4 puntos (pesos iguales, exacta para
    tension lineal) es la media de volumen del elemento.
    """
    feb = Path(feb)
    u = ultimo_registro(feb.with_name(feb.stem + "_u.txt"), 3)
    s = ultimo_registro(feb.with_name(feb.stem + "_s.txt"), 6)
    return {"u": u, "sigma": s * 1e6}


# ---------------------------------------------------------------------------
# Magnitudes derivadas del campo de tension
# ---------------------------------------------------------------------------

def von_mises(s):
    """von Mises = sqrt(3 J2) de tensiones Voigt [xx yy zz yz xz xy]."""
    s = np.asarray(s, float)
    return np.sqrt(np.maximum(
        0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2
               + (s[:, 2] - s[:, 0]) ** 2)
        + 3.0 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2), 0.0))


def eps_eff(sig, E, nu):
    """Deformacion efectiva de Pistoia a partir de la tension (Voigt, Pa).

    spinpy la calcula de sigma y eps; aqui eps sale de sigma por la
    flexibilidad isotropa, con distorsiones de ingenieria. Es la misma
    energia: U = sigma . S sigma / 2, eps_eff = sqrt(2 U / E).
    """
    S = np.zeros((6, 6))
    S[:3, :3] = -nu / E
    np.fill_diagonal(S[:3, :3], 1.0 / E)
    S[3, 3] = S[4, 4] = S[5, 5] = 2.0 * (1.0 + nu) / E
    sig = np.asarray(sig, float)
    eps = sig @ S.T
    U = 0.5 * np.einsum("ij,ij->i", sig, eps)
    return np.sqrt(np.maximum(2.0 * U / E, 0.0))


def escribir_prescrito(nodos, elems, u, ruta, E_s=20e9, nu_s=0.30):
    """Malla con el desplazamiento de TODOS sus nodos impuesto (validacion).

    Es la prueba del orden de nodos: con un campo cuadratico exacto, un TET10
    en el orden que FEBio espera reproduce la deformacion lineal exacta; en
    otro orden, FEBio interpreta los nodos intermedios en aristas que no son
    las suyas, la geometria del elemento se curva y la tension sale mal SIN
    ningun error. Un campo lineal NO sirve para detectarlo: un elemento
    isoparametrico reproduce campos lineales con cualquier geometria.
    Pocos nodos: una condicion por nodo y grado de libertad.
    """
    ruta = Path(ruta)
    nodos = np.asarray(nodos, float)
    elems = np.asarray(elems, dtype=np.int64)
    u = np.asarray(u, float).reshape(-1, 3)
    tipo = {8: "hex8", 10: "tet10", 4: "tet4"}[elems.shape[1]]
    regla = ' elem_type="TET10G8"' if tipo == "tet10" else ""
    L = ['<?xml version="1.0" encoding="ISO-8859-1"?>',
         '<febio_spec version="4.0">',
         '\t<Module type="solid">\n\t\t<units>mm-N-s</units>\n\t</Module>',
         "\t<Control>\n\t\t<analysis>STATIC</analysis>\n"
         "\t\t<time_steps>1</time_steps>\n\t\t<step_size>1</step_size>\n"
         "\t\t<solver>\n\t\t\t<max_refs>50</max_refs>\n"
         '\t\t\t<qn_method type="BFGS">\n\t\t\t\t<max_ups>0</max_ups>\n'
         "\t\t\t</qn_method>\n\t\t\t<dtol>1e-9</dtol>\n"
         "\t\t\t<etol>1e-12</etol>\n\t\t\t<rtol>1e-12</rtol>\n"
         '\t\t\t<linear_solver type="pardiso"/>\n\t\t</solver>\n\t</Control>',
         '\t<Material>\n\t\t<material id="1" name="m" type="isotropic elastic">'
         f"\n\t\t\t<E>{E_s / 1e6:.17g}</E>\n\t\t\t<v>{nu_s:.17g}</v>\n"
         "\t\t</material>\n\t</Material>",
         '\t<Mesh>\n\t\t<Nodes name="todos">']
    L += [f'\t\t\t<node id="{i}">{p[0]:.17g},{p[1]:.17g},{p[2]:.17g}</node>'
          for i, p in enumerate(nodos, start=1)]
    L += [f'\t\t</Nodes>\n\t\t<Elements type="{tipo}" name="s">']
    L += [f'\t\t\t<elem id="{i}">' + ",".join(map(str, c)) + "</elem>"
          for i, c in enumerate((elems + 1).tolist(), start=1)]
    L += ["\t\t</Elements>"]
    L += [f'\t\t<NodeSet name="n{i}">{i}</NodeSet>'
          for i in range(1, nodos.shape[0] + 1)]
    L += ["\t</Mesh>",
          f'\t<MeshDomains>\n\t\t<SolidDomain name="s" mat="m"{regla}/>\n'
          "\t</MeshDomains>", "\t<Boundary>"]
    for i in range(1, nodos.shape[0] + 1):
        for c, d in enumerate("xyz"):
            L.append(f'\t\t<bc type="prescribed displacement" node_set="n{i}">'
                     f'<dof>{d}</dof><value lc="1">{u[i - 1, c]:.17g}</value>'
                     "<relative>0</relative></bc>")
    L += ["\t</Boundary>",
          '\t<LoadData>\n\t\t<load_controller id="1" type="loadcurve">\n'
          "\t\t\t<interpolate>LINEAR</interpolate>\n\t\t\t<points>\n"
          "\t\t\t\t<point>0,0</point>\n\t\t\t\t<point>1,1</point>\n"
          "\t\t\t</points>\n\t\t</load_controller>\n\t</LoadData>",
          '\t<Output>\n\t\t<plotfile type="febio">\n'
          '\t\t\t<var type="displacement"/>\n\t\t</plotfile>\n\t\t<logfile>\n'
          f'\t\t\t<node_data data="ux;uy;uz" delim=" " file="{ruta.stem}_u.txt"/>'
          f'\n\t\t\t<element_data data="sx;sy;sz;syz;sxz;sxy" delim=" " '
          f'file="{ruta.stem}_s.txt"/>\n\t\t</logfile>\n\t</Output>',
          "</febio_spec>"]
    ruta.write_text("\n".join(L) + "\n", encoding="utf-8")
    return {"ruta": str(ruta), "tipo": tipo, "n_nodos": int(nodos.shape[0])}


def huella(nodos, elems):
    """SHA-256 de la malla (coordenadas float64 y conectividad int64, en C)."""
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(np.asarray(nodos, np.float64)).tobytes())
    h.update(np.ascontiguousarray(np.asarray(elems, np.int64)).tobytes())
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Mallas
# ---------------------------------------------------------------------------

MALLAS = ("hex8", "tet10")

# Caras de un tetraedro por sus esquinas locales, cada una con la esquina
# opuesta, y aristas C3D10 (ranuras 4..9) para montar las tri6.
_CARAS_TET = ((1, 2, 3, 0), (0, 3, 2, 1), (0, 1, 3, 2), (0, 2, 1, 3))
_ARISTA_C3D10 = {(0, 1): 4, (1, 2): 5, (0, 2): 6, (0, 3): 7, (1, 3): 8,
                 (2, 3): 9}

#: Perdida de volumen de la malla suave, en % del volumen de VOXELES del hueso
#: portante, a partir de la cual el registro va con reservas. DECLARADO ANTES
#: DE MEDIR (2026-09-24), a partir de lo que ya estaba medido en `solido.py`
#: (-5 % frente a voxeles, -1.9 % frente a la superficie cruda, spinodoide a
#: 48^3): la rigidez escala ~ rho^2, asi que un 3 % de volumen son ~6 % de
#: rigidez que no vienen de los escalones. Es NUESTRO criterio.
PERDIDA_VOLUMEN_MAX_PCT = 3.0

#: Hueso descartado al exigir que la malla una base y techo, en % del volumen
#: de la malla, a partir del cual hay reservas. Mismo umbral que el 1 % de
#: `informe.DESCONEXION_RESERVAS`, por coherencia. NUESTRO criterio.
DESCARTE_MAX_PCT = 1.0


def _volumen_tet(nodos, elems):
    P = nodos[elems[:, :4]]
    return np.einsum("ij,ij->i", P[:, 1] - P[:, 0],
                     np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 0])) / 6.0


def caras_borde_tet10(elems):
    """Caras de borde de una malla TET10: (elemento, tri6 saliente).

    Una cara de esquinas que aparece en un solo tetraedro es de borde. Se
    orienta con la normal hacia FUERA del tetraedro (la esquina opuesta queda
    detras), que es lo que necesita la presion de FEBio, y se completa con sus
    nodos intermedios en el orden tri6 (esquinas a, b, c; aristas ab, bc, ca).
    """
    elems = np.asarray(elems, dtype=np.int64)
    n = elems.shape[0]
    tri, dueno = [], []
    for a, b, c, _ in _CARAS_TET:
        tri.append(elems[:, [a, b, c]])
        dueno.append(np.arange(n))
    tri = np.concatenate(tri)
    dueno = np.concatenate(dueno)
    loc = np.concatenate([np.tile(np.array(f[:3]), (n, 1)) for f in _CARAS_TET])
    clave = np.sort(tri, axis=1)
    _, inv, cuenta = np.unique(clave, axis=0, return_inverse=True,
                               return_counts=True)
    borde = cuenta[inv.ravel()] == 1
    tri, dueno, loc = tri[borde], dueno[borde], loc[borde]

    def medio(i, j):
        k = np.array([_ARISTA_C3D10[tuple(sorted((int(p), int(q))))]
                      for p, q in zip(i, j)])
        return elems[dueno, k]

    m_ab = medio(loc[:, 0], loc[:, 1])
    m_bc = medio(loc[:, 1], loc[:, 2])
    m_ca = medio(loc[:, 2], loc[:, 0])
    return dueno, np.column_stack([tri, m_ab, m_bc, m_ca])


def _orientar_saliente(nodos, elems, dueno, caras):
    """Invierte las caras cuya normal apunta hacia dentro de su elemento."""
    P = nodos[caras[:, :3]]
    n = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    c_el = nodos[elems[dueno, :4]].mean(axis=1)
    dentro = np.einsum("ij,ij->i", n, c_el - P[:, 0]) > 0
    caras = caras.copy()
    # (a, b, c, ab, bc, ca) -> (a, c, b, ca, bc, ab)
    caras[dentro] = caras[dentro][:, [0, 2, 1, 5, 4, 3]]
    return caras


def _componentes_portantes(nodos, elems, z0, z1, tol):
    """Elementos de los componentes que tocan base Y techo."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    e4 = elems[:, :4]
    n = nodos.shape[0]
    filas = np.repeat(e4[:, 0], 3)
    cols = e4[:, 1:].ravel()
    G = coo_matrix((np.ones(filas.size), (filas, cols)), shape=(n, n))
    ncomp, lab = connected_components(G, directed=False)
    z = nodos[:, 2]
    en_base = set(np.unique(lab[z <= z0 + tol]))
    en_techo = set(np.unique(lab[z >= z1 - tol]))
    buenos = np.array(sorted(en_base & en_techo), dtype=np.int64)
    usados = np.unique(e4)
    n_comp = int(np.unique(lab[usados]).size)
    return np.isin(lab[e4[:, 0]], buenos), n_comp


def _compactar(nodos, elems):
    usados, inv = np.unique(elems, return_inverse=True)
    return nodos[usados], inv.reshape(elems.shape)


def opciones_malla(**kw):
    """Opciones de la malla suave con sus valores por omision, validadas."""
    from .solido import DECIMADO, SUAVIZADO_BANDA, SUAVIZADO_ITER
    o = {"suavizado": SUAVIZADO_ITER, "banda": SUAVIZADO_BANDA,
         "decimado": DECIMADO, "tam_max_mm": None, "minratio": None,
         "corregir_volumen": True}
    for k, v in kw.items():
        if k not in o:
            raise ValueError(f"opcion de malla desconocida: {k}")
        o[k] = v
    return o


def maxvolume_de_tamano(L):
    """Volumen del tetraedro regular de arista L: L^3 / (6 sqrt 2)."""
    return None if not L else float(L) ** 3 / (6.0 * np.sqrt(2.0))


def preparar_tet10(nodos, elems, forma, spacing):
    """Post-proceso de una malla TET10 dentro del cubo de la mascara.

    El cubo es el de `solido.superficie_cerrada`: planos en -h/2 y
    (n - 1/2) h en cada eje. Se exige que la MALLA una base y techo (se
    descartan los componentes que no toquen los dos planos z, y se reporta el
    volumen perdido), se extraen las caras de borde tri6 salientes, las del
    techo (cargadas) y la capa superficial: elementos con una cara de borde
    que NO este sobre un plano del cubo.

    Devuelve (nodos, elems, vol, caras_techo, superficie, z0, z1, informe).
    Sirve igual para una malla que no salga de voxeles (la validacion con la
    esfera analitica usa esta misma funcion).
    """
    spacing = np.asarray(spacing, float)
    z0, z1 = -0.5 * spacing[2], (forma[2] - 0.5) * spacing[2]
    tol = 1e-6 * float(spacing[2])
    # Los nodos vienen de marching cubes en FLOAT32: a z ~ 5 mm el techo queda
    # a ~4e-7 mm de su plano (medido en el VOI proximal de H4 a 48^3), mas que
    # cualquier tolerancia razonable en unidades de h. Se devuelven a su plano
    # EXACTO, en float64, los nodos a menos de 1e-3 h de el. Una arista con
    # los dos extremos en el plano tiene alli tambien su nodo intermedio.
    nodos = np.array(nodos, dtype=np.float64, copy=True)
    for e in range(3):
        for v in (-0.5 * spacing[e], (forma[e] - 0.5) * spacing[e]):
            cerca = np.abs(nodos[:, e] - v) <= 1e-3 * spacing[e]
            nodos[cerca, e] = v
    vol = _volumen_tet(nodos, elems)
    if (vol <= 0).any():
        raise RuntimeError(f"{int((vol <= 0).sum())} tetraedros con "
                           "volumen no positivo.")
    V_total = float(vol.sum())
    port, n_comp = _componentes_portantes(nodos, elems, z0, z1, tol)
    elems, vol = elems[port], vol[port]
    nodos, elems = _compactar(nodos, elems)
    V = float(vol.sum())
    dueno, caras_b = caras_borde_tet10(elems)
    caras_b = _orientar_saliente(nodos, elems, dueno, caras_b)
    # Caras sobre los planos del cubo: todas sus esquinas en el plano.
    en_plano = np.zeros(caras_b.shape[0], bool)
    techo = np.zeros(caras_b.shape[0], bool)
    for e in range(3):
        a, b = -0.5 * spacing[e], (forma[e] - 0.5) * spacing[e]
        c = nodos[caras_b[:, :3], e]
        t = 1e-6 * float(spacing[e])
        en_a = np.all(np.abs(c - a) <= t, axis=1)
        en_b = np.all(np.abs(c - b) <= t, axis=1)
        en_plano |= en_a | en_b
        if e == 2:
            techo = en_b
    caras = caras_b[techo]
    sup = np.zeros(elems.shape[0], bool)
    sup[dueno[~en_plano]] = True
    V_caja = float(np.prod(np.asarray(forma) * spacing))
    inf = {"n_componentes": n_comp, "V_malla_total_mm3": V_total,
           "V_malla_mm3": V,
           "descartado_pct": 100.0 * (V_total - V) / V_total,
           "BVTV_malla": V / V_caja,
           "vol_elem_min_mm3": float(vol.min()),
           "vol_elem_mediana_mm3": float(np.median(vol))}
    return nodos, elems, vol, caras, sup, z0, z1, inf


def malla_de_tet10(nodos, elems, forma, spacing, informe=None):
    """Dict de malla (como el de `mallar`) a partir de un TET10 cualquiera."""
    spacing = np.asarray(spacing, float)
    nodos, elems, vol, caras, sup, z0, z1, inf = preparar_tet10(
        np.asarray(nodos, float), np.asarray(elems, np.int64), forma, spacing)
    inf = {**(informe or {}), **inf, "tipo": "tet10",
           "n_nodos": int(nodos.shape[0]), "n_elems": int(elems.shape[0]),
           "n_superficie": int(sup.sum()), "n_caras_techo": int(caras.shape[0])}
    return {"tipo": "tet10", "nodos": nodos, "elems": elems,
            "caras_techo": caras, "vol_elem": vol, "superficie": sup,
            "A_bruta": float(forma[0] * spacing[0] * forma[1] * spacing[1]),
            "H": float(z1 - z0), "z0": z0, "z1": z1, "forma": tuple(forma),
            "spacing": spacing, "eje": 2, "huella": huella(nodos, elems),
            "informe": inf}


def mallar(BW, spacing, tipo="hex8", eje=2, progreso=None, **opciones):
    """Mascara -> malla lista para FEBio, con el eje de carga en z.

    Siempre se filtra PRIMERO a lo portante (`resistencia._solo_portante`,
    el mismo filtro del ensayo de la app): un fragmento que no une base y
    techo solo aporta una submatriz casi singular. Para cargar en X o Y se
    PERMUTA la mascara como `ensayo_compresion_eje` y se malla el resultado:
    la carga siempre va en z.

    Devuelve un dict con
      tipo, nodos (N, 3) mm, elems (M, 8|10) base 0,
      caras_techo   quad4 / tri6 cargadas, normal +z,
      vol_elem      volumen de cada elemento (mm^3),
      superficie    elementos de la capa superficial (bool, M): los que
                    tienen una cara en el borde libre del hueso. Las caras
                    sobre los seis planos del cubo NO cuentan —la regla de
                    las isocaps, la misma de `resistencia.capa_superficie`—.
      A_bruta, H, z0, z1, forma, spacing, eje, huella (SHA-256), informe.

    Con TET10, `informe` declara la perdida de volumen frente a los voxeles
    del hueso portante y frente a la superficie cruda de marching cubes, y el
    volumen descartado al exigir que la MALLA —no la mascara— una base y
    techo: el suavizado puede cortar un puntal de un voxel.
    """
    from .resistencia import _PERM, _solo_portante, capa_superficie
    from .solido import malla_hex, malla_tet10

    if tipo not in MALLAS:
        raise ValueError(f"malla desconocida: {tipo!r}")
    BW = np.asarray(BW, dtype=bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    p = _PERM[int(eje)]
    BW = np.transpose(BW, p)
    spacing = spacing[list(p)]
    nx, ny, nz = BW.shape

    BW_res = _solo_portante(BW)
    if not BW_res[:, :, 0].any() or not BW_res[:, :, -1].any():
        raise ValueError("No hay hueso que una la base con el techo.")
    frac_port = float(BW_res.sum() / max(BW.sum(), 1))
    V_vox = float(BW_res.sum() * np.prod(spacing))
    A_bruta = float(nx * spacing[0] * ny * spacing[1])
    inf = {"tipo": tipo, "eje": int(eje), "forma": [nx, ny, nz],
           "spacing_mm": spacing.tolist(), "BVTV_voxel": float(BW.mean()),
           "frac_portante_voxel": frac_port, "V_voxel_portante_mm3": V_vox}

    if tipo == "hex8":
        nodos, elems, _ = malla_hex(BW_res, spacing)
        from .escribe import _caras_techo_hex
        caras = _caras_techo_hex(nodos, elems)
        vol = np.full(elems.shape[0], float(np.prod(spacing)))
        sup = capa_superficie(BW_res)[np.nonzero(BW_res)]
        z0, z1 = 0.0, float(nz * spacing[2])
        inf.update(V_malla_mm3=float(vol.sum()), vol_pct_voxel=100.0,
                   descartado_pct=0.0)
    else:
        o = opciones_malla(**opciones)
        objetivo = V_vox if o["corregir_volumen"] else None
        nodos, elems, _s, it = malla_tet10(
            BW_res, spacing, suavizado=o["suavizado"], banda=o["banda"],
            decimado=o["decimado"], progreso=progreso,
            maxvolume=maxvolume_de_tamano(o["tam_max_mm"]),
            minratio=o["minratio"], ejes_planos=(0, 1, 2),
            volumen_objetivo=objetivo)
        if not it["orden_ok"]:
            raise RuntimeError("Orden de nodos TET10 incorrecto: "
                               + it["orden"]["msg"])
        nodos, elems, vol, caras, sup, z0, z1, extra = preparar_tet10(
            nodos, elems, BW.shape, spacing)
        V = extra["V_malla_mm3"]
        inf.update(extra)
        inf.update(
            vol_pct_voxel=100.0 * V / V_vox,
            vol_pct_MC=100.0 * V / it["V_mc"] if it["V_mc"] > 0 else None,
            perdida_pct=100.0 - 100.0 * V / V_vox,
            opciones={k: v for k, v in o.items()},
            maxvolume_mm3=maxvolume_de_tamano(o["tam_max_mm"]),
            tetgen={k: it[k] for k in ("estanca", "bordes_abiertos",
                                       "tri_superficie", "orden",
                                       "correccion_volumen", "tiempo_s")})
        inf["reservas_volumen"] = abs(inf["perdida_pct"]) > PERDIDA_VOLUMEN_MAX_PCT
        inf["reservas_descarte"] = inf["descartado_pct"] > DESCARTE_MAX_PCT

    inf.update(n_nodos=int(nodos.shape[0]), n_elems=int(elems.shape[0]),
               n_gdl=int(3 * nodos.shape[0]),
               n_superficie=int(sup.sum()), n_caras_techo=int(caras.shape[0]))
    return {"tipo": tipo, "nodos": nodos, "elems": elems,
            "caras_techo": caras, "vol_elem": vol, "superficie": sup,
            "A_bruta": A_bruta, "H": float(z1 - z0), "z0": z0, "z1": z1,
            "forma": (nx, ny, nz), "spacing": spacing, "eje": int(eje),
            "huella": huella(nodos, elems), "informe": inf}


def tamano_previsto(BW, spacing, tipo, n=None):
    """Elementos, GDL y memoria de FEBio previstos, sin mallar.

    hex8: EXACTO (voxeles portantes y sus nodos a la resolucion `n`). tet10:
    `tiempos.n_tet10` a partir de los voxeles de superficie, +-4 % en hueso
    trabecular (ver `tiempos`). Cuesta menos de un segundo a 48^3.
    """
    from . import tiempos
    from .resistencia import _solo_portante, capa_superficie
    if n is None:
        n = N_HEX_DEF if tipo == "hex8" else N_TET_DEF
    B0, _sp = _remuestrear(BW, spacing, n)
    B = _solo_portante(B0)
    frac = float(B.sum() / max(B0.sum(), 1))
    if tipo == "hex8":
        i, j, k = np.nonzero(B)
        nn = np.array(B.shape) + 1
        esq = np.array([[a, b, c] for a in (0, 1) for b in (0, 1)
                        for c in (0, 1)])
        ids = np.unique(((i[:, None] + esq[:, 0]) + nn[0] * (
            (j[:, None] + esq[:, 1]) + nn[1] * (k[:, None] + esq[:, 2]))))
        gdl = 3 * int(ids.size)
        n_el = int(i.size)
    else:
        n_el = int(round(tiempos.n_tet10(capa_superficie(B).sum(), B.sum())))
        gdl = int(round(tiempos.gdl_tet10(n_el)))
    return {"tipo": tipo, "n": int(max(B.shape)), "n_elems": n_el,
            "gdl": gdl, "memoria_MB": tiempos.memoria_febio(tipo, gdl),
            "exacto": tipo == "hex8", "frac_portante": frac}


def memoria_equipo_MB():
    """RAM fisica total del equipo en MB (None si no se puede leer)."""
    try:
        import ctypes

        class MEMORYSTATUSEX(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong),
                        ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong),
                        ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong),
                        ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong),
                        ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("sullAvailExtendedVirtual", ctypes.c_ulonglong)]
        m = MEMORYSTATUSEX()
        m.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
        return m.ullTotalPhys / 2 ** 20
    except Exception:
        try:
            return (os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
                    / 2 ** 20)
        except (ValueError, OSError, AttributeError):
            return None


class ErrorMalla(RuntimeError):
    """El mallado TET10 fallo (tetgen abortado o proceso caido)."""


def _mallar_hijo(BW, spacing, tipo, eje, opciones):
    return mallar(BW, spacing, tipo, eje=eje, **opciones)


def mallar_aislado(BW, spacing, tipo="tet10", eje=2, **opciones):
    """`mallar` en un PROCESO HIJO: si tetgen se cae, no se lleva la app.

    tetgen es codigo C sin red: sobre el VOI proximal de H4 remuestreado a
    20^3 (puntales de ~1 voxel, que el suavizado deja casi degenerados) mato
    el interprete con una violacion de segmento. En la GUI eso cerraba la
    ventana entera, con horas de sesion. Aqui el hijo muere solo y se lanza
    `ErrorMalla`, que la etapa registra como fallo y la cadena sigue.
    El hijo cuesta unos segundos de arranque (importar numpy y pyvista);
    con hex8 no hay tetgen y se malla en el mismo proceso.
    """
    if tipo == "hex8":
        return mallar(BW, spacing, tipo, eje=eje, **opciones)
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    from concurrent.futures.process import BrokenProcessPool
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=1, mp_context=ctx) as ex:
        fut = ex.submit(_mallar_hijo, np.asarray(BW, bool),
                        np.asarray(spacing, float), tipo, int(eje),
                        dict(opciones))
        try:
            return fut.result()
        except BrokenProcessPool as e:
            raise ErrorMalla(
                "El mallado TET10 se cayo (tetgen abortado). Suele pasar "
                "cuando la resolucion deja puntales de uno o dos voxeles: "
                "subir la resolucion de la mascara o usar Ladrillos.") from e


# ---------------------------------------------------------------------------
# Protocolos y analisis
# ---------------------------------------------------------------------------

#: Una entrada por protocolo; agregar uno aqui lo agrega a los dialogos.
#: 'app' es el ensayo por omision de la app (`resistencia`). 'tapia2026'
#: reproduce los cuatro valores que definen el ensayo de Tapia, Gonzalez,
#: Vidal & Salinas, Biology 2026;15:722 —los mismos que `visor.PAPER_*`, lo
#: comprueba el bloque 27—; alli el mallado fue SOLID187 de 0.05 mm en
#: ANSYS, de ahi `tam_elem_mm`. Si se edita cualquiera de los valores
#: publicados, el registro deja de llamarse Tapia (`nombre_protocolo`).
#: 'homogeneizacion' da el tensor elastico: periodico con hex8 (el mismo
#: problema discreto que `elastic.homogeneizar`) y cotas KUBC/SUBC con tet10.
PROTOCOLOS_FEBIO = {
    "app": {"tipo": "compresion", "etiqueta": "Ensayo de la app",
            "E_s": 20e9, "nu": 0.30, "sigma_app": 1e6, "carga_N": None,
            "apoyo": "deslizante"},
    "tapia2026": {"tipo": "compresion",
                  "etiqueta": "Protocolo de Tapia et al. (2026)",
                  "E_s": 18e9, "nu": 0.30, "sigma_app": None,
                  "carga_N": 100.0, "apoyo": "empotrado",
                  "tam_elem_mm": 0.05},
    "homogeneizacion": {"tipo": "homogeneizacion",
                        "etiqueta": "Homogeneización (tensor elástico)",
                        "E_s": 20e9, "nu": 0.30, "amplitud": 1e-5,
                        "escala_vacio": 1e-6},
}

#: Claves que definen un protocolo publicado: si alguna cambia, el registro
#: se renombra (p. ej. «Tapia (modificado)»).
CLAVES_PROTOCOLO = ("E_s", "nu", "sigma_app", "carga_N", "apoyo", "amplitud",
                    "escala_vacio")

#: Analisis de un protocolo de compresion. 'nl_plato' y 'nl_pistoia'
#: necesitan el lineal (la deformacion equivalente y la carga de Pistoia
#: salen de el): si no se pidio, se corre igual y se declara.
ANALISIS = ("lineal", "nl_fuerza", "nl_plato", "nl_pistoia")

#: Cargas del analisis lineal, como fraccion de la del protocolo: dos cargas
#: pequenas y extrapolacion a carga nula, u = 2 u(s1) - u(s2) reescalado.
#: Son las de `comparativa_febio` (1 y 2 kPa frente a 1 MPa).
FRACCIONES_LINEAL = (1e-3, 2e-3)


def protocolo(clave, **cambios):
    """Copia de un protocolo con `cambios` aplicados y su nombre de registro."""
    base = PROTOCOLOS_FEBIO[clave]
    p = dict(base)
    p.update(cambios)
    p["clave"] = clave
    p["modificado"] = any(
        k in base and not _iguales(p.get(k), base.get(k))
        for k in CLAVES_PROTOCOLO)
    p["nombre"] = nombre_protocolo(clave, p["modificado"])
    return p


def _iguales(a, b):
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return float(a) == float(b)
    return a == b


def nombre_protocolo(clave, modificado=False):
    base = {"app": "App", "tapia2026": "Tapia",
            "homogeneizacion": "Homogeneización"}.get(clave, clave)
    return f"{base} (modificado)" if modificado else base


def _sigma_ref(prot, A_bruta):
    """Tension aparente del protocolo en Pa (A_bruta en mm^2).

    Con `carga_N` se deduce de la seccion BRUTA, como `ensayo_compresion`
    con `unidad='mm'`: sigma[Pa] = 1e6 F[N] / A[mm^2].
    """
    if prot.get("carga_N") is not None:
        return 1e6 * float(prot["carga_N"]) / float(A_bruta)
    return float(prot["sigma_app"])


def _uz_medio(malla, u):
    """Desplazamiento vertical medio del techo, ponderado por AREA.

    Es la integral de uz sobre las caras cargadas entre su area. Con quad4
    bilineales cada nodo de una cara pesa A/4; con tri6 las esquinas pesan
    CERO y cada nodo intermedio A/3 (integral de las funciones de forma
    cuadraticas). La media simple de los nodos del techo —la definicion de
    la app con voxeles— pesaria de mas las esquinas de las caras.
    """
    from .escribe import area_caras
    caras = malla["caras_techo"]
    A = area_caras(malla["nodos"], caras)
    uz = u[:, 2]
    if caras.shape[1] == 4:
        m = uz[caras].mean(axis=1)
    else:
        m = uz[caras[:, 3:]].mean(axis=1)
    return float((A * m).sum() / A.sum())


def _uz_medio_nodal(malla, u):
    """Media simple de uz en los nodos del techo: la definicion de la app."""
    z = malla["nodos"][:, 2]
    tolz = 1e-9 * max(float(z.max() - z.min()), 1.0)
    return float(u[z >= z.max() - tolz, 2].mean())


def _fuerza(malla, sigma):
    """F = -sum(sigma_zz V_e) / H (N), con sigma en Pa y mm.

    FEBio 4.5 escribe reacciones nulas en los `zero displacement`, asi que la
    fuerza se toma de la integral de volumen de la tension, exacta en el
    problema discreto lineal (desplazamiento virtual v = z e_z). En no lineal
    es tension de Cauchy sobre el volumen de referencia: error del orden de
    la deformacion del tejido (~1e-3).
    """
    return float(-(sigma[:, 2] * malla["vol_elem"]).sum() / malla["H"]) * 1e-6


def _escribir(malla, ruta, prot, sigma=None, eps_plato=None, material="svk",
              pasos=1, comentario="", adaptativo=False):
    from .escribe import escribir_febio_ensayo
    return escribir_febio_ensayo(
        malla["nodos"], malla["elems"], ruta, E_s=prot["E_s"],
        nu_s=prot["nu"], sigma_app=(sigma if sigma is not None else 0.0),
        A_bruta=malla["A_bruta"], apoyo=prot["apoyo"], eps_plato=eps_plato,
        caras_techo=malla["caras_techo"], material=material, pasos=pasos,
        comentario=comentario, adaptativo=adaptativo)


def estadisticos(malla, sigma, prot, sigma_ref, F_total_Pa_mm2):
    """p99 de superficie y Pistoia de un campo de FEBio, ponderados por volumen.

    Se usan las MISMAS funciones que la app (`criterio_pistoia`,
    `estadisticos_vm`), a las que se pasa `vol_solido`: con hexaedros de un
    voxel los pesos son iguales y el resultado es bit a bit el de siempre.
    Con TET10 la tension de cada elemento es la media de sus puntos de Gauss.
    """
    from .resistencia import criterio_pistoia
    res = {"ok": True, "sigma_app": float(sigma_ref),
           "F_total": float(F_total_Pa_mm2),
           "eps_eff_solido": eps_eff(sigma, prot["E_s"], prot["nu"]),
           "vm_solido": von_mises(sigma),
           "superficie_solido": np.asarray(malla["superficie"], bool),
           "vol_solido": np.asarray(malla["vol_elem"], float)}
    return criterio_pistoia(res)


def ensayo(malla, prot, carpeta, analisis=("lineal",), exe=None, hilos=None,
           cancelar=None, material="svk", pasos=1, progreso=None,
           conservar=True, prefijo="ensayo"):
    """Ensayo de compresion de un protocolo en FEBio sobre `malla`.

    Corre los `analisis` pedidos (ver `ANALISIS`) y devuelve un registro con
    procedencia. Cada corrida de FEBio que falla se registra y NO detiene las
    demas; solo `Cancelado` sale hacia arriba.

      lineal      FEBio a s1 y s2 (FRACCIONES_LINEAL de la carga) y
                  extrapolacion a carga nula: es comparable 1:1 con la app.
                  Da E_app, la fuerza (comprobacion de equilibrio), el p99 de
                  von Mises en la capa superficial y Pistoia.
      nl_fuerza   no lineal a la carga del protocolo, con fuerza impuesta.
      nl_plato    plato rigido sin friccion: lineal extrapolado del plato y no
                  lineal a la deformacion equivalente a la carga del
                  protocolo (la del ensayo lineal con fuerza), como
                  `comparativa_febio/plato.py`.
      nl_pistoia  no lineal a la carga de fallo de Pistoia del lineal.

    E_app se define con uz medio del techo PONDERADO POR AREA (`_uz_medio`).
    Con hex8 se da ademas `E_app_nodal`, la definicion de la app (media simple
    de nodos), que es la que se compara 1:1 con `ensayo_compresion`.
    """
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    analisis = [a for a in ANALISIS if a in set(analisis)]
    if not analisis:
        raise ValueError("No se pidio ningun analisis.")
    necesita_lineal = "lineal" in analisis or any(
        a in analisis for a in ("nl_plato", "nl_pistoia"))
    sref = _sigma_ref(prot, malla["A_bruta"])
    F_ref_N = sref * malla["A_bruta"] * 1e-6
    H = malla["H"]
    exe = Path(exe) if exe else localizar()
    reg = {"protocolo": prot.get("clave"), "nombre": prot.get("nombre"),
           "modificado": bool(prot.get("modificado")),
           "parametros": {k: prot.get(k) for k in CLAVES_PROTOCOLO
                          if k in prot},
           "malla": malla["tipo"], "huella_malla": malla["huella"],
           "informe_malla": malla["informe"], "eje": malla["eje"],
           "sigma_ref_Pa": sref, "F_ref_N": F_ref_N,
           "material_no_lineal": material, "pasos": int(pasos),
           "febio": {"version": None, "ruta": str(exe) if exe else None,
                     "origen": origen(exe)},
           "corridas": [], "fallos": [], "analisis_pedidos": list(analisis),
           "carpeta": str(carpeta)}
    n_tot = (2 if necesita_lineal else 0) + analisis.count("nl_fuerza") \
        + 3 * analisis.count("nl_plato") + analisis.count("nl_pistoia")
    hechas = [0]

    def run(nombre, **kw):
        """Escribe y corre un .feb; devuelve el campo leido o None."""
        if progreso:
            progreso(hechas[0], n_tot, nombre)
        feb = carpeta / f"{prefijo}_{nombre}.feb"
        info = _escribir(malla, feb, prot, **kw)
        try:
            r = correr(feb, exe=exe, hilos=hilos, cancelar=cancelar)
        except ErrorFEBio as e:
            reg["fallos"].append({"corrida": nombre, "msg": str(e),
                                  "log": e.log[-2000:]})
            hechas[0] += 1
            return None
        hechas[0] += 1
        reg["febio"]["version"] = r["version"]
        reg["corridas"].append({"corrida": nombre, "feb": feb.name,
                                "tiempo_s": r["tiempo_s"],
                                "memoria_MB": r["memoria_MB"],
                                "iteraciones": r["iteraciones"],
                                "F_N_escrita": info["F_N"],
                                "control": info["control"]})
        out = leer(feb)
        if not conservar:
            for suf in (".xplt",):
                f = feb.with_suffix(suf)
                if f.exists():
                    f.unlink()
        return out

    lin = None
    if necesita_lineal:
        campos = []
        for f in FRACCIONES_LINEAL:
            c = run(f"lineal_{f:g}", sigma=f * sref)
            campos.append(None if c is None else
                          {k: v / f for k, v in c.items()})
        if all(c is not None for c in campos):
            u = 2 * campos[0]["u"] - campos[1]["u"]
            s = 2 * campos[0]["sigma"] - campos[1]["sigma"]
            uz = _uz_medio(malla, u)
            eps_app = abs(uz) / H
            F = _fuerza(malla, s)
            lin = {"E_app": sref / eps_app, "eps_app": eps_app,
                   "F_N": F, "dF_rel": abs(F - F_ref_N) / F_ref_N,
                   "desvio_nl_s1": float(np.abs(campos[0]["u"] - u).max()
                                         / np.abs(u).max()),
                   "u_max_mm": float(np.abs(u).max())}
            if malla["tipo"] == "hex8":
                lin["E_app_nodal"] = sref / (abs(_uz_medio_nodal(malla, u)) / H)
            p = estadisticos(malla, s, prot, sref, sref * malla["A_bruta"])
            lin["pistoia"] = {k: v for k, v in p.items()
                              if isinstance(v, (int, float, bool, str))}
            reg["_campos_lineal"] = {"u": u, "sigma": s}
        reg["lineal"] = lin

    if "nl_fuerza" in analisis:
        c = run("nl_fuerza", sigma=sref, material=material, pasos=pasos,
                adaptativo=True)
        if c is not None:
            e = abs(_uz_medio(malla, c["u"])) / H
            reg["nl_fuerza"] = {"E_app": sref / e, "eps_app": e,
                                "F_N": _fuerza(malla, c["sigma"])}
            if lin:
                reg["nl_fuerza"]["dE_rel"] = reg["nl_fuerza"]["E_app"] \
                    / lin["E_app"] - 1.0
            reg["_campos_nl_fuerza"] = c

    if "nl_plato" in analisis and lin:
        eps_de = lambda sig: sig / lin["E_app"]           # noqa: E731
        Ep = []
        for f in FRACCIONES_LINEAL:
            c = run(f"plato_lineal_{f:g}", eps_plato=eps_de(f * sref))
            Ep.append(None if c is None else
                      _fuerza(malla, c["sigma"]) * 1e6
                      / malla["A_bruta"] / eps_de(f * sref))
        c = run("nl_plato", eps_plato=eps_de(sref), material=material,
                pasos=pasos, adaptativo=True)
        d = {"eps_plato": eps_de(sref)}
        if all(x is not None for x in Ep):
            d["E_app_lineal"] = 2 * Ep[0] - Ep[1]
            d["cociente_plato_fuerza_lineal"] = d["E_app_lineal"] / lin["E_app"]
        if c is not None:
            d["F_N"] = _fuerza(malla, c["sigma"])
            d["E_app"] = d["F_N"] * 1e6 / malla["A_bruta"] / d["eps_plato"]
            if "E_app_lineal" in d:
                d["dE_rel"] = d["E_app"] / d["E_app_lineal"] - 1.0
        reg["nl_plato"] = d

    if "nl_pistoia" in analisis and lin and lin["pistoia"].get("ok"):
        k = lin["pistoia"]["factor"]
        c = run("nl_pistoia", sigma=k * sref, material=material, pasos=pasos,
                adaptativo=True)
        d = {"factor": k, "sigma_Pa": k * sref}
        if c is not None:
            e = abs(_uz_medio(malla, c["u"])) / H
            d.update(E_app=k * sref / e, eps_app=e,
                     dE_rel=(k * sref / e) / lin["E_app"] - 1.0,
                     F_N=_fuerza(malla, c["sigma"]))
        reg["nl_pistoia"] = d

    if progreso:
        progreso(n_tot, n_tot, "listo")
    reg["ok"] = bool(reg["corridas"]) and not reg["fallos"]
    return reg


# ---------------------------------------------------------------------------
# Un analisis completo: el UNICO camino de calculo de las dos puertas
# ---------------------------------------------------------------------------

#: Resolucion por omision de la malla de ladrillos: la de los ensayos que se
#: reportan (Estudio_Convergencia: 40-48; a 32^3 la rigidez del VOI proximal pierde un
#: ~10 %).
N_HEX_DEF = 40
#: Resolucion por omision de la mascara de la que sale la malla suave. Algo
#: mas fina que la de ladrillos (el suavizado de Taubin se come puntales de
#: uno o dos voxeles, y a 40^3 Tb.Th/h ronda 2), pero no 64: el VOI proximal
#: de H4 a 64^3 son 350 000 TET10 y 1.9 millones de GDL (medido), que con
#: Pardiso no caben en un equipo de 16 GB. A 48^3: 179 000 y 1.0 millones.
N_TET_DEF = 48


def _remuestrear(BW, spacing, n):
    from .elastic import remuestrear_bw
    if n is None or int(n) == max(BW.shape):
        return np.asarray(BW, bool), np.asarray(spacing, float)
    return remuestrear_bw(BW, spacing, int(n))


def ensayo_app(BW, spacing, prot, eje=2):
    """El ensayo LINEAL de la app con los valores del protocolo.

    Es la columna «app (hex8, lineal)» de la tabla. Mismo `resistencia` que
    el panel de Analisis mecanico; nada propio.
    """
    from .resistencia import criterio_pistoia, ensayo_compresion_eje
    kw = {"E_s": prot["E_s"], "nu_s": prot["nu"], "apoyo": prot["apoyo"]}
    if prot.get("carga_N") is not None:
        kw["carga_N"] = prot["carga_N"]
    else:
        kw["sigma0"] = prot["sigma_app"]
    t0 = time.perf_counter()
    r = ensayo_compresion_eje(BW, spacing, eje=eje, **kw)
    out = {"ok": bool(r.get("ok")), "msg": r.get("msg", ""),
           "tiempo_s": time.perf_counter() - t0,
           "solver": r.get("solver"), "residuo_rel": r.get("residuo_rel")}
    if r.get("ok"):
        p = criterio_pistoia(r)
        out.update(E_app=r["E_app"], n_elem=r["n_elem"],
                   frac_portante=r["frac_portante"],
                   pistoia={k: v for k, v in p.items()
                            if isinstance(v, (int, float, bool, str))})
    return out


def analizar(BW, spacing, prot, malla="hex8", analisis=("lineal",), eje=2,
             carpeta=".", n=None, exe=None, hilos=None, cancelar=None,
             material="svk", pasos=1, progreso=None, conservar=True,
             comparar_app=True, opciones_malla=None, etiqueta="estructura",
             conv_malla=False):
    """Una estructura, un protocolo de compresion, una malla, un eje.

    Es lo que corren las dos puertas de la GUI («FEM automatico (FEBio)» y
    «Analizar con FEBio…») y la CLI: remuestrea la mascara a `n` (vecino mas
    proximo, como la app; None = la resolucion de siempre de cada malla),
    malla, corre FEBio (`ensayo`) y, si `comparar_app`, el ensayo lineal de
    la app sobre la MISMA mascara remuestreada a la resolucion de ladrillos.

    `conv_malla` (solo TET10): repite el lineal con la mitad del tamano de
    elemento y sin decimar, y guarda `convergencia_malla` (dE_rel, dp99_rel),
    que `informe.comprobar` califica.

    Devuelve el registro de `ensayo` con `app`, `estructura`, `n`,
    `tiempo_malla_s` y los objetos en memoria `_malla` y `_campos_*` (que
    `registro_json` quita antes de guardar).
    """
    if n is None:
        n = N_HEX_DEF if malla == "hex8" else N_TET_DEF
    t0 = time.perf_counter()
    B, sp = _remuestrear(BW, spacing, n)
    op = dict(opciones_malla or {})
    if malla == "tet10" and op.get("tam_max_mm") is None and \
            prot.get("tam_elem_mm"):
        op["tam_max_mm"] = prot["tam_elem_mm"]
    m = mallar_aislado(B, sp, malla, eje=eje,
                       **(op if malla == "tet10" else {}))
    t_malla = time.perf_counter() - t0
    nombre = f"{etiqueta}_{prot.get('clave')}_{malla}_{'XYZ'[eje]}"
    reg = ensayo(m, prot, Path(carpeta) / nombre, analisis=analisis,
                 exe=exe, hilos=hilos, cancelar=cancelar, material=material,
                 pasos=pasos, progreso=progreso, conservar=conservar,
                 prefijo=nombre)
    reg.update(estructura=etiqueta, n=int(max(B.shape)), eje_nombre="XYZ"[eje],
               tiempo_malla_s=t_malla, opciones_malla=op, _malla=m)
    if conv_malla and malla == "tet10" and reg.get("lineal"):
        # Segunda malla: la MITAD del tamano de elemento (de la mediana del
        # volumen de la primera si no se fijo) y sin decimar la superficie,
        # sobre la MISMA mascara. Solo el lineal: es lo que se compara.
        V_med = float(np.median(m["vol_elem"]))
        tam1 = op.get("tam_max_mm") or (6.0 * np.sqrt(2.0) * V_med) ** (1 / 3)
        op2 = dict(op, tam_max_mm=0.5 * float(tam1), decimado=0)
        try:
            m2 = mallar_aislado(B, sp, "tet10", eje=eje, **op2)
            r2 = ensayo(m2, prot, Path(carpeta) / (nombre + "_fina"),
                        analisis=("lineal",), exe=exe, hilos=hilos,
                        cancelar=cancelar, conservar=conservar,
                        prefijo=nombre + "_fina")
            l1, l2 = reg["lineal"], r2.get("lineal")
            d = {"tam_mm": [float(tam1), 0.5 * float(tam1)],
                 "n_elems": [int(m["elems"].shape[0]),
                             int(m2["elems"].shape[0])],
                 "huella_fina": m2["huella"], "fallos": r2["fallos"]}
            if l2:
                d["dE_rel"] = l2["E_app"] / l1["E_app"] - 1.0
                p1 = l1["pistoia"].get("vm_p99_superficie")
                p2 = l2["pistoia"].get("vm_p99_superficie")
                if p1 and p2:
                    d["dp99_rel"] = p2 / p1 - 1.0
            reg["convergencia_malla"] = d
        except ErrorMalla as e:
            reg["convergencia_malla"] = {"fallos": [{"corrida": "malla_fina",
                                                     "msg": str(e)}]}
    if comparar_app and "lineal" in reg["analisis_pedidos"]:
        Ba, spa = _remuestrear(BW, spacing, N_HEX_DEF if malla == "tet10"
                               else n)
        reg["app"] = ensayo_app(Ba, spa, prot, eje=eje)
        reg["app"]["n"] = int(max(Ba.shape))
    return reg


def registro_json(reg):
    """Copia serializable del registro: sin campos ni mallas en memoria."""
    def limpio(v):
        if isinstance(v, dict):
            return {k: limpio(x) for k, x in v.items()
                    if not str(k).startswith("_")}
        if isinstance(v, (list, tuple)):
            return [limpio(x) for x in v]
        if isinstance(v, np.generic):
            return v.item()
        if isinstance(v, np.ndarray):
            return v.tolist()
        return v
    return limpio(reg)


def fila_tabla(reg):
    """Numeros de la tabla comparativa de un registro (MPa y adimensionales).

    `dE_app_implementacion` (solo hex8): FEBio frente a la app con la
    DEFINICION de la app (media simple de nodos del techo); debe ser ~0.
    Las demas diferencias frente a la app mezclan malla y definicion y se
    rotulan asi en la GUI.
    """
    f = {"estructura": reg.get("estructura"), "protocolo": reg.get("nombre"),
         "malla": reg.get("malla"), "eje": reg.get("eje_nombre"),
         "n": reg.get("n"), "ok": reg.get("ok"),
         "fallos": len(reg.get("fallos", []))}
    lin = reg.get("lineal") or {}
    app = reg.get("app") or {}
    inf = reg.get("informe_malla") or {}
    f["BVTV_malla"] = inf.get("BVTV_malla", inf.get("BVTV_voxel"))
    f["perdida_volumen_pct"] = inf.get("perdida_pct")
    f["descartado_pct"] = inf.get("descartado_pct")
    if lin:
        p = lin.get("pistoia", {})
        f.update(E_app_MPa=lin["E_app"] / 1e6,
                 vm_p99_sup_MPa=(p.get("vm_p99_superficie", np.nan) / 1e6),
                 sigma_fallo_MPa=(p.get("sigma_fallo", np.nan) / 1e6),
                 dF_rel=lin.get("dF_rel"))
        if "E_app_nodal" in lin and app.get("ok"):
            f["dE_app_implementacion"] = lin["E_app_nodal"] / app["E_app"] - 1
    if app.get("ok"):
        pa = app.get("pistoia", {})
        f.update(E_app_app_MPa=app["E_app"] / 1e6,
                 vm_p99_sup_app_MPa=pa.get("vm_p99_superficie", np.nan) / 1e6,
                 sigma_fallo_app_MPa=pa.get("sigma_fallo", np.nan) / 1e6)
    for k in ("nl_fuerza", "nl_plato", "nl_pistoia"):
        d = reg.get(k) or {}
        if "dE_rel" in d:
            f[f"dE_{k}"] = d["dE_rel"]
    if (reg.get("nl_plato") or {}).get("cociente_plato_fuerza_lineal"):
        f["plato_fuerza_lineal"] = reg["nl_plato"]["cociente_plato_fuerza_lineal"]
    return f


# ---------------------------------------------------------------------------
# Homogeneizacion
# ---------------------------------------------------------------------------
#
# hex8   PERIODICA: el MISMO problema discreto que `elastic.homogeneizar`
#        (celda completa, vacio a 1e-6 E_s, nodo imagen = nodo + salto E.dx)
#        escrito con restricciones lineales de FEBio. Validado a 1e-6 del
#        maximo de C en el bloque 25 y en `comparativa_febio/porcino/`.
# tet10  COTAS. Una malla suave no tiene nodos emparejados en caras opuestas,
#        y el VOI real no es periodico, asi que no hay tensor periodico que
#        calcular. Se dan dos condiciones de borde uniformes sobre el hueso que
#        toca las caras del cubo:
#          KUBC  desplazamiento u = E.x impuesto; C[:, j] = <sigma> sobre la
#                celda (los poros no tienen tension). Cota SUPERIOR.
#          SUBC  traccion uniforme en cada cara, con la fuerza de la seccion
#                BRUTA (Sigma.n A_bruta) repartida sobre el hueso de esa cara
#                —misma convencion que el ensayo de la app—; la deformacion
#                media sale de los desplazamientos medios de caras opuestas.
#                Cota INFERIOR aproximada: con poros abiertos en las caras
#                <eps> no puede integrarse sobre todo el borde, y esta es la
#                variante practica en micro-EF de hueso (Pahr y Zysset 2008).
#        La tabla dice siempre «cotas KUBC/SUBC», nunca «tensor periodico».

VOIGT = ("xx", "yy", "zz", "yz", "xz", "xy")


def tensor_de_voigt(e):
    """Voigt [xx yy zz yz xz xy] con distorsiones de INGENIERIA -> tensor."""
    return np.array([[e[0], e[5] / 2, e[4] / 2],
                     [e[5] / 2, e[1], e[3] / 2],
                     [e[4] / 2, e[3] / 2, e[2]]], float)


def _cabecera(w, comentario, n_materiales, rtol=1e-12):
    w('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
    w(f"<!-- {comentario} -->\n")
    w('<febio_spec version="4.0">\n')
    w('\t<Module type="solid">\n\t\t<units>mm-N-s</units>\n\t</Module>\n')
    w("\t<Control>\n\t\t<analysis>STATIC</analysis>\n")
    w("\t\t<time_steps>1</time_steps>\n\t\t<step_size>1</step_size>\n")
    w("\t\t<solver>\n\t\t\t<max_refs>50</max_refs>\n")
    w('\t\t\t<qn_method type="BFGS">\n\t\t\t\t<max_ups>0</max_ups>\n'
      "\t\t\t</qn_method>\n\t\t\t<dtol>1e-9</dtol>\n")
    w(f"\t\t\t<etol>1e-12</etol>\n\t\t\t<rtol>{rtol:g}</rtol>\n")
    w("\t\t\t<lstol>0.9</lstol>\n")
    w('\t\t\t<linear_solver type="pardiso"/>\n\t\t</solver>\n')
    w("\t</Control>\n")


def _pie(w, stem, nodos=True):
    w('\t<LoadData>\n\t\t<load_controller id="1" name="rampa" '
      'type="loadcurve">\n')
    w("\t\t\t<interpolate>LINEAR</interpolate>\n\t\t\t<points>\n")
    w("\t\t\t\t<point>0,0</point>\n\t\t\t\t<point>1,1</point>\n")
    w("\t\t\t</points>\n\t\t</load_controller>\n\t</LoadData>\n")
    # El plotfile NO es opcional: FEBio 4.5.0 termina con una violacion de
    # acceso (0xC0000005) si Output solo tiene logfile.
    w('\t<Output>\n\t\t<plotfile type="febio">\n'
      '\t\t\t<var type="displacement"/>\n\t\t</plotfile>\n\t\t<logfile>\n')
    w(f'\t\t\t<element_data data="sx;sy;sz;syz;sxz;sxy" delim=" " '
      f'file="{stem}_s.txt"/>\n')
    if nodos:
        w(f'\t\t\t<node_data data="ux;uy;uz" delim=" " file="{stem}_u.txt"/>\n')
    w("\t\t</logfile>\n\t</Output>\n</febio_spec>\n")


def _nodos_elems(w, nodos, grupos, tipo):
    w('\t<Mesh>\n\t\t<Nodes name="todos">\n')
    w("".join(f'\t\t\t<node id="{i}">{p[0]:.12g},{p[1]:.12g},{p[2]:.12g}'
              "</node>\n" for i, p in enumerate(nodos, start=1)))
    w("\t\t</Nodes>\n")
    for nom, ids, conn in grupos:
        w(f'\t\t<Elements type="{tipo}" name="{nom}">\n')
        w("".join(f'\t\t\t<elem id="{i}">' + ",".join(map(str, c))
                  + "</elem>\n"
                  for i, c in zip(ids.tolist(), (conn + 1).tolist())))
        w("\t\t</Elements>\n")


def escribir_periodico(BW, spacing, ruta, E_voigt, E_s=20e9, nu_s=0.30,
                       escala_vacio=1e-6):
    """Celda periodica hex8 bajo la deformacion macroscopica E_voigt.

    Copia de `comparativa_febio/porcino/homog_febio.escribir_periodico`
    (validada a 1e-6 contra `elastic.homogeneizar`), que sigue alli igual.
    Cada nodo con algun indice igual a n es DEPENDIENTE de su imagen modulo n
    con el salto E.(x' - x); la traslacion se quita fijando el nodo (0,0,0).
    """
    ruta = Path(ruta)
    BW = np.asarray(BW, bool)
    nx, ny, nz = BW.shape
    sp = np.asarray(spacing, float)
    I, J, K = np.meshgrid(np.arange(nx + 1), np.arange(ny + 1),
                          np.arange(nz + 1), indexing="ij")
    I, J, K = (a.ravel(order="F") for a in (I, J, K))
    nodos = np.stack([I * sp[0], J * sp[1], K * sp[2]], axis=1).astype(float)
    ijk = np.stack([I, J, K], axis=1)

    def nid(a, b, c):
        return a + (nx + 1) * b + (nx + 1) * (ny + 1) * c

    ex, ey, ez = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz),
                             indexing="ij")
    ex, ey, ez = (a.ravel(order="F") for a in (ex, ey, ez))
    elems = np.stack([nid(ex, ey, ez), nid(ex + 1, ey, ez),
                      nid(ex + 1, ey + 1, ez), nid(ex, ey + 1, ez),
                      nid(ex, ey, ez + 1), nid(ex + 1, ey, ez + 1),
                      nid(ex + 1, ey + 1, ez + 1), nid(ex, ey + 1, ez + 1)],
                     axis=1)
    solido = BW.ravel(order="F")
    Emac = tensor_de_voigt(np.asarray(E_voigt, float))
    E_MPa = float(E_s) / 1e6
    n_lim = np.array([nx, ny, nz])
    dep = np.where((ijk == n_lim).any(axis=1))[0]
    img = ijk[dep] % n_lim
    ind = img[:, 0] + (nx + 1) * img[:, 1] + (nx + 1) * (ny + 1) * img[:, 2]
    offset = (nodos[dep] - nodos[ind]) @ Emac.T

    with open(ruta, "w", encoding="utf-8") as fh:
        w = fh.write
        _cabecera(w, "Celda periodica hex8 (spinpy/febio.py). E_macro (Voigt, "
                  "ingenieria) = " + ", ".join(f"{v:.6g}" for v in E_voigt), 2)
        w("\t<Material>\n")
        for i, (nom, E) in enumerate((("hueso", E_MPa),
                                      ("vacio", E_MPa * escala_vacio)), 1):
            w(f'\t\t<material id="{i}" name="{nom}" type="isotropic elastic">'
              f"\n\t\t\t<E>{E:.17g}</E>\n\t\t\t<v>{nu_s:.9g}</v>\n"
              "\t\t</material>\n")
        w("\t</Material>\n")
        ids = np.arange(1, elems.shape[0] + 1)
        grupos = [(nom, ids[m], elems[m]) for nom, m in
                  (("hueso", solido), ("vacio", ~solido)) if m.any()]
        _nodos_elems(w, nodos, grupos, "hex8")
        w('\t\t<NodeSet name="origen">1</NodeSet>\n\t</Mesh>\n')
        w("\t<MeshDomains>\n")
        for nom, _, _ in grupos:
            w(f'\t\t<SolidDomain name="{nom}" mat="{nom}"/>\n')
        w("\t</MeshDomains>\n\t<Boundary>\n")
        w('\t\t<bc name="origen" type="zero displacement" node_set="origen">'
          "\n\t\t\t<x_dof>1</x_dof>\n\t\t\t<y_dof>1</y_dof>\n"
          "\t\t\t<z_dof>1</z_dof>\n\t\t</bc>\n")
        partes = []
        for d_, i_, off in zip((dep + 1).tolist(), (ind + 1).tolist(),
                               offset.tolist()):
            for c, dn in enumerate("xyz"):
                partes.append(
                    f'\t\t<bc type="linear constraint"><node>{d_}</node>'
                    f"<dof>{dn}</dof><offset>{off[c]:.17g}</offset>"
                    f"<child_dof><node>{i_}</node><dof>{dn}</dof>"
                    "<value>1</value></child_dof></bc>\n")
        w("".join(partes))
        w("\t</Boundary>\n")
        _pie(w, ruta.stem, nodos=False)
    return {"ruta": str(ruta), "n_elems": int(elems.shape[0]),
            "n_restricciones": int(3 * dep.size)}


def _caras_planos(malla):
    """Caras tri6 de borde sobre cada uno de los seis planos del cubo."""
    nodos, elems = malla["nodos"], malla["elems"]
    dueno, cb = caras_borde_tet10(elems)
    cb = _orientar_saliente(nodos, elems, dueno, cb)
    sp, forma = malla["spacing"], malla["forma"]
    out = {}
    for e in range(3):
        a, b = -0.5 * sp[e], (forma[e] - 0.5) * sp[e]
        c = nodos[cb[:, :3], e]
        t = 1e-6 * float(sp[e])
        out[(e, 0)] = (cb[np.all(np.abs(c - a) <= t, axis=1)], a)
        out[(e, 1)] = (cb[np.all(np.abs(c - b) <= t, axis=1)], b)
    return out


def _mayor_componente(malla):
    """La malla TET10 reducida a su mayor componente conexo."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    nodos, elems = malla["nodos"], malla["elems"]
    e4 = elems[:, :4]
    n = nodos.shape[0]
    G = coo_matrix((np.ones(3 * e4.shape[0]),
                    (np.repeat(e4[:, 0], 3), e4[:, 1:].ravel())), shape=(n, n))
    _, lab = connected_components(G, directed=False)
    le = lab[e4[:, 0]]
    vol = malla["vol_elem"]
    tot = np.bincount(le, weights=vol)
    k = int(np.argmax(tot))
    keep = le == k
    m = malla_de_tet10(nodos, elems[keep], malla["forma"], malla["spacing"],
                       informe=malla["informe"])
    m["informe"]["descartado_homog_pct"] = float(
        100.0 * (vol.sum() - vol[keep].sum()) / vol.sum())
    return m


def _escribir_tet_uniforme(malla, ruta, E_s, nu_s, E_voigt=None, S_voigt=None):
    """KUBC (E_voigt, deformacion) o SUBC (S_voigt, tension Pa) sobre TET10."""
    ruta = Path(ruta)
    nodos, elems = malla["nodos"], malla["elems"]
    planos = _caras_planos(malla)
    E_MPa = float(E_s) / 1e6
    with open(ruta, "w", encoding="utf-8") as fh:
        w = fh.write
        tipo = "KUBC" if E_voigt is not None else "SUBC"
        from .escribe import RTOL
        _cabecera(w, f"Homogeneizacion {tipo} TET10 (spinpy/febio.py)", 1,
                  rtol=RTOL["tet10"])
        w('\t<Material>\n\t\t<material id="1" name="hueso" '
          'type="isotropic elastic">\n')
        w(f"\t\t\t<E>{E_MPa:.17g}</E>\n\t\t\t<v>{nu_s:.9g}</v>\n"
          "\t\t</material>\n\t</Material>\n")
        _nodos_elems(w, nodos, [("hueso", np.arange(1, elems.shape[0] + 1),
                                 elems)], "tet10")
        if E_voigt is not None:
            borde = np.unique(np.concatenate(
                [c.ravel() for c, _ in planos.values() if c.size]))
            w('\t\t<NodeSet name="caras">\n\t\t\t'
              + ",".join(map(str, (borde + 1).tolist())) + "\n\t\t</NodeSet>\n")
        else:
            sp, forma = malla["spacing"], malla["forma"]
            lo = -0.5 * sp
            hi = (np.asarray(forma) - 0.5) * sp
            anclas = []
            for objetivo in ([lo[0], lo[1], lo[2]], [hi[0], lo[1], lo[2]],
                             [lo[0], hi[1], lo[2]]):
                anclas.append(int(np.argmin(np.linalg.norm(
                    nodos - np.asarray(objetivo), axis=1))))
            for k, a in enumerate(anclas):
                w(f'\t\t<NodeSet name="ancla{k}">{a + 1}</NodeSet>\n')
            for (e, lado), (c, _) in planos.items():
                if c.size:
                    w(f'\t\t<Surface name="cara{e}{lado}">\n')
                    w("".join(f'\t\t\t<tri6 id="{i}">'
                              + ",".join(map(str, f)) + "</tri6>\n"
                              for i, f in enumerate((c + 1).tolist(), 1)))
                    w("\t\t</Surface>\n")
        w("\t</Mesh>\n")
        w('\t<MeshDomains>\n\t\t<SolidDomain name="hueso" mat="hueso" '
          f'elem_type="{_regla_tet10()}"/>\n\t</MeshDomains>\n')
        w("\t<Boundary>\n")
        if E_voigt is not None:
            F = np.eye(3) + tensor_de_voigt(np.asarray(E_voigt, float))
            w('\t\t<bc name="kubc" type="prescribed deformation" '
              'node_set="caras">\n\t\t\t<scale lc="1">1</scale>\n'
              "\t\t\t<F>" + ",".join(f"{v:.17g}" for v in F.ravel())
              + "</F>\n\t\t</bc>\n")
        else:
            for k, dofs in enumerate(((1, 1, 1), (0, 1, 1), (0, 0, 1))):
                w(f'\t\t<bc name="ancla{k}" type="zero displacement" '
                  f'node_set="ancla{k}">\n'
                  f"\t\t\t<x_dof>{dofs[0]}</x_dof>\n\t\t\t<y_dof>{dofs[1]}"
                  f"</y_dof>\n\t\t\t<z_dof>{dofs[2]}</z_dof>\n\t\t</bc>\n")
        w("\t</Boundary>\n")
        if S_voigt is not None:
            from .escribe import area_caras
            Sg = tensor_de_voigt(np.asarray(S_voigt, float))
            # tensor_de_voigt divide las componentes 3..5 entre 2 (ingenieria);
            # para TENSIONES no hay factor: se deshace.
            for i, j in ((1, 2), (0, 2), (0, 1)):
                Sg[i, j] *= 2.0
                Sg[j, i] *= 2.0
            sp, forma = malla["spacing"], malla["forma"]
            L = np.asarray(forma) * sp
            w("\t<Loads>\n")
            for (e, lado), (c, _) in planos.items():
                if not c.size:
                    continue
                nvec = np.zeros(3)
                nvec[e] = 1.0 if lado else -1.0
                A_bruta = float(np.prod(np.delete(L, e)))
                A_osea = float(area_caras(nodos, c).sum())
                t = Sg @ nvec * (A_bruta / A_osea) / 1e6       # MPa
                w(f'\t\t<surface_load name="t{e}{lado}" type="traction" '
                  f'surface="cara{e}{lado}">\n\t\t\t<scale lc="1">1</scale>\n'
                  "\t\t\t<traction>" + ",".join(f"{v:.17g}" for v in t)
                  + "</traction>\n\t\t</surface_load>\n")
            w("\t</Loads>\n")
        _pie(w, ruta.stem, nodos=True)
    return planos


def _regla_tet10():
    from .escribe import REGLA_TET10
    return REGLA_TET10


def homogeneizar(BW, spacing, prot, malla="hex8", carpeta=".", n=None,
                 exe=None, hilos=None, cancelar=None, progreso=None,
                 comparar_app=True, opciones_malla=None, etiqueta="estructura",
                 conservar=True):
    """Tensor elastico en FEBio: periodico (hex8) o cotas KUBC/SUBC (tet10).

    Cada columna se obtiene a dos amplitudes y se extrapola a amplitud nula
    (FEBio es no lineal geometricamente). Devuelve un registro con C en Pa
    (Voigt, ingenieria) —`C_periodico`, o `C_KUBC` y `C_SUBC`— y, para hex8
    con `comparar_app`, `C_app` de `elastic.homogeneizar` sobre la misma
    mascara (el mismo problema discreto: su diferencia mide la
    implementacion).
    """
    from .escribe import area_caras
    if n is None:
        n = N_HEX_DEF if malla == "hex8" else N_TET_DEF
    B, sp = _remuestrear(BW, spacing, n)
    carpeta = Path(carpeta) / f"{etiqueta}_homog_{malla}"
    carpeta.mkdir(parents=True, exist_ok=True)
    exe = Path(exe) if exe else localizar()
    E_s, nu = prot["E_s"], prot["nu"]
    amp = float(prot.get("amplitud", 1e-5))
    reg = {"protocolo": prot.get("clave"), "nombre": prot.get("nombre"),
           "modificado": bool(prot.get("modificado")), "tipo": "homogeneizacion",
           "malla": malla, "estructura": etiqueta, "n": int(max(B.shape)),
           "parametros": {k: prot.get(k) for k in CLAVES_PROTOCOLO if k in prot},
           "febio": {"version": None, "ruta": str(exe) if exe else None,
                     "origen": origen(exe)},
           "corridas": [], "fallos": [], "carpeta": str(carpeta)}
    n_tot = 12 if malla == "hex8" else 24
    hechas = [0]

    def run(feb):
        if progreso:
            progreso(hechas[0], n_tot, feb.stem)
        try:
            r = correr(feb, exe=exe, hilos=hilos, cancelar=cancelar)
        except ErrorFEBio as e:
            reg["fallos"].append({"corrida": feb.stem, "msg": str(e),
                                  "log": e.log[-2000:]})
            return None
        finally:
            hechas[0] += 1
        reg["febio"]["version"] = r["version"]
        reg["corridas"].append({"corrida": feb.stem, "tiempo_s": r["tiempo_s"],
                                "memoria_MB": r["memoria_MB"]})
        return r

    if malla == "hex8":
        V = float(np.prod(np.asarray(B.shape) * sp))
        C = np.full((6, 6), np.nan)
        for j in range(6):
            cols = []
            for a in (amp, 2 * amp):
                e = np.zeros(6)
                e[j] = a
                feb = carpeta / f"per_{VOIGT[j]}_{a:g}.feb"
                escribir_periodico(B, sp, feb, e, E_s=E_s, nu_s=nu,
                                   escala_vacio=prot.get("escala_vacio", 1e-6))
                if run(feb) is None:
                    break
                s = ultimo_registro(feb.with_name(feb.stem + "_s.txt"), 6)
                cols.append(s.mean(axis=0) * 1e6 / a)
            if len(cols) == 2:
                C[:, j] = 2 * cols[0] - cols[1]
        reg["C_periodico"] = C
        reg["condicion"] = "periodica"
        if comparar_app:
            from .elastic import homogeneizar as homog_app
            t0 = time.perf_counter()
            Ca, info = homog_app(B, E_s, nu, vox_size=sp)
            reg["C_app"] = np.asarray(Ca)
            reg["app"] = {"ok": bool(info.get("ok")),
                          "tiempo_s": time.perf_counter() - t0,
                          "residuo_rel": info.get("residuo_rel")}
            if np.isfinite(C).all():
                reg["dC_app_rel"] = float(np.abs(C - Ca).max()
                                          / np.abs(Ca).max())
    else:
        op = dict(opciones_malla or {})
        m = _mayor_componente(mallar_aislado(B, sp, "tet10", **op))
        V = float(np.prod(np.asarray(m["forma"]) * m["spacing"]))
        vol = m["vol_elem"]
        CK = np.full((6, 6), np.nan)
        SS = np.full((6, 6), np.nan)
        for j in range(6):
            colK, colS = [], []
            for a in (amp, 2 * amp):
                e = np.zeros(6)
                e[j] = a
                feb = carpeta / f"kubc_{VOIGT[j]}_{a:g}.feb"
                _escribir_tet_uniforme(m, feb, E_s, nu, E_voigt=e)
                if run(feb) is not None:
                    s = leer(feb)["sigma"]
                    colK.append((s * vol[:, None]).sum(axis=0) / V / a)
                sig = np.zeros(6)
                sig[j] = a * E_s
                feb = carpeta / f"subc_{VOIGT[j]}_{a:g}.feb"
                planos = _escribir_tet_uniforme(m, feb, E_s, nu, S_voigt=sig)
                if run(feb) is not None:
                    u = leer(feb)["u"]
                    G = np.zeros((3, 3))
                    L = np.asarray(m["forma"]) * m["spacing"]
                    for k in range(3):
                        med = []
                        for lado in (0, 1):
                            c, _ = planos[(k, lado)]
                            Af = area_caras(m["nodos"], c)
                            um = u[c[:, 3:]].mean(axis=1)
                            med.append((Af[:, None] * um).sum(0) / Af.sum())
                        G[:, k] = (med[1] - med[0]) / L[k]
                    eps = np.array([G[0, 0], G[1, 1], G[2, 2], G[1, 2] + G[2, 1],
                                    G[0, 2] + G[2, 0], G[0, 1] + G[1, 0]])
                    colS.append(eps / (a * E_s))
            if len(colK) == 2:
                CK[:, j] = 2 * colK[0] - colK[1]
            if len(colS) == 2:
                SS[:, j] = 2 * colS[0] - colS[1]
        reg["C_KUBC"] = CK
        reg["S_SUBC"] = SS
        reg["C_SUBC"] = (np.linalg.inv(0.5 * (SS + SS.T))
                         if np.isfinite(SS).all() else SS)
        reg["condicion"] = "cotas KUBC/SUBC"
        reg["informe_malla"] = m["informe"]
        reg["huella_malla"] = m["huella"]
        reg["_malla"] = m
    reg["V_celda_mm3"] = V
    reg["ok"] = not reg["fallos"]
    return reg


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    """python -m spinpy.febio VOI.vtk [...] --protocolo app tapia2026 ..."""
    import argparse
    import json

    from .io import leer_voi
    ap = argparse.ArgumentParser(
        prog="python -m spinpy.febio",
        description="Ensayos de spinpy resueltos en FEBio (hex8 y/o TET10).")
    ap.add_argument("vois", nargs="+", help="VOI(s): .vtk, .mat, carpeta o .tif")
    ap.add_argument("--protocolo", nargs="+", default=["app"],
                    choices=[k for k, v in PROTOCOLOS_FEBIO.items()
                             if v["tipo"] == "compresion"])
    ap.add_argument("--malla", nargs="+", default=["hex8", "tet10"],
                    choices=list(MALLAS))
    ap.add_argument("--analisis", nargs="+", default=["lineal"],
                    choices=list(ANALISIS))
    ap.add_argument("--ejes", default="Z", help="Z, X, Y o XYZ")
    ap.add_argument("--n-hex", type=int, default=N_HEX_DEF)
    ap.add_argument("--n-tet", type=int, default=N_TET_DEF)
    ap.add_argument("--tam", type=float, default=None,
                    help="tamano maximo de tetraedro (mm)")
    ap.add_argument("--sin-correccion", action="store_true",
                    help="no corregir la perdida de volumen del suavizado")
    ap.add_argument("--conv-malla", action="store_true",
                    help="TET10: repetir el lineal con la mitad del tamano")
    ap.add_argument("--material", default="svk", choices=["svk",
                                                          "neohookeano"])
    ap.add_argument("--pasos", type=int, default=1)
    ap.add_argument("--hilos", type=int, default=None)
    ap.add_argument("--febio", default=None, help="ruta a febio4.exe")
    ap.add_argument("--salida", default="./febio_resultados")
    a = ap.parse_args(argv)

    exe = localizar(a.febio)
    if exe is None:
        raise SystemExit("No se encontro FEBio (febio4.exe).")
    salida = Path(a.salida)
    salida.mkdir(parents=True, exist_ok=True)
    print(f"FEBio {version(exe)} ({origen(exe)}): {exe}")
    op = {"tam_max_mm": a.tam, "corregir_volumen": not a.sin_correccion}
    filas = []
    for ruta in a.vois:
        BW, spc = leer_voi(ruta)
        et = Path(ruta).stem
        for pk in a.protocolo:
            prot = protocolo(pk)
            for mt in a.malla:
                for e in a.ejes.upper():
                    eje = "XYZ".index(e)
                    print(f"== {et} · {prot['nombre']} · {mt} · {e}", flush=True)
                    reg = analizar(
                        BW, spc, prot, malla=mt, analisis=a.analisis, eje=eje,
                        carpeta=salida, n=a.n_hex if mt == "hex8" else a.n_tet,
                        exe=exe, hilos=a.hilos, material=a.material,
                        pasos=a.pasos, opciones_malla=op, etiqueta=et,
                        conv_malla=a.conv_malla)
                    f = fila_tabla(reg)
                    filas.append(f)
                    print("   " + ", ".join(
                        f"{k}={v:.5g}" if isinstance(v, float) else f"{k}={v}"
                        for k, v in f.items()), flush=True)
                    with open(salida / "resultados_fem.jsonl", "a",
                              encoding="utf-8") as fh:
                        fh.write(json.dumps(registro_json(reg),
                                            ensure_ascii=False) + "\n")
    return filas


if __name__ == "__main__":
    main()
