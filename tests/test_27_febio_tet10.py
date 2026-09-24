"""
Bloque 27 — FEBio dentro de spinpy: malla suave TET10 y estadisticos por volumen.

REFERENCIA
  `PLAN_Analizar_con_FEBio.md` §9 y `spinpy/febio.py`. El camino hex8 ya esta
  validado contra la app en el bloque 25 y en `comparativa_febio/`.

LO QUE SE VERIFICA
  (1) Orden de nodos TET10 de FEBio (con FEBio): un solo elemento con un campo
      de desplazamientos CUADRATICO impuesto en sus diez nodos. En el orden
      C3D10 —el de `solido.malla_tet10`— la tension media del elemento es la
      exacta en el centroide; con los nodos intermedios permutados FEBio
      aborta o da otra tension. Un campo lineal no lo detectaria: un elemento
      isoparametrico reproduce campos lineales con cualquier geometria.
  (2) Bloque macizo (con FEBio): E_app = E_s con hex8 y con TET10, y la
      fuerza de la integral de volumen de sigma_zz igual a la aplicada.
  (3) Estadisticos ponderados (sin FEBio): con pesos iguales, `criterio_
      pistoia` y `estadisticos_vm` dan BIT A BIT lo mismo con y sin
      `vol_solido`; la formula ponderada con pesos iguales coincide con numpy
      (hazen y tipo 7) a redondeo; y un caso con pesos distintos calculado a
      mano.
  (4) Escritor (sin FEBio): `escribir_febio_ensayo` con hex8, un paso y St.
      Venant-Kirchhoff escribe el MISMO problema que `escribir_febio` (todas
      las lineas salvo el comentario de cabecera), asi que hereda su
      validacion; con TET10 las caras tri6 del techo miran a +z, presion por
      area osea = sigma_app * A_bruta, y la regla va como ATRIBUTO.
  (5) Lectura sin FEBio: `leer` toma el ULTIMO bloque de un logfile de
      ejemplo, ordenado por identificador, y pasa MPa a Pa.
  (6) Malla (sin FEBio): una cavidad cerrada NO se rellena (antes PyMeshFix
      borraba su pared y tetgen la mallaba maciza); la correccion de volumen
      lleva la malla al volumen de voxeles; la capa superficial excluye las
      caras del cubo (un bloque macizo no tiene capa).
  (7) El protocolo de Tapia de `febio.PROTOCOLOS_FEBIO` es el de
      `visor.PAPER_*` (leido del fuente, sin importar Qt), y editar un valor
      lo renombra «Tapia (modificado)».
  (8) Homogeneizacion (con FEBio): en un bloque macizo las cotas KUBC y SUBC
      con TET10 coinciden con la matriz D isotropa.
  (9) Calificacion (sin FEBio): `informe.comprobar` gradua los registros de
      FEBio con los motivos declarados (fallo, equilibrio, perdida de volumen,
      cotas, pico con malla suave, protocolo modificado).
  (10) Modelo de tiempo y memoria (sin FEBio): reproduce sus puntos de
      calibracion (`tiempos`, FEBio) con el error del propio ajuste.

TOLERANCIAS DECLARADAS ANTES DE MEDIR
  (1) 1e-6 relativa al maximo de la tension exacta (extrapolada a amplitud
      nula); con orden erroneo: aborta o error > 1e-2.
  (2) E_app 1e-6 relativa; fuerza 1e-6 relativa.
  (3) bit a bit (==) con pesos iguales; 1e-12 de max|v| la formula frente a
      numpy (se declaro «relativa» al valor, lo que no tiene sentido cerca
      de cero: ver la prueba); el caso a mano exacto a 1e-15.
  (4) igualdad de lineas; area x presion 1e-12 relativa.
  (5) exacto.
  (6) cavidad: volumen de la malla a 1 % del de voxeles (sin la cavidad
      seria +3.4 %); correccion: 0.5 %.
  (7) exacto.
  (8) 1e-6 del maximo de D.
  (9) exacto (estados y motivos).
  (10) memoria: 10 % hex8, 30 % TET10 (residuo maximo del ajuste: 5 y 29 %);
      tamano de la malla TET10: 5 % en estructuras trabeculares (residuo 4 %).
      Los TIEMPOS no se comprueban: se midieron con la CPU al 100 % y son
      provisionales (ver `tiempos`).
"""
import ast
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

from spinpy import febio
from spinpy.escribe import area_caras, escribir_febio, escribir_febio_ensayo
from spinpy.resistencia import (_matriz_D, criterio_pistoia, estadisticos_vm,
                                cuantiles_vm_superficie, percentil_ponderado)
from spinpy.solido import malla_hex, malla_tet10

BLOQUE = "27 FEBio TET10"
REF = "PLAN_Analizar_con_FEBio.md §9"
EXE = febio.localizar()
con_febio = pytest.mark.skipif(EXE is None, reason="FEBio 4 no esta instalado")
VISOR = Path(__file__).resolve().parent.parent / "visor.py"


def _anotar(registro, prueba, esp, obt, err, crit, ok, nota=""):
    registro.anotar(BLOQUE, prueba, REF, esp, obt, err, crit, bool(ok),
                    nota=nota)


# ---------------------------------------------------------------------------
# (1) Orden de nodos
# ---------------------------------------------------------------------------

_C = np.array([[0.1, 0.05, 0.0], [1.3, 0.2, 0.1], [0.3, 1.1, -0.1],
               [0.2, 0.3, 0.9]])
_AR = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
_X = np.vstack([_C] + [0.5 * (_C[a] + _C[b]) for a, b in _AR])


def _campo():
    rng = np.random.default_rng(1)
    L = rng.normal(size=(3, 3))
    Q = rng.normal(size=(3, 3, 3))
    Q = (Q + Q.transpose(0, 2, 1)) / 2
    u = lambda x: x @ L.T + np.einsum("ijk,nj,nk->ni", Q, x, x)  # noqa: E731
    G = L + 2 * np.einsum("ijk,k->ij", Q, _C.mean(0))
    eps = np.array([G[0, 0], G[1, 1], G[2, 2], G[1, 2] + G[2, 1],
                    G[0, 2] + G[2, 0], G[0, 1] + G[1, 0]])
    return u, _matriz_D(20e9, 0.3) @ eps


def _sigma_elemento(orden, d):
    u, _ = _campo()
    s = {}
    for a in (1e-6, 2e-6):
        f = d / f"t_{a:g}.feb"
        febio.escribir_prescrito(_X, np.array([orden]), a * u(_X), f)
        febio.correr(f, exe=EXE)
        s[a] = febio.leer(f)["sigma"][0] / a
    return 2 * s[1e-6] - s[2e-6]


@con_febio
def test_orden_nodos_tet10(registro, tmp_path):
    _, sig = _campo()
    err = float(np.abs(_sigma_elemento(list(range(10)), tmp_path) - sig).max()
                / np.abs(sig).max())
    _anotar(registro, "TET10 en orden C3D10: tension exacta con campo "
            "cuadratico", 0.0, err, err, "1e-6", err < 1e-6)
    assert err < 1e-6
    # Orden erroneo: intermedios de las aristas (1,2) y (0,2) intercambiados
    d = tmp_path / "mal"
    d.mkdir()
    try:
        e2 = float(np.abs(_sigma_elemento([0, 1, 2, 3, 4, 6, 5, 7, 8, 9], d)
                          - sig).max() / np.abs(sig).max())
        detectado = e2 > 1e-2
    except febio.ErrorFEBio:
        e2, detectado = float("inf"), True
    _anotar(registro, "orden erroneo detectado (aborta o error > 1e-2)", 1,
            int(detectado), 0.0, "aborta o > 1e-2", detectado,
            nota=f"error {e2:.3g}")
    assert detectado


# ---------------------------------------------------------------------------
# (2) Bloque macizo
# ---------------------------------------------------------------------------

@con_febio
@pytest.mark.parametrize("tipo", ["hex8", "tet10"])
def test_bloque_macizo(registro, tmp_path, tipo):
    m = febio.mallar(np.ones((5, 5, 7), bool), np.full(3, 0.1), tipo)
    r = febio.ensayo(m, febio.protocolo("app"), tmp_path, ("lineal",),
                     exe=EXE)
    lin = r["lineal"]
    dE = abs(lin["E_app"] / 20e9 - 1.0)
    _anotar(registro, f"bloque macizo {tipo}: E_app = E_s", 20e9,
            lin["E_app"], dE, "1e-6", dE < 1e-6)
    _anotar(registro, f"bloque macizo {tipo}: fuerza = sigma A", 0.0,
            lin["dF_rel"], lin["dF_rel"], "1e-6", lin["dF_rel"] < 1e-6)
    assert dE < 1e-6 and lin["dF_rel"] < 1e-6
    assert m["informe"]["n_superficie"] == 0


# ---------------------------------------------------------------------------
# (3) Estadisticos ponderados
# ---------------------------------------------------------------------------

def _res(n=5000, semilla=3):
    rng = np.random.default_rng(semilla)
    return {"ok": True, "sigma_app": 1e6, "F_total": 2.5e7,
            "eps_eff_solido": rng.gamma(2.0, 1e-4, n),
            "vm_solido": rng.gamma(3.0, 2e6, n),
            "superficie_solido": rng.random(n) < 0.3}


def test_pesos_iguales_bit_a_bit(registro):
    r = _res()
    base_p, base_s = criterio_pistoia(r), estadisticos_vm(r)
    base_q = cuantiles_vm_superficie(r)
    for w in (1.0, 0.0625, 3.7e-4):
        rw = dict(r, vol_solido=np.full(r["vm_solido"].size, w))
        p, s, q = (criterio_pistoia(rw), estadisticos_vm(rw),
                   cuantiles_vm_superficie(rw))
        iguales = (all(p[k] == base_p[k] for k in base_p if k != "ponderado_volumen")
                   and s == base_s and q == base_q)
        _anotar(registro, f"pesos iguales ({w:g}): Pistoia, p99 sup. y "
                "cuantiles bit a bit", 1, int(iguales), 0.0, "==", iguales)
        assert iguales
        assert p["ponderado_volumen"] is False


def test_formula_ponderada_con_pesos_iguales(registro):
    from spinpy.resistencia import _percentil_formula
    v = np.random.default_rng(5).normal(size=777)
    w = np.full(v.size, 0.3)
    # Relativo a la ESCALA de los datos (max |v|), no al valor: la primera
    # version dividia por el valor y en q=50 (valor ~1.5e-5 de una normal
    # centrada) daba 3e-9 con un error absoluto de 4e-14. El criterio mal
    # planteado se corrigio; la tolerancia, 1e-12, es la declarada.
    esc = float(np.abs(v).max())
    peor = 0.0
    for metodo in ("hazen", "linear"):
        for q in (0.0, 1.0, 2.5, 50.0, 98.0, 99.0, 99.9, 100.0):
            a = _percentil_formula(v, w, q, metodo)
            b = float(np.percentile(v, q, method=metodo))
            peor = max(peor, abs(a - b) / esc)
    _anotar(registro, "formula ponderada = numpy con pesos iguales (Q5, Q7)",
            0.0, peor, peor, "1e-12 de max|v|", peor < 1e-12)
    assert peor < 1e-12


def test_ponderado_a_mano(registro):
    # valores [1, 2], pesos [3, 1]: puntos medios 1.5 y 3.5 de W = 4.
    # hazen: posiciones 0.375 y 0.875 -> q=50 da 1 + 0.125/0.5 = 1.25.
    # tipo 7: posiciones 0 y 1 -> q=50 da 1.5; q=25 da 1.25.
    v, w = np.array([2.0, 1.0]), np.array([1.0, 3.0])
    casos = [(percentil_ponderado(v, w, 50, "hazen"), 1.25),
             (percentil_ponderado(v, w, 50, "linear"), 1.5),
             (percentil_ponderado(v, w, 25, "linear"), 1.25),
             (percentil_ponderado(v, w, 10, "hazen"), 1.0)]
    err = max(abs(a - b) for a, b in casos)
    _anotar(registro, "percentil ponderado calculado a mano", 0.0, err, err,
            "1e-15", err < 1e-15)
    assert err < 1e-15


def test_pistoia_cuenta_volumen(registro):
    """Con pesos distintos, la fraccion de VOLUMEN sobre el umbral es ~2 %."""
    rng = np.random.default_rng(9)
    n = 20000
    e = rng.gamma(2.0, 1e-4, n)
    w = rng.uniform(0.1, 3.0, n)
    r = {"ok": True, "sigma_app": 1e6, "F_total": 1.0,
         "eps_eff_solido": e, "vm_solido": e * 1e10, "vol_solido": w,
         "superficie_solido": np.ones(n, bool)}
    p = criterio_pistoia(r)
    frac = float(w[e > p["eps_eff_p"]].sum() / w.sum())
    err = abs(frac - 0.02)
    # Un elemento pesa como mucho 3/(0.1 n): la fraccion es exacta a eso.
    tol = 3.0 / w.sum() * 2
    _anotar(registro, "Pistoia ponderado: 2 % del VOLUMEN sobre el umbral",
            0.02, frac, err, "peso de un elemento", err < tol)
    assert err < tol and p["ponderado_volumen"] is True


# ---------------------------------------------------------------------------
# (4) Escritores
# ---------------------------------------------------------------------------

def _giroide(n, nivel=-0.6):
    x = np.linspace(0.0, 1.0, n)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    f = (np.sin(4 * np.pi * X) * np.cos(4 * np.pi * Y)
         + np.sin(4 * np.pi * Y) * np.cos(4 * np.pi * Z)
         + np.sin(4 * np.pi * Z) * np.cos(4 * np.pi * X))
    return f > nivel


def _sin_comentario(ruta):
    txt = Path(ruta).read_text(encoding="utf-8")
    i, j = txt.index("<!--"), txt.index("-->") + 3
    return (txt[:i] + txt[j:]).splitlines()


@pytest.mark.parametrize("apoyo,plato", [("deslizante", None),
                                         ("empotrado", None),
                                         ("deslizante", 1e-4)])
def test_escritor_hex8_igual_al_validado(registro, tmp_path, apoyo, plato):
    from spinpy.resistencia import _solo_portante
    BW = _solo_portante(_giroide(10))
    nodos, elems, _ = malla_hex(BW, np.full(3, 0.1))
    # Mismo nombre en dos carpetas: el nombre del .feb entra en los logfile.
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    a, b = tmp_path / "a" / "g.feb", tmp_path / "b" / "g.feb"
    escribir_febio(nodos, elems, a, sigma_app=1e6, A_bruta=1.0, apoyo=apoyo,
                   eps_plato=plato)
    escribir_febio_ensayo(nodos, elems, b, sigma_app=1e6, A_bruta=1.0,
                          apoyo=apoyo, eps_plato=plato)
    iguales = _sin_comentario(a) == _sin_comentario(b)
    _anotar(registro, f"escritor general = validado (hex8, {apoyo}, "
            f"{'plato' if plato else 'fuerza'})", 1, int(iguales), 0.0,
            "lineas iguales", iguales)
    assert iguales


def test_escritor_tet10(registro, tmp_path):
    m = febio.mallar(np.ones((4, 4, 5), bool), np.full(3, 0.1), "tet10")
    ruta = tmp_path / "t.feb"
    inf = escribir_febio_ensayo(m["nodos"], m["elems"], ruta, sigma_app=1e6,
                                A_bruta=m["A_bruta"],
                                caras_techo=m["caras_techo"])
    raiz = ET.parse(ruta).getroot()
    caras = raiz.findall("Mesh/Surface/tri6")
    p = float(raiz.find("Loads/surface_load/pressure").text)
    c = m["caras_techo"]
    P = m["nodos"][c[:, :3]]
    nz = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])[:, 2]
    F = p * float(area_caras(m["nodos"], c).sum())
    errF = abs(F - m["A_bruta"]) / m["A_bruta"]          # 1 MPa * A (N)
    dom = raiz.find("MeshDomains/SolidDomain")
    ok = (len(caras) == c.shape[0] and bool((nz > 0).all()) and errF < 1e-12
          and dom.get("elem_type") == "TET10G8"
          and raiz.find("Output/plotfile") is not None
          # piso de redondeo del residuo con TET10 (ver escribe.RTOL)
          and float(raiz.find("Control/solver/rtol").text) == 1e-10
          and raiz.find("Mesh/Elements").get("type") == "tet10")
    _anotar(registro, "escritor TET10: tri6 hacia +z, F = sigma A, regla "
            "como atributo, plotfile", 1, int(ok), errF, "1e-12 / exacto", ok)
    assert ok and inf["tipo"] == "tet10"


# ---------------------------------------------------------------------------
# (5) Lectura
# ---------------------------------------------------------------------------

def test_leer_logfile(registro, tmp_path):
    u = ("*Step  = 0\n*Time  = 0\n*Data  = ux;uy;uz\n1 0 0 0\n2 0 0 0\n"
         "*Step  = 1\n*Time  = 1\n*Data  = ux;uy;uz\n2 4 5 6\n1 1 2 3\n")
    s = ("*Step  = 0\n*Time  = 0\n*Data  = sx;sy;sz;syz;sxz;sxy\n"
         "1 0 0 0 0 0 0\n*Step  = 1\n*Time  = 1\n"
         "*Data  = sx;sy;sz;syz;sxz;sxy\n1 1 2 3 4 5 6\n")
    (tmp_path / "m_u.txt").write_text(u)
    (tmp_path / "m_s.txt").write_text(s)
    r = febio.leer(tmp_path / "m.feb")
    ok = (np.array_equal(r["u"], [[1, 2, 3], [4, 5, 6]])
          and np.array_equal(r["sigma"], [[1e6, 2e6, 3e6, 4e6, 5e6, 6e6]]))
    _anotar(registro, "leer: ultimo paso, ordenado, MPa -> Pa", 1, int(ok),
            0.0, "exacto", ok)
    assert ok


# ---------------------------------------------------------------------------
# (6) Malla
# ---------------------------------------------------------------------------

def _cavidad(n=16, r_rel=0.2):
    h = 1.0 / n
    c = np.arange(n) * h
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    o = (n / 2.0 - 0.5) * h
    return np.sqrt((X - o) ** 2 + (Y - o) ** 2 + (Z - o) ** 2) > r_rel, \
        np.full(3, h)


def test_cavidad_cerrada_no_se_rellena(registro):
    BW, spc = _cavidad()
    V = float(BW.sum() * np.prod(spc))
    _, _, _, inf = malla_tet10(BW, spc)
    err = abs(inf["volumen"] / V - 1.0)
    _anotar(registro, "cavidad cerrada: la malla no la rellena", V,
            inf["volumen"], err, "1 % (rellena: +3.4 %)", err < 0.01,
            nota=f"poros quitados {inf['poros']['n_poros_quitados']}")
    assert err < 0.01 and inf["poros"]["n_poros_quitados"] == 1


def test_correccion_de_volumen(registro):
    from spinpy import generar_mascara
    B, _, _ = generar_mascara(resolution=24, wave_number=6 * np.pi,
                              num_waves=200, thetas=[90, 90, 90], rho=0.4,
                              seed=2)
    B = np.asarray(B, bool)
    m = febio.mallar(B, np.full(3, 0.1), "tet10")
    err = abs(m["informe"]["vol_pct_voxel"] - 100.0) / 100.0
    _anotar(registro, "correccion de volumen: malla = voxeles portantes",
            100.0, m["informe"]["vol_pct_voxel"], err, "0.5 %", err < 0.005)
    assert err < 0.005
    assert m["informe"]["n_superficie"] > 0


# ---------------------------------------------------------------------------
# (7) Protocolo de Tapia
# ---------------------------------------------------------------------------

def _constantes_visor():
    arbol = ast.parse(VISOR.read_text(encoding="utf-8"))
    out = {}
    for n in arbol.body:
        if isinstance(n, ast.Assign) and len(n.targets) == 1 and \
                isinstance(n.targets[0], ast.Name) and \
                n.targets[0].id.startswith("PAPER_"):
            try:
                out[n.targets[0].id] = ast.literal_eval(n.value)
            except ValueError:
                pass
    return out


@con_febio
def test_homogeneizacion_bloque(registro, tmp_path):
    D = _matriz_D(20e9, 0.3)
    r = febio.homogeneizar(np.ones((3, 3, 3), bool), np.full(3, 0.1),
                           febio.protocolo("homogeneizacion"), malla="tet10",
                           carpeta=tmp_path, n=3, exe=EXE, comparar_app=False)
    for k in ("C_KUBC", "C_SUBC"):
        err = float(np.abs(r[k] - D).max() / np.abs(D).max())
        _anotar(registro, f"bloque macizo TET10: {k} = D isotropa", 0.0, err,
                err, "1e-6 del max de D", err < 1e-6)
        assert err < 1e-6
    assert not r["fallos"]


def _doc(*registros):
    return {"resultados": {"febio": {"registros": list(registros)}},
            "procedencia": {"version_formato": 99}}


def test_calificacion_febio(registro):
    from spinpy import informe
    base = {"estructura": "voi", "nombre": "App", "eje_nombre": "Z",
            "fallos": []}
    lin = {"E_app": 1e9, "dF_rel": 1e-9,
           "pistoia": {"ok": True, "sigma_fallo": 5e6, "vm_max": 1e8,
                       "vm_p99_superficie": 5e7, "vm_n_superficie": 5000}}
    hex_ok = dict(base, malla="hex8", n=40, lineal=lin,
                  informe_malla={"frac_portante_voxel": 1.0})
    tet_perd = dict(base, malla="tet10", lineal=lin,
                    informe_malla={"frac_portante_voxel": 1.0,
                                   "perdida_pct": 5.0, "descartado_pct": 0.0})
    sin_lin = dict(base, malla="hex8", n=40, lineal=None,
                   fallos=[{"corrida": "lineal_0.001"}],
                   informe_malla={"frac_portante_voxel": 1.0})
    desequil = dict(hex_ok, lineal=dict(lin, dF_rel=1e-3))
    homog = dict(base, tipo="homogeneizacion", malla="tet10",
                 informe_malla={"frac_portante_voxel": 1.0})
    modif = dict(hex_ok, modificado=True, nombre="Tapia (modificado)")

    def estados(r):
        return {(i["magnitud"]): (i["estado"], set(i["motivos"]))
                for i in informe.comprobar(_doc(r))}

    e = estados(hex_ok)
    casos = [
        ("hex8 limpio: E_app citable", e["febio_E_app"][0] == "citable"),
        ("hex8 limpio: vm_max no citable", e["vm_max"][0] == "no_citable"),
        ("hex8 limpio: p99 superficie citable",
         e["vm_p99_superficie"][0] == "citable"),
        ("TET10 con perdida 5 %: reservas",
         "febio_perdida_volumen" in estados(tet_perd)["febio_E_app"][1]),
        ("TET10: p99 superficie con reservas por el suavizado",
         "febio_suavizado_pico" in estados(tet_perd)["vm_p99_superficie"][1]),
        ("sin lineal y con fallos: no citable",
         estados(sin_lin)["febio_E_app"][0] == "no_citable"),
        ("fuerza desequilibrada: no citable",
         estados(desequil)["febio_E_app"][0] == "no_citable"),
        ("homogeneizacion TET10: cotas con reservas",
         estados(homog)["febio_homog"] == ("reservas", {"febio_cotas"})),
        ("protocolo modificado: reservas",
         "febio_protocolo_modificado"
         in estados(modif)["febio_E_app"][1]),
    ]
    for nombre, ok in casos:
        _anotar(registro, "calificacion: " + nombre, 1, int(ok), 0.0,
                "exacto", ok)
    assert all(ok for _n, ok in casos), [n for n, ok in casos if not ok]
    # Las etiquetas distinguen protocolo y malla.
    it = informe.comprobar(_doc(tet_perd))[0]
    assert "FEBio App, tet10" in informe.etiqueta_item(it)


# (gdl, memoria MB) medidos: cavidad.jsonl, 2026-09-24.
_MEM_HEX8 = [(14568, 133.4), (46032, 464.2), (105558, 1236.6),
             (201966, 2739.6)]
_MEM_TET10 = [(204492, 2099.9), (468441, 4731.8), (33267, 221.2),
              (79722, 526.0), (155202, 966.4), (227211, 1446.7),
              (120747, 877.7), (128127, 924.9), (186525, 1407.4),
              (190398, 1484.4), (234570, 1604.3), (269409, 1800.9),
              (262008, 1753.7), (336450, 2523.4), (397155, 2733.8),
              (367623, 2628.5)]
# (voxeles de superficie, voxeles de hueso, TET10 medidos): VOI H4 48 y 64,
# spinodoide 32^3.
_TET_TRAB = [(24509, 31129, 178929), (49284, 73381, 350378),
             (6574, 11627, 46708)]


def test_modelo_memoria_y_malla(registro):
    from spinpy import tiempos
    peor = {}
    for tipo, pts, tol in (("hex8", _MEM_HEX8, 0.10),
                           ("tet10", _MEM_TET10, 0.30)):
        e = max(abs(tiempos.memoria_febio(tipo, g) / m - 1) for g, m in pts)
        peor[tipo] = e
        _anotar(registro, f"memoria FEBio {tipo} en sus puntos medidos", 0.0,
                e, e, f"{tol:.0%}", e < tol)
    e = max(abs(tiempos.n_tet10(s, v) / t - 1) for s, v, t in _TET_TRAB)
    _anotar(registro, "numero de TET10 previsto (trabecular)", 0.0, e, e,
            "5 %", e < 0.05)
    assert peor["hex8"] < 0.10 and peor["tet10"] < 0.30 and e < 0.05


def test_tapia_igual_al_visor(registro):
    c = _constantes_visor()
    t = febio.PROTOCOLOS_FEBIO["tapia2026"]
    ok = (t["E_s"] == c["PAPER_E_S"] and t["nu"] == c["PAPER_NU"]
          and t["carga_N"] == c["PAPER_CARGA_N"]
          and t["apoyo"] == c["PAPER_APOYO"])
    p0 = febio.protocolo("tapia2026")
    p1 = febio.protocolo("tapia2026", E_s=17e9)
    ok = ok and p0["nombre"] == "Tapia" and p1["nombre"] == "Tapia (modificado)"
    _anotar(registro, "protocolo Tapia = visor.PAPER_*, editado se renombra",
            1, int(ok), 0.0, "exacto", ok)
    assert ok
