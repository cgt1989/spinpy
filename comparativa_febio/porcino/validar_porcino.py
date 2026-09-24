"""validar_porcino.py — Los analisis mecanicos de la app, contra FEBio 4.

    cd Port_Python
    python comparativa_febio/porcino/validar_porcino.py              # todo
    python comparativa_febio/porcino/validar_porcino.py compresion comparado
    python comparativa_febio/porcino/validar_porcino.py homogeneizacion \
        --estructuras voi

Sobre el VOI porcino V1 y los DOS candidatos ajustados a el en la sesion del
informe automatico del 23-09-2026 (`Informes_Porcino/VOI_V1_ambas_2026-09-23`),
repite con spinpy cada analisis mecanico que la app propone y resuelve el MISMO
problema discreto con FEBio 4.5:

  compresion       ensayo de la app: 40^3, apoyo deslizante, E_s 20 GPa,
                   1 MPa sobre la seccion bruta, eje Z; Pistoia y von Mises
  comparado        protocolo de Tapia et al. (2026): 40^3, empotrado,
                   E_s 18 GPa, 100 N
  homogeneizacion  tensor periodico a 32^3 (restricciones lineales en FEBio)
  convergencia     E_app a 22, 28, 34 y 40^3, las resoluciones de la sesion
  fallo            fallo progresivo (Simulaciones in silico), 32^3, 10 pasos,
                   con el bucle entero repetido usando FEBio como solver

LAS ESTRUCTURAS NO SE REDIBUJAN: el VOI se lee del .mat de la sesion y se
comprueba su SHA-256; los candidatos se regeneran desde `reproduccion.json` por
`informe._regenerar` (la unica traduccion parametros -> generador) y se
comprueba su huella bit a bit. Si algo no casa, el script se detiene.

LECTURA DE LOS NUMEROS — lo mismo que en `comparar_febio.py`:

* FEBio es no lineal geometricamente. El problema LINEAL se contrasta con la
  extrapolacion a carga nula u_lin = 2 u(s1) - u(2 s1), reescalada a la carga
  del ensayo; el termino no lineal de primer orden se cancela.
* Ademas se resuelve FEBio a la carga del ensayo, y en el ensayo de la app
  tambien a la carga de fallo de Pistoia, con fuerza impuesta (como la app) y
  con un plato rigido: eso no valida ni invalida a spinpy, mide hasta donde
  se sostiene la hipotesis lineal.
* Los estadisticos (p99 de superficie, Pistoia) se calculan con las MISMAS
  funciones de spinpy sobre los dos campos.

Salidas: `resultados/<analisis>.jsonl` (una linea por estructura, o por punto),
`campos/<analisis>_<estructura>.npz` (campos por elemento para los mapas de
color) y `corridas/` (los .feb y las salidas de FEBio; no se versionan).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
PORT = AQUI.parent.parent
sys.path.insert(0, str(PORT))
sys.path.insert(0, str(AQUI.parent))
sys.path.insert(0, str(AQUI))

from comparar_febio import FEBIO, _rel, _ultimo_registro, _von_mises, correr_febio  # noqa: E402
from homog_febio import escribir_periodico                            # noqa: E402
from spinpy import __version__ as VERSION_SPINPY                      # noqa: E402
from spinpy import leer_voi, procedencia                              # noqa: E402
from spinpy.elastic import (constantes_ingenieria, distancia_log_euclidea,  # noqa: E402
                            homogeneizar, remuestrear_bw)
from spinpy.escribe import escribir_febio                             # noqa: E402
from spinpy.informe import _regenerar, huella_mascara, sha256_archivo  # noqa: E402
from spinpy.resistencia import (FRAC_CRITICA, EPS_CRITICA,           # noqa: E402
                                _solo_portante, criterio_pistoia,
                                cuantiles_vm_superficie, ensayo_compresion)
from spinpy.simulacion import RIGIDEZ_DANADA, fallo_progresivo       # noqa: E402
from spinpy.solido import malla_hex                                  # noqa: E402

# Copia de `resultados_sesion.json` y `reproduccion.json` de la sesion del
# informe automatico (Informes_Porcino/VOI_V1_ambas_2026-09-23), para que el
# estudio viaje entero. El VOI no se copia: es de Koria, Mengoni y Brockett
# (2020, Research Data Leeds, doi:10.5518/787, CC BY 4.0) y se busca en la ruta
# que guardo la sesion, o en la que se pase con --voi; su SHA-256 se comprueba.
SESION = AQUI / "sesion"
RUTA_VOI = None
ESTRUCTURAS = ("voi", "spinodoide", "dual-lattice")
# Donde guarda la sesion cada estructura: (sufijo de la clave, clave interna)
EN_SESION = {"voi": ("", "voi"), "spinodoide": ("", "spin"),
             "dual-lattice": ("_dual", "spin")}

NU = 0.30
APP = dict(n=40, apoyo="deslizante", E_s=20e9, sigma=1e6)
TAPIA = dict(n=40, apoyo="empotrado", E_s=18e9, carga_N=100.0)
N_HOMOG = 32
N_CONV = (22, 28, 34, 40)
FALLO = dict(n_mec=32, pasos=10)
SIGMA_LIN = (1e3, 2e3)        # Pa: dos cargas pequenas para extrapolar a cero
EPS_PLATO_LIN = (1e-6, 2e-6)  # idem, en deformacion de plato
AMP_HOMOG = (1e-5, 2e-5)      # idem, en deformacion macroscopica

DIR_RES = AQUI / "resultados"
DIR_CAMPOS = AQUI / "campos"
DIR_CORR = AQUI / "corridas"
DIR_CACHE = AQUI / "cache"


# ---------------------------------------------------------------------------
# Estructuras y sesion
# ---------------------------------------------------------------------------

def sesion():
    return json.load(open(SESION / "resultados_sesion.json", encoding="utf-8"))


def reproduccion():
    return json.load(open(SESION / "reproduccion.json", encoding="utf-8"))


def cargar(est):
    """(BW, spacing) de la estructura, con la huella comprobada."""
    rep = reproduccion()
    DIR_CACHE.mkdir(parents=True, exist_ok=True)
    if est == "voi":
        ruta = RUTA_VOI or rep["voi"]["ruta"]
        if sha256_archivo(ruta) != rep["voi"]["sha256"]:
            raise RuntimeError(f"El VOI {ruta} no es el de la sesion")
        BW, sp = leer_voi(ruta)
        return np.asarray(BW, bool), np.asarray(sp, float)
    a = rep["ajustes"][est]
    cache = DIR_CACHE / f"{est}.npy"
    if cache.exists():
        BW = np.load(cache)
    else:
        BW = np.asarray(_regenerar(a["parametros"]), bool)
    if huella_mascara(BW) != a["sha256_mascara"]:
        raise RuntimeError(f"La mascara regenerada de {est} no es la de la "
                           "sesion (huella distinta)")
    np.save(cache, BW)
    return BW, np.full(3, float(a["parametros"]["spacing_mm"]))


def _registro_sesion(doc, analisis, est):
    suf, k = EN_SESION[est]
    r = doc["resultados"].get(analisis + suf) or {}
    return (r.get("por_estructura") or {}).get(k)


def _escribir(analisis, fila):
    DIR_RES.mkdir(parents=True, exist_ok=True)
    fila["procedencia"] = procedencia.bloque(
        fila.get("estructura", "?"), semilla=20260720)
    fila["procedencia"]["febio"] = fila.get("version_febio")
    fila["fecha"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    with open(DIR_RES / f"{analisis}.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(fila, ensure_ascii=False, default=float) + "\n")


# ---------------------------------------------------------------------------
# FEBio
# ---------------------------------------------------------------------------

def _dos_materiales(feb, modulo_elem, E_s):
    """Reescribe un .feb de `escribir_febio` con el tejido ablandado aparte.

    `modulo_elem` (n_elem,) es la rigidez relativa de cada elemento en el
    orden de la malla; los elementos con modulo < 1 pasan a un segundo
    material de E = modulo * E_s. Solo se admite UN valor de ablandamiento,
    que es lo que hace `fallo_progresivo`.
    """
    blandos = modulo_elem < 1.0
    if not blandos.any():
        return
    vals = np.unique(modulo_elem[blandos])
    if vals.size != 1:
        raise ValueError("solo se admite un nivel de ablandamiento")
    txt = feb.read_text(encoding="utf-8")
    ini = txt.index('\t\t<Elements type="hex8" name="solido">\n')
    fin = txt.index("\t\t</Elements>\n", ini) + len("\t\t</Elements>\n")
    lineas = txt[ini:fin].splitlines(keepends=True)[1:-1]
    sano = "".join(l for l, b in zip(lineas, blandos) if not b)
    blando = "".join(l for l, b in zip(lineas, blandos) if b)
    bloque = ('\t\t<Elements type="hex8" name="solido">\n' + sano
              + "\t\t</Elements>\n"
              + '\t\t<Elements type="hex8" name="danado">\n' + blando
              + "\t\t</Elements>\n")
    txt = txt[:ini] + bloque + txt[fin:]
    E2 = float(vals[0]) * E_s / 1e6
    txt = txt.replace(
        "\t\t</material>\n\t</Material>\n",
        "\t\t</material>\n\t\t<material id=\"2\" name=\"danado\" "
        "type=\"isotropic elastic\">\n"
        f"\t\t\t<E>{E2:.17g}</E>\n\t\t\t<v>{NU:.9g}</v>\n"
        "\t\t</material>\n\t</Material>\n", 1)
    txt = txt.replace(
        '\t\t<SolidDomain name="solido" mat="hueso"/>\n',
        '\t\t<SolidDomain name="solido" mat="hueso"/>\n'
        '\t\t<SolidDomain name="danado" mat="danado"/>\n', 1)
    feb.write_text(txt, encoding="utf-8")


class Caso:
    """Una malla de spinpy lista para resolverse en FEBio con varias cargas."""

    def __init__(self, BW, spc, apoyo, E_s, carpeta, nombre, modulo=None):
        self.BW = np.asarray(BW, bool)
        self.spc = np.asarray(spc, float)
        self.apoyo, self.E_s = apoyo, float(E_s)
        self.dir = Path(carpeta)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.nombre = nombre
        port = _solo_portante(self.BW)
        self.nodos, self.elems, _ = malla_hex(port, self.spc)
        self.ijk = np.stack(np.nonzero(port), axis=1)
        self.modulo = (np.ones(self.elems.shape[0]) if modulo is None
                       else np.asarray(modulo, float)[tuple(self.ijk.T)])
        nx, ny, _nz = self.BW.shape
        self.A_bruta = float(nx * self.spc[0] * ny * self.spc[1])
        z = self.nodos[:, 2]
        self.H = float(z.max() - z.min())
        self.techo = z >= z.max() - 1e-9 * max(self.H, 1.0)
        self.V_e = float(np.prod(self.spc))
        self.version = None
        self.tiempos = {}

    def correr(self, etiqueta, sigma=None, eps_plato=None):
        """Corre FEBio. Devuelve (u mm, sigma Pa) sin reescalar, o None."""
        feb = self.dir / f"{self.nombre}_{etiqueta}.feb"
        escribir_febio(self.nodos, self.elems, feb, E_s=self.E_s, nu_s=NU,
                       sigma_app=1e6 if sigma is None else sigma,
                       A_bruta=self.A_bruta, apoyo=self.apoyo,
                       eps_plato=eps_plato)
        _dos_materiales(feb, self.modulo, self.E_s)
        try:
            t, ver = correr_febio(feb)
        except RuntimeError as e:
            print(f"   FEBio {etiqueta}: {e}")
            return None
        self.version = ver
        self.tiempos[etiqueta] = t
        u = _ultimo_registro(feb.with_name(feb.stem + "_u.txt"), 3).ravel()
        s = _ultimo_registro(feb.with_name(feb.stem + "_s.txt"), 6) * 1e6
        return u, s

    # --- magnitudes derivadas, con las definiciones de spinpy ---------------
    def E_app(self, u, sigma):
        eps = abs(float(np.mean(u[2::3][self.techo]))) / self.H
        return sigma / eps

    def fuerza(self, s):
        """Fuerza en Pa*mm^2 por la integral de volumen de sigma_zz."""
        return float(-s[:, 2].sum() * self.V_e / self.H)

    def eps_eff(self, s):
        """Deformacion efectiva de Pistoia desde la tension, elemento a
        elemento con SU modulo, y dividida por E_s como en spinpy."""
        E_e = self.E_s * self.modulo
        S = np.zeros((self.modulo.size, 6, 6))
        S[:, :3, :3] = (-NU / E_e)[:, None, None]
        for k in range(3):
            S[:, k, k] = 1.0 / E_e
            S[:, 3 + k, 3 + k] = 2.0 * (1.0 + NU) / E_e
        eps = np.einsum("nij,nj->ni", S, s)
        U = 0.5 * np.einsum("ni,ni->n", s, eps)
        return np.sqrt(np.maximum(2.0 * U / self.E_s, 0.0))

    def lineal(self, sigma_ref):
        """Campo lineal extrapolado a carga nula, reescalado a sigma_ref."""
        r1 = self.correr("lin1", sigma=SIGMA_LIN[0])
        r2 = self.correr("lin2", sigma=SIGMA_LIN[1])
        if r1 is None or r2 is None:
            raise RuntimeError(f"FEBio no convergio a carga pequena en "
                               f"{self.nombre}")
        (u1, s1), (u2, s2) = r1, r2
        e1 = sigma_ref / SIGMA_LIN[0]
        u = (2 * u1 - u2 / 2) * e1
        s = (2 * s1 - s2 / 2) * e1
        desvio = float(np.abs(u1 * e1 - u).max() / np.abs(u).max())
        return u, s, desvio

    def plato_lineal(self):
        r1 = self.correr("plato1", eps_plato=EPS_PLATO_LIN[0])
        r2 = self.correr("plato2", eps_plato=EPS_PLATO_LIN[1])
        if r1 is None or r2 is None:
            return None
        (u1, s1), (u2, s2) = r1, r2
        s = (2 * s1 - s2 / 2) / EPS_PLATO_LIN[0]      # tension por unidad eps
        return self.fuerza(s) / self.A_bruta          # = E_app del plato


def res_febio(caso, res_s, u, s):
    """Diccionario tipo `ensayo_compresion` a partir del campo de FEBio, para
    pasarlo por las MISMAS funciones de spinpy (Pistoia, percentiles)."""
    solido = caso.modulo >= 1.0
    vm = _von_mises(s)
    ee = caso.eps_eff(s)
    return {"ok": True, "sigma_app": res_s["sigma_app"],
            "F_total": res_s["F_total"],
            "E_app": caso.E_app(u, res_s["sigma_app"]),
            "eps_eff_solido": ee[solido], "vm_solido": vm[solido],
            "superficie_solido": res_s["superficie_solido"],
            "_vm": vm, "_ee": ee}


def _desp_elem(caso, u):
    return np.linalg.norm(u.reshape(-1, 3), axis=1)[caso.elems].mean(axis=1)


# ---------------------------------------------------------------------------
# Ensayo de compresion (app) y analisis comparado (Tapia)
# ---------------------------------------------------------------------------

CLAVES_PISTOIA = ("vm_p99_superficie", "factor", "sigma_fallo", "vm_max",
                  "vm_p99", "eps_eff_p", "vm_media")


def ensayo(est, analisis, BW, spc, cfg, doc):
    """Un ensayo entero: spinpy, FEBio lineal, FEBio no lineal, plato."""
    n = cfg["n"]
    bw, spr = remuestrear_bw(BW, spc, n)
    fila = {"estructura": est, "analisis": analisis, "n": n,
            "apoyo": cfg["apoyo"], "E_s_Pa": cfg["E_s"],
            "BVTV": float(bw.mean()), "spinpy": VERSION_SPINPY}
    print(f"\n== {analisis} / {est}: {bw.shape}, BV/TV {bw.mean():.4f}")

    t0 = time.perf_counter()
    res = ensayo_compresion(bw, spr, E_s=cfg["E_s"], nu_s=NU,
                            sigma0=cfg.get("sigma", 1e6), apoyo=cfg["apoyo"],
                            carga_N=cfg.get("carga_N"), unidad="mm")
    fila["t_spinpy_s"] = time.perf_counter() - t0
    if not res["ok"]:
        raise RuntimeError(f"spinpy no resolvio {est}: {res['msg']}")
    sigma = float(res["sigma_app"])
    fila.update(sigma_app_Pa=sigma, F_N=res["F_total"] * 1e-6,
                n_elem=res["n_elem"], n_dof=res["n_dof"],
                frac_portante=res["frac_portante"],
                solver_spinpy=res.get("solver"),
                residuo_spinpy=res.get("residuo_rel"))

    caso = Caso(bw, spr, cfg["apoyo"], cfg["E_s"],
                DIR_CORR / analisis / est, f"{analisis}_{est}")
    if 3 * caso.nodos.shape[0] != res["n_dof"] or \
            caso.elems.shape[0] != res["n_elem"]:
        raise RuntimeError("La malla exportada no es la que resolvio spinpy")

    # --- problema lineal ---------------------------------------------------
    u_f, s_f, desvio = caso.lineal(sigma)
    fila["desvio_nl_carga_pequena"] = desvio
    fila["version_febio"] = caso.version
    E_s_app = res["E_app"]
    E_f = caso.E_app(u_f, sigma)
    fila.update(E_app_spinpy_MPa=E_s_app / 1e6, E_app_febio_MPa=E_f / 1e6,
                dE_rel=_rel(E_f, E_s_app))
    fila["F_febio_N"] = caso.fuerza(s_f) * 1e-6
    fila["dF_rel"] = _rel(caso.fuerza(s_f), res["F_total"])
    fila["du_max_rel"] = float(np.abs(u_f - res["u"]).max()
                               / np.abs(res["u"]).max())
    fila["dsig_max_rel"] = float(np.abs(s_f - res["sig"]).max()
                                 / res["vm"].max())
    rf = res_febio(caso, res, u_f, s_f)
    fila["dvm_max_elem_rel"] = float(np.abs(rf["_vm"] - res["vm"]).max()
                                     / res["vm"].max())
    d_s, d_f = _desp_elem(caso, res["u"]), _desp_elem(caso, u_f)
    fila["ddesp_max_rel"] = float(np.abs(d_f - d_s).max() / d_s.max())
    fila["desp_max_spinpy_mm"] = float(res["desp_max"])
    fila["desp_max_febio_mm"] = float(np.linalg.norm(u_f.reshape(-1, 3),
                                                     axis=1).max())
    p_s, p_f = criterio_pistoia(res), criterio_pistoia(rf)
    for k in CLAVES_PISTOIA:
        fila[f"{k}_spinpy"] = p_s.get(k)
        fila[f"{k}_febio"] = p_f.get(k)
        fila[f"d_{k}_rel"] = _rel(p_f[k], p_s[k])
    fila["vm_n_superficie"] = p_s.get("vm_n_superficie")
    fila["F_fallo_spinpy_N"] = p_s["F_fallo"] * 1e-6
    fila["F_fallo_febio_N"] = p_f["F_fallo"] * 1e-6
    q_s = cuantiles_vm_superficie(res, escala=1e-6)
    q_f = cuantiles_vm_superficie(rf, escala=1e-6)
    fila["cuantiles_vm_superficie"] = {"p": q_s.get("p"),
                                       "spinpy_MPa": q_s.get("valor"),
                                       "febio_MPa": q_f.get("valor")}

    # --- lo que guardo la sesion: ¿la estructura regenerada es la misma? ----
    reg = _registro_sesion(doc, "resistencia" if analisis == "compresion"
                           else "analisis_comparado", est)
    if reg is not None:
        if analisis == "compresion":
            z = reg["ejes"]["Z"]
            fila["sesion"] = {"E_app_MPa": z["E_app"] / 1e6,
                              "factor": z["factor"],
                              "vm_p99_superficie": z["vm_p99_superficie"]}
            fila["d_sesion_E_rel"] = _rel(E_s_app, z["E_app"])
            fila["d_sesion_factor_rel"] = _rel(p_s["factor"], z["factor"])
        else:
            fila["sesion"] = {"E_app_MPa": reg["E_app"] / 1e6,
                              "vm_p99_superficie_MPa":
                                  reg["vm_p99_superficie"]}
            fila["d_sesion_E_rel"] = _rel(E_s_app, reg["E_app"])
            fila["d_sesion_p99_rel"] = _rel(p_s["vm_p99_superficie"] / 1e6,
                                            reg["vm_p99_superficie"])

    campos = {"ijk": caso.ijk.astype(np.int16), "forma": np.array(bw.shape),
              "spacing": spr, "sigma_app": sigma,
              "vm_spinpy": res["vm"].astype(np.float32),
              "vm_febio": rf["_vm"].astype(np.float32),
              "ee_spinpy": res["eps_eff"].astype(np.float32),
              "ee_febio": rf["_ee"].astype(np.float32),
              "desp_spinpy": d_s.astype(np.float32),
              "desp_febio": d_f.astype(np.float32),
              "superficie": np.asarray(res["superficie_solido"], bool)}

    # --- no lineal: a la carga del ensayo y a la de fallo de Pistoia -------
    cargas = {"ensayo": sigma}
    if analisis == "compresion":
        cargas["fallo"] = float(p_s["sigma_fallo"])
    for nom, sg in cargas.items():
        r = caso.correr(f"nl_{nom}", sigma=sg)
        if r is None:
            fila[f"nl_{nom}"] = {"sigma_Pa": sg, "converge": False}
            continue
        u_nl, s_nl = r
        E_nl = caso.E_app(u_nl, sg)
        rn = res_febio(caso, res, u_nl * sigma / sg, s_nl * sigma / sg)
        ee_nl = caso.eps_eff(s_nl)[caso.modulo >= 1.0]
        pn = criterio_pistoia(rn)
        fila[f"nl_{nom}"] = {
            "sigma_Pa": sg, "converge": True,
            "E_app_MPa": E_nl / 1e6, "dE_rel": E_nl / E_s_app - 1.0,
            "du_max_rel": float(np.abs(u_nl * sigma / sg - res["u"]).max()
                                / np.abs(res["u"]).max()),
            "d_p99_sup_rel": (pn["vm_p99_superficie"]
                              / p_s["vm_p99_superficie"] - 1.0),
            # a la carga de fallo, ¿que fraccion del tejido pasa de 0.7 %?
            # En lineal es EXACTAMENTE el 2 % por construccion.
            "frac_sobre_eps_crit": float(np.mean(ee_nl >= EPS_CRITICA)),
            "desp_max_mm": float(np.linalg.norm(u_nl.reshape(-1, 3),
                                                axis=1).max()),
        }
        if nom == "ensayo":
            campos["vm_febio_nl"] = (rn["_vm"]).astype(np.float32)
            campos["desp_febio_nl"] = (_desp_elem(caso, u_nl)
                                       ).astype(np.float32)
        print(f"   no lineal a {sg / 1e6:.3g} MPa: dE "
              f"{fila[f'nl_{nom}']['dE_rel']:+.3e}  frac>0.7% "
              f"{fila[f'nl_{nom}']['frac_sobre_eps_crit']:.4f}")

    # --- plato rigido (solo el ensayo de la app) ---------------------------
    if analisis == "compresion":
        E_pl = caso.plato_lineal()
        if E_pl is not None:
            fila["E_app_plato_MPa"] = E_pl / 1e6
            fila["plato_sobre_fuerza"] = E_pl / E_s_app
            eps_f = float(p_s["sigma_fallo"]) / E_s_app
            r = caso.correr("plato_nl_fallo", eps_plato=eps_f)
            if r is not None:
                E_pn = caso.fuerza(r[1]) / caso.A_bruta / eps_f
                fila["plato_nl_fallo"] = {"eps": eps_f,
                                          "E_app_MPa": E_pn / 1e6,
                                          "dE_rel": E_pn / E_pl - 1.0}

    fila["t_febio_s"] = caso.tiempos
    DIR_CAMPOS.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(DIR_CAMPOS / f"{analisis}_{est}.npz", **campos)

    print(f"   E_app spinpy {fila['E_app_spinpy_MPa']:.4f}  FEBio "
          f"{fila['E_app_febio_MPa']:.4f} MPa  dE {fila['dE_rel']:.1e}  "
          f"du {fila['du_max_rel']:.1e}  dsig {fila['dsig_max_rel']:.1e}")
    print(f"   p99sup d {fila['d_vm_p99_superficie_rel']:.1e}  Pistoia d "
          f"{fila['d_factor_rel']:.1e}  sesion dE "
          f"{fila.get('d_sesion_E_rel', float('nan')):.1e}")
    return fila


# ---------------------------------------------------------------------------
# Homogeneizacion periodica
# ---------------------------------------------------------------------------

def homogeneizacion(est, BW, spc, doc):
    bw, spr = remuestrear_bw(BW, spc, N_HOMOG)
    print(f"\n== homogeneizacion / {est}: {bw.shape}, BV/TV {bw.mean():.4f}")
    t0 = time.perf_counter()
    Cs, info = homogeneizar(bw, 20e9, NU, vox_size=spr)
    t_s = time.perf_counter() - t0
    if not info["ok"]:
        raise RuntimeError(info["msg"])
    d = DIR_CORR / "homogeneizacion" / est
    d.mkdir(parents=True, exist_ok=True)
    C = np.zeros((6, 6))
    tiempos, ver = {}, None
    for j in range(6):
        cols = []
        for a in AMP_HOMOG:
            e = np.zeros(6)
            e[j] = a
            feb = d / f"homog_{est}_{j}_{a:g}.feb"
            escribir_periodico(bw, spr, feb, e, E_s=20e9, nu_s=NU)
            t, ver = correr_febio(feb)
            tiempos[feb.stem] = t
            s = _ultimo_registro(feb.with_name(feb.stem + "_s.txt"), 6)
            cols.append(s.mean(axis=0) * 1e6 / a)
        C[:, j] = 2 * cols[0] - cols[1]
        print(f"   caso {j}: {sum(list(tiempos.values())[-2:]):.0f} s")
    Csym = 0.5 * (C + C.T)
    ces, cef = constantes_ingenieria(Cs), constantes_ingenieria(Csym)
    fila = {"estructura": est, "analisis": "homogeneizacion", "n": N_HOMOG,
            "BVTV": float(bw.mean()), "version_febio": ver,
            "C_spinpy_Pa": Cs.tolist(), "C_febio_Pa": C.tolist(),
            "dC_max_rel": float(np.abs(C - Cs).max() / np.abs(Cs).max()),
            "asimetria_febio": float(np.abs(C - C.T).max() / np.abs(C).max()),
            "d_log_euclidea": float(distancia_log_euclidea(Cs, Csym)),
            "constantes_spinpy": {k: float(v) for k, v in ces.items()
                                  if np.isscalar(v)},
            "constantes_febio": {k: float(v) for k, v in cef.items()
                                 if np.isscalar(v)},
            "solver_spinpy": info.get("solver"),
            "residuo_spinpy": info.get("residuo_rel"),
            "t_spinpy_s": t_s, "t_febio_s": tiempos}
    fila["d_constantes_rel"] = {k: _rel(fila["constantes_febio"][k], v)
                                for k, v in fila["constantes_spinpy"].items()
                                if k in fila["constantes_febio"] and v}
    reg = _registro_sesion(doc, "elastico", est)
    if reg is not None:
        Cses = np.asarray(reg["C_Pa"], float)
        fila["d_sesion_C_rel"] = float(np.abs(Cs - Cses).max()
                                       / np.abs(Cses).max())
    print(f"   dC {fila['dC_max_rel']:.1e}  d_LE {fila['d_log_euclidea']:.1e}"
          f"  sesion {fila.get('d_sesion_C_rel', float('nan')):.1e}")
    return fila


# ---------------------------------------------------------------------------
# Convergencia en malla
# ---------------------------------------------------------------------------

def convergencia(est, BW, spc, doc):
    filas = []
    ses = doc["resultados"].get("convergencia") or {}
    ses_pts = {p["n"]: p for p in ses.get("puntos", [])} \
        if ses.get("estructura_codigo") == est else {}
    for n in N_CONV:
        bw, spr = remuestrear_bw(BW, spc, n)
        res = ensayo_compresion(bw, spr, E_s=20e9, nu_s=NU, sigma0=1e6,
                                apoyo="deslizante")
        caso = Caso(bw, spr, "deslizante", 20e9,
                    DIR_CORR / "convergencia" / est, f"conv_{est}_{n}")
        u_f, s_f, _ = caso.lineal(1e6)
        E_f = caso.E_app(u_f, 1e6)
        fila = {"estructura": est, "analisis": "convergencia", "n": n,
                "h_mm": float(spr[0]), "BVTV": float(bw.mean()),
                "frac_portante": res["frac_portante"],
                "n_elem": res["n_elem"], "version_febio": caso.version,
                "E_app_spinpy_MPa": res["E_app"] / 1e6,
                "E_app_febio_MPa": E_f / 1e6,
                "dE_rel": _rel(E_f, res["E_app"]),
                "t_febio_s": caso.tiempos}
        if n in ses_pts:
            fila["E_app_sesion_MPa"] = ses_pts[n]["E_app"] / 1e6
            fila["d_sesion_E_rel"] = _rel(res["E_app"], ses_pts[n]["E_app"])
        print(f"   conv {est} n={n}: spinpy {fila['E_app_spinpy_MPa']:.3f} "
              f"FEBio {fila['E_app_febio_MPa']:.3f}  dE {fila['dE_rel']:.1e}")
        filas.append(fila)
    return filas


# ---------------------------------------------------------------------------
# Fallo progresivo: el bucle entero, dos veces
# ---------------------------------------------------------------------------

def fallo(est, BW, spc, doc):
    """`simulacion.fallo_progresivo` de spinpy, y el MISMO algoritmo con
    FEBio resolviendo cada paso. Los dos bucles son independientes: cada uno
    decide que se rompe a partir de su propio campo."""
    print(f"\n== fallo progresivo / {est}")
    t0 = time.perf_counter()
    sp_res = fallo_progresivo(BW, spc, pasos=FALLO["pasos"],
                              n_mec=FALLO["n_mec"], guardar_mascaras=True)
    t_s = time.perf_counter() - t0
    actual, spr = remuestrear_bw(BW, spc, FALLO["n_mec"])
    modulo = np.ones(actual.shape)
    danado = np.zeros(actual.shape, bool)
    filas = []
    for i in range(FALLO["pasos"] + 1):
        caso = Caso(actual, spr, "deslizante", 20e9,
                    DIR_CORR / "fallo" / est, f"fallo_{est}_{i}",
                    modulo=modulo)
        u, s, _ = caso.lineal(1e6)
        E_f = caso.E_app(u, 1e6)
        ee = caso.eps_eff(s)
        campo = np.full(actual.shape, np.nan)
        campo[tuple(caso.ijk.T)] = ee
        sano = actual & ~danado & np.isfinite(campo)
        e = campo[sano]
        umbral = float(np.percentile(e, 100.0 * (1.0 - FRAC_CRITICA)))
        k = EPS_CRITICA / umbral
        F_f = k * 1e6 * caso.A_bruta * 1e-6
        rotos = sano & (k * np.nan_to_num(campo, nan=0.0)
                        >= EPS_CRITICA * (1.0 - 1e-12))
        p = sp_res["pasos"][i]
        # ¿el tejido danado de FEBio en este paso es el mismo que el de spinpy?
        m_s = sp_res["mascaras"][i] if i < len(sp_res["mascaras"]) else None
        dif = (int(np.sum(m_s != (actual & ~danado)))
               if m_s is not None else None)
        fila = {"estructura": est, "analisis": "fallo", "paso": i,
                "n_mec": FALLO["n_mec"], "version_febio": caso.version,
                "E_app_spinpy_MPa": p["E_app"] / 1e6,
                "E_app_febio_MPa": E_f / 1e6,
                "dE_rel": _rel(E_f, p["E_app"]),
                "F_fallo_spinpy_N": p["F_fallo"],
                "F_fallo_febio_N": F_f,
                "dF_rel": _rel(F_f, p["F_fallo"]),
                "rotos_spinpy": p["rotos"], "rotos_febio": int(rotos.sum()),
                "voxeles_intactos_distintos": dif,
                "dano_acumulado": p["dano_acumulado"],
                "E_rel_spinpy": p.get("E_rel"),
                "tras_colapso": p.get("tras_colapso")}
        print(f"   paso {i}: E spinpy {fila['E_app_spinpy_MPa']:.3f} FEBio "
              f"{fila['E_app_febio_MPa']:.3f} (dE {fila['dE_rel']:.1e}); "
              f"rotos {p['rotos']}/{int(rotos.sum())}; intactos distintos "
              f"{dif}")
        filas.append(fila)
        if i < FALLO["pasos"]:
            danado |= rotos
            modulo[rotos] = RIGIDEZ_DANADA
    filas[0]["resumen_spinpy"] = sp_res.get("resumen")
    filas[0]["t_spinpy_s"] = t_s
    return filas


# ---------------------------------------------------------------------------

ANALISIS = ("compresion", "comparado", "homogeneizacion", "convergencia",
            "fallo")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("analisis", nargs="*", default=list(ANALISIS))
    ap.add_argument("--estructuras", nargs="*", default=list(ESTRUCTURAS))
    ap.add_argument("--voi", default=None,
                    help="ruta de VOI_V1.mat si no esta donde la guardo la "
                         "sesion (se comprueba su SHA-256)")
    a = ap.parse_args()
    global RUTA_VOI
    RUTA_VOI = a.voi
    if not FEBIO.exists():
        sys.exit(f"No se encuentra FEBio en {FEBIO}")
    doc = sesion()
    print(f"spinpy {VERSION_SPINPY}; sesion {SESION}")
    estructuras = {e: cargar(e) for e in a.estructuras}
    for an in a.analisis:
        for est, (BW, spc) in estructuras.items():
            t0 = time.perf_counter()
            if an == "compresion":
                filas = [ensayo(est, an, BW, spc, APP, doc)]
            elif an == "comparado":
                filas = [ensayo(est, an, BW, spc, TAPIA, doc)]
            elif an == "homogeneizacion":
                filas = [homogeneizacion(est, BW, spc, doc)]
            elif an == "convergencia":
                filas = convergencia(est, BW, spc, doc)
            elif an == "fallo":
                filas = fallo(est, BW, spc, doc)
            else:
                sys.exit(f"analisis desconocido: {an}")
            for f in filas:
                f["t_total_s"] = time.perf_counter() - t0
                _escribir(an, f)


if __name__ == "__main__":
    main()
