"""
estudio.py: Las correcciones de borde (V2.1.0 y V2.1.1) sobre VOIs reales.

    cd comparativa_motores/vois_reales
    python estudio.py CARPETA_DE_VOIS [validacion] [practico] [morfometria] [campos]

Los VOIs no viajan con el repositorio; su huella SHA-256 queda en el JSON.
Dos estudios por VOI:

VALIDACION (configuracion embebida en hueso real). El VOI completo,
remuestreado a 96^3 voxeles, hace de BLOQUE; su cubo central de 64^3 (2/3
del lado, margen de 16 voxeles) hace de "VOI de ensayo". El cubo central se
ensaya como lo haria la app (`fem.ensayo` con 'lineal' y 'lineal_plato') y
se compara con el bloque ensayado con plato y medido solo en la region del
cubo central: alli no hay techo cargado ni caras cortadas. Es el mismo
diseno que se valido en el espinodoide (comparativa_motores/correcciones/),
ahora con hueso de micro-CT. Sensibilidad: el bloque tambien con traccion.

PRACTICO. El VOI completo por `fem.analizar` (el camino de la GUI), con
ladrillos a 40^3 (la resolucion por omision) y malla suave a 40^3, y
NGSolve, en procesos hijo como la GUI: lo que el usuario ve antes y despues
de las correcciones. La malla suave a 48^3 (su omision) no cabe en 15 GB
con estos VOIs: `fem.tamano_previsto` da 7,8 a 10 GB para el directo y el
primer intento (C1, BDDC) agoto la memoria.

CAMPOS. Tension vertical y de von Mises del bloque (con plato) y del cubo
aislado (traccion y plato) en la rejilla del cubo de ensayo: un corte
vertical y dos perfiles (por distancia a la cara lateral y por altura).

Escribe resultados/validacion.json, practico.json, morfometria.json,
campos.json y campos.npz (cortes).
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent.parent))
from spinpy import fem, motores                                 # noqa: E402
from spinpy.elastic import remuestrear_bw                       # noqa: E402
from spinpy.io import leer_voi                                  # noqa: E402

RES = AQUI / "resultados"
E_PA, NU = 20e9, 0.30
N_BLOQUE, N_INTERIOR = 96, 64
#: Margen alternativo del nucleo (mm) para VOIs pequenos, EXPLORATORIO: no es
#: el de la app (0,625 mm) ni esta validado.
MARGEN_EXPLORATORIO = 0.25
#: Resolucion del estudio practico (ladrillos y malla suave).
N_PRACTICO = 40

VOIS = [("VOI_C1.mat", "C1", "porcino"),
        ("VOI_C2.mat", "C2", "porcino"),
        ("VOI_proximal_cubico.vtk", "Prox. cúbico", "equino"),
        ("VOI_proximal_PCAaligned.vtk", "Prox. PCA", "equino"),
        ("VOI_medio_PCAaligned.vtk", "Medio PCA", "equino")]


def huella(ruta):
    return hashlib.sha256(Path(ruta).read_bytes()).hexdigest()


def guardar(nombre, datos):
    RES.mkdir(exist_ok=True)
    (RES / f"{nombre}.json").write_text(json.dumps(datos, indent=1,
                                                   default=float))


def cargar(nombre):
    f = RES / f"{nombre}.json"
    return json.loads(f.read_text()) if f.exists() else {}


def limpio(d):
    return {k: v for k, v in d.items() if not str(k).startswith("_")}


# ---------------------------------------------------------------------------
# Validacion embebida
# ---------------------------------------------------------------------------

def resolver_bloque(malla, control):
    """Bloque con plato o traccion; devuelve u (mm) y sigma (Pa)."""
    carga = (1e-3,) if control == "plato" else (1.0,)
    p = motores.problema_de_malla(malla, E=E_PA / 1e6, nu=NU, control=control,
                                  cargas=carga)
    t0 = time.perf_counter()
    out = motores.resolver(p, "ngsolve")
    t = time.perf_counter() - t0
    s = fem.tension_elemental(p["nodos"], malla["elems"], out["u"], E_PA, NU)
    return out["u"], s, t, out["meta"].get("solver"), int(3 * p["nodos"].shape[0])


def region_rel(nuc):
    """E (MPa) y p99 relativo a la tension media de la region."""
    if not nuc.get("ok"):
        return {"ok": False}
    return {"ok": True, "E": nuc["E_app"] / 1e6,
            "p99": nuc["vm_p99_superficie"] / nuc["sigma_Pa"],
            "margen_mm": nuc["margen_mm"], "frac_volumen": nuc["frac_volumen"]}


def validacion(carpeta):
    out = cargar("validacion")
    for archivo, nombre, especie in VOIS:
        if nombre in out:
            continue
        ruta = Path(carpeta) / archivo
        BW0, sp0 = leer_voi(ruta)
        L = float(BW0.shape[0] * sp0[0])
        B, sp = remuestrear_bw(BW0, sp0, N_BLOQUE)
        h = float(sp[0])
        m = (N_BLOQUE - N_INTERIOR) // 2
        Bi = B[m:m + N_INTERIOR, m:m + N_INTERIOR, m:m + N_INTERIOR]
        d = {"especie": especie, "archivo": archivo, "sha256": huella(ruta),
             "lado_voi_mm": L, "h_mm": h, "lado_interior_mm": N_INTERIOR * h,
             "margen_embebido_mm": m * h, "BVTV_bloque": float(B.mean()),
             "BVTV_interior": float(Bi.mean())}
        # Referencia: el bloque, medido en la region del cubo interior
        mb = fem.mallar(B, sp, "hex8")
        ref = {}
        for control in ("plato", "fuerza"):
            u, s, t, solver, gdl = resolver_bloque(mb, control)
            ref[control] = {
                "t_s": t, "solver": solver, "gdl": gdl,
                "interior": region_rel(fem.magnitudes_nucleo(
                    mb, u, s, margen=m * h)),
                "nucleo": region_rel(fem.magnitudes_nucleo(
                    mb, u, s, margen=m * h + fem.NUCLEO_MARGEN_MM)),
                "nucleo_exp": region_rel(fem.magnitudes_nucleo(
                    mb, u, s, margen=m * h + MARGEN_EXPLORATORIO))}
            print(nombre, "bloque", control, f"{gdl} GDL {t:.0f}s",
                  ref[control]["interior"], flush=True)
        d["referencia"] = ref
        # Metodos: el cubo interior ensayado como un VOI, por la app
        mi = fem.mallar(Bi, sp, "hex8")
        reg = fem.ensayo(mi, fem.protocolo("app", E_s=E_PA),
                         analisis=["lineal", "lineal_plato"],
                         motor="ngsolve", aislado=False)
        sref = reg["sigma_ref_Pa"]
        lin, lp, cor = reg["lineal"], reg["lineal_plato"], reg["corregido"]
        met = {"B0": {"E": lin["E_app"] / 1e6,
                      "p99": lin["pistoia"]["vm_p99_superficie"] / sref},
               "F1": {"E": lp["E_app"] / 1e6,
                      "p99": lp["pistoia"]["vm_p99_superficie"] / sref}}
        for k, n in (("F2", lin["nucleo"]), ("F12", lp["nucleo"])):
            met[k] = ({"E": n["E_app"] / 1e6,
                       "p99": n["vm_p99_superficie"] / sref}
                      if n.get("ok") else {"no_disponible": n.get("msg")})
        # Nucleo exploratorio (margen menor), con los campos del ensayo
        for k, campos in (("F2_exp", reg["_campos_lineal"]),
                          ("F12_exp", reg["_campos_lineal_plato"])):
            n = fem.magnitudes_nucleo(mi, campos["u"], campos["sigma"],
                                      sigma_ref=sref,
                                      margen=MARGEN_EXPLORATORIO)
            met[k] = ({"E": n["E_app"] / 1e6,
                       "p99": n["vm_p99_superficie"] / sref}
                      if n.get("ok") else {"no_disponible": n.get("msg")})
        d["metodos"] = met
        d["corregido"] = {k: (v / 1e6 if k == "E_app" else
                              v / sref if k == "vm_p99_superficie" else v)
                          for k, v in cor.items()}
        d["gdl_interior"] = int(3 * mi["nodos"].shape[0])
        d["sigma_fallo_MPa"] = {
            "traccion": lin["pistoia"]["sigma_fallo"] / 1e6,
            "plato": lp["pistoia"]["sigma_fallo"] / 1e6}
        out[nombre] = d
        guardar("validacion", out)
        print(nombre, "metodos", {k: v.get("E") for k, v in met.items()},
              flush=True)
    return out


# ---------------------------------------------------------------------------
# Uso practico: el VOI completo por la app
# ---------------------------------------------------------------------------

def practico(carpeta):
    out = cargar("practico")
    prot = fem.protocolo("app", E_s=E_PA)
    for archivo, nombre, especie in VOIS:
        BW, sp = leer_voi(Path(carpeta) / archivo)
        for malla in ("hex8", "tet10"):
            k = f"{nombre}|{malla}"
            if k in out:
                continue
            t0 = time.perf_counter()
            try:
                r = fem.analizar(BW, sp, prot, malla=malla,
                                 analisis=["lineal", "lineal_plato"],
                                 motores_fem=["ngsolve"], aislado=True,
                                 comparar_app=(malla == "hex8"), n=N_PRACTICO,
                                 hilos=4)[0]
            except Exception as e:                      # noqa: BLE001
                out[k] = {"error": f"{type(e).__name__}: {e}"}
                guardar("practico", out)
                continue
            t = time.perf_counter() - t0
            lin, lp = r.get("lineal") or {}, r.get("lineal_plato") or {}
            cor = r.get("corregido") or {}
            inf = r.get("informe_malla") or {}
            out[k] = {
                "especie": especie, "ok": r["ok"], "fallos": r["fallos"],
                "n": r.get("n"), "t_s": t, "gdl": inf.get("n_gdl"),
                "perdida_volumen_pct": inf.get("perdida_pct"),
                "E_app": lin.get("E_app", np.nan) / 1e6,
                "E_corregido": cor.get("E_app", np.nan) / 1e6,
                "metodo_E": cor.get("metodo_E"),
                "E_nucleo_traccion": (lin.get("nucleo") or {}).get(
                    "E_app", np.nan) / 1e6,
                "E_plato": lp.get("E_app", np.nan) / 1e6,
                "p99": (lin.get("pistoia") or {}).get(
                    "vm_p99_superficie", np.nan) / 1e6,
                "p99_corregido": cor.get("vm_p99_superficie", np.nan) / 1e6,
                "sigma_fallo": (lin.get("pistoia") or {}).get(
                    "sigma_fallo", np.nan) / 1e6,
                "sigma_fallo_plato": (lp.get("pistoia") or {}).get(
                    "sigma_fallo", np.nan) / 1e6,
                "E_app_app": ((r.get("app") or {}).get("E_app", np.nan)
                              / 1e6)}
            guardar("practico", out)
            print(k, {x: out[k][x] for x in ("ok", "E_app", "E_corregido",
                                             "metodo_E", "t_s")}, flush=True)
    return out


# ---------------------------------------------------------------------------
# Campos: donde actuan los artefactos y las correcciones
# ---------------------------------------------------------------------------

def _rejilla(malla, sp, n, desplaz, campo):
    """Campo por elemento hex8 -> rejilla n^3 del cubo de ensayo (NaN = poro).

    `desplaz`: voxeles entre el origen de la malla y el del cubo (16 en el
    bloque, 0 en el cubo aislado)."""
    c = malla["nodos"][malla["elems"]].mean(axis=1)
    ijk = np.floor(c / np.asarray(sp)).astype(int) - desplaz
    dentro = np.all((ijk >= 0) & (ijk < n), axis=1)
    G = np.full((n, n, n), np.nan)
    i, j, k = ijk[dentro].T
    G[i, j, k] = campo[dentro]
    return G


def _perfiles(Szz, Svm, h):
    """Perfiles normalizados por la tension media del cubo (poros incluidos).

    lateral: tension vertical media en capas a distancia d de la cara
    lateral mas proxima (poros cuentan 0); vertical: von Mises medio en el
    hueso de cada capa horizontal."""
    n = Szz.shape[0]
    S0 = np.nan_to_num(Szz)
    sm = S0.mean()
    c = np.arange(n) + 0.5
    dx = np.minimum(c, n - c)
    d = np.minimum(dx[:, None, None], dx[None, :, None]) * np.ones((1, 1, n))
    bordes = np.arange(0, n // 2 + 1, 2)
    lat = [float(S0[(d >= a) & (d < b)].mean() / sm)
           for a, b in zip(bordes[:-1], bordes[1:])]
    vert = [float(np.nanmean(Svm[:, :, k]) / sm) for k in range(n)]
    return {"d_mm": list((bordes[:-1] + 1) * h), "lateral": lat,
            "z_mm": list(c * h), "vertical": vert, "sigma_media": float(sm)}


def campos(carpeta):
    """Bloque con plato y cubo aislado con traccion y con plato: tension
    vertical y de von Mises en la rejilla del cubo de ensayo."""
    out = cargar("campos")
    f_npz = RES / "campos.npz"
    cortes = dict(np.load(f_npz)) if f_npz.exists() else {}
    for archivo, nombre, especie in VOIS:
        if nombre in out:
            continue
        BW0, sp0 = leer_voi(Path(carpeta) / archivo)
        B, sp = remuestrear_bw(BW0, sp0, N_BLOQUE)
        h = float(sp[0])
        m = (N_BLOQUE - N_INTERIOR) // 2
        Bi = B[m:m + N_INTERIOR, m:m + N_INTERIOR, m:m + N_INTERIOR]
        d = {"especie": especie, "h_mm": h, "voxeles_hueso_cubo": int(Bi.sum())}
        casos = (("referencia", B, m, "plato"), ("traccion", Bi, 0, "fuerza"),
                 ("plato", Bi, 0, "plato"))
        for caso, M, desplaz, control in casos:
            malla = fem.mallar(M, sp, "hex8")
            u, s, t, _, gdl = resolver_bloque(malla, control)
            Szz = _rejilla(malla, sp, N_INTERIOR, desplaz, -s[:, 2])
            Svm = _rejilla(malla, sp, N_INTERIOR, desplaz, fem.von_mises(s))
            p = _perfiles(Szz, Svm, h)
            p["frac_portante"] = float(np.isfinite(Szz).sum() / Bi.sum())
            p["t_s"], p["gdl"] = t, gdl
            d[caso] = p
            sm = p["sigma_media"]
            clave = nombre.replace(" ", "_").replace(".", "")
            cortes[f"{clave}|{caso}|zz"] = (Szz[:, N_INTERIOR // 2, :]
                                            / sm).astype(np.float32)
            cortes[f"{clave}|{caso}|vm"] = (Svm[:, N_INTERIOR // 2, :]
                                            / sm).astype(np.float32)
            print(nombre, caso, f"{gdl} GDL {t:.0f}s portante "
                  f"{p['frac_portante']:.3f}", flush=True)
        out[nombre] = d
        guardar("campos", out)
        np.savez_compressed(f_npz, **cortes)
    return out


def morfometria_vois(carpeta):
    from spinpy.morphometry import morfometria
    out = cargar("morfometria")
    for archivo, nombre, especie in VOIS:
        if nombre in out and "perfil_z" in out[nombre]:
            continue
        BW, sp = leer_voi(Path(carpeta) / archivo)
        if nombre in out:
            # Perfil de BV/TV por capa horizontal (eje de carga), anadido
            # despues: no hace falta repetir la morfometria.
            out[nombre]["perfil_z"] = BW.mean(axis=(0, 1)).tolist()
            guardar("morfometria", out)
            continue
        m = morfometria(BW, sp, do_mil=True)
        out[nombre] = {"especie": especie, "forma": list(BW.shape),
                       "h_mm": float(sp[0]),
                       **{k: float(m[k]) for k in ("BVTV", "TbTh", "TbSp",
                                                   "TbN", "DA") if k in m},
                       "perfil_z": BW.mean(axis=(0, 1)).tolist()}
        guardar("morfometria", out)
        print(nombre, out[nombre], flush=True)
    return out


if __name__ == "__main__":
    carpeta = sys.argv[1]
    etapas = sys.argv[2:] or ["morfometria", "validacion", "practico"]
    if "morfometria" in etapas:
        morfometria_vois(carpeta)
    if "validacion" in etapas:
        validacion(carpeta)
    if "practico" in etapas:
        practico(carpeta)
    if "campos" in etapas:
        campos(carpeta)
