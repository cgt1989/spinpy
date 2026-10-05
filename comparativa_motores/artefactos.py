"""
artefactos.py: Diagnostico de artefactos numericos de los ensayos FEM.

    cd comparativa_motores
    python artefactos.py            # corre todo y escribe resultados/artefactos.json

Sobre el espinodoide de referencia de la comparativa (`casos.ESPINODOIDE`,
5 mm de lado) se resuelve el ensayo de compresion lineal de la app con
NGSolve en varias mallas (hex8 a 32^3, 48^3 y 64^3; TET10 a 32^3 y 48^3) y
condiciones de contorno (traccion uniforme o plato rigido; base deslizante o
empotrada), y se localiza de donde salen los valores altos de von Mises.

Cada elemento recibe unas ETIQUETAS de origen posible (no excluyentes):

  borde_lateral   centroide a <= BANDA mm de una cara lateral del VOI
                  (trabeculas cortadas por el recorte: superficie libre
                  artificial).
  techo, base     centroide a <= BANDA mm del plano cargado o del apoyo.
  ancla           a <= BANDA mm de uno de los dos nodos que fijan el solido
                  rigido en el apoyo deslizante.
  union_singular  (hex8) voxel que toca un vertice o una arista donde el
                  solido se une sin compartir cara: una bisagra que la malla
                  de voxeles crea y que no existe en el hueso.
  arista_entrante (hex8) voxel en una arista concava de la escalera (tres
                  de los cuatro voxeles alrededor de la arista son hueso).
  mala_calidad    (TET10) diedro minimo < 5 grados o radio-arista > 3.
  diminuto        (TET10) volumen < 1e-3 de la mediana.

Para cada etiqueta se mide, en la capa superficial (la que se cita): que
fraccion del volumen ocupa, que fraccion de la cola por encima del p99 cae en
ella (enriquecimiento = cociente de ambas) y cuanto cambia el p99 citado si
sus elementos se excluyen. Ademas:

  * un indicador de error por recuperacion de tensiones (Zienkiewicz-Zhu):
    distancia entre la tension del elemento y la promediada en sus nodos;
  * la persistencia de las zonas calientes al cambiar la resolucion y el tipo
    de elemento (la estructura es la misma: misma semilla);
  * la rigidez aparente por las dos definiciones (media de nodos del techo y
    media por area) y su cambio con la malla y las condiciones de contorno;
  * el error del iterativo CG + pyamg de scikit-fem frente al directo en la
    malla TET10 a 32^3, y donde se localiza.
"""

from __future__ import annotations

import json
import sys
import time
from itertools import product
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(AQUI.parent))
import casos                                                    # noqa: E402
from spinpy import fem, motores                                 # noqa: E402
from spinpy.resistencia import _solo_portante, percentil_ponderado  # noqa: E402

E_MPA, NU = 20000.0, 0.30
LADO = casos.LADO_MM
BANDA = 0.30          # mm: ~1/3 de la longitud de onda del espinodoide (0,83 mm)
SALIDA = AQUI / "resultados" / "artefactos.json"
CAMPOS = AQUI / "resultados" / "artefactos_campos"


# ---------------------------------------------------------------------------
# Etiquetas geometricas de la malla de voxeles
# ---------------------------------------------------------------------------

def _componentes_2x2x2():
    """Componentes 6-conexas del solido en cada uno de los 256 bloques 2x2x2."""
    tabla = np.zeros(256, np.int8)
    for codigo in range(256):
        b = np.array([(codigo >> k) & 1 for k in range(8)], bool)
        tabla[codigo] = ndimage.label(b.reshape(2, 2, 2))[1]
    return tabla


_TABLA = _componentes_2x2x2()


def etiquetas_voxel(B):
    """(union_singular, arista_entrante, n_vertices_singulares, n_aristas_singulares)
    por voxel de `B` (alineadas con np.nonzero(B))."""
    nx, ny, nz = B.shape
    P = np.pad(B, 1).astype(np.uint8)
    # VERTICES: bloque de los 8 voxeles que lo rodean. Mas de una componente
    # 6-conexa = el solido se une por el vertice o por una arista.
    cod = np.zeros((nx + 1, ny + 1, nz + 1), np.int16)
    for k, (a, b, c) in enumerate(product((0, 1), repeat=3)):
        cod |= P[a:a + nx + 1, b:b + ny + 1, c:c + nz + 1].astype(np.int16) << k
    sing = _TABLA[cod] > 1
    union = np.zeros(B.shape, bool)
    for a, b, c in product((0, 1), repeat=3):
        # el voxel (i, j, k) toca los vertices (i+a, j+b, k+c)
        union |= sing[a:a + nx, b:b + ny, c:c + nz]
    union &= B
    # ARISTAS: los cuatro voxeles alrededor de cada arista, por eje.
    entr = np.zeros(B.shape, bool)
    n_ar_sing = 0
    for eje in range(3):
        Q = np.moveaxis(P, eje, 2)              # la arista va a lo largo de z
        mx, my, mz = Q.shape[0] - 2, Q.shape[1] - 2, Q.shape[2] - 2
        c00 = Q[0:mx + 1, 0:my + 1, 1:mz + 1]
        c10 = Q[1:mx + 2, 0:my + 1, 1:mz + 1]
        c01 = Q[0:mx + 1, 1:my + 2, 1:mz + 1]
        c11 = Q[1:mx + 2, 1:my + 2, 1:mz + 1]
        s = c00 + c10 + c01 + c11
        tres = s == 3
        diag = (s == 2) & (c00 == c11)
        n_ar_sing += int(diag.sum())
        marca = np.zeros((mx, my, mz), bool)
        for (a, b) in product((0, 1), repeat=2):
            # voxel (i, j) toca las aristas (i+a, j+b)
            marca |= tres[a:a + mx, b:b + my, :]
        entr |= np.moveaxis(marca, 2, eje)
    entr &= B
    idx = np.nonzero(B)
    return union[idx], entr[idx], int(sing.sum()), n_ar_sing


# ---------------------------------------------------------------------------
# Calidad de los tetraedros (esquinas)
# ---------------------------------------------------------------------------

def calidad_tet(nodos, elems):
    """(radio-arista, diedro minimo en grados) por elemento."""
    P = nodos[elems[:, :4]]
    a, b, c, d = P[:, 0], P[:, 1], P[:, 2], P[:, 3]
    L = np.linalg.norm(np.stack([b - a, c - a, d - a, c - b, d - b, d - c], 1),
                       axis=2)
    u, v, w = b - a, c - a, d - a
    V = np.abs(np.einsum("ij,ij->i", u, np.cross(v, w))) / 6
    uu, vv, ww = (u * u).sum(1), (v * v).sum(1), (w * w).sum(1)
    cr = np.linalg.norm(uu[:, None] * np.cross(v, w)
                        + vv[:, None] * np.cross(w, u)
                        + ww[:, None] * np.cross(u, v), axis=1) / (12 * V)
    n = []
    for i, j, k, o in ((0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 3, 1), (1, 2, 3, 0)):
        nn = np.cross(P[:, j] - P[:, i], P[:, k] - P[:, i])
        s = np.sign(np.einsum("ij,ij->i", nn, P[:, o] - P[:, i]))
        n.append(-nn * s[:, None] / np.linalg.norm(nn, axis=1)[:, None])
    dmin = np.full(len(elems), 180.0)
    for x in range(4):
        for y in range(x + 1, 4):
            ang = 180 - np.degrees(np.arccos(np.clip((n[x] * n[y]).sum(1),
                                                     -1, 1)))
            dmin = np.minimum(dmin, ang)
    return cr / L.min(1), dmin


# ---------------------------------------------------------------------------
# Indicador de error por recuperacion (Zienkiewicz-Zhu)
# ---------------------------------------------------------------------------

def indicador_zz(nodos, elems, sig, vol):
    """|sigma_e - sigma*_e| (norma de Frobenius del tensor) por elemento.

    sigma* es la tension promediada en las ESQUINAS, ponderada por volumen,
    y devuelta al elemento como la media de sus esquinas. En una solucion
    exacta suave ambas coinciden; la diferencia es grande donde el campo da
    saltos entre elementos vecinos (error de discretizacion).
    """
    esq = elems[:, :8] if elems.shape[1] == 8 else elems[:, :4]
    k = esq.shape[1]
    acc = np.zeros((nodos.shape[0], 6))
    pes = np.zeros(nodos.shape[0])
    for j in range(k):
        np.add.at(acc, esq[:, j], sig * vol[:, None])
        np.add.at(pes, esq[:, j], vol)
    s_nod = acc / np.maximum(pes, 1e-300)[:, None]
    s_rec = s_nod[esq].mean(axis=1)
    d = sig - s_rec
    w = np.array([1, 1, 1, 2, 2, 2], float)
    return np.sqrt((d * d * w).sum(1))


# ---------------------------------------------------------------------------
# Un ensayo
# ---------------------------------------------------------------------------

_MALLAS = {}


def malla(tipo, n):
    if (tipo, n) not in _MALLAS:
        BW, sp = casos.espinodoide(n)
        t0 = time.perf_counter()
        m = fem.mallar(BW, sp, tipo)
        m["_t_malla"] = time.perf_counter() - t0
        m["_BW"] = BW
        _MALLAS[(tipo, n)] = m
    return _MALLAS[(tipo, n)]


def ensayo(tipo, n, control="fuerza", apoyo="deslizante", motor="ngsolve",
           solver="auto"):
    m = malla(tipo, n)
    nodos, elems = m["nodos"], np.asarray(m["elems"], np.int64)
    cargas = (1.0,) if control == "fuerza" else (0.001,)
    p = motores.problema_de_malla(m, E=E_MPA, nu=NU, apoyo=apoyo,
                                  control=control, cargas=cargas,
                                  solver=solver)
    t0 = time.perf_counter()
    out = motores.resolver(p, motor)
    t_sol = time.perf_counter() - t0
    u = out["u"]
    nodos = p["nodos"]                              # aristas rectificadas
    sig = fem.tension_elemental(nodos, elems, u, E_MPA, NU)
    vm = fem.von_mises(sig)
    F = float(out["F_reac"][-1])
    s_app = F / m["A_bruta"]                        # MPa, seccion bruta
    H = m["H"]
    uz_area = fem._uz_medio(m, u)
    uz_nod = fem._uz_medio_nodal(m, u)
    vol = np.asarray(m["vol_elem"], float)
    sup = np.asarray(m["superficie"], bool)
    # Coordenadas fisicas comunes: hex8 ocupa [0, L]; la malla suave tiene
    # los planos del cubo en -h/2 y (n - 1/2) h.
    desp = 0.0 if tipo == "hex8" else 0.5 * float(m["spacing"][0])
    cen = nodos[elems[:, :8] if tipo == "hex8" else elems[:, :4]].mean(1) + desp
    anc = nodos[[p["meta"]["ancla_xy"], p["meta"]["ancla_y"]]] + desp
    et = {
        "borde_lateral": np.minimum(cen[:, :2], LADO - cen[:, :2]).min(1) <= BANDA,
        "techo": cen[:, 2] >= LADO - BANDA,
        "base": cen[:, 2] <= BANDA,
    }
    et["ancla"] = (cKDTree(anc).query(cen)[0] <= BANDA) if apoyo == "deslizante" \
        else np.zeros(len(cen), bool)
    extra = {}
    if tipo == "hex8":
        B = _solo_portante(m["_BW"])
        assert int(B.sum()) == elems.shape[0]
        un, en, nvs, nas = etiquetas_voxel(B)
        et["union_singular"], et["arista_entrante"] = un, en
        extra.update(vertices_singulares=nvs, aristas_singulares=nas)
    else:
        ra, dm = calidad_tet(nodos, elems)
        et["mala_calidad"] = (dm < 5.0) | (ra > 3.0)
        et["diminuto"] = vol < 1e-3 * np.median(vol)
        extra.update(perdida_volumen_pct=m["informe"].get("perdida_pct"),
                     descartado_pct=m["informe"].get("descartado_pct"),
                     rectificacion_max_mm=p["meta"]["rectificacion_max_mm"])
    eta = indicador_zz(nodos, elems, sig, vol)
    return {"tipo": tipo, "n": n, "control": control, "apoyo": apoyo,
            "motor": motor, "solver": out["meta"].get("solver"),
            "residuo_rel": out["meta"].get("residuo_rel"),
            "n_gdl": int(3 * nodos.shape[0]) if tipo == "hex8" else
            int(out["meta"].get("n_gdl", 3 * nodos.shape[0])),
            "n_elems": int(elems.shape[0]), "t_malla_s": m["_t_malla"],
            "t_sol_s": t_sol, "F_N": F, "sigma_app_MPa": s_app,
            "E_app_area_MPa": s_app / (abs(uz_area) / H),
            "E_app_nodal_MPa": s_app / (abs(uz_nod) / H),
            "V_mm3": float(vol.sum()), **extra,
            "_u": u, "_vm": vm / s_app, "_eta": eta / s_app, "_sup": sup,
            "_vol": vol, "_cen": cen, "_et": et}


# ---------------------------------------------------------------------------
# Atribucion de la cola
# ---------------------------------------------------------------------------

def p99(v, w):
    return float(percentil_ponderado(v, w, 99))


def atribucion(r):
    vm, w, sup = r["_vm"], r["_vol"], r["_sup"]
    vs, ws = vm[sup], w[sup]
    q = p99(vs, ws)
    cola = vs >= q
    fil = {"p99_sup": q, "max_sup": float(vs.max()),
           "max_global": float(vm.max()), "p50_sup": float(
               percentil_ponderado(vs, ws, 50))}
    todas = np.zeros(sup.sum(), bool)
    clases = {}
    for k, mask in r["_et"].items():
        ms = mask[sup]
        todas |= ms
        fv = float(ws[ms].sum() / ws.sum())
        fc = float(ws[cola & ms].sum() / ws[cola].sum())
        q_sin = p99(vs[~ms], ws[~ms]) if (~ms).any() else np.nan
        clases[k] = {"frac_vol_sup": fv, "frac_cola": fc,
                     "enriquecimiento": fc / fv if fv > 0 else np.nan,
                     "p99_sin": q_sin, "dp99_pct": 100 * (q_sin / q - 1),
                     "max_en_clase": float(vs[ms].max()) if ms.any() else np.nan}
    # Lo que queda: superficie interior, lejos de bordes, sin defectos de malla
    lim = ~todas
    clases["ninguna"] = {"frac_vol_sup": float(ws[lim].sum() / ws.sum()),
                         "frac_cola": float(ws[cola & lim].sum() / ws[cola].sum())}
    clases["ninguna"]["enriquecimiento"] = (clases["ninguna"]["frac_cola"]
                                            / clases["ninguna"]["frac_vol_sup"])
    clases["ninguna"]["p99_solo"] = p99(vs[lim], ws[lim])
    fil["clases"] = clases
    # Indicador ZZ: global y en la cola
    eta, vmr = r["_eta"], vm
    fil["zz_global"] = float(np.sqrt((eta ** 2 * w).sum()
                                     / ((vmr ** 2) * w).sum()))
    fil["zz_rel_mediana_sup"] = float(np.median(eta[sup] / np.maximum(vmr[sup], 1e-12)))
    fil["zz_rel_mediana_cola"] = float(np.median(
        (eta[sup] / np.maximum(vmr[sup], 1e-12))[cola]))
    return fil


def persistencia(ra, rb):
    """Concordancia espacial de los campos de von Mises de dos mallas.

    Se lleva el campo de `rb` a los centroides superficiales de `ra` (vecino
    mas proximo entre los superficiales de rb) y se mide la correlacion de
    rangos y que fraccion de la cola (>= p99) de `ra` tiene, a <= BANDA mm,
    algun elemento de la cola de `rb`.
    """
    from scipy.stats import spearmanr
    sa, sb = ra["_sup"], rb["_sup"]
    ca, cb = ra["_cen"][sa], rb["_cen"][sb]
    va, vb = ra["_vm"][sa], rb["_vm"][sb]
    qa, qb = p99(va, ra["_vol"][sa]), p99(vb, rb["_vol"][sb])
    d, i = cKDTree(cb).query(ca)
    rho = float(spearmanr(va, vb[i]).statistic)
    cola_b = cb[vb >= qb]
    da = cKDTree(cola_b).query(ca[va >= qa])[0]
    return {"spearman": rho, "cola_cerca": float((da <= BANDA).mean()),
            "dist_vecino_mediana_mm": float(np.median(d))}


# ---------------------------------------------------------------------------
# Campana
# ---------------------------------------------------------------------------

CORRIDAS = [
    ("hex8", 32, "fuerza", "deslizante"),
    ("hex8", 48, "fuerza", "deslizante"),
    ("hex8", 64, "fuerza", "deslizante"),
    ("hex8", 48, "plato", "deslizante"),
    ("hex8", 48, "fuerza", "empotrado"),
    ("tet10", 32, "fuerza", "deslizante"),
    ("tet10", 48, "fuerza", "deslizante"),
    ("tet10", 48, "plato", "deslizante"),
]


def clave(r):
    return f"{r['tipo']}_n{r['n']}_{r['control']}_{r['apoyo']}" + (
        "" if r["motor"] == "ngsolve" else f"_{r['motor']}")


def _publico(r):
    return {k: v for k, v in r.items() if not k.startswith("_")}


def main():
    CAMPOS.mkdir(parents=True, exist_ok=True)
    res, todo = {}, {}
    for tipo, n, control, apoyo in CORRIDAS:
        r = ensayo(tipo, n, control, apoyo)
        k = clave(r)
        r.update(atribucion(r))
        todo[k] = r
        res[k] = _publico(r)
        np.savez_compressed(CAMPOS / f"{k}.npz", vm=r["_vm"], eta=r["_eta"],
                            sup=r["_sup"], vol=r["_vol"], cen=r["_cen"],
                            **{f"et_{a}": b for a, b in r["_et"].items()})
        print(f"{k}: {r['n_gdl']} GDL, {r['t_sol_s']:.1f} s, "
              f"E_area {r['E_app_area_MPa']:.2f}, p99 {r['p99_sup']:.3f}",
              flush=True)
    # Error del iterativo de scikit-fem en TET10 (seccion 4.7 del informe)
    try:
        r = ensayo("tet10", 32, motor="skfem", solver="iterativo")
        ref = todo["tet10_n32_fuerza_deslizante"]
        k = clave(r)
        r.update(atribucion(r))
        err = np.abs(r["_vm"] - ref["_vm"]) / np.maximum(ref["_vm"], 1e-12)
        ra, dm = calidad_tet(malla("tet10", 32)["nodos"],
                             np.asarray(malla("tet10", 32)["elems"]))
        malo = ref["_et"]["mala_calidad"] | ref["_et"]["diminuto"]
        top = np.argsort(err)[-int(0.001 * err.size):]
        r["error_iterativo"] = {
            "du_rel": float(np.abs(r["_u"] - ref["_u"]).max()
                            / np.abs(ref["_u"]).max()),
            "vm_err_max": float(err.max()), "vm_err_p99": float(
                np.percentile(err, 99)),
            "vm_err_mediana": float(np.median(err)),
            "frac_malo_en_top01": float(malo[top].mean()),
            "frac_malo_total": float(malo.mean()),
            "dp99_pct": 100 * (r["p99_sup"] / ref["p99_sup"] - 1),
            "dE_area_pct": 100 * (r["E_app_area_MPa"]
                                  / ref["E_app_area_MPa"] - 1)}
        res[k] = _publico(r)
        print(k, r["error_iterativo"], flush=True)
    except Exception as e:                                # noqa: BLE001
        res["tet10_n32_skfem"] = {"error": f"{type(e).__name__}: {e}"}
    # Persistencia entre mallas
    pares = [("hex8_n32_fuerza_deslizante", "hex8_n64_fuerza_deslizante"),
             ("hex8_n48_fuerza_deslizante", "hex8_n64_fuerza_deslizante"),
             ("tet10_n32_fuerza_deslizante", "tet10_n48_fuerza_deslizante"),
             ("hex8_n48_fuerza_deslizante", "tet10_n48_fuerza_deslizante"),
             ("hex8_n64_fuerza_deslizante", "tet10_n48_fuerza_deslizante"),
             ("hex8_n48_fuerza_deslizante", "hex8_n48_plato_deslizante"),
             ("tet10_n48_fuerza_deslizante", "tet10_n48_plato_deslizante"),
             ("hex8_n48_fuerza_deslizante", "hex8_n48_fuerza_empotrado")]
    res["persistencia"] = {f"{a}|{b}": persistencia(todo[a], todo[b])
                           for a, b in pares}
    for k, v in res["persistencia"].items():
        print(k, v, flush=True)
    SALIDA.write_text(json.dumps(res, indent=1, default=float))
    print("->", SALIDA)


if __name__ == "__main__":
    main()
