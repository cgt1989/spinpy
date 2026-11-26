"""Comparacion de familias generadoras frente a un VOI, con la maquinaria de spinpy.

Uso: python exp.py VOI FAMILIA
FAMILIA: spinodoide | dual | open_foam | levelset | levelset_tpl

Protocolo comun (el de spinpy):
  * candidato a res = min(shape, 96) sobre el cubo fisico del VOI (lc/res mm/voxel)
  * morfometria = spinpy.morphometry.morfometria, misma para VOI y candidatos
  * error = spinpy.error.error_morfometrico (con y sin el termino DA2)
  * busqueda escalonada A (densidad x tamano) + B (23 presets de estiramiento) + C (refinado)
  * el ganador se reevalua con 10 semillas NUEVAS (1001..1010), distintas de la de busqueda
"""
import sys, json, time, math, os
sys.path.insert(0, '../../..')
import numpy as np

from spinpy.io import leer_voi
from spinpy.morphometry import morfometria, metricas_forma, conectividad, indice_smi
from spinpy.error import error_morfometrico
from spinpy.fit import (ajustar_spinodoide, longitud_caracteristica,
                        resolucion_comparacion, fraccion_poro_conexa,
                        SEMILLA_POR_DEFECTO, RHO_MAX, RHO_MIN)
from spinpy.fit_dual import ajustar_dual_lattice, PRESETS_COMPLETO, PESO_DA2
from spinpy.grf import generar_mascara
from spinpy.dual_lattice import generar_dual_lattice

OUT = os.environ.get('EXP_OUT', './out')
os.makedirs(OUT, exist_ok=True)
SEMILLAS_NUEVAS = list(range(1001, 1011))
ESCALARES = ['BVTV', 'TbTh', 'TbSp', 'TbN', 'BSBV', 'PoTot', 'DA', 'DA2']


def errores(m_voi, m):
    e0, _ = error_morfometrico(m_voi, m)
    e1, _ = error_morfometrico(m_voi, m, peso_da2=PESO_DA2)
    return float(e0), float(e1)


def escalares(m):
    return {k: float(m[k]) for k in ESCALARES if k in m and np.ndim(m[k]) == 0}


# ---------------------------------------------------------------------------
# generadores de morphoGen
# ---------------------------------------------------------------------------
def _dominio():
    import morphogen as mg
    return mg, mg.Domain.cube(1.0, periodic=True)


def gen_open_foam(res, d, N, e, seed):
    mg, dom = _dominio()
    model = mg.models.OpenFoam(n_cells=int(round(N)), volume_fraction=float(d),
                               anisotropy={'stretch': [float(x) for x in e]},
                               relaxation=10, node_ratio=1.3)  # plantilla trabecular_bone
    micro = model.generate(dom, seed=int(seed))
    return np.asarray(micro.voxelize(shape=res).data).astype(bool)


def gen_levelset(res, d, N, e, seed, plantilla=False):
    mg, dom = _dominio()
    pf = 0.6
    rg = (pf / (4.0 / 3.0 * math.pi * float(N))) ** (1.0 / 3.0)   # radio geometrico de la semilla
    ex, ey, ez = (float(x) for x in e)
    b, c = ey / ex, ez / ex
    a = rg / (b * c) ** (1.0 / 3.0)
    kw = dict(curvature=-0.7, c1=50.0, c2=1.0) if plantilla else {}
    model = mg.models.LevelSetFoam(seed_radius=a, aspect=(b, c), packing='dense',
                                   packing_fraction=pf, volume_fraction=float(d), **kw)
    micro = model.generate(dom, seed=int(seed))
    return np.asarray(micro.voxelize(shape=res).data).astype(bool)


GENERADORES = {
    'open_foam': gen_open_foam,
    'levelset': gen_levelset,
    'levelset_tpl': lambda *a: gen_levelset(*a, plantilla=True),
}


def buscar(familia, m_voi, res, sp_cand, bvtv, tiempo_max):
    """Busqueda escalonada A+B+C con el mismo presupuesto que ajustar_dual_lattice."""
    gen = GENERADORES[familia]
    dens = np.linspace(max(RHO_MIN, bvtv * 0.7), min(RHO_MAX, bvtv * 1.3), 5)
    tamanos = np.unique(np.round(np.geomspace(50, 800, 6)))
    traza, t0 = [], time.time()

    def evaluar(d, N, e, etapa):
        try:
            BW = gen(res, d, N, e, SEMILLA_POR_DEFECTO)
        except ValueError as exc:   # combinacion inviable para el modelo (p. ej. no cabe el empaquetado)
            traza.append(dict(etapa=etapa, d=float(d), N=float(N), e=[float(x) for x in e],
                              err=1e6, err_da2=1e6, inviable=str(exc)[:80]))
            return 1e6
        m = morfometria(BW, sp_cand)
        e0, e1 = errores(m_voi, m)
        traza.append(dict(etapa=etapa, d=float(d), N=float(N), e=[float(x) for x in e],
                          err=e0, err_da2=e1))
        return e1  # el objetivo de busqueda incluye DA2, como el dual-lattice de spinpy

    best = dict(err=np.inf, d=dens[0], N=tamanos[0], e=(1.0, 1.0, 1.0))
    for d in dens:
        for N in tamanos:
            er = evaluar(d, N, best['e'], 'A')
            if er < best['err']:
                best.update(err=er, d=d, N=N)
    eb = best['e']
    for e in PRESETS_COMPLETO:
        er = evaluar(best['d'], best['N'], e, 'B')
        if er < best['err']:
            best.update(err=er, e=e)
    d_ref = np.unique(np.clip(best['d'] * np.array([0.92, 1.0, 1.08]), RHO_MIN, RHO_MAX))
    N_ref = np.unique(np.round(best['N'] * np.array([0.8, 1.0, 1.25])))
    for d in d_ref:
        for N in N_ref:
            er = evaluar(d, N, best['e'], 'C')
            if er < best['err']:
                best.update(err=er, d=d, N=N)
    return best, traza, time.time() - t0


def main(voi, familia):
    BW_voi, spacing = leer_voi(f'../../../datos/porcino_koria2020/VOI_{voi}.mat')
    spacing = np.asarray(spacing, float).ravel()
    m_voi = morfometria(BW_voi, spacing)
    lc = longitud_caracteristica(BW_voi.shape, spacing)
    res = resolucion_comparacion(BW_voi.shape)
    sp_cand = np.full(3, lc / res)
    bvtv = float(m_voi['BVTV'])
    res_json = dict(voi=voi, familia=familia, res=res, lc=lc, voi_escalares=escalares(m_voi))
    t0 = time.time()

    if familia == 'spinodoide':
        r = ajustar_spinodoide(BW_voi, spacing, modo='completo', m_voi=m_voi)
        p = r['parametros']
        th = np.clip(np.round(p['thetas']), 0, 90)
        R = np.asarray(p['R'])
        gen = lambda s: generar_mascara(resolution=res, wave_number=p['wave_number_rad'],
                                        num_waves=p['num_waves'], thetas=th, rho=p['densidad'],
                                        R=R, esquema=p['esquema'], seed=int(s), precision='f32')[0]
        res_json['parametros'] = {k: p[k] for k in ('densidad', 'wave_number_pi', 'num_waves', 'thetas', 'esquema')}
        res_json['n_eval'] = r['n_evaluaciones']
        res_json['error_busqueda'] = r['error_busqueda']
        res_json['tiempo_busqueda_s'] = r['tiempo_s']
    elif familia == 'dual':
        r = ajustar_dual_lattice(BW_voi, spacing, modo='completo', m_voi=m_voi)
        p = r['parametros']
        R = np.asarray(p['R'])
        gen = lambda s: generar_dual_lattice(res, p['celdas'], p['densidad'],
                                             estiramiento=p['estiramiento'],
                                             irregularidad=p['irregularidad'], R=R, seed=int(s))[0]
        res_json['parametros'] = {k: p[k] for k in ('densidad', 'celdas', 'estiramiento', 'irregularidad')}
        res_json['n_eval'] = r['n_evaluaciones']
        res_json['error_busqueda'] = r['error_busqueda']
        res_json['tiempo_busqueda_s'] = r['tiempo_s']
    else:
        best, traza, dt = buscar(familia, m_voi, res, sp_cand, bvtv, None)
        g = GENERADORES[familia]
        gen = lambda s: g(res, best['d'], best['N'], best['e'], s)
        res_json['parametros'] = dict(densidad=float(best['d']), N=float(best['N']),
                                      estiramiento=[float(x) for x in best['e']])
        res_json['n_eval'] = len(traza)
        res_json['error_busqueda'] = float(best['err'])
        res_json['tiempo_busqueda_s'] = dt
        res_json['traza'] = traza
        res_json['n_inviables'] = sum(1 for t in traza if 'inviable' in t)

    # ---- reevaluacion del ganador con semillas nuevas ----
    reps, masks = [], []
    for s in SEMILLAS_NUEVAS:
        BW = gen(s)
        m = morfometria(BW, sp_cand)
        e0, e1 = errores(m_voi, m)
        reps.append(dict(semilla=s, err=e0, err_da2=e1, **escalares(m)))
        masks.append(BW)
    np.savez_compressed(f'{OUT}/masks_{voi}_{familia}.npz', BW=np.stack(masks))
    res_json['replicas'] = reps
    e0 = np.array([r['err'] for r in reps]); e1 = np.array([r['err_da2'] for r in reps])
    res_json['err_media'] = float(e0.mean()); res_json['err_sd'] = float(e0.std(ddof=1))
    res_json['err_da2_media'] = float(e1.mean()); res_json['err_da2_sd'] = float(e1.std(ddof=1))
    res_json['tiempo_total_s'] = time.time() - t0
    json.dump(res_json, open(f'{OUT}/res_{voi}_{familia}.json', 'w'), indent=1)
    print(voi, familia, 'err %.4f +- %.4f | con DA2 %.4f +- %.4f | n_eval %d | %.0fs' % (
        res_json['err_media'], res_json['err_sd'], res_json['err_da2_media'],
        res_json['err_da2_sd'], res_json['n_eval'], res_json['tiempo_total_s']))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
