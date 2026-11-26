"""Experimento 4: level_set_foam con sus PARAMETROS DE FORMA ajustados (curvatura K_c y engrosamiento de nodos c1, c2),
densidad fija al BV/TV medido, objetivo ampliado de objetivo.py.

Uso: python exp4.py VOI          (la familia se llama 'levelset_shape')
Busqueda (~60 evaluaciones, el presupuesto de los ajustes de spinpy):
  A  tamano N (6 valores)          B  23 presets de estiramiento        C  refinado de N (3)
  D  rejilla de forma: K_c x nodos (5 x 4 = 20)                          E  reajuste de N (4)
  F  afinado de K_c (2) y de nodos (2)
"""
import os, sys, json, time, math
sys.path.insert(0, '.')
os.environ.setdefault('EXP_OUT', './out4')
import numpy as np
import exp3 as X                      # parches en memoria, objetivo ampliado, utilidades
from exp3 import OBJ, PES, errores, escalares, SEMILLAS_NUEVAS, PRESETS_COMPLETO
from spinpy.io import leer_voi
from spinpy.morphometry import morfometria
from spinpy.fit import longitud_caracteristica, resolucion_comparacion

OUT = os.environ['EXP_OUT']
os.makedirs(OUT, exist_ok=True)
FAM = 'levelset_shape'
KC = [-1.0, -0.7, -0.35, 0.0, 0.3]
NODOS = [(0.0, 0.0), (10.0, 0.5), (50.0, 1.0), (150.0, 2.0)]


def gen_ls(res, d, N, e, seed, kc, c1, c2):
    import morphogen as mg
    dom = mg.Domain.cube(1.0, periodic=True)
    pf = 0.6
    rg = (pf / (4.0 / 3.0 * math.pi * float(N))) ** (1.0 / 3.0)
    ex, ey, ez = (float(x) for x in e)
    b, c = ey / ex, ez / ex
    a = rg / (b * c) ** (1.0 / 3.0)
    model = mg.models.LevelSetFoam(seed_radius=a, aspect=(b, c), packing='dense', packing_fraction=pf,
                                   volume_fraction=float(d), curvature=float(kc), c1=float(c1), c2=float(c2))
    return np.asarray(model.generate(dom, seed=int(seed)).voxelize(shape=res).data).astype(bool)


def buscar(m_voi_ext, res, sp_cand, bvtv):
    traza, t0 = [], time.time()
    best = dict(err=np.inf, N=None, e=(1.0, 1.0, 1.0), kc=0.0, c1=0.0, c2=0.0)

    def ev(N, e, kc, c1, c2, etapa):
        par = dict(etapa=etapa, N=float(N), e=[float(x) for x in e], kc=float(kc), c1=float(c1), c2=float(c2))
        try:
            BW = gen_ls(res, bvtv, N, e, SEMILLA, kc, c1, c2)
        except (ValueError, RuntimeError, ArithmeticError) as exc:
            traza.append(dict(par, err_ext=1e6, inviable=str(exc)[:80])); return
        m = OBJ.medir_ext(BW, sp_cand)
        e_b, e_x = errores(m_voi_ext, m)
        traza.append(dict(par, err_base=e_b, err_ext=e_x))
        if e_x < best['err']:
            best.update(err=e_x, N=N, e=tuple(e), kc=kc, c1=c1, c2=c2)

    for N in np.unique(np.round(np.geomspace(50, 800, 6))):                           # A
        ev(N, best['e'], 0.0, 0.0, 0.0, 'A')
    for e in PRESETS_COMPLETO:                                                          # B
        ev(best['N'], e, 0.0, 0.0, 0.0, 'B')
    for f in (0.8, 1.0, 1.25):                                                          # C
        ev(round(best['N'] * f), best['e'], best['kc'], best['c1'], best['c2'], 'C')
    base = (best['N'], best['e'])
    for kc in KC:                                                                       # D
        for c1, c2 in NODOS:
            ev(base[0], base[1], kc, c1, c2, 'D')
    for f in (0.75, 0.9, 1.1, 1.3):                                                     # E
        ev(round(best['N'] * f), best['e'], best['kc'], best['c1'], best['c2'], 'E')
    kc0, c10, c20 = best['kc'], best['c1'], best['c2']
    for dk in (-0.15, 0.15):                                                            # F
        if kc0 + dk < 0.95:
            ev(best['N'], best['e'], kc0 + dk, c10, c20, 'F')
    for fn in (0.6, 1.6):
        if c10 > 0:
            ev(best['N'], best['e'], best['kc'], c10 * fn, c20, 'F')
        else:
            ev(best['N'], best['e'], best['kc'], 5.0 * fn, 0.5, 'F')
    return best, traza, time.time() - t0


SEMILLA = 20260720


def main(voi):
    BW_voi, spacing = leer_voi(f'../../../datos/porcino_koria2020/VOI_{voi}.mat')
    spacing = np.asarray(spacing, float).ravel()
    m_voi = morfometria(BW_voi, spacing)
    m_voi_ext = dict(m_voi)
    m_voi_ext.update(OBJ.metricas_forma_ext(BW_voi, spacing, tbth=True, poro=True, curvatura=True, conectividad=True))
    lc = longitud_caracteristica(BW_voi.shape, spacing)
    res = resolucion_comparacion(BW_voi.shape)
    sp_cand = np.full(3, lc / res)
    bvtv = float(m_voi['BVTV'])
    best, traza, dt = buscar(m_voi_ext, res, sp_cand, bvtv)
    gen = lambda s: gen_ls(res, bvtv, best['N'], best['e'], s, best['kc'], best['c1'], best['c2'])
    reps, masks = [], []
    for s in SEMILLAS_NUEVAS:
        BW = gen(s); m = OBJ.medir_ext(BW, sp_cand); e_b, e_x = errores(m_voi_ext, m)
        reps.append(dict(semilla=s, err_base=e_b, err_ext=e_x, **escalares(m))); masks.append(BW)
    np.savez_compressed(f'{OUT}/masks_{voi}_{FAM}.npz', BW=np.stack(masks))
    r = dict(voi=voi, familia=FAM, res=res, lc=lc, voi_escalares=escalares(m_voi_ext), pesos=PES,
             parametros=dict(densidad=bvtv, N=float(best['N']), estiramiento=[float(x) for x in best['e']],
                             curvatura=float(best['kc']), c1=float(best['c1']), c2=float(best['c2'])),
             n_eval=len(traza), n_inviables=sum(1 for t in traza if 'inviable' in t),
             error_busqueda=float(best['err']), tiempo_busqueda_s=dt, traza=traza, replicas=reps)
    for k in ('err_base', 'err_ext'):
        v = np.array([x[k] for x in reps]); r[k + '_media'] = float(v.mean()); r[k + '_sd'] = float(v.std(ddof=1))
    json.dump(r, open(f'{OUT}/res_{voi}_{FAM}.json', 'w'), indent=1)
    print(voi, FAM, 'err_ext %.4f +- %.4f | %s | n_eval %d (inviables %d) | %.0fs' % (
        r['err_ext_media'], r['err_ext_sd'], r['parametros'], r['n_eval'], r['n_inviables'], dt))


if __name__ == '__main__':
    main(sys.argv[1])
