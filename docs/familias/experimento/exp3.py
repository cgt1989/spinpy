"""Experimento 3: objetivo AMPLIADO (objetivo.py) con la DENSIDAD FIJADA al BV/TV medido del VOI, en las 5 familias.

Uso: python exp2.py VOI FAMILIA     (FAMILIA: spinodoide | dual | open_foam | levelset | levelset_tpl)
Objetivo = error base (con DA2) + PoDm + TbTh_CV + H_desv + ConnD (peso 1 cada uno, TbSp a 0).
Las propiedades NO incluidas en el objetivo (SMI, n. de poros, mecanica) se miden en stage2.py.
"""
import sys, json, time, math, os
sys.path.insert(0, '.')
sys.path.insert(0, '../../..')
import numpy as np
import objetivo as OBJ
PES = OBJ.activar()                      # parches en memoria ANTES de usar los ajustes de spinpy

import spinpy.error as E
import spinpy.fit as _F, spinpy.fit_dual as _FD
PIN = {'rho': None}                      # BV/TV medido del VOI en curso (se fija en main)
_gm, _gdl = _F.generar_mascara, _FD.generar_dual_lattice
def _gm_pin(*a, **kw):
    kw['rho'] = PIN['rho']; return _gm(*a, **kw)
def _gdl_pin(resolution, celdas, rho, *a, **kw):
    return _gdl(resolution, celdas, PIN['rho'], *a, **kw)
_F.generar_mascara = _gm_pin; _FD.generar_dual_lattice = _gdl_pin
from spinpy.io import leer_voi
from spinpy.morphometry import morfometria
from spinpy.fit import (ajustar_spinodoide, longitud_caracteristica, resolucion_comparacion,
                        SEMILLA_POR_DEFECTO, RHO_MAX, RHO_MIN)
from spinpy.fit_dual import ajustar_dual_lattice, PRESETS_COMPLETO
from spinpy.grf import generar_mascara
from spinpy.dual_lattice import generar_dual_lattice
import exp as X1                         # reutiliza generadores de morphoGen (open_foam, level set)

OUT = os.environ.get('EXP_OUT', './out3')
os.makedirs(OUT, exist_ok=True)
SEMILLAS_NUEVAS = list(range(1001, 1011))
ESCALARES = ['BVTV', 'TbTh', 'TbSp', 'TbN', 'BSBV', 'PoTot', 'DA', 'DA2', 'TbTh_CV', 'H_desv', 'PoDm', 'ConnD']


def errores(m_voi, m):
    e_b, _ = E.error_morfometrico(m_voi, m, peso_da2=1.0)                 # objetivo del piloto (con DA2)
    e_x, _ = E.error_morfometrico(m_voi, m, peso_da2=1.0, pesos=PES)     # objetivo ampliado
    return float(e_b), float(e_x)


def escalares(m):
    return {k: float(m[k]) for k in ESCALARES if k in m and np.ndim(m[k]) == 0}


def buscar(familia, m_voi_ext, res, sp_cand, bvtv):
    gen = X1.GENERADORES[familia]
    dens = np.array([bvtv])           # densidad fijada
    tamanos = np.unique(np.round(np.geomspace(50, 800, 6)))
    traza, t0 = [], time.time()

    def evaluar(d, N, e, etapa):
        try:
            BW = gen(res, d, N, e, SEMILLA_POR_DEFECTO)
        except ValueError as exc:
            traza.append(dict(etapa=etapa, d=float(d), N=float(N), e=[float(x) for x in e],
                              err_ext=1e6, inviable=str(exc)[:80]))
            return 1e6
        m = OBJ.medir_ext(BW, sp_cand)
        e_b, e_x = errores(m_voi_ext, m)
        traza.append(dict(etapa=etapa, d=float(d), N=float(N), e=[float(x) for x in e],
                          err_base=e_b, err_ext=e_x))
        return e_x

    best = dict(err=np.inf, d=dens[0], N=tamanos[0], e=(1.0, 1.0, 1.0))
    for d in dens:
        for N in tamanos:
            er = evaluar(d, N, best['e'], 'A')
            if er < best['err']:
                best.update(err=er, d=d, N=N)
    for e in PRESETS_COMPLETO:
        er = evaluar(best['d'], best['N'], e, 'B')
        if er < best['err']:
            best.update(err=er, e=e)
    d_ref = np.array([bvtv])
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
    m_voi = morfometria(BW_voi, spacing)                   # base: VOI a resolucion nativa (como spinpy)
    m_voi_ext = dict(m_voi)
    m_voi_ext.update(OBJ.metricas_forma_ext(BW_voi, spacing, tbth=True, poro=True, curvatura=True,
                                            conectividad=True))   # terminos nuevos: VOI al voxel del candidato
    lc = longitud_caracteristica(BW_voi.shape, spacing)
    res = resolucion_comparacion(BW_voi.shape)
    sp_cand = np.full(3, lc / res)
    bvtv = float(m_voi['BVTV']); PIN['rho'] = bvtv
    res_json = dict(voi=voi, familia=familia, res=res, lc=lc, voi_escalares=escalares(m_voi_ext), pesos=PES)
    t0 = time.time()

    if familia == 'spinodoide':
        r = ajustar_spinodoide(BW_voi, spacing, modo='completo', m_voi=m_voi, pesos=PES)
        p = r['parametros']
        th = np.clip(np.round(p['thetas']), 0, 90)
        R = np.asarray(p['R'])
        gen = lambda s: generar_mascara(resolution=res, wave_number=p['wave_number_rad'],
                                        num_waves=p['num_waves'], thetas=th, rho=bvtv,
                                        R=R, esquema=p['esquema'], seed=int(s), precision='f32')[0]
        res_json['parametros'] = {k: p[k] for k in ('densidad', 'wave_number_pi', 'num_waves', 'thetas', 'esquema')}
        res_json['n_eval'] = r['n_evaluaciones']; res_json['error_busqueda'] = r['error_busqueda']
        res_json['tiempo_busqueda_s'] = r['tiempo_s']
    elif familia == 'dual':
        r = ajustar_dual_lattice(BW_voi, spacing, modo='completo', m_voi=m_voi, pesos=PES)
        p = r['parametros']
        R = np.asarray(p['R'])
        gen = lambda s: generar_dual_lattice(res, p['celdas'], bvtv, estiramiento=p['estiramiento'],
                                             irregularidad=p['irregularidad'], R=R, seed=int(s))[0]
        res_json['parametros'] = {k: p[k] for k in ('densidad', 'celdas', 'estiramiento', 'irregularidad')}
        res_json['n_eval'] = r['n_evaluaciones']; res_json['error_busqueda'] = r['error_busqueda']
        res_json['tiempo_busqueda_s'] = r['tiempo_s']
    else:
        best, traza, dt = buscar(familia, m_voi_ext, res, sp_cand, bvtv)
        g = X1.GENERADORES[familia]
        gen = lambda s: g(res, best['d'], best['N'], best['e'], s)
        res_json['parametros'] = dict(densidad=float(best['d']), N=float(best['N']),
                                      estiramiento=[float(x) for x in best['e']])
        res_json['n_eval'] = len(traza); res_json['error_busqueda'] = float(best['err'])
        res_json['tiempo_busqueda_s'] = dt
        res_json['traza'] = traza
        res_json['n_inviables'] = sum(1 for t in traza if 'inviable' in t)

    reps, masks = [], []
    for s in SEMILLAS_NUEVAS:
        BW = gen(s)
        m = OBJ.medir_ext(BW, sp_cand)
        e_b, e_x = errores(m_voi_ext, m)
        reps.append(dict(semilla=s, err_base=e_b, err_ext=e_x, **escalares(m)))
        masks.append(BW)
    np.savez_compressed(f'{OUT}/masks_{voi}_{familia}.npz', BW=np.stack(masks))
    res_json['replicas'] = reps
    for k in ('err_base', 'err_ext'):
        v = np.array([r_[k] for r_ in reps])
        res_json[k + '_media'] = float(v.mean()); res_json[k + '_sd'] = float(v.std(ddof=1))
    res_json['tiempo_total_s'] = time.time() - t0
    json.dump(res_json, open(f'{OUT}/res_{voi}_{familia}.json', 'w'), indent=1)
    print(voi, familia, 'err_ext %.4f +- %.4f | err_base %.4f +- %.4f | n_eval %d | %.0fs' % (
        res_json['err_ext_media'], res_json['err_ext_sd'], res_json['err_base_media'],
        res_json['err_base_sd'], res_json['n_eval'], res_json['tiempo_total_s']))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
