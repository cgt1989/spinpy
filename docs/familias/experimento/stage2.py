"""Fase 2: propiedades FUERA de la funcion de error (no se optimizaron) para el ganador de cada familia.

Uso: python stage2.py VOI FAMILIA     (FAMILIA puede ser 'voi')
Lee out/masks_{VOI}_{FAMILIA}.npz (10 replicas del ganador a 96^3) y calcula, para las 3 primeras:
  * Conn.D, SMI, TbTh_CV, H_desv, fraccion de poro conexa (morfologia no ajustada)
  * homogeneizacion elastica periodica (hex8 + pyamg, la de spinpy) a res_mec 16 (3 replicas) y 32 (1 replica),
    en los EJES MATERIALES (E1>=E2>=E3, G, nu) y la distancia log-euclidea del tensor al del VOI.
La rigidez NO es citable como valor absoluto (README de spinpy); sirve para comparar VOI y candidatos al mismo res_mec.
"""
import sys, json, time, os
sys.path.insert(0, '../../..')
import numpy as np
from spinpy.io import leer_voi
from spinpy.morphometry import morfometria, metricas_forma
from spinpy.fit import fraccion_poro_conexa, metricas_mecanicas
from spinpy.elastic import remuestrear_bw, distancia_log_euclidea

OUT = os.environ.get('EXP_OUT', './out')


def fuera_de_muestra(BW, sp):
    m = morfometria(BW, sp, extra=True)
    f = metricas_forma(BW, sp, tbth=True, curvatura=True)
    fp, npor = fraccion_poro_conexa(BW)
    return dict(ConnD=float(m['ConnD']), SMI=float(m['SMI']), TbTh_CV=float(f['TbTh_CV']),
                H_desv=float(f['H_desv']), poro_conexo=fp, n_poros=int(npor),
                BVTV=float(m['BVTV']), TbTh=float(m['TbTh']), DA=float(m['DA']))


def mecanica(BW, sp, res_mec):
    t = time.time()
    r = metricas_mecanicas(BW, sp, res_mec)
    r['t_s'] = time.time() - t
    return r


def main(voi, fam):
    BW_voi, spacing = leer_voi(f'../../../datos/porcino_koria2020/VOI_{voi}.mat')
    spacing = np.asarray(spacing, float).ravel()
    res = {'voi': voi, 'familia': fam}
    if fam == 'voi':
        sp96 = None
        b96, s96 = remuestrear_bw(BW_voi, spacing, 96)
        res['morfo_fina'] = fuera_de_muestra(BW_voi, spacing)
        res['morfo_96'] = fuera_de_muestra(b96, s96)
        mecs = {16: [mecanica(BW_voi, spacing, 16)], 32: [mecanica(BW_voi, spacing, 32)]}
        sp_list = [spacing]
        res['morfo'] = [res['morfo_96']]
    else:
        masks = np.load(f'{OUT}/masks_{voi}_{fam}.npz')['BW']
        lc = float(np.mean(np.array(BW_voi.shape) * spacing))
        sp = np.full(3, lc / masks.shape[1])
        res['morfo'] = [fuera_de_muestra(masks[i], sp) for i in range(3)]
        mecs = {16: [mecanica(masks[i], sp, 16) for i in range(3)], 32: [mecanica(masks[0], sp, 32)]}
    res['mec'] = {str(k): v for k, v in mecs.items()}
    json.dump(res, open(f'{OUT}/stage2_{voi}_{fam}.json', 'w'), indent=1, default=float)
    print(voi, fam, 'listo')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
