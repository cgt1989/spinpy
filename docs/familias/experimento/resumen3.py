import sys, json, os, numpy as np
sys.path.insert(0, '../../..')
from spinpy.elastic import distancia_log_euclidea
S = './'
VOIS = ['C1', 'C5', 'V3', 'V2']
F = ['spinodoide', 'dual', 'open_foam', 'levelset', 'levelset_tpl']
NUEVAS = ['open_foam', 'levelset', 'levelset_tpl']
R = {'R1': 'out', 'R2': 'out2', 'R3': 'out3'}


def J(carpeta, pref, v, f):
    p = S + f'{carpeta}/{pref}_{v}_{f}.json'
    return json.load(open(p)) if os.path.exists(p) else None


def dist_mec(ref, L):
    ds = [distancia_log_euclidea(np.array(ref['C_principal_rel']), np.array(x['C_principal_rel']))
          for x in L if 'C_principal_rel' in x]
    return float(np.mean(ds)) if ds else np.nan


def medidas(carpeta, v, f):
    s2 = J(carpeta, 'stage2', v, f)
    vo = J('out3', 'stage2', v, 'voi')
    if s2 is None:
        return None
    mv = vo['morfo_96']; m = s2['morfo']
    d = {k: np.mean([x[k] for x in m]) for k in ['ConnD', 'SMI', 'TbTh_CV', 'H_desv', 'n_poros', 'BVTV']}
    return dict(ConnD=abs(d['ConnD'] - mv['ConnD']) / mv['ConnD'],
                CV=abs(d['TbTh_CV'] - mv['TbTh_CV']) / mv['TbTh_CV'],
                H=abs(d['H_desv'] - mv['H_desv']) / mv['H_desv'],
                SMI=abs(d['SMI'] - mv['SMI']),
                poros=abs(np.log((d['n_poros'] + 1) / (mv['n_poros'] + 1))),
                d16=dist_mec(vo['mec']['16'][0], s2['mec']['16']),
                d32=dist_mec(vo['mec']['32'][0], s2['mec']['32']),
                _abs=dict(ConnD=d['ConnD'], SMI=d['SMI'], CV=d['TbTh_CV'], H=d['H_desv'], poros=d['n_poros']))


def dBV(carpeta, v, f):
    r = J(carpeta, 'res', v, f)
    if r is None:
        return np.nan
    vb = r['voi_escalares']['BVTV']
    return 100 * (np.mean([x['BVTV'] for x in r['replicas']]) - vb) / vb


def err_ext(carpeta, v, f):
    r = J(carpeta, 'res', v, f)
    return None if r is None else r['err_ext_media']


T = {}
for v in VOIS:
    vo = J('out3', 'stage2', v, 'voi')['morfo_96']
    print(f"\n=== {v} (ronda 3, densidad fija)  VOI: ConnD {vo['ConnD']:.1f}  CV {vo['TbTh_CV']:.2f}  H {vo['H_desv']:.2f}  SMI {vo['SMI']:.2f}  poros {vo['n_poros']}")
    print('%-13s %-15s %-8s | %6s %5s %5s | %6s %6s %7s %7s' % ('familia', 'err_ext', 'dBVTV%', 'ConnD', 'CV', 'H', 'SMI', 'poros', 'dL16', 'dL32'))
    for f in F:
        r = J('out3', 'res', v, f); m = medidas('out3', v, f)
        if r is None or m is None:
            print('%-13s (pendiente)' % f); continue
        T[(v, f)] = dict(r=r, m=m)
        a = m['_abs']
        print('%-13s %.4f±%.4f   %+6.1f   | %6.1f %5.2f %5.2f | %6.2f %6.0f %7.3f %7.3f' % (
            f, r['err_ext_media'], r['err_ext_sd'], dBV('out3', v, f), a['ConnD'], a['CV'], a['H'], a['SMI'], a['poros'], m['d16'], m['d32']))

if len(T) == len(VOIS) * len(F):
    print('\n=== Desviacion de BV/TV del ganador respecto al VOI (%), media de |.| y maximo, por ronda')
    for rn, c in R.items():
        print(rn, {f: '%.1f (max %.1f)' % (np.mean([abs(dBV(c, v, f)) for v in VOIS]), max(abs(dBV(c, v, f)) for v in VOIS)) for f in F})

    print('\n=== Ronda 3, error ampliado por VOI (mejor en negrita en el informe)')
    print('%-13s ' % '' + ' '.join('%8s' % v for v in VOIS) + '   media')
    for f in F:
        e = [T[(v, f)]['r']['err_ext_media'] for v in VOIS]
        print('%-13s ' % f + ' '.join('%8.4f' % x for x in e) + '   %.4f' % np.mean(e))

    ks = ['ConnD', 'CV', 'H', 'SMI', 'poros', 'd16', 'd32']
    print('\n=== Ronda 3: desviacion media respecto al VOI (4 VOI). EN el objetivo: ConnD, CV, H. FUERA: SMI, poros, d16, d32')
    print('%-13s %8s ' % ('', 'err_ext') + ' '.join('%7s' % k for k in ks))
    for f in F:
        print('%-13s %8.4f ' % (f, np.mean([T[(v, f)]['r']['err_ext_media'] for v in VOIS])) +
              ' '.join('%7.3f' % np.nanmean([T[(v, f)]['m'][k] for v in VOIS]) for k in ks))

    print('\n=== Comparacion de rondas en lo que NUNCA estuvo en el objetivo (media de 4 VOI): SMI | poros | d16 | d32')
    M = {rn: {f: {v: medidas(c, v, f) for v in VOIS} for f in F} for rn, c in R.items()}
    for f in F:
        cols = []
        for k in ['SMI', 'poros', 'd16', 'd32']:
            cols.append(' / '.join('%.2f' % np.nanmean([M[rn][f][v][k] for v in VOIS]) for rn in R))
        print('%-13s  ' % f + '   |   '.join(cols))
    print('(cada grupo: ronda1 / ronda2 / ronda3)')

    cr = ['err_ext', 'ConnD', 'CV', 'H', 'SMI', 'poros', 'd32']
    def val(v, f, c):
        return T[(v, f)]['r']['err_ext_media'] if c == 'err_ext' else T[(v, f)]['m'][c]
    print('\n=== Ronda 3: en cuantos de los 4 VOI es MEJOR cada variante de morphoGen (menor desviacion)')
    print('%-13s %-10s ' % ('', 'frente a') + ' '.join('%7s' % c for c in cr))
    for f in NUEVAS:
        for ref in ['spinodoide', 'dual']:
            print('%-13s %-10s ' % (f, ref) + ' '.join('%5d/4' % sum(val(v, f, c) < val(v, ref, c) for v in VOIS) for c in cr))
    print('\n=== Ronda 3: la mas cercana de las 5 en cada criterio (VOI)')
    print('%-13s ' % '' + ' '.join('%-12s' % c for c in cr))
    for f in F:
        print('%-13s ' % f + ' '.join('%-12s' % (','.join(v for v in VOIS if min(F, key=lambda g: val(v, g, c)) == f) or '-') for c in cr))
    print('\n=== Parametros ganadores, ronda 3')
    for v in VOIS:
        print(v, {f: T[(v, f)]['r']['parametros'] for f in ('dual', 'levelset', 'open_foam')})
