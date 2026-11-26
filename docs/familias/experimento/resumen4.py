import sys, json, os, numpy as np
sys.path.insert(0, '../../..')
from spinpy.elastic import distancia_log_euclidea
S = './'
VOIS = ['C1', 'C5', 'V3', 'V2']
FAMS = ['spinodoide', 'dual', 'open_foam', 'levelset', 'levelset_tpl', 'levelset_shape']
CARP = {f: 'out3' for f in FAMS}; CARP['levelset_shape'] = 'out4'


def J(carpeta, pref, v, f):
    return json.load(open(S + f'{carpeta}/{pref}_{v}_{f}.json'))


def dist_mec(ref, L):
    ds = [distancia_log_euclidea(np.array(ref['C_principal_rel']), np.array(x['C_principal_rel']))
          for x in L if 'C_principal_rel' in x]
    return float(np.mean(ds)) if ds else np.nan


def med(v, f):
    s2 = J(CARP[f], 'stage2', v, f); vo = J('out3', 'stage2', v, 'voi'); mv = vo['morfo_96']; m = s2['morfo']
    d = {k: np.mean([x[k] for x in m]) for k in ['ConnD', 'SMI', 'TbTh_CV', 'H_desv', 'n_poros']}
    return dict(ConnD=abs(d['ConnD'] - mv['ConnD']) / mv['ConnD'], CV=abs(d['TbTh_CV'] - mv['TbTh_CV']) / mv['TbTh_CV'],
                H=abs(d['H_desv'] - mv['H_desv']) / mv['H_desv'], SMI=abs(d['SMI'] - mv['SMI']),
                poros=abs(np.log((d['n_poros'] + 1) / (mv['n_poros'] + 1))),
                d16=dist_mec(vo['mec']['16'][0], s2['mec']['16']), d32=dist_mec(vo['mec']['32'][0], s2['mec']['32']),
                abs_=d)


T = {(v, f): dict(r=J(CARP[f], 'res', v, f), m=med(v, f)) for v in VOIS for f in FAMS}

print('=== levelset_shape: parametros ganadores y rendimiento por VOI')
for v in VOIS:
    r = T[(v, 'levelset_shape')]['r']; p = r['parametros']; m = T[(v, 'levelset_shape')]['m']; vo = J('out3', 'stage2', v, 'voi')['morfo_96']
    print(f"{v}: N={p['N']:.0f} estir={p['estiramiento']} K_c={p['curvatura']:+.2f} c1={p['c1']:.0f} c2={p['c2']:.1f} | err_ext {r['err_ext_media']:.4f}±{r['err_ext_sd']:.4f} "
          f"(busqueda {r['error_busqueda']:.4f}) | inviables {r['n_inviables']}/{r['n_eval']} | "
          f"VOI SMI {vo['SMI']:.2f} poros {vo['n_poros']} -> SMI {m['abs_']['SMI']:.2f} poros {m['abs_']['n_poros']:.0f} d16 {m['d16']:.3f} d32 {m['d32']:.3f}")

print('\n=== Error ampliado por VOI (ronda con densidad fija; levelset_shape = forma ajustada)')
print('%-15s ' % '' + ' '.join('%8s' % v for v in VOIS) + '   media')
for f in FAMS:
    e = [T[(v, f)]['r']['err_ext_media'] for v in VOIS]
    print('%-15s ' % f + ' '.join('%8.4f' % x for x in e) + '   %.4f' % np.mean(e))

ks = ['ConnD', 'CV', 'H', 'SMI', 'poros', 'd16', 'd32']
print('\n=== Desviacion media respecto al VOI (4 VOI). EN el objetivo: ConnD, CV, H. FUERA: SMI, poros, d16, d32')
print('%-15s %8s ' % ('', 'err_ext') + ' '.join('%7s' % k for k in ks))
for f in FAMS:
    print('%-15s %8.4f ' % (f, np.mean([T[(v, f)]['r']['err_ext_media'] for v in VOIS])) +
          ' '.join('%7.3f' % np.nanmean([T[(v, f)]['m'][k] for v in VOIS]) for k in ks))

print('\n=== levelset_shape frente a levelset (forma por defecto) y a tpl (forma de la plantilla), por VOI: diferencia (shape - otra); negativo = shape mejor')
for ref in ('levelset', 'levelset_tpl'):
    print('-- frente a', ref)
    print('%-5s ' % '' + ' '.join('%8s' % k for k in ['err_ext'] + ks))
    for v in VOIS:
        a, b = T[(v, 'levelset_shape')], T[(v, ref)]
        print('%-5s ' % v + ' %+8.4f' % (a['r']['err_ext_media'] - b['r']['err_ext_media']) + ' '.join('%+8.3f' % (a['m'][k] - b['m'][k]) for k in ks))

cr = ['err_ext', 'ConnD', 'CV', 'H', 'SMI', 'poros', 'd32']
def val(v, f, c):
    return T[(v, f)]['r']['err_ext_media'] if c == 'err_ext' else T[(v, f)]['m'][c]
print('\n=== En cuantos de los 4 VOI es MEJOR levelset_shape que cada otra familia')
print('%-15s ' % 'frente a' + ' '.join('%7s' % c for c in cr))
for ref in FAMS[:-1]:
    print('%-15s ' % ref + ' '.join('%5d/4' % sum(val(v, 'levelset_shape', c) < val(v, ref, c) for v in VOIS) for c in cr))

print('\n=== Control de sobreajuste: error de la semilla de busqueda vs semillas nuevas (media de 4 VOI)')
for f in ('levelset', 'levelset_tpl', 'levelset_shape'):
    bs = np.mean([T[(v, f)]['r']['error_busqueda'] for v in VOIS]); nv = np.mean([T[(v, f)]['r']['err_ext_media'] for v in VOIS])
    print('%-15s busqueda %.4f  semillas nuevas %.4f  cociente %.2f' % (f, bs, nv, nv / bs))
print('\n=== Aportacion de cada etapa en levelset_shape (mejor err_ext acumulado, por VOI)')
for v in VOIS:
    tr = T[(v, 'levelset_shape')]['r']['traza']; row = []
    for et in 'ABCDEF':
        sub = [t['err_ext'] for t in tr if t['etapa'] in 'ABCDEF'[:'ABCDEF'.index(et) + 1] and t['err_ext'] < 1e5]
        row.append('%s %.4f' % (et, min(sub)))
    print(v, ' | '.join(row))
