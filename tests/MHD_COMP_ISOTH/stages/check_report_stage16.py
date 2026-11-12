#!/usr/bin/env python3
"""The numbers of REPORT_stage16.md against runs/stage16/{judge_tm,judge1_hf,compare_conv,conv_table}.json: each
number is formatted as in the report and must appear there (the costs: the rows of chapter 8 add up to the totals).

    python3 tests/MHD_COMP_ISOTH/stages/check_report_stage16.py [root]
"""
import os
import sys
import json
import paths

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(paths.RUNS, 'stage16')
REP = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'REPORT_stage16.md')
KEYS = ('E_B W1', 'E_B W2', 'E_m W1', 'E_m W2')
fails = []


def has(s, what):
    if s not in text:
        fails.append('%s: "%s" not in the report' % (what, s))


def e1(x):
    return ('%.1e' % x).replace('e-0', 'e-').replace('e+0', 'e+')


def f3(x):
    return ('%.3f' % x).replace('-', '−')


def sg(x):
    return ('%+.3f' % x).replace('-', '−')


def row(lab, cells):
    """A table row: label then the cells joined by ' | '."""
    has('| %s | %s |' % (lab, ' | '.join(cells)), 'row ' + lab)


text = open(REP).read()
J = json.load(open(os.path.join(ROOT, 'judge_tm.json')))
A, B, M = J['A'], J['B'], J['ms05']
labs = dict(zip(KEYS, ('E_B、W1', 'E_B、W2', 'E_m、W1', 'E_m、W2')))

# 4.1 checks
has(e1(J['checks']['mean']), 'mean check')
near = J['checks']['near']
has(e1(max(abs(x['t'] - x['t_near']) for x in near)), 'near dt')
has(e1(max(x['dlnE'] for x in near)), 'near dlnE')
over = [x['dlnE'] - x['bound'] for x in near if x['dlnE'] > x['bound']]
has('%d 件は上限以下、%d 件は上限を %s 以下' % (len(near) - len(over), len(over), e1(1.3e-15)), 'near counts')
assert max(over) <= 1.3e-15
has('%.4f' % J['delta_rerun_tm'], 'delta_rerun_tm')
has('%.4f' % J['delta_rerun_15'], 'delta_rerun_15')
has('%.4f' % J['delta_interp'], 'delta_interp')
has('%.4f' % max(v for k, v in J['delta_interp_all'].items() if k.startswith('ES ')), 'delta_interp ES')

# 4.2 variant A
for lab, v in sorted(A['cand'].items()):
    has(f3(v), 'cand A ' + lab)
assert A['B'] == 'H_r8' and A['verdict'] == 'ES better'
c = A['counts']
assert (c['plus'], c['minus'], c['d_es'], c['need']) == (3, 0, 4, 3)
for k in KEYS:
    t = A['table'][k + ' ES256']
    pe, pb = A['pairs']['ES'][k + ' ES256'], A['pairs']['H_r8'][k + ' ES256']
    row(labs[k], [f3(t['kES']), f3(t['kB']), sg(t['d']), f3(t['delta']), f3(pe['sigma']), f3(pb['sigma']),
                  f3(t['DES']), f3(t['DB'])])
has('ES で %s、H(8) で %s' % ('、'.join('%.1f' % A['pairs']['ES'][k + ' ES256']['k'] for k in KEYS),
                            '、'.join('%.1f' % A['pairs']['H_r8'][k + ' ES256']['k'] for k in KEYS)), 'k_x A')
cost = J['stage15']['cost']
has('ES が %.2f core-h、H(8) が %.2f core-h で、比は %.2f' % (cost['ES'], cost['B'], A['cost']['c_ratio']), 'cost')
has('k̄ は ES で %.2f、H(8) で %.2f で、(k̄_ES/k̄_B*)⁴ = %.2f' % (A['kbar']['ES'], A['kbar']['B'],
                                                            A['cost']['k_ratio4']), 'kbar')
has('%.2f > %.2f' % (A['cost']['c_ratio'], A['cost']['k_ratio4']), 'cost rule A')
assert not A['cost']['ok']

# 4.3 variant B
S15 = J['stage15']['ms1']
assert B['B'] == 'H_r8' and B['verdict'] == 'ES better' and S15['verdict'] == 'undecided'
has('H(8)（κ の平均 %s）' % f3(B['cand']['H_r8']), 'cand B')
for k in KEYS:
    row(labs[k], [sg(S15['table'][k + ' ES256']['d']), f3(B['table'][k + ' ES256']['delta']),
                  f3(S15['table'][k + ' ES256']['delta'])])
    assert B['table'][k + ' ES256']['delta'] == A['table'][k + ' ES256']['delta']
has('%.2f > %.2f' % (B['cost']['c_ratio'], B['cost']['k_ratio4']), 'cost rule B')

# 4.5 M_s = 0.5 and self-convergence
assert M['B'] == 'H_r0.25' and M['verdict'] == 'undecided'
assert (M['counts']['plus'], M['counts']['minus'], M['counts']['d_es'], M['counts']['d_b']) == (0, 3, 2, 2)
S05 = J['stage15']['ms05']['table']
for k in KEYS:
    t = M['table'][k]
    row(labs[k], [f3(t['kES']), f3(t['kB']), sg(t['d']), f3(t['delta']), f3(S05[k]['delta']), f3(t['DES']),
                  f3(t['DB'])])
has('ES が %.4f、H(1/4) が %.4f' % (M['table']['E_B W2']['DES'], M['table']['E_B W2']['DB']), 'D ms05 E_B W2')
sp = J['self_ms1_64']['pairs']
for lab, name in (('ES', 'ES'), ('H_r8', 'H(8)'), ('L_r1', 'L(1)')):
    has('%s で %s' % (name, '、'.join(f3(sp[lab][k]['kappa']) for k in KEYS)), 'self64 ' + lab)
has('L(2) で 4 組とも %s' % f3(sp['L_r2']['E_B W1']['kappa']), 'self64 L_r2')
assert len({round(sp['L_r2'][k]['kappa'], 3) for k in KEYS}) == 1
has('（%s）' % '、'.join(f3(J['self_ms1_128']['pairs']['ES'][k]['kappa']) for k in KEYS), 'self128 ES')
for k in KEYS:
    assert abs(J['self_ms1_128']['pairs']['ES'][k]['kappa'] - A['table'][k + ' ES256']['kES']) < 1e-12

# 4.4 E_m W2 kappa_j
for lab, name in (('ES', 'ES'), ('H_r8', 'H(8)')):
    has('%s で %s' % (name, '、'.join(f3(x) for x in A['pairs'][lab]['E_m W2 ES256']['kappas'])), 'kappa_j ' + lab)

# 5.3 judge1
J1 = json.load(open(os.path.join(ROOT, 'judge1_hf.json')))['A']
assert J1['B'] == 'H_r8' and J1['verdict'] == 'ES better' and J1['table'] == A['table']
has('H(8) が %s、HF(8) が %s、L(1) が %s、L(2) が %s' % tuple(f3(J1['cand'][x]) for x in ('H_r8', 'HF_r8', 'L_r1',
                                                                                   'L_r2')), 'cand judge1')
hf = [J1['pairs']['HF_r8'][k + ' ES256']['kappa'] for k in KEYS]
h = [A['table'][k + ' ES256']['kB'] for k in KEYS]
has('%s で、H(8) の %s と %.3f 以内' % ('、'.join(f3(x) for x in hf), '、'.join(f3(x) for x in h), 0.007),
    'judge1 kappa')
assert max(abs(a - b) for a, b in zip(hf, h)) <= 0.007
di = json.load(open(os.path.join(ROOT, 'judge1_hf.json')))['delta_interp_all']
has('%.4f' % max(v for k, v in di.items() if k.startswith('HF_r8')), 'judge1 delta_interp')
has('κ_j は %s' % '、'.join(f3(x) for x in J1['pairs']['HF_r8']['E_m W2 ES256']['kappas']), 'kappa_j HF')


# 3 conv
def e2(x):
    return ('%.2e' % x).replace('e-0', 'e-').replace('e+0', 'e+')


def o2(x):
    return ('%.2f' % x).replace('-', '−')


def rng(xs, f, sep='–'):
    return f(min(xs)) + sep + f(max(xs))


C = json.load(open(os.path.join(ROOT, 'compare_conv.json')))
assert C['ok'] == {'CV1': True, 'CV2': True}
for name, c in C['cases'].items():
    assert c['complete'] and c['lines'] == [21, 21] and c['old_max'] == c['old_end'] == c['dt_max'] == 0.0, name
    assert c['last_same_values'] and c['last_same_body'], name
    assert c['last_same_bytes'] == (c['runs'] == 1) == (name in ('I1_N4096', 'I1_N4096_eta2')), name
T = json.load(open(os.path.join(ROOT, 'conv_table.json')))
CASES = (('I1', 'I1、bulk 粘性'), ('I1_eta2', 'I1、C_η = 2'), ('RDpc', 'RDpc、bulk 粘性'),
         ('RDpc_eta2', 'RDpc、C_η = 2'))
for key, lab in CASES:
    new = T[key]['new']
    for i, N in enumerate(sorted(new, key=int)):
        r = new[N]
        cells = [N, e2(r['fv']), o2(r['order_fv']) if 'order_fv' in r else '–', e2(r['self']) if 'self' in r else '–',
                 o2(r['order_self']) if 'order_self' in r else '–']
        has(('| %s | ' % lab if i == 0 else '| | ') + ' | '.join(cells) + ' |', 'conv row %s %s' % (key, N))
assert e2(T['I1_fv']['I1_N4096.npz']) == '9.50e-4' and e2(T['RDpc_fv']['RDpc_N4096.npz']) == '1.53e-3'
has('I1 で 9.50e-4、RDpc で 1.53e-3', 'FV self')
old_N = ('256', '512', '1024', '2048')
md = [v for key, _ in CASES for N, v in T[key]['maxdiff'].items() if N in old_N]
has('9.8e-5–4.1e-3', 'maxdiff range')
assert e1(min(md)) + '–' + e1(max(md)) == '9.8e-5–4.1e-3'
assert all(v == 0.0 for key, _ in CASES for N, v in T[key]['maxdiff'].items() if N not in old_N)


def reldiff(q, keys):
    return max(abs(T[k]['new'][N][q] - T[k]['old'][N][q])/abs(T[k]['old'][N][q])
               for k in keys for N in old_N if N in T[k]['new'] and N in T[k]['old'])


has('FV との ρ の L1 で %.2f%% 以下、I1 の理想解との L1 で %.2f%% 以下、領域ごとの L1 で %.2f%% 以下' % (
    100*reldiff('fv', [k for k, _ in CASES]), 100*reldiff('ideal', ('I1', 'I1_eta2')),
    100*max(abs(T[k]['new'][N]['regions'][g] - T[k]['old'][N]['regions'][g])/T[k]['old'][N]['regions'][g]
            for k, _ in CASES for N in old_N if N in T[k]['new'] and N in T[k]['old']
            for g in T[k]['new'][N]['regions'])), 'reldiff')
has('前方の ρ の差は %s 以下' % e1(max(abs(T[k]['new'][N]['ahead'][0] - T[k]['old'][N]['ahead'][0])
                                  for k in ('I1', 'I1_eta2') for N in old_N)), 'ahead diff')
has('%s から %s' % (e2(T['I1']['old']['2048']['fv']), e2(T['I1']['new']['2048']['fv'])), 'I1 2048 old new')
nst = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stage3_nstep.json')))
n1 = [nst['I1_N%d' % N][0] for N in (256, 512, 1024, 2048)]
has('、'.join('%d' % n for n in n1 + [2*n1[-1], 4*n1[-1]]), 'nstep I1')
I1, I1e, RD, RDe = (T[k]['new'] for k, _ in CASES)
has(rng([I1[N]['fv'] for N in ('1024', '2048', '4096', '8192')], e2), 'I1 fv range')
has(rng([RD[N]['fv'] for N in RD], e2), 'RDpc fv range')
has(rng([I1[N]['order_fv'] for N in ('1024', '2048', '4096')] + [RD[N]['order_fv'] for N in RD if 'order_fv' in RD[N]],
        o2, ' から '), 'bulk order_fv range')
has('次数 %s で縮む' % rng([r['order_self'] for d in (I1, RD) for r in d.values() if 'order_self' in r], o2),
    'bulk order_self range')
has('約 %.1f 倍（I1）と約 %.1f 倍（RDpc）' % (I1['8192']['fv']/I1['4096']['self'], RD['8192']['fv']/RD['4096']['self']),
    'bulk ratio')
has('%s、%s、%s から %s、%s' % tuple(o2(I1e[N]['order_fv']) for N in ('256', '512', '1024', '2048', '4096')),
    'I1 eta2 order_fv')
has('次数は %s で' % rng([r['order_self'] for r in I1e.values() if 'order_self' in r], o2), 'I1 eta2 order_self')
has('距離 %s は、N = 4096 と 8192 の距離 %s と同じ' % (e2(I1e['8192']['fv']), e2(I1e['4096']['self'])), 'I1 eta2 8192')
has('%s、%s と下がり、N = 4096 → 8192 では %s' % (o2(RDe['256']['order_fv']), o2(RDe['512']['order_fv']),
                                                 o2(RDe['4096']['order_fv'])), 'RDpc eta2 order_fv')
has('距離 %s は、N = 4096 と 8192 の距離 %s の %.1f 倍' % (e2(RDe['8192']['fv']), e2(RDe['4096']['self']),
                                                      RDe['8192']['fv']/RDe['4096']['self']), 'RDpc eta2 ratio')
has('%s しか求められない' % o2(RDe['256']['order_self']), 'RDpc eta2 order_self')

# 8 cost: the rows add up to the total and to the numbers of the header
cost_rows = [r for r in (ln.split('|') for ln in text.splitlines() if ln.startswith(('| conv ', '| F')))
             if len(r) == 7 and r[2].strip().isdigit()]
conv_h = sum(float(r[5]) for r in cost_rows if r[1].startswith(' conv '))
hf_h = sum(float(r[5]) for r in cost_rows if r[1].startswith(' F'))
has('| 計 | | | | %.2f |' % (conv_h + hf_h), 'cost total')
has('ジョブ %d 本、合計 %.2f core-h（conv %.2f、HF %.2f。' % (len(cost_rows), conv_h + hf_h, conv_h, hf_h), 'cost header')
has('HF の 3 本の合計は %.2f core-h' % hf_h, 'cost HF')

print('\n'.join(fails) if fails else 'all numbers found')
sys.exit(1 if fails else 0)
