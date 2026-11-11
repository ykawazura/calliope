#!/usr/bin/env python3
"""Stage 16: the 1D conv runs of Stage 3 again with the final binary (gen_stage16.py), and K of Stage 15 again with
time-matched spectra (a post hoc judgement). The criteria and the definitions were fixed in tasks/todo.md (Stage 16)
before the runs and before the judgement.

    python3 analyze_stage16.py selftest        the functions below on synthetic data
    python3 analyze_stage16.py conv [root]     CV1 (completion) and CV2 (all 21 lines against runs/stage3_f/conv)
                                               of <root>/conv: <root>/compare_conv.json and the last states in
                                               <root>/conv_last (ps1)
    python3 analyze_stage16.py table [root]    the conv table of s13d (N <= 2048: runs/stage14/mhd1d/*_last; 4096,
                                               8192: <root>/conv_last) beside the earlier one (stage3_aux.conv), with
                                               the rho L1 between successive N and the orders: <root>/conv_table.json
                                               and <root>/fig_s16_conv.png
    python3 analyze_stage16.py rejudge [root]  K of Stage 15 from runs/stage15/spec.json with the time-matched pairs
                                               (variant A: kappa, D and sigma; variant B: sigma only), the reference
                                               judgement of M_s = 0.5 and the self-convergence reports:
                                               <root>/judge_tm.json and <root>/fig_s16_kappa_t.png
    python3 analyze_stage16.py hf <batch> [root]
                                               C and E of the HF runs of <root>/<batch> (analyze_stage15.batch_runs:
                                               <root>/<batch>/summary.json) and the gate of the batch in
                                               <root>/<batch>.json: hf64 r_star (the smallest r from which every
                                               larger rung is complete), hf128 r_C (the smallest complete rung); no
                                               gate while a run is pending or not run (ps1; the stop rule of the cost is
                                               analyze_stage15.py watch <batch> <root>)
    python3 analyze_stage16.py spec_hf [root]  per-record spectra of the complete HF runs of <root>/hf128 and
                                               <root>/hf256 (analyze_stage15.run_spec): <root>/spec_hf.json (ps1)
    python3 analyze_stage16.py judge2 [root]   K at M_s = 1 on 128^3 against ES-256 and HF-256 (8 pairs, need = 6)
                                               from runs/stage15/spec.json and <root>/spec_hf.json: variant A,
                                               variant B and the definitions of Stage 15 (reference only), and the
                                               (128, 256) self-convergence of ES and HF(r_C): <root>/judge_hf.json
    python3 analyze_stage16.py judge1 [root]   variant A of rejudge (ES-256 the only reference) with the complete HF
                                               runs of <root>/spec_hf.json among the 128^3 runs (a sensitivity fixed
                                               after H-F was stopped; the verdict of judge_tm.json stays):
                                               <root>/judge1_hf.json

Time matching. The records of a run X in a window are sorted by t, and a record within DUP of the previous kept one is
left out; ln E of each shell is linear in t between the kept records. E_X is taken at the record times t_j of the
reference R with the weights w_j of R (analyze_stage15.window_weights): the window means are sum_j w_j E(t_j), kappa_j
comes from E_X(t_j)/E_R(t_j), and sigma is the weighted standard deviation of kappa_j. The default root is
runs/stage16.
"""
import os
import re
import sys
import glob
import json
import math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_stage3 as a3
import analyze_stage4 as a4
import analyze_stage14 as a14
import analyze_stage15 as a15
import gen_stage15 as g15
import gen_stage16 as g16
from analyze_stage1 import run_dirs

DUP = 1e-3          # tau0: a record of X within DUP of the previous kept one is left out
S15 = os.path.join('runs', 'stage15')
S14_LAST = {'': os.path.join('runs', 'stage14', 'mhd1d', 'suite_last'),
            '_eta2': os.path.join('runs', 'stage14', 'mhd1d', 'eta_suite_last')}
N_SUITE = (256, 512, 1024, 2048)
N_CONV = (4096, 8192)


# ---------------------------------------------------------------- time-matched pairs
def series(run, W, q):
    """(t, E) of the window records of run: sorted by t, a record within DUP of the previous kept one left out."""
    d = run['win'][W]
    keep = []
    for i in sorted(range(len(d['t'])), key=lambda i: d['t'][i]):
        if not keep or d['t'][i] - d['t'][keep[-1]] > DUP:
            keep.append(i)
    return np.array([d['t'][i] for i in keep], float), np.array([d['recs'][q][i] for i in keep], float)


def half(t, E):
    """Every other record, both ends kept (Delta_interp)."""
    i = list(range(0, len(t), 2))
    if i[-1] != len(t) - 1:
        i.append(len(t) - 1)
    return t[i], E[i]


def at(t, E, tj):
    """E (records x shells) at the times tj: the record itself at its own time and the nearest record within WIN_TOL
    outside [t[0], t[-1]]; between two records ln E linear in t (nan where E <= 0); an error further out."""
    out = []
    for x in tj:
        if x < t[0] - a15.WIN_TOL or x > t[-1] + a15.WIN_TOL:
            raise ValueError('t = %.6f is outside the records [%.6f, %.6f] by more than WIN_TOL' % (x, t[0], t[-1]))
        i = int(np.searchsorted(t, x))
        if x <= t[0]:
            out.append(E[0])
        elif x >= t[-1]:
            out.append(E[-1])
        elif x == t[i]:
            out.append(E[i])
        else:
            a = (x - t[i - 1])/(t[i] - t[i - 1])
            with np.errstate(divide='ignore', invalid='ignore'):
                lo, hi = (np.where(E[k] > 0, np.log(E[k]), np.nan) for k in (i - 1, i))
            out.append(np.exp((1 - a)*lo + a*hi))
    return np.array(out)


def wsigma(k, w):
    """sqrt(sum w (k - kbar)^2/(1 - sum w^2)), w normalized to 1 and kbar = sum w k (ddof = 1 for equal weights);
    0 for one record."""
    k, w = np.asarray(k, float), np.asarray(w, float)
    w = w/np.sum(w)
    s2 = 1.0 - np.sum(w**2)
    if len(k) < 2 or s2 <= 0:
        return 0.0
    return float(math.sqrt(np.sum(w*(k - np.sum(w*k))**2)/s2))


def pair_tm(run, ref, q, W, cap, d_hi, kc, every_other=False):
    """analyze_stage15.pair with the time-matched spectra (variant A): k_x, kappa, censoring and D of the window
    means, sigma and kappa_j at the record times t_j of ref; every_other takes every other record of run (half)."""
    R = ref['win'][W]
    tj, w = np.array(R['t'], float), np.array(R['w'], float)
    Er = np.array(R['recs'][q], float)
    t, E = series(run, W, q)
    if every_other:
        t, E = half(t, E)
    Ex = at(t, E, tj)
    n = min(Ex.shape[1], Er.shape[1])
    T = a15.ratio(w @ Ex[:, :n], w @ Er[:, :n])
    k, cens = a15.k_cross(T, cap)
    ks = [a15.k_cross(a15.ratio(Ex[j, :n], Er[j, :n]), cap)[0]/kc for j in range(len(tj))]
    return dict(k=k, kappa=k/kc, cens=cens, D=a15.dev(T, d_hi), sigma=wsigma(ks, w), kappas=ks, t=tj.tolist())


def pair_sig(run, ref, q, W, cap, d_hi, kc):
    """Variant B: analyze_stage15.pair (k_x, kappa, D of Stage 15) with sigma of pair_tm."""
    p, tm = a15.pair(run, ref, q, W, cap, d_hi, kc), pair_tm(run, ref, q, W, cap, d_hi, kc)
    return dict(p, sigma15=p['sigma'], sigma=tm['sigma'], kappas_tm=tm['kappas'], t=tm['t'])


def mean_check(sp):
    """max relative difference of sum_j w_j E(t_j) from the stored window mean over the runs, windows, QJ and the
    shells where both are positive."""
    worst = 0.0
    for v in sp.values():
        for d in v['win'].values():
            for q in a15.QJ:
                m, s = np.array(d['mean'][q]), np.array(d['w']) @ np.array(d['recs'][q], float)
                ok = (m > 0) & (s > 0)
                worst = max(worst, float(np.max(np.abs(s[ok] - m[ok])/m[ok])))
    return worst


def near_check(runs, ref, n_hi=a15.D_HI):
    """At each record time t_j of ref within DUP of a kept record t_i of a run (t_j != t_i): max |ln E_X(t_j) -
    ln E_X(t_i)| over the shells N0..n_hi and the bound |t_j - t_i| max |d ln E/dt| of the intervals next to t_i."""
    out = []
    for lab, run in sorted(runs.items()):
        for W in a15.WIN:
            if W not in run['win'] or W not in ref['win']:
                continue
            for q in a15.QJ:
                t, E = series(run, W, q)
                L = np.log(E[:, a15.N0:n_hi + 1])
                for x in ref['win'][W]['t']:
                    i = int(np.argmin(np.abs(t - x)))
                    if abs(t[i] - x) >= DUP or t[i] == x:
                        continue
                    d = float(np.max(np.abs(np.log(at(t, E, [x])[0, a15.N0:n_hi + 1]) - L[i])))
                    sl = max(float(np.max(np.abs(L[j + 1] - L[j])))/(t[j + 1] - t[j]) for j in (i - 1, i)
                             if 0 <= j < len(t) - 1)
                    out.append(dict(run=lab, W=W, q=q, t=x, t_near=float(t[i]), dlnE=d, bound=abs(x - t[i])*sl))
    return out


def pairs_of(sp, prob):
    """{label: (64^3 run, 128^3 run)} of the labels on both grids (as analyze_stage15.judge)."""
    out = {}
    for v in sp.values():
        i = v['info']
        if i['problem'] == prob and i['N'] == 64:
            w = a15.find(sp, prob, 128, i['label'])
            if w:
                out[i['label']] = (v, w)
    return out


def show(name, J):
    print('== %s: candidates %s' % (name, ', '.join('%s %.3f' % x for x in sorted(J.get('cand', {}).items()))))
    if 'B' not in J:
        print(J['verdict'])
        return
    print('B* = %s: %s %s' % (J['B'], J['verdict'], J['counts']))
    for k, t in J['table'].items():
        print('  %-20s kES %.3f%s kB %.3f%s d %+.3f delta %.3f  D_ES %.3f D_B %.3f  sigma_ES %.3f sigma_B %.3f' % (
            k, t['kES'], '+' if t['censES'] else ' ', t['kB'], '+' if t['censB'] else ' ', t['d'], t['delta'],
            t['DES'], t['DB'], J['pairs']['ES'][k]['sigma'], J['pairs'][J['B']][k]['sigma']))
    if 'cost' in J:
        print('  cost', json.dumps(J['cost']))


def rejudge(root):
    sp = json.load(open(os.path.join(S15, 'spec.json')))
    j15 = json.load(open(os.path.join(S15, 'judge.json')))
    kc128, kc64 = 128/3.0, 64/3.0
    es256, s8 = sp[a15.ES256], sp[a15.S8]
    refs = {'ES256': es256}
    runs = {v['info']['label']: v for v in sp.values() if v['info']['problem'] == 'turb1' and v['info']['N'] == 128
            and v['info']['role'] == 'b'}
    es128 = runs['ES']
    out = dict(note='post hoc: defined in tasks/todo.md (Stage 16) before this judgement; the pre-fixed verdict is '
                    'that of Stage 15', flags=['H-256 missing: ES-256 only, need = 3 of 4 (as Stage 15)'],
               checks=dict(mean=mean_check(sp), near=near_check(runs, es256)))
    args = (a15.CAP_K, a15.D_HI, kc128)
    out['delta_rerun_tm'] = max(abs(pair_tm(es128, es256, q, 'W1', *args)['kappa']
                                    - pair_tm(s8, es256, q, 'W1', *args)['kappa']) for q in a15.QJ)
    out['delta_rerun_15'] = max(abs(a15.pair(es128, es256, q, 'W1', *args)['kappa']
                                    - a15.pair(s8, es256, q, 'W1', *args)['kappa']) for q in a15.QJ)
    di = {}
    for lab, run in sorted(runs.items()):
        if all(W in run['win'] for W in a15.WIN):
            for q in a15.QJ:
                for W in a15.WIN:
                    a, b = (pair_tm(run, es256, q, W, *args, every_other=e)['kappa'] for e in (False, True))
                    di['%s %s %s' % (lab, q, W)] = abs(a - b)
    out['delta_interp_all'] = di
    out['delta_interp'] = max(di.values())
    for key, fn, extra in (('A', pair_tm, max(out['delta_rerun_tm'], out['delta_interp'])),
                           ('B', pair_sig, out['delta_rerun_15'])):
        J = a15.judge_set(runs, refs, {}, kc128, a15.CAP_K, a15.D_HI, extra, pairfn=fn)
        if 'B' in J:
            J['cost'] = a15.cost_ok(runs['ES']['cost_tstop'], runs[J['B']]['cost_tstop'],
                                    [J['pairs']['ES'][k]['k'] for k in J['table']],
                                    [J['pairs'][J['B']][k]['k'] for k in J['table']])
        out[key] = J
    out['ms05'] = a15.selfconv_set(pairs_of(sp, 'turb05'), kc64, a15.CAP_K64, a15.D_HI64, pairfn=pair_tm)
    out['self_ms1_64'] = a15.selfconv_set(pairs_of(sp, 'turb1'), kc64, a15.CAP_K64, a15.D_HI64, judge=False,
                                          pairfn=pair_tm)
    out['self_ms1_128'] = a15.selfconv_set({'ES': (es128, es256)}, kc128, a15.CAP_K, a15.D_HI, judge=False,
                                           pairfn=pair_tm)
    out['stage15'] = dict(ms1=dict((k, j15['ms1'][k]) for k in ('B', 'verdict', 'counts', 'table', 'cand')),
                          ms05=dict((k, j15['ms05'].get(k)) for k in ('B', 'verdict', 'counts', 'table')),
                          delta_rerun=j15['delta_rerun'], cost=j15.get('cost'))
    os.makedirs(root, exist_ok=True)
    a14.dump(os.path.join(root, 'judge_tm.json'), out)
    print('checks: mean %.2e, near %d records, max dlnE %.2e, max dlnE/bound %.3f' % (
        out['checks']['mean'], len(out['checks']['near']), max(x['dlnE'] for x in out['checks']['near']),
        max(x['dlnE']/x['bound'] for x in out['checks']['near'] if x['bound'] > 0)))
    print('Delta_rerun: time-matched %.4f, Stage 15 %.4f; Delta_interp %.4f (%s)' % (
        out['delta_rerun_tm'], out['delta_rerun_15'], out['delta_interp'], max(di, key=di.get)))
    show('variant A (kappa, D, sigma time-matched), M_s = 1, 128^3', out['A'])
    show('variant B (sigma time-matched), M_s = 1, 128^3', out['B'])
    show('M_s = 0.5, (64, 128) reference, time-matched', out['ms05'])
    print('Stage 15: B* = %s, %s; M_s = 0.5: %s' % (j15['ms1']['B'], j15['ms1']['verdict'], j15['ms05']['verdict']))
    for key in ('self_ms1_64', 'self_ms1_128'):
        for lab, v in out[key]['pairs'].items():
            print('%s %-10s %s' % (key, lab, ' '.join('%s %.3f%s/%.3f' % (k, x['kappa'], '+' if x['cens'] else '',
                                                                         x['D']) for k, x in v.items())))
    fig_kappa_t(root, runs, es256, out['A'])
    return out


def fig_kappa_t(root, runs, ref, J):
    """kappa_j against t of ES and B*: time-matched (kappa_j at the records of ref, window kappa as a segment) and
    Stage 15 (each record of the run against the window mean of ref)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if 'B' not in J:
        return
    fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    for i, q in enumerate(a15.QJ):
        for lab, col in (('ES', a15.COLORS['ES']), (J['B'], a15.COLORS['B'])):
            for j, W in enumerate(a15.WIN):
                p = J['pairs'][lab]['%s %s ES256' % (q, W)]
                p15 = a15.pair(runs[lab], ref, q, W, a15.CAP_K, a15.D_HI, 128/3.0)
                first = j == 0
                ax[i].plot(p['t'], p['kappas'], 'o-', color=col, lw=1.5, ms=4,
                           label='%s, time-matched' % lab if first else None)
                ax[i].plot(runs[lab]['win'][W]['t'], p15['kappas'], 'o--', color=col, lw=1.0, ms=5, mfc='none',
                           alpha=0.7, label='%s, Stage 15 (record / window mean)' % lab if first else None)
                ax[i].hlines(p['kappa'], *a15.WIN[W], color=col, lw=3, alpha=0.4,
                             label='%s, window kappa (time-matched)' % lab if first else None)
        ax[i].set_ylabel('kappa = k_x / k_c')
        ax[i].set_title('%s against ES-256 (128^3, M_s = 1)' % q)
    ax[-1].set_xlabel('t / tau0')
    fig.legend(*ax[0].get_legend_handles_labels(), loc='lower center', ncol=3, fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(root, 'fig_s16_kappa_t.png'), dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------- conv: CV1, CV2 (ps1)
def last_file(case):
    """(path, t) of the out1d file with the largest t over the runs of the case."""
    best = (None, -1.0)
    for d in run_dirs(case):
        for fn in glob.glob(os.path.join(d, 'out1d', 'line_*.dat')):
            with open(fn) as f:
                t = float(re.search(r't =\s*([-+0-9.eE]+)', f.readline()).group(1))
            if t > best[1]:
                best = (fn, t)
    return best


def body(fn):
    with open(fn) as f:
        return [s for s in f if not s.startswith('#')]


def walls(case):
    """Wall times (s) of the runs of the case (/usr/bin/time in calliope.out.std)."""
    w = []
    for d in run_dirs(case):
        fn = os.path.join(d, 'calliope.out.std')
        if os.path.exists(fn):
            w += [float(x) for x in re.findall(r'^wall ([0-9.]+) s$', open(fn, errors='replace').read(), re.M)]
    return w


def conv(root):
    res = {}
    for name in g16.CONV:
        case, old = os.path.join(root, 'conv', name), os.path.join(g16.OLD_CONV, name)
        m, s = a14.lines_1d(case)
        _, so = a14.lines_1d(old)
        d, dt = a14.against(s, so)
        comp = bool(s) and a14.complete_1d(case, m, s)
        (fa, ta), (fb, tb) = last_file(case), last_file(old)
        w = walls(case)
        r = dict(status=a4.status(case), runs=len(run_dirs(case)), lines=(len(s), len(so)), complete=comp,
                 old_max=max(d) if d else None, old_end=d[-1] if d else None, dt_max=dt,
                 last_t=(ta, tb), last_same_values=bool(s and so and np.array_equal(s[-1][1], so[-1][1])),
                 last_same_bytes=bool(fa and fb and a4.same_bytes(fa, fb)),
                 last_same_body=bool(fa and fb and body(fa) == body(fb)),
                 wall=w, core_h=sum(w)*int(open(os.path.join(case, 'nproc')).read())/3600)
        r['ok_old'] = bool(comp and len(s) == len(so) and max(d) <= a14.TOL_OLD1D and dt == 0.0)
        res[name] = r
        print('%-16s %-7s runs %d lines %s complete %-5s max diff %s end %s same last (values %s, data lines %s, '
              'file %s) %.1f core-h' % (name, r['status'], r['runs'], r['lines'], comp, r['old_max'], r['old_end'],
                                        r['last_same_values'], r['last_same_body'], r['last_same_bytes'], r['core_h']))
    ok = dict(CV1=len(res) == len(g16.CONV) and all(r['complete'] for r in res.values()),
              CV2=len(res) == len(g16.CONV) and all(r['ok_old'] for r in res.values()))
    print(' '.join('%s: %s' % (k, 'PASS' if v else 'FAIL') for k, v in ok.items()))
    a3.last_states(os.path.join(root, 'conv'), os.path.join(root, 'conv_last'))
    a14.dump(os.path.join(root, 'compare_conv.json'), dict(cases=res, ok=ok))
    return ok


# ---------------------------------------------------------------- conv table (local)
def rows(prob, Ws, ref):
    """{N: L1 against the FV solution (rho, the 6 variables), I1 against the ideal shock and the state ahead, the rho
    L1 by regions, the rho L1 to 2N (self) and the orders log2(e_N/e_2N) of fv and self}."""
    import stage3_aux as aux
    out = {}
    Ns = sorted(Ws)
    for N in Ns:
        W = Ws[N]
        r = dict(fv=aux.l1(W, ref), fv_all=sum(aux.l1(W, ref, i) for i in range(6)))
        if prob == 'I1':
            e, up = aux.ideal_I1(W)[:2]
            r.update(ideal=e, ahead=[float(x) for x in up])
        a = np.abs(W[0] - ref.reshape(6, N, -1).mean(2)[0])/N
        r['regions'] = {p: float(a[m].sum()) for p, m in aux.regions(prob, N)}
        if 2*N in Ws:
            r['self'] = aux.l1(W, Ws[2*N])
        out[N] = r
    for N in Ns:
        for key in ('fv', 'self'):
            if key in out[N] and 2*N in out and key in out[2*N]:
                out[N]['order_' + key] = math.log2(out[N][key]/out[2*N][key])
    return out


def table(root):
    import stage3_aux as aux
    out = {}
    for prob in ('I1', 'RDpc'):
        ref = aux.last(np.load(os.path.join(aux.FVWIN, aux.FV_CONV[prob]))['W'])
        out[prob + '_fv'] = {f: aux.l1(aux.last(np.load(os.path.join(aux.FVWIN, f))['W']), ref)
                             for f in aux.FV_SELF[prob]}
        for tag, combo in (('', 'e0_n0'), ('_eta2', 'e2_n0')):
            new, old = {}, {}
            for N in N_SUITE:
                fn = os.path.join(S14_LAST[tag], '%s_N%d%s.npz' % (prob, N, tag))
                fo = os.path.join(aux.NUT, '%s_N%d_%s.npz' % (prob, N, combo))
                if os.path.exists(fn):
                    new[N] = np.load(fn)['W']
                if os.path.exists(fo):
                    old[N] = aux.last(np.load(fo)['W'])
            for N in N_CONV:
                new[N] = np.load(os.path.join(root, 'conv_last', '%s_N%d%s.npz' % (prob, N, tag)))['W']
                old[N] = np.load(os.path.join(aux.CONV_LAST, '%s_N%d%s.npz' % (prob, N, tag)))['W']
            key = prob + tag
            out[key] = dict(new=rows(prob, new, ref), old=rows(prob, old, ref),
                            maxdiff={N: a3.maxdiff(new[N], old[N]) for N in sorted(new) if N in old})
            for N, r in sorted(out[key]['new'].items()):
                o = out[key]['old'].get(N, {})
                print('%-10s N=%5d  FV rho %.3e (earlier %s)  order %s  self %s order %s  all %.3e  maxdiff %s' % (
                    key, N, r['fv'], '%.3e' % o['fv'] if o else '-',
                    '%.2f' % r['order_fv'] if 'order_fv' in r else '-',
                    '%.3e' % r['self'] if 'self' in r else '-', '%.2f' % r['order_self'] if 'order_self' in r else '-',
                    r['fv_all'], '%.1e' % out[key]['maxdiff'][N] if N in out[key]['maxdiff'] else '-'))
    a14.dump(os.path.join(root, 'conv_table.json'), out)
    fig_conv(root, out)
    return out


def broken(N, y):
    """N and y with nan where N is not twice the previous one, so that no line crosses a missing N."""
    xs, ys = [], []
    for i, n in enumerate(N):
        if i and n != 2*N[i - 1]:
            xs.append(np.nan)
            ys.append(np.nan)
        xs.append(n)
        ys.append(y[i])
    return xs, ys


def fig_conv(root, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for i, prob in enumerate(('I1', 'RDpc')):
        for tag, col, lab in (('', 'C0', 'C_eta = 0'), ('_eta2', 'C1', 'C_eta = 2')):
            T = out[prob + tag]
            N = sorted(T['new'])
            ax[i].loglog(*broken(N, [T['new'][n]['fv'] for n in N]), 'o-', color=col, lw=1.5, ms=5,
                         label='%s: against FV N = 16384' % lab)
            Ns = [n for n in N if 'self' in T['new'][n]]
            ax[i].loglog(*broken(Ns, [T['new'][n]['self'] for n in Ns]), 's--', color=col, lw=1.0, ms=5,
                         label='%s: N against 2N' % lab)
            No = sorted(T['old'])
            ax[i].loglog(No, [T['old'][n]['fv'] for n in No], 'x', color=col, ms=9, mew=1.0,
                         label='%s: earlier table (prototype, s3g)' % lab)
        ax[i].set_title('%s: rho L1 at t_end' % prob)
        ax[i].set_xlabel('N')
    ax[0].set_ylabel('L1 of rho')
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(root, 'fig_s16_conv.png'), dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------- HF gates (ps1)
def hf_gate(batch, runs):
    """The gate of an HF batch from {name: run_summary} (None while a run is pending or not run): hf64 r_star,
    hf128 r_C."""
    if any(r['kind'] in ('pending', 'not run') for r in runs.values()):
        return None
    rungs = sorted((r['info']['r'], r['complete']) for r in runs.values())
    out = dict(rungs=rungs, kind={k: r['kind'] for k, r in sorted(runs.items())})
    if batch == 'hf64':
        out['r_star'], out['nonmono'] = a15.r_star(rungs)
    elif batch == 'hf128':
        out['r_C'] = min((r for r, ok in rungs if ok), default=None)
    return out


def hf(batch, root):
    out = hf_gate(batch, a15.batch_runs(root, batch))
    if out is None:
        print('%s: a run is pending, no gate' % batch)
        return None
    with open(os.path.join(root, batch + '.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != 'kind'}))
    return out


# ---------------------------------------------------------------- two references: ES-256 and HF-256
def spec_hf(root):
    runs = os.path.dirname(os.path.normpath(root))
    out = {}
    for batch in ('hf128', 'hf256'):
        for c in a15.cases_of(root, batch):
            if a15.run_summary(c)['complete']:
                out[os.path.relpath(c, runs)] = a15.run_spec(c, batch)
    with open(os.path.join(root, 'spec_hf.json'), 'w') as f:
        json.dump(out, f)
    for k, v in out.items():
        print('%-40s %-10s %s' % (k, v['info']['label'], ' '.join(sorted(v['win']))))
    return out


def kref_of(hf, es, pairfn):
    """{(q, W): pair of HF-256 against ES-256} with pairfn (cap CAP_REF; k is k_ref,W)."""
    return {(q, W): pairfn(hf, es, q, W, a15.CAP_REF, a15.D_HI, 1.0) for q in a15.QJ for W in a15.WIN}


def judge2(root):
    """The judgement of tasks/todo.md (Stage 16, section 3), fixed before the HF runs."""
    sp = json.load(open(os.path.join(S15, 'spec.json')))
    sp.update(json.load(open(os.path.join(root, 'spec_hf.json'))))
    with open(os.path.join(root, 'hf128.json')) as f:
        rC = json.load(f)['r_C']
    kc128 = 128/3.0
    args = (a15.CAP_K, a15.D_HI, kc128)
    es256, s8 = sp[a15.ES256], sp[a15.S8]
    hf256 = a15.find(sp, 'turb1', 256, g15.label('HF', rC), roles=('hf256',))
    if hf256 is None or not all(W in hf256['win'] for W in a15.WIN):
        raise SystemExit('HF-256 (HF_r%g) with both windows is not in spec_hf.json' % rC)
    refs = dict(ES256=es256, HF256=hf256)
    runs = {v['info']['label']: v for v in sp.values() if v['info']['problem'] == 'turb1' and v['info']['N'] == 128
            and v['info']['role'] in ('b', 'hf128') and all(W in v['win'] for W in a15.WIN)}
    es128 = runs['ES']
    out = dict(note='post hoc (time matching, as judge_tm.json); 8 pairs against ES-256 and HF-256, need = 6; '
                    'fixed in tasks/todo.md (Stage 16, section 3) before the HF runs', r_C=rC, runs=sorted(runs),
               flags=[])
    out['dt_rec'] = max(min(abs(tj - t) for t in series(run, W, 'E_B')[0]) for run in runs.values()
                        for W in a15.WIN for tj in hf256['win'][W]['t'])
    out['delta_rerun_tm'] = max(abs(pair_tm(es128, es256, q, 'W1', *args)['kappa']
                                    - pair_tm(s8, es256, q, 'W1', *args)['kappa']) for q in a15.QJ)
    di = {'%s %s %s' % (lab, q, W): abs(pair_tm(run, es256, q, W, *args)['kappa']
                                         - pair_tm(run, es256, q, W, *args, every_other=True)['kappa'])
          for lab, run in sorted(runs.items()) for q in a15.QJ for W in a15.WIN}
    out['delta_interp_all'], out['delta_interp'] = di, max(di.values())
    k15 = kref_of(hf256, es256, a15.pair)
    out['delta_rerun_15'] = max(abs(a15.pair(es128, ref, q, 'W1', min(a15.CAP_K, k15[(q, 'W1')]['k']), a15.D_HI,
                                             kc128)['kappa']
                                    - a15.pair(s8, ref, q, 'W1', min(a15.CAP_K, k15[(q, 'W1')]['k']), a15.D_HI,
                                               kc128)['kappa'])
                                for q in a15.QJ for ref in refs.values())
    for key, fn, kfn, extra in (('A', pair_tm, pair_tm, max(out['delta_rerun_tm'], out['delta_interp'])),
                                ('B', pair_sig, a15.pair, out['delta_rerun_15']),
                                ('S15', a15.pair, a15.pair, out['delta_rerun_15'])):
        kp = kref_of(hf256, es256, kfn)
        J = a15.judge_set(runs, refs, {k: v['k'] for k, v in kp.items()}, kc128, a15.CAP_K, a15.D_HI, extra,
                          pairfn=fn)
        J['k_ref'] = {'%s %s' % k: dict(k=v['k'], cens=v['cens']) for k, v in kp.items()}
        out['flags'] += ['%s: k_ref %s = %.2f, the references differ at the large scales' % (key, k, v['k'])
                         for k, v in J['k_ref'].items() if v['k'] < a15.N0 + 1]
        if 'B' in J:
            J['cost'] = dict(ES=runs['ES']['cost_tstop'], B=runs[J['B']]['cost_tstop'],
                             **a15.cost_ok(runs['ES']['cost_tstop'], runs[J['B']]['cost_tstop'],
                                           [J['pairs']['ES'][k]['k'] for k in J['table']],
                                           [J['pairs'][J['B']][k]['k'] for k in J['table']]))
        out[key] = J
    sc = {'ES': (es128, es256)}
    if g15.label('HF', rC) in runs:
        sc[g15.label('HF', rC)] = (runs[g15.label('HF', rC)], hf256)
    out['self_ms1_128'] = a15.selfconv_set(sc, kc128, a15.CAP_K, a15.D_HI, judge=False, pairfn=pair_tm)
    out['cost_ref'] = dict(ES256=es256['cost_tstop'], HF256=hf256['cost_tstop'])
    a14.dump(os.path.join(root, 'judge_hf.json'), out)
    print('runs %s; HF-256 = HF_r%g; dt_rec %.2e tau0' % (' '.join(out['runs']), rC, out['dt_rec']))
    print('Delta_rerun: time-matched %.4f, Stage 15 %.4f; Delta_interp %.4f (%s)' % (
        out['delta_rerun_tm'], out['delta_rerun_15'], out['delta_interp'], max(di, key=di.get)))
    for key, name in (('A', 'variant A (kappa, D, sigma time-matched)'), ('B', 'variant B (sigma time-matched)'),
                      ('S15', 'definitions of Stage 15 (reference only)')):
        print('k_ref %s: %s' % (key, ' '.join('%s %.2f%s' % (k, v['k'], '+' if v['cens'] else '')
                                             for k, v in out[key]['k_ref'].items())))
        show('%s, M_s = 1, 128^3, ES-256 and HF-256' % name, out[key])
    for lab, v in out['self_ms1_128']['pairs'].items():
        print('self (128, 256) %-8s %s' % (lab, ' '.join('%s %.3f%s/%.3f' % (k, x['kappa'], '+' if x['cens'] else '',
                                                                              x['D']) for k, x in v.items())))
    for f in out['flags']:
        print('FLAG', f)
    return out


def judge1(root):
    """The sensitivity of tasks/todo.md (Stage 16, after H-F was stopped): variant A of rejudge with the complete HF
    runs on 128^3 added to the candidates of B*."""
    sp = json.load(open(os.path.join(S15, 'spec.json')))
    hf = json.load(open(os.path.join(root, 'spec_hf.json')))
    kc128 = 128/3.0
    args = (a15.CAP_K, a15.D_HI, kc128)
    es256, s8 = sp[a15.ES256], sp[a15.S8]
    runs = {}
    for v in list(sp.values()) + list(hf.values()):
        i = v['info']
        if i['problem'] == 'turb1' and i['N'] == 128 and i['role'] in ('b', 'hf128'):
            if i['label'] in runs:
                raise SystemExit('two 128^3 runs labelled %s' % i['label'])
            runs[i['label']] = v
    hfs = sorted(lab for lab in runs if runs[lab]['info']['role'] == 'hf128')
    if not hfs:
        raise SystemExit('%s: no complete HF run on 128^3' % os.path.join(root, 'spec_hf.json'))
    es128 = runs['ES']
    out = dict(note='sensitivity of variant A (judge_tm.json), fixed in tasks/todo.md (Stage 16) after H-F was stopped '
                    'and before the HF kappa; the verdict of judge_tm.json is not replaced', hf=hfs,
               flags=['H-256 and HF-256 missing: ES-256 only, need = 3 of 4 (as Stage 15)'],
               t_rec={lab: {W: series(runs[lab], W, 'E_B')[0].tolist() for W in a15.WIN} for lab in hfs})
    out['delta_rerun_tm'] = max(abs(pair_tm(es128, es256, q, 'W1', *args)['kappa']
                                    - pair_tm(s8, es256, q, 'W1', *args)['kappa']) for q in a15.QJ)
    di = {'%s %s %s' % (lab, q, W): abs(pair_tm(run, es256, q, W, *args)['kappa']
                                         - pair_tm(run, es256, q, W, *args, every_other=True)['kappa'])
          for lab, run in sorted(runs.items()) if all(W in run['win'] for W in a15.WIN)
          for q in a15.QJ for W in a15.WIN}
    out['delta_interp_all'], out['delta_interp'] = di, max(di.values())
    J = a15.judge_set(runs, {'ES256': es256}, {}, kc128, a15.CAP_K, a15.D_HI,
                      max(out['delta_rerun_tm'], out['delta_interp']), pairfn=pair_tm)
    if 'B' in J:
        J['cost'] = dict(ES=runs['ES']['cost_tstop'], B=runs[J['B']]['cost_tstop'],
                         **a15.cost_ok(runs['ES']['cost_tstop'], runs[J['B']]['cost_tstop'],
                                       [J['pairs']['ES'][k]['k'] for k in J['table']],
                                       [J['pairs'][J['B']][k]['k'] for k in J['table']]))
    out['A'] = J
    a14.dump(os.path.join(root, 'judge1_hf.json'), out)
    for lab in hfs:
        print('%s records: %s' % (lab, '; '.join('%s %s' % (W, ' '.join('%.4f' % t for t in ts))
                                                  for W, ts in out['t_rec'][lab].items())))
    print('Delta_rerun %.4f; Delta_interp %.4f (%s)' % (out['delta_rerun_tm'], out['delta_interp'],
                                                         max(di, key=di.get)))
    show('variant A with %s, M_s = 1, 128^3' % ', '.join(hfs), J)
    for lab in hfs + ['H_r8']:
        print('%-6s %s' % (lab, ' '.join('%s %.3f%s/%.3f/%.3f' % (k, x['kappa'], '+' if x['cens'] else '', x['D'],
                                                                 x['sigma']) for k, x in J['pairs'][lab].items())))
    return out


# ---------------------------------------------------------------- selftest
def selftest():
    import tempfile
    import stage3_aux as aux
    bad = []

    def check(name, ok):
        print('%-60s %s' % (name, 'PASS' if ok else 'FAIL'))
        if not ok:
            bad.append(name)

    n = np.arange(50, dtype=float)

    def spec(t, b=0.05):
        # ln E linear in t; shell 0 empty
        return np.concatenate(([0.0], n[1:]**-2.0))*np.exp(-(0.2 + b*n)*t)

    def run_of(ts, f, W='W1'):
        ww = a15.window_weights(list(ts), *a15.WIN[W])
        recs = np.array([f(ts[i]) for i, _ in ww])
        w = np.array([x for _, x in ww])
        return {'win': {W: dict(t=[float(ts[i]) for i, _ in ww], w=w.tolist(),
                                mean={q: (w @ recs).tolist() for q in a15.QJ},
                                recs={q: recs.tolist() for q in a15.QJ})}}

    tX = [1.0, 0.75, 0.875, 1.0004, 1.125, 1.25]
    tR = [0.81756, 0.89309, 0.90646, 0.92113, 0.96499, 0.99999, 1.02022, 1.15584]
    X, R = run_of(tX, spec), run_of(tR, lambda t: spec(t, 0.06))
    t, E = series(X, 'W1', 'E_B')
    check('series: sorted, record within DUP left out', list(t) == [0.75, 0.875, 1.0, 1.125, 1.25])
    check('half: both ends', list(half(t, E)[0]) == [0.75, 1.0, 1.25]
          and list(half(t[:4], E[:4])[0]) == [0.75, 1.0, 1.125])
    Ex = at(t, E, tR)
    ex = np.array([spec(x) for x in tR])
    check('at: exact for ln E linear in t', np.max(np.abs(Ex[:, 1:]/ex[:, 1:] - 1)) < 1e-12
          and np.all(np.isnan(Ex[:, 0])))
    check('at: a record at its own time', np.array_equal(at(t, E, [t[2]])[0], E[2]))
    check('at: nearest record within WIN_TOL', np.array_equal(at(t, E, [t[-1] + 0.005])[0], E[-1])
          and np.array_equal(at(t, E, [t[0] - 0.005])[0], E[0]))
    try:
        at(t, E, [t[-1] + 0.02])
        check('at: error beyond WIN_TOL', False)
    except ValueError:
        check('at: error beyond WIN_TOL', True)
    kc = 128/3.0
    p, ph = (pair_tm(X, R, 'E_B', 'W1', a15.CAP_K, a15.D_HI, kc, every_other=e) for e in (False, True))
    w = np.array(R['win']['W1']['w'])
    Tex = a15.ratio(w @ ex, w @ np.array(R['win']['W1']['recs']['E_B']))
    check('pair_tm: window kappa of the exact spectra', abs(p['k'] - a15.k_cross(Tex, a15.CAP_K)[0]) < 1e-9
          and not p['cens'] and 3 < p['k'] < 42)
    check('pair_tm: Delta_interp = 0 for ln E linear in t', abs(p['kappa'] - ph['kappa']) < 1e-10)
    X2 = run_of(tR, spec)
    p2, q2 = pair_tm(X2, R, 'E_B', 'W1', a15.CAP_K, a15.D_HI, kc), a15.pair(X2, R, 'E_B', 'W1', a15.CAP_K, a15.D_HI, kc)
    kj = [a15.k_cross(a15.ratio(spec(x), spec(x, 0.06)), a15.CAP_K)[0]/kc for x in tR]
    check('pair_tm = pair at the same record times (k, D)', p2['k'] == q2['k'] and abs(p2['D'] - q2['D']) < 1e-14)
    check('pair_tm: kappa_j of each record', np.allclose(p2['kappas'], kj, rtol=0, atol=1e-12))
    ps = pair_sig(X, R, 'E_B', 'W1', a15.CAP_K, a15.D_HI, kc)
    check('pair_sig: kappa of pair, sigma of pair_tm', ps['kappa'] == a15.pair(X, R, 'E_B', 'W1', a15.CAP_K,
                                                                             a15.D_HI, kc)['kappa']
          and ps['sigma'] == p['sigma'])
    k = np.array([0.3, 0.5, 0.4, 0.7])
    check('wsigma: ddof = 1 for equal weights', abs(wsigma(k, [2, 2, 2, 2]) - np.std(k, ddof=1)) < 1e-15
          and wsigma([0.4], [1.0]) == 0.0)
    check('wsigma: weights normalized', abs(wsigma(k, [1, 2, 3, 4]) - wsigma(k, [0.1, 0.2, 0.3, 0.4])) < 1e-15)
    check('mean_check', mean_check({'x': X, 'r': R}) < 1e-14)
    nc = near_check({'x': X}, R)
    check('near_check: within the bound', len(nc) > 0 and all(0 < x['dlnE'] <= x['bound']*(1 + 1e-9) for x in nc))
    with tempfile.TemporaryDirectory() as tmp:
        for d, fn, tt in (('chunk_1', 'line_00000.dat', 0.5), ('', 'line_00001.dat', 1.0)):
            os.makedirs(os.path.join(tmp, d, 'out1d'), exist_ok=True)
            with open(os.path.join(tmp, d, 'out1d', fn), 'w') as f:
                f.write('# t =  %.16e  step = 7\n  1.0 2.0\n' % tt)
        f, tt = last_file(tmp)
        check('last_file: largest t over the runs', tt == 1.0 and f.endswith('line_00001.dat')
              and body(f) == ['  1.0 2.0\n'])
    W = np.random.default_rng(1).standard_normal((6, 64))
    check('l1 of the same state', aux.l1(W, W) == 0.0 and aux.l1(W[:, ::2], W) == 0.0)
    runs = {'HF_r%g' % r: dict(info=dict(r=r), complete=ok, kind='complete' if ok else 'invalid')
            for r, ok in ((0.25, True), (0.5, False), (1.0, True), (2.0, True), (4.0, True))}
    g64, g128 = hf_gate('hf64', runs), hf_gate('hf128', runs)
    check('hf_gate: r_star, r_C', g64['r_star'] == 1.0 and g64['nonmono'] and g128['r_C'] == 0.25
          and 'r_star' not in g128)
    runs['HF_r8'] = dict(info=dict(r=8.0), complete=False, kind='pending')
    check('hf_gate: none while pending', hf_gate('hf64', runs) is None)
    runs['HF_r8']['kind'] = 'not run'
    check('hf_gate: none while not run', hf_gate('hf64', runs) is None)
    print('selftest: %s' % ('all PASS' if not bad else 'FAIL: ' + ', '.join(bad)))
    return not bad


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'selftest':
        sys.exit(0 if selftest() else 1)
    if cmd in ('spec_hf', 'judge2', 'judge1'):
        dict(spec_hf=spec_hf, judge2=judge2, judge1=judge1)[cmd](sys.argv[2] if len(sys.argv) > 2 else g16.ROOT)
        sys.exit(0)
    if cmd == 'hf' and len(sys.argv) > 2 and sys.argv[2] in g16.BATCHES:
        hf(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else g16.ROOT)
        sys.exit(0)
    if cmd not in ('conv', 'table', 'rejudge'):
        sys.exit(__doc__)
    dict(conv=conv, table=table, rejudge=rejudge)[cmd](sys.argv[2] if len(sys.argv) > 2 else os.path.join('runs',
                                                                                                         'stage16'))
