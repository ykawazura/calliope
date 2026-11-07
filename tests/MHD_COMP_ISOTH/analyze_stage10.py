#!/usr/bin/env python3
"""Stage 10: the mass budget of the lowest cells of rho in turb1 (M_s = 1) on 128^3, from the dense rerun of the
control of Stage 8 (out3d every 0.01 tau0).

    python3 analyze_stage10.py selftest      synthetic checks of the tools below
    python3 analyze_stage10.py repro ROOT    gate R: rerun against the control      -> ROOT/stage10/repro.json
    python3 analyze_stage10.py run ROOT      budgets of the records of 0.60-1.05 tau0 -> ROOT/stage10/budget.json
    python3 analyze_stage10.py gate ROOT     gates 1a and 1b                        -> ROOT/stage10/gate.json
    python3 analyze_stage10.py judge ROOT    windows, labels, verdict, figures      -> ROOT/stage10/summary.json, fig/*.png

ROOT is the runs directory. judge stops unless repro and gate passed.

m_K(t) is the sum of rho over the K lowest cells. It is the min over the sets of K cells of the sum over the set, so
where the set is unique dm_K/dt = sum over the set of d rho/dt; the set changes where two values are equal, so m_K
is continuous and there is no term from the change of the set. With d rho/dt = -div m + T (filter_lnrho = F and
mass_correction = F, as in Stage 9) and S_T = sum over the set of T,
    I_T = int S_T/m_K dt (trapezoid),  ln m_K(end) - ln m_K(start) = I_adv + I_svv

Criteria (tasks/todo.md, Stage 10, fixed before the job):
    gate R   max |ln rho_min - ln rho_min of the control| <= 0.01 over 0-1.05 tau0 (control interpolated in t)
    gate 1a  |rho_min of the record/rho_min of the row - 1| <= 1e-10 at all records of 0.60-1.05 tau0
    gate 1b  res = |g_c - g_c,f|/A <= 0.05 at the same records, none left out, where g_c is the central difference
             of ln rho_min of the rows, g_c,f the same from analyze_stage9.field_min, A = max(a(x_min), a(p_b), a(p_f)),
             a = (|div m| + |T|)/rho and p_b, p_f the cells of the predicted min of the rows before and after
    windows  (K = 8, 64; t_* = time of the min of m_K over 0.60-1.05 tau0) W_p = [0.60 tau0, t_*], W_s = [0.75 tau0, t_*]
    checks   closed: |Delta ln m_K - I_adv - I_svv| <= 0.5 max(|I_svv|, 0.1)
             stable: |I_svv - I_svv of every second record| <= 0.5 max(|I_svv|, 0.1)
    label    undecided: not closed, not stable, or Delta ln m_K > -0.1
             dig: I_svv <= -0.1 and I_svv/|Delta ln m_K| <= -0.5;  fill: I_svv >= 0.1;  mid otherwise
    verdict  A: the four labels of K = 8, 64 on W_p, W_s are dig; B: they are fill; C otherwise
Reported, not used: K = 1, 512; W_a = [0.75, 1.05] tau0; theta = 0.2, 0.1 and the core with the formulas of Stage 9;
I_svv from records thinned to 0.02, 0.05, 0.1 tau0; sign changes of S_adv/m_K and S_svv/m_K; the terms at x_min.
"""
import contextlib
import io
import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_ic  # noqa: E402
import analyze_stage4 as a4  # noqa: E402
import svv3d  # noqa: E402
import analyze_stage7 as a7  # noqa: E402
import analyze_stage9 as a9  # noqa: E402

TAU0 = a9.TAU0
CONTROL = os.path.join('stage8', 'coef', 'e0.5_C2')
DENSE = os.path.join('stage10', 'dense', 'e0.5_C2')
T0_P, T0_S, T1 = 0.60, 0.75, 1.05  # tau0
T_SLACK = 0.001    # tau0: a record is on the step before the multiple of the interval
NREC = 46          # records of 0.60, 0.61, ..., 1.05 tau0
KS, KS_MAIN = (1, 8, 64, 512), (8, 64)
WINS = ('p', 's')
THIN = (1, 2, 5, 10)  # 0.01, 0.02, 0.05, 0.1 tau0
LEVELS = ('0.2', '0.1', 'core')
REPRO_TOL, REC_TOL = 0.01, 1e-10
GATE1_TOL = a9.GATE1_TOL
CHECK_TOL, I_MIN, DROP_MIN, DIG_FRAC = 0.5, 0.1, 0.1, 0.5
COLOR = {'adv': a9.COLOR['adv'], 'svv': a9.COLOR['svv'], 'net': '#0b0b0b', 'sum': '#898781'}


def write_json(root, name, obj):
    d = os.path.join(root, 'stage10')
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), 'w') as f:
        json.dump(obj, f, indent=1, default=lambda x: x.tolist() if hasattr(x, 'tolist') else str(x))


def read_json(root, name):
    p = os.path.join(root, 'stage10', name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


# ---------------------------------------------------------------------------------------------- gate R

def repro_rows(t, lr, t0, lr0, tmax):
    """Differences of ln rho_min against the control interpolated to t, up to tmax; first index where t differs."""
    m = t <= tmax
    d = lr[m] - np.interp(t[m], t0, lr0)
    i = int(np.argmax(np.abs(d)))
    n = min(len(t), len(t0))
    split = np.nonzero(np.abs(t[:n] - t0[:n]) > 1e-12)[0]
    big = np.nonzero(np.abs(d) > 1e-10)[0]
    return dict(rows=int(m.sum()), max_abs=float(abs(d[i])), t_tau_max=float(t[m][i]/TAU0),
                t_tau_split=None if len(split) == 0 else float(t[split[0]]/TAU0),
                t_tau_1e10=None if len(big) == 0 else float(t[m][big[0]]/TAU0), PASS=bool(abs(d[i]) <= REPRO_TOL))


def repro(root):
    t, lr, gaps, nrow, dm = a9.cons_rows(os.path.join(root, DENSE))
    t0, lr0, _, nrow0, _ = a9.cons_rows(os.path.join(root, CONTROL))
    r = repro_rows(t, lr, t0, lr0, T1*TAU0)
    r.update(rows_all=nrow, rows_control=nrow0, gaps=gaps, max_dmass_corr=dm, t_end_tau=float(t[-1]/TAU0))
    write_json(root, 'repro.json', r)
    print('gate R: %s  max|d ln rho_min| %.3e at %.4f tau0 (limit %.2f), rows %d of %d (control %d), t differs from %s, '
          'd > 1e-10 from %s' % ('PASS' if r['PASS'] else 'FAIL', r['max_abs'], r['t_tau_max'], REPRO_TOL, r['rows'],
                                 nrow, nrow0, r['t_tau_split'], r['t_tau_1e10']))
    print('chain gaps %d, max|dmass_corr| %.1e, last row %.4f tau0' % (gaps, dm, r['t_end_tau']))
    return r


# ------------------------------------------------------------------------------------------ the budget

def lowest(rho, terms, Ks, p0=None):
    """For each K of Ks: m_K, S_T and the largest rho of the K lowest cells; near = fraction of them within 4 cells
    of p0."""
    r = rho.ravel()
    Kmax = max(Ks)
    idx = np.argpartition(r, Kmax - 1)[:Kmax]
    idx = idx[np.argsort(r[idx], kind='stable')]
    rc = np.cumsum(r[idx])
    pos = None if p0 is None else np.array(np.unravel_index(idx, rho.shape)).T
    out = {}
    for K in Ks:
        o = dict(m=float(rc[K - 1]), rho_top=float(r[idx[K - 1]]),
                 S={k: float(T.ravel()[idx[:K]].sum()) for k, T in terms.items()})
        if pos is not None:
            d = np.abs(pos[:K] - np.array(p0))
            d = np.minimum(d, np.array(rho.shape) - d)
            o['near'] = float(np.mean(np.sum(d**2, axis=1) <= 16))
        out[str(K)] = o
    return out


def gate_row(t, lr, t0, s, rate, a):
    """analyze_stage9.fieldmin_row with the scale A = max a over the min of the record and the predicted min of the
    rows before and after (a: field)."""
    i = int(np.argmin(np.abs(t - t0)))
    _, p0 = a9.field_min(s, rate, 0.0)
    cells = [p0, a9.field_min(s, rate, t[i-1] - t[i])[1], a9.field_min(s, rate, t[i+1] - t[i])[1]]
    A = max(float(a[p]) for p in cells)
    return dict(a9.fieldmin_row(t, lr, t0, s, rate, A), A=A, a_min=float(a[p0]))


def terms_of(g, op, sn):
    """The four parts of T and adv = -div m on the grid (formulas of analyze_stage9.record_budget)."""
    T = a9.svv_split(op, sn.s, sn.u, sn.rho, op.eps(sn.rho, sn.u, sn.b))
    T['adv'] = -op.ifft(1j*(g.kxd*sn.mk[0] + g.kyd*sn.mk[1] + g.kzd*sn.mk[2])*op.mL)
    return T


def record_budget(g, op, sn, t, lr, t0):
    """The budget and the row of gate 1b of one record (Snap); the terms as real arrays are not kept."""
    dV = g.Lx*g.Ly*g.Lz/g.ntot
    T = terms_of(g, op, sn)
    svv = sum(T[k] for k in a9.TERMS[1:])
    rho_min = float(sn.rho.min())
    p0 = tuple(int(x) for x in np.unravel_index(int(np.argmin(sn.rho)), sn.rho.shape))
    gate = gate_row(t, lr, t0, sn.s, (svv + T['adv'])/sn.rho, (np.abs(T['adv']) + np.abs(svv))/sn.rho)
    return dict(rho_min=rho_min, p=list(p0), at_min={k: float(T[k][p0]) for k in T}, gate=gate,
                low=lowest(sn.rho, T, KS, p0), dV=dV, sum_abs={k: float(np.abs(T[k]).sum()*dV) for k in T},
                sum_all={k: float(T[k].sum()*dV) for k in T}, **a9.sublevel(sn.rho, T, (0.2, 0.1, 2*rho_min), dV))


def snap_of(g, op, sn, cons, t, lr, tt, path, rec):
    rb = record_budget(g, op, sn, t, lr, tt)
    return dict(t=tt, t_tau=tt/TAU0, path=path, rec=rec, cons=a7.cons_check(cons, tt, rb['rho_min']), **rb)


def records(case):
    """(t, out3d directory, record) of 0.60-1.05 tau0, by time."""
    recs = [r for r in a4.out3d_records(case) if T0_P - T_SLACK <= r[0]/TAU0 <= T1 + T_SLACK]
    return sorted(recs, key=lambda r: r[0])


def run(root):
    case = os.path.join(root, DENSE)
    g, op = svv3d.from_case(case)
    cons = a7.cons_table(case)
    t, lr = a9.cons_rows(case)[:2]
    recs = records(case)
    print('records %d (expected %d), rows %d, all records %d' % (len(recs), NREC, len(t), len(a4.out3d_records(case))),
          flush=True)
    res = dict(case=case, N=g.nlx, eps0=op.eps0, cs2=op.cs2, snaps=[], M0=float(cons[0, 2]),
               drift_M=float(np.max(np.abs(cons[:, 2] - cons[0, 2]))/cons[0, 2]))
    for tt, path, rec in recs:
        x = snap_of(g, op, a7.Snap(g, op, path, rec), cons, t, lr, tt, path, rec)
        res['snaps'].append(x)
        print('t/tau0 %.4f  rho_min %.4f  rel %+.1e  res %.4f' % (tt/TAU0, x['rho_min'], x['cons']['rel'], x['gate']['res']),
              flush=True)
    write_json(root, 'budget.json', res)
    return res


# ------------------------------------------------------------------------------------------ gates 1a, 1b

def gate(root):
    b = read_json(root, 'budget.json')
    if b is None:
        raise SystemExit('gate: no budget.json; the step "run" did not finish, so nothing is decided')
    print('   t/tau0   rho_min      rel        A      a_min    d_b   d_f   g_c      g_c(f)    res  status')
    for s in b['snaps']:
        r = s['gate']
        print('  %7.4f  %.5f  %+.1e  %8.3f %8.3f %5.1f %5.1f %8.3f %8.3f %.4f %s'
              % (s['t_tau'], s['rho_min'], s['cons']['rel'], r['A'], r['a_min'], r['d_b'], r['d_f'], r['g_c'], r['g_c_f'],
                 r['res'], r['status']))
    rel = [abs(s['cons']['rel']) for s in b['snaps']]
    res = [s['gate']['res'] for s in b['snaps']]
    nrec = len(b['snaps'])
    out = dict(n=nrec, nrec_ok=nrec == NREC, max_rel=max(rel), max_res=max(res),
               t_tau_max_res=b['snaps'][int(np.argmax(res))]['t_tau'],
               PASS_1a=max(rel) <= REC_TOL, PASS_1b=max(res) <= GATE1_TOL,
               res_a_min=[abs(s['gate']['g_c'] - s['gate']['g_c_f'])/s['gate']['a_min'] for s in b['snaps']],
               sum_rel={k: max(abs(s['sum_all'][k])/s['sum_abs'][k] for s in b['snaps'] if s['sum_abs'][k] > 0)
                        for k in b['snaps'][0]['sum_all']},
               m_K_min={K: min(s['low'][K]['m'] for s in b['snaps']) for K in b['snaps'][0]['low']},
               drift_M=b['drift_M'])
    out['PASS'] = bool(out['PASS_1a'] and out['PASS_1b'])
    write_json(root, 'gate.json', out)
    print('records %d (expected %d); gate 1a %s: max rel %.2e (limit %.0e); gate 1b %s: max res %.4f at %.4f tau0 (limit %.2f)'
          % (nrec, NREC, 'PASS' if out['PASS_1a'] else 'FAIL', out['max_rel'], REC_TOL,
             'PASS' if out['PASS_1b'] else 'FAIL', out['max_res'], out['t_tau_max_res'], GATE1_TOL))
    print('max|sum over the domain|/sum|.|: %s' % ' '.join('%s %.1e' % kv for kv in out['sum_rel'].items()))
    print('min m_K: %s; drift of the mass %.2e' % (out['m_K_min'], out['drift_M']))
    return out


# -------------------------------------------------------------------------------------- windows and labels

def thin(i0, i1, step):
    """Every step-th index of i0..i1, with i1 always kept."""
    idx = list(range(i0, i1 + 1, step))
    if idx[-1] != i1:
        idx.append(i1)
    return idx


def integ(t, m, Sa, Ss):
    """I_adv, I_svv, Delta ln m_K and the error of the closure over the records given."""
    w = a9.trap_weights(t)
    m = np.asarray(m, dtype=float)
    Ia, Is = float(np.sum(w*np.asarray(Sa)/m)), float(np.sum(w*np.asarray(Ss)/m))
    dl = float(np.log(m[-1]/m[0]))
    return dict(I_adv=Ia, I_svv=Is, dl=dl, eps=abs(dl - Ia - Is))


def k_window(t, m, Sa, Ss, i0, i1):
    """integ on records i0..i1, the checks of the criteria and I_svv of the thinned records."""
    if i1 <= i0:
        return dict(empty=True, i0=i0, i1=i1, t0=float(t[i0]/TAU0), t1=float(t[i1]/TAU0))

    def on(step):
        idx = thin(i0, i1, step)
        return integ(t[idx], m[idx], Sa[idx], Ss[idx])
    w = on(1)
    lim = CHECK_TOL*max(abs(w['I_svv']), I_MIN)
    I2 = on(2)['I_svv']
    w.update(empty=False, i0=i0, i1=i1, t0=float(t[i0]/TAU0), t1=float(t[i1]/TAU0), lim=lim, I_svv2=I2,
             closed=w['eps'] <= lim, stable=abs(w['I_svv'] - I2) <= lim,
             thin={int(st): on(st)['I_svv'] for st in THIN}, thin_eps={int(st): on(st)['eps'] for st in THIN})
    return w


def k_label(w):
    if w['empty'] or not (w['closed'] and w['stable']) or w['dl'] > -DROP_MIN:
        return 'undecided'
    if w['I_svv'] <= -I_MIN and w['I_svv']/abs(w['dl']) <= -DIG_FRAC:
        return 'dig'
    if w['I_svv'] >= I_MIN:
        return 'fill'
    return 'mid'


def overall(lab):
    """The verdict from lab[K][win]."""
    vals = [lab[K][w] for K in KS_MAIN for w in WINS]
    if all(v == 'dig' for v in vals):
        return 'A'
    if all(v == 'fill' for v in vals):
        return 'B'
    return 'C'


def sign_stats(t, x):
    """Sign changes of x over the records and the fraction of the time (trapezoid) with x < 0."""
    w = a9.trap_weights(t)
    x = np.asarray(x, dtype=float)
    return dict(changes=int(np.sum(np.sign(x[1:]) != np.sign(x[:-1]))), frac_neg=float(np.sum(w*(x < 0))/np.sum(w)))


def series(snaps, K):
    k = str(K)
    t = np.array([s['t'] for s in snaps])
    m = np.array([s['low'][k]['m'] for s in snaps])
    Sa = np.array([s['low'][k]['S']['adv'] for s in snaps])
    Ss = np.array([sum(s['low'][k]['S'][x] for x in a9.TERMS[1:]) for s in snaps])
    return t, m, Sa, Ss


def analyse_k(snaps, K):
    t, m, Sa, Ss = series(snaps, K)
    tau = t/TAU0
    n = len(t)
    istar = int(np.argmin(m))
    i_s = int(np.argmin(np.abs(tau - T0_S)))
    wins = dict(p=k_window(t, m, Sa, Ss, 0, istar), s=k_window(t, m, Sa, Ss, i_s, istar),
                a=k_window(t, m, Sa, Ss, i_s, n - 1))
    lab = {w: k_label(wins[w]) for w in WINS}
    out = dict(K=K, t_star_tau=float(tau[istar]), m_star=float(m[istar]), wins=wins, label=lab,
               rho_top=[s['low'][str(K)]['rho_top'] for s in snaps], near=[s['low'][str(K)]['near'] for s in snaps])
    for w, (a, b) in dict(p=(0, istar), s=(i_s, istar), a=(i_s, n - 1)).items():
        if b > a:
            out['sign_' + w] = dict(adv=sign_stats(t[a:b + 1], Sa[a:b + 1]), svv=sign_stats(t[a:b + 1], Ss[a:b + 1]))
    # time series for the figures, per tau0
    out['rate'] = dict(t_tau=tau.tolist(), adv=(TAU0*Sa/m).tolist(), svv=(TAU0*Ss/m).tolist(), ln_m=np.log(m).tolist())
    return out


def theta_tables(snaps):
    """The formulas of Stage 9 on theta = 0.2, 0.1 and the core over 0.60-1.05 and 0.75-1.05 tau0, on all the records
    and on those of 0.1 tau0."""
    t = np.array([s['t'] for s in snaps])
    i_s = int(np.argmin(np.abs(t/TAU0 - T0_S)))
    out = {}
    for name, a in (('0.60-1.05', 0), ('0.75-1.05', i_s)):
        for step in (1, 10):
            idx = thin(a, len(t) - 1, step)
            for j, reg in enumerate(LEVELS):
                S = {k: [snaps[i]['S'][k][j] for i in idx] for k in a9.TERMS}
                M = [snaps[i]['M'][j] for i in idx]
                D = None if reg == 'core' else [snaps[i]['D'][j] for i in idx]
                w = a9.window(t[idx], S, M, D)
                out['%s/%d/%s' % (name, step, reg)] = dict(
                    empty=w['empty'], label=a9.label(w), n_empty=sum(x == 0.0 for x in M),
                    **({} if w['empty'] else dict(Y_adv=w['Y']['adv'], Y_svv=w['Y']['svv'], r=w['r'], E_D=w.get('E_D'))))
    return out


def judge(root):
    rp, gt = read_json(root, 'repro.json'), read_json(root, 'gate.json')
    if rp is None or gt is None or not (rp['PASS'] and gt['PASS']):
        raise SystemExit('judge: gate R or gate 1 did not pass; no labels (criteria of Stage 10)')
    snaps = read_json(root, 'budget.json')['snaps']
    k = {str(K): analyse_k(snaps, K) for K in KS}
    s = dict(k=k, verdict=overall({K: k[str(K)]['label'] for K in KS_MAIN}), theta=theta_tables(snaps),
             at_min=[dict(t_tau=x['t_tau'], rho=x['rho_min'], p=x['p'], adv=x['at_min']['adv'],
                          svv=sum(x['at_min'][y] for y in a9.TERMS[1:])) for x in snaps],
             repro=rp, gate={x: gt[x] for x in gt if x != 'res_a_min'})
    write_json(root, 'summary.json', s)
    print_judge(s)
    figures(root, s)
    return s


def print_judge(s):
    for K in KS:
        o = s['k'][str(K)]
        i = o['wins']['p']['i1']
        print('K = %d: t_* %.4f tau0, m_K/K at t_* %.4f, rho_top at 0.60 tau0 %.4f, at t_* %.4f; near the min: %.2f, %.2f'
              % (K, o['t_star_tau'], o['m_star']/K, o['rho_top'][0], o['rho_top'][i], o['near'][0], o['near'][i]))
        print('  win  t0      t1      dln m    I_adv    I_svv    eps     lim    I_svv2  closed stable  label     '
              'I_svv(0.02, 0.05, 0.1 tau0)')
        for w in ('p', 's', 'a'):
            x = o['wins'][w]
            if x['empty']:
                print('  %s    empty' % w)
                continue
            print('  %s  %.4f  %.4f %+8.3f %+8.3f %+8.3f  %.3f  %.3f  %+7.3f  %-5s  %-5s  %-9s %+.3f %+.3f %+.3f'
                  % (w, x['t0'], x['t1'], x['dl'], x['I_adv'], x['I_svv'], x['eps'], x['lim'], x['I_svv2'], x['closed'],
                     x['stable'], o['label'].get(w, '-'), x['thin'][2], x['thin'][5], x['thin'][10]))
        for w in ('p', 'a'):
            if 'sign_' + w in o:
                print('  signs on %s: adv %s, svv %s' % (w, o['sign_' + w]['adv'], o['sign_' + w]['svv']))
    print('verdict: (%s)' % s['verdict'])
    for key, v in s['theta'].items():
        if v['empty']:
            print('  %-20s empty' % key)
        else:
            print('  %-20s Y_adv %+8.3f Y_svv %+8.3f r %+7.3f E_D %s label %s (empty records %d)'
                  % (key, v['Y_adv'], v['Y_svv'], v['r'], 'none' if v['E_D'] is None else '%.3f' % v['E_D'], v['label'],
                     v['n_empty']))


def figures(root, s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = os.path.join(root, 'stage10', 'fig')
    os.makedirs(d, exist_ok=True)
    # cumulative I_T from 0.60 tau0
    fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.8), sharey=True)
    for ax, K in zip(axs, KS_MAIN):
        r = s['k'][str(K)]['rate']
        tau = np.array(r['t_tau'])
        ca = np.concatenate(([0.0], np.cumsum(0.5*np.diff(tau)*(np.array(r['adv'])[1:] + np.array(r['adv'])[:-1]))))
        cs = np.concatenate(([0.0], np.cumsum(0.5*np.diff(tau)*(np.array(r['svv'])[1:] + np.array(r['svv'])[:-1]))))
        ax.plot(tau, ca, color=COLOR['adv'], lw=2.0, label='I_adv')
        ax.plot(tau, cs, color=COLOR['svv'], lw=2.0, label='I_svv')
        ax.plot(tau, ca + cs, color=COLOR['sum'], lw=2.0, ls='--', label='I_adv + I_svv')
        ax.plot(tau, np.array(r['ln_m']) - r['ln_m'][0], color=COLOR['net'], lw=2.0, label='ln m_K(t)/m_K(0.60 tau0)')
        ax.axvline(s['k'][str(K)]['t_star_tau'], color='#898781', lw=0.8, ls=':')
        ax.axhline(0, color='#898781', lw=0.8)
        ax.set_xlabel('t / tau0')
        ax.set_title('K = %d' % K, fontsize=10)
        ax.grid(alpha=0.35, lw=0.5, color='#e1e0d9')
    axs[0].set_ylabel('e-folds from 0.60 tau0')
    axs[0].legend(frameon=False, fontsize=8)
    fig.savefig(os.path.join(d, 'stage10_cumulative.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    # the rates
    fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.8), sharey=True)
    for ax, K in zip(axs, KS_MAIN):
        r = s['k'][str(K)]['rate']
        for key, name in (('adv', 'tau0 S_adv/m_K'), ('svv', 'tau0 S_svv/m_K')):
            ax.plot(r['t_tau'], r[key], color=COLOR[key], lw=2.0, label=name)
        ax.axhline(0, color='#898781', lw=0.8)
        ax.set_xlabel('t / tau0')
        ax.set_title('K = %d' % K, fontsize=10)
        ax.grid(alpha=0.35, lw=0.5, color='#e1e0d9')
    axs[0].set_ylabel('rate per tau0')
    axs[0].legend(frameon=False, fontsize=8)
    fig.savefig(os.path.join(d, 'stage10_rates.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


# -------------------------------------------------------------------------------------------- tests

def moving_hole(g, t, kap=4.0, c=1.0):
    """rho = 1 - a(t) h(x - c t, y, z), a = 0.5 + 0.3 sin t, h = exp(kap (cos X + cos Y + cos Z - 3)), with its rate."""
    X = 2*np.pi*g.X/g.Lx - c*t
    h = np.exp(kap*(np.cos(X) + np.cos(2*np.pi*g.Y/g.Ly) + np.cos(2*np.pi*g.Z/g.Lz) - 3))
    a, da = 0.5 + 0.3*np.sin(t), 0.3*np.cos(t)
    return 1 - a*h, -da*h - a*c*kap*np.sin(X)*h


class State:
    """The attributes of analyze_stage7.Snap from real m and b."""

    def __init__(self, op, rho, m, b):
        self.rho, self.s = rho, np.log(rho)
        self.mk = [op.fft(x) for x in m]
        self.u = [x/rho for x in m]
        self.b = b


def synthetic_budget(nsub=4):
    """budget.json of a hole on 24^3 with rho stepped by Euler with -div m + T at frozen m and b, nsub steps per
    0.01 tau0: one row per step, so that d rho/dt holds as in the code up to the step."""
    g = gen_ic.Grid(24, 24, 24)
    op = svv3d.SVV(g, 1.0, 0.5, 0.5, True, True)
    kl = 2*np.pi/g.Lx
    rho = 1 - 0.6*np.exp(2*(np.cos(kl*g.X) + np.cos(kl*g.Y) + np.cos(kl*g.Z) - 3))
    m = [0.1*np.sin(kl*g.X), 0.05*np.sin(kl*g.Z), 0.05*np.sin(kl*g.Y)]
    b = [0*rho, 0*rho, 1 + 0*rho]
    t0, h = T0_P*TAU0, 0.01*TAU0/nsub
    nstep = (NREC - 1)*nsub + nsub
    t, rmin, recs = [], [], []
    for i in range(nstep + 1):
        sn = State(op, rho, m, b)
        t.append(t0 + i*h)
        rmin.append(float(rho.min()))
        if i % nsub == 0 and len(recs) < NREC:
            recs.append((t[-1], sn))
        rho = rho + h*sum(terms_of(g, op, sn).values())
    t, lr = np.array(t), np.log(rmin)
    cons = np.zeros((len(t), 17))
    cons[:, 0], cons[:, 11] = t, rmin
    return dict(snaps=[snap_of(g, op, sn, cons, t, lr, tt, 'none', j) for j, (tt, sn) in enumerate(recs)], drift_M=0.0)


def selftest():
    ok = True

    def rep(name, good, msg=''):
        nonlocal ok
        ok &= bool(good)
        print('%-44s %s  %s' % (name, msg, 'PASS' if good else 'FAIL'))

    # 1. lowest against direct sums; K = all cells gives the total
    rng = np.random.default_rng(3)
    r = rng.uniform(0.1, 2.0, (12, 12, 12))
    x = rng.standard_normal((12, 12, 12))
    lo = lowest(r, dict(x=x), (1, 8, 64), (0, 0, 0))
    o = np.argsort(r.ravel())
    rep('lowest = direct sums', all(abs(lo[str(K)]['m'] - r.ravel()[o[:K]].sum()) < 1e-13 and
                                    abs(lo[str(K)]['S']['x'] - x.ravel()[o[:K]].sum()) < 1e-13 and
                                    lo[str(K)]['rho_top'] == r.ravel()[o[K - 1]] for K in (1, 8, 64)))
    pos = np.array(np.unravel_index(o[:8], r.shape)).T
    d = np.minimum(pos, 12 - pos)
    rep('lowest: near', abs(lo['8']['near'] - np.mean(np.sum(d**2, axis=1) <= 16)) < 1e-14)
    y = np.sin(2*np.pi*np.arange(12)/12)[:, None, None]*np.ones((12, 12, 12))
    all_ = lowest(r, dict(y=y), (12**3,))[str(12**3)]
    rep('K = all: total mass, sum of T = 0', abs(all_['m'] - r.sum()) < 1e-11 and abs(all_['S']['y']) < 1e-11,
        'sum %.1e' % abs(all_['S']['y']))
    # 2. the identity on a moving hole: the trapezoid converges to Delta ln m_K
    g = gen_ic.Grid(48, 48, 48)
    t0, t1 = 0.2, 1.0
    for K in (8, 64):
        err = []
        for n in (9, 33, 129, 513):
            ts = np.linspace(t0, t1, n)
            m, S = [], []
            for tt in ts:
                rho, rate = moving_hole(g, tt)
                lo = lowest(rho, dict(r=rate), (K,))[str(K)]
                m.append(lo['m'])
                S.append(lo['S']['r'])
            w = integ(ts, m, S, np.zeros(n))
            err.append(w['eps']/abs(w['dl']))
        rep('identity of m_K, K = %d' % K, err[-1] < 2e-3 and err[-1] < 0.1*err[0],
            'rel err ' + ' '.join('%.1e' % x for x in err))
    # 3. integ by hand: S_adv/m = 2, 2, 0; S_svv/m = -1, 0, -1; weights 0.5, 1.5, 1
    w = integ(np.array([0.0, 1.0, 3.0]), [1.0, 2.0, 4.0], [2.0, 4.0, 0.0], [-1.0, 0.0, -4.0])
    rep('integ by hand', abs(w['I_adv'] - 4) < 1e-14 and abs(w['I_svv'] + 1.5) < 1e-14 and
        abs(w['dl'] - np.log(4)) < 1e-14 and abs(w['eps'] - abs(np.log(4) - 2.5)) < 1e-14)
    rep('thin keeps the end', thin(0, 5, 2) == [0, 2, 4, 5] and thin(0, 4, 2) == [0, 2, 4] and thin(3, 4, 10) == [3, 4])
    # 4. k_window on a set made to close: m = exp(-t), S_adv = 0, S_svv = -m; a stable, closed dig
    ts = np.linspace(0.0, 1.0, 21)
    m = np.exp(-ts)
    w = k_window(ts, m, np.zeros(21), -m, 0, 20)
    rep('k_window: closed, stable', w['closed'] and w['stable'] and abs(w['I_svv'] + 1) < 1e-12 and k_label(w) == 'dig')
    rep('k_window: empty when t_* is before the start', k_window(ts, m, 0*m, -m, 5, 3)['empty'] and
        k_label(k_window(ts, m, 0*m, -m, 5, 3)) == 'undecided')
    # 5. labels and verdict by hand

    def W(dl, Is, closed=True, stable=True):
        return dict(empty=False, dl=dl, I_svv=Is, closed=closed, stable=stable)
    rep('labels', k_label(W(-1, -0.6)) == 'dig' and k_label(W(-1, -0.4)) == 'mid' and k_label(W(-1, 0.3)) == 'fill' and
        k_label(W(-1, -0.05)) == 'mid' and k_label(W(-1, 0.05)) == 'mid' and k_label(W(-0.05, -1)) == 'undecided' and
        k_label(W(-1, -0.6, closed=False)) == 'undecided' and k_label(W(-1, -0.6, stable=False)) == 'undecided' and
        k_label(W(-0.15, -0.1)) == 'dig' and k_label(W(-0.15, -0.09)) == 'mid')

    def L(a, b, c, d):
        return {8: dict(p=a, s=b), 64: dict(p=c, s=d)}
    rep('verdict', overall(L('dig', 'dig', 'dig', 'dig')) == 'A' and overall(L('fill', 'fill', 'fill', 'fill')) == 'B' and
        overall(L('dig', 'dig', 'dig', 'mid')) == 'C' and overall(L('fill', 'dig', 'fill', 'fill')) == 'C' and
        overall(L('undecided', 'dig', 'dig', 'dig')) == 'C')
    rep('sign_stats', sign_stats([0.0, 1.0, 3.0], [1.0, -1.0, -1.0]) == dict(changes=1, frac_neg=2.5/3))
    # 6. gate_row: rows from the predicted min pass, a wrong rate fails, A is the largest a of the three cells
    g = gen_ic.Grid(16, 16, 16)
    s = 0.5*np.cos(2*np.pi*g.X/g.Lx)*np.cos(2*np.pi*g.Y/g.Ly) + 0.2*np.sin(2*np.pi*g.Z/g.Lz)
    rate = 0.3 + 0.1*np.sin(2*np.pi*g.X/g.Lx)  # 0.3 at both minima of s
    tr = np.array([0.99, 1.0, 1.01])
    lr = np.array([a9.field_min(s, rate, x - 1.0)[0] for x in tr])
    a = np.ones_like(s)
    r = gate_row(tr, lr, 1.0, s, rate, a)
    rep('gate_row: consistent rows', r['res'] < 1e-12 and r['status'] == 'ok' and r['A'] == 1.0)
    rep('gate_row: wrong rate fails', gate_row(tr, lr, 1.0, s, 3*rate, a)['status'] == 'FAIL')
    a2 = a.copy()
    a2[a9.field_min(s, rate, -0.01)[1]] = 5.0
    r2 = gate_row(tr, lr + np.array([0.0, 0.0, 0.01]), 1.0, s, rate, a2)
    rep('gate_row: scale of the predicted cell', r2['A'] == 5.0 and r2['a_min'] <= 5.0)
    # 7. repro_rows: the same rows, a shifted time base, a change over the limit
    t = np.linspace(0.0, 4.0, 41)
    lr = np.sin(t)
    a = repro_rows(t, lr, t, lr, 3.5)
    b = repro_rows(t, lr + 0.02*(t > 2), t, lr, 3.5)
    rep('repro_rows', a['PASS'] and a['max_abs'] == 0 and a['t_tau_split'] is None and not b['PASS'] and
        abs(b['t_tau_1e10'] - t[t > 2][0]/TAU0) < 1e-14)
    # 8. the records are read from the time.dat of the real layout, the gate stops without a budget
    with tempfile.TemporaryDirectory() as case:
        tq = TAU0*np.array([T0_P - 2*T_SLACK, T0_P - 1e-4, 0.8, T1 - 1e-4, T1 + 2*T_SLACK])
        for d, x in ((os.path.join(case, 'chunk_1'), tq[:3]), (case, tq[3:])):
            os.makedirs(os.path.join(d, 'out3d'))
            np.savetxt(os.path.join(d, 'out3d', 'time.dat'), np.column_stack((x, 0*x)))
        rc = records(case)
        try:
            gate(case)
            nobudget = False
        except SystemExit:
            nobudget = True
    rep('records: window, time order, chunks', [r[0] for r in rc] == sorted(tq[1:4].tolist()) and
        rc[0][1].endswith(os.path.join('chunk_1', 'out3d')) and rc[0][2] == 1 and rc[-1][2] == 0, '')
    rep('gate stops without the budget', nobudget, '')
    # 9. the whole chain on a state stepped with -div m + T: the table of the rows, gates, windows, figures, the stop
    b = synthetic_budget()
    with tempfile.TemporaryDirectory() as root, contextlib.redirect_stdout(io.StringIO()):
        write_json(root, 'budget.json', b)
        write_json(root, 'repro.json', dict(PASS=True))
        gt = gate(root)
        s = judge(root)
        figs = sorted(os.listdir(os.path.join(root, 'stage10', 'fig')))
        back = read_json(root, 'summary.json')
        write_json(root, 'repro.json', dict(PASS=False))
        try:
            judge(root)
            stopped = False
        except SystemExit:
            stopped = True
    rep('chain: gates pass', gt['PASS'] and gt['n'] == NREC and gt['nrec_ok'] and gt['max_rel'] == 0,
        'max res %.4f' % gt['max_res'])
    rep('chain: the sums vanish', max(gt['sum_rel'].values()) <= 1e-12, 'max %.1e' % max(gt['sum_rel'].values()))
    pw = [s['k'][str(K)]['wins']['p'] for K in KS_MAIN]
    rep('chain: closed on W_p, K = 8, 64', all(x['eps'] <= 0.02*abs(x['dl']) and x['closed'] and x['stable'] for x in pw),
        'eps/|dln m| ' + ' '.join('%.1e' % (x['eps']/abs(x['dl'])) for x in pw))
    a = s['k']['8']
    tau = np.array(a['rate']['t_tau'])
    rep('chain: W_p ends at the min of m_K', abs(a['wins']['p']['t1'] - a['t_star_tau']) < 1e-14 and
        a['t_star_tau'] == tau[np.argmin(a['rate']['ln_m'])] and abs(a['wins']['p']['t0'] - T0_P) < 1e-3, '')
    rep('chain: table of theta, core present', s['theta']['0.60-1.05/1/core']['empty'] is False and
        s['theta']['0.60-1.05/1/0.1']['empty'] is True, '')
    rep('chain: json, figures, stop without gates', back['verdict'] == s['verdict'] and len(figs) == 2 and stopped, s['verdict'])
    print('selftest:', 'all PASS' if ok else 'FAIL')
    return ok


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'selftest':
        sys.exit(0 if selftest() else 1)
    {'repro': repro, 'run': run, 'gate': gate, 'judge': judge}[cmd](sys.argv[2])
