#!/usr/bin/env python3
"""Stage 11: the dependence on svv_eps0 of the mass budget of the 64 lowest cells of rho in turb1 (M_s = 1, 128^3),
from the existing out3d records (every 0.1 tau0) of the cases e0.25_C2, e0.5_C2 (control) and e1.0_C2 of Stage 8.

    python3 analyze_stage11.py selftest     synthetic checks and the error estimate from the dense rerun of Stage 10
    python3 analyze_stage11.py run ROOT     budgets of the records of 0.60-1.00 tau0 -> ROOT/stage11/budget_<eps0>.json
    python3 analyze_stage11.py gate ROOT    gates N, C, 1a and 1b                    -> ROOT/stage11/gate.json
    python3 analyze_stage11.py judge ROOT   F, u, label, reported values, figure     -> ROOT/stage11/summary.json, fig/*.png

ROOT is the runs directory (it holds stage8/coef and stage10/budget.json). judge stops unless the gates passed.

The budget is that of analyze_stage10: m_K = sum of rho over the K lowest cells, I_T = int S_T/m_K dt (trapezoid) and
Delta ln m_K = I_adv + I_svv. At fixed fields T_svv is linear in svv_eps0 (svv3d.SVV.eps, analyze_stage9.svv_split), so
a departure of I_svv from a value proportional to eps0 comes from the fields.

Criteria (tasks/todo.md, Stage 11, fixed before the values of K = 64 of e1.0 and e0.25 were computed):
    gate N   5 records of 0.60-1.00 tau0 in each case and eps0 of calliope.in 0.25, 0.5, 1.0
    gate C   |I_svv(K = 64) of the 5 records of the control - the same of the 5 records of the dense rerun at the same
             times| <= 0.05 (a record of the dense rerun is at the same time if within T_SLACK)
    gate 1a  |rho_min of the record/rho_min of the row - 1| <= 1e-10 at the 15 records
    gate 1b  res of analyze_stage10.gate_row <= 0.05 at the 15 records, none left out
    test     (K = 64) eps_e = |Delta ln m_64 - I_adv - I_e|, delta_e = max(0.187, eps_e),
             F = (I_1.0 - I_0.25)/(1.5 I_0.5),  u = sqrt(delta_1.0^2 + delta_0.25^2)/|1.5 I_0.5| + |F| delta_0.5/|I_0.5|;
             [F - u, F + u] holds 1 and not 0: prop; 0 and not 1: const; neither and 0 < F < 1: mid;
             neither otherwise: out; both: undecided
0.187 is the largest |error| of I_svv(K = 64) of 5 records 0.1 tau0 apart against the dense integral over the same
window, over the 6 starts 0.60, 0.61, ..., 0.65 tau0 of the dense rerun (selftest recomputes it).
Reported, not used: K = 8 and 512 (delta 1.480, 0.156); f_hi, f_lo; I_e/eps0; dl max(|u| + c_f); the rates of the
records; the differences of I_adv and Delta ln m_64 from the control; R_64; ||T_svv||_1/eps0; R of the rows.
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
import analyze_stage8 as a8  # noqa: E402
import analyze_stage9 as a9  # noqa: E402
import analyze_stage10 as a10  # noqa: E402
from analyze_stage1 import read_nml  # noqa: E402

TAU0 = a9.TAU0
RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'runs')
CASES = {'0.25': 'e0.25_C2', '0.5': 'e0.5_C2', '1.0': 'e1.0_C2'}  # control: 0.5
EPS = {e: float(e) for e in CASES}
T0, T1 = 0.60, 1.00  # tau0
T_SLACK = a10.T_SLACK
NREC = 5
KS, K_MAIN = (8, 64, 512), 64
DELTA = {8: 1.480, 64: 0.187, 512: 0.156}
DENSE_STEP, PHASES = 10, 6  # records of the dense rerun per 0.1 tau0; starts 0.60, ..., 0.65 tau0
DENSE_I64 = -0.796          # I_svv(K = 64) of the dense records of 0.60, 0.70, ..., 1.00 tau0 (plan)
I_ADV_ERR = 0.405           # the same estimate as DELTA for I_adv, K = 64
GATE_C_TOL = 0.05
REC_TOL, GATE1_TOL = a10.REC_TOL, a10.GATE1_TOL
R_STAGE8 = {'0.25': 0.0500, '1.0': -0.0864}  # R of 0.6-1.0 tau0 (REPORT_stage8.md)
LABELS = ('prop', 'const', 'mid', 'out', 'undecided')
COLOR = {'0.25': '#eb6834', '0.5': '#0b0b0b', '1.0': '#2a78d6', 'ref': '#898781'}  # the colors of Stage 8


def case_dir(root, e):
    return os.path.join(root, 'stage8', 'coef', CASES[e])


def write_json(root, name, obj):
    d = os.path.join(root, 'stage11')
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), 'w') as f:
        json.dump(obj, f, indent=1, default=lambda x: x.tolist() if hasattr(x, 'tolist') else str(x))


def read_json(root, name):
    p = os.path.join(root, 'stage11', name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def tmean(t, x):
    w = a9.trap_weights(t)
    return float(np.sum(w*np.asarray(x, dtype=float))/np.sum(w))


# ------------------------------------------------------------------------------------------ the budget

def record_budget(g, op, sn, t, lr, t0):
    """The budget of the K lowest cells, the row of gate 1b, the domain sums and dl_d max(|u_d| + c_f,d) of one record."""
    dV = g.Lx*g.Ly*g.Lz/g.ntot
    T = a10.terms_of(g, op, sn)
    svv = sum(T[k] for k in a9.TERMS[1:])
    gate = a10.gate_row(t, lr, t0, sn.s, (svv + T['adv'])/sn.rho, (np.abs(T['adv']) + np.abs(svv))/sn.rho)
    both = dict(adv=T['adv'], svv=svv)
    return dict(rho_min=float(sn.rho.min()), gate=gate, low=a10.lowest(sn.rho, T, KS), dV=dV,
                sum_abs={k: float(np.abs(v).sum()*dV) for k, v in both.items()},
                sum_all={k: float(v.sum()*dV) for k, v in both.items()},
                c_dl=[x/op.eps0 for x in op.eps(sn.rho, sn.u, sn.b)])


def records(case):
    """(t, out3d directory, record) of 0.60-1.00 tau0, by time."""
    recs = [r for r in a4.out3d_records(case) if T0 - T_SLACK <= r[0]/TAU0 <= T1 + T_SLACK]
    return sorted(recs, key=lambda r: r[0])


def run(root):
    for e in CASES:
        case = case_dir(root, e)
        g, op = svv3d.from_case(case)
        cons = a7.cons_table(case)
        t, lr = a9.cons_rows(case)[:2]
        recs = records(case)
        print('%s: eps0 %g, records %d (expected %d), rows %d' % (CASES[e], op.eps0, len(recs), NREC, len(t)), flush=True)
        res = dict(case=case, N=g.nlx, eps0=op.eps0, cs2=op.cs2, nml=read_nml(os.path.join(case, 'calliope.in')),
                   snaps=[], M0=float(cons[0, 2]), drift_M=float(np.max(np.abs(cons[:, 2] - cons[0, 2]))/cons[0, 2]),
                   L_rows=a8.mean_log(t/TAU0, np.exp(lr), T0, T1))
        for tt, path, rec in recs:
            rb = record_budget(g, op, a7.Snap(g, op, path, rec), t, lr, tt)
            res['snaps'].append(dict(t=tt, t_tau=tt/TAU0, path=path, rec=rec,
                                     cons=a7.cons_check(cons, tt, rb['rho_min']), **rb))
            print('  t/tau0 %.4f  rho_min %.4f  rel %+.1e  res %.4f' % (tt/TAU0, rb['rho_min'],
                                                                       res['snaps'][-1]['cons']['rel'], rb['gate']['res']),
                  flush=True)
        write_json(root, 'budget_%s.json' % e, res)


# ------------------------------------------------------------------------------------------------ gates

def load(root):
    b = {e: read_json(root, 'budget_%s.json' % e) for e in CASES}
    if any(v is None for v in b.values()):
        raise SystemExit('gate: a budget_<eps0>.json is missing; the step "run" did not finish, so nothing is decided')
    return b


def same_time(ts, td):
    """Index of the record of td nearest to each time of ts and the largest distance in tau0."""
    idx = [int(np.argmin(np.abs(td - x))) for x in ts]
    return idx, float(np.max(np.abs(td[idx] - ts))/TAU0)


def gate(root):
    b = load(root)
    dense = a10.read_json(root, 'budget.json')
    if dense is None:
        raise SystemExit('gate: no stage10/budget.json; gate C cannot be evaluated, so nothing is decided')
    n = {e: len(b[e]['snaps']) for e in CASES}
    pass_n = all(n[e] == NREC and abs(b[e]['eps0'] - EPS[e]) <= 1e-12 for e in CASES)
    # gate C
    c = b['0.5']['snaps']
    idx, dt = same_time(np.array([s['t'] for s in c]), np.array([s['t'] for s in dense['snaps']]))
    I_c = a10.integ(*a10.series(c, K_MAIN))['I_svv']
    I_d = a10.integ(*a10.series([dense['snaps'][i] for i in idx], K_MAIN))['I_svv']
    pass_c = bool(dt <= T_SLACK and abs(I_c - I_d) <= GATE_C_TOL)
    print('gate N %s: records %s (expected %d), eps0 %s' % ('PASS' if pass_n else 'FAIL', n, NREC,
                                                           {e: b[e]['eps0'] for e in CASES}))
    print('gate C %s: I_svv(K = 64) control %+.4f, dense %+.4f (plan %+.3f), |diff| %.4f (limit %.2f), '
          'max |t - t_dense| %.1e tau0 (limit %.0e)' % ('PASS' if pass_c else 'FAIL', I_c, I_d, DENSE_I64,
                                                       abs(I_c - I_d), GATE_C_TOL, dt, T_SLACK))
    print('   case  t/tau0   rho_min      rel        A      a_min    g_c      g_c(f)    res  res(a_min) status')
    rel, res, rows = [], [], {}
    for e in CASES:
        rows[e] = []
        for s in b[e]['snaps']:
            r = s['gate']
            old = abs(r['g_c'] - r['g_c_f'])/r['a_min']
            rows[e].append(dict(t_tau=s['t_tau'], rel=s['cons']['rel'], res=r['res'], res_a_min=old, A=r['A'],
                                a_min=r['a_min']))
            rel.append(abs(s['cons']['rel']))
            res.append(r['res'])
            print('  %5s  %.4f  %.5f  %+.1e  %8.3f %8.3f %8.3f %8.3f  %.4f  %.4f  %s'
                  % (e, s['t_tau'], s['rho_min'], s['cons']['rel'], r['A'], r['a_min'], r['g_c'], r['g_c_f'], r['res'],
                     old, r['status']))
    out = dict(n=n, PASS_N=pass_n, I_c=I_c, I_d=I_d, dt_tau=dt, dense_idx=idx, PASS_C=pass_c, max_rel=max(rel),
               max_res=max(res), PASS_1a=max(rel) <= REC_TOL, PASS_1b=max(res) <= GATE1_TOL, rows=rows,
               sum_rel={e: {k: max(abs(s['sum_all'][k])/s['sum_abs'][k] for s in b[e]['snaps']) for k in ('adv', 'svv')}
                        for e in CASES},
               m_K_min={e: {K: min(s['low'][K]['m'] for s in b[e]['snaps']) for K in b[e]['snaps'][0]['low']}
                        for e in CASES},
               drift_M={e: b[e]['drift_M'] for e in CASES},
               nml_diff=sorted(k for k in set(b['0.5']['nml']) | set(b['0.25']['nml']) | set(b['1.0']['nml'])
                               if len({b[e]['nml'].get(k) for e in CASES}) > 1))
    out['PASS'] = bool(out['PASS_N'] and out['PASS_C'] and out['PASS_1a'] and out['PASS_1b'])
    write_json(root, 'gate.json', out)
    print('gate 1a %s: max rel %.2e (limit %.0e); gate 1b %s: max res %.4f (limit %.2f)'
          % ('PASS' if out['PASS_1a'] else 'FAIL', out['max_rel'], REC_TOL, 'PASS' if out['PASS_1b'] else 'FAIL',
             out['max_res'], GATE1_TOL))
    for e in CASES:
        print('%s: max|sum over the domain|/sum|.| adv %.1e svv %.1e; min m_K %s; drift of the mass %.2e'
              % (CASES[e], out['sum_rel'][e]['adv'], out['sum_rel'][e]['svv'], out['m_K_min'][e], out['drift_M'][e]))
    print('entries of calliope.in that differ: %s' % out['nml_diff'])
    print('gates: %s' % ('PASS' if out['PASS'] else 'FAIL'))
    return out


# ------------------------------------------------------------------------------------------------- test

def label(F, u):
    has1, has0 = F - u <= 1 <= F + u, F - u <= 0 <= F + u
    if has1 and has0:
        return 'undecided'
    if has1:
        return 'prop'
    if has0:
        return 'const'
    return 'mid' if 0 < F < 1 else 'out'


def test(I, eps, delta):
    """F, u and the label from I[e] = I_svv and eps[e] = the closure error of the three cases."""
    d = {e: max(delta, eps[e]) for e in CASES}
    F = (I['1.0'] - I['0.25'])/(1.5*I['0.5'])
    u = float(np.hypot(d['1.0'], d['0.25'])/abs(1.5*I['0.5']) + abs(F)*d['0.5']/abs(I['0.5']))
    return dict(F=float(F), u=u, lo=float(F - u), hi=float(F + u), delta=d, label=label(F, u))


def one_sided(I, d):
    """f_hi = (I_1.0 - I_0.5)/I_0.5, f_lo = (I_0.25 - I_0.5)/(-0.5 I_0.5): 1 if proportional, 0 if constant; the
    widths are the linear bounds with the deltas d."""
    a = abs(I['0.5'])
    f_hi = (I['1.0'] - I['0.5'])/I['0.5']
    f_lo = (I['0.25'] - I['0.5'])/(-0.5*I['0.5'])
    return dict(f_hi=float(f_hi), u_hi=float((d['1.0'] + abs(f_hi + 1)*d['0.5'])/a),
                f_lo=float(f_lo), u_lo=float((2*d['0.25'] + abs(2 - f_lo)*d['0.5'])/a))


def judge(root):
    gt = read_json(root, 'gate.json')
    if gt is None or not gt['PASS']:
        raise SystemExit('judge: the gates did not pass; no label (criteria of Stage 11)')
    b = load(root)
    k = {}
    for K in KS:
        w = {e: a10.integ(*a10.series(b[e]['snaps'], K)) for e in CASES}
        I = {e: w[e]['I_svv'] for e in CASES}
        tst = test(I, {e: w[e]['eps'] for e in CASES}, DELTA[K])
        k[str(K)] = dict(w=w, test=tst, one=one_sided(I, tst['delta']), I_over_eps={e: I[e]/EPS[e] for e in CASES},
                         D={e: {x: w[e][x] - w['0.5'][x] for x in ('I_svv', 'I_adv', 'dl')} for e in CASES})
    rates, R64, c_dl, l1 = {}, {}, {}, {}
    for e in CASES:
        sn = b[e]['snaps']
        t, m, Sa, Ss = a10.series(sn, K_MAIN)
        rates[e] = dict(t_tau=(t/TAU0).tolist(), adv=(TAU0*Sa/m).tolist(), svv=(TAU0*Ss/m).tolist(), m=m.tolist())
        c_dl[e] = [tmean(t, [s['c_dl'][d] for s in sn]) for d in range(3)]
        l1[e] = tmean(t, [s['sum_abs']['svv'] for s in sn])/EPS[e]
    lm = {e: (np.array(rates[e]['t_tau']), np.log(rates[e]['m'])) for e in CASES}
    for e in CASES:
        R64[e] = dict(mean=tmean(*lm[e]) - tmean(*lm['0.5']), first=float(lm[e][1][0] - lm['0.5'][1][0]),
                      R_rows=b[e]['L_rows'] - b['0.5']['L_rows'], R_stage8=R_STAGE8.get(e))
    s = dict(k=k, label=k[str(K_MAIN)]['test']['label'], rates=rates, R64=R64, c_dl=c_dl, l1_over_eps=l1,
             gate={x: gt[x] for x in gt if x != 'rows'})
    write_json(root, 'summary.json', s)
    print_judge(s)
    figure(root, s)
    return s


def print_judge(s):
    for K in KS:
        o = s['k'][str(K)]
        print('K = %d%s: case  dln m    I_adv    I_svv    eps    delta  I_svv/eps0  D_svv    D_adv    D_dl'
              % (K, '' if K == K_MAIN else ' (reported)'))
        for e in CASES:
            w = o['w'][e]
            print('          %5s %+8.4f %+8.4f %+8.4f  %.4f  %.3f  %+8.4f  %+8.4f %+8.4f %+8.4f'
                  % (e, w['dl'], w['I_adv'], w['I_svv'], w['eps'], o['test']['delta'][e], o['I_over_eps'][e],
                     o['D'][e]['I_svv'], o['D'][e]['I_adv'], o['D'][e]['dl']))
        t, f = o['test'], o['one']
        print('  F %+.4f  u %.4f  [%+.4f, %+.4f]  label %s;  f_hi %+.4f (u %.4f), f_lo %+.4f (u %.4f)'
              % (t['F'], t['u'], t['lo'], t['hi'], t['label'], f['f_hi'], f['u_hi'], f['f_lo'], f['u_lo']))
    print('label (K = %d): %s' % (K_MAIN, s['label']))
    for e in CASES:
        r = s['R64'][e]
        print('%5s: R_64 %+.4f, at 0.60 tau0 %+.4f; R of the rows %+.4f (Stage 8 %s); dl max(|u|+c_f) %s; '
              '||T_svv||_1/eps0 %.4e' % (e, r['mean'], r['first'], r['R_rows'], r['R_stage8'],
                                        ' '.join('%.5f' % x for x in s['c_dl'][e]), s['l1_over_eps'][e]))
    print('rates of K = %d per tau0 (S_adv/m, S_svv/m):' % K_MAIN)
    for e in CASES:
        r = s['rates'][e]
        print('  %5s ' % e + '  '.join('%.4f: %+.3f %+.3f' % x for x in zip(r['t_tau'], r['adv'], r['svv'])))


def figure(root, s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = os.path.join(root, 'stage11', 'fig')
    os.makedirs(d, exist_ok=True)
    o = s['k'][str(K_MAIN)]
    fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.8))
    ax = axs[0]
    x = np.linspace(0.0, 1.05, 2)
    I5 = o['w']['0.5']['I_svv']
    ax.plot(x, I5*x/0.5, color=COLOR['ref'], lw=1.4, ls='--')
    ax.plot(x, I5 + 0*x, color=COLOR['ref'], lw=1.4, ls=':')
    ax.annotate('proportional', (1.0, 2*I5), xytext=(-6, 0), textcoords='offset points', ha='right', va='center',
                fontsize=8, color='#55534e')
    ax.annotate('constant', (1.0, I5), xytext=(0, 6), textcoords='offset points', ha='right', fontsize=8,
                color='#55534e')
    for e in CASES:  # the abscissa names the case, so this panel has no legend; error bars: delta
        ax.errorbar(EPS[e], o['w'][e]['I_svv'], yerr=o['test']['delta'][e], color=COLOR[e], marker='o', ms=8, lw=2.0,
                    capsize=4, mec='white', mew=2)
    ax.axhline(0, color='#898781', lw=0.8)
    ax.set_xlabel('svv_eps0')
    ax.set_ylabel('I_svv (K = %d, 0.60-1.00 tau0)' % K_MAIN)
    ax.set_title('F = %.2f, u = %.2f' % (o['test']['F'], o['test']['u']), fontsize=10)
    ax.grid(alpha=0.35, lw=0.5, color='#e1e0d9')
    ax = axs[1]
    for e in CASES:
        r = s['rates'][e]
        ax.plot(r['t_tau'], r['svv'], color=COLOR[e], lw=2.0, marker='o', ms=6, mec='white', mew=1.5,
                label='eps0 = %s' % e)
    ax.axhline(0, color='#898781', lw=0.8)
    ax.set_xlabel('t / tau0')
    ax.set_ylabel('tau0 S_svv/m_%d' % K_MAIN)
    ax.grid(alpha=0.35, lw=0.5, color='#e1e0d9')
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(os.path.join(d, 'stage11_eps0.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


# ------------------------------------------------------------------------------------------------- tests

def dense_errors(dense, K):
    """For the starts 0.60, ..., 0.65 tau0: I of the 5 records 0.1 tau0 apart minus I of all the records between."""
    out = []
    for p in range(PHASES):
        i5 = a10.integ(*a10.series([dense['snaps'][p + DENSE_STEP*j] for j in range(NREC)], K))
        ia = a10.integ(*a10.series(dense['snaps'][p:p + DENSE_STEP*(NREC - 1) + 1], K))
        out.append(dict(svv=i5['I_svv'] - ia['I_svv'], adv=i5['I_adv'] - ia['I_adv'], I5=i5['I_svv'], eps5=i5['eps']))
    return out


def synthetic_snaps(ts, e, a=0.3, c=1.2, shift=0.0):
    """Records with constant rates S_adv/m = a and S_svv/m = -c eps0 + shift, so that the trapezoid is exact."""
    out = []
    for i, t in enumerate(ts):
        r = a - c*EPS[e] + shift
        m = np.exp(r*(t - ts[0]))
        low = {str(K): dict(m=K/64*m, rho_top=0.1, S=dict(adv=a*K/64*m, Ls=(r - a)*K/64*m, Lu=0.0, Hs=0.0, Hu=0.0))
               for K in (1,) + KS}
        out.append(dict(t=float(t), t_tau=float(t/TAU0), rec=i, cons=dict(rel=0.0), rho_min=0.1, low=low,
                        gate=dict(res=0.01, A=2.0, a_min=1.0, g_c=0.0, g_c_f=0.02, status='ok'),
                        sum_abs=dict(adv=1.0, svv=EPS[e]), sum_all=dict(adv=0.0, svv=0.0), c_dl=[0.02, 0.02, 0.02]))
    return out


def synthetic_root(root, nrec=None, shift=0.0):
    ts = TAU0*np.array([0.5999, 0.6998, 0.7999, 0.8999, 0.9999])
    for e in CASES:
        sn = synthetic_snaps(ts if nrec is None or e != '1.0' else ts[:nrec], e)
        write_json(root, 'budget_%s.json' % e, dict(eps0=EPS[e], snaps=sn, drift_M=0.0, L_rows=-1.0 - EPS[e],
                                                     nml=dict(svv_eps0=e, nlx='128')))
    td = TAU0*(0.5999 + 0.01*np.arange(46))
    a10.write_json(root, 'budget.json', dict(snaps=synthetic_snaps(td, '0.5', shift=shift)))


def selftest():
    ok = True

    def rep(name, cond, info=''):
        nonlocal ok
        ok &= bool(cond)
        print('  %-58s %s  %s' % (name, 'PASS' if cond else 'FAIL', info))
    # 1. delta and the dense value of the plan from the dense rerun of Stage 10
    dense = a10.read_json(RUNS, 'budget.json')
    rep('stage10/budget.json present', dense is not None and len(dense['snaps']) == 46, RUNS)
    if dense is not None:
        for K in KS:
            er = dense_errors(dense, K)
            mx = max(abs(x['svv']) for x in er)
            rep('delta of K = %d' % K, round(mx, 3) == DELTA[K], 'max|err| %.4f, errors %s'
                % (mx, ' '.join('%+.3f' % x['svv'] for x in er)))
        er = dense_errors(dense, K_MAIN)
        rep('dense I_svv(K = 64) at 0.60 tau0 = plan', round(er[0]['I5'], 3) == DENSE_I64, '%.4f' % er[0]['I5'])
        rep('error of I_adv (K = 64)', round(max(abs(x['adv']) for x in er), 3) == I_ADV_ERR,
            '%.4f' % max(abs(x['adv']) for x in er))
    # 2. the labels by hand and the widths of the plan (I_0.5 = -0.8)

    def I(a, b, c):
        return {'0.25': a, '0.5': b, '1.0': c}
    z = I(0.0, 0.0, 0.0)
    p, cst = test(I(-0.4, -0.8, -1.6), z, 0.187), test(I(-0.8, -0.8, -0.8), z, 0.187)
    rep('u of the plan, delta 0.187: F = 0 -> 0.22, F = 1 -> 0.45', abs(p['F'] - 1) < 1e-14 and cst['F'] == 0 and
        round(cst['u'], 2) == 0.22 and round(p['u'], 2) == 0.45, '%.4f %.4f' % (cst['u'], p['u']))
    p1, c1 = test(I(-0.4, -0.8, -1.6), I(0.0, 0.0, 0.563), 0.187), test(I(-0.8, -0.8, -0.8), I(0.0, 0.0, 0.563), 0.187)
    rep('u of the plan, eps 0.563 in one case: 0.49, 0.73', round(c1['u'], 2) == 0.49 and round(p1['u'], 2) == 0.73,
        '%.4f %.4f' % (c1['u'], p1['u']))
    big = I(0.563, 0.563, 0.563)
    cases = [(p, 'prop'), (cst, 'const'), (test(I(-0.6, -0.8, -1.2), z, 0.01), 'mid'),
             (test(I(-0.2, -0.8, -3.0), z, 0.01), 'out'), (test(I(-1.6, -0.8, -0.4), z, 0.01), 'out'),
             (test(I(-0.6, -0.8, -1.2), big, 0.187), 'undecided')]
    rep('labels of the five kinds', all(x['label'] == y for x, y in cases) and {y for _, y in cases} == set(LABELS),
        ' '.join('%s(F %+.2f u %.2f)' % (x['label'], x['F'], x['u']) for x, _ in cases))
    rep('label at the ends of the interval', label(0.5, 0.5) == 'undecided' and label(1.5, 0.5) == 'prop' and
        label(-0.5, 0.5) == 'const' and label(1.0 + 1e-9, 1e-10) == 'out')
    o = one_sided(I(-0.4, -0.8, -1.6), dict(p['delta']))
    oc = one_sided(I(-0.8, -0.8, -0.8), dict(p['delta']))
    rep('f_hi, f_lo: 1 if proportional, 0 if constant', abs(o['f_hi'] - 1) < 1e-14 and abs(o['f_lo'] - 1) < 1e-14 and
        oc['f_hi'] == 0 and oc['f_lo'] == 0)
    # 3. T_svv is linear in eps0 at fixed fields; record_budget is consistent on a smooth state
    g = gen_ic.Grid(16, 16, 16)
    kl = 2*np.pi/g.Lx
    rho = 1 - 0.6*np.exp(2*(np.cos(kl*g.X) + np.cos(kl*g.Y) + np.cos(kl*g.Z) - 3))
    m = [0.1*np.sin(kl*g.X), 0.05*np.sin(kl*g.Z), 0.05*np.sin(kl*g.Y)]
    b = [0.1*np.cos(kl*g.Y), 0*rho, 1 + 0*rho]
    T = {}
    for e0 in (0.5, 1.0):
        op = svv3d.SVV(g, 1.0, e0, 0.5, True, True)
        T[e0] = a10.terms_of(g, op, a10.State(op, rho, m, b))
    s5, s1 = (sum(T[x][k] for k in a9.TERMS[1:]) for x in (0.5, 1.0))
    rep('T_svv doubles with eps0, adv unchanged', np.abs(s1 - 2*s5).max() <= 1e-13*np.abs(s5).max() and
        np.array_equal(T[0.5]['adv'], T[1.0]['adv']), 'max|T_svv| %.2e' % np.abs(s5).max())
    op = svv3d.SVV(g, 1.0, 0.5, 0.5, True, True)
    sn = a10.State(op, rho, m, b)
    rate = (s5 + T[0.5]['adv'])/rho
    tr = np.array([0.99, 1.0, 1.01])
    lr = np.array([a9.field_min(sn.s, rate, x - 1.0)[0] for x in tr])
    rb = record_budget(g, op, sn, tr, lr, 1.0)
    rep('record_budget: rows from the state pass gate 1b', rb['gate']['res'] < 1e-12, 'res %.1e' % rb['gate']['res'])
    rep('record_budget: domain sums vanish, c_dl = eps/eps0',
        max(abs(rb['sum_all'][k])/rb['sum_abs'][k] for k in ('adv', 'svv')) < 1e-12 and
        np.allclose(rb['c_dl'], [x/0.5 for x in op.eps(rho, sn.u, b)], rtol=1e-15, atol=0),
        'sums %.1e %.1e' % tuple(abs(rb['sum_all'][k])/rb['sum_abs'][k] for k in ('adv', 'svv')))
    rep('record_budget: m_K of the K lowest cells', all(
        abs(rb['low'][str(K)]['m'] - np.sort(rho.ravel())[:K].sum()) < 1e-12 for K in KS))
    # 4. the records of the window from the time.dat of the real layout
    with tempfile.TemporaryDirectory() as case:
        tq = TAU0*np.array([T0 - 2*T_SLACK, T0 - 1e-4, 0.7, 0.8, 0.9, T1 - 1e-4, T1 + 2*T_SLACK])
        os.makedirs(os.path.join(case, 'out3d'))
        np.savetxt(os.path.join(case, 'out3d', 'time.dat'), np.column_stack((tq, 0*tq)))
        rc = records(case)
    rep('records: the 5 of 0.60-1.00 tau0', [r[0] for r in rc] == tq[1:6].tolist() and rc[0][2] == 1)
    # 5. the chain on synthetic budgets: proportional rates, gate C, the stops
    with tempfile.TemporaryDirectory() as root, contextlib.redirect_stdout(io.StringIO()):
        synthetic_root(root)
        gt = gate(root)
        sm = judge(root)
        figs = os.listdir(os.path.join(root, 'stage11', 'fig'))
        synthetic_root(root, shift=0.2)
        gtc = gate(root)
        synthetic_root(root, nrec=4)
        gtn = gate(root)
        try:
            judge(root)
            stopped = False
        except SystemExit:
            stopped = True
    o = sm['k'][str(K_MAIN)]
    rep('chain: gates pass, closure exact', gt['PASS'] and gt['I_c'] == gt['I_d'] and
        max(o['w'][e]['eps'] for e in CASES) < 1e-14, 'dt %.1e' % gt['dt_tau'])
    rep('chain: F = 1, label prop, I/eps0 equal', abs(o['test']['F'] - 1) < 1e-12 and sm['label'] == 'prop' and
        np.ptp(list(o['I_over_eps'].values())) < 1e-12, 'F %.6f u %.3f' % (o['test']['F'], o['test']['u']))
    # ln m_64 = (0.3 - 1.2 eps0)(t - t_0): the difference of the time means is -0.6 (t_end - t_0)/2 = -0.12 pi for 1.0
    rep('chain: R_64, R of the rows, figure', abs(sm['R64']['1.0']['mean'] + 0.12*np.pi) < 1e-12 and
        sm['R64']['1.0']['first'] == 0 and abs(sm['R64']['1.0']['R_rows'] + 0.5) < 1e-14 and
        figs == ['stage11_eps0.png'], 'R_64 %+.5f' % sm['R64']['1.0']['mean'])
    rep('chain: gate C fails on a shifted dense rerun', not gtc['PASS_C'] and not gtc['PASS'],
        '|diff| %.3f' % abs(gtc['I_c'] - gtc['I_d']))
    rep('chain: gate N fails on 4 records, judge stops', not gtn['PASS_N'] and stopped)
    print('selftest:', 'all PASS' if ok else 'FAIL')
    return ok


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'selftest':
        sys.exit(0 if selftest() else 1)
    {'run': run, 'gate': gate, 'judge': judge}[cmd](sys.argv[2])
