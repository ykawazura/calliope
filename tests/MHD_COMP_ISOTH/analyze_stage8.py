#!/usr/bin/env python3
"""Stage 8: the dependence of the holes of rho in turb1 (M_s = 1, 128^3) on svv_eps0 and C_shock, from the runs
of gen_stage8.py (stage8/coef) and the run their control repeats (gen_stage8.REF).

    python3 analyze_stage8.py rhomin ROOT   rho_min(t), L, m and the checks of the runs  -> ROOT/stage8/rhomin.json
    python3 analyze_stage8.py tail ROOT     analyze_stage7.tail_case at T_OUT             -> ROOT/stage8/tail_<case>.json
    python3 analyze_stage8.py hole ROOT     analyze_stage7.hole_case at T_WIN             -> ROOT/stage8/hole_<case>.json
    python3 analyze_stage8.py judge ROOT    criteria, figures         -> ROOT/stage8/summary.json, ROOT/stage8/fig/*.png

ROOT is the runs directory. judge needs rhomin and tail; without hole the auxiliary hole quantities are left out.

Criteria (tasks/todo.md, Stage 8, fixed before the runs), for each case c of the coef batch:
    L_c  time mean over WIN of ln rho_min(t) (conservation.dat of every step, trapezoids in t; the same end for
         every run, analyze_stage8.windows)
    m_c  min of rho_min(t) over 0 <= t <= T_STOP tau0
    H_c  sum over the out3d records at T_WIN of the number of local minima with rho < H_THR
    V_c  mean over the same records of F_V(rho < V_THR)
    R_c  L_c - L_control
    noise: |L_control - L_ref| <= NOISE, otherwise both axes are 'undecidable'
    axis (doubled, halved coefficient): deeper if R_doubled <= -DELTA; fills if R_doubled >= DELTA and
    R_halved <= -DELTA; insensitive if both |R| < DELTA; non-monotone otherwise
A case that did not reach T_STOP tau0 (analyze_stage8.reached) is compared with the control over WIN_SHORT if it reached its end, and left
out otherwise; an axis with one case left is judged from that case alone (partial; fixed with the code, before
the runs): doubled R <= -DELTA or halved R >= DELTA deeper, doubled R >= DELTA or halved R <= -DELTA fills,
|R| < DELTA insensitive.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_stage7 as a7  # noqa: E402
from analyze_stage6 import dt_history, dt_stats  # noqa: E402
from gen_stage8 import COEFS, REF, T_STOP  # noqa: E402

TAU0 = np.pi  # pi/M_s (case.json)
CASES = {name: os.path.join('stage8', 'coef', name) for name, _, _ in COEFS}
COEF = {name: dict(svv_eps0=e, C_shock=C) for name, e, C in COEFS}
CONTROL = 'e0.5_C2'
AXES = {'svv_eps0': dict(doubled='e1.0_C2', halved='e0.25_C2'), 'C_shock': dict(doubled='e0.5_C4', halved='e0.5_C1')}
WIN, WIN_SHORT = (0.6, T_STOP), (0.6, 1.0)  # tau0
DELTA = 0.1
NOISE = DELTA/10
T_OUT = tuple(round(0.5 + 0.1*j, 10) for j in range(8))  # 0.5, ..., 1.2 tau0
T_WIN = T_OUT[1:]
H_THR, V_THR = '0.1', '0.2'  # keys of nmin (analyze_stage7.NHOLE_THR) and F_V (analyze_stage7.THETA)
COLOR = {'e1.0_C2': '#2a78d6', 'e0.25_C2': '#eb6834', 'e0.5_C4': '#1baf7a', 'e0.5_C1': '#eda100',
         CONTROL: '#0b0b0b', 'ref': '#898781'}  # categorical slots 1-4 (dataviz), primary and muted ink


def write_json(root, name, obj):
    d = os.path.join(root, 'stage8')
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), 'w') as f:
        json.dump(obj, f, indent=1, default=lambda x: x.tolist() if hasattr(x, 'tolist') else str(x))


def read_json(root, name):
    p = os.path.join(root, 'stage8', name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def runs(root):
    """{label: case directory} of the coef cases and of REF ('ref')."""
    out = {k: os.path.join(root, v) for k, v in CASES.items()}
    out['ref'] = os.path.join(root, REF)
    return out


# ------------------------------------------------------------------------------------------- rho_min

def reached(t_final, dt_last, t):
    """True if a run that ended at t_final with the last step dt_last reached t: the rule at which the code stops
    at t_end (time.t_end_reached), so the last t of a run that stopped there may be up to dt/2 short of t_end."""
    return t_final >= t - 0.5*dt_last


def mean_log(t, y, t0, t1):
    """Time mean of ln y over t0 <= t <= t1 (trapezoids, the ends interpolated)."""
    inside = (t > t0) & (t < t1)
    tt = np.concatenate(([t0], t[inside], [t1]))
    ly = np.interp(tt, t, np.log(y))
    return float(np.sum(0.5*(ly[1:] + ly[:-1])*np.diff(tt))/(t1 - t0))


def windows(res):
    """{key of w: dict(t1, L={run: L or None})} for w = WIN, WIN_SHORT and the runs of res (run_stats). A run covers
    w = (t0, t1) if its rows start by t0 and it reached t1 (reached). The means of all the runs covering w
    end at the same time: t1, or the earliest last row (up to T_STOP tau0) of these runs if it precedes t1 (by less
    than a step of that run)."""
    out = {}
    for w in (WIN, WIN_SHORT):
        cov = [k for k, r in res.items() if r['series']['t_tau'][0] <= w[0] and reached(r['t_final'], r['dt_last'], w[1])]
        t1 = min([w[1]] + [res[k]['series']['t_tau'][-1] for k in cov])
        out['%g-%g' % w] = dict(t1=t1, L={k: mean_log(np.array(res[k]['series']['t_tau']),
                                                      np.array(res[k]['series']['rho_min']), w[0], t1)
                                          if k in cov else None for k in res})
    return out


def run_stats(case):
    """exit_status, the final t (t + dt of the last row of cfl.dat), m, the mass and the steps of the run, and
    rho_min(t) up to T_STOP tau0."""
    cons = a7.cons_table(case)
    t, i = np.unique(cons[:, 0], return_index=True)
    c = cons[i]
    h = dt_history(case)
    st, _ = dt_stats(h, TAU0)
    last = h['cfls'][-1][-1]
    t_final = float(last[0] + last[1])
    sel = t <= T_STOP*TAU0*(1 + 1e-12)
    k = int(np.argmin(c[sel, 11]))
    status = open(os.path.join(case, 'exit_status')).read().strip() if os.path.exists(
        os.path.join(case, 'exit_status')) else None
    return dict(case=case, exit_status=status, t_final=t_final/TAU0, dt_last=float(last[1])/TAU0,
                reached=reached(t_final/TAU0, float(last[1])/TAU0, T_STOP),
                m=float(c[sel, 11][k]), t_m=float(t[sel][k]/TAU0),
                mass=dict(rel_max=float(np.abs(c[:, 3]).max()), corr_tot=float(c[-1, 17]),
                          M_last_over_first=float(c[-1, 2]/c[0, 2] - 1)),
                finite=bool(np.all(np.isfinite(c))) and st['finite'], steps=st['steps'],
                changes=st['changes'], max_dt_over_cfl=st['max_dt_over_cfl'], dt_min=st['dt_min'],
                dt_max=st['dt_max'], series=dict(t_tau=(t[sel]/TAU0).tolist(), rho_min=c[sel, 11].tolist()))


def rhomin(root):
    res = {name: run_stats(case) for name, case in runs(root).items()}
    for wk, v in windows(res).items():
        print('L over %s tau0: end %.6f tau0' % (wk, v['t1']))
        for k, r in res.items():
            r.setdefault('L', {})[wk] = v['L'][k]
            r.setdefault('L_end', {})[wk] = v['t1']
    for name, r in res.items():
        print('%-9s %-14s t_final/tau0 %.6f  L(0.6-1.2) %s  L(0.6-1.0) %s  m %.4e at %.3f  steps %d (%d changes)  '
              'max dt/dt_cfl %.4f  mass rel %.1e corr %.1e  finite %s'
              % (name, r['exit_status'], r['t_final'], fmt(r['L']['%g-%g' % WIN]), fmt(r['L']['%g-%g' % WIN_SHORT]),
                 r['m'], r['t_m'], r['steps'], r['changes'], r['max_dt_over_cfl'], r['mass']['rel_max'],
                 r['mass']['corr_tot'], r['finite']))
    write_json(root, 'rhomin.json', res)
    return res


def fmt(x, f='%+.4f'):
    return '-' if x is None else f % x


# ----------------------------------------------------------------------------------------- tail, hole

def tail(root):
    for name, case in runs(root).items():
        if name != 'ref':
            print('== %s' % name)
            write_json(root, 'tail_%s.json' % name, a7.tail_case(case, T_OUT))


def hole(root):
    for name, case in runs(root).items():
        if name != 'ref':
            print('== %s' % name)
            write_json(root, 'hole_%s.json' % name, a7.hole_case(case, T_WIN))


# --------------------------------------------------------------------------------------------- judge

def axis_verdict(Rd, Rh):
    """Verdict of an axis from R of the doubled and the halved coefficient (None: case left out)."""
    if Rd is None and Rh is None:
        return dict(verdict='undecidable', partial=True)
    if Rd is None or Rh is None:
        # one case left: the sign of R as the effect of more dissipation
        R = Rd if Rh is None else -Rh
        v = 'deeper' if R <= -DELTA else 'fills' if R >= DELTA else 'insensitive'
        return dict(verdict=v, partial=True)
    if Rd <= -DELTA:
        v = 'deeper'
    elif Rd >= DELTA and Rh <= -DELTA:
        v = 'fills'
    elif abs(Rd) < DELTA and abs(Rh) < DELTA:
        v = 'insensitive'
    else:
        v = 'non-monotone'
    return dict(verdict=v, partial=False)


def next_step(svv, C):
    if svv == 'deeper' and C in ('fills', 'insensitive'):
        return 'N: revise the minimum principle of ES-SVV (ask), 3b with e1.0'
    if svv == 'fills' and C == 'fills':
        return 'P: 3b with e1.0, then 4 (DNS reference, ask)'
    if svv == 'insensitive':
        return 'dx: skip 3b, test the SVV mass term at the minima on the snapshots'
    return 'consult'


def hole_aux(hd):
    """Over the holes at T_WIN: the fraction with SVV < 0 at the point, the median of SVV/|-div m| and the fraction
    of suspect holes (analyze_stage7.classify)."""
    hs = [h for s in hd['snaps'] if s['t_req'] in T_WIN for h in s['holes']]
    svv = np.array([h['budget']['svv'] for h in hs])
    adv = np.array([h['budget']['adv'] for h in hs])
    return dict(n=len(hs), frac_svv_neg=float(np.mean(svv < 0)), med_svv_over_adv=float(np.median(svv/np.abs(adv))),
                frac_suspect=float(np.mean([h['suspect'] for h in hs])))


def judge(root):
    rm = read_json(root, 'rhomin.json')
    tails = {k: read_json(root, 'tail_%s.json' % k) for k in CASES}
    holes = {k: read_json(root, 'hole_%s.json' % k) for k in CASES}
    wl, ws = '%g-%g' % WIN, '%g-%g' % WIN_SHORT
    ctrl = rm[CONTROL]
    # the window of each case against the control
    win = {}
    for k in CASES:
        if rm[k]['reached'] and ctrl['reached']:
            win[k] = wl
        elif rm[k]['L'][ws] is not None and ctrl['L'][ws] is not None:
            win[k] = ws
        else:
            win[k] = None
    nw = win[CONTROL]
    noise = dict(window=nw, L_control=ctrl['L'][nw] if nw else None, L_ref=rm['ref']['L'][nw] if nw else None)
    noise['dL'] = noise['L_control'] - noise['L_ref'] if nw else None
    noise['PASS'] = nw is not None and abs(noise['dL']) <= NOISE
    # H and V on the records at T_WIN present in every case
    common = [t for t in T_WIN if all(any(s['t_req'] == t for s in tails[k]['snaps']) for k in CASES)]
    rows = {}
    for k in CASES:
        sn = [s for s in tails[k]['snaps'] if s['t_req'] in common]
        w = win[k]
        rows[k] = dict(COEF[k], window=w, L=rm[k]['L'][w] if w else None, m=rm[k]['m'], t_m=rm[k]['t_m'],
                       H=int(sum(s['nmin'][H_THR] for s in sn)), V=float(np.mean([s['F_V'][V_THR] for s in sn])),
                       reached=rm[k]['reached'], t_final=rm[k]['t_final'], exit_status=rm[k]['exit_status'])
    r0 = rows[CONTROL]
    for k, r in rows.items():
        Lc = ctrl['L'][r['window']] if r['window'] else None
        r['R'] = r['L'] - Lc if r['window'] else None
        r['ln_m'] = float(np.log(r['m']/r0['m']))
        r['H_ratio'] = r['H']/r0['H'] if r0['H'] else None
        r['V_ratio'] = r['V']/r0['V'] if r0['V'] else None
        # direction of each quantity: -1 deeper holes, +1 shallower, 0 within DELTA (R, ln m) or 10% (H, V)
        r['dir'] = dict(R=0 if r['R'] is None or abs(r['R']) < DELTA else int(np.sign(r['R'])),
                        ln_m=0 if abs(r['ln_m']) < DELTA else int(np.sign(r['ln_m'])),
                        H=0 if r['H_ratio'] is None or abs(r['H_ratio'] - 1) < 0.1 else -int(np.sign(r['H_ratio'] - 1)),
                        V=0 if r['V_ratio'] is None or abs(r['V_ratio'] - 1) < 0.1 else -int(np.sign(r['V_ratio'] - 1)))
        if holes[k] is not None:
            r['hole'] = hole_aux(holes[k])
    axes = {}
    for a, d in AXES.items():
        if not noise['PASS']:
            axes[a] = dict(verdict='undecidable', partial=False, reason='noise')
        else:
            axes[a] = axis_verdict(rows[d['doubled']]['R'], rows[d['halved']]['R'])
        axes[a].update(R_doubled=rows[d['doubled']]['R'], R_halved=rows[d['halved']]['R'])
    nxt = next_step(axes['svv_eps0']['verdict'], axes['C_shock']['verdict']) if noise['PASS'] else 'consult'
    summ = dict(criteria=dict(WIN=WIN, WIN_SHORT=WIN_SHORT, L_end=ctrl['L_end'], DELTA=DELTA, NOISE=NOISE,
                              H_THR=H_THR, V_THR=V_THR, T_WIN=T_WIN, common=common),
                noise=noise, rows=rows, axes=axes, next=nxt,
                checks={k: dict(exit_status=v['exit_status'], t_final=v['t_final'], finite=v['finite'],
                                max_dt_over_cfl=v['max_dt_over_cfl'], mass=v['mass'], steps=v['steps'])
                        for k, v in rm.items()})
    write_json(root, 'summary.json', summ)
    print_judge(summ)
    figures(root, rm, tails)
    return summ


def print_judge(s):
    n = s['noise']
    print('noise (control - ref over %s): %s  limit %.3f  %s'
          % (n['window'], fmt(n['dL'], '%+.2e'), NOISE, 'PASS' if n['PASS'] else 'FAIL'))
    print('%-9s %8s %7s %-8s %9s %8s %10s %8s %5s %8s %8s %6s  dir(R,m,H,V)  hole(n, SVV<0, med SVV/|adv|, suspect)'
          % ('case', 'svv_eps0', 'C_shock', 'window', 'L', 'R', 'm', 'ln(m/m0)', 'H', 'H/H0', 'V', 'V/V0'))
    for k, r in s['rows'].items():
        hh = r.get('hole')
        print('%-9s %8.2f %7.1f %-8s %9s %8s %10.4e %+8.3f %5d %8s %8.2e %6s  %-13s %s'
              % (k, r['svv_eps0'], r['C_shock'], r['window'] or '-', fmt(r['L']), fmt(r['R']), r['m'], r['ln_m'],
                 r['H'], fmt(r['H_ratio'], '%.3f'), r['V'], fmt(r['V_ratio'], '%.3f'),
                 ','.join('%+d' % r['dir'][q] for q in ('R', 'ln_m', 'H', 'V')),
                 '-' if hh is None else '%d, %.2f, %+.3f, %.2f' % (hh['n'], hh['frac_svv_neg'], hh['med_svv_over_adv'],
                                                                   hh['frac_suspect'])))
    for a, v in s['axes'].items():
        print('axis %-8s: %s%s  (R doubled %s, halved %s)' % (a, v['verdict'], ' (partial)' if v['partial'] else '',
                                                             fmt(v['R_doubled']), fmt(v['R_halved'])))
    print('next:', s['next'])


def figures(root, rm, tails):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = os.path.join(root, 'stage8', 'fig')
    os.makedirs(d, exist_ok=True)
    label = {k: 'svv_eps0 %.2f, C_shock %g%s' % (COEF[k]['svv_eps0'], COEF[k]['C_shock'],
                                                  ' (control)' if k == CONTROL else '') for k in CASES}
    label['ref'] = 'Stage 6 run (svv_eps0 0.50, C_shock 2)'
    order = [AXES['svv_eps0']['doubled'], AXES['svv_eps0']['halved'], AXES['C_shock']['doubled'],
             AXES['C_shock']['halved'], CONTROL, 'ref']
    # rho_min(t)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.axvspan(WIN[0], WIN[1], color='#e1e0d9', alpha=0.5, lw=0)
    for k in order:
        s = rm[k]['series']
        ax.plot(s['t_tau'], s['rho_min'], color=COLOR[k], lw=2.0 if k != 'ref' else 1.4,
                ls='--' if k == 'ref' else '-', label=label[k], zorder=3 if k == CONTROL else 2)
    ax.set_yscale('log')
    ax.set_xlim(0, T_STOP)
    ax.set_xlabel('t / tau0')
    ax.set_ylabel('rho_min')
    ax.set_title('turb1, 128^3: rho_min(t) (shaded: window of L)', fontsize=10)
    ax.grid(alpha=0.35, lw=0.5, color='#e1e0d9')
    ax.legend(frameon=False, fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=2)
    fig.savefig(os.path.join(d, 'stage8_rhomin.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    # number of local minima below H_THR at T_OUT
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for k in order[:-1]:
        sn = sorted(tails[k]['snaps'], key=lambda s: s['t_req'])
        ax.plot([s['t_tau'] for s in sn], [s['nmin'][H_THR] for s in sn], color=COLOR[k], lw=2.0, marker='o', ms=5,
                label=label[k], zorder=3 if k == CONTROL else 2)
    ax.axvspan(WIN[0], WIN[1], color='#e1e0d9', alpha=0.5, lw=0)
    ax.set_xlabel('t / tau0')
    ax.set_ylabel('local minima with rho < %s' % H_THR)
    ax.set_title('turb1, 128^3: number of holes (shaded: records of H)', fontsize=10)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.35, lw=0.5, color='#e1e0d9')
    ax.legend(frameon=False, fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=2)
    fig.savefig(os.path.join(d, 'stage8_holes.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    a = sys.argv[1:]
    cmd = {'rhomin': rhomin, 'tail': tail, 'hole': hole, 'judge': judge}
    if len(a) != 2 or a[0] not in cmd:
        sys.exit(__doc__)
    cmd[a[0]](a[1])
