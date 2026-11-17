"""Markdown tables and figures of the Stage 2 runs (python3 stage2.py table).

Reads runs/stage2/<case>.json and .npz of stage2.py; writes runs/stage2/tables.md and fig/stage2_*.png
next to this file. The values of a case are those of its last attempt (the dt/2, dt/4 reruns of a failed run).
"""
import os
import json
import glob
import math
import numpy as np
import riemann_exact as rx
from gen_stage1 import RIEMANN
from stage2 import SHOCK_VARIANTS, NS
from stage3_table import COLOR, MARK, GRAY, e, stats, status, gmax, order, table

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig')


def load(out):
    R = {}
    for f in glob.glob(os.path.join(out, '*.json')):
        with open(f) as fh:
            d = json.load(fh)
        R[d['case']['name']] = d
    return R


def shock(d):
    """Worst values over the shocks and star states of a Riemann run (None if it did not reach t_end)."""
    if d is None or not d['ok']:
        return None
    sh, st = d.get('shocks', []), d.get('stars', [])
    return dict(speed=gmax([abs(s['speed_rel_err']) for s in sh if 'speed_rel_err' in s]),
                star=gmax([abs(s['rho_s_rel_err']) for s in st if 'rho_s_rel_err' in s]),
                pos=gmax([abs(s['dx_cells']) for s in sh]),
                over=gmax([max(s['overshoot'], s['undershoot']) for s in sh]),
                rh_mass=gmax([abs(s['rh_mass']) for s in sh if 'rh_mass' in s]),
                rh_mom=gmax([abs(s['rh_mom']) for s in sh if 'rh_mom' in s]),
                L1=d['final']['L1_rho'] if d.get('final') else None)


def star_table(R):
    """rho* error of each interface at the largest N (interface 1 at L/4, 2 at 3L/4)."""
    N = NS[-1]
    head = ['方式'] + ['%s 界面%d' % (prob, i) for prob in RIEMANN for i in (1, 2)]
    rows = []
    for v in SHOCK_VARIANTS:
        r = [v]
        for prob in RIEMANN:
            d = R.get('%s_N%d_%s' % (prob, N, v))
            st = {x['interface']: x['rho_s_rel_err'] for x in d.get('stars', [])} if d and d['ok'] else {}
            r += [e(st[i], '%+.1e') if i in st else '-' for i in (1, 2)]
        rows.append(r)
    return table(head, rows)


def conservation(R):
    head = ['方式', '本数', 'dM/dt', 'dP/dt', 'R_ζ', 'R_svv', 'ΔM/M0', 'ΔP', 'E 増加', '∫ε_ideal / ∫visc', '時間/step']
    rows = []
    wall_C = {}
    for v in SHOCK_VARIANTS + ['D4', 'lin', 'lin_own']:
        runs = [d for d in R.values() if d['case']['kind'] == 'riemann' and d['case']['var'] == v]
        if not runs:
            continue
        S = [stats(d) for d in runs]
        pr = lambda k: e(gmax([s.get(k) for s in S])) + ' (' + e(gmax([s.get(k + '_raw') for s in S])) + ')'
        nE = sum(s['n_E_up'] for s in S)
        ratio = [s['int_ideal_eps'] / s['int_visc'] for s in S if s.get('int_visc')]
        # wall time per step relative to C of the same (problem, N)
        rel = []
        for d in runs:
            a = d['attempts'][-1]
            c = R.get('%s_N%d_C' % (d['case']['prob'], d['case']['N']))
            if c is not None and a['steps'] > 0:
                ac = c['attempts'][-1]
                rel.append((a['wall'] / a['steps']) / (ac['wall'] / ac['steps']))
        rows.append([v, '%d/%d' % (sum(d['ok'] for d in runs), len(runs)), pr('max_dM'), pr('max_dP'), pr('max_R_zeta'),
                     pr('max_R_svv') if v not in ('lin', 'lin_own') else '-',
                     e(gmax([s['max_dM_over_M0'] for s in S])), e(gmax([s.get('max_dP_over_scale') for s in S])),
                     '%d (%s)' % (nE, e(gmax([s.get('max_E_up_over_E0') for s in S]))),
                     '%s – %s' % (e(min(ratio), '%+.1e'), e(max(ratio), '%+.1e')) if ratio else '-',
                     e(float(np.median(rel)), '%.2f') if rel else '-'])
    return table(head, rows)


def accuracy(R, prob):
    head = (['方式'] + ['速度 N=%d' % N for N in NS] + ['ρ* N=%d' % N for N in NS] +
            ['L1(ρ) N=%d' % NS[-1], 'L1 次数', 'over', '位置', '状態'])
    rows = []
    for v in SHOCK_VARIANTS + (['D4'] if prob == 'U2' else []):
        ds = [R.get('%s_N%d_%s' % (prob, N, v)) for N in NS]
        if all(d is None for d in ds):
            continue
        m = [shock(d) for d in ds]
        g = lambda k: [x[k] if x else None for x in m]
        bad = ['N=%d: %s' % (N, status(d)) for N, d in zip(NS, ds) if d is not None and status(d) != 'ok']
        rows.append([v] + [e(x) for x in g('speed')] + [e(x) for x in g('star')] +
                    [e(g('L1')[-1], '%.4f'), e(order(NS, g('L1')), '%.2f'), e(gmax(g('over')), '%.3f'),
                     e(gmax(g('pos')), '%.2f'), '; '.join(bad) or 'ok'])
    return table(head, rows)


def rh(R):
    head = ['問題', '方式'] + ['質量 N=%d' % N for N in NS] + ['運動量 N=%d' % N for N in NS]
    rows = []
    for prob in RIEMANN:
        for v in ('C', 'B', 'A', 'Cexp', 'Aexp'):
            m = [shock(R.get('%s_N%d_%s' % (prob, N, v))) for N in NS]
            rows.append([prob, v] + [e(x['rh_mass'] if x else None) for x in m] + [e(x['rh_mom'] if x else None) for x in m])
    return table(head, rows)


def svv(R):
    head = ['問題', 'N', 'SVV', 'dt', 'ρ*', '速度', 'over', 'L1(ρ)', '∫ε_svv / ∫visc', '状態']
    rows = []
    for prob in ('U2', 'U4'):
        for N in (1024, 2048):
            for v, lab in (('C', 'entropy'), ('lin', 'linear（ES の dt）'), ('lin_own', 'linear（自身の dt）')):
                d = R.get('%s_N%d_%s' % (prob, N, v))
                if d is None:
                    continue
                m, s, a = shock(d), stats(d), d['attempts'][-1]
                rows.append([prob, str(N), lab, '%.3e' % a['dt'], e(m['star'] if m else None), e(m['speed'] if m else None),
                             e(m['over'] if m else None, '%.3f'), e(m['L1'] if m else None, '%.4f'),
                             e(s['int_svv_eps'] / s['int_visc'], '%+.1e'), status(d)])
    return table(head, rows)


def same_operator(R, out):
    """max |rho - rho_C| / max rho_C of the final profiles of D3, D4 against C (D1)."""
    head = ['問題', 'N', '方式', 'max|Δρ| / max ρ']
    rows = []
    for prob in RIEMANN:
        for N in NS:
            pc = os.path.join(out, '%s_N%d_C.npz' % (prob, N))
            for v in ('D3', 'D4'):
                pv = os.path.join(out, '%s_N%d_%s.npz' % (prob, N, v))
                if os.path.exists(pc) and os.path.exists(pv):
                    a, b = np.load(pc), np.load(pv)
                    rows.append([prob, str(N), v, e(float(np.max(np.abs(b['rho'] - a['rho'])) / np.max(a['rho'])))])
    return table(head, rows)


def smooth(R):
    rows = []
    for N in (32, 64):
        for tag in ('C0', 'C2'):
            r = ['%d' % N, tag]
            for q in ('max_abs_amp_err', 'max_abs_phase_err'):
                for f in 'CBA':
                    d = R.get('sound_N%d_%s_%s' % (N, tag, f))
                    r.append(e(d['sound'][q] if d and 'sound' in d else None))
            for f in 'CBA':
                d = R.get('sound_N%d_%s_%s' % (N, tag, f))
                r.append(e(stats(d)['max_dM_over_M0'] if d else None))
            rows.append(r)
    t1 = table(['N', 'C_shock', '振幅 C', '振幅 B', '振幅 A', '位相 C', '位相 B', '位相 A', 'ΔM/M0 C', 'ΔM/M0 B', 'ΔM/M0 A'],
               rows)
    rows = []
    Ns = (32, 64, 128, 256, 512)
    for tag in ('C0_nosvv', 'C0_svv', 'C2_svv'):
        for f in 'CBA':
            ys = []
            for N in Ns:
                d = R.get('simple_N%d_%s_%s' % (N, tag, f))
                ys.append(d['simple']['final']['Linf_u'] if d and 'simple' in d else None)
            ords = ['-' if a is None or b is None or a <= 0 or b <= 0 else '%.2f' % (math.log(a / b) / math.log(2))
                    for a, b in zip(ys[:-1], ys[1:])]
            dM = gmax([stats(R['simple_N%d_%s_%s' % (N, tag, f)])['max_dM_over_M0']
                       for N in Ns if 'simple_N%d_%s_%s' % (N, tag, f) in R])
            rows.append([tag, f] + [e(y) for y in ys] + ['/'.join(ords), e(dM)])
    t2 = table(['設定', '方式'] + ['L∞(u) N=%d' % N for N in Ns] + ['次数（隣どうし）', 'max ΔM/M0'], rows)
    return t1, t2


def failures(R):
    rows = [[n, status(d)] for n, d in sorted(R.items()) if status(d) != 'ok']
    return table(['ケース', '試行'], rows) if rows else '（なし）\n'


def fig_profiles(R, out, prob='U4', N=2048):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rL, rR, uL, uR, te = RIEMANN[prob]
    L = 2.0 * math.pi
    fig, ax = plt.subplots(1, 3, figsize=(12.0, 3.6), constrained_layout=True)
    x_ex = np.linspace(0.0, L, 8192, endpoint=False)
    rho_ex, _, _ = rx.periodic_solution(x_ex, te, rL, uL, rR, uR, L, 1.0)
    wins = [(0.0, L, 'full domain'), (0.25 * L - 0.25, 0.25 * L + 0.25, 'collision shocks, x = L/4'),
            (0.75 * L - 0.7, 0.75 * L + 0.7, 'rarefaction star, x = 3L/4')]
    for a, (x0, x1, title) in zip(ax, wins):
        sel = (x_ex >= x0) & (x_ex <= x1)
        a.plot(x_ex[sel], rho_ex[sel], color=GRAY, lw=1.0, label='exact')
        for f in 'CBA':
            p = os.path.join(out, '%s_N%d_%s.npz' % (prob, N, f))
            if not os.path.exists(p):
                continue
            z = np.load(p)
            s = (z['x'] >= x0) & (z['x'] <= x1)
            a.plot(z['x'][s], z['rho'][s], color=COLOR[f], lw=1.4, marker=MARK[f] if x1 - x0 < 2 else None, ms=3,
                   markevery=max(1, int(s.sum() // 40)), label=f)
        a.set_title(title, fontsize=10)
        a.set_xlabel('x')
        a.set_xlim(x0, x1)
        a.grid(color='#e5e4df', lw=0.6)
    ax[0].set_yscale('log')
    ax[0].set_ylabel(r'$\rho$')
    ax[2].set_ylim(0.0, 0.05)
    ax[0].legend(frameon=False, fontsize=9)
    fig.suptitle('%s, N = %d, t = %g: formulations A, B, C (same dt, ES-SVV, C_shock = 2)' % (prob, N, te), fontsize=10)
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, 'stage2_%s_profiles.png' % prob), dpi=130)
    plt.close(fig)


def fig_convergence(R):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    probs = list(RIEMANN)
    fig, ax = plt.subplots(2, len(probs), figsize=(12.0, 5.6), sharex=True, constrained_layout=True)
    for j, prob in enumerate(probs):
        for f in 'CBA':
            for v, ls in ((f, '-'), (f + 'exp', '--')):
                if v not in SHOCK_VARIANTS:
                    continue
                m = [shock(R.get('%s_N%d_%s' % (prob, N, v))) for N in NS]
                ds = [R.get('%s_N%d_%s' % (prob, N, v)) for N in NS]
                star1 = [next((abs(x['rho_s_rel_err']) for x in d['stars'] if x['interface'] == 1), None)
                         if d and d['ok'] else None for d in ds]
                for i, ys in enumerate(([x['speed'] if x else None for x in m], star1)):
                    pts = [(N, y) for N, y in zip(NS, ys) if y]
                    if pts:
                        ax[i, j].loglog(*zip(*pts), color=COLOR[f], ls=ls, lw=1.4, marker=MARK[f], ms=5)
        ax[0, j].set_title(prob, fontsize=10)
        ax[1, j].set_xlabel('N')
        for i in range(2):
            ax[i, j].grid(color='#e5e4df', lw=0.6, which='both')
            ax[i, j].set_xticks(NS)
            ax[i, j].set_xticklabels([str(N) for N in NS])
            ax[i, j].minorticks_off()
    ax[0, 0].set_ylabel('|shock speed error|')
    ax[1, 0].set_ylabel(r'|$\rho_*$ error| at x = L/4')
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], color=COLOR[v[0]], ls='--' if v.endswith('exp') else '-', marker=MARK[v[0]], ms=5, lw=1.4,
                 label=v + (' (expanded viscosity)' if v.endswith('exp') else '')) for v in ('C', 'Cexp', 'B', 'A', 'Aexp')]
    fig.legend(handles=hs, loc='outside right center', frameon=False, fontsize=9)
    fig.suptitle('Relative shock speed error (top, worst shock) and relative error of the star density between '
                 'the shocks at x = L/4 (bottom)', fontsize=10)
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, 'stage2_shock_convergence.png'), dpi=130)
    plt.close(fig)


def main(out):
    R = load(out)
    t_sound, t_simple = smooth(R)
    parts = ['# Stage 2 tables (%d cases)\n' % len(R),
             '## 保存量（衝撃波、全 (問題, N) の最大値。段のフィルタを掛けた後 (掛ける前)）\n', conservation(R)]
    for prob in RIEMANN:
        parts += ['## 衝撃波の精度: %s\n' % prob, accuracy(R, prob)]
    parts += ['## 星状態 ρ* の相対誤差（N = %d、界面ごと）\n' % NS[-1], star_table(R)]
    parts += ['## RH 残差（当てはめた速度 / 測った状態から作る速度 − 1）\n', rh(R),
              '## SVV の形\n', svv(R), '## D3, D4 と C（D1）の最終分布の差\n', same_operator(R, out),
              '## 音波\n', t_sound, '## 単純波（t = t_b / 2）\n', t_simple, '## 破綻と再実行\n', failures(R)]
    with open(os.path.join(out, 'tables.md'), 'w') as f:
        f.write('\n'.join(parts))
    fig_profiles(R, out)
    fig_convergence(R)
    print('wrote', os.path.join(out, 'tables.md'), 'and', FIG)
