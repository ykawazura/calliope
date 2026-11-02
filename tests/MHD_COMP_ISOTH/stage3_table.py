"""Markdown tables and figures of the Stage 3 runs (python3 stage3.py table).

Reads runs/stage3/<case>.json and .npz of stage3.py; writes runs/stage3/tables.md and fig/stage3_*.png next to
this file. The values of a case are those of its last attempt (the dt/2, dt/4 reruns of a failed run).
"""
import os
import json
import glob
import math
import numpy as np
import mhd_exact as mx
import fv_mhd1d as fv
import stage3 as s3
from stage2_table import e, stats, status, gmax, order, table, COLOR, MARK, GRAY

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig')
FAMS = ('perp', 'coplanar', 'td', 'rd', 'compound')
NAMES = ['rho', 'u_x', 'u_y', 'u_z', 'B_y', 'B_z']


def load(out):
    R = {}
    for f in glob.glob(os.path.join(out, '*.json')):
        with open(f) as fh:
            d = json.load(fh)
        R[d['case']['name']] = d
    return R


def base(c):
    """Case name without _N<N>."""
    return c['name'].replace('_N%d' % c['N'], '')


def groups(R, fams, Ns=s3.NS):
    """[(base name, {N: run})] in the order of stage3.cases() (C_eta = 0 first, then each C_eta of R)."""
    out, seen = [], {}
    etas = sorted({d['case']['C_eta'] for d in R.values()} | {0.0})
    for c in s3.cases(etas):
        if c['fam'] in fams and c['N'] in Ns and c['name'] in R:
            b = base(c)
            if b not in seen:
                seen[b] = {}
                out.append((b, seen[b]))
            seen[b][c['N']] = R[c['name']]
    return out


def ok_an(d):
    return d['analysis'] if d is not None and d['ok'] else None


def runs(R):
    rows = []
    for f in FAMS:
        ds = [d for d in R.values() if d['case']['fam'] == f]
        bad = ['%s (%s)' % (d['case']['name'], status(d)) for d in sorted(ds, key=lambda d: d['case']['name'])
               if status(d) != 'ok']
        rows.append([f, str(len(ds)), str(sum(d['ok'] for d in ds)), '; '.join(bad) or '-'])
    return table(['family', 'cases', 'reached t_end', 'not at the first dt'], rows)


def conservation(R):
    """Worst values over the cases of each family."""
    rows = []
    for f in FAMS:
        S = [stats(d) for d in R.values() if d['case']['fam'] == f and d['ok']]
        if not S:
            continue
        g = lambda ks: gmax([gmax([s.get(k) for k in ks]) for s in S])
        mn = lambda k: min(s[k] for s in S if s.get(k) is not None)
        rows.append([f, e(g(['drift_M'])), e(g(['drift_Px', 'drift_Py', 'drift_Pz'])), e(g(['drift_By', 'drift_Bz'])),
                     e(g(['rate_dM_raw'])), e(g(['rate_dPx_raw', 'rate_dPy_raw', 'rate_dPz_raw'])),
                     e(g(['rate_dBy_raw', 'rate_dBz_raw'])), e(g(['R_zeta'])), e(g(['R_svv'])), e(g(['R_eta'])),
                     '%d' % sum(s['n_E_up'] for s in S), '%d' % sum(s['n_svv_pos'] for s in S),
                     e(g(['max_svv_b_eps']), '%+.1e'), '%.3g' % mn('rho_min'), '%.2f' % g(['max_dt_over_limit']),
                     '%.2f' % g(['max_dt_lam_svv'])])
    return table(['family', 'ΔM/M', '\\|ΔP\\|', '\\|Δ∫B\\|', 'dM/dt raw', 'dP/dt raw', 'd∫B/dt raw', 'R_ζ', 'R_svv',
                  'R_η', 'E up (steps)', 'ε_svv > 0', 'max ε_svv,B', 'ρ_min', 'max dt/dt_lim', 'max dt λ_svv'], rows)


def front_rows(R, fams, key, fmt='%.1e', signed=False):
    """One row per front (case, kind) of the families, one column per N."""
    rows = []
    for b, byN in groups(R, fams):
        fr0 = next((ok_an(d)['fronts'] for d in byN.values() if ok_an(d) and ok_an(d).get('fronts')), None)
        if fr0 is None:
            continue
        for k, f0 in enumerate(fr0):
            r = ['%s %d %s' % (b, k, f0['kind'])]
            for N in s3.NS:
                an = ok_an(byN.get(N))
                v = an['fronts'][k].get(key) if an and an.get('fronts') else None
                if isinstance(v, list):
                    v = gmax([x for x in v if x is not None])
                if v is not None and not signed:
                    v = abs(v)
                r.append(e(v, '%+.1e' if signed else fmt))
            rows.append(r)
    return table(['front'] + ['N = %d' % N for N in s3.NS], rows)


def l1_rows(R, fams, comps=range(6)):
    rows = []
    for b, byN in groups(R, fams):
        ys = []
        for N in s3.NS:
            an = ok_an(byN.get(N))
            ys.append(sum(an['l1'][i] for i in comps) if an and 'l1' in an else None)
        if all(y is None for y in ys):
            continue
        p = order(s3.NS, ys)
        rows.append([b] + [e(y) for y in ys] + [e(p, '%.2f')])
    return table(['case'] + ['N = %d' % N for N in s3.NS] + ['order'], rows)


def rd(R):
    rows = []
    for b, byN in groups(R, ('rd',)):
        for N in s3.NS:
            an = ok_an(byN.get(N))
            if an is None:
                rows.append([b, str(N), status(byN[N]) if N in byN else '-'] + [''] * 9)
                continue
            dv = an['dev_max']
            rows.append([b, str(N), e(an['speed_err'], '%+.1e'), e(an['shift_max_dl'], '%.2f'), e(dv['bt']),
                         e(dv['rho']), e(dv['ux']), e(dv['alfv']), e(dv['D']), e(dv['zeta']), e(dv.get('z')),
                         e(an['l1'][4] + an['l1'][5])])
    return table(['case', 'N', 'speed err', 'max shift/Δ', '\\|B_t\\| - 1', 'ρ - 1', 'u_x', 'u_t + B_t', 'D', 'ζ',
                  'u_z, B_z', 'L1 B_t'], rows)


def td(R):
    rows = []
    for b, byN in groups(R, ('td',)):
        for N in s3.NS:
            an = ok_an(byN.get(N))
            if an is None:
                continue
            d0, dv, w = an['dev_initial'], an['dev_max'], an['width_dl']
            rows.append([b, str(N), e(an['speed_err'], '%+.1e'), e(d0['ptot']), e(dv['ptot']), e(dv['ux']),
                         e(dv['ut']), e(dv['D']), e(dv['zeta']), '%.2f / %.2f / %.2f' % (w['data'], w['initial'],
                                                                                       w['final']),
                         e(an['l1'][0])])
    return table(['case', 'N', 'speed err', 'P_tot t=0', 'max P_tot', 'max \\|u_x - U\\|', 'max \\|u_t\\|',
                  'max \\|D\\|', 'max ζ', 'width/Δ (data / t=0 / t_end)', 'L1 ρ'], rows)


def compound(R):
    """L1 distance (sum over the six variables) from the FV reference at each T_REF."""
    rows = []
    Ns = s3.NS + (4096,)
    for b, byN in groups(R, ('compound',), Ns):
        an0 = next((ok_an(d) for d in byN.values() if ok_an(d) and ok_an(d).get('l1_fv')), None)
        if an0 is None:
            continue
        for k, x in enumerate(an0['l1_fv']):
            ys = []
            for N in Ns:
                an = ok_an(byN.get(N))
                ys.append(sum(an['l1_fv'][k]['l1']) if an and an.get('l1_fv') else None)
            rows.append(['%s (FV N = %d)' % (b, an0['fv_N']), '%.3g' % x['t']] + [e(y) for y in ys] +
                        [e(order(Ns[-3:], ys[-3:]), '%.2f'), e(x.get('fv_self'))])
    return table(['case', 't'] + ['N = %d' % N for N in Ns] + ['order (3 finest)', 'FV N/2 - N'], rows)


def compound_vars(R, N=2048):
    """L1 per variable from the FV reference at the last T_REF."""
    rows = []
    for b, byN in groups(R, ('compound',), (N,)):
        an = ok_an(byN.get(N))
        if an and an.get('l1_fv'):
            rows.append([b] + [e(v) for v in an['l1_fv'][-1]['l1']])
    return table(['case (N = %d)' % N] + NAMES, rows)


# the finest FV run of the coplanar cases whose dissipative limit is not the ideal solution (runs/stage3/fvwin/):
# I1 (stage3.py fvwin and fv_i1.py) and RDpc (fv_rdpc.py, one period), with the next coarser one for the self-check
FV_NONUNIQUE = {'I1': ('I1_N16384_conv.npz', 'I1_N4096.npz'), 'RDpc': ('RDpc_N16384.npz', 'RDpc_N4096.npz')}


def fv_nonunique(R, out):
    """L1 distance of rho (and the sum over the six variables) at t_end from the finest FV solution, for I1 and RDpc.
    The ideal solution (intermediate shock, 180 deg RD) is not the limit of the dissipative coplanar problem, so the
    converged FV solution is the reference; the last column is the FV N/4 - N distance on the same samples."""
    rows = []
    for b, byN in groups(R, ('coplanar', 'rd')):
        key = b.split('_')[0]
        if key not in FV_NONUNIQUE:
            continue
        fns = [os.path.join(out, 'fvwin', f) for f in FV_NONUNIQUE[key]]
        if not os.path.exists(fns[0]):
            fns = fns[1:]
        if not os.path.exists(fns[0]):
            continue
        Wf = np.load(fns[0])['W']
        r_rho, r_all = [], []
        for N in s3.NS:
            d = npz(out, '%s_N%d%s' % (key, N, b[len(key):])) if ok_an(byN.get(N)) else None
            if d is None:
                r_rho.append(None)
                r_all.append(None)
                continue
            Ws = Wf[:, ::Wf.shape[1] // N]
            l1 = [float(np.mean(np.abs(d['W'][-1][i] - Ws[i]))) for i in range(6)]
            r_rho.append(l1[0])
            r_all.append(sum(l1))
        fv_self = None
        if len(fns) > 1 and os.path.exists(fns[1]):
            Wc = np.load(fns[1])['W']
            fv_self = float(np.mean(np.abs(Wc[0] - Wf[0, ::Wf.shape[1] // Wc.shape[1]])))
        rows.append([b, 'FV N = %d' % Wf.shape[1]] + ['%s / %s' % (e(a), e(s)) for a, s in zip(r_rho, r_all)] +
                    [e(order(s3.NS, r_rho), '%.2f'), e(fv_self)])
    return table(['case', 'reference'] + ['N = %d (ρ / all)' % N for N in s3.NS] + ['order ρ', 'FV ρ: N/4 - N'],
                 rows)


def energy(R, N=1024):
    """Energy changes over the run relative to E0 (instruction 6): shock viscosity, resistivity, the SVV of (rho, m)
    and B, the ideal terms and the time integration (dE - sum dt rk_eps)."""
    rows = []
    for c in s3.cases(sorted({d['case']['C_eta'] for d in R.values()} | {0.0})):
        d = R.get(c['name'])
        if c['N'] != N or d is None or not d['ok']:
            continue
        s = stats(d)
        E0 = s['E0']
        rows.append([base(c), '%.4g' % E0, e(s['dE'] / E0, '%+.1e')] +
                    [e(s['int_rk_%s_eps' % k] / E0, '%+.1e') for k in ('visc', 'eta', 'svv', 'svv_b', 'ideal')] +
                    [e(s['time_err'] / E0, '%+.1e')])
    return table(['case (N = %d)' % N, 'E0', 'ΔE', 'shock visc.', 'η_shock', 'SVV (ρ, m)', 'SVV B', 'ideal', 'time'],
                 rows)


# ------------------------------------------------------------ figures
def npz(out, name):
    fn = os.path.join(out, name + '.npz')
    return np.load(fn) if os.path.exists(fn) else None


def fig_coplanar(R, out, N=1024, figdir=FIG):
    """rho and B_y of the five coplanar shocks at t_end with the exact solution on the valid points."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    probs = [p for p in mx.COPLANAR if npz(out, '%s_N%d' % (p, N)) is not None]
    if not probs:
        return
    fig, ax = plt.subplots(2, len(probs), figsize=(3.2 * len(probs), 5.0), squeeze=False)
    x = s3.grid(N)
    for j, p in enumerate(probs):
        d = npz(out, '%s_N%d' % (p, N))
        c = dict(fam='coplanar', prob=p, N=N, C_eta=0.0)
        t, W = float(d['t'][-1]), d['W'][-1]
        We, valid = s3.exact(c, x, t)
        S = mx.coplanar_case(p)[0]
        xs = s3.cop_interval(p)[0] + S * t
        dx = s3.pdist(x, xs)
        sel = np.abs(dx) < 0.08
        o = np.argsort(dx[sel])
        for i, (k, nm) in enumerate(((0, 'ρ'), (4, 'B_y'))):
            a = ax[i, j]
            v = valid[sel][o]
            a.plot(dx[sel][o], np.where(v, We[k][sel][o], np.nan), color=GRAY, lw=1.0, label='exact')
            a.plot(dx[sel][o], W[k][sel][o], color=COLOR['C'], lw=1.5, marker='.', ms=3, label='proto')
            a.set_xlabel('x − x_s')
            a.set_ylabel(nm)
            if i == 0:
                a.set_title('%s (%s), t = %.3g' % (p, mx.coplanar_case(p)[3], t))
    ax[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    os.makedirs(figdir, exist_ok=True)
    fig.savefig(os.path.join(figdir, 'stage3_coplanar_N%d.png' % N), dpi=150)
    plt.close(fig)


def fig_compound(R, out, Ns=(512, 2048), figdir=FIG):
    """rho, u_y and B_y of C1-C3 at the last T_REF with the FV reference."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 3, figsize=(10.5, 7.5))
    shade = {Ns[0]: '#9ec5f0', Ns[-1]: COLOR['C']}
    for j, p in enumerate(mx.COMPOUND):
        ref = s3.fv_reference(p)
        for i, (k, nm) in enumerate(((0, 'ρ'), (2, 'u_y'), (4, 'B_y'))):
            a = ax[i, j]
            if ref is not None:
                Nf, tr, Wr, _ = ref
                a.plot(np.arange(Nf) / Nf, Wr[-1][k], color=GRAY, lw=1.0, label='FV N = %d' % Nf)
            for N in Ns:
                d = npz(out, '%s_N%d' % (p, N))
                if d is not None:
                    a.plot(s3.grid(N), d['W'][-1][k], color=shade[N], lw=1.0, label='proto N = %d' % N)
            a.set_xlim(0.0, 0.5)
            a.set_ylabel(nm)
            if i == 0:
                a.set_title('%s, t = %.3g' % (p, fv.T_REF[-1]))
            if i == 2:
                a.set_xlabel('x')
    ax[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    os.makedirs(figdir, exist_ok=True)
    fig.savefig(os.path.join(figdir, 'stage3_compound.png'), dpi=150)
    plt.close(fig)


def fig_rd(R, out, N=512, figdir=FIG):
    """B_y, B_z and |B_t| of the RDs after one period with the exact (initial) profile."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names = [b for b in ('RDh_w0', 'RDh_w8', 'RDp_w0', 'RDp_w8', 'RDpc') if npz(out, '%s_N%d' % (b, N)) is not None]
    if not names:
        return
    fig, ax = plt.subplots(2, len(names), figsize=(3.0 * len(names), 4.8), squeeze=False)
    x = s3.grid(N)
    for j, b in enumerate(names):
        d = npz(out, '%s_N%d' % (b, N))
        W0, W = d['W0'], d['W'][-1]
        for i, (f, nm) in enumerate(((lambda W: W[4], 'B_y'), (lambda W: np.hypot(W[4], W[5]), '|B_t|'))):
            a = ax[i, j]
            a.plot(x, f(W0), color=GRAY, lw=1.0, label='exact')
            a.plot(x, f(W), color=COLOR['C'], lw=1.0, label='proto')
            a.set_xlim(0.15, 0.35)
            a.set_ylabel(nm)
            if i == 0:
                a.set_title('%s, N = %d, t = 1' % (b, N))
            else:
                a.set_xlabel('x')
    ax[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    os.makedirs(figdir, exist_ok=True)
    fig.savefig(os.path.join(figdir, 'stage3_rd_N%d.png' % N), dpi=150)
    plt.close(fig)


def main(out, figdir=FIG):
    R = load(out)
    parts = [
        ('Runs', runs(R)),
        ('Conservation and semi-discrete identities (worst case of each family)', conservation(R)),
        ('Front speed error (S_fit - S)/|S|', front_rows(R, ('perp', 'coplanar', 'compound'), 'dS_rel', signed=True)),
        ('Plateau error next to each front: max |mean - exact| / max jump', front_rows(R, ('perp', 'coplanar', 'compound'),
                                                                                       'plat_err')),
        ('RH residual with the fitted speed (max over the equations, relative)',
         front_rows(R, ('perp', 'coplanar', 'compound'), 'rh')),
        ('rho overshoot / jump', front_rows(R, ('perp', 'coplanar', 'compound'), 'over_rho')),
        ('Front offset at t_end / Δ', front_rows(R, ('perp', 'coplanar', 'compound'), 'off_final_dl', signed=True)),
        ('L1 error on the valid points (sum over the six variables)', l1_rows(R, ('perp', 'coplanar', 'compound'))),
        ('Rotational discontinuities (one period; max over the snapshots)', rd(R)),
        ('Tangential discontinuities (one period; max over the snapshots)', td(R)),
        ('Compound Riemann problems: L1 distance from the FV reference', compound(R)),
        ('Compound Riemann problems: L1 per variable at the last T_REF', compound_vars(R)),
        ('I1 and RDpc: L1 distance from the converged FV solution at t_end', fv_nonunique(R, out)),
        ('Energy decomposition (relative to E0)', energy(R)),
    ]
    with open(os.path.join(out, 'tables.md'), 'w') as f:
        for h, t in parts:
            f.write('## %s\n\n%s\n' % (h, t))
    fig_coplanar(R, out, figdir=figdir)
    fig_compound(R, out, figdir=figdir)
    fig_rd(R, out, figdir=figdir)
    print('wrote %s and %s/stage3_*.png' % (os.path.join(out, 'tables.md'), figdir))
