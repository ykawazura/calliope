#!/usr/bin/env python3
"""Stage 15 analysis of the runs of gen_stage15.py (ES-SVV against the standard dissipation in 3D decaying
turbulence). The criteria were fixed before the runs (tasks/todo.md, Stage 15).

    python3 analyze_stage15.py selftest              checks of the functions below on synthetic data
    python3 analyze_stage15.py pre [base]            phase 0: P0a (budget error with dt and dt/2, reported), P0b
                                                     (planes of P = 8 against P = 64 after 0 and 10 steps,
                                                     TOL_DEC10), eps_ref per M_s from <prob>_N64_ES; writes
                                                     <base>/pre.json (ok: P0b passed and both ES runs complete
                                                     with R <= TOL_RES)
    python3 analyze_stage15.py runs <batch> [base]   C (completion, failure) and E (budget) of every case of the
                                                     batch: <base>/<batch>/summary.json; batch a also writes
                                                     <base>/a.json (r_star, r0_ok, diag), batch b <base>/b.json
                                                     (r_C)
    python3 analyze_stage15.py watch <batch> [base]  the stop rule: an 'incomplete' case whose cost projected to
                                                     its t_end exceeds COST_CAP gets the exit_status stopped_cost
                                                     and cost_projection.json. Run it only while no job runs the
                                                     cases of the batch (the job script moves exit_status away
                                                     while it runs a case)
    python3 analyze_stage15.py spec [base]           per-record and window-mean spectra of the complete runs of
                                                     pre, a, b, c and of the reused runs ES256, S8: <base>/spec.json
    python3 analyze_stage15.py judge [base]          K (k_x, D, sigma, delta, B*, verdict), the cost rule, the
                                                     reference judgement of M_s = 0.5 and the self-convergence
                                                     reports from spec.json: <base>/judge.json and fig_*.png;
                                                     the H reference on 256^3 is the rung of gen_stage15.r_256

The default base is runs/stage15; the reused runs are read from runs/stage6q and runs/stage8 next to it.
"""
import os
import re
import sys
import json
import glob
import math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_ic
import analyze_stage4 as a4
from analyze_stage1 import load, load_runs, run_dirs
from analyze_stage6 import chunk_costs
import gen_stage15 as g15

WIN = {'W1': (0.75, 1.25), 'W2': (1.5, 2.0)}   # tau0
WIN_TOL = 0.01      # tau0: a record belongs to a window when t/tau0 is in [lo - WIN_TOL, hi + WIN_TOL]
WIN_END = 0.05      # tau0: a window of a run is used when its last out3d record is at t >= hi - WIN_END - WIN_TOL
                    # (half the record spacing of S8; a run stops up to dt/2 before t_end)
T_TOL = 1.25        # |ln T| <= ln T_TOL: the spectra agree
N0 = 3              # first shell of the scan of k_x and of D
CAP_K = 42.0        # cap of k_x on 128^3 (k_c(128) = 42.7)
CAP_REF = 84.0      # cap of k_ref,W on 256^3
D_HI = 32           # upper end of D on 128^3
CAP_K64, D_HI64 = 21.0, 16    # the same for the self-convergence (64, 128)
COST_CAP = {64: 16.0, 128: 96.0, 256: 450.0}   # core-h per run
PROJ_ROWS = 200     # rows of cfl.dat for the mean dt of the projection
QJ = ('E_B', 'E_m')                 # components of the judgement
QS = ('E_B', 'E_m', 'P_s', 'E_sru')  # components stored
ES256 = os.path.join('stage6q', 'turb', 'turb1_N256')
S8 = os.path.join('stage8', 'coef', 'e0.5_C2')
SHARES = dict(ideal=(9,), visc=(10,), svv=(20, -21), svv_b=(21,), eta=(26,), nu=(16,), res=(15,), lmd=(17,),
              force=(14,))
INVALID = re.compile(r'\[proc\s+(\d+)\] invalid state at t =\s+(\S+), step\s+(\d+), stage\s+(\d+): (\w+)'
                     r'\(iy, iz, ix\) = \(\s*(\d+),\s*(\d+),\s*(\d+)\) =\s+(\S+)')
ABORT = re.compile(r'\[proc\s+(\d+)\] ABORT: (.*)')


# ---------------------------------------------------------------- the quantities of K
def window_weights(ts, lo, hi):
    """[(i, w)] of the times ts (tau0) in [lo - WIN_TOL, hi + WIN_TOL]: w is the length of the part of [lo, hi]
    nearest to ts[i] over the window length (nearest-record weights; trapezoid weights for equal spacing).
    A record within 1e-9 t of an earlier one is left out."""
    sel = []
    for i, t in enumerate(ts):
        if lo - WIN_TOL <= t <= hi + WIN_TOL and not any(abs(t - ts[j]) < 1e-9*max(abs(t), 1e-300) for j in sel):
            sel.append(i)
    sel.sort(key=lambda i: ts[i])
    if not sel:
        return []
    t = np.array([ts[i] for i in sel])
    edge = np.concatenate(([lo], 0.5*(t[1:] + t[:-1]), [hi]))
    w = np.clip(edge[1:], lo, hi) - np.clip(edge[:-1], lo, hi)
    return [(i, float(x/(hi - lo))) for i, x in zip(sel, w)]


def ratio(a, b):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.asarray(a, float)/np.asarray(b, float)


def k_cross(T, cap, n0=N0):
    """(k_x, censored): the first shell n >= n0 with |ln T(n)| > ln T_TOL, interpolated linearly in |ln T|
    between n - 1 and n (k_x = n0 when n0 already exceeds; a non-finite T(n) gives n - 1); cap when no shell up
    to floor(cap) exceeds (censored)."""
    thr = math.log(T_TOL)
    with np.errstate(divide='ignore', invalid='ignore'):
        y = np.abs(np.log(np.asarray(T, float)))
    for n in range(n0, min(int(math.floor(cap)), len(y) - 1) + 1):
        if not np.isfinite(y[n]) or y[n] > thr:
            if n == n0:
                return float(n0), False
            if not np.isfinite(y[n]):
                return float(n - 1), False
            return float(n - 1 + (thr - y[n - 1])/(y[n] - y[n - 1])), False
    return float(cap), True


def dev(T, n_hi, n0=N0):
    """D: the trapezoid integral of |ln T(n)| over ln n on n = n0, ..., n_hi, over ln(n_hi/n0) (nan for
    n_hi < n0)."""
    if n_hi < n0:
        return float('nan')
    n = np.arange(n0, int(n_hi) + 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        y = np.abs(np.log(np.asarray(T, float)[n]))
    if not np.all(np.isfinite(y)):
        return float('inf')
    if len(n) == 1:
        return float(y[0])
    x = np.log(n)
    return float(np.sum(0.5*(y[1:] + y[:-1])*np.diff(x))/(x[-1] - x[0]))


def sigma(kappas):
    return float(np.std(kappas, ddof=1)) if len(kappas) > 1 else 0.0


def verdict(kES, kB, dES, dB, delta, need):
    """'ES better', 'ES worse' or 'undecided' over the pairs, with the counts: plus (kES - kB > delta), minus
    (kES - kB < -delta), D_ES < D_B and D_B < D_ES."""
    diff = np.asarray(kES) - np.asarray(kB)
    delta = np.asarray(delta)
    c = dict(n=len(diff), need=need, plus=int(np.sum(diff > delta)), minus=int(np.sum(diff < -delta)),
             d_es=int(np.sum(np.asarray(dES) < np.asarray(dB))), d_b=int(np.sum(np.asarray(dB) < np.asarray(dES))))
    if c['plus'] >= need and c['minus'] == 0 and c['d_es'] >= need:
        v = 'ES better'
    elif c['minus'] >= need and c['plus'] == 0 and c['d_b'] >= need:
        v = 'ES worse'
    else:
        v = 'undecided'
    return v, c


def need_of(npairs):
    return 3*npairs//4


def gmean(x):
    return float(np.exp(np.mean(np.log(x))))


def cost_ok(cES, cB, kES, kB):
    """c_ES/c_B* < (kbar_ES/kbar_B*)^4 with kbar the geometric mean of k_x over the pairs."""
    lhs, rhs = cES/cB, (gmean(kES)/gmean(kB))**4
    return dict(c_ratio=lhs, k_ratio4=rhs, ok=bool(lhs < rhs))


# ---------------------------------------------------------------- spectra
_GRIDS = {}


def grid(N):
    if N not in _GRIDS:
        _GRIDS[N] = gen_ic.Grid(N, N, N)
    return _GRIDS[N]


def spectra(g, fk):
    """E_B, E_m (m = rho u, the coefficients of out3d), P_s (s = ln rho) and E_sru (m/sqrt(rho)) of a record of
    analyze_stage4.out3d_at, with the shells of analyze_stage4.shell_spectra_3d."""
    K = np.sqrt(g.kx**2 + g.ky**2 + g.kz**2)
    n = np.rint(K).astype(int)
    w = np.broadcast_to(np.where((g.iky == 0) | (g.iky == g.nly//2), 1.0, 2.0).reshape(1, 1, -1), K.shape)
    nb = g.nlx//2 + 1
    sel = n < nb

    def shell(p):
        return np.bincount(n[sel], weights=(w*p)[sel], minlength=nb)[:nb]

    def e(vk):
        return shell(sum(0.5*np.abs(v)**2 for v in vk))

    rho = g.ifft(fk['rho'])
    sr = np.sqrt(rho)
    return dict(E_B=e([fk[q] for q in ('bx', 'by', 'bz')]), E_m=e([fk[q] for q in ('mx', 'my', 'mz')]),
                P_s=shell(np.abs(g.fft(np.log(rho)))**2*(K > 0)),
                E_sru=e([g.fft(g.ifft(fk[q])/sr) for q in ('mx', 'my', 'mz')]))


# ---------------------------------------------------------------- runs: cost, completion, failure, budget
def nproc_of(case):
    try:
        return int(open(os.path.join(case, 'nproc')).read())
    except OSError:
        return int(a4.case_info(case)['nproc'])


def t_end_of(case):
    """t_end of the job script (file t_end), else nrun dt of case.json."""
    try:
        return float(open(os.path.join(case, 't_end')).read())
    except OSError:
        m = a4.case_info(case)
        return m['nrun']*m['dt']


def cost_to(case, T=None):
    """core-h (wall of mpirun x ranks) of the runs of the case; with T, a run counts the fraction of its rows
    with t + dt/2 < T. (core-h, number of runs without a wall line)."""
    np_, c, miss = nproc_of(case), 0.0, 0
    cc = chunk_costs(case)
    for x in cc:
        if x['wall'] is None:
            miss += 1
            continue
        f = 1.0
        if T is not None:
            cfl = load(os.path.join(case, x['run'], 'cfl.dat'))
            f = float(np.mean(cfl[:, 0] + 0.5*cfl[:, 1] < T))
        c += x['wall']*np_/3600.0*f
    return c, miss + len(run_dirs(case)) - len(cc)


def projection(case, t_stop):
    """Cost of the case projected to t_stop: spent + steps left x wall per step of its last run x ranks, the
    steps left from the mean dt of the last PROJ_ROWS rows of cfl.dat and t + dt of the last row."""
    spent, miss = cost_to(case)
    cc = [x for x in chunk_costs(case) if x['wall'] is not None]
    cfl = np.concatenate(load_runs(case, 'cfl.dat'))
    t_now = float(cfl[-1, 0] + cfl[-1, 1])
    dtm = float(np.mean(cfl[-PROJ_ROWS:, 1]))
    left = max(0.0, (t_stop - t_now)/dtm)
    wps = cc[-1]['wall']/cc[-1]['steps']
    return dict(spent=spent, projected=spent + left*wps*nproc_of(case)/3600.0, steps_left=left, dt_mean=dtm,
                t_now=t_now, t_stop=t_stop, wall_per_step=wps, runs_without_wall=miss)


def failure(case):
    """The invalid-state lines (check_state of advance.F90) and ABORT lines (abort_run) of the .std of all runs:
    the first invalid state (earliest t, then step, stage, rank), the fields named, and the abort messages."""
    inv, ab = [], []
    for d in run_dirs(case):
        for f in glob.glob(os.path.join(d, '*.std')):
            for line in open(f, errors='replace'):
                m = INVALID.search(line)
                if m:
                    p, t, st, sg, q, iy, iz, ix, v = m.groups()
                    inv.append(dict(proc=int(p), t=float(t), step=int(st), stage=int(sg), field=q,
                                    iyzx=(int(iy), int(iz), int(ix)), value=v, run=os.path.relpath(d, case)))
                m = ABORT.search(line)
                if m:
                    ab.append(dict(proc=int(m.group(1)), msg=m.group(2).strip()))
    first = min(inv, key=lambda x: (x['t'], x['step'], x['stage'], x['proc'])) if inv else None
    return dict(first=first, fields=sorted({x['field'] for x in inv}), n_invalid=len(inv),
                abort=sorted({x['msg'] for x in ab}), n_abort=len(ab))


def decimate(a, n=400, tail=20):
    i = sorted(set(range(0, len(a), max(1, len(a)//n))) | set(range(max(0, len(a) - tail), len(a))))
    return a[i]


def pileup(E, kc):
    """E(floor(k_c)) over min E(n) on n in [ceil(k_c/2), floor(k_c)] (> 1: a rise before k_c)."""
    lo, hi = int(math.ceil(kc/2)), int(math.floor(kc))
    return float(E[hi]/np.min(E[lo:hi + 1]))


def run_summary(case):
    """C and E of one run."""
    import analyze_stage3 as ast3
    m = a4.case_info(case)
    N, tau0 = m['N'], m['tau0']
    r = dict(info={k: m.get(k) for k in ('problem', 'N', 'tau0', 'scheme', 'r', 'gamma_c', 'eps_ref', 'k_c', 'coef',
                                         'C_shock', 'role')}, status=a4.status(case))
    if r['status'] == 'not run' or not load_runs(case, 'cfl.dat'):
        r['complete'] = False
        r['kind'] = 'not run'
        return r
    cons = np.concatenate(load_runs(case, 'conservation.dat'))
    bud = np.concatenate(load_runs(case, 'budget.dat'))
    cfl = np.concatenate(load_runs(case, 'cfl.dat'))
    te, t_last = t_end_of(case), float(cons[-1, 0])
    r['complete'] = bool(r['status'] == 'exit 0' and t_last >= te - 0.5*cfl[-1, 1] - 1e-12*te)
    r.update(t_end=te/tau0, t_last=t_last/tau0, steps=len(cfl), dt_min=float(cfl[:, 1].min()),
             dt_max=float(cfl[:, 1].max()))
    if r['complete']:
        r['kind'] = 'complete'
    elif r['status'] == 'incomplete':
        r['kind'] = 'pending'
    elif r['status'] == 'stopped_cost':
        r['kind'] = 'stopped_cost'
    else:
        fl = failure(case)
        r['failure'] = fl
        r['kind'] = 'invalid' if fl['n_invalid'] else ('abort' if fl['n_abort'] else 'other')
        if r['kind'] == 'other':
            r['std_tail'] = open(os.path.join(case, 'calliope.out.std'), errors='replace').read().splitlines()[-5:]
    if not r['complete']:
        rc, rf = decimate(cons[:, (0, 11)]), decimate(cfl[:, :2])
        r['series'] = dict(t=(rc[:, 0]/tau0).tolist(), rho_min=rc[:, 1].tolist(), tf=(rf[:, 0]/tau0).tolist(),
                           dt=rf[:, 1].tolist())
        recs = a4.out3d_records(case)
        if recs:
            t, path, rec = recs[-1]
            g = grid(N)
            sp = spectra(g, a4.out3d_at(g, path, rec))
            r['last_out3d'] = dict(t=t/tau0, E_B=sp['E_B'].tolist(), E_m=sp['E_m'].tolist(),
                                   pileup_B=pileup(sp['E_B'], N/3.0), pileup_m=pileup(sp['E_m'], N/3.0))
    st = ast3.stats(case, dict(dir='x', case=dict(N=N)), [])
    dt = bud[:, 1]
    tot = float(np.sum(dt*bud[:, 7]))
    nu = bud[:, 16]
    r['E'] = dict(drift_M=st['drift_M'], dM_end=float((cons[-1, 2] - cons[0, 2])/cons[0, 2]),
                  R_zeta=st['R_zeta'], R_svv_A=a4.svv_res_rel(case), R_eta=st['R_eta'],
                  err_rel=abs(float(np.sum(dt*bud[:, 8])))/abs(tot) if tot else None,
                  nu_pos_frac=float(np.mean(nu > 0)), nu_pos_int=float(np.sum(dt*np.maximum(nu, 0.0))),
                  nu_int=float(np.sum(dt*nu)),
                  nu_pos_rel=float(np.sum(dt*np.maximum(nu, 0.0)))/abs(tot) if tot else None,
                  rho_min=st['rho_min'], rho_max=st['rho_max'], zeta_max=st['zeta_max'],
                  divb_max=float(cons[:, 19].max()))
    r['share'] = {k: (float(np.sum(dt*sum(np.sign(c)*bud[:, abs(c)] for c in cs)))/tot if tot else None)
                  for k, cs in SHARES.items()}
    r['share_sum'] = sum(r['share'].values()) if tot else None
    r['cost'], r['runs_without_wall'] = cost_to(case)
    r['cost_tstop'] = cost_to(case, g15.T_STOP*tau0)[0]
    res = (('R_zeta', 'R_svv_A', 'R_eta') if m.get('scheme') == 'ES' else
           ('R_zeta', 'R_eta') if m.get('scheme') == 'ESf' else ())
    if res:
        r['E']['es_ok'] = all(r['E'][k] is None or r['E'][k] <= a4.TOL_RES for k in res)
    return r


def cases_of(base, batch):
    return sorted(os.path.dirname(f) for f in glob.glob(os.path.join(base, batch, '*', 'case.json')))


def print_runs(runs):
    print('%-28s %-14s %-12s %6s %8s %9s %9s %9s %9s %9s %8s %7s' % ('case', 'status', 'kind', 't/tau0', 'rho_min',
          'dt_min', 'dM/M0', 'R_max', 'err_rel', 'nu>0', 'share', 'core-h'))
    for name, r in runs.items():
        if r['kind'] == 'not run':
            print('%-28s %-14s %-12s' % (name, r['status'][:14], r['kind']))
            continue
        E = r['E']
        rmax = max(v for v in (E['R_zeta'], E['R_svv_A'], E['R_eta'], 0.0) if v is not None)
        print('%-28s %-14s %-12s %6.3f %8.2e %9.2e %9.2e %9.2e %9.2e %9.3f %8.4f %7.2f' % (
            name, r['status'][:14], r['kind'], r['t_last'], E['rho_min'], r['dt_min'], E['dM_end'], rmax,
            E['err_rel'] or 0.0, E['nu_pos_frac'], r['share_sum'] or 0.0, r['cost']))
        if 'failure' in r and r['failure']['first']:
            f = r['failure']['first']
            print('    first invalid: t/tau0 = %.4f step %d stage %d %s%s = %s (fields %s)' % (
                f['t']/r['info']['tau0'], f['step'], f['stage'], f['field'], f['iyzx'], f['value'],
                ' '.join(r['failure']['fields'])))
        if 'last_out3d' in r:
            L = r['last_out3d']
            print('    last out3d t/tau0 = %.3f, pile-up E_B %.2f E_m %.2f' % (L['t'], L['pileup_B'], L['pileup_m']))
    print('total core-h %.2f' % sum(r.get('cost', 0.0) for r in runs.values()))


def r_star(rungs):
    """The smallest r from which every larger rung is complete ([(r, complete)]), or None; and whether the
    complete rungs are not a top streak (a complete rung below an incomplete one)."""
    rs = sorted(rungs, reverse=True)
    best = None
    for r, ok in rs:
        if not ok:
            break
        best = r
    nonmono = any(ok for r, ok in rs if best is None or r < best)
    return best, nonmono


def batch_runs(base, batch):
    runs = {os.path.basename(c): run_summary(c) for c in cases_of(base, batch)}
    with open(os.path.join(base, batch, 'summary.json'), 'w') as f:
        json.dump(runs, f, indent=1)
    print_runs(runs)
    if batch == 'a':
        out = dict(r_star={}, r0_ok={}, diag=[], rungs={}, nonmono={})
        for prob in g15.PROBS:
            out['r_star'][prob], out['rungs'][prob], out['nonmono'][prob] = {}, {}, {}
            for scheme, rs in (('L', g15.R_L), ('H', g15.R_H)):
                rungs = [(r, runs.get('%s_N64_%s' % (prob, g15.label(scheme, r)), {}).get('complete', False))
                         for r in rs]
                out['r_star'][prob][scheme], out['nonmono'][prob][scheme] = r_star(rungs)
                out['rungs'][prob][scheme] = rungs
                if not any(ok for _, ok in rungs):
                    out['diag'].append([prob, scheme])
            out['r0_ok'][prob] = runs.get('%s_N64_R0' % prob, {}).get('complete', False)
        with open(os.path.join(base, 'a.json'), 'w') as f:
            json.dump(out, f, indent=1)
        print(json.dumps({k: out[k] for k in ('r_star', 'r0_ok', 'diag', 'nonmono')}))
    elif batch == 'b':
        out = dict(r_C=None, r_star_128={}, rungs={})
        for prob in g15.PROBS:
            out['rungs'][prob], out['r_star_128'][prob] = {}, {}
            for scheme in ('L', 'H'):
                rungs = sorted((r['info']['r'], r['complete']) for r in runs.values()
                               if r['info']['problem'] == prob and r['info']['scheme'] == scheme)
                out['rungs'][prob][scheme] = rungs
                out['r_star_128'][prob][scheme] = r_star(rungs)[0]
        for scheme in ('H', 'L'):
            ok = [r for r, c in out['rungs']['turb1'][scheme] if c]
            if ok:
                out['r_C'] = [scheme, min(ok)]
                break
        with open(os.path.join(base, 'b.json'), 'w') as f:
            json.dump(out, f, indent=1)
        print(json.dumps({k: out[k] for k in ('r_C', 'r_star_128')}))
    return runs


def watch(base, batch):
    """The stop rule on the 'incomplete' cases of the batch."""
    for c in cases_of(base, batch):
        st = a4.status(c)
        if st != 'incomplete':
            continue
        m = a4.case_info(c)
        p = projection(c, t_end_of(c))
        p['cap'] = COST_CAP[m['N']]
        p['stop'] = p['projected'] > p['cap']
        print('%-28s spent %7.2f projected %8.2f cap %6.1f t/tau0 %.3f steps left %.0f%s' % (
            os.path.basename(c), p['spent'], p['projected'], p['cap'], p['t_now']/m['tau0'], p['steps_left'],
            '  -> stopped_cost' if p['stop'] else ''))
        if p['stop']:
            with open(os.path.join(c, 'cost_projection.json'), 'w') as f:
                json.dump(p, f, indent=1)
            with open(os.path.join(c, 'exit_status'), 'w') as f:
                f.write('stopped_cost\n')


# ---------------------------------------------------------------- phase 0
def pre(base):
    root = os.path.join(base, 'pre')
    out = dict(p0a={}, p0b={}, eps_ref={}, eps_peak_t={}, es={}, coefs={})
    for lab in ('ES', 'L_r1', 'H_r1'):
        d = {}
        for k in (1, 2):
            c = os.path.join(root, 'dt%d_turb1_N64_%s' % (k, lab))
            bud = load_runs(c, 'budget.dat')
            d['status%d' % k] = a4.status(c)
            if bud:
                bud = np.concatenate(bud)
                d['e%d' % k] = abs(float(np.sum(bud[:, 1]*bud[:, 8])))/abs(float(np.sum(bud[:, 1]*bud[:, 7])))
                d['t%d' % k] = float(bud[-1, 0] + bud[-1, 1])/a4.case_info(c)['tau0']
        if d.get('e2'):
            d['ratio'] = d['e1']/d['e2']
        out['p0a'][lab] = d
    g = grid(64)
    ok = True
    for lab in ('L_r1', 'H_r1'):
        ca, cb = (os.path.join(root, 'dec_turb1_N64_%s_np%d' % (lab, P)) for P in g15.P0B_NP)
        d = dict(status=[a4.status(ca), a4.status(cb)], planes={})
        for pl in a4.TURB_PLANES:
            for (ta, sa, pa, ia), (tb, sb, pb, ib) in zip(a4.plane_records(ca, pl), a4.plane_records(cb, pl)):
                d['planes']['%s %d' % (pl, sb)] = dict(step=(sa, sb), **a4.plane_diff(
                    a4.plane_at(g, pa, ia, pl), a4.plane_at(g, pb, ib, pl), 2*np.pi/64))
        v = d['planes'].values()
        d['max'] = max((x['max'] for x in v), default=None)
        d['ok'] = bool(all(s in ('exit 0', 'dt_changed') for s in d['status']) and v
                       and all(x['step'][0] == x['step'][1] for x in v)
                       and {x['step'][1] for x in v} == {0, g15.P0B_NSTEP} and d['max'] <= a4.TOL_DEC10)
        ok &= d['ok']
        out['p0b'][lab] = d
    for prob in g15.PROBS:
        c = os.path.join(root, '%s_N64_ES' % prob)
        r = run_summary(c)
        out['es'][prob] = r
        ok &= r['complete'] and r.get('E', {}).get('es_ok') is True
        if r['complete']:
            bud = np.concatenate(load_runs(c, 'budget.dat'))
            i = int(np.argmax(-bud[:, 6]))
            out['eps_ref'][prob] = g15.eps_peak_ref(c)
            out['eps_peak_t'][prob] = float(bud[i, 0] + 0.5*bud[i, 1])/r['info']['tau0']
            out['coefs'][prob] = {N: {s: g15.coefs(s, 1.0, N, out['eps_ref'][prob])[3] for s in ('L', 'H')}
                                  for N in (64, 128, 256)}
    out['ok'] = bool(ok)
    with open(os.path.join(base, 'pre.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print('P0a', json.dumps(out['p0a']))
    for lab, d in out['p0b'].items():
        print('P0b %s: status %s, max %s, ok %s' % (lab, d['status'], d['max'], d['ok']))
    for prob, r in out['es'].items():
        print('ES %s: %s complete %s, eps_ref %s at t/tau0 %s, R ok %s' % (
            prob, r['status'], r['complete'], out['eps_ref'].get(prob), out['eps_peak_t'].get(prob),
            r.get('E', {}).get('es_ok')))
    print('ok', out['ok'])
    return out


# ---------------------------------------------------------------- spectra of the runs
def run_spec(case, role, ext=False):
    """Per-record spectra (QS) of the out3d records in the windows and the window means, for the windows the
    run reached (last out3d record at t >= hi - WIN_END - WIN_TOL)."""
    m = a4.case_info(case)
    N, tau0 = m['N'], m['tau0']
    g = grid(N)
    cons = np.concatenate(load_runs(case, 'conservation.dat'))
    t_last = float(cons[-1, 0])/tau0
    recs = a4.out3d_records(case)
    ts = [t/tau0 for t, _, _ in recs]
    t_rec = max(ts, default=-1.0)
    info = dict(problem=m['problem'], N=N, tau0=tau0, scheme='ES' if ext else m.get('scheme'),
                r=None if ext else m.get('r'), role=role, t_last=t_last, t_rec=t_rec)
    info['label'] = g15.label(info['scheme'], info['r'])
    out = dict(info=info, win={})
    for W, (lo, hi) in WIN.items():
        if t_rec < hi - WIN_END - WIN_TOL:
            continue
        ww = window_weights(ts, lo, hi)
        sp = [spectra(g, a4.out3d_at(g, recs[i][1], recs[i][2])) for i, _ in ww]
        out['win'][W] = dict(t=[ts[i] for i, _ in ww], w=[w for _, w in ww],
                             mean={q: sum(w*s[q] for (_, w), s in zip(ww, sp)).tolist() for q in QS},
                             recs={q: [s[q].tolist() for s in sp] for q in QS})
        print('# %s %s: %d records %s' % (case, W, len(ww), ' '.join('%.3f' % ts[i] for i, _ in ww)),
              file=sys.stderr)
    out['cost'] = cost_to(case)[0]
    out['cost_tstop'] = cost_to(case, g15.T_STOP*tau0)[0]
    return out


def spec(base):
    runs = os.path.dirname(os.path.normpath(base))
    out = {}
    for batch in ('pre', 'a', 'b', 'c'):
        for c in cases_of(base, batch):
            name = os.path.basename(c)
            if batch == 'pre' and name not in ['%s_N64_ES' % p for p in g15.PROBS]:
                continue
            if run_summary(c)['complete']:
                out[os.path.relpath(c, runs)] = run_spec(c, batch)
    out[ES256] = run_spec(os.path.join(runs, ES256), 'ES256', ext=True)
    out[S8] = run_spec(os.path.join(runs, S8), 'S8', ext=True)
    with open(os.path.join(base, 'spec.json'), 'w') as f:
        json.dump(out, f)
    for k, v in out.items():
        print('%-40s %-8s %s' % (k, v['info']['label'], ' '.join(sorted(v['win']))))
    return out


# ---------------------------------------------------------------- the judgement
def find(sp, prob, N, label, roles=('pre', 'a', 'b', 'c', 'ES256')):
    for k, v in sp.items():
        i = v['info']
        if i['problem'] == prob and i['N'] == N and i['label'] == label and i['role'] in roles:
            return v
    return None


def pair(run, ref, q, W, cap, d_hi, kc):
    """k_x, kappa, censoring, D and sigma of run against ref (window means; per record of run for sigma)."""
    E, Er = np.array(run['win'][W]['mean'][q]), np.array(ref['win'][W]['mean'][q])
    n = min(len(E), len(Er))
    T = ratio(E[:n], Er[:n])
    k, cens = k_cross(T, cap)
    ks = [k_cross(ratio(np.array(e)[:n], Er[:n]), cap)[0]/kc for e in run['win'][W]['recs'][q]]
    return dict(k=k, kappa=k/kc, cens=cens, D=dev(T, d_hi), sigma=sigma(ks), kappas=ks)


def judge_set(runs, refs, kref, kc, cap, d_hi, extra_delta, need=None):
    """The pairs (QJ x WIN x refs) of every run against the refs; B* among the standard dissipations; delta;
    verdict of ES against B*. runs: {label: spec}; refs: {name: spec}; kref: {(q, W): k_ref or None}."""
    keys = [(q, W, rn) for q in QJ for W in WIN for rn in refs if W in refs[rn]['win']]
    P = {}
    for lab, run in runs.items():
        if all(W in run['win'] for W in WIN):
            P[lab] = {}
            for q, W, rn in keys:
                kr = kref.get((q, W))
                P[lab]['%s %s %s' % (q, W, rn)] = pair(run, refs[rn], q, W, min(cap, kr) if kr else cap,
                                                        int(math.floor(min(d_hi, kr))) if kr else d_hi, kc)
    cand = {lab: float(np.mean([x['kappa'] for x in v.values()])) for lab, v in P.items()
            if lab.split('_')[0] in ('L', 'H', 'R0')}
    out = dict(pairs=P, cand=cand, npairs=len(keys))
    if 'ES' not in P or not cand:
        out['verdict'] = 'no ES or no candidate'
        return out
    B = max(cand, key=cand.get)
    ks = ['%s %s %s' % k for k in keys]
    es, b = P['ES'], P[B]
    delta = [max(2.0/kc, 2*math.hypot(es[k]['sigma'], b[k]['sigma']), extra_delta) for k in ks]
    need = need or need_of(len(ks))
    v, cnt = verdict([es[k]['kappa'] for k in ks], [b[k]['kappa'] for k in ks], [es[k]['D'] for k in ks],
                     [b[k]['D'] for k in ks], delta, need)
    out.update(B=B, delta=dict(zip(ks, delta)), verdict=v, counts=cnt,
               table={k: dict(kES=es[k]['kappa'], kB=b[k]['kappa'], d=es[k]['kappa'] - b[k]['kappa'], delta=dl,
                              DES=es[k]['D'], DB=b[k]['D'], censES=es[k]['cens'], censB=b[k]['cens'])
                      for k, dl in zip(ks, delta)},
               kbar=dict(ES=gmean([es[k]['k'] for k in ks]), B=gmean([b[k]['k'] for k in ks])))
    return out


def judge(base):
    sp = json.load(open(os.path.join(base, 'spec.json')))
    kc128, kc64 = 128/3.0, 64/3.0
    out = dict(flags=[])
    # references and k_ref,W
    es256 = sp[ES256]
    rC = None
    try:
        rC = g15.r_256(os.path.join(base, 'c'))
    except OSError:
        pass
    h256 = find(sp, 'turb1', 256, g15.label(*rC), roles=('c',)) if rC else None
    refs = {'ES256': es256}
    if h256 and all(W in h256['win'] for W in WIN):
        refs['H256'] = h256
    else:
        out['flags'].append('H-256 missing: ES-256 only, need = 3 npairs//4 (consult the user)')
    kref = {}
    if 'H256' in refs:
        for q in QJ:
            for W in WIN:
                T = ratio(es256['win'][W]['mean'][q], h256['win'][W]['mean'][q])
                k, cens = k_cross(T, CAP_REF)
                kref[(q, W)] = k
                out.setdefault('k_ref', {})['%s %s' % (q, W)] = dict(k=k, cens=cens)
                if k < N0 + 1:
                    out['flags'].append('k_ref %s %s = %.2f: the references differ at the large scales' % (q, W, k))
    # M_s = 1 on 128^3
    runs = {v['info']['label']: v for v in sp.values() if v['info']['problem'] == 'turb1' and v['info']['N'] == 128
            and v['info']['role'] == 'b'}
    es128, s8 = runs.get('ES'), sp.get(S8)
    dr = []
    if es128 and s8 and 'W1' in s8['win']:
        for q in QJ:
            for rn, ref in refs.items():
                kr = kref.get((q, 'W1'))
                cap = min(CAP_K, kr) if kr else CAP_K
                a, b = (pair(x, ref, q, 'W1', cap, D_HI, kc128)['kappa'] for x in (es128, s8))
                dr.append(abs(a - b))
    out['delta_rerun'] = max(dr) if dr else None
    if out['delta_rerun'] is None:
        out['flags'].append('Delta_rerun not available')
    J = judge_set(runs, refs, kref, kc128, CAP_K, D_HI, out['delta_rerun'] or 0.0)
    out['ms1'] = J
    if 'B' in J:
        out['cost'] = dict(ES=runs['ES']['cost_tstop'], B=runs[J['B']]['cost_tstop'],
                           ES256=es256['cost_tstop'], H256=h256['cost_tstop'] if h256 else None,
                           **cost_ok(runs['ES']['cost_tstop'], runs[J['B']]['cost_tstop'],
                                     [J['pairs']['ES'][k]['k'] for k in J['table']],
                                     [J['pairs'][J['B']][k]['k'] for k in J['table']]))
    # M_s = 0.5: (64, 128) self-convergence, the reference judgement
    sc05 = {}
    for v in sp.values():
        i = v['info']
        if i['problem'] == 'turb05' and i['N'] == 64:
            w = find(sp, 'turb05', 128, i['label'])
            if w:
                sc05[i['label']] = (v, w)
    J05 = selfconv_set(sc05, kc64, CAP_K64, D_HI64)
    out['ms05'] = J05
    # self-convergence reports
    sc1 = {}
    for v in sp.values():
        i = v['info']
        if i['problem'] == 'turb1' and i['N'] == 64:
            w = find(sp, 'turb1', 128, i['label'])
            if w:
                sc1[i['label']] = (v, w)
    out['self_ms1_64'] = selfconv_set(sc1, kc64, CAP_K64, D_HI64, judge=False)
    sc2 = {}
    if es128:
        sc2['ES'] = (es128, es256)
    if h256 and rC and g15.label(*rC) in runs:
        sc2[g15.label(*rC)] = (runs[g15.label(*rC)], h256)
    out['self_ms1_128'] = selfconv_set(sc2, kc128, CAP_K, D_HI, judge=False)
    with open(os.path.join(base, 'judge.json'), 'w') as f:
        json.dump(out, f, indent=1)
    report(out)
    figures(base, sp, runs, refs, J, sc05, J05)
    return out


def selfconv_set(pairs, kc, cap, d_hi, judge=True):
    """(N, 2N) self-convergence of each label: T = E_N/E_2N with the 2N run as the reference (QJ x WIN); with
    judge, B* among the standard dissipations and the verdict as judge_set with delta = max(2/k_c, 2 sigma)."""
    P = {}
    for lab, (lo, hi) in pairs.items():
        if all(W in lo['win'] and W in hi['win'] for W in WIN):
            P[lab] = {'%s %s' % (q, W): pair(lo, hi, q, W, cap, d_hi, kc) for q in QJ for W in WIN}
    out = dict(pairs=P)
    if not judge:
        return out
    cand = {lab: float(np.mean([x['kappa'] for x in v.values()])) for lab, v in P.items()
            if lab.split('_')[0] in ('L', 'H', 'R0')}
    out['cand'] = cand
    if 'ES' not in P or not cand:
        out['verdict'] = 'no ES or no candidate'
        return out
    B = max(cand, key=cand.get)
    ks = sorted(P['ES'])
    es, b = P['ES'], P[B]
    delta = [max(2.0/kc, 2*math.hypot(es[k]['sigma'], b[k]['sigma'])) for k in ks]
    v, cnt = verdict([es[k]['kappa'] for k in ks], [b[k]['kappa'] for k in ks], [es[k]['D'] for k in ks],
                     [b[k]['D'] for k in ks], delta, need_of(len(ks)))
    out.update(B=B, verdict=v, counts=cnt,
               table={k: dict(kES=es[k]['kappa'], kB=b[k]['kappa'], d=es[k]['kappa'] - b[k]['kappa'], delta=dl,
                              DES=es[k]['D'], DB=b[k]['D'], censES=es[k]['cens'], censB=b[k]['cens'])
                      for k, dl in zip(ks, delta)})
    return out


def report(out):
    for f in out['flags']:
        print('FLAG', f)
    if 'k_ref' in out:
        print('k_ref,W', ', '.join('%s %.2f%s' % (k, v['k'], '+' if v['cens'] else '')
                                   for k, v in out['k_ref'].items()))
    print('Delta_rerun', out['delta_rerun'])
    for key, name in (('ms1', 'M_s = 1, 128^3'), ('ms05', 'M_s = 0.5, (64, 128) reference')):
        J = out[key]
        print('== %s: candidates %s' % (name, ', '.join('%s %.3f' % x for x in sorted(J.get('cand', {}).items()))))
        if 'B' not in J:
            print(J['verdict'])
            continue
        print('B* = %s: %s %s' % (J['B'], J['verdict'], J['counts']))
        for k, t in J['table'].items():
            print('  %-20s kES %.3f%s kB %.3f%s d %+.3f delta %.3f  D_ES %.3f D_B %.3f' % (
                k, t['kES'], '+' if t['censES'] else ' ', t['kB'], '+' if t['censB'] else ' ', t['d'], t['delta'],
                t['DES'], t['DB']))
    if 'cost' in out:
        print('cost', json.dumps(out['cost']))
    for key in ('self_ms1_64', 'self_ms1_128'):
        for lab, v in out[key]['pairs'].items():
            print('%s %-10s %s' % (key, lab, ' '.join('%s %.3f%s/%.3f' % (k, x['kappa'], '+' if x['cens'] else '',
                                                                         x['D']) for k, x in v.items())))


COLORS = {'ES': 'C0', 'B': 'C1', 'ESf': 'C2', 'ES256': 'k', 'H256': '0.5'}


def figures(base, sp, runs, refs, J, sc05, J05):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    lines = [('ES', runs.get('ES'), COLORS['ES'], '-')]
    if 'B' in J:
        lines.append(('B* = ' + J['B'], runs[J['B']], COLORS['B'], '-'))
    if 'ESf' in runs:
        lines.append(('ESf', runs['ESf'], COLORS['ESf'], '-'))
    lines += [(rn + ' (256^3)', refs[rn], COLORS[rn], '--') for rn in refs]
    fig, ax = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    for i, q in enumerate(QJ):
        for j, W in enumerate(WIN):
            for name, run, col, ls in lines:
                if run is None or W not in run['win']:
                    continue
                E = np.array(run['win'][W]['mean'][q])
                n = np.arange(len(E))
                ax[i, j].loglog(n[1:], E[1:], color=col, ls=ls, lw=1.2, label=name)
            ax[i, j].axvline(128/3.0, color='0.7', lw=0.8)
            ax[i, j].set_title('%s, %s = [%g, %g] tau0' % (q, W, *WIN[W]))
            ax[i, j].set_xlabel('n')
    ax[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(base, 'fig_spec_ms1.png'), dpi=120)
    plt.close(fig)
    fig, ax = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
    for i, q in enumerate(QJ):
        for j, W in enumerate(WIN):
            for name, run, col, ls in lines:
                if run is None or W not in run['win'] or run['info']['N'] != 128:
                    continue
                for rn, ref, rls in zip(refs, refs.values(), ('-', ':')):
                    E, Er = np.array(run['win'][W]['mean'][q]), np.array(ref['win'][W]['mean'][q])
                    n = min(len(E), len(Er))
                    ax[i, j].semilogx(np.arange(1, n), np.log(ratio(E[1:n], Er[1:n])), color=col, ls=rls, lw=1.2,
                                      label='%s / %s' % (name, rn))
            for s in (1, -1):
                ax[i, j].axhline(s*math.log(T_TOL), color='0.6', lw=0.8)
            ax[i, j].axvline(128/3.0, color='0.7', lw=0.8)
            ax[i, j].set_ylim(-1.0, 1.0)
            ax[i, j].set_title('ln T, %s, %s' % (q, W))
            ax[i, j].set_xlabel('n')
    ax[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(base, 'fig_T_ms1.png'), dpi=120)
    plt.close(fig)
    if 'B' in J05:
        fig, ax = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
        for i, q in enumerate(QJ):
            for j, W in enumerate(WIN):
                for lab, col in (('ES', COLORS['ES']), (J05['B'], COLORS['B'])):
                    lo, hi = sc05[lab]
                    E, Er = np.array(lo['win'][W]['mean'][q]), np.array(hi['win'][W]['mean'][q])
                    n = min(len(E), len(Er))
                    ax[i, j].semilogx(np.arange(1, n), np.log(ratio(E[1:n], Er[1:n])), color=col, lw=1.2,
                                      label='%s: 64 / 128' % lab)
                for s in (1, -1):
                    ax[i, j].axhline(s*math.log(T_TOL), color='0.6', lw=0.8)
                ax[i, j].axvline(64/3.0, color='0.7', lw=0.8)
                ax[i, j].set_ylim(-1.0, 1.0)
                ax[i, j].set_title('ln T, M_s = 0.5, %s, %s' % (q, W))
                ax[i, j].set_xlabel('n')
        ax[0, 0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(base, 'fig_T_ms05.png'), dpi=120)
        plt.close(fig)


# ---------------------------------------------------------------- selftest
def selftest():
    import tempfile
    bad = []

    def check(name, ok):
        print('%-60s %s' % (name, 'PASS' if ok else 'FAIL'))
        if not ok:
            bad.append(name)

    n = np.arange(60, dtype=float)
    thr = math.log(T_TOL)
    # k_x: ln T linear in n crosses ln 1.25 at 3 + thr/a exactly
    a = 0.02
    T = np.exp(a*(n - 3))
    check('k_cross linear', abs(k_cross(T, CAP_K)[0] - (3 + thr/a)) < 1e-12 and not k_cross(T, CAP_K)[1])
    check('k_cross negative ln T', abs(k_cross(1/T, CAP_K)[0] - (3 + thr/a)) < 1e-12)
    check('k_cross censored', k_cross(np.ones(60), CAP_K) == (CAP_K, True))
    check('k_cross censored at non-integer cap', k_cross(np.ones(60), 30.5) == (30.5, True))
    check('k_cross exceeds at n0', k_cross(np.full(60, 2.0), CAP_K) == (3.0, False))
    T2 = np.ones(60)
    T2[10] = np.nan
    check('k_cross non-finite', k_cross(T2, CAP_K) == (9.0, False))
    T2[10], T2[11] = 0.0, 1.0
    check('k_cross zero', k_cross(T2, CAP_K) == (9.0, False))
    T3 = np.ones(60)
    T3[40] = 2.0
    check('k_cross beyond cap', k_cross(T3, 39.9) == (39.9, True) and k_cross(T3, 40.0)[0] > 39.0)
    # D
    check('D constant', abs(dev(np.full(60, math.exp(-0.1)), 32) - 0.1) < 1e-12)
    Tl = np.exp(np.log(np.maximum(n, 1)))    # ln T = ln n
    x = np.log(np.arange(3, 33))
    want = np.sum(0.5*(x[1:] + x[:-1])*np.diff(x))/(x[-1] - x[0])
    check('D ln n (trapezoid)', abs(dev(Tl, 32) - want) < 1e-12 and abs(want - 0.5*(x[0] + x[-1])) < 1e-12)
    check('D single shell', abs(dev(np.full(60, 1.5), 3) - math.log(1.5)) < 1e-15)
    check('D non-finite', dev(T2, 32) == float('inf') and math.isnan(dev(np.ones(60), 2.9)))
    # sigma
    check('sigma ddof 1', abs(sigma([1.0, 2.0, 3.0]) - 1.0) < 1e-15 and sigma([0.4]) == 0.0)
    # verdict
    d8 = [0.05]*8
    check('verdict better', verdict([0.8]*8, [0.6]*8, [0.1]*8, [0.2]*8, d8, 6)[0] == 'ES better')
    check('verdict worse', verdict([0.6]*8, [0.8]*8, [0.2]*8, [0.1]*8, d8, 6)[0] == 'ES worse')
    check('verdict 5 of 8 undecided', verdict([0.8]*5 + [0.6]*3, [0.6]*8, [0.1]*8, [0.2]*8, d8, 6)[0] == 'undecided')
    check('verdict one minus undecided',
          verdict([0.8]*7 + [0.5], [0.6]*8, [0.1]*8, [0.2]*8, d8, 6)[0] == 'undecided')
    check('verdict D fails undecided', verdict([0.8]*8, [0.6]*8, [0.3]*8, [0.2]*8, d8, 6)[0] == 'undecided')
    check('verdict 6 of 8 within delta', verdict([0.8]*6 + [0.62]*2, [0.6]*8, [0.1]*8, [0.2]*8, d8, 6)[0]
          == 'ES better')
    check('verdict per-pair delta', verdict([0.8]*8, [0.6]*8, [0.1]*8, [0.2]*8, [0.3]*3 + [0.05]*5, 6)[0]
          == 'undecided')
    check('verdict 3 of 4', verdict([0.8]*3 + [0.6], [0.6]*4, [0.1]*4, [0.2]*4, d8[:4], need_of(4))[0]
          == 'ES better' and need_of(8) == 6 and need_of(4) == 3 and need_of(6) == 4)
    # window weights
    ww = window_weights([0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375], 0.75, 1.25)
    check('window trapezoid', [i for i, _ in ww] == [2, 3, 4, 5, 6]
          and np.allclose([w for _, w in ww], [1/8, 1/4, 1/4, 1/4, 1/8], rtol=0, atol=1e-15))
    e256 = [0.5, 0.818, 0.893, 0.906, 0.921, 0.965, 1.0, 1.0, 1.020, 1.156, 1.5, 1.503, 1.763, 2.0]
    ww = window_weights(e256, 0.75, 1.25)
    t = np.array([e256[i] for i, _ in ww])
    edge = np.concatenate(([0.75], 0.5*(t[1:] + t[:-1]), [1.25]))
    check('window ES-256 W1', [i for i, _ in ww] == [1, 2, 3, 4, 5, 6, 8, 9]
          and np.allclose([w for _, w in ww], np.diff(edge)/0.5, rtol=0, atol=1e-15)
          and abs(sum(w for _, w in ww) - 1) < 1e-15)
    ww = window_weights(e256, 1.5, 2.0)
    check('window ES-256 W2', [i for i, _ in ww] == [10, 11, 12, 13] and abs(sum(w for _, w in ww) - 1) < 1e-15
          and abs(ww[0][1] - 0.0015/0.5) < 1e-12)
    ww = window_weights([0.7395, 0.7405, 1.255, 1.255*(1 + 1e-12), 1.265], 0.75, 1.25)
    check('window tolerance and duplicates', [i for i, _ in ww] == [1, 2] and abs(sum(w for _, w in ww) - 1) < 1e-15)
    # cost rule
    c = cost_ok(10.0, 5.0, [30.0]*8, [25.0]*8)
    check('cost rule', c['ok'] == ((30/25)**4 > 2.0) and abs(c['k_ratio4'] - 1.2**4) < 1e-12)
    check('cost rule fails', not cost_ok(10.0, 2.0, [30.0]*8, [25.0]*8)['ok'])
    check('gmean', abs(gmean([2.0, 8.0]) - 4.0) < 1e-14)
    # r_star
    check('r_star', r_star([(0.25, False), (0.5, True), (1.0, True), (2.0, True)]) == (0.5, False)
          and r_star([(0.25, True), (0.5, False), (1.0, True)]) == (1.0, True)
          and r_star([(1.0, False), (2.0, False)]) == (None, False)
          and r_star([(1.0, True), (2.0, False)]) == (None, True))
    # projection, cost_to and failure on a synthetic case
    with tempfile.TemporaryDirectory() as d:
        c = os.path.join(d, 'case')
        os.makedirs(os.path.join(c, 'chunk_1'))
        json.dump(dict(N=64, tau0=math.pi, nrun=0, dt=1.0), open(os.path.join(c, 'case.json'), 'w'))
        open(os.path.join(c, 'nproc'), 'w').write('64\n')
        open(os.path.join(c, 't_end'), 'w').write('%.17g\n' % (2*math.pi))
        dt1, dt2 = 1e-3, 5e-4
        t1 = np.arange(1000)*dt1
        t2 = t1[-1] + dt1 + np.arange(400)*dt2
        for run, t, dt, wall in (('chunk_1', t1, dt1, 1000.0), ('', t2, dt2, 600.0)):
            cfl = np.zeros((len(t), 17))
            cfl[:, 0], cfl[:, 1] = t, dt
            np.savetxt(os.path.join(c, run, 'cfl.dat'), cfl)
            std = 'wall %.1f s\n' % wall
            if run:
                std += ('[proc      3] invalid state at t =  1.2000000000000000E+00, step       1201, stage  2: '
                        'my(iy, iz, ix) = (       5,       6,       7) =                      NaN\n'
                        '[proc      1] invalid state at t =  1.2000000000000000E+00, step       1201, stage  2: '
                        'lnrho(iy, iz, ix) = (       1,       2,       3) =  -2.1000000000000000E+01\n'
                        '[proc      1] ABORT: non-finite value or |ln rho| > s_abort\n')
            open(os.path.join(c, run, 'calliope.out.std'), 'w').write(std)
        p = projection(c, 2*math.pi)
        spent = (1000.0 + 600.0)*64/3600
        t_now = t2[-1] + dt2
        want = spent + (2*math.pi - t_now)/dt2*(600.0/400)*64/3600
        check('projection', abs(p['spent'] - spent) < 1e-12 and abs(p['projected'] - want) < 1e-9*want
              and p['runs_without_wall'] == 0)
        T = 0.5
        f1 = float(np.mean(t1 + 0.5*dt1 < T))
        check('cost_to T', abs(cost_to(c, T)[0] - 1000.0*64/3600*f1) < 1e-12 and 0 < f1 < 1)
        fl = failure(c)
        check('failure parse', fl['n_invalid'] == 2 and fl['first']['proc'] == 1 and fl['first']['field'] == 'lnrho'
              and fl['first']['iyzx'] == (1, 2, 3) and fl['fields'] == ['lnrho', 'my']
              and fl['abort'] == ['non-finite value or |ln rho| > s_abort'] and fl['first']['run'] == 'chunk_1')
    # spectra against analyze_stage4.shell_spectra_3d
    g = gen_ic.Grid(16, 16, 16)
    rng = np.random.default_rng(1)
    fr = dict(rho=1.0 + 0.3*np.sin(g.X + 2*g.Y)*np.cos(g.Z), **{q: rng.standard_normal(g.X.shape) for q in
                                                                 ('mx', 'my', 'mz', 'bx', 'by', 'bz')})
    fk = {q: g.fft(v)*(g.mask if q != 'rho' else 1.0) for q, v in fr.items()}
    s, s4 = spectra(g, fk), a4.shell_spectra_3d(g, fk)[0]
    check('spectra = shell_spectra_3d (E_B, E_sru, P_lnrho)',
          all(np.allclose(s[q], s4[p], rtol=1e-13, atol=0) for q, p in (('E_B', 'E_B'), ('E_sru', 'E_sru'),
                                                                         ('P_s', 'P_lnrho'))))
    K = np.sqrt(g.kx**2 + g.ky**2 + g.kz**2)
    m = [g.ifft(fk[q]*(K < 6.5)) for q in ('mx', 'my', 'mz')]
    sb = spectra(g, dict(fk, **{q: g.fft(v) for q, v in zip(('mx', 'my', 'mz'), m)}))
    check('E_m Parseval (|k| < 6.5)', abs(np.sum(sb['E_m']) - 0.5*np.mean(sum(v**2 for v in m)))
          < 1e-12*np.sum(sb['E_m']) and np.all(sb['E_m'][7:] < 1e-28))
    E = np.ones(17)
    E[4] = 0.5
    check('pileup', abs(pileup(E, 16/3.0) - 2.0) < 1e-15 and pileup(np.ones(17), 16/3.0) == 1.0)
    print('selftest: %s' % ('all PASS' if not bad else 'FAIL: ' + ', '.join(bad)))
    return not bad


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'selftest':
        sys.exit(0 if selftest() else 1)
    if cmd in ('runs', 'watch'):
        base = sys.argv[3] if len(sys.argv) > 3 else os.path.join('runs', 'stage15')
        (batch_runs if cmd == 'runs' else watch)(base, sys.argv[2])
    else:
        base = sys.argv[2] if len(sys.argv) > 2 else os.path.join('runs', 'stage15')
        dict(pre=pre, spec=spec, judge=judge)[cmd](base)
