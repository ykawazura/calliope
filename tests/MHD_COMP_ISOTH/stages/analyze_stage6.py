#!/usr/bin/env python3
"""Stage 6 analysis of the runs of gen_stage6.py (reset_method = 'track' of src/time.F90, turb1 on 256^3).

    python3 analyze_stage6.py <batch> [base]   reads <base>/<batch>/ (default base: runs/stage6), prints the
                                               checks and writes <base>/<batch>/summary.json (and PNG files)

The runs of Stage 4 are read from runs/stage4 next to base. The criteria were fixed before the runs
(tasks/todo.md, Stage 6).

Batches
    regress   ISOTH pairs <case>_old, <case>_new as analyze_stage4.regress_pair (dec_turb1_N64 must also change
              dt, or it would not test 'decrement', and run nstep steps); MHD_INCOMP pairs (incomp_*): every file
              except exe byte for byte, *.std without the lines of the date and the timers (STD_VOLATILE) and with
              the lines of the FFTW shim of all ranks (SHIM_LINE) as a multiset, calliope_run.in without
              max_wall_time; equal statuses
    track     decomposition: dec_N64_np1, np8 against dec_N64_np64 (TDEC_*; each with the status dt_changed of the
              job script and nstep steps); turb1_N64_t, turb1_N128_t against turb1_N64, turb1_N128_a of
              runs/stage4/turb (TRACK_*); dt/dt_cfl in fig_track_dt.png
    turb      turb1 on 64^3, 128^3 (runs/stage4/turb) and 256^3 (turb1_N256): analyze_stage4.turb_convergence
              (S4.3; on [0, t_c] by turb_convergence_upto when the chain of N256 was cut at t_c < t_end) and
              analyze_stage4.turb_figures in <base>/turb, dt of 'track' and the cost per step
"""
import os
import re
import sys
import json
import glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_stage4 as a4
from analyze_stage1 import load, load_runs, read_nml, num, strv, run_dirs

TDEC_REF = 'dec_N64_np64'
TDEC_P = (1, 8, 64)
TDEC_TOL_DT = 0.0      # max|dt(P) - dt(64)| over the rows/max dt(64): bitwise, as dt of 'track' is rounded to
                       # track_digits digits (runs/stage6, without the rounding, gave 3.1e-9 and 4.0e-9)
TDEC_TOL_CONS = 1e-10  # the same for E_K, E_M, rho_min of conservation.dat
TDEC_COLS = dict(E_K=7, E_M=8, rho_min=11)
TRACK_RUNS = {64: ('turb1_N64_t', 'turb1_N64'), 128: ('turb1_N128_t', 'turb1_N128_a')}  # 'track', Stage 4
TRACK_PRED = {64: 5019, 128: 12194}  # steps of the simulation of 'track' (3 digits) on the cfl.dat of the Stage 4
                                     # run (without the rounding: 5055, 12251, runs/stage6 5024, 12186)
TRACK_TAU = (1, 2, 3)  # t/tau0 of the state
TRACK_COLS = dict(E_K=7, E_M=8, rho_min=11, rho_max=12)
TRACK_RES_FRAC = 0.1   # d_track <= max(TRACK_RES_FRAC d_res, TRACK_FLOOR); d_res: 64^3 against 128^3 of Stage 4
TRACK_FLOOR = 1e-3
RK3_LIM = 2.51         # check steps: max(dt_prev, dt) lam_svv (above it the run stops, advance.F90)
COST_TIMERS = ('get_svv_es', 'svv_power', 'stage_budget')  # ISOTH timers whose share chunk_costs reports
N256 = 'turb1_N256'
N256_CUT = 'incomplete'  # status of the chain of N256 stopped at the cap of 20 jobs before t_end (agreed with the
                         # user, tasks/todo.md): judged on [0, t_c], t_c its last t (turb_convergence_upto)
TURB_PEAK_DROP = a4.OT_CONV_EPS  # a run cut at t_c has passed the peak of eps when eps_w at the end of [0, t_c]
                                 # is <= (1 - drop) of its maximum on [0, t_c] (fixed when the cut became likely,
                                 # before t_c was known)
TURB_EPS_WIN = 0.02  # [tau0] eps_w(t) = (E(t) - E(t + w))/w, w = TURB_EPS_WIN tau0, on t = 0, w/4, ...: eps per step
                     # of N256 varies by 25% within 0.05 while rho_min ~ 0.01 (26.4-33.4 on t = 2.80-2.85)
STD_VOLATILE = re.compile(r'executed on|wall time exceeds|\d min +[\d.]+ %|total from timer is|init \+ adv|^wall ')
RUN_IN_VOLATILE = re.compile(r'^\s*max_wall_time\s*=')
SHIM_LINE = re.compile(r'^fftw_shim: ')  # written by every rank to stderr: its place among the lines of rank 0 varies


def base4_of(root):
    return os.path.join(os.path.dirname(os.path.dirname(os.path.normpath(root))), 'stage4')


def dt_history(case):
    """cfl.dat of each run of the case (chunk_1, ..., the case directory), dt0, reset_method and svv_pi_check
    (default 100) of calliope.in, and the changes of dt: (run, row, dt before, dt after) where dt of a row
    differs from that of the row before (row 0 of the first run: from dt0; of a continued run: from the last row
    of the run before, whose dt the job script passes on). Row j of a run is its step istep_adv = j; the power
    method converges at the rows j > 0 that are multiples of svv_pi_check (check steps)."""
    nml = read_nml(os.path.join(case, 'calliope.in'))
    dt0 = num(nml['dt'])
    cfls = load_runs(case, 'cfl.dat')
    prev, ch = dt0, []
    for k, c in enumerate(cfls):
        dt = np.concatenate(([prev], c[:, 1]))
        for j in np.nonzero(dt[1:] != dt[:-1])[0]:
            ch.append((k, int(j), float(dt[j]), float(dt[j + 1])))
        prev = c[-1, 1]
    return dict(cfls=cfls, dt0=dt0, reset=strv(nml.get('reset_method', "'default'")),
                pi_check=int(num(nml.get('svv_pi_check', '100'))), changes=ch)


def dt_stats(h, tau0):
    """Steps, changes of dt (raises only at check steps), max dt/dt_cfl over the rows, max(dt_prev, dt) lam_svv
    over the check steps, and the mean of dt/dt_cfl over the rows in [j, j + 1) tau0 (j = 0, 1, 2)."""
    rows, lam = [], 0.0
    for c in h['cfls']:
        dt = c[:, 1]
        j = np.arange(len(c))
        chk = (j > 0) & (j % h['pi_check'] == 0)
        if chk.any():
            lam = max(lam, float(np.max(np.maximum(dt[np.nonzero(chk)[0] - 1], dt[chk])*c[chk, 12])))
        rows.append(np.column_stack((c[:, 0], dt, np.minimum(c[:, 2], c[:, 3]))))
    a = np.vstack(rows)
    r = a[:, 1]/a[:, 2]
    up = [x for x in h['changes'] if x[3] > x[2]]
    mean = {}
    for j in range(3):
        m = (a[:, 0] >= j*tau0) & (a[:, 0] < (j + 1)*tau0)
        mean['%d-%d' % (j, j + 1)] = float(np.mean(r[m])) if m.any() else None
    return dict(reset=h['reset'], steps=len(a), changes=len(h['changes']), raises=len(up),
                decreases=len(h['changes']) - len(up), raise_rows=sorted({x[1] for x in up}),
                raise_at_check=all(x[1] > 0 and x[1] % h['pi_check'] == 0 for x in up),
                max_dt_over_cfl=float(r.max()), max_dt_lam_check=lam, mean_ratio=mean,
                dt_min=float(a[:, 1].min()), dt_max=float(a[:, 1].max()),
                finite=bool(all(np.all(np.isfinite(c)) for c in h['cfls']))), a


def steps_advanced(case):
    """'# of steps advanced' of calliope.out.std (written at the normal end of a run), or None."""
    m = re.search(r'# of steps advanced\s+(\d+)', open(os.path.join(case, 'calliope.out.std'), errors='replace').read())
    return int(m.group(1)) if m else None


def ran_nstep(case, steps):
    """A run without the file t_end whose dt changed: the job script gives it the status 'dt_changed' (without the
    exit code of mpirun), so the end of the run is checked by the steps: nstep of calliope.in in calliope.out.std
    and in cfl.dat (steps)."""
    n = int(num(read_nml(os.path.join(case, 'calliope.in'))['nstep']))
    return a4.status(case) == 'dt_changed' and steps_advanced(case) == n == steps


def rel_max(a, b):
    """max|a - b|/max|b|."""
    return float(np.max(np.abs(a - b))/max(float(np.max(np.abs(b))), 1e-300))


# ------------------------------------------------------------ regress
def files_of(case):
    return {os.path.relpath(os.path.join(d, f), case) for d, _, fs in os.walk(case) for f in fs}


def lines_of(p):
    return open(p, 'rb').read().decode(errors='replace').splitlines()


def same_lines(a, b, drop):
    def kept(p):
        return [ln for ln in lines_of(p) if not drop.search(ln)]
    return kept(a) == kept(b)


def same_std(a, b):
    """*.std: the lines without STD_VOLATILE and SHIM_LINE in order, the lines SHIM_LINE as a multiset."""
    def shim(p):
        return sorted(ln for ln in lines_of(p) if SHIM_LINE.search(ln))
    return same_lines(a, b, re.compile(STD_VOLATILE.pattern + '|' + SHIM_LINE.pattern)) and shim(a) == shim(b)


def regress_all_files(old, new):
    """Every file of the runs old and new except exe: *.std as same_std, calliope_run.in without RUN_IN_VOLATILE,
    the others byte for byte; the file sets equal."""
    fo, fn = files_of(old) - {'exe'}, files_of(new) - {'exe'}
    d = {}
    for f in sorted(fo & fn):
        a, b = os.path.join(old, f), os.path.join(new, f)
        if f.endswith('.std'):
            d[f] = same_std(a, b)
        elif os.path.basename(f) == 'calliope_run.in':
            d[f] = same_lines(a, b, RUN_IN_VOLATILE)
        else:
            d[f] = a4.same_bytes(a, b)
    r = dict(status=(a4.status(old), a4.status(new)), same=d, only_old=sorted(fo - fn), only_new=sorted(fn - fo))
    r['ok'] = r['status'][0] == r['status'][1] and all(d.values()) and not r['only_old'] and not r['only_new']
    return r


def regress(root):
    out, ok_all = {}, True
    for old in sorted(glob.glob(os.path.join(root, '*_old'))):
        name = os.path.basename(old)[:-4]
        new = old[:-4] + '_new'
        if name.startswith('incomp'):
            r = regress_all_files(old, new)
        else:
            r = a4.regress_pair(old, new)
        if name.startswith('dec_'):
            # the pair tests 'decrement' only if reset_dt changed dt
            hs = [dt_history(c) for c in (old, new)]
            r['dt_changes'] = [len(h['changes']) for h in hs]
            r['ran_nstep'] = [ran_nstep(c, sum(len(x) for x in h['cfls'])) for c, h in zip((old, new), hs)]
            r['ok'] &= min(r['dt_changes']) > 0 and all(r['ran_nstep'])
        out[name] = r
        ok_all &= r['ok']
        bad = [k for k, v in r['same'].items() if v is False] + r.get('only_old', []) + r.get('only_new', [])
        print('%-28s %-24s ok=%s  files %d  differ: %s%s' % (name, '/'.join(r['status']), r['ok'], len(r['same']),
                                                           bad or '-', '  dt changes %s' % r['dt_changes']
                                                           if 'dt_changes' in r else ''))
    ok_all &= bool(out)
    print('regress:', 'PASS' if ok_all else 'FAIL')
    return dict(cases=out, ok=ok_all)


# ------------------------------------------------------------ track
def track_decomp(root):
    """dec_N64_np<P> against TDEC_REF: rows, the rows where dt changes, max|d dt|/max dt (TDEC_TOL_DT) and
    E_K, E_M, rho_min of conservation.dat (TDEC_TOL_CONS); each P: dt raised at least once, only at check steps."""
    out = {}
    ref = os.path.join(root, TDEC_REF)
    hb = dt_history(ref)
    dtb = np.concatenate([c[:, 1] for c in hb['cfls']])
    cb = np.concatenate(load_runs(ref, 'conservation.dat'))
    for P in TDEC_P:
        case = os.path.join(root, 'dec_N64_np%d' % P)
        h = dt_history(case)
        s, _ = dt_stats(h, a4.case_info(case)['tau0'])
        d = dict(status=a4.status(case), stats=s)
        d['ok'] = ran_nstep(case, s['steps']) and s['finite'] and s['raises'] > 0 and s['raise_at_check']
        if case != ref:
            dta = np.concatenate([c[:, 1] for c in h['cfls']])
            ca = np.concatenate(load_runs(case, 'conservation.dat'))
            d['rows'] = (len(dta), len(dtb))
            d['same_change_rows'] = [x[:2] for x in h['changes']] == [x[:2] for x in hb['changes']]
            if len(dta) == len(dtb) and ca.shape == cb.shape:
                d['dt'] = rel_max(dta, dtb)
                d['cons'] = {q: rel_max(ca[:, c], cb[:, c]) for q, c in TDEC_COLS.items()}
                d['ok'] &= (d['same_change_rows'] and d['dt'] <= TDEC_TOL_DT
                            and max(d['cons'].values()) <= TDEC_TOL_CONS)
            else:
                d['ok'] = False
        out['np%d' % P] = d
        print('decomp P = %-3d %-8s ok=%s  steps %d  changes %d (raises at rows %s)  dt %s  cons %s'
              % (P, d['status'], d['ok'], s['steps'], s['changes'], s['raise_rows'],
                 '%.1e' % d['dt'] if 'dt' in d else '-',
                 {q: '%.1e' % v for q, v in d['cons'].items()} if 'cons' in d else '-'))
    return out


def track_state(case, m, r):
    """eps peak and its time (analyze_stage4.turb_run r) and E_K, E_M, rho_min, rho_max at t = j tau0
    (j in TRACK_TAU) from conservation.dat of every step, interpolated in t."""
    cons = np.concatenate(load_runs(case, 'conservation.dat'))
    t, i = np.unique(cons[:, 0], return_index=True)
    q = dict(eps_peak=r['eps_peak']['eps'], t_peak=r['eps_peak']['t'])
    for j in TRACK_TAU:
        for name, c in TRACK_COLS.items():
            q['%s %g' % (name, j)] = float(np.interp(j*m['tau0'], t, cons[i, c]))
    return q


def track_figure(root, rows):
    """dt/dt_cfl against t/tau0 of the Stage 4 run and the 'track' run, one panel per N, as fig_track_dt.png."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, len(rows), figsize=(4.6*len(rows), 3.4), squeeze=False, sharey=True)
    for a, (N, runs) in zip(ax[0], sorted(rows.items())):
        for (label, tau0, arr, s), c in zip(runs, ('C0', 'C1')):
            a.plot(arr[:, 0]/tau0, arr[:, 1]/arr[:, 2], color=c, lw=1.0,
                   label="%s, %s (%d steps)" % (label, s['reset'], s['steps']))
        for y in (0.8, 0.7):
            a.axhline(y, color='0.5', lw=0.6, ls=':')
        a.set_title('turb1, N = %d (dotted: track_safety, track_raise)' % N, fontsize=9)
        a.set_xlabel('t/tau0')
        a.set_ylim(0, 1.05)
        a.legend(frameon=False, fontsize=8, loc='lower left')
        a.grid(alpha=0.3, lw=0.5)
    ax[0, 0].set_ylabel('dt/dt_cfl')
    fig.tight_layout()
    fig.savefig(os.path.join(root, 'fig_track_dt.png'), dpi=110)
    plt.close(fig)


def track(root):
    base4 = os.path.join(base4_of(root), 'turb')
    dec = track_decomp(root)
    runs, q, dts, fig = {}, {}, {}, {}
    for N, names in sorted(TRACK_RUNS.items()):
        fig[N] = []
        for name, case in ((names[1], os.path.join(base4, names[1])), (names[0], os.path.join(root, names[0]))):
            m = a4.case_info(case)
            r = a4.turb_run(case, m)
            r['R'] = a4.turb_residual(r)
            r['status'] = a4.status(case)
            runs[name], q[name] = r, track_state(case, m, r)
            dts[name], arr = dt_stats(dt_history(case), m['tau0'])
            fig[N].append(('stage 4' if case.startswith(base4) else 'stage 6', m['tau0'], arr, dts[name]))
            a4.turb_print(name, r['status'], r)
    ref64, ref128 = q[TRACK_RUNS[64][1]], q[TRACK_RUNS[128][1]]
    d_res = {k: a4.rel(ref64[k], ref128[k]) for k in ref128}
    lim = {k: max(TRACK_RES_FRAC*v, TRACK_FLOOR) for k, v in d_res.items()}
    checks = {'decomp ' + k: dict(ok=v['ok']) for k, v in dec.items()}
    valid = {}
    for N, (tname, rname) in sorted(TRACK_RUNS.items()):
        d = {k: a4.rel(q[tname][k], q[rname][k]) for k in q[rname]}
        s, r = dts[tname], runs[tname]
        valid[N] = dict(d_track=d, bad=[k for k in d if d[k] > lim[k]], steps=(s['steps'], dts[rname]['steps']),
                        pred=TRACK_PRED[N], ratio_to_ref=s['steps']/dts[rname]['steps'])
        ok_run = (r['status'] == 'exit 0' and r['complete'] and s['finite'] and s['max_dt_over_cfl'] <= 1.0
                  and s['max_dt_lam_check'] <= RK3_LIM and s['raise_at_check'])
        checks['track N%d' % N] = dict(ok=bool(ok_run and not valid[N]['bad']), run=ok_run, bad=valid[N]['bad'])
    for name, s in dts.items():
        print('%-14s %-10s steps %6d  changes %3d (raises %3d, at check steps %s)  max dt/dt_cfl %.4f  '
              'max dt*lam (check) %.3f  dt in [%.3e, %.3e]  mean dt/dt_cfl %s'
              % (name, s['reset'], s['steps'], s['changes'], s['raises'], s['raise_at_check'], s['max_dt_over_cfl'],
                 s['max_dt_lam_check'], s['dt_min'], s['dt_max'],
                 {k: (round(v, 3) if v is not None else None) for k, v in s['mean_ratio'].items()}))
    print('%-12s %10s %10s %10s %10s' % ('quantity', 'd_res', 'limit', 'd_track 64', 'd_track 128'))
    for k in d_res:
        print('%-12s %10.2e %10.2e %10.2e %10.2e' % (k, d_res[k], lim[k], valid[64]['d_track'][k],
                                                   valid[128]['d_track'][k]))
    for N, v in sorted(valid.items()):
        print('N = %d: steps %d (Stage 4 %d, ratio %.3f; predicted %d)  exceeds the limit: %s'
              % (N, v['steps'][0], v['steps'][1], v['ratio_to_ref'], v['pred'], v['bad'] or '-'))
    track_figure(root, fig)
    ok_all = all(c['ok'] for c in checks.values())
    for k, c in checks.items():
        if not c['ok']:
            print('%-20s FAIL %s' % (k, c))
    print('track:', 'PASS' if ok_all else 'FAIL')
    return dict(decomp=dec, dt=dts, state=q, d_res=d_res, limit=lim, valid=valid, checks=checks,
                runs={k: {x: r[x] for x in ('info', 'stats', 'share', 'end', 'eps_peak', 'dt_lim', 'complete',
                                            't_last', 't_end', 'R', 'status')} for k, r in runs.items()}, ok=ok_all)


# ------------------------------------------------------------ turb
def chunk_costs(case):
    """Each run of the chain: t of its first and last row and steps (rows of cfl.dat), the mean number of power
    iterations over the rows with more than one (the check steps; column 16 of cfl.dat), the timer 'step (solve)'
    (max over the ranks), the shares of COST_TIMERS in it [%] and the wall time of mpirun [s] from its .std."""
    out = []
    for d in run_dirs(case):
        c = load(os.path.join(d, 'cfl.dat'))
        fs = glob.glob(os.path.join(d, '*.std'))
        if c is None or len(fs) != 1:
            continue
        txt = open(fs[0], errors='replace').read()
        m = re.search(r'^\s+step \(solve\)\s+(\S+)\s+(\S+)\s+\S+ %', txt, re.M)
        w = re.search(r'^wall (\S+) s', txt, re.M)
        sh = {}
        for k in COST_TIMERS:
            s = re.search(r'^\s+%s\s+\S+\s+\S+\s+(\S+) %%' % re.escape(k), txt, re.M)
            sh[k] = float(s.group(1)) if s else None
        it = c[c[:, 15] > 1, 15]
        out.append(dict(run=os.path.relpath(d, case), t0=float(c[0, 0]), t1=float(c[-1, 0]), steps=len(c),
                        iter_check=float(it.mean()) if len(it) else None, share=sh,
                        solve=float(m.group(2)) if m else None, wall=float(w.group(1)) if w else None))
    return out


def per_step(costs, key):
    c = [x for x in costs if x[key] is not None]
    return sum(x[key] for x in c)/sum(x['steps'] for x in c) if c else None


def upto(case, tc, w):
    """E_K, E_M, rho_min, rho_max at tc (conservation.dat, interpolated), the peak of eps_w on [0, tc] (window w,
    TURB_EPS_WIN) and eps_w of the last window, and the shares of the dissipation on [0, tc] (budget.dat, the
    steps that end by tc) as analyze_stage4.turb_run."""
    cons = np.concatenate(load_runs(case, 'conservation.dat'))
    t, k = np.unique(cons[:, 0], return_index=True)
    tg = np.arange(0.0, tc - w + 1e-12*tc, 0.25*w)
    ew = (np.interp(tg, t, cons[k, 10]) - np.interp(tg + w, t, cons[k, 10]))/w
    i = int(np.argmax(ew))
    bud = np.concatenate(load_runs(case, 'budget.dat'))
    b = bud[bud[:, 0] + bud[:, 1] <= tc*(1 + 1e-12)]
    dt = b[:, 1]
    tot = float(np.sum(dt*b[:, 7]))
    share = dict(ideal=float(np.sum(dt*b[:, 9]))/tot, visc=float(np.sum(dt*b[:, 10]))/tot,
                 svv=float(np.sum(dt*(b[:, 20] - b[:, 21])))/tot, svv_b=float(np.sum(dt*b[:, 21]))/tot,
                 eta=float(np.sum(dt*b[:, 26]))/tot)
    return dict(state={q: float(np.interp(tc, t, cons[k, j])) for q, j in TRACK_COLS.items()},
                peak=dict(t=float(tg[i] + 0.5*w), eps=float(ew[i])), eps_last=float(ew[-1]), share=share)


def turb_convergence_upto(runs, paths, tc):
    """S4.3 of analyze_stage4.turb_convergence for turb1_N128_a against N256 on [0, tc] (the chain of N256 cut
    before t_end): the state at tc in place of t_end, the shares on [0, tc], the peak of eps only when both runs
    have passed it on [0, tc] (eps_w, TURB_PEAK_DROP; otherwise not judged), and k_eff and a90 at the records t <= tc
    (turb_run keeps no later record of N256). Items that cannot be judged are listed in 'unjudged'."""
    c = a4.turb_convergence({k: dict(r, complete=True) for k, r in runs.items()})['turb1']  # pairs, ratios
    tau0 = runs[N256]['info']['tau0']
    a, b = (upto(paths[k], tc, TURB_EPS_WIN*tau0) for k in ('turb1_N128_a', N256))
    c['t_c'] = tc
    c['state_t'] = {k: v for k, v in c['state_t'].items() if float(k)*tau0 <= tc}
    c['state'] = {q: a4.rel(a['state'][q], b['state'][q]) for q in ('rho_min', 'rho_max', 'E_K', 'E_M')}
    c['share_visc'] = a4.rel(a['share']['visc'], b['share']['visc'])
    c['upto'] = {'turb1_N128_a': a, N256: b}
    ok = max(c['state'].values()) <= a4.OT_CONV_STATE and c['share_visc'] <= a4.OT_CONV_SHARE
    unjudged = []
    c['peak_passed'] = {k: x['eps_last'] <= (1 - TURB_PEAK_DROP)*x['peak']['eps']
                        for k, x in (('turb1_N128_a', a), (N256, b))}
    if all(c['peak_passed'].values()):
        c['eps_peak'] = a4.rel(a['peak']['eps'], b['peak']['eps'])
        c['t_peak'] = a4.rel(a['peak']['t'], b['peak']['t'])
        ok &= c['eps_peak'] <= a4.OT_CONV_EPS and c['t_peak'] <= a4.OT_CONV_TPEAK
    else:
        c.pop('eps_peak', None)
        c.pop('t_peak', None)
        unjudged += ['eps_peak', 't_peak']
    if c.get('keff_ratio'):
        ok &= all(a4.OT_KEFF_RATIO[0] <= v <= a4.OT_KEFF_RATIO[1] for v in c['keff_ratio'].values())
    else:
        unjudged.append('keff_ratio')
    if c['a90_last']:
        ok &= all(a4.OT_ZFRAC_RATIO[0] <= v <= a4.OT_ZFRAC_RATIO[1] for v in c['a90_last'].values())
    else:
        unjudged.append('a90')
    c['unjudged'] = unjudged
    c['ok'] = bool(ok)
    return c


def turb(root):
    base4 = os.path.join(base4_of(root), 'turb')
    paths = {'turb1_N64': os.path.join(base4, 'turb1_N64'), 'turb1_N128_a': os.path.join(base4, 'turb1_N128_a'),
             N256: os.path.join(root, N256)}
    runs, checks, dts, costs = {}, {}, {}, {}
    for name, case in paths.items():
        m = a4.case_info(case)
        r = runs[name] = a4.turb_run(case, m)
        r['R'] = a4.turb_residual(r)
        r['status'] = a4.status(case)
        dts[name] = dt_stats(dt_history(case), m['tau0'])[0]
        c = chunk_costs(case)
        costs[name] = dict(runs=c, solve_per_step=per_step(c, 'solve'), wall_per_step=per_step(c, 'wall'),
                           nproc=m['nproc'])
        s = dts[name]
        cut = name == N256 and r['status'] == N256_CUT and not r['complete']
        checks['run ' + name] = dict(ok=bool((r['status'] == 'exit 0' and r['complete'] or cut) and s['finite']
                                             and r['R'] <= a4.TOL_RES), complete=r['complete'], cut=cut, R=r['R'])
        a4.turb_print(name, r['status'], r)
        print('   dt: %s steps %d changes %d (raises %d) max dt/dt_cfl %.4f max dt*lam (check) %.3f mean dt/dt_cfl %s'
              % (s['reset'], s['steps'], s['changes'], s['raises'], s['max_dt_over_cfl'], s['max_dt_lam_check'],
                 {k: (round(v, 3) if v is not None else None) for k, v in s['mean_ratio'].items()}))
        cs = costs[name]
        print('   cost: %d runs, step (solve) %s s/step, wall %s s/step on %d ranks'
              % (len(c), '%.4f' % cs['solve_per_step'] if cs['solve_per_step'] else '-',
                 '%.4f' % cs['wall_per_step'] if cs['wall_per_step'] else '-', m['nproc']))
        for x in c if len(c) > 1 else ():
            f = lambda v, p: p % v if v is not None else '-'
            print('     %-9s t %.4f-%.4f steps %5d wall %s solve %s s/step %s mean iter (check) %s'
                  % (x['run'], x['t0'], x['t1'], x['steps'], f(x['wall'] and x['wall']/x['steps'], '%.3f'),
                     f(x['solve'] and x['solve']/x['steps'], '%.3f'),
                     ' '.join('%s %s%%' % (k, f(v, '%.1f')) for k, v in x['share'].items()),
                     f(x['iter_check'], '%.1f')))
    if runs[N256]['complete']:
        conv = a4.turb_convergence(runs)
    else:
        conv = {'turb1': turb_convergence_upto(runs, paths, runs[N256]['t_last'])}
    a4.turb_print_conv(conv)
    for c in conv.values():
        for k in ('t_c', 'peak_passed', 'unjudged', 'upto'):
            if k in c:
                print('   %s %s' % (k, c[k]))
    for name, c in conv.items():
        if 'ok' in c:
            checks['conv ' + name] = dict(ok=c['ok'])
    for k, c in checks.items():
        if not c['ok']:
            print('%-24s FAIL %s' % (k, c))
    a4.turb_figures(root, runs)
    ok_all = bool(checks) and all(c['ok'] for c in checks.values())
    print('turb:', 'PASS' if ok_all else 'FAIL')
    return dict(runs=runs, dt=dts, cost=costs, conv=conv, checks=checks, ok=ok_all)


BATCHES = {'regress': regress, 'track': track, 'turb': turb}


if __name__ == '__main__':
    batch = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.join('runs', 'stage6')
    root = os.path.join(base, batch)
    res = BATCHES[batch](root)
    with open(os.path.join(root, 'summary.json'), 'w') as f:
        json.dump(res, f, indent=1, default=str)
