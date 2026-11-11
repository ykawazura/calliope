#!/usr/bin/env python3
"""Stage 13: ES-SVV without the mass fluxes (svv_mass = F) on turb1 (M_s = 1, 128^3) against the control of Stage 8
(stage8/coef/e0.5_C2), from the runs of gen_stage13.py.

    python3 analyze_stage13.py selftest      synthetic checks of the tools below
    python3 analyze_stage13.py verify ROOT   gates V1-V5 of the verification job -> verify.json
    python3 analyze_stage13.py rhomin ROOT   rho_min(t), L, m and the checks of N and C -> rhomin.json
    python3 analyze_stage13.py fields ROOT   H, V and the shell spectrum of s at T_OUT -> fields_<N|C>.json
    python3 analyze_stage13.py judge ROOT    R, the label, the reported quantities -> summary.json, fig/*.png
    python3 analyze_stage13.py instab ROOT   the stop of N (dt, lambda_svv, rho_min, err_t; s at the out3d records of N
                                             and C) -> instab.json, fig/stage13_instab_*.png (ps1)

(the files in ROOT/ST/; ROOT/stage13/ and ROOT/stage13b/ hold the failed verifications (jobs 486708 and 487523)
and are not written; stage13c failed before the job and has no runs)

ROOT is the runs directory: ST/regress and ST/track (verify), ST/nomass/e0.5_C2 (N) and stage8/coef/e0.5_C2 (C; only
read, nothing is written into it). verify, rhomin and fields run on ps1 (fields reads the out3d records on the login
node); judge reads the JSON files only.

Gates (tasks/todo.md, Stage 13, revised with the approval of the user before the second verification job; if one
fails, the main run is not submitted)
    V1  analyze_stage6.regress on ST/regress: the 8 pairs equal (ISOTH bitwise, incomp_ot16 byte for byte)
    V2  analyze_stage6.track_decomp on ST/track (svv_mass = F): P = 1, 8 against P = 64
    V3  the 3 runs of track pass judge_stage1.checks_common with the ES-SVV checks (dMdt_raw <= 1e-12, dP <= 1e-11,
        R_zeta <= 1e-10, R_svv <= 1e-10, min D_svv >= 0), with completed and R_zeta of v3_case; completed is False
        for the control C (ended by t_end), which shows that this check can fail
    V4  calliope.in of the 3 runs has svv_mass = F; on ic/, one power iteration from the initial vector of the code
        (mass = False) gives lambda within V4_SEED for h = 1e-7 and 1e-5 (J v0 is not round-off, which would scale as
        1/h); lambda_svv of the first row of cfl.dat differs by at most V4_TOL from svv3d (mass = False, the initial
        vector of the code, svv_pi_init iterations, no tol) on ic/, and by at least V4_SEP from the same with
        mass = True
    V5  .std of the 3 runs has the line of svv_mass = F and not "Reading svv_parameters failed"; the 3 runs used
        gen_stage13.EXE_NEW (file exe), whose md5 is recorded
Judgment (the same section)
    R = L_N - L_C (analyze_stage8.windows on N and C only) over WIN if both reached T_STOP tau0, else over
    WIN_SHORT; no label if N did not reach 1.0 tau0. Label (unrounded R): R >= +DELTA 'shallower', R <= -DELTA
    'deeper', else 'unchanged within DELTA'.
Reported, not used: ln(m_N/m_C) (direction by DELTA), H and V (by 10%) as in Stage 8; R over the other window;
steps, the medians of dt and lambda_svv over each run (cfl.dat) and their ratios N/C; max R_svv (rounding scale),
min and max D_svv (budget.dat); the shell spectrum E_s(n) of s at T_OUT and the fractions
    f_hi = sum_{n >= k_c} E_s/sum_{n >= 1} E_s,  f_mid = sum_{k_c/2 <= n < k_c} E_s/sum_{n >= 1} E_s
(n the integer shell of |k| (shells), k_c = 2/3 of the Nyquist number, n = 0 (the mean of s) left out),
their ratios N/C at each record and the flag 'more high wavenumbers of s' if a ratio exceeds HI_RATIO at any record.
"""
import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_ic  # noqa: E402
import svv3d  # noqa: E402
import analyze_stage4 as a4  # noqa: E402
import analyze_stage6 as a6  # noqa: E402
import analyze_stage7 as a7  # noqa: E402
import analyze_stage8 as a8  # noqa: E402
import judge_stage1 as j1  # noqa: E402
from analyze_stage1 import analyze, load_runs, read_nml, run_dirs  # noqa: E402
from gen_stage13 import EXE_NEW  # noqa: E402

TAU0 = a8.TAU0
ST = 'stage13d'   # the initial vector for svv_mass = F set by svv_seed_nomass (stage13, 13b, 13c: failed)
RUNS = dict(N=os.path.join(ST, 'nomass', 'e0.5_C2'), C=os.path.join('stage8', 'coef', 'e0.5_C2'))
DELTA, WIN, WIN_SHORT, T_OUT, T_WIN = a8.DELTA, a8.WIN, a8.WIN_SHORT, a8.T_OUT, a8.T_WIN
WL, WS = '%g-%g' % WIN, '%g-%g' % WIN_SHORT
H_THR, V_THR = a8.H_THR, a8.V_THR
DIR_TOL = 0.1      # H, V: a direction if |ratio - 1| >= DIR_TOL (Stage 8)
HI_RATIO = 10.0    # f_hi, f_mid: the flag if N/C > HI_RATIO at a record
V4_TOL, V4_SEP, V4_SEED = 1e-6, 1e-4, 1e-6
H_SEED = (1e-7, 1e-5)  # V4: the finite-difference steps of the check of the initial vector
MASS_LINE = 'svv_mass = F: no mass fluxes'  # shock_viscosity.F90, init_svv
READ_FAIL = 'Reading svv_parameters failed'
LABELS = ('shallower', 'deeper', 'unchanged within DELTA')
COLOR = dict(N='#2a78d6', C='#0b0b0b')  # categorical slot 1 and primary ink (Stage 8)
NAME = dict(N='svv_mass = F (Stage 13)', C='control (Stage 8)')


def write_json(root, name, obj):
    d = os.path.join(root, ST)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), 'w') as f:
        json.dump(obj, f, indent=1, default=lambda x: x.tolist() if hasattr(x, 'tolist') else str(x))


def read_json(root, name):
    p = os.path.join(root, ST, name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


# ------------------------------------------------------------------------------------------------ verify

def steps_of(case):
    return sum(len(c) for c in load_runs(case, 'cfl.dat'))


def v3_case(case):
    """judge_stage1.common of a run of 'track' with 2 items replaced, then checks_common; ok needs the ES-SVV checks
    (es, es_identity) among the 6. completed: analyze_stage6.ran_nstep, as in track_decomp (the status of such a run
    is 'dt_changed', not 'exit 0'). R_zeta: the ratio of the RK-combined budget, max_abs_R_zeta_over_visc; the max
    over the stages (R_zeta_stage, max_R_zeta_ratio of judge_stage1) is O(0.3) at step 0, where for the solenoidal
    u of the IC numerator and denominator are both round-off."""
    o = analyze(case)
    r = j1.common(o)
    steps = steps_of(case)
    r.update(completed=a6.ran_nstep(case, steps), steps=steps, R_zeta_stage=r.pop('R_zeta', None))
    if o.get('budget'):
        r['R_zeta'] = o['budget']['max_abs_R_zeta_over_visc']
    c = j1.checks_common(r)
    return dict(ok=bool(r['es'] and r['es_identity'] and len(c) == 6 and all(c)), checks=[bool(x) for x in c],
                **{k: r.get(k) for k in ('exit', 'completed', 'steps', 'dMdt_raw', 'dP', 'R_zeta', 'R_zeta_stage',
                                         'R_svv', 'min_D_svv', 'es', 'es_identity')})


def v4_case(case):
    """On case/ic: lambda of one power iteration (mass = False) for the steps H_SEED, and lambda_svv of the first
    step (svv_pi_init iterations from the initial vector of the code, no tol) with (mass = True) or without the mass
    fluxes, against the first row of cfl.dat."""
    g, op = svv3d.from_case(case)
    mass_in = op.mass
    niter = int(read_nml(os.path.join(case, 'calliope.in')).get('svv_pi_init', '50'))
    lk, mk, bk, _ = gen_ic.read_restart(os.path.join(case, 'ic'), g)
    s, m, rho, u, b = svv3d.real_state(g, lk, mk, bk)
    eps = op.eps(rho, u, b)
    lam = {}
    for mass in (True, False):
        op.mass = mass
        lam[mass] = op.power(s, m, eps, 'rho', niter, 0.0)[0]
    l1 = [op.lam_first(s, m, eps, 'rho', h) for h in H_SEED]
    lam_code = float(load_runs(case, 'cfl.dat')[0][0, 12])
    r = dict(mass_in=mass_in, lam1=l1, seed_rel=abs(l1[0]/l1[1] - 1), lam_code=lam_code, lam_F=lam[False],
             lam_T=lam[True], rel_F=abs(lam_code/lam[False] - 1), rel_T=abs(lam_code/lam[True] - 1))
    r['ok'] = bool(mass_in is False and r['seed_rel'] <= V4_SEED and r['rel_F'] <= V4_TOL and r['rel_T'] >= V4_SEP)
    return r


def std_text(case):
    return ''.join(open(os.path.join(d, 'calliope.out.std'), errors='replace').read() for d in run_dirs(case)
                   if os.path.exists(os.path.join(d, 'calliope.out.std')))


def v5_case(case, repo):
    txt = std_text(case)
    exe = open(os.path.join(case, 'exe')).read().strip()
    with open(os.path.join(repo, exe), 'rb') as f:
        md5 = hashlib.md5(f.read()).hexdigest()
    r = dict(line=MASS_LINE in txt, read_fail=READ_FAIL in txt, exe=exe, md5=md5)
    r['ok'] = bool(r['line'] and not r['read_fail'] and exe == EXE_NEW)
    return r


def verify(root):
    repo = os.path.dirname(os.path.abspath(root))
    st13 = os.path.join(root, ST)
    out = dict(V1=a6.regress(os.path.join(st13, 'regress')))
    out['V1']['PASS'] = bool(out['V1']['ok'] and len(out['V1']['cases']) == 8)
    track = os.path.join(st13, 'track')
    dec = a6.track_decomp(track)
    out['V2'] = dict(cases=dec, PASS=bool(len(dec) == 3 and all(d['ok'] for d in dec.values())))
    cases = {'np%d' % P: os.path.join(track, 'dec_N64_np%d' % P) for P in a6.TDEC_P}
    for gate, fn in (('V3', v3_case), ('V4', v4_case), ('V5', lambda c: v5_case(c, repo))):
        res = {k: fn(c) for k, c in cases.items()}
        out[gate] = dict(cases=res, PASS=bool(all(r['ok'] for r in res.values())))
    out['V5']['md5'] = sorted(set(r['md5'] for r in out['V5']['cases'].values()))
    out['V5']['PASS'] &= len(out['V5']['md5']) == 1
    # the completion check of V3 can fail: the control ended by t_end
    ctrl = os.path.join(root, RUNS['C'])
    out['V3']['control_completed'] = bool(a6.ran_nstep(ctrl, steps_of(ctrl)))
    out['V3']['PASS'] &= not out['V3']['control_completed']
    for k, r in out['V3']['cases'].items():
        print('V3 %-5s ok=%s  %s, ran nstep %s (%d steps)  dMdt_raw %.1e  dP %.1e  R_zeta %.1e (max over the stages '
              '%.3f)  R_svv %.1e  min D_svv %.3e  es %s/%s'
              % (k, r['ok'], r['exit'], r['completed'], r['steps'], r['dMdt_raw'], r['dP'], np.nan
                 if r['R_zeta'] is None else r['R_zeta'], np.nan if r['R_zeta_stage'] is None else r['R_zeta_stage'],
                 r['R_svv'], r['min_D_svv'], r['es'], r['es_identity']))
    print('V3 control %s: ran nstep %s (must be False)' % (RUNS['C'], out['V3']['control_completed']))
    for k, r in out['V4']['cases'].items():
        print('V4 %-5s ok=%s  svv_mass in calliope.in %s  lambda_1(h=%g, %g) = %.10e, %.10e (rel %.1e)  lambda: code '
              '%.12e  mass=F %.12e (rel %.1e)  mass=T %.12e (rel %.2e)'
              % (k, r['ok'], r['mass_in'], H_SEED[0], H_SEED[1], r['lam1'][0], r['lam1'][1], r['seed_rel'],
                 r['lam_code'], r['lam_F'], r['rel_F'], r['lam_T'], r['rel_T']))
    for k, r in out['V5']['cases'].items():
        print('V5 %-5s ok=%s  line %s  read failed %s  exe %s  md5 %s' % (k, r['ok'], r['line'], r['read_fail'],
                                                                         r['exe'], r['md5']))
    out['PASS'] = all(out[g]['PASS'] for g in ('V1', 'V2', 'V3', 'V4', 'V5'))
    print('gates:', '  '.join('%s %s' % (g, 'PASS' if out[g]['PASS'] else 'FAIL')
                              for g in ('V1', 'V2', 'V3', 'V4', 'V5')), ' ->', 'PASS' if out['PASS'] else 'FAIL')
    write_json(root, 'verify.json', out)
    return out


# ------------------------------------------------------------------------------------------------ rhomin

def run_extras(case):
    """Steps, dt and lambda_svv of cfl.dat, R_svv and D_svv of budget.dat over the whole run; reads only."""
    cfl = np.concatenate(load_runs(case, 'cfl.dat'))
    bud = np.concatenate(load_runs(case, 'budget.dat'))
    D = -bud[:, 22]
    return dict(steps=len(cfl), dt_med=float(np.median(cfl[:, 1])), dt_min=float(cfl[:, 1].min()),
                dt_max=float(cfl[:, 1].max()), lam_med=float(np.median(cfl[:, 12])), lam_min=float(cfl[:, 12].min()),
                lam_max=float(cfl[:, 12].max()), max_dt_lam=float(cfl[:, 14].max()),
                R_svv_max=float(bud[:, 25].max()), R_svv_over_D_max=float(bud[:, 24].max()),
                D_min=float(D.min()), D_max=float(D.max()),
                series=dict(t_tau=(cfl[:, 0]/TAU0).tolist(), dt=cfl[:, 1].tolist(), lam=cfl[:, 12].tolist()))


def rhomin(root):
    res = {k: a8.run_stats(os.path.join(root, v)) for k, v in RUNS.items()}
    for wk, v in a8.windows(res).items():
        print('L over %s tau0: end %.6f tau0' % (wk, v['t1']))
        for k, r in res.items():
            r.setdefault('L', {})[wk] = v['L'][k]
            r.setdefault('L_end', {})[wk] = v['t1']
    for k, r in res.items():
        r['extras'] = run_extras(r['case'])
    # the common checks of N (analyze writes summary.json into N only)
    res['N']['common'] = j1.common(analyze(res['N']['case']))
    for k, r in res.items():
        e = r['extras']
        print('%s %-14s t_final/tau0 %.6f  L(%s) %s  L(%s) %s  m %.4e at %.3f  steps %d  dt med %.4e  lambda med %.4e  '
              'max dt lambda %.3f  R_svv max %.1e  D_svv %.3e .. %.3e  finite %s'
              % (k, r['exit_status'], r['t_final'], WL, a8.fmt(r['L'][WL]), WS, a8.fmt(r['L'][WS]), r['m'], r['t_m'],
                 e['steps'], e['dt_med'], e['lam_med'], e['max_dt_lam'], e['R_svv_max'], e['D_min'], e['D_max'],
                 r['finite']))
    c = res['N']['common']
    print('N common: exit %s  dMdt_raw %.1e  dP %.1e  R_zeta %.1e  R_svv %.1e  min D_svv %.3e  checks %s'
          % (c['exit'], c.get('dMdt_raw', np.nan), c.get('dP', np.nan), c.get('R_zeta', np.nan),
             c.get('R_svv', np.nan), c.get('min_D_svv', np.nan), j1.checks_common(c)))
    write_json(root, 'rhomin.json', res)
    return res


# ------------------------------------------------------------------------------------------------ fields

def band_fractions(E, kc):
    n = np.arange(len(E))
    tot = float(np.sum(E[n >= 1]))
    return dict(hi=float(np.sum(E[n >= kc]))/tot, mid=float(np.sum(E[(n >= 0.5*kc) & (n < kc)]))/tot)


def fields(root):
    for k, rel in RUNS.items():
        case = os.path.join(root, rel)
        g, op = svv3d.from_case(case)
        if not g.Lx == g.Ly == g.Lz:
            raise SystemExit('fields: the shells assume a cubic box')
        sh = shells(g)
        kc = (g.nlx//2)*2/3
        cons = a7.cons_table(case)
        res = dict(run=k, case=case, N=g.nlx, kc=kc, mass=op.mass, snaps=[])
        print('%s: %s, svv_mass %s' % (k, case, op.mass), flush=True)
        for tt in T_OUT:
            r = a7.pick(case, tt)
            if r is None:
                print('  t/tau0 %.4f: no record' % tt)
                continue
            sn = a7.Snap(g, op, r[1], r[2], full=False)
            st = a7.tail_stats(sn.rho, sn.s)
            E = spectrum(op, sn.s, sh)
            st.update(t=r[0], t_tau=r[0]/TAU0, t_req=tt, path=r[1], rec=r[2], spec_s=E, frac=band_fractions(E, kc),
                      cons=a7.cons_check(cons, r[0], st['rho_min']))
            res['snaps'].append(st)
            print('  t/tau0 %.4f  rho_min %.4e (cons rel %.1e)  minima < %s: %d  F_V(< %s) %.3e  f_hi %.3e  f_mid %.3e'
                  % (st['t_tau'], st['rho_min'], st['cons']['rel'], H_THR, st['nmin'][H_THR], V_THR,
                     st['F_V'][V_THR], st['frac']['hi'], st['frac']['mid']), flush=True)
        write_json(root, 'fields_%s.json' % k, res)


# ------------------------------------------------------------------------------------------------- judge

def label(R):
    return LABELS[0] if R >= DELTA else LABELS[1] if R <= -DELTA else LABELS[2]


def window(rm):
    if rm['N']['reached'] and rm['C']['reached']:
        return WL
    if rm['N']['L'][WS] is not None and rm['C']['L'][WS] is not None:
        return WS
    return None


def direction(x, tol):
    """-1 deeper holes, +1 shallower, 0 within tol (x: R or ln m; for H and V pass -(ratio - 1))."""
    return 0 if x is None or abs(x) < tol else int(np.sign(x))


def decide(rm, fl):
    """The summary from rhomin.json (rm) and fields_<N|C>.json (fl)."""
    w = window(rm)
    R = {wk: (rm['N']['L'][wk] - rm['C']['L'][wk] if rm['N']['L'][wk] is not None and rm['C']['L'][wk] is not None
              else None) for wk in (WL, WS)}
    out = dict(window=w, R=R[w] if w else None, label=label(R[w]) if w else None,
               R_other={wk: v for wk, v in R.items() if wk != w}, L={k: rm[k]['L'] for k in RUNS},
               L_end=rm['C']['L_end'])
    out['ln_m'] = float(np.log(rm['N']['m']/rm['C']['m']))
    common = [t for t in T_WIN if all(any(s['t_req'] == t for s in fl[k]['snaps']) for k in RUNS)]
    sn = {k: {s['t_req']: s for s in fl[k]['snaps']} for k in RUNS}
    H = {k: int(sum(sn[k][t]['nmin'][H_THR] for t in common)) for k in RUNS}
    V = {k: float(np.mean([sn[k][t]['F_V'][V_THR] for t in common])) if common else None for k in RUNS}
    out.update(common=common, H=H, V=V, H_ratio=H['N']/H['C'] if H['C'] else None,
               V_ratio=V['N']/V['C'] if V['C'] else None)
    out['dir'] = dict(R=direction(out['R'], DELTA), ln_m=direction(out['ln_m'], DELTA),
                      H=direction(None if out['H_ratio'] is None else -(out['H_ratio'] - 1), DIR_TOL),
                      V=direction(None if out['V_ratio'] is None else -(out['V_ratio'] - 1), DIR_TOL))
    e = {k: rm[k]['extras'] for k in RUNS}
    out['run'] = {k: dict(exit_status=rm[k]['exit_status'], t_final=rm[k]['t_final'], reached=rm[k]['reached'],
                          m=rm[k]['m'], t_m=rm[k]['t_m'], finite=rm[k]['finite'], mass=rm[k]['mass'],
                          **{q: v for q, v in e[k].items() if q != 'series'}) for k in RUNS}
    out['ratio'] = {q: e['N'][q]/e['C'][q] for q in ('steps', 'dt_med', 'lam_med')}
    out['common_checks_N'] = rm['N'].get('common')
    # the spectrum of s: fractions at the records of T_OUT present in both runs, their ratios and the flag
    ts = [t for t in T_OUT if t in sn['N'] and t in sn['C']]
    frac = {b: {k: [sn[k][t]['frac'][b] for t in ts] for k in RUNS} for b in ('hi', 'mid')}
    for b, v in frac.items():
        v['ratio'] = [x/y for x, y in zip(v['N'], v['C'])]
        v['max_ratio'] = max(v['ratio']) if ts else None
        v['flag'] = bool(ts) and v['max_ratio'] > HI_RATIO
    out['spec'] = dict(t=ts, frac=frac, flag=any(v['flag'] for v in frac.values()), kc=fl['C']['kc'],
                       mean={k: np.mean([sn[k][t]['spec_s'] for t in common], axis=0).tolist() if common else None
                             for k in RUNS})
    return out


def judge(root):
    rm = read_json(root, 'rhomin.json')
    fl = {k: read_json(root, 'fields_%s.json' % k) for k in RUNS}
    if rm is None or None in fl.values():
        raise SystemExit('judge: needs rhomin.json and fields_N.json, fields_C.json')
    s = decide(rm, fl)
    s['criteria'] = dict(WIN=WIN, WIN_SHORT=WIN_SHORT, DELTA=DELTA, DIR_TOL=DIR_TOL, HI_RATIO=HI_RATIO, H_THR=H_THR,
                         V_THR=V_THR, T_WIN=T_WIN, T_OUT=T_OUT)
    write_json(root, 'summary.json', s)
    print_judge(s)
    figures(root, rm, fl, s)
    return s


def print_judge(s):
    print('window %s (end %s tau0)  R = L_N - L_C = %s  label: %s' % (
        s['window'], a8.fmt(s['L_end'].get(s['window']) if s['window'] else None, '%.6f'), a8.fmt(s['R']), s['label']))
    print('L: N %s  C %s;  R over the other window %s' % (
        {wk: a8.fmt(v) for wk, v in s['L']['N'].items()}, {wk: a8.fmt(v) for wk, v in s['L']['C'].items()},
        {wk: a8.fmt(v) for wk, v in s['R_other'].items()}))
    print('ln(m_N/m_C) %+.4f  H %d/%d (%s)  V %.3e/%.3e (%s)  dir(R, ln m, H, V) %s  records %s' % (
        s['ln_m'], s['H']['N'], s['H']['C'], a8.fmt(s['H_ratio'], '%.3f'), s['V']['N'] or np.nan,
        s['V']['C'] or np.nan, a8.fmt(s['V_ratio'], '%.3f'),
        ','.join('%+d' % s['dir'][q] for q in ('R', 'ln_m', 'H', 'V')), s['common']))
    for k, r in s['run'].items():
        print('%s %-14s t_final/tau0 %.6f  m %.4e at %.3f  steps %d  dt med %.4e  lambda med %.4e  R_svv max %.1e  '
              'D_svv %.3e .. %.3e' % (k, r['exit_status'], r['t_final'], r['m'], r['t_m'], r['steps'], r['dt_med'],
                                     r['lam_med'], r['R_svv_max'], r['D_min'], r['D_max']))
    print('ratios N/C: %s' % {q: '%.4f' % v for q, v in s['ratio'].items()})
    for b, v in s['spec']['frac'].items():
        print('f_%-3s N %s  C %s  N/C %s  max %s  flag %s' % (
            b, ['%.2e' % x for x in v['N']], ['%.2e' % x for x in v['C']], ['%.2f' % x for x in v['ratio']],
            a8.fmt(v['max_ratio'], '%.3f'), v['flag']))
    print('more high wavenumbers of s:', s['spec']['flag'])


def figures(root, rm, fl, s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = os.path.join(root, ST, 'fig')
    os.makedirs(d, exist_ok=True)
    grid = dict(color='#e4e3df', lw=0.6)
    # rho_min, dt and lambda_svv against t
    fig, axs = plt.subplots(3, 1, figsize=(7.2, 7.6), sharex=True)
    for k in RUNS:
        ser, ex = rm[k]['series'], rm[k]['extras']['series']
        axs[0].semilogy(ser['t_tau'], ser['rho_min'], color=COLOR[k], lw=1.4, label=NAME[k])
        axs[1].plot(ex['t_tau'], ex['dt'], color=COLOR[k], lw=1.4, label=NAME[k])
        axs[2].plot(ex['t_tau'], ex['lam'], color=COLOR[k], lw=1.4, label=NAME[k])
    if s['window']:
        w = WIN if s['window'] == WL else WIN_SHORT
        axs[0].axvspan(w[0], s['L_end'][s['window']], color='#f0efec', zorder=0)
    for ax, yl in zip(axs, ('rho_min', 'dt', 'lambda_svv')):
        ax.set_ylabel(yl)
        ax.grid(**grid)
    axs[0].legend(frameon=False, fontsize=8)
    axs[-1].set_xlabel('t/tau0')
    fig.tight_layout()
    fig.savefig(os.path.join(d, 'stage13_rhomin.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    # the spectrum of s (mean over the common records of T_WIN), its ratio and the fractions against t
    kc = s['spec']['kc']
    fig, axs = plt.subplots(2, 2, figsize=(9.0, 6.6))
    if s['spec']['mean']['N'] is not None:
        E = {k: np.asarray(s['spec']['mean'][k]) for k in RUNS}
        n = np.arange(len(E['C']))
        sel = n >= 1
        for k in RUNS:
            axs[0, 0].loglog(n[sel], E[k][sel], color=COLOR[k], lw=1.4, label=NAME[k])
        axs[0, 1].loglog(n[sel], E['N'][sel]/E['C'][sel], color=COLOR['N'], lw=1.4)
        axs[0, 1].axhline(1, color='#898781', lw=0.8)
        for ax in axs[0]:
            for kk in (0.5*kc, kc):
                ax.axvline(kk, color='#898781', lw=0.8, ls=':')
            ax.set_xlabel('shell n = |k| L/2 pi')
        axs[0, 0].set_ylabel('E_s(n), mean over %s tau0' % ', '.join('%g' % t for t in s['common']))
        axs[0, 1].set_ylabel('E_s of N / E_s of C')
        axs[0, 0].legend(frameon=False, fontsize=8)
    for j, b in enumerate(('hi', 'mid')):
        v = s['spec']['frac'][b]
        for k in RUNS:
            axs[1, j].semilogy(s['spec']['t'], v[k], 'o-', color=COLOR[k], lw=1.4, ms=5, label=NAME[k])
        axs[1, j].set_xlabel('t/tau0')
        axs[1, j].set_ylabel('f_hi (n >= k_c)' if b == 'hi' else 'f_mid (k_c/2 <= n < k_c)')
        axs[1, j].legend(frameon=False, fontsize=8)
    for ax in axs.ravel():
        ax.grid(**grid)
    fig.tight_layout()
    fig.savefig(os.path.join(d, 'stage13_spec.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ instab

DT_FRAC = (0.5, 0.1, 1e-2, 1e-4, 1e-8)  # the first rows of N at t >= T_REF with dt below these fractions of DT_REF
T_REF = 0.3       # tau0; DT_REF = the median dt of N over t < T_REF
RK3_LIM = 2.51    # the parameter rk3_lim of advance.F90 (abort if dt lambda_svv > RK3_LIM at a power-method step)
T_TAB = (0.1, 0.2, 0.3, 0.4, 0.42, 0.44, 0.45, 0.46, 0.465)  # tau0; the rows of N and C at these times
DIP = (0.43, 0.45)  # tau0; the first drop of rho_min of N
ERR_X = 10.0      # the first row of N with |err_t| > ERR_X max|err_t| of C
L_CUT = 5         # s along x, y, z through the min of rho, L_CUT cells each side
L_RATE = 200      # the last L_RATE rows of conservation.dat: the slope of rho_min against t


def series(case):
    """The columns of cfl.dat (rows: steps), conservation.dat (unique t) and budget.dat used by instab."""
    cfl = np.concatenate(load_runs(case, 'cfl.dat'))
    cons = a7.cons_table(case)
    _, i = np.unique(cons[:, 0], return_index=True)
    cons = cons[i]
    bud = np.concatenate(load_runs(case, 'budget.dat'))
    return dict(t=cfl[:, 0]/TAU0, dt=cfl[:, 1], dt_adv=cfl[:, 2], dt_visc=cfl[:, 3], cmax=cfl[:, 4:7].max(1),
                eps=cfl[:, 9:12].max(1), lam=cfl[:, 12], dt_lam=cfl[:, 14], t_c=cons[:, 0]/TAU0, rho_min=cons[:, 11],
                E=cons[:, 10], t_b=bud[:, 0]/TAU0, dEdt=bud[:, 6], err_t=bud[:, 8], D=-bud[:, 22])


def near(t, x):
    return int(min(np.searchsorted(t, x), len(t) - 1))


def row(S, i):
    """Row i of cfl.dat and the row of conservation.dat at its t."""
    j = near(S['t_c'], S['t'][i])
    return dict(row=int(i), t_tau=float(S['t'][i]), rho_min=float(S['rho_min'][j]), E=float(S['E'][j]),
                **{q: float(S[q][i]) for q in ('dt', 'dt_adv', 'dt_visc', 'cmax', 'eps', 'lam', 'dt_lam')})


def shells(g):
    """Integer shells of |k| in units of 2 pi/L (a cubic box) and the weights of the r2c half spectrum."""
    n = np.rint(np.sqrt(g.ikx.reshape(-1, 1, 1)**2 + g.iky.reshape(1, 1, -1)**2 + g.ikz.reshape(1, -1, 1)**2))
    w = np.where((g.iky == 0) | (2*g.iky == g.nly), 1.0, 2.0).reshape(1, 1, -1)
    return n.astype(int), np.broadcast_to(w, n.shape)


def spectrum(op, f, sh):
    """sum over the shell of |f_k|^2 (the sum over the shells is <f^2>)."""
    n, w = sh
    return np.bincount(n.ravel(), weights=(w*np.abs(op.fft(f))**2).ravel(), minlength=int(n.max()) + 1)


def cube_fractions(g, op, f, sh):
    """The fractions of sum_{n >= 1} w|f_k|^2 (n, w of shells) outside the cube max_d |k_d| <=
    (nl_d - 1)//3 that the 2/3 filter keeps (H) and on the Nyquist planes |k_d| = nl_d/2 (nyq)."""
    n, w = sh
    P = w*np.abs(op.fft(f))**2
    H = np.zeros(n.shape, bool)
    nyq = np.zeros(n.shape, bool)
    for ik, nl in ((np.abs(g.ikx).reshape(-1, 1, 1), g.nlx), (np.abs(g.ikz).reshape(1, -1, 1), g.nlz),
                   (np.abs(g.iky).reshape(1, 1, -1), g.nly)):
        H |= ik > (nl - 1)//3
        nyq |= 2*ik == nl
    tot = P[n >= 1].sum()
    return dict(H=float(P[H].sum()/tot), nyq=float(P[nyq].sum()/tot))


def div_n(g, op, u):
    """Spectral div of u = (u_x, u_y, u_z) with k_d = 0 at Nyquist (div_N of advance.F90 without the filter)."""
    return sum(op.ifft(1j*k*op.fft(x)) for k, x in zip((g.kxd, g.kyd, g.kzd), u))


def spike(s, p):
    """s along x, y, z through p (L_CUT cells each side) and the depth mean(6 neighbours) - s(p)."""
    cuts = {lab: a7.cut(s, p, d, L_CUT) for lab, d in (('x', (1, 0, 0)), ('y', (0, 1, 0)), ('z', (0, 0, 1)))}
    nb = [c[L_CUT + o] for c in cuts.values() for o in (-1, 1)]
    return cuts, float(np.mean(nb) - s[p])


def snap_instab(g, op, sh, kc, path, rec):
    """At the min of rho of one record: s, the depth of spike, |u|, |B|, v_A = |B|/sqrt(rho), div u (spectral,
    k_d = 0 at Nyquist as div_N, no filter); over the box: max|u|, min div u, the cells with rho < 0.1 and 0.01, the
    shell spectrum of s and its fractions hi (band_fractions), H and nyq (cube_fractions). div_N m at the min of rho
    is -d rho/dt there (ds/dt = -(div_N m)/rho, m is 2/3-filtered after each stage: m_out = max|m_k| outside the
    mask / max|m_k| checks it)."""
    sn = a7.Snap(g, op, path, rec, full=True)
    p = tuple(int(a) for a in np.unravel_index(int(np.argmin(sn.rho)), sn.rho.shape))
    u2 = sum(x**2 for x in sn.u)
    b2 = sum(x**2 for x in sn.b)
    div = div_n(g, op, sn.u)
    m = [op.ifft(x) for x in sn.mk]
    divm = div_n(g, op, m)
    m_out = max(float(np.abs(x*(1 - g.mask)).max()/np.abs(x).max()) for x in sn.mk)
    E = spectrum(op, sn.s, sh)
    cuts, depth = spike(sn.s, p)
    return dict(rho_min=float(sn.rho[p]), s_min=float(sn.s[p]), cell_ikj=p, cell_xyz=(p[0], p[2], p[1]), depth=depth,
                cuts=cuts, u_max=float(np.sqrt(u2.max())), u_p=float(np.sqrt(u2[p])), B_p=float(np.sqrt(b2[p])),
                vA_p=float(np.sqrt(b2[p]/sn.rho[p])), div_p=float(div[p]), div_min=float(div.min()),
                divm_p=float(divm[p]), m_p=float(np.sqrt(sum(x[p]**2 for x in m))), m_out=m_out,
                n_01=int(np.sum(sn.rho < 0.1)), n_001=int(np.sum(sn.rho < 0.01)), sigma_s=float(sn.s.std()),
                frac=dict(hi=band_fractions(E, kc)['hi'], **cube_fractions(g, op, sn.s, sh)), spec_s=E)


def instab(root):
    """The stop of N: the reason (exit_status, ABORT lines of .std), the first rows with dt below DT_FRAC DT_REF, the
    rows at T_TAB, the drop of rho_min in DIP, E, err_t; at the out3d records of N (and of C at the same times):
    snap_instab. -> instab.json, fig/stage13_instab_t.png, fig/stage13_instab_s.png."""
    cases = {k: os.path.join(root, v) for k, v in RUNS.items()}
    S = {k: series(c) for k, c in cases.items()}
    N, C = S['N'], S['C']
    nml = read_nml(os.path.join(cases['N'], 'calliope.in'))
    st = os.path.join(cases['N'], 'exit_status')
    out = dict(criteria=dict(DT_FRAC=DT_FRAC, T_REF=T_REF, RK3_LIM=RK3_LIM, T_TAB=T_TAB, DIP=DIP, ERR_X=ERR_X,
                             L_CUT=L_CUT, L_RATE=L_RATE))
    out['stop'] = dict(exit_status=open(st).read().strip() if os.path.exists(st) else None,
                       abort=[x for x in std_text(cases['N']).splitlines() if 'ABORT' in x], steps=len(N['t']),
                       t_final=float(N['t'][-1] + N['dt'][-1]/TAU0), dt_last=float(N['dt'][-1]),
                       finite=bool(all(np.all(np.isfinite(v)) for v in N.values())),
                       dt_lam_max=float(N['dt_lam'].max()), dt_lam_max_t=float(N['t'][np.argmax(N['dt_lam'])]),
                       s_abort=nml.get('s_abort'), rk3_lim=RK3_LIM)
    ref = float(np.median(N['dt'][N['t'] < T_REF]))
    out['DT_REF'] = ref
    out['below'] = {}
    for f in DT_FRAC:
        i = np.nonzero((N['dt'] < f*ref) & (N['t'] >= T_REF))[0]
        out['below']['%g' % f] = row(N, int(i[0])) if len(i) else None
    out['table'] = {k: [row(S[k], near(S[k]['t'], x)) for x in T_TAB if x <= S[k]['t'][-1]] for k in RUNS}
    out['table']['N'].append(row(N, len(N['t']) - 1))
    # the drop of rho_min in DIP and the largest rho_min after it until dt < 0.1 DT_REF
    t01 = out['below']['0.1']['t_tau']
    sel = np.nonzero((N['t_c'] >= DIP[0]) & (N['t_c'] <= DIP[1]))[0]
    j = int(sel[np.argmin(N['rho_min'][sel])])
    aft = np.nonzero((N['t_c'] > N['t_c'][j]) & (N['t_c'] < t01))[0]
    k = int(aft[np.argmax(N['rho_min'][aft])])
    selC = np.nonzero((C['t_c'] >= DIP[0]) & (C['t_c'] <= DIP[1]))[0]
    out['dip'] = dict(t_tau=float(N['t_c'][j]), rho_min=float(N['rho_min'][j]), row=row(N, near(N['t'], N['t_c'][j])),
                      rho_min_C=float(C['rho_min'][near(C['t_c'], N['t_c'][j])]),
                      rho_min_C_min=float(C['rho_min'][selC].min()), t_back=float(N['t_c'][k]),
                      rho_back=float(N['rho_min'][k]))
    out['E'] = {q: dict(first=float(S[q]['E'][0]), last=float(S[q]['E'][-1]), n_inc=int(np.sum(np.diff(S[q]['E']) > 0)),
                        max_inc=float(np.diff(S[q]['E']).max())) for q in RUNS}
    errC = float(np.abs(C['err_t']).max())
    big = np.nonzero(np.abs(N['err_t']) > ERR_X*errC)[0]
    out['err_t'] = dict(C_max=errC, N_max=float(np.abs(N['err_t']).max()),
                        N_max_t=float(N['t_b'][np.argmax(np.abs(N['err_t']))]),
                        before={q: float(np.abs(S[q]['err_t'][S[q]['t_b'] < t01]).max()) for q in RUNS},
                        first=None if not len(big) else dict(t_tau=float(N['t_b'][big[0]]),
                                                              err_t=float(N['err_t'][big[0]]),
                                                              dEdt=float(N['dEdt'][big[0]]), D=float(N['D'][big[0]])))
    # the out3d records of N after t = 0 and the records of C at the same times
    op_g = {q: svv3d.from_case(c) for q, c in cases.items()}
    g = op_g['N'][0]
    sh = shells(g)
    kc = (g.nlx//2)*2/3
    out['records'] = []
    for t, path, rec in a4.out3d_records(cases['N']):
        if t/TAU0 < a7.T_TOL:
            continue
        r = dict(t_tau=t/TAU0, N=snap_instab(*op_g['N'], sh, kc, path, rec), cfl_N=row(N, near(N['t'], t/TAU0)))
        rc = a7.pick(cases['C'], t/TAU0)
        r['C'] = None if rc is None else dict(t_tau=rc[0]/TAU0, **snap_instab(*op_g['C'], sh, kc, rc[1], rc[2]))
        out['records'].append(r)
    out['kc'] = kc
    # d rho_min/dt over the last L_RATE rows against -div_N m at the min of rho of the last record; ln dt and
    # ln lambda_svv against ln rho_min from the first row with dt < 1e-4 DT_REF to the end
    tc = N['t_c'][-L_RATE:]*TAU0
    i0 = out['below']['0.0001']['row']
    rho = N['rho_min'][np.minimum(np.searchsorted(N['t_c'], N['t'][i0:]), len(N['t_c']) - 1)]
    out['rate'] = dict(L_RATE=L_RATE, t_tau=(float(tc[0]/TAU0), float(tc[-1]/TAU0)),
                       drho_dt=float(np.polyfit(tc - tc[0], N['rho_min'][-L_RATE:], 1)[0]),
                       minus_divm=-out['records'][-1]['N']['divm_p'], t_rec=out['records'][-1]['t_tau'],
                       row0=i0, p_dt=float(np.polyfit(np.log(rho), np.log(N['dt'][i0:]), 1)[0]),
                       p_lam=float(np.polyfit(np.log(rho), np.log(N['lam'][i0:]), 1)[0]),
                       p_cmax=float(np.polyfit(np.log(rho), np.log(N['cmax'][i0:]), 1)[0]))
    print_instab(out)
    write_json(root, 'instab.json', out)
    instab_figures(root, S, out)
    return out


def print_instab(o):
    s = o['stop']
    print('stop: exit_status %s  ABORT lines %d  steps %d  t_final/tau0 %.6f  dt last %.3e  finite %s  max dt lambda '
          '%.3f at %.6f tau0  (s_abort %s, rk3_lim %s)'
          % (s['exit_status'], len(s['abort']), s['steps'], s['t_final'], s['dt_last'], s['finite'], s['dt_lam_max'],
             s['dt_lam_max_t'], s['s_abort'], s['rk3_lim']))
    fmt = ('  %-6s row %5d  t/tau0 %.6f  dt %.3e  dt_adv %.3e  dt_visc %.3e  cmax %.4e  eps %.3e  lambda %.4e  '
           'dt lambda %.3f  rho_min %.4e')
    print('DT_REF (median dt of N over t < %g tau0) %.4e; first rows at t >= %g tau0 below f DT_REF:'
          % (T_REF, o['DT_REF'], T_REF))
    for f, r in o['below'].items():
        print(fmt % ((f, r['row'], r['t_tau'], r['dt'], r['dt_adv'], r['dt_visc'], r['cmax'], r['eps'], r['lam'],
                      r['dt_lam'], r['rho_min'])) if r else '  %s: none' % f)
    for k, rows in o['table'].items():
        print('rows of %s:' % k)
        for r in rows:
            print(fmt % (('', r['row'], r['t_tau'], r['dt'], r['dt_adv'], r['dt_visc'], r['cmax'], r['eps'], r['lam'],
                          r['dt_lam'], r['rho_min'])))
    d = o['dip']
    print('dip: N rho_min %.4e at %.5f tau0 (C %.4e there, min over DIP %.4e); back to %.4e at %.5f tau0; dt/DT_REF '
          '%.3f, cmax %.3e there' % (d['rho_min'], d['t_tau'], d['rho_min_C'], d['rho_min_C_min'], d['rho_back'],
                                     d['t_back'], d['row']['dt']/o['DT_REF'], d['row']['cmax']))
    for k, e in o['E'].items():
        print('E %s: %.6f -> %.6f  increases %d (max %.2e)' % (k, e['first'], e['last'], e['n_inc'], e['max_inc']))
    e = o['err_t']
    print('err_t: max|.| C %.2e, N %.2e at %.6f tau0; before dt < 0.1 DT_REF: N %.2e, C %.2e; first N > %g C max: %s'
          % (e['C_max'], e['N_max'], e['N_max_t'], e['before']['N'], e['before']['C'], ERR_X,
             None if e['first'] is None else '%.6f tau0 (err_t %.3e, dE/dt %.3e, D_svv %.3e)'
             % (e['first']['t_tau'], e['first']['err_t'], e['first']['dEdt'], e['first']['D'])))
    for r in o['records']:
        for k in RUNS:
            x = r[k]
            if x is None:
                continue
            print('%.5f %s  rho_min %.4e at (x,y,z) %s  s %.3f  depth %.3f  |u| %.3e (max %.3e)  |B| %.3f  v_A %.4e  '
                  'div u %.3e (min %.3e)  n(<0.1) %d n(<0.01) %d  f_hi %.3e f_H %.3e f_nyq %.3e'
                  % (r['t_tau'], k, x['rho_min'], x['cell_xyz'], x['s_min'], x['depth'], x['u_p'], x['u_max'],
                     x['B_p'], x['vA_p'], x['div_p'], x['div_min'], x['n_01'], x['n_001'], x['frac']['hi'],
                     x['frac']['H'], x['frac']['nyq']))
        print('         cmax of cfl.dat at the record (N) %.4e;  N: div_N m %.4e, |m| %.3e at the min of rho, '
              'm_out %.1e' % (r['cfl_N']['cmax'], r['N']['divm_p'], r['N']['m_p'], r['N']['m_out']))
    r = o['rate']
    print('d rho_min/dt over the last %d rows (%.7f-%.7f tau0) %.4e;  -div_N m at the min of rho (%.5f tau0) %.4e'
          % (r['L_RATE'], r['t_tau'][0], r['t_tau'][1], r['drho_dt'], r['t_rec'], r['minus_divm']))
    print('from row %d: d ln dt / d ln rho_min %.3f, d ln lambda / d ln rho_min %.3f, d ln cmax / d ln rho_min %.3f'
          % (r['row0'], r['p_dt'], r['p_lam'], r['p_cmax']))


def instab_figures(root, S, o):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = os.path.join(root, ST, 'fig')
    os.makedirs(d, exist_ok=True)
    grid = dict(color='#e4e3df', lw=0.6)
    N = S['N']
    tN = o['stop']['t_final']
    r0 = max(0, o['below']['0.1']['row'] - 300)
    rows = np.arange(r0, len(N['t']))
    rho_rows = N['rho_min'][np.minimum(np.searchsorted(N['t_c'], N['t'][rows]), len(N['t_c']) - 1)]
    # left: rho_min, dt, cmax and lambda_svv against t (N and C up to the end of N); right: N against the step
    fig, axs = plt.subplots(4, 2, figsize=(10.0, 9.6), sharex='col')
    for k in RUNS:
        c = S[k]
        m = c['t'] <= tN
        mc = c['t_c'] <= tN
        axs[0, 0].semilogy(c['t_c'][mc], c['rho_min'][mc], color=COLOR[k], lw=1.4, label=NAME[k])
        axs[1, 0].semilogy(c['t'][m], c['dt'][m], color=COLOR[k], lw=1.4, label='dt, ' + NAME[k])
        axs[2, 0].semilogy(c['t'][m], c['cmax'][m], color=COLOR[k], lw=1.4, label=NAME[k])
        axs[3, 0].semilogy(c['t'][m], c['lam'][m], color=COLOR[k], lw=1.4, label=NAME[k])
    for ax, x, sl in ((axs[1, 0], N['t'], slice(None)), (axs[1, 1], rows, rows)):
        ax.semilogy(x, N['dt_adv'][sl], color=COLOR['N'], lw=1.0, ls='--', label='dt_adv, N')
        ax.semilogy(x, N['dt_visc'][sl], color=COLOR['N'], lw=1.0, ls=':', label='dt_visc, N')
    axs[0, 1].semilogy(rows, rho_rows, color=COLOR['N'], lw=1.4)
    axs[1, 1].semilogy(rows, N['dt'][rows], color=COLOR['N'], lw=1.4, label='dt, N')
    axs[2, 1].semilogy(rows, N['cmax'][rows], color=COLOR['N'], lw=1.4)
    axs[3, 1].semilogy(rows, N['lam'][rows], color=COLOR['N'], lw=1.4)
    for f in ('0.1', '1e-08'):
        if o['below'].get(f):
            for ax in axs[:, 1]:
                ax.axvline(o['below'][f]['row'], color='#898781', lw=0.8, ls=':')
    for i, yl in enumerate(('rho_min', 'dt', 'max(|u_d| + c_f,d)', 'lambda_svv')):
        for ax in axs[i]:
            ax.set_ylabel(yl)
            ax.grid(**grid)
    axs[0, 0].legend(frameon=False, fontsize=8)
    axs[1, 0].legend(frameon=False, fontsize=7, ncol=2)
    axs[1, 1].legend(frameon=False, fontsize=7)
    axs[-1, 0].set_xlabel('t/tau0')
    axs[-1, 1].set_xlabel('step of N (dotted: dt < 0.1 and 1e-8 DT_REF)')
    axs[0, 1].set_title('N from step %d to the end' % r0, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(d, 'stage13_instab_t.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    # s: the fractions against t, their ratios N/C, the spectra and the cuts through the min of rho at the end of N
    rec = o['records']
    both = [r for r in rec if r['C'] is not None]
    ls = dict(hi=('-', 'o'), H=('--', 's'), nyq=(':', '^'))
    lab = dict(hi='f_hi (n >= k_c)', H='f_H (outside the 2/3 cube)', nyq='f_nyq (Nyquist planes)')
    fig, axs = plt.subplots(2, 2, figsize=(10.0, 7.4))
    for b, (l, mk) in ls.items():
        axs[0, 0].semilogy([r['t_tau'] for r in rec], [r['N']['frac'][b] for r in rec], color=COLOR['N'], ls=l,
                           marker=mk, ms=5, lw=1.4, label='N, ' + lab[b])
        axs[0, 0].semilogy([r['t_tau'] for r in both], [r['C']['frac'][b] for r in both], color=COLOR['C'], ls=l,
                           marker=mk, ms=5, lw=1.4, label='C, ' + lab[b])
        axs[0, 1].plot([r['t_tau'] for r in both], [r['N']['frac'][b]/r['C']['frac'][b] for r in both],
                       color=COLOR['N'], ls=l, marker=mk, ms=5, lw=1.4, label=lab[b])
    axs[0, 1].axhline(HI_RATIO, color='#898781', lw=0.8, ls='--')
    axs[0, 1].axhline(1, color='#898781', lw=0.8)
    axs[0, 1].set_ylim(0, HI_RATIO*1.1)
    axs[0, 0].set_ylabel('fraction of the variance of s')
    axs[0, 1].set_ylabel('N / C (dashed: HI_RATIO)')
    for ax in axs[0]:
        ax.set_xlabel('t/tau0')
        ax.legend(frameon=False, fontsize=7)
    last = both[-1]
    n = np.arange(len(last['N']['spec_s']))
    for k, x, c, nm in (('N', last['N'], COLOR['N'], 'N %.4f' % last['t_tau']),
                        ('C', last['C'], COLOR['C'], 'C %.4f' % last['C']['t_tau']),
                        ('N', rec[-1]['N'], '#eb6834', 'N %.4f (last record)' % rec[-1]['t_tau'])):
        axs[1, 0].loglog(n[1:], np.asarray(x['spec_s'])[1:], color=c, lw=1.4, label=nm + ' tau0')
    axs[1, 0].axvline(o['kc'], color='#898781', lw=0.8, ls=':')
    axs[1, 0].set_xlabel('shell n = |k| L/2 pi (dotted: k_c)')
    axs[1, 0].set_ylabel('E_s(n)')
    axs[1, 0].legend(frameon=False, fontsize=7)
    off = np.arange(-L_CUT, L_CUT + 1)
    for lab_, c in (('x', COLOR['N']), ('y', '#eb6834'), ('z', COLOR['C'])):
        axs[1, 1].plot(off, rec[-1]['N']['cuts'][lab_], color=c, lw=1.4, marker='o', ms=4, label='along ' + lab_)
    axs[1, 1].set_xlabel('cells from the min of rho of N\nat %.4f tau0, (x, y, z) = %s'
                         % (rec[-1]['t_tau'], tuple(rec[-1]['N']['cell_xyz'])))
    axs[1, 1].set_ylabel('s = ln rho')
    axs[1, 1].legend(frameon=False, fontsize=7)
    for ax in axs.ravel():
        ax.grid(**grid)
    fig.tight_layout()
    fig.savefig(os.path.join(d, 'stage13_instab_s.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------- selftest

def synthetic(LN, LC, reachedN=(True, True), mN=0.05, mC=0.0566, hN=(5, 5), hC=(5, 5), vN=0.01, vC=0.01,
              fN=(1e-3, 1e-2), fC=(1e-3, 1e-2)):
    """rm and fl for decide(): L of N and C over (WL, WS), N reached (1.2, 1.0) tau0, m, the nmin of H_THR at the
    first and the other records of T_WIN, F_V of V_THR, (f_hi, f_mid) at the records of T_OUT."""
    def run(L, reached, m):
        ex = dict(steps=100, dt_med=1e-3, dt_min=1e-3, dt_max=1e-3, lam_med=1.0, lam_min=1.0, lam_max=1.0,
                  max_dt_lam=1e-3, R_svv_max=1e-14, R_svv_over_D_max=1e-12, D_min=1.0, D_max=2.0)
        return dict(exit_status='exit 0', t_final=1.2, reached=reached[0], m=m, t_m=0.9, finite=True, mass={},
                    L={WL: L[0] if reached[0] else None, WS: L[1] if reached[1] else None}, L_end={WL: 1.2, WS: 1.0},
                    extras=ex)
    rm = dict(N=run(LN, reachedN, mN), C=run(LC, (True, True), mC))

    def snaps(h, v, f):
        return [dict(t_req=t, nmin={H_THR: h[0] if t == T_WIN[0] else h[1]}, F_V={V_THR: v},
                     frac=dict(hi=f[0], mid=f[1]), spec_s=[0.0, 1.0, 0.5]) for t in T_OUT]
    fl = dict(N=dict(kc=1.5, snaps=snaps(hN, vN, fN)), C=dict(kc=1.5, snaps=snaps(hC, vC, fC)))
    return rm, fl


def selftest():
    ok = True

    def check(name, good):
        nonlocal ok
        ok &= bool(good)
        print('%-60s %s' % (name, 'PASS' if good else 'FAIL'))

    # 1. the label at the bounds, unrounded
    for R, want in ((DELTA, LABELS[0]), (np.nextafter(DELTA, 0), LABELS[2]), (-DELTA, LABELS[1]),
                    (np.nextafter(-DELTA, 0), LABELS[2]), (0.0, LABELS[2])):
        check('label(%r) = %s' % (R, want), label(R) == want)
    # 2. windows, R, the other window, ln m, H, V and their directions
    s = decide(*synthetic((-2.0, -2.1), (-2.3747, -2.4)))
    check('both reached 1.2: WIN, R = +0.3747, shallower',
          s['window'] == WL and abs(s['R'] - 0.3747) < 1e-12 and s['label'] == LABELS[0]
          and abs(s['R_other'][WS] - 0.3) < 1e-12)
    s = decide(*synthetic((-2.0, -2.45), (-2.3747, -2.4), reachedN=(False, True)))
    check('N reached 1.0 only: WIN_SHORT, R = -0.05, unchanged', s['window'] == WS and abs(s['R'] + 0.05) < 1e-12
          and s['label'] == LABELS[2] and s['R_other'][WL] is None)
    s = decide(*synthetic((-2.0, -2.45), (-2.3747, -2.4), reachedN=(False, False)))
    check('N short of 1.0: no label', s['window'] is None and s['label'] is None and s['R'] is None)
    s = decide(*synthetic((-2.6, -2.6), (-2.3747, -2.4), mN=0.0566*np.exp(-0.2), hN=(6, 6), vN=0.0089))
    H_want = 5 + 5*(len(T_WIN) - 1)
    check('deeper; ln m = -0.2; H, V directions',
          s['label'] == LABELS[1] and abs(s['ln_m'] + 0.2) < 1e-12 and s['H']['C'] == H_want
          and s['H']['N'] == 6*len(T_WIN) and s['dir'] == dict(R=-1, ln_m=-1, H=-1, V=+1))
    # 3. the flag of the fractions: > HI_RATIO at one record, not = HI_RATIO
    s = decide(*synthetic((-2.0, -2.1), (-2.0, -2.1), fN=(1e-2, 1e-2)))
    check('f_hi ratio = HI_RATIO: no flag',
          not s['spec']['flag'] and abs(s['spec']['frac']['hi']['max_ratio'] - 10) < 1e-12)
    s = decide(*synthetic((-2.0, -2.1), (-2.0, -2.1), fN=(1e-3, 0.1000001)))
    check('f_mid ratio > HI_RATIO: flag', s['spec']['flag'] and s['spec']['frac']['mid']['flag']
          and not s['spec']['frac']['hi']['flag'])
    # 4. the fractions of known modes on 32^3 (k_c = 10.67; shells 4: none, 6 and 8: mid, 12 and (8, 8, 0): hi)
    g = gen_ic.Grid(32, 32, 32)
    op = svv3d.SVV(g, 1.0)
    sh = shells(g)
    kc = (g.nlx//2)*2/3
    i, k, j = np.indices(op.shape)  # (x, z, y)
    th = lambda n, a: 2*np.pi*n*a/32  # noqa: E731
    amps = dict(x4=0.3, x6=0.2, y8=0.25, z12=0.1, xy8=0.15)
    f = (-1.0 + amps['x4']*np.cos(th(4, i)) + amps['x6']*np.cos(th(6, i)) + amps['y8']*np.sin(th(8, j))
         + amps['z12']*np.cos(th(12, k)) + amps['xy8']*np.cos(th(8, i) + th(8, j)))
    E = spectrum(op, f, sh)
    p = {q: a**2/2 for q, a in amps.items()}
    tot = sum(p.values())
    fr = band_fractions(E, kc)
    check('shells: sum over n >= 1 = variance', abs(np.sum(E[1:])/np.var(f) - 1) < 1e-12)
    check('fractions: hi = (z12 + xy8)/all, mid = (x6 + y8)/all',
          abs(fr['hi'] - (p['z12'] + p['xy8'])/tot) < 1e-12 and abs(fr['mid'] - (p['x6'] + p['y8'])/tot) < 1e-12)
    # 5. instab: H (max|k_d| > (32 - 1)//3 = 10: z12 and the Nyquist mode x16), nyq (x16; its variance is a^2), the
    # depth of a one-cell spike
    an = 0.05
    cf = cube_fractions(g, op, f + an*np.cos(th(16, i)), sh)
    check('instab: H = (z12 + x16)/all, nyq = x16/all',
          abs(cf['H'] - (p['z12'] + an**2)/(tot + an**2)) < 1e-12 and abs(cf['nyq'] - an**2/(tot + an**2)) < 1e-12)
    sp = 0.1*np.cos(th(1, i))
    q = (3, 5, 7)
    sp[q] -= 4.0
    cuts, depth = spike(sp, q)
    c = [np.cos(th(1, x)) for x in (2, 3, 4)]
    check('instab: spike depth = 4 + 0.1 ((c2 + c4 + 4 c3)/6 - c3), the cuts along x and y',
          abs(depth - (4.0 + 0.1*((c[0] + c[2] + 4*c[1])/6 - c[1]))) < 1e-12 and cuts['x'][L_CUT] == sp[q]
          and cuts['y'][L_CUT + 1] == sp[3, 5, 8] and cuts['z'][L_CUT - 1] == sp[3, 4, 7])
    # div_n on L = 2 pi: u = (0.3 sin 3x + 0.05 cos 16x, 0.2 sin 2y, 0.1 cos 5z); the Nyquist mode x16 drops out
    u = (0.3*np.sin(th(3, i)) + 0.05*np.cos(th(16, i)), 0.2*np.sin(th(2, j)), 0.1*np.cos(th(5, k)))
    dv = 0.9*np.cos(th(3, i)) + 0.4*np.cos(th(2, j)) - 0.5*np.sin(th(5, k))
    check('instab: div_n = 0.9 cos 3x + 0.4 cos 2y - 0.5 sin 5z', np.abs(div_n(g, op, u) - dv).max() < 1e-12)
    print('all PASS' if ok else 'some FAIL')
    return ok


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'selftest'
    if cmd == 'selftest':
        sys.exit(0 if selftest() else 1)
    fn = dict(verify=verify, rhomin=rhomin, fields=fields, judge=judge, instab=instab)[cmd]
    r = fn(sys.argv[2])
    if cmd == 'verify':
        sys.exit(0 if r['PASS'] else 1)
