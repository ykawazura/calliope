#!/usr/bin/env python3
"""Stage 1 pass/fail table from the batches of gen_stage1.py.

    python3 judge_stage1.py <base> [batch ...]      (base: e.g. runs/stage1_svv)

Runs analyze_stage1.analyze on every case of the batches found under <base> and checks
the criteria fixed in tasks/todo.md before the runs ("Stage 1 の合格基準"). Prints one
Markdown table per batch and writes <base>/judge.json.
"""
import sys
import os
import glob
import re
import json
import math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_stage1 import analyze, read_lines, run_dirs

CRIT = dict(dMdt_raw=1e-12, dP=1e-11, R_zeta=1e-10, R_svv=1e-10, D_svv_min=0.0, sound=1e-5, star=1e-2,
            speed=1e-2, pos_cells=2.0, overshoot=0.05, L1_order=0.8, simple_order=1.8, simple_floor=1e-9,
            mass_dt_order=2.7, decomp=1e-10, dt_star=1e-3)
# the shock criteria apply for N >= N_SHOCK (U4: the star state rho* = e^-4 spans 33 cells at N = 512)
N_SHOCK = {'U4': 1024}


def ok(flag):
    return 'pass' if flag else '**FAIL**'


def orders(Ns, errs):
    return [math.log(errs[i] / errs[i + 1]) / math.log(Ns[i + 1] / Ns[i]) for i in range(len(Ns) - 1)]


def common(o):
    """Checks that apply to every run: completion, mass, momentum, R_zeta."""
    r = dict(exit=o.get('exit'))
    r['completed'] = (o.get('exit') == 'exit 0')
    b, c = o.get('budget'), o.get('cons')
    if b:
        r['dMdt_raw'] = b['max_abs_dMdt_raw_over_abs']
        r['R_zeta'] = b['max_R_zeta_ratio']
        r['int_eps_svv'] = b.get('int_eps_svv', 0.0)
        r['int_eps_shock'] = b['int_eps_shock']
        if 'min_D_svv' in b:
            r['min_D_svv'] = b['min_D_svv']
            r['int_D_svv'] = b['int_D_svv']
            # |R_svv| relative to its rounding scale A_svv = sum |dE/dU.U_svv| dV (|R_svv|/D_svv when
            # A_svv is not written); |R_svv|/D_svv is kept for the report
            r['R_svv'] = b.get('max_R_svv_scale', b['max_R_svv_ratio'])
            r['R_svv_over_D'] = b['max_R_svv_ratio']
    if c:
        r['dP'] = c['max_abs_dP_over_scale']
        r['dM'] = c['max_abs_dM_over_M0']
        r['rho_min'] = c['rho_min']
    v = o.get('svv', {})
    r['es'] = v.get('form') == 'entropy' and v.get('eps0', 0.0) > 0.0
    r['es_identity'] = r['es'] and not v.get('filter_lnrho')
    return r


def checks_common(r):
    """completion, mass, momentum, R_zeta; R_svv and D_svv >= 0 for the entropy-stable SVV"""
    c = [r['completed'], r.get('dMdt_raw', 1) <= CRIT['dMdt_raw'], r.get('dP', 1) <= CRIT['dP'],
         r.get('R_zeta', 1) <= CRIT['R_zeta']]
    if r['es']:
        c.append(r.get('min_D_svv', -1) >= CRIT['D_svv_min'])
    if r['es_identity']:
        c.append(r.get('R_svv', 1) <= CRIT['R_svv'])
    return c


def judge_shock(cases):
    rows, by_problem = [], {}
    for case in cases:
        o = analyze(case)
        r = common(o)
        r['case'] = os.path.basename(case)
        f = o.get('final')
        if f:
            r['L1'] = f['L1_rho']
            r['t'] = f['t']
            r['valid'] = f['valid']
            r['star'] = max(abs(s['rho_s_rel_err']) for s in o['stars'] if 'rho_s' in s)
            sh = o['shocks']
            r['speed'] = max((abs(m['speed_rel_err']) for m in sh if 'speed_rel_err' in m), default=float('nan'))
            r['pos'] = max((abs(m['dx_cells']) for m in sh), default=float('nan'))
            r['over'] = max((max(m['overshoot'], m['undershoot']) for m in sh), default=float('nan'))
        name, N = r['case'].split('_N')[0], int(r['case'].split('_N')[1].split('_')[0])
        r['N'], r['problem'] = N, name
        by_problem.setdefault(name, []).append(r)
        rows.append(r)
    return rows, by_problem


def fmt(x, f='%.2e'):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return '-'
    return f % x


def table_shock(rows, conv=True):
    out = ['| case | exit | ρ_min | L1(ρ) | ρ* err | speed err | pos (cells) | over/under | dMdt_raw | dP | R_ζ '
           '| R_svv | min D_svv | ∫ε_svv | 判定 |',
           '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        checks = checks_common(r)
        shock = r['N'] >= N_SHOCK.get(r['problem'], 512) and 'L1' in r
        if shock:
            checks += [r['star'] <= CRIT['star'], r['speed'] <= CRIT['speed'], r['pos'] <= CRIT['pos_cells'],
                       r['over'] <= CRIT['overshoot']]
        r['pass'] = all(checks)
        verdict = ok(r['pass']) if shock or not r['pass'] else 'pass（衝撃波の項目は参考）'
        out.append('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
            r['case'], r['exit'], fmt(r.get('rho_min'), '%.4f'), fmt(r.get('L1'), '%.4f'), fmt(r.get('star')),
            fmt(r.get('speed')), fmt(r.get('pos'), '%.2f'), fmt(r.get('over'), '%.3f'), fmt(r.get('dMdt_raw')),
            fmt(r.get('dP')), fmt(r.get('R_zeta')), fmt(r.get('R_svv')) if r['es_identity'] else '-',
            fmt(r.get('min_D_svv')) if r['es'] else '-', fmt(r.get('int_eps_svv')), verdict))
    return out


def judge_batch_shock(cases):
    rows, byp = judge_shock(cases)
    out = table_shock(rows)
    out.append('')
    out.append('| problem | N | L1(ρ) | observed order (N → 2N) |')
    out.append('|---|---|---|---|')
    res = {}
    for name, rs in byp.items():
        rs = sorted([r for r in rs if 'L1' in r], key=lambda r: r['N'])
        Ns, L1 = [r['N'] for r in rs], [r['L1'] for r in rs]
        od = orders(Ns, L1) if len(Ns) > 1 else []
        sel = [i for i in range(len(Ns)) if Ns[i] >= 512]
        fit = -np.polyfit(np.log([Ns[i] for i in sel]), np.log([L1[i] for i in sel]), 1)[0] if len(sel) > 1 else float('nan')
        res[name] = dict(N=Ns, L1=L1, orders=od, fit_order_512=fit, fit_pass=bool(fit >= CRIT['L1_order']))
        for i, r in enumerate(rs):
            out.append('| %s | %d | %.4e | %s |' % (name, r['N'], r['L1'], fmt(od[i], '%.3f') if i < len(od) else '-'))
        out.append('| %s | fit N ≥ 512 | | %.3f (%s) |' % (name, fit, ok(res[name]['fit_pass'])))
    return out, dict(rows=rows, convergence=res)


def es_summary(rows):
    """R_svv and min D_svv over the rows of a batch (entropy-stable SVV)"""
    es = [r for r in rows if r['es']]
    if not es:
        return []
    idn = [r for r in es if r['es_identity']]
    Dmin = min(r.get('min_D_svv', -1.0) for r in es)
    s = ['', 'ES-SVV: min D_svv = %.2e over %d runs（基準 ≥ 0: %s）' % (Dmin, len(es), ok(Dmin >= CRIT['D_svv_min']))]
    if idn:
        R = max(r.get('R_svv', 1.0) for r in idn)
        RD = max(r.get('R_svv_over_D', 1.0) for r in idn)
        s.append('ES-SVV: max |R_svv|/A_svv = %.2e over %d runs with filter_lnrho = F（基準 ≤ %.0e: %s）、'
                 'max |R_svv|/D_svv = %.2e（参考）' % (R, len(idn), CRIT['R_svv'], ok(R <= CRIT['R_svv']), RD))
    return s


def judge_batch_scan(cases):
    rows, _ = judge_shock(cases)
    return table_shock(rows), dict(rows=rows)


def judge_batch_smooth(cases):
    out = ['| case | exit | 量 | 値 | 判定 |', '|---|---|---|---|---|']
    simple = {}
    res = dict(rows=[])
    for case in cases:
        o = analyze(case)
        r = common(o)
        name = os.path.basename(case)
        if 'sound' in o:
            w = o['sound']
            r.update(amp=w['max_abs_amp_err'], phase=w['max_abs_phase_err'])
            crit = name.endswith('C0')
            flag = r['completed'] and (not crit or (w['max_abs_amp_err'] <= CRIT['sound'] and w['max_abs_phase_err'] <= CRIT['sound']))
            out.append('| %s | %s | max\\|amp err\\|, max\\|phase err\\| | %.2e, %.2e | %s |'
                       % (name, r['exit'], w['max_abs_amp_err'], w['max_abs_phase_err'], ok(flag) if crit else '参考'))
        if 'simple' in o:
            fin = o['simple']['final']
            r.update(Linf_s=fin['Linf_s'], Linf_u=fin['Linf_u'], t_over_tb=fin['t_over_tb'])
            tag = name.split('_', 2)[2]
            N = int(name.split('_N')[1].split('_')[0])
            simple.setdefault(tag, []).append((N, fin['Linf_s'], fin['Linf_u']))
        r['case'] = name
        res['rows'].append(r)
    out.append('')
    # C_shock = 0: Linf(u) drops by >= 10 per doubling of N until the floor (rounding, time error) is reached;
    # C_shock = 2: observed order >= 1.8 from N = 64 (zeta ~ dl^2; N = 32 -> 64 is reported)
    out.append('| simple wave | N | Linf(s) | Linf(u) | order (u) | 判定 (N → 2N) |')
    out.append('|---|---|---|---|---|---|')
    res['simple'] = {}
    for tag, v in sorted(simple.items()):
        v.sort()
        Ns, es, eu = [a[0] for a in v], [a[1] for a in v], [a[2] for a in v]
        od = orders(Ns, eu)
        flags = []
        for i in range(len(od)):
            if tag.startswith('C0'):
                flags.append(eu[i] / eu[i + 1] >= 10.0 or eu[i] <= CRIT['simple_floor'])
            else:
                flags.append(od[i] >= CRIT['simple_order'] if Ns[i] >= 64 else None)
        res['simple'][tag] = dict(N=Ns, Linf_s=es, Linf_u=eu, orders=od, pass_=[f for f in flags])
        for i in range(len(Ns)):
            f = flags[i] if i < len(flags) else None
            out.append('| %s | %d | %.3e | %.3e | %s | %s |' % (tag, Ns[i], es[i], eu[i],
                       fmt(od[i], '%.2f') if i < len(od) else '-', '-' if i >= len(flags) else ('参考' if f is None else ok(f))))
    out += es_summary(res['rows'])
    return out, res


def judge_batch_dtconv(cases):
    series = sorted([c for c in cases if '_dt' in os.path.basename(c)], key=lambda c: int(c.split('_dt')[-1]))
    pairs = [c for c in cases if os.path.basename(c).endswith('_cflsvv1')]
    out, res = judge_dt_series(series) if series else ([], {})
    for c in pairs:
        o, p = judge_cflsvv_pair(c)
        out += [''] + o
        res[os.path.basename(c)] = p
    return out, res


def judge_cflsvv_pair(case):
    """case (cfl_svv = 1) against shock/<problem>_N<N> (cfl_svv = 2): the dt of the ES-SVV limits the error"""
    ref = os.path.join(os.path.dirname(os.path.dirname(case)), 'shock', os.path.basename(case).rsplit('_cflsvv1', 1)[0])
    o1, o2 = analyze(case), analyze(ref)
    r1, r2 = common(o1), common(o2)
    d = dict(case=case, ref=ref, exit=o1.get('exit'), exit_ref=o2.get('exit'))
    ok_run = r1['completed'] and r2['completed'] and o1.get('final') and o2.get('final')
    if ok_run:
        d['d_star'] = max(abs(a['rho_s_rel_err'] - b['rho_s_rel_err'])
                          for a, b in zip(o1['stars'], o2['stars']) if 'rho_s' in a and 'rho_s' in b)
        l1, l2 = read_lines(case)[-1], read_lines(ref)[-1]
        d['t'] = (l1[0], l2[0])
        d['d_rho_max'] = float(np.max(np.abs(l1[2][:, 1] - l2[2][:, 1])))
        d['L1'] = (o1['final']['L1_rho'], o2['final']['L1_rho'])
        d['star_err'] = (max(abs(a['rho_s_rel_err']) for a in o1['stars'] if 'rho_s' in a),
                         max(abs(a['rho_s_rel_err']) for a in o2['stars'] if 'rho_s' in a))
    d['pass'] = bool(ok_run and d['d_star'] <= CRIT['dt_star'])
    out = ['| case | exit | dt | t_end | ρ* err | L1(ρ) | 判定 |', '|---|---|---|---|---|---|---|']
    for o, i in ((o1, 0), (o2, 1)):
        out.append('| %s | %s | %.4e | %s | %s | %s | |' % (os.path.relpath(o['case'], os.path.dirname(os.path.dirname(case))),
                   o.get('exit'), o['cfl']['dt_min'], fmt(o['cons']['t_end'], '%.4f'),
                   fmt(d['star_err'][i]) if ok_run else '-', fmt(d['L1'][i], '%.4f') if ok_run else '-'))
    out.append('| 差 | | | | Δ(ρ*/ρ*_exact) = %s | max\\|Δρ\\| = %s | %s（基準 ≤ %.0e） |'
               % (fmt(d.get('d_star')), fmt(d.get('d_rho_max')), ok(d['pass']), CRIT['dt_star']))
    return out, d


def judge_dt_series(cases):
    out = ['| case | dt | max\\|ΔM/M0\\| | order | max\\|ρ − ρ_finest\\| | order |', '|---|---|---|---|---|---|']
    data = []
    for case in cases:
        o = analyze(case)
        lines = read_lines(case)
        data.append((float(o['cfl']['dt_min']), o['cons']['max_abs_dM_over_M0'], lines[-1][2][:, 1], o))
    ref = data[-1][2]
    dts = [d[0] for d in data]
    dM = [d[1] for d in data]
    dr = [float(np.max(np.abs(d[2] - ref))) for d in data]
    oM = [math.log(dM[i] / dM[i + 1]) / math.log(dts[i] / dts[i + 1]) for i in range(len(dts) - 1)]
    orr = [math.log(dr[i] / dr[i + 1]) / math.log(dts[i] / dts[i + 1]) for i in range(len(dts) - 2)]
    for i, d in enumerate(data):
        out.append('| %s | %.4e | %.3e | %s | %s | %s |' % (os.path.basename(cases[i]),
                   dts[i], dM[i], fmt(oM[i], '%.2f') if i < len(oM) else '-', fmt(dr[i]) if i < len(data) - 1 else '(ref)',
                   fmt(orr[i], '%.2f') if i < len(orr) else '-'))
    out.append('')
    out.append('質量誤差の次数（最小）: %.2f（基準 ≥ %.1f: %s）' % (min(oM), CRIT['mass_dt_order'], ok(min(oM) >= CRIT['mass_dt_order'])))
    out += es_summary([common(d[3]) for d in data])
    return out, dict(dt=dts, dM=dM, d_rho=dr, order_M=oM, order_rho=orr, pass_=bool(min(oM) >= CRIT['mass_dt_order']))


def judge_batch_decomp(cases):
    out = ['| case | exit | 出力数 | max\\|Δρ\\| vs x_np1 | max\\|Δu\\| vs x_np1 | 判定 |', '|---|---|---|---|---|---|']
    ref_case = [c for c in cases if c.endswith('_x_np1')][0]
    ref = read_lines(ref_case)
    res = []
    for case in cases:
        o = analyze(case)
        d = os.path.basename(case).split('_')[2]
        idx = 'xyz'.index(d)
        lines = read_lines(case)
        dr = du = 0.0
        for (t0, s0, a0), (t1, s1, a1) in zip(ref, lines):
            assert s0 == s1, (case, s0, s1)
            dr = max(dr, float(np.max(np.abs(a1[:, 1] - a0[:, 1]))))
            du = max(du, float(np.max(np.abs(a1[:, 2 + idx] - a0[:, 2]))))
        flag = o.get('exit') == 'exit 0' and len(lines) == len(ref) and max(dr, du) <= CRIT['decomp']
        out.append('| %s | %s | %d | %.2e | %.2e | %s |' % (os.path.basename(case), o.get('exit'), len(lines), dr, du, ok(flag)))
        res.append(dict(case=case, d_rho=dr, d_u=du, n=len(lines), pass_=bool(flag), common=common(o)))
    out += es_summary([r['common'] for r in res])
    return out, res


def estimate_only(case):
    """True if every run of the case logged its FFT plans (fftw_shim.c) and all of them use FFTW_ESTIMATE."""
    flags = []
    for d in run_dirs(case):
        p = os.path.join(d, 'calliope.out.std')
        if not os.path.exists(p):
            return False
        f = [int(x) for x in re.findall(r'^fftw_shim: .* -> (\d+)$', open(p, errors='replace').read(), re.M)]
        if not f:
            return False
        flags += f
    return bool(flags) and all(x & 64 for x in flags)


def judge_batch_smoke(cases):
    """NaN injection stops with exit != 0; zeta of the self test bit-identical for 1, 2, 4 ranks;
    the colliding flows pass the common checks; the run continued from restart/ after one step
    (collide_U2_N512_est_chain) is bit-identical to collide_U2_N512_est, both with FFTW_ESTIMATE."""
    out = ['| case | exit | 検査 | 判定 |', '|---|---|---|---|']
    res = []
    byname = {os.path.basename(c): c for c in cases}
    for case in cases:
        name = os.path.basename(case)
        o = analyze(case)
        r = dict(case=case, exit=o.get('exit'))
        if name == 'nan':
            txt = open(os.path.join(case, 'calliope.out.std'), errors='replace').read()
            flag = o.get('exit') not in (None, 'exit 0') and 'non-finite value' in txt
            what = 'exit != 0 と non-finite の停止メッセージ'
        elif name.startswith('selftest_np'):
            z = [np.loadtxt(os.path.join(byname[n], 'shock_selftest.dat')) for n in ('selftest_np1', name)]
            r['max_diff'] = float(np.max(np.abs(z[1] - z[0]))) if z[0].shape == z[1].shape else float('inf')
            flag = o.get('exit') == 'exit 0' and r['max_diff'] == 0.0
            what = 'ζ の max|差| vs np1 = %.1e' % r['max_diff']
        elif name.endswith('_chain'):
            ref = byname[name.rsplit('_chain', 1)[0]]
            # the continued run also writes out1d at its start; compare at the output times of ref
            l1 = {t: a for t, _, a in read_lines(case)}
            l0 = read_lines(ref)
            o0 = analyze(ref)
            r['runs'] = o.get('runs')
            r['n_cons'] = (o['cons']['n'], o0['cons']['n'])
            r['max_diff'] = max(float(np.max(np.abs(l1[t] - a0))) if t in l1 and l1[t].shape == a0.shape
                                else float('inf') for t, _, a0 in l0)
            # both runs with FFTW_ESTIMATE in every plan (fftw_shim.c logs the flags of each plan)
            r['estimate'] = all(estimate_only(c) for c in (ref, case))
            flag = (o.get('exit') == 'exit 0' and r['runs'] == 2 and r['n_cons'][0] == r['n_cons'][1]
                    and r['max_diff'] == 0.0 and o['cons']['max_abs_dM_over_M0'] == o0['cons']['max_abs_dM_over_M0']
                    and r['estimate'])
            what = '%d 回の実行、1D 出力 %d 時刻の max|差| vs %s = %.1e、conservation.dat の行数 %d/%d、FFTW_ESTIMATE %s' % (
                r['runs'], len(l0), os.path.basename(ref), r['max_diff'], r['n_cons'][0], r['n_cons'][1],
                '両方' if r['estimate'] else 'なし')
        elif name.endswith('_nosvv_fs'):
            # configuration of the first smoke run (C = 1, no SVV, filtered ln rho): reproduces its stop
            r.update(common(o))
            flag = None
            what = '初回 smoke の設定（C_shock = 1、SVV なし、ln ρ にフィルタ）の再現、t_end = %.4f、ρ_min = %.2e' % (
                o['cons']['t_end'], o['cons']['rho_min'])
        else:
            r.update(common(o))
            flag = all(checks_common(r))
            what = '共通検査（質量・運動量・R_ζ%s）' % ('・R_svv・D_svv' if r['es'] else '')
        r['pass_'] = flag if flag is None else bool(flag)
        out.append('| %s | %s | %s | %s |' % (name, o.get('exit'), what, '参考' if flag is None else ok(flag)))
        res.append(r)
    return out, res


def judge_batch_cflmax(cases):
    out = ['| case | dt | exit | t_end | max dt/limit | max dt·λ_svv | dt·λ_svv（検査ステップ） | ρ_min | L1(ρ) |',
           '|---|---|---|---|---|---|---|---|---|']
    res = []
    for case in sorted(cases):
        o = analyze(case)
        c, f = o.get('cfl', {}), o.get('final') or {}
        out.append('| %s | %.3e | %s | %s | %s | %s | %s | %s | %s |' % (os.path.basename(case), c.get('dt_min', float('nan')),
                   o.get('exit'), fmt(o.get('cons', {}).get('t_end'), '%.4f'), fmt(c.get('max_dt_over_limit'), '%.3f'),
                   fmt(c.get('max_dt_lam'), '%.3f'), fmt(c.get('max_dt_lam_check'), '%.3f'),
                   fmt(o.get('cons', {}).get('rho_min'), '%.4f'), fmt(f.get('L1_rho'), '%.4f')))
        res.append(dict(case=case, exit=o.get('exit'), cfl=c, final=f))
    return out, res


JUDGE = {'smoke': judge_batch_smoke, 'shock': judge_batch_shock, 'scan': judge_batch_scan, 'smooth': judge_batch_smooth,
         'dtconv': judge_batch_dtconv, 'decomp': judge_batch_decomp, 'cflmax': judge_batch_cflmax}

if __name__ == '__main__':
    base = sys.argv[1]
    batches = sys.argv[2:] or [b for b in JUDGE if os.path.isdir(os.path.join(base, b))]
    allres = {}
    for b in batches:
        cases = sorted(d for d in glob.glob(os.path.join(base, b, '*')) if os.path.isdir(d))
        if not cases:
            continue
        lines, res = JUDGE[b](cases)
        allres[b] = res
        print('### %s\n' % b)
        print('\n'.join(lines))
        print()
    with open(os.path.join(base, 'judge.json'), 'w') as f:
        json.dump(allres, f, indent=1, default=float)
