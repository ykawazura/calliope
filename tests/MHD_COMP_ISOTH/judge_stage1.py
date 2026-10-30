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
import json
import math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_stage1 import analyze, read_lines

CRIT = dict(dMdt_raw=1e-12, dP=1e-11, R_zeta=1e-10, sound=1e-5, star=1e-2, speed=1e-2, pos_cells=2.0,
            overshoot=0.05, L1_order=0.8, simple_order=1.8, mass_dt_order=2.7, decomp=1e-10)


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
    if c:
        r['dP'] = c['max_abs_dP_over_scale']
        r['dM'] = c['max_abs_dM_over_M0']
        r['rho_min'] = c['rho_min']
    return r


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
        r['N'] = N
        by_problem.setdefault(name, []).append(r)
        rows.append(r)
    return rows, by_problem


def fmt(x, f='%.2e'):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return '-'
    return f % x


def table_shock(rows, conv=True):
    out = ['| case | exit | ρ_min | L1(ρ) | ρ* err | speed err | pos (cells) | over/under | dMdt_raw | dP | R_ζ | ∫ε_svv | 判定 |',
           '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        checks = [r['completed'], r.get('dMdt_raw', 1) <= CRIT['dMdt_raw'], r.get('dP', 1) <= CRIT['dP'],
                  r.get('R_zeta', 1) <= CRIT['R_zeta']]
        if r['N'] >= 512 and 'L1' in r:
            checks += [r['star'] <= CRIT['star'], r['speed'] <= CRIT['speed'], r['pos'] <= CRIT['pos_cells'],
                       r['over'] <= CRIT['overshoot']]
        out.append('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
            r['case'], r['exit'], fmt(r.get('rho_min'), '%.4f'), fmt(r.get('L1'), '%.4f'), fmt(r.get('star')),
            fmt(r.get('speed')), fmt(r.get('pos'), '%.2f'), fmt(r.get('over'), '%.3f'), fmt(r.get('dMdt_raw')),
            fmt(r.get('dP')), fmt(r.get('R_zeta')), fmt(r.get('int_eps_svv')), ok(all(checks))))
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
        res[name] = dict(N=Ns, L1=L1, orders=od, fit_order_512=fit)
        for i, r in enumerate(rs):
            out.append('| %s | %d | %.4e | %s |' % (name, r['N'], r['L1'], fmt(od[i], '%.3f') if i < len(od) else '-'))
        out.append('| %s | fit N ≥ 512 | | %.3f (%s) |' % (name, fit, ok(fit >= CRIT['L1_order'])))
    return out, dict(rows=rows, convergence=res)


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
    out.append('| simple wave | N | Linf(s) | Linf(u) | order (u) |')
    out.append('|---|---|---|---|---|')
    res['simple'] = {}
    for tag, v in sorted(simple.items()):
        v.sort()
        Ns, es, eu = [a[0] for a in v], [a[1] for a in v], [a[2] for a in v]
        od = orders(Ns, eu)
        res['simple'][tag] = dict(N=Ns, Linf_s=es, Linf_u=eu, orders=od)
        for i in range(len(Ns)):
            out.append('| %s | %d | %.3e | %.3e | %s |' % (tag, Ns[i], es[i], eu[i], fmt(od[i], '%.2f') if i < len(od) else '-'))
    return out, res


def judge_batch_dtconv(cases):
    out = ['| case | dt | max\\|ΔM/M0\\| | order | max\\|ρ − ρ_finest\\| | order |', '|---|---|---|---|---|---|']
    data = []
    for case in sorted(cases, key=lambda c: int(c.split('_dt')[-1])):
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
        out.append('| %s | %.4e | %.3e | %s | %s | %s |' % (os.path.basename(cases[0]).rsplit('_dt', 1)[0] + '_dt%d' % 2 ** i,
                   dts[i], dM[i], fmt(oM[i], '%.2f') if i < len(oM) else '-', fmt(dr[i]) if i < len(data) - 1 else '(ref)',
                   fmt(orr[i], '%.2f') if i < len(orr) else '-'))
    out.append('')
    out.append('質量誤差の次数（最小）: %.2f（基準 ≥ %.1f: %s）' % (min(oM), CRIT['mass_dt_order'], ok(min(oM) >= CRIT['mass_dt_order'])))
    return out, dict(dt=dts, dM=dM, d_rho=dr, order_M=oM, order_rho=orr)


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
        res.append(dict(case=case, d_rho=dr, d_u=du, n=len(lines)))
    return out, res


def judge_batch_cflmax(cases):
    out = ['| case | dt | exit | t_end | max dt/limit | ρ_min | L1(ρ) |', '|---|---|---|---|---|---|---|']
    res = []
    for case in sorted(cases):
        o = analyze(case)
        c, f = o.get('cfl', {}), o.get('final') or {}
        out.append('| %s | %.3e | %s | %s | %s | %s | %s |' % (os.path.basename(case), c.get('dt_min', float('nan')), o.get('exit'),
                   fmt(o.get('cons', {}).get('t_end'), '%.4f'), fmt(c.get('max_dt_over_limit'), '%.3f'),
                   fmt(o.get('cons', {}).get('rho_min'), '%.4f'), fmt(f.get('L1_rho'), '%.4f')))
        res.append(dict(case=case, exit=o.get('exit'), cfl=c, final=f))
    return out, res


JUDGE = {'shock': judge_batch_shock, 'scan': judge_batch_scan, 'smooth': judge_batch_smooth,
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
