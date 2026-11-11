#!/usr/bin/env python3
"""Analysis of the Stage 3 Fortran runs (cases of gen_stage3.py), in the format of the prototype runs of stage3.py.

    python3 analyze_stage3.py case <case> ...        summary of each case (stage3.json in the case directory)
    python3 analyze_stage3.py table <batch dir> <out> <name>.json/.npz of each case in <out> (the format of
                                                      runs/stage3/), and <out>/tables.md, <out>/fig by stage3_table
    python3 analyze_stage3.py agree <batch dir>       the cases against the prototype runs runs/stage3/<case>.npz
    python3 analyze_stage3.py decomp <batch dir>      the decompositions and directions of each case against _x_np1
    python3 analyze_stage3.py last <batch dir> <out>  <out>/<case>.npz (t, W) of the last line of each case (conv)

The snapshots are the out1d lines at the steps of the prototype snapshots, rotated back to (rho, u_n, u_t1, u_t2, B_t1,
B_t2) (gen_stage3.rotate), with eta from column 11 of the lines (C_eta > 0; 0 otherwise). The statistics are those of
stage3.Stats from conservation.dat, budget.dat and cfl.dat, per step: drift_* from conservation.dat (the drift of the
integrals of B from the out1d lines), rate_dM_raw from budget.dat (no dP/dt or dB/dt columns there), R_zeta, R_svv,
R_eta relative to -visc_diss, D_svv and eta_diss, the energy changes from the stage-weighted rates of budget.dat,
eta_max from cfl.dat.
"""
import sys
import os
import json
import glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_stage1 import load_runs, read_lines, analyze as analyze1
import stage3 as s3


def meta(case):
    with open(os.path.join(case, 'case.json')) as f:
        return json.load(f)


def snaps_of(case, m):
    """[(t, W (6, N), zeta, eta, D)] of the out1d lines at the steps of the prototype snapshots (multiples of
    nstep // NWRITE and nstep: a restart writes one more line at its first step), the steps of the lines (the step
    of the run that wrote the line; equal to the step from t = 0 for a single run) and max |B_n - B_x| over them."""
    c, i = m['case'], 'xyz'.index(m['dir'])
    x = s3.grid(c['N'])
    m_w = m['nstep'] // s3.NWRITE
    out, steps, dbn = [], [], 0.0
    for t, step, a in read_lines(case):
        n = int(round(t / m['dt']))
        if n % m_w != 0 and n != m['nstep']:
            continue
        ncol = 11 if c['C_eta'] > 0.0 else 10
        assert a.shape == (c['N'], ncol) and np.max(np.abs(a[:, 0] - x)) < 1e-12, (case, t, a.shape)
        u, b = a[:, 2:5], a[:, 5:8]
        W = np.array([a[:, 1], u[:, i], u[:, (i + 1) % 3], u[:, (i + 2) % 3], b[:, (i + 1) % 3], b[:, (i + 2) % 3]])
        out.append((t, W, a[:, 8], a[:, 10] if ncol == 11 else np.zeros(c['N']), a[:, 9]))
        steps.append(step)
        dbn = max(dbn, float(np.max(np.abs(b[:, i] - s3.bx_of(c)))))
    return out, steps, dbn


def stats(case, m, snaps):
    """stage3.Stats keys from the per-step files of the run (None where the files have no counterpart)."""
    d = {}
    cons, bud, cfl = (np.concatenate(r) if r else None
                      for r in (load_runs(case, n) for n in ('conservation.dat', 'budget.dat', 'cfl.dat')))
    i = 'xyz'.index(m['dir'])
    if cons is not None:
        M0 = cons[0, 2]
        d.update(M0=M0, drift_M=float(np.max(np.abs(cons[:, 2] - M0)) / M0),
                 rho_min=float(cons[:, 11].min()), rho_max=float(cons[:, 12].max()),
                 zeta_max=float(cons[:, 13].max()), max_E_kin=float(cons[:, 7].max()))
        for k, j in (('Px', i), ('Py', (i + 1) % 3), ('Pz', (i + 2) % 3)):
            d['drift_' + k] = float(np.max(np.abs(cons[:, 4 + j] - cons[0, 4 + j])))
    dl = s3.L / m['case']['N']
    for k, j in (('By', 4), ('Bz', 5)):
        B = [float(np.sum(W[j]) * dl) for t, W, *_ in snaps]
        d['drift_' + k] = float(np.max(np.abs(np.array(B) - B[0]))) if B else None
    if bud is not None:
        dt = bud[:, 1]
        pos = lambda a: a > 0.0
        visc, Dsvv, etad = -bud[:, 11], -bud[:, 22], -bud[:, 27]
        d.update(E0=float(bud[0, 18]), E_end=float(bud[-1, 19]), dE=float(bud[-1, 19] - bud[0, 18]),
                 n_E_up=int(np.sum(bud[:, 19] > bud[:, 18])),
                 max_E_up=float(np.max(bud[:, 19] - bud[:, 18])) if np.any(bud[:, 19] > bud[:, 18]) else None,
                 rate_dM_raw=float(np.max(np.abs(bud[:, 3])[pos(bud[:, 5])] / bud[pos(bud[:, 5]), 5]))
                 if np.any(pos(bud[:, 5])) else None,
                 R_zeta=float(np.max(np.abs(bud[pos(visc), 12]) / visc[pos(visc)])) if np.any(pos(visc)) else None,
                 R_svv=float(np.max(np.abs(bud[Dsvv != 0, 23]) / np.abs(Dsvv[Dsvv != 0])))
                 if np.any(Dsvv != 0) else None,
                 R_eta=float(np.max(np.abs(bud[pos(etad), 28]) / etad[pos(etad)])) if np.any(pos(etad)) else None,
                 min_D_svv=float(np.min(Dsvv[Dsvv != 0])) if np.any(Dsvv != 0) else None,
                 n_svv_pos=int(np.sum(bud[:, 20] - bud[:, 21] > 0.0)),
                 max_svv_b_eps=float(np.max(bud[:, 21])), n_svv_b_pos=int(np.sum(bud[:, 21] > 0.0)),
                 int_rk_ideal_eps=float(np.sum(dt * bud[:, 9])), int_rk_visc_eps=float(np.sum(dt * bud[:, 10])),
                 int_rk_svv_eps=float(np.sum(dt * (bud[:, 20] - bud[:, 21]))),
                 int_rk_svv_b_eps=float(np.sum(dt * bud[:, 21])), int_rk_eta_eps=float(np.sum(dt * bud[:, 26])),
                 int_rk_res_eps=float(np.sum(dt * bud[:, 15])), int_rk_nu_eps=float(np.sum(dt * bud[:, 16])),
                 int_rk_eps=float(np.sum(dt * bud[:, 7])), time_err=float(np.sum(dt * bud[:, 8])),
                 int_rk_visc=float(np.sum(dt * visc)), int_rk_D_svv=float(np.sum(dt * Dsvv)),
                 int_rk_eta_diss=float(np.sum(dt * etad)))
    if cfl is not None:
        d.update(max_dt_over_limit=float(np.max(cfl[:, 1] / np.minimum(cfl[:, 2], cfl[:, 3]))),
                 max_dt_lam_svv=float(np.max(cfl[:, 14])), dt_min=float(cfl[:, 1].min()), dt_max=float(cfl[:, 1].max()),
                 eta_max=float(cfl[:, 16].max()))
    return d


def analyze_case(case):
    """stage3.json of the case: the stage3 record (case, t_end, attempts, ok, analysis) and the Fortran summary."""
    m = meta(case)
    c = m['case']
    s1 = analyze1(case)                     # exit status, conservation/budget/cfl summaries (summary.json)
    snaps, steps, dbn = snaps_of(case, m)
    m_w = m['nstep'] // s3.NWRITE
    t_stop = snaps[-1][0] if snaps else 0.0
    ok = (bool(snaps) and abs(t_stop - m['t_end']) < 1e-9 * m['t_end'] and s1.get('exit') == 'exit 0')
    status = 'ok' if ok else 'stopped (%s)' % s1.get('exit')
    st = stats(case, m, snaps)
    with np.errstate(all='ignore'):
        an = s3.analyze(c, snaps, s3.initial(c)) if snaps else {}
    att = dict(dt_factor=1, dt=m['dt'], nstep=m['nstep'], status=status, t_stop=t_stop,
               steps=int(round(t_stop / m['dt'])), stats=st)
    out = dict(case=c, t_end=m['t_end'], bx=s3.bx_of(c), dir=m['dir'], attempts=[att], ok=ok, analysis=an,
               fortran=dict(exit=s1.get('exit'), runs=s1.get('runs'), line_steps=steps,
                            line_steps_expected=list(range(0, m['nstep'], m_w)) + [m['nstep']],
                            max_dBn=dbn, cons=s1.get('cons'), budget=s1.get('budget'), cfl=s1.get('cfl')))
    with open(os.path.join(case, 'stage3.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)
    return out, snaps


def cases_in(batch):
    return sorted(os.path.dirname(f) for f in glob.glob(os.path.join(batch, '*', 'case.json')))


def table(batch, out):
    """The cases of the batch as runs/stage3-style records in out, and the tables and figures of stage3_table."""
    import stage3_table
    os.makedirs(out, exist_ok=True)
    for case in cases_in(batch):
        d, snaps = analyze_case(case)
        if not snaps:
            continue
        name = d['case']['name']
        np.savez_compressed(os.path.join(out, name + '.npz'), t=np.array([s[0] for s in snaps]),
                            W=np.array([s[1] for s in snaps]), zeta=np.array([s[2] for s in snaps]),
                            eta=np.array([s[3] for s in snaps]), W0=s3.initial(d['case']))
        with open(os.path.join(out, name + '.json'), 'w') as f:
            json.dump(d, f, indent=1, default=float)
    stage3_table.main(out, figdir=os.path.join(out, 'fig'))


def last_states(batch, out):
    """<out>/<case>.npz (t, W) of the last out1d line of each case of the batch that reached t_end (conv: no prototype
    runs)."""
    os.makedirs(out, exist_ok=True)
    for case in cases_in(batch):
        m = meta(case)
        snaps = snaps_of(case, m)[0]
        if not snaps or abs(snaps[-1][0] - m['t_end']) > 1e-9 * m['t_end']:
            print('%s: not at t_end (t = %s), skipped' % (os.path.basename(case), snaps[-1][0] if snaps else None))
            continue
        np.savez_compressed(os.path.join(out, os.path.basename(case) + '.npz'), t=snaps[-1][0], W=snaps[-1][1])


def maxdiff(Wa, Wb):
    """max over the six variables of max |a - b| / max(1, max |b|)."""
    return max(float(np.max(np.abs(Wa[i] - Wb[i])) / max(1.0, float(np.max(np.abs(Wb[i]))))) for i in range(6))


def agree(batch):
    """The Fortran lines against the prototype snapshots of runs/stage3/<case>.npz at the same steps."""
    print('%-18s %6s %6s %10s %10s %10s %10s %10s' % ('case', 'lines', 'steps', 'max|dt|', 'diff t0', 'diff max',
                                                     'diff end', 'zeta end'))
    for case in cases_in(batch):
        m = meta(case)
        name = m['case']['name']
        fn = os.path.join(s3.OUT, name + '.npz')
        pj = os.path.join(s3.OUT, name + '.json')
        if not os.path.exists(fn) or not os.path.exists(pj):
            print('%-18s no prototype run' % name)
            continue
        with open(pj) as f:
            att = [a for a in json.load(f)['attempts'] if a['status'] == 'ok']
        assert att and att[0]['nstep'] == m['nstep'], (name, att, m['nstep'])
        p = np.load(fn)
        snaps, steps, _ = snaps_of(case, m)
        n = min(len(snaps), len(p['t']))
        exp = list(range(0, m['nstep'], m['nstep'] // s3.NWRITE)) + [m['nstep']]
        dts = [abs(snaps[k][0] - p['t'][k]) for k in range(n)]
        d = [maxdiff(snaps[k][1], p['W'][k]) for k in range(n)]
        dz = float(np.max(np.abs(snaps[n - 1][2] - p['zeta'][n - 1]))) if n else None
        print('%-18s %3d/%-2d %6s %10.2e %10.2e %10.2e %10.2e %10.2e' %
              (name, n, len(p['t']), 'ok' if steps[:n] == exp[:n] else 'DIFF', max(dts), d[0], max(d), d[-1], dz))


def decomp(batch):
    """Each run against the run of the same case along x with one rank (tag _x_np1), at every line."""
    ref = {}
    for case in cases_in(batch):
        if case.endswith('_x_np1'):
            m = meta(case)
            ref[m['case']['name']] = snaps_of(case, m)[0]
    print('%-26s %6s %10s %10s' % ('case', 'lines', 'diff max', 'diff end'))
    for case in cases_in(batch):
        m = meta(case)
        r = ref.get(m['case']['name'])
        s = snaps_of(case, m)[0]
        if r is None or not s:
            print('%-26s missing' % os.path.basename(case))
            continue
        n = min(len(r), len(s))
        d = [maxdiff(s[k][1], r[k][1]) for k in range(n)]
        print('%-26s %3d/%-2d %10.2e %10.2e' % (os.path.basename(case), n, len(r), max(d), d[-1]))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'case':
        for case in sys.argv[2:]:
            d, _ = analyze_case(case)
            print(os.path.basename(case.rstrip('/')), d['attempts'][0]['status'], d['fortran']['exit'])
    elif cmd == 'table':
        table(sys.argv[2], sys.argv[3])
    elif cmd == 'agree':
        agree(sys.argv[2])
    elif cmd == 'decomp':
        decomp(sys.argv[2])
    elif cmd == 'last':
        last_states(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
