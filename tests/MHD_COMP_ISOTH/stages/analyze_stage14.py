#!/usr/bin/env python3
"""Stage 14: the reruns of gen_stage14.py with the final binary against the earlier runs.

    python3 analyze_stage14.py selftest                 the comparison functions on synthetic data
    python3 analyze_stage14.py pairs <root>             P: <root>/pairs, <case>_old against <case>_new (bytes)
    python3 analyze_stage14.py mhd1d <root>             D1, D3, D4: <root>/mhd1d against runs/stage3_f
    python3 analyze_stage14.py proto <out> <dir> ...    D2: <dir>/<case>.npz (analyze_stage3.py last) against the
                                                        last snapshot of the prototype runs/stage3/<case>.npz
    python3 analyze_stage14.py records <out> <new> <old>  D5: <new>/<case>.json (analyze_stage3.py table) against
                                                        <old>/<case>.json (runs/stage3_f/all_json)
    python3 analyze_stage14.py ot <root>                O1-O4: <root>/ot/summary.json against <root>/old/ot/summary.json
                                                        (analyze_stage4.py ot of both) and the data of the runs
    python3 analyze_stage14.py seams <root>             the continued runs of <root>/mhd1d/{suite,eta_suite}
                                                        with the rows of the earlier run (seam_copy), and records
    python3 analyze_stage14.py lam <root>               lam_svv (cfl.dat) of <root>/seams against the earlier run,
                                                        step by step, before and after the first restart
    python3 analyze_stage14.py estimate <root>          <root>/estimate (gen_stage14.py estimate): the pairs as
                                                        in P, and the records of the runs
    python3 analyze_stage14.py explain <out> <root> <old>  the D5 differences: <root>/seams_rec and
                                                        <root>/estimate_rec against <old> (runs/stage3_f/all_json)
                                                        and <root>/all_json

pairs, mhd1d, ot, seams, lam and estimate write <root>/compare_<command>.json (ot also
<root>/fig_s14_ot_diff.png), proto, records and explain the file <out>. The tolerances are those of the Stage 14 plan
(tasks/todo.md), fixed before the runs.
"""
import sys
import os
import re
import glob
import json
import collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_ic
import stage3 as s3
import analyze_stage3 as a3
import analyze_stage4 as a4
import gen_stage14 as g14
from analyze_stage1 import load_runs, run_dirs

TOL_PROTO = 1e-11              # D2: final state against the prototype
TOL_OLD1D = 2e-11              # D3: all 21 lines against the earlier Fortran run
TOL_DEC1D = 1e-13              # D4: against _x_np1
REC_RTOL, REC_ATOL = 1e-9, 1e-11   # D5: |a - b| <= REC_RTOL max(|a|, |b|) + REC_ATOL
REC_SKIP = ('fortran.runs',)
# depend on where a run was continued (analyze_stage1.analyze): compared only for the same number of runs
REC_CHUNK = ('fortran.line_steps', 'fortran.cons.n', 'fortran.cons.dmass_corr_tot', 'fortran.budget.n')
TOL_OT = 1e-6                  # O4: relative difference up to t = OT_T
OT_T = np.pi
OT_COLS = dict(M=2, E_K=7, E_M=8, E=10, rho_min=11, rho_max=12, zeta_max=13)
OT_WINDOWS = (1.0, 2.0, np.pi, 4.0, 5.0, 6.0, 2*np.pi)
ESTIMATE = 64                  # FFTW_ESTIMATE
NEG = (('ot_N128_eta0_old', 'ot_N128_eta2_new'), ('I1_N256_old', 'I1_N256_eta2_new'))
N_PAIRS, N_1D, N_REC = 14, 207, 182


def dump(path, d):
    with open(path, 'w') as f:
        json.dump(d, f, indent=1, default=float)
    print('wrote', path)


# ------------------------------------------------------------ P: ESTIMATE pairs
def cmp_tree(a, b, sub):
    """{relative path: same bytes} of the files under a/sub and b/sub present in both, and those in one only."""
    def files(c):
        return {os.path.relpath(p, c) for p in glob.glob(os.path.join(c, sub, '**'), recursive=True)
                if os.path.isfile(p)}
    fa, fb = files(a), files(b)
    return {f: a4.same_bytes(os.path.join(a, f), os.path.join(b, f)) for f in sorted(fa & fb)}, sorted(fa ^ fb)


def estimate_used(case):
    """(number of plans logged by fftw_shim in the .std of all runs, whether all have FFTW_ESTIMATE)."""
    flags = []
    for d in run_dirs(case):
        fn = os.path.join(d, 'calliope.out.std')
        if os.path.exists(fn):
            flags += [int(f) for f in re.findall(r'^fftw_shim: .* -> (\d+)$', open(fn, errors='replace').read(), re.M)]
    return len(flags), bool(flags) and all(f & ESTIMATE for f in flags)


def pair(old, new):
    """analyze_stage4.regress_pair, and the bytes of out1d/, out2d/ and restart/ (files in one run only fail)."""
    r = a4.regress_pair(old, new)
    same, only = {}, []
    for sub in ('out1d', 'out2d', 'restart'):
        s, o = cmp_tree(old, new, sub)
        same.update(s)
        only += o
    r['same'].update(same)
    r['only'] = only
    r['estimate'] = [estimate_used(c) for c in (old, new)]
    r['differ'] = sorted(k for k, v in r['same'].items() if k != 'conservation.dat ncol' and not v)
    r['ok'] = bool(r['ok'] and tuple(r['status']) == ('exit 0', 'exit 0') and not r['differ'] and not only
                   and all(n > 0 and e for n, e in r['estimate']))
    return r


def pairs(root):
    d = os.path.join(root, 'pairs')
    out = {os.path.basename(o)[:-4]: pair(o, o[:-4] + '_new') for o in sorted(glob.glob(os.path.join(d, '*_old')))}
    neg = {'%s vs %s' % n: pair(*(os.path.join(d, c) for c in n)) for n in NEG}
    for name, r in sorted(out.items()):
        print('%-18s %-5s status %s plans %s files %d differ %s only %s' % (
            name, r['ok'], r['status'], r['estimate'], len(r['same']), r['differ'] or '-', r['only'] or '-'))
    for name, r in neg.items():
        print('negative %-40s differ %d of %d' % (name, len(r['differ']), len(r['same'])))
    ok = (len(out) == N_PAIRS and all(r['ok'] for r in out.values())
          and all(r['differ'] and tuple(r['status']) == ('exit 0', 'exit 0') for r in neg.values()))
    print('P:', 'PASS' if ok else 'FAIL')
    dump(os.path.join(root, 'compare_pairs.json'), dict(pairs=out, negative=neg, ok=ok))
    return ok


# ------------------------------------------------------------ D1, D3, D4: 1D against the earlier Fortran runs
def lines_1d(case):
    m = a3.meta(case)
    return m, a3.snaps_of(case, m)[0]


def complete_1d(case, m, s):
    """D1: exit 0, the lines at the steps j nstep/NWRITE and nstep, the last at t_end."""
    exp = list(range(0, m['nstep'], m['nstep'] // s3.NWRITE)) + [m['nstep']]
    steps = [int(round(x[0]/m['dt'])) for x in s]
    return (a4.status(case) == 'exit 0' and steps == exp and abs(s[-1][0] - m['t_end']) <= 1e-9*m['t_end'])


def against(s, so):
    """maxdiff of the lines of two runs at the same steps: [diff], max |t - t_old|."""
    k = min(len(s), len(so))
    return [a3.maxdiff(s[i][1], so[i][1]) for i in range(k)], max((abs(s[i][0] - so[i][0]) for i in range(k)),
                                                                   default=None)


def decomp_1d(batch):
    """{case: max over the lines of maxdiff against _x_np1 of the same problem} of a batch directory."""
    runs = {os.path.basename(c): lines_1d(c) for c in a3.cases_in(batch)}
    ref = {m['case']['name']: s for c, (m, s) in runs.items() if c.endswith('_x_np1')}
    out = {}
    for c, (m, s) in runs.items():
        r = ref.get(m['case']['name'])
        d, _ = against(s, r) if r and s else ([], None)
        out[c] = dict(max=max(d) if d else None, lines=(len(s), len(r) if r else 0))
    return out


def mhd1d(root):
    res = {}
    for batch, n in g14.MHD1D_BATCHES:
        cases = a3.cases_in(os.path.join(root, 'mhd1d', batch))
        assert len(cases) == n, (batch, len(cases))
        for case in cases:
            name = os.path.basename(case)
            m, s = lines_1d(case)
            mo, so = lines_1d(os.path.join(g14.OLD3, batch, name))
            d, dt = against(s, so)
            comp = bool(s) and complete_1d(case, m, s)
            r = dict(status=a4.status(case), runs=len(run_dirs(case)) - 1 or 1, lines=(len(s), len(so)),
                     complete=comp, old_max=max(d) if d else None, old_end=d[-1] if d else None, dt_max=dt)
            r['ok_old'] = bool(comp and len(s) == len(so) and max(d) <= TOL_OLD1D and dt == 0.0)
            res['%s/%s' % (batch, name)] = r
        worst = {c: np.inf if r['old_max'] is None else r['old_max']
                 for c, r in res.items() if c.startswith(batch + '/')}
        k = max(worst, key=worst.get)
        print('%-11s %3d cases: complete %3d, max diff against the earlier run %.2e (%s)' % (
            batch, len(cases), sum(res[c]['complete'] for c in worst), worst[k], k))
    for batch in ('decomp', 'eta_decomp'):
        new, old = decomp_1d(os.path.join(root, 'mhd1d', batch)), decomp_1d(os.path.join(g14.OLD3, batch))
        for c, d in new.items():
            r = res['%s/%s' % (batch, c)]
            r['dec'], r['dec_old'] = d['max'], old[c]['max']
            r['ok_dec'] = bool(d['max'] is not None and d['lines'][0] == d['lines'][1] and d['max'] <= TOL_DEC1D)
        print('%-11s against _x_np1: max %.2e (earlier runs %.2e)' % (
            batch, *(max(np.inf if d['max'] is None else d['max'] for d in x.values()) for x in (new, old))))
    bad = {k: r for k, r in res.items() if not (r['complete'] and r['ok_old'] and r.get('ok_dec', True))}
    for k, r in sorted(bad.items()):
        print('FAIL %-34s %s' % (k, r))
    ok = dict(D1=len(res) == N_1D and all(r['complete'] for r in res.values()),
              D3=len(res) == N_1D and all(r['ok_old'] for r in res.values()),
              D4=all(r['ok_dec'] for r in res.values() if 'ok_dec' in r))
    print(' '.join('%s: %s' % (k, 'PASS' if v else 'FAIL') for k, v in ok.items()))
    dump(os.path.join(root, 'compare_mhd1d.json'), dict(cases=res, ok=ok))
    return ok


# ------------------------------------------------------------ D2: final states against the prototype
def proto(out, dirs):
    res = {}
    for d in dirs:
        for fn in sorted(glob.glob(os.path.join(d, '*.npz'))):
            name = os.path.basename(fn)[:-4]
            a, p = np.load(fn), np.load(os.path.join(s3.OUT, name + '.npz'))
            res[name] = dict(diff=a3.maxdiff(a['W'], p['W'][-1]), dt=float(a['t']) - float(p['t'][-1]))
    ok = len(res) == N_REC and all(r['diff'] <= TOL_PROTO and r['dt'] == 0.0 for r in res.values())
    for eta in (False, True):
        sel = {k: r for k, r in res.items() if k.endswith('_eta2') == eta}
        if sel:
            k = max(sel, key=lambda k: sel[k]['diff'])
            print('%-9s %3d cases: max diff %.2e (%s), max |dt| %.1e' % (
                'eta_suite' if eta else 'suite', len(sel), sel[k]['diff'], k, max(abs(r['dt']) for r in sel.values())))
    print('D2:', 'PASS' if ok else 'FAIL')
    dump(out, dict(cases=res, ok=ok))
    return ok


# ------------------------------------------------------------ D5: records
def close(a, b):
    if a == b or (a != a and b != b):
        return True
    return abs(a - b) <= REC_RTOL*max(abs(a), abs(b)) + REC_ATOL


def cmp_values(a, b, path='', skip=(), out=None):
    """[(path, a, b)] of the items of two JSON values that differ: floats by close(), bool, int, str and None
    exactly; the paths in skip are left out."""
    out = [] if out is None else out
    if path in skip:
        return out
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            p = '%s.%s' % (path, k) if path else k
            if p in skip:
                continue
            if k in a and k in b:
                cmp_values(a[k], b[k], p, skip, out)
            else:
                out.append((p, a.get(k, '<missing>'), b.get(k, '<missing>')))
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append((path + '[len]', len(a), len(b)))
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                cmp_values(x, y, '%s[%d]' % (path, i), skip, out)
    elif any(isinstance(v, (bool, str)) or v is None for v in (a, b)) or (isinstance(a, int) and isinstance(b, int)):
        if a != b or type(a) is not type(b):
            out.append((path, a, b))
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if not close(float(a), float(b)):
            out.append((path, a, b))
    else:
        out.append((path, a, b))
    return out


def records(out, new, old):
    res = {}
    names = sorted({os.path.basename(f)[:-5] for d in (new, old) for f in glob.glob(os.path.join(d, '*.json'))})
    for name in names:
        fa, fb = os.path.join(new, name + '.json'), os.path.join(old, name + '.json')
        if not (os.path.exists(fa) and os.path.exists(fb)):
            res[name] = dict(missing=[f for f in (fa, fb) if not os.path.exists(f)], ok=False)
            continue
        a, b = (json.load(open(f)) for f in (fa, fb))
        ra, rb = a['fortran']['runs'], b['fortran']['runs']
        d = cmp_values(a, b, skip=REC_SKIP + (REC_CHUNK if ra != rb else ()))
        res[name] = dict(runs=(ra, rb), n=len(d), differ=[list(x) for x in d], ok=not d)
    bad = {k: r for k, r in res.items() if not r['ok']}
    for k, r in sorted(bad.items()):
        print('%-18s %s' % (k, r.get('missing') or ['%s: %s / %s' % tuple(x) for x in r['differ'][:6]]))
    ok = len(res) == N_REC and not bad
    print('D5: %d records, %d differ: %s' % (len(res), len(bad), 'PASS' if ok else 'FAIL'))
    dump(out, dict(cases=res, ok=ok))
    return ok


# ------------------------------------------------------------ D5: the records that differ
SEAM_FILES = ('conservation.dat', 'budget.dat', 'cfl.dat')
SEAM_LINKS = ('calliope.in', 'case.json', 'exit_status', 'calliope.out.std')
LAM_COL = 12                   # lam_svv in cfl.dat


def steps_in(case, name, dt):
    return collections.Counter(int(round(t/dt)) for r in load_runs(case, name) for t in r[:, 0])


def seam_copy(new, old, dst):
    """dst: the runs of new (chunk_k and the case directory) with the rows of SEAM_FILES at a step that old has
    fewer times left out (the later rows), and links to SEAM_LINKS and out1d/ of new. Returns {file: [steps
    left out, steps that old has more often]}."""
    dt = a3.meta(new)['dt']
    res = {}
    dirs = run_dirs(new)
    os.makedirs(dst)
    for d in dirs:
        sub = os.path.join(dst, os.path.relpath(d, new))
        os.makedirs(sub, exist_ok=True)
        if os.path.isdir(os.path.join(d, 'out1d')):
            os.symlink(os.path.abspath(os.path.join(d, 'out1d')), os.path.join(sub, 'out1d'))
    for f in SEAM_LINKS:
        if os.path.exists(os.path.join(new, f)):
            os.symlink(os.path.abspath(os.path.join(new, f)), os.path.join(dst, f))
    for name in SEAM_FILES:
        cn, co = steps_in(new, name, dt), steps_in(old, name, dt)
        drop = cn - co
        res[name] = [sorted(drop.elements()), sorted((co - cn).elements())]
        for d in reversed(dirs):
            fn = os.path.join(d, name)
            if not os.path.exists(fn):
                continue
            keep = []
            for line in reversed(open(fn).readlines()):
                k = None if line.lstrip().startswith('#') or not line.strip() else int(round(float(line.split()[0])/dt))
                if k is not None and drop[k] > 0:
                    drop[k] -= 1
                else:
                    keep.append(line)
            with open(os.path.join(dst, os.path.relpath(d, new), name), 'w') as f:
                f.writelines(reversed(keep))
    return res


def seams(root):
    """seam_copy of the continued runs (new or earlier) of suite and eta_suite in <root>/seams/<batch>/<case>, and
    the records of the copies (analyze_stage3.analyze_case) in <root>/seams_rec, for explain."""
    rec = os.path.join(root, 'seams_rec')
    os.makedirs(rec)
    res = {}
    for batch in ('suite', 'eta_suite'):
        for case in a3.cases_in(os.path.join(root, 'mhd1d', batch)):
            name = os.path.basename(case)
            old = os.path.join(g14.OLD3, batch, name)
            runs = (len(run_dirs(case)), len(run_dirs(old)))
            if runs == (1, 1):
                continue
            dst = os.path.join(root, 'seams', batch, name)
            res[name] = dict(runs=runs, rows=seam_copy(case, old, dst))
            dump(os.path.join(rec, name + '.json'), a3.analyze_case(dst)[0])
            print('%-18s runs %s left out %s, more in the earlier run %s' % (
                name, runs, {f[:4]: r[0] for f, r in res[name]['rows'].items()},
                {f[:4]: r[1] for f, r in res[name]['rows'].items() if r[1]} or '-'))
    dump(os.path.join(root, 'compare_seams.json'), res)
    return res


def lam_by_step(case, dt):
    """({step: [lam_svv of each row at the step]}, [first step of each run]) of cfl.dat."""
    out, starts = collections.defaultdict(list), []
    for r in load_runs(case, 'cfl.dat'):
        starts.append(int(round(r[0, 0]/dt)))
        for x in r:
            out[int(round(x[0]/dt))].append(x[LAM_COL])
    return out, starts


def lam(root):
    """lam_svv of the seam copies (<root>/seams) against the earlier run at the steps with one row in each, before
    the first restart of either run and after it."""
    res = {}
    for dst in sorted(glob.glob(os.path.join(root, 'seams', '*', '*'))):
        name = os.path.basename(dst)
        dt = a3.meta(dst)['dt']
        (a, sa), (o, so) = lam_by_step(dst, dt), lam_by_step(os.path.join(g14.OLD3, os.path.basename(
            os.path.dirname(dst)), name), dt)
        k1 = min(sa[1:] + so[1:])
        r = {k: rel(a[k][0], o[k][0]) for k in a if len(a[k]) == 1 and len(o.get(k, ())) == 1}
        res[name] = dict(first_restart=k1, starts=[sa, so],
                         before=[sum(k < k1 for k in r), max([v for k, v in r.items() if k < k1], default=0.0)],
                         after=[sum(k >= k1 for k in r), max([v for k, v in r.items() if k >= k1], default=0.0)])
        print('%-18s first restart %6d  before: %6d steps max rel %.1e  after: %6d steps max rel %.1e' % (
            name, k1, *res[name]['before'], *res[name]['after']))
    dump(os.path.join(root, 'compare_lam.json'), res)
    return res


def estimate(root):
    """<root>/estimate (gen_stage14.py estimate): the pairs as in P, and the records of the _new runs
    (analyze_stage3.analyze_case) in <root>/estimate_rec, for explain."""
    rec = os.path.join(root, 'estimate_rec')
    os.makedirs(rec)
    res = {}
    for o in sorted(glob.glob(os.path.join(root, 'estimate', '*_old'))):
        name = os.path.basename(o)[:-4]
        r = res[name] = pair(o, o[:-4] + '_new')
        print('%-18s %-5s status %s plans %s files %d differ %s only %s' % (
            name, r['ok'], r['status'], r['estimate'], len(r['same']), r['differ'] or '-', r['only'] or '-'))
        dump(os.path.join(rec, name + '.json'), a3.analyze_case(o[:-4] + '_new')[0])
    ok = sorted(res) == sorted(g14.EST_CASES) and all(r['ok'] for r in res.values())
    print('ESTIMATE pairs:', 'PASS' if ok else 'FAIL')
    dump(os.path.join(root, 'compare_estimate.json'), dict(pairs=res, ok=ok))
    return ok


def item(path):
    return re.sub(r'\[\d+\]', '', path)


def rel(a, b):
    """|a - b|/max(|a|, |b|) of two numbers (None for other values)."""
    if not all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in (a, b)):
        return None
    return abs(a - b)/max(abs(a), abs(b)) if a != b else 0.0


def by_item(diffs):
    """{item: [number of differing values, max rel]} of [(path, a, b)]."""
    out = {}
    for p, a, b in diffs:
        n, m = out.get(item(p), (0, 0.0))
        r = rel(a, b)
        out[item(p)] = [n + 1, m if r is None else max(m, r)]
    return out


def explain(out, root, old):
    """The D5 differences of the continued runs after seam_copy (<root>/seams_rec against <old>, line_steps left
    out), and the records of the ESTIMATE runs (<root>/estimate_rec, E) against the MEASURE runs of s13d
    (<root>/all_json, N) and s3g (<old>, O)."""
    def rec(d, n):
        with open(os.path.join(d, n + '.json')) as f:
            return json.load(f)
    d5 = json.load(open(os.path.join(root, 'compare_records.json')))['cases']
    res = dict(seams={}, estimate={})
    for fn in sorted(glob.glob(os.path.join(root, 'seams_rec', '*.json'))):
        n = os.path.basename(fn)[:-5]
        d = cmp_values(rec(os.path.dirname(fn), n), rec(old, n), skip=REC_SKIP + ('fortran.line_steps',))
        res['seams'][n] = dict(before=by_item(d5[n]['differ']), after=by_item(d))
        print('seams %-18s %3d -> %2d values: %s' % (n, d5[n]['n'], len(d), ', '.join(
            '%s %.0e' % (k.split('.')[-1], v[1]) for k, v in sorted(res['seams'][n]['after'].items())) or '-'))
    tot = {k: collections.Counter() for k in ('NO', 'EN', 'EO')}
    for fn in sorted(glob.glob(os.path.join(root, 'estimate_rec', '*.json'))):
        n = os.path.basename(fn)[:-5]
        e, rn, ro = rec(os.path.dirname(fn), n), rec(os.path.join(root, 'all_json'), n), rec(old, n)
        r = res['estimate'][n] = {k: by_item(cmp_values(a, b, skip=REC_SKIP))
                                  for k, (a, b) in dict(NO=(rn, ro), EN=(e, rn), EO=(e, ro)).items()}
        for k, v in r.items():
            tot[k].update(v.keys())
        print('estimate %-18s %s' % (n, '  '.join('%s %d' % (k, sum(x[0] for x in v.values())) for k, v in r.items())))
    items = sorted(set().union(*tot.values()))
    print('%-58s %s' % ('item (cases with a difference, max rel)', '   '.join('%-14s' % k for k in tot)))
    for it in items:
        print('%-58s %s' % (it, '   '.join('%3d %.1e' % (tot[k][it], max(
            [r[k][it][1] for r in res['estimate'].values() if it in r[k]], default=0.0)) for k in tot)))
    res['items'] = {it: {k: tot[k][it] for k in tot} for it in items}
    dump(out, res)
    return res


# ------------------------------------------------------------ O1-O4: OT
def by_step(rows, dt):
    """{step from t = 0: row} (a continued run repeats its first row; the later one is kept)."""
    return {int(round(r[0]/dt)): r for r in rows}


def cons_diff(a, b, dt, tmax):
    """{column: max|a - b|/max|b|} of conservation.dat rows a, b at the common steps with t <= tmax."""
    ra, rb = by_step(a, dt), by_step(b, dt)
    n = sorted(k for k in set(ra) & set(rb) if rb[k][0] <= tmax + 1e-9)
    A, B = np.array([ra[k] for k in n]), np.array([rb[k] for k in n])
    return {q: float(np.max(np.abs(A[:, j] - B[:, j]))/max(float(np.max(np.abs(B[:, j]))), 1e-300))
            for q, j in OT_COLS.items()}


PLANE_VEC = dict(rho=('rho',), ux=('ux', 'uy', 'uz'), bx=('bx', 'by', 'bz'))


def plane_rel(fa, fb, dl):
    """max|fa - fb| of rho, u, B relative to max|rho|, max|u_i|, max|B_i| of fb (a component that vanishes, uy and
    by of OT, against the scale of its vector) and the 'max' of them; zeta, D, eta as plane_diff."""
    d = {}
    for qs in PLANE_VEC.values():
        sc = max(float(np.abs(fb[q]).max()) for q in qs)
        for q in qs:
            e = float(np.abs(fa[q] - fb[q]).max())
            d[q] = e/sc if sc > 0 else (0.0 if e == 0 else np.inf)
    d['max'] = max(d[q] for q in a4.PLANE_SCALE1)
    pd = a4.plane_diff(fa, fb, dl)
    d.update({q: pd[q] for q in ('zeta', 'D', 'eta')})
    return d


def conv_items(c):
    """{item: passed} of the checks of ot_convergence (their conjunction is c['ok'])."""
    it = {'state %s' % q: v <= a4.OT_CONV_STATE for q, v in c.get('state', {}).items()}
    for q, tol in (('eps_peak', a4.OT_CONV_EPS), ('t_peak', a4.OT_CONV_TPEAK), ('share_visc', a4.OT_CONV_SHARE)):
        if q in c:
            it[q] = c[q] is not None and c[q] <= tol
    lo, hi = a4.OT_KEFF_RATIO
    it.update({'keff_ratio %s' % k: lo <= v <= hi for k, v in c.get('keff_ratio', {}).items()})
    lo, hi = a4.OT_ZFRAC_RATIO
    it['a90_last'] = bool(c.get('a90_last'))
    it.update({'a90_last %s' % k: lo <= v <= hi for k, v in c.get('a90_last', {}).items()})
    return it


def share_at(r, t):
    return float(np.interp(t, r['cum']['t'], r['cum']['visc'])/np.interp(t, r['cum']['t'], r['cum']['eps']))


def ot_scalars(r, tj):
    """The scalars of S4.3 up to tj: the peak of eps and its t, the share of the shock viscosity at tj and the
    state at t = 1, 2, tj."""
    s = dict(eps_peak=r['eps_peak']['eps'], t_peak=r['eps_peak']['t'], share_visc=share_at(r, tj))
    for t in (1.0, 2.0, tj):
        s.update({'%s(%.4f)' % (q, t): a4.series_at(r, q, t) for q in ('E_K', 'E_M', 'rho_min', 'rho_max')})
    return s


def ot_pair(new, old, m, rn, ro):
    """new against old up to tj = min(OT_T, t_end) (N = 1024: t_end = nstep//2 dt, below pi for an odd nstep);
    the conservation.dat columns also up to the later OT_WINDOWS."""
    N, dt = m['N'], m['dt']
    tj = min(OT_T, m['nrun']*dt)
    ca, cb = (np.concatenate(load_runs(c, 'conservation.dat')) for c in (new, old))
    cons = {'%.4f' % t: cons_diff(ca, cb, dt, t) for t in sorted({min(t, m['nrun']*dt) for t in OT_WINDOWS})}
    g = gen_ic.Grid(N, 2, N)
    pa, pb = (by_step(a4.plane_records(c), dt) for c in (new, old))
    planes = []
    for k in sorted(set(pa) & set(pb)):
        (ta, _, fa, ia), (tb, _, fb, ib) = pa[k], pb[k]
        planes.append(dict(t=tb, dt=ta - tb, **plane_rel(a4.plane_at(g, fa, ia), a4.plane_at(g, fb, ib), 2*np.pi/N)))
    sa, sb = ot_scalars(rn, tj), ot_scalars(ro, tj)
    scal = {k: a4.rel(sa[k], sb[k]) for k in sa}
    zeta = {k: dict(frac=(z['frac'], ro['zeta'][k]['frac']), a90=(z.get('a90'), ro['zeta'][k].get('a90')),
                    same=z['frac'] == ro['zeta'][k]['frac'] and z.get('a90') == ro['zeta'][k].get('a90'))
            for k, z in rn['zeta'].items() if float(k) <= tj + 1e-3 and k in ro['zeta']}
    early = [p for p in planes if p['t'] <= tj + 1e-9]
    r = dict(t_judge=tj, cons=cons, planes=planes, scalars=dict(new=sa, old=sb, rel=scal), zeta=zeta,
             max_cons=max(cons['%.4f' % tj].values()), max_planes=max(p['max'] for p in early),
             max_scalars=max(scal.values()), planes_n=(len(early), len(planes)))
    r['ok'] = bool(r['max_cons'] <= TOL_OT and r['max_planes'] <= TOL_OT and r['max_scalars'] <= TOL_OT
                   and all(p['dt'] == 0.0 for p in planes))
    return r


def ot_figure(root, res):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    Ns = sorted({r['N'] for r in res.values()})
    shade = {N: plt.cm.Blues(0.4 + 0.6*i/max(1, len(Ns) - 1)) for i, N in enumerate(Ns)}
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for name, r in sorted(res.items(), key=lambda x: (x[1]['C_eta'], x[1]['N'])):
        t = [p['t'] for p in r['planes']]
        ax.semilogy(t, [max(p['max'], 1e-17) for p in r['planes']], color=shade[r['N']], lw=1.5,
                    ls='-' if r['C_eta'] == 0 else '--', marker='o', ms=3,
                    label='N = %d, C_eta = %g' % (r['N'], r['C_eta']))
    ax.axvline(OT_T, color='0.5', lw=0.8, ls=':')
    ax.axhline(TOL_OT, color='0.5', lw=0.8, ls=':')
    ax.set_xlabel('t')
    ax.set_ylabel('max |new - old| / max |old| (rho, u, B)')
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(root, 'fig_s14_ot_diff.png'), dpi=120)
    plt.close(fig)


def ot(root):
    new, old = (json.load(open(os.path.join(root, d, 'summary.json'))) for d in ('ot', os.path.join('old', 'ot')))
    runs = list(g14.OT_RUNS)
    o1 = {k: c for k, c in new['checks'].items() if k.startswith('run ')}
    o2 = {k: c for k, c in new['checks'].items() if k.startswith('decomp ')}
    ok = dict(O1=sorted(o1) == sorted('run ' + c for c in runs) and all(c['ok'] for c in o1.values()),
              O2=len(o2) == 3 and all(c['ok'] for c in o2.values()))
    for k in sorted(o1) + sorted(o2):
        c, co = new['checks'][k], old['checks'].get(k, {})
        print('%-28s %-5s %s   (earlier: %s %s)' % (k, c['ok'], ' '.join('%s %.2e' % (q, c[q]) for q in ('R', 'max')
                                                                       if q in c), co.get('ok'),
                                                  ' '.join('%s %.2e' % (q, co[q]) for q in ('R', 'max') if q in co)))
    conv = {}
    for e in sorted(set(new['conv']) | set(old['conv'])):
        cn, co = new['conv'].get(e, {}), old['conv'].get(e, {})
        itn, ito = conv_items(cn), conv_items(co)
        conv[e] = dict(ok=(cn.get('ok'), co.get('ok')), fail_new=sorted(k for k, v in itn.items() if not v),
                       fail_old=sorted(k for k, v in ito.items() if not v),
                       consistent=(all(itn.values()) == cn.get('ok'), all(ito.values()) == co.get('ok')),
                       values={k: (cn.get(k), co.get(k)) for k in ('state', 'eps_peak', 't_peak', 'share_visc',
                                                                    'keff_ratio', 'a90_last')},
                       keff={p: (v['keff'], co.get('pairs', {}).get(p, {}).get('keff'))
                             for p, v in cn.get('pairs', {}).items()})
        print('conv %-5s ok %s fail new %s earlier %s' % (e, conv[e]['ok'], conv[e]['fail_new'], conv[e]['fail_old']))
    ok['O3'] = all(c['ok'][0] == c['ok'][1] and c['fail_new'] == c['fail_old'] and all(c['consistent'])
                   for c in conv.values()) and len(conv) == 2
    res = {}
    for name in runs:
        m = a4.case_info(os.path.join(root, 'ot', name))
        r = ot_pair(os.path.join(root, 'ot', name), os.path.join(root, 'old', 'ot', name), m,
                    new['runs'][name], old['runs'][name])
        r.update(N=m['N'], C_eta=m['C_eta'])
        res[name] = r
        late = [p for p in r['planes'] if p['t'] > r['t_judge'] + 1e-9]
        print('%-20s %-5s t <= pi: cons %.1e planes %.1e (%d) scalars %.1e zeta same %s; later planes %s' % (
            name, r['ok'], r['max_cons'], r['max_planes'], r['planes_n'][0], r['max_scalars'],
            all(z['same'] for z in r['zeta'].values()), ' '.join('%.2f:%.0e' % (p['t'], p['max']) for p in late[::4])))
    ok['O4'] = len(res) == len(runs) and all(r['ok'] for r in res.values())
    print(' '.join('%s: %s' % (k, 'PASS' if v else 'FAIL') for k, v in ok.items()))
    ot_figure(root, res)
    dump(os.path.join(root, 'compare_ot.json'), dict(checks=dict(o1=o1, o2=o2), conv=conv, runs=res, ok=ok))
    return ok


# ------------------------------------------------------------ selftest
def selftest():
    import tempfile
    ok = []

    def check(name, cond):
        ok.append(bool(cond))
        print('%-52s %s' % (name, 'ok' if cond else 'FAILED'))

    rec = dict(case=dict(name='X', N=256), fortran=dict(runs=1, line_steps=[0, 10], cons=dict(n=21, dM=1e-6)),
               analysis=dict(l1=[0.1, 0.2], ok=True, n=3, s=None, nan=float('nan')))
    dup = json.loads(json.dumps(rec))
    check('records: identical', cmp_values(rec, dup) == [])
    dup['analysis']['l1'][1] *= 1 + 1e-12
    check('records: relative 1e-12 passes', cmp_values(rec, dup) == [])
    dup['analysis']['l1'][1] *= 1 + 1e-8
    check('records: relative 1e-8 fails', [x[0] for x in cmp_values(rec, dup)] == ['analysis.l1[1]'])
    dup = json.loads(json.dumps(rec))
    dup['fortran']['cons']['dM'] += 5e-12
    check('records: absolute 5e-12 on 1e-6 passes', cmp_values(rec, dup) == [])
    dup['analysis']['n'] = 4
    dup['case']['name'] = 'Y'
    dup['analysis']['s'] = 0.0
    check('records: int, str, None differ', len(cmp_values(rec, dup)) == 3)
    dup = json.loads(json.dumps(rec))
    dup['fortran']['runs'], dup['fortran']['line_steps'] = 2, [0, 3]
    del dup['analysis']['ok']
    check('records: missing key, list differ', len(cmp_values(rec, dup, skip=REC_SKIP)) == 2)
    check('records: skipped chunk keys', len(cmp_values(rec, dup, skip=REC_SKIP + REC_CHUNK)) == 1)
    dup['analysis']['l1'] = [0.1]
    check('records: list length', ('analysis.l1[len]', 2, 1) in cmp_values(rec, dup))

    with tempfile.TemporaryDirectory() as tmp:
        rng = np.random.default_rng(1)
        cases = {}
        for c in ('a_old', 'a_new', 'b_new'):
            d = os.path.join(tmp, c)
            for sub in ('out1d', 'out2d', 'restart'):
                os.makedirs(os.path.join(d, sub))
            with open(os.path.join(d, 'calliope.in'), 'w') as f:
                f.write('&time_parameters\n nstep = 10\n/\n')
            with open(os.path.join(d, 'exit_status'), 'w') as f:
                f.write('exit 0\n')
            with open(os.path.join(d, 'calliope.out.std'), 'w') as f:
                f.write('fftw_shim: r2c  n = 128 howmany = 1 flags 0 -> 64\nfftw_shim: c2r  n = 128 howmany = 1 '
                        'flags 1 -> 65\n')
            cases[c] = d
        cons = rng.random((5, 20))
        for c, d in cases.items():
            np.savetxt(os.path.join(d, 'conservation.dat'), cons)
            for f in ('cfl.dat', 'budget.dat'):
                np.savetxt(os.path.join(d, f), cons[:, :8])
            (cons[:2] if c != 'b_new' else cons[:2] + 1e-15).tofile(os.path.join(d, 'out2d', 'plane_y0.dat'))
            cons[0].tofile(os.path.join(d, 'restart', 'lnrho.dat'))
            np.savetxt(os.path.join(d, 'out1d', 'line_x.dat'), cons[:3])
        r = pair(cases['a_old'], cases['a_new'])
        check('pairs: identical runs pass', r['ok'] and r['estimate'] == [(2, True), (2, True)])
        r = pair(cases['a_old'], cases['b_new'])
        check('pairs: a plane differing by 1e-15 fails', not r['ok'] and r['differ'] == ['out2d/plane_y0.dat'])
        with open(os.path.join(cases['b_new'], 'restart', 'extra.dat'), 'w') as f:
            f.write('1\n')
        check('pairs: a file in one run only fails',
              pair(cases['a_old'], cases['b_new'])['only'] == ['restart/extra.dat'])
        with open(os.path.join(cases['a_new'], 'calliope.out.std'), 'w') as f:
            f.write('fftw_shim: r2c  n = 128 howmany = 1 flags 0 -> 0\n')
        check('pairs: a plan without FFTW_ESTIMATE fails', not pair(cases['a_old'], cases['a_new'])['ok'])
        os.remove(os.path.join(cases['a_new'], 'calliope.out.std'))
        check('pairs: no shim log fails', not pair(cases['a_old'], cases['a_new'])['ok'])

    dt = 0.01
    a = np.zeros((6, 20))
    a[:, 0] = [0.0, 0.01, 0.02, 0.02, 0.03, 0.04]       # a continued run repeats the row at t = 0.02
    a[:, 2:] = 1.0 + np.arange(6)[:, None]*0.1
    b = np.delete(a, 2, axis=0)
    check('OT: conservation rows matched by step', max(cons_diff(a, b, dt, 1.0).values()) == 0.0)
    b2 = b.copy()
    b2[3, 7] *= 1 + 2e-6
    d = cons_diff(a, b2, dt, 1.0)
    check('OT: conservation difference 2e-6 found', d['E_K'] > TOL_OT and d['M'] == 0.0)
    check('OT: rows after tmax left out', cons_diff(a, b2, dt, 0.02)['E_K'] == 0.0)
    g = {q: rng.random((8, 8)) + 1.0 for q in gen_ic.PLANE_Q}
    g['uy'] = np.zeros((8, 8))
    h = {q: v.copy() for q, v in g.items()}
    check('OT: identical planes', plane_rel(g, h, 0.1)['max'] == 0.0)
    h['bx'] = h['bx']*(1 + 3e-6)
    h['uy'] = h['uy'] + 1e-20
    p = plane_rel(h, g, 0.1)
    check('OT: plane difference relative to max|old|', TOL_OT < p['bx'] <= 3e-6 and p['max'] == p['bx'])
    check('OT: vanishing uy against max|u|', 0 < p['uy'] < 1e-19)
    c = dict(state=dict(rho_min=0.03, E_K=0.01), eps_peak=0.01, t_peak=0.001, share_visc=0.05,
             keff_ratio={'1.0000': 1.5, '3.0000': 2.0}, a90_last={'1.0000': 0.5})
    it = conv_items(c)
    check('OT: convergence items', sorted(k for k, v in it.items() if not v) == ['keff_ratio 1.0000',
                                                                                  'state rho_min'])
    with tempfile.TemporaryDirectory() as tmp:
        new, old = os.path.join(tmp, 'new'), os.path.join(tmp, 'old')
        rows = lambda steps: np.array([[k*0.5, 0.5] + [k + 1.0]*4 for k in steps])
        # earlier run: a seam after step 3 without the budget row of step 3; new: a seam after step 4 whose second
        # run repeats the conservation row of step 5
        for d, cons, bud in ((old, ([0, 1, 2, 3], [4, 5, 6]), ([0, 1, 2], [4, 5])),
                             (new, ([0, 1, 2, 3, 4, 5], [5, 6]), ([0, 1, 2, 3, 4], [5]))):
            os.makedirs(os.path.join(d, 'chunk_1', 'out1d'))
            os.makedirs(os.path.join(d, 'out1d'))
            with open(os.path.join(d, 'case.json'), 'w') as f:
                json.dump(dict(dt=0.5, nstep=6), f)
            for k, sub in enumerate(('chunk_1', '')):
                np.savetxt(os.path.join(d, sub, 'conservation.dat'), rows(cons[k]), header='1:t')
                np.savetxt(os.path.join(d, sub, 'budget.dat'), rows(bud[k]))
        r = seam_copy(new, old, os.path.join(tmp, 'copy'))
        check('seams: rows left out', r['budget.dat'] == [[3], []] and r['conservation.dat'] == [[5], []])
        cp = os.path.join(tmp, 'copy')
        check('seams: copy has the rows of the earlier run', all(
            steps_in(cp, f, 0.5) == steps_in(old, f, 0.5) for f in ('budget.dat', 'conservation.dat')))
        check('seams: the later repeated row left out', len(load_runs(cp, 'conservation.dat')[0]) == 6
              and os.path.islink(os.path.join(cp, 'chunk_1', 'out1d'))
              and os.path.islink(os.path.join(cp, 'case.json')))
        r = seam_copy(old, new, os.path.join(tmp, 'copy2'))
        check('seams: rows missing in new reported', r['budget.dat'] == [[], [3]])

    check('records: items without indices', by_item([('a.b[1].c', 1.0, 2.0), ('a.b[3].c', 1.0, 1.5),
                                                     ('x', 'p', 'q')]) == {'a.b.c': [2, 0.5], 'x': [1, 0.0]})
    one = np.ones((6, 4))
    check('maxdiff of analyze_stage3: relative to max(1, max|b|)',
          a3.maxdiff(2*one, 4*one) == 0.5 and np.isclose(a3.maxdiff(0.6*one, 0.5*one), 0.1))
    print('selftest:', 'PASS' if all(ok) else 'FAIL')
    return all(ok)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'selftest':
        sys.exit(0 if selftest() else 1)
    elif cmd in ('pairs', 'mhd1d', 'ot', 'seams', 'lam', 'estimate'):
        globals()[cmd](sys.argv[2])
    elif cmd == 'proto':
        proto(sys.argv[2], sys.argv[3:])
    elif cmd == 'records':
        records(*sys.argv[2:5])
    elif cmd == 'explain':
        explain(*sys.argv[2:5])
    else:
        sys.exit(__doc__)
