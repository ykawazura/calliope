#!/usr/bin/env python3
"""The cases of stage1_es run with the Stage 1 binary es (2126274) and with the final binary s13d (07716cc), both with
FFTW_ESTIMATE, and compared. The regression chain s3 -> s4 -> s4b -> s6 -> s13d is bitwise (Stage 14); es -> s3 had
not been compared in Fortran.

    python3 tests/MHD_COMP_ISOTH/Euler1D/check_s13d.py gen             runs/euler1d_s13d/<batch>_<case>_{old,new} and
                                                                       runs/cases_euler1d_s13d.txt
    python3 tests/MHD_COMP_ISOTH/Euler1D/check_s13d.py judge [<root>]  the pairs of <root> (runs/euler1d_s13d) and the
                                                                       negative control -> <root>/compare.json
    python3 tests/MHD_COMP_ISOTH/Euler1D/check_s13d.py pair <old> <new>  one pair, printed (a dry run on earlier output)
    python3 tests/MHD_COMP_ISOTH/Euler1D/check_s13d.py gen1 <batch/case>   a CHAINED case run in one piece:
                                                                       runs/euler1d_s13d/unchained/<batch>_<case>_{old,new}
                                                                       and runs/cases_euler1d_s13d_unchained.txt
    python3 tests/MHD_COMP_ISOTH/Euler1D/check_s13d.py values <root> [<batch/case> ...]
                                                                       the shock values of the pairs of <root> against the
                                                                       Stage 1 judge.json -> <root>/values.json

Run in the top directory of the code (ps1: ~/work/calliope_isoth) and submit the list with
    qsub -N e1s13d -l select=1:ncpus=48:mem=96gb:mpiprocs=48:ompthreads=1 \\
         -v CASEFILE=tests/MHD_COMP_ISOTH/runs/cases_euler1d_s13d.txt,PAR=48 tests/MHD_COMP_ISOTH/stages/job.pbs-ps-fes-dev
and the list of gen1 (shock/U4_N4096: 4392 s x 4 ranks in Stage 1, with MEASURE) on A_S, whose wall time limit is 24 h:
    qsub -q A_S -N e1u4 -l select=1:ncpus=8:mem=16gb:mpiprocs=8:ompthreads=1 -l walltime=02:30:00 \\
         -v CASEFILE=tests/MHD_COMP_ISOTH/runs/cases_euler1d_s13d_unchained.txt,PAR=8,BUDGET=147,I_MPI_ADJUST_ALLREDUCE=1 \\
         tests/MHD_COMP_ISOTH/stages/job.pbs-ps-fes-dev

gen copies calliope.in and nproc of each case of stage1_es/ to the two cases (md5 checked), writes the file exe and
an empty fftw_estimate in each, <root>/manifest.json and the list, by decreasing wall time of the Stage 1 run. It
leaves out the cases that ran as a restart chain (CHAINED: where a chunk ends depends on the wall time). It writes
nothing in stage1_es/ and stops if a case directory, the list or the manifest exists.

A pair passes when (fixed before the runs, tasks/MHD_COMP_ISOTH/todo.md, 2026-11-17)
- the exit status of both runs is that of the Stage 1 run (as the job script now labels it: expected(); changed after
  the runs, since the three Stage 1 runs stopped by the dt check had 'exit 0' and both of their reruns 'dt_stop'),
- every FFT plan of both runs has FFTW_ESTIMATE,
- out1d/, out2d/, out3d/, restart/ and shock_selftest.dat are the same bytes, with no file in one run only,
- the columns that es writes (conservation.dat 1-18, cfl.dat 1-16, budget.dat 1-26) have the same names in the
  headers and the same text in each row, except column 10 of budget.dat: eps_ideal of s13d also subtracts eps_eta
  (0 for C_eta = 0), which may change its rounding (1 ulp between s3 and s3e, REPORT_stage3). Its largest difference
  is reported, not judged.

Of the CHAINED cases only shock/U4_N4096 gives values of the paper (rho* error, overshoot, a point of the L1 fit). gen1
makes its pair to run in one piece (on A_dev the chunks of es and s13d would end at different steps), with the
Allreduce fixed (I_MPI_ADJUST_ALLREDUCE=1: with four ranks the sums of budget.dat depend on the order of the additions,
job 497729). The pair is judged as above (judge <root>), and values checks (fixed before the run, todo.md, 2026-11-17)
that judge_stage1.judge_shock gives the same L1(rho), rho* error, overshoot, speed and position errors for both runs. It
reports their relative differences from the Stage 1 judge.json (the run of the paper, chained, with MEASURE) and the L1
fits with L1 of N = 4096 replaced.
"""
import sys
import os
import re
import glob
import json
import fnmatch
import hashlib
import numpy as np
TEST = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(TEST, '..', 'stages'))
import paths
import analyze_stage4 as a4
import analyze_stage14 as a14
import judge_stage1 as j1

SRC = os.path.relpath(os.path.join(TEST, 'stage1_es'))
ROOT = os.path.join(paths.RUNS, 'euler1d_s13d')
LIST = os.path.join(paths.RUNS, 'cases_euler1d_s13d.txt')
UNCHAINED = os.path.join(ROOT, 'unchained')
LIST1 = os.path.join(paths.RUNS, 'cases_euler1d_s13d_unchained.txt')
JUDGE1 = os.path.join(SRC, 'judge.json')
# shock values of judge_stage1.judge_shock that must be the same in the two runs
VALUES = ('L1', 'star', 'over', 'speed', 'pos')
EXE = (('old', 'bin/calliope_isoth_es'), ('new', 'bin/calliope_isoth_s13d'))
NCASES = 102
CHAINED = ('dtconv/U4_N2048_cflsvv1', 'shock/U4_N4096', 'smoke/collide_U2_N512_est_chain')
INPUTS = ('calliope.in', 'nproc')
# output of the Stage 1 runs (job.pbs-ps-fes-dev, judge_stage1.py) and fftw_estimate, which gen writes anew
KNOWN = INPUTS + ('fftw_estimate', 'conservation.dat', 'budget.dat', 'cfl.dat', 'calliope.out.std', 'calliope_run.in',
                  'exit_status', 'out1d', 'out2d', 'out3d', 'restart', 'shock_selftest.dat', 'summary.json',
                  'chunk_*', '*.out.nc', '*.modes.out')
SOLUTION = ('out1d', 'out2d', 'out3d', 'restart')
TABLES = (('conservation.dat', 18), ('cfl.dat', 16), ('budget.dat', 26))
REPORTED = ('budget.dat', 10)
NEG = ('scan_U4_N1024_C1_old', 'scan_U4_N1024_C4_new')


def md5(path):
    with open(path, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()


def wall(case):
    """Wall time in s of the run (last 'wall X s' of calliope.out.std, from /usr/bin/time in the job)."""
    w = re.findall(r'^wall ([0-9.]+) s$', open(os.path.join(case, 'calliope.out.std'), errors='replace').read(), re.M)
    if not w:
        raise SystemExit('%s: no wall time in calliope.out.std' % case)
    return float(w[-1])


def wall_chained(case):
    """Wall time in s of a run continued from restart/: that of chunk_1, chunk_2, ... and the case directory."""
    return sum(wall(d) for d in glob.glob(os.path.join(case, 'chunk_*')) + [case])


def sources(chained=False):
    """{batch/case: directory} of stage1_es less CHAINED (checked: exactly the cases with chunk_*), or of CHAINED."""
    cases = {os.path.relpath(os.path.dirname(f), SRC): os.path.dirname(f)
             for f in glob.glob(os.path.join(SRC, '*', '*', 'calliope.in'))}
    found = sorted(k for k, d in cases.items() if glob.glob(os.path.join(d, 'chunk_*')))
    if len(cases) != NCASES or tuple(found) != CHAINED:
        raise SystemExit('%s: %d cases (expected %d), chained %s (expected %s)' % (SRC, len(cases), NCASES, found,
                                                                                  list(CHAINED)))
    for k, d in cases.items():
        unknown = [n for n in os.listdir(d) if not any(fnmatch.fnmatch(n, p) for p in KNOWN)]
        missing = [n for n in INPUTS + ('exit_status', 'calliope.out.std') if not os.path.exists(os.path.join(d, n))]
        if unknown or missing:
            raise SystemExit('%s: unknown files %s, missing %s' % (d, unknown, missing))
    return {k: d for k, d in cases.items() if (k in CHAINED) == chained}


def gen(src=None, root=ROOT, lst=LIST, walltime=wall):
    src = sources() if src is None else src
    for _, exe in EXE:
        if not os.path.isfile(exe):
            raise SystemExit('%s: no such binary' % exe)
    todo = {k: {tag: os.path.join(root, k.replace('/', '_') + '_' + tag) for tag, _ in EXE} for k in src}
    dsts = [d for v in todo.values() for d in v.values()]
    manifest = os.path.join(root, 'manifest.json')
    for f in dsts + [lst, manifest]:
        if os.path.lexists(f):
            raise SystemExit('%s exists: not written again' % f)
    if len(set(dsts)) != len(dsts):
        raise SystemExit('duplicate case names')
    man = {}
    for k, d in sorted(src.items()):
        sums = {n: md5(os.path.join(d, n)) for n in INPUTS}
        for tag, exe in EXE:
            dst = todo[k][tag]
            os.makedirs(dst)
            for n in INPUTS:
                with open(os.path.join(d, n), 'rb') as f, open(os.path.join(dst, n), 'wb') as g:
                    g.write(f.read())
                if md5(os.path.join(dst, n)) != sums[n]:
                    raise SystemExit('%s: md5 of %s differs from %s' % (dst, n, d))
            with open(os.path.join(dst, 'exe'), 'w') as f:
                f.write(exe + '\n')
            open(os.path.join(dst, 'fftw_estimate'), 'w').close()
        man[k] = dict(src=d, status=a4.status(d), wall=walltime(d), nproc=int(open(os.path.join(d, 'nproc')).read()),
                      md5=sums, **todo[k])
    order = sorted(man, key=lambda k: (-man[k]['wall'], k))
    with open(lst, 'w') as f:
        f.write('\n'.join(man[k][tag] for k in order for tag, _ in EXE) + '\n')
    with open(manifest, 'w') as f:
        json.dump(dict(cases=man, exe=dict(EXE), list=lst), f, indent=1)
    core_s = sum(2 * m['wall'] * m['nproc'] for m in man.values())
    print('%s: %d cases (%d pairs), Stage 1 wall time of the pairs %.2f core-h, longest %s %.0f s x %d ranks' % (
        lst, 2 * len(man), len(man), core_s / 3600, order[0], man[order[0]]['wall'], man[order[0]]['nproc']))


def gen1(key):
    src = sources(chained=True)
    if key not in src:
        raise SystemExit('%s: not one of %s' % (key, list(CHAINED)))
    gen({key: src[key]}, UNCHAINED, LIST1, wall_chained)


def header(path):
    """{column: name} of the header line '# 1:t 2:dt ...' of a table."""
    with open(path, errors='replace') as f:
        for line in f:
            if not line.startswith('#'):
                break
            if line.startswith('# 1:'):
                p = re.split(r'\s(\d+):', ' ' + line[1:].strip())
                return {int(c): n.strip() for c, n in zip(p[1::2], p[2::2])}
    return {}


def rows(path):
    with open(path, errors='replace') as f:
        return [line.split() for line in f if not line.startswith('#') and line.strip()]


def table(old, new, name, ncol):
    """Columns 1..ncol of <name> in the two runs: names, number of rows and text of each column (REPORTED: the
    largest |difference| and that relative to the largest |value| of the old run instead)."""
    fo, fn = os.path.join(old, name), os.path.join(new, name)
    there = [os.path.exists(fo), os.path.exists(fn)]
    if not any(there):
        return dict(present=False, ok=True)
    if not all(there):
        return dict(present=there, ok=False)
    ho, hn = header(fo), header(fn)
    ro, rn = rows(fo), rows(fn)
    rep = REPORTED[1] if name == REPORTED[0] else None
    cols = [c for c in range(1, ncol + 1) if c != rep]
    r = dict(present=True, nrow=[len(ro), len(rn)], ncol=[sorted({len(x) for x in ro}), sorted({len(x) for x in rn})],
             names_differ=['%d: %s | %s' % (c, ho.get(c), hn.get(c)) for c in cols
                           if ho.get(c) is None or ho.get(c) != hn.get(c)])
    short = [i for i, x in enumerate(ro + rn) if len(x) < ncol]
    r['short_rows'] = len(short)
    if len(ro) != len(rn) or short:
        r['ok'] = False
        return r
    first = {}
    for c in cols:
        d = [i for i, (a, b) in enumerate(zip(ro, rn)) if a[c - 1] != b[c - 1]]
        if d:
            first[c] = dict(n=len(d), row=d[0] + 1, old=ro[d[0]][c - 1], new=rn[d[0]][c - 1])
    r['differ'] = first
    if rep:
        r['reported'] = dict(column=rep, name=[ho.get(rep), hn.get(rep)])
        try:
            vo, vn = [float(x[rep - 1]) for x in ro], [float(x[rep - 1]) for x in rn]
            dmax = max((abs(a - b) for a, b in zip(vo, vn)), default=0.0)
            scale = max((abs(a) for a in vo), default=0.0)
            r['reported'].update(nrow_differ=sum(a != b for a, b in zip(vo, vn)), max_abs=dmax,
                                 max_rel=dmax / scale if scale > 0 else None, max_value=scale)
        except ValueError as e:
            r['reported']['error'] = str(e)
    r['ok'] = not r['names_differ'] and not first
    return r


def expected(src):
    """Status that job.pbs-ps-fes-dev (status_of) now gives a run whose calliope.out.std is that of the Stage 1 run
    <src>. The Stage 1 job wrote 'exit N' only; dt_stop and dt_changed came with 3b784b1 (Stage 4). es has no t_end."""
    std = open(os.path.join(src, 'calliope.out.std'), errors='replace').read()
    if 'This run is stopped' in std:
        return 'dt_stop'
    if re.search(r'dt is (de|in)creased', std):
        return 'dt_changed'
    return a4.status(src)


def pair(old, new, expect=None):
    st = [a4.status(c) for c in (old, new)]
    est = [a14.estimate_used(c) for c in (old, new)]
    same, only = {}, []
    for sub in SOLUTION:
        s, o = a14.cmp_tree(old, new, sub)
        same.update(s)
        only += o
    there = [os.path.exists(os.path.join(c, 'shock_selftest.dat')) for c in (old, new)]
    if all(there):
        same['shock_selftest.dat'] = a4.same_bytes(*(os.path.join(c, 'shock_selftest.dat') for c in (old, new)))
    elif any(there):
        only.append('shock_selftest.dat')
    tabs = {name: table(old, new, name, n) for name, n in TABLES}
    differ = sorted(k for k, v in same.items() if not v)
    ok = (st[0] == st[1] and (st[0].startswith('exit') or st[0] in ('dt_stop', 'dt_changed'))
          and (expect is None or st[0] == expect)
          and all(n > 0 and e for n, e in est) and len(same) > 0 and not differ and not only
          and all(t['ok'] for t in tabs.values()))
    return dict(old=old, new=new, ok=bool(ok), status=st, expect=expect, estimate=est, nfiles=len(same),
                differ=differ, only=only, tables=tabs)


def problem(t):
    if t['present'] is not True:
        return 'in one run only'
    if t['nrow'][0] != t['nrow'][1] or t['short_rows']:
        return 'rows %s, %d short' % (t['nrow'], t['short_rows'])
    if t['names_differ']:
        return 'names %s' % t['names_differ']
    return 'columns %s' % sorted(t['differ'])


def line(name, r):
    rep = r['tables'][REPORTED[0]].get('reported', {})
    bad = ['%s %s' % (n, problem(t)) for n, t in r['tables'].items() if not t['ok']]
    return '%-36s %-4s %-8s plans %s files %3d differ %d only %d %s col10 %s' % (
        name, 'ok' if r['ok'] else 'FAIL', r['status'][1], [n for n, _ in r['estimate']], r['nfiles'],
        len(r['differ']), len(r['only']), '; '.join(bad),
        '%d rows, max %.3g (rel %s)' % (rep['nrow_differ'], rep['max_abs'],
                                         '%.3g' % rep['max_rel'] if rep['max_rel'] is not None else '-')
        if 'max_abs' in rep else '-')


def judge(root):
    man = json.load(open(os.path.join(root, 'manifest.json')))['cases']
    out = {k: pair(m['old'], m['new'], expected(m['src'])) for k, m in sorted(man.items())}
    for k, r in out.items():
        print(line(k, r))
    neg = pair(*(os.path.join(ROOT, c) for c in NEG))
    neg['detected'] = bool(neg['differ'] or any(not t['ok'] for t in neg['tables'].values()))
    print(line('negative control', neg), 'detected' if neg['detected'] else 'NOT DETECTED')
    reps = [r['tables'][REPORTED[0]].get('reported', {}) for r in out.values()]
    reps = [x for x in reps if 'max_abs' in x]
    summ = dict(pairs=len(out), passed=sum(r['ok'] for r in out.values()),
                failed=sorted(k for k, r in out.items() if not r['ok']),
                status=sorted({r['status'][1] for r in out.values()}),
                col10=dict(pairs=len(reps), pairs_differ=sum(x['nrow_differ'] > 0 for x in reps),
                           max_abs=max((x['max_abs'] for x in reps), default=None),
                           max_rel=max((x['max_rel'] for x in reps if x['max_rel'] is not None), default=None)),
                negative_detected=neg['detected'])
    print(json.dumps(summ, indent=1))
    with open(os.path.join(root, 'compare.json'), 'w') as f:
        json.dump(dict(summary=summ, pairs=out, negative=neg), f, indent=1)


def fits(N, L1):
    """Fitted orders of L1 for N >= 512 (judge_stage1.judge_batch_shock) and N >= 1024 (plot_euler1d.euler_l1)."""
    N, L1 = np.array(N, float), np.array(L1)
    return {m: -np.polyfit(np.log(N[N >= m]), np.log(L1[N >= m]), 1)[0] for m in (512, 1024)}


def values(root, keys):
    """judge_stage1.judge_shock of the pairs of <root> (it writes summary.json in each run) against the rows of the
    Stage 1 judge.json; <root>/values.json when no case is listed."""
    man = json.load(open(os.path.join(root, 'manifest.json')))['cases']
    ref = json.load(open(JUDGE1))['shock']
    out = {}
    for k in keys or sorted(man):
        batch, case = k.split('/')
        if batch != 'shock':
            raise SystemExit('%s: values is for the batch shock' % k)
        rows, _ = j1.judge_shock([man[k]['old'], man[k]['new']])
        r0 = [x for x in ref['rows'] if x['case'] == case][0]
        num = sorted(n for n, v in r0.items() if isinstance(v, float) and isinstance(rows[0].get(n), float))
        same = {n: rows[0][n] == rows[1][n] for n in num}
        rel = {n: (rows[1][n] - r0[n]) / abs(r0[n]) if r0[n] else None for n in num}
        conv = ref['convergence'][case.split('_N')[0]]
        N = int(case.split('_N')[1])
        L1 = [rows[1]['L1'] if n == N else v for n, v in zip(conv['N'], conv['L1'])]
        out[k] = dict(old=rows[0], new=rows[1], stage1=r0, same=same, ok=all(same[n] for n in VALUES), rel=rel,
                      fit_stage1=fits(conv['N'], conv['L1']), fit_new=fits(conv['N'], L1))
        print('%s: %s (%s the same in old and new; differ: %s)' % (
            k, 'ok' if out[k]['ok'] else 'FAIL', '/'.join(VALUES), sorted(n for n in num if not same[n]) or 'none'))
        for n in num:
            print('  %-14s new %.10e  Stage 1 %.10e  rel %s' % (n, rows[1][n], r0[n],
                                                                '%.3e' % rel[n] if rel[n] is not None else '-'))
        print('  fit N >= 512: %.6f (Stage 1 %.6f), N >= 1024: %.6f (Stage 1 %.6f)' % (
            out[k]['fit_new'][512], out[k]['fit_stage1'][512], out[k]['fit_new'][1024], out[k]['fit_stage1'][1024]))
    if not keys:
        with open(os.path.join(root, 'values.json'), 'w') as f:
            json.dump(out, f, indent=1, default=float)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'gen':
        gen()
    elif cmd == 'judge':
        judge(sys.argv[2] if len(sys.argv) > 2 else ROOT)
    elif cmd == 'pair':
        r = pair(sys.argv[2], sys.argv[3])
        print(line(os.path.basename(sys.argv[3]), r))
        print(json.dumps(r, indent=1))
    elif cmd == 'gen1':
        gen1(sys.argv[2])
    elif cmd == 'values':
        values(sys.argv[2], sys.argv[3:])
    else:
        raise SystemExit(__doc__)
