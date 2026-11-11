#!/usr/bin/env python3
"""Stage 14: the Stage 3 1D MHD Riemann problems and the Stage 4 OT vortex again with the final binary EXE_NEW.

The cases are copies of the input files of the earlier runs, not generated again: the defaults of isoth_case
changed after Stage 4, and gen_stage4.ot_nstep_from takes nstep from the runs on the next lower N.

    python3 tests/MHD_COMP_ISOTH/stages/gen_stage14.py pairs   runs/stage14/pairs: <case>_old (EXE_OLD_1D or EXE_OLD_OT)
                                                                and <case>_new (EXE_NEW) of the 1D smoke and eta_smoke
                                                                cases and of ot_N128_eta0, eta2, all with FFTW_ESTIMATE
    python3 tests/MHD_COMP_ISOTH/stages/gen_stage14.py mhd1d   runs/stage14/mhd1d/<batch> of MHD1D_BATCHES (EXE_NEW;
                                                                fftw_estimate as in the earlier case)
    python3 tests/MHD_COMP_ISOTH/stages/gen_stage14.py ot      runs/stage14/ot of OT_JOBS (EXE_NEW, as above), and
                                                                runs/stage14/old/ot: links to the earlier cases, for
                                                                analyze_stage4.py ot runs/stage14/old
    python3 tests/MHD_COMP_ISOTH/stages/gen_stage14.py estimate  runs/stage14/estimate: <case>_old (EXE_OLD_1D) and
                                                                <case>_new (EXE_NEW) of EST_CASES, with FFTW_ESTIMATE

Run in the top directory of the code (ps1: ~/work/calliope_isoth). Each command checks all its cases before it
writes anything and stops if a case directory or a list exists or an earlier case holds a file it does not know.
It copies the input files (INPUTS, ic/), compares their md5 with the originals, writes the file "exe", the lists
runs/cases_s14_<name>.txt (cases by decreasing cost, as gen_stage3.py) and runs/stage14/manifest_<command>.json.
Nothing is written in runs/stage3_f and runs/stage4.
"""
import sys
import os
import glob
import json
import shutil
import fnmatch
import hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_stage3 import cost

OLD3 = os.path.join('runs', 'stage3_f')
OLD4 = os.path.join('runs', 'stage4', 'ot')
ROOT = os.path.join('runs', 'stage14')
EXE_OLD_1D = 'bin/calliope_isoth_s3g'
EXE_OLD_OT = 'bin/calliope_isoth_s4'
EXE_NEW = 'bin/calliope_isoth_s13d'
INPUTS = ('calliope.in', 'nproc', 'case.json', 'profile1d.dat', 't_end', 'fftw_estimate')
REQUIRED = ('calliope.in', 'nproc', 'case.json')
# output of a run (job.pbs-ps-fes-dev, analyze_stage3.analyze_case, analyze_stage1) and the file exe, which is
# written anew
OUTPUTS = ('conservation.dat', 'budget.dat', 'cfl.dat', 'calliope.out.std', 'calliope_run.in', 'exit_status',
           'out1d', 'out2d', 'out3d', 'restart', 'stage3.json', 'summary.json', 'superseded_by', 'shock_selftest.dat',
           'exe', 'chunk_*', '*.out.nc', '*.modes.out')
PAIRS_1D = (('smoke', 8), ('eta_smoke', 4))
PAIRS_OT = ('ot_N128_eta0', 'ot_N128_eta2')
MHD1D_BATCHES = (('suite', 131), ('eta_suite', 51), ('decomp', 15), ('eta_decomp', 10))
OT_RUNS = ('ot_N128_eta0', 'ot_N128_eta2', 'ot_N256_eta0', 'ot_N256_eta2', 'ot_N512_eta0', 'ot_N512_eta2',
           'ot_N1024_eta0_p256')
OT_DEC = tuple('dec_N256_np%d%s' % (P, s) for P in (4, 16, 64, 256) for s in ('_n10', ''))
# one job each (A_dev: 30 min): N = 128, 256 and the decomposition test; N = 512; N = 1024 (continued)
OT_JOBS = (('ot1', OT_RUNS[:4] + OT_DEC), ('ot2', OT_RUNS[4:6]), ('ot3', OT_RUNS[6:]))
# the cases of suite and eta_suite run once (no restart) whose records D5 fails, less TDS_w1pre_N2048 and
# TDS_w4pre_N2048 (over 20 min each: a continued run would stop at a different step in the two runs)
EST_CASES = ('C1_N1024_eta2', 'C1_N512', 'C2_N1024', 'C2_N256', 'P1b_N256_eta2', 'P1b_N512', 'P1b_N512_eta2',
             'P1c_N1024', 'P1c_N256', 'RDh_w4_N1024', 'RDh_w4_N1024_eta2', 'RDh_w4_N256', 'RDh_w4_N256_eta2',
             'RDh_w4_N512', 'RDh_w8_N256', 'RDh_w8_N512', 'RDp_w4_N1024', 'RDp_w4_N2048', 'RDp_w4_N256',
             'RDp_w4_N512', 'RDp_w8_N1024', 'RDp_w8_N256', 'RDp_w8_N512', 'S1_N2048_eta2', 'TDA_w0pre_N1024',
             'TDA_w1pre_N512', 'TDA_w4_N1024', 'TDA_w4_N256', 'TDA_w4_N512', 'TDA_w4pre_N256', 'TDS_w1_N1024',
             'TDS_w1_N512_eta2', 'TDS_w1pre_N1024', 'TDS_w1pre_N256', 'TDS_w1pre_N512', 'TDS_w4pre_N1024',
             'TDS_w4pre_N256', 'TDS_w4pre_N512')


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def cases_in(batch):
    return sorted(os.path.dirname(f) for f in glob.glob(os.path.join(batch, '*', 'case.json')))


def check_src(src):
    """Input files of the earlier case src (INPUTS present, and 'ic' for a directory ic/)."""
    names = os.listdir(src)
    unknown = [n for n in names if n not in INPUTS + ('ic',) and not any(fnmatch.fnmatch(n, p) for p in OUTPUTS)]
    missing = [n for n in REQUIRED if n not in names]
    if unknown or missing:
        raise SystemExit('%s: unknown files %s, missing %s' % (src, unknown, missing))
    return [n for n in INPUTS + ('ic',) if n in names]


def copy_case(src, dst, exe, estimate):
    """dst with the input files of src (md5 checked) and the file exe; estimate True writes fftw_estimate, None
    copies it if src has it. Returns the manifest entry."""
    names = check_src(src)
    os.makedirs(dst)
    files = []
    for n in names:
        if n == 'fftw_estimate' and estimate is not None:
            continue
        if n == 'ic':
            shutil.copytree(os.path.join(src, n), os.path.join(dst, n))
            files += sorted(os.path.join(n, f) for f in os.listdir(os.path.join(src, n)))
        else:
            shutil.copy2(os.path.join(src, n), os.path.join(dst, n))
            files.append(n)
    sums = {}
    for f in files:
        a, b = md5(os.path.join(src, f)), md5(os.path.join(dst, f))
        if a != b:
            raise SystemExit('%s: md5 of %s differs from %s' % (dst, f, src))
        sums[f] = a
    with open(os.path.join(dst, 'exe'), 'w') as f:
        f.write(exe + '\n')
    if estimate:
        open(os.path.join(dst, 'fftw_estimate'), 'w').close()
    return dict(src=src, exe=exe, fftw_estimate=os.path.exists(os.path.join(dst, 'fftw_estimate')), md5=sums)


def plan_pairs():
    todo = []
    for batch, n in PAIRS_1D:
        srcs = cases_in(os.path.join(OLD3, batch))
        if len(srcs) != n:
            raise SystemExit('%s: %d cases, expected %d' % (batch, len(srcs), n))
        todo += [(s, EXE_OLD_1D) for s in srcs]
    todo += [(os.path.join(OLD4, c), EXE_OLD_OT) for c in PAIRS_OT]
    out = [(src, os.path.join(ROOT, 'pairs', os.path.basename(src) + '_' + tag), exe, True)
           for src, old in todo for tag, exe in (('old', old), ('new', EXE_NEW))]
    return {'pairs': [o[1] for o in out]}, out


def plan_mhd1d():
    out = []
    for batch, n in MHD1D_BATCHES:
        srcs = cases_in(os.path.join(OLD3, batch))
        if len(srcs) != n:
            raise SystemExit('%s: %d cases, expected %d' % (batch, len(srcs), n))
        out += [(s, os.path.join(ROOT, 'mhd1d', batch, os.path.basename(s)), EXE_NEW, None) for s in srcs]
    return {'mhd1d': [o[1] for o in out]}, out


def plan_ot():
    lists, out = {}, []
    for job, names in OT_JOBS:
        out += [(os.path.join(OLD4, c), os.path.join(ROOT, 'ot', c), EXE_NEW, None) for c in names]
        lists[job] = [os.path.join(ROOT, 'ot', c) for c in names]
    return lists, out


def plan_estimate():
    srcs = [os.path.join(OLD3, 'eta_suite' if c.endswith('_eta2') else 'suite', c) for c in EST_CASES]
    out = [(src, os.path.join(ROOT, 'estimate', os.path.basename(src) + '_' + tag), exe, True)
           for src in srcs for tag, exe in (('old', EXE_OLD_1D), ('new', EXE_NEW))]
    return {'estimate': [o[1] for o in out]}, out


def old_links():
    """runs/stage14/old/ot/<case> -> the earlier OT cases (OT_RUNS, OT_DEC)."""
    d = os.path.join(ROOT, 'old', 'ot')
    links = [(os.path.abspath(os.path.join(OLD4, c)), os.path.join(d, c)) for c in OT_RUNS + OT_DEC]
    return links


PLANS = {'pairs': plan_pairs, 'mhd1d': plan_mhd1d, 'ot': plan_ot, 'estimate': plan_estimate}


def main(cmd):
    lists, todo = PLANS[cmd]()
    links = old_links() if cmd == 'ot' else []
    files = {name: os.path.join('runs', 'cases_s14_%s.txt' % name) for name in lists}
    manifest = os.path.join(ROOT, 'manifest_%s.json' % cmd)
    for src, dst, exe, est in todo:
        check_src(src)
        if os.path.lexists(dst):
            raise SystemExit('%s exists: not written again' % dst)
        if os.path.abspath(dst).startswith(os.path.abspath(OLD3)) or os.path.abspath(dst).startswith(
                os.path.abspath(os.path.dirname(OLD4))):
            raise SystemExit('%s: inside the earlier runs' % dst)
    for target, link in links:
        if not os.path.isdir(target) or os.path.lexists(link):
            raise SystemExit('%s: no target %s or the link exists' % (link, target))
    for f in list(files.values()) + [manifest]:
        if os.path.exists(f):
            raise SystemExit('%s exists: not written again' % f)
    if len({d for _, d, _, _ in todo}) != len(todo):
        raise SystemExit('duplicate case names')
    man = {dst: copy_case(src, dst, exe, est) for src, dst, exe, est in todo}
    for target, link in links:
        os.makedirs(os.path.dirname(link), exist_ok=True)
        os.symlink(target, link)
    for name, cases in lists.items():
        with open(files[name], 'w') as f:
            f.write('\n'.join(sorted(cases, key=cost, reverse=True)) + '\n')
    with open(manifest, 'w') as f:
        json.dump(dict(cases=man, links=dict((l, t) for t, l in links), lists=files), f, indent=1)
    for name, cases in lists.items():
        print('%s: %d cases' % (files[name], len(cases)))
    if links:
        print('%s: %d links' % (os.path.dirname(links[0][1]), len(links)))


if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in PLANS:
        sys.exit(__doc__)
    main(sys.argv[1])
