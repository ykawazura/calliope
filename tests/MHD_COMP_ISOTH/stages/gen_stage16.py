#!/usr/bin/env python3
"""Stage 16: the 1D conv runs of Stage 3 (runs/stage3_f/conv) again with the final binary gen_stage14.EXE_NEW, and
the HF runs of the second 256^3 reference (H of gen_stage15 with filter_lnrho = F as ES, turb1).

    python3 tests/MHD_COMP_ISOTH/stages/gen_stage16.py conv    runs/stage16/conv/<case> of CONV (fftw_estimate as in the
                                                                earlier case)
    python3 tests/MHD_COMP_ISOTH/stages/gen_stage16.py hf64    F1: runs/stage16/hf64/turb1_N64_HF_r<r>, r of
                                                                gen_stage15.R_H
    python3 tests/MHD_COMP_ISOTH/stages/gen_stage16.py hf128   F2: HF(r), HF(2r), HF(4r) on 128^3 with r = r_star of
                                                                runs/stage16/hf64.json (analyze_stage16.py hf hf64)
    python3 tests/MHD_COMP_ISOTH/stages/gen_stage16.py hf256   F3: HF(r) on 256^3 with r = r_C of
                                                                runs/stage16/hf128.json (analyze_stage16.py hf hf128)

Run in the top directory of the code (ps1: ~/work/calliope_isoth). As gen_stage14.py, conv checks all its cases
before it writes anything and stops if a case directory, the list or the manifest exists or an earlier case holds a
file it does not know. It copies the input files with gen_stage14.copy_case (md5 checked) and writes the list
runs/cases_s16_conv.txt (cases by decreasing cost) and runs/stage16/manifest_conv.json. Nothing is written in
runs/stage3_f. The hf batches write their cases with gen_stage15.run_case (eps_ref of runs/stage15/pre.json, which
is only read), check each calliope.in against the ES run of gen_stage15.REF on the same N (it may differ in t_end,
write_intvl_3d, svv_eps0 and the hyperviscosity only; ic/ the same) and write the list runs/cases_s16_<batch>.txt
(largest r first: the small rungs, which may stop at the end of the budget, start last).
"""
import sys
import os
import json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_stage14 as g14
import gen_stage15 as g15
from gen_stage3 import cost

ROOT = os.path.join('runs', 'stage16')
S15 = os.path.join('runs', 'stage15')
OLD_CONV = os.path.join(g14.OLD3, 'conv')
CONV = tuple('%s_N%d%s' % (p, N, t) for p in ('I1', 'RDpc') for N in (4096, 8192) for t in ('', '_eta2'))


def plan_conv():
    srcs = g14.cases_in(OLD_CONV)
    if sorted(os.path.basename(s) for s in srcs) != sorted(CONV):
        raise SystemExit('%s: cases %s, expected %s' % (OLD_CONV, [os.path.basename(s) for s in srcs], CONV))
    return [(s, os.path.join(ROOT, 'conv', os.path.basename(s))) for s in srcs]


PLANS = dict(conv=plan_conv)


def gate(name, key):
    """key of runs/stage16/<name>, written by analyze_stage16.py hf."""
    with open(os.path.join(ROOT, name)) as f:
        return json.load(f)[key]


def hf_cases(N, rs):
    """The HF cases of turb1 on N^3 for the rungs rs in runs/stage16/hf<N>."""
    eps = g15.eps_of(os.path.join(S15, 'a'))['turb1']     # runs/stage15/pre.json, gates of phase 0 passed
    root = os.path.join(ROOT, 'hf%d' % N)
    cases = []
    for r in rs:
        c = g15.run_case(root, 'turb1', N, 'HF', r, eps, dict(role='hf%d' % N, stage=16))
        g15.check_ref(c, N, extra={'svv_eps0'} | {k.lower() for k in g15.coefs('HF', r, N, eps)[1]})
        cases.append(c)
    return cases


def batch_hf64():
    return hf_cases(64, g15.R_H)


def batch_hf128():
    r = gate('hf64.json', 'r_star')
    if r is None:
        raise SystemExit('hf64.json: no r*_HF on 64^3 (consult the user)')
    return hf_cases(128, (r, 2*r, 4*r))


def batch_hf256():
    r = gate('hf128.json', 'r_C')
    if r is None:
        raise SystemExit('hf128.json: no complete rung on 128^3 (stop and consult the user)')
    return hf_cases(256, (r,))


BATCHES = dict(hf64=batch_hf64, hf128=batch_hf128, hf256=batch_hf256)


def main_hf(cmd):
    lst = os.path.join('runs', 'cases_s16_%s.txt' % cmd)
    if os.path.exists(lst):
        raise SystemExit('%s exists: not written again' % lst)
    cases = BATCHES[cmd]()
    with open(lst, 'w') as f:
        f.write('\n'.join(sorted(cases, key=lambda c: json.load(open(os.path.join(c, 'case.json')))['r'],
                                 reverse=True)) + '\n')
    print('%s: %d cases' % (lst, len(cases)))


def main(cmd):
    todo = PLANS[cmd]()
    lst = os.path.join('runs', 'cases_s16_%s.txt' % cmd)
    manifest = os.path.join(ROOT, 'manifest_%s.json' % cmd)
    for src, dst in todo:
        g14.check_src(src)
        if os.path.lexists(dst):
            raise SystemExit('%s exists: not written again' % dst)
        if os.path.abspath(dst).startswith(os.path.abspath(g14.OLD3) + os.sep):
            raise SystemExit('%s: inside the earlier runs' % dst)
    for f in (lst, manifest):
        if os.path.exists(f):
            raise SystemExit('%s exists: not written again' % f)
    man = {dst: g14.copy_case(src, dst, g14.EXE_NEW, None) for src, dst in todo}
    with open(lst, 'w') as f:
        f.write('\n'.join(sorted((d for _, d in todo), key=cost, reverse=True)) + '\n')
    with open(manifest, 'w') as f:
        json.dump(dict(cases=man, list=lst), f, indent=1)
    print('%s: %d cases' % (lst, len(todo)))


if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in set(PLANS) | set(BATCHES):
        sys.exit(__doc__)
    (main if sys.argv[1] in PLANS else main_hf)(sys.argv[1])
