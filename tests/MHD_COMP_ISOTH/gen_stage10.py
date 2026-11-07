#!/usr/bin/env python3
"""Stage 10 (dense rerun of the 128^3 control of Stage 8) case generator for the Fortran code.

    python3 gen_stage10.py <batch> [base]   writes <base>/<batch>/<case>/ (default base: runs/stage10) as
                                            gen_stage4.py and prints the CASES string for job.pbs-ps-fes-dev

A case directory that exists is not written again (the batch stops with an error). The control of Stage 8 is
read from runs/stage8 next to base and is only read.

Batches
    dense     turb1 on 128^3 with 'track' with the arguments of the control (stage8/coef/e0.5_C2, CONTROL) except
              out3d every OUT3D tau0 and the stop at T_STOP tau0; stops with an error if calliope.in differs from
              that of the control in anything else or if the files of ic/ differ
"""
import hashlib
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_stage4 as s4
from gen_stage6 import new_case
from analyze_stage1 import read_nml

CONTROL = os.path.join('stage8', 'coef', 'e0.5_C2')
T_STOP = 1.051    # tau0: the record of 1.05 tau0 keeps a following row
OUT3D = 0.01      # tau0
NAME, EPS0, C_SHOCK = 'e0.5_C2', 0.5, 2.0
ALLOWED = {'t_end', 'write_intvl_3d'}  # what may differ from the control (read_nml: lower case)


def md5s(d):
    out = {}
    for name in sorted(os.listdir(d)):
        with open(os.path.join(d, name), 'rb') as f:
            out[name] = hashlib.md5(f.read()).hexdigest()
    return out


def check_control(case, control):
    """Stops unless calliope.in differs from the control in ALLOWED only and ic/ is the same; returns the changes."""
    a, b = read_nml(os.path.join(case, 'calliope.in')), read_nml(os.path.join(control, 'calliope.in'))
    diff = {k: (a.get(k), b.get(k)) for k in sorted(set(a) | set(b)) if a.get(k) != b.get(k)}
    if set(diff) != ALLOWED:
        raise SystemExit('%s: differs from the control in %s' % (case, sorted(diff)))
    ma, mb = md5s(os.path.join(case, 'ic')), md5s(os.path.join(control, 'ic'))
    if ma != mb:
        raise SystemExit('%s: ic/ differs from the control' % case)
    return diff, ma


def batch_dense(root):
    runs = os.path.dirname(os.path.dirname(root))
    nstep = s4.turb_nstep_from(os.path.join(runs, 'stage4', 'turb'), 'turb1', 128, 64)
    new_case(os.path.join(root, NAME))
    c = s4.turb_case(root, NAME, 'turb1', 128, nstep, s4.TURB_NP[128],
                     dict(role='dense', stage=10, ref=CONTROL, svv_eps0=EPS0, C_shock=C_SHOCK, t_stop=T_STOP),
                     cons=0.0, out3d=OUT3D, track=True, C_shock=C_SHOCK, svv_eps0=EPS0, t_stop=T_STOP)
    diff, ic = check_control(c, os.path.join(runs, CONTROL))
    for k, v in diff.items():
        print('# %s: %s (control %s)' % (k, v[0], v[1]), file=sys.stderr)
    print('# ic/ same as the control: %d files' % len(ic), file=sys.stderr)
    return [c]


BATCHES = {'dense': batch_dense}


if __name__ == '__main__':
    batch = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.join('runs', 'stage10')
    print(' '.join(BATCHES[batch](os.path.join(base, batch))))
