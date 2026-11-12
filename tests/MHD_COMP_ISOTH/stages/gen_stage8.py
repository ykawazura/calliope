#!/usr/bin/env python3
"""Stage 8 (dependence of the holes of rho on the dissipation coefficients) case generator for the Fortran code.

    python3 gen_stage8.py <batch> [base]   writes <base>/<batch>/<case>/ (default base: runs/stage8) as
                                           gen_stage4.py and prints the CASES string for job.pbs-ps-fes-dev

A case directory that exists is not written again (the batch stops with an error). The runs of Stage 4 and 6
are read from runs/stage4 and runs/stage6q next to base.

Batches
    coef      turb1 on 128^3 with 'track' up to T_STOP tau0, with the initial state, dt0, ranks and output of
              stage6q/track/turb1_N128_t (REF) except out3d every OUT3D tau0, for (svv_eps0, C_shock) of COEFS:
              the first four change one coefficient by a factor of 2 from the default (0.5, 2), the last
              (the default) is the control
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths
import gen_stage4 as s4
from gen_stage6 import new_case
from analyze_stage1 import read_nml

REF = os.path.join('stage6q', 'track', 'turb1_N128_t')
T_STOP = 1.2      # tau0: past the min of rho_min on 128^3 (about 1.0 tau0) and the second dip on 256^3 (1.18 tau0)
OUT3D = 0.1       # tau0
# (name, svv_eps0, C_shock), by decreasing expected cost (eps0 = 1 halves dt_visc of the SVV); the control last
COEFS = (('e1.0_C2', 1.0, 2.0), ('e0.5_C4', 0.5, 4.0), ('e0.5_C1', 0.5, 1.0), ('e0.25_C2', 0.25, 2.0),
         ('e0.5_C2', 0.5, 2.0))


def batch_coef(root):
    """The cases of COEFS (module docstring); dt0 is that of REF (the nstep of gen_stage6.batch_track), checked
    against its calliope.in."""
    runs = os.path.dirname(os.path.dirname(root))
    nstep = s4.turb_nstep_from(os.path.join(runs, 'stage4', 'turb'), 'turb1', 128, 64)
    cases = []
    for name, eps0, C in COEFS:
        new_case(os.path.join(root, name))
        c = s4.turb_case(root, name, 'turb1', 128, nstep, s4.TURB_NP[128],
                         dict(role='coef', stage=8, ref=REF, svv_eps0=eps0, C_shock=C, t_stop=T_STOP),
                         cons=0.0, out3d=OUT3D, track=True, C_shock=C, svv_eps0=eps0, t_stop=T_STOP)
        ref = read_nml(os.path.join(runs, REF, 'calliope.in'))
        if read_nml(os.path.join(c, 'calliope.in'))['dt'] != ref['dt']:
            raise SystemExit('%s: dt differs from %s' % (c, REF))
        cases.append(c)
    return cases


BATCHES = {'coef': batch_coef}


if __name__ == '__main__':
    batch = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.join(paths.RUNS, 'stage8')
    print(' '.join(sorted(BATCHES[batch](os.path.join(base, batch)), key=s4.cost, reverse=True)))
