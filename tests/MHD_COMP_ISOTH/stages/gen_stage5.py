#!/usr/bin/env python3
"""Stage 5 (the hole of rho at M_s = 4) case generator for the Fortran code.

    python3 gen_stage5.py <batch> [base]   writes <base>/<batch>/<case>/ (default base: runs/stage5) as
                                           gen_stage4.py and prints the CASES string for job.pbs-ps-fes-dev

A case directory that exists is not written again (the batch stops with an error).

Batches
    res       turb4 and turb4b (gen_stage4.TURB_PROBS) on 128^3 with the settings of turb4_N64_a and
              turb4b_N64 (runs/stage4/turb): is the hole of rho in one cell along B0 an effect of the
              resolution? Adaptive dt from the start, conservation.dat, budget.dat, cfl.dat at every step and
              out3d every RES_3D in t (see batch_res)
    calib     100 steps of turb2b and turb3b on 64^3 (8 ranks) from t = 0 as gen_stage4.py calib: dt_adv and
              dt_visc at the start for the dt of ms
    ms        turb2b and turb3b (M_s = 2, 3 with rms|dB| = B0 as turb1 and turb4b) on 64^3: the upper limit
              of M_s that runs to t_end on 64^3 with the settings of turb4b_N64 (see batch_ms)
"""
import os
import sys
import math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths
import gen_stage4 as s4

RES_N = 128
RES_NP = 64          # one host of 64 cores per case (job with select=2:ncpus=64, place=scatter, PAR=128)
RES_CALIB = 'turb4_N128_np64'   # runs/stage4/calib: 100 steps of turb4 on 128^3 from t = 0
RES_3D = 0.01        # out3d interval in t (about 120 MB per record on 128^3); the hole formed at t = 0.10-0.18 on 64^3
MS_PROBS = ('turb2b', 'turb3b')
MS_N = 64
MS_CALIB_NP = 8      # as calib/turb1_N64_np8 and turb4_N64_np8 of runs/stage4
MS_3D = 0.1          # out3d interval in tau0 (about 17 MB per record on 64^3, 31 records up to t_end)


def new_case(path):
    """path, which must not exist (a case is not written again)."""
    if os.path.exists(path):
        raise SystemExit('%s exists: not written again' % path)
    return path


def calib_dt(case):
    """TURB_PILOT x the min over the calib run of min(dt_adv, dt_visc) (cfl.dat), as gen_stage4.turb_pilot_nstep.
    The calib run on 128^3 ended with dt_changed (its dt of 3e-4 was above dt_visc), which leaves the limits
    in cfl.dat valid."""
    from analyze_stage1 import load_runs
    cfl = np.vstack(load_runs(case, 'cfl.dat'))
    return s4.TURB_PILOT*float(np.minimum(cfl[:, 2], cfl[:, 3]).min())


def batch_res(root):
    """turb4 and turb4b on RES_N^3 (RES_NP ranks) with reset_method = 'decrement' and the chain up to t_end.
    Both start with the dt of the turb4 calib run on RES_N^3 (on 64^3 turb4b_N64 also took the dt of turb4),
    not with gen_stage4.turb_nstep_from: the min dt_cfl of the runs on 64^3 is set by the hole."""
    base4 = os.path.join(os.path.dirname(os.path.dirname(root)), 'stage4')
    dt = calib_dt(os.path.join(base4, 'calib', RES_CALIB))
    cases = []
    for prob, name in (('turb4', 'turb4_N%d_a' % RES_N), ('turb4b', 'turb4b_N%d' % RES_N)):
        new_case(os.path.join(root, name))
        tau0 = math.pi/s4.TURB_PROBS[prob][0]
        nstep = int(math.ceil(s4.TURB_NTAU*tau0/dt))
        cases.append(s4.turb_case(root, name, prob, RES_N, nstep, RES_NP,
                                  dict(role='resolution', stage=5, dt_from=RES_CALIB, ref64=name.replace('128', '64')),
                                  cons=0.0, out3d=RES_3D/tau0, adaptive=True))
    return cases


def batch_calib(root):
    """gen_stage4.batch_calib for MS_PROBS on MS_N^3 with MS_CALIB_NP ranks (runs/stage5/calib)."""
    probs = tuple((p, MS_N, MS_CALIB_NP) for p in MS_PROBS)
    for p, N, P in probs:
        new_case(os.path.join(root, '%s_N%d_np%d' % (p, N, P)))
    return s4.batch_calib(root, probs, header='Stage 5 calibration of turb2b, turb3b (gen_stage5.py calib)')


def batch_ms(root):
    """MS_PROBS on MS_N^3 with the settings of turb4b_N64 (runs/stage4/turb): s4.TURB_NP[64] ranks,
    reset_method = 'decrement' and the chain up to t_end = 3 tau0, conservation.dat, budget.dat, cfl.dat at
    every step; dt0 = TURB_PILOT x the min dt_cfl of the calib run of the problem (gen_stage4.turb_pilot_nstep
    on runs/stage5/calib). out3d every MS_3D tau0 for the shape of a hole."""
    base = os.path.dirname(root)
    cases = []
    for prob in MS_PROBS:
        name = '%s_N%d' % (prob, MS_N)
        new_case(os.path.join(root, name))
        cases.append(s4.turb_case(root, name, prob, MS_N, s4.turb_pilot_nstep(prob, base), s4.TURB_NP[MS_N],
                                  dict(role='ms_limit', stage=5, dt_from='calib/%s_N%d_np%d' % (prob, MS_N, MS_CALIB_NP)),
                                  cons=0.0, out3d=MS_3D, adaptive=True))
    return cases


BATCHES = {'res': batch_res, 'calib': batch_calib, 'ms': batch_ms}


if __name__ == '__main__':
    batch = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.join(paths.RUNS, 'stage5')
    print(' '.join(sorted(BATCHES[batch](os.path.join(base, batch)), key=s4.cost, reverse=True)))
