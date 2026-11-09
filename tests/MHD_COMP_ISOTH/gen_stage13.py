#!/usr/bin/env python3
"""Stage 13 (ES-SVV without the mass fluxes, svv_mass = F) case generator for the Fortran code.

    python3 gen_stage13.py <batch> [base]   writes <base>/<batch>/<case>/ (default base: runs/stage13d) as
                                            gen_stage4.py and prints the CASES string for job.pbs-ps-fes-dev

A case directory that exists is not written again (the batch stops with an error). The runs of Stage 4 and the
control of Stage 8 are read from runs/stage4 and runs/stage8 next to base and are only read. Every case has the
file "exe" (EXE_NEW, or the pair of gen_stage6.batch_regress). runs/stage13 and bin/*_s13 are the first
verification (job 486708, failed), before the m part of the initial vector of svv_power for svv_mass = F;
runs/stage13b and bin/*_s13b the second (job 487523, failed in V1), with that part inside the loop of the default
path; bin/*_s13c the third (failed in the check of svv_power's machine code before the job; no runs), with that part
in a loop of its own in svv_power. Here (s13d) the initial vector for svv_mass = F is set by svv_seed_nomass, and
svv_power is as in Stage 6.

Batches
    regress   gen_stage6.batch_regress (the 7 ISOTH inputs and the MHD_INCOMP input, FFTW_ESTIMATE) with the
              binaries of Stage 6 (EXE_OLD, INCOMP_OLD) and of this stage (EXE_NEW, INCOMP_NEW); svv_mass is not
              set, so the default path must be bitwise the same
    track     dec_N64_np<P>: gen_stage6.batch_track (turb1 on 64^3 with 'track' from TDEC_DT0 dt_cfl(0) over
              TDEC_NSTEP steps on P = 1, 8, 64) with svv_mass = F (decomposition)
    nomass    turb1 on 128^3 with the arguments of the control (stage8/coef/e0.5_C2, CONTROL) and svv_mass = F;
              stops with an error if calliope.in differs from that of the control in anything else or if the
              files of ic/ differ
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_stage4 as s4
import gen_stage6 as s6
import gen_stage8 as s8
from gen_stage10 import check_control

EXE_OLD = 'bin/calliope_isoth_s6'     # src of Stage 6 (0c374aa), the same as that of 35c6ed0
EXE_NEW = 'bin/calliope_isoth_s13d'
INCOMP_OLD = 'bin/calliope_incomp_s6'
INCOMP_NEW = 'bin/calliope_incomp_s13d'
CONTROL = os.path.join('stage8', 'coef', 'e0.5_C2')
NAME, EPS0, C_SHOCK = 'e0.5_C2', 0.5, 2.0


def set_exe(c, exe=EXE_NEW):
    with open(os.path.join(c, 'exe'), 'w') as f:
        f.write(exe + '\n')
    return c


def batch_regress(root):
    s6.EXE_OLD, s6.EXE_NEW, s6.INCOMP_OLD, s6.INCOMP_NEW = EXE_OLD, EXE_NEW, INCOMP_OLD, INCOMP_NEW
    return s6.batch_regress(root)


def batch_track(root):
    base4 = os.path.join(os.path.dirname(os.path.dirname(root)), 'stage4')
    nstep = s6.turb_nstep('turb1', s6.TDEC_DT0*s6.dt_cfl0(os.path.join(base4, 'calib', 'turb1_N64_np8')))
    cases = []
    for P in s6.TDEC_P:
        name = 'dec_N64_np%d' % P
        s6.new_case(os.path.join(root, name))
        cases.append(set_exe(s4.turb_case(root, name, 'turb1', 64, nstep, P, dict(role='decomp', stage=13),
                                          nrun=s6.TDEC_NSTEP, planes=0.0, cons=0.0, out3d=0.0, track=True,
                                          svv_mass=False)))
    return cases


def batch_nomass(root):
    runs = os.path.dirname(os.path.dirname(root))
    nstep = s4.turb_nstep_from(os.path.join(runs, 'stage4', 'turb'), 'turb1', 128, 64)
    s6.new_case(os.path.join(root, NAME))
    c = s4.turb_case(root, NAME, 'turb1', 128, nstep, s4.TURB_NP[128],
                     dict(role='nomass', stage=13, ref=CONTROL, svv_eps0=EPS0, C_shock=C_SHOCK, t_stop=s8.T_STOP),
                     cons=0.0, out3d=s8.OUT3D, track=True, C_shock=C_SHOCK, svv_eps0=EPS0, t_stop=s8.T_STOP,
                     svv_mass=False)
    diff, ic = check_control(c, os.path.join(runs, CONTROL), allowed={'svv_mass'})
    for k, v in diff.items():
        print('# %s: %s (control %s)' % (k, v[0], v[1]), file=sys.stderr)
    print('# ic/ same as the control: %d files' % len(ic), file=sys.stderr)
    return [set_exe(c)]


BATCHES = {'regress': batch_regress, 'track': batch_track, 'nomass': batch_nomass}


if __name__ == '__main__':
    batch = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.join('runs', 'stage13d')
    print(' '.join(sorted(BATCHES[batch](os.path.join(base, batch)), key=s4.cost, reverse=True)))
