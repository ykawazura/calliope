#!/usr/bin/env python3
"""Stage 6 (dt that follows dt_cfl, turb1 on 256^3) case generator for the Fortran code.

    python3 gen_stage6.py <batch> [base]   writes <base>/<batch>/<case>/ (default base: runs/stage6) as
                                           gen_stage4.py and prints the CASES string for job.pbs-ps-fes-dev

A case directory that exists is not written again (the batch stops with an error). The runs of Stage 4 are
read from runs/stage4 next to base.

Batches
    regress   bitwise regression of the change of src/time.F90 (reset_method = 'track', t_end) and
              src/calliope.F90: the cases of gen_stage4.batch_regress run by EXE_OLD and EXE_NEW, turb1 on 64^3
              in which 'decrement' reduces dt (DEC_*), and the MHD_INCOMP input incomp_regress.in (next to this
              file) run by INCOMP_OLD and INCOMP_NEW; FFTW_ESTIMATE
    track     turb1 on 64^3 with 'track' from dt0 = TDEC_DT0 dt_cfl(0) over TDEC_NSTEP steps on P = TDEC_P
              (decomposition; dt is raised at step 100), and turb1 with 'track' up to t_end on 64^3 and 128^3
              with the dt0, ranks and output of turb1_N64 and turb1_N128_a (runs/stage4/turb)
    turb      turb1 on 256^3 (N256_NP ranks) with 'track' up to t_end, output as turb1_N128_a
"""
import os
import sys
import math
import shutil
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_stage4 as s4

EXE_OLD = 'bin/calliope_isoth_s4b'     # src of Stage 5b (HEAD 0f3c805)
EXE_NEW = 'bin/calliope_isoth_s6'
INCOMP_OLD = 'bin/calliope_incomp_base'
INCOMP_NEW = 'bin/calliope_incomp_s6'
INCOMP_INPUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'incomp_regress.in')
INCOMP_NP = 4
DEC_DT0 = 2.0       # 'decrement' regression: dt0 = DEC_DT0 dt_cfl(0), so that reset_dt reduces dt at the first step
DEC_NSTEP = 100
DEC_NP = 8
TDEC_DT0 = 0.3      # 'track' decomposition: dt0 = TDEC_DT0 dt_cfl(0) < track_raise dt_cfl(0)
TDEC_NSTEP = 400
TDEC_P = (1, 8, 64)
N256_NP = 256       # one node
# dt_cfl(128^3)/dt_cfl(256^3) assumed equal to the min dt_visc of turb1 on 64^3 over that on 128^3 (REPORT_stage4 8.3)
N256_DT_RATIO = 2.58
N256_DT_SAFETY = 0.8   # track_safety of src/time.F90


def new_case(path):
    """path, which must not exist (a case is not written again)."""
    if os.path.exists(path):
        raise SystemExit('%s exists: not written again' % path)
    return path


def dt_cfl0(case):
    """min(dt_adv, dt_visc) of the first row of cfl.dat of a run (its first step, from t = 0)."""
    from analyze_stage1 import load_runs
    cfl = load_runs(case, 'cfl.dat')[0]
    return float(min(cfl[0, 2], cfl[0, 3]))


def turb_nstep(prob, dt):
    """nstep with t_end/nstep <= dt (t_end = TURB_NTAU tau0 of prob)."""
    return int(math.ceil(s4.TURB_NTAU*math.pi/s4.TURB_PROBS[prob][0]/dt))


def pair(root, c, old, new):
    """Copies <root>/<c>_old and <root>/<c>_new of the case directory c with the file "exe" (old, new) and
    FFTW_ESTIMATE."""
    out = []
    for tag, exe in (('old', old), ('new', new)):
        dst = new_case(os.path.join(root, os.path.basename(c) + '_' + tag))
        shutil.copytree(c, dst)
        with open(os.path.join(dst, 'exe'), 'w') as f:
            f.write(exe + '\n')
        open(os.path.join(dst, 'fftw_estimate'), 'w').close()
        out.append(dst)
    return out


def batch_regress(root):
    """gen_stage4.batch_regress with EXE_OLD, EXE_NEW; turb1 on 64^3 with 'decrement' from DEC_DT0 dt_cfl(0)
    (DEC_NSTEP steps, DEC_NP ranks: reset_dt reduces dt before the first step, and later when dt_cfl falls below
    dt); the MHD_INCOMP regression input (INCOMP_NP ranks) with INCOMP_OLD, INCOMP_NEW."""
    base4 = os.path.join(os.path.dirname(os.path.dirname(root)), 'stage4')
    new_case(root)
    s4.EXE_OLD, s4.EXE_NEW = EXE_OLD, EXE_NEW
    cases = s4.batch_regress(root)
    src = os.path.join(root, 'src')
    nstep = turb_nstep('turb1', DEC_DT0*dt_cfl0(os.path.join(base4, 'calib', 'turb1_N64_np8')))
    c = s4.turb_case(src, 'dec_turb1_N64', 'turb1', 64, nstep, DEC_NP, dict(role='regress', stage=6),
                     nrun=DEC_NSTEP, planes=0.0, cons=0.0, out3d=0.0, adaptive=True)
    # the run ends after DEC_NSTEP steps: without the file t_end the job script does not continue it
    os.remove(os.path.join(c, 't_end'))
    cases += pair(root, c, EXE_OLD, EXE_NEW)
    c = os.path.join(src, 'incomp_ot16')
    os.makedirs(c)
    shutil.copy(INCOMP_INPUT, os.path.join(c, 'calliope.in'))
    with open(os.path.join(c, 'nproc'), 'w') as f:
        f.write('%d\n' % INCOMP_NP)
    cases += pair(root, c, INCOMP_OLD, INCOMP_NEW)
    return cases


def batch_track(root):
    """dec_N64_np<P>: turb1 on 64^3 with 'track' from TDEC_DT0 dt_cfl(0) over TDEC_NSTEP steps on P = TDEC_P
    ranks (conservation.dat and cfl.dat at every step, no planes, no out3d). turb1_N64_t, turb1_N128_t: turb1
    with 'track' up to t_end with the dt0, ranks and output of turb1_N64 and turb1_N128_a, cons at every step."""
    base4 = os.path.join(os.path.dirname(os.path.dirname(root)), 'stage4')
    nstep = turb_nstep('turb1', TDEC_DT0*dt_cfl0(os.path.join(base4, 'calib', 'turb1_N64_np8')))
    cases = []
    for P in TDEC_P:
        name = 'dec_N64_np%d' % P
        new_case(os.path.join(root, name))
        cases.append(s4.turb_case(root, name, 'turb1', 64, nstep, P, dict(role='decomp', stage=6), nrun=TDEC_NSTEP,
                                  planes=0.0, cons=0.0, out3d=0.0, track=True))
    for N, nstep, ref in ((64, s4.turb_pilot_nstep('turb1', base4), 'turb1_N64'),
                          (128, s4.turb_nstep_from(os.path.join(base4, 'turb'), 'turb1', 128, 64), 'turb1_N128_a')):
        name = 'turb1_N%d_t' % N
        new_case(os.path.join(root, name))
        cases.append(s4.turb_case(root, name, 'turb1', N, nstep, s4.TURB_NP[N],
                                  dict(role='track', stage=6, ref='stage4/turb/' + ref), cons=0.0, track=True))
    return cases


def batch_turb(root):
    """turb1_N256: turb1 on 256^3 (N256_NP ranks) with 'track' up to t_end, cons at every step, planes and out3d
    as turb1_N128_a; dt0 = N256_DT_SAFETY dt_cfl(0) of turb1_N128_a/N256_DT_RATIO ('track' corrects it at the
    first step if it exceeds dt_cfl, and raises it at step 100 if it is below track_raise dt_cfl)."""
    base4 = os.path.join(os.path.dirname(os.path.dirname(root)), 'stage4')
    name = 'turb1_N256'
    new_case(os.path.join(root, name))
    dt0 = N256_DT_SAFETY*dt_cfl0(os.path.join(base4, 'turb', 'turb1_N128_a'))/N256_DT_RATIO
    return [s4.turb_case(root, name, 'turb1', 256, turb_nstep('turb1', dt0), N256_NP,
                         dict(role='production', stage=6, dt_from='stage4/turb/turb1_N128_a'), cons=0.0, track=True)]


BATCHES = {'regress': batch_regress, 'track': batch_track, 'turb': batch_turb}


if __name__ == '__main__':
    batch = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.join('runs', 'stage6')
    print(' '.join(sorted(BATCHES[batch](os.path.join(base, batch)), key=s4.cost, reverse=True)))
