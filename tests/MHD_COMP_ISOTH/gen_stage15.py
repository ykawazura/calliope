#!/usr/bin/env python3
"""Stage 15 (ES-SVV against the standard dissipation in 3D decaying turbulence) case generator for the Fortran code.

    python3 gen_stage15.py <batch> [base]   writes <base>/<batch>/<case>/ (default base: runs/stage15) as
                                            gen_stage4.py and prints the CASES string for job.pbs-ps-fes-dev

A case directory that exists is not written again (the batch stops with an error). The runs of Stages 4, 6 and 8
are read from runs/stage4, runs/stage6q and runs/stage8 next to base and are only read. Every case has the file
"exe" (EXE). No Fortran change: the schemes differ in the namelists only.

Schemes (C_shock = 2, C_eta = 0, 'track' from the dt0 of turb1 on the same N, t_end = T_STOP tau0, out3d every
OUT3D tau0, conservation.dat, budget.dat and cfl.dat at every step)
    ES        svv_eps0 = 0.5, filter_lnrho = F (the production scheme)
    ESf       ES with filter_lnrho = T (separates the filter from the SVV)
    ES25      svv_eps0 = 0.25 (64^3 only)
    R0        svv_eps0 = 0, nu = eta = 0, filter_lnrho = T
    L_r<r>    svv_eps0 = 0, filter_lnrho = T, nu = eta = (27/4) gamma_c, lmd = 0 (Laplacian of m and B)
    H_r<r>    svv_eps0 = 0, filter_lnrho = T, nu_h = eta_h = lmd_h = (27/4)^4 gamma_c, exponents 4 (k^8)
    HF_r<r>   H_r<r> with filter_lnrho = F as ES (Stage 16, gen_stage16.py)
The decay rate of advance.F90 is nu k^2/k2_max + nu_h (k^2/k2_max)^4 with k2_max = 3 (N/2)^2 (src/grid.F90), so
that at k_c = N/3 (k_c^2/k2_max = 4/27) both families decay at gamma_c = r eps_ref^(1/3) k_c^(2/3). eps_ref is the
peak of -dE/dt of ES on 64^3 over (2 pi)^3, per M_s (pre.json, analyze_stage15.py pre).

Batches
    pre   phase 0 on 64^3 (M_s = 1 unless noted): dt<k>_*: ES, L_r1, H_r1 with the fixed dt0/k of turb1_N64
          (k = 1, 2; reset_method = 'stop') up to P0A_T tau0 on P0A_NP ranks; dec_*_np<P>: L_r1, H_r1 over
          P0B_NSTEP steps with dt0 on P = P0B_NP ranks, FFTW_ESTIMATE, planes at steps 0 and P0B_NSTEP;
          <prob>_N64_ES for both M_s (also the ES of phase A). L_r1 and H_r1 here take eps_ref of the Stage 6 run
          with the arguments of turb1_N64_ES (EPS_PROV_REF): eps_ref of pre.json is fixed after this batch.
    a     phase A on 64^3 (64 ranks), both M_s: ESf, ES25, R0, L_r<r> (R_L), H_r<r> (R_H)
    a4    the largest rung of a family with no complete rung on 64^3 with C_shock = 4 (a.json: diag); a
          diagnostic only
    b     phase B on 128^3 (128 ranks), per M_s: ES, L(max(1, r*_L)), L(2 max(1, r*_L)), H(r*_H), H(2 r*_H),
          H(4 r*_H) (a.json: r_star), ESf for M_s = 1, R0 where it completed on 64^3 (a.json: r0_ok)
    c     phase C: turb1 on 256^3 (256 ranks) with the scheme of r_256() (r_C of b.json, or r_256 of c.json), up
          to T_STOP tau0 (the job script chains it; analyze_stage15.py watch projects the cost after each job)
"""
import json
import math
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_stage4 as s4
import gen_stage6 as s6
from gen_stage6 import new_case
from gen_stage10 import md5s
from gen_stage13 import set_exe
from analyze_stage1 import read_nml, num

EXE = 'bin/calliope_isoth_s13d'
PROBS = ('turb05', 'turb1')     # M_s = 0.5, 1
T_STOP = 2.0                    # tau0
OUT3D = 0.125                   # tau0
NP = {64: 64, 128: 128, 256: 256}
R_L = (0.25, 0.5, 1.0, 2.0, 4.0)
R_H = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
KMAX2_KC2 = 27.0/4.0            # k2_max/k_c^2 = 3 (N/2)^2/(N/3)^2
P0A_T = 0.2                     # tau0
P0A_NP = 16
P0B_NP = (8, 64)
P0B_NSTEP = 10
C_SHOCK_DIAG = 4.0
# the runs whose calliope.in each new ES case must match (dt0, ic/; t_end and out3d differ), relative to runs/
REF = {64: os.path.join('stage6q', 'track', 'turb1_N64_t'), 128: os.path.join('stage8', 'coef', 'e0.5_C2'),
       256: os.path.join('stage6q', 'turb', 'turb1_N256')}
EPS_PROV_REF = REF[64]
ALLOWED_ES = {'t_end', 'write_intvl_3d'}


def label(scheme, r=None):
    return scheme if r is None else '%s_r%g' % (scheme, r)


def coefs(scheme, r, N, eps_ref):
    """(svv_eps0, phys, num, info) of a scheme (module docstring) on N^3."""
    info = dict(scheme=scheme, r=r)
    if scheme in ('ES', 'ESf', 'ES25'):
        return (0.25 if scheme == 'ES25' else 0.5), None, (dict(filter_lnrho=True) if scheme == 'ESf' else None), info
    if scheme == 'R0':
        return 0.0, None, dict(filter_lnrho=True), info
    kc = N/3.0
    gc = r*eps_ref**(1.0/3.0)*kc**(2.0/3.0)
    if scheme == 'L':
        c = KMAX2_KC2*gc
        phys = dict(nu=c, eta=c, lmd=0.0)
    else:
        c = KMAX2_KC2**4*gc
        phys = dict(nu_h=c, eta_h=c, lmd_h=c, nu_h_exp=4, eta_h_exp=4, lmd_h_exp=4)
    info.update(gamma_c=gc, eps_ref=eps_ref, k_c=kc, coef=c)
    return 0.0, phys, (None if scheme == 'HF' else dict(filter_lnrho=True)), info


def check_coefs(c, svv_eps0, phys, nm):
    """Stops unless calliope.in of c holds svv_eps0, phys and nm exactly."""
    d = read_nml(os.path.join(c, 'calliope.in'))
    want = dict(svv_eps0=svv_eps0, **(phys or {}), filter_lnrho=bool(nm and nm.get('filter_lnrho')))
    for k, v in want.items():
        got = d[k.lower()]
        ok = (got == ('T' if v else 'F')) if isinstance(v, bool) else num(got) == v
        if not ok:
            raise SystemExit('%s: %s = %s, not %r' % (c, k, got, v))


def runs_of(root):
    return os.path.dirname(os.path.dirname(root))


def nstep_turb1(runs, N):
    """nstep of turb1 on N^3 as the reference runs (REF): dt0 = TURB_NTAU tau0/nstep."""
    base4 = os.path.join(runs, 'stage4')
    if N == 64:
        return s4.turb_pilot_nstep('turb1', base4)
    if N == 128:
        return s4.turb_nstep_from(os.path.join(base4, 'turb'), 'turb1', 128, 64)
    return s6.turb_nstep('turb1', s6.N256_DT_SAFETY*s6.dt_cfl0(os.path.join(base4, 'turb', 'turb1_N128_a'))
                         / s6.N256_DT_RATIO)


def nstep_of(runs, prob, N):
    """nstep with the dt0 of turb1 on N^3: tau0 = pi/M_s, so nstep scales as 1/M_s."""
    return int(round(nstep_turb1(runs, N)*s4.TURB_PROBS['turb1'][0]/s4.TURB_PROBS[prob][0]))


def eps_peak_ref(case):
    """max(-dE/dt) of budget.dat over (2 pi)^3."""
    from analyze_stage1 import load_runs
    bud = np.concatenate(load_runs(case, 'budget.dat'))
    return float((-bud[:, 6]).max())/(2.0*math.pi)**3


def run_case(root, prob, N, scheme, r, eps_ref, info, nproc=None, C_shock=2.0, name=None, **kw):
    """One case <prob>_N<N>_<label> ('track' up to T_STOP tau0 unless kw says otherwise); returns its path. The
    stage of case.json is 15 unless info has another."""
    name = name or '%s_N%d_%s' % (prob, N, label(scheme, r))
    new_case(os.path.join(root, name))
    svv_eps0, phys, nm, ci = coefs(scheme, r, N, eps_ref)
    args = dict(cons=0.0, out3d=OUT3D, track=True, t_stop=T_STOP)
    args.update(kw)
    c = s4.turb_case(root, name, prob, N, nstep_of(runs_of(root), prob, N), nproc or NP[N],
                     dict(info, stage=info.get('stage', 15), **ci, C_shock=C_shock), C_shock=C_shock,
                     svv_eps0=svv_eps0, phys=phys, num=nm, **args)
    check_coefs(c, svv_eps0, phys, nm)
    return set_exe(c, EXE)


def check_ref(c, N, extra=()):
    """Stops unless calliope.in of c differs from that of the reference run on N^3 (REF) in ALLOWED_ES | extra
    only and ic/ is the same; prints the differences."""
    ref = os.path.join(runs_of(os.path.dirname(c)), REF[N])
    a, b = read_nml(os.path.join(c, 'calliope.in')), read_nml(os.path.join(ref, 'calliope.in'))
    diff = {k: (a.get(k), b.get(k)) for k in sorted(set(a) | set(b)) if a.get(k) != b.get(k)}
    if not set(diff) <= ALLOWED_ES | set(extra):
        raise SystemExit('%s: differs from %s in %s' % (c, REF[N], sorted(set(diff) - ALLOWED_ES - set(extra))))
    ic = md5s(os.path.join(c, 'ic'))
    if ic != md5s(os.path.join(ref, 'ic')):
        raise SystemExit('%s: ic/ differs from that of %s' % (c, REF[N]))
    print('# %s: %s; ic/ same as %s (%d files)' % (os.path.basename(c), ', '.join(
        '%s %s (ref %s)' % (k, *v) for k, v in diff.items()), REF[N], len(ic)), file=sys.stderr)


def load(root, name):
    with open(os.path.join(os.path.dirname(root), name)) as f:
        return json.load(f)


def eps_of(root):
    """eps_ref per M_s of pre.json, after the gates of phase 0 (P0b, both ES runs complete) passed."""
    pre = load(root, 'pre.json')
    if not pre.get('ok'):
        raise SystemExit('pre.json: the gates of phase 0 did not pass')
    return pre['eps_ref']


def batch_pre(root):
    runs = runs_of(root)
    eps = eps_peak_ref(os.path.join(runs, EPS_PROV_REF))
    print('# provisional eps_ref (M_s = 1) %.6e from %s' % (eps, EPS_PROV_REF), file=sys.stderr)
    n1 = nstep_turb1(runs, 64)
    nrun = int(math.ceil(P0A_T/s4.TURB_NTAU*n1))
    cases = []
    for scheme, r in (('ES', None), ('L', 1.0), ('H', 1.0)):
        for k in (1, 2):
            name = 'dt%d_turb1_N64_%s' % (k, label(scheme, r))
            new_case(os.path.join(root, name))
            svv_eps0, phys, nm, ci = coefs(scheme, r, 64, eps)
            c = s4.turb_case(root, name, 'turb1', 64, k*n1, P0A_NP, dict(ci, role='dt', stage=15, dt_div=k,
                             eps_ref_from=EPS_PROV_REF), nrun=k*nrun, planes=0.0, cons=0.0, out3d=0.0,
                             svv_eps0=svv_eps0, phys=phys, num=nm)
            check_coefs(c, svv_eps0, phys, nm)
            cases.append(set_exe(c, EXE))
    dt = s4.TURB_NTAU*math.pi/n1
    for scheme in ('L', 'H'):
        for P in P0B_NP:
            name = 'dec_turb1_N64_%s_np%d' % (label(scheme, 1.0), P)
            new_case(os.path.join(root, name))
            svv_eps0, phys, nm, ci = coefs(scheme, 1.0, 64, eps)
            c = s4.turb_case(root, name, 'turb1', 64, n1, P, dict(ci, role='decomp', stage=15,
                             eps_ref_from=EPS_PROV_REF), nrun=P0B_NSTEP, planes=s4.line_intvl(P0B_NSTEP, dt, 1),
                             cons=0.0, out3d=0.0, svv_eps0=svv_eps0, phys=phys, num=nm)
            check_coefs(c, svv_eps0, phys, nm)
            open(os.path.join(c, 'fftw_estimate'), 'w').close()
            cases.append(set_exe(c, EXE))
    for prob in PROBS:
        c = run_case(root, prob, 64, 'ES', None, None, dict(role='es'))
        if prob == 'turb1':
            check_ref(c, 64)
        cases.append(c)
    return cases


def batch_a(root):
    eps = eps_of(root)
    cases = []
    for prob in PROBS:
        for scheme, rs in (('ESf', (None,)), ('ES25', (None,)), ('R0', (None,)), ('L', R_L), ('H', R_H)):
            cases += [run_case(root, prob, 64, scheme, r, eps[prob], dict(role='screen')) for r in rs]
    return cases


def batch_a4(root):
    eps, a = eps_of(root), load(root, 'a.json')
    cases = []
    for prob, scheme in a['diag']:
        r = max(R_L if scheme == 'L' else R_H)
        cases.append(run_case(root, prob, 64, scheme, r, eps[prob], dict(role='diag'), C_shock=C_SHOCK_DIAG,
                              name='%s_N64_%s_C4' % (prob, label(scheme, r))))
    return cases


def b_schemes(a, prob):
    """(scheme, r) of phase B for prob from a.json."""
    out = [('ES', None)] + ([('ESf', None)] if prob == 'turb1' else [])
    out += [('R0', None)] if a['r0_ok'][prob] else []
    rL, rH = a['r_star'][prob]['L'], a['r_star'][prob]['H']
    if rL is not None:
        out += [('L', max(1.0, rL)), ('L', 2*max(1.0, rL))]
    if rH is not None:
        out += [('H', rH), ('H', 2*rH), ('H', 4*rH)]
    return out


def batch_b(root):
    eps, a = eps_of(root), load(root, 'a.json')
    cases = []
    for prob in PROBS:
        for scheme, r in b_schemes(a, prob):
            c = run_case(root, prob, 128, scheme, r, eps[prob], dict(role='compare'))
            if prob == 'turb1' and scheme == 'ES':
                check_ref(c, 128)
            cases.append(c)
    return cases


def r_256(root):
    """(scheme, r) of the 256^3 reference: r_256 of c.json, written by hand with the user's choice after a rung
    stopped on 256^3 (H_r8 at 0.614 tau0, 2026-11-11: H_r32), if it exists; else r_C of b.json."""
    if os.path.exists(os.path.join(os.path.dirname(root), 'c.json')):
        return load(root, 'c.json')['r_256']
    return load(root, 'b.json')['r_C']


def batch_c(root):
    eps = eps_of(root)
    scheme, r = r_256(root)
    c = run_case(root, 'turb1', 256, scheme, r, eps['turb1'], dict(role='reference'))
    check_ref(c, 256, extra={'svv_eps0', 'filter_lnrho'} | {k.lower() for k in coefs(scheme, r, 256, 1.0)[1]})
    return [c]


BATCHES = {'pre': batch_pre, 'a': batch_a, 'a4': batch_a4, 'b': batch_b, 'c': batch_c}


if __name__ == '__main__':
    batch = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.join('runs', 'stage15')
    print(' '.join(sorted(BATCHES[batch](os.path.join(base, batch)), key=s4.cost, reverse=True)))
