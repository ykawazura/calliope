#!/usr/bin/env python3
"""Stage 1 (1D isothermal Euler, B = 0) case generator.

    python3 gen_stage1.py <batch> [base]

writes <base>/<batch>/<case>/{calliope.in,nproc} (default base: runs/stage1)
and prints the CASES string for job.pbs-ps-fes-dev.
cs = 1 (beta0 = 2), lx = 1 throughout, i.e. L = 2 pi (grid.F90 multiplies lx by 2 pi).
Every run has a fixed dt with reset_method = 'stop' (time_parameters), so a run that
violates min(dt_adv, dt_visc) stops instead of silently changing dt. dt is set from
DT512, the minimum of min(dt_adv, dt_visc) over the run of proto1d.py at N = 512 with
the default parameters, times SAFETY, scaled by 512/N (both limits scale as 1/N).

Batches
    smoke   sound wave, NaN injection, colliding flow (default and C = 1 without SVV), ζ self test
    smooth  sound wave (N = 32, 64) and simple wave before breaking (N = 32 ... 512)
    shock   U1, U2, U4, Sod 10:1 at N = 256 ... 4096 with the defaults
    scan    U2 and Sod at N = 1024: C_shock, SVV off / Q(>k_c) = 0, filter_lnrho, mass_correction, pruned
    dtconv  U2 at N = 512 with dt / 1, 2, 4, 8 (time-integration error, mass error ~ dt^3)
    cflmax  U2 at N = 512 with dt = f * DT512 (f up to 2.4) and the CFL check relaxed
    decomp  U2 at N = 512 along x, y, z with 1, 2, 4 ranks
"""
import sys
import os
import math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from isoth_case import write_case

# Riemann problems: (rho_L, rho_R, u_L, u_R, t_end). The R state is on (L/4, 3L/4);
# t_end < the time the fans of the two interfaces meet (riemann_exact.fan_edges):
# U1 1.20, U2 0.92, U4 0.60, Sod 0.90.
RIEMANN = {'U1': (1.0, 1.0, 1.0, -1.0, 0.8), 'U2': (1.0, 1.0, 2.0, -2.0, 0.8),
           'U4': (1.0, 1.0, 4.0, -4.0, 0.4), 'Sod': (10.0, 1.0, 0.0, 0.0, 0.6)}
# min(dt_adv, dt_visc) over the run at N = 512 (proto1d.py, defaults: C = 2, SVV eps0 = 0.5)
DT512 = {'U1': 1.030e-3, 'U2': 5.665e-4, 'U4': 2.982e-4, 'Sod': 1.124e-3}
SAFETY = 0.7


def transverse(nproc):
    """(nly, nlz) for a 1D run along x: the real-space split (dim_decomp = 2)
    is z over dims(1) and x over dims(2); nlz must be divisible by dims(1)."""
    return {1: (2, 2), 2: (2, 2), 4: (2, 2), 8: (2, 4)}[nproc]


def box_1d(n, dir='x', nproc=1, nt=None):
    nt1, nt2 = transverse(nproc) if nt is None else (nt, nt)
    if dir == 'x':
        return dict(nlx=n, nly=nt1, nlz=nt2)
    if dir == 'y':
        return dict(nlx=nt2, nly=n, nlz=nt1)
    return dict(nlx=nt2, nly=nt1, nlz=n)


def riemann_ic(dir, rho_L, rho_R, uL, uR, width=1.0):
    idx = 'xyz'.index(dir)
    u_L = [0.0, 0.0, 0.0]; u_L[idx] = float(uL)
    u_R = [0.0, 0.0, 0.0]; u_R[idx] = float(uR)
    return ('riemann', dict(dir="'%s'" % dir, rho_L=float(rho_L), rho_R=float(rho_R),
                            u_L=u_L, u_R=u_R, b_L=[0.0] * 3, b_R=[0.0] * 3, width=float(width)))


def nstep_for(t_end, dt):
    n = int(round(t_end / dt))
    assert abs(n * dt - t_end) < 1e-9 * t_end, (t_end, dt)
    return n


def fixed_dt(t_end, dt_max):
    """dt <= dt_max with t_end = nstep*dt exactly."""
    n = int(math.ceil(t_end / dt_max * (1.0 - 1e-12)))
    return dict(dt=t_end / n, nstep=n)


def riemann_case(path, name, N, dir='x', nproc=1, t_end=None, dt_factor=1.0, nwrite=20, nt=None, **kw):
    rho_L, rho_R, uL, uR, te = RIEMANN[name]
    te = te if t_end is None else t_end
    time = fixed_dt(te, SAFETY * dt_factor * DT512[name] * 512.0 / N)
    time.update(kw.pop('time', {}))
    diag = dict(write_intvl_1D=te / nwrite)
    diag.update(kw.pop('isoth_diag', {}))
    return write_case(path, nproc=nproc, box=dict(box_1d(N, dir, nproc, nt), **kw.pop('box', {})),
                      time=time, isoth_diag=diag, ic=riemann_ic(dir, rho_L, rho_R, uL, uR), **kw)


def sound_ic(mode=1, amp=1e-3):
    return ('sound_wave', dict(dir="'x'", mode=mode, amp=amp, rho0=1.0))


def batch_smoke(root):
    cases = []
    # sound wave, one period (omega = cs k = 1)
    cases.append(write_case(os.path.join(root, 'sound64'), nproc=1, box=box_1d(64),
                            time=dict(dt=2.0 * math.pi / 2000, nstep=2000),
                            isoth_diag=dict(write_intvl_1D=0.5 * math.pi), ic=sound_ic()))
    # NaN injection at step 5: must abort with a non-zero exit code
    cases.append(write_case(os.path.join(root, 'nan'), nproc=2, box=box_1d(64, nproc=2),
                            time=dict(dt=1e-3, nstep=20), num=dict(debug_nan_step=5), ic=sound_ic()))
    # colliding flow U = 2: defaults, and the configuration of the first smoke run (C = 1, no SVV,
    # filtered ln rho); both are cross-checked with proto1d.py step by step
    cases.append(riemann_case(os.path.join(root, 'collide_U2_N512'), 'U2', 512, t_end=0.1,
                              time=dict(dt=1e-4, nstep=1000), nwrite=2))
    cases.append(riemann_case(os.path.join(root, 'collide_U2_N512_C1_nosvv_fs'), 'U2', 512, t_end=0.1,
                              time=dict(dt=1e-4, nstep=1000), nwrite=2, shock=dict(C_shock=1.0),
                              svv=dict(svv_eps0=0.0), num=dict(filter_lnrho=True)))
    # shock-viscosity self test across decompositions
    for np_ in (1, 2, 4):
        cases.append(write_case(os.path.join(root, 'selftest_np%d' % np_), nproc=np_,
                                box=dict(nlx=32, nly=16, nlz=16), time=dict(dt=1e-4, nstep=1),
                                shock=dict(C_shock=1.0, shock_selftest=True, shock_nmax=2, shock_nsmooth=2),
                                ic=sound_ic()))
    return cases


def batch_smooth(root):
    cases = []
    # sound wave over one period: C = 0 (inviscid; SVV does not act on k = 1 < m_N) and the default C = 2
    for N in (32, 64):
        for tag, shock in (('C0', dict(C_shock=0.0)), ('C2', {})):
            cases.append(write_case(os.path.join(root, 'sound_N%d_%s' % (N, tag)), nproc=1, box=box_1d(N),
                                    time=dict(dt=2.0 * math.pi / 4000, nstep=4000), shock=shock,
                                    isoth_diag=dict(write_intvl_1D=0.25 * math.pi), ic=sound_ic()))
    # simple wave amp = 0.5, k = 1: t_b = 1/(cs amp k) = 2, run to t = 1 = t_b/2;
    # dt = 2.5e-4 for every N so that the time error (~dt^3) stays below the spatial error
    for N in (32, 64, 128, 256, 512):
        for tag, shock, svv in (('C0_nosvv', dict(C_shock=0.0), dict(svv_eps0=0.0)),
                                ('C0_svv', dict(C_shock=0.0), {}), ('C2_svv', {}, {})):
            cases.append(write_case(os.path.join(root, 'simple_N%d_%s' % (N, tag)), nproc=1, box=box_1d(N),
                                    time=dict(dt=2.5e-4, nstep=4000), shock=shock, svv=svv,
                                    isoth_diag=dict(write_intvl_1D=0.25),
                                    ic=('simple_wave', dict(dir="'x'", mode=1, amp=0.5, rho0=1.0))))
    return cases


def batch_shock(root):
    cases = []
    for name in ('U1', 'U2', 'U4', 'Sod'):
        for N in (256, 512, 1024, 2048, 4096):
            cases.append(riemann_case(os.path.join(root, '%s_N%d' % (name, N)), name, N))
    return cases


def batch_scan(root):
    cases = []
    N = 1024
    variants = [
        ('C1', dict(shock=dict(C_shock=1.0))),
        ('C4', dict(shock=dict(C_shock=4.0), dt_factor=0.5)),
        ('C0', dict(shock=dict(C_shock=0.0))),
        ('nosvv', dict(svv=dict(svv_eps0=0.0))),
        ('svvQ0', dict(svv=dict(svv_above_kc=False))),
        ('eps0.25', dict(svv=dict(svv_eps0=0.25))),
        ('eps1', dict(svv=dict(svv_eps0=1.0), dt_factor=0.5)),
        ('fs', dict(num=dict(filter_lnrho=True))),
        ('fs_mc', dict(num=dict(filter_lnrho=True, mass_correction=True))),
        ('mc', dict(num=dict(mass_correction=True))),
        # 3/2 rule: pruned = T truncates every field (filter_lnrho = T); transverse nl = 6 (nlc = 6)
        ('pruned', dict(num=dict(filter_lnrho=True), box=dict(pruned=True), nt=6)),
        ('fs_nt6', dict(num=dict(filter_lnrho=True), nt=6)),
    ]
    for name in ('U2', 'Sod', 'U4'):
        for tag, kw in variants:
            if name == 'U4' and tag not in ('C1', 'C4', 'nosvv', 'svvQ0'):
                continue
            cases.append(riemann_case(os.path.join(root, '%s_N%d_%s' % (name, N, tag)), name, N, **kw))
    return cases


def batch_dtconv(root):
    cases = []
    for f in (1, 2, 4, 8):
        cases.append(riemann_case(os.path.join(root, 'U2_N512_dt%d' % f), 'U2', 512, t_end=0.2,
                                  dt_factor=1.0 / f, nwrite=4))
    return cases


def batch_cflmax(root):
    # the CFL check is relaxed (cfl = 2, cfl_visc = 4) so that dt = f*DT512 runs until it fails
    cases = []
    for f in (1.0, 1.2, 1.6, 2.0, 2.4, 3.0):
        cases.append(riemann_case(os.path.join(root, 'U2_N512_f%.1f' % f), 'U2', 512, t_end=0.4,
                                  time=dict(fixed_dt(0.4, f * DT512['U2']), cfl=2.0), num=dict(cfl_visc=4.0),
                                  nwrite=8))
    return cases


def batch_decomp(root):
    cases = []
    for d in 'xyz':
        for np_ in (1, 2, 4):
            cases.append(riemann_case(os.path.join(root, 'U2_N512_%s_np%d' % (d, np_)), 'U2', 512, dir=d,
                                      nproc=np_, t_end=0.2, nwrite=4))
    return cases


BATCHES = {'smoke': batch_smoke, 'smooth': batch_smooth, 'shock': batch_shock, 'scan': batch_scan,
           'dtconv': batch_dtconv, 'cflmax': batch_cflmax, 'decomp': batch_decomp}

def cost(case):
    """nstep * (grid points) * log2(grid points), for listing the cases by decreasing cost."""
    from analyze_stage1 import read_nml
    d = read_nml(os.path.join(case, 'calliope.in'))
    n = int(d['nlx']) * int(d['nly']) * int(d['nlz'])
    return int(d['nstep']) * n * math.log2(n)


if __name__ == '__main__':
    batch = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.join('runs', 'stage1')
    print(' '.join(sorted(BATCHES[batch](os.path.join(base, batch)), key=cost, reverse=True)))
