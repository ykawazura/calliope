#!/usr/bin/env python3
"""Stage 4 (2D / 3D isothermal MHD) case generator for the Fortran code.

    python3 gen_stage4.py <batch> [base]   writes <base>/<batch>/<case>/{calliope.in, nproc, ic/, case.json, ...}
                                           (default base: runs/stage4) and prints the CASES string for
                                           job.pbs-ps-fes-dev, by decreasing cost

The initial conditions are written by gen_ic.py in the restart format (ic/, copied to restart/ by the job
script) and read with init_type = 'restart'. 2D runs are in the x-z plane (nly = 2, dim_decomp = 1).

Batches
    readback  smooth asymmetric fields on 16^2 (x-z) with 1, 4 ranks and 16^3 with 1, 8 ranks;
              nstep = 0 (restart/ written back must equal ic/) and nstep = 1 (out2d planes at t = 0
              against the fields of gen_ic.py; conservation.dat columns 19-20 for div B)
    regress   Stage 1 and 3 cases run by the Stage 3 binary (bin/calliope_isoth_s3h, built from 60d78cc)
              and the Stage 4 binary (bin/calliope_isoth_s4) with FFTW_ESTIMATE: the outputs must agree
              bit for bit (conservation.dat in its columns 1-18)
    calib     100 steps of the OT vortex (256^2, 1024^2, 2048^2) and of the turbulence (64^3-256^3) for the
              cost per step and the dt limits at the start (S4.0c)
    alfven    CP Alfven waves along k = (1,0,1), (1,0,2) (x-z) and (1,1,2) with its permutations (3D), ideal
              against default settings, three dt and several N and ranks (S4.1)
    planar    planar problems of Stage 3 along an oblique normal, N = 256 (S4.2: ctl, eqv, phy; see
              batch_planar)
    planar2   the phy cases of planar on N = 512 with dt from the N = 256 runs and the compressed 1D runs with
              the bandwidth of the 2D run along k = (1, 2) on N = 256 (written next to them)
    planar3   the phy cases of planar on N = 1024 (as planar2)
    ot        the isothermal OT vortex (x-z) on N = 128, 256 with C_eta = 0, 2 and the MPI decomposition test
              on N = 256 (S4.3; see batch_ot)
    ot2, ot3  the OT vortex on N = 512 (C_eta = 0, 2) and 1024 (C_eta = 0) with dt from the runs on N/2
              (written next to them)
    turb      decaying turbulence (M_s = 1, 4) on 64^3 with dt from the calib runs, and the MPI decomposition
              test on 64^3 (S4.4; see batch_turb)
    turb2     the turbulence on 128^3 with dt from the runs on 64^3 (written next to them); turb2a, turb4a:
              turb1_N128 and turb4_N64 again with an adaptive dt; turb4b: M_s = 4 with rms|dB| = B0 on 64^3;
              turb4b2: turb4b on 128^3 (after turb4b_N64)
    scale2d   strong and weak scaling of the OT vortex (x-z): 1024^2 on 16-1024 ranks and 256^2-2048^2 with
              4096 points per rank, each over SCALE_NSTEP steps (S4.5; see batch_scale2d)
    scale3d   the same for the turbulence (M_s = 4): 256^3 on 64-1024 ranks and 64^3-256^3 with 128^3/16
              points per rank (written next to the scale2d cases)
"""
import sys
import os
import json
import math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from isoth_case import write_case
from gen_stage1 import box_1d
import gen_ic
import stage3 as s3

EXE_OLD = 'bin/calliope_isoth_s3h'
EXE_NEW = 'bin/calliope_isoth_s4'


def box_2d(N, lx=1.0, lz=1.0):
    return dict(nlx=N, nly=2, nlz=N, lx=lx, ly=1.0, lz=lz)


def box_3d(N, l=1.0):
    return dict(nlx=N, nly=N, nlz=N, lx=l, ly=l, lz=l)


def mpi(box):
    """x-z plane: 1D decomposition (the real-space pencils are split in x only)."""
    return dict(dim_decomp=1) if box['nly'] == 2 else dict(dim_decomp=2)


def ic_case(path, g, state, nproc, box, header='', info=None, estimate=True, **kw):
    """Case with the initial state (lnrho_k, m_k, b_k) in ic/ and init_type = 'restart'."""
    write_case(path, nproc=nproc, header=header, mpi=mpi(box), box=box, ic=('restart', {}), **kw)
    gen_ic.write_restart(os.path.join(path, 'ic'), g, *state)
    if estimate:
        open(os.path.join(path, 'fftw_estimate'), 'w').close()
    with open(os.path.join(path, 'case.json'), 'w') as f:
        json.dump(dict(info or {}, nl=g.nl, nproc=nproc), f, indent=1)
    return path


def smooth_fields(g):
    """Asymmetric smooth lnrho, u, B (every mode inside the 2/3 mask for N >= 16). In the x-z plane
    (nly = 2) Y is set to 0: y = pi of the second grid point would give a y-Nyquist mode, which the mask
    removes from m and B but not from ln rho (filter_lnrho = F)."""
    X, Z = g.X, g.Z
    Y = g.Y if g.nly > 2 else 0*g.Y
    lnrho = 0.1*np.sin(X + 2*Z) + 0.05*np.cos(2*X - Z) + 0.03*np.cos(X + 2*Y - 3*Z)
    u = (0.2*np.sin(Z) + 0.1*np.cos(X + Z) + 0.05*np.sin(Y),
         0.1*np.sin(X - Z) + 0.07*np.cos(2*Y - X),
         0.15*np.cos(2*X + Z) + 0.05 + 0.04*np.sin(3*Y + Z))
    b = (0.3 + 0.2*np.cos(X + 2*Z) + 0.05*np.sin(Y - Z),
         0.1*np.sin(2*X - Z) + 0.08*np.cos(Y + 2*X),
         1.0 + 0.1*np.sin(X + 3*Z) + 0.06*np.cos(2*Y))
    return lnrho, u, b


CS2 = 1.0   # c_s^2 = cs2va2 = beta0/2 of the Stage 4 runs (beta0 = 2: c_s = v_A0 = 1)


def cfl_dt(g, lnrho, u, b, cfl=0.5, cs2=CS2):
    """dt_adv of advance.F90 for real fields: cfl/sum_d max(|u_d| + c_f,d)/dl_d over the active directions
    (nl > 2), c_f,d the fast speed along d."""
    rho = np.exp(lnrho)
    va2 = (b[0]**2 + b[1]**2 + b[2]**2)/rho
    rate = 0.0
    for d, (n, L) in enumerate(((g.nlx, g.Lx), (g.nly, g.Ly), (g.nlz, g.Lz))):
        if n > 2:
            cf2 = 0.5*(cs2 + va2 + np.sqrt(np.maximum(0.0, (cs2 + va2)**2 - 4.0*cs2*b[d]**2/rho)))
            rate += (np.abs(u[d]) + np.sqrt(cf2)).max()/(L/n)
    return cfl/rate


def cp_basis(kvec):
    """(n, e1, e2) of a CP Alfven wave with the wave vector kvec: n = k/|k|; e1 = y for k_y = 0 (x-z plane),
    else n x z normalized; e2 = n x e1."""
    k = np.asarray(kvec, float)
    n = k/np.linalg.norm(k)
    e1 = np.array([0.0, 1.0, 0.0]) if k[1] == 0 else np.cross(n, [0.0, 0.0, 1.0])
    e1 /= np.linalg.norm(e1)
    return n, e1, np.cross(n, e1)


def cp_fields(g, kvec, amp, basis, t=0.0, B0=1.0):
    """Circularly polarized Alfven wave along B0 n, an exact solution for any amplitude (|B| and rho
    uniform, u.grad = 0): rho = 1, B = B0 n + dB, u = -dB, dB = amp B0 (e1 cos phi + e2 sin phi),
    phi = k.x - |k| B0 t (propagation along +n at v_A = B0)."""
    n, e1, e2 = basis
    phi = kvec[0]*g.X + kvec[1]*g.Y + kvec[2]*g.Z - np.linalg.norm(kvec)*B0*t
    dB = [amp*B0*(e1[d]*np.cos(phi) + e2[d]*np.sin(phi)) for d in range(3)]
    return g.zeros(), tuple(-x for x in dB), tuple(B0*n[d] + dB[d] for d in range(3))


OT_B0 = math.sqrt(3.0/5.0)


def ot_fields(g, B0=OT_B0):
    """Isothermal Orszag-Tang vortex of Mignone (2007, J. Comput. Phys. 225, 1427, sec. 5.5) scaled to c_s = 1
    (his a = 2: u and B divided by a, t multiplied by a, so that his t = 1.5, 3 are t = 3, 6), in the x-z
    plane (his y is z): rho = 1, u = (-sin z, 0, sin x), B = B0 (-sin z, 0, sin 2x), B0 = sqrt(3/5)
    (beta = 2 c_s^2/<B^2> = 10/3)."""
    X, Z = g.X, g.Z
    return (g.zeros(), (-np.sin(Z), 0*X, np.sin(X)), (-B0*np.sin(Z), 0*X, B0*np.sin(2*X)))


TURB_SEED = 20261103
TURB_GRID = 16   # the random modes (|k| <= 4) are drawn on this grid, so that they do not depend on N
# turbulence problems: name -> (M_s = rms|u|, rms|dB|); rms|dB| = M_s is E_M = E_K. turb4b keeps the dB of turb1
# (rms|dB| = B0): with E_M = E_K at M_s = 4 the mean p_B is 8.5 against p_th = 1, and turb4_N64_a stopped at
# t = 0.13 on a hole of rho in one cell (REPORT_stage4 8.3). turb2b and turb3b lie between turb1 and turb4b
# (same seed and dB, u scaled) for the upper limit of M_s on 64^3 (gen_stage5.py ms)
TURB_PROBS = {'turb1': (1.0, 1.0), 'turb4': (4.0, 4.0), 'turb4b': (4.0, 1.0), 'turb2b': (2.0, 1.0), 'turb3b': (3.0, 1.0)}


def turb_label(prob):
    Ms, dB = TURB_PROBS[prob]
    return 'M_s = %g' % Ms + ('' if dB == Ms else ', rms dB = %g B0' % dB)


def solenoidal_k(g0, rng, kp=2.0, kmin=1.0, kmax=4.0):
    """Fourier coefficients on g0 of a real solenoidal random vector field with E(k) ~ k^4 exp(-2 k^2/kp^2)
    on kmin <= |k| <= kmax and unit rms of |v|: white noise of rng in real space, shaped by
    |v_k| ~ k exp(-k^2/kp^2) and projected on k.v_k = 0."""
    k2 = g0.kx**2 + g0.ky**2 + g0.kz**2
    kk = np.sqrt(k2)
    amp = np.where((kk >= kmin) & (kk <= kmax), kk*np.exp(-k2/kp**2), 0.0)
    vk = g0.project_divfree(tuple(g0.fft(rng.standard_normal(g0.zeros().shape))*amp for _ in range(3)))
    rms = np.sqrt(np.mean(sum(g0.ifft(c)**2 for c in vk)))
    return tuple(c/rms for c in vk)


def embed(g0, f0, g):
    """Fourier coefficients f0 of grid g0 at the same wavenumbers on the grid g (same box, g at least as
    fine as g0; f0 must vanish at the Nyquist wavenumbers of g0)."""
    f = np.zeros((g.nkx, g.nkz, g.nky), complex)
    f[np.ix_(g0.ikx % g.nkx, g0.ikz % g.nkz, g0.iky)] = f0
    return f


def turb_fields(g, Ms, dB=None, seed=TURB_SEED, B0=1.0):
    """Decaying turbulence: rho = 1, u and dB solenoidal random fields (solenoidal_k, u first, then dB, from
    default_rng(seed)) with rms|u| = Ms and rms|dB| = dB (default Ms: E_M = E_K; c_s = 1), and B = B0 z + dB. The modes are
    drawn on TURB_GRID^3, so that every N sees the same fields. Returns the fields and a dict with
    sigma_c = 2<u.dB>/<u^2 + dB^2> and H_m = <A.dB> (A = i k x dB_k/k^2)."""
    g0 = gen_ic.Grid(TURB_GRID, TURB_GRID, TURB_GRID, g.Lx/(2*np.pi), g.Ly/(2*np.pi), g.Lz/(2*np.pi))
    rng = np.random.default_rng(seed)
    uk = tuple(Ms*c for c in solenoidal_k(g0, rng))
    dB = Ms if dB is None else dB
    bk = tuple(dB*c for c in solenoidal_k(g0, rng))
    u0, b0 = [g0.ifft(c) for c in uk], [g0.ifft(c) for c in bk]
    k2 = g0.kx**2 + g0.ky**2 + g0.kz**2
    k2inv = np.where(k2 > 0, 1/np.where(k2 > 0, k2, 1), 0.0)
    ak = (1j*(g0.ky*bk[2] - g0.kz*bk[1])*k2inv, 1j*(g0.kz*bk[0] - g0.kx*bk[2])*k2inv,
          1j*(g0.kx*bk[1] - g0.ky*bk[0])*k2inv)
    info = dict(Ms=Ms, dB=dB, seed=seed, B0=B0,
                sigma_c=float(2*np.mean(sum(a*b for a, b in zip(u0, b0)))/np.mean(sum(a**2 + b**2 for a, b in zip(u0, b0)))),
                H_m=float(np.mean(sum(g0.ifft(a)*b for a, b in zip(ak, b0)))))
    u = tuple(g.ifft(embed(g0, c, g)) for c in uk)
    b = tuple(g.ifft(embed(g0, c, g)) for c in bk)
    return (g.zeros(), u, (b[0], b[1], b[2] + B0)), info


def problem_fields(prob, g):
    """(fields, info) of 'ot' and of the turbulence problems of TURB_PROBS."""
    if prob == 'ot':
        return ot_fields(g), dict(problem='ot', B0=OT_B0)
    return turb_fields(g, *TURB_PROBS[prob])


CALIB_PROBS = (('ot', 256, 16), ('ot', 1024, 256), ('ot', 2048, 256),
               ('turb4', 64, 8), ('turb4', 128, 64), ('turb4', 256, 256), ('turb1', 64, 8))


def batch_calib(root, probs=CALIB_PROBS, header='Stage 4 cost calibration (gen_stage4.py calib)'):
    """100 steps of the OT vortex (x-z) and of the turbulence for the cost per step (timers of
    calliope.out.std) and dt_adv, dt_visc at the start; dt = 0.5 dt_adv(t = 0) of the initial fields,
    reset_method = 'multiply', FFTW_MEASURE (as the production runs). Run one case at a time (no PAR).
    probs: (problem, N, ranks) of the cases."""
    cases = []
    for prob, N, P in probs:
        box = box_2d(N) if prob == 'ot' else box_3d(N)
        g = gen_ic.Grid(box['nlx'], box['nly'], box['nlz'])
        fields, info = problem_fields(prob, g)
        dt = float('%.1e' % (0.5*cfl_dt(g, *fields)))
        cases.append(ic_case(os.path.join(root, '%s_N%d_np%d' % (prob, N, P)), g,
                             gen_ic.state_from_real(g, *fields), P, box, estimate=False,
                             header=header, info=dict(info, batch='calib', dt=dt),
                             time=dict(dt=dt, nstep=100, reset_method="'multiply'"),
                             isoth_diag=dict(write_intvl_cons=0.0, budget=True)))
    return cases


ALFVEN_NPER = 10
IDEAL = dict(shock=dict(C_shock=0.0, C_eta=0.0), svv=dict(svv_eps0=0.0))
DEFAULT = dict(shock=dict(C_shock=2.0, C_eta=2.0), svv=dict(svv_eps0=0.5))
# the dissipation terms of DEFAULT one at a time (cpdiag)
SETTINGS = dict(ideal=IDEAL, default=DEFAULT,
                noeta=dict(shock=dict(C_shock=2.0, C_eta=0.0), svv=dict(svv_eps0=0.5)),
                nosvv=dict(shock=dict(C_shock=2.0, C_eta=2.0), svv=dict(svv_eps0=0.0)),
                etaonly=dict(shock=dict(C_shock=0.0, C_eta=2.0), svv=dict(svv_eps0=0.0)),
                svvonly=dict(shock=dict(C_shock=0.0, C_eta=0.0), svv=dict(svv_eps0=0.5)))


def roll_basis(basis, s):
    """basis of the problem rotated by the cyclic permutation of the axes x -> y -> z (s times)."""
    return tuple(np.roll(v, s) for v in basis)


def alfven_case(root, dim, kvec, amp, N, setting, n_p, nproc=1, perm=0):
    """CP Alfven wave over ALFVEN_NPER periods T = 2 pi/|k| (v_A = B0 = 1) with dt = T/n_p, out3d every
    period; perm > 0: the problem of kvec0 = roll(kvec, -perm) with the basis rolled by perm."""
    box = box_2d(N) if dim == '2d' else box_3d(N)
    g = gen_ic.Grid(box['nlx'], box['nly'], box['nlz'])
    basis = roll_basis(cp_basis(np.roll(kvec, -perm)), perm)
    T = 2*np.pi/np.linalg.norm(kvec)
    name = '%s_k%s_a%s_N%d_%s_p%d' % (dim, ''.join('%d' % k for k in kvec), {0.1: '01', 1.0: '1'}[amp], N,
                                      setting, n_p)
    if nproc > 1:
        name += '_np%d' % nproc
    intvl = T*(1 + 1e-9)
    return ic_case(os.path.join(root, name), g, gen_ic.state_from_real(g, *cp_fields(g, kvec, amp, basis)),
                   nproc, box, header='Stage 4 CP Alfven wave (gen_stage4.py alfven)',
                   info=dict(batch='alfven', dim=dim, kvec=[int(k) for k in kvec], amp=amp,
                             basis=[[float(x) for x in v] for v in basis],
                             setting=setting, n_p=n_p, T=T, nper=ALFVEN_NPER, perm=perm),
                   time=dict(dt=T/n_p, nstep=ALFVEN_NPER*n_p),
                   diag=dict(write_intvl_3D=intvl),
                   isoth_diag=dict(write_intvl_cons=T/20*(1 + 1e-9), budget=True),
                   **SETTINGS[setting])


def alfven_np(dim, kvec, amp, Ns):
    """Smallest n_p = 100 2^m with T/n_p <= 0.8 dt_adv (cfl = 0.5) of the wave on the finest grid of Ns."""
    N = max(Ns)
    g = gen_ic.Grid(N, 2 if dim == '2d' else N, N)
    dt = 0.8*cfl_dt(g, *cp_fields(g, kvec, amp, cp_basis(kvec)))
    n_p = 100
    while 2*np.pi/np.linalg.norm(kvec)/n_p > dt:
        n_p *= 2
    return n_p


def batch_alfven(root):
    """CP Alfven waves (S4.1). For each (k, amplitude): the ideal (C_shock = C_eta = 0, no SVV) and the
    default (C_shock = 2, C_eta = 2, ES-SVV) run on every N with a common dt = T/n_p (spatially exact
    single mode: the error must not depend on N), dt/2 and dt/4 on the coarsest N (error ~ dt^3);
    3D: the cyclic permutations of k = (1,1,2); MPI: 2D N = 32 with 4 and 16 ranks, 3D N = 16 with 8."""
    cases = []
    groups = [('2d', (1, 0, 1), (16, 32, 64)), ('2d', (1, 0, 2), (16, 32, 64)), ('3d', (1, 1, 2), (16, 32))]
    for dim, kvec, Ns in groups:
        for amp in (0.1, 1.0):
            n_p = alfven_np(dim, kvec, amp, Ns)
            for setting in ('ideal', 'default'):
                for N in Ns:
                    cases.append(alfven_case(root, dim, kvec, amp, N, setting, n_p, nproc=8 if N == 32 and dim == '3d' else 1))
                for m in (2, 4):
                    cases.append(alfven_case(root, dim, kvec, amp, Ns[0], setting, m*n_p))
                if dim == '3d' and amp == 1.0:
                    for perm in (1, 2):
                        cases.append(alfven_case(root, dim, tuple(int(k) for k in np.roll(kvec, perm)), amp, 16, setting, n_p,
                                                 perm=perm))
            if amp == 1.0:
                if dim == '2d' and kvec == (1, 0, 1):
                    for P in (4, 16):
                        cases.append(alfven_case(root, dim, kvec, amp, 32, 'default', n_p, nproc=P))
                if dim == '3d':
                    cases.append(alfven_case(root, dim, kvec, amp, 16, 'default', n_p, nproc=8))
    return cases


def batch_cpdiag(root):
    """The A = 1 CP waves of S4.1 that grow faster with the default settings than the bound of alfven, with
    the dissipation terms of DEFAULT one at a time (SETTINGS): 2D N = 16 and 3D N = 16, 32 with the n_p of
    batch_alfven, so that the default and ideal runs of runs/stage4/alfven are the references."""
    cases = []
    for dim, kvec, Ns, runs in (('2d', (1, 0, 1), (16, 32, 64), ((16, 1),)), ('2d', (1, 0, 2), (16, 32, 64), ((16, 1),)),
                                ('3d', (1, 1, 2), (16, 32), ((16, 1), (32, 8)))):
        n_p = alfven_np(dim, kvec, 1.0, Ns)
        for N, P in runs:
            for setting in ('noeta', 'nosvv', 'etaonly', 'svvonly'):
                cases.append(alfven_case(root, dim, kvec, 1.0, N, setting, n_p, nproc=P))
    return cases


def batch_readback(root):
    cases = []
    for dim, N, nps in (('2d', 16, (1, 4)), ('3d', 16, (1, 8))):
        box = box_2d(N) if dim == '2d' else box_3d(N)
        g = gen_ic.Grid(box['nlx'], box['nly'], box['nlz'])
        state = gen_ic.state_from_real(g, *smooth_fields(g))
        for np_ in nps:
            for nstep in (0, 1):
                cases.append(ic_case(os.path.join(root, '%s_N%d_np%d_n%d' % (dim, N, np_, nstep)), g, state, np_,
                                     box, header='Stage 4 read-back test (gen_stage4.py readback)',
                                     info=dict(batch='readback', dim=dim),
                                     time=dict(dt=1e-3, nstep=nstep),
                                     isoth_diag=dict(write_intvl_1D=1.0, write_intvl_planes=1.0,
                                                     write_intvl_cons=0.0, budget=True),
                                     shock=dict(C_eta=2.0)))
    return cases


def batch_regress(root):
    """The same input run by EXE_OLD and EXE_NEW (file "exe"), FFTW_ESTIMATE."""
    import gen_stage1 as g1
    import gen_stage3 as g3
    tmp = []
    src = os.path.join(root, 'src')
    # Stage 1: sound wave, colliding flow U2 (shock viscosity + SVV + budget), shock self test with 1, 2, 4 ranks
    tmp.append(write_case(os.path.join(src, 'sound64'), nproc=1, box=g1.box_1d(64),
                          time=dict(dt=2.0*math.pi/2000, nstep=2000),
                          isoth_diag=dict(write_intvl_1D=0.5*math.pi), ic=g1.sound_ic()))
    tmp.append(g1.riemann_case(os.path.join(src, 'collide_U2_N512'), 'U2', 512, t_end=0.1,
                               time=dict(dt=1e-4, nstep=1000), nwrite=2))
    for np_ in (1, 2, 4):
        tmp.append(write_case(os.path.join(src, 'selftest_np%d' % np_), nproc=np_,
                              box=dict(nlx=32, nly=16, nlz=16), time=dict(dt=1e-4, nstep=1),
                              shock=dict(C_shock=1.0, shock_selftest=True, shock_nmax=2, shock_nsmooth=2),
                              ic=g1.sound_ic()))
    # Stage 3: MHD with the shock resistivity, along z with 4 ranks
    tmp.append(g3.case_dir(src, g3.find(g3.eta_name('C3_N256')), dir='z', nproc=4, tag='_z_np4'))
    cases = []
    for c in tmp:
        for tag, exe in (('old', EXE_OLD), ('new', EXE_NEW)):
            dst = os.path.join(root, os.path.basename(c) + '_' + tag)
            os.makedirs(dst, exist_ok=True)
            for f in os.listdir(c):
                with open(os.path.join(c, f), 'rb') as fi, open(os.path.join(dst, f), 'wb') as fo:
                    fo.write(fi.read())
            with open(os.path.join(dst, 'exe'), 'w') as f:
                f.write(exe + '\n')
            open(os.path.join(dst, 'fftw_estimate'), 'w').close()
            cases.append(dst)
    return cases


# ------------------------------------------------------------------ S4.2 planar problems along an oblique normal
PLANAR_PROBS = ('F1', 'F1_eta2', 'C2', 'C2_eta2', 'TDS_w1', 'RDh_w4', 'RDpc')
PLANAR_K = ((1, 1), (1, 2))
PLANAR_CTL = ('F1', 'C2_eta2', 'TDS_w1', 'RDh_w4')
PLANAR_EQV = ('F1', 'F1_eta2', 'C2_eta2', 'TDS_w1', 'RDh_w4')
PLANAR_EQV_NSTEP, PLANAR_EQV_NLINE = 200, 10
# nstep of the first runs (N = 256) = PLANAR_MARGIN nstep of the Stage 3 run at the same N: the rate of the shock
# viscosity of the initial profile is about 2 |k| times the 1D one in the 2D run (kd2_max of two directions, gradients
# |k| times steeper on the compressed profile) and the run lasts t_end/|k|; the compressed 1D run has |k|^2 kd2_max
# of the 1D one (spacing L/(|k| N)), so its rate is |k|^3 times and its nstep PLANAR_MARGIN max(1, |k|^2/2) times
PLANAR_MARGIN = 2.4
PLANAR_NP = {256: 16, 512: 64, 1024: 256}
PLANAR_PLANES = {256: 4, 512: 2, 1024: 1}


def s3_case(prob, N):
    """The stage3 case of prob ('F1', 'F1_eta2', 'TDS_w1', ...) at N."""
    import gen_stage3 as g3
    base, _, eta = prob.partition('_eta')
    return g3.find('%s_N%d' % (base, N) + ('_eta' + eta if eta else ''))


def s3_nstep(c):
    import gen_stage3 as g3
    return json.load(open(g3.NSTEP_FILE))[c['name']][0]


def planar_frame(kvec):
    """(n, t1, t2) of the planar problem with the normal along kvec = (a, b) in the x-z plane: n = (a, 0, b)/|k|,
    t1 = (b, 0, -a)/|k| (in the plane), t2 = y; n x t1 = t2, as x, y, z of the 1D problem along x."""
    a, b = kvec
    kn = math.hypot(a, b)
    return np.array([a, 0.0, b])/kn, np.array([b, 0.0, -a])/kn, np.array([0.0, 1.0, 0.0])


def planar_index(N, kvec):
    """p(i, l) = mod(a i + b l, N), the point of the 1D profile at the point (i, l) of the x-z grid: the state
    f(i, l) = F(p(i, l)) stays a function of p, since every operator of the scheme commutes with the shifts of
    the grid (i, l) -> (i + 1, l), (i, l + 1), which shift p by a and b."""
    a, b = kvec
    return np.mod(a*np.arange(N)[:, None] + b*np.arange(N)[None, :], N)


def planar_fields(W, bx, N, kvec):
    """lnrho, u, B on the x-z grid N^2 (nly = 2) of the 1D profile W = (rho, u_n, u_t1, u_t2, B_t1, B_t2), B_n = bx,
    given on N points of its period L: f(i, l) = W(p(i, l)) with the components along (n, t1, t2) of planar_frame.
    On the box L^2 this is the 1D problem compressed by 1/|k| (period L/|k| along n)."""
    n, t1, t2 = planar_frame(kvec)
    Wp = np.asarray(W)[:, planar_index(N, kvec)][..., None]*np.ones(2)     # (6, nlx, nlz, nly)
    u = tuple(Wp[1]*n[d] + Wp[2]*t1[d] + Wp[3]*t2[d] for d in range(3))
    b = tuple(bx*n[d] + Wp[4]*t1[d] + Wp[5]*t2[d] for d in range(3))
    return np.log(Wp[0]), u, b


def line_intvl(nstep, dt, nline):
    """write_intvl_1D (or _planes) for nline + 1 records at the steps j nstep/nline (as gen_stage3.case_dir):
    write_due fires in the step whose (tt, tt + dt] holds j intvl = j (nstep/nline) dt + j dt/(2 nstep) and the
    record holds the state at tt."""
    assert nstep % nline == 0, (nstep, nline)
    return nstep*dt/nline*(1.0 + 0.5/nstep)


def planar_2d_case(root, name, c, kvec, nstep, nproc, info, nrun=None, nline=None, planes=0, estimate=False,
                   shock=None, svv=None):
    """2D run (x-z, N^2) of the stage3 case c along the normal kvec over t_end/|k| with dt = (t_end/|k|)/nstep
    (nrun steps: default nstep), reset_method = 'stop', out1d (the line along x at z = 0: p = i for a = 1) and
    out2d planes."""
    N, kn = c['N'], math.hypot(*kvec)
    te = s3.t_end(c)/kn
    dt = te/nstep
    nrun = nstep if nrun is None else nrun
    nline = nline or s3.NWRITE
    box = box_2d(N, s3.LX, s3.LX)
    g = gen_ic.Grid(N, 2, N, s3.LX, 1.0, s3.LX)
    state = gen_ic.state_from_real(g, *planar_fields(s3.initial(c), s3.bx_of(c), N, kvec))
    return ic_case(os.path.join(root, name), g, state, nproc, box, estimate=estimate,
                   header='Stage 4 planar problem %s along k = %s (gen_stage4.py planar)' % (c['name'], kvec),
                   info=dict(info, batch='planar', dim='2d', case=c, kvec=list(kvec), kn=kn, t_end=te, nstep=nstep,
                             nrun=nrun, dt=dt, nline=nline, planes=planes),
                   time=dict(dt=dt, nstep=nrun),
                   shock=dict(dict(C_eta=c['C_eta']), **(shock or {})), svv=svv or {},
                   isoth_diag=dict(write_intvl_1D=line_intvl(nrun, dt, nline),
                                   write_intvl_planes=line_intvl(nrun, dt, planes) if planes else 0.0,
                                   write_intvl_cons=0.0, budget=True))


def planar_1d_case(root, name, c, kvec, nstep, info, nrun=None, nline=None, shock=None, svv=None, estimate=False):
    """1D run along x of the stage3 case c compressed by 1/|k| (box L/|k|, t_end/|k|, the same dt and steps as
    planar_2d_case) with Delta = L/N of the x-z grid (shock_dl): the problem that the 2D run solves on p."""
    import gen_stage3 as g3
    N, kn = c['N'], math.hypot(*kvec)
    te = s3.t_end(c)/kn
    dt = te/nstep
    nrun = nstep if nrun is None else nrun
    nline = nline or s3.NWRITE
    path = os.path.join(root, name)
    write_case(path, nproc=1, header='Stage 4 1D reference of %s compressed by 1/|k|, |k| = %.6f (gen_stage4.py '
                                     'planar)' % (c['name'], kn),
               box=dict(box_1d(N, 'x', 1), lx=s3.LX/kn, ly=s3.LX/kn, lz=s3.LX/kn),
               time=dict(dt=dt, nstep=nrun),
               shock=dict(dict(C_eta=c['C_eta'], shock_dl=s3.L/N), **(shock or {})), svv=svv or {},
               isoth_diag=dict(write_intvl_1D=line_intvl(nrun, dt, nline), write_intvl_cons=0.0, budget=True),
               ic=('profile1d', dict(dir="'x'", table="'profile1d.dat'")))
    W = s3.initial(c)
    u, b = g3.rotate(W, s3.bx_of(c), 'x')
    np.savetxt(os.path.join(path, 'profile1d.dat'), np.column_stack([np.arange(N)*(s3.L/kn/N), W[0], *u, *b]),
               fmt='%.17e', header='%s compressed by 1/|k| = 1/%.6f along x: x, rho, u_x, u_y, u_z, b_x, b_y, b_z'
                                   % (c['name'], kn))
    if estimate:
        open(os.path.join(path, 'fftw_estimate'), 'w').close()
    with open(os.path.join(path, 'case.json'), 'w') as f:
        json.dump(dict(info, batch='planar', dim='1d', case=c, kvec=list(kvec), kn=kn, t_end=te, nstep=nstep,
                       nrun=nrun, dt=dt, nline=nline, nl=[N, 2, 2], nproc=1), f, indent=1)
    return path


def planar_nstep(c, kvec, factor=PLANAR_MARGIN):
    """Multiple of NWRITE above factor x the Stage 3 nstep of c (see PLANAR_MARGIN)."""
    return s3.NWRITE*int(math.ceil(factor*s3_nstep(c)/s3.NWRITE))


def batch_planar(root):
    """S4.2 first job (N = 256).
    ctl  k = (1, 0): the 2D run against the 1D run along x with the same dt (2x the Stage 3 steps: dt_adv and
         dt_visc of the 2D run are not below half of the 1D ones for a state that depends on x only)
    eqv  k = (1, 1), no SVV, 200 steps: the 2D run against the compressed 1D run with shock_nmax = shock_nsmooth = 2
         (on p the z and x passes of get_zeta are the same 3-point operators: max over p +- 2 and
         [1, 4, 6, 4, 1]/16) and shock_dl = L/N; the 2/3 mask (|m| < N/3) and k_d are those of the 1D run
    phy  k = (1, 1), (1, 2), default settings: the 2D run and the compressed 1D run (default stencil) over
         t_end/|k|; their L1 errors against the 1D run of Stage 3 (same N) separate the effect of Delta = L/N
         (|k| times the spacing of the lattice p along n) from that of the 2D stencil, mask and SVV; along
         k = (1, 2) the mask keeps half the modes of the compressed 1D run on N, whose run on N/2 (_1dh, planar2)
         has the bandwidth of the 2D run."""
    cases = []
    N = 256
    for prob in PLANAR_CTL:
        c = s3_case(prob, N)
        ns = 2*s3_nstep(c)
        info = dict(role='ctl', prob=prob, nstep_s3=s3_nstep(c))
        cases.append(planar_2d_case(root, 'ctl_%s_N%d_2d' % (prob, N), c, (1, 0), ns, PLANAR_NP[N], info,
                                    planes=PLANAR_PLANES[N], estimate=True))
        cases.append(planar_1d_case(root, 'ctl_%s_N%d_1d' % (prob, N), c, (1, 0), ns, info, estimate=True))
    for prob in PLANAR_EQV:
        c = s3_case(prob, N)
        ns = planar_nstep(c, (1, 1))
        info = dict(role='eqv', prob=prob, nstep_s3=s3_nstep(c))
        kw = dict(nrun=PLANAR_EQV_NSTEP, nline=PLANAR_EQV_NLINE, svv=dict(svv_eps0=0.0), estimate=True)
        cases.append(planar_2d_case(root, 'eqv_%s_k11_N%d_2d' % (prob, N), c, (1, 1), ns, PLANAR_NP[N], info,
                                    planes=1, **kw))
        cases.append(planar_1d_case(root, 'eqv_%s_k11_N%d_1d' % (prob, N), c, (1, 1), ns, info,
                                    shock=dict(shock_nmax=2, shock_nsmooth=2), **kw))
    cases += batch_planar_phy(root, (N,))
    return cases


def batch_planar_phy(root, Ns, dt_from=None, bw_only=False):
    """phy cases on Ns; dt_from: directory of the N = 256 runs whose min(dt_adv, dt_visc) (cfl.dat) gives dt on N
    (PLANAR_DT_SAFETY (256/N)/1.05^log2(N/256), from the Stage 3 steps: x 2.02-2.08 per doubling of N). For
    max(|a|, |b|) > 1 also the compressed 1D run with the bandwidth of the 2D run (_1dh, planar_bw_case; the dt of
    the compressed 1D run on N, whose dt limits are the lower); bw_only: only these."""
    cases = []
    for prob in PLANAR_PROBS:
        for kvec in PLANAR_K:
            for N in Ns:
                c = s3_case(prob, N)
                k = '%d%d' % kvec
                info = dict(role='phy', prob=prob, nstep_s3=s3_nstep(c))
                ns = {}
                for dim, factor in (('2d', PLANAR_MARGIN), ('1d', PLANAR_MARGIN*max(1.0, (kvec[0]**2 + kvec[1]**2)/2))):
                    if dt_from is None:
                        ns[dim] = planar_nstep(c, kvec, factor)
                    else:
                        ns[dim] = planar_nstep_from(os.path.join(dt_from, 'phy_%s_k%s_N256_%s' % (prob, k, dim)), c,
                                                    kvec)
                if not bw_only:
                    cases.append(planar_2d_case(root, 'phy_%s_k%s_N%d_2d' % (prob, k, N), c, kvec, ns['2d'],
                                                PLANAR_NP[N], info, planes=PLANAR_PLANES[N]))
                    cases.append(planar_1d_case(root, 'phy_%s_k%s_N%d_1d' % (prob, k, N), c, kvec, ns['1d'], info))
                ch = planar_bw_case(c, kvec)
                if ch is not None:
                    cases.append(planar_1d_case(root, 'phy_%s_k%s_N%d_1dh' % (prob, k, N), ch, kvec, ns['1d'], info,
                                                shock=dict(shock_dl=s3.L/N)))
    return cases


def planar_bw_case(c, kvec):
    """The stage3 case c on N/max(|a|, |b|) points with the same profile (width w in cells halved), or None for
    max(|a|, |b|) = 1. Mode m of the profile sits at (a m, b m) on the x-z grid, so the 2/3 mask keeps
    |m| < N/(3 max(|a|, |b|)) (k = (1, 2): N/6, half the bandwidth of the compressed 1D run on N points);
    the compressed 1D run of this case on N/max(|a|, |b|) points with Delta = L/N has the bandwidth and Delta of
    the 2D run."""
    r = max(abs(kvec[0]), abs(kvec[1]))
    if r == 1:
        return None
    assert c['N'] % r == 0, (c, kvec)
    ch = dict(c, N=c['N']//r, name='%s_bw%d' % (c['name'], r))
    if 'w' in c:
        ch['w'] = c['w']/r     # stage3.rd_W, td_W: width w L/N
    return ch


PLANAR_DT_SAFETY = 0.85


def planar_nstep_from(case256, c, kvec):
    """nstep on c['N'] from the N = 256 run case256 (see batch_planar_phy)."""
    status = open(os.path.join(case256, 'exit_status')).read().split()
    assert status == ['exit', '0'], (case256, status)
    cfl = np.loadtxt(os.path.join(case256, 'cfl.dat'), ndmin=2)
    r = c['N']/256
    dt = PLANAR_DT_SAFETY*np.minimum(cfl[:, 2], cfl[:, 3]).min()/r/1.05**math.log2(r)
    return s3.NWRITE*int(math.ceil(s3.t_end(c)/math.hypot(*kvec)/dt/s3.NWRITE))


def batch_planar2(root):
    """S4.2 second job: the phy cases on N = 512 with dt from the N = 256 runs and the _1dh runs on N = 256
    (added after the first job)."""
    return batch_planar_phy(root, (512,), dt_from=root) + batch_planar_phy(root, (256,), bw_only=True)


def batch_planar3(root):
    """S4.2 third job (after planar2): the phy cases on N = 1024 with dt from the N = 256 runs."""
    return batch_planar_phy(root, (1024,), dt_from=root)


# S4.3: t_end = 2 pi (extended from pi of the plan) covers t = 3 and 6, Mignone's t = 1.5 and 3 (see ot_fields)
OT_TEND = 2*math.pi
OT_NP = {128: 4, 256: 16, 512: 128, 1024: 512, 2048: 512}
OT_PLANES = 0.25    # out2d planes (the whole field of the x-z plane) at t = 0, 0.25, ..., 6.25 and t_end
OT_CONS = 0.01      # conservation.dat
OT_PILOT = 4        # N = 128, 256: dt = dt_adv(t = 0)/OT_PILOT (dt_cfl(t) is not known before these runs)
OT_DT_SAFETY = 0.85
OT_DECOMP_P = (4, 16, 64, 256)
OT_DECOMP_T = 1.0   # the decomposition runs (N = 256) stop at about t = 1


def ot_case(root, name, N, nstep, nproc, C_eta, info, nrun=None, planes=OT_PLANES, cons=OT_CONS):
    """OT vortex (x-z, N^2) over OT_TEND with dt = OT_TEND/nstep (nrun steps: default nstep), C_shock = 2,
    ES-SVV and C_eta, reset_method = 'stop', FFTW_MEASURE."""
    box = box_2d(N)
    g = gen_ic.Grid(N, 2, N)
    dt = OT_TEND/nstep
    nrun = nstep if nrun is None else nrun
    return ic_case(os.path.join(root, name), g, gen_ic.state_from_real(g, *ot_fields(g)), nproc, box,
                   estimate=False, header='Stage 4 isothermal OT vortex (gen_stage4.py ot)',
                   info=dict(info, batch='ot', problem='ot', B0=OT_B0, N=N, C_eta=C_eta, t_end=OT_TEND,
                             nstep=nstep, nrun=nrun, dt=dt, planes=planes),
                   time=dict(dt=dt, nstep=nrun),
                   shock=dict(C_shock=2.0, C_eta=C_eta), svv=dict(svv_eps0=0.5),
                   isoth_diag=dict(write_intvl_planes=planes, write_intvl_cons=cons, budget=True))


def ot_pilot_nstep(N):
    g = gen_ic.Grid(N, 2, N)
    return int(math.ceil(OT_PILOT*OT_TEND/cfl_dt(g, *ot_fields(g))))


def min_dt_cfl(case):
    """min over the run (all chunks) of min(dt_adv, dt_visc) of cfl.dat; the run must have ended with exit 0."""
    from analyze_stage1 import load_runs
    status = open(os.path.join(case, 'exit_status')).read().split()
    assert status == ['exit', '0'], (case, status)
    cfl = np.vstack(load_runs(case, 'cfl.dat'))
    return float(np.minimum(cfl[:, 2], cfl[:, 3]).min())


def ot_nstep_from(root, N, Nsrc):
    """nstep on N from the runs on Nsrc (C_eta = 0, 2): OT_DT_SAFETY (Nsrc/N)/1.05^log2(N/Nsrc) of their min
    dt_cfl (dt_adv, dt_visc and cfl_svv/lambda_svv all scale as 1/N for a resolved shock, as for planar)."""
    dt = min(min_dt_cfl(os.path.join(root, 'ot_N%d_eta%d' % (Nsrc, e))) for e in (0, 2))
    r = N/Nsrc
    return int(math.ceil(OT_TEND/(OT_DT_SAFETY*dt/r/1.05**math.log2(r))))


def batch_ot(root):
    """S4.3 first job: the OT vortex on N = 128, 256 (C_eta = 0, 2; dt = dt_adv(0)/OT_PILOT) and the MPI
    decomposition test on N = 256 (C_eta = 2, the dt of ot_N256_eta2): P = 4, 16, 64, 256 over 10 steps
    (planes at steps 0 and 10) and up to t = OT_DECOMP_T (conservation.dat at every step, planes as the
    production runs; P = 256 holds one x plane per rank)."""
    cases = []
    for N in (128, 256):
        ns = ot_pilot_nstep(N)
        for e in (0, 2):
            cases.append(ot_case(root, 'ot_N%d_eta%d' % (N, e), N, ns, OT_NP[N], float(e), dict(role='pilot')))
    N, ns = 256, ot_pilot_nstep(256)
    dt = OT_TEND/ns
    nlong = int(round(OT_DECOMP_T/dt))
    for P in OT_DECOMP_P:
        cases.append(ot_case(root, 'dec_N%d_np%d_n10' % (N, P), N, ns, P, 2.0, dict(role='decomp'), nrun=10,
                             planes=line_intvl(10, dt, 1), cons=0.0))
        cases.append(ot_case(root, 'dec_N%d_np%d' % (N, P), N, ns, P, 2.0, dict(role='decomp'), nrun=nlong,
                             cons=0.0))
    return cases


def batch_ot2(root):
    """S4.3 second job: N = 512 (C_eta = 0, 2) with dt from the N = 256 runs."""
    return [ot_case(root, 'ot_N512_eta%d' % e, 512, ot_nstep_from(root, 512, 256), OT_NP[512], float(e),
                    dict(role='production', dt_from=256)) for e in (0, 2)]


def batch_ot3(root):
    """S4.3 third job: N = 1024 with dt from the N = 512 runs, C_eta = 0 only (the default after S4.2; C_eta = 2
    converges at first order over N = 128-512, and a run on 1024^2 costs about 3 chained jobs on 512 ranks)."""
    return [ot_case(root, 'ot_N1024_eta0', 1024, ot_nstep_from(root, 1024, 512), OT_NP[1024], 0.0,
                    dict(role='production', dt_from=512))]


def batch_ot3b(root):
    """S4.3: ot_N1024_eta0 again, on the 256 ranks of one node and up to t = pi (the t_end of the plan; the peak of
    eps is at t = 1.68, and after t = 4 the runs on N = 256, 512 do not converge pointwise: plasmoids). On two
    nodes (512 ranks) a step cost 1.85 s against 0.136 s on one node (job 481249, (T(110) - T(10))/100), nearly all
    of it in the FFTs: the transposes of the 1D decomposition of the x-z plane are all-to-all among all ranks,
    with messages of about 64 B through the one NIC of a node. So ot_N1024_eta0 was killed; its file
    superseded_by names this case. dt as for ot_N1024_eta0."""
    nstep = ot_nstep_from(root, 1024, 512)
    assert nstep % 2 == 0
    return [ot_case(root, 'ot_N1024_eta0_p256', 1024, nstep, 256, 0.0, dict(role='production', dt_from=512),
                    nrun=nstep//2)]


# S4.4: decaying turbulence (turb_fields) over TURB_NTAU tau0, tau0 = pi/u_rms (u_rms = M_s), C_shock = 2,
# ES-SVV, C_eta = 0 (the default kept after S4.2). Times of the outputs in units of tau0.
TURB_NTAU = 3.0
TURB_NP = {64: 64, 128: 128}
TURB_PLANES = 0.125   # plane_y0, plane_z0 (out2d)
TURB_3D = 0.5         # out3d (Fourier, ~1.9 GB at 256^3 with the record gaps of mpiio_write_var); the end of
                      # every chunk adds one
TURB_CONS = 0.005     # conservation.dat, budget.dat, cfl.dat (0: every step)
TURB_PILOT = 0.5      # N = 64: dt = TURB_PILOT x the min dt_cfl of the calib run (100 steps, runs/stage4/calib),
                      # cfl.dat at every step; the N = 64 runs give dt_cfl(t) for the larger N
TURB_DT_SAFETY = 0.85
TURB_DECOMP_P = (1, 8, 64)
TURB_DECOMP_NSTEP = 1000   # the long decomposition runs (M_s = 4, against the run on 64 ranks)
TRACK_NMAX = 10  # upper bound of nstep in units of t_end/dt0 for turb_case(track=True) up to t_end


def turb_case(root, name, prob, N, nstep, nproc, info, nrun=None, planes=None, cons=TURB_CONS, out3d=TURB_3D,
              adaptive=False, track=False):
    """Decaying turbulence prob (TURB_PROBS) on N^3 over TURB_NTAU tau0 with dt = t_end/nstep (nrun
    steps: default nstep), reset_method = 'stop', FFTW_MEASURE. planes, cons, out3d in units of tau0
    (planes: default TURB_PLANES tau0; an explicit value is taken as is, e.g. from line_intvl).
    adaptive: reset_method = 'decrement' and the file t_end, with which the job script chains the run up to
    t_end with the dt of the last step (dt = t_end/nstep at the start).
    track (Stage 6): reset_method = 'track' (dt follows dt_cfl, src/time.F90) from dt = t_end/nstep; with nrun the
    run takes nrun steps, without nrun it stops at t_end of calliope.in (also in the file t_end for the chain) and
    nstep = TRACK_NMAX nstep is only an upper bound."""
    box = box_3d(N)
    g = gen_ic.Grid(N, N, N)
    fields, tinfo = problem_fields(prob, g)
    tau0 = math.pi/tinfo['Ms']
    t_end = TURB_NTAU*tau0
    dt = t_end/nstep
    to_t_end = track and nrun is None
    nrun = (TRACK_NMAX*nstep if to_t_end else nstep) if nrun is None else nrun
    planes = TURB_PLANES*tau0 if planes is None else planes
    time = dict(dt=dt, nstep=nrun, reset_method="'decrement'") if adaptive else dict(dt=dt, nstep=nrun)
    if track:
        time['reset_method'] = "'track'"
        if to_t_end:
            time['t_end'] = t_end
    path = ic_case(os.path.join(root, name), g, gen_ic.state_from_real(g, *fields), nproc, box, estimate=False,
                   header='Stage 4 decaying isothermal MHD turbulence (gen_stage4.py turb)',
                   info=dict(info, **tinfo, batch='turb', problem=prob, N=N, tau0=tau0, t_end=t_end,
                             nstep=nstep, nrun=nrun, dt=dt, planes=planes, adaptive=adaptive,
                             **(dict(track=True, to_t_end=to_t_end) if track else {})),
                   time=time,
                   shock=dict(C_shock=2.0, C_eta=0.0), svv=dict(svv_eps0=0.5),
                   isoth_diag=dict(write_intvl_planes=planes, write_intvl_cons=cons*tau0, budget=True),
                   diag=dict(write_intvl_3D=out3d*tau0 if out3d > 0 else 0.0))
    if adaptive or to_t_end:
        with open(os.path.join(path, 't_end'), 'w') as f:
            f.write(repr(t_end) + '\n')
    return path


def turb_pilot_nstep(prob, base):
    """nstep on N = 64 from the calib run of prob (calib/<prob>_N64_np8, 100 steps from t = 0)."""
    from analyze_stage1 import load_runs
    cfl = np.vstack(load_runs(os.path.join(base, 'calib', '%s_N64_np8' % prob), 'cfl.dat'))
    dt = TURB_PILOT*float(np.minimum(cfl[:, 2], cfl[:, 3]).min())
    return int(math.ceil(TURB_NTAU*math.pi/TURB_PROBS[prob][0]/dt))


def batch_turb(root):
    """S4.4 first job: M_s = 1, 4 on 64^3 (64 ranks, dt from the calib runs, cfl.dat at every step) and the
    MPI decomposition test of M_s = 4 on 64^3 with its dt: P = 1, 8, 64 over 10 steps (planes at steps 0 and
    10) and P = 1, 8 over TURB_DECOMP_NSTEP steps (conservation.dat at every step, against turb4_N64)."""
    base = os.path.dirname(root)
    cases = [turb_case(root, '%s_N64' % p, p, 64, turb_pilot_nstep(p, base), TURB_NP[64],
                       dict(role='pilot'), cons=0.0) for p in ('turb1', 'turb4')]
    ns = turb_pilot_nstep('turb4', base)
    dt = TURB_NTAU*math.pi/4.0/ns
    for P in TURB_DECOMP_P:
        cases.append(turb_case(root, 'dec_N64_np%d_n10' % P, 'turb4', 64, ns, P, dict(role='decomp'), nrun=10,
                               planes=line_intvl(10, dt, 1), cons=0.0, out3d=0.0))
        if P != TURB_NP[64]:
            cases.append(turb_case(root, 'dec_N64_np%d' % P, 'turb4', 64, ns, P, dict(role='decomp'),
                                   nrun=TURB_DECOMP_NSTEP, planes=0.0, cons=0.0, out3d=0.0))
    return cases


def turb_nstep_from(root, prob, N, Nsrc):
    """nstep on N from the run of prob on Nsrc, as ot_nstep_from."""
    dt = min_dt_cfl(os.path.join(root, '%s_N%d' % (prob, Nsrc)))
    r = N/Nsrc
    return int(math.ceil(TURB_NTAU*math.pi/TURB_PROBS[prob][0]/(TURB_DT_SAFETY*dt/r/1.05**math.log2(r))))


def batch_turb2(root):
    """S4.4 second job: M_s = 1 on 128^3 (128 ranks) with dt from turb1_N64. M_s = 4 on 128^3 is decided after
    turb4_N64_a (batch turb4a): its dt is set by SVV and falls as the fast speed grows (turb4_N64 stopped)."""
    return [turb_case(root, 'turb1_N128', 'turb1', 128, turb_nstep_from(root, 'turb1', 128, 64), TURB_NP[128],
                      dict(role='production', dt_from=64))]


def batch_turb4a(root):
    """S4.4: M_s = 4 on 64^3 with an adaptive dt. turb4_N64 stopped at t = 0.10 (step 1471, reset_method =
    'stop'): lambda_svv grows as cmax^2 (eps_d and F_svv both grow with the fast speed) while rho_min falls, and
    dt_visc fell below its dt. Same initial state, nstep, dt and ranks as turb4_N64 (its first 1470 steps are
    reproduced), with reset_method = 'decrement' and the chain up to t_end (file t_end)."""
    return [turb_case(root, 'turb4_N64_a', 'turb4', 64, turb_pilot_nstep('turb4', os.path.dirname(root)),
                      TURB_NP[64], dict(role='production'), cons=0.0, adaptive=True)]


def batch_turb2a(root):
    """S4.4: turb1_N128 with an adaptive dt. turb1_N128 stopped at t = 2.22 (step 4057, reset_method = 'stop'):
    its dt_visc (SVV) fell to the dt from turb1_N64, as cmax and F_svv grew more on 128^3 (dt lambda_svv 1.61).
    Same initial state, nstep, dt and ranks, with reset_method = 'decrement' and the chain up to t_end."""
    return [turb_case(root, 'turb1_N128_a', 'turb1', 128, turb_nstep_from(root, 'turb1', 128, 64), TURB_NP[128],
                      dict(role='production', dt_from=64), cons=0.0, adaptive=True)]


def batch_turb4b(root):
    """S4.4: M_s = 4 with rms|dB| = B0 (turb4b) on 64^3, adaptive dt from the dt of turb4_N64 (the min dt_cfl
    of the calib run of turb4 halved; the u of turb4b is that of turb4 and its B is weaker)."""
    return [turb_case(root, 'turb4b_N64', 'turb4b', 64, turb_pilot_nstep('turb4', os.path.dirname(root)),
                      TURB_NP[64], dict(role='production'), cons=0.0, adaptive=True)]


def batch_turb4b2(root):
    """S4.4: turb4b on 128^3 (128 ranks) after turb4b_N64, adaptive from the start (the fixed dt from the run on
    64^3 was too large for turb1_N128); the first dt is that of turb_nstep_from."""
    return [turb_case(root, 'turb4b_N128', 'turb4b', 128, turb_nstep_from(root, 'turb4b', 128, 64), TURB_NP[128],
                      dict(role='production', dt_from=64), cons=0.0, adaptive=True)]


SCALE_NSTEP = (10, 110)  # cost per step = (T(110) - T(10))/100: the initialization and the first steps cancel
SCALE_2D = dict(strong=(1024, (16, 32, 64, 128, 256, 512, 1024)), weak=((256, 16), (512, 64), (2048, 1024)))
SCALE_3D = dict(strong=(256, (64, 128, 256, 512, 1024)), weak=((64, 2), (128, 16)))


def scale_case(root, prob, N, P, n, dt):
    """n steps of prob from t = 0 on P ranks with the settings of the production runs (C_shock = 2, C_eta = 0,
    ES-SVV, budget at every step, FFTW_MEASURE) and no output but conservation.dat at the end."""
    box = box_2d(N) if prob == 'ot' else box_3d(N)
    g = gen_ic.Grid(box['nlx'], box['nly'], box['nlz'])
    fields, info = problem_fields(prob, g)
    return ic_case(os.path.join(root, '%s_N%d_np%d_n%d' % (prob, N, P, n)), g, gen_ic.state_from_real(g, *fields),
                   P, box, estimate=False, header='Stage 4 parallel performance (gen_stage4.py scale)',
                   info=dict(info, batch='scale', role='timing', problem=prob, N=N, nrun=n, dt=dt),
                   time=dict(dt=dt, nstep=n), shock=dict(C_shock=2.0, C_eta=0.0), svv=dict(svv_eps0=0.5),
                   isoth_diag=dict(write_intvl_cons=0.0, budget=True))


def scale_cases(root, prob, table, dt_of):
    pts = [(table['strong'][0], P) for P in table['strong'][1]] + list(table['weak'])
    return [scale_case(root, prob, N, P, n, dt_of(N)) for N, P in pts for n in SCALE_NSTEP]


def batch_scale2d(root):
    """S4.5, one job without PAR (select=4: the run on 1024 ranks; the others leave cores idle so that
    no two runs share a node): the OT vortex with dt = 0.5 dt_adv(t = 0) (as calib, which ran without a
    change of dt on 256^2-2048^2)."""
    def dt_of(N):
        g = gen_ic.Grid(N, 2, N)
        return float('%.1e' % (0.5*cfl_dt(g, *problem_fields('ot', g)[0])))
    return scale_cases(root, 'ot', SCALE_2D, dt_of)


def batch_scale3d(root):
    """S4.5, one job without PAR (select=4): the turbulence M_s = 4 with dt = 0.8 min(dt_adv, dt_visc) of
    the calib run on the same N (100 steps from t = 0; dt_visc of the SVV is below 0.5 dt_adv)."""
    from analyze_stage1 import load_runs
    import glob

    def dt_of(N):
        cal, = glob.glob(os.path.join(os.path.dirname(root), 'calib', 'turb4_N%d_np*' % N))
        cfl = np.vstack(load_runs(cal, 'cfl.dat'))
        return float('%.1e' % (0.8*float(np.minimum(cfl[:, 2], cfl[:, 3]).min())))
    return scale_cases(root, 'turb4', SCALE_3D, dt_of)


BATCHES = {'readback': batch_readback, 'regress': batch_regress, 'calib': batch_calib, 'alfven': batch_alfven,
           'cpdiag': batch_cpdiag, 'planar': batch_planar, 'planar2': batch_planar2,
           'planar3': batch_planar3, 'ot': batch_ot, 'ot2': batch_ot2, 'ot3': batch_ot3, 'ot3b': batch_ot3b,
           'turb': batch_turb, 'turb2': batch_turb2, 'turb2a': batch_turb2a, 'turb4a': batch_turb4a,
           'turb4b': batch_turb4b, 'turb4b2': batch_turb4b2, 'scale2d': batch_scale2d, 'scale3d': batch_scale3d}
# batch -> directory of its cases under base (default: the batch name)
ROOT_OF = {'planar2': 'planar', 'planar3': 'planar', 'ot2': 'ot', 'ot3': 'ot', 'ot3b': 'ot', 'turb2': 'turb', 'turb2a': 'turb', 'turb4a': 'turb', 'turb4b': 'turb', 'turb4b2': 'turb', 'scale2d': 'scale', 'scale3d': 'scale'}


def cost(case):
    """nstep * (grid points) * log2(grid points), for listing the cases by decreasing cost."""
    from analyze_stage1 import read_nml
    d = read_nml(os.path.join(case, 'calliope.in'))
    n = int(d['nlx'])*int(d['nly'])*int(d['nlz'])
    return max(int(d['nstep']), 1)*n*math.log2(n)


if __name__ == '__main__':
    batch = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.join('runs', 'stage4')
    print(' '.join(sorted(BATCHES[batch](os.path.join(base, ROOT_OF.get(batch, batch))), key=cost, reverse=True)))
