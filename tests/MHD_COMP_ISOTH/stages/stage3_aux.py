#!/usr/bin/env python3
"""Auxiliary Stage 3 studies of the structures in which B_t reverses (I1, RDpc; REPORT_stage3.md, section 10).

    python3 stage3_aux.py structure              fixed points of the steady structure ODE of F1, S1, I1 and the
                                                 number of positive eigenvalues of its Jacobian
    python3 stage3_aux.py fv PROB N ...          FV runs (fv_mhd1d) of I1 (the setup of stage3.py fvwin) or RDpc
                                                 at N cells -> runs/stage3/fvwin/I1_N<N>_conv.npz, RDpc_N<N>.npz
    python3 stage3_aux.py nut NAMES COMBOS [n]   prototype runs with eta and the tangential viscosity nu_t, e.g.
                                                 nut I1_N256,RDpc_N256 0:0,2:0,0:2,2:2 4 -> runs/stage3/nut
    python3 stage3_aux.py compare                rho L1 of the nut runs against the converged FV solutions and
                                                 between successive N
    python3 stage3_aux.py conv                   the same for the Fortran conv runs at N = 4096, 8192
                                                 (runs/stage3_f/conv_last), with the rho L1 by regions
    python3 stage3_aux.py fig                    fig/stage3_reversal.png (I1, RDpc at N = 1024: C_eta = 0, 2, FV, ideal)

Structure ODE. In the frame of the shock, with the flow along +x (upstream at x = -inf), w = u_x - S > 0 and the
mass flux j = rho w, the steady coplanar equations integrate to D X' = F(X), X = (w, u_y, B_y), D = diag(rho zeta,
rho nu, eta),
    F = (j w + c^2 j/w + B_y^2/2 - Pi,  j u_y - B_x B_y - T,  w B_y - B_x u_y - E)
with the constants Pi, T, E of the upstream state. The Jacobian J of F at a fixed point is symmetric, so the
dimension of the unstable manifold is the number n+(J) of positive eigenvalues of J for every positive D (Sylvester's
law of inertia; checked here with random D). nu = 0 makes u_y algebraic and replaces J by the Schur complement
J2 = [[j (1 - c^2/w^2), B_y], [B_y, w - w_A]], w_A = B_x^2/j, with n+(J) = 1 + n+(J2). With zeta only (nu = eta = 0)
B_y = K/(w - w_A) (K constant): a continuous structure cannot cross w = w_A. The non-coplanar components (u_z, B_z)
add the block [[j, -B_x], [-B_x, w]] (with nu) or w - w_A (without nu). A connection i -> k has dimension
n+(i) - n+(k) generically (1: isolated orbit, 0: none).

nu_t (nut): dm_t/dt += d(rho nu_t d_N u_t), nu_t = C_nu Delta^2 S[max(0, -D)] (the zeta detector), energy
-sum rho nu_t (d_N u_t)^2 dV; fixed dt = t_end/(2 nstep), nstep of the C_eta = 2 prototype run (stage3_nstep.json).
"""
import sys
import os
import json
import time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stage3 as s3
import mhd_exact as mx
from proto_mhd1d import MHD1D

FVWIN = os.path.join(s3.OUT, 'fvwin')
NUT = os.path.join(s3.OUT, 'nut')
# converged FV solutions of the comparisons (compare)
FV_CONV = {'I1': 'I1_N16384_conv.npz', 'RDpc': 'RDpc_N16384.npz'}
# coarser FV solutions (the error of the reference itself) and the final states of the Fortran conv runs (conv)
FV_SELF = {'I1': ('I1_N1024_conv.npz', 'I1_N4096.npz'), 'RDpc': ('RDpc_N256.npz', 'RDpc_N1024.npz', 'RDpc_N4096.npz')}
CONV_LAST = os.path.normpath(os.path.join(s3.OUT, '..', 'stage3_f', 'conv_last'))


def last(W):
    return W[-1] if W.ndim == 3 else W


# ---------------------------------------------------------------------------------------------------------------
# structure
def wave_speeds(rho, by, bx, c):
    ca2, v2 = bx ** 2 / rho, (bx ** 2 + by ** 2) / rho
    d = np.sqrt((c ** 2 + v2) ** 2 - 4 * c ** 2 * ca2)
    return np.sqrt((c ** 2 + v2 - d) / 2), np.sqrt(ca2), np.sqrt((c ** 2 + v2 + d) / 2)


def npos(A):
    return int(np.sum(np.linalg.eigvalsh(A) > 0.0))


def structure(name, c=1.0, rng=np.random.default_rng(1)):
    W1, bx = np.array(mx.COPLANAR_UP), mx.COPLANAR_BX
    r, kind = mx.COPLANAR[name]
    X, S, W2, k = mx.pick(mx.coplanar_shocks(W1, bx, r, c, side=-1), kind)
    rho1, uy1, by1 = W1[0], W1[2], W1[4]
    w1 = W1[1] - S
    j = rho1 * w1
    Pi, T, E = j * w1 + c ** 2 * rho1 + by1 ** 2 / 2, j * uy1 - bx * by1, w1 * by1 - bx * uy1
    wA, K = bx ** 2 / j, E + bx * T / j
    f = lambda w: j * w + c ** 2 * j / w + (K / (w - wA)) ** 2 / 2 - Pi
    ws = w1 * np.logspace(-3, 2, 400001)
    fw = f(ws)
    roots = []
    for i in np.nonzero(np.sign(fw[:-1]) != np.sign(fw[1:]))[0]:
        a, b = ws[i], ws[i + 1]
        if (a - wA) * (b - wA) <= 0.0:          # the pole of B_y
            continue
        for _ in range(200):
            m = 0.5 * (a + b)
            a, b = (m, b) if np.sign(f(m)) == np.sign(f(a)) else (a, m)
        if not np.any(np.isclose(0.5 * (a + b), roots, rtol=1e-9)):   # a root on a sample point is found twice
            roots.append(0.5 * (a + b))
    print('== %s (%s, r = %g): S = %.6f, w_1 = %.6f, w_A = %.6f' % (name, k, r, S, w1, wA))
    print('   type        w      rho      B_y      u_y   |w|/c_s  |w|/c_A  |w|/c_f   n+: coplanar zeta+eta, '
          'zeta+nu+eta; non-coplanar zeta+eta, zeta+nu+eta')
    for w in sorted(roots, reverse=True):
        rho, by = j / w, K / (w - wA)
        uy = (T + bx * by) / j
        cs, ca, cf = wave_speeds(rho, by, bx, c)
        typ = 1 if w > cf else 2 if w > ca else 3 if w > cs else 4
        J = np.array([[j * (1 - c ** 2 / w ** 2), 0.0, by], [0.0, j, -bx], [by, -bx, w]])
        J2 = np.array([[J[0, 0], by], [by, w - wA]])
        Jz = np.array([[j, -bx], [-bx, w]])
        n = (npos(J2), npos(J), npos(J2) + int(w > wA), npos(J) + npos(Jz))
        assert n[1] == 1 + n[0]
        for _ in range(20):                      # Sylvester: n+(D^-1 J) = n+(J) for D > 0
            D = np.diag(rng.uniform(0.1, 10.0, 3))
            assert int(np.sum(np.linalg.eigvals(np.linalg.solve(D, J)).real > 0.0)) == n[1]
        print('   %d    %9.5f %8.5f %+8.5f %+8.5f %8.3f %8.3f %8.3f   %d  %d   %d  %d'
              % (typ, w, rho, by, uy, w / cs, w / ca, w / cf, *n))


# ---------------------------------------------------------------------------------------------------------------
# FV runs
def fv_run(prob, N):
    import fv_mhd1d as fv
    os.makedirs(FVWIN, exist_ok=True)
    t0 = time.time()
    if prob == 'RDpc':
        c = dict(fam='rd', prob='RDp', N=N, w=0, cop=True, C_eta=0.0)
        W, n = fv.run(s3.initial(c), 1.0, s3.L, 1.0)
        np.savez_compressed(os.path.join(FVWIN, 'RDpc_N%d.npz' % N), x=s3.grid(N), W=W,
                            We=s3.rd_W(c, s3.grid(N), 1.0))
    else:
        c = dict(fam='coplanar', prob=prob, N=N, C_eta=0.0)
        W, n = fv.run(s3.initial(c), mx.COPLANAR_BX, s3.L, s3.t_end(c))
        np.savez_compressed(os.path.join(FVWIN, '%s_N%d_conv.npz' % (prob, N)), W=W)
    print(prob, N, n, '%.0fs' % (time.time() - t0), flush=True)


# ---------------------------------------------------------------------------------------------------------------
# nu_t
class NuT(MHD1D):
    def __init__(self, N, C_nu=0.0, **kw):
        super().__init__(N, **kw)
        self.C_nu = C_nu

    def fields(self, sh, mh, bh):
        f = super().fields(sh, mh, bh)
        f['nu'] = self.C_nu * self.Delta ** 2 * f['zeta'] / (self.C * self.Delta ** 2)
        return f

    def rhs(self, sh, mh, bh):
        fs, fm, fb, f = super().rhs(sh, mh, bh)
        if self.C_nu > 0.0:
            for i in (1, 2):
                fm[i] = fm[i] + 1j * self.kd * self.fwd(f['rho'] * f['nu'] * self.ddx(f['uh'][i], self.Md))
        return fs, fm, fb, f


def nut_run(args):
    name, C_eta, C_nu = args
    base = next(c for c in s3.cases() if c['name'] == name)
    nst = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stage3_nstep.json')))
    ref = name + '_eta2' if name + '_eta2' in nst else name
    te = s3.t_end(base)
    nstep = 2 * (nst[ref][0] if ref in nst else s3.nstep_of(base, s3.pilot_dt(base)['dtmin']))
    p = NuT(base['N'], C_nu=C_nu, lx=s3.LX, bx=s3.bx_of(base), C_eta=C_eta)
    p.set_W(s3.initial(base))
    t0 = time.time()
    with np.errstate(all='ignore'):
        st = p.run(te / nstep, nstep)
    np.savez_compressed(os.path.join(NUT, '%s_e%g_n%g.npz' % (name, C_eta, C_nu)), W=p.prim())
    return name, C_eta, C_nu, st, nstep, time.time() - t0


def l1(a, b, i=0):
    """rho L1 between two runs, the finer one sampled at the points of the coarser one."""
    if a.shape[1] > b.shape[1]:
        a, b = b, a
    return float(np.mean(np.abs(a[i] - b[i, ::b.shape[1] // a.shape[1]])))


def compare():
    for prob in ('I1', 'RDpc', 'C2'):
        if prob == 'C2':
            ref = last(s3.fv_reference(prob)[2])
        else:
            ref = last(np.load(os.path.join(FVWIN, FV_CONV[prob]))['W'])
        for combo in ('e0_n0', 'e0.5_n0', 'e2_n0', 'e8_n0', 'e0_n2', 'e2_n2'):
            Ws = {N: np.load(os.path.join(NUT, '%s_N%d_%s.npz' % (prob, N, combo)))['W'] for N in (256, 512, 1024, 2048)
                  if os.path.exists(os.path.join(NUT, '%s_N%d_%s.npz' % (prob, N, combo)))}
            if not Ws:
                continue
            Ns = sorted(Ws)
            print('%-4s %-7s FV %s | self %s' % (prob, combo, ' '.join('%d:%.2e' % (N, l1(Ws[N], ref)) for N in Ns),
                  ' '.join('%d-%d:%.2e' % (a, b, l1(Ws[a], Ws[b])) for a, b in zip(Ns[:-1], Ns[1:]))))
            if prob == 'I1':
                print('              ideal I1 %s' % ' '.join('%d:%.2e' % (N, ideal_I1(Ws[N])[0]) for N in Ns))
    for lab, W in [('FV N16384', last(np.load(os.path.join(FVWIN, FV_CONV['I1']))['W']))] + \
            [('%s N2048' % combo, np.load(os.path.join(NUT, 'I1_N2048_%s.npz' % combo))['W'])
             for combo in ('e0_n0', 'e2_n0', 'e0_n2')]:
        e, up, peak, plat = ideal_I1(W)
        print('I1 %-12s ideal %.2e  ahead (rho, u_x, B_y) = (%.4f, %+.4f, %+.4f)  max rho behind %.4f  '
              'plateau (rho, B_y) = (%.4f, %+.4f)' % (lab, e, *up, peak, *plat))


def ideal_I1(W):
    """rho L1 of W (the final state of I1) against the ideal shock on the valid points (per unit length), the mean
    (rho, u_x, B_y) on 0.01 < x - x_s < 0.08 ahead of the shock, max rho on -0.03 < x - x_s < 0 and the mean
    (rho, B_y) on -0.08 < x - x_s < -0.04."""
    N = W.shape[1]
    c = dict(fam='coplanar', prob='I1', N=N, C_eta=0.0)
    x, te = s3.grid(N), s3.t_end(c)
    We, valid = s3.exact(c, x, te)
    S = mx.coplanar_case('I1')[0]
    d = (x - (s3.cop_interval('I1')[0] + S * te) + 0.5) % 1.0 - 0.5
    up, dn, beh = (d > 0.01) & (d < 0.08), (d > -0.08) & (d < -0.04), (d > -0.03) & (d < 0.0)
    return (float(np.sum(np.abs(W[0] - We[0])[valid]) / N), W[[0, 1, 4]][:, up].mean(1), W[0, beh].max(),
            W[[0, 4]][:, dn].mean(1))


def fig(N=1024):
    """rho, u_y, B_y of I1 (t_end, near the shock) and RDpc (t = 1) at N: C_eta = 0 and 2 with the converged FV
    solution and the ideal one (I1: on the valid points only) -> fig/stage3_reversal.png."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from stage3_table import COLOR, GRAY
    fig, ax = plt.subplots(2, 3, figsize=(10.5, 5.6))
    for r, prob in enumerate(('I1', 'RDpc')):
        fv = np.load(os.path.join(FVWIN, FV_CONV[prob]))
        Wf = last(fv['W'])
        xf, x = s3.grid(Wf.shape[1]), s3.grid(N)
        if prob == 'I1':
            c = dict(fam='coplanar', prob='I1', N=N, C_eta=0.0)
            We, valid = s3.exact(c, x, s3.t_end(c))
            We = np.where(valid, We, np.nan)
            xs = s3.cop_interval('I1')[0] + mx.coplanar_case('I1')[0] * s3.t_end(c)
            lim = (xs - 0.12, xs + 0.1)
        else:
            We, lim = fv['We'][:, ::Wf.shape[1] // N], (0.0, 0.5)
        for j, (k, nm) in enumerate(((0, 'ρ'), (2, 'u_y'), (4, 'B_y'))):
            a = ax[r, j]
            a.plot(x, We[k], color=GRAY, lw=1.0, ls=':', label='ideal')
            a.plot(xf, Wf[k], color=GRAY, lw=1.5, label='FV N = %d' % Wf.shape[1])
            for tag, col, lab in (('', COLOR['C'], 'bulk'), ('_eta2', COLOR['B'], 'C_eta = 2')):
                a.plot(x, np.load(os.path.join(s3.OUT, '%s_N%d%s.npz' % (prob, N, tag)))['W'][-1][k], color=col,
                       lw=1.0, label='%s, N = %d' % (lab, N))
            a.set_xlim(*lim)
            a.set_ylabel(nm)
            if r == 1:
                a.set_xlabel('x')
            if j == 0:
                a.set_title('%s, t = %.3g' % (prob, 0.24 if prob == 'I1' else 1.0), loc='left')
    ax[0, 1].legend(frameon=False, fontsize=7, loc='center left')
    fig.tight_layout()
    fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig', 'stage3_reversal.png'), dpi=150)
    plt.close(fig)


def regions(prob, N):
    """Regions of the rho L1 of conv: I1 by d = x - x_s (shock, ahead, behind as in ideal_I1, the rest), RDpc by the
    distance from the RDs at x = 0.25, 0.75."""
    x = s3.grid(N)
    if prob == 'I1':
        c = dict(fam='coplanar', prob='I1', N=N, C_eta=0.0)
        d = (x - (s3.cop_interval('I1')[0] + mx.coplanar_case('I1')[0] * s3.t_end(c)) + 0.5) % 1.0 - 0.5
        return [('shock', np.abs(d) < 0.01), ('ahead', (d > 0.01) & (d < 0.08)),
                ('behind', (d > -0.08) & (d < -0.01)), ('rest', np.abs(d) >= 0.08)]
    dist = np.min(np.abs((x[:, None] - np.array([0.25, 0.75])[None, :] + 0.5) % 1.0 - 0.5), axis=1)
    return [('[%g,%g)' % (lo, hi), (dist >= lo) & (dist < hi)) for lo, hi in ((0, 0.01), (0.01, 0.05), (0.05, 0.15),
                                                                              (0.15, 1.0))]


def conv():
    """The Fortran conv runs (final states in CONV_LAST: analyze_stage3.py last) after the nut runs at N <= 2048:
    L1 against the converged FV solution (rho and the sum over the 6 variables; I1 also against the ideal shock) and
    the rho L1 split by regions (sum |rho - rho_FV|/N over each region, rho_FV the mean over each cell), with the
    coarser FV solutions for the error of the reference itself."""
    for prob in ('I1', 'RDpc'):
        ref = last(np.load(os.path.join(FVWIN, FV_CONV[prob]))['W'])
        rows = [('FV', last(np.load(os.path.join(FVWIN, f))['W'])) for f in FV_SELF[prob]]
        for tag, combo in (('', 'e0_n0'), ('_eta2', 'e2_n0')):
            fns = [('proto', os.path.join(NUT, '%s_N%d_%s.npz' % (prob, N, combo))) for N in (256, 512, 1024, 2048)]
            fns += [('fort', os.path.join(CONV_LAST, '%s_N%d%s.npz' % (prob, N, tag))) for N in (4096, 8192)]
            rows += [(src + tag, np.load(fn)['W']) for src, fn in fns if os.path.exists(fn)]
        for lab, W in rows:
            N = W.shape[1]
            line = '%-4s %-11s N=%5d  FV rho %.2e  all %.2e' % (prob, lab, N, l1(W, ref),
                                                              sum(l1(W, ref, i) for i in range(6)))
            if prob == 'I1':
                e, up = ideal_I1(W)[:2]
                line += '  ideal %.2e  ahead (%.4f %+.4f %+.4f)' % (e, *up)
            a = np.abs(W[0] - ref.reshape(6, N, -1).mean(2)[0]) / N
            print(line + '  |  ' + '  '.join('%s %.2e' % (p, a[m].sum()) for p, m in regions(prob, N)))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'fig':
        fig()
    if cmd == 'structure':
        for name in ('F1', 'S1', 'I1'):
            structure(name)
    elif cmd == 'fv':
        for N in map(int, sys.argv[3:]):
            fv_run(sys.argv[2], N)
    elif cmd == 'nut':
        import multiprocessing as mp
        os.makedirs(NUT, exist_ok=True)
        jobs = [(n, *map(float, s.split(':'))) for n in sys.argv[2].split(',') for s in sys.argv[3].split(',')]
        jobs.sort(key=lambda j: -int(j[0].split('_N')[1]))
        with mp.Pool(int(sys.argv[4]) if len(sys.argv) > 4 else 2) as pool:
            for r in pool.imap_unordered(nut_run, jobs):
                print('%-12s eta %-4g nu %-4g %-6s nstep %5d  %.0fs' % r, flush=True)
    elif cmd == 'compare':
        compare()
    elif cmd == 'conv':
        conv()
