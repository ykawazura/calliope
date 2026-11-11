#!/usr/bin/env python3
"""Entropy-stable SVV of rho and m (svv_form = 'entropy') on a 3D state, for the static analysis of Stage 5.

Port of the eps_d of advance.F90 and get_svv_eps, and of svv_es_core, svv_G and svv_power
(shock_viscosity.F90), with variants of the coefficients:

    eps     'global' : eps_d = svv_eps0 dl_d max(|u_d| + c_f,d)                 (the code)
            'local'  : eps_d(x) = svv_eps0 dl_d (|u_d| + c_f,d)(x)               (svv_eps_local = T)
    weight  'rho'    : rho                                                       (the code)
            'min<n>' : -S_max^n(-rho), S_max = 3-point max over y, z, x          (svv_rho_nmin = n)
            'const<r>': the constant r (static test only)

The weight W replaces rho in a = eps W/cs^2, b = eps W and kap W; kap = eps/(cs^2 + max|u|^2).
For any W >= 0 and eps >= 0 the matrix K = eps (W/rho) dU/dw is positive semidefinite at each point, so
sum (w_r d rho/dt + u.dm/dt) = -D_svv <= 0, and sum d rho/dt = sum dm/dt = 0 (sig_d = 0 at k_d = 0).
eps and the constant W are fixed over the power iterations; 'min<n>' is computed from the rho of each
evaluation, as in svv_G.

F = lambda/sum_d eps_d,global lambda_s,d is the ratio of the largest |eigenvalue| of the Jacobian of
G = ((d rho/dt)/rho, dm/dt) to the linear bound of the code (column 14 of cfl.dat); F = 1 for a uniform state.

    python3 svv3d.py selftest                      uniform state (F = 1), identities on a smooth state
    python3 svv3d.py synth [N]                     F against U/c_s for 1D velocity profiles (rho = 1)
    python3 svv3d.py hole [N]                      F against the depth of a hole at rest, and its power of rho_h
    python3 svv3d.py kernel [N ...]                share of the L and H parts of the SVV near a point (32^3, 64^3)
    python3 svv3d.py diag [N | <case> [sub]] [--json f]
                                                   the diagonal of the Jacobian at the centre of that hole (or at
                                                   the min of rho of a restart), against F
    python3 svv3d.py bound [N | <case> [sub]] [--json f]
                                                   at the same points, lower bounds of lambda from the frozen SVV:
                                                   Q (trial vector w_rho at the point), Q4 (Rayleigh-Ritz over the
                                                   4 point perturbations of rho, m) and the fill R
    python3 svv3d.py state <case> [sub] [--json f] [--eig]
                                                   the void of a restart (sub = restart, ic or chunk_k/restart):
                                                   d rho/dt at the min of rho, F of the variants, the predictor
    python3 svv3d.py lam <case> <eps> <weight> [niter [sub]]
                                                   lambda after niter (default svv_pi_init = 50) iterations from
                                                   the initial vector of the code, as in the first cfl.dat row
                                                   of a run started from <case>/<sub>
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_ic  # noqa: E402
from analyze_stage1 import read_nml, num, strv  # noqa: E402

# threads of the FFTs (SVV3D_WORKERS; few by default, for a shared login node)
WORKERS = int(os.environ.get('SVV3D_WORKERS', '4'))
try:
    import scipy.fft as _sfft

    def _rfftn(a):
        return _sfft.rfftn(a, axes=(0, 1, 2), workers=WORKERS)

    def _irfftn(a, s):
        return _sfft.irfftn(a, s=s, axes=(0, 1, 2), workers=WORKERS)
except ImportError:  # pragma: no cover
    def _rfftn(a):
        return np.fft.rfftn(a, axes=(0, 1, 2))

    def _irfftn(a, s):
        return np.fft.irfftn(a, s=s, axes=(0, 1, 2))

H_FD = 1e-7  # svv_power


class SVV:
    """Operators on the grid g (gen_ic.Grid; real arrays [i, k, j] = (x, z, y), Fourier (nkx, nkz, nky))."""

    def __init__(self, g, cs2, eps0=0.5, mN=0.5, above_kc=True, filter_mom=True, mass=True):
        self.g, self.cs2, self.eps0, self.filter_mom, self.mass = g, cs2, eps0, filter_mom, mass
        self.shape = (g.nlx, g.nlz, g.nly)
        self.ntot = g.ntot
        self.active = [g.nlx > 2, g.nly > 2, g.nlz > 2]
        self.dl = [g.Lx/g.nlx, g.Ly/g.nly, g.Lz/g.nlz]
        # (k, ik, nl, L, shape) of x, y, z
        dirs = ((g.kx, g.ikx, g.nlx, g.Lx, (-1, 1, 1)), (g.ky, g.iky, g.nly, g.Ly, (1, 1, -1)),
                (g.kz, g.ikz, g.nlz, g.Lz, (1, -1, 1)))
        self.kq2, self.sig, self.lam_s = [], [], []
        for d, (k, ik, n, L, shp) in enumerate(dirs):
            k = k.ravel()
            kc = 2*np.pi*(n//2)/L*2/3  # svv_kernel: kmax_noprune*2/3
            m = mN*kc
            a = np.abs(k)
            with np.errstate(divide='ignore', invalid='ignore'):
                r = (kc - a)/(a - m)
                q = np.where(r < 30, np.exp(-np.where(r < 30, r, 0)**2), 0.0)
            Q = np.where(a >= kc, 1.0 if above_kc else 0.0, np.where(a > m, q, 0.0))
            kq2 = k**2*Q if self.active[d] else np.zeros_like(k)
            sg = np.where(np.abs(ik) == n//2, np.sqrt(kq2) + 0j, 1j*np.where(k >= 0, 1, -1)*np.sqrt(kq2))
            self.kq2.append(kq2.reshape(shp))
            self.sig.append(sg.reshape(shp))
            self.lam_s.append(float(kq2.max()))  # filter_lnrho = F: every mode is kept in ln rho
        self.mL = g.mask if filter_mom else np.ones_like(g.mask)
        self.mH = 1 - self.mL

    # FFTs normalized as the code: f_k = FFT(f)/N, f = sum_k f_k e^{ikx}
    def fft(self, f):
        return _rfftn(f)/self.ntot

    def ifft(self, fk):
        return _irfftn(fk*self.ntot, self.shape)

    def cmax(self, rho, u, b):
        """(|u_d| + c_f,d)(x) of advance.F90 for d = x, y, z."""
        va2 = sum(bi**2 for bi in b)/rho
        out = []
        for d in range(3):
            vad2 = b[d]**2/rho
            cf2 = 0.5*(self.cs2 + va2 + np.sqrt(np.maximum(0.0, (self.cs2 + va2)**2 - 4*self.cs2*vad2)))
            out.append(np.abs(u[d]) + np.sqrt(cf2))
        return out

    def eps(self, rho, u, b, local=False):
        c = self.cmax(rho, u, b)
        return [self.eps0*self.dl[d]*(c[d] if local else float(c[d].max())) if self.active[d] else 0.0
                for d in range(3)]

    def lam_lin(self, eps_g):
        return sum(e*l for e, l in zip(eps_g, self.lam_s))

    def smin(self, f, n):
        """-S_max^n(-f): n passes of the 3-point min over y, z, x (pass_1d order) in the active directions."""
        for _ in range(n):
            for ax, d in ((2, 1), (1, 2), (0, 0)):
                if self.active[d]:
                    f = np.minimum(np.minimum(np.roll(f, 1, ax), f), np.roll(f, -1, ax))
        return f

    def weight(self, rho, kind):
        if kind == 'rho':
            return rho
        if kind.startswith('min'):
            return self.smin(rho, int(kind[3:]))
        if kind.startswith('const'):
            return np.full(self.shape, float(kind[5:]))
        raise ValueError(kind)

    def _grad(self, xk, d, low):
        return self.ifft(self.sig[d]*(self.mL if low else self.mH)*xk)

    def _adj(self, y, d, low):
        return np.conj(self.sig[d])*(self.mL if low else self.mH)*self.fft(y)

    def tend(self, s, u, W, eps, parts=False):
        """svv_es_core: (d rho/dt, dm/dt (3 real arrays), D_svv as the grid sum); with parts also the
        A_L and A_H parts of d rho/dt. mass = False (svv_mass = F): f_r = f_h = 0 and f_m = b g_m."""
        cs2 = self.cs2
        umax2 = float(np.max(sum(x**2 for x in u)))
        wk = self.fft(cs2*s - 0.5*sum(x**2 for x in u)) if self.mass else None
        uk = [self.fft(x) for x in u]
        rL = np.zeros_like(uk[0])
        rH = np.zeros_like(uk[0])
        mk = [np.zeros_like(uk[0]) for _ in range(3)]
        dsum = 0.0
        for d in range(3):
            e = eps[d]
            if np.all(np.asarray(e) == 0):
                continue
            gm = [self._grad(x, d, True) for x in uk]
            b = e*W
            if not self.mass:
                dsum += float(np.sum(b*(gm[0]**2 + gm[1]**2 + gm[2]**2)))
                for c in range(3):
                    mk[c] -= self._adj(b*gm[c], d, True)
                continue
            gr = self._grad(wk, d, True)
            a = e*W/cs2
            t = gr + u[0]*gm[0] + u[1]*gm[1] + u[2]*gm[2]
            dsum += float(np.sum(a*t**2 + b*(gm[0]**2 + gm[1]**2 + gm[2]**2)))
            fr = a*t
            rL -= self._adj(fr, d, True)
            for c in range(3):
                mk[c] -= self._adj(u[c]*fr + b*gm[c], d, True)
            if self.filter_mom:
                gh = self._grad(wk, d, False)
                kap = e/(cs2 + umax2)
                dsum += float(np.sum(kap*W*gh**2))
                rH -= self._adj(kap*W*gh, d, False)
        rdot = self.ifft(rL + rH)
        mdot = [self.ifft(x) for x in mk]
        if parts:
            return rdot, mdot, dsum, self.ifft(rL), self.ifft(rH)
        return rdot, mdot, dsum

    def G(self, s, m, eps, wkind, wconst=None):
        """svv_G: ((d rho/dt)/rho, dm/dt) of the state (s, m) with eps fixed."""
        rho = np.exp(s)
        u = [x/rho for x in m]
        W = wconst if wconst is not None else self.weight(rho, wkind)
        rdot, mdot, _ = self.tend(s, u, W, eps)
        return [rdot/rho] + mdot

    def v0(self):
        """Initial vector of svv_power: the s part, and m = 0 except for mass = False, where each m component is
        prod_d cos(2 pi K_d i_d/nl_d + 1) with K_d = (nl_d - 1)//3, the top index kept by the 2/3 filter on m."""
        g = self.g
        amp = np.ones(self.shape)
        pm = np.ones(self.shape)
        nsg = np.zeros(self.shape, dtype=int)
        idx = np.indices(self.shape)  # (i, k, j)
        for d, ax, n in ((0, 0, g.nlx), (1, 2, g.nly), (2, 1, g.nlz)):
            if self.active[d]:
                amp += 0.3*np.cos(2*np.pi*idx[ax]/n + 1.0)
                nsg += idx[ax]
                pm *= np.cos(2*np.pi*(((n - 1)//3)*idx[ax])/n + 1.0)
        v = np.zeros((4,) + self.shape)
        v[0] = np.where(nsg % 2 == 0, 1.0, -1.0)*amp
        if not self.mass:
            v[1:] = pm
        return v

    def norm(self, v):
        return float(np.sqrt(np.sum(v**2)/self.ntot))

    def power(self, s, m, eps, wkind, niter=50, tol=0.0, v=None, wconst=None):
        """svv_power: (lambda, iterations, converged, vector, history of lambda)."""
        G0 = np.array(self.G(s, m, eps, wkind, wconst))
        v = self.v0() if v is None else v.copy()
        lam_old, lam, hist, conv, it = 0.0, 0.0, [], False, 0
        for it in range(1, niter + 1):
            nv = self.norm(v)
            if nv == 0:
                break
            v /= nv
            Gp = np.array(self.G(s + H_FD*v[0], [m[c] + H_FD*v[c + 1] for c in range(3)], eps, wkind, wconst))
            v = (Gp - G0)/H_FD
            lam = self.norm(v)
            hist.append(lam)
            conv = tol > 0 and it > 1 and abs(lam - lam_old) <= tol*lam
            lam_old = lam
            if conv:
                break
        return lam, it, conv, v, hist

    def lam_first(self, s, m, eps, wkind, h, v=None, wconst=None):
        """lambda of one power iteration from v (default v0) with the finite-difference step h instead of H_FD:
        it scales as 1/h if J v is round-off (v in the null space of J)."""
        global H_FD
        h0, H_FD = H_FD, h
        try:
            return self.power(s, m, eps, wkind, niter=1, v=v, wconst=wconst)[0]
        finally:
            H_FD = h0

    def eig(self, s, m, eps, wkind, wconst=None, k=3, tol=1e-6, ncv=24, maxiter=2000):
        """Largest |eigenvalue| of the finite-difference Jacobian (the J v of svv_power) by ARPACK
        (k eigenvalues of largest magnitude). Returns (|lambda_1|, the k eigenvalues, number of J v,
        the eigenvector of lambda_1 as (4,) + shape)."""
        from scipy.sparse.linalg import LinearOperator, eigs
        G0 = np.array(self.G(s, m, eps, wkind, wconst))
        shp = G0.shape
        cnt = [0]

        def mv(x):
            cnt[0] += 1
            v = np.asarray(x).reshape(shp)
            nv = self.norm(v)
            if nv == 0:
                return np.zeros(G0.size)
            v = v/nv
            Gp = np.array(self.G(s + H_FD*v[0], [m[c] + H_FD*v[c + 1] for c in range(3)], eps, wkind, wconst))
            return ((Gp - G0)/H_FD*nv).ravel()

        A = LinearOperator((G0.size, G0.size), matvec=mv, dtype=float)
        vals, vecs = eigs(A, k=k, which='LM', tol=tol, v0=self.v0().ravel(), ncv=ncv, maxiter=maxiter)
        i = int(np.argmax(np.abs(vals)))
        return float(np.abs(vals[i])), vals, cnt[0], np.abs(vecs[:, i]).reshape(shp)

    def predictor(self, rho, u, W, eps, eps_g, R=2):
        """Lambda = max over x, y = x +- r e_d (r = 0..R, d active) of
        (eps_d(x)/eps_d,global) (W(x)/rho(y)) (q/2 + sqrt(1 + q^2/4))^2, q = |u(x) - u(y)|/c_s:
        the largest eigenvalue of K(x)/eps_d,global times H(y) = dw/dU(y) for one pair of points,
        since dU/dw(x) (dU/dw(y))^-1 is similar to (rho(x)/rho(y)) [[1, q^T], [q, q q^T + I]].
        Returns (Lambda, its x, its y, max of the ratio factor alone, max of the u factor alone)."""
        cs = np.sqrt(self.cs2)
        best = (0.0, None, None)
        rmax = umax = 0.0
        for d, ax in ((0, 0), (1, 2), (2, 1)):
            if not self.active[d]:
                continue
            er = (np.asarray(eps[d])/eps_g[d])*np.ones(self.shape)
            for r in range(-R, R + 1):
                # values at y = x + r e_d, on the grid of x
                rho_y = np.roll(rho, -r, ax)
                q = np.sqrt(sum((u[c] - np.roll(u[c], -r, ax))**2 for c in range(3)))/cs
                uf = (q/2 + np.sqrt(1 + q**2/4))**2
                rat = er*W/rho_y
                lam = rat*uf
                i = int(np.argmax(lam))
                rmax = max(rmax, float(rat.max()))
                umax = max(umax, float(uf.max()))
                if lam.flat[i] > best[0]:
                    x = np.unravel_index(i, self.shape)
                    y = list(x)
                    y[ax] = (y[ax] + r) % self.shape[ax]
                    best = (float(lam.flat[i]), tuple(int(a) for a in x), tuple(int(a) for a in y))
        return best + (rmax, umax)


def from_case(case):
    """(Grid, SVV) with the parameters of <case>/calliope.in (defaults of params.F90)."""
    d = read_nml(os.path.join(case, 'calliope.in'))
    g = gen_ic.Grid(int(d['nlx']), int(d['nly']), int(d['nlz']),
                    num(d.get('lx', '1.0')), num(d.get('ly', '1.0')), num(d.get('lz', '1.0')))
    lg = {'.true.': True, '.false.': False, 't': True, 'f': False}
    op = SVV(g, num(d.get('beta0', '1.0'))/2, num(d.get('svv_eps0', '0.5')), num(d.get('svv_mn', '0.5')),
             lg[d.get('svv_above_kc', '.true.').lower()], lg[d.get('filter_mom', '.true.').lower()],
             lg[d.get('svv_mass', '.true.').lower()])
    assert strv(d.get('svv_form', "'entropy'")) == 'entropy'
    assert lg[d.get('filter_lnrho', '.false.').lower()] is False
    return g, op


def real_state(g, lk, mk, bk):
    """s, m, rho, u, B in real space as the code has them (m and B masked in Fourier space by the state)."""
    s = g.ifft(lk)
    rho = np.exp(s)
    m = [g.ifft(x) for x in mk]
    u = [x/rho for x in m]
    b = [g.ifft(x) for x in bk]
    return s, m, rho, u, b


VARIANTS = [('global', 'rho'), ('local', 'rho'), ('global', 'const1'), ('global', 'constmean'),
            ('global', 'min1'), ('global', 'min2'), ('global', 'min3'), ('local', 'min1'), ('local', 'min2'),
            ('local', 'min3')]


def lam_variant(op, s, m, rho, u, b, ek, wk, niter=400, tol=1e-5, method='power'):
    """F and lambda of the variant (ek, wk) at the state; the constant weight 'constmean' is <rho>.
    method = 'power' (svv_power, niter iterations at most) or 'eig' (ARPACK, converged)."""
    eg = op.eps(rho, u, b)
    eps = op.eps(rho, u, b, local=(ek == 'local'))
    wconst = None
    if wk == 'constmean':
        wconst = np.full(op.shape, float(rho.mean()))
    elif wk.startswith('const'):
        wconst = np.full(op.shape, float(wk[5:]))
    t0 = time.time()
    if method == 'eig':
        lam, vals, nit, v = op.eig(s, m, eps, wk, wconst=wconst)
        conv, hist = True, [complex(x) for x in vals]
        hist = [[x.real, x.imag] for x in hist]
    else:
        lam, nit, conv, v, hist = op.power(s, m, eps, wk, niter, tol, wconst=wconst)
    W = wconst if wconst is not None else op.weight(rho, wk)
    P = op.predictor(rho, u, W, eps, eg)
    P1 = op.predictor(rho, u, W, eps, eg, R=1)
    # where the eigenvector lives: the cell of max |v_s| + |v_m|/rho... reported as max of the s part and of |v_m|
    vs = np.abs(v[0])
    vm = np.sqrt(v[1]**2 + v[2]**2 + v[3]**2)
    return dict(eps=ek, weight=wk, lam=lam, F=lam/op.lam_lin(eg), nit=nit, conv=bool(conv), sec=time.time() - t0,
                lam_hist_tail=hist[-3:], Lambda=P[0], Lambda_x=P[1], Lambda_y=P[2], ratio_max=P[3], ufac_max=P[4],
                Lambda_R1=P1[0],
                vs_argmax=tuple(int(a) for a in np.unravel_index(np.argmax(vs), vs.shape)),
                vm_argmax=tuple(int(a) for a in np.unravel_index(np.argmax(vm), vm.shape)),
                W_at_min=float(W.flat[int(np.argmin(rho))]))


# ---------------------------------------------------------------------------------------------- tests

def smooth_state(g, amp_s=0.4, amp_u=0.8, amp_b=0.6):
    X, Y, Z = g.X, g.Y, g.Z
    lnrho = amp_s*(np.sin(X + 2*Z + 0.3) + 0.5*np.cos(2*X - Z + Y) + 0.3*np.sin(Y - Z))
    u = (amp_u*np.sin(Z + 0.2*Y), amp_u*0.7*np.cos(X - Z), amp_u*0.5*np.sin(X + Y + 1))
    b = (amp_b*np.cos(Y - 2*Z), amp_b*np.sin(X + Z), 1.0 + amp_b*0.5*np.cos(X - Y))
    return gen_ic.state_from_real(g, lnrho, u, b)


def selftest():
    ok = True
    # 1. uniform state with u != 0 and B != 0: F = 1 for each variant
    g = gen_ic.Grid(16, 16, 16)
    op = SVV(g, 1.0)
    rho0, u0, b0 = 1.3, (0.7, -0.4, 1.1), (0.3, 0.2, 1.0)
    lk, mk, bk = gen_ic.state_from_real(g, np.full(op.shape, np.log(rho0)), tuple(np.full(op.shape, x) for x in u0),
                                        tuple(np.full(op.shape, x) for x in b0))
    s, m, rho, u, b = real_state(g, lk, mk, bk)
    for ek, wk in [('global', 'rho'), ('local', 'rho'), ('global', 'min2'), ('local', 'min2'),
                   ('global', 'const%.17g' % rho0)]:
        r = lam_variant(op, s, m, rho, u, b, ek, wk, niter=300, tol=1e-9)
        good = abs(r['F'] - 1) < 1e-6
        ok &= good
        print('uniform 16^3  eps %-6s weight %-8s F - 1 = %+.2e (%d it)  Lambda %.6f  %s'
              % (ek, wk, r['F'] - 1, r['nit'], r['Lambda'], 'PASS' if good else 'FAIL'))
    # 2. smooth state: sum (w_r rdot + u.mdot) = -D_svv, sum rdot = sum mdot = 0
    for N in (16, 24):
        g = gen_ic.Grid(N, N, N)
        op = SVV(g, 1.0)
        s, m, rho, u, b = real_state(g, *smooth_state(g))
        for ek, wk in [('global', 'rho'), ('local', 'rho'), ('global', 'min1'), ('global', 'min3'), ('local', 'min2'),
                       ('global', 'const1')]:
            eps = op.eps(rho, u, b, local=(ek == 'local'))
            W = op.weight(rho, wk)
            rdot, mdot, dsum = op.tend(s, u, W, eps)
            wr = op.cs2*s - 0.5*sum(x**2 for x in u)
            prod = wr*rdot + sum(u[c]*mdot[c] for c in range(3))
            res = (np.sum(prod) + dsum)/np.sum(np.abs(prod))
            scale = np.sum(np.abs(rdot)) + sum(np.sum(np.abs(x)) for x in mdot)
            cons = max(abs(np.sum(rdot)), *(abs(np.sum(x)) for x in mdot))/scale
            good = abs(res) < 1e-12 and cons < 1e-13 and dsum > 0
            ok &= good
            print('smooth %2d^3  eps %-6s weight %-6s  (sum w.G + D)/sum|w.G| = %+.1e  max|sum G|/sum|G| = %.1e  '
                  'D = %.4e  %s' % (N, ek, wk, res, cons, dsum, 'PASS' if good else 'FAIL'))
    # 3. the m part of the tendency lies in the mask (A_L^* carries M)
    g = gen_ic.Grid(16, 16, 16)
    op = SVV(g, 1.0)
    s, m, rho, u, b = real_state(g, *smooth_state(g))
    _, mdot, _ = op.tend(s, u, rho, op.eps(rho, u, b))
    out = max(np.max(np.abs(op.fft(x)*op.mH)) for x in mdot)/max(np.max(np.abs(op.fft(x))) for x in mdot)
    good = out < 1e-14
    ok &= good
    print('mask of dm/dt: max outside / max = %.1e  %s' % (out, 'PASS' if good else 'FAIL'))
    ok &= selftest_nomass()
    print('all PASS' if ok else 'some FAIL')
    return ok


def selftest_nomass():
    """mass = False (svv_mass = F): uniform state lambda = max over the mask of m of sum_d eps_d k_d^2 Q_d (not
    F = 1: lam_lin includes the modes |k_d| >= k_c of ln rho, which have no SVV here; Stage 13 corrected this
    criterion with the approval of the user, F - 1 is printed as well); on the smooth state rdot = 0,
    sum u.mdot = -D_svv, D_svv = sum_d sum b |A_L u|^2 > 0, sum mdot = 0; lambda of 50 power iterations from v0
    differs from mass = True; on a state with s = 0 and u, B below the SVV band, J v = 0 for the s part of v0 (one
    iteration gives lambda ~ 1/h), not for v0 (its m part), and 50 iterations reach lambda_mask."""
    ok = True
    # 1. uniform state with u != 0 and B != 0: J_mm = sum_d eps_d A_L^* A_L, whose largest eigenvalue is the max
    # over the mask of m (2/3 rule) of sum_d eps_d k_d^2 Q_d
    g = gen_ic.Grid(16, 16, 16)
    op = SVV(g, 1.0, mass=False)
    rho0, u0, b0 = 1.3, (0.7, -0.4, 1.1), (0.3, 0.2, 1.0)
    lk, mk, bk = gen_ic.state_from_real(g, np.full(op.shape, np.log(rho0)), tuple(np.full(op.shape, x) for x in u0),
                                        tuple(np.full(op.shape, x) for x in b0))
    s, m, rho, u, b = real_state(g, lk, mk, bk)
    for ek, wk in [('global', 'rho'), ('local', 'rho'), ('global', 'min2'), ('local', 'min2'),
                   ('global', 'const%.17g' % rho0)]:
        r = lam_variant(op, s, m, rho, u, b, ek, wk, niter=300, tol=1e-9)
        eps = op.eps(rho, u, b, local=(ek == 'local'))
        lam_mask = float(np.max(op.mL*sum(np.max(np.asarray(eps[d]))*op.kq2[d] for d in range(3))))
        res = r['lam']/lam_mask - 1
        good = abs(res) < 1e-6
        ok &= good
        print('mass=F uniform 16^3  eps %-6s weight %-8s lambda/lambda_mask - 1 = %+.2e (%d it)  F - 1 = %+.4f  %s'
              % (ek, wk, res, r['nit'], r['F'] - 1, 'PASS' if good else 'FAIL'))
    # 2. smooth state
    for N in (16, 24):
        g = gen_ic.Grid(N, N, N)
        op = SVV(g, 1.0, mass=False)
        s, m, rho, u, b = real_state(g, *smooth_state(g))
        for ek, wk in [('global', 'rho'), ('local', 'rho'), ('global', 'min1'), ('global', 'min3'), ('local', 'min2'),
                       ('global', 'const1')]:
            eps = op.eps(rho, u, b, local=(ek == 'local'))
            W = op.weight(rho, wk)
            rdot, mdot, dsum = op.tend(s, u, W, eps)
            prod = sum(u[c]*mdot[c] for c in range(3))
            res = (np.sum(prod) + dsum)/np.sum(np.abs(prod))
            cons = max(abs(np.sum(x)) for x in mdot)/sum(np.sum(np.abs(x)) for x in mdot)
            dref = sum(float(np.sum(eps[d]*W*sum(op._grad(op.fft(x), d, True)**2 for x in u))) for d in range(3))
            dd = abs(dsum - dref)/dref
            good = bool(np.all(rdot == 0)) and abs(res) < 1e-12 and cons < 1e-13 and dsum > 0 and dd < 1e-12
            ok &= good
            print('mass=F smooth %2d^3  eps %-6s weight %-6s  max|rdot| = %.0e  (sum u.G + D)/sum|u.G| = %+.1e  '
                  'max|sum G|/sum|G| = %.1e  |D/D_ref - 1| = %.1e  D = %.4e  %s'
                  % (N, ek, wk, np.max(np.abs(rdot)), res, cons, dd, dsum, 'PASS' if good else 'FAIL'))
    # 3. lambda after 50 iterations from v0 (as at the first step of the code) with and without the mass fluxes
    g = gen_ic.Grid(16, 16, 16)
    s, m, rho, u, b = real_state(g, *smooth_state(g))
    lam = {}
    for mass in (True, False):
        op = SVV(g, 1.0, mass=mass)
        lm, _, _, _, hist = op.power(s, m, op.eps(rho, u, b), 'rho', niter=50)
        lam[mass] = lm
        conv = abs(hist[-1] - hist[-2])/hist[-1]
        good = conv <= 1e-3
        ok &= good
        print('smooth 16^3 mass=%s: lambda_50 = %.10f  |lambda_50/lambda_49 - 1| = %.1e  %s'
              % ('TF'[not mass], lm, conv, 'PASS' if good else 'FAIL'))
    diff = abs(lam[False] - lam[True])/lam[True]
    good = diff >= 1e-4
    ok &= good
    print('smooth 16^3 |lambda(F)/lambda(T) - 1| = %.3e (>= 1e-4)  %s' % (diff, 'PASS' if good else 'FAIL'))
    # 4. s = 0, u and B at |k_d| <= 2 < svv_mN k_c = 4 (24^3): A_L(u v_s) = 0, so J v = 0 for the s part of v0;
    # (a) shows that the test (b) of v0 can fail
    g = gen_ic.Grid(24, 24, 24)
    X, Y, Z = g.X, g.Y, g.Z
    u = (0.8*np.sin(Z + 2*Y + 0.2), 0.5*np.cos(2*X - Z), 0.4*np.sin(X + Y + 1))
    b = (0.6*np.cos(Y - 2*Z), 0.6*np.sin(X + Z), 1.0 + 0.3*np.cos(X - Y))
    s, m, rho, u, b = real_state(g, *gen_ic.state_from_real(g, np.zeros(X.shape), u, b))
    op = SVV(g, 1.0, mass=False)
    eps = op.eps(rho, u, b)
    vs = op.v0()
    vs[1:] = 0.0
    l7, l5 = op.lam_first(s, m, eps, 'rho', 1e-7, v=vs), op.lam_first(s, m, eps, 'rho', 1e-5, v=vs)
    good = l7 > 10*l5
    ok &= good
    print('low-k 24^3 mass=F, s part of v0: lambda_1(h=1e-7) = %.3e  lambda_1(h=1e-5) = %.3e  ratio %.3g (> 10)  %s'
          % (l7, l5, l7/l5, 'PASS' if good else 'FAIL'))
    l7, l5 = op.lam_first(s, m, eps, 'rho', 1e-7), op.lam_first(s, m, eps, 'rho', 1e-5)
    good = abs(l7/l5 - 1) <= 1e-6
    ok &= good
    print('low-k 24^3 mass=F, v0: lambda_1(h=1e-7) = %.10f  |lambda_1(1e-7)/lambda_1(1e-5) - 1| = %.1e (<= 1e-6)  %s'
          % (l7, abs(l7/l5 - 1), 'PASS' if good else 'FAIL'))
    lm = op.power(s, m, eps, 'rho', niter=50)[0]
    lam_mask = float(np.max(op.mL*sum(np.max(np.asarray(eps[d]))*op.kq2[d] for d in range(3))))
    good = abs(lm/lam_mask - 1) <= 1e-6
    ok &= good
    print('low-k 24^3 mass=F, v0: lambda_50 = %.10f  lambda_50/lambda_mask - 1 = %+.1e (|.| <= 1e-6)  %s'
          % (lm, lm/lam_mask - 1, 'PASS' if good else 'FAIL'))
    return ok


def synth(N=32):
    """F against U/c_s: rho = 1, B = B0 z (v_A = c_s), u_c = U f(x) with f = sin x (smooth) or
    tanh(sin(x)/dl) (a jump of about one cell), c = x (compressive) or y (shear); m masked as the state."""
    g = gen_ic.Grid(N, N, N)
    op = SVV(g, 1.0)
    out = []
    dl = g.Lx/N
    for prof, f in (('sin', np.sin(g.X)), ('step', np.tanh(np.sin(g.X)/dl))):
        for comp in (0, 1):
            for U in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
                u = [np.zeros(op.shape) for _ in range(3)]
                u[comp] = U*f
                st = gen_ic.state_from_real(g, np.zeros(op.shape), tuple(u),
                                            (np.zeros(op.shape), np.zeros(op.shape), np.ones(op.shape)))
                s, m, rho, uu, b = real_state(g, *st)
                r = lam_variant(op, s, m, rho, uu, b, 'global', 'rho', niter=600, tol=1e-6)
                dumax = max(float(np.max(np.abs(uu[comp] - np.roll(uu[comp], -1, 0)))), 0.0)
                out.append(dict(profile=prof, comp='xy'[comp], U=U, F=r['F'], nit=r['nit'], conv=r['conv'],
                                Lambda=r['Lambda'], Lambda_R1=r['Lambda_R1'], du1_max=dumax))
                print('%-4s u_%s U/c_s = %4.2f  F = %8.4f (%3d it, conv %d)  Lambda(R=1, 2) = %8.4f %8.4f  '
                      'max|du| over 1 cell = %.3f'
                      % (prof, 'xy'[comp], U, r['F'], r['nit'], r['conv'], r['Lambda_R1'], r['Lambda'], dumax))
    return out


HOLE_VARIANTS = [('global', 'rho'), ('local', 'rho'), ('global', 'const1'), ('global', 'min1'), ('global', 'min2'),
                 ('global', 'min3'), ('local', 'min1'), ('local', 'min2'), ('local', 'min3')]


def hole_state(g, op, rh, sig=0.6):
    """The hole of hole() on g: (s, m, rho, u, b) as real_state."""
    dl = g.Lx/g.nlx
    c = np.pi  # centre on a grid point
    r2 = ((g.X - c)**2 + (g.Y - c)**2 + (g.Z - c)**2)/(sig*dl)**2
    st = gen_ic.state_from_real(g, np.log(rh)*np.exp(-r2/2), tuple(np.zeros(op.shape) for _ in range(3)),
                                (np.zeros(op.shape), np.zeros(op.shape), np.ones(op.shape)))
    return real_state(g, *st)


def hole(N=32, sig=0.6, rhos=(1e-1, 1e-2, 1e-3, 1e-4), variants=HOLE_VARIANTS):
    """F against the depth rho_h of a hole at rest: ln rho = ln(rho_h) exp(-r^2/(2 (sig dl)^2)), u = 0,
    B = z (v_A = c_s = 1 outside). With sig = 0.6 the 6 neighbours have rho = rho_h^0.25.
    lambda by ARPACK; p = d ln F/d ln rho_h between the last two depths."""
    g = gen_ic.Grid(N, N, N)
    op = SVV(g, 1.0)
    out = []
    for ek, wk in variants:
        Fs = []
        for rh in rhos:
            s, m, rho, u, b = hole_state(g, op, rh, sig)
            r = lam_variant(op, s, m, rho, u, b, ek, wk, method='eig')
            Fs.append(r['F'])
            out.append(dict(eps=ek, weight=wk, rho_h=rh, rho_min=float(rho.min()), F=r['F'], lam=r['lam'], nJv=r['nit'],
                            eps_g=op.eps(rho, u, b), vs_argmax=r['vs_argmax']))
        p = np.log(Fs[-1]/Fs[-2])/np.log(rhos[-1]/rhos[-2])
        lam_p = [o['lam'] for o in out[-2:]]
        q = np.log(lam_p[1]/lam_p[0])/np.log(rhos[-1]/rhos[-2])
        out.append(dict(eps=ek, weight=wk, p_F=p, p_lam=q))
        print('hole  eps %-6s weight %-6s  F(rho_h = %s) = %s  p_F = %+.3f  p_lambda = %+.3f'
              % (ek, wk, ' '.join('%g' % x for x in rhos), ' '.join('%9.3f' % x for x in Fs), p, q))
        sys.stdout.flush()
    return out


def kernel(Ns=(32, 64), radii=(0, 1, 2, 3, 4, 6, 8, 16)):
    """Spread of the L (2/3 mask) and H (outside it) parts of the SVV: K(y) = sum_d (A_d delta)(y)^2 for a delta at
    the origin, A_d = sig_d m_L or sig_d m_H (the operators of tend). The share of sum_y K on the 6 neighbours
    and within the Chebyshev distance r (cells) of the origin."""
    out = []
    for N in Ns:
        g = gen_ic.Grid(N, N, N)
        op = SVV(g, 1.0)
        d0 = np.zeros(op.shape)
        d0[0, 0, 0] = 1.0
        dk = op.fft(d0)
        i = np.arange(N)
        i = np.minimum(i, N - i)
        cheb = np.maximum.reduce(np.meshgrid(i, i, i, indexing='ij'))
        taxi = sum(np.meshgrid(i, i, i, indexing='ij'))
        for part, low in (('L', True), ('H', False)):
            K = sum(op._grad(dk, d, low)**2 for d in range(3))
            sh = {r: float(K[cheb <= r].sum()/K.sum()) for r in radii if r < N//2}
            nb6 = float(K[taxi == 1].sum()/K.sum())
            out.append(dict(N=N, part=part, sum=float(K.sum()), nb6=nb6, share=sh))
            print('kernel %d^3 %s  sum %.4g  6 neighbours %.3f  within r:' % (N, part, K.sum(), nb6),
                  ' '.join('r%d %.3f' % x for x in sh.items()))
    return out


def jlocal(op, s, m, rho, eps, wk, wconst=None, q=10.0, nmax=27):
    """The Jacobian of G (the J v of svv_power) restricted to (s, m_x, m_y, m_z) on the set P of the nmax points
    nearest the min x0 of rho with rho < q rho_min (x0 first), by J_ij = (G_i(U + h e_j) - G_i(U))/h, h = H_FD.
    Returns J_ss and max_c J_mc,mc at x0, G_s(x0) = (d rho/dt)/rho (J_ss = (d rho/dt)'/rho - G_s: the second part
    is the fill by the SVV over rho), max |eigenvalue| of the 4 x 4 block at x0 and of the block on P, and |P|."""
    i0 = np.unravel_index(int(np.argmin(rho)), rho.shape)
    cand = np.argwhere(rho < q*rho[i0])
    n = np.array(rho.shape)
    dd = np.abs(cand - np.array(i0))
    dd = np.minimum(dd, n - dd)
    P = [tuple(int(x) for x in cand[k]) for k in np.argsort((dd**2).sum(1), kind='stable')[:nmax]]
    G0 = op.G(s, m, eps, wk, wconst)
    J = np.zeros((4*len(P), 4*len(P)))
    for jp, p in enumerate(P):
        for c in range(4):
            e = np.zeros(op.shape)
            e[p] = H_FD
            if c == 0:
                Gp = op.G(s + e, m, eps, wk, wconst)
            else:
                Gp = op.G(s, [m[k] + (e if k == c - 1 else 0.0) for k in range(3)], eps, wk, wconst)
            for ip, pp in enumerate(P):
                for c2 in range(4):
                    J[4*ip + c2, 4*jp + c] = (Gp[c2][pp] - G0[c2][pp])/H_FD
    return (J[0, 0], max(J[1, 1], J[2, 2], J[3, 3]), float(G0[0][i0]), float(max(abs(np.linalg.eigvals(J[:4, :4])))),
            float(max(abs(np.linalg.eigvals(J)))), len(P))


def diag(case=None, sub='restart', N=32, sig=0.6, rhos=(1e-2, 1e-4), variants=HOLE_VARIANTS):
    """jlocal at the min of rho over the linear bound of the code, against F (ARPACK): at the centre of the hole
    of hole() on N^3 (case = None) or at the min of a restart."""
    if case is None:
        g = gen_ic.Grid(N, N, N)
        op = SVV(g, 1.0)
        states = [(rh, hole_state(g, op, rh, sig)) for rh in rhos]
    else:
        g, op = from_case(case)
        lk, mk, bk, tt = gen_ic.read_restart(os.path.join(case, sub), g)
        st = real_state(g, lk, mk, bk)
        states = [(float(st[2].min()), st)]
    out = []
    for rh, (s, m, rho, u, b) in states:
        ll = op.lam_lin(op.eps(rho, u, b))
        for ek, wk in variants:
            r = lam_variant(op, s, m, rho, u, b, ek, wk, method='eig')
            wconst = float(rho.mean()) if wk == 'constmean' else float(wk[5:]) if wk.startswith('const') else None
            Js, Jm, Gs, B4, BP, nP = jlocal(op, s, m, rho, op.eps(rho, u, b, local=(ek == 'local')), wk, wconst)
            out.append(dict(rho_min=rh, eps=ek, weight=wk, F=r['F'], J_ss=Js/ll, J_mm=Jm/ll, G_s=Gs/ll, B4=B4/ll,
                            BP=BP/ll, nP=nP))
            print('diag  rho_min %.3e  eps %-6s weight %-6s  F = %10.4g  J_ss = %10.4g  max J_mm = %10.4g  '
                  '-G_s = %10.4g  B4 = %10.4g  B_P = %10.4g (|P| %d)  (over lam_lin; J_ss/F %+.3f, -G_s/J_ss %.3f, '
                  'B4/F %.3f, B_P/F %.3f)' % (rh, ek, wk, r['F'], Js/ll, Jm/ll, -Gs/ll, B4/ll, BP/ll, nP, Js/ll/r['F'],
                                             -Gs/Js, B4/ll/r['F'], BP/ll/r['F']))
            sys.stdout.flush()
    return out


def bound(case=None, sub='restart', N=32, sig=0.6, rhos=(1e-2, 1e-4), variants=HOLE_VARIANTS):
    """Lower bounds of lambda at the min x0 of rho (states as diag) from the frozen SVV, dU/dt = -A* K A H dU on
    U = (rho, m), H = dw/dU. -A* K A H is similar to the symmetric -H^1/2 A* K A H^1/2 (and through ds = d rho/rho
    to the frozen operator on (s, m); the Jacobian of G adds the derivatives of the coefficients and
    -(d rho/dt)/rho on the diagonal), so for dw = H dU the Rayleigh quotient sum dw.K dw/(dU.H dU) bounds its
    largest |eigenvalue| from below. With e = delta(x - x0) and P_d = A_L,d e, P_H,d = A_H,d e:
    - Q: dw = (e, 0). Q = cs^2 S/rho(x0), S = sum_d sum_y [a_d P_d^2 + kap_d W P_H,d^2](y).
    - Q4: dw = c e for any 4-vector c (the 4 point perturbations of (rho, m) at x0): the largest eigenvalue of
      N c = q H^-1(x0) c, N_ij = sum_d sum_y [a_d P_d^2 tau_i tau_j + b_d P_d^2 delta_ij (i, j > 0)
      + kap_d W P_H,d^2 delta_i0 delta_j0], tau = (1, u(y)).
    Since b_d = cs^2 a_d, Q4/Q is (up to the kap part) the largest eigenvalue of [[1, d^T], [d, I + E]],
    d = (<u> - u(x0))/cs, E = <(u - u(x0))(u - u(x0))^T>/cs^2, <.> the mean over y with the weights sum_d a_d P_d^2:
    1 + lambda_max(E) <= Q4/Q <= 2 + Mrel^2 (Mrel = tr(E)^1/2), so Q4 > Q needs the velocity on the kernel to differ
    from that at x0.
    The fill R = (d rho/dt)_svv(x0) = -sum_d sum_y [P_d a t + P_H,d kap W g_h](y) has the coefficients and the
    kernel of S, so mu = R/(cs^2 S) is the mean of -t/(cs^2 P_d) and -g_h/(cs^2 P_H,d) with the weights of S,
    and Q = R/(mu rho(x0)). The Jacobian of G has -(d rho/dt)/rho = -R/rho(x0) on the (s, s) diagonal at x0."""
    if case is None:
        g = gen_ic.Grid(N, N, N)
        op = SVV(g, 1.0)
        states = [(rh, hole_state(g, op, rh, sig)) for rh in rhos]
    else:
        g, op = from_case(case)
        lk, mk, bk, tt = gen_ic.read_restart(os.path.join(case, sub), g)
        st = real_state(g, lk, mk, bk)
        states = [(float(st[2].min()), st)]
    cs2 = op.cs2
    out = []
    for rh, (s, m, rho, u, b) in states:
        ll = op.lam_lin(op.eps(rho, u, b))
        i0 = np.unravel_index(int(np.argmin(rho)), rho.shape)
        r0 = float(rho[i0])
        u0 = np.array([float(x[i0]) for x in u])
        Hinv = r0*np.block([[np.ones((1, 1))/cs2, u0[None, :]/cs2], [u0[:, None]/cs2, np.eye(3) + np.outer(u0, u0)/cs2]])
        e0 = np.zeros(op.shape)
        e0[i0] = 1.0
        e0k = op.fft(e0)
        P2 = [op._grad(e0k, d, True)**2 for d in range(3)]
        PH2 = [op._grad(e0k, d, False)**2 for d in range(3)]
        umax2 = float(np.max(sum(x**2 for x in u)))
        tau = [np.ones(op.shape)] + list(u)
        for ek, wk in variants:
            eps = op.eps(rho, u, b, local=(ek == 'local'))
            W = (np.full(op.shape, float(rho.mean()) if wk == 'constmean' else float(wk[5:])) if wk.startswith('const')
                 else op.weight(rho, wk))
            Nm = np.zeros((4, 4))
            wL = np.zeros(op.shape)
            for d in range(3):
                if np.all(np.asarray(eps[d]) == 0):
                    continue
                aP2 = eps[d]*W/cs2*P2[d]
                wL += aP2
                for i in range(4):
                    for j in range(i, 4):
                        Nm[i, j] += float(np.sum(aP2*tau[i]*tau[j]))
                bP2 = float(np.sum(eps[d]*W*P2[d]))
                for i in range(1, 4):
                    Nm[i, i] += bP2
                if op.filter_mom:
                    Nm[0, 0] += float(np.sum(eps[d]/(cs2 + umax2)*W*PH2[d]))
            Nm = np.triu(Nm) + np.triu(Nm, 1).T
            S = Nm[0, 0]
            Q = cs2*S/r0
            Q4 = float(np.max(np.linalg.eigvals(np.linalg.solve(Hinv, Nm)).real))
            Mrel = float(np.sqrt(np.sum(wL*sum((x - x0)**2 for x, x0 in zip(u, u0)))/np.sum(wL)/cs2))
            R = float(op.tend(s, u, W, eps)[0][i0])
            out.append(dict(rho_min=rh, eps=ek, weight=wk, Q=Q/ll, Q4=Q4/ll, R=R, R_rho=R/r0/ll, mu=R/(cs2*S),
                            M0=float(np.sqrt(u0 @ u0/cs2)), Mrel=Mrel))
            print('bound  rho_min %.3e  |u|/cs %.2f  eps %-6s weight %-6s  Q = %10.4g  Q4 = %10.4g  R/rho = %+10.4g '
                  '(over lam_lin)  Mrel %.2f  R = %+10.4g  mu = %+.4g'
                  % (rh, np.sqrt(u0 @ u0/cs2), ek, wk, Q/ll, Q4/ll, R/r0/ll, Mrel, R, R/(cs2*S)))
            sys.stdout.flush()
    return out


def state(case, sub='restart', jout=None, variants=VARIANTS, method='power'):
    g, op = from_case(case)
    lk, mk, bk, tt = gen_ic.read_restart(os.path.join(case, sub), g)
    s, m, rho, u, b = real_state(g, lk, mk, bk)
    i0 = np.unravel_index(int(np.argmin(rho)), rho.shape)
    res = dict(case=case, sub=sub, t=tt, N=list(g.nl), cs2=op.cs2, rho_min=float(rho[i0]), rho_mean=float(rho.mean()),
               imin=[int(a) for a in i0])
    # neighbours of the min along x, z, y (array axes 0, 1, 2)
    nb = {}
    for ax, name in ((0, 'x'), (1, 'z'), (2, 'y')):
        vals = []
        for o in (-1, 1):
            j = list(i0)
            j[ax] = (j[ax] + o) % rho.shape[ax]
            vals.append(float(rho[tuple(j)]))
        nb[name] = vals
    res['rho_nb'] = nb
    eg = op.eps(rho, u, b)
    el = op.eps(rho, u, b, local=True)
    res['eps_global'] = eg
    res['eps_local_at_min'] = [float(np.asarray(e)[i0]) for e in el]
    nbs = [tuple((np.array(i0) + o*np.eye(3, dtype=int)[ax]) % np.array(rho.shape)) for ax in range(3) for o in (-1, 1)]
    res['eps_local_max_nb'] = [float(max(np.asarray(e)[j] for j in nbs)) for e in el]
    res['lam_lin'] = op.lam_lin(eg)
    c = op.cmax(rho, u, b)
    res['cmax_at_min'] = [float(x[i0]) for x in c]
    res['vA_at_min'] = float(np.sqrt(sum(x[i0]**2 for x in b)/rho[i0]))
    # d rho/dt at the min: -div_N (M m) and the SVV (A_L and A_H parts), with the eps and weight of the code
    divm = g.ifft(1j*(g.kxd*mk[0] + g.kyd*mk[1] + g.kzd*mk[2])*op.mL)
    rdot, mdot, dsum, rL, rH = op.tend(s, u, rho, eg, parts=True)
    res['budget_at_min'] = dict(adv=float(-divm[i0]), svv=float(rdot[i0]), svv_L=float(rL[i0]), svv_H=float(rH[i0]))
    res['D_svv'] = dsum*g.Lx*g.Ly*g.Lz/g.ntot
    # the same with the variants of the weight (one evaluation)
    bud = {}
    for ek, wk in variants:
        if wk == 'constmean':
            W = np.full(op.shape, float(rho.mean()))
        else:
            W = op.weight(rho, wk)
        r2, _, _, l2, h2 = op.tend(s, u, W, op.eps(rho, u, b, local=(ek == 'local')), parts=True)
        bud['%s/%s' % (ek, wk)] = dict(svv=float(r2[i0]), svv_L=float(l2[i0]), svv_H=float(h2[i0]))
    res['budget_svv_variants'] = bud
    print('%s/%s t = %.6f  rho_min %.3e at %s  neighbours %s' % (case, sub, tt, rho[i0], tuple(i0), nb))
    print('  eps_global %s  lam_lin %.4e  v_A at min %.1f  c_max at min %s'
          % (['%.4g' % e for e in eg], res['lam_lin'], res['vA_at_min'], ['%.4g' % x for x in res['cmax_at_min']]))
    print('  d rho/dt at min: -div m %+.4e  SVV %+.4e (A_L %+.4e, A_H %+.4e)'
          % tuple(res['budget_at_min'][k] for k in ('adv', 'svv', 'svv_L', 'svv_H')))
    res['variants'] = []
    for ek, wk in variants:
        r = lam_variant(op, s, m, rho, u, b, ek, wk, method=method)
        r['F_rho_min'] = r['F']*res['rho_min']
        res['variants'].append(r)
        print('  eps %-6s weight %-9s F = %10.4f  F rho_min = %.3e  (%3d it, conv %d, %.0f s)  Lambda(R=1, 2) %10.4f '
              '%10.4f  ratio %9.3f  u-factor %6.3f  W(min) %.3e  v_s max at %s, v_m max at %s'
              % (ek, wk, r['F'], r['F_rho_min'], r['nit'], r['conv'], r['sec'], r['Lambda_R1'], r['Lambda'], r['ratio_max'],
                 r['ufac_max'], r['W_at_min'], r['vs_argmax'], r['vm_argmax']))
        sys.stdout.flush()
    if jout:
        with open(jout, 'w') as f:
            json.dump(res, f, indent=1)
    return res


def lam(case, ek, wk, niter=None, sub='restart'):
    g, op = from_case(case)
    d = read_nml(os.path.join(case, 'calliope.in'))
    niter = int(niter) if niter else int(d.get('svv_pi_init', '50'))
    lk, mk, bk, tt = gen_ic.read_restart(os.path.join(case, sub), g)
    s, m, rho, u, b = real_state(g, lk, mk, bk)
    eps = op.eps(rho, u, b, local=(ek == 'local'))
    wconst = np.full(op.shape, float(wk[5:])) if wk.startswith('const') else None
    lam_, nit, conv, v, hist = op.power(s, m, eps, wk, niter, 0.0, wconst=wconst)
    eg = op.eps(rho, u, b)
    print('%s/%s t = %.6f eps %s weight %s: lambda after %d iterations %.12e  F %.8f  eps_global %s'
          % (case, sub, tt, ek, wk, nit, lam_, lam_/op.lam_lin(eg), ' '.join('%.12e' % e for e in eg)))
    return lam_


if __name__ == '__main__':
    a = sys.argv[1:]
    if not a or a[0] == 'selftest':
        sys.exit(0 if selftest() else 1)
    elif a[0] == 'synth':
        r = synth(int(a[1]) if len(a) > 1 else 32)
        if '--json' in a:
            with open(a[a.index('--json') + 1], 'w') as f:
                json.dump(r, f, indent=1)
    elif a[0] == 'hole':
        r = hole(int(a[1]) if len(a) > 1 and a[1][0] != '-' else 32)
        if '--json' in a:
            with open(a[a.index('--json') + 1], 'w') as f:
                json.dump(r, f, indent=1)
    elif a[0] == 'kernel':
        kernel(tuple(int(x) for x in a[1:]) if len(a) > 1 else (32, 64))
    elif a[0] in ('diag', 'bound'):
        fn = diag if a[0] == 'diag' else bound
        rest = [x for x in a[1:] if x != '--json' and (a.index(x) == 0 or a[a.index(x) - 1] != '--json')]
        if rest and not rest[0].isdigit():
            r = fn(rest[0], rest[1] if len(rest) > 1 else 'restart')
        else:
            r = fn(N=int(rest[0]) if rest else 32)
        if '--json' in a:
            with open(a[a.index('--json') + 1], 'w') as f:
                json.dump(r, f, indent=1)
    elif a[0] == 'state':
        jout = a[a.index('--json') + 1] if '--json' in a else None
        meth = 'eig' if '--eig' in a else 'power'
        rest = [x for x in a[1:] if x not in ('--json', '--eig') and x != jout]
        state(rest[0], rest[1] if len(rest) > 1 else 'restart', jout, method=meth)
    elif a[0] == 'lam':
        lam(*a[1:])
