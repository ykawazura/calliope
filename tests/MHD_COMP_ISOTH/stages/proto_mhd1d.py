#!/usr/bin/env python3
"""1D numpy prototype of the MHD_COMP_ISOTH discretization with B (formulation C), for Stage 3.

It mirrors advance.F90 / shock_viscosity.F90 for a run along x with nly = nlz = 2 (proto1d.py is the B = 0 case):
    state s^ = F[ln rho], m^ = F[rho u] (3 components), b^ = F[B_y, B_z] (rfft coefficients); B_x is a constant
    (div_free removes the k != 0 modes of B_x); m and B are filtered by the 2/3 mask after every stage and at the
    start (filter_mom, filter_mag), s is not (filter_lnrho = F),
    ds/dt = -(d_N m_x)/rho,  d_N = F^-1 i k_d M F (k_d = 0 at Nyquist, M = the m filter),
    dm/dt = -d [rho u_x u - B_x B + (cs^2 rho + |B|^2/2 - q) e_x],  q = rho zeta D,  D = d_N u_x,  d = F^-1 i k_d F,
    dB_y/dt = d(B_x u_y - B_y u_x),  dB_z/dt = d(B_x u_z - B_z u_x)   (get_ext_terms: curl (b x u)),
    zeta = C Delta^2 S[max(0, -D)] (proto1d.zeta_of),
    SVV: eps = svv_eps0 dl max(|u_x| + c_f), c_f the fast speed along x, fixed at stage 1 of each step;
      svv_form = 'entropy': the ES-SVV of proto1d in (rho, m) with w = (cs^2 s - |u|^2/2, u) and the three
        components of u, kappa_h = eps cs^2/(cs^2 + max |u|^2); lam_svv by the power method on (s, m_x, m_y, m_z);
      'linear': lambda = eps k^2 Q on m and on rho, as in proto1d;
      both: dB^/dt -= lambda B^ (the linear SVV of B is used with either form),
    dt_adv = cfl dl/max(|u_x| + c_f); dt_visc as in advance.F90 with the B term eps lam_b,
    time stepping: eSSPIFRK3 of advance_common.F90 without linear terms.
The shock resistivity (C_eta of shock_parameters in the Fortran code, Stage 3; off by default):
    eta = C_eta Delta^2 S[max(0, -D)] (the zeta detector), dB/dt|_eta = -curl_N(eta J), J = curl_N B,
    i.e. dB_t/dt|_eta = d(eta d B_t); after the B filter dE_B/dt|_eta = -sum eta |J|^2 dV, and eta_max kd2_max
    is added to the B term of dt_visc.
The defaults are those of the Fortran code (params.F90).
track_budget = True adds budget() of the stage-1 state to each hist entry; 'stages' also adds rk_<key>, the
stage-weighted rates of the step (dE = dt rk_eps up to the time-integration error).

    from proto_mhd1d import MHD1D
    p = MHD1D(512, lx=1/(2*math.pi), bx=1.0); p.set_prim(rho, ux, uy, uz, by, bz); p.run_until(0.1, 1e-3)
"""
import math
import numpy as np
from proto1d import Isoth1D


class MHD1D(Isoth1D):
    def __init__(self, N, bx=0.0, filter_b=True, C_eta=0.0, **kw):
        for k in ('formulation', 'visc_form', 'deriv_mask', 'filter_rhs', 'pad'):
            assert k not in kw, k          # the variants of proto1d are not carried over
        super().__init__(N, **kw)
        self.bx = bx
        self.Fb = self.mask if filter_b else np.ones_like(self.mask)
        self.lam_b = np.max(self.kq2 * (self.Fb > 0))
        self.C_eta = C_eta
        self.nk = N // 2 + 1

    # ------------------------------------------------------------ initial conditions
    def set_su(self, s, u):
        """proto1d's set_su (used by riemann and wave1d of proto1d): B_t = 0, u_t = 0."""
        z = np.zeros(self.Ng)
        self.set_s_prim(s, u, z, z, z, z)

    def set_s_prim(self, s, ux, uy, uz, by, bz):
        z = np.zeros(self.Ng)
        rho = np.exp(s + z)
        self.sh = self.fwd(s + z) * self.Fs
        self.mh = np.array([self.fwd(rho * (u + z)) for u in (ux, uy, uz)]) * self.Fm
        self.bh = np.array([self.fwd(b + z) for b in (by, bz)]) * self.Fb

    def set_prim(self, rho, ux=0.0, uy=0.0, uz=0.0, by=0.0, bz=0.0):
        """State from primitive fields on the grid (arrays or scalars), filtered as at the start of the Fortran run."""
        self.set_s_prim(np.log(rho + np.zeros(self.Ng)), ux, uy, uz, by, bz)

    def set_W(self, W):
        """W = (rho, u_x, u_y, u_z, B_y, B_z) on the grid (the order of mhd_exact)."""
        self.set_prim(*W)

    def riemann_mhd(self, WL, WR, width=1.0):
        """init_riemann of fields.F90 with B: every primitive variable mixed by the tanh steps at L/4 and 3L/4."""
        x, L = self.x, self.L

        def step(x0):
            if width > 0.0:
                return 0.5 * (1.0 + np.tanh((x - x0) / (width * self.dl)))
            return np.where(x > x0, 1.0, np.where(x == x0, 0.5, 0.0))
        h = step(0.25 * L) - step(0.75 * L)
        self.set_W([(1.0 - h) * a + h * b for a, b in zip(WL, WR)])

    def prim(self):
        """(rho, u_x, u_y, u_z, B_y, B_z) of the current state on the grid."""
        f = self.fields(self.sh, self.mh, self.bh)
        return np.array([f['rho'], *f['u'], *f['B']])

    # ------------------------------------------------------------ operators
    def detector(self, D):
        """S[max(0, -D)] of zeta_of, without the factor C Delta^2."""
        z = np.maximum(0.0, -D)
        for _ in range(self.nmax):
            z = np.maximum(np.maximum(np.roll(z, 1), z), np.roll(z, -1))
        for _ in range(self.nsmooth):
            z = 0.5 * (z + 0.5 * (np.roll(z, 1) + np.roll(z, -1)))
        return z

    def cfast(self, rho, B):
        """Fast speed along x (advance.F90, maxs(6))."""
        c2 = self.cs ** 2
        va2 = (self.bx ** 2 + B[0] ** 2 + B[1] ** 2) / rho
        return np.sqrt(0.5 * (c2 + va2 + np.sqrt(np.maximum(0.0, (c2 + va2) ** 2 - 4.0 * c2 * self.bx ** 2 / rho))))

    def fields(self, sh, mh, bh):
        s = self.bwd(sh)
        rho = np.exp(s)
        m = np.array([self.bwd(a) for a in mh])
        u = m / rho
        uh = np.array([self.fwd(a) for a in u])
        D = self.ddx(uh[0], self.Md)
        z = self.detector(D) if self.C > 0.0 or self.C_eta > 0.0 else np.zeros_like(D)
        zeta = self.C * self.Delta ** 2 * z
        eta = self.C_eta * self.Delta ** 2 * z
        B = np.array([self.bwd(a) for a in bh])
        return dict(s=s, rho=rho, m=m, u=u, uh=uh, D=D, zeta=zeta, eta=eta, B=B, q=rho * zeta * D)

    def vmax(self, rho, u, B):
        return np.max(np.abs(u[0]) + self.cfast(rho, B))

    def svv_eps(self, v):
        return self.svv_eps0 * self.dl * v if self.svv_eps0 > 0.0 else 0.0

    def svv_entropy3(self, s, rho, u, uh):
        """ES-SVV with the three components of u (get_svv_es): (d rho/dt, dm^/dt (3), D_svv density)."""
        c2 = self.cs ** 2
        u2 = np.sum(u * u, axis=0)
        wh = self.fwd(c2 * s - 0.5 * u2)
        g_r = self.bwd(self.sigL * wh)
        g_m = np.array([self.bwd(self.sigL * a) for a in uh])
        g_h = self.bwd(self.sigH * wh)
        f_r = self.eps * rho / c2 * (g_r + np.sum(u * g_m, axis=0))
        f_m = u * f_r + self.eps * rho * g_m
        f_h = self.eps / (c2 + np.max(u2)) * rho * g_h
        rdot = self.bwd(-np.conj(self.sigL) * self.fwd(f_r) - np.conj(self.sigH) * self.fwd(f_h))
        mdh = np.array([-np.conj(self.sigL) * self.fwd(a) for a in f_m])
        return rdot, mdh, f_r * g_r + np.sum(f_m * g_m, axis=0) + f_h * g_h

    def svv_G(self, s, m):
        """SVV tendencies (ds/dt, dm/dt (3)) in real space for the frozen eps (power method)."""
        rho = np.exp(s)
        u = m / rho
        rdot, mdh, _ = self.svv_entropy3(s, rho, u, np.array([self.fwd(a) for a in u]))
        return rdot / rho, np.array([self.bwd(a) for a in mdh])

    def svv_power(self, s, m, niter, tol=0.0, h=1e-7):
        """proto1d.svv_power on (s, m_x, m_y, m_z); |v|^2 = mean(v_s^2) + sum_c mean(v_mc^2)."""
        if self.pi_v is None:
            i = np.arange(self.Ng)
            vs = (-1.0) ** i * (1.0 + 0.3 * np.cos(2.0 * math.pi * i / self.Ng + 1.0))
            self.pi_v = (vs, np.zeros((3, self.Ng)))
        G0s, G0m = self.svv_G(s, m)
        vs, vm = self.pi_v
        lam_old = None
        for it in range(1, niter + 1):
            nv = math.sqrt(np.mean(vs * vs) + np.sum(np.mean(vm * vm, axis=1)))
            vs, vm = vs / nv, vm / nv
            Gs, Gm = self.svv_G(s + h * vs, m + h * vm)
            vs, vm = (Gs - G0s) / h, (Gm - G0m) / h
            lam = math.sqrt(np.mean(vs * vs) + np.sum(np.mean(vm * vm, axis=1)))
            done = tol > 0.0 and lam_old is not None and abs(lam - lam_old) <= tol * lam
            lam_old = lam
            if done:
                break
        self.lam_svv, self.pi_v = lam, (vs, vm)
        return it

    def prepare(self):
        """eps and lam_svv of the stage-1 state (proto1d.prepare)."""
        if self.prepared == self.step:
            return
        s0 = self.bwd(self.sh)
        rho0 = np.exp(s0)
        m0 = np.array([self.bwd(a) for a in self.mh])
        B0 = np.array([self.bwd(a) for a in self.bh])
        self.eps = self.svv_eps(self.vmax(rho0, m0 / rho0, B0))
        self.pi_conv = None
        if self.svv_form == 'entropy' and self.eps > 0.0:
            if self.pi_v is None:
                self.svv_power(s0, m0, self.svv_pi_init)
            elif self.svv_pi_check > 0 and self.step % self.svv_pi_check == 0:
                n = self.svv_power(s0, m0, 200, tol=1e-3)
                self.pi_conv = (self.lam_svv, n)
            else:
                self.svv_power(s0, m0, 1)
        self.prepared = self.step

    # ------------------------------------------------------------ tendencies
    def rhs(self, sh, mh, bh):
        """(ds^/dt, dm^/dt, dB^/dt, fields) before the stage filter, with the fluxes of get_nonlinear_terms."""
        f = self.fields(sh, mh, bh)
        s, rho, u, B, q = f['s'], f['rho'], f['u'], f['B'], f['q']
        bx, ikd = self.bx, 1j * self.kd
        f['divm'] = self.ddx(mh[0], self.Fm)
        fs = self.fwd(-f['divm'] / rho)
        ptot = self.cs ** 2 * rho + 0.5 * (bx * bx + B[0] ** 2 + B[1] ** 2) - q
        fm = np.array([-ikd * self.fwd(rho * u[0] * u[0] - bx * bx + ptot),
                       -ikd * self.fwd(rho * u[0] * u[1] - bx * B[0]),
                       -ikd * self.fwd(rho * u[0] * u[2] - bx * B[1])])
        fb = np.array([ikd * self.fwd(bx * u[1] - B[0] * u[0]),
                       ikd * self.fwd(bx * u[2] - B[1] * u[0])])
        f['d_svv'] = None
        if self.eps > 0.0:
            lam = self.eps * self.kq2
            if self.svv_form == 'linear':
                fm = fm - lam * mh
                fs = fs + self.fwd(self.bwd(-lam * self.fwd(rho)) / rho)
            else:
                rdot, mdh, f['d_svv'] = self.svv_entropy3(s, rho, u, f['uh'])
                fs = fs + self.fwd(rdot / rho)
                fm = fm + mdh
            fb = fb - lam * bh
        if self.C_eta > 0.0:
            g = np.array([self.bwd(ikd * a) for a in bh])
            f['J'] = g
            fb = fb + np.array([ikd * self.fwd(f['eta'] * a) for a in g])
        return fs, fm, fb, f

    def parts(self, sh, mh, bh):
        """Tendencies split into ideal, visc (shock viscosity), svv (rho, m), svv_b (B) and eta, before the filter."""
        f = self.fields(sh, mh, bh)
        s, rho, u, B, q = f['s'], f['rho'], f['u'], f['B'], f['q']
        bx, ikd, c2 = self.bx, 1j * self.kd, self.cs ** 2
        zs, zm, zb = np.zeros(self.nk, complex), np.zeros((3, self.nk), complex), np.zeros((2, self.nk), complex)
        divm = self.ddx(mh[0], self.Fm)
        pt = c2 * rho + 0.5 * (bx * bx + B[0] ** 2 + B[1] ** 2)
        ideal = (self.fwd(-divm / rho),
                 np.array([-ikd * self.fwd(rho * u[0] * u[0] - bx * bx + pt),
                           -ikd * self.fwd(rho * u[0] * u[1] - bx * B[0]),
                           -ikd * self.fwd(rho * u[0] * u[2] - bx * B[1])]),
                 np.array([ikd * self.fwd(bx * u[1] - B[0] * u[0]), ikd * self.fwd(bx * u[2] - B[1] * u[0])]))
        visc = (zs, np.array([ikd * self.fwd(q), zs, zs]), zb)
        svv, svv_b, eta = (zs, zm, zb), (zs, zm, zb), (zs, zm, zb)
        f['d_svv'], f['J'] = None, None
        if self.eps > 0.0:
            lam = self.eps * self.kq2
            if self.svv_form == 'linear':
                svv = (self.fwd(self.bwd(-lam * self.fwd(rho)) / rho), -lam * mh, zb)
            else:
                rdot, mdh, f['d_svv'] = self.svv_entropy3(s, rho, u, f['uh'])
                svv = (self.fwd(rdot / rho), mdh, zb)
            svv_b = (zs, zm, -lam * bh)
        if self.C_eta > 0.0:
            g = np.array([self.bwd(ikd * a) for a in bh])
            f['J'] = g
            eta = (zs, zm, np.array([ikd * self.fwd(f['eta'] * a) for a in g]))
        f['divm'] = divm
        return dict(ideal=ideal, visc=visc, svv=svv, svv_b=svv_b, eta=eta), f

    def budget(self, sh=None, mh=None, bh=None):
        """Semi-discrete budget at (sh, mh, bh) (default: the current state), as proto1d.budget with B:
        for each part rho dot = rho F^-1(s dot), m dot, B dot on the grid; dM, dPx..dPz, dBy, dBz (sums times dV);
        eps = sum(w_rho rho dot + u.m dot + B.B dot) dV, w_rho = cs^2 s - |u|^2/2.
        Keys <part>_<rate> after the stage filter and <part>_<rate>_raw before it; visc = sum rho zeta D^2 dV,
        D_svv = sum d_svv dV, eta_diss = sum eta |J|^2 dV; R_zeta = visc_eps + visc, R_svv = svv_eps + D_svv,
        R_eta = eta_eps + eta_diss (all of them should be round-off); A_M, A_P = sum|ideal rho dot|, sum|ideal m dot|."""
        sh = self.sh if sh is None else sh
        mh = self.mh if mh is None else mh
        bh = self.bh if bh is None else bh
        P, f = self.parts(sh, mh, bh)
        s, rho, u, B = f['s'], f['rho'], f['u'], f['B']
        dV = self.L / self.Ng
        w = self.cs ** 2 * s - 0.5 * np.sum(u * u, axis=0)
        names = ('dM', 'dPx', 'dPy', 'dPz', 'dBy', 'dBz', 'eps')
        out = {}
        for tag, Fs, Fm, Fb in (('_raw', 1.0, 1.0, 1.0), ('', self.Fs, self.Fm, self.Fb)):
            for k, (a, b, c) in P.items():
                rdot = rho * self.bwd(a * Fs)
                mdot = np.array([self.bwd(x * Fm) for x in b])
                bdot = np.array([self.bwd(x * Fb) for x in c])
                e = w * rdot + np.sum(u * mdot, axis=0) + np.sum(B * bdot, axis=0)
                vals = (np.sum(rdot), *np.sum(mdot, axis=1), *np.sum(bdot, axis=1), np.sum(e))
                for r, v in zip(names, vals):
                    out[k + '_' + r + tag] = v * dV
                if tag == '_raw' and k == 'ideal':
                    out['A_M'] = np.sum(np.abs(rdot)) * dV
                    out['A_P'] = np.sum(np.abs(mdot)) * dV
                    out['A_B'] = np.sum(np.abs(bdot)) * dV
            for r in names:
                out[r + tag] = sum(out[k + '_' + r + tag] for k in P)
        out['visc'] = np.sum(rho * f['zeta'] * f['D'] ** 2) * dV
        out['D_svv'] = np.sum(f['d_svv']) * dV if f['d_svv'] is not None else 0.0
        out['eta_diss'] = np.sum(f['eta'] * f['J'] ** 2) * dV if f['J'] is not None else 0.0
        for tag in ('_raw', ''):
            out['R_zeta' + tag] = out['visc_eps' + tag] + out['visc']
            out['R_svv' + tag] = out['svv_eps' + tag] + out['D_svv'] if f['d_svv'] is not None else 0.0
            out['R_eta' + tag] = out['eta_eps' + tag] + out['eta_diss']
        return out

    def monitor(self, f):
        rho, u, B = f['rho'], f['u'], f['B']
        dV = self.L / self.Ng
        eta_rate = f['eta'].max() * self.kd2_max
        if self.svv_form == 'entropy':
            rate = max(f['zeta'].max() * self.kd2_max / self.cfl_visc + self.lam_svv / self.cfl_svv,
                       (self.eps * self.lam_b + eta_rate) / self.cfl_visc)
            dt_visc = 1.0 / rate if rate > 0 else np.inf
        else:
            rate = max(f['zeta'].max() * self.kd2_max + self.eps * self.lam_m, self.eps * self.lam_s,
                       self.eps * self.lam_b + eta_rate)
            dt_visc = self.cfl_visc / rate if rate > 0 else np.inf
        E_kin = np.sum(0.5 * rho * np.sum(u * u, axis=0)) * dV
        E_mag = np.sum(0.5 * (self.bx ** 2 + B[0] ** 2 + B[1] ** 2)) * dV
        E_int = np.sum(self.cs ** 2 * rho * (f['s'] - 1.0)) * dV
        P = np.sum(rho * u, axis=1) * dV
        return dict(t=self.t, M=np.sum(rho) * dV, Px=P[0], Py=P[1], Pz=P[2],
                    By=np.sum(B[0]) * dV, Bz=np.sum(B[1]) * dV,
                    E=E_kin + E_mag + E_int, E_kin=E_kin, E_mag=E_mag,
                    rho_min=rho.min(), rho_max=rho.max(), zeta_max=f['zeta'].max(), eta_max=f['eta'].max(),
                    D_min=f['D'].min(), visc=np.sum(rho * f['zeta'] * f['D'] ** 2) * dV, eps_svv=self.eps,
                    lam_svv=self.lam_svv, D_svv=np.sum(f['d_svv']) * dV if f['d_svv'] is not None else 0.0,
                    dt_adv=self.cfl * self.dl / self.vmax(rho, u, B), dt_visc=dt_visc)

    # ------------------------------------------------------------ time stepping
    def advance(self, dt):
        self.prepare()
        s0, m0, b0 = self.sh, self.mh, self.bh
        f0s, f0m, f0b, f = self.rhs(s0, m0, b0)
        mon = self.monitor(f)
        if self.track_budget:
            mon.update(self.budget(s0, m0, b0))
        if self.pi_conv is not None:
            mon['pi_conv'] = self.pi_conv
            if dt * self.pi_conv[0] > 2.51:
                raise FloatingPointError('svv_power: dt*lam_svv = %.3f > 2.51 at t = %.6g' % (dt * self.pi_conv[0], self.t))
        Fs, Fm, Fb = self.Fs, self.Fm, self.Fb
        s1 = (s0 + 2.0 / 3.0 * dt * f0s) * Fs
        m1 = (m0 + 2.0 / 3.0 * dt * f0m) * Fm
        b1 = (b0 + 2.0 / 3.0 * dt * f0b) * Fb
        f1s, f1m, f1b, _ = self.rhs(s1, m1, b1)
        s2 = (2.0 / 3.0 * s0 + 1.0 / 3.0 * s1 + 4.0 / 9.0 * dt * f1s) * Fs
        m2 = (2.0 / 3.0 * m0 + 1.0 / 3.0 * m1 + 4.0 / 9.0 * dt * f1m) * Fm
        b2 = (2.0 / 3.0 * b0 + 1.0 / 3.0 * b1 + 4.0 / 9.0 * dt * f1b) * Fb
        f2s, f2m, f2b, _ = self.rhs(s2, m2, b2)
        if self.track_budget == 'stages':
            # the rates of the step, sum_i b_i X_i over the stages with b = (1/4, 3/16, 9/16) (budget.dat)
            B1, B2 = self.budget(s1, m1, b1), self.budget(s2, m2, b2)
            mon.update({'rk_' + k: 0.25 * mon[k] + 0.1875 * B1[k] + 0.5625 * B2[k] for k in B1})
        self.sh = (37.0 / 64.0 * s0 + 5.0 / 32.0 * dt * f0s + 27.0 / 64.0 * s2 + 9.0 / 16.0 * dt * f2s) * Fs
        self.mh = (37.0 / 64.0 * m0 + 5.0 / 32.0 * dt * f0m + 27.0 / 64.0 * m2 + 9.0 / 16.0 * dt * f2m) * Fm
        self.bh = (37.0 / 64.0 * b0 + 5.0 / 32.0 * dt * f0b + 27.0 / 64.0 * b2 + 9.0 / 16.0 * dt * f2b) * Fb
        self.t += dt
        self.step += 1
        return mon, f

    def finite(self):
        return bool(np.all(np.isfinite(self.sh)) and np.all(np.isfinite(self.mh)) and np.all(np.isfinite(self.bh)))

    def dt_limit(self):
        self.prepare()
        mon = self.monitor(self.rhs(self.sh, self.mh, self.bh)[3])
        return min(mon['dt_adv'], mon['dt_visc'])

    def run_until(self, t_end, dt_max, every=0, callback=None, maxstep=10 ** 7):
        """Adaptive dt = min(dt_max, dt_limit); returns 'ok', 'svv_unstable' or 'nonfinite'."""
        while self.t < t_end * (1.0 - 1e-12) and self.step < maxstep:
            dt = min(dt_max, self.dt_limit(), t_end - self.t)
            try:
                mon, f = self.advance(dt)
            except FloatingPointError:
                return 'svv_unstable'
            mon['dt'] = dt
            self.hist.append(mon)
            if not self.finite():
                return 'nonfinite'
            if callback is not None and every > 0 and self.step % every == 0:
                callback(self, f)
        return 'ok'

    def run(self, dt, nstep, every=0, callback=None, stop_on_cfl=True):
        for n in range(nstep):
            try:
                mon, f = self.advance(dt)
            except FloatingPointError:
                return 'svv_unstable'
            mon['dt'] = dt
            self.hist.append(mon)
            if not self.finite():
                return 'nonfinite'
            if stop_on_cfl and dt > min(mon['dt_adv'], mon['dt_visc']):
                return 'cfl'
            if callback is not None and every > 0 and n % every == 0:
                callback(self, f)
        return 'ok'

    def state(self):
        f = self.fields(self.sh, self.mh, self.bh)
        return dict(x=self.x, s=f['s'], rho=f['rho'], u=f['u'], B=f['B'], D=f['D'], zeta=f['zeta'], eta=f['eta'])
