#!/usr/bin/env python3
"""1D numpy prototype of the MHD_COMP_ISOTH discretization (B = 0, formulation C).

It mirrors advance.F90 / shock_viscosity.F90 for a run along x with nly = nlz = 2:
    state s^ = F[ln rho], m^ (rfft coefficients), filtered by the 2/3 mask after every stage,
    ds/dt = -(div_N m)/rho,  dm/dt = -div_N [rho u u + cs^2 rho - q],  q = rho zeta D,
    D = div_N u,  div_N = F^-1 i k_d M F  (k_d = 0 at Nyquist, M = 2/3 mask),
    zeta = C Delta^2 S[max(0, -D)],  S = (3-point max)^nmax then ([1,2,1]/4)^nsmooth,
    SVV (shock_viscosity.F90), eps = svv_eps0 dl max(|u| + cs) fixed at stage 1 of each step:
      svv_form = 'entropy' (default): flux form in the entropy variables w = (cs^2 s - u^2/2, u),
        A_L = F^-1 sig M F, A_H = F^-1 sig (1 - M) F, sig = i k sqrt(Q) (k sqrt(Q) at Nyquist),
        f_rho = (eps rho/cs^2)(A_L w_rho + u A_L u), f_m = u f_rho + eps rho A_L u,
        f_h = (kappa_h rho/cs^2) A_H w_rho, kappa_h = eps cs^2/(cs^2 + max u^2),
        d rho/dt = -A_L^* f_rho - A_H^* f_h, dm/dt = -A_L^* f_m, ds/dt = (d rho/dt)/rho:
        dE/dt = -sum(f_rho A_L w_rho + f_m A_L u + f_h A_H w_rho) <= 0 for any state;
        dt_visc = 1/(zeta_max kd2_max/cfl_visc + lam_svv/cfl_svv), lam_svv by the power method
        on the Jacobian of the SVV terms (finite difference, eps frozen; svv_power)
      svv_form = 'linear': lambda(k) = eps k^2 Q(|k|), dm^/dt -= lambda m^,
        ds/dt += (1/rho) F^-1(-lambda F rho), dt_visc = cfl_visc / max(zeta_max kd2_max + eps lam_m, eps lam_s),
    time stepping: eSSPIFRK3 of advance_common.F90 without linear terms.
The defaults are those of the Fortran code (params.F90).
Prototype-only variants (REPORT_stage2.md), d_N = F^-1 i k_d M_d F with M_d = deriv_mask (default: the m filter):
    formulation = 'A': state (s, u), ds/dt = -u d_N s - d_N u, du/dt = -u d_N u - cs^2 d_N s;
                  'B': state (s, u), ds/dt = -d_N(rho u)/rho (rho u on the grid), du/dt as in A;
                  'C' (default): the Fortran formulation above. filter_m filters u in A and B.
      The dissipative terms are those of C in (rho, m), mapped pointwise for A and B:
      ds/dt = (d rho/dt)/rho, du/dt = (dm/dt - u d rho/dt)/rho.
    visc_form = 'expanded': dm/dt = rho zeta d_N D + D d_N(rho zeta) instead of d_N q;
    filter_rhs = True: the m (u) filter acts on the initial state and on the tendency instead of the state;
    budget(): mass, momentum and energy rates of the ideal, viscous and SVV parts.
It is used to cross-check the Fortran output (conservation.dat) and to study the
discretization quickly; it is not a replacement for the Fortran code.

    from proto1d import Isoth1D
    p = Isoth1D(512); p.riemann(1.0, 1.0, 2.0, -2.0); p.run(dt=1e-4, nstep=100)
"""
import math
import numpy as np


class Isoth1D:
    def __init__(self, N, lx=1.0, cs=1.0, C=2.0, nmax=1, nsmooth=1, filter_s=False, filter_m=True,
                 shock_dl=0.0, cfl=0.5, cfl_visc=1.0, pad=1, svv_eps0=0.5, svv_mN=0.5, svv_above_kc=True,
                 svv_form='entropy', cfl_svv=2.0, svv_pi_init=50, svv_pi_check=100,
                 formulation='C', visc_form='flux', deriv_mask=None, filter_rhs=False):
        # pad > 1 (prototype only): nonlinear terms on a pad*N grid (over-integration), to test aliasing
        assert formulation in ('A', 'B', 'C') and visc_form in ('flux', 'expanded'), (formulation, visc_form)
        assert deriv_mask in (None, 'mask', 'none') and not (filter_rhs and not filter_m), (deriv_mask, filter_rhs)
        self.formulation, self.visc_form, self.filter_rhs = formulation, visc_form, filter_rhs
        # the Stage 1 code path (rhs) for C, flux form and the default masks; the other variants use parts()
        self.legacy = formulation == 'C' and visc_form == 'flux' and deriv_mask is None and not filter_rhs
        self.N, self.L, self.cs = N, 2.0 * math.pi * lx, cs
        self.pad, self.Ng = pad, pad * N
        self.dl = self.L / N
        self.x = np.arange(self.Ng) * self.L / self.Ng
        n = np.arange(N // 2 + 1)
        self.k = 2.0 * math.pi / self.L * n
        self.kd = self.k.copy()
        self.kd[-1] = 0.0                                  # Nyquist
        self.mask = (n < N / 3.0).astype(float)            # dealias.F90: |k| >= kmax*2/3 -> 0
        self.Fs = self.mask if filter_s else np.ones_like(self.mask)
        self.Fm = self.mask if filter_m else np.ones_like(self.mask)
        # Md: mask in d_N and the low/high split of the SVV; F0m: filter of the initial m (u); Fr: filter of dm/dt
        self.Md = self.Fm if deriv_mask is None else (self.mask if deriv_mask == 'mask' else np.ones_like(self.mask))
        self.F0m = self.Fm
        if filter_rhs:
            self.Md, self.F0m, self.Fr = self.mask, self.mask, self.mask
            self.Fm = np.ones_like(self.mask)
        self.C, self.nmax, self.nsmooth = C, nmax, nsmooth
        self.Delta = shock_dl if shock_dl > 0.0 else self.dl
        self.cfl, self.cfl_visc = cfl, cfl_visc
        self.kd2_max = np.max((self.kd * self.Md) ** 2)
        # SVV kernel (shock_viscosity.F90: svv_kernel); k itself, so the Nyquist mode is damped
        kc = (2.0 * math.pi / self.L * (N // 2)) * 2.0 / 3.0
        mN = svv_mN * kc
        a = self.k
        with np.errstate(divide='ignore', invalid='ignore'):
            r = np.where(a > mN, (kc - a) / np.where(a > mN, a - mN, 1.0), np.inf)
        Q = np.where(a >= kc, 1.0 if svv_above_kc else 0.0, np.where((a > mN) & (r < 30.0), np.exp(-np.minimum(r, 30.0) ** 2), 0.0))
        self.svv_eps0 = svv_eps0
        self.kq2 = a ** 2 * Q
        self.lam_m = np.max(self.kq2 * (self.F0m > 0))
        self.lam_s = np.max(self.kq2 * (self.Fs > 0))
        self.eps = 0.0
        assert svv_form in ('entropy', 'linear'), svv_form
        self.svv_form, self.cfl_svv = svv_form, cfl_svv
        self.svv_pi_init, self.svv_pi_check = svv_pi_init, svv_pi_check
        sig = 1j * a * np.sqrt(Q)
        if N % 2 == 0:
            sig[-1] = a[-1] * np.sqrt(Q[-1])               # real at Nyquist, so that the Nyquist mode is damped
        self.sigL, self.sigH = sig * self.Md, sig * (1.0 - self.Md)
        self.lam_svv, self.pi_v, self.pi_conv = 0.0, None, None
        self.prepared = -1
        self.t, self.step = 0.0, 0
        self.hist = []
        self.track_budget = False          # budget() of the stage-1 state in each hist entry

    # ------------------------------------------------------------ transforms
    def fwd(self, f):
        return np.fft.rfft(f)[:self.N // 2 + 1] / self.Ng

    def bwd(self, fh):
        if self.pad > 1:
            fh = np.concatenate([fh, np.zeros(self.Ng // 2 - self.N // 2, dtype=complex)])
        return np.fft.irfft(fh * self.Ng, n=self.Ng)

    # ------------------------------------------------------------ initial conditions
    def set_state(self, rho, u):
        self.set_su(np.log(rho), u)

    def set_su(self, s, u):
        self.sh = self.fwd(s) * self.Fs
        self.mh = self.fwd(np.exp(s) * u if self.formulation == 'C' else u) * self.F0m

    def wave1d(self, kind, amp=1e-3, mode=1, rho0=1.0):
        """init_wave1d of fields.F90 along x: 'sound_wave' s = ln rho0 + amp cos(k x), u = 0;
        'simple_wave' s = ln rho0 + amp sin(k x), u = cs amp sin(k x)."""
        kx = 2.0 * math.pi * mode / self.L * self.x
        if kind == 'sound_wave':
            self.set_su(math.log(rho0) + amp * np.cos(kx), np.zeros_like(kx))
        else:
            self.set_su(math.log(rho0) + amp * np.sin(kx), self.cs * amp * np.sin(kx))

    def riemann(self, rho_L, rho_R, u_L, u_R, width=1.0, width2=None):
        """init_riemann of fields.F90; width2 (prototype only) sets the tanh width at 3L/4."""
        x, L = self.x, self.L

        def step(x0, w):
            if w > 0.0 and self.pad > 1:
                xc = np.arange(self.N) * self.dl
                return np.interp(x, xc, 0.5 * (1.0 + np.tanh((xc - x0) / (w * self.dl))), period=L)
            if w > 0.0:
                return 0.5 * (1.0 + np.tanh((x - x0) / (w * self.dl)))
            return np.where(x > x0, 1.0, np.where(x == x0, 0.5, 0.0))
        h = step(0.25 * L, width) - step(0.75 * L, width if width2 is None else width2)
        rho = (1.0 - h) * rho_L + h * rho_R
        self.set_state(rho, (1.0 - h) * u_L + h * u_R)

    # ------------------------------------------------------------ operators
    def ddx(self, fh, M):
        return self.bwd(1j * self.kd * M * fh)

    def zeta_of(self, D):
        z = np.maximum(0.0, -D)
        for _ in range(self.nmax * self.pad):
            z = np.maximum(np.maximum(np.roll(z, 1), z), np.roll(z, -1))
        for _ in range(self.nsmooth):
            z = 0.5 * (z + 0.5 * (np.roll(z, self.pad) + np.roll(z, -self.pad)))
        return self.C * self.Delta ** 2 * z

    def fields(self, sh, mh):
        s = self.bwd(sh)
        rho = np.exp(s)
        if self.formulation == 'C':
            m = self.bwd(mh)
            u = m / rho
            D = self.ddx(self.fwd(u), self.Md)
        else:
            u = self.bwd(mh)
            m = rho * u
            D = self.ddx(mh, self.Md)
        zeta = self.zeta_of(D) if self.C > 0.0 else np.zeros_like(D)
        return s, m, rho, u, D, zeta

    def svv_eps(self, u):
        return self.svv_eps0 * self.dl * np.max(np.abs(u) + self.cs) if self.svv_eps0 > 0.0 else 0.0

    def svv_entropy(self, s, rho, u, uh=None):
        """ES-SVV: (d rho/dt, dm^/dt, D_svv density) for the frozen eps; uh = F[u] if already known."""
        c2 = self.cs ** 2
        wh = self.fwd(c2 * s - 0.5 * u * u)
        uh = self.fwd(u) if uh is None else uh
        g_r = self.bwd(self.sigL * wh)
        g_m = self.bwd(self.sigL * uh)
        g_h = self.bwd(self.sigH * wh)
        f_r = self.eps * rho / c2 * (g_r + u * g_m)
        f_m = u * f_r + self.eps * rho * g_m
        f_h = self.eps / (c2 + np.max(u * u)) * rho * g_h
        rdot = self.bwd(-np.conj(self.sigL) * self.fwd(f_r) - np.conj(self.sigH) * self.fwd(f_h))
        mdh = -np.conj(self.sigL) * self.fwd(f_m)
        return rdot, mdh, f_r * g_r + f_m * g_m + f_h * g_h

    def svv_G(self, s, m):
        """SVV tendencies (ds/dt, dm/dt) in real space for the frozen eps (power method)."""
        rho = np.exp(s)
        rdot, mdh, _ = self.svv_entropy(s, rho, m / rho)
        return rdot / rho, self.bwd(mdh)

    def svv_power(self, s, m, niter, tol=0.0, h=1e-7):
        """Power method on the Jacobian J of svv_G at (s, m); J v by finite differences.
        The estimate |J v|/|v| is kept in lam_svv and v in pi_v (warm start for the next step).
        Returns the number of iterations (stops early when the relative change is < tol)."""
        if self.pi_v is None:
            # deterministic start (function of the global grid index only): Nyquist-like in s, with a modulation
            # (shock_viscosity.F90: svv_power, whose sign and modulation run over the active directions only)
            i = np.arange(self.Ng)
            vs = (-1.0) ** i * (1.0 + 0.3 * np.cos(2.0 * math.pi * i / self.Ng + 1.0))
            self.pi_v = (vs, np.zeros(self.Ng))
        G0s, G0m = self.svv_G(s, m)
        vs, vm = self.pi_v
        lam_old = None
        for it in range(1, niter + 1):
            nv = math.sqrt(np.mean(vs * vs) + np.mean(vm * vm))
            vs, vm = vs / nv, vm / nv
            Gs, Gm = self.svv_G(s + h * vs, m + h * vm)
            vs, vm = (Gs - G0s) / h, (Gm - G0m) / h
            lam = math.sqrt(np.mean(vs * vs) + np.mean(vm * vm))
            # converged: two successive iterates of this call agree to tol
            done = tol > 0.0 and lam_old is not None and abs(lam - lam_old) <= tol * lam
            lam_old = lam
            if done:
                break
        self.lam_svv, self.pi_v = lam, (vs, vm)
        return it

    def prepare(self):
        """Stage-1 quantities fixed over the step (Fortran: get_nonlinear_terms, istage = 1):
        eps and, for svv_form = 'entropy', the power-method estimate lam_svv
        (svv_pi_init iterations at the first step, then one per step; every svv_pi_check steps
        until converged, kept in pi_conv for the check dt*lam <= 2.51 in advance)."""
        if self.prepared == self.step:
            return
        s0 = self.bwd(self.sh)
        m0 = self.bwd(self.mh)
        if self.formulation != 'C':
            m0 = np.exp(s0) * m0          # A, B: the power method runs on the SVV of C in (s, m)
        self.eps = self.svv_eps(m0 / np.exp(s0))
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

    def rhs(self, sh, mh):
        if not self.legacy:
            P, f = self.parts(sh, mh)
            return P['ideal'][0] + P['visc'][0] + P['svv'][0], P['ideal'][1] + P['visc'][1] + P['svv'][1], f
        s, m, rho, u, D, zeta = self.fields(sh, mh)
        q = rho * zeta * D
        divm = self.ddx(mh, self.Fm)
        fs = self.fwd(-divm / rho)
        fm = -1j * self.kd * self.fwd(rho * u * u + self.cs ** 2 * rho - q)
        d_svv = None
        if self.eps > 0.0 and self.svv_form == 'linear':
            lam = self.eps * self.kq2
            fm = fm - lam * mh
            fs = fs + self.fwd(self.bwd(-lam * self.fwd(rho)) / rho)
        elif self.eps > 0.0:
            rdot, mdh, d_svv = self.svv_entropy(s, rho, u)
            fs = fs + self.fwd(rdot / rho)
            fm = fm + mdh
        return fs, fm, dict(s=s, rho=rho, u=u, D=D, zeta=zeta, q=q, divm=divm, d_svv=d_svv)

    def parts(self, sh, mh):
        """Spectral tendencies (ds/dt, dm/dt) for C or (ds/dt, du/dt) for A and B, split into the ideal,
        shock-viscosity and SVV parts, before the stage filter (Fr included for filter_rhs); and the fields."""
        s, m, rho, u, D, zeta = self.fields(sh, mh)
        c2, Md, ikd = self.cs ** 2, self.Md, 1j * self.kd * self.Md
        C = self.formulation == 'C'
        zero = np.zeros_like(sh)
        divm = self.ddx(self.fwd(m), Md) if not C else self.ddx(mh, Md)
        if C:
            ideal = (self.fwd(-divm / rho), -ikd * self.fwd(rho * u * u + c2 * rho))
        else:
            ds = self.ddx(sh, Md)
            sdot = -divm / rho if self.formulation == 'B' else -u * ds - D
            ideal = (self.fwd(sdot), self.fwd(-u * D - c2 * ds))
        # dissipative parts: d rho/dt (grid) and dm/dt (spectral) of C, then mapped to (s, u) for A and B
        q = rho * zeta * D
        if self.visc_form == 'flux':
            mv = ikd * self.fwd(q)
        else:
            rz = rho * zeta
            mv = self.fwd(rz * self.ddx(self.fwd(D), Md) + D * self.ddx(self.fwd(rz), Md))
        visc = (zero, mv if C else self.fwd(self.bwd(mv) / rho))
        d_svv = None
        svv = (zero, zero)
        if self.eps > 0.0:
            if self.svv_form == 'linear':
                lam = self.eps * self.kq2
                rdot = self.bwd(-lam * self.fwd(rho))
                mdh = -lam * (mh if C else self.fwd(m))
            else:
                rdot, mdh, d_svv = self.svv_entropy(s, rho, u, None if C else mh)
            svv = (self.fwd(rdot / rho), mdh if C else self.fwd((self.bwd(mdh) - u * rdot) / rho))
        P = dict(ideal=ideal, visc=visc, svv=svv)
        if self.filter_rhs:
            P = {k: (a, b * self.Fr) for k, (a, b) in P.items()}
        return P, dict(s=s, rho=rho, u=u, D=D, zeta=zeta, q=q, divm=divm, d_svv=d_svv)

    def budget(self, sh=None, mh=None):
        """Semi-discrete budget at (sh, mh) (default: the current state). The tendencies of each part are mapped
        to (d rho/dt, dm/dt) on the grid (rho ds/dt, rho du/dt + u rho ds/dt for A and B); for each part
        dM = sum(rho dot) dV, dP = sum(m dot) dV, eps = sum(w_rho rho dot + u m dot) dV with
        w_rho = cs^2 s - u^2/2, the dE/drho at fixed m of E = sum[rho u^2/2 + cs^2 rho (s - 1)] dV.
        Keys <part>_<rate> after the stage filter (Fs, Fm) and <part>_<rate>_raw before it; the scales
        A_M = sum|rho dot_ideal| dV, A_P = sum|m dot_ideal| dV, A_svv = sum|w_rho rho dot + u m dot|_svv dV (raw),
        visc = sum(rho zeta D^2) dV, D_svv = sum(d_svv) dV; R_zeta = eps_visc + visc, R_svv = eps_svv + D_svv."""
        sh = self.sh if sh is None else sh
        mh = self.mh if mh is None else mh
        P, f = self.parts(sh, mh)
        s, rho, u = f['s'], f['rho'], f['u']
        dV = self.L / self.Ng
        w = self.cs ** 2 * s - 0.5 * u * u
        out = {}
        for tag, Fs, Fm in (('_raw', 1.0, 1.0), ('', self.Fs, self.Fm)):
            for k, (a, b) in P.items():
                rdot = rho * self.bwd(a * Fs)
                b = self.bwd(b * Fm)
                mdot = b if self.formulation == 'C' else rho * b + u * rdot
                e = w * rdot + u * mdot
                out[k + '_dM' + tag], out[k + '_dP' + tag] = np.sum(rdot) * dV, np.sum(mdot) * dV
                out[k + '_eps' + tag] = np.sum(e) * dV
                if tag == '_raw' and k == 'ideal':
                    out['A_M'], out['A_P'] = np.sum(np.abs(rdot)) * dV, np.sum(np.abs(mdot)) * dV
                if tag == '_raw' and k == 'svv':
                    out['A_svv'] = np.sum(np.abs(e)) * dV
            for r in ('dM', 'dP', 'eps'):
                out[r + tag] = sum(out[k + '_' + r + tag] for k in P)
        out['visc'] = np.sum(rho * f['zeta'] * f['D'] ** 2) * dV
        out['D_svv'] = np.sum(f['d_svv']) * dV if f['d_svv'] is not None else 0.0
        for tag in ('_raw', ''):
            out['R_zeta' + tag] = out['visc_eps' + tag] + out['visc']
            out['R_svv' + tag] = out['svv_eps' + tag] + out['D_svv'] if f['d_svv'] is not None else 0.0
        return out

    def monitor(self, f):
        rho, u = f['rho'], f['u']
        dV = self.L / self.Ng
        if self.svv_form == 'entropy':
            rate = f['zeta'].max() * self.kd2_max / self.cfl_visc + self.lam_svv / self.cfl_svv
            dt_visc = 1.0 / rate if rate > 0 else np.inf
        else:
            e = self.svv_eps(u)
            rate = max(f['zeta'].max() * self.kd2_max + e * self.lam_m, e * self.lam_s)
            dt_visc = self.cfl_visc / rate if rate > 0 else np.inf
        E = np.sum(0.5 * rho * u * u + self.cs ** 2 * rho * (f['s'] - 1.0)) * dV
        return dict(t=self.t, M=np.sum(rho) * dV, P=np.sum(rho * u) * dV, E=E, E_kin=np.sum(0.5 * rho * u * u) * dV,
                    rho_min=rho.min(), rho_max=rho.max(), zeta_max=f['zeta'].max(), D_min=f['D'].min(),
                    visc=np.sum(rho * f['zeta'] * f['D'] ** 2) * dV, eps_svv=self.eps,
                    lam_svv=self.lam_svv, F_svv=self.lam_svv / (self.eps * self.k[-1] ** 2) if self.eps > 0 else 0.0,
                    D_svv=np.sum(f['d_svv']) * dV if f['d_svv'] is not None else 0.0,
                    D_svv_min=f['d_svv'].min() if f['d_svv'] is not None else 0.0,
                    dt_adv=self.cfl / (np.max(np.abs(u) + self.cs) / self.dl), dt_visc=dt_visc)

    def svv_budget(self, f):
        """(energy rate of the SVV tendencies, D_svv) at the state of f (svv_form = 'entropy');
        the rate sum[(dE/ds) ds/dt + u dm/dt] dV should equal -D_svv up to round-off."""
        s, rho, u = f['s'], f['rho'], f['u']
        rdot, mdh, d = self.svv_entropy(s, rho, u)
        dV = self.L / self.Ng
        dEds = rho * (self.cs ** 2 * (s + 1.0) - 0.5 * u * u)
        rate = np.sum(dEds * (rdot / rho) + u * self.bwd(mdh)) * dV
        return rate, np.sum(d) * dV

    # ------------------------------------------------------------ time stepping
    def advance(self, dt):
        self.prepare()
        s0, m0 = self.sh, self.mh
        f0s, f0m, f = self.rhs(s0, m0)
        mon = self.monitor(f)
        if self.track_budget:
            mon.update(self.budget(s0, m0))
        if self.pi_conv is not None:
            mon['pi_conv'] = self.pi_conv
            if dt * self.pi_conv[0] > 2.51:
                raise FloatingPointError('svv_power: dt*lam_svv = %.3f > 2.51 at t = %.6g' % (dt * self.pi_conv[0], self.t))
        s1 = (s0 + 2.0 / 3.0 * dt * f0s) * self.Fs
        m1 = (m0 + 2.0 / 3.0 * dt * f0m) * self.Fm
        f1s, f1m, _ = self.rhs(s1, m1)
        s2 = (2.0 / 3.0 * s0 + 1.0 / 3.0 * s1 + 4.0 / 9.0 * dt * f1s) * self.Fs
        m2 = (2.0 / 3.0 * m0 + 1.0 / 3.0 * m1 + 4.0 / 9.0 * dt * f1m) * self.Fm
        f2s, f2m, _ = self.rhs(s2, m2)
        self.sh = (37.0 / 64.0 * s0 + 5.0 / 32.0 * dt * f0s + 27.0 / 64.0 * s2 + 9.0 / 16.0 * dt * f2s) * self.Fs
        self.mh = (37.0 / 64.0 * m0 + 5.0 / 32.0 * dt * f0m + 27.0 / 64.0 * m2 + 9.0 / 16.0 * dt * f2m) * self.Fm
        self.t += dt
        self.step += 1
        return mon, f

    def dt_limit(self):
        """min(dt_adv, dt_visc) of the current state (the CFL of advance.F90)."""
        self.prepare()
        mon = self.monitor(self.rhs(self.sh, self.mh)[2])
        return min(mon['dt_adv'], mon['dt_visc'])

    def run_until(self, t_end, dt_max, every=0, callback=None, maxstep=10 ** 7):
        """Adaptive dt = min(dt_max, dt_limit) (the CFL numbers are in dt_adv, dt_visc)."""
        while self.t < t_end * (1.0 - 1e-12) and self.step < maxstep:
            dt = min(dt_max, self.dt_limit(), t_end - self.t)
            try:
                mon, f = self.advance(dt)
            except FloatingPointError:
                return 'svv_unstable'
            mon['dt'] = dt
            self.hist.append(mon)
            if not np.all(np.isfinite(f['rho'])) or not np.all(np.isfinite(self.sh)):
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
            if not np.all(np.isfinite(f['rho'])):
                return 'nonfinite'
            if stop_on_cfl and dt > min(mon['dt_adv'], mon['dt_visc']):
                return 'cfl'
            if callback is not None and every > 0 and n % every == 0:
                callback(self, f)
        return 'ok'

    def state(self):
        s, m, rho, u, D, zeta = self.fields(self.sh, self.mh)
        return dict(x=self.x, s=s, rho=rho, u=u, D=D, zeta=zeta)


def from_case(case):
    """Isoth1D configured from <case>/calliope.in (riemann along x, B = 0); returns (proto, dt, nstep)."""
    import os
    from analyze_stage1 import read_nml, num, nums, strv
    d = read_nml(os.path.join(case, 'calliope.in'))
    T = lambda k, dflt: strv(d.get(k, dflt)).upper().startswith(('T', '.T'))
    assert strv(d['init_type']) == 'riemann' and strv(d.get('dir', "'x'")) == 'x'
    p = Isoth1D(int(d['nlx']), lx=num(d['lx']), cs=math.sqrt(num(d['beta0']) / 2.0),
                C=num(d['c_shock']), nmax=int(d['shock_nmax']), nsmooth=int(d['shock_nsmooth']),
                filter_s=T('filter_lnrho', 'F'), filter_m=T('filter_mom', 'T'), shock_dl=num(d['shock_dl']),
                cfl=num(d['cfl']), cfl_visc=num(d['cfl_visc']), svv_eps0=num(d.get('svv_eps0', '0.0')),
                svv_mN=num(d.get('svv_mn', '0.5')), svv_above_kc=T('svv_above_kc', 'T'),
                svv_form=strv(d.get('svv_form', "'entropy'")).strip("'\""), cfl_svv=num(d.get('cfl_svv', '2.0')),
                svv_pi_init=int(num(d.get('svv_pi_init', '50'))), svv_pi_check=int(num(d.get('svv_pi_check', '100'))))
    p.riemann(num(d['rho_l']), num(d['rho_r']), nums(d['u_l'])[0], nums(d['u_r'])[0], width=num(d['width']))
    # the Fortran volume integrals include the transverse area (2 pi ly)(2 pi lz)
    p.area = (2.0 * math.pi) ** 2 * num(d['ly']) * num(d['lz'])
    return p, num(d['dt']), int(d['nstep'])


if __name__ == '__main__':
    # cross-check with a Fortran run:  python3 proto1d.py <case_dir> [nstep]
    # conservation.dat row i and cfl.dat row i are the state after i steps (stage 1 of step i+1)
    import sys, os
    from analyze_stage1 import load
    case = sys.argv[1]
    p, dt, nstep = from_case(case)
    nstep = min(nstep, int(sys.argv[2])) if len(sys.argv) > 2 else nstep
    st = p.run(dt, nstep)
    print('status', st, 'steps', p.step)
    c = load(os.path.join(case, 'conservation.dat'))
    f = load(os.path.join(case, 'cfl.dat'))
    cols = [('M', 2), ('Px', 4), ('E', 10), ('rho_min', 11), ('rho_max', 12), ('zeta_max', 13), ('D_min', 14),
            ('visc', 15)]
    worst = {}
    for i in range(min(len(c), len(p.hist))):
        h = p.hist[i]
        for name, j in cols:
            key = 'P' if name == 'Px' else name
            hv = h[key] * p.area if name in ('M', 'Px', 'E', 'visc') else h[key]
            scale = abs(c[i, j]) if name != 'Px' else abs(c[0, 2])
            e = abs(hv - c[i, j]) / max(scale, 1e-300)
            worst[name] = max(worst.get(name, 0.0), e)
        if f is not None and f.shape[1] > 9 and i < len(f):
            worst['eps_svv'] = max(worst.get('eps_svv', 0.0), abs(h['eps_svv'] - f[i, 9]) / max(f[i, 9], 1e-300))
            worst['dt_visc'] = max(worst.get('dt_visc', 0.0), abs(h['dt_visc'] / f[i, 3] - 1.0))
            worst['dt_adv'] = max(worst.get('dt_adv', 0.0), abs(h['dt_adv'] / f[i, 2] - 1.0))
    print('rows compared:', min(len(c), len(p.hist)))
    print('max relative difference (Px relative to M0):')
    for k, v in worst.items():
        print('  %-9s %.2e' % (k, v))
