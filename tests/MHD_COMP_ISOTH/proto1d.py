#!/usr/bin/env python3
"""1D numpy prototype of the MHD_COMP_ISOTH discretization (B = 0, formulation C).

It mirrors advance.F90 / shock_viscosity.F90 for a run along x with nly = nlz = 2:
    state s^ = F[ln rho], m^ (rfft coefficients), filtered by the 2/3 mask after every stage,
    ds/dt = -(div_N m)/rho,  dm/dt = -div_N [rho u u + cs^2 rho - q],  q = rho zeta D,
    D = div_N u,  div_N = F^-1 i k_d M F  (k_d = 0 at Nyquist, M = 2/3 mask),
    zeta = C Delta^2 S[max(0, -D)],  S = (3-point max)^nmax then ([1,2,1]/4)^nsmooth,
    SVV (shock_viscosity.F90): lambda(k) = eps k^2 Q(|k|), eps = svv_eps0 dl max(|u| + cs)
        fixed at stage 1 of each step; dm^/dt -= lambda m^, ds/dt += (1/rho) F^-1(-lambda F rho),
    dt_visc = cfl_visc / max(zeta_max kd2_max + eps lam_m, eps lam_s),
    time stepping: eSSPIFRK3 of advance_common.F90 without linear terms.
The defaults are those of the Fortran code (params.F90).
It is used to cross-check the Fortran output (conservation.dat) and to study the
discretization quickly; it is not a replacement for the Fortran code.

    from proto1d import Isoth1D
    p = Isoth1D(512); p.riemann(1.0, 1.0, 2.0, -2.0); p.run(dt=1e-4, nstep=100)
"""
import math
import numpy as np


class Isoth1D:
    def __init__(self, N, lx=1.0, cs=1.0, C=2.0, nmax=1, nsmooth=1, filter_s=False, filter_m=True,
                 shock_dl=0.0, cfl=0.5, cfl_visc=1.0, pad=1, svv_eps0=0.5, svv_mN=0.5, svv_above_kc=True):
        # pad > 1 (prototype only): nonlinear terms on a pad*N grid (over-integration), to test aliasing
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
        self.C, self.nmax, self.nsmooth = C, nmax, nsmooth
        self.Delta = shock_dl if shock_dl > 0.0 else self.dl
        self.cfl, self.cfl_visc = cfl, cfl_visc
        self.kd2_max = np.max((self.kd * self.Fm) ** 2)
        # SVV kernel (shock_viscosity.F90: svv_kernel); k itself, so the Nyquist mode is damped
        kc = (2.0 * math.pi / self.L * (N // 2)) * 2.0 / 3.0
        mN = svv_mN * kc
        a = self.k
        with np.errstate(divide='ignore', invalid='ignore'):
            r = np.where(a > mN, (kc - a) / np.where(a > mN, a - mN, 1.0), np.inf)
        Q = np.where(a >= kc, 1.0 if svv_above_kc else 0.0, np.where((a > mN) & (r < 30.0), np.exp(-np.minimum(r, 30.0) ** 2), 0.0))
        self.svv_eps0 = svv_eps0
        self.kq2 = a ** 2 * Q
        self.lam_m = np.max(self.kq2 * (self.Fm > 0))
        self.lam_s = np.max(self.kq2 * (self.Fs > 0))
        self.eps = 0.0
        self.t, self.step = 0.0, 0
        self.hist = []

    # ------------------------------------------------------------ transforms
    def fwd(self, f):
        return np.fft.rfft(f)[:self.N // 2 + 1] / self.Ng

    def bwd(self, fh):
        if self.pad > 1:
            fh = np.concatenate([fh, np.zeros(self.Ng // 2 - self.N // 2, dtype=complex)])
        return np.fft.irfft(fh * self.Ng, n=self.Ng)

    # ------------------------------------------------------------ initial conditions
    def set_state(self, rho, u):
        self.sh = self.fwd(np.log(rho)) * self.Fs
        self.mh = self.fwd(rho * u) * self.Fm

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
        m = self.bwd(mh)
        rho = np.exp(s)
        u = m / rho
        D = self.ddx(self.fwd(u), self.Fm)
        zeta = self.zeta_of(D) if self.C > 0.0 else np.zeros_like(D)
        return s, m, rho, u, D, zeta

    def svv_eps(self, u):
        return self.svv_eps0 * self.dl * np.max(np.abs(u) + self.cs) if self.svv_eps0 > 0.0 else 0.0

    def rhs(self, sh, mh):
        s, m, rho, u, D, zeta = self.fields(sh, mh)
        q = rho * zeta * D
        divm = self.ddx(mh, self.Fm)
        fs = self.fwd(-divm / rho)
        fm = -1j * self.kd * self.fwd(rho * u * u + self.cs ** 2 * rho - q)
        if self.eps > 0.0:
            lam = self.eps * self.kq2
            fm = fm - lam * mh
            fs = fs + self.fwd(self.bwd(-lam * self.fwd(rho)) / rho)
        return fs, fm, dict(s=s, rho=rho, u=u, D=D, zeta=zeta, q=q, divm=divm)

    def monitor(self, f):
        rho, u = f['rho'], f['u']
        e = self.svv_eps(u)
        rate = max(f['zeta'].max() * self.kd2_max + e * self.lam_m, e * self.lam_s)
        dV = self.L / self.Ng
        E = np.sum(0.5 * rho * u * u + self.cs ** 2 * rho * (f['s'] - 1.0)) * dV
        return dict(t=self.t, M=np.sum(rho) * dV, P=np.sum(rho * u) * dV, E=E,
                    rho_min=rho.min(), rho_max=rho.max(), zeta_max=f['zeta'].max(), D_min=f['D'].min(),
                    visc=np.sum(rho * f['zeta'] * f['D'] ** 2) * dV, eps_svv=self.eps,
                    dt_adv=self.cfl / (np.max(np.abs(u) + self.cs) / self.dl),
                    dt_visc=self.cfl_visc / rate if rate > 0 else np.inf)

    # ------------------------------------------------------------ time stepping
    def advance(self, dt):
        s0, m0 = self.sh, self.mh
        rho0 = np.exp(self.bwd(s0))
        self.eps = self.svv_eps(self.bwd(m0) / rho0)       # fixed over the step (Fortran: stage 1)
        f0s, f0m, f = self.rhs(s0, m0)
        mon = self.monitor(f)
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
        mon = self.monitor(self.rhs(self.sh, self.mh)[2])
        return min(mon['dt_adv'], mon['dt_visc'])

    def run_until(self, t_end, dt_max, every=0, callback=None, maxstep=10 ** 7):
        """Adaptive dt = min(dt_max, dt_limit) (the CFL numbers are in dt_adv, dt_visc)."""
        while self.t < t_end * (1.0 - 1e-12) and self.step < maxstep:
            dt = min(dt_max, self.dt_limit(), t_end - self.t)
            mon, f = self.advance(dt)
            mon['dt'] = dt
            self.hist.append(mon)
            if not np.all(np.isfinite(f['rho'])) or not np.all(np.isfinite(self.sh)):
                return 'nonfinite'
            if callback is not None and every > 0 and self.step % every == 0:
                callback(self, f)
        return 'ok'

    def run(self, dt, nstep, every=0, callback=None, stop_on_cfl=True):
        for n in range(nstep):
            mon, f = self.advance(dt)
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
                svv_mN=num(d.get('svv_mn', '0.5')), svv_above_kc=T('svv_above_kc', 'T'))
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
