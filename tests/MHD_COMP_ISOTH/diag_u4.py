#!/usr/bin/env python3
"""Diagnostics of the U4 colliding flow (rho = 1, u = +-4) with proto1d.py (REPORT_stage1.md, sec. U4).

    python3 diag_u4.py cause  N [default zmask smoothcoll]   star noise with/without zeta in the rarefaction half
    python3 diag_u4.py remote N [N ...]                      grid-scale forcing of the star from far away (t = 0.3)
    python3 diag_u4.py onset  N [default|taper]              smooth collision: star noise vs the time the shocks form
    python3 diag_u4.py taper  N [taper default]              raised-cosine taper of i k and the m-SVV near k_c
    python3 diag_u4.py fsmass Q1|Q0                          filter_lnrho = T: mass error with Q(k >= k_c) = 1 or 0 (U2)
    python3 diag_u4.py svvsplit                              dE/dt of the SVV split into its m and ln rho terms (U4, U2, Sod)

The star region is |x - 3L/4| < 0.5 t (rarefaction star, rho* = e^-4); grid-scale noise is measured with
the local high-pass hp (4th difference, unit gain at k = k_c), not with a spectral band (sinc tails of the shocks).
The SVV is svv_form = 'linear' (the form before the entropy-stable SVV, whose U4 results these diagnose),
and dt is that of its Fortran shock batch (gen_stage1.py: SAFETY*DT512_LINEAR*512/N)."""
import sys
import os
import math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from proto1d import Isoth1D

DT_U4_512 = 0.7 * 2.982e-4
RHO_STAR = math.exp(-4.0)
LINEAR = dict(svv_form='linear')


def hp(f):
    return (-np.roll(f, 2) + 4 * np.roll(f, 1) - 6 * f + 4 * np.roll(f, -1) - np.roll(f, -2)) / 9.0


class ZMask(Isoth1D):
    """zeta = 0 for x > L/2 (no shock viscosity in the rarefaction half; the shocks keep theirs)."""
    def zeta_of(self, D):
        return super().zeta_of(D) * (self.x < 0.5 * self.L)


class Taper(Isoth1D):
    """i k and the m-SVV multiplied by sigma = 1 (k <= k1 k_c) -> 0 (k = k_c), raised cosine;
    D = d_N u uses the same i k, so the R_zeta identity still holds. The SVV of ln rho is unchanged."""
    def __init__(self, N, k1=0.8, **kw):
        super().__init__(N, **kw)
        n = np.arange(N // 2 + 1)
        kc = N / 3.0
        r = np.clip((n - k1 * kc) / ((1.0 - k1) * kc), 0.0, 1.0)
        sig = np.where(n < kc, 0.5 * (1.0 + np.cos(math.pi * r)), 0.0)
        self.kd = self.kd * sig
        self.kq2_s = self.kq2.copy()
        self.kq2 = self.kq2 * np.where(n < kc, sig, 1.0)

    def rhs(self, sh, mh):
        s, m, rho, u, D, zeta = self.fields(sh, mh)
        q = rho * zeta * D
        divm = self.ddx(mh, self.Fm)
        fs = self.fwd(-divm / rho)
        fm = -1j * self.kd * self.fwd(rho * u * u + self.cs ** 2 * rho - q)
        if self.eps > 0.0:
            fm = fm - self.eps * self.kq2 * mh
            fs = fs + self.fwd(self.bwd(-self.eps * self.kq2_s * self.fwd(rho)) / rho)
        return fs, fm, dict(s=s, rho=rho, u=u, D=D, zeta=zeta, q=q, divm=divm)


def star_line(N, tag, p):
    st = p.state()
    x = st['x'] - 0.75 * p.L
    star = np.abs(x) < 0.5 * p.t
    fan = (np.abs(x) > 1.2 * p.t) & (np.abs(x) < 4.8 * p.t)
    m = st['rho'] * st['u']
    return ('N=%5d %-10s t=%.3f  star: hp m %.1e s %.1e u %.1e  rho/rho*-1 %+.2e  zeta mean %.1e | fan zeta mean %.1e | rho_min %.4f' % (
        N, tag, p.t, hp(m)[star].std(), hp(st['s'])[star].std(), hp(st['u'])[star].std(),
        st['rho'][star].mean() / RHO_STAR - 1, st['zeta'][star].mean(), st['zeta'][fan].mean(), st['rho'].min()))


def advance_to(p, t_out, dt):
    p.run(dt, int(round((t_out - p.t) / dt)), stop_on_cfl=False)
    return np.all(np.isfinite(p.sh)) and np.all(np.isfinite(p.mh))


def cause(N, tags):
    dt = DT_U4_512 * 512 / N
    for tag in tags or ['default', 'zmask', 'smoothcoll']:
        p = (ZMask if tag == 'zmask' else Isoth1D)(N, **LINEAR)
        # smoothcoll: tanh width N/(2 pi) cells = 1 in x at the collision (no shock before t ~ 0.25)
        p.riemann(1.0, 1.0, 4.0, -4.0, width=(N / (2 * math.pi) if tag == 'smoothcoll' else 1.0), width2=1.0)
        for t_out in ((0.2,) if tag == 'smoothcoll' else (0.2, 0.4)):
            with np.errstate(all='ignore'):
                ok = advance_to(p, t_out, dt)
            if not ok:
                print('N=%5d %-10s non-finite by t=%.3f' % (N, tag, p.t))
                break
            print(star_line(N, tag, p), flush=True)


def remote(Ns):
    """The RHS of m and s is linear in the field it acts on once (rho, u, q) are fixed:
    RHS[W f] + RHS[(1 - W) f] = RHS[f]. The far part (1 - W) reaches the star only through spectral nonlocality."""
    for N in Ns:
        p = Isoth1D(N, **LINEAR)
        p.riemann(1.0, 1.0, 4.0, -4.0)
        L, dt = p.L, DT_U4_512 * 512 / N
        x0 = 0.75 * L
        advance_to(p, 0.3, dt)
        p.eps = p.svv_eps(p.bwd(p.mh) / np.exp(p.bwd(p.sh)))
        _, _, f = p.rhs(p.sh, p.mh)
        rho, u, q, x = f['rho'], f['u'], f['q'], p.x
        m = rho * u
        star = np.abs(x - x0) < 0.2 * p.t
        F = rho * u * u + p.cs ** 2 * rho - q
        lam = p.eps * p.kq2
        lam_c = p.eps * (N / 3.0) ** 2

        def m_rhs(w):
            return p.bwd((-1j * p.kd * p.fwd(w * F) - lam * p.fwd(w * m)) * p.Fm)

        for a in (1.0, 2.0, 2.8):
            W = 0.5 * (np.tanh((x - (x0 - a)) / 0.1) - np.tanh((x - (x0 + a)) / 0.1))
            far = hp(m_rhs(1.0 - W))[star].std()
            print('N=%5d t=%.2f a=%.1f: far hp(dm/dt) %.2e  far/lambda(k_c) %.2e | measured hp m %.2e  lambda(k_c) %.0f' % (
                N, p.t, a, far, far / lam_c, hp(m)[star].std(), lam_c), flush=True)
        Wsh = 0.5 * (np.tanh((x - (0.25 * L - 0.3)) / 0.05) - np.tanh((x - (0.25 * L + 0.3)) / 0.05))
        print('          shocks only (|x - L/4| < 0.3): hp(dm/dt) %.2e' % hp(m_rhs(Wsh))[star].std(), flush=True)


def onset(N, tag):
    p = Taper(N, **LINEAR) if tag == 'taper' else Isoth1D(N, **LINEAR)
    p.riemann(1.0, 1.0, 4.0, -4.0, width=N / (2 * math.pi), width2=1.0)
    dt = DT_U4_512 * 512 / N
    for j in range(1, 21):
        advance_to(p, 0.025 * j, dt)
        st = p.state()
        x = st['x']
        star = np.abs(x - 0.75 * p.L) < 0.5 * p.t
        coll = np.abs(x - 0.25 * p.L) < 1.0
        print('N=%5d %-7s t=%.3f  star hp(u) %.1e hp(m) %.1e  rho/rho*-1 %+.2e | zeta_max collision %.1e  rarefaction half %.1e' % (
            N, tag, p.t, hp(st['u'])[star].std(), hp(st['rho'] * st['u'])[star].std(),
            st['rho'][star].mean() / RHO_STAR - 1, st['zeta'][coll].max(), st['zeta'][x > 0.5 * p.L].max()), flush=True)


def taper(N, tags):
    dt = DT_U4_512 * 512 / N
    for tag in tags or ['taper']:
        p = Taper(N, **LINEAR) if tag == 'taper' else Isoth1D(N, **LINEAR)
        p.riemann(1.0, 1.0, 4.0, -4.0)
        for t_out in (0.2, 0.4):
            if not advance_to(p, t_out, dt):
                print('N=%5d %-10s non-finite by t=%.3f' % (N, tag, p.t))
                break
            print(star_line(N, tag, p) + '  dP %.1e' % abs(p.hist[-1]['P'] - p.hist[0]['P']), flush=True)


def fsmass(which):
    q1 = which == 'Q1'
    p = Isoth1D(1024, filter_s=True, svv_above_kc=q1, **LINEAR)
    p.riemann(1.0, 1.0, 2.0, -2.0)
    dt = 0.7 * 5.665e-4 * 0.5
    status = p.run(dt, int(round(0.4 / dt)), stop_on_cfl=False)
    M = np.array([h['M'] for h in p.hist])
    print('U2 N=1024 filter_lnrho=T Q(k>=k_c)=%d: %s, t=%.3f, max|dM/M0| %.3e' % (q1, status, p.t, np.abs(M / M[0] - 1).max()))


def svvsplit():
    """E = sum m^2/(2 rho) + cs^2 rho (s - 1): dE/dt = sum u mdot + sum (cs^2 s - u^2/2) rho sdot.
    Per unit transverse area (the Fortran budget is (2 pi)^2 times larger)."""
    cases = {'U4': (1.0, 1.0, 4.0, -4.0, 2.982e-4), 'U2': (1.0, 1.0, 2.0, -2.0, 5.665e-4),
             'Sod': (10.0, 1.0, 0.0, 0.0, 1.124e-3)}
    N = 1024
    for name, (rL, rR, uL, uR, dt512) in cases.items():
        p = Isoth1D(N, **LINEAR)
        p.riemann(rL, rR, uL, uR)
        dt = 0.7 * dt512 * 512 / N
        for t_out in (0.05, 0.2, 0.4):
            advance_to(p, t_out, dt)
            s = p.bwd(p.sh)
            rho = np.exp(s)
            u = p.bwd(p.mh) / rho
            lam = p.svv_eps(u) * p.kq2
            dV = p.L / p.Ng
            mdot = p.bwd(-lam * p.mh * p.Fm)
            sdot = p.bwd(p.Fs * p.fwd(p.bwd(-lam * p.fwd(rho)) / rho))
            e_m = np.sum(u * mdot) * dV
            e_kin = np.sum(-0.5 * u * u * rho * sdot) * dV
            e_int = np.sum(p.cs ** 2 * s * rho * sdot) * dV
            _, _, f = p.rhs(p.sh, p.mh)
            visc = np.sum(f['rho'] * f['zeta'] * f['D'] ** 2) * dV
            print('%-3s t=%.2f  eps_svv %+.3e = m %+.3e + ln rho (kinetic) %+.3e + ln rho (internal) %+.3e | eps_shock %+.3e' % (
                name, p.t, e_m + e_kin + e_int, e_m, e_kin, e_int, -visc), flush=True)


if __name__ == '__main__':
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == 'cause':
        cause(int(args[0]), args[1:])
    elif cmd == 'remote':
        remote([int(a) for a in args])
    elif cmd == 'onset':
        onset(int(args[0]), args[1] if len(args) > 1 else 'default')
    elif cmd == 'taper':
        taper(int(args[0]), args[1:])
    elif cmd == 'fsmass':
        fsmass(args[0])
    elif cmd == 'svvsplit':
        svvsplit()
    else:
        sys.exit(__doc__)
