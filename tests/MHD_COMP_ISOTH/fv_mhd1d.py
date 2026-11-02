#!/usr/bin/env python3
"""Second-order finite-volume solver of 1D isothermal MHD: the reference solutions of Stage 3.

Cells i = 0..N-1 of width dx = L/N centred at x_i = i dx, so the grid points j L/N_s of a spectral run with
N_s | N are cell centres. Conserved U = (rho, rho u_x, rho u_y, rho u_z, B_y, B_z), B_x constant, sound speed c.
Interface flux: isothermal HLLD after Mignone (2007, J. Comput. Phys. 225, 1427). rho and u_x are constant in
the fan between the outer waves S_L, S_R (Davis estimates): rho* is the HLL average and u* = F_rho^hll/rho*, and
the fluxes of rho and m_x are the HLL ones. The tangential components (u_t, B_t) of the states *L, *R follow
from the jump conditions across S_L, S_R, and the state ** between the rotational waves u* -+ |B_x|/sqrt(rho*)
solves the (linear) tangential Riemann problem between *L and *R; the tangential flux is the physical flux of
the state at x/t = 0. An isolated rotational discontinuity is therefore exact at first order.
Reconstruction: primitive variables with the MC limiter. Time: SSP-RK2 (Heun), dt = cfl dx / max(|u_x| + c_f).
Boundaries: 'periodic' or 'outflow' (zero gradient).

python3 fv_mhd1d.py [quick]        self-tests against the exact single waves of mhd_exact.py
python3 fv_mhd1d.py ref NAME LOG2N  reference of mhd_exact.COMPOUND[NAME] on 2^LOG2N cells: R on (L/4, 3L/4)
                                    as in riemann_mhd(width = 0), periodic, snapshots T_REF -> runs/stage3/fvref
"""
import math
import os
import sys
import time

import numpy as np

import mhd_exact as mx


def cfast(rho, by, bz, bx, c):
    a = c * c + (bx * bx + by * by + bz * bz) / rho
    return np.sqrt(0.5 * (a + np.sqrt(np.maximum(a * a - 4.0 * c * c * bx * bx / rho, 0.0))))


def hlld(WL, WR, bx, c):
    """Interface fluxes (6, n) from the left and right primitive states (6, n)."""
    rL, uL, vL, wL, yL, zL = WL
    rR, uR, vR, wR, yR, zR = WR
    bx2 = bx * bx
    mL, mR = rL * uL, rR * uR
    pL = c * c * rL + 0.5 * (bx2 + yL * yL + zL * zL)
    pR = c * c * rR + 0.5 * (bx2 + yR * yR + zR * zR)
    FL = (mL, mL * uL + pL - bx2, mL * vL - bx * yL, mL * wL - bx * zL, yL * uL - bx * vL, zL * uL - bx * wL)
    FR = (mR, mR * uR + pR - bx2, mR * vR - bx * yR, mR * wR - bx * zR, yR * uR - bx * vR, zR * uR - bx * wR)
    cL, cR = cfast(rL, yL, zL, bx, c), cfast(rR, yR, zR, bx, c)
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)
    idS = 1.0 / (SR - SL)
    rs = (SR * rR - SL * rL - (FR[0] - FL[0])) * idS
    F0 = (SR * FL[0] - SL * FR[0] + SL * SR * (rR - rL)) * idS
    F1 = (SR * FL[1] - SL * FR[1] + SL * SR * (mR - mL)) * idS
    us = F0 / rs
    sq = np.sqrt(rs)
    ca = abs(bx) / sq
    SaL = np.maximum(us - ca, SL)
    SaR = np.minimum(us + ca, SR)

    def star(r, u, v, w, y, z, S):
        # rho (S - u) = rho* (S - u*) holds for the HLL state, so the 2 x 2 tangential system has the
        # determinant rho* (S - u*)^2 - B_x^2; it vanishes when S coincides with a rotational wave
        R = r * (S - u)
        den = R * (S - us) - bx2
        ok = np.abs(den) > 1e-10 * (np.abs(R * (S - us)) + bx2)
        dd = np.where(ok, den, 1.0)
        fac = np.where(ok, bx * (us - u) / dd, 0.0)
        g = np.where(ok, (R * (S - u) - bx2) / dd, 1.0)
        return v - y * fac, w - z * fac, y * g, z * g

    vl, wl, yl, zl = star(rL, uL, vL, wL, yL, zL, SL)
    vr, wr, yr, zr = star(rR, uR, vR, wR, yR, zR, SR)
    s = math.copysign(1.0, bx) if bx != 0.0 else 0.0
    v2 = 0.5 * (vl + vr) + 0.5 * s * (yr - yl) / sq
    w2 = 0.5 * (wl + wr) + 0.5 * s * (zr - zl) / sq
    y2 = 0.5 * (yl + yr) + 0.5 * s * sq * (vr - vl)
    z2 = 0.5 * (zl + zr) + 0.5 * s * sq * (wr - wl)
    inL, inR = SaL > 0.0, SaR < 0.0
    v = np.where(inL, vl, np.where(inR, vr, v2))
    w = np.where(inL, wl, np.where(inR, wr, w2))
    y = np.where(inL, yl, np.where(inR, yr, y2))
    z = np.where(inL, zl, np.where(inR, zr, z2))
    ms = rs * us
    Ff = (F0, F1, ms * v - bx * y, ms * w - bx * z, y * us - bx * v, z * us - bx * w)
    left, right = SL >= 0.0, SR <= 0.0
    return np.array([np.where(left, FL[k], np.where(right, FR[k], Ff[k])) for k in range(6)])


def mc_slope(W, bc):
    if bc == 'periodic':
        dm = W - np.roll(W, 1, axis=1)
        dp = np.roll(W, -1, axis=1) - W
    else:
        Wg = np.concatenate([W[:, :1], W, W[:, -1:]], axis=1)
        dm = Wg[:, 1:-1] - Wg[:, :-2]
        dp = Wg[:, 2:] - Wg[:, 1:-1]
    a = np.minimum(np.minimum(np.abs(2.0 * dm), np.abs(2.0 * dp)), 0.5 * np.abs(dm + dp))
    return np.where(dm * dp > 0.0, np.sign(dm) * a, 0.0)


def prim(U):
    r = U[0]
    return np.array([r, U[1] / r, U[2] / r, U[3] / r, U[4], U[5]])


def hll(WL, WR, bx, c):
    """HLL fluxes (6, n) with the same outer speeds as hlld(): the comparison of the numerical dissipation."""
    UL, UR = mx.cons(WL), mx.cons(WR)
    FL, FR = np.asarray(mx.flux(WL, bx, c)), np.asarray(mx.flux(WR, bx, c))
    cL, cR = cfast(WL[0], WL[4], WL[5], bx, c), cfast(WR[0], WR[4], WR[5], bx, c)
    SL = np.minimum(WL[1] - cL, WR[1] - cR)
    SR = np.maximum(WL[1] + cL, WR[1] + cR)
    F = (SR * FL - SL * FR + SL * SR * (UR - UL)) / (SR - SL)
    return np.where(SL >= 0.0, FL, np.where(SR <= 0.0, FR, F))


SOLVERS = {'hlld': hlld, 'hll': hll}


def rhs(U, dx, bx, c, bc, solver='hlld'):
    W = prim(U)
    sl = mc_slope(W, bc)
    Wm, Wp = W - 0.5 * sl, W + 0.5 * sl
    flux = SOLVERS[solver]
    if bc == 'periodic':
        F = flux(Wp, np.roll(Wm, -1, axis=1), bx, c)           # F[:, i] at x_{i+1/2}
        return -(F - np.roll(F, 1, axis=1)) / dx
    # outflow: ghost states equal to the boundary cells
    WL = np.concatenate([W[:, :1], Wp], axis=1)
    WR = np.concatenate([Wm, W[:, -1:]], axis=1)
    F = flux(WL, WR, bx, c)                                     # F[:, i] at x_{i-1/2}, i = 0..N
    return -(F[:, 1:] - F[:, :-1]) / dx


def run(W0, bx, L, t_end, c=1.0, cfl=0.4, bc='periodic', solver='hlld', every=0, callback=None):
    """Advance the primitive state W0 (6, N) to t_end; returns (W (6, N), number of steps)."""
    W0 = np.asarray(W0, dtype=float)
    N = W0.shape[1]
    dx = L / N
    U = mx.cons(W0)
    t, n = 0.0, 0
    while t < t_end * (1.0 - 1e-14):
        W = prim(U)
        if not np.all(np.isfinite(U)) or np.min(W[0]) <= 0.0:
            raise FloatingPointError('fv_mhd1d: non-finite or non-positive rho at t = %g (step %d)' % (t, n))
        vmax = np.max(np.abs(W[1]) + cfast(W[0], W[4], W[5], bx, c))
        dt = min(cfl * dx / vmax, t_end - t)
        U1 = U + dt * rhs(U, dx, bx, c, bc, solver)
        U = 0.5 * (U + U1 + dt * rhs(U1, dx, bx, c, bc, solver))
        t += dt
        n += 1
        if every and callback is not None and n % every == 0:
            callback(t, prim(U))
    return prim(U), n


# ------------------------------------------------------------ references
T_REF = (0.025, 0.05, 0.075, 0.1)
REF_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'runs', 'stage3',
                                        'fvref'))


def compound_ic(name, N, L=1.0):
    """riemann_mhd(width = 0) of proto_mhd1d on x_i = i L/N: H = 1/2 at L/4 and 3L/4."""
    WL, WR, bx = mx.COMPOUND[name]
    x = np.arange(N) * (L / N)
    h = mx.step_H(x, 0.25 * L, 0.75 * L, 0.0, L)
    return x, np.array([(1.0 - h) * a + h * b for a, b in zip(WL, WR)]), bx


def reference(name, log2n):
    N = 2 ** log2n
    x, W, bx = compound_ic(name, N)
    os.makedirs(REF_DIR, exist_ok=True)
    t, steps, snaps, t0 = 0.0, 0, [], time.time()
    for tn in T_REF:
        W, n = run(W, bx, 1.0, tn - t)
        t, steps = tn, steps + n
        snaps.append(W.copy())
        print('%s N=%d t=%.3f steps %d  min rho %.4f  %.0f s' % (name, N, t, steps, W[0].min(), time.time() - t0),
              flush=True)
    np.savez_compressed(os.path.join(REF_DIR, '%s_N%d.npz' % (name, N)), t=np.array(T_REF), W=np.array(snaps),
                        bx=bx, steps=steps, wall=time.time() - t0)


# ------------------------------------------------------------ self-tests
def _order(errs, Ns):
    return [math.log(errs[i] / errs[i + 1]) / math.log(Ns[i + 1] / Ns[i]) for i in range(len(errs) - 1)]


def _l1(W, We, dx, comps=range(6)):
    return sum(np.sum(np.abs(W[i] - We[i])) for i in comps) * dx


def _riemann_outflow(WL, WR, N, x0):
    x = (np.arange(N) + 0.5) / N
    H = np.where(x < x0, 0.0, 1.0)
    return x, np.array([WL[i] + (WR[i] - WL[i]) * H for i in range(6)])


if __name__ == '__main__' and len(sys.argv) == 4 and sys.argv[1] == 'ref':
    reference(sys.argv[2], int(sys.argv[3]))
elif __name__ == '__main__':
    t0 = time.time()
    Ns = (256, 512, 1024, 2048) if 'quick' in sys.argv else (512, 1024, 2048, 4096)
    # 1. first-order exactness of an isolated RD: one HLLD flux of a jump that is a rotational discontinuity
    for psi in (0.5 * math.pi, math.pi):
        rho, ux, bxr = 1.3, 0.2, 0.9
        V = ux + bxr / math.sqrt(rho)
        WL = np.array([rho, ux, -math.cos(0.3) / math.sqrt(rho), -math.sin(0.3) / math.sqrt(rho),
                       math.cos(0.3), math.sin(0.3)])[:, None]
        WR = np.array([rho, ux, -math.cos(0.3 + psi) / math.sqrt(rho), -math.sin(0.3 + psi) / math.sqrt(rho),
                       math.cos(0.3 + psi), math.sin(0.3 + psi)])[:, None]
        F = hlld(WL, WR, bxr, 1.0)[:, 0]
        Fe = mx.flux(WR[:, 0], bxr) if V < 0 else mx.flux(WL[:, 0], bxr)
        assert np.max(np.abs(F - np.asarray(Fe))) < 1e-13, (F, Fe)
        assert np.max(np.abs(rhs(np.repeat(mx.cons(WL), 8, axis=1), 0.1, bxr, 1.0, 'periodic'))) < 1e-13
    print('HLLD: an isolated RD (dpsi = pi/2, pi) gets the upwind flux to 1e-13; a uniform state has zero rhs')

    # 2. single coplanar shocks F1, S1, S2 and the intermediate I1 (outflow, the jump at the cell face x0 = 5/16)
    print('\nsingle coplanar shocks (outflow, x0 = 5/16): L1 error of all variables against the exact step')
    x0 = 0.3125
    for name, t_end in (('F1', 0.15), ('S1', 0.3), ('S2', 0.3), ('I1', 0.25)):
        S, Wd, Wu, kind = mx.coplanar_case(name)
        errs = []
        for N in Ns:
            x, W0 = _riemann_outflow(Wd, Wu, N, x0)
            W, n = run(W0, mx.COPLANAR_BX, 1.0, t_end, bc='outflow')
            We = mx.single_wave_sample(S, Wd, Wu)((x - x0) / t_end)
            errs.append(_l1(W, We, 1.0 / N))
        print('  %-3s %-13s S = %.4f  L1 %s  order %s'
              % (name, kind, S, ' '.join('%.2e' % e for e in errs), ' '.join('%.2f' % o for o in _order(errs, Ns))))
    print('  (I1, a 2 -> 3 shock, needs a saddle-saddle connection: the FV solution emits waves and does not '
          'converge to it)')

    # 3. perpendicular Riemann problems against the exact solution (outflow, x0 = 0.5)
    print('\nperpendicular Riemann problems (outflow, x0 = 0.5, t = 0.1)')
    for name in ('P1b', 'P2', 'P3'):
        WL, WR = mx.PERP[name]
        w = mx.perp_waves(WL, WR)
        errs = []
        for N in Ns:
            x, W0 = _riemann_outflow(WL, WR, N, 0.5)
            W, n = run(W0, 0.0, 1.0, 0.1, bc='outflow')
            We = mx.perp_sample((x - 0.5) / 0.1, WL, WR, w=w)
            errs.append(_l1(W, We, 1.0 / N))
        print('  %-4s L1 %s  order %s' % (name, ' '.join('%.2e' % e for e in errs),
                                          ' '.join('%.2f' % o for o in _order(errs, Ns))))

    # 4. C3: two fast and two slow shocks (outflow, x0 = 0.5, t = 0.15): states between the shocks
    print('\nC3 colliding flows (outflow, x0 = 0.5, t = 0.15): plateau states against fast_slow_pair')
    WL, WR = (1.0, 1.0, 0.0, 0.0, 1.0, 0.0), (1.0, -1.0, 0.0, 0.0, 1.0, 0.0)
    pr = mx.fast_slow_pair(WL, WR, 1.0)
    errs = []
    for N in Ns:
        x, W0 = _riemann_outflow(WL, WR, N, 0.5)
        W, n = run(W0, 1.0, 1.0, 0.15, bc='outflow')
        We = mx.chain_sample(pr['S'], pr['W'])((x - 0.5) / 0.15)
        errs.append(_l1(W, We, 1.0 / N))
        xi = (x - 0.5) / 0.15
        Sf, Ss = pr['S'][3], pr['S'][2]
        mid = (xi > Ss + 0.25 * (Sf - Ss)) & (xi < Sf - 0.25 * (Sf - Ss))
        cen = np.abs(xi) < 0.5 * Ss
        line = '  N=%5d  fast-down plateau err %.1e  centre err %.1e' % (
            N, max(np.max(np.abs(W[i][mid] - pr['W'][3][i])) for i in (0, 1, 2, 4)),
            max(np.max(np.abs(W[i][cen] - pr['W'][2][i])) for i in (0, 1, 2, 4)))
        print(line)
    print('  L1 %s  order %s' % (' '.join('%.2e' % e for e in errs), ' '.join('%.2f' % o for o in _order(errs, Ns))))

    # 5. RDs transported one period (periodic, rho = B_x = |B_t| = 1, u_x = 0, t = 1): L1 of (u_t, B_t)
    print('\nRD (periodic, rho = B_x = |B_t| = 1, u_x = 0, t = 1): L1 of (u_t, B_t)')
    for dpsi, w in ((0.5 * math.pi, 0.0), (math.pi, 0.0), (math.pi, 0.02)):
        errs = []
        for N in Ns[:3]:
            x = np.arange(N) / N
            W0 = mx.rd_profile(x, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, dpsi, 0.25, 0.75, w, 1.0)
            W, n = run(W0, 1.0, 1.0, 1.0)
            errs.append(_l1(W, W0, 1.0 / N, comps=(2, 3, 4, 5)))
            dr = np.max(np.abs(W[0] - 1.0))
        print('  dpsi = %.2f width %.2f  L1 %s  order %s  (N = %d: max|rho - 1| %.1e)'
              % (dpsi, w, ' '.join('%.2e' % e for e in errs), ' '.join('%.2f' % o for o in _order(errs, Ns[:3])),
                 Ns[2], dr))
    print('  (the sharp coplanar reversal keeps B_z = u_z = 0 exactly and cannot rotate: it splits into other waves)')
    print('\nelapsed %.0f s' % (time.time() - t0))
