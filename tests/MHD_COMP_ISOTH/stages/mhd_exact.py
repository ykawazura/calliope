#!/usr/bin/env python3
"""Exact solutions of 1D isothermal MHD (sound speed c, B_x constant) for the Stage 3 tests.

State W = (rho, u_x, u_y, u_z, B_y, B_z) (first axis), conserved U = (rho, rho u, B_y, B_z), flux along x
    F = (rho u_x, rho u_x^2 + c^2 rho + B^2/2 - B_x^2, rho u_x u_t - B_x B_t, u_x B_t - B_x u_t)   (t = y, z).
Characteristic speeds: c_A = |B_x|/sqrt(rho), c_f,s^2 = [c^2 + v_A^2 +- sqrt((c^2 + v_A^2)^2 - 4 c^2 c_A^2)]/2.

Perpendicular MHD (B_x = 0): b = |B_t|/rho is advected and continuous across fast waves, so each side is a
barotropic gas with P = c^2 rho + b^2 rho^2/2, c_f^2 = c^2 + b^2 rho; u_t and the direction of B_t are
constant on each side, and the tangential discontinuity (TD) at u* separates the two sides (P, u_x continuous):
    shock        f_K(P) = sqrt((P - P_K)(1/rho_K - 1/rho)),   mass flux j_K = sqrt((P - P_K)/(1/rho_K - 1/rho))
    rarefaction  f_K(P) = Phi(rho; rho_K) = int_{rho_K}^{rho} c_f/rho' drho'
                        = 2(s - s_K) + c ln(rho/rho_K) - 2c ln((s + c)/(s_K + c)),   s = c_f(rho)
    u* = u_L - f_L(P*) = u_R + f_R(P*) (bisection in ln P*),  S_L = u_L - j_L/rho_L,  S_R = u_R + j_R/rho_R,
    fans: u = u_L - Phi(rho; rho_L), x/t = u - c_f (left);  u = u_R + Phi(rho; rho_R), x/t = u + c_f (right).
Coplanar shocks (u_z = B_z = 0, B_y1 != 0): in the shock frame w = u_x - S, j = rho w; with X = w_1^2,
a^2 = B_x^2/rho_1, b_1^2 = B_y1^2/rho_1 and the compression r = rho_2/rho_1,
    B_y2/B_y1 = r (X - a^2)/(X - r a^2),   u_y2 = u_y1 + B_x (B_y2 - B_y1)/j,
    (X - r c^2)(X - r a^2)^2 = (r b_1^2/2) X ((r + 1) X - 2 r a^2)   (normal momentum);
the roots X > 0 are the slow, intermediate and fast shocks of compression r (classify()).
Alfven waves: with rho, u_x and |B_t| uniform the tangential equations are linear, so any profile with
u_t = -sigma B_t/sqrt(rho) + const translates at V = u_x + sigma B_x/sqrt(rho) (circularly polarized wave,
rotational discontinuity (RD)). A TD with c^2 rho + B_t^2/2 uniform and u = 0 is steady.
"""
import math
import numpy as np


# ------------------------------------------------------------ fluxes and speeds
def cons(W):
    rho, ux, uy, uz, by, bz = (np.asarray(v, dtype=float) for v in W)
    return np.array([rho, rho * ux, rho * uy, rho * uz, by, bz])


def prim(U):
    rho = np.asarray(U[0], dtype=float)
    return np.array([rho, U[1] / rho, U[2] / rho, U[3] / rho, U[4], U[5]])


def flux(W, bx, c=1.0):
    rho, ux, uy, uz, by, bz = (np.asarray(v, dtype=float) for v in W)
    pt = c * c * rho + 0.5 * (bx * bx + by * by + bz * bz)
    return np.array([rho * ux, rho * ux * ux + pt - bx * bx, rho * ux * uy - bx * by, rho * ux * uz - bx * bz,
                     ux * by - bx * uy, ux * bz - bx * uz])


def speeds(W, bx, c=1.0):
    """(c_s, c_A, c_f) along x; c_s = c c_A/c_f avoids the cancellation of the difference form."""
    rho, by, bz = (np.asarray(W[i], dtype=float) for i in (0, 4, 5))
    ca2 = bx * bx / rho
    va2 = ca2 + (by * by + bz * bz) / rho
    cf = np.sqrt(0.5 * (c * c + va2 + np.sqrt(np.maximum(0.0, (c * c + va2) ** 2 - 4.0 * c * c * ca2))))
    return c * np.sqrt(ca2) / cf, np.sqrt(ca2), cf


def rh_residual(W1, W2, S, bx, c=1.0):
    """F(W2) - F(W1) - S (U2 - U1): zero for a discontinuity of speed S."""
    return flux(W2, bx, c) - flux(W1, bx, c) - S * (cons(W2) - cons(W1))


# ------------------------------------------------------------ perpendicular MHD (B_x = 0)
def perp_P(rho, b, c=1.0):
    return c * c * rho + 0.5 * b * b * rho * rho


def perp_rho(P, b, c=1.0):
    return 2.0 * P / (c * c + np.sqrt(c ** 4 + 2.0 * b * b * P))


def perp_cf(rho, b, c=1.0):
    return np.sqrt(c * c + b * b * rho)


def perp_Phi(rho, rho_K, b, c=1.0):
    """int_{rho_K}^{rho} c_f/rho' drho' in the difference form (c ln(rho/rho_K) for b -> 0)."""
    s, sK = perp_cf(rho, b, c), perp_cf(rho_K, b, c)
    ds = b * b * (rho - rho_K) / (s + sK)
    return 2.0 * ds + c * np.log(rho / rho_K) - 2.0 * c * np.log1p(ds / (sK + c))


def perp_f(P, rho_K, b, c=1.0):
    """(f_K(P), mass flux j_K) of one side; j_K = 0 for a rarefaction."""
    PK = perp_P(rho_K, b, c)
    rho = perp_rho(P, b, c)
    if P > PK:
        dtau = (rho - rho_K) / (rho * rho_K)
        return math.sqrt((P - PK) * dtau), math.sqrt((P - PK) / dtau)
    return float(perp_Phi(rho, rho_K, b, c)), 0.0


def _b(W):
    return math.hypot(W[4], W[5]) / W[0]


def perp_waves(WL, WR, c=1.0):
    """Wave structure of the perpendicular Riemann problem: dict with P*, u*, the star states W*_L, W*_R
    and, for each side ('L', 'R'), the type and speeds (as riemann_exact.waves)."""
    WL, WR = [float(v) for v in WL], [float(v) for v in WR]
    bL, bR = _b(WL), _b(WR)
    PL, PR = perp_P(WL[0], bL, c), perp_P(WR[0], bR, c)
    g = lambda q: perp_f(math.exp(q), WL[0], bL, c)[0] + perp_f(math.exp(q), WR[0], bR, c)[0] + WR[1] - WL[1]
    lo, hi = math.log(min(PL, PR)) - 50.0, math.log(max(PL, PR)) + 50.0
    assert g(lo) < 0.0 < g(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    Ps = math.exp(0.5 * (lo + hi))
    fL, jL = perp_f(Ps, WL[0], bL, c)
    fR, jR = perp_f(Ps, WR[0], bR, c)
    us = 0.5 * ((WL[1] - fL) + (WR[1] + fR))
    w = dict(P_s=Ps, u_s=us, bL=bL, bR=bR)
    for side, W, b, j, sgn in (('L', WL, bL, jL, -1.0), ('R', WR, bR, jR, 1.0)):
        rs = float(perp_rho(Ps, b, c))
        q = rs / W[0]
        w['W' + side] = (rs, us, W[2], W[3], W[4] * q, W[5] * q)
        if Ps > perp_P(W[0], b, c):
            w[side] = dict(type='shock', speed=W[1] + sgn * j / W[0], j=j, r=q)
        else:
            w[side] = dict(type='rarefaction', head=W[1] + sgn * float(perp_cf(W[0], b, c)),
                           tail=us + sgn * float(perp_cf(rs, b, c)))
    return w


def fan_edges(w):
    """Leftmost and rightmost wave speeds of a Riemann fan."""
    a = w['L']['speed'] if w['L']['type'] == 'shock' else w['L']['head']
    b = w['R']['speed'] if w['R']['type'] == 'shock' else w['R']['head']
    return a, b


def perp_sample(xi, WL, WR, c=1.0, w=None):
    """W (6, n) at the similarity coordinates xi = x/t (rarefactions by bisection in ln rho)."""
    if w is None:
        w = perp_waves(WL, WR, c)
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    out = np.empty((6,) + xi.shape)
    us = w['u_s']
    for side, W, b, sgn in (('L', WL, w['bL'], -1.0), ('R', WR, w['bR'], 1.0)):
        sel = xi <= us if side == 'L' else xi > us
        far = (lambda e: xi < e) if side == 'L' else (lambda e: xi > e)
        Ws = np.asarray(w['W' + side], dtype=float)
        wv = w[side]
        outer = sel & far(wv['speed'] if wv['type'] == 'shock' else wv['head'])
        out[:, outer] = np.asarray(W, dtype=float)[:, None]
        if wv['type'] == 'shock':
            out[:, sel & ~outer] = Ws[:, None]
            continue
        fan = sel & ~outer & far(wv['tail'])
        star = sel & ~outer & ~fan
        out[:, star] = Ws[:, None]
        x = xi[fan]
        lo, hi = np.full(x.shape, math.log(Ws[0])), np.full(x.shape, math.log(W[0]))
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            r = np.exp(mid)
            h = W[1] + sgn * perp_Phi(r, W[0], b, c) + sgn * perp_cf(r, b, c) - x
            # left fan: h = u - c_f - xi decreases with rho; right fan: h = u + c_f - xi increases
            up = (h > 0.0) if side == 'L' else (h < 0.0)
            lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
        r = np.exp(0.5 * (lo + hi))
        q = r / W[0]
        out[:, fan] = np.array([r, W[1] + sgn * perp_Phi(r, W[0], b, c), np.full_like(r, W[2]),
                                np.full_like(r, W[3]), W[4] * q, W[5] * q])
    return out


# ------------------------------------------------------------ coplanar shocks
def classify(W1, W2, S, bx, c=1.0, tol=1e-9):
    """(i1, i2, kind): the states 1: |w| > c_f, 2: c_A < |w| < c_f, 3: c_s < |w| < c_A, 4: |w| < c_s
    upstream and downstream in the shock frame; kind = 'fast' (1->2), 'slow' (3->4), 'intermediate',
    with '(switch-off)' or '(switch-on)' appended when |w| = c_A on one side within tol."""
    idx = []
    for W in (W1, W2):
        cs, ca, cf = (float(v) for v in speeds(W, bx, c))
        a = abs(W[1] - S)
        idx.append(1 if a > cf else 2 if a > ca else 3 if a > cs else 4)
    i1, i2 = idx
    kind = 'fast' if (i1, i2) == (1, 2) else 'slow' if (i1, i2) == (3, 4) else 'intermediate'
    ca1, ca2 = float(speeds(W1, bx, c)[1]), float(speeds(W2, bx, c)[1])
    if abs(abs(W1[1] - S) - ca1) <= tol * ca1:
        kind += ' (switch-off)' if abs(W2[4]) <= tol * (abs(W1[4]) + abs(bx)) else ' (|w1| = c_A1)'
    if abs(abs(W2[1] - S) - ca2) <= tol * ca2:
        kind += ' (switch-on)' if abs(W1[4]) <= tol * (abs(W2[4]) + abs(bx)) else ' (|w2| = c_A2)'
    return i1, i2, kind


def coplanar_shocks(W1, bx, r, c=1.0, side=1):
    """Shocks of compression r = rho_2/rho_1 > 1 into the upstream state W1 (u_z = B_z = 0, B_y1 != 0).
    side = +1: the shock moves to +x into W1 (upstream on the right, w_1 = u_x1 - S < 0);
    side = -1: to -x (upstream on the left, w_1 > 0).
    Returns [(X, S, W2, kind)] sorted by X = w_1^2."""
    rho1, ux1, uy1, uz1, by1, bz1 = (float(v) for v in W1)
    assert by1 != 0.0 and uz1 == 0.0 and bz1 == 0.0 and r > 1.0
    a2, b2 = bx * bx / rho1, by1 * by1 / rho1
    A, Bq = r * c * c, r * a2
    p = np.array([1.0, -(2.0 * Bq + A) - 0.5 * r * (r + 1.0) * b2, Bq * Bq + 2.0 * A * Bq + r * r * a2 * b2,
                  -A * Bq * Bq])
    dp = np.polyder(p)
    out = []
    for z in np.roots(p):
        if abs(z.imag) > 1e-8 * abs(z) or z.real <= 0.0:
            continue
        X = z.real
        for _ in range(4):
            d = np.polyval(dp, X)
            if d != 0.0:
                X -= np.polyval(p, X) / d
        w1 = -side * math.sqrt(X)
        S = ux1 - w1
        j = rho1 * w1
        by2 = by1 * r * (X - a2) / (X - Bq)
        W2 = (r * rho1, S + w1 / r, uy1 + bx * (by2 - by1) / j, 0.0, by2, 0.0)
        out.append((X, S, W2, classify(W1, W2, S, bx, c)[2]))
    return sorted(out, key=lambda e: e[0])


def pick(shocks, kind):
    """The shock of the given kind ('fast', 'slow', 'intermediate') from coplanar_shocks()."""
    sel = [e for e in shocks if e[3].split()[0] == kind]
    assert len(sel) == 1, (kind, [e[3] for e in shocks])
    return sel[0]


def fast_slow_pair(WL, WR, bx, c=1.0, guess=(1.5, 1.4, 1.5, 1.4), tol=1e-14):
    """Riemann problem with a fast and a slow shock on each side and no rotation (u_z = B_z = 0):
    Newton on (r_fL, r_sL, r_fR, r_sR) so that the two central states agree in (rho, u_x, u_y, B_y).
    Returns dict(r=..., S=(S_fL, S_sL, S_sR, S_fR), W=(W_L, W_fL, W_c, W_fR, W_R)) (W_c from the left)."""
    def chain(q):
        rfL, rsL, rfR, rsR = q
        _, SfL, WfL, _ = pick(coplanar_shocks(WL, bx, rfL, c, -1), 'fast')
        _, SsL, WcL, _ = pick(coplanar_shocks(WfL, bx, rsL, c, -1), 'slow')
        _, SfR, WfR, _ = pick(coplanar_shocks(WR, bx, rfR, c, +1), 'fast')
        _, SsR, WcR, _ = pick(coplanar_shocks(WfR, bx, rsR, c, +1), 'slow')
        res = np.array([WcL[i] - WcR[i] for i in (0, 1, 2, 4)])
        return res, (SfL, SsL, SsR, SfR), (WfL, WcL, WcR, WfR)
    q = np.array(guess, dtype=float)
    for it in range(60):
        res, S, Ws = chain(q)
        if np.max(np.abs(res)) < tol:
            break
        J = np.empty((4, 4))
        for k in range(4):
            h = 1e-7 * q[k]
            dq = q.copy()
            dq[k] += h
            J[:, k] = (chain(dq)[0] - res) / h
        q = q - np.linalg.solve(J, res)
    res, S, (WfL, WcL, WcR, WfR) = chain(q)
    assert np.max(np.abs(res)) < 1e-12, res
    assert S[0] < S[1] < S[2] < S[3], S
    return dict(r=tuple(q), S=S, W=(tuple(WL), WfL, WcL, WfR, tuple(WR)), res=np.max(np.abs(res)), WcR=WcR)


def single_wave_sample(S, WL, WR):
    """sample(xi) of a single discontinuity of speed S between WL and WR."""
    WL, WR = np.asarray(WL, dtype=float), np.asarray(WR, dtype=float)
    return lambda xi: np.where(np.atleast_1d(xi)[None, :] < S, WL[:, None], WR[:, None])


def chain_sample(S, Ws):
    """sample(xi) of discontinuities at the speeds S[0] < S[1] < ... separating the states Ws[0], Ws[1], ..."""
    Ws = np.asarray(Ws, dtype=float)
    return lambda xi: Ws[np.searchsorted(np.asarray(S), np.atleast_1d(xi), side='right')].T


# ------------------------------------------------------------ periodic superposition
def periodic(x, t, xa, xb, L, WL, WR, sol_a, sol_b):
    """Solution of the periodic step, R on (xa, xb) and L elsewhere, at time t > 0: the solutions of the
    interfaces (L | R) at xa and (R | L) at xb superposed. sol_a, sol_b = (sample, (a, b)): sample(xi) -> W
    (6, n) at xi = (x - x0)/t, a <= b the fan edges; sample = None marks an unknown fan bounded by (a, b),
    whose points are invalid. Returns (W, valid): valid is False everywhere once the fans overlap."""
    x = np.mod(np.asarray(x, dtype=float), L)
    inside = (x > xa) & (x < xb)
    W = np.where(inside[None, :], np.asarray(WR, dtype=float)[:, None], np.asarray(WL, dtype=float)[:, None])
    valid = np.ones(x.shape, dtype=bool)
    edges = []
    for x0, (smp, (a, b)) in ((xa, sol_a), (xb, sol_b)):
        assert max(-a, b) * t < 0.5 * L, 'fan wider than L/2 on one side'
        xi = (np.mod(x - x0 + 0.5 * L, L) - 0.5 * L) / t
        sel = (xi >= a) & (xi <= b)
        if smp is None:
            valid[sel] = False
        elif np.any(sel):
            W[:, sel] = smp(xi[sel])
        edges.append((x0 + a * t, x0 + b * t))
    (a1, b1), (a2, b2) = edges
    if not ((b1 < a2) and (b2 < a1 + L)):
        valid[:] = False
    return W, valid


def char_bounds(WL, WR, bx, c=1.0):
    """(min(u_x - c_f), max(u_x + c_f)) over the two states: the fan edges of the Riemann problem when its
    outermost waves are rarefactions (a shock moves faster; check with the FV reference)."""
    lam = []
    for W in (WL, WR):
        cf = float(speeds(W, bx, c)[2])
        lam += [W[1] - cf, W[1] + cf]
    return min(lam), max(lam)


# ------------------------------------------------------------ smooth profiles and linearly degenerate waves
def step_H(x, xa, xb, w, L):
    """Periodic smoothed step: 1 on (xa, xb), 0 outside, tanh of width w (length) at both interfaces;
    H = [1 + tanh(g/(w (2 pi/L) sin(pi l/L)))]/2, g = cos(2 pi (x - x_mid)/L) - cos(pi l/L), l = xb - xa,
    so that H ~ [1 + tanh((x - xa)/w)]/2 near xa. w = 0: the step with H = 1/2 at xa and xb."""
    l, xm = xb - xa, 0.5 * (xa + xb)
    g = np.cos(2.0 * math.pi * (np.asarray(x, dtype=float) - xm) / L) - math.cos(math.pi * l / L)
    if w == 0.0:
        return np.where(np.abs(g) < 1e-14, 0.5, (g > 0.0).astype(float))
    return 0.5 * (1.0 + np.tanh(g / (w * 2.0 * math.pi / L * math.sin(math.pi * l / L))))


def rd_profile(x, t, rho0, ux0, bx, bt, psi_L, dpsi, xa, xb, w, L, sigma=1.0):
    """Alfven waves of angle psi = psi_L + dpsi H(x - V t) (step_H), B_t = bt (cos psi, sin psi),
    u_t = -sigma B_t/sqrt(rho0), V = ux0 + sigma bx/sqrt(rho0): exact for any w (two RDs for w = 0)."""
    V = ux0 + sigma * bx / math.sqrt(rho0)
    psi = psi_L + dpsi * step_H(np.asarray(x, dtype=float) - V * t, xa, xb, w, L)
    by, bz = bt * np.cos(psi), bt * np.sin(psi)
    one = np.ones_like(psi)
    return np.array([rho0 * one, ux0 * one, -sigma * by / math.sqrt(rho0), -sigma * bz / math.sqrt(rho0), by, bz])


def cp_alfven(x, t, rho0, ux0, bx, bt, mode, L, sigma=1.0, psi0=0.0):
    """Circularly polarized Alfven wave psi = 2 pi mode (x - V t)/L + psi0 (exact at any amplitude bt)."""
    V = ux0 + sigma * bx / math.sqrt(rho0)
    psi = 2.0 * math.pi * mode * (np.asarray(x, dtype=float) - V * t) / L + psi0
    by, bz = bt * np.cos(psi), bt * np.sin(psi)
    one = np.ones_like(psi)
    return np.array([rho0 * one, ux0 * one, -sigma * by / math.sqrt(rho0), -sigma * bz / math.sqrt(rho0), by, bz])


def linear_mode(x, t, rho0, bx, by0, family, amp, mode, L, c=1.0):
    """Right-going fast or slow eigenmode on (rho0, u = 0, B_x, B_y0): delta rho = amp rho0 cos(k x - omega t),
    delta u_x = v delta rho/rho0, delta B_y = (v^2 - c^2) delta rho/B_y0, delta u_y = -B_x delta B_y/(rho0 v),
    v = c_f or c_s. Returns (W, omega)."""
    cs, ca, cf = (float(v) for v in speeds((rho0, 0, 0, 0, by0, 0), bx, c))
    v = cf if family == 'fast' else cs
    k = 2.0 * math.pi * mode / L
    dr = amp * rho0 * np.cos(k * np.asarray(x, dtype=float) - v * k * t)
    dby = (v * v - c * c) * dr / by0
    one = np.ones_like(dr)
    return np.array([rho0 + dr, v * dr / rho0, -bx * dby / (rho0 * v), 0.0 * one, by0 + dby, 0.0 * one]), v * k


def td_density(by, bz, P0, c=1.0):
    """rho of the steady TD c^2 rho + B_t^2/2 = P0 (u = 0)."""
    rho = (P0 - 0.5 * (np.asarray(by) ** 2 + np.asarray(bz) ** 2)) / (c * c)
    assert np.all(rho > 0.0)
    return rho


# ------------------------------------------------------------ test problems (c = 1)
# perpendicular Riemann problems (B_x = 0): W = (rho, u_x, u_y, u_z, B_y, B_z)
PERP = {'P1a': ((1.0, 1.0, 0, 0, 1.0, 0), (1.0, -1.0, 0, 0, 1.0, 0)),
        'P1b': ((1.0, 2.0, 0, 0, 1.0, 0), (1.0, -2.0, 0, 0, 1.0, 0)),
        'P1c': ((1.0, 4.0, 0, 0, 1.0, 0), (1.0, -4.0, 0, 0, 1.0, 0)),
        'P2': ((10.0, 0.0, 0, 0, 1.0, 0), (1.0, 0.0, 0, 0, 1.0, 0)),
        'P3': ((1.0, 0.0, 0, 0, 1.0, 0), (0.125, 0.0, 0, 0, -1.0, 0))}
# coplanar single shocks into the upstream state (1, 0, 0, 0, 1, 0), B_x = 1: (r, kind)
COPLANAR_UP, COPLANAR_BX = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0), 1.0
COPLANAR = {'F1': (2.0, 'fast'), 'F2': (4.0, 'fast'), 'S1': (1.5, 'slow'), 'S2': (1.9, 'slow'),
            'I1': (1.5, 'intermediate')}
# Riemann problems with fast and slow waves (W_L, W_R, B_x); no exact fan: references from fv_mhd1d.py.
# C2 is C1 with B_y reversed (the isothermal analogue of Brio & Wu 1988); C3 is exact on the colliding side
COMPOUND = {'C1': ((1.0, 0.0, 0.0, 0.0, 1.0, 0.0), (0.125, 0.0, 0.0, 0.0, 1.0, 0.0), 0.75),
            'C2': ((1.0, 0.0, 0.0, 0.0, 1.0, 0.0), (0.125, 0.0, 0.0, 0.0, -1.0, 0.0), 0.75),
            'C3': ((1.0, 1.0, 0.0, 0.0, 1.0, 0.0), (1.0, -1.0, 0.0, 0.0, 1.0, 0.0), 1.0)}


def coplanar_case(name, side=1, c=1.0):
    """(S, W_down, W_up, kind) of a coplanar test shock moving to +x (side = 1)."""
    r, kind = COPLANAR[name]
    X, S, W2, k = pick(coplanar_shocks(COPLANAR_UP, COPLANAR_BX, r, c, side), kind)
    return S, W2, COPLANAR_UP, k


if __name__ == '__main__':
    import riemann_exact as re
    ok = lambda a, b, tol: abs(a - b) <= tol * max(1.0, abs(b))
    # 1. b -> 0: the isothermal Euler solution
    for name, (rL, uL, rR, uR) in {'U1': (1, 1, 1, -1), 'U2': (1, 2, 1, -2), 'U4': (1, 4, 1, -4),
                                   'Sod': (10, 0, 1, 0), 'recede': (1, -2, 1, 2)}.items():
        we = re.waves(rL, uL, rR, uR)
        for bt in (0.0, 1e-9):
            w = perp_waves((rL, uL, 0, 0, bt * rL, 0), (rR, uR, 0, 0, bt * rR, 0))
            assert ok(w['WL'][0], we['rho_s'], 1e-9) and ok(w['u_s'], we['u_s'], 1e-9), (name, bt)
            assert ok(fan_edges(w)[1], fan_edges(we)[1], 1e-9) and ok(fan_edges(w)[0], fan_edges(we)[0], 1e-9)
            xi = np.linspace(-5.0, 5.0, 2001)
            rho_e, u_e = re.sample(xi, rL, uL, rR, uR, w=we)
            Wp = perp_sample(xi, (rL, uL, 0, 0, bt * rL, 0), (rR, uR, 0, 0, bt * rR, 0), w=w)
            assert np.max(np.abs(Wp[0] - rho_e)) < 1e-8 and np.max(np.abs(Wp[1] - u_e)) < 1e-8, name
    print('b -> 0: perp_waves/perp_sample = riemann_exact (U1, U2, U4, Sod, receding) to 1e-9')
    # 2. Phi against Gauss-Legendre quadrature and its b -> 0 limit
    xg, wg = np.polynomial.legendre.leggauss(64)
    for b, r0, r1 in ((1.0, 1.0, 0.3), (2.0, 5.0, 0.01), (1e-6, 1.0, 0.5), (0.5, 0.2, 3.0)):
        # int c_f dln(rho): smooth in ln(rho)
        q0, q1 = math.log(r0), math.log(r1)
        quad = 0.5 * (q1 - q0) * np.sum(wg * perp_cf(np.exp(0.5 * (q1 - q0) * xg + 0.5 * (q1 + q0)), b))
        assert abs(perp_Phi(r1, r0, b) - quad) < 1e-12 * max(1.0, abs(quad)), (b, r0, r1)
    print('Phi(rho; rho_K) = quadrature of c_f/rho to 1e-12')
    # 3. P1b closed form: rho* = 3, P* = 7.5, j = 3, S = +-1, u* = 0
    w = perp_waves(*PERP['P1b'])
    assert ok(w['WL'][0], 3.0, 1e-13) and ok(w['P_s'], 7.5, 1e-13) and ok(w['R']['j'], 3.0, 1e-12)
    assert ok(w['R']['speed'], 1.0, 1e-12) and ok(w['L']['speed'], -1.0, 1e-12) and abs(w['u_s']) < 1e-13
    print('P1b: rho* = 3, P* = 7.5, j = 3, S = +-1 to 1e-12')
    # 4. perpendicular problems: RH residuals, fan ODE, TD balance
    print('\nperpendicular (B_x = 0)')
    for name, (WL, WR) in PERP.items():
        w = perp_waves(WL, WR)
        res = []
        for side, W in (('L', WL), ('R', WR)):
            if w[side]['type'] == 'shock':
                res.append(np.max(np.abs(rh_residual(W, w['W' + side], w[side]['speed'], 0.0))))
            else:
                # self-similar fan: dF/dxi = xi dU/dxi (central differences, O(h^2))
                a, b = sorted((w[side]['head'], w[side]['tail']))
                xi = np.linspace(a, b, 41)[1:-1]
                h = 1e-5
                Wp, Wm = perp_sample(xi + h, WL, WR, w=w), perp_sample(xi - h, WL, WR, w=w)
                ode = (flux(Wp, 0.0) - flux(Wm, 0.0) - xi * (cons(Wp) - cons(Wm))) / (2 * h)
                res.append(np.max(np.abs(ode)))
        WsL, WsR = w['WL'], w['WR']
        td = rh_residual(WsL, WsR, w['u_s'], 0.0)
        res.append(np.max(np.abs(td)))
        assert max(res) < 1e-8, (name, res)
        line = '%-4s P* = %.6f u* = %+.6f  rho*_L = %.6f B*_L = %+.6f  rho*_R = %.6f B*_R = %+.6f' % (
            name, w['P_s'], w['u_s'], WsL[0], WsL[4], WsR[0], WsR[4])
        for side in ('L', 'R'):
            v = w[side]
            line += ('  %s shock S = %+.6f' % (side, v['speed']) if v['type'] == 'shock' else
                     '  %s fan [%+.6f, %+.6f]' % (side, v['head'], v['tail']))
        print(line + '  (max residual %.1e)' % max(res))
    # 5. coplanar shocks: F1 closed form, RH residuals, classification
    S, W2, W1, k = coplanar_case('F1')
    X = (S - W1[1]) ** 2
    assert ok(X, 4.0 + 2.0 * math.sqrt(2.0), 1e-13) and ok(W2[4], 1.0 + math.sqrt(2.0), 1e-12)
    roots = [e[0] for e in coplanar_shocks(COPLANAR_UP, 1.0, 2.0)]
    assert np.allclose(roots, [1.0, 4.0 - 2.0 * math.sqrt(2.0), 4.0 + 2.0 * math.sqrt(2.0)], rtol=1e-13)
    print('\nF1: X^3 - 9X^2 + 16X - 8 = (X - 1)(X^2 - 8X + 8), X_fast = 4 + 2 sqrt(2), B_y2 = 1 + sqrt(2) to 1e-12')
    print('coplanar shocks into (rho, u_x, u_y, B_y) = (1, 0, 0, 1), B_x = 1 (c_s, c_A, c_f = %.6f, %.6f, %.6f)'
          % tuple(float(v) for v in speeds(COPLANAR_UP, 1.0)))
    for name in COPLANAR:
        S, W2, W1, k = coplanar_case(name)
        res = np.max(np.abs(rh_residual(W1, W2, S, COPLANAR_BX)))
        assert res < 1e-12, (name, res)
        i1, i2, _ = classify(W1, W2, S, COPLANAR_BX)
        print('%-3s r = %.2f  S = %.6f  down (%.6f, %.6f, %+.6f, %+.6f)  %d->%d %s  (RH residual %.1e)'
              % (name, COPLANAR[name][0], S, W2[0], W2[1], W2[2], W2[4], i1, i2, k, res))
    for r in (2.0,):
        for X, S, W2, k in coplanar_shocks(COPLANAR_UP, 1.0, r):
            print('    r = %.1f: X = %.6f S = %.6f down (%.6f, %.6f, %+.6f, %+.6f) %s'
                  % (r, X, S, W2[0], W2[1], W2[2], W2[4], k))
    # mirror symmetry (x, u_x, B_y) -> (-x, -u_x, -B_y): the left-moving shock into the mirrored upstream state
    S, W2, W1, k = coplanar_case('F1')
    _, Sm, W2m, km = pick(coplanar_shocks((1.0, 0.0, 0.0, 0.0, -1.0, 0.0), 1.0, 2.0, side=-1), 'fast')
    assert ok(Sm, -S, 1e-13) and np.allclose(W2m, (W2[0], -W2[1], W2[2], 0.0, -W2[4], 0.0), atol=1e-13)
    # 6. C3: colliding (1, +-1, 0, 1), B_x = 1: fast and slow shocks on each side
    WL, WR = (1.0, 1.0, 0.0, 0.0, 1.0, 0.0), (1.0, -1.0, 0.0, 0.0, 1.0, 0.0)
    p = fast_slow_pair(WL, WR, 1.0)
    for i, (Wa, Wb) in enumerate(zip(p['W'][:2], p['W'][1:3])):
        assert np.max(np.abs(rh_residual(Wa, Wb, p['S'][i], 1.0))) < 1e-12
    for i, (Wa, Wb) in enumerate(zip((p['WcR'], p['W'][3]), p['W'][3:])):
        assert np.max(np.abs(rh_residual(Wa, Wb, p['S'][i + 2], 1.0))) < 1e-12
    # parity (x, u) -> (-x, -u): S_fR = -S_fL, u_c = 0
    assert ok(p['S'][3], -p['S'][0], 1e-11) and abs(p['W'][2][1]) < 1e-12 and abs(p['W'][2][2]) < 1e-12
    print('\nC3: r_f = %.6f r_s = %.6f  S_f = %+.6f S_s = %+.6f  fast down (%.6f, %+.6f, %+.6f, %.6f)'
          '  centre (%.6f, %+.1e, %+.1e, %.6f)  (matching residual %.1e)'
          % (p['r'][2], p['r'][3], p['S'][3], p['S'][2], *[p['W'][3][i] for i in (0, 1, 2, 4)],
             *[p['W'][2][i] for i in (0, 1, 2, 4)], p['res']))
    # 7. Alfven waves and the linear modes: residual of the PDE with spectral derivatives
    L, N = 2.0 * math.pi, 256
    x = np.arange(N) * L / N
    k = np.fft.rfftfreq(N, d=L / N) * 2.0 * math.pi

    def pde_residual(fun, bx, h=1e-4):
        """max |dU/dt + dF/dx| at t = 0 (dU/dt by central differences of fun(t), dF/dx spectrally)."""
        dUdt = (cons(fun(h)) - cons(fun(-h))) / (2.0 * h)
        F = flux(fun(0.0), bx)
        dFdx = np.fft.irfft(1j * k * np.fft.rfft(F, axis=1), n=N, axis=1)
        return np.max(np.abs(dUdt + dFdx))
    r_cp = pde_residual(lambda t: cp_alfven(x, t, 2.0, 0.3, 1.2, 0.8, 2, L, sigma=-1.0), 1.2)
    r_rd = pde_residual(lambda t: rd_profile(x, t, 1.0, 0.0, 1.0, 1.0, 0.0, math.pi, 0.25 * L, 0.75 * L,
                                             8 * L / N, L), 1.0, h=1e-5)
    assert r_cp < 1e-7 and r_rd < 1e-6, (r_cp, r_rd)
    print('Alfven waves (CP mode 2, RD dpsi = pi, width 8 cells): PDE residual %.1e, %.1e (dt by central differences)'
          % (r_cp, r_rd))
    for fam in ('fast', 'slow'):
        amp = 1e-6
        Wl = lambda t: linear_mode(x, t, 1.0, 1.0, 0.5, fam, amp, 1, L)[0]
        r_lin = pde_residual(Wl, 1.0, h=1e-3)
        assert r_lin < 1e-10, (fam, r_lin)
        print('linear %s mode (B_x = 1, B_y0 = 0.5, amp 1e-6): omega = %.6f, PDE residual %.1e (O(amp^2))'
              % (fam, linear_mode(x, 0.0, 1.0, 1.0, 0.5, fam, amp, 1, L)[1], r_lin))
    # 8. step_H: periodic, 1/2 at the interfaces, the tanh profile near xa
    H = step_H(x, 0.2 * L, 0.7 * L, 0.0, L)
    assert H[np.argmin(np.abs(x - 0.25 * L))] == 1.0 and step_H(np.array([0.2 * L]), 0.2 * L, 0.7 * L, 0.0, L)[0] == 0.5
    xs = 0.2 * L + np.linspace(-0.05, 0.05, 11)
    assert np.max(np.abs(step_H(xs, 0.2 * L, 0.7 * L, 0.02, L) - 0.5 * (1 + np.tanh((xs - 0.2 * L) / 0.02)))) < 2e-3
    # 9. periodic superposition of a perpendicular problem
    WL, WR = PERP['P2']
    w = perp_waves(WL, WR)
    w2 = perp_waves(WR, WL)
    tt = 0.3
    Wp, valid = periodic(x, tt, 0.25 * L, 0.75 * L, L, WL, WR,
                         (lambda xi: perp_sample(xi, WL, WR, w=w), fan_edges(w)),
                         (lambda xi: perp_sample(xi, WR, WL, w=w2), fan_edges(w2)))
    assert np.all(valid)
    near = np.abs(x - 0.25 * L) < np.abs(x - 0.75 * L)
    W1 = perp_sample((x - 0.25 * L) / tt, WL, WR, w=w)
    assert np.allclose(Wp[:, near], W1[:, near])
    print('step_H and the periodic superposition: ok')
