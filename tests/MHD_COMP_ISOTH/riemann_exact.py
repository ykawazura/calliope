#!/usr/bin/env python3
"""Exact Riemann solution of 1D isothermal Euler (sound speed c).

Middle state rho*, u*:
    u* = u_L - f(rho*; rho_L) = u_R + f(rho*; rho_R),
    f(rho; rho_K) = c (rho - rho_K)/sqrt(rho rho_K)   (shock, rho > rho_K)
                  = c ln(rho/rho_K)                   (rarefaction)
Shock speeds:  S_L = u_L - c sqrt(rho*/rho_L),  S_R = u_R + c sqrt(rho*/rho_R)
               (upstream Mach number M1 = sqrt(rho2/rho1), i.e. rho2/rho1 = M1^2).
Rarefaction fans use the Riemann invariants u +- c ln rho.

The periodic test of fields.F90 (init_riemann) has the R state on (L/4, 3L/4):
interface (L | R) at L/4 and its mirror (R | L) at 3L/4. periodic_solution()
superposes the two solutions, valid until waves from the two interfaces meet.
"""
import math
import numpy as np


def f_K(rho, rho_K, c):
    if rho > rho_K:
        return c * (rho - rho_K) / math.sqrt(rho * rho_K)
    return c * math.log(rho / rho_K)


def middle_state(rho_L, u_L, rho_R, u_R, c=1.0):
    """rho*, u* by bisection in ln(rho) (g is monotonically increasing)."""
    g = lambda s: f_K(math.exp(s), rho_L, c) + f_K(math.exp(s), rho_R, c) + (u_R - u_L)
    lo, hi = math.log(min(rho_L, rho_R)) - 50.0, math.log(max(rho_L, rho_R)) + 50.0
    assert g(lo) < 0.0 < g(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    rs = math.exp(0.5 * (lo + hi))
    us = u_L - f_K(rs, rho_L, c)
    return rs, us


def waves(rho_L, u_L, rho_R, u_R, c=1.0):
    """Wave structure: dict with rho*, u*, and for each side ('L', 'R') the type and speeds."""
    rs, us = middle_state(rho_L, u_L, rho_R, u_R, c)
    w = dict(rho_s=rs, u_s=us)
    if rs > rho_L:
        w['L'] = dict(type='shock', speed=u_L - c * math.sqrt(rs / rho_L), mach=math.sqrt(rs / rho_L))
    else:
        w['L'] = dict(type='rarefaction', head=u_L - c, tail=us - c)
    if rs > rho_R:
        w['R'] = dict(type='shock', speed=u_R + c * math.sqrt(rs / rho_R), mach=math.sqrt(rs / rho_R))
    else:
        w['R'] = dict(type='rarefaction', head=u_R + c, tail=us + c)
    return w


def sample(xi, rho_L, u_L, rho_R, u_R, c=1.0, w=None):
    """(rho, u) at similarity coordinate xi = x/t (array)."""
    if w is None:
        w = waves(rho_L, u_L, rho_R, u_R, c)
    xi = np.asarray(xi, dtype=float)
    rho = np.empty_like(xi)
    u = np.empty_like(xi)
    rs, us = w['rho_s'], w['u_s']
    for n, x in np.ndenumerate(xi):
        if x <= us:  # left of the contact (there is none in isothermal flow, u* separates L and R waves)
            if w['L']['type'] == 'shock':
                r, v = (rho_L, u_L) if x < w['L']['speed'] else (rs, us)
            else:
                if x < w['L']['head']:
                    r, v = rho_L, u_L
                elif x > w['L']['tail']:
                    r, v = rs, us
                else:
                    v = x + c
                    r = rho_L * math.exp((u_L - v) / c)
        else:
            if w['R']['type'] == 'shock':
                r, v = (rho_R, u_R) if x > w['R']['speed'] else (rs, us)
            else:
                if x > w['R']['head']:
                    r, v = rho_R, u_R
                elif x < w['R']['tail']:
                    r, v = rs, us
                else:
                    v = x - c
                    r = rho_R * math.exp((v - u_R) / c)
        rho[n], u[n] = r, v
    return rho, u


def max_speed(w):
    s = []
    for side in ('L', 'R'):
        if w[side]['type'] == 'shock':
            s.append(abs(w[side]['speed']))
        else:
            s += [abs(w[side]['head']), abs(w[side]['tail'])]
    return max(s)


def fan_edges(w):
    """Leftmost and rightmost wave speeds of a Riemann fan."""
    a = w['L']['speed'] if w['L']['type'] == 'shock' else w['L']['head']
    b = w['R']['speed'] if w['R']['type'] == 'shock' else w['R']['head']
    return a, b


def periodic_solution(x, t, rho_L, u_L, rho_R, u_R, L=1.0, c=1.0):
    """Superposed solution of the two interfaces of init_riemann at time t > 0.
    A point belongs to interface n if it lies inside that fan, x0_n + [a_n, b_n] t;
    between the fans the state is R on the interval that started as (L/4, 3L/4) and L
    elsewhere. Returns (rho, u, valid); valid is False once the two fans overlap."""
    x = np.mod(np.asarray(x, dtype=float), L)
    ifaces = [(0.25 * L, (rho_L, u_L, rho_R, u_R)), (0.75 * L, (rho_R, u_R, rho_L, u_L))]
    rho = np.where((x > 0.25 * L) & (x < 0.75 * L), rho_R, rho_L).astype(float)
    u = np.where((x > 0.25 * L) & (x < 0.75 * L), u_R, u_L).astype(float)
    edges = []
    for x0, st in ifaces:
        w = waves(*st, c)
        a, b = fan_edges(w)
        xi = (np.mod(x - x0 + 0.5 * L, L) - 0.5 * L) / t
        sel = (xi >= a) & (xi <= b)
        rho[sel], u[sel] = sample(xi[sel], *st, c, w)
        edges.append((x0 + a * t, x0 + b * t))
    (a1, b1), (a2, b2) = edges
    valid = (b1 < a2) and (b2 < a1 + L)
    return rho, u, valid


if __name__ == '__main__':
    # checks against the closed forms of the colliding / receding flows
    for U in (0.5, 1.0, 2.0, 4.0):
        w = waves(1.0, U, 1.0, -U)
        sq = 0.5 * (U + math.sqrt(U * U + 4.0))
        r = sq * sq
        Vs = U / (r - 1.0)
        assert abs(w['rho_s'] / r - 1) < 1e-12 and abs(w['R']['speed'] / Vs - 1) < 1e-12
        assert abs(w['R']['mach'] ** 2 - w['rho_s']) < 1e-12  # rho2/rho1 = M1^2
        wr = waves(1.0, -U, 1.0, U)
        assert abs(wr['rho_s'] / math.exp(-U) - 1) < 1e-12
        # Rankine-Hugoniot in the lab frame (mass and momentum fluxes)
        S, rs, us = w['R']['speed'], w['rho_s'], w['u_s']
        assert abs((rs * (us - S)) - (1.0 * (-U - S))) < 1e-12
        assert abs((rs * us * (us - S) + rs) - (-U * (-U - S) + 1.0)) < 1e-11
        # periodic superposition: the receding pair at 3L/4 reaches beyond L/2 at late times
        xg = np.linspace(0.0, 2.0 * math.pi, 4096, endpoint=False)
        tt = 0.8 * 0.5 * math.pi / (max_speed(w) + max_speed(wr))
        rp, up, ok = periodic_solution(xg, tt, 1.0, U, 1.0, -U, L=2.0 * math.pi)
        assert ok
        r1, u1 = sample((xg - 0.25 * 2.0 * math.pi) / tt, 1.0, U, 1.0, -U)
        r2, u2 = sample((xg - 0.75 * 2.0 * math.pi) / tt, 1.0, -U, 1.0, U)
        near1 = np.abs(xg - 0.5 * math.pi) < np.abs(xg - 1.5 * math.pi)
        assert np.allclose(np.where(near1, r1, r2), rp) and np.allclose(np.where(near1, u1, u2), up)
        print('U = %4.1f: rho2/rho1 = %.6f, V_s = %.6f, M1 = %.6f, rho*(rarefaction) = %.6f'
              % (U, w['rho_s'], w['R']['speed'], w['R']['mach'], wr['rho_s']))
    # Sod-type rho_L/rho_R = 10
    w = waves(10.0, 0.0, 1.0, 0.0)
    print('Sod-type 10:1: rho* = %.6f, u* = %.6f, S_R = %.6f, M1 = %.6f, fan = [%.6f, %.6f]'
          % (w['rho_s'], w['u_s'], w['R']['speed'], w['R']['mach'], w['L']['head'], w['L']['tail']))
