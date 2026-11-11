#!/usr/bin/env python3
"""Linear stability of the CP Alfven wave of S4.1 (gen_stage4.alfven_case: rho = 1, B0 = 1, c_s = 1,
dB/B0 = eta) by Hill matrices in the wave frame, where the pump is the steady flow u0 = -B(phi), phi = k0.x.

A perturbation of the Bloch family q is sum_n X_n exp(i (q + n k0).x); the pump couples n to n +- 1 (n +- 2
through B B). Two operators:
    physical  ideal isothermal MHD in (rho, u, b), members |n - n_c| <= nmax (the continuum, truncated)
    galerkin  the semi-discrete code without dissipation (formulation C) linearized in its own variables
              (s, m, B), with the members restricted to the 2/3 mask of m and B:
                  d_t ds = -i K.dm,
                  d_t dm = -i K.F,  F = -(dm B_m + B_m dm) - (dB B_m + B_m dB) + (B_m.dB) I - T_m ds
                                        + [m = 0] cs2 ds I,   T_m = sum_{a+b=m} B_a B_b,
                  d_t dB = -i K x (B_m x (dm + dB));
              the s members outside the mask are passive at linear order (div_N m0 = 0) and are left out.
Both give the same rates without truncation (check). Rates are in units of w0 = |k0| (v_A = 1). Not modelled:
the SVV, and eta_shock, whose term curl(eta J0) with the current J0 of the pump is of first order in the
perturbation but not linear (eta ~ S[max(0, -dD)]); the shock viscosity (q = rho zeta D) is of second order.
The runs of S4.1 (analyze_stage4.py cpdiag) follow the Galerkin rate of their dominant family without eta_shock
and exceed it with eta_shock (3D k0 = (1,1,2), default: 0.30 against 0.06 at N = 16, 0.13 against 0.095 at
N = 32); the excess falls with N.

    python3 cp_stability.py check
    python3 cp_stability.py <k0, e.g. 102 or 112> N [N ...]   (2D x-z if k0_y = 0)
"""
import sys

import numpy as np

CS2 = 1.0
NMAX = 40   # physical: |n - n_c| <= NMAX. The truncation in n is itself unstable if too short at large k_perp
            # (k0 = (1,0,1), q = (-21,0,6): 0.062 for NMAX = 20, 0 for 40-80); physical_max checks the largest
            # rate against 1.5 NMAX


def cross_mat(a):
    return np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]], dtype=complex)


def basis(kvec):
    """n, e1, e2 of gen_stage4.cp_basis."""
    k = np.asarray(kvec, float)
    n = k/np.linalg.norm(k)
    e1 = np.array([0.0, 1.0, 0.0]) if k[1] == 0 else np.cross(n, [0.0, 0.0, 1.0])
    e1 /= np.linalg.norm(e1)
    return n, e1, np.cross(n, e1)


def perp_basis(K):
    K = np.asarray(K, float)
    a = np.array([1.0, 0, 0]) if abs(K[0]) < 0.9*np.linalg.norm(K) else np.array([0, 1.0, 0])
    p1 = np.cross(K, a)
    p1 /= np.linalg.norm(p1)
    p2 = np.cross(K/np.linalg.norm(K), p1)
    return np.stack([p1, p2], axis=1).astype(complex)  # 3 x 2


def pump(k0, eta):
    n_hat, e1, e2 = basis(k0)
    c = 0.5*eta*(e1 - 1j*e2)
    return {0: n_hat.astype(complex), 1: c, -1: np.conj(c)}


def physical(k0, eta, q, members, cs2=CS2):
    """Hill matrix of ideal isothermal MHD; unknowns per member rho, u (3), b in the plane normal to K (2)."""
    k0 = np.asarray(k0, float)
    Bh = pump(k0, eta)
    uh = {m: -v for m, v in Bh.items()}
    Jh = {m: np.cross(1j*m*k0, v) for m, v in Bh.items()}
    Fh = {}
    for a in Jh:
        for b in Bh:
            Fh[a + b] = Fh.get(a + b, 0) + np.cross(Jh[a], Bh[b])
    idx = {n: i for i, n in enumerate(members)}
    K = {n: np.asarray(q, float) + n*k0 for n in members}
    Q = {n: perp_basis(K[n]) for n in members}
    L = np.zeros((6*len(members), 6*len(members)), dtype=complex)
    for n in members:
        o = 6*idx[n]
        iK = 1j*K[n]
        curl = Q[n].conj().T @ cross_mat(iK)
        for m in (-1, 0, 1):
            p = n - m
            if p not in idx:
                continue
            i = 6*idx[p]
            iKp, imk = 1j*K[p], 1j*m*k0
            B, u0, F, Qp = Bh[m], uh[m], Fh.get(m, np.zeros(3)), Q[p]
            L[o, i] += -(u0 @ iKp)
            L[o+1:o+4, i] += -F
            L[o+1:o+4, i+1:i+4] += -(u0 @ iKp)*np.eye(3) - np.outer(u0, imk)
            L[o+1:o+4, i+4:i+6] += (-cross_mat(B) @ cross_mat(iKp) + cross_mat(np.cross(imk, B))) @ Qp
            L[o+4:o+6, i+4:i+6] += curl @ cross_mat(u0) @ Qp
            L[o+4:o+6, i+1:i+4] += -curl @ cross_mat(B)
            if m == 0:
                L[o, i+1:i+4] += -iK
                L[o+1:o+4, i] += -cs2*iK
    return L


def galerkin(k0, eta, q, members, cs2=CS2):
    """Hill matrix of the code (see the module header); unknowns per member s, m (3), B normal to K (2)."""
    k0 = np.asarray(k0, float)
    Bh = pump(k0, eta)
    T = {}
    for a in Bh:
        for b in Bh:
            T[a + b] = T.get(a + b, 0) + np.outer(Bh[a], Bh[b])
    idx = {n: i for i, n in enumerate(members)}
    K = {n: np.asarray(q, float) + n*k0 for n in members}
    Q = {n: perp_basis(K[n]) for n in members}
    L = np.zeros((6*len(members), 6*len(members)), dtype=complex)
    I3 = np.eye(3)
    for n in members:
        o = 6*idx[n]
        Kn = K[n]
        L[o, o+1:o+4] += -1j*Kn
        for m in range(-2, 3):
            p = n - m
            if p not in idx:
                continue
            i = 6*idx[p]
            if m in Bh:
                Bm = Bh[m]
                # (K.F)_i = sum_j K_j F_ji: K.(dm B) = (K.dm) B, K.(B dm) = (K.B) dm
                Adm = -(np.outer(Bm, Kn) + (Kn @ Bm)*I3)
                AdB = (Adm + np.outer(Kn, Bm)) @ Q[p]
                L[o+1:o+4, i+1:i+4] += -1j*Adm
                L[o+1:o+4, i+4:i+6] += -1j*AdB
                M = -Q[n].conj().T @ (1j*cross_mat(Kn) @ cross_mat(Bm))
                L[o+4:o+6, i+1:i+4] += M
                L[o+4:o+6, i+4:i+6] += M @ Q[p]
            L[o+1:o+4, i] += -1j*(-(Kn @ T[m]) + (cs2*Kn if m == 0 else 0))
    return L


def rate(L, k0):
    return float(np.linalg.eigvals(L).real.max())/np.linalg.norm(k0)


def mask_members(k0, q, N, active=(True, True, True), nmax=200):
    """n with |q + n k0|_d < N/3 in the active directions and 0 in the others."""
    k0, q = np.asarray(k0), np.asarray(q)
    out = []
    for n in range(-nmax, nmax + 1):
        K = q + n*k0
        if all(abs(K[d]) < N/3 if active[d] else K[d] == 0 for d in range(3)):
            out.append(n)
    return out


def families(k0, N, active=(True, True, True)):
    """One q per family (q ~ q + k0 ~ -q) with members inside the mask, the family of k0 itself (q = 0) excluded:
    its members n != 0 give its rates (the member K = 0 is constant: mass, momentum and mean B), zero (check)."""
    k0 = np.asarray(k0)
    lim = int(np.ceil(N/3))
    rng = [range(-lim, lim + 1) if a else range(1) for a in active]
    seen, out = set(), []
    for qx in rng[0]:
        for qy in rng[1]:
            for qz in rng[2]:
                q = np.array([qx, qy, qz])
                if not all(abs(q[d]) < N/3 for d in range(3)):
                    continue
                key = tuple(int(x) for x in np.cross(q, k0))
                if not any(key) or key in seen or tuple(-x for x in key) in seen:
                    continue
                seen.add(key)
                out.append(q)
    return out


def physical_rate(k0, eta, q, nmax=NMAX):
    k0 = np.asarray(k0)
    nc = -int(round(np.asarray(q) @ k0/(k0 @ k0)))
    return rate(physical(k0, eta, q, list(range(nc - nmax, nc + nmax + 1))), k0)


def galerkin_rate(k0, eta, q, N, active=(True, True, True)):
    return rate(galerkin(k0, eta, q, mask_members(k0, q, N, active)), k0)


def k_perp(k0, q):
    k0 = np.asarray(k0, float)
    return float(np.linalg.norm(np.cross(q, k0))/np.linalg.norm(k0))


def physical_max(k0, eta, N, active=(True, True, True)):
    """Largest physical rate over the families of the grid: (rate, q, rate of q with 1.5 NMAX members)."""
    g, q = max(((physical_rate(k0, eta, q), tuple(int(x) for x in q)) for q in families(k0, N, active)),
               key=lambda t: t[0])
    return g, q, physical_rate(k0, eta, q, int(1.5*NMAX))


def galerkin_max(k0, eta, N, active=(True, True, True)):
    return max(((galerkin_rate(k0, eta, q, N, active), tuple(int(x) for x in q)) for q in families(k0, N, active)),
               key=lambda t: t[0])


if __name__ == '__main__':
    if sys.argv[1] == 'check':
        k0 = np.array([1, 0, 1])
        print('eta = 0: max Re/w0 %.1e' % rate(galerkin(k0, 0.0, (-5, 0, 4), list(range(-20, 21))), k0))
        for q in ((-5, 0, 4), (2, 0, -1)):
            mem = list(range(-40, 41))
            print('no truncation q = %s: galerkin %.4f physical %.4f'
                  % (q, rate(galerkin(k0, 1.0, q, mem), k0), rate(physical(k0, 1.0, q, mem), k0)))
        for k0 in ((1, 0, 1), (1, 0, 2), (1, 1, 2)):
            mem = [n for n in range(-NMAX, NMAX + 1) if n != 0]
            print('family of the wave k0 = %s: max Re/w0 %.1e' % (k0, rate(galerkin(k0, 1.0, (0, 0, 0), mem), k0)))
        for k0, q in (((1, 0, 2), (1, 0, -1)), ((1, 1, 2), (-5, -3, -1))):
            print('physical k0 = %s q = %s: NMAX 20 %.4f 30 %.4f'
                  % (k0, q, physical_rate(k0, 1.0, q, 20), physical_rate(k0, 1.0, q, 30)))
        sys.exit(0)
    k0 = np.array([int(c) for c in sys.argv[1]])
    active = (True, k0[1] != 0, True)
    for N in [int(a) for a in sys.argv[2:]]:
        g, qg = galerkin_max(k0, 1.0, N, active)
        p, qp, p15 = physical_max(k0, 1.0, N, active)
        print('k0 = %s N = %3d  galerkin max %.3f (q = %s, k_perp %.2f, physical %.3f)  physical max %.3f '
              '(%.3f with 1.5 NMAX; q = %s, k_perp %.2f)' % (tuple(k0), N, g, qg, k_perp(k0, qg),
                                                              physical_rate(k0, 1.0, qg), p, p15, qp, k_perp(k0, qp)))
