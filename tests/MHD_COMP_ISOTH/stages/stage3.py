#!/usr/bin/env python3
"""Stage 3: 1D isothermal MHD with the prototype proto_mhd1d.py (formulation C) against the exact solutions of
mhd_exact.py and the FV references of fv_mhd1d.py (REPORT_stage3.md).

    python3 stage3.py check                  proto_mhd1d against proto1d (B = 0), smooth waves, semi-discrete identities
    python3 stage3.py fvwin [N] [name ...]   FV runs (N = 4096) of the coplanar single shocks: is the window of the
                                             exact solution valid?  -> runs/stage3/fvwin/<name>_N<N>.json, .npz
    python3 stage3.py run [-f] [-j n] [-eta C,...] [pat ...]  run the cases (n = 6 processes) into
                                             runs/stage3/<case>.json, .npz (cases with a .json are skipped; -f reruns
                                             them); -eta: the cases with C_eta = C (name suffix _eta<C>) in place of
                                             the bulk-viscosity cases (C_eta = 0)
    python3 stage3.py table                  Markdown tables runs/stage3/tables.md
    python3 stage3.py list                   the case names

Units and grid: c = 1 (beta0 = 2), L = 1 (lx = 1/(2 pi)), x_j = j L/N; W = (rho, u_x, u_y, u_z, B_y, B_z).
Each case runs with a fixed dt = t_end/nstep <= SAFETY*DTMIN, DTMIN the minimum of min(dt_adv, dt_visc) over an
adaptive pilot run with the same parameters (cached in runs/stage3/pilot); nstep is a multiple of NWRITE, so the
snapshots fall on t_end j/NWRITE. A failed run (non-finite values, CFL stop, dt lambda_svv > 2.51) is repeated with
dt/2 and dt/4; all attempts are kept. budget() is evaluated at every stage (track_budget = 'stages'): the
semi-discrete identities use the stage-1 state of each step, the energy decomposition (instruction 6) the
stage-weighted rates rk_* of the step, and the time-integration error is dE - sum dt rk_eps.

Families (initial data on the grid, steps with H = 1/2 on a grid point unless stated):
    perp      P1a-c, P2, P3 (B_x = 0), R on (L/4, 3L/4); exact Riemann fans (perp_waves) superposed
    coplanar  F1, F2, S1, S2, I1: downstream outside, upstream on (xa, xb); exact single shock at xa, the fan of
              xb unknown (bounded by char_bounds; checked with fvwin)
    td        TDS (u = 0) and TDA (u_x = 1, one period): c^2 rho + B_y^2/2 = 2 with B_y = sqrt(3) outside, 1 inside,
              tanh width w cells; 'pre': B_y filtered by the 2/3 mask before rho is computed (balanced on the grid)
    rd        rho = B_x = |B_t| = 1, u_x = 0, u_t = -B_t (V = 1), dpsi = pi/2 (RDh) or pi (RDp), width w cells,
              one period; RDpc: dpsi = pi with the jumps between grid points, B_z = u_z = 0 exactly (coplanar)
    compound  C1-C3 of mhd_exact, compared at T_REF with the FV reference (finest N available); C3 also with the
              exact fast-slow shock pair on the colliding side
"""
import sys
import os
import json
import math
import time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mhd_exact as mx
import fv_mhd1d as fv
from proto_mhd1d import MHD1D
from gen_stage1 import SAFETY

OUT = os.path.normpath(os.path.join(HERE, '..', '..', '..', 'runs', 'stage3'))
LX, L = 1.0 / (2.0 * math.pi), 1.0
NS = (256, 512, 1024, 2048)
NWRITE = 20
PERP_TEND = {'P1a': 0.11, 'P1b': 0.09, 'P1c': 0.06, 'P2': 0.09, 'P3': 0.055}
COP_XA = 0.125
COP_TEND = {'F1': 0.16, 'F2': 0.105, 'S1': 0.25, 'S2': 0.25, 'I1': 0.24}
TD_P0 = 2.0
TD_W = ((0.5, 0.0, 0.0, 0.0, math.sqrt(3.0), 0.0), (1.5, 0.0, 0.0, 0.0, 1.0, 0.0))      # outside, inside
RD_DPSI = {'RDh': 0.5 * math.pi, 'RDp': math.pi}


# ------------------------------------------------------------ cases
def case_name(c):
    s = c['prob']
    if c['fam'] == 'td':
        s += '_w%d%s' % (c['w'], 'pre' if c['pre'] else '')
    elif c['fam'] == 'rd':
        s += 'c' if c['cop'] else '_w%d' % c['w']
    s += '_N%d' % c['N']
    if c['C_eta'] > 0.0:
        s += '_eta%g' % c['C_eta']
    return s


def cases(etas=(0.0,)):
    cs = []
    for C_eta in etas:
        for N in NS:
            cs += [dict(fam='perp', prob=p, N=N) for p in mx.PERP]
            cs += [dict(fam='coplanar', prob=p, N=N) for p in mx.COPLANAR]
            cs += [dict(fam='td', prob=p, N=N, w=w, pre=pre) for p in ('TDS', 'TDA') for w in (0, 1, 4)
                   for pre in (False, True)]
            cs += [dict(fam='rd', prob=p, N=N, w=w, cop=False) for p in RD_DPSI for w in (0, 4, 8)]
            cs.append(dict(fam='rd', prob='RDp', N=N, w=0, cop=True))
        for N in NS + (4096,):
            cs += [dict(fam='compound', prob=p, N=N) for p in mx.COMPOUND]
        for c in cs:
            c.setdefault('C_eta', C_eta)
    for c in cs:
        c['name'] = case_name(c)
    return cs


def t_end(c):
    f = c['fam']
    if f == 'perp':
        return PERP_TEND[c['prob']]
    if f == 'coplanar':
        return COP_TEND[c['prob']]
    if f == 'compound':
        return fv.T_REF[-1]
    return 1.0


def bx_of(c):
    if c['fam'] == 'compound':
        return mx.COMPOUND[c['prob']][2]
    return {'perp': 0.0, 'td': 0.0, 'rd': 1.0, 'coplanar': mx.COPLANAR_BX}[c['fam']]


def grid(N):
    return np.arange(N) * (L / N)


def make(c):
    return MHD1D(c['N'], lx=LX, bx=bx_of(c), C_eta=c['C_eta'])


def mix(WA, WB, h):
    return np.array([(1.0 - h) * a + h * b for a, b in zip(WA, WB)])


def cop_interval(prob):
    """(xa, xb) of a coplanar case: at t_end the valid upstream (xa + S t, xb + a t) and downstream
    (xb + b t - L, xa + S t) intervals are equally long, l = xb - xa = (L + (2 S - a - b) t_end)/2, rounded to a
    multiple of 1/256 (grid points for N >= 256)."""
    S, Wd, Wu, kind = mx.coplanar_case(prob)
    a, b = mx.char_bounds(Wu, Wd, mx.COPLANAR_BX)
    l = 0.5 * (L + (2.0 * S - a - b) * COP_TEND[prob])
    return COP_XA, COP_XA + round(l * 256) / 256.0


def cop_window(prob):
    """Last time at which the exact solution of the coplanar case is valid somewhere: the shock from xa meets the
    left edge of the fan of xb, or the right edge of that fan comes round the period and reaches the shock."""
    S, Wd, Wu, kind = mx.coplanar_case(prob)
    xa, xb = cop_interval(prob)
    a, b = mx.char_bounds(Wu, Wd, mx.COPLANAR_BX)
    t1 = (xb - xa) / (S - a)
    t2 = (L - (xb - xa)) / (b - S) if b > S else math.inf
    return min(t1, t2)


def rd_geometry(c):
    dl = L / c['N']
    if c['cop']:
        return 0.25 * L + 0.5 * dl, 0.75 * L + 0.5 * dl, 0.0
    return 0.25 * L, 0.75 * L, c['w'] * dl


def rd_W(c, x, t):
    xa, xb, w = rd_geometry(c)
    W = mx.rd_profile(x, t, 1.0, 0.0, 1.0, 1.0, 0.0, RD_DPSI[c['prob']], xa, xb, w, L)
    if c['cop']:
        W[3] = 0.0
        W[5] = 0.0
    return W


_PROTO = {}


def _proto(N):
    if N not in _PROTO:
        _PROTO[N] = MHD1D(N, lx=LX)
    return _PROTO[N]


def td_U(c):
    return 1.0 if c['prob'] == 'TDA' else 0.0


def td_W(c, x, t):
    N = len(x)
    U = td_U(c)
    (_, _, _, _, bo, _), (_, _, _, _, bi, _) = TD_W
    if c['pre']:
        p = _proto(N)
        by = bo + (bi - bo) * mx.step_H(x, 0.25 * L, 0.75 * L, c['w'] * L / N, L)
        by = p.bwd(p.fwd(by) * p.Fb * np.exp(-1j * p.k * U * t))
    else:
        by = bo + (bi - bo) * mx.step_H(x - U * t, 0.25 * L, 0.75 * L, c['w'] * L / N, L)
    z = np.zeros(N)
    return np.array([mx.td_density(by, 0.0, TD_P0), U + z, z, z, by, z])


def c3_pair():
    WL, WR, bx = mx.COMPOUND['C3']
    return mx.fast_slow_pair(WL, WR, bx)


def initial(c):
    """W (6, N) of the case on the grid (also the initial-data file of the Fortran runs)."""
    x = grid(c['N'])
    f = c['fam']
    if f in ('perp', 'compound'):
        WL, WR = mx.PERP[c['prob']] if f == 'perp' else mx.COMPOUND[c['prob']][:2]
        return mix(WL, WR, mx.step_H(x, 0.25 * L, 0.75 * L, 0.0, L))
    if f == 'coplanar':
        S, Wd, Wu, kind = mx.coplanar_case(c['prob'])
        xa, xb = cop_interval(c['prob'])
        return mix(Wd, Wu, mx.step_H(x, xa, xb, 0.0, L))
    if f == 'rd':
        return rd_W(c, x, 0.0)
    return td_W(c, x, 0.0)


def exact(c, x, t):
    """(W, valid) at t > 0; W = None where no exact solution is known (C1, C2)."""
    f = c['fam']
    allv = np.ones(len(x), dtype=bool)
    if f == 'rd':
        return rd_W(c, x, t), allv
    if f == 'td':
        return td_W(c, x, t), allv
    if f == 'perp':
        WL, WR = mx.PERP[c['prob']]
        w1, w2 = mx.perp_waves(WL, WR), mx.perp_waves(WR, WL)
        return mx.periodic(x, t, 0.25 * L, 0.75 * L, L, WL, WR,
                           (lambda xi: mx.perp_sample(xi, WL, WR, w=w1), mx.fan_edges(w1)),
                           (lambda xi: mx.perp_sample(xi, WR, WL, w=w2), mx.fan_edges(w2)))
    if f == 'coplanar':
        # Wu on (xa + S t, xb), Wd elsewhere (the downstream side of the shock is the outer state); unknown in the
        # fan of xb; invalid everywhere once the fan meets the shock on either side (mx.periodic would stop at the
        # left edge xi = 0 of the sample region (0, S) instead)
        S, Wd, Wu, kind = mx.coplanar_case(c['prob'])
        xa, xb = cop_interval(c['prob'])
        a, b = mx.char_bounds(Wu, Wd, mx.COPLANAR_BX)
        xs = xa + S * t
        up = np.mod(x - xs, L) < np.mod(xb - xs, L)
        W = np.where(up[None, :], np.asarray(Wu)[:, None], np.asarray(Wd)[:, None])
        valid = ~((np.mod(x - xb, L) <= b * t) | (np.mod(xb - x, L) <= -a * t))     # the fan may exceed L/2
        if not (xs < xb + a * t and xb + b * t < xs + L):
            valid[:] = False
        return W, valid
    if c['prob'] == 'C3':
        WL, WR, bx = mx.COMPOUND['C3']
        pr = c3_pair()
        return mx.periodic(x, t, 0.25 * L, 0.75 * L, L, WL, WR,
                           (mx.chain_sample(pr['S'], pr['W']), (pr['S'][0], pr['S'][-1])),
                           (None, mx.char_bounds(WR, WL, bx)))
    return None, ~allv


def features(c):
    """(fronts, obstacles): the exact discontinuities dict(kind, x0, S, Wm, Wp) (Wm on the left) with a density
    jump, and the other moving points (x0, speed) of the exact solution (fan edges)."""
    f, F, ob = c['fam'], [], []
    if f == 'perp':
        WL, WR = mx.PERP[c['prob']]
        for x0, A, B in ((0.25 * L, WL, WR), (0.75 * L, WR, WL)):
            w = mx.perp_waves(A, B)
            for side, Wm, Wp in (('L', A, w['WL']), ('R', w['WR'], B)):
                if w[side]['type'] == 'shock':
                    F.append(dict(kind='fast', x0=x0, S=w[side]['speed'], Wm=Wm, Wp=Wp))
                else:
                    ob += [(x0, w[side]['head']), (x0, w[side]['tail'])]
                if side == 'L':
                    F.append(dict(kind='td', x0=x0, S=w['u_s'], Wm=w['WL'], Wp=w['WR']))
    elif f == 'coplanar':
        S, Wd, Wu, kind = mx.coplanar_case(c['prob'])
        xa, xb = cop_interval(c['prob'])
        F.append(dict(kind=kind, x0=xa, S=S, Wm=Wd, Wp=Wu))
        ob += [(xb, v) for v in mx.char_bounds(Wu, Wd, mx.COPLANAR_BX)]
    elif f == 'compound' and c['prob'] == 'C3':
        WL, WR, bx = mx.COMPOUND['C3']
        pr = c3_pair()
        for i, kind in enumerate(('fast', 'slow', 'slow', 'fast')):
            F.append(dict(kind=kind, x0=0.25 * L, S=pr['S'][i], Wm=pr['W'][i], Wp=pr['W'][i + 1]))
        ob += [(0.75 * L, v) for v in mx.char_bounds(WR, WL, bx)]
    for e in F:
        e['Wm'], e['Wp'], e['S'] = [float(v) for v in e['Wm']], [float(v) for v in e['Wp']], float(e['S'])
    return [e for e in F if abs(e['Wm'][0] - e['Wp'][0]) > 1e-9], ob


# ------------------------------------------------------------ measurements
def pdist(x, x0):
    """Signed periodic distance x - x0 in [-L/2, L/2)."""
    return np.mod(np.asarray(x) - x0 + 0.5 * L, L) - 0.5 * L


def crossing(x, v, xe, vm, half):
    """Offset from xe (|offset| <= half) of the crossing of v through vm closest to xe (linear interpolation); nan if
    there is none."""
    d = pdist(x, xe)
    idx = np.where(np.abs(d) <= half)[0]
    idx = idx[np.argsort(d[idx])]
    if len(idx) < 2:
        return math.nan
    g, dd = v[idx] - vm, d[idx]
    cand = list(dd[g == 0.0])
    k = np.where(g[:-1] * g[1:] < 0.0)[0]
    cand += list(dd[k] + (dd[k + 1] - dd[k]) * g[k] / (g[k] - g[k + 1]))
    return float(min(cand, key=abs)) if cand else math.nan


def erode(m, n):
    m0 = m.copy()
    for k in range(1, n + 1):
        m &= np.roll(m0, k) & np.roll(m0, -k)
    return m


def window(F, ob, k, t, dl):
    """Half-width of the search window of front k at t: 0.45 times the distance to the nearest other moving point,
    at most 0.05 and at least 3 cells."""
    xk = F[k]['x0'] + F[k]['S'] * t
    others = [e['x0'] + e['S'] * t for i, e in enumerate(F) if i != k] + [x0 + s * t for x0, s in ob]
    gap = min([abs(float(pdist(xo, xk))) for xo in others] + [L])
    return max(3.0 * dl, min(0.05, 0.45 * gap))


def fft1(W):
    return np.fft.fft(W[4] + 1j * W[5])[1]


def analyze(c, snaps, W0):
    """Measurements of a run from snaps = [(t, W (6, N), zeta, eta, D)] (snaps[0] at t = 0) and the initial data W0."""
    N = c['N']
    x, dl, te, bx = grid(N), L / N, t_end(c), bx_of(c)
    out = dict(zeta_max=max(float(s[2].max()) for s in snaps), eta_max=max(float(s[3].max()) for s in snaps))
    if len(snaps) < 2:
        return out
    tF, WF = snaps[-1][0], snaps[-1][1]
    out['t_final'] = tF
    We, valid = exact(c, x, tF)
    if We is not None:
        out['l1'] = [float(np.sum(np.abs(WF[i] - We[i])[valid]) * dl) for i in range(6)]
        out['valid_frac'] = float(valid.mean())
    F, ob = features(c)
    if F:
        n_er = max(8, int(math.ceil(0.02 / dl)))
        states, plats = [], []
        for e in F:
            for key in ('Wm', 'Wp'):
                Ws = np.array(e[key])
                i = next((j for j, s in enumerate(states) if np.max(np.abs(s - Ws)) < 1e-12), None)
                if i is None:
                    states.append(Ws)
                    i = len(states) - 1
                    m = valid & np.all(np.abs(We - Ws[:, None]) <= 1e-12 * (1.0 + np.abs(Ws[:, None])), axis=0)
                    m = erode(m, n_er)
                    p = dict(W=list(Ws), npts=int(m.sum()))
                    if m.sum() > 0:
                        p['mean_err'] = list(WF[:, m].mean(axis=1) - Ws)
                        p['max_err'] = list(np.max(np.abs(WF[:, m] - Ws[:, None]), axis=1))
                        p['mean'] = list(WF[:, m].mean(axis=1))
                    plats.append(p)
                e['i' + key[1]] = i
        out['plateaus'] = plats
        fr = []
        for k, e in enumerate(F):
            vm = 0.5 * (e['Wm'][0] + e['Wp'][0])
            ts, ds = [], []
            for t, W, *_ in snaps[1:]:
                if t >= 0.25 * te - 1e-12:
                    ts.append(t)
                    ds.append(crossing(x, W[0], e['x0'] + e['S'] * t, vm, window(F, ob, k, t, dl)))
            r = dict(kind=e['kind'], S=e['S'], off_final_dl=ds[-1] / dl if ds else None, im=e['im'], ip=e['ip'])
            # plateau error on the two sides relative to the largest jump of the front
            pe = [plats[i] for i in (e['im'], e['ip'])]
            jump = float(np.max(np.abs(np.array(e['Wp']) - np.array(e['Wm']))))
            r['plat_npts'] = min(p['npts'] for p in pe)
            if r['plat_npts'] > 0:
                r['plat_err'] = max(max(abs(v) for v in p['mean_err']) for p in pe) / jump
            ok = np.isfinite(ds)
            if ok.sum() >= 2:
                slope = float(np.polyfit(np.array(ts)[ok], np.array(ds)[ok], 1)[0])
                r.update(S_fit=e['S'] + slope, dS=slope, dS_rel=slope / abs(e['S']) if abs(e['S']) > 1e-9 else None)
            # overshoots of rho and B_y next to the front at the final time
            if ds and math.isfinite(ds[-1]):
                xm = e['x0'] + e['S'] * tF + ds[-1]
                band = np.abs(pdist(x, xm)) <= min(max(20 * dl, 0.03), window(F, ob, k, tF, dl))
                for i, nm in ((0, 'rho'), (4, 'By')):
                    lo, hi = sorted((e['Wm'][i], e['Wp'][i]))
                    if hi - lo > 1e-9:
                        r['over_' + nm] = float(max(0.0, WF[i][band].max() - hi) / (hi - lo))
                        r['under_' + nm] = float(max(0.0, lo - WF[i][band].min()) / (hi - lo))
            # Rankine-Hugoniot residual of the measured plateaus on both sides with the fitted (and the exact) speed
            pm, pp = plats[e['im']], plats[e['ip']]
            if 'mean' in pm and 'mean' in pp:
                Wm, Wp = np.array(e['Wm']), np.array(e['Wp'])
                sc = np.abs(mx.flux(Wp, bx) - mx.flux(Wm, bx)) + abs(e['S']) * np.abs(mx.cons(Wp) - mx.cons(Wm))
                big = sc > 1e-9 * sc.max()
                for tag, S in (('', r.get('S_fit')), ('_Sexact', e['S'])):
                    if S is not None:
                        res = mx.rh_residual(np.array(pm['mean']), np.array(pp['mean']), S, bx)
                        r['rh' + tag] = [float(abs(a) / s) if b else None for a, s, b in zip(res, sc, big)]
            fr.append(r)
        out['fronts'] = fr
    if c['fam'] in ('rd', 'td'):
        V = 1.0 if c['fam'] == 'rd' else td_U(c)
        c10 = fft1(W0)
        sh = [float(np.mod(-np.angle(fft1(W) / c10) / (2.0 * math.pi) * L - V * t + 0.5 * L, L) - 0.5 * L)
              for t, W, *_ in snaps]
        out.update(shift_final_dl=sh[-1] / dl, shift_max_dl=max(abs(s) for s in sh) / dl,
                   speed_err=sh[-1] / tF)
        rest = snaps[1:]
        if c['fam'] == 'rd':
            # alfv: the Alfvenic relation u_t = -B_t; z: u_z, B_z of the coplanar RDpc (zero in the data)
            dev = lambda W, D, z: dict(rho=np.max(np.abs(W[0] - 1.0)), ux=np.max(np.abs(W[1])),
                                       bt=np.max(np.abs(np.hypot(W[4], W[5]) - 1.0)), D=np.max(np.abs(D)),
                                       zeta=np.max(z), alfv=np.max(np.abs(W[2:4] + W[4:6])),
                                       **(dict(z=np.max(np.abs(W[[3, 5]]))) if c['cop'] else {}))
        else:
            dev = lambda W, D, z: dict(ux=np.max(np.abs(W[1] - V)), ut=np.max(np.abs(W[2:4])),
                                       ptot=np.max(np.abs(W[0] + 0.5 * (W[4] ** 2 + W[5] ** 2) - TD_P0)),
                                       D=np.max(np.abs(D)), zeta=np.max(z))
        d0 = dev(snaps[0][1], snaps[0][4], snaps[0][2])
        out['dev_initial'] = {k: float(v) for k, v in d0.items()}
        out['dev_max'] = {k: float(max(dev(W, D, z)[k] for t, W, z, eta, D in rest)) for k in d0}
        if c['fam'] == 'td':
            def width(W, t):
                xe = 0.25 * L + V * t
                lo, hi = TD_W[0][0], TD_W[1][0]
                a = crossing(x, W[0], xe, lo + 0.1 * (hi - lo), 0.1)
                b = crossing(x, W[0], xe, lo + 0.9 * (hi - lo), 0.1)
                return (b - a) / dl
            out['width_dl'] = dict(data=width(W0, 0.0), initial=width(snaps[0][1], 0.0), final=width(WF, tF))
    if c['fam'] == 'compound':
        ref = fv_reference(c['prob'])
        if ref is not None:
            Nf, tr, Wr, Wo = ref
            out['fv_N'] = Nf
            l1 = []
            for k, tk in enumerate(tr):
                j = [s for s in snaps if abs(s[0] - tk) < 1e-9]
                if not j:
                    continue
                Wf = Wr[k][:, ::Nf // N]
                e = dict(t=float(tk), l1=[float(np.sum(np.abs(j[0][1][i] - Wf[i])) * dl) for i in range(6)])
                if Wo is not None:
                    e['fv_self'] = float(sum(np.sum(np.abs(Wf[i] - Wo[k][i, ::Wo.shape[2] // N])) for i in range(6)) * dl)
                l1.append(e)
            out['l1_fv'] = l1
    return out


def fv_reference(prob):
    """(N, t, W (len(t), 6, N), W at N/2 or None) of the finest FV reference of a compound problem."""
    for n in (65536, 32768):
        fn = os.path.join(fv.REF_DIR, '%s_N%d.npz' % (prob, n))
        if os.path.exists(fn):
            d = np.load(fn)
            fo = os.path.join(fv.REF_DIR, '%s_N%d.npz' % (prob, n // 2))
            return n, d['t'], d['W'], np.load(fo)['W'] if os.path.exists(fo) else None
    return None


class Stats:
    """Running maxima, minima and time integrals over the hist entries."""
    PARTS = ('ideal', 'visc', 'svv', 'svv_b', 'eta')

    def __init__(self):
        self.d, self.h0, self.E_prev = {}, None, None

    def mx(self, k, v):
        self.d[k] = max(self.d.get(k, -np.inf), float(v))

    def mn(self, k, v):
        self.d[k] = min(self.d.get(k, np.inf), float(v))

    def add(self, k, v):
        self.d[k] = self.d.get(k, 0.0) + float(v)

    def update(self, hist):
        for h in hist:
            if self.h0 is None:
                self.h0 = h
                self.d.update(M0=h['M'], E0=h['E'], n_E_up=0, n_svv_pos=0, n_svv_b_pos=0)
            h0, dt = self.h0, h['dt']
            self.mx('drift_M', abs(h['M'] - h0['M']) / h0['M'])
            for k in ('Px', 'Py', 'Pz', 'By', 'Bz'):
                self.mx('drift_' + k, abs(h[k] - h0[k]))
            self.mx('max_E_kin', h['E_kin'])
            self.mn('rho_min', h['rho_min'])
            self.mx('rho_max', h['rho_max'])
            self.mx('zeta_max', h['zeta_max'])
            self.mx('eta_max', h['eta_max'])
            self.mx('max_dt_over_limit', dt / min(h['dt_adv'], h['dt_visc']))
            self.mx('max_dt_lam_svv', dt * h['lam_svv'])
            if self.E_prev is not None and h['E'] > self.E_prev:
                self.d['n_E_up'] += 1
                self.mx('max_E_up', h['E'] - self.E_prev)
            self.E_prev = h['E']
            for tag in ('', '_raw'):
                for k, A in (('dM', 'A_M'), ('dPx', 'A_P'), ('dPy', 'A_P'), ('dPz', 'A_P'), ('dBy', 'A_B'),
                             ('dBz', 'A_B')):
                    if h[A] > 0.0:
                        self.mx('rate_' + k + tag, abs(h[k + tag]) / h[A])
                if h['visc'] > 0.0:
                    self.mx('R_zeta' + tag, abs(h['R_zeta' + tag]) / h['visc'])
                if h['D_svv'] != 0.0:
                    self.mx('R_svv' + tag, abs(h['R_svv' + tag]) / abs(h['D_svv']))
                if h['eta_diss'] > 0.0:
                    self.mx('R_eta' + tag, abs(h['R_eta' + tag]) / h['eta_diss'])
            if h['D_svv'] != 0.0:
                self.mn('min_D_svv', h['D_svv'])
                self.d['n_svv_pos'] += h['svv_eps'] > 0.0
            self.mx('max_svv_b_eps', h['svv_b_eps'])
            self.d['n_svv_b_pos'] += h['svv_b_eps'] > 0.0
            # energy decomposition with the stage-weighted rates of the step
            for k in self.PARTS:
                self.add('int_rk_%s_eps' % k, dt * h['rk_%s_eps' % k])
            for k in ('eps', 'visc', 'D_svv', 'eta_diss'):
                self.add('int_rk_' + k, dt * h['rk_' + k])

    def result(self):
        return {k: (None if isinstance(v, float) and not math.isfinite(v) else v) for k, v in self.d.items()}


# ------------------------------------------------------------ runs
def pilot_dt(c):
    """DTMIN of the case from an adaptive run (dt = min(dt_adv, dt_visc)), cached in runs/stage3/pilot."""
    fn = os.path.join(OUT, 'pilot', c['name'] + '.json')
    if os.path.exists(fn):
        with open(fn) as f:
            return json.load(f)
    p = make(c)
    p.set_W(initial(c))
    t0 = time.time()
    with np.errstate(all='ignore'):
        st = p.run_until(t_end(c), np.inf)
    lim = [min(h['dt_adv'], h['dt_visc']) for h in p.hist]
    d = dict(status=st, t=p.t, steps=p.step, dtmin=float(min(lim)), wall=time.time() - t0)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with open(fn, 'w') as f:
        json.dump(d, f, indent=1)
    return d


def nstep_of(c, dtmin):
    te = t_end(c)
    return NWRITE * int(math.ceil(te / (SAFETY * dtmin) / NWRITE * (1.0 - 1e-12)))


def snap(p):
    f = p.fields(p.sh, p.mh, p.bh)
    return p.t, np.array([f['rho'], *f['u'], *f['B']]), f['zeta'], f['eta'], f['D']


def run_case(c):
    t_wall = time.time()
    te = t_end(c)
    pilot = pilot_dt(c)
    n0 = nstep_of(c, pilot['dtmin'])
    W0 = initial(c)
    attempts = []
    for fac in (1, 2, 4):
        p = make(c)
        p.set_W(W0)
        p.track_budget = 'stages'
        nstep = n0 * fac
        dt = te / nstep
        stats, snaps, status, t0 = Stats(), [snap(p)], 'ok', time.time()
        with np.errstate(all='ignore'):
            for j in range(1, NWRITE + 1):
                status = p.run(dt, j * nstep // NWRITE - p.step)
                stats.update(p.hist)
                p.hist = []
                if status != 'ok':
                    break
                snaps.append(snap(p))
            res = stats.result()
            if p.finite():
                E = p.monitor(p.rhs(p.sh, p.mh, p.bh)[3])['E']
                res.update(E_end=E, dE=E - res['E0'], time_err=E - res['E0'] - res['int_rk_eps'])
        attempts.append(dict(dt_factor=fac, dt=dt, nstep=nstep, status=status, t_stop=p.t, steps=p.step,
                             wall=time.time() - t0, stats=res))
        if status == 'ok':
            break
    with np.errstate(all='ignore'):
        an = analyze(c, snaps, W0)
    out = dict(case=c, t_end=te, bx=bx_of(c), pilot=pilot, attempts=attempts, ok=status == 'ok', analysis=an)
    np.savez_compressed(os.path.join(OUT, c['name'] + '.npz'), t=np.array([s[0] for s in snaps]),
                        W=np.array([s[1] for s in snaps]), zeta=np.array([s[2] for s in snaps]),
                        eta=np.array([s[3] for s in snaps]), W0=W0)
    with open(os.path.join(OUT, c['name'] + '.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)
    return c['name'], [(a['dt_factor'], a['status'], round(a['t_stop'], 4)) for a in attempts], time.time() - t_wall


def cost(c):
    """Rough cost N^2 t_end vmax (the steps scale as N t_end vmax, a step as N)."""
    W = initial(c)
    rho, B = W[0], W[4:6]
    bx = bx_of(c)
    va2 = (bx * bx + B[0] ** 2 + B[1] ** 2) / rho
    cf = np.sqrt(0.5 * (1 + va2 + np.sqrt((1 + va2) ** 2 - 4 * bx * bx / rho)))
    return c['N'] ** 2 * t_end(c) * np.max(np.abs(W[1]) + cf)


def run(args):
    import multiprocessing as mp
    force = '-f' in args
    nproc = int(args[args.index('-j') + 1]) if '-j' in args else 6
    etas = tuple(float(v) for v in args[args.index('-eta') + 1].split(',')) if '-eta' in args else (0.0,)
    pats = [a for i, a in enumerate(args) if a not in ('-f', '-j', '-eta') and (i == 0 or args[i - 1] not in ('-j', '-eta'))]
    os.makedirs(OUT, exist_ok=True)
    todo = [c for c in cases(etas) if (not pats or any(s in c['name'] for s in pats))
            and (force or not os.path.exists(os.path.join(OUT, c['name'] + '.json')))]
    todo.sort(key=cost, reverse=True)
    print('%d cases' % len(todo), flush=True)
    with mp.Pool(nproc) as pool:
        for name, att, wall in pool.imap_unordered(run_case, todo):
            print('%-24s %s  %.0fs' % (name, att, wall), flush=True)


def fvwin(args):
    """FV runs of the coplanar cases with the initial data of the prototype: the exact solution must agree with the
    FV solution on the valid points away from the shock (> 0.02); the fan of xb must not have entered them."""
    N = int(args[0]) if args and args[0].isdigit() else 4096
    names = [a for a in args if not a.isdigit()] or list(mx.COPLANAR)
    d = os.path.join(OUT, 'fvwin')
    os.makedirs(d, exist_ok=True)
    x = grid(N)
    for prob in names:
        c = dict(fam='coplanar', prob=prob, N=N, C_eta=0.0)
        te, t0 = t_end(c), time.time()
        W, n = fv.run(initial(c), mx.COPLANAR_BX, L, te)
        We, valid = exact(c, x, te)
        S, Wd, Wu, kind = mx.coplanar_case(prob)
        xa, xb = cop_interval(prob)
        far = valid & (np.abs(pdist(x, xa + S * te)) > 0.02)
        err = [float(np.max(np.abs(W[i][far] - We[i][far]))) if far.any() else math.nan for i in range(6)]
        l1 = [float(np.sum(np.abs(W[i] - We[i])[valid]) / N) for i in range(6)]
        r = dict(prob=prob, N=N, t_end=te, t_window=cop_window(prob), steps=n, valid_frac=float(valid.mean()),
                 max_err_far=err, l1=l1, wall=time.time() - t0)
        np.savez_compressed(os.path.join(d, '%s_N%d.npz' % (prob, N)), x=x, W=W, We=We, valid=valid)
        with open(os.path.join(d, '%s_N%d.json' % (prob, N)), 'w') as f:
            json.dump(r, f, indent=1)
        print('%-3s %-12s N=%d t=%.3f (window %.3f) valid %.2f  max|W - We| away from the shock %s  L1 %s  %.0fs'
              % (prob, kind, N, te, r['t_window'], r['valid_frac'], ' '.join('%.1e' % e for e in err),
                 ' '.join('%.1e' % e for e in l1), r['wall']), flush=True)


# ------------------------------------------------------------ checks
def check():
    from proto1d import Isoth1D
    from gen_stage1 import RIEMANN, DTMIN
    # 1. B = 0, u_t = 0: proto_mhd1d is proto1d
    print('1. B = 0, u_t = 0: MHD1D against proto1d.Isoth1D (fixed dt of Stage 1)')
    for prob, N in (('U2', 512), ('U4', 1024)):
        rL, rR, uL, uR, te = RIEMANN[prob]
        dt = SAFETY * DTMIN[(prob, N)]
        nstep = int(round(te / dt))
        a, b = Isoth1D(N), MHD1D(N)
        a.riemann(rL, rR, uL, uR)
        b.riemann(rL, rR, uL, uR)
        sa, sb = a.run(dt, nstep), b.run(dt, nstep)
        A, B = a.state(), b.state()
        print('   %s N=%d %d steps %s %s: |d rho|/max %.1e  |d u|/max %.1e  u_t %.0e B_t %.0e  lam_svv %.1e  '
              'dt_visc %.1e (max relative)' % (
                  prob, N, nstep, sa, sb, np.max(np.abs(A['rho'] - B['rho'])) / np.max(A['rho']),
                  np.max(np.abs(A['u'] - B['u'][0])) / np.max(np.abs(A['u'])), np.max(np.abs(B['u'][1:])),
                  np.max(np.abs(B['B'])), max(abs(x['lam_svv'] / y['lam_svv'] - 1) for x, y in zip(a.hist, b.hist)),
                  max(abs(x['dt_visc'] / y['dt_visc'] - 1) for x, y in zip(a.hist, b.hist))))
    # 2. smooth waves: n steps per period from dt_limit of the finest grid, the same n on every grid
    mode1 = lambda a: np.fft.rfft(a)[1]
    print('2. linear fast and slow modes (rho0 = B_x = B_y0 = 1, amp 1e-6), one period T, dt = T/n, '
          'n = ceil(T/(f dt_limit(N = 64)))')
    for fam in ('fast', 'slow'):
        q = MHD1D(64, lx=LX, bx=1.0)
        W0, om = mx.linear_mode(q.x, 0.0, 1.0, 1.0, 1.0, fam, 1e-6, 1, L)
        q.set_W(W0)
        T = 2.0 * math.pi / om
        for f in (1.0, 0.5, 0.25):
            nper = int(math.ceil(T / (f * q.dt_limit())))
            for N in (16, 32, 64):
                p = MHD1D(N, lx=LX, bx=1.0)
                W0, _ = mx.linear_mode(p.x, 0.0, 1.0, 1.0, 1.0, fam, 1e-6, 1, L)
                p.set_W(W0)
                st = p.run(T / nper, nper)
                W = p.prim()
                We, _ = mx.linear_mode(p.x, T, 1.0, 1.0, 1.0, fam, 1e-6, 1, L)
                a0, a1 = mode1(W0[0] - 1.0), mode1(W[0] - 1.0)
                print('   %-4s f=%.2f n=%4d N=%2d %s  |a|/|a0| - 1 %+.3e (RK3: %+.3e)  phase %+.2e  max|W - We|/amp %.3e'
                      % (fam, f, nper, N, st, abs(a1 / a0) - 1, -(2.0 * math.pi) ** 4 / (24.0 * nper ** 3),
                         np.angle(a1 / a0), max(np.max(np.abs(W[i] - We[i])) for i in range(6)) / 1e-6))
    print('   circularly polarized Alfven wave (rho0 = B_x = |B_t| = 1, mode 1), t = 1, n = ceil(1/(f dt_limit(N = 64)))')
    q = MHD1D(64, lx=LX, bx=1.0)
    q.set_W(mx.cp_alfven(q.x, 0.0, 1.0, 0.0, 1.0, 1.0, 1, L))
    for f in (1.0, 0.5):
        nper = int(math.ceil(1.0 / (f * q.dt_limit())))
        for N in (16, 64):
            p = MHD1D(N, lx=LX, bx=1.0)
            W0 = mx.cp_alfven(p.x, 0.0, 1.0, 0.0, 1.0, 1.0, 1, L)
            p.set_W(W0)
            st = p.run(1.0 / nper, nper)
            W = p.prim()
            print('   f=%.2f n=%4d N=%2d %s  max|W - W0| %.2e  |rho - 1| %.0e  |u_x| %.0e  ||B_t| - 1| %.1e  max zeta %.0e'
                  % (f, nper, N, st, max(np.max(np.abs(W[i] - W0[i])) for i in range(6)), np.max(np.abs(W[0] - 1)),
                     np.max(np.abs(W[1])), np.max(np.abs(np.hypot(W[4], W[5]) - 1)),
                     max(h['zeta_max'] for h in p.hist)))
    # 3. semi-discrete identities at rough states
    print('3. identities at t = 0.08 of tanh Riemann problems (N = 512), relative to A_M, A_P, A_B or the dissipation')
    rough = {'C2': ((1, 0, 0, 0, 1, 0), (0.125, 0, 0, 0, -1, 0), 0.75),
             'rot': ((1, 0.3, 0.2, -0.1, 1, 0.5), (0.125, -0.2, -0.3, 0.4, -1, -0.7), 0.75)}
    for nm, (WL, WR, bx) in rough.items():
        for form in ('entropy', 'linear'):
            for C_eta in (0.0, 2.0):
                p = MHD1D(512, lx=LX, bx=bx, C_eta=C_eta, svv_form=form)
                p.riemann_mhd(WL, WR, width=1.0)
                st = p.run_until(0.08, 1.0)
                b = p.budget()
                P, f = p.parts(p.sh, p.mh, p.bh)
                fs, fm, fb, _ = p.rhs(p.sh, p.mh, p.bh)
                sm = lambda i: sum(v[i] for v in P.values())
                rel = lambda a, r: np.max(np.abs(a - r)) / np.max(np.abs(r))
                q = lambda k, s: abs(b[k]) / s if s else 0.0
                print('   %-3s %-7s C_eta=%.0f %s %4d steps | rhs - sum(parts) %.0e | dM %.0e dP %.0e %.0e %.0e '
                      'dB %.0e %.0e (raw) | R_zeta %.0e R_svv %.0e R_eta %.0e | eps_svv_b %+.1e eta_eps %+.1e'
                      % (nm, form, C_eta, st, p.step, max(rel(fs, sm(0)), rel(fm, sm(1)), rel(fb, sm(2))),
                         q('dM_raw', b['A_M']), q('dPx_raw', b['A_P']), q('dPy_raw', b['A_P']),
                         q('dPz_raw', b['A_P']), q('dBy_raw', b['A_B']), q('dBz_raw', b['A_B']),
                         q('R_zeta', b['visc']), q('R_svv', abs(b['D_svv'])), q('R_eta', b['eta_diss']),
                         b['svv_b_eps'], b['eta_eps']))
    # 4. the stage-weighted rates: dE = dt rk_eps + O(dt^k) on one step of a rough state
    p = MHD1D(512, lx=LX, bx=0.75)
    p.riemann_mhd(*rough['rot'][:2], width=1.0)
    p.run_until(0.04, 1.0)
    dt0 = p.dt_limit()
    print('4. one step of the rot state at t = 0.04: (E1 - E0 - dt rk_eps)/(dt rk_eps)')
    for f in (1.0, 0.5, 0.25):
        q = MHD1D(512, lx=LX, bx=0.75)
        q.sh, q.mh, q.bh, q.pi_v, q.step, q.t = p.sh.copy(), p.mh.copy(), p.bh.copy(), p.pi_v, p.step, p.t
        q.track_budget = 'stages'
        mon, _ = q.advance(f * dt0)
        E1 = q.monitor(q.rhs(q.sh, q.mh, q.bh)[3])['E']
        print('   dt = %.2f dt_limit: %.2e' % (f, (E1 - mon['E'] - f * dt0 * mon['rk_eps']) / (f * dt0 * mon['rk_eps'])))


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'run':
        run(sys.argv[2:])
    elif cmd == 'check':
        check()
    elif cmd == 'fvwin':
        fvwin(sys.argv[2:])
    elif cmd == 'list':
        for c in cases():
            print(c['name'])
    elif cmd == 'table':
        import stage3_table
        stage3_table.main(OUT)
    else:
        sys.exit(__doc__)


if __name__ == '__main__':
    main()
