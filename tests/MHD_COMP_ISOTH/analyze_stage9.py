#!/usr/bin/env python3
"""Stage 9: the mass budget of the holes of rho in turb1 (M_s = 1) and the split of the SVV term of the mass
equation, from the out3d records and conservation.dat of the existing runs (no job, no change of the code).

    python3 analyze_stage9.py selftest          synthetic checks of the tools below
    python3 analyze_stage9.py run ROOT [label]  budgets of the records          -> ROOT/stage9/budget_<label>.json
    python3 analyze_stage9.py closure ROOT      gate 1: d ln rho_min/dt          -> ROOT/stage9/closure.json
    python3 analyze_stage9.py fieldmin ROOT     diagnostic of gate 1: min over the field -> ROOT/stage9/fieldmin.json
    python3 analyze_stage9.py judge ROOT        gate 2, criteria, figures -> ROOT/stage9/summary.json, fig/*.png
    python3 analyze_stage9.py cancel ROOT       diagnostic of the s and u parts -> ROOT/stage9/cancel.json

ROOT is the runs directory, label a key of CASES (default: all). judge stops unless the cases of the verdict
(JUDGED) passed closure.

With filter_lnrho = F and mass_correction = F (checked on the runs) d rho/dt = -div m + T, and T = d rho/dt of
svv_es_core splits exactly into (svv_split; G_d = sig_d m_L, G_H,d = sig_d m_H, kap_d = e_d/(c^2 + max|u|^2))
    T_Ls = -sum_d G_d^*(e_d rho G_d s)
    T_Lu = -sum_d G_d^*((e_d rho/c^2)[u.G_d u - G_d(u^2/2)])
    T_Hs = -sum_d G_H,d^*(kap_d rho c^2 G_H,d s)
    T_Hu = +sum_d G_H,d^*(kap_d rho G_H,d(u^2/2))
With Q = 1 and no filter, T_Lu = 0 and T_Ls = e lap rho, which can only fill a set {rho < theta}.

Budget of a hole: D(theta, t) = int (theta - rho)_+ dV vanishes continuously on the boundary of the set, so
    dD/dt = -S(theta, t),  S = int_{rho < theta} d rho/dt dV = S_adv + S_svv
and S_svv < 0 means that T raises D at that time (digs). Over a window with the trapezoid weights w_t, for
theta = 0.2, 0.1 and the core {rho < 2 rho_min(t)} (theta moves: S keeps its meaning at each time, D does not
have the identity over the window):
    A_T = sum w_t S_T / sum w_t M,  M = int_{rho < theta} rho dV,  Y_T = tau0 A_T,  r = A_svv/|A_adv|

Criteria (tasks/todo.md, Stage 9, fixed before any budget was computed):
    gate 1  at each record, g_c = d ln rho_min/dt of conservation.dat (central difference on unequal steps) and
            g_o = (-div m + T)/rho at the min, a = (|div m| + |T|)/rho; records with |forward - backward| > 0.2 a
            (the min may have jumped) or without a following row are left out; pass if |g_c - g_o| <= 0.05 a at
            all the others and they are at least half of the records
    gate 2  E_D = |D(end) - D(start) + sum w_t (S_adv + S_svv)|/|sum w_t S_adv| <= 0.5 (fixed theta), otherwise
            the label of the region is 'undecided'
    label   dig: r <= -0.5 and |Y_svv| >= 0.1; fill: A_svv >= 0; mid otherwise; undecided: gate 2, or empty
    verdict (first that holds; 128^3 control and 256^3)
            a  both N: theta = 0.1 and core dig
            b  both N: all three regions fill
            c  both N: core digs, and theta = 0.2 or 0.1 fills
            d  the core digs at one N and fills at the other
            e  otherwise
    mechanism  u: A_u < 0 and (A_s >= 0 or |A_u| >= 2|A_s|); s: A_s < 0 and (A_u >= 0 or |A_s| >= 2|A_u|); both
Sensitivity (reported, not used): RTH, YTH and the windows of WIN_ALT.

Gate 1 was fixed for all six cases; closure failed on one row of e0.5_C4, so before any label was seen it was
narrowed to the two cases of the verdict, with the threshold unchanged (REPORT_stage9.md). fieldmin is the
diagnostic of that row and not a criterion: the min of a neighbouring row is the min over all cells of
s(x) + h ds/dt(x), so it shows a change of the cell of the min; its a is taken at the min of the record only.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_ic  # noqa: E402
import svv3d  # noqa: E402
import analyze_stage7 as a7  # noqa: E402
import analyze_stage8 as a8  # noqa: E402

TAU0 = np.pi  # pi/M_s (case.json)
LABEL256 = 'N256'
CASES = dict(a8.CASES)
CASES[LABEL256] = a7.CASES[256]
JUDGED = {128: a8.CONTROL, 256: LABEL256}
T128 = a8.T_OUT  # 0.5, ..., 1.2 tau0
T256 = a7.HOLE_T[256] + (1.0202, 1.1558)
WIN = {128: a8.T_WIN, 256: a7.HOLE_T[256][1:]}  # t_req of the records of the window; 0.5 tau0 only for gate 1
WIN_ALT = {128: a8.T_WIN[:5], 256: T256[1:]}
TERMS = ('adv', 'Ls', 'Lu', 'Hs', 'Hu')
# levels of the budget: the regions of the criteria on a log grid, for the figure
THETA = tuple(sorted(set(np.round(np.geomspace(0.02, 1.0, 25), 4).tolist()) | {0.1, 0.2}))
REGIONS = ('0.2', '0.1', 'core')
GATE1_KINK, GATE1_TOL, GATE2_TOL = 0.2, 0.05, 0.5
RTH, YTH = -0.5, 0.1
RTH_ALT, YTH_ALT = (-0.25, -1.0), (0.05, 0.2)
R_DOUBLED = -0.083  # Stage 8: R of svv_eps0 = 1 (REPORT_stage8.md)
COLOR = {'adv': '#2a78d6', 'svv': '#eb6834', '0.2': '#1baf7a', '0.1': '#eda100', 'core': '#0b0b0b'}


def write_json(root, name, obj):
    d = os.path.join(root, 'stage9')
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), 'w') as f:
        json.dump(obj, f, indent=1, default=lambda x: x.tolist() if hasattr(x, 'tolist') else str(x))


def read_json(root, name):
    p = os.path.join(root, 'stage9', name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


# ---------------------------------------------------------------------------------------- the budget

def svv_split(op, s, u, rho, eps):
    """The four parts (module docstring) of d rho/dt of svv3d.SVV.tend with the weight rho, on the real grid."""
    cs2 = op.cs2
    q = 0.5*sum(x**2 for x in u)
    umax2 = 2*float(q.max())
    sk, qk = op.fft(s), op.fft(q)
    uk = [op.fft(x) for x in u]
    out = {k: np.zeros_like(sk) for k in TERMS[1:]}
    for d in range(3):
        e = eps[d]
        if e == 0:
            continue
        com = u[0]*op._grad(uk[0], d, True) + u[1]*op._grad(uk[1], d, True) + u[2]*op._grad(uk[2], d, True) \
            - op._grad(qk, d, True)
        out['Ls'] -= op._adj(e*rho*op._grad(sk, d, True), d, True)
        out['Lu'] -= op._adj(e*rho/cs2*com, d, True)
        if op.filter_mom:
            kap = e/(cs2 + umax2)
            out['Hs'] -= op._adj(kap*rho*cs2*op._grad(sk, d, False), d, False)
            out['Hu'] += op._adj(kap*rho*op._grad(qk, d, False), d, False)
    return {k: op.ifft(v) for k, v in out.items()}


def sublevel(rho, terms, levels, dV):
    """For each level: n = cells with rho < level, M, D and S_T of the set (cumulative sums over the sorted rho)."""
    r = rho.ravel()
    idx = np.argsort(r, kind='stable')
    r = r[idx]
    levels = np.asarray(levels, dtype=float)
    n = np.searchsorted(r, levels, side='left')
    M = np.concatenate(([0.0], np.cumsum(r)))[n]*dV
    S = {k: np.concatenate(([0.0], np.cumsum(T.ravel()[idx])))[n]*dV for k, T in terms.items()}
    return dict(n=n.tolist(), M=M.tolist(), D=(levels*n*dV - M).tolist(), S={k: v.tolist() for k, v in S.items()})


def record_budget(g, op, sn):
    """The budget of one record (Snap); the terms as real arrays are not kept."""
    dV = g.Lx*g.Ly*g.Lz/g.ntot
    eg = op.eps(sn.rho, sn.u, sn.b)
    T = svv_split(op, sn.s, sn.u, sn.rho, eg)
    T['adv'] = -op.ifft(1j*(g.kxd*sn.mk[0] + g.kyd*sn.mk[1] + g.kzd*sn.mk[2])*op.mL)
    rho_min = float(sn.rho.min())
    res = sublevel(sn.rho, T, THETA + (2*rho_min,), dV)
    pts = a7.pick_holes(sn.rho, a7.local_min(sn.rho))
    holes = [dict(p=list(p), rho=float(sn.rho[p]), **{k: float(T[k][p]) for k in TERMS}) for p in pts]
    return dict(rho_min=rho_min, rho_max=float(sn.rho.max()), eps=eg, dV=dV, holes=holes,
                sum_abs={k: float(np.abs(T[k]).sum()*dV) for k in TERMS},
                sum_all={k: float(T[k].sum()*dV) for k in TERMS}, **res)


def run_case(root, label):
    """budget_<label>.json of the records of T128 or T256; the comparison with CHECK at the records it has."""
    case = os.path.join(root, CASES[label])
    g, op = svv3d.from_case(case)
    N = g.nlx
    cons = a7.cons_table(case)
    print('%s: N %d svv_eps0 %g conservation.dat rows %d' % (label, N, op.eps0, len(cons)), flush=True)
    res = dict(label=label, N=N, case=case, eps0=op.eps0, cs2=op.cs2, theta=list(THETA), snaps=[])
    for tt in (T256 if N == 256 else T128):
        r = a7.pick(case, tt)
        if r is None:
            print('%s t/tau0 %.4f: no record' % (label, tt))
            continue
        sn = a7.Snap(g, op, r[1], r[2])
        rb = record_budget(g, op, sn)
        snap = dict(t=r[0], t_tau=r[0]/TAU0, t_req=tt, path=r[1], rec=r[2],
                    cons=a7.cons_check(cons, r[0], rb['rho_min']), **rb)
        h = rb['holes'][0]
        snap['check'] = a7.check_ref(a7.CHECK.get((N, tt)), dict(t=r[0], holes=[dict(
            p=tuple(h['p']), budget=dict(adv=h['adv'], svv=sum(h[k] for k in TERMS[1:])))]))
        res['snaps'].append(snap)
        print('%s t/tau0 %.4f record done (rho_min of the row: rel %+.1e, check %s)'
              % (label, r[0]/TAU0, snap['cons']['rel'], None if snap['check'] is None else snap['check']['PASS']),
              flush=True)
    write_json(root, 'budget_%s.json' % label, res)
    return res


def run(root, labels=None):
    for label in labels or CASES:
        run_case(root, label)


# ------------------------------------------------------------- the cancellation of the s and u parts

def pair_stats(ss, uu, su):
    """From |a|^2, |b|^2 and a.b: |a + b|/|a|, |b|/|a| and the correlation of a and b."""
    return dict(net_over_s=float(np.sqrt(max(ss + uu + 2*su, 0.0)/ss)), u_over_s=float(np.sqrt(uu/ss)),
                corr=float(su/np.sqrt(ss*uu)) if uu > 0 else None)


def cancel_flux(op, s, u, rho, theta):
    """The flux of each band without e_d, rho*(a + b) with L: a = c^2 G_d s, b = u.G_d u - G_d(u^2/2) and
    H: a = c^2 G_H,d s, b = -G_H,d(u^2/2), over the field and over {rho < theta}; summed over d."""
    q = 0.5*sum(x**2 for x in u)
    sk, qk = op.fft(s), op.fft(q)
    uk = [op.fft(x) for x in u]
    regions = dict(all=rho, hole=rho*(rho < theta))
    acc = {b + '_' + r: np.zeros(3) for b in 'LH' for r in regions}

    def add(band, a, b):
        for r, w in regions.items():
            acc[band + '_' + r] += [float(((w*a)**2).sum()), float(((w*b)**2).sum()), float((w*w*a*b).sum())]

    for d in range(3):
        add('L', op.cs2*op._grad(sk, d, True),
            sum(u[c]*op._grad(uk[c], d, True) for c in range(3)) - op._grad(qk, d, True))
        add('H', op.cs2*op._grad(sk, d, False), -op._grad(qk, d, False))
    return {k: pair_stats(*v) for k, v in acc.items()}


def cancel_shells(g, op, T, adv, edges=(0.0, 0.5, 0.75, 1.0, 1.51)):
    """|T_k|/|adv_k| and their correlation in shells of max_d |k_d|/kc_d, kc_d = the 2/3 limit of direction d."""
    kq = 0.0
    for k, n, L in ((g.kx, g.nlx, g.Lx), (g.ky, g.nly, g.Ly), (g.kz, g.nlz, g.Lz)):
        kq = np.maximum(kq, np.abs(k)/(2*np.pi*(n//2)/L*2/3))
    w = np.where(g.ky == 0, 1.0, 2.0)  # the half spectrum of the r2c transform
    Tk, Ak = op.fft(T), op.fft(adv)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (kq > lo) & (kq <= hi)
        pt, pa = float((w*np.abs(Tk)**2)[m].sum()), float((w*np.abs(Ak)**2)[m].sum())
        c = float((w*(Tk*np.conj(Ak)).real)[m].sum())
        out.append(dict(lo=lo, hi=hi, T_over_adv=float(np.sqrt(pt/pa)) if pa > 0 else None,
                        corr=c/np.sqrt(pt*pa) if pa > 0 and pt > 0 else None))
    return out


def cancel(root, theta=0.2):
    """Diagnostic, not a criterion: does the cancellation of the s and u parts already hold for the fluxes?"""
    res = {}
    for N, tts in ((128, (T128[3], T128[5])), (256, (a7.HOLE_T[256][3], a7.HOLE_T[256][4]))):
        case = os.path.join(root, CASES[JUDGED[N]])
        g, op = svv3d.from_case(case)
        res[N] = []
        for tt in tts:
            r = a7.pick(case, tt)
            sn = a7.Snap(g, op, r[1], r[2])
            T = svv_split(op, sn.s, sn.u, sn.rho, op.eps(sn.rho, sn.u, sn.b))
            Ts, Tu = T['Ls'] + T['Hs'], T['Lu'] + T['Hu']
            adv = -op.ifft(1j*(g.kxd*sn.mk[0] + g.kyd*sn.mk[1] + g.kzd*sn.mk[2])*op.mL)
            hole = sn.rho < theta
            tend = {r_: pair_stats(float((Ts[m]**2).sum()), float((Tu[m]**2).sum()), float((Ts[m]*Tu[m]).sum()))
                    for r_, m in (('all', np.ones_like(hole)), ('hole', hole))}
            rec = dict(t_tau=r[0]/TAU0, theta=theta, n_hole=int(hole.sum()),
                       flux=cancel_flux(op, sn.s, sn.u, sn.rho, theta), tend=tend,
                       shells=cancel_shells(g, op, Ts + Tu, adv))
            res[N].append(rec)
            print('N %d t/tau0 %.4f n(rho < %g) %d' % (N, rec['t_tau'], theta, rec['n_hole']))
            for k, v in rec['flux'].items():
                print('  flux %-7s net/s %.3f u/s %.3f corr %+.3f' % (k, v['net_over_s'], v['u_over_s'], v['corr']))
            for k, v in tend.items():
                print('  T    %-7s net/s %.3f u/s %.3f corr %+.3f' % (k, v['net_over_s'], v['u_over_s'], v['corr']))
            for v in rec['shells']:
                print('  shell %.2f-%.2f T/adv %s corr %s' % (v['lo'], v['hi'], v['T_over_adv'], v['corr']), flush=True)
    write_json(root, 'cancel.json', res)
    return res


# ---------------------------------------------------------------------------------------------- gates

def closure_row(t, lr, t0, g_o, a):
    """Gate 1 at time t0 from the times t (ascending) and ln rho_min lr of conservation.dat."""
    i = int(np.argmin(np.abs(t - t0)))
    out = dict(dt_row=float(t[i] - t0), g_o=float(g_o), a=float(a))
    if i == 0 or i == len(t) - 1:
        return dict(out, status='end')
    hb, hf = t[i] - t[i-1], t[i+1] - t[i]
    gb, gf = (lr[i] - lr[i-1])/hb, (lr[i+1] - lr[i])/hf
    gc = (hb*gf + hf*gb)/(hb + hf)
    out.update(g_c=float(gc), g_b=float(gb), g_f=float(gf), res=float(abs(gc - g_o)/a))
    if abs(gf - gb) > GATE1_KINK*a:
        return dict(out, status='kink')
    return dict(out, status='ok' if abs(gc - g_o) <= GATE1_TOL*a else 'FAIL')


def gate1(rows):
    """pass of gate 1 from the statuses of the records of a case."""
    st = [r['status'] for r in rows]
    n_ok = st.count('ok')
    return dict(n=len(st), n_ok=n_ok, n_fail=st.count('FAIL'), n_kink=st.count('kink'), n_end=st.count('end'),
                PASS=st.count('FAIL') == 0 and 2*n_ok >= len(st))


def cons_rows(case):
    """(t, ln rho_min) of the chain of conservation.dat, by time; the repeated row of a continued run once. The
    number of rows and of steps that differ from the dt column tell whether the chain has gaps."""
    c = a7.cons_table(case)
    c = c[np.argsort(c[:, 0], kind='stable')]
    c = c[np.append(np.diff(c[:, 0]) > 0, True)]
    return c[:, 0], np.log(c[:, 11]), int(np.sum(np.abs(np.diff(c[:, 0]) - c[:-1, 1]) > 1e-9)), len(c), \
        float(np.abs(c[:, 16]).max())


def closure(root):
    res = {}
    for label in CASES:
        b = read_json(root, 'budget_%s.json' % label)
        if b is None:
            print('%s: no budget' % label)
            continue
        t, lr, ngap, nrow, corr = cons_rows(os.path.join(root, CASES[label]))
        rows = []
        print('%s: rows %d, steps that differ from the dt column %d, max|dmass_corr| %.1e'
              % (label, nrow, ngap, corr))
        print('   t/tau0   rho_min      g_c        g_o     |g_c-g_o|/a      a     status   g_b       g_f')
        for sn in b['snaps']:
            h = sn['holes'][0]
            T = sum(h[k] for k in TERMS[1:])
            r = closure_row(t, lr, sn['t'], (h['adv'] + T)/h['rho'], (abs(h['adv']) + abs(T))/h['rho'])
            r.update(t_tau=sn['t_tau'], rho_min=sn['rho_min'], rel_rho_min=sn['cons']['rel'])
            rows.append(r)
            print('  %7.4f %9.4e %s %9.3f %10s %9.3f %-6s %s %s'
                  % (sn['t_tau'], sn['rho_min'], a8.fmt(r.get('g_c'), '%9.3f'), r['g_o'],
                     a8.fmt(r.get('res'), '%.4f'), r['a'], r['status'], a8.fmt(r.get('g_b'), '%9.3f'),
                     a8.fmt(r.get('g_f'), '%9.3f')))
        g = gate1(rows)
        print('   gate 1: %s (%d checked, %d kink, %d without a following row, %d over %.2f)'
              % ('PASS' if g['PASS'] else 'FAIL', g['n_ok'] + g['n_fail'], g['n_kink'], g['n_end'], g['n_fail'],
                 GATE1_TOL))
        res[label] = dict(rows=rows, nrow=nrow, ngap=ngap, max_dmass_corr=corr, **g)
    res['PASS'] = all(res[k]['PASS'] for k in CASES if k in res) and all(k in res for k in CASES)
    write_json(root, 'closure.json', res)
    print('gate 1, all cases:', 'PASS' if res['PASS'] else 'FAIL')
    return res


def field_min(s, rate, h):
    """The min over the field of s + h rate (first order in h) and its cell."""
    f = s + h*rate
    i = int(np.argmin(f))
    return float(f.flat[i]), tuple(int(x) for x in np.unravel_index(i, f.shape))


def fieldmin_row(t, lr, t0, s, rate, a):
    """The differences of ln rho_min of the rows around t0 and the same from field_min; d_b, d_f: cells between
    the min at t0 and the predicted min of the row."""
    i = int(np.argmin(np.abs(t - t0)))
    s0, p0 = field_min(s, rate, 0.0)
    out = dict(a=float(a), lr_row=float(lr[i]), s0=s0)
    hb = t[i] - t[i-1]
    fb, pb = field_min(s, rate, -hb)
    gb, gb_f = (lr[i] - lr[i-1])/hb, (s0 - fb)/hb
    out.update(g_b=float(gb), g_b_f=float(gb_f), d_b=a7.pdist(pb, p0, s.shape[0]))
    if i == len(t) - 1:
        res = abs(gb - gb_f)/a
    else:
        hf = t[i+1] - t[i]
        ff, pf = field_min(s, rate, hf)
        gf, gf_f = (lr[i+1] - lr[i])/hf, (ff - s0)/hf
        gc, gc_f = (hb*gf + hf*gb)/(hb + hf), (hb*gf_f + hf*gb_f)/(hb + hf)
        res = abs(gc - gc_f)/a
        out.update(g_f=float(gf), g_f_f=float(gf_f), d_f=a7.pdist(pf, p0, s.shape[0]), g_c=float(gc), g_c_f=float(gc_f))
    return dict(out, res=float(res), status='ok' if res <= GATE1_TOL else 'FAIL')


def kink_rate(t, lr, t0, t1):
    """Rows of [t0, t1] where the step-wise d ln rho_min/dt jumps by more than 5 median; the min of two branches
    can only lower the slope (concave)."""
    g = np.diff(lr)/np.diff(t)
    dg = np.diff(g)
    m = (t[1:-1] >= t0) & (t[1:-1] <= t1)
    thr = 5*float(np.median(np.abs(dg[m])))
    return dict(rows=int(m.sum()), thr=thr, n=int(np.sum(np.abs(dg[m]) > thr)), n_concave=int(np.sum(dg[m] < -thr)))


def fieldmin(root):
    res = {}
    for label in CASES:
        b = read_json(root, 'budget_%s.json' % label)
        case = os.path.join(root, CASES[label])
        g, op = svv3d.from_case(case)
        t, lr, _, nrow, _ = cons_rows(case)
        rows = []
        print('%s: rows %d' % (label, nrow))
        print('   t/tau0     a      g_b    g_b(f)  d_b     g_f    g_f(f)  d_f     g_c    g_c(f)  res/a  status', flush=True)
        for sb in b['snaps']:
            sn = a7.Snap(g, op, sb['path'], sb['rec'])
            T = svv_split(op, sn.s, sn.u, sn.rho, op.eps(sn.rho, sn.u, sn.b))
            div = op.ifft(1j*(g.kxd*sn.mk[0] + g.kyd*sn.mk[1] + g.kzd*sn.mk[2])*op.mL)
            rate = (sum(T.values()) - div)/sn.rho
            p = tuple(sb['holes'][0]['p'])
            r = fieldmin_row(t, lr, sb['t'], sn.s, rate, (abs(div[p]) + abs(sum(v[p] for v in T.values())))/sn.rho[p])
            r.update(t_tau=sb['t_tau'])
            rows.append(r)
            print('  %7.4f %8.3f %8.3f %8.3f %4.1f %s %s %s %s %s %.4f %s'
                  % (sb['t_tau'], r['a'], r['g_b'], r['g_b_f'], r['d_b'], a8.fmt(r.get('g_f'), '%8.3f'),
                     a8.fmt(r.get('g_f_f'), '%8.3f'), a8.fmt(r.get('d_f'), '%4.1f'), a8.fmt(r.get('g_c'), '%8.3f'),
                     a8.fmt(r.get('g_c_f'), '%8.3f'), r['res'], r['status']), flush=True)
        n_fail = sum(r['status'] == 'FAIL' for r in rows)
        kr = kink_rate(t, lr, 0.5*TAU0, 1.2*TAU0)
        print('   fieldmin: %s (%d records, %d over %.2f); slope jumps in 0.5-1.2 tau0: %s'
              % ('PASS' if n_fail == 0 else 'FAIL', len(rows), n_fail, GATE1_TOL, kr), flush=True)
        res[label] = dict(rows=rows, n=len(rows), n_fail=n_fail, kink=kr, PASS=n_fail == 0)
    res['PASS'] = all(res[k]['PASS'] for k in CASES)
    write_json(root, 'fieldmin.json', res)
    print('fieldmin, all cases:', 'PASS' if res['PASS'] else 'FAIL')
    return res


# ---------------------------------------------------------------------------------------------- criteria

def trap_weights(t):
    t = np.asarray(t, dtype=float)
    w = np.zeros_like(t)
    w[:-1] += 0.5*np.diff(t)
    w[1:] += 0.5*np.diff(t)
    return w


def region_level(b, reg):
    """Index of the level of a region in the lists of a record: the core is the last one."""
    return len(b['theta']) if reg == 'core' else int(np.argmin(np.abs(np.array(b['theta']) - float(reg))))


def window(t, S, M, D=None):
    """A_T, Y_T, r and E_D of one region over records at t; S[k], M, D are per record (D: fixed theta only)."""
    w = trap_weights(t)
    Mb = float(np.sum(w*np.asarray(M)))
    I = {k: float(np.sum(w*np.asarray(v))) for k, v in S.items()}
    I['svv'] = sum(I[k] for k in TERMS[1:])
    I.update(s=I['Ls'] + I['Hs'], u=I['Lu'] + I['Hu'], L=I['Ls'] + I['Lu'], H=I['Hs'] + I['Hu'])
    out = dict(I=I, Mbar=Mb, empty=Mb == 0.0)
    if Mb == 0.0:
        return out
    out['A'] = {k: v/Mb for k, v in I.items()}
    out['Y'] = {k: TAU0*v/Mb for k, v in I.items()}
    out['r'] = I['svv']/abs(I['adv']) if I['adv'] != 0.0 else float('nan')
    if D is not None:
        out['dD'] = float(D[-1] - D[0])
        out['E_D'] = abs(out['dD'] + I['adv'] + I['svv'])/abs(I['adv']) if I['adv'] != 0.0 else float('inf')
    return out


def label(w, rth=RTH, yth=YTH):
    """The label of a region from window()."""
    if w['empty'] or w.get('E_D', 0.0) > GATE2_TOL:
        return 'undecided'
    if w['A']['svv'] >= 0:
        return 'fill'
    return 'dig' if w['r'] <= rth and abs(w['Y']['svv']) >= yth else 'mid'


def mechanism(As, Au):
    if Au < 0 and (As >= 0 or abs(Au) >= 2*abs(As)):
        return 'u'
    if As < 0 and (Au >= 0 or abs(As) >= 2*abs(Au)):
        return 's'
    return 'both'


def verdict(lab):
    """The verdict from lab[N][region], N = 128, 256."""
    Ns = (128, 256)
    if all(lab[N]['0.1'] == 'dig' and lab[N]['core'] == 'dig' for N in Ns):
        return 'a'
    if all(lab[N][r] == 'fill' for N in Ns for r in REGIONS):
        return 'b'
    if all(lab[N]['core'] == 'dig' and 'fill' in (lab[N]['0.2'], lab[N]['0.1']) for N in Ns):
        return 'c'
    if {lab[128]['core'], lab[256]['core']} == {'dig', 'fill'}:
        return 'd'
    return 'e'


def regions_of(b, times):
    """{region: window()} of budget b over the records with t_req in times."""
    snaps = [s for s in b['snaps'] if s['t_req'] in times]
    t = [s['t'] for s in snaps]
    out = {}
    for reg in REGIONS:
        i = region_level(b, reg)
        S = {k: [s['S'][k][i] for s in snaps] for k in TERMS}
        out[reg] = window(t, S, [s['M'][i] for s in snaps], None if reg == 'core' else [s['D'][i] for s in snaps])
        svv = np.sum([S[k] for k in TERMS[1:]], axis=0)
        ns = [s['n'][i] for s in snaps]
        out[reg].update(mech=None if out[reg]['empty'] else mechanism(out[reg]['A']['s'], out[reg]['A']['u']),
                        n_neg=int(np.sum(svv < 0)), n_rec=len(snaps), n_empty=int(sum(n == 0 for n in ns)),
                        n_cells=ns, svv_t=svv.tolist(), t_tau=[s['t_tau'] for s in snaps],
                        M_t=[s['M'][i] for s in snaps], adv_t=S['adv'])
    return out


def curve_of(b, times):
    """Y_adv and Y_svv against theta over the window (figure); None where the set is empty over the window."""
    snaps = [s for s in b['snaps'] if s['t_req'] in times]
    t = [s['t'] for s in snaps]
    out = dict(theta=b['theta'], adv=[], svv=[])
    for i in range(len(b['theta'])):
        S = {k: [s['S'][k][i] for s in snaps] for k in TERMS}
        w = window(t, S, [s['M'][i] for s in snaps])
        out['adv'].append(None if w['empty'] else w['Y']['adv'])
        out['svv'].append(None if w['empty'] else w['Y']['svv'])
    return out


def hole_aux(b, times):
    """Over the holes of the records of the window: the fraction with T < 0 at the point, and of each sign of T^(s),
    T^(u); the median of T/|adv|."""
    hs = [h for s in b['snaps'] if s['t_req'] in times for h in s['holes']]
    T = np.array([sum(h[k] for k in TERMS[1:]) for h in hs])
    Ts = np.array([h['Ls'] + h['Hs'] for h in hs])
    Tu = np.array([h['Lu'] + h['Hu'] for h in hs])
    adv = np.array([h['adv'] for h in hs])
    return dict(n=len(hs), frac_svv_neg=float(np.mean(T < 0)), frac_s_neg=float(np.mean(Ts < 0)),
                frac_u_neg=float(np.mean(Tu < 0)), med_svv_over_adv=float(np.median(T/np.abs(adv))))


def judge(root):
    cl, fm = read_json(root, 'closure.json'), read_json(root, 'fieldmin.json')
    if cl is None or not all(cl[JUDGED[N]]['PASS'] for N in JUDGED):
        sys.exit('judge: gate 1 did not pass for a case of the verdict (or closure was not run); see closure.json')
    B = {k: read_json(root, 'budget_%s.json' % k) for k in CASES}
    main = {N: regions_of(B[JUDGED[N]], WIN[N]) for N in JUDGED}
    lab = {N: {r: label(main[N][r]) for r in REGIONS} for N in JUDGED}
    v = verdict(lab)
    # sensitivity: thresholds, one at a time, then the windows
    sens = {}
    for name, kw in [('r %g' % x, dict(rth=x)) for x in RTH_ALT] + [('Y %g' % x, dict(yth=x)) for x in YTH_ALT]:
        lb = {N: {r: label(main[N][r], **kw) for r in REGIONS} for N in JUDGED}
        sens[name] = dict(verdict=verdict(lb), lab=lb)
    alt = {N: regions_of(B[JUDGED[N]], WIN_ALT[N]) for N in JUDGED}
    lb = {N: {r: label(alt[N][r]) for r in REGIONS} for N in JUDGED}
    sens['window'] = dict(verdict=verdict(lb), lab=lb, r={N: {r: alt[N][r].get('r') for r in REGIONS} for N in JUDGED},
                          E_D={N: {r: alt[N][r].get('E_D') for r in REGIONS} for N in JUDGED})
    # auxiliary: the 128^3 cases, the holes, the first order estimate of doubling svv_eps0
    cases128 = {k: regions_of(B[k], WIN[128]) for k in a8.CASES}
    aux = dict(A_over_e={k: {r: (None if w[r]['empty'] else w[r]['A']['svv']/B[k]['eps0']) for r in REGIONS}
                         for k, w in cases128.items()},
               r_cases={k: {r: w[r].get('r') for r in REGIONS} for k, w in cases128.items()},
               holes={N: hole_aux(B[JUDGED[N]], WIN[N]) for N in JUDGED},
               R_doubled=R_DOUBLED)
    span = (B[JUDGED[128]]['snaps'][-1]['t'] - [s for s in B[JUDGED[128]]['snaps'] if s['t_req'] == WIN[128][0]][0]['t'])
    aux['first_order'] = {r: None if main[128][r]['empty'] else 0.5*span*main[128][r]['A']['svv'] for r in REGIONS}
    summ = dict(criteria=dict(RTH=RTH, YTH=YTH, GATE2_TOL=GATE2_TOL), verdict=v, label=lab, main=main, sens=sens,
                aux=aux, curve={N: curve_of(B[JUDGED[N]], WIN[N]) for N in JUDGED},
                gate1=dict(closure=cl['PASS'], closure_cases={k: cl[k]['PASS'] for k in CASES},
                           fieldmin_cases=None if fm is None else {k: fm[k]['PASS'] for k in CASES}))
    write_json(root, 'summary.json', summ)
    print_judge(summ)
    figures(root, summ)
    return summ


def print_judge(s):
    for N in (128, 256):
        print('N = %d  (Y: ln mass per tau0)' % N)
        print('  region  label      r      Y_adv    Y_svv | Y_Ls    Y_Lu    Y_Hs    Y_Hu  | Y_s     Y_u   mech  E_D  '
              'S_svv<0  empty/rec')
        for reg in REGIONS:
            w = s['main'][N][reg]
            if w['empty']:
                print('  %-6s  %-9s (empty)' % (reg, s['label'][N][reg]))
                continue
            Y = w['Y']
            print('  %-6s  %-9s %+6.3f %+8.3f %+8.3f | %+6.3f %+7.3f %+7.3f %+7.3f | %+6.3f %+6.3f %-4s %s %3d/%-3d %3d/%-3d'
                  % (reg, s['label'][N][reg], w['r'], Y['adv'], Y['svv'], Y['Ls'], Y['Lu'], Y['Hs'], Y['Hu'], Y['s'],
                     Y['u'], w['mech'], a8.fmt(w.get('E_D'), '%.2f'), w['n_neg'], w['n_rec'], w['n_empty'], w['n_rec']))
    print('verdict:', s['verdict'])
    print('sensitivity:', {k: v['verdict'] for k, v in s['sens'].items()})
    print('aux holes:', s['aux']['holes'])
    print('aux A_svv/e:', s['aux']['A_over_e'])
    print('aux r of the 128^3 cases:', s['aux']['r_cases'])
    print('aux first order estimate of doubling svv_eps0 (128^3):', s['aux']['first_order'], 'against R =', R_DOUBLED)


def figures(root, s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = os.path.join(root, 'stage9', 'fig')
    os.makedirs(d, exist_ok=True)
    # Y_adv and Y_svv against theta
    fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.8), sharex=True, sharey=True)
    for ax, N in zip(axs, (128, 256)):
        c = s['curve'][N]
        for k, name in (('adv', 'Y_adv (-div m)'), ('svv', 'Y_svv (SVV)')):
            y = [np.nan if x is None else x for x in c[k]]
            ax.plot(c['theta'], y, color=COLOR[k], lw=2.0, marker='o', ms=4, label=name)
        ax.axhline(0, color='#898781', lw=0.8)
        ax.set_xscale('log')
        ax.set_yscale('symlog', linthresh=1.0)
        ax.set_xlabel('theta')
        ax.set_title('%d^3, window %.3f-%.3f tau0' % (N, WIN[N][0], WIN[N][-1]), fontsize=10)
        ax.grid(alpha=0.35, lw=0.5, color='#e1e0d9')
    axs[0].set_ylabel('tau0 x (rate of the mass of {rho < theta}) / mass')
    axs[0].legend(frameon=False, fontsize=8)
    fig.savefig(os.path.join(d, 'stage9_theta.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    # S_svv/M at the records
    fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.8), sharey=True)
    for ax, N in zip(axs, (128, 256)):
        for reg in REGIONS:
            w = s['main'][N][reg]
            y = [TAU0*a/m if m > 0 else np.nan for a, m in zip(w['svv_t'], w['M_t'])]
            ax.plot(w['t_tau'], y, color=COLOR[reg], lw=2.0, marker='o', ms=4,
                    label='core (rho < 2 rho_min)' if reg == 'core' else 'rho < %s' % reg)
        ax.axhline(0, color='#898781', lw=0.8)
        ax.set_yscale('symlog', linthresh=1.0)
        ax.set_xlabel('t / tau0')
        ax.set_title('%d^3' % N, fontsize=10)
        ax.grid(alpha=0.35, lw=0.5, color='#e1e0d9')
    axs[0].set_ylabel('tau0 S_svv / M')
    axs[0].legend(frameon=False, fontsize=8)
    fig.savefig(os.path.join(d, 'stage9_time.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


# -------------------------------------------------------------------------------------------- tests

def smooth_rough(g, seed=1):
    """A state with content in all bands, for the identities that hold on any state."""
    rng = np.random.default_rng(seed)
    op0 = svv3d.SVV(g, 1.0)

    def rough(amp):
        f = rng.standard_normal(op0.shape)
        fk = op0.fft(f)
        k2 = g.kx**2 + g.ky**2 + g.kz**2
        return amp*op0.ifft(fk*np.exp(-k2/(0.5*float(k2.max()))))
    return rough(0.3), [rough(0.2), rough(0.2), rough(0.2)]


def selftest():
    ok = True

    def rep(name, good, msg):
        nonlocal ok
        ok &= bool(good)
        print('%-44s %s  %s' % (name, msg, 'PASS' if good else 'FAIL'))

    # 1. the four parts sum to tend, on a state with content in all bands; each sums to 0 over the domain
    g = gen_ic.Grid(32, 32, 32)
    op = svv3d.SVV(g, 1.0, 0.5, 0.5, True, True)
    s, u = smooth_rough(g)
    rho = np.exp(s)
    eps = [0.07, 0.11, 0.05]
    T = svv_split(op, s, u, rho, eps)
    rdot, _, _, rL, rH = op.tend(s, u, rho, eps, parts=True)
    scale = float(np.abs(rdot).max())
    rep('sum of the parts = tend', np.abs(sum(T.values()) - rdot).max() <= 1e-12*scale,
        'max diff/max|rdot| %.1e' % (np.abs(sum(T.values()) - rdot).max()/scale))
    rep('Ls + Lu = rL, Hs + Hu = rH', max(np.abs(T['Ls'] + T['Lu'] - rL).max(),
                                          np.abs(T['Hs'] + T['Hu'] - rH).max()) <= 1e-12*scale,
        'H band max %.3e, L band max %.3e' % (np.abs(rH).max(), np.abs(rL).max()))
    rep('H band present in the test', np.abs(rH).max() > 1e-3*np.abs(rL).max() and np.abs(T['Lu']).max() > 0, '')
    sums = max(abs(float(v.sum()))/float(np.abs(v).sum()) for v in T.values())
    rep('each part sums to 0', sums <= 1e-12, 'max %.1e' % sums)
    # cancel_flux: u = 0 leaves the s part alone; s = u^2/(2 c^2) makes the flux of the H band vanish
    cf0 = cancel_flux(op, s, [np.zeros(op.shape)]*3, rho, 2.0)
    rep('cancel_flux, u = 0', abs(cf0['L_all']['net_over_s'] - 1) < 1e-14 and abs(cf0['H_all']['net_over_s'] - 1) < 1e-14,
        '')
    cf1 = cancel_flux(op, 0.5*sum(x**2 for x in u)/op.cs2, u, rho, 2.0)
    rep('cancel_flux, flat w in the H band', cf1['H_all']['net_over_s'] < 1e-10 and cf1['H_all']['corr'] < -1 + 1e-10,
        'net/s %.1e corr %+.6f' % (cf1['H_all']['net_over_s'], cf1['H_all']['corr']))
    # 2. u = 0: T^(u) = 0 and the rest is unchanged in form (kap = e/c^2)
    zero = np.zeros(op.shape)
    T0 = svv_split(op, s, [zero, zero, zero], rho, eps)
    r0, _, _ = op.tend(s, [zero, zero, zero], rho, eps)
    rep('u = 0: T^(u) = 0, T^(s) = tend', np.abs(T0['Lu']).max() == 0 and np.abs(T0['Hu']).max() == 0 and
        np.abs(T0['Ls'] + T0['Hs'] - r0).max() <= 1e-12*np.abs(r0).max(), '')
    # 3. Q = 1 and no filter, smooth state: G_d = d/dx_d, T_Lu = 0, T_Ls = e lap rho, S(theta) >= 0
    N = 64
    g = gen_ic.Grid(N, N, N)
    op = svv3d.SVV(g, 1.0, 0.5, 0.5, True, False)
    for d, (k, shp) in enumerate(((g.kx, (-1, 1, 1)), (g.ky, (1, 1, -1)), (g.kz, (1, -1, 1)))):
        k = k.ravel()
        op.sig[d] = (1j*k).reshape(shp)
    kl = [2*np.pi/g.Lx, 2*np.pi/g.Ly, 2*np.pi/g.Lz]
    s = 0.5*np.cos(kl[0]*g.X)*np.cos(kl[1]*g.Y) + 0.3*np.sin(kl[2]*g.Z)
    u = [0.2*np.sin(kl[1]*g.Y), 0.15*np.cos(kl[2]*g.Z + kl[0]*g.X), 0.1*np.sin(kl[0]*g.X)]
    rho = np.exp(s)
    e = 0.03
    T = svv_split(op, s, u, rho, [e, e, e])
    rk = op.fft(rho)
    lap = op.ifft(-(g.kx**2 + g.ky**2 + g.kz**2)*rk)
    d1 = float(np.abs(T['Ls'] - e*lap).max()/np.abs(e*lap).max())
    d2 = float(np.abs(T['Lu']).max()/np.abs(e*lap).max())
    rep('Q = 1: T_Ls = e lap rho, T_Lu = 0', d1 < 1e-10 and d2 < 1e-10, 'T_Ls %.1e  T_Lu %.1e' % (d1, d2))
    levels = np.linspace(float(rho.min())*1.01, float(rho.max())*0.99, 60)
    S = np.array(sublevel(rho, dict(T=T['Ls']), levels, g.Lx*g.Ly*g.Lz/g.ntot)['S']['T'])
    rel = float(S.min()/np.abs(S).max())
    rep('Q = 1: S(theta) >= 0', rel >= -1e-3, 'min S/max|S| %+.1e' % rel)
    # 4. sublevel against direct sums
    rng = np.random.default_rng(2)
    r = rng.uniform(0.1, 2.0, (16, 16, 16))
    x = rng.standard_normal((16, 16, 16))
    sl = sublevel(r, dict(x=x), [0.5, 1.0], 0.25)
    good = True
    for j, th in enumerate((0.5, 1.0)):
        m = r < th
        good &= sl['n'][j] == int(m.sum()) and abs(sl['M'][j] - float(r[m].sum())*0.25) < 1e-12 \
            and abs(sl['S']['x'][j] - float(x[m].sum())*0.25) < 1e-12 \
            and abs(sl['D'][j] - float((th - r[m]).sum())*0.25) < 1e-12
    rep('sublevel = direct sums', good, '')
    # 5. dD/dt = -S on rho = 1 - 0.7 cos(t) cos(x) sin(y) with theta fixed: trapezoid converges to Delta D
    g = gen_ic.Grid(64, 64, 64)
    dV = g.Lx*g.Ly*g.Lz/g.ntot
    sh = np.cos(2*np.pi*g.X/g.Lx)*np.sin(2*np.pi*g.Y/g.Ly)
    t0, t1, th = 0.3, 1.1, 0.5
    err = []
    for n in (5, 9, 17, 33):
        ts = np.linspace(t0, t1, n)
        S, D = [], []
        for t in ts:
            sl = sublevel(1 - 0.7*np.cos(t)*sh, dict(a=0.7*np.sin(t)*sh), [th], dV)
            S.append(sl['S']['a'][0])
            D.append(sl['D'][0])
        err.append(abs(D[-1] - D[0] + float(np.sum(trap_weights(ts)*np.array(S))))/abs(D[-1] - D[0]))
    rep('dD/dt = -S (trapezoid)', err[-1] < 1e-2 and err[-1] < 0.25*err[0],
        'rel err %s' % np.array2string(np.array(err), precision=1))
    # 6. weights, window, labels, verdict, mechanism by hand
    rep('trapezoid weights', np.allclose(trap_weights([0.0, 1.0, 3.0]), [0.5, 1.5, 1.0]), '')
    S = dict(adv=[2.0, 4.0, 0.0], Ls=[-1.0, 0.0, 0.0], Lu=[0.0, 1.0, 0.0], Hs=[0.0, 0.0, -2.0], Hu=[0.0, 0.0, 0.0])
    w = window([0.0, 1.0, 3.0], S, [1.0, 2.0, 0.0], D=[1.0, 2.0, 3.0])
    # I_adv = 1 + 6 = 7, I_Ls = -0.5, I_Lu = 1.5, I_Hs = -2, I_svv = -1; Mbar = 0.5 + 3 = 3.5; dD = 2
    good = abs(w['I']['adv'] - 7) < 1e-14 and abs(w['I']['svv'] + 1) < 1e-14 and abs(w['Mbar'] - 3.5) < 1e-14 and \
        abs(w['A']['svv'] + 1/3.5) < 1e-14 and abs(w['Y']['adv'] - TAU0*2) < 1e-14 and \
        abs(w['r'] + 1/7) < 1e-14 and abs(w['E_D'] - 8/7) < 1e-14 and \
        abs(w['I']['s'] + 2.5) < 1e-14 and abs(w['I']['u'] - 1.5) < 1e-14
    rep('window by hand', good, '')
    wl = lambda r, Y, ED=0.1: dict(empty=False, r=r, A=dict(svv=np.sign(r)*abs(Y)), Y=dict(svv=Y), E_D=ED)
    good = label(wl(-0.5, -0.1)) == 'dig' and label(wl(-0.4, -0.1)) == 'mid' and label(wl(-0.6, -0.05)) == 'mid' \
        and label(wl(0.2, 0.3)) == 'fill' and label(wl(-0.9, -1.0, 0.6)) == 'undecided' \
        and label(dict(empty=True)) == 'undecided' and label(wl(-0.4, -0.1), rth=-0.25) == 'dig'
    rep('labels', good, '')

    def lb(a, b, c, a2=None, b2=None, c2=None):
        return {128: {'0.2': a, '0.1': b, 'core': c},
                256: {'0.2': a2 or a, '0.1': b2 or b, 'core': c2 or c}}
    good = verdict(lb('fill', 'dig', 'dig')) == 'a' and verdict(lb('fill', 'fill', 'fill')) == 'b' and \
        verdict(lb('fill', 'mid', 'dig')) == 'c' and verdict(lb('mid', 'fill', 'dig')) == 'c' and \
        verdict(lb('mid', 'mid', 'dig', c2='fill')) == 'd' and verdict(lb('mid', 'mid', 'mid')) == 'e' and \
        verdict(lb('fill', 'fill', 'dig', 'mid', 'mid', 'dig')) == 'e' and \
        verdict(lb('undecided', 'dig', 'dig')) == 'a' and verdict(lb('fill', 'undecided', 'fill')) == 'e'
    rep('verdict table', good, '')
    good = mechanism(1.0, -0.1) == 'u' and mechanism(-1.0, -2.0) == 'u' and mechanism(-1.0, -1.5) == 'both' and \
        mechanism(-0.1, 1.0) == 's' and mechanism(-2.0, -1.0) == 's' and mechanism(1.0, 1.0) == 'both'
    rep('mechanism', good, '')
    # 7. gate 1: a quadratic on unequal steps, a kink, the last row
    t = np.array([0.0, 0.1, 0.3, 0.4])
    lr = 1.0 + 2.0*t - 3.0*t**2
    r = closure_row(t, lr, 0.1, 2.0 - 6.0*0.1, 10.0)
    good = r['status'] == 'ok' and abs(r['g_c'] - 1.4) < 1e-12
    good &= closure_row(t, lr, 0.1, 2.0, 10.0)['status'] == 'FAIL'
    good &= closure_row(t, lr, 0.4, 0.0, 10.0)['status'] == 'end'
    lk = np.array([0.0, 0.1, -1.0, -1.1])
    good &= closure_row(t, lk, 0.1, 0.0, 10.0)['status'] == 'kink'
    rep('closure_row', good, '')
    rows = [dict(status=x) for x in ('ok', 'ok', 'kink', 'end')]
    good = gate1(rows)['PASS'] and not gate1(rows + [dict(status='kink')])['PASS'] and \
        not gate1(rows + [dict(status='FAIL')])['PASS']
    rep('gate 1', good, '')
    # 8. the min moves to another cell between the rows: field_min gives both slopes of the kink
    sf, rf = np.ones((8, 8, 8)), np.zeros((8, 8, 8))
    sf[1, 2, 3], rf[1, 2, 3] = 0.0, 1.0
    sf[6, 2, 3], rf[6, 2, 3] = 0.001, -1.0
    tf = np.array([-0.01, 0.0, 0.01])
    r = fieldmin_row(tf, np.array([-0.01, 0.0, -0.009]), 0.0, sf, rf, 2.0)
    good = abs(r['g_b_f'] - 1.0) < 1e-12 and abs(r['g_f_f'] + 0.9) < 1e-12 and r['d_b'] == 0 and r['d_f'] == 3 and \
        r['res'] < 1e-12 and r['status'] == 'ok'
    r = fieldmin_row(tf, np.array([-0.01, 0.0, 0.01]), 0.0, sf, rf, 2.0)
    good &= r['status'] == 'FAIL'
    r = fieldmin_row(tf[:2], np.array([-0.01, 0.0]), 0.0, sf, rf, 2.0)
    good &= r['status'] == 'ok' and 'g_f' not in r
    rep('fieldmin_row', good, '')
    kr = kink_rate(np.arange(0.0, 1.0, 0.01), np.minimum(np.arange(0.0, 1.0, 0.01), 0.5 - 0.5*np.arange(0.0, 1.0, 0.01)) +
                   1e-6*np.sin(40*np.arange(0.0, 1.0, 0.01)), 0.0, 1.0)
    rep('kink_rate', kr['n'] >= 1 and kr['n_concave'] == kr['n'], str(kr))
    print('selftest:', 'all PASS' if ok else 'FAILED')
    return ok


if __name__ == '__main__':
    a = sys.argv[1:]
    if a[:1] == ['selftest']:
        sys.exit(0 if selftest() else 1)
    if len(a) >= 2 and a[0] == 'run':
        run(a[1], a[2:])
    elif len(a) == 2 and a[0] == 'closure':
        closure(a[1])
    elif len(a) == 2 and a[0] == 'fieldmin':
        fieldmin(a[1])
    elif len(a) == 2 and a[0] == 'judge':
        judge(a[1])
    elif len(a) == 2 and a[0] == 'cancel':
        cancel(a[1])
    else:
        sys.exit(__doc__)
