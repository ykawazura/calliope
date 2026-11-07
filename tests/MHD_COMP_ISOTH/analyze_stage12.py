#!/usr/bin/env python3
"""Stage 12: which of the four parts of T_svv (Ls, Lu, Hs, Hu) is not proportional to svv_eps0 at the 64 lowest
cells of rho (Q1), and does the amplitude of their arguments fall with eps0 so that the size of the parts does not
change (Q2, the mechanism of REPORT_stage11.md, section 10)? Cases e0.25_C2, e0.5_C2 (control) and e1.0_C2 of
Stage 8, the 15 records of 0.60-1.00 tau0 of Stage 11.

    python3 analyze_stage12.py selftest      synthetic checks and the errors delta_T from the dense rerun of Stage 10
    python3 analyze_stage12.py terms ROOT    Q1: gates C_T, Sigma and the labels   -> ROOT/stage12/terms.json
    python3 analyze_stage12.py fields ROOT   Q2: amplitudes, norms, D_svv, spectra -> ROOT/stage12/fields_<eps0>.json
    python3 analyze_stage12.py judge ROOT    Q2: classes, conclusions, figures     -> ROOT/stage12/summary.json, fig/*.png

ROOT is the runs directory (stage8/coef, stage10/budget.json, stage11/budget_<eps0>.json and gate.json). terms reads
no out3d; fields reads the 15 out3d records (login node, no job); judge stops unless gates B, R and Sigma2 passed.

The parts of analyze_stage9.svv_split as T_X = -+ sum_d G_d^*(coef_d rho a_X,d) (+ for Hu only), q = u^2/2:
    a_Ls,d = G_d s,    a_Lu,d = (u.G_d u - G_d q)/c^2   (low band,  coef_d = e_d)
    a_Hs,d = G_H,d s,  a_Hu,d = G_H,d q/c^2             (high band, coef_d = kappa_d c^2, kappa_d = e_d/(c^2 + max u^2))
A_X = sqrt(<sum_d a_X,d^2>) and ||T_X||_2 = sqrt(<T_X^2>), < > the mean over the grid.

Criteria (tasks/todo.md, Stage 12, fixed before the values of the parts of e1.0 and e0.25 were computed):
 Q1 (K = 64; Ls, Lu, Hs, Hu)
    gate C_T  |I_T of the 5 records of the control - the same of the dense rerun at the same times| <= 0.25 delta_T
    gate Sig  |Ls + Lu + Hs + Hu - S_svv| <= 1e-10 sum|S_T| at the 15 records (true by construction: S_svv of
              analyze_stage10.series is the sum of the parts; hence gate Sigma2)
    test      delta_T,e = delta_T max(1, r_T,e), r_T,e = tmean|S_T/m_64|_e/tmean|S_T/m_64|_0.5; F, u and the label
              of analyze_stage11.test; undecided if delta_T/|I_T,0.5| >= 0.5
 Q2 (the domain)
    gate B    the parts built from the arguments = analyze_stage9.svv_split (relative L2 difference <= 1e-12)
    gate R    t and record as in Stage 11, |S_T(K = 64) - the same of stage11/budget_<eps0>.json| <= 1e-8 sum|S_T|
    gate Sig2 |sum over the 64 cells of rdot of svv3d.SVV.tend - the sum of the parts of Stage 11| <= 1e-10 sum|S_T|
    p_A,X = ln(A_X,1.0/A_X,0.25)/ln 4 and q_X the same of ||T_X||_2 (window means by the trapezoid; also per record);
    classes of p_A: (-inf,-1.25] out, (-1.25,-0.75] -1, (-0.75,-0.25] -1/2, (-0.25,+0.25) 0, [+0.25,inf) out, of q
    the same shifted by +1; a class only if the mean and the 5 records agree, else undecided; (p_A, q) = (-1, 0):
    mechanism, (-1/2, +1/2): amplitude falls but T not constant, (0, +1): no mechanism, else the pair as it is.
delta_T is the largest |error| of I_T of 5 records 0.1 tau0 apart against the dense integral over the same window,
over the 6 starts 0.60-0.65 tau0 of the dense rerun of Stage 10 (selftest recomputes it).
Reported, not used: the sums s, u, L, H, svv and K = 8, 512 (Q1); f_hi, f_lo; r_T,e; ||T_X||_1, ||T_svv||_1,2, the
amplitudes on the 64 cells, the exponents against the control, D_svv, c_dl, the coefficients, the shell spectra (Q2).
"""
import contextlib
import io
import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_ic  # noqa: E402
import svv3d  # noqa: E402
import analyze_stage7 as a7  # noqa: E402
import analyze_stage9 as a9  # noqa: E402
import analyze_stage10 as a10  # noqa: E402
import analyze_stage11 as a11  # noqa: E402
from analyze_stage1 import read_nml  # noqa: E402

TAU0, RUNS = a9.TAU0, a11.RUNS
CASES, EPS, COLOR = a11.CASES, a11.EPS, a11.COLOR
NREC, KS, K_MAIN, T_SLACK = a11.NREC, a11.KS, a11.K_MAIN, a11.T_SLACK
PARTS = a9.TERMS[1:]  # Ls, Lu, Hs, Hu: the judged parts
SUMS = {'s': ('Ls', 'Hs'), 'u': ('Lu', 'Hu'), 'L': ('Ls', 'Lu'), 'H': ('Hs', 'Hu'), 'svv': PARTS}
ALL = PARTS + tuple(SUMS)
# delta_T, I_T of the 5 records from 0.60 tau0 and I_T of the dense records of the same window (K = 64), as in the plan
DELTA_TXT = {
    8: dict(Ls='43.86', Lu='45.03', Hs='11.38', Hu='11.20', s='50.16', u='51.41', L='1.6115', H='0.4407',
            svv='1.4799'),
    64: dict(Ls='4.8605', Lu='4.7495', Hs='2.2871', Hu='2.2322', s='4.9419', u='4.8965', L='0.1420', H='0.3214',
             svv='0.1872'),
    512: dict(Ls='1.4838', Lu='1.5692', Hs='0.2300', Hu='0.1983', s='1.5758', u='1.6410', L='0.2032', H='0.0476',
              svv='0.1556')}
I5_TXT = dict(Ls='+53.5584', Lu='-54.5174', Hs='+6.8366', Hu='-6.6737', s='+60.3950', u='-61.1911', L='-0.9590',
              H='+0.1629', svv='-0.7961')
DENSE_TXT = dict(Ls='+51.2799', Lu='-52.1788', Hs='+8.1475', Hu='-8.1646', s='+59.4273', u='-60.3434', L='-0.8989',
                 H='-0.0172', svv='-0.9161')
DELTA = {K: {T: float(x) for T, x in d.items()} for K, d in DELTA_TXT.items()}
GATE_C_FRAC = 0.25  # 0.05/0.187 of gate C of Stage 11
SIG_TOL, B_TOL, R_TOL = 1e-10, 1e-12, 1e-8
DEN_MAX = 0.5
CONSIST_TOL = 0.25
CLASS_P = ('out (< -1)', '-1', '-1/2', '0', 'out (> 0)')
CLASS_Q = ('out (< 0)', '0', '+1/2', '+1', 'out (> +1)')
VERDICT = {(1, 1): 'mechanism', (2, 2): 'amplitude falls, T not constant', (3, 3): 'no mechanism'}
PART_COLOR = {'Ls': '#2a78d6', 'Lu': '#eb6834', 'Hs': '#1baf7a', 'Hu': '#eda100'}  # categorical order of Stage 9


def write_json(root, name, obj):
    d = os.path.join(root, 'stage12')
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), 'w') as f:
        json.dump(obj, f, indent=1, default=lambda x: x.tolist() if hasattr(x, 'tolist') else str(x))


def read_json(root, name):
    p = os.path.join(root, 'stage12', name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def same_txt(x, txt):
    """x rounded as the text of the plan (its number of decimals and sign) is the text."""
    n = len(txt.split('.')[1])
    return ('%+.*f' if txt[0] in '+-' else '%.*f') % (n, x) == txt


# ----------------------------------------------------------------------------------------------- Q1: terms

def members(T):
    return SUMS.get(T, (T,))


def part_series(snaps, K, T):
    k = str(K)
    t = np.array([s['t'] for s in snaps])
    m = np.array([s['low'][k]['m'] for s in snaps])
    S = np.array([sum(s['low'][k]['S'][x] for x in members(T)) for s in snaps])
    return t, m, S


def part_I(snaps, K, T):
    """I_T = int S_T/m_K dt by the trapezoid over the records given (analyze_stage10.integ for one part)."""
    t, m, S = part_series(snaps, K, T)
    return float(np.sum(a9.trap_weights(t)*S/m))


def part_rate(snaps, K, T):
    t, m, S = part_series(snaps, K, T)
    return a11.tmean(t, np.abs(S/m))


def dense_part_errors(dense, K, T):
    """For the starts 0.60, ..., 0.65 tau0: I_T of the 5 records 0.1 tau0 apart minus I_T of all the records between."""
    out = []
    for p in range(a11.PHASES):
        i5 = part_I([dense['snaps'][p + a11.DENSE_STEP*j] for j in range(NREC)], K, T)
        ia = part_I(dense['snaps'][p:p + a11.DENSE_STEP*(NREC - 1) + 1], K, T)
        out.append(dict(err=i5 - ia, I5=i5, Ia=ia))
    return out


def part_test(I, r, delta):
    """analyze_stage11.test with delta_T,e = delta max(1, r_e); undecided if delta/|I_0.5| >= DEN_MAX."""
    out = a11.test(I, {e: delta*max(1.0, r[e]) for e in CASES}, delta)
    out['den'] = float(delta/abs(I['0.5'])) if I['0.5'] != 0 else float('inf')
    out['den_rule'] = bool(out['den'] >= DEN_MAX)
    if out['den_rule']:
        out['label'] = 'undecided'
    out['r'] = {e: float(r[e]) for e in CASES}
    return out


def terms(root):
    g11 = a11.read_json(root, 'gate.json')
    if g11 is None or not g11['PASS']:
        raise SystemExit('terms: no passed gates of Stage 11 for these records; nothing is decided')
    b = a11.load(root)
    dense = a10.read_json(root, 'budget.json')
    if dense is None:
        raise SystemExit('terms: no stage10/budget.json; gate C_T cannot be evaluated, so nothing is decided')
    c = b['0.5']['snaps']
    idx, dt = a11.same_time(np.array([s['t'] for s in c]), np.array([s['t'] for s in dense['snaps']]))
    gc = {}
    for T in PARTS:
        Ic, Id = part_I(c, K_MAIN, T), part_I([dense['snaps'][i] for i in idx], K_MAIN, T)
        gc[T] = dict(I_c=Ic, I_d=Id, diff=abs(Ic - Id), limit=GATE_C_FRAC*DELTA[K_MAIN][T])
        gc[T]['PASS'] = bool(gc[T]['diff'] <= gc[T]['limit'])
    sig = 0.0
    for e in CASES:
        for s in b[e]['snaps']:
            S = s['low'][str(K_MAIN)]['S']
            sig = max(sig, abs(sum(S[x] for x in PARTS) - a10.series([s], K_MAIN)[3][0])/sum(abs(S[x]) for x in PARTS))
    gate = dict(C=gc, dt_tau=dt, dense_idx=idx, PASS_C=bool(dt <= T_SLACK and all(gc[T]['PASS'] for T in PARTS)),
                sigma=sig, PASS_Sigma=bool(sig <= SIG_TOL))
    gate['PASS'] = gate['PASS_C'] and gate['PASS_Sigma']
    out = dict(gate=gate)
    print('gate C_T %s: max |t - t_dense| %.1e tau0 (limit %.0e)' % ('PASS' if gate['PASS_C'] else 'FAIL', dt, T_SLACK))
    for T in PARTS:
        print('  %-3s control %+9.4f  dense %+9.4f  |diff| %.4f  (limit %.4f)  %s'
              % (T, gc[T]['I_c'], gc[T]['I_d'], gc[T]['diff'], gc[T]['limit'], 'PASS' if gc[T]['PASS'] else 'FAIL'))
    print('gate Sigma %s: max |sum of the parts - S_svv|/sum|S_T| %.1e (limit %.0e; true by construction)'
          % ('PASS' if gate['PASS_Sigma'] else 'FAIL', sig, SIG_TOL))
    if not gate['PASS']:
        write_json(root, 'terms.json', out)
        print('gates: FAIL; no label (criteria of Stage 12)')
        return out
    k = {}
    for K in KS:
        k[str(K)] = {}
        for T in ALL:
            I = {e: part_I(b[e]['snaps'], K, T) for e in CASES}
            r0 = part_rate(b['0.5']['snaps'], K, T)
            tst = part_test(I, {e: part_rate(b[e]['snaps'], K, T)/r0 for e in CASES}, DELTA[K][T])
            k[str(K)][T] = dict(I=I, test=tst, one=a11.one_sided(I, tst['delta']), judged=K == K_MAIN and T in PARTS)
    out['k'] = k
    out['labels'] = {T: k[str(K_MAIN)][T]['test']['label'] for T in PARTS}
    out['sign_same'] = {e: {T: bool(np.sign(k[str(K_MAIN)][T]['I'][e]) == np.sign(k[str(K_MAIN)][T]['I']['0.5']))
                            for T in PARTS} for e in CASES}
    write_json(root, 'terms.json', out)
    print_terms(out)
    return out


def print_terms(out):
    for K in KS:
        print('K = %d%s:  part   I_0.25     I_0.5      I_1.0    r_0.25  r_1.0  delta_T     F       u       '
              '[lo, hi]        den   label' % (K, '' if K == K_MAIN else ' (reported)'))
        for T in ALL:
            o = out['k'][str(K)][T]
            t = o['test']
            print('        %-4s %+9.4f %+9.4f %+9.4f  %6.3f %6.3f %8.4f  %+7.3f %7.3f  [%+7.3f, %+7.3f]  %5.3f  %s%s'
                  % (T, o['I']['0.25'], o['I']['0.5'], o['I']['1.0'], t['r']['0.25'], t['r']['1.0'], t['delta']['0.5'],
                     t['F'], t['u'], t['lo'], t['hi'], t['den'], t['label'], '' if o['judged'] else ' (reported)'))
    print('labels (K = %d): %s' % (K_MAIN, out['labels']))
    print('same sign as the control (K = %d): %s' % (K_MAIN, out['sign_same']))


# --------------------------------------------------------------------------------------------- Q2: fields

def bands(op, s, u):
    """The arguments a_X,d of the four parts (module docstring) by active direction, and max u^2 of svv_split."""
    cs2 = op.cs2
    q = 0.5*sum(x**2 for x in u)
    sk, qk = op.fft(s), op.fft(q)
    uk = [op.fft(x) for x in u]
    a = {X: [] for X in PARTS}
    for d in range(3):
        if not op.active[d]:
            continue
        com = u[0]*op._grad(uk[0], d, True) + u[1]*op._grad(uk[1], d, True) + u[2]*op._grad(uk[2], d, True) \
            - op._grad(qk, d, True)
        a['Ls'].append((d, op._grad(sk, d, True)))
        a['Lu'].append((d, com/cs2))
        if op.filter_mom:
            a['Hs'].append((d, op._grad(sk, d, False)))
            a['Hu'].append((d, op._grad(qk, d, False)/cs2))
    return a, 2*float(q.max())


def assemble(op, a, rho, eps, umax2):
    """T_X = -+ sum_d G_d^*(coef_d rho a_X,d) on the real grid (zero if X has no argument)."""
    out = {}
    for X in PARTS:
        low, sgn = X[0] == 'L', 1.0 if X == 'Hu' else -1.0
        acc = np.zeros(op.mL.shape, dtype=complex)
        for d, ad in a[X]:
            if eps[d] == 0:
                continue
            coef = eps[d] if low else eps[d]/(op.cs2 + umax2)*op.cs2
            acc += sgn*op._adj(coef*rho*ad, d, low)
        out[X] = op.ifft(acc)
    return out


def shells(g):
    """Integer shells of |k| in units of 2 pi/L (a cubic box) and the weights of the r2c half spectrum."""
    n = np.rint(np.sqrt(g.ikx.reshape(-1, 1, 1)**2 + g.iky.reshape(1, 1, -1)**2 + g.ikz.reshape(1, -1, 1)**2))
    w = np.where((g.iky == 0) | (2*g.iky == g.nly), 1.0, 2.0).reshape(1, 1, -1)
    return n.astype(int), np.broadcast_to(w, n.shape)


def spectrum(op, f, sh):
    """sum over the shell of |f_k|^2 (the sum over the shells is <f^2>)."""
    n, w = sh
    return np.bincount(n.ravel(), weights=(w*np.abs(op.fft(f))**2).ravel(), minlength=int(n.max()) + 1)


def record_fields(g, op, sn, sh):
    """Amplitudes, norms, gate B, the sums of the 64 lowest cells, D_svv and the spectra of one record (Snap)."""
    eps = op.eps(sn.rho, sn.u, sn.b)
    T = a9.svv_split(op, sn.s, sn.u, sn.rho, eps)
    a, umax2 = bands(op, sn.s, sn.u)
    Tb = assemble(op, a, sn.rho, eps, umax2)
    B = {X: float(np.sqrt(np.sum((Tb[X] - T[X])**2)/np.sum(T[X]**2))) for X in PARTS}
    del Tb
    rdot, _, dsum = op.tend(sn.s, sn.u, sn.rho, eps)
    full = dict(T, svv=sum(T[X] for X in PARTS))
    low = a10.lowest(sn.rho, dict(T, rdot=rdot), (K_MAIN,))[str(K_MAIN)]
    idx = np.argpartition(sn.rho.ravel(), K_MAIN - 1)[:K_MAIN]
    amp2 = {X: sum(ad**2 for _, ad in a[X]) for X in PARTS}
    act = [d for d in range(3) if op.active[d]]
    return dict(rho_min=float(sn.rho.min()), B=B, low=low, umax2=umax2, eps=eps, c_dl=[x/op.eps0 for x in eps],
                A={X: float(np.sqrt(np.mean(amp2[X]))) for X in PARTS},
                A64={X: float(np.sqrt(np.mean(amp2[X].ravel()[idx]))) for X in PARTS},
                L2={X: float(np.sqrt(np.mean(v**2))) for X, v in full.items()},
                L1={X: float(np.mean(np.abs(v))) for X, v in full.items()},
                sum_rel={X: float(abs(v.sum())/np.abs(v).sum()) for X, v in full.items()},
                D=float(dsum*g.Lx*g.Ly*g.Lz/g.ntot),
                coef=dict(L=float(np.mean([eps[d] for d in act])),
                          H=float(np.mean([eps[d]*op.cs2/(op.cs2 + umax2) for d in act]))),
                spec=dict(s=spectrum(op, sn.s, sh), u=sum(spectrum(op, x, sh) for x in sn.u)))


def repro(ref, t, rec, low):
    """Gates R and Sigma2 of one record against the record of stage11/budget_<eps0>.json."""
    if ref is None:
        return dict(same=False, rel=float('inf'), sigma2=float('inf'))
    S11 = ref['low'][str(K_MAIN)]['S']
    tot = sum(abs(S11[x]) for x in PARTS)
    return dict(same=bool(ref['rec'] == rec and abs(ref['t'] - t) <= 1e-12*TAU0),
                rel=max(abs(low['S'][x] - S11[x]) for x in PARTS)/tot,
                sigma2=abs(low['S']['rdot'] - sum(S11[x] for x in PARTS))/tot)


def fields(root):
    for e in CASES:
        b11 = a11.read_json(root, 'budget_%s.json' % e)
        if b11 is None:
            raise SystemExit('fields: no stage11/budget_%s.json; gate R cannot be evaluated' % e)
        case = a11.case_dir(root, e)
        g, op = svv3d.from_case(case)
        if not g.Lx == g.Ly == g.Lz:
            raise SystemExit('fields: the shells assume a cubic box')
        sh = shells(g)
        kc = (g.nlx//2)*2/3
        recs = a11.records(case)
        print('%s: eps0 %g, records %d (expected %d)' % (CASES[e], op.eps0, len(recs), NREC), flush=True)
        res = dict(case=case, N=g.nlx, eps0=op.eps0, cs2=op.cs2, kc=kc,
                   m=svv3d.num(read_nml(os.path.join(case, 'calliope.in')).get('svv_mn', '0.5'))*kc, snaps=[])
        for i, (tt, path, rec) in enumerate(recs):
            f = record_fields(g, op, a7.Snap(g, op, path, rec), sh)
            f['R'] = repro(b11['snaps'][i] if i < len(b11['snaps']) else None, tt, rec, f['low'])
            res['snaps'].append(dict(t=tt, t_tau=tt/TAU0, path=path, rec=rec, **f))
            print('  t/tau0 %.4f  rho_min %.4f  B %.1e  R %s %.1e  Sigma2 %.1e  D %.4e  max sum_rel %.1e'
                  % (tt/TAU0, f['rho_min'], max(f['B'].values()), f['R']['same'], f['R']['rel'], f['R']['sigma2'],
                     f['D'], max(f['sum_rel'].values())), flush=True)
        write_json(root, 'fields_%s.json' % e, res)


# ---------------------------------------------------------------------------------------------- Q2: judge

def expo(x1, x0, base=4.0):
    return float(np.log(x1/x0)/np.log(base))


def klass(p, shift=0.0):
    """Index of CLASS_P (shift 0, p_A) or CLASS_Q (shift 1, q)."""
    x = p - shift
    if x <= -1.25:
        return 0
    if x <= -0.75:
        return 1
    if x <= -0.25:
        return 2
    if x < 0.25:
        return 3
    return 4


def decide(mean, recs, shift=0.0):
    c = klass(mean, shift)
    return dict(mean=float(mean), recs=[float(x) for x in recs],
                cls=c if all(klass(x, shift) == c for x in recs) else None)


def judge(root):
    f = {e: read_json(root, 'fields_%s.json' % e) for e in CASES}
    if any(v is None for v in f.values()):
        raise SystemExit('judge: a fields_<eps0>.json is missing; the step "fields" did not finish, so nothing is decided')
    sn = {e: f[e]['snaps'] for e in CASES}
    gate = dict(n={e: len(sn[e]) for e in CASES},
                B=max(max(s['B'].values()) for e in CASES for s in sn[e]),
                R_same=all(s['R']['same'] for e in CASES for s in sn[e]),
                R=max(s['R']['rel'] for e in CASES for s in sn[e]),
                Sigma2=max(s['R']['sigma2'] for e in CASES for s in sn[e]))
    gate['PASS'] = bool(all(n == NREC for n in gate['n'].values()) and gate['B'] <= B_TOL and gate['R_same'] and
                        gate['R'] <= R_TOL and gate['Sigma2'] <= SIG_TOL)
    print('gates %s: records %s (expected %d); B %.1e (limit %.0e); R same %s, %.1e (limit %.0e); Sigma2 %.1e '
          '(limit %.0e)' % ('PASS' if gate['PASS'] else 'FAIL', gate['n'], NREC, gate['B'], B_TOL, gate['R_same'],
                            gate['R'], R_TOL, gate['Sigma2'], SIG_TOL))
    if not gate['PASS']:
        write_json(root, 'summary.json', dict(gate=gate))
        raise SystemExit('judge: the gates did not pass; no class (criteria of Stage 12)')
    t = {e: np.array([s['t'] for s in sn[e]]) for e in CASES}

    def wm(e, key, X=None):
        return a11.tmean(t[e], [s[key] if X is None else s[key][X] for s in sn[e]])

    def ex(key, X=None, e1='1.0', e0='0.25', base=4.0):
        """The exponent of the window means and per record."""
        mean = expo(wm(e1, key, X), wm(e0, key, X), base)
        recs = [expo(a[key] if X is None else a[key][X], b[key] if X is None else b[key][X], base)
                for a, b in zip(sn[e1], sn[e0])]
        return mean, recs
    q2 = {}
    for X in PARTS:
        pA, qT, pc = ex('A', X), ex('L2', X), ex('coef', 'L' if X[0] == 'L' else 'H')
        dP, dQ = decide(*pA), decide(*qT, shift=1.0)
        if dP['cls'] is None or dQ['cls'] is None:
            verdict = 'undecided'
        else:
            verdict = VERDICT.get((dP['cls'], dQ['cls']), 'pair (%s, %s)' % (CLASS_P[dP['cls']], CLASS_Q[dQ['cls']]))
        phi = (qT[0] - pc[0], [x - y for x, y in zip(qT[1], pc[1])])
        cons = abs(phi[0] - pA[0]) <= CONSIST_TOL and all(abs(x - y) <= CONSIST_TOL for x, y in zip(phi[1], pA[1]))
        q2[X] = dict(p_A=dP, q=dQ, verdict=verdict, p_coef=dict(mean=pc[0], recs=pc[1]),
                     phi=dict(mean=phi[0], recs=phi[1]), consistent=bool(cons),
                     cls_p=None if dP['cls'] is None else CLASS_P[dP['cls']],
                     cls_q=None if dQ['cls'] is None else CLASS_Q[dQ['cls']],
                     rep=dict(L1=ex('L1', X)[0], A64=ex('A64', X)[0],
                              A_hi=ex('A', X, '1.0', '0.5', 2.0)[0], A_lo=ex('A', X, '0.5', '0.25', 2.0)[0],
                              L2_hi=ex('L2', X, '1.0', '0.5', 2.0)[0], L2_lo=ex('L2', X, '0.5', '0.25', 2.0)[0]))
    rep = dict(svv_L1=ex('L1', 'svv')[0], svv_L2=ex('L2', 'svv')[0], D=ex('D')[0], umax2=ex('umax2')[0],
               coef_L=ex('coef', 'L')[0], coef_H=ex('coef', 'H')[0],
               c_dl=expo(np.mean([wm('1.0', 'c_dl', d) for d in range(3)]),
                         np.mean([wm('0.25', 'c_dl', d) for d in range(3)])),
               means={e: dict(A={X: wm(e, 'A', X) for X in PARTS}, L2={X: wm(e, 'L2', X) for X in ALL if X in sn[e][0]['L2']},
                              L1={X: wm(e, 'L1', X) for X in PARTS + ('svv',)}, A64={X: wm(e, 'A64', X) for X in PARTS},
                              D=wm(e, 'D'), coef={b: wm(e, 'coef', b) for b in ('L', 'H')},
                              c_dl=[wm(e, 'c_dl', d) for d in range(3)], umax2=wm(e, 'umax2')) for e in CASES})
    spec = {e: {v: np.sum(a9.trap_weights(t[e])[:, None]*np.array([s['spec'][v] for s in sn[e]]), axis=0)
                / np.sum(a9.trap_weights(t[e])) for v in ('s', 'u')} for e in CASES}
    s = dict(gate=gate, q2=q2, rep=rep, spec=spec, kc=f['0.5']['kc'], m=f['0.5']['m'])
    write_json(root, 'summary.json', s)
    print_judge(s)
    tm = read_json(root, 'terms.json')
    figures(root, s, tm if tm is not None and tm['gate']['PASS'] else None)
    return s


def print_judge(s):
    print('part  p_A (mean; records)                         class      q (mean; records)                         '
          'class     verdict')
    for X in PARTS:
        o = s['q2'][X]
        print('%-4s %+.3f; %s  %-9s %+.3f; %s  %-9s %s'
              % (X, o['p_A']['mean'], ' '.join('%+.3f' % x for x in o['p_A']['recs']), o['cls_p'], o['q']['mean'],
                 ' '.join('%+.3f' % x for x in o['q']['recs']), o['cls_q'], o['verdict']))
    for X in PARTS:
        o = s['q2'][X]
        r = o['rep']
        print('%-4s p_coef %+.3f  phi %+.3f  |phi - p_A| %.3f  (records max %.3f)  consistent %s;  L1 %+.3f  A64 %+.3f  '
              'A hi/lo %+.3f %+.3f  L2 hi/lo %+.3f %+.3f'
              % (X, o['p_coef']['mean'], o['phi']['mean'], abs(o['phi']['mean'] - o['p_A']['mean']),
                 max(abs(x - y) for x, y in zip(o['phi']['recs'], o['p_A']['recs'])), o['consistent'], r['L1'],
                 r['A64'], r['A_hi'], r['A_lo'], r['L2_hi'], r['L2_lo']))
    r = s['rep']
    print('reported: svv L1 %+.3f L2 %+.3f; D_svv %+.3f; max u^2 %+.3f; coef L %+.3f H %+.3f; c_dl %+.3f'
          % (r['svv_L1'], r['svv_L2'], r['D'], r['umax2'], r['coef_L'], r['coef_H'], r['c_dl']))
    for e in CASES:
        m = r['means'][e]
        print('  %5s A %s  L2 %s  D %.4e' % (e, ' '.join('%s %.4e' % (X, m['A'][X]) for X in PARTS),
                                            ' '.join('%s %.4e' % (X, v) for X, v in m['L2'].items()), m['D']))


def figures(root, s, tm):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = os.path.join(root, 'stage12', 'fig')
    os.makedirs(d, exist_ok=True)
    grid = dict(alpha=0.35, lw=0.5, color='#e1e0d9')
    if tm is not None:  # Q1: I_T against eps0 with the lines of proportional and constant
        fig, axs = plt.subplots(2, 2, figsize=(8.4, 6.6))
        x = np.linspace(0.0, 1.05, 2)
        for ax, T in zip(axs.ravel(), PARTS):
            o = tm['k'][str(K_MAIN)][T]
            I5 = o['I']['0.5']
            ax.plot(x, I5*x/0.5, color=COLOR['ref'], lw=1.4, ls='--')
            ax.plot(x, I5 + 0*x, color=COLOR['ref'], lw=1.4, ls=':')
            for e in CASES:  # the abscissa names the case; error bars: delta_T,e
                ax.errorbar(EPS[e], o['I'][e], yerr=o['test']['delta'][e], color=COLOR[e], marker='o', ms=7, lw=2.0,
                            capsize=4, mec='white', mew=1.5)
            ax.axhline(0, color='#898781', lw=0.8)
            ax.set_title('%s: F = %.2f, u = %.2f, %s' % (T, o['test']['F'], o['test']['u'], o['test']['label']),
                         fontsize=10)
            ax.set_xlabel('svv_eps0')
            ax.set_ylabel('I_%s (K = %d)' % (T, K_MAIN))
            ax.grid(**grid)
        axs[0, 0].annotate('proportional (dashed), constant (dotted)', (0.02, 0.96), xycoords='axes fraction',
                           va='top', fontsize=8, color='#55534e')
        fig.tight_layout()
        fig.savefig(os.path.join(d, 'stage12_terms.png'), dpi=150, bbox_inches='tight')
        plt.close(fig)
    # Q2: amplitudes and norms relative to the control, with the slopes of the classes
    fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.9))
    ex = np.array([0.25, 1.0])
    for ax, key, slopes, lab in ((axs[0], 'A', (-1.0, -0.5, 0.0), 'A_X / A_X(0.5)'),
                                 (axs[1], 'L2', (0.0, 0.5, 1.0), '||T_X||_2 / ||T_X||_2(0.5)')):
        for p, ls in zip(slopes, ('--', '-.', ':')):
            ax.plot(ex, (ex/0.5)**p, color=COLOR['ref'], lw=1.2, ls=ls)
            ax.annotate('%+g' % p if p else '0', (1.0, (1.0/0.5)**p), xytext=(4, 0), textcoords='offset points',
                        va='center', fontsize=8, color='#55534e')
        for X in PARTS:
            y = [s['rep']['means'][e][key][X]/s['rep']['means']['0.5'][key][X] for e in CASES]
            ax.plot([EPS[e] for e in CASES], y, color=PART_COLOR[X], lw=2.0, marker='o', ms=7, mec='white', mew=1.5,
                    label=X)
        ax.set_xscale('log', base=2)
        ax.set_yscale('log', base=2)
        ax.set_xlabel('svv_eps0')
        ax.set_ylabel(lab)
        ax.grid(**grid)
        ax.legend(frameon=False, fontsize=8, loc='best')
    fig.tight_layout()
    fig.savefig(os.path.join(d, 'stage12_bands.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    # spectra of s and u (window means) and their ratios to the control; dotted: m and kc of the SVV
    fig, axs = plt.subplots(2, 2, figsize=(9.0, 6.6), sharex=True)
    for j, v in enumerate(('s', 'u')):
        k = np.arange(len(s['spec']['0.5'][v]))
        sel = k >= 1
        for e in CASES:
            axs[0, j].loglog(k[sel], np.asarray(s['spec'][e][v])[sel], color=COLOR[e], lw=1.6, label='eps0 = %s' % e)
            if e != '0.5':
                axs[1, j].semilogx(k[sel], np.asarray(s['spec'][e][v])[sel]/np.asarray(s['spec']['0.5'][v])[sel],
                                   color=COLOR[e], lw=1.6, label='eps0 = %s' % e)
        axs[1, j].axhline(1, color='#898781', lw=0.8)
        for ax in axs[:, j]:
            for kk in (s['m'], s['kc']):
                ax.axvline(kk, color='#898781', lw=0.8, ls=':')
            ax.grid(**grid)
        axs[0, j].set_ylabel('E_%s(k)' % v)
        axs[1, j].set_ylabel('E_%s / E_%s(0.5)' % (v, v))
        axs[1, j].set_xlabel('|k| L/2 pi')
        axs[0, j].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(d, 'stage12_spectra.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


# ------------------------------------------------------------------------------------------------- tests

def synthetic_budgets(root, c, shift=0.0):
    """stage11/budget_<eps0>.json, gate.json and stage10/budget.json with constant m and S, so that I_T,e = c[T](e)
    (the window is 0.4 tau0); shift adds to S_Ls of the dense records."""
    W = 0.4*TAU0

    def snaps(ts, e, sh=0.0):
        out = []
        for i, tt in enumerate(ts):
            S = {T: (c[T](EPS[e]) + (sh if T == 'Ls' else 0.0))/W for T in PARTS}
            S['adv'] = 0.0
            out.append(dict(t=float(tt), t_tau=float(tt/TAU0), rec=i,
                            low={str(K): dict(m=1.0, rho_top=0.1, S=dict(S)) for K in KS}))
        return out
    ts = TAU0*np.array([0.5999, 0.6998, 0.7999, 0.8999, 0.9999])
    for e in CASES:
        a11.write_json(root, 'budget_%s.json' % e, dict(eps0=EPS[e], snaps=snaps(ts, e)))
    a11.write_json(root, 'gate.json', dict(PASS=True))
    a10.write_json(root, 'budget.json', dict(snaps=snaps(TAU0*(0.5999 + 0.01*np.arange(46)), '0.5', shift)))


def synthetic_fields(root, p, q, bump=None, bad_B=False):
    """fields_<eps0>.json with A_X = (eps0/0.5)^p_X and ||T_X||_2 = (eps0/0.5)^q_X (times a slow factor of t);
    bump = (X, record, factor) multiplies A_X of e1.0 at one record."""
    for e in CASES:
        out = []
        for i in range(NREC):
            tt = TAU0*(0.6 + 0.1*i)
            fac = 1.0 + 0.1*i
            A = {X: fac*(EPS[e]/0.5)**p[X] for X in PARTS}
            if bump is not None and e == '1.0' and i == bump[1]:
                A[bump[0]] *= bump[2]
            L2 = {X: fac*(EPS[e]/0.5)**q[X] for X in PARTS}
            L2['svv'] = 0.1*fac
            out.append(dict(t=tt, t_tau=tt/TAU0, rec=i, A=A, A64=dict(A), L2=L2, L1=dict(L2), D=EPS[e]**0.0,
                            B={X: 1e-6 if bad_B and i == 2 else 0.0 for X in PARTS},
                            R=dict(same=True, rel=0.0, sigma2=0.0), coef=dict(L=0.02*EPS[e], H=0.01*EPS[e]),
                            c_dl=[0.02, 0.02, 0.02], umax2=1.0, spec=dict(s=[0.0, 1.0, 0.5], u=[0.0, 2.0, 1.0])))
        write_json(root, 'fields_%s.json' % e, dict(eps0=EPS[e], snaps=out, kc=1.5, m=0.75))


def selftest():
    ok = True

    def rep(name, cond, info=''):
        nonlocal ok
        ok &= bool(cond)
        print('  %-58s %s  %s' % (name, 'PASS' if cond else 'FAIL', info))
    # 1. delta_T, I5 and the dense integral of the plan from the dense rerun of Stage 10
    dense = a10.read_json(RUNS, 'budget.json')
    rep('stage10/budget.json present', dense is not None and len(dense['snaps']) == 46, RUNS)
    if dense is not None:
        for K in KS:
            bad = []
            for T in ALL:
                er = dense_part_errors(dense, K, T)
                mx = max(abs(x['err']) for x in er)
                if not same_txt(mx, DELTA_TXT[K][T]) or (K == K_MAIN and not (
                        same_txt(er[0]['I5'], I5_TXT[T]) and same_txt(er[0]['Ia'], DENSE_TXT[T]))):
                    bad.append('%s %.4f %+.4f %+.4f' % (T, mx, er[0]['I5'], er[0]['Ia']))
            rep('delta_T of K = %d (and I5, dense of K = 64), 9 kinds' % K, not bad, ' '.join(bad))
        er = dense_part_errors(dense, K_MAIN, 'svv')
        rep('delta of svv = 0.187 of Stage 11', round(max(abs(x['err']) for x in er), 3) == a11.DELTA[K_MAIN])
    # 2. the labels of the plan (I_0.5 and delta of the table, r as given)

    def I(a, b, c):
        return {'0.25': a, '0.5': b, '1.0': c}
    r1 = I(1.0, 1.0, 1.0)
    d64 = DELTA[K_MAIN]
    p = part_test(I(26.8, 53.5584, 107.1), I(0.5, 1.0, 2.0), d64['Ls'])
    c = part_test(I(53.5584, 53.5584, 53.5584), r1, d64['Ls'])
    rep('plan, Ls: F = 1 -> u 0.226 [0.77, 1.23] prop; F = 0 -> 0.086',
        round(p['u'], 3) == 0.226 and round(p['lo'], 2) == 0.77 and round(p['hi'], 2) == 1.23 and p['label'] == 'prop'
        and round(c['u'], 3) == 0.086 and c['label'] == 'const', 'F %.5f u %.4f; u %.4f' % (p['F'], p['u'], c['u']))
    p = part_test(I(3.4183, 6.8366, 13.67), I(0.5, 1.0, 2.0), d64['Hs'])
    c = part_test(I(6.8366, 6.8366, 6.8366), r1, d64['Hs'])
    rep('plan, Hs: F = 1 -> u 0.83 [0.17, 1.83] prop; F = 0 -> 0.32',
        round(p['u'], 2) == 0.83 and round(p['lo'], 2) == 0.17 and round(p['hi'], 2) == 1.83 and p['label'] == 'prop'
        and round(c['u'], 2) == 0.32 and c['label'] == 'const', 'F %.5f u %.4f; u %.4f' % (p['F'], p['u'], c['u']))
    h = part_test(I(0.1629/2, 0.1629, 0.1629*2), r1, d64['H'])
    rep('denominator rule: H (0.3214/0.1629 = 1.97) undecided', h['den_rule'] and h['label'] == 'undecided' and
        abs(h['F'] - 1) < 1e-12, 'den %.3f' % h['den'])
    rep('denominators of the four parts < 0.5', all(d64[T]/abs(float(I5_TXT[T])) < DEN_MAX for T in PARTS),
        ' '.join('%s %.3f' % (T, d64[T]/abs(float(I5_TXT[T]))) for T in PARTS))
    cases = [(part_test(I(25.0, 50.0, 75.0), r1, 1.0), 'mid'), (part_test(I(10.0, 50.0, 200.0), I(1, 1, 4), 1.0), 'out'),
             (part_test(I(60.0, 50.0, 40.0), r1, 1.0), 'out'), (part_test(I(30.0, 50.0, 70.0), r1, 20.0), 'undecided')]
    rep('labels mid, out (two ways), undecided', all(x['label'] == y for x, y in cases),
        ' '.join('%s(F %+.2f u %.2f)' % (x['label'], x['F'], x['u']) for x, _ in cases))
    rep('delta_T,e = delta_T max(1, r)', p['delta'] == {'0.25': d64['Hs'], '0.5': d64['Hs'], '1.0': 2*d64['Hs']})
    # 3. the classes at their edges and the rule of the 5 records
    kp = [klass(x) for x in (-1.25, -1.2499, -0.75, -0.7499, -0.25, -0.2499, 0.2499, 0.25)]
    kq = [klass(x, 1.0) for x in (-0.25, -0.2499, 0.25, 0.2501, 0.75, 0.7501, 1.2499, 1.25)]
    rep('classes of p_A and q at the edges', kp == [0, 1, 1, 2, 2, 3, 3, 4] and kq == [0, 1, 1, 2, 2, 3, 3, 4],
        '%s %s' % (kp, kq))
    rep('5 records: one other class -> undecided', decide(-1.0, [-1.0]*5)['cls'] == 1 and
        decide(-1.0, [-1.0]*4 + [-0.7])['cls'] is None and decide(0.0, [0.1]*5, 1.0)['cls'] == 1)
    # 4. a synthetic state on 32^3: the arguments rebuild svv_split, scaling, Parseval, record_fields
    g = gen_ic.Grid(32, 32, 32)
    kl = 2*np.pi/g.Lx
    rng = np.random.default_rng(12)
    rho = 1 - 0.6*np.exp(2*(np.cos(kl*g.X) + np.cos(kl*g.Y) + np.cos(kl*g.Z) - 3)) + 0.01*rng.standard_normal(g.X.shape)
    m = [0.1*np.sin(kl*g.X) + 0.01*rng.standard_normal(g.X.shape), 0.05*np.sin(kl*g.Z), 0.05*np.sin(kl*g.Y)]
    b = [0.1*np.cos(kl*g.Y), 0*rho, 1 + 0*rho]
    op = svv3d.SVV(g, 0.7, 0.5, 0.5, True, True)  # cs2 != 1: the order of the products differs
    sn = a10.State(op, rho, m, b)
    eps = op.eps(rho, sn.u, b)
    T = a9.svv_split(op, sn.s, sn.u, rho, eps)
    a, umax2 = bands(op, sn.s, sn.u)
    Tb = assemble(op, a, rho, eps, umax2)
    dif = {X: float(np.sqrt(np.sum((Tb[X] - T[X])**2)/np.sum(T[X]**2))) for X in PARTS}
    rep('arguments rebuild svv_split (relative L2)', max(dif.values()) <= B_TOL,
        ' '.join('%s %.1e' % x for x in dif.items()))
    T2 = assemble(op, a, rho, [2*x for x in eps], umax2)
    rep('T_X doubles with eps0 (arguments unchanged)',
        all(np.abs(T2[X] - 2*Tb[X]).max() <= 1e-13*np.abs(Tb[X]).max() for X in PARTS))

    def amp(aa):
        return {X: float(np.sqrt(np.mean(sum(ad**2 for _, ad in aa[X])))) for X in PARTS}
    A0, As, Au = amp(a), amp(bands(op, 2*sn.s, sn.u)[0]), amp(bands(op, sn.s, [2*x for x in sn.u])[0])
    rep('s x 2: A_Ls, A_Hs x 2; u x 2: A_Lu, A_Hu x 4', all(
        abs(As[X]/A0[X] - f1) < 1e-12 and abs(Au[X]/A0[X] - f2) < 1e-12
        for X, f1, f2 in (('Ls', 2, 1), ('Hs', 2, 1), ('Lu', 1, 4), ('Hu', 1, 4))), ' '.join('%s %.3e' % x for x in A0.items()))
    sh = shells(g)
    es, eu = spectrum(op, sn.s, sh), sum(spectrum(op, x, sh) for x in sn.u)
    rep('Parseval of the shell spectra', abs(es.sum()/np.mean(sn.s**2) - 1) < 1e-12 and
        abs(eu.sum()/np.mean(sum(x**2 for x in sn.u)) - 1) < 1e-12, 'shells %d' % len(es))
    rf = record_fields(g, op, sn, sh)
    S = rf['low']['S']
    rep('record_fields: gate B, Sigma2 on the state, sums vanish, D > 0',
        max(rf['B'].values()) <= B_TOL and abs(S['rdot'] - sum(S[X] for X in PARTS)) <= SIG_TOL*sum(abs(S[X]) for X in PARTS)
        and max(rf['sum_rel'].values()) < 1e-12 and rf['D'] > 0 and rf['A'] == A0,
        'Sigma2 %.1e, max sum_rel %.1e, D %.3e' % (abs(S['rdot'] - sum(S[X] for X in PARTS))/sum(abs(S[X]) for X in PARTS),
                                                    max(rf['sum_rel'].values()), rf['D']))
    ref = dict(rec=3, t=1.0, low={str(K_MAIN): dict(S={X: S[X] for X in PARTS})})
    rr = repro(ref, 1.0, 3, rf['low'])
    rep('repro: the same record passes R and Sigma2', rr['same'] and rr['rel'] == 0 and rr['sigma2'] <= SIG_TOL)
    # 5. the chains on synthetic files
    cfun = dict(Ls=lambda x: 100*x, Lu=lambda x: -45.0 + 0*x, Hs=lambda x: 6.8*np.sqrt(x/0.5), Hu=lambda x: -20*(x/0.5)**2)
    pp = dict(Ls=-1.0, Lu=-0.5, Hs=0.0, Hu=-1.0)
    qq = dict(Ls=0.0, Lu=0.5, Hs=1.0, Hu=0.6)
    with tempfile.TemporaryDirectory() as root, contextlib.redirect_stdout(io.StringIO()):
        synthetic_budgets(root, cfun)
        tm = terms(root)
        synthetic_budgets(root, cfun, shift=0.25*DELTA[K_MAIN]['Ls'] + 0.01)
        tmc = terms(root)
        synthetic_budgets(root, cfun)
        terms(root)
        synthetic_fields(root, pp, qq)
        sm = judge(root)
        figs = sorted(os.listdir(os.path.join(root, 'stage12', 'fig')))
        synthetic_fields(root, pp, qq, bump=('Ls', 2, 4**0.6))
        smb = judge(root)
        synthetic_fields(root, pp, qq, bad_B=True)
        try:
            judge(root)
            stopped = False
        except SystemExit:
            stopped = True
    rep('chain terms: labels prop, const, undecided, out', tm['gate']['PASS'] and tm['labels'] ==
        dict(Ls='prop', Lu='const', Hs='undecided', Hu='out'), '%s' % tm.get('labels'))
    rep('chain terms: gate C_T fails on a shifted dense rerun, no label', not tmc['gate']['PASS_C'] and 'k' not in tmc)
    v = {X: sm['q2'][X]['verdict'] for X in PARTS}
    rep('chain judge: verdicts', v == dict(Ls='mechanism', Lu='amplitude falls, T not constant', Hs='no mechanism',
                                           Hu='pair (-1, +1/2)') and all(sm['q2'][X]['consistent'] for X in PARTS[:3])
        and not sm['q2']['Hu']['consistent'] and figs == ['stage12_bands.png', 'stage12_spectra.png',
                                                          'stage12_terms.png'], '%s' % v)
    rep('chain judge: one record of another class -> undecided', smb['q2']['Ls']['verdict'] == 'undecided' and
        smb['q2']['Lu']['verdict'] == v['Lu'])
    rep('chain judge: gate B fails, judge stops', stopped)
    print('selftest:', 'all PASS' if ok else 'FAIL')
    return ok


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'selftest':
        sys.exit(0 if selftest() else 1)
    {'terms': terms, 'fields': fields, 'judge': judge}[cmd](sys.argv[2])
