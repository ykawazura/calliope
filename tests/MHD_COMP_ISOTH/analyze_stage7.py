#!/usr/bin/env python3
"""Stage 7: the holes of rho in turb1 (M_s = 1) on 64^3, 128^3 and 256^3, from the out3d records (and
conservation.dat) of the existing runs only.

    python3 analyze_stage7.py selftest          synthetic checks of the tools below
    python3 analyze_stage7.py tail ROOT N       density PDF and its tail          -> ROOT/stage7/tail_N<N>.json
    python3 analyze_stage7.py hole ROOT N       structure of the deepest holes     -> ROOT/stage7/hole_N<N>.json
    python3 analyze_stage7.py judge ROOT        criteria, comparison across N, figures
                                                                -> ROOT/stage7/summary.json, ROOT/stage7/fig/*.png

ROOT is the runs directory (CASES below). out3d rho is the FFT of exp(s) without the mask (fields.F90
lnrho_to_rho), so the grid rho is exact and s = ln rho; m and B are the (masked) state.

tail: rho_min, rho_max, sigma_s, the volume PDF of s (bins of 0.1), F_V(rho < theta), F_M(rho < theta) and the
number of cells n(rho < theta) for THETA, and the number of local minima (26 neighbours, periodic) below
NHOLE_THR, at TAIL_T (all N) and TAIL_T256 (256^3).

hole: at HOLE_T, the min of rho and the NHOLE deepest local minima with rho < 0.1 at least SEP cells from
each other. At each hole x0 (b = B/|B| at x0, cells as the unit of length):
    point values; D = div u, D_par = b b : grad u, D_perp = D - D_par (spectral, k_d of the code);
    background: medians over the shell SHELL_R of cells around x0;
    the curvature widths w_i = sqrt(2/lambda_i) of the Hessian of s (central differences), the thin axis n;
    the region rho < sqrt(rho_min rho_bg) connected to x0 (26 neighbours, +-REG_R cells): cells, semi-axes
    sqrt(5 var) (a filled ellipsoid), the angle of its long axis to b;
    cuts (+-CUT_L steps) along the lattice direction (of DIRS) nearest to b and to n;
    Pi = (P_tot,bg - P_tot(x0))/(P_th,bg - P_th(x0)), Lambda = (|B|/rho)(x0)/(|B|/rho)_bg,
    chi_par = ln(rho_bg,par/rho(x0))/(du_par/2c_s) over +-CHI_L steps of the cut along b (u_par = u.b(x0));
    the band split of the depth of s (L: |k_d| <= m_N in every direction, M: the rest of the 2/3 mask,
    H: outside the mask; as svv3d.SVV); the odd-even check of s along the 3 axes;
    d rho/dt = -div_N m + SVV (svv3d.SVV.eps and tend with the weight rho, as the code).
256^3 also: the move of the min between neighbouring snapshots against |u| dt.

judge: the criteria fixed before the analysis (tasks/todo.md, Stage 7; REPORT_stage7.md):
    tail, for each pair (coarse, fine) and time, on the theta with n >= TAIL_NMIN on the coarse grid:
        converged   all |F_V,fine/F_V,coarse - 1| <= 0.25 and n_fine(rho < rho_min,coarse) <= 100
        fattening   F_V,fine/F_V,coarse >= 1.25 for 2 or more theta <= 0.2
        intermediate otherwise; over the times: fattening if any time is, converged if all are
    hole, each: balance |Pi| <= 0.3; drain |D_par| >= |D_perp|, Lambda in [0.5, 2] rho_bg/rho_min and
        chi_par in [0.5, 2]; perp |D_perp| > |D_par| and Lambda in [0.5, 2]; suspect (numerical origin) if
        none of them holds or s oscillates odd-even; 'many' suspect holes: half of the 256^3 holes or more
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_ic  # noqa: E402
import svv3d  # noqa: E402
from analyze_stage1 import load_runs, num, read_nml  # noqa: E402
import analyze_stage4 as s4  # noqa: E402

CASES = {64: 'stage4/turb/turb1_N64', 128: 'stage4/turb/turb1_N128_a', 256: 'stage6q/turb/turb1_N256'}
TAU0 = np.pi  # pi/M_s (case.json)
T_TOL = 0.003  # tau0; the record nearest to a requested time
TAIL_T = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
TAIL_T256 = (0.578, 0.700, 0.818, 0.893, 0.9065, 0.921, 0.965)
HOLE_T = {64: (0.5, 1.0), 128: (0.5, 1.0), 256: (0.5, 0.578, 0.700, 0.818, 0.893, 0.9065, 0.921, 0.965, 1.0)}
THETA = (0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02)
NHOLE_THR = (0.2, 0.1, 0.05)
PDF_EDGES = np.round(np.arange(-8.0, 3.0001, 0.1), 10)
NLOW = 2000  # smallest values of rho kept for n(rho < x) at any x
TAIL_NMIN = 100
HOLE_THR, NHOLE, SEP = 0.1, 4, 8.0
SHELL_R = (4.0, 8.0)
REG_R, CUT_L, CHI_L, PROF_L = 16, 8, 4, 16
OE_FRAC, OE_RUN = 0.05, 4
# 13 lattice directions (x, y, z): axes, face and body diagonals
DIRS = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1),
        (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]
# the min of rho of the restart copies (runs/stage6q/ana/state256.log, state256b.log; 5 significant digits)
CHECK = {(256, 0.578): dict(t=1.816833, imin=(168, 55, 240), adv=1.0099, svv=-0.63503),
         (256, 0.893): dict(t=2.805711, imin=(34, 22, 7), adv=-2.9441, svv=3.3421)}
COLOR_N = {64: '#2a78d6', 128: '#eb6834', 256: '#1baf7a'}


# ------------------------------------------------------------------------------------------- reading

def field_at(g, path, rec, q):
    """One field of a record of out3d/ (the layout of analyze_stage4.out3d_at)."""
    n = g.nky*g.nkz*g.nkx
    nrec = len(np.atleast_2d(np.loadtxt(os.path.join(path, 'time.dat'))))
    a = np.memmap(os.path.join(path, q + '.dat'), dtype='<c16', mode='r')
    assert len(a) in (nrec*n, (2*nrec - 1)*n), '%s: %d values for %d records' % (path, len(a), nrec)
    j = rec*(2 if len(a) == (2*nrec - 1)*n else 1)
    return np.array(a[j*n:(j + 1)*n]).reshape(g.nky, g.nkz, g.nkx).transpose(2, 1, 0)


def pick(case, t_tau):
    """(t, out3d directory, record) nearest to t_tau tau0 (within T_TOL), or None."""
    recs = s4.out3d_records(case)
    best = min(recs, key=lambda r: abs(r[0]/TAU0 - t_tau))
    return best if abs(best[0]/TAU0 - t_tau) <= T_TOL else None


class Snap:
    """The state of one out3d record in real space (svv3d.real_state from rho instead of s)."""

    def __init__(self, g, op, path, rec, full=True):
        self.rho = op.ifft(field_at(g, path, rec, 'rho'))
        if not (np.all(np.isfinite(self.rho)) and self.rho.min() > 0):
            raise ValueError('%s record %d: rho not finite or not positive' % (path, rec))
        self.s = np.log(self.rho)
        if full:
            self.mk = [field_at(g, path, rec, q) for q in ('mx', 'my', 'mz')]
            m = [op.ifft(x) for x in self.mk]
            self.u = [x/self.rho for x in m]
            self.b = [op.ifft(field_at(g, path, rec, q)) for q in ('bx', 'by', 'bz')]


def cons_table(case):
    return np.concatenate(load_runs(case, 'conservation.dat'))


def cons_check(cons, t, rho_min):
    """rho_min of the conservation.dat row nearest to t."""
    i = int(np.argmin(np.abs(cons[:, 0] - t)))
    return dict(dt=float(cons[i, 0] - t), rho_min=float(cons[i, 11]), rel=float(rho_min/cons[i, 11] - 1))


# ------------------------------------------------------------------------------------ grid helpers
# array index p = (i, k, j) = (x, z, y); offsets o = (dx, dy, dz) in cells

def at(f, p, o):
    n = f.shape
    return f[(p[0] + o[0]) % n[0], (p[1] + o[2]) % n[1], (p[2] + o[1]) % n[2]]


def cut(f, p, d, L):
    """f at p + k d, k = -L..L (d a lattice direction (dx, dy, dz))."""
    k = np.arange(-L, L + 1)
    n = f.shape
    return f[(p[0] + k*d[0]) % n[0], (p[1] + k*d[2]) % n[1], (p[2] + k*d[1]) % n[2]]


def shell_offsets(r1=SHELL_R[0], r2=SHELL_R[1]):
    r = np.arange(-int(r2), int(r2) + 1)
    o = np.array(np.meshgrid(r, r, r, indexing='ij')).reshape(3, -1).T
    d = np.sqrt((o**2).sum(1))
    return o[(d >= r1) & (d <= r2)]


SHELL = shell_offsets()


def shell_vals(f, p, off=SHELL):
    n = f.shape
    return f[(p[0] + off[:, 0]) % n[0], (p[1] + off[:, 2]) % n[1], (p[2] + off[:, 1]) % n[2]]


def pdist(p, q, n):
    """Periodic distance in cells between array indices on the same grid."""
    d = np.abs(np.array(p) - np.array(q))
    d = np.minimum(d, n - d)
    return float(np.sqrt((d**2).sum()))


def phys(p, dl):
    """(x, y, z) of the array index p."""
    return [p[0]*dl, p[2]*dl, p[1]*dl]


def local_min(f):
    """Strict local minima over the 26 neighbours (periodic)."""
    from scipy.ndimage import minimum_filter
    fp = np.ones((3, 3, 3), dtype=bool)
    fp[1, 1, 1] = False
    return f < minimum_filter(f, footprint=fp, mode='wrap')


def nearest_dir(v):
    """The lattice direction of DIRS nearest to the unit vector v (x, y, z), oriented along v, and the angle."""
    best, c = None, -1.0
    for d in DIRS:
        cd = float(np.dot(d, v))/np.sqrt(np.dot(d, d))
        if abs(cd) > c:
            best, c = (d if cd >= 0 else tuple(-x for x in d)), abs(cd)
    return tuple(int(x) for x in best), float(np.degrees(np.arccos(min(1.0, c))))


def hess(f, p):
    """Hessian of f at p by central differences (cells), (x, y, z)."""
    e = np.eye(3, dtype=int)
    c = at(f, p, (0, 0, 0))
    H = np.zeros((3, 3))
    for a in range(3):
        H[a, a] = at(f, p, e[a]) - 2*c + at(f, p, -e[a])
        for b in range(a + 1, 3):
            H[a, b] = H[b, a] = (at(f, p, e[a] + e[b]) - at(f, p, e[a] - e[b]) - at(f, p, e[b] - e[a])
                                 + at(f, p, -e[a] - e[b]))/4
    return H


def region(rho, p, thr, R=REG_R):
    """The cells with rho < thr connected to p (26 neighbours) within +-R: count, semi-axes sqrt(5 var)
    (descending), their axes (x, y, z) and whether the region reaches the edge of the cube."""
    from scipy.ndimage import label
    n = rho.shape
    idx = [(p[a] + np.arange(-R, R + 1)) % n[a] for a in range(3)]
    sub = rho[np.ix_(*idx)]
    lab, _ = label(sub < thr, structure=np.ones((3, 3, 3), dtype=int))
    c = np.argwhere(lab == lab[R, R, R]) - R
    xyz = c[:, [0, 2, 1]].astype(float)
    edge = bool(np.any(np.abs(c) == R))
    if len(xyz) < 2:
        return dict(cells=int(len(xyz)), semi=[0.0, 0.0, 0.0], axes=np.eye(3).tolist(), edge=edge)
    ev, evec = np.linalg.eigh(np.cov(xyz.T, bias=True))
    o = np.argsort(ev)[::-1]
    return dict(cells=int(len(xyz)), semi=[float(np.sqrt(5*max(ev[i], 0.0))) for i in o],
                axes=evec[:, o].T.tolist(), edge=edge)


def oddeven(s, p, L=CUT_L, frac=OE_FRAC):
    """The longest run of sign alternations of the second difference of s along each axis (+-L), counting only
    the points with |d2| >= frac |d2(x0)|."""
    runs = []
    for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        c = cut(s, p, d, L + 1)
        d2 = c[2:] - 2*c[1:-1] + c[:-2]
        thr = frac*abs(d2[L])
        run = best = 0
        for i in range(len(d2) - 1):
            if abs(d2[i]) >= thr and abs(d2[i + 1]) >= thr and d2[i]*d2[i + 1] < 0:
                run += 1
                best = max(best, run)
            else:
                run = 0
        runs.append(best)
    return runs


def band_masks(g, mN):
    """L (|k_d| <= m_d = m_N k_c,d in every direction, where the kernel of svv3d.SVV is 0), M (the rest of the
    2/3 mask) and H (outside it)."""
    L = np.ones(g.mask.shape, dtype=bool)
    for k, n, Lb in ((g.kx, g.nlx, g.Lx), (g.ky, g.nly, g.Ly), (g.kz, g.nlz, g.Lz)):
        kc = 2*np.pi*(n//2)/Lb*2/3
        L = L & (np.abs(k) <= mN*kc)
    L = L.astype(float)
    return dict(L=L, M=g.mask - L, H=1 - g.mask)


def grad_at(op, u, pts):
    """G[p, i, j] = d u_i/d x_j at the points (spectral, k_d with the Nyquist mode set to 0)."""
    g = op.g
    kd = (g.kxd, g.kyd, g.kzd)
    ix = tuple(np.array(pts).T)
    G = np.zeros((len(pts), 3, 3))
    for i in range(3):
        uk = op.fft(u[i])
        for j in range(3):
            G[:, i, j] = op.ifft(1j*kd[j]*uk)[ix]
    return G


# --------------------------------------------------------------------------------------------- tail

def tail_stats(rho, s, nb_lm=None):
    n = rho.size
    hist = np.histogram(s, bins=PDF_EDGES)[0]
    tot = float(rho.sum())
    out = dict(rho_min=float(rho.min()), rho_max=float(rho.max()), mean_s=float(s.mean()), sigma_s=float(s.std()),
               pdf=(hist/(n*0.1)).tolist(), pdf_out=[int((s < PDF_EDGES[0]).sum()), int((s >= PDF_EDGES[-1]).sum())],
               F_V={}, F_M={}, n={})
    for th in THETA:
        sel = rho < th
        out['n'][str(th)] = int(sel.sum())
        out['F_V'][str(th)] = float(sel.sum()/n)
        out['F_M'][str(th)] = float(rho[sel].sum()/tot)
    lm = local_min(rho) if nb_lm is None else nb_lm
    v = rho[lm]
    out['nmin'] = {str(th): int((v < th).sum()) for th in NHOLE_THR}
    k = min(NLOW, n - 1)
    out['lowest'] = np.sort(np.partition(rho.ravel(), k)[:k]).tolist()
    return out


def tail(root, N):
    res = tail_case(os.path.join(root, CASES[N]), TAIL_T + (TAIL_T256 if N == 256 else ()))
    write_json(root, 'tail_N%d.json' % N, res)
    return res


def tail_case(case, times):
    """tail_stats of the out3d records of case nearest to times (tau0)."""
    g, op = svv3d.from_case(case)
    N = g.nlx
    cons = cons_table(case)
    res = dict(N=N, case=case, snaps=[])
    for tt in times:
        r = pick(case, tt)
        if r is None:
            print('N %d t/tau0 %.4f: no record' % (N, tt))
            continue
        sn = Snap(g, op, r[1], r[2], full=False)
        st = tail_stats(sn.rho, sn.s)
        st.update(t=r[0], t_tau=r[0]/TAU0, t_req=tt, path=r[1], rec=r[2], cons=cons_check(cons, r[0], st['rho_min']))
        res['snaps'].append(st)
        print('N %3d t/tau0 %.4f  rho %.4e .. %.3f  sigma_s %.4f  F_V(<0.1) %.3e  n(<0.05) %d  minima <0.2/0.1/0.05 '
              '%s  cons: dt %.1e rel %.1e'
              % (N, st['t_tau'], st['rho_min'], st['rho_max'], st['sigma_s'], st['F_V']['0.1'], st['n']['0.05'],
                 '/'.join(str(st['nmin'][str(x)]) for x in NHOLE_THR), st['cons']['dt'], st['cons']['rel']))
        sys.stdout.flush()
    return res


# --------------------------------------------------------------------------------------------- hole

def pick_holes(rho, lm, thr=HOLE_THR, nmax=NHOLE, sep=SEP):
    """The min of rho and up to nmax local minima below thr, deepest first, sep cells apart."""
    i0 = tuple(int(a) for a in np.unravel_index(int(np.argmin(rho)), rho.shape))
    cand = np.argwhere(lm & (rho < thr))
    cand = cand[np.argsort(rho[tuple(cand.T)], kind='stable')]
    out = [i0]
    for c in cand:
        c = tuple(int(a) for a in c)
        if all(pdist(c, q, rho.shape[0]) >= sep for q in out):
            out.append(c)
            if len(out) == nmax + 1:
                break
    return out


def hole_metrics(op, F, u, b, bands, p, Gp, budget=None, rho_bg=None):
    """The quantities of a hole at p (module docstring); F holds rho, s, Pth, PB, Ptot, Bor."""
    cs2 = op.cs2
    cs = np.sqrt(cs2)
    val = {k: float(F[k][p]) for k in ('rho', 's', 'Pth', 'PB', 'Ptot', 'Bor')}
    bv = np.array([float(x[p]) for x in b])
    uv = np.array([float(x[p]) for x in u])
    Bm = float(np.linalg.norm(bv))
    bh = bv/Bm
    D = float(np.trace(Gp))
    Dpar = float(bh @ Gp @ bh)
    sh = {k: float(np.median(shell_vals(F[k], p))) for k in ('rho', 's', 'Pth', 'PB', 'Ptot', 'Bor')}
    rbg = sh['rho'] if rho_bg is None else rho_bg
    out = dict(p=list(p), rho=val['rho'], B=Bm, Pth=val['Pth'], PB=val['PB'], Ptot=val['Ptot'],
               beta=2*val['Pth']/Bm**2, vA=Bm/np.sqrt(val['rho']), u=float(np.linalg.norm(uv)),
               upar=float(uv @ bh), bhat=bh.tolist(), D=D, Dpar=Dpar, Dperp=D - Dpar, bg=sh)
    out['Pi'] = (sh['Ptot'] - val['Ptot'])/(sh['Pth'] - val['Pth'])
    out['Lambda'] = val['Bor']/sh['Bor']
    out['depth'] = rbg/val['rho']
    # shape
    H = hess(F['s'], p)
    lam, vec = np.linalg.eigh(H)
    out['hess_eig'] = lam.tolist()
    out['w'] = [float(np.sqrt(2/x)) if x > 0 else None for x in lam[::-1]]  # thin first
    nh = vec[:, int(np.argmax(lam))]
    out['nhat'] = nh.tolist()
    out['angle_n_b'] = float(np.degrees(np.arccos(min(1.0, abs(float(nh @ bh))))))
    reg = region(F['rho'], p, np.sqrt(val['rho']*rbg))
    reg['angle_long_b'] = float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(reg['axes'][0], bh)))))))
    out['region'] = reg
    # cuts along the lattice directions nearest to b and n
    out['cut'] = {}
    for name, v in (('b', bh), ('n', nh)):
        d, ang = nearest_dir(v)
        c = {q: cut(F[q], p, d, CUT_L).tolist() for q in ('rho', 'Pth', 'PB', 'Ptot', 'Bor')}
        c['upar'] = sum(cut(u[k], p, d, CUT_L)*bh[k] for k in range(3)).tolist()
        c.update(dir=list(d), angle=ang, step=float(np.sqrt(np.dot(d, d))))
        out['cut'][name] = c
    cb = out['cut']['b']
    du = cb['upar'][CUT_L + CHI_L] - cb['upar'][CUT_L - CHI_L]
    rpar = 0.5*(cb['rho'][CUT_L + CHI_L] + cb['rho'][CUT_L - CHI_L])
    out['du_par'], out['rho_bg_par'] = du, rpar
    out['chi_par'] = float(np.log(rpar/val['rho'])/(du/(2*cs))) if du != 0 else None
    # band split of the depth of s
    if bands is not None:
        ds = val['s'] - sh['s']
        out['band'] = {k: float((v[p] - np.median(shell_vals(v, p)))/ds) for k, v in bands.items()}
    out['oddeven_run'] = oddeven(F['s'], p)
    if budget is not None:
        out['budget'] = budget
    out.update(classify(out))
    return out


def classify(h):
    ratio = h['depth']
    tags = []
    if abs(h['Pi']) <= 0.3:
        tags.append('balance')
    chi = h['chi_par']
    if (abs(h['Dpar']) >= abs(h['Dperp']) and 0.5*ratio <= h['Lambda'] <= 2*ratio
            and chi is not None and 0.5 <= chi <= 2):
        tags.append('drain')
    if abs(h['Dperp']) > abs(h['Dpar']) and 0.5 <= h['Lambda'] <= 2:
        tags.append('perp')
    oe = max(h['oddeven_run']) >= OE_RUN
    return dict(tags=tags, oddeven=oe, suspect=(not tags) or oe)


def fields(op, sn):
    B2 = sum(x**2 for x in sn.b)
    return dict(rho=sn.rho, s=sn.s, Pth=op.cs2*sn.rho, PB=0.5*B2, Ptot=op.cs2*sn.rho + 0.5*B2, Bor=np.sqrt(B2)/sn.rho)


def check_ref(ref, snap):
    """The comparison with ref (a value of CHECK, 5 significant digits), or None without ref."""
    if ref is None:
        return None
    h = snap['holes'][0]

    def same(x, r):
        return abs(x - r) <= 0.5*10.0**(np.floor(np.log10(abs(r))) - 4)
    ok = dict(t=abs(snap['t'] - ref['t']) <= 5e-7, imin=tuple(h['p']) == ref['imin'],
              adv=same(h['budget']['adv'], ref['adv']), svv=same(h['budget']['svv'], ref['svv']))
    return dict(ref=ref, got=dict(t=snap['t'], imin=h['p'], adv=h['budget']['adv'], svv=h['budget']['svv']),
                ok=ok, PASS=all(ok.values()))


def hole(root, N):
    res = hole_case(os.path.join(root, CASES[N]), HOLE_T[N], {t: v for (n, t), v in CHECK.items() if n == N})
    write_json(root, 'hole_N%d.json' % N, res)
    return res


def hole_case(case, times, refs=None):
    """The holes of the out3d records of case nearest to times (tau0); refs: {t_req: value of CHECK}."""
    g, op = svv3d.from_case(case)
    N = g.nlx
    refs = refs or {}
    cons = cons_table(case)
    dl = g.Lx/g.nlx
    bands_k = band_masks(g, num(read_nml(os.path.join(case, 'calliope.in')).get('svv_mn', '0.5')))
    res = dict(N=N, case=case, cs2=op.cs2, dl=dl, snaps=[])
    prev = None
    for tt in times:
        r = pick(case, tt)
        if r is None:
            print('N %d t/tau0 %.4f: no record' % (N, tt))
            continue
        sn = Snap(g, op, r[1], r[2])
        F = fields(op, sn)
        lm = local_min(sn.rho)
        pts = pick_holes(sn.rho, lm)
        G = grad_at(op, sn.u, pts)
        sk = op.fft(sn.s)
        bands = {k: op.ifft(sk*v) for k, v in bands_k.items()}
        eg = op.eps(sn.rho, sn.u, sn.b)
        divm = op.ifft(1j*(g.kxd*sn.mk[0] + g.kyd*sn.mk[1] + g.kzd*sn.mk[2])*op.mL)
        rdot, _, dsum, rL, rH = op.tend(sn.s, sn.u, sn.rho, eg, parts=True)
        snap = dict(t=r[0], t_tau=r[0]/TAU0, t_req=tt, path=r[1], rec=r[2], rho_min=float(sn.rho.min()),
                    cons=cons_check(cons, r[0], float(sn.rho.min())), eps_global=eg, lam_lin=op.lam_lin(eg),
                    n_lm={str(x): int((sn.rho[lm] < x).sum()) for x in NHOLE_THR}, holes=[])
        for i, p in enumerate(pts):
            bud = dict(adv=float(-divm[p]), svv=float(rdot[p]), svv_L=float(rL[p]), svv_H=float(rH[p]))
            h = hole_metrics(op, F, sn.u, sn.b, bands, p, G[i], bud)
            h['x'] = phys(p, dl)
            snap['holes'].append(h)
        # profile of rho along the thin axis of the min (+-PROF_L steps)
        h0 = snap['holes'][0]
        d = tuple(h0['cut']['n']['dir'])
        snap['profile'] = dict(dir=list(d), x0=h0['x'], step=float(np.sqrt(np.dot(d, d)))*dl,
                               rho=cut(sn.rho, tuple(h0['p']), d, PROF_L).tolist())
        if prev is not None:
            p0, p1 = tuple(prev['p']), tuple(h0['p'])
            dtt = r[0] - prev['t']
            snap['track'] = dict(dt=dtt, move=pdist(p0, p1, g.nlx)*dl, u_prev_dt=prev['u']*dtt, u_now_dt=h0['u']*dtt,
                                 rho_prev_at_now=float(prev['rho_arr'][p1]), rho_now_at_prev=float(sn.rho[p0]))
        prev = dict(p=h0['p'], t=r[0], u=h0['u'], rho_arr=sn.rho)
        snap['check'] = check_ref(refs.get(tt), snap)
        res['snaps'].append(snap)
        print_snap(N, snap)
        del F, G, bands, divm, rdot, rL, rH, sk
    return res


def print_snap(N, sn):
    print('N %3d t/tau0 %.4f  rho_min %.4e (cons: dt %.1e rel %.1e)  minima <0.2/0.1/0.05 %s  lam_lin %.4e'
          % (N, sn['t_tau'], sn['rho_min'], sn['cons']['dt'], sn['cons']['rel'],
             '/'.join(str(sn['n_lm'][str(x)]) for x in NHOLE_THR), sn['lam_lin']))
    for h in sn['holes']:
        print('   %-15s rho %.3e depth %6.1f beta %.2e vA %5.2f |u| %5.2f  D %+8.2f D_par %+8.2f  Pi %+7.3f '
              'Lambda %7.2f chi %s  w %s  cells %d semi %s  ang(long,b) %4.1f  band L/M/H %s  oe %s  '
              'adv %+.4e svv %+.4e  %s%s'
              % (tuple(h['p']), h['rho'], h['depth'], h['beta'], h['vA'], h['u'], h['D'], h['Dpar'], h['Pi'],
                 h['Lambda'], '%.2f' % h['chi_par'] if h['chi_par'] is not None else 'n/a',
                 '/'.join('%.2f' % w if w is not None else 'inf' for w in h['w']), h['region']['cells'],
                 '/'.join('%.1f' % x for x in h['region']['semi']), h['region']['angle_long_b'],
                 '/'.join('%.2f' % h['band'][k] for k in 'LMH'), h['oddeven_run'], h['budget']['adv'],
                 h['budget']['svv'], ','.join(h['tags']) or '-', ' SUSPECT' if h['suspect'] else ''))
    if 'track' in sn:
        tr = sn['track']
        print('   track: move %.4f  |u|dt prev %.4f now %.4f  rho_prev(x_now) %.3e  rho_now(x_prev) %.3e'
              % (tr['move'], tr['u_prev_dt'], tr['u_now_dt'], tr['rho_prev_at_now'], tr['rho_now_at_prev']))
    if sn['check'] is not None:
        print('   check against the restart copy: %s %s'
              % ('PASS' if sn['check']['PASS'] else 'FAIL', sn['check']['got']))
    sys.stdout.flush()


# -------------------------------------------------------------------------------------------- judge

def tail_verdict(co, fi):
    """The tail criteria for one pair of snapshots (coarse, fine) at the same time."""
    use = [th for th in THETA if co['n'][str(th)] >= TAIL_NMIN]
    rV = {str(th): fi['F_V'][str(th)]/co['F_V'][str(th)] for th in use}
    rM = {str(th): fi['F_M'][str(th)]/co['F_M'][str(th)] for th in use}
    nbelow = int(np.searchsorted(fi['lowest'], co['rho_min']))
    nbelow = nbelow if nbelow < len(fi['lowest']) else None  # None: more than NLOW
    conv = bool(use) and all(abs(x - 1) <= 0.25 for x in rV.values()) and nbelow is not None and nbelow <= 100
    fat = sum(1 for th in use if th <= 0.2 and rV[str(th)] >= 1.25) >= 2
    return dict(theta_used=use, ratio_V=rV, ratio_M=rM, n_fine_below_coarse_min=nbelow,
                verdict='fattening' if fat else 'converged' if conv else 'intermediate')


def cross(root, holes, tails):
    """For t = 0.5, 1 tau0: each hole of 256^3 against 64^3 and 128^3 (min of rho within 2 coarse cells of the
    same position, the nearest analysed coarse hole) and the profiles along the thin axis of the 256^3 min."""
    out = {}
    for tt in (0.5, 1.0):
        fine = next((s for s in holes[256]['snaps'] if s['t_req'] == tt), None)
        if fine is None:
            continue
        ent = dict(t_tau=fine['t_tau'], holes=[], profiles={256: fine['profile']})
        coarse = {}
        for N in (64, 128):
            case = os.path.join(root, CASES[N])
            g, op = svv3d.from_case(case)
            r = pick(case, tt)
            coarse[N] = (g, Snap(g, op, r[1], r[2], full=False).rho, next(s for s in holes[N]['snaps']
                                                                          if s['t_req'] == tt))
        for h in fine['holes']:
            e = dict(p=h['p'], x=h['x'], rho=h['rho'], tags=h['tags'], suspect=h['suspect'])
            for N, (g, rho, cs) in coarse.items():
                dl = g.Lx/g.nlx
                q = tuple(int(round(h['x'][a]/dl)) % g.nlx for a in (0, 2, 1))  # array index (x, z, y)
                r2 = np.arange(-2, 3)
                sub = rho[np.ix_((q[0] + r2) % g.nlx, (q[1] + r2) % g.nlz, (q[2] + r2) % g.nly)]
                dists = []
                for ch in cs['holes']:
                    d = np.abs(np.array(ch['x']) - np.array(h['x']))
                    d = np.minimum(d, g.Lx - d)
                    dists.append((float(np.sqrt((d**2).sum())), ch['rho']))
                near = min(dists)
                e['N%d' % N] = dict(rho_min_2cells=float(sub.min()), nearest_hole_dist=near[0]/(2*np.pi/256),
                                    nearest_hole_rho=near[1])
            ent['holes'].append(e)
        pf = fine['profile']
        for N, (g, rho, _) in coarse.items():
            dl = g.Lx/g.nlx
            d = tuple(pf['dir'])
            q = tuple(int(round(pf['x0'][a]/dl)) % g.nlx for a in (0, 2, 1))
            L = max(1, int(round(PROF_L*N/256)))  # the same length as +-PROF_L steps of 256^3
            ent['profiles'][N] = dict(dir=list(d), step=float(np.sqrt(np.dot(d, d)))*dl,
                                      x0=[q[0]*dl, q[2]*dl, q[1]*dl], rho=cut(rho, q, d, L).tolist())
        out[str(tt)] = ent
    return out


def judge(root):
    tails = {N: read_json(root, 'tail_N%d.json' % N) for N in CASES}
    holes = {N: read_json(root, 'hole_N%d.json' % N) for N in CASES}
    summ = dict(tail={}, hole={}, checks={})
    # tail criteria
    for co, fi in ((64, 128), (128, 256)):
        key = '%d-%d' % (co, fi)
        summ['tail'][key] = {}
        for tt in TAIL_T:
            a = next((s for s in tails[co]['snaps'] if s['t_req'] == tt), None)
            b = next((s for s in tails[fi]['snaps'] if s['t_req'] == tt), None)
            if a and b:
                summ['tail'][key][str(tt)] = tail_verdict(a, b)
        vs = [v['verdict'] for v in summ['tail'][key].values()]
        summ['tail'][key]['overall'] = ('fattening' if 'fattening' in vs else
                                        'converged' if vs and all(v == 'converged' for v in vs) else 'intermediate')
    # holes
    for N in CASES:
        hs = [(s['t_tau'], h) for s in holes[N]['snaps'] for h in s['holes']]
        summ['hole'][N] = dict(n=len(hs), suspect=sum(h['suspect'] for _, h in hs),
                               tags={k: sum(k in h['tags'] for _, h in hs) for k in ('balance', 'drain', 'perp')},
                               oddeven=sum(h['oddeven'] for _, h in hs))
    many = summ['hole'][256]['suspect'] >= 0.5*summ['hole'][256]['n']
    ov = summ['tail']['128-256']['overall']
    if ov == 'fattening' or many:
        nxt = 'step 3 (svv_eps0, C_shock) and a revision of the minimum principle of ES-SVV'
    elif ov == 'converged':
        nxt = 'step 4 (DNS reference, needs approval) or 256^3 as the upper limit'
    else:
        nxt = 'step 3 (the dependence on the coefficients)'
    summ['verdict'] = dict(tail_128_256=ov, suspect_many_256=bool(many), next=nxt)
    # reading checks
    for N in CASES:
        for kind, d in (('tail', tails[N]), ('hole', holes[N])):
            for s in d['snaps']:
                summ['checks'].setdefault(str(N), []).append(dict(kind=kind, t_tau=s['t_tau'], cons=s['cons']))
    summ['checks']['restart'] = [s['check'] for s in holes[256]['snaps'] if s['check'] is not None]
    summ['cross'] = cross(root, holes, tails)
    write_json(root, 'summary.json', summ)
    figures(root, tails, holes, summ)
    print_judge(summ)
    return summ


def print_judge(summ):
    for key, d in summ['tail'].items():
        for tt, v in d.items():
            if tt == 'overall':
                continue
            print('tail %s t/tau0 %s: %-12s theta %s  F_V ratio %s  n_fine(rho < rho_min,coarse) %s'
                  % (key, tt, v['verdict'], v['theta_used'], ' '.join('%.2f' % x for x in v['ratio_V'].values()),
                     v['n_fine_below_coarse_min']))
        print('tail %s overall: %s' % (key, d['overall']))
    for N, h in summ['hole'].items():
        print('holes N %s: %d analysed, suspect %d, tags %s, odd-even %d' % (N, h['n'], h['suspect'], h['tags'],
                                                                             h['oddeven']))
    for c in summ['checks']['restart']:
        print('restart check: %s  got %s  ref %s' % ('PASS' if c['PASS'] else 'FAIL', c['got'], c['ref']))
    bad = [(N, c) for N, cs in summ['checks'].items() if N != 'restart' for c in cs
           if abs(c['cons']['dt']) < 1e-9 and abs(c['cons']['rel']) > 1e-9]
    print('rho_min against conservation.dat at the same t: %s' % ('all within 1e-9' if not bad else bad))
    for tt, e in summ['cross'].items():
        for h in e['holes']:
            print('cross t/tau0 %s: 256^3 hole %s rho %.3e %s | 128^3: min(2 cells) %.3e nearest hole %.1f fine '
                  'cells (rho %.3e) | 64^3: min(2 cells) %.3e nearest hole %.1f (rho %.3e)'
                  % (tt, tuple(h['p']), h['rho'], ','.join(h['tags']) or '-', h['N128']['rho_min_2cells'],
                     h['N128']['nearest_hole_dist'], h['N128']['nearest_hole_rho'], h['N64']['rho_min_2cells'],
                     h['N64']['nearest_hole_dist'], h['N64']['nearest_hole_rho']))
    print('verdict: %s' % summ['verdict'])


def figures(root, tails, holes, summ):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fd = os.path.join(root, 'stage7', 'fig')
    os.makedirs(fd, exist_ok=True)
    sc = 0.5*(PDF_EDGES[1:] + PDF_EDGES[:-1])
    # 1. PDF of s and F_V(rho < theta) at 0.5 and 1 tau0
    fig, ax = plt.subplots(2, 2, figsize=(9.2, 6.4))
    for i, tt in enumerate((0.5, 1.0)):
        for N in CASES:
            s = next((x for x in tails[N]['snaps'] if x['t_req'] == tt), None)
            if s is None:
                continue
            pdf = np.array(s['pdf'])
            ax[i, 0].semilogy(sc[pdf > 0], pdf[pdf > 0], color=COLOR_N[N], lw=1.5, label='N = %d' % N)
            ax[i, 1].loglog(THETA, [max(s['F_V'][str(th)], 1e-12) for th in THETA], color=COLOR_N[N], lw=1.5,
                            marker='o', ms=4, label='N = %d' % N)
        ax[i, 0].set_title('volume PDF of s = ln rho, t = %.1f tau0' % tt, fontsize=9)
        ax[i, 0].set_xlabel('s')
        ax[i, 0].set_ylabel('PDF')
        ax[i, 0].set_xlim(-6, 2)
        ax[i, 1].set_title('F_V(rho < theta), t = %.1f tau0' % tt, fontsize=9)
        ax[i, 1].set_xlabel('theta')
        ax[i, 1].set_ylabel('volume fraction')
        ax[i, 1].set_ylim(1e-8, 1)
        for a in ax[i]:
            a.legend(frameon=False, fontsize=8)
            a.grid(alpha=0.3, lw=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(fd, 'stage7_tail.png'), dpi=110)
    plt.close(fig)
    # 2. rho_min(t) of conservation.dat
    fig, a = plt.subplots(figsize=(6.4, 3.6))
    for N in CASES:
        c = cons_table(os.path.join(root, CASES[N]))
        a.semilogy(c[:, 0]/TAU0, c[:, 11], color=COLOR_N[N], lw=1.2, label='N = %d' % N)
    a.set_xlabel('t/tau0')
    a.set_ylabel('rho_min')
    a.set_title('turb1 (M_s = 1): min of rho over the grid', fontsize=9)
    a.legend(frameon=False, fontsize=8)
    a.grid(alpha=0.3, lw=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(fd, 'stage7_rhomin.png'), dpi=110)
    plt.close(fig)
    # 3. rho along the thin axis of the min of 256^3, three N
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4), sharey=True)
    for a, (tt, e) in zip(ax, sorted(summ['cross'].items())):
        for N in CASES:
            pf = e['profiles'][N] if N in e['profiles'] else e['profiles'][str(N)]
            r = np.array(pf['rho'])
            L = (len(r) - 1)//2
            a.semilogy(np.arange(-L, L + 1)*pf['step'], r, color=COLOR_N[N], lw=1.5, marker='o', ms=3,
                       label='N = %d' % N)
        a.set_title('t = %s tau0, along %s through the min of 256^3'
                    % (tt, tuple(e['profiles'][256 if 256 in e['profiles'] else '256']['dir'])), fontsize=9)
        a.set_xlabel('distance from the min of 256^3')
        a.legend(frameon=False, fontsize=8)
        a.grid(alpha=0.3, lw=0.5)
    ax[0].set_ylabel('rho')
    fig.tight_layout()
    fig.savefig(os.path.join(fd, 'stage7_hole.png'), dpi=110)
    plt.close(fig)
    # 4. pressures around the min of 256^3 at 0.893 tau0, along the cuts near n and b
    sn = next((s for s in holes[256]['snaps'] if s['t_req'] == 0.893), None)
    if sn is not None:
        h = sn['holes'][0]
        fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4), sharey=True)
        for a, key in zip(ax, ('n', 'b')):
            c = h['cut'][key]
            ell = np.arange(-CUT_L, CUT_L + 1)*c['step']
            # red-green CVD distance 7.2: told apart by line style and marker as well
            for q, col, ls, mk, lab in (('Pth', '#6250d6', '-', 'o', 'P_th'), ('PB', '#e34948', '--', 's', 'P_B'),
                                        ('Ptot', '#008300', '-', '^', 'P_tot')):
                v = np.array(c[q])
                a.plot(ell, v - v[CUT_L], color=col, ls=ls, lw=1.5, marker=mk, ms=4, label='%s - %s(x0)' % (lab, lab))
            a.set_title('256^3, t = %.3f tau0, min of rho: along %s (%.0f deg from %s)'
                        % (sn['t_tau'], tuple(c['dir']), c['angle'], 'n' if key == 'n' else 'b'), fontsize=8)
            a.set_xlabel('distance (cells)')
            a.legend(frameon=False, fontsize=8)
            a.grid(alpha=0.3, lw=0.5)
        ax[0].set_ylabel('pressure - value at the min')
        fig.tight_layout()
        fig.savefig(os.path.join(fd, 'stage7_balance.png'), dpi=110)
        plt.close(fig)


# ---------------------------------------------------------------------------------------------- I/O

def write_json(root, name, obj):
    d = os.path.join(root, 'stage7')
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), 'w') as f:
        json.dump(obj, f, indent=1, default=lambda x: x.tolist() if hasattr(x, 'tolist') else str(x))


def read_json(root, name):
    with open(os.path.join(root, 'stage7', name)) as f:
        d = json.load(f)
    return d


# -------------------------------------------------------------------------------------------- tests

def selftest():
    ok = True

    def rep(name, good, msg):
        nonlocal ok
        ok &= bool(good)
        print('%-34s %s  %s' % (name, msg, 'PASS' if good else 'FAIL'))

    # 1. rotated anisotropic Gaussian hole s = -2 exp(-r.M.r/2): Hessian widths = sigma_i; region semi-axes
    N = 48
    g = gen_ic.Grid(N, N, N)
    dl = g.Lx/N
    p = (24, 24, 24)
    r = [g.X/dl - 24, g.Y/dl - 24, g.Z/dl - 24]
    a1, a2, a3 = 0.4, 0.7, 1.1
    Rz = np.array([[np.cos(a1), -np.sin(a1), 0], [np.sin(a1), np.cos(a1), 0], [0, 0, 1]])
    Ry = np.array([[np.cos(a2), 0, np.sin(a2)], [0, 1, 0], [-np.sin(a2), 0, np.cos(a2)]])
    Rx = np.array([[1, 0, 0], [0, np.cos(a3), -np.sin(a3)], [0, np.sin(a3), np.cos(a3)]])
    R = Rz @ Ry @ Rx
    sig = np.array([3.0, 5.0, 8.0])
    M = R @ np.diag(1/sig**2) @ R.T
    q = sum(M[i, j]*r[i]*r[j] for i in range(3) for j in range(3))
    s = -2*np.exp(-q/2)
    lam, vec = np.linalg.eigh(hess(s, p))
    w = np.sqrt(2/lam[::-1])
    nh = vec[:, int(np.argmax(lam))]
    ang = np.degrees(np.arccos(abs(nh @ R[:, 0])))
    rep('Hessian widths (3, 5, 8 cells)', np.all(np.abs(w/sig - 1) < 0.05) and ang < 3,
        'w %s  angle(thin, sigma=3 axis) %.2f deg' % (np.round(w, 3), ang))
    reg = region(np.exp(s), p, np.exp(-1.0))
    ex = sig[::-1]*np.sqrt(2*np.log(2))
    ang = np.degrees(np.arccos(abs(np.dot(reg['axes'][0], R[:, 2]))))
    rep('region semi-axes', np.all(np.abs(np.array(reg['semi'])/ex - 1) < 0.05) and ang < 3 and not reg['edge'],
        'semi %s (expected %s)  angle(long) %.2f deg  cells %d' % (np.round(reg['semi'], 2), np.round(ex, 2), ang,
                                                                   reg['cells']))
    # 2. Pi and Lambda: rho = 1 - 0.5 exp(-r^2/8) at rest; B_z balancing P_tot, flat, or prop. to rho
    N = 32
    g = gen_ic.Grid(N, N, N)
    op = svv3d.SVV(g, 1.0)
    dl = g.Lx/N
    p = (16, 16, 16)
    r2 = ((g.X/dl - 16)**2 + (g.Y/dl - 16)**2 + (g.Z/dl - 16)**2)
    rho = 1 - 0.5*np.exp(-r2/8)
    zero = np.zeros(op.shape)
    u = [zero, zero, zero]
    G0 = np.zeros((3, 3))
    for name, bz, want in (('Pi, P_tot balanced', np.sqrt(1 + 2*(1 - rho)), ('Pi', 0.0)),
                           ('Pi, flat B', np.ones(op.shape), ('Pi', 1.0)),
                           ('Lambda, B prop. to rho', rho, ('Lambda', 1.0))):
        b = [zero, zero, bz]
        sn = type('S', (), dict(rho=rho, s=np.log(rho), b=b))
        h = hole_metrics(op, fields(op, sn), u, b, None, p, G0)
        rep(name, abs(h[want[0]] - want[1]) < 1e-10, '%s = %.3e' % (want[0], h[want[0]]))
    # 3. chi_par of a symmetric isothermal rarefaction along b = z: u_z = U tanh(z/w), rho = exp(-U(1 - |tanh|)/c_s)
    U = 0.8
    th = np.tanh(g.Z/dl - 16)
    rho = np.exp(-U*(1 - np.abs(th)))
    b = [zero, zero, np.ones(op.shape)]
    sn = type('S', (), dict(rho=rho, s=np.log(rho), b=b))
    h = hole_metrics(op, fields(op, sn), [zero, zero, U*th], b, None, p, G0)
    rep('chi_par of a rarefaction', abs(h['chi_par'] - 1) < 0.01 and h['cut']['b']['dir'] == [0, 0, 1],
        'chi_par = %.4f' % h['chi_par'])
    # 4. D_par and D_perp: u = (0.3 sin x, 0, 0.2 sin z), b = z
    u = [0.3*np.sin(g.X), zero, 0.2*np.sin(g.Z)]
    pts = [(5, 9, 3)]
    Gp = grad_at(op, u, pts)[0]
    bh = np.array([0, 0, 1.0])
    Dpar = bh @ Gp @ bh
    Dperp = np.trace(Gp) - Dpar
    x, z = g.x[5], g.z[9]
    rep('D_par, D_perp', abs(Dpar - 0.2*np.cos(z)) < 1e-12 and abs(Dperp - 0.3*np.cos(x)) < 1e-12,
        'D_par %.6f (%.6f)  D_perp %.6f (%.6f)' % (Dpar, 0.2*np.cos(z), Dperp, 0.3*np.cos(x)))
    # 5. local minima across the periodic boundary
    f = np.ones((16, 16, 16))
    f[0, 5, 5], f[15, 15, 15], f[8, 8, 8] = 0.5, 0.4, 0.6
    f[0, 0, 0], f[15, 0, 0] = 0.7, 0.65  # both next to (15, 15, 15) across the boundary
    lm = local_min(f)
    got = sorted(tuple(int(a) for a in x) for x in np.argwhere(lm & (f < 0.9)))
    rep('local minima (periodic)', got == [(0, 5, 5), (8, 8, 8), (15, 15, 15)], '%s' % got)
    hs = pick_holes(f, lm, thr=0.9, sep=9)  # (0, 5, 5) is 8.5 cells from (15, 15, 15) across the boundary
    rep('pick_holes (sep 9, periodic)', hs == [(15, 15, 15), (8, 8, 8)], '%s' % hs)
    # 6. lattice directions and cuts: f[i, k, j] = x + 10 y + 100 z (cells)
    tests = [((1, 1, 0.1), (1, 1, 0)), ((0.1, 0, -1), (0, 0, -1)), ((1, -1, 1.1), (1, -1, 1)),
             ((-0.2, 0.9, 1), (0, 1, 1))]
    good = all(nearest_dir(np.array(v)/np.linalg.norm(v))[0] == d for v, d in tests)
    i, k, j = np.indices((32, 32, 32))
    f = (i + 10*j + 100*k).astype(float)
    steps = [float(np.diff(cut(f, (10, 10, 10), d, 2)).mean()) for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, -1, 1))]
    good &= steps == [1.0, 10.0, 100.0, 91.0]
    rep('13 directions, cuts', good, 'steps %s' % steps)
    # 7. odd-even: a Gaussian dip with and without a grid-scale oscillation along x
    s = -np.exp(-r2/8)
    sx = s + 0.05*np.where((g.X/dl).round().astype(int) % 2 == 0, 1.0, -1.0)*np.exp(-r2/32)
    a, b_ = oddeven(s, p), oddeven(sx, p)
    rep('odd-even check', max(a) < OE_RUN and b_[0] >= OE_RUN, 'smooth %s, oscillating %s' % (a, b_))
    # 8. band masks: the parts add up, L inside the mask
    bm = band_masks(g, 0.5)
    good = np.all(bm['L'] + bm['M'] + bm['H'] == 1) and np.all(bm['M'] >= 0)
    rep('band masks', good, 'modes L %d M %d H %d' % tuple(int(bm[x].sum()) for x in 'LMH'))
    print('all PASS' if ok else 'some FAIL')
    return ok


if __name__ == '__main__':
    a = sys.argv[1:]
    if not a or a[0] == 'selftest':
        sys.exit(0 if selftest() else 1)
    elif a[0] == 'tail':
        tail(a[1], int(a[2]))
    elif a[0] == 'hole':
        hole(a[1], int(a[2]))
    elif a[0] == 'judge':
        judge(a[1])
    else:
        sys.exit(__doc__)
