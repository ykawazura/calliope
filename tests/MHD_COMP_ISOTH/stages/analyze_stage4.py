#!/usr/bin/env python3
"""Analysis of the Stage 4 runs (gen_stage4.py) of the Fortran code.

    python3 analyze_stage4.py <batch> [base]   (default base: runs/stage4) prints the checks of the batch and
                                               writes <base>/<batch>/summary.json

Batches
    readback  restart/ of the nstep = 0 runs equals ic/ byte for byte; the out2d planes at t = 0 equal rho,
              u, B and D = div_N u of gen_ic.py (relative error <= TOL_READ); planes, out1d and
              conservation.dat do not depend on the number of ranks (byte for byte); div B columns <= TOL_DIVB
    regress   the outputs of the Stage 3 and Stage 4 binaries agree bit for bit (conservation.dat columns
              1-18, the other files whole)
    alfven    CP Alfven waves of S4.1 (see alfven)
    cpdiag    growth of the A = 1 CP waves with the dissipation terms one at a time (see cpdiag)
    planar    planar problems of S4.2 along an oblique normal: k = (1, 0) against 1D (TOL_CTL), k = (1, 1)
              without SVV against the compressed 1D run (TOL_EQV), the translation symmetry f(i, l) = F(p) on
              the planes (TOL_SYM), and the 2D runs with the default settings against the compressed 1D runs
              and the 1D runs of Stage 3 (runs/stage3_f next to base)
"""
import sys
import os
import re
import json
import glob
import math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths
import gen_ic
from analyze_stage1 import load, load_runs, read_lines, run_dirs
import stage3 as s3
import cp_stability as cp

TOL_READ = 1e-13
TOL_DIVB = 1e-13
TOL_RED = 1e-13


def reduction_diff(a, b, name):
    """max over the columns of |a - b|/scale for two outputs of reductions; the scale is the larger of the
    column's own max|a| and the scale of the file (conservation.dat: E; budget.dat: E/dt; cfl.dat: none),
    so that the columns at the rounding level of a larger total (M - M0, R_X, err_t) are measured
    against that total."""
    x, y = load(a), load(b)
    if x is None or y is None or x.shape != y.shape:
        return float('inf')
    if name == 'conservation.dat':
        ref = np.abs(x[:, 10]).max()
    elif name == 'budget.dat':
        ref = (np.abs(x[:, 18])/x[:, 1]).max()
    else:
        ref = 0.0
    scale = np.maximum(np.abs(x).max(axis=0), ref)
    scale[scale == 0] = 1.0
    return float((np.abs(x - y).max(axis=0)/scale).max())


def case_info(case):
    return json.load(open(os.path.join(case, 'case.json')))


def grid_of(case):
    from analyze_stage1 import read_nml, num
    d = read_nml(os.path.join(case, 'calliope.in'))
    return gen_ic.Grid(int(d['nlx']), int(d['nly']), int(d['nlz']),
                       num(d.get('lx', '1.0')), num(d.get('ly', '1.0')), num(d.get('lz', '1.0')))


def nml_failed(case):
    """Namelist groups of calliope.in whose reading failed according to the .std of the runs (the code
    prints 'Reading <group> failed' and keeps the defaults)."""
    groups = set(re.findall(r'^&(\w+)', open(os.path.join(case, 'calliope.in')).read(), re.M))
    bad = set()
    for path in run_dirs(case):
        for f in glob.glob(os.path.join(path, '*.std')):
            bad |= set(re.findall(r'Reading (\w+) failed', open(f, errors='replace').read()))
    return sorted(bad & groups)


def status(case):
    try:
        s = open(os.path.join(case, 'exit_status')).read().strip()
    except OSError:
        return 'not run'
    bad = nml_failed(case)
    return '%s, namelist read failed: %s' % (s, ' '.join(bad)) if bad else s


def superseded(case):
    """Name of the case that replaces this one (file superseded_by: a run stopped or killed for a reason recorded
    in the report and run again under a new name), or None. Such a case is listed but not judged."""
    try:
        return open(os.path.join(case, 'superseded_by')).read().strip()
    except OSError:
        return None


def same_bytes(a, b):
    try:
        return open(a, 'rb').read() == open(b, 'rb').read()
    except OSError:
        return False


def readback(root):
    out = {}
    ok_all = True
    for case in sorted(glob.glob(os.path.join(root, '*_n[01]'))):
        name = os.path.basename(case)
        g = grid_of(case)
        r = dict(status=status(case))
        ic = os.path.join(case, 'ic')
        if name.endswith('_n0'):
            r['restart_eq_ic'] = all(same_bytes(os.path.join(ic, n + '.dat'), os.path.join(case, 'restart', n + '.dat'))
                                     for n in gen_ic.NAMES)
            r['ok'] = r['status'] == 'exit 0' and r['restart_eq_ic']
        else:
            lk, mk, bk, _ = gen_ic.read_restart(ic, g)
            rho, u, b = gen_ic.real_from_state(g, lk, mk, bk)
            D = g.div_d(tuple(g.fft(ui)*g.mask for ui in u))
            ref = dict(rho=rho, ux=u[0], uy=u[1], uz=u[2], bx=b[0], by=b[1], bz=b[2], D=D)
            planes = [('y0', (slice(None), slice(None), 0))]
            if g.nly > 2:
                planes.append(('z0', (slice(None), 0, slice(None))))
            err = {}
            for pl, sl in planes:
                p = gen_ic.read_planes(os.path.join(case, 'out2d'), g, pl)
                r['nrec_' + pl] = len(p['t'])
                for q, f in ref.items():
                    scale = max(abs(f).max(), abs(ref['ux']).max())
                    err['%s_%s' % (pl, q)] = float(abs(p[q][0] - f[sl]).max()/scale)
            r['err_read_max'] = max(err.values())
            r['err_read'] = err
            c = load(os.path.join(case, 'conservation.dat'))
            r['divb_rms_max'] = float(c[:, 18].max())
            r['divb_max_max'] = float(c[:, 19].max())
            r['ok'] = (r['status'] == 'exit 0' and r['err_read_max'] <= TOL_READ and r['nrec_y0'] == 2
                       and max(r['divb_rms_max'], r['divb_max_max']) <= TOL_DIVB)
        out[name] = r
        ok_all &= r['ok']
    # the same input with more ranks: the state (restart/ after the step, planes, out1d) byte-identical; the
    # reductions (sums over the ranks by MPI_Allreduce, whose order depends on the ranks) within TOL_RED of
    # the scale of their file
    for case in sorted(glob.glob(os.path.join(root, '*_np[2-9]*_n1'))):
        ref = case.replace('_np%s_' % case.rsplit('_np', 1)[1].split('_')[0], '_np1_')
        files = (['restart/%s.dat' % n for n in gen_ic.NAMES]
                 + ['out2d/plane_y0.dat', 'out2d/plane_z0.dat', 'out1d/line_00000.dat', 'out1d/line_00001.dat'])
        files = [f for f in files if os.path.exists(os.path.join(ref, f))]
        d = {f: same_bytes(os.path.join(ref, f), os.path.join(case, f)) for f in files}
        red = {f: reduction_diff(os.path.join(ref, f), os.path.join(case, f), f)
               for f in ('conservation.dat', 'cfl.dat', 'budget.dat')}
        r = out[os.path.basename(case)]
        r['same_as_np1'] = d
        r['reduction_diff'] = red
        r['ok'] &= all(d.values()) and max(red.values()) <= TOL_RED
        ok_all &= r['ok']
    for name, r in out.items():
        line = '%-18s %-8s ok=%s' % (name, r['status'], r['ok'])
        if 'restart_eq_ic' in r:
            line += '  restart==ic %s' % r['restart_eq_ic']
        if 'err_read_max' in r:
            worst = max(r['err_read'], key=r['err_read'].get)
            line += '  read err %.2e (%s)  divB rms %.1e max %.1e' % (r['err_read_max'], worst, r['divb_rms_max'],
                                                                       r['divb_max_max'])
        if 'same_as_np1' in r:
            line += '  state same as np1: %s' % ('all' if all(r['same_as_np1'].values()) else
                                                 [f for f, v in r['same_as_np1'].items() if not v])
            line += '  reductions %.1e' % max(r['reduction_diff'].values())
        print(line)
    print('readback:', 'PASS' if ok_all else 'FAIL')
    return dict(cases=out, ok=ok_all)


def regress_pair(old, new):
    """Statuses and the files of the run old compared with those of new: columns 1-18 of conservation.dat,
    cfl.dat, budget.dat, shock_selftest.dat and out1d/* byte for byte."""
    r = dict(status=(status(old), status(new)))
    d = {}
    a, b = load(os.path.join(old, 'conservation.dat')), load(os.path.join(new, 'conservation.dat'))
    if a is not None and b is not None:
        d['conservation.dat[1:18]'] = a.shape[0] == b.shape[0] and np.array_equal(a[:, :18], b[:, :18])
        d['conservation.dat ncol'] = (a.shape[1], b.shape[1])
    for f in ['cfl.dat', 'budget.dat', 'shock_selftest.dat'] + sorted(os.path.relpath(p, old) for p in glob.glob(os.path.join(old, 'out1d', '*'))):
        if os.path.exists(os.path.join(old, f)):
            d[f] = same_bytes(os.path.join(old, f), os.path.join(new, f))
    r['same'] = d
    r['ok'] = r['status'][0] == r['status'][1] and all(v for k, v in d.items() if k != 'conservation.dat ncol')
    return r


def regress(root):
    out = {}
    ok_all = True
    for old in sorted(glob.glob(os.path.join(root, '*_old'))):
        name = os.path.basename(old)[:-4]
        r = regress_pair(old, old[:-4] + '_new')
        d = r['same']
        out[name] = r
        ok_all &= r['ok']
        bad = [k for k, v in d.items() if v is False]
        print('%-28s %-16s ok=%s  files %d  differ: %s' % (name, '/'.join(r['status']), r['ok'], len(d) - 1, bad or '-'))
    print('regress:', 'PASS' if ok_all else 'FAIL')
    return dict(cases=out, ok=ok_all)


TOL_INV = 1e-12    # |rho - 1|, |u.n|/A
TOL_NUM = 1e-12    # zeta, eta <= TOL_NUM A dl (dl = L/N); default - ideal, N, axis permutation
TOL_DISS = 1e-14   # energy removed by the shock viscosity, SVV and eta over a period / E, above the eta floor:
                   # the round-off eta (<= TOL_NUM A dl) on the current of the wave |J0| = |k0| A removes up to
                   # eta_max |k0|^2 A^2 V T/E (2e-14 - 7e-14 at A = 1)
TOL_DEF = 1e-10    # |default - ideal|/A
TOL_MPI = 1e-13
DT3 = (7.0, 9.0)   # error(dt)/error(dt/2) of the RK3
TOL_PRED = 1e-11  # max|(u, B) - RK3 solution of the eigenmode|/A
# A = 1: the wave is linearly unstable (cp_stability.py), so that the round-off grows; the field checks use
# t = T only and the growth rate gamma of the perturbation (offpump_growth) is checked instead
GROWTH_AMAX = 1e-4   # upper end of the fit window of the perturbation amplitude (linear regime)
GROWTH_TOL = 0.25    # ideal: |gamma/gamma_galerkin - 1| for the dominant family
GROWTH_PHYS = 1.25   # default, N >= GROWTH_NMIN: gamma <= GROWTH_PHYS x the largest physical rate on the grid
GROWTH_NMIN = 32     # below: the default runs grow faster than this bound, recorded, not judged. The shock
                     # resistivity drives it (cpdiag): curl(eta J0) with eta ~ S[max(0, -dD)] is of first order
                     # in the perturbation, and the runs without it follow the Galerkin rate of their family


def snapshots(case, g, T):
    """out3d snapshots at t = m T: dict m -> (t, rho, u, B) in real space (the first record of each m)."""
    o = gen_ic.read_out3d(os.path.join(case, 'out3d'), g)
    out = {}
    for r, t in enumerate(o['t']):
        m = int(round(t/T))
        if m in out or abs(t - m*T) > 1e-9*T:
            continue
        rho = g.ifft(o['rho'][r])
        out[m] = (t, rho, tuple(g.ifft(o[n][r])/rho for n in ('mx', 'my', 'mz')),
                  tuple(g.ifft(o[n][r]) for n in ('bx', 'by', 'bz')))
    return out


def roll_fields(v, g, s):
    """Vector field v (components on [i, k, j] = [x, z, y] arrays) of a cubic grid rotated by the cyclic
    permutation x -> y -> z of the axes, s times: v'_{d+1}(roll(x)) = v_d(x)."""
    for _ in range(s):
        std = [c.transpose(0, 2, 1) for c in v]                       # [x, y, z]
        rot = [np.transpose(std[(d - 1) % 3], (2, 0, 1)) for d in range(3)]
        v = [c.transpose(0, 2, 1) for c in rot]
    return v


def offpump_growth(case, g, info):
    """Growth of the perturbation of the CP wave from the out3d records (one per period): the dominant k (the
    Fourier mode k != 0, +-k0 of largest |(m, B)_k|/a at the last record with max <= GROWTH_AMAX) and A(t), the
    amplitude of its family (the members k_dom + n k0 and their conjugates, root of the sum of |(m, B)_k|^2/a^2);
    gamma/w0 fitted to ln A on the records t >= T with 10 A(T) <= A <= GROWTH_AMAX (None with fewer than 3).
    The max over all modes is not used for the fit: the round-off of the other families grows with the number of
    steps and hides the slow growth (k0 = (1, 0, 1), N = 16: 4e-14 at 10 T with 3200 steps per period)."""
    o = gen_ic.read_out3d(os.path.join(case, 'out3d'), g)
    t, first = np.unique(o['t'], return_index=True)   # the last record may be written twice
    k0, amp = np.array(info['kvec']), info['amp']
    IX, IZ, IY = np.meshgrid(g.ikx, g.ikz, g.iky, indexing='ij')
    keep = np.ones(IX.shape, bool)
    for sgn in (-1, 0, 1):
        keep &= ~((IX == sgn*k0[0]) & (IY == sgn*k0[1]) & (IZ == sgn*k0[2]))
    a = np.array([np.where(keep, np.sqrt(sum(abs(o[n][r])**2 for n in ('mx', 'my', 'mz', 'bx', 'by', 'bz'))), 0.0)
                  for r in first])/amp
    A_all = a.reshape(len(t), -1).max(axis=1)
    rd = max(r for r in range(len(t)) if A_all[r] <= GROWTH_AMAX)
    i, k, j = np.unravel_index(int(np.argmax(a[rd])), a[rd].shape)
    kd = (int(g.ikx[i]), int(g.iky[j]), int(g.ikz[k]))
    cr, c0 = np.cross(np.stack([IX, IY, IZ], axis=-1), k0), np.cross(kd, k0)
    fam = keep & (np.all(cr == c0, axis=-1) | np.all(cr == -c0, axis=-1))
    A = np.sqrt((a[:, fam]**2).sum(axis=1))
    i1 = int(np.argmin(abs(t - info['T'])))
    fit = [r for r in range(i1, len(t)) if 10*A[i1] <= A[r] <= GROWTH_AMAX and A[r] > 0]
    gamma = (float(np.polyfit(t[fit], np.log(A[fit]), 1)[0])/np.linalg.norm(k0) if len(fit) >= 3 else None)
    return dict(gamma=gamma, nfit=len(fit), A_T=float(A[i1]), A_end=float(A[-1]), A_all_end=float(A_all[-1]),
                t_end=float(t[-1]/info['T']), k_dom=kd, k_perp=cp.k_perp(k0, kd))


def physical_max_cached(root, k0, amp, N, active):
    """cp_stability.physical_max over the families of the grid, cached in <root>/cp_rates.json."""
    path = os.path.join(root, 'cp_rates.json')
    try:
        cache = json.load(open(path))
    except (OSError, ValueError):
        cache = {}
    key = '%s a=%g N=%d active=%s' % (tuple(int(x) for x in k0), amp, N, ''.join('1' if x else '0' for x in active))
    if key not in cache:
        g, q, g15 = cp.physical_max(k0, amp, N, active)
        cache[key] = dict(rate=g, q=q, rate_15nmax=g15)
        json.dump(cache, open(path, 'w'), indent=1)
    return cache[key]


def alfven(root):
    """CP Alfven waves of S4.1 (gen_stage4.alfven_cases), rho = 1, B = n + dB, u = -dB, |dB| = A, over 10 periods.
    Per run: the exit, the snapshots at every period, |rho - 1|, |u.n|/A <= TOL_INV and the RK3 solution of the
    eigenmode (TOL_PRED) at all periods (A = 1: at t = T), div B <= TOL_DIVB, and over the run (A = 1: the first
    period) zeta, eta <= TOL_NUM A dl and the energy removed by the dissipation over a period <= TOL_DISS E
    above the eta floor (see TOL_DISS).
    A = 1 is linearly unstable, and the round-off seeds the instability (cp_stability.py): per run, the growth
    rate gamma of the perturbation (offpump_growth) with, for the ideal runs, |gamma/gamma_galerkin - 1| <=
    GROWTH_TOL for the dominant family (the instability of the Galerkin truncation of the code), and for the
    default runs with N >= GROWTH_NMIN, gamma <= GROWTH_PHYS x the largest physical rate of the grid (or no
    measurable growth); the ideal runs may be stopped by the CFL guard (exit_status dt_stop) after t = T.
    Across runs (A = 1: at t = T only; the last period is recorded): default - ideal <= TOL_DEF, the MPI
    decomposition (all snapshots) <= TOL_MPI, the axis permutations and the N independence <= TOL_NUM, and
    error(dt)/error(dt/2) in DT3."""
    import gen_stage4 as gs
    cases, snap = {}, {}
    for case in sorted(glob.glob(os.path.join(root, '*_p[0-9]*'))):
        name = os.path.basename(case)
        info = case_info(case)
        g = grid_of(case)
        T, amp, kvec = info['T'], info['amp'], info['kvec']
        r = dict(status=status(case), info=info)
        cases[name] = r
        stopped = amp >= 1.0 and info['setting'] == 'ideal' and r['status'] == 'dt_stop'
        if r['status'] != 'exit 0' and not stopped:
            r['ok'] = False
            continue
        basis = tuple(np.array(v) for v in info['basis'])
        n = basis[0]
        s = snapshots(case, g, T)
        snap[name] = s
        err = {}
        for m, (t, rho, u, b) in sorted(s.items()):
            _, uex, bex = gs.cp_fields(g, kvec, amp, basis, t=t)
            err[m] = dict(B=max(abs(b[d] - bex[d]).max() for d in range(3))/amp,
                          u=max(abs(u[d] - uex[d]).max() for d in range(3))/amp,
                          rho=float(abs(rho - 1).max()),
                          upar=float(abs(sum(u[d]*n[d] for d in range(3))).max())/amp)
        # RK3 solution of the eigenmode (exp(-i w t), w = |k| B0): the CP subspace is invariant under the
        # stages and the RHS is linear in it, so that the scheme gives the exact wave with exp(-i w t) replaced
        # by R(-i w dt)^n, R(z) = 1 + z + z^2/2 + z^3/6: dB_num = |R^n| dB_ex(t'), t' = -arg(R^n)/w.
        # pred = max|B - B_num|/A (round-off if the scheme is RK3 on the exact spatial operator), and
        # B_pred = max|B_num - B_ex|/A on the grid (the error the scheme should have)
        w, dt = np.linalg.norm(kvec), T/info['n_p']
        z = -1j*w*dt
        R = 1 + z + z**2/2 + z**3/6
        for m, e in err.items():
            t, _, u, b = s[m]
            Rn = R**int(round(t/dt))
            _, un, bn = gs.cp_fields(g, kvec, amp, basis, t=-np.angle(Rn)/w)
            bn = tuple(n[d] + abs(Rn)*(bn[d] - n[d]) for d in range(3))
            un = tuple(abs(Rn)*x for x in un)
            _, _, bex = gs.cp_fields(g, kvec, amp, basis, t=t)
            e['pred'] = max(max(abs(b[d] - bn[d]).max(), abs(u[d] - un[d]).max()) for d in range(3))/amp
            e['B_pred'] = max(abs(bn[d] - bex[d]).max() for d in range(3))/amp
        r['err'] = err
        r['nsnap'] = len(s)
        r['pred_T'] = err[1]['pred'] if 1 in err else float('inf')
        r['pred_all'] = max(e['pred'] for e in err.values())
        c, f, bg = load(os.path.join(case, 'conservation.dat')), load(os.path.join(case, 'cfl.dat')), \
            load(os.path.join(case, 'budget.dat'))
        dl = 2*np.pi/max(g.nl)
        E0, V = c[0, 10], c[0, 2]           # V = M0 (rho = 1)
        r['divb'] = float(c[:, 18:20].max())
        r['dM'] = float(abs(c[:, 3]).max())
        r['dP'] = float(abs(c[:, 4:7]).max())
        r['dE'] = float(abs(c[:, 10]/E0 - 1).max())
        # cfl.dat (0-based) 7: zeta_max, 16: eta_max. A = 1: judged on the first period, before the
        # instability has grown (D = 0 holds for the wave only); the whole run in *_all
        for sfx, wf, wb in (('_all', slice(None), slice(None)),
                            ('', f[:, 0] <= T*(1 + 1e-9) if amp >= 1.0 else slice(None),
                             bg[:, 0] <= T*(1 + 1e-9) if amp >= 1.0 else slice(None))):
            r['zeta' + sfx] = float(f[wf, 7].max())/(amp*dl)
            r['eta' + sfx] = float(f[wf, 16].max())/(amp*dl)
            r['diss' + sfx] = float(abs(bg[wb][:, [10, 20, 21, 22, 26]]).max())*T/E0
            r['diss_floor' + sfx] = float(f[wf, 16].max())*float(np.dot(kvec, kvec))*amp**2*V*T/E0
        r['dt_margin'] = float((np.minimum(f[:, 2], f[:, 3])/f[:, 1]).min())
        r['R'] = float(abs(bg[:, [12, 23, 28]]).max())*T/E0
        r['err_t'] = float(abs(bg[:, 8]).max())*T/E0
        last = max(err)
        r['inv_T'] = max(err[1]['rho'], err[1]['upar'])
        r['inv_all'] = max(max(e['rho'], e['upar']) for e in err.values())
        # growth of |rho - 1| (round-off seeded instability of the finite-amplitude wave): rate per period
        ms = [m for m in sorted(err) if m >= 1 and err[m]['rho'] > 0]
        r['rho_growth'] = (float(np.polyfit(ms, np.log(np.array([err[m]['rho'] for m in ms])), 1)[0])
                           if len(ms) >= 3 else None)
        complete = r['nsnap'] == info['nper'] + 1 and last == info['nper']
        r['ok'] = ((complete or (stopped and 1 in err)) and r['inv_T'] <= TOL_INV and r['divb'] <= TOL_DIVB
                   and r['pred_T'] <= TOL_PRED and (amp >= 1.0 or r['pred_all'] <= TOL_PRED)
                   and r['zeta'] <= TOL_NUM and r['eta'] <= TOL_NUM and r['diss'] <= TOL_DISS + r['diss_floor']
                   and (amp >= 1.0 or r['inv_all'] <= TOL_INV))
        gr = offpump_growth(case, g, info)
        N, active = info['nl'][0], tuple(n > 2 for n in g.nl)
        # the rates of the family of k_dom (k_dom along k0: the family of the wave itself, neutral, see
        # cp_stability.py check)
        fam = any(np.cross(gr['k_dom'], kvec))
        gr['galerkin'] = cp.galerkin_rate(kvec, amp, gr['k_dom'], N, active) if fam else 0.0
        gr['physical'] = cp.physical_rate(kvec, amp, gr['k_dom']) if fam else 0.0
        r['growth'] = gr
        if amp >= 1.0:
            if info['setting'] == 'ideal':
                gr['judged'] = 'galerkin'
                gr['ok'] = (gr['gamma'] is not None and gr['galerkin'] > 0
                            and abs(gr['gamma']/gr['galerkin'] - 1) <= GROWTH_TOL)
            elif N >= GROWTH_NMIN:
                pm = physical_max_cached(root, kvec, amp, N, active)
                gr['physical_max'] = pm
                gr['judged'] = 'physical_max'
                gr['ok'] = gr['gamma'] is None or gr['gamma'] <= GROWTH_PHYS*pm['rate']
            else:
                gr['judged'] = 'no'
                gr['ok'] = True
            r['ok'] &= gr['ok']

    def key(i, **kw):
        d = dict(dim=i['dim'], kvec=tuple(i['kvec']), amp=i['amp'], setting=i['setting'], n_p=i['n_p'],
                 N=i['nl'][0], nproc=i['nproc'], perm=i['perm'])
        d.update(kw)
        return tuple(sorted(d.items()))
    by_key = {key(r['info']): name for name, r in cases.items()}

    def field_diff(a, b, amp, ms, maps=lambda v: v, sub=None):
        """max over the snapshots ms and the components of |a - b|/amp for (u, B), b mapped by maps and
        sampled on a's grid by sub."""
        out = 0.0
        for m in ms:
            if m not in snap[a] or m not in snap[b]:
                return float('inf')
            for q in (2, 3):
                va, vb = snap[a][m][q], maps(snap[b][m][q])
                for d in range(3):
                    x = vb[d] if sub is None else vb[d][sub]
                    out = max(out, float(abs(va[d] - x).max())/amp)
        return out

    checks = {}
    ok_all = all(r['ok'] for r in cases.values())
    for name, r in cases.items():
        if name not in snap:
            continue
        i = r['info']
        amp, nper = i['amp'], i['nper']
        # default vs ideal (t = T and the last period)
        if i['setting'] == 'default':
            ref = by_key.get(key(i, setting='ideal'))
            if ref in snap:
                d1, dn = field_diff(name, ref, amp, [1]), field_diff(name, ref, amp, [nper])
                checks['def_vs_ideal ' + name] = dict(T=d1, last=dn, ok=d1 <= TOL_DEF and (amp >= 1.0 or dn <= TOL_DEF))
        # MPI decomposition
        if i['nproc'] > 1:
            ref = by_key.get(key(i, nproc=1))
            if ref in snap:
                d = field_diff(name, ref, amp, sorted(snap[name]))
                checks['mpi ' + name] = dict(diff=d, ok=d <= TOL_MPI)
        # axis permutation: the problem rolled back to kvec0
        if i['perm'] > 0:
            k0 = tuple(int(k) for k in np.roll(i['kvec'], -i['perm']))
            ref = by_key.get(key(i, kvec=k0, perm=0))
            if ref in snap:
                g = grid_of(os.path.join(root, name))
                mp = lambda v, s=i['perm']: roll_fields(v, g, s)
                d1, dn = field_diff(name, ref, amp, [1], maps=mp), field_diff(name, ref, amp, [nper], maps=mp)
                checks['perm ' + name] = dict(T=d1, last=dn, ok=d1 <= TOL_NUM and (amp >= 1.0 or dn <= TOL_NUM))
        # N independence: the error field on a finer grid sampled on the coarsest one
        Ns = sorted({c['info']['nl'][0] for c in cases.values() if key(c['info'], N=0, nproc=0) == key(i, N=0, nproc=0)})
        if i['nl'][0] == Ns[0] and len(Ns) > 1 and i['nproc'] == 1:
            for N in Ns[1:]:
                fine = [nm for nm, c in cases.items()
                        if key(c['info'], nproc=0) == key(i, N=N, nproc=0)]
                if not fine or fine[0] not in snap:
                    continue
                f = fine[0]
                st = N//Ns[0]
                g = grid_of(os.path.join(root, name))
                sub = (slice(None, None, st), slice(None, None, st), slice(None, None, st if g.nly > 2 else 1))
                d1, dn = field_diff(name, f, amp, [1], sub=sub), field_diff(name, f, amp, [nper], sub=sub)
                checks['N%d vs N%d %s' % (N, Ns[0], name)] = dict(T=d1, last=dn,
                                                                  ok=d1 <= TOL_NUM and (amp >= 1.0 or dn <= TOL_NUM))
    # dt^3: error at T and at the last period for n_p, 2 n_p, 4 n_p
    for name, r in cases.items():
        if 'err' not in r:
            continue
        i = r['info']
        seq = [by_key.get(key(i, n_p=m*i['n_p'])) for m in (1, 2, 4)]
        if all(s in cases and 'err' in cases[s] for s in seq) and not any(
                key(i, n_p=i['n_p']//2) == key(c['info']) for c in cases.values()):
            ratio = {}
            for m in (1, i['nper']):
                if not all(m in cases[s]['err'] for s in seq):
                    continue
                e = [cases[s]['err'][m]['B'] for s in seq]
                ratio[m] = (e[0]/e[1], e[1]/e[2], e)
            judged = [1] if i['amp'] >= 1.0 else [1, i['nper']]
            okr = all(m in ratio and DT3[0] <= x <= DT3[1] for m in judged for x in ratio.get(m, (0, 0))[:2])
            checks['dt3 ' + name] = dict(ratio={m: v[:2] for m, v in ratio.items()},
                                         err={m: v[2] for m, v in ratio.items()}, ok=okr)
    ok_all &= all(c['ok'] for c in checks.values())
    for name, r in cases.items():
        if 'err' not in r:
            print('%-36s %-8s ok=False' % (name, r['status']))
            continue
        e1, en = r['err'][1], r['err'][max(r['err'])]
        gr = r['growth']
        print('%-36s ok=%-5s errB(T) %.2e errB(end) %.2e |rho-1| %.1e/%.1e u.n %.1e divB %.1e zeta %.1e/%.1e '
              'eta %.1e/%.1e diss %.1e/%.1e (floor %.1e) R %.1e dE %.1e dt_margin %.2f pred %.1e/%.1e'
              % (name, r['ok'], e1['B'], en['B'], r['inv_T'], r['inv_all'], en['upar'], r['divb'], r['zeta'],
                 r['zeta_all'], r['eta'], r['eta_all'], r['diss'], r['diss_all'], r['diss_floor'], r['R'], r['dE'],
                 r['dt_margin'],
                 r['pred_T'], r['pred_all']))
        print('%-36s %s  gamma/w0 %s (%d records; family A(T) %.1e A(%.2fT) %.1e, all modes %.1e)  k_dom %s '
              'k_perp %.2f: galerkin %.3f physical %.3f%s  judged %s ok=%s'
              % ('', r['status'], '-' if gr['gamma'] is None else '%.3f' % gr['gamma'], gr['nfit'], gr['A_T'],
                 gr['t_end'], gr['A_end'], gr['A_all_end'], gr['k_dom'], gr['k_perp'], gr['galerkin'], gr['physical'],
                 '  physical max %.3f (q %s)' % (gr['physical_max']['rate'], tuple(gr['physical_max']['q']))
                 if 'physical_max' in gr else '', gr.get('judged', '-'), gr.get('ok', '-')))
    for k, c in checks.items():
        print('%-60s ok=%-5s %s' % (k, c['ok'], {a: b for a, b in c.items() if a != 'ok'}))
    print('alfven:', 'PASS' if ok_all else 'FAIL')
    return dict(cases=cases, checks=checks, ok=ok_all)


def cpdiag(root):
    """Growth of the A = 1 CP waves with the dissipation terms of the default settings one at a time
    (gen_stage4.batch_cpdiag) next to the ideal and default runs of the same (dim, k, N, n_p, nproc) in
    <root>/../alfven: gamma/w0 and the dominant k of offpump_growth, the Galerkin and physical rates of its
    family, and max eta, zeta/(A dl). A diagnostic: recorded, not judged (ok: every run finished). Per group:
    gamma/gamma_galerkin of noeta and etaonly (dominant family), and |A_end(nosvv)/A_end(etaonly) - 1| (zeta
    inert at linear order: q = rho zeta D is of second order)."""
    import gen_stage4 as gs
    ref_root = os.path.join(os.path.dirname(os.path.normpath(root)), 'alfven')
    groups = {}
    for case in sorted(glob.glob(os.path.join(root, '*_p[0-9]*'))):
        i = case_info(case)
        k = (i['dim'], tuple(i['kvec']), i['nl'][0], i['n_p'], i['nproc'])
        groups.setdefault(k, {})[i['setting']] = case
        for setting in ('ideal', 'default'):
            ref = case.replace(root.rstrip('/'), ref_root).replace('_%s_' % i['setting'], '_%s_' % setting)
            if os.path.isdir(ref):
                groups[k][setting] = ref
    res, ok = {}, True
    for k, runs in sorted(groups.items()):
        dim, kvec, N, n_p, P = k
        print('== %s k = %s N = %d n_p = %d nproc = %d' % (dim, kvec, N, n_p, P))
        for setting in [x for x in gs.SETTINGS if x in runs]:
            case = runs[setting]
            name = os.path.basename(case)
            st = status(case)
            r = dict(status=st)
            res[name] = r
            if st not in ('exit 0', 'dt_stop'):
                ok = False
                print('   %-8s %s' % (setting, st))
                continue
            i, g = case_info(case), grid_of(case)
            gr = offpump_growth(case, g, i)
            active = tuple(n > 2 for n in g.nl)
            fam = any(np.cross(gr['k_dom'], kvec))
            gr['galerkin'] = cp.galerkin_rate(kvec, i['amp'], gr['k_dom'], N, active) if fam else 0.0
            gr['physical'] = cp.physical_rate(kvec, i['amp'], gr['k_dom']) if fam else 0.0
            f = load(os.path.join(case, 'cfl.dat'))
            dl = 2*np.pi/max(g.nl)
            r.update(growth=gr, zeta=float(f[:, 7].max())/(i['amp']*dl), eta=float(f[:, 16].max())/(i['amp']*dl))
            print('   %-8s %-8s gamma/w0 %s (%d records)  A(%.2fT) %.1e  k_dom %-12s k_perp %5.2f: galerkin %.3f '
                  'physical %.3f  max zeta %.1e eta %.1e'
                  % (setting, st, '-' if gr['gamma'] is None else '%.3f' % gr['gamma'], gr['nfit'], gr['t_end'],
                     gr['A_end'], gr['k_dom'], gr['k_perp'], gr['galerkin'], gr['physical'], r['zeta'], r['eta']))
        g_of = {x: res[os.path.basename(runs[x])].get('growth') for x in runs}
        rel = lambda x: ('-' if not g_of.get(x) or g_of[x]['gamma'] is None or g_of[x]['galerkin'] <= 0
                         else '%.2f' % (g_of[x]['gamma']/g_of[x]['galerkin']))
        za = ('%.1e' % abs(g_of['nosvv']['A_end']/g_of['etaonly']['A_end'] - 1)
              if g_of.get('nosvv') and g_of.get('etaonly') else '-')
        print('   gamma/galerkin: noeta %s etaonly %s default %s  |A_end(nosvv)/A_end(etaonly) - 1| %s'
              % (rel('noeta'), rel('etaonly'), rel('default'), za))
    print('cpdiag:', 'all runs finished' if ok else 'runs missing or failed')
    return dict(runs=res, ok=ok)


# ------------------------------------------------------------------ S4.2 planar problems along an oblique normal
TOL_CTL = 1e-11    # k = (1, 0): 2D against 1D lines (W absolute, the states being O(1); zeta, D, eta relative to
                   # the larger of their max and their rounding scale, see rounding_floor)
TOL_EQV = 1e-10    # k = (1, 1) without SVV: 2D against the compressed 1D run, as TOL_CTL
TOL_SYM = 1e-12    # max|f(i, l) - f(p(i, l), 0)| on the planes, as TOL_CTL
TOL_RES = 1e-10    # R_zeta, R_eta of budget.dat relative to the dissipation, R_svv to its rounding scale (svv_res_rel)
PHY_L1 = 1.3       # L1(2D)/L1(compressed 1D with the bandwidth and Delta of the 2D run: on N points for k = (1, 1),
                   # N/2 for k = (1, 2), see gen_stage4.planar_bw_case) of the summed L1 error of W
PHY_L1_STOP = 2.0  # above: the anisotropy of the mask and the stencil is investigated before S4.3
PHY_TD_UN = 1.5    # max|u_n - U| of TDS against the 1D runs with the bandwidth of the 2D run (as PHY_L1)
PHY_RD_DM = 2e-5   # |M - M0|/M0 of RDpc
PLANE_SCALE1 = ('rho', 'ux', 'uy', 'uz', 'bx', 'by', 'bz')


def svv_res_rel(case):
    """max|R_svv| relative to its rounding scale sum|dE/dU.U_svv| dV (budget.dat column 26). |R_svv|/D_svv of
    analyze_stage3.stats has no meaning where D_svv is itself at the rounding level (the smooth initial data of
    OT: D_svv ~ 1e-25 at t = 0, the ratio up to 1e-3 at N = 256)."""
    return float(np.abs(np.concatenate(load_runs(case, 'budget.dat'))[:, 25]).max())


def rounding_floor(rho, u2, b2, dl):
    """Rounding scales of D, zeta and eta, the floor of the scale their differences are divided by. Where D = 0
    holds analytically (RD, TD, OT at t = 0), max|D| is itself round-off and the difference relative to it O(1).
    The round-off of u is set by the forces, not by |u| (TDS: |u| ~ 1e-3 next to a magnetic pressure of O(1)), so
    the speed is U = max(|u| + c_f), c_f^2 = c_s^2 + |B|^2/rho; D: U/dl (the spectral derivative), zeta, eta: U dl
    (Delta^2 D without the coefficient)."""
    import gen_stage4 as gs
    U = float(np.max(np.sqrt(u2) + np.sqrt(gs.CS2 + b2/rho)))
    return dict(D=U/dl, zeta=U*dl, eta=U*dl)


def plane_floor(f, dl):
    """rounding_floor of a plane record f (rho, ux, ..., bz)."""
    return rounding_floor(f['rho'], f['ux']**2 + f['uy']**2 + f['uz']**2, f['bx']**2 + f['by']**2 + f['bz']**2, dl)


def planar_snaps(case, m):
    """[(t, W (6, N), zeta, eta, D)] of the out1d lines of a planar run at the steps j nrun/nline (2D: the line
    along x at z = 0, where p = i; compressed 1D: the run along x), in the units of the 1D problem of Stage 3 (t,
    zeta, eta times |k|, D over |k|) with W along (n, t1, t2) of gen_stage4.planar_frame (1D: x, y, z) on the
    grid of Stage 3, the steps of the lines and max|B_n - b_x| over them."""
    import gen_stage4 as gs
    c, kn = m['case'], m['kn']
    N = c['N']
    if m['dim'] == '2d':
        frame, x = np.array(gs.planar_frame(m['kvec'])), np.arange(N)*(s3.L/N)
    else:
        frame, x = np.eye(3), np.arange(N)*(s3.L/kn/N)
    m_w = m['nrun']//m['nline']
    out, steps, dbn = [], [], 0.0
    for t, step, a in read_lines(case):
        n = int(round(t/m['dt']))
        if n % m_w != 0 and n != m['nrun']:
            continue
        ncol = 11 if c['C_eta'] > 0.0 else 10
        assert a.shape == (N, ncol) and np.max(np.abs(a[:, 0] - x)) < 1e-12, (case, t, a.shape)
        u, b = a[:, 2:5] @ frame.T, a[:, 5:8] @ frame.T
        W = np.array([a[:, 1], u[:, 0], u[:, 1], u[:, 2], b[:, 1], b[:, 2]])
        out.append((kn*t, W, kn*a[:, 8], kn*a[:, 10] if ncol == 11 else np.zeros(N), a[:, 9]/kn))
        steps.append(step)
        dbn = max(dbn, float(np.max(np.abs(b[:, 0] - s3.bx_of(c)))))
    return out, steps, dbn


def lines_diff(sa, sb):
    """max over the lines of |a - b| for two runs with the same steps: W absolute, zeta, eta, D relative to the
    larger of their max in b and rounding_floor (in the units of the 1D problem, dl = L/N, B_n left out of c_f);
    max: the largest of these."""
    if not sa or len(sa) != len(sb):
        return dict(n=(len(sa), len(sb)), max=float('inf'))
    d = dict(n=len(sa), t=max(abs(a[0] - b[0]) for a, b in zip(sa, sb)),
             W=float(max(np.abs(a[1] - b[1]).max() for a, b in zip(sa, sb))))
    W = sb[0][1]
    floor = rounding_floor(W[0], (W[1:4]**2).sum(axis=0), (W[4:6]**2).sum(axis=0), s3.L/W.shape[1])
    for k, q in ((2, 'zeta'), (3, 'eta'), (4, 'D')):
        sc = max(max(float(np.abs(b[k]).max()) for b in sb), floor[q])
        d[q] = float(max(np.abs(a[k] - b[k]).max() for a, b in zip(sa, sb)))/sc
    d['max'] = max(d['W'], d['zeta'], d['eta'], d['D'])
    return d


def planar_symmetry(case, m):
    """max over the records of out2d/plane_y0.dat of |f(i, l) - f(p(i, l), 0)| (rho, u, B absolute; zeta, D,
    eta relative to the larger of their max in the record and rounding_floor); max: the largest of these."""
    import gen_stage4 as gs
    g = grid_of(case)
    P = gs.planar_index(g.nlx, m['kvec'])
    d, nrec = {q: 0.0 for q in gen_ic.PLANE_Q}, 0
    for path in run_dirs(case):
        if not os.path.exists(os.path.join(path, 'out2d', 'plane_time.dat')):
            continue
        pl = gen_ic.read_planes(os.path.join(path, 'out2d'), g, 'y0')
        for r in range(len(pl['t'])):
            nrec += 1
            floor = plane_floor({q: pl[q][r] for q in PLANE_SCALE1}, g.Lx/g.nlx)
            for q in gen_ic.PLANE_Q:
                f = pl[q][r]
                sc = 1.0 if q in PLANE_SCALE1 else max(float(np.abs(f).max()), floor[q])
                d[q] = max(d[q], float(np.abs(f - f[P, 0]).max())/sc)
    d['nrec'], d['max'] = nrec, max(d[q] for q in gen_ic.PLANE_Q)
    return d


def s3_result(base, name):
    """stage3.json of the Stage 3 run of the case (runs/stage3_f/{suite, eta_suite} next to base)."""
    for b in ('suite', 'eta_suite'):
        f = os.path.join(base, '..', 'stage3_f', b, name, 'stage3.json')
        if os.path.exists(f):
            return json.load(open(f))
    return None


def l1_sum(an):
    """Summed L1 error of W at the final time against the exact solution (compound: the FV reference) and of
    B_t1, B_t2."""
    if 'l1' in an:
        l1 = an['l1']
    elif an.get('l1_fv'):
        l1 = an['l1_fv'][-1]['l1']
    else:
        return None, None
    return float(sum(l1)), float(l1[4] + l1[5])


def front_over(an):
    """Largest overshoot of rho and B_t1 next to the fronts (stage3.analyze, relative to the jump)."""
    fr = an.get('fronts', [])
    return {q: max([f.get(q, 0.0) or 0.0 for f in fr], default=None) for q in ('over_rho', 'over_By')}


def planar(root):
    import analyze_stage3 as ast3
    import gen_stage4 as gs
    base = os.path.dirname(os.path.abspath(root))
    cases = {}
    for f in sorted(glob.glob(os.path.join(root, '*', 'case.json'))):
        case = os.path.dirname(f)
        m = case_info(case)
        r = dict(status=status(case), info={k: m[k] for k in ('role', 'prob', 'dim', 'kvec', 'kn', 'nstep',
                                                               'nrun', 'dt')})
        r['info']['N'] = m['case']['N']
        cases[os.path.basename(case)] = r
        if r['status'] != 'exit 0':
            continue
        snaps, steps, dbn = planar_snaps(case, m)
        r['snaps'] = snaps
        r['dBn'] = dbn
        r['complete'] = bool(snaps) and abs(snaps[-1][0] - m['kn']*m['nrun']*m['dt']) <= 1e-9*snaps[-1][0]
        st = ast3.stats(case, dict(m, dir='x'), snaps)
        st['R_svv_A'] = svv_res_rel(case)
        cons = np.concatenate(load_runs(case, 'conservation.dat'))
        frame = np.array(gs.planar_frame(m['kvec'])) if m['dim'] == '2d' else np.eye(3)
        P = cons[:, 4:7] @ frame.T
        st.update({'drift_P' + q: float(np.abs(P[:, j] - P[0, j]).max()/cons[0, 2]) for j, q in
                   enumerate(('n', 't1', 't2'))})
        r['stats'] = {k: st.get(k) for k in ('drift_M', 'drift_Pn', 'drift_Pt1', 'drift_Pt2', 'drift_By', 'drift_Bz',
                                             'rho_min', 'rho_max', 'R_zeta', 'R_svv', 'R_svv_A', 'R_eta',
                                             'max_dt_over_limit',
                                             'max_dt_lam_svv', 'zeta_max', 'eta_max', 'dE', 'n_E_up', 'time_err')}
        if m['role'] == 'phy':
            with np.errstate(all='ignore'):
                r['analysis'] = s3.analyze(m['case'], snaps, s3.initial(m['case']))
        if m['dim'] == '2d':
            r['sym'] = planar_symmetry(case, m)
    checks, ok_all = {}, True

    def need(name):
        r = cases.get(name)
        return r if r is not None and r['status'] == 'exit 0' and r.get('complete') else None

    for name, r in cases.items():
        if r['status'] != 'exit 0' or not r.get('complete'):
            checks['run ' + name] = dict(ok=False, status=r['status'], complete=r.get('complete'))
            continue
        st = r['stats']
        res = max(v for v in (st['R_zeta'], st['R_svv_A'], st['R_eta'], 0.0) if v is not None)
        c = dict(R=res, drift_M=st['drift_M'], dBn=r['dBn'], ok=res <= TOL_RES)
        if 'sym' in r:
            c['sym'] = r['sym']['max']
            c['sym_nrec'] = r['sym']['nrec']
            c['ok'] = c['ok'] and r['sym']['nrec'] >= 2 and r['sym']['max'] <= TOL_SYM
        checks['run ' + name] = c
    for name in sorted(cases):
        if not name.endswith('_2d') or cases[name]['info']['role'] not in ('ctl', 'eqv'):
            continue
        a, b = need(name), need(name[:-3] + '_1d')
        tol = TOL_CTL if cases[name]['info']['role'] == 'ctl' else TOL_EQV
        d = lines_diff(a['snaps'], b['snaps']) if a and b else dict(max=float('inf'))
        checks['lines ' + name[:-3]] = dict(d, ok=d['max'] <= tol)
    # phy: 2D, compressed 1D and the 1D run of Stage 3
    phy = {}
    for name in sorted(cases):
        i = cases[name]['info']
        if i['role'] != 'phy' or i['dim'] != '2d':
            continue
        a, b, h = need(name), need(name[:-3] + '_1d'), need(name[:-3] + '_1dh')
        c3 = gs.s3_case(i['prob'], i['N'])
        o = s3_result(base, c3['name'])
        row = dict(prob=i['prob'], kvec=i['kvec'], N=i['N'])
        for tag, an, st in (('2d', a and a['analysis'], a and a['stats']), ('c1d', b and b['analysis'], b and b['stats']),
                            ('c1dh', h and h['analysis'], h and h['stats']),
                            ('1d', o and o['analysis'], o and o['attempts'][-1]['stats'])):
            if an is None:
                continue
            row['L1_' + tag], row['L1B_' + tag] = l1_sum(an)
            row.update({k + '_' + tag: v for k, v in front_over(an).items()})
            row['drift_M_' + tag] = st.get('drift_M')
            if c3['fam'] == 'td':
                row['un_' + tag] = an['dev_max']['ux']
                row['width_' + tag] = an['width_dl']['final']
            if c3['fam'] == 'rd':
                row['speed_err_' + tag] = an['speed_err']
        # reference with the bandwidth of the 2D run along n (gen_stage4.planar_bw_case)
        refs = ('c1dh',) if gs.planar_bw_case(c3, i['kvec']) is not None else ('c1d', '1d')
        ok = True
        if row.get('L1_2d') and row.get('L1_c1d'):
            row['r_2d_c1d'] = row['L1_2d']/row['L1_c1d']
        if row.get('L1_2d') and row.get('L1_' + refs[0]):
            row['ref'] = refs[0]
            row['r_2d_ref'] = row['L1_2d']/row['L1_' + refs[0]]
            ok &= row['r_2d_ref'] <= PHY_L1
            row['investigate'] = row['r_2d_ref'] > PHY_L1_STOP
        else:
            ok = False
        if row.get('L1_1d'):
            for tag in ('2d', 'c1d', 'c1dh'):
                if row.get('L1_' + tag):
                    row['r_%s_1d' % tag] = row['L1_' + tag]/row['L1_1d']
        if c3['fam'] == 'td' and 'un_2d' in row:
            ref = max(row.get('un_' + t, 0.0) for t in refs)
            row['r_un'] = row['un_2d']/ref if ref > 0 else None
            ok &= row['r_un'] is not None and row['r_un'] <= PHY_TD_UN
        if c3['name'].startswith('RDpc') and 'drift_M_2d' in row:
            ok &= row['drift_M_2d'] <= PHY_RD_DM
        row['ok'] = ok
        phy[name[:-3]] = row
        checks['phy ' + name[:-3]] = row
    # C_eta = 0 against 2 (fast shock and compound): L1(B) and the overshoot of B_t1
    eta = {}
    for name, row in phy.items():
        if '_eta2_' not in name:
            continue
        r0 = phy.get(name.replace('_eta2', ''))
        if r0 is None:
            continue
        eta[name] = {tag: dict(L1B=(r0.get('L1B_' + tag), row.get('L1B_' + tag)),
                               over_By=(r0.get('over_By_' + tag), row.get('over_By_' + tag)))
                     for tag in ('2d', 'c1d', '1d')}
    ok_all = all(c['ok'] for c in checks.values())
    for name, r in sorted(cases.items()):
        st = r.get('stats', {})
        print('%-28s %-8s complete=%-5s dM %s R_zeta %s R_svv/D %s R_svv/A %s R_eta %s sym %s dt/lim %s'
              % (name, r['status'], r.get('complete'), st.get('drift_M'), st.get('R_zeta'), st.get('R_svv'),
                 st.get('R_svv_A'), st.get('R_eta'), r.get('sym', {}).get('max'), st.get('max_dt_over_limit')))
    for k, c in checks.items():
        if not k.startswith('run ') or not c['ok']:
            print('%-40s ok=%-5s %s' % (k, c['ok'], {a: b for a, b in c.items() if a != 'ok'}))
    for k, e in eta.items():
        print('eta %-30s %s' % (k, e))
    print('planar:', 'PASS' if ok_all else 'FAIL')
    for r in cases.values():
        r.pop('snaps', None)
    return dict(cases=cases, checks=checks, eta=eta, ok=ok_all)


# ------------------------------------------------------------ S4.3 OT vortex
OT_TIMES = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)  # spectra and zeta statistics at the plane records nearest to these t
OT_KEFF_TOL = 0.25    # k_eff: largest k with |E_N/E_2N - 1| <= OT_KEFF_TOL on 1..k (our threshold)
OT_ZETA_FRAC = 0.01   # zeta > OT_ZETA_FRAC zeta_max(t): the region of zeta (reported; it includes the smooth
                      # compressions, whose zeta/zeta_max ~ w Delta |grad u|/du is above 0.01 for N up to about 2000)
OT_A90 = 0.9          # a90: the smallest area fraction holding OT_A90 of the shock-viscous dissipation rho zeta D^2;
                      # that of a shock is independent of Delta on an area ~ Delta, that of the smooth compressions ~ Delta^2
OT_CONV_STATE = 0.02  # top two N: rho_min, rho_max, E_K, E_M at t_end
OT_STATE_TIMES = (1.0, 2.0, np.pi, 4.0, 5.0, 6.0)  # the same, reported at these t (pi: t_end of the plan)
OT_CONV_EPS = 0.05    # peak of eps(t)
OT_CONV_TPEAK = 0.02  # time of the peak
OT_CONV_SHARE = 0.10  # share of the shock viscosity in the dissipation over [0, t_end]
OT_KEFF_RATIO = (1.7, 2.3)
OT_ZFRAC_RATIO = (0.4, 0.65)  # a90(2N)/a90(N)
TOL_DEC10 = 1e-12     # decomposition: planes after 10 steps against P = 16 (rho, u, B absolute; zeta, D, eta / max)
SPEC_Q = ('E_u', 'E_sru', 'E_B', 'E_tot', 'P_drho', 'P_lnrho', 'chi', 'E_zp', 'E_zm')


def plane_records(case, plane='y0'):
    """[(t, step, path of plane_<plane>.dat, record)] over all runs of the case (a continued run repeats the
    state at its start)."""
    out = []
    for d in run_dirs(case):
        fn = os.path.join(d, 'out2d', 'plane_time.dat')
        if os.path.exists(fn):
            for row in np.atleast_2d(np.loadtxt(fn)):
                out.append((float(row[1]), int(row[2]), os.path.join(d, 'out2d', 'plane_%s.dat' % plane),
                            int(row[0])))
    return out


def plane_at(g, path, rec, plane='y0'):
    """dict q -> (nlx, nlz) of one record of plane_y0.dat, (nlx, nly) of plane_z0.dat (memory map: a record of
    1024^2 is 84 MB)."""
    a = np.memmap(path, dtype='<f8', mode='r').reshape(-1, len(gen_ic.PLANE_Q), g.nlz if plane == 'y0' else g.nly,
                                                       g.nlx)
    return {q: np.array(a[rec, iq].T) for iq, q in enumerate(gen_ic.PLANE_Q)}


def shell_spectra(f):
    """Shell-summed spectra (bins |k| in [n - 1/2, n + 1/2), n = 0..N/2; L = 2 pi) of a plane f of plane_at:
    E_u = |u_k|^2/2, E_sru of sqrt(rho) u, E_B, E_tot = E_sru + E_B, P_drho = |rho_k|^2 (k > 0), P_lnrho,
    chi = compressive/total of u_k (k.u_k/|k|), E_zp, E_zm of the Elsaesser fields u +- B/sqrt(rho)."""
    N = f['rho'].shape[0]
    k1 = np.fft.fftfreq(N, 1.0/N)
    KX, KZ = np.meshgrid(k1, k1, indexing='ij')
    K = np.hypot(KX, KZ)
    n = np.rint(K).astype(int)
    nb = N//2 + 1
    sel = n < nb

    def shell(p):
        return np.bincount(n[sel], weights=p[sel], minlength=nb)[:nb]

    def fk(a):
        return np.fft.fft2(a)/a.size

    def e(vs):
        return shell(sum(0.5*np.abs(fk(v))**2 for v in vs))

    rho = f['rho']
    u = [f['ux'], f['uy'], f['uz']]
    b = [f['bx'], f['by'], f['bz']]
    sr, isr = np.sqrt(rho), 1/np.sqrt(rho)
    s = dict(E_u=e(u), E_sru=e([sr*v for v in u]), E_B=e(b))
    s['E_tot'] = s['E_sru'] + s['E_B']
    s['P_drho'] = shell(np.abs(fk(rho))**2*(K > 0))
    s['P_lnrho'] = shell(np.abs(fk(np.log(rho)))**2*(K > 0))
    uxk, uzk = fk(u[0]), fk(u[2])
    with np.errstate(invalid='ignore', divide='ignore'):
        comp = np.where(K > 0, np.abs(KX*uxk + KZ*uzk)**2/np.where(K > 0, K, 1)**2, 0.0)
        s['chi'] = shell(0.5*comp)/s['E_u']
    s['E_zp'] = e([v + w*isr for v, w in zip(u, b)])
    s['E_zm'] = e([v - w*isr for v, w in zip(u, b)])
    return {q: s[q] for q in SPEC_Q}


def zeta_stats(f):
    """Region zeta > OT_ZETA_FRAC zeta_max: its area fraction, its share of the shock-viscous dissipation density
    rho zeta D^2, and the fraction of it with D < 0; a90 (OT_A90)."""
    z = f['zeta']
    zmax = float(z.max())
    if zmax <= 0.0:
        return dict(zeta_max=zmax, frac=0.0, diss_share=None, compressive=None, a90=None)
    m = z > OT_ZETA_FRAC*zmax
    w = f['rho']*z*f['D']**2
    ws = np.sort(w.ravel())[::-1]
    a90 = float((np.searchsorted(np.cumsum(ws), OT_A90*ws.sum()) + 1)/ws.size) if ws.sum() > 0 else None
    return dict(zeta_max=zmax, frac=float(m.mean()), diss_share=float(w[m].sum()/w.sum()) if w.sum() > 0 else None,
                compressive=float((f['D'][m] < 0).mean()), a90=a90)


def k_eff(Ea, Eb):
    """Largest k with |Ea/Eb - 1| <= OT_KEFF_TOL for 1..k (Ea on N, Eb on 2N; the 2/3 mask of N removes Ea
    above about N/3)."""
    kmax = min(len(Ea), len(Eb))
    k = 0
    for j in range(1, kmax):
        if Eb[j] <= 0 or abs(Ea[j]/Eb[j] - 1) > OT_KEFF_TOL:
            break
        k = j
    return k


def ot_run(case, m):
    """Time series, budget, spectra and zeta statistics of one OT run."""
    import analyze_stage3 as ast3
    N = m['N']
    g = gen_ic.Grid(N, 2, N)
    r = dict(info={k: m.get(k) for k in ('role', 'N', 'C_eta', 'nstep', 'nrun', 'dt', 'nproc')})
    cons = np.concatenate(load_runs(case, 'conservation.dat'))
    bud = np.concatenate(load_runs(case, 'budget.dat'))
    te = m['nrun']*m['dt']
    r['t_last'] = float(cons[-1, 0])
    r['complete'] = abs(r['t_last'] - te) <= 1e-9*te
    st = ast3.stats(case, dict(dir='x', case=dict(N=N)), [])
    st['R_svv_A'] = svv_res_rel(case)
    r['stats'] = {k: st.get(k) for k in ('drift_M', 'drift_Px', 'drift_Py', 'drift_Pz', 'rho_min', 'rho_max',
                                         'zeta_max', 'eta_max', 'R_zeta', 'R_svv', 'R_svv_A', 'R_eta',
                                         'max_dt_over_limit',
                                         'max_dt_lam_svv', 'dE', 'n_E_up', 'time_err', 'int_rk_eps',
                                         'int_rk_ideal_eps', 'int_rk_visc_eps', 'int_rk_svv_eps',
                                         'int_rk_svv_b_eps', 'int_rk_eta_eps')}
    r['stats']['divb_rms'] = float(cons[:, 18].max())
    r['stats']['divb_max'] = float(cons[:, 19].max())
    tot = st.get('int_rk_eps')
    r['share'] = {k: (st['int_rk_%s_eps' % k]/tot if tot else None) for k in ('ideal', 'visc', 'svv', 'svv_b', 'eta')}
    last = cons[-1]
    r['end'] = dict(t=float(last[0]), E_K=float(last[7]), E_M=float(last[8]), E=float(last[10]),
                    rho_min=float(last[11]), rho_max=float(last[12]))
    # eps(t) = -dE/dt of each step and the dissipation rates of the mechanisms
    eps = -bud[:, 6]
    i = int(np.argmax(eps))
    r['eps_peak'] = dict(t=float(bud[i, 0] + 0.5*bud[i, 1]), eps=float(eps[i]))
    stride = max(1, len(bud)//2000)
    # the integrals of r['share'] over [0, t] at the end of the steps (and of the last step), for the shares at
    # a common t of runs that end at different t
    keep = np.unique(np.r_[np.arange(0, len(bud), stride), len(bud) - 1])
    r['cum'] = dict(t=(bud[keep, 0] + bud[keep, 1]).tolist(),
                    **{k: np.cumsum(bud[:, 1]*bud[:, c])[keep].tolist() for k, c in (('eps', 7), ('visc', 10))})
    r['series'] = dict(t=(bud[::stride, 0] + 0.5*bud[::stride, 1]).tolist(), eps=eps[::stride].tolist(),
                       visc=(-bud[::stride, 10]).tolist(), svv=(-bud[::stride, 20]).tolist(),
                       eta=(-bud[::stride, 26]).tolist(), ideal=(-bud[::stride, 9]).tolist(),
                       tc=cons[:, 0].tolist(), E_K=cons[:, 7].tolist(), E_M=cons[:, 8].tolist(),
                       rho_min=cons[:, 11].tolist(), rho_max=cons[:, 12].tolist(), zeta_max=cons[:, 13].tolist())
    recs = plane_records(case)
    r['records'] = len(recs)
    r['spec'], r['zeta'] = {}, {}
    for t in OT_TIMES + (te,):
        if not recs or t > te + 1e-9:
            continue
        tr, step, path, rec = min(recs, key=lambda x: abs(x[0] - t))
        if abs(tr - t) > 2*m['dt']:
            continue
        f = plane_at(g, path, rec)
        key = '%.4f' % t
        r['spec'][key] = dict(t=tr, **{q: v.tolist() for q, v in shell_spectra(f).items()})
        r['zeta'][key] = dict(t=tr, **zeta_stats(f))
    return r


def rel(a, b):
    return abs(a - b)/abs(b) if b else float('inf')


def plane_diff(fa, fb, dl):
    """max|fa - fb| of two planes: rho, u, B absolute; zeta, D, eta relative to the larger of max|fb| and
    rounding_floor (max|fb| alone is round-off where D = 0, OT at t = 0); max: the largest of these."""
    floor = plane_floor(fb, dl)
    d = {}
    for q in gen_ic.PLANE_Q:
        sc = 1.0 if q in PLANE_SCALE1 else max(float(np.abs(fb[q]).max()), floor[q])
        d[q] = float(np.abs(fa[q] - fb[q]).max())/(sc if sc > 0 else 1.0)
    d['max'] = max(d[q] for q in gen_ic.PLANE_Q)
    return d


def ot_decomp(root, cases):
    """Planes after 10 steps and the per-step conservation.dat up to OT_DECOMP_T against P = 16."""
    out = {}
    ref10, ref = os.path.join(root, 'dec_N256_np16_n10'), os.path.join(root, 'dec_N256_np16')
    if cases.get('dec_N256_np16_n10', {}).get('status') != 'exit 0':
        return out
    g = gen_ic.Grid(256, 2, 256)
    for P in (4, 64, 256):
        name = 'dec_N256_np%d' % P
        d = {}
        if cases.get(name + '_n10', {}).get('status') == 'exit 0':
            ra, rb = plane_records(os.path.join(root, name + '_n10')), plane_records(ref10)
            d['n10'] = {}
            for (ta, sa, pa, ia), (tb, sb, pb, ib) in zip(ra, rb):
                diff = plane_diff(plane_at(g, pa, ia), plane_at(g, pb, ib), 2*np.pi/256)
                d['n10'][str(sb)] = dict(step=(sa, sb), **diff)
            last = d['n10'][max(d['n10'], key=int)]
            d['ok10'] = last['step'][0] == last['step'][1] == 10 and last['max'] <= TOL_DEC10
        if cases.get(name, {}).get('status') == 'exit 0' and cases.get('dec_N256_np16', {}).get('status') == 'exit 0':
            a = np.concatenate(load_runs(os.path.join(root, name), 'conservation.dat'))
            b = np.concatenate(load_runs(ref, 'conservation.dat'))
            n = min(len(a), len(b))
            cols = dict(E_K=7, E_M=8, E=10, rho_min=11, rho_max=12, zeta_max=13)
            growth = {}
            for tmax in (0.01, 0.1, 0.5, 1.0):
                sel = b[:n, 0] <= tmax + 1e-12
                # relative to the largest |value| of the column up to tmax (zeta_max = 0 at t = 0)
                growth['%g' % tmax] = {q: float(np.max(np.abs(a[:n][sel, j] - b[:n][sel, j]))
                                                / max(float(np.max(np.abs(b[:n][sel, j]))), 1e-300))
                                       for q, j in cols.items()}
            d['cons'] = growth
            ra, rb = plane_records(os.path.join(root, name)), plane_records(ref)
            d['planes'] = [dict(t=tb, step=(sa, sb), **plane_diff(plane_at(g, pa, ia), plane_at(g, pb, ib),
                                                                  2*np.pi/256))
                           for (ta, sa, pa, ia), (tb, sb, pb, ib) in zip(ra, rb)]
        out[name] = d
    return out


def ot_convergence(runs):
    """Top two N of each C_eta: the checks of S4.3; k_eff, the ratio of the zeta area fractions (zfrac, reported)
    and of a90 (judged) for each pair (N, 2N)."""
    conv = {}
    for e in sorted({r['info']['C_eta'] for r in runs.values()}):
        rs = sorted((r for r in runs.values() if r['info']['C_eta'] == e and r.get('complete')),
                    key=lambda r: r['info']['N'])
        Ns = [r['info']['N'] for r in rs]
        c = dict(Ns=Ns, pairs={})
        for ra, rb in zip(rs, rs[1:]):
            p = dict(keff={}, zfrac={}, a90={})
            for key, sa in ra['spec'].items():
                sb = rb['spec'].get(key)
                if sb:
                    p['keff'][key] = k_eff(sa['E_tot'], sb['E_tot'])
                za, zb = ra['zeta'].get(key), rb['zeta'].get(key)
                if za and zb and za['frac'] > 0:
                    p['zfrac'][key] = zb['frac']/za['frac']
                if za and zb and za.get('a90') and zb.get('a90'):
                    p['a90'][key] = zb['a90']/za['a90']
            c['pairs']['%d-%d' % (ra['info']['N'], rb['info']['N'])] = p
        if len(rs) >= 2:
            a, b = rs[-2], rs[-1]
            # judged at the earlier end of the two (N = 1024 runs up to t = pi, the t_end of the plan)
            tj = c['t_judge'] = min(a['t_last'], b['t_last'])
            at = lambda r, q: r['end'][q] if abs(r['t_last'] - tj) <= 1e-9*tj else series_at(r, q, tj)
            share = lambda r: float(np.interp(tj, r['cum']['t'], r['cum']['visc'])
                                    / np.interp(tj, r['cum']['t'], r['cum']['eps']))
            c['state'] = {q: rel(at(a, q), at(b, q)) for q in ('rho_min', 'rho_max', 'E_K', 'E_M')}
            c['state_t'] = {'%.4f' % t: {q: rel(series_at(a, q, t), series_at(b, q, t))
                                         for q in ('rho_min', 'rho_max', 'E_K', 'E_M')}
                            for t in OT_STATE_TIMES if t <= tj + 1e-9}
            # the peak of eps (t = 1.68 for N = 128-512) must lie in [0, t_judge] for both
            assert max(a['eps_peak']['t'], b['eps_peak']['t']) <= tj
            c['eps_peak'] = rel(a['eps_peak']['eps'], b['eps_peak']['eps'])
            c['t_peak'] = rel(a['eps_peak']['t'], b['eps_peak']['t'])
            c['share_visc'] = rel(share(a), share(b)) if b['share']['visc'] else None
            pr = list(c['pairs'].values())
            ok = (max(c['state'].values()) <= OT_CONV_STATE and c['eps_peak'] <= OT_CONV_EPS
                  and c['t_peak'] <= OT_CONV_TPEAK and c['share_visc'] is not None
                  and c['share_visc'] <= OT_CONV_SHARE)
            # band and shock width: the ratio over successive pairs needs three N
            if len(pr) >= 2:
                c['keff_ratio'] = {k: pr[-1]['keff'][k]/pr[-2]['keff'][k] for k in pr[-1]['keff']
                                   if pr[-2]['keff'].get(k)}
                ok &= all(OT_KEFF_RATIO[0] <= v <= OT_KEFF_RATIO[1] for v in c['keff_ratio'].values())
            c['zfrac_last'] = pr[-1]['zfrac']
            c['a90_last'] = pr[-1]['a90']
            ok &= bool(pr[-1]['a90']) and all(OT_ZFRAC_RATIO[0] <= v <= OT_ZFRAC_RATIO[1] for v in pr[-1]['a90'].values())
            c['ok'] = bool(ok)
        conv['eta%g' % e] = c
    return conv


def ot_figures(root, runs):
    """eps(t), E_tot(k) at t = 1 and 3 and the zeta maps at t = 3 of the top two N (C_eta = 0, the default, which
    alone has N = 1024) as PNG files in root."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    Ns = sorted({r['info']['N'] for r in runs.values()})
    shade = {N: plt.cm.Blues(0.35 + 0.65*i/max(1, len(Ns) - 1)) for i, N in enumerate(Ns)}
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    for name, r in sorted(runs.items(), key=lambda x: x[1]['info']['N']):
        e = r['info']['C_eta']
        a = ax[0 if e == 0 else 1]
        a.plot(r['series']['t'], r['series']['eps'], color=shade[r['info']['N']], lw=1.5, label='N = %d' % r['info']['N'])
        a.set_title('C_eta = %g' % e)
        a.set_xlabel('t')
    ax[0].set_ylabel('eps = -dE/dt')
    for a in ax:
        a.legend(frameon=False, fontsize=8)
        a.grid(alpha=0.3, lw=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(root, 'fig_ot_eps.png'), dpi=110)
    plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    for j, key in enumerate(('1.0000', '3.0000')):
        for name, r in sorted(runs.items(), key=lambda x: x[1]['info']['N']):
            s = r['spec'].get(key)
            if r['info']['C_eta'] != 0 or not s:
                continue
            k = np.arange(len(s['E_tot']))
            ax[j].loglog(k[1:], np.array(s['E_tot'][1:]), color=shade[r['info']['N']], lw=1.5,
                         label='N = %d' % r['info']['N'])
            ax[j].axvline(r['info']['N']/3, color=shade[r['info']['N']], lw=0.6, ls=':')
        ax[j].set_title('E_tot(k), t = %s, C_eta = 0' % key[:3])
        ax[j].set_xlabel('k')
        ax[j].legend(frameon=False, fontsize=8)
        ax[j].grid(alpha=0.3, lw=0.5, which='both')
    fig.tight_layout()
    fig.savefig(os.path.join(root, 'fig_ot_spec.png'), dpi=110)
    plt.close(fig)
    top = [n for n, r in sorted(runs.items(), key=lambda x: x[1]['info']['N']) if r['info']['C_eta'] == 0][-2:]
    if top:
        fig, ax = plt.subplots(1, len(top), figsize=(4.4*len(top), 4))
        ax = np.atleast_1d(ax)
        for a, name in zip(ax, top):
            r = runs[name]
            N = r['info']['N']
            recs = plane_records(os.path.join(root, name))
            tr, step, path, rec = min(recs, key=lambda x: abs(x[0] - 3.0))
            f = plane_at(gen_ic.Grid(N, 2, N), path, rec)
            im = a.imshow(np.log10(np.maximum(f['zeta'], 1e-12)).T, origin='lower', cmap='Greys',
                          extent=(0, 2*np.pi, 0, 2*np.pi), vmin=np.log10(f['zeta'].max()) - 4)
            a.set_title('log10 zeta, N = %d, t = %.2f' % (N, tr))
            a.set_xlabel('x')
            a.set_ylabel('z')
            fig.colorbar(im, ax=a, shrink=0.8)
        fig.tight_layout()
        fig.savefig(os.path.join(root, 'fig_ot_zeta.png'), dpi=110)
        plt.close(fig)


OT_LATE_TIMES = (3.0, np.pi, 5.0, 6.0)


def ot_late(root, runs):
    """J_y = d_z B_x - d_x B_z of the runs with C_eta = 0 on the plane y0 at OT_LATE_TIMES: max|J| and the range
    of rho for each (N, t), and the maps (fig_ot_late_J.png)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rs = sorted(((r['info']['N'], name) for name, r in runs.items() if r['info']['C_eta'] == 0 and r.get('complete')))
    if not rs:
        return {}
    out = {}
    fig, ax = plt.subplots(len(OT_LATE_TIMES), len(rs), figsize=(3.3*len(rs), 3.1*len(OT_LATE_TIMES)), squeeze=False)
    for j, (N, name) in enumerate(rs):
        recs = plane_records(os.path.join(root, name))
        k = np.fft.fftfreq(N, 1.0/N)
        KX, KZ = np.meshgrid(k, k, indexing='ij')
        for i, t in enumerate(OT_LATE_TIMES):
            tr, step, path, rec = min(recs, key=lambda x: abs(x[0] - t))
            if abs(tr - t) > 0.01:  # no plane at t (planes every 0.25; N = 1024 ends at t = pi)
                ax[i, j].axis('off')
                continue
            f = plane_at(gen_ic.Grid(N, 2, N), path, rec)
            J = np.real(np.fft.ifft2(1j*KZ*np.fft.fft2(f['bx']) - 1j*KX*np.fft.fft2(f['bz'])))
            out['%s %.4f' % (name, tr)] = dict(N=N, t=float(tr), J_max=float(np.abs(J).max()),
                                               rho_min=float(f['rho'].min()), rho_max=float(f['rho'].max()))
            im = ax[i, j].imshow(J.T, origin='lower', cmap='RdBu_r', vmin=-40, vmax=40,
                                 extent=(0, 2*np.pi, 0, 2*np.pi))
            ax[i, j].set_title('J_y, N = %d, t = %.2f' % (N, tr), fontsize=9)
            ax[i, j].set_xticks([])
            ax[i, j].set_yticks([])
    fig.colorbar(im, ax=ax, shrink=0.6)
    fig.savefig(os.path.join(root, 'fig_ot_late_J.png'), dpi=100)
    plt.close(fig)
    return out


def ot(root):
    cases, runs = {}, {}
    for f in sorted(glob.glob(os.path.join(root, '*', 'case.json'))):
        case = os.path.dirname(f)
        name = os.path.basename(case)
        m = case_info(case)
        cases[name] = dict(status=status(case), role=m['role'], superseded_by=superseded(case))
        if name.startswith('ot_') and cases[name]['status'] == 'exit 0' and not cases[name]['superseded_by']:
            runs[name] = ot_run(case, m)
    checks = {}
    for name, r in runs.items():
        st = r['stats']
        res = max(v for v in (st['R_zeta'], st['R_svv_A'], st['R_eta'], 0.0) if v is not None)
        checks['run ' + name] = dict(complete=r['complete'], R=res, drift_M=st['drift_M'],
                                     ok=bool(r['complete'] and res <= TOL_RES))
    for name, c in cases.items():
        if name.startswith('ot_') and c['status'] != 'exit 0' and not c['superseded_by']:
            checks['run ' + name] = dict(ok=False, status=c['status'])
    conv = ot_convergence(runs)
    dec = ot_decomp(root, cases)
    for name, d in dec.items():
        if 'ok10' in d:
            checks['decomp ' + name] = dict(ok=d['ok10'], max=d['n10'][max(d['n10'], key=int)]['max'])
    for name, c in conv.items():
        if 'ok' in c:
            checks['conv ' + name] = dict(ok=c['ok'])
    for name, r in sorted(runs.items(), key=lambda x: (x[1]['info']['C_eta'], x[1]['info']['N'])):
        st, sh = r['stats'], r['share']
        print('%-14s complete=%-5s dM %.1e R %.1e dt/lim %.2f dt*lam %.2f rho [%.4f, %.4f] divB %.1e '
              'eps_peak %.4f at %.3f share visc %.3f svv %.3f svv_b %.3f eta %.3f ideal %.3f err_t %.1e'
              % (name, r['complete'], st['drift_M'], checks['run ' + name]['R'], st['max_dt_over_limit'],
                 st['max_dt_lam_svv'], st['rho_min'], st['rho_max'], st['divb_max'], r['eps_peak']['eps'],
                 r['eps_peak']['t'], sh['visc'], sh['svv'], sh['svv_b'], sh['eta'], sh['ideal'], st['time_err']))
        print('   zeta frac %s' % {k: round(z['frac'], 5) for k, z in r['zeta'].items()})
        print('   a90 %s' % {k: (round(z['a90'], 5) if z.get('a90') else None) for k, z in r['zeta'].items()})
    for name, c in conv.items():
        print('conv %s: Ns %s pairs %s' % (name, c['Ns'], c['pairs']))
        for k in ('state', 'eps_peak', 't_peak', 'share_visc', 'keff_ratio', 'zfrac_last', 'a90_last', 'ok'):
            if k in c:
                print('   %s %s' % (k, c[k]))
        for t, v in c.get('state_t', {}).items():
            print('   state at t = %s: %s' % (t, ' '.join('%s %.1e' % kv for kv in v.items())))
    for name, d in dec.items():
        print('decomp %s: n10 %s' % (name, {k: '%.1e' % v['max'] for k, v in d.get('n10', {}).items()}))
        if 'cons' in d:
            print('   cons %s' % {t: '%.1e' % max(v.values()) for t, v in d['cons'].items()})
            print('   planes %s' % ['%.2f: %.1e' % (p['t'], p['max']) for p in d['planes']])
    for k, c in checks.items():
        if not c['ok']:
            print('%-30s FAIL %s' % (k, c))
    late = {}
    if runs:
        ot_figures(root, runs)
        late = ot_late(root, runs)
        for k, v in late.items():
            print('late %-26s max|J| %.1f rho [%.3f, %.3f]' % (k, v['J_max'], v['rho_min'], v['rho_max']))
    ok_all = bool(checks) and all(c['ok'] for c in checks.values())
    print('ot:', 'PASS' if ok_all else 'FAIL')
    return dict(cases=cases, runs=runs, checks=checks, conv=conv, decomp=dec, late=late, ok=ok_all)


# ------------------------------------------------------------ S4.4 decaying turbulence
TURB_TAU = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)  # t/tau0 of the spectra (out3d) and zeta statistics (planes)
TURB_DEC_STEPS = (10, 100, 1000)  # the long decomposition runs against turb4_N64 at these steps
TURB_PLANES = ('y0', 'z0')        # y0 holds B0 (along z), z0 is normal to it


def out3d_records(case):
    """[(t, out3d directory, record)] over all runs of the case (the end of every chunk adds a record, and a
    continued run repeats the state at its start)."""
    out = []
    for d in run_dirs(case):
        fn = os.path.join(d, 'out3d', 'time.dat')
        if os.path.exists(fn):
            for rec, t in enumerate(np.atleast_2d(np.loadtxt(fn))[:, 0]):
                out.append((float(t), os.path.join(d, 'out3d'), rec))
    return out


def out3d_at(g, path, rec):
    """dict name -> Fourier array (nkx, nkz, nky) of one record of out3d/ (memory map; the layout with the gaps
    of mpiio_write_var as gen_ic._read_out3d_k)."""
    n = g.nky*g.nkz*g.nkx
    nrec = len(np.atleast_2d(np.loadtxt(os.path.join(path, 'time.dat'))))
    out = {}
    for q in ('rho', 'mx', 'my', 'mz', 'bx', 'by', 'bz'):
        a = np.memmap(os.path.join(path, q + '.dat'), dtype='<c16', mode='r')
        assert len(a) in (nrec*n, (2*nrec - 1)*n), '%s: %d values for %d records' % (path, len(a), nrec)
        j = rec*(2 if len(a) == (2*nrec - 1)*n else 1)
        out[q] = np.array(a[j*n:(j + 1)*n]).reshape(g.nky, g.nkz, g.nkx).transpose(2, 1, 0)
    return out


def shell_spectra_3d(g, fk):
    """Shell-summed spectra of a record of out3d_at as shell_spectra (bins |k| in [n - 1/2, n + 1/2), L = 2 pi;
    the modes with 0 < k_y < Nyquist count twice, r2c), and s_out: the fraction of the variance of
    s = ln rho outside the 2/3 mask (s itself is not masked)."""
    K = np.sqrt(g.kx**2 + g.ky**2 + g.kz**2)
    n = np.rint(K).astype(int)
    w = np.broadcast_to(np.where((g.iky == 0) | (g.iky == g.nly//2), 1.0, 2.0).reshape(1, 1, -1), K.shape)
    nb = g.nlx//2 + 1
    sel = n < nb

    def shell(p):
        return np.bincount(n[sel], weights=(w*p)[sel], minlength=nb)[:nb]

    def e(vk):
        return shell(sum(0.5*np.abs(v)**2 for v in vk))

    rho = g.ifft(fk['rho'])
    sr = np.sqrt(rho)
    m = [g.ifft(fk[q]) for q in ('mx', 'my', 'mz')]
    b = [g.ifft(fk[q]) for q in ('bx', 'by', 'bz')]
    uk = [g.fft(v/rho) for v in m]
    s = dict(E_u=e(uk), E_sru=e([g.fft(v/sr) for v in m]), E_B=e([fk[q] for q in ('bx', 'by', 'bz')]))
    s['E_tot'] = s['E_sru'] + s['E_B']
    s['P_drho'] = shell(np.abs(fk['rho'])**2*(K > 0))
    sk = g.fft(np.log(rho))
    p = np.abs(sk)**2*(K > 0)
    s['P_lnrho'] = shell(p)
    s_out = float(np.sum(w*p*(1 - g.mask))/np.sum(w*p)) if np.sum(p) > 0 else 0.0
    with np.errstate(invalid='ignore', divide='ignore'):
        comp = np.abs(g.kx*uk[0] + g.ky*uk[1] + g.kz*uk[2])**2/np.where(K > 0, K, 1)**2
        s['chi'] = shell(0.5*comp)/s['E_u']
    s['E_zp'] = e([g.fft(v/r + w_/sr) for v, r, w_ in zip(m, (rho,)*3, b)])
    s['E_zm'] = e([g.fft(v/r - w_/sr) for v, r, w_ in zip(m, (rho,)*3, b)])
    return {q: s[q] for q in SPEC_Q}, s_out


def turb_run(case, m):
    """Time series, budget, dt limits, spectra (out3d), s_out and zeta statistics (planes y0, z0) of one
    turbulence run (also a run that stopped: complete = False)."""
    import analyze_stage3 as ast3
    N, tau0 = m['N'], m['tau0']
    g = gen_ic.Grid(N, N, N)
    r = dict(info={k: m.get(k) for k in ('role', 'problem', 'Ms', 'N', 'nstep', 'nrun', 'dt', 'nproc', 'tau0',
                                         'sigma_c', 'H_m')})
    cons = np.concatenate(load_runs(case, 'conservation.dat'))
    bud = np.concatenate(load_runs(case, 'budget.dat'))
    cfl = np.concatenate(load_runs(case, 'cfl.dat'))
    r['t_last'] = float(cons[-1, 0])
    if m.get('adaptive') or m.get('to_t_end'):
        # dt changed by reset_dt: the chain ends within dt/2 of t_end (file t_end), as in the job script
        te = float(open(os.path.join(case, 't_end')).read())
        r['complete'] = abs(r['t_last'] - te) <= 0.5*cfl[-1, 1] + 1e-12*te
    else:
        te = m['nrun']*m['dt']
        r['complete'] = abs(r['t_last'] - te) <= 1e-9*te
    r['t_end'] = te
    # largest dt of the run, for the match of output records to the times a tau0 ('track' raises dt above dt0)
    dtmax = max(m['dt'], float(cfl[:, 1].max()))
    st = ast3.stats(case, dict(dir='x', case=dict(N=N)), [])
    st['R_svv_A'] = svv_res_rel(case)
    r['stats'] = {k: st.get(k) for k in ('drift_M', 'drift_Px', 'drift_Py', 'drift_Pz', 'rho_min', 'rho_max',
                                         'zeta_max', 'R_zeta', 'R_svv', 'R_svv_A', 'R_eta', 'max_dt_over_limit',
                                         'max_dt_lam_svv', 'dE', 'n_E_up', 'time_err', 'int_rk_eps',
                                         'int_rk_ideal_eps', 'int_rk_visc_eps', 'int_rk_svv_eps',
                                         'int_rk_svv_b_eps', 'int_rk_eta_eps')}
    r['stats']['divb_max'] = float(cons[:, 19].max())
    tot = st.get('int_rk_eps')
    r['share'] = {k: (st['int_rk_%s_eps' % k]/tot if tot else None) for k in ('ideal', 'visc', 'svv', 'svv_b', 'eta')}
    last = cons[-1]
    r['end'] = dict(t=float(last[0]), E_K=float(last[7]), E_M=float(last[8]), E=float(last[10]),
                    rho_min=float(last[11]), rho_max=float(last[12]))
    # the limit of dt (min of dt_adv and dt_visc of cfl.dat), where it is reached, and its parts
    lim = np.minimum(cfl[:, 2], cfl[:, 3])
    i = int(np.argmin(lim))
    r['dt_lim'] = dict(min=float(lim[i]), t=float(cfl[i, 0]), adv=float(cfl[i, 2]), visc=float(cfl[i, 3]),
                       lam_svv=float(cfl[i, 12]), max_dt_lam=float(cfl[:, 14].max()))
    eps = -bud[:, 6]
    i = int(np.argmax(eps))
    r['eps_peak'] = dict(t=float(bud[i, 0] + 0.5*bud[i, 1]), eps=float(eps[i]))
    sb, sc, sf = max(1, len(bud)//2000), max(1, len(cons)//2000), max(1, len(cfl)//2000)
    r['series'] = dict(t=(bud[::sb, 0] + 0.5*bud[::sb, 1]).tolist(), eps=eps[::sb].tolist(),
                       visc=(-bud[::sb, 10]).tolist(), svv=(-bud[::sb, 20]).tolist(), eta=(-bud[::sb, 26]).tolist(),
                       ideal=(-bud[::sb, 9]).tolist(), tc=cons[::sc, 0].tolist(), E_K=cons[::sc, 7].tolist(),
                       E_M=cons[::sc, 8].tolist(), rho_min=cons[::sc, 11].tolist(),
                       rho_max=cons[::sc, 12].tolist(), zeta_max=cons[::sc, 13].tolist(),
                       tf=cfl[::sf, 0].tolist(), dt_adv=cfl[::sf, 2].tolist(), dt_visc=cfl[::sf, 3].tolist(),
                       dt_lam=cfl[::sf, 14].tolist())
    r['spec'], r['s_out'], r['zeta'] = {}, {}, {}
    recs = out3d_records(case)
    for a in TURB_TAU:
        t = a*tau0
        if not recs or t > r['t_last'] + 1e-9:
            continue
        tr, path, rec = min(recs, key=lambda x: abs(x[0] - t))
        if abs(tr - t) > 2*dtmax:
            continue
        sp, so = shell_spectra_3d(g, out3d_at(g, path, rec))
        key = '%.2f' % a
        r['spec'][key] = dict(t=tr, **{q: v.tolist() for q, v in sp.items()})
        r['s_out'][key] = so
    for pl in TURB_PLANES:
        recs = plane_records(case, pl)
        r['zeta'][pl] = {}
        for a in TURB_TAU:
            t = a*tau0
            if not recs or t > r['t_last'] + 1e-9:
                continue
            tr, step, path, rec = min(recs, key=lambda x: abs(x[0] - t))
            if abs(tr - t) <= 2*dtmax:
                r['zeta'][pl]['%.2f' % a] = dict(t=tr, **zeta_stats(plane_at(g, path, rec, pl)))
    return r


def series_at(r, q, t):
    s = r['series']
    return float(np.interp(t, s['tc'], s[q]))


def turb_convergence(runs):
    """Each problem: the checks of S4.3 on the top two N (complete runs), the state at t/tau0 = 1, 2, 3 (judged
    at t_end), k_eff, zeta fraction and a90 ratios (planes y0, z0) for each pair (N, 2N)."""
    conv = {}
    for prob in sorted({r['info']['problem'] for r in runs.values()}):
        rs = sorted((r for r in runs.values() if r['info']['problem'] == prob and r['complete']),
                    key=lambda r: r['info']['N'])
        c = dict(Ns=[r['info']['N'] for r in rs], pairs={})
        for ra, rb in zip(rs, rs[1:]):
            p = dict(keff={}, zfrac={}, a90={})
            for key, sa in ra['spec'].items():
                if key in rb['spec'] and key != '0.00':
                    p['keff'][key] = k_eff(sa['E_tot'], rb['spec'][key]['E_tot'])
            for pl in TURB_PLANES:
                for key, za in ra['zeta'][pl].items():
                    zb = rb['zeta'][pl].get(key)
                    if zb and za['frac'] > 0:
                        p['zfrac'][pl + ' ' + key] = zb['frac']/za['frac']
                    if zb and za.get('a90') and zb.get('a90'):
                        p['a90'][pl + ' ' + key] = zb['a90']/za['a90']
            c['pairs']['%d-%d' % (ra['info']['N'], rb['info']['N'])] = p
        if len(rs) >= 2:
            a, b = rs[-2], rs[-1]
            tau0 = b['info']['tau0']
            c['state_t'] = {'%g' % j: {q: rel(series_at(a, q, j*tau0), series_at(b, q, j*tau0))
                                       for q in ('rho_min', 'rho_max', 'E_K', 'E_M')} for j in (1, 2)}
            c['state'] = {q: rel(a['end'][q], b['end'][q]) for q in ('rho_min', 'rho_max', 'E_K', 'E_M')}
            c['eps_peak'] = rel(a['eps_peak']['eps'], b['eps_peak']['eps'])
            c['t_peak'] = rel(a['eps_peak']['t'], b['eps_peak']['t'])
            c['share_visc'] = rel(a['share']['visc'], b['share']['visc']) if b['share']['visc'] else None
            pr = list(c['pairs'].values())
            ok = (max(c['state'].values()) <= OT_CONV_STATE and c['eps_peak'] <= OT_CONV_EPS
                  and c['t_peak'] <= OT_CONV_TPEAK and c['share_visc'] is not None
                  and c['share_visc'] <= OT_CONV_SHARE)
            if len(pr) >= 2:
                c['keff_ratio'] = {k: pr[-1]['keff'][k]/pr[-2]['keff'][k] for k in pr[-1]['keff']
                                   if pr[-2]['keff'].get(k)}
                ok &= all(OT_KEFF_RATIO[0] <= v <= OT_KEFF_RATIO[1] for v in c['keff_ratio'].values())
            c['a90_last'] = pr[-1]['a90']
            ok &= bool(pr[-1]['a90']) and all(OT_ZFRAC_RATIO[0] <= v <= OT_ZFRAC_RATIO[1]
                                              for v in pr[-1]['a90'].values())
            c['ok'] = bool(ok)
        conv[prob] = c
    return conv


def turb_decomp(root, cases):
    """Planes y0, z0 after 10 steps against P = 64 (TOL_DEC10), and conservation.dat of the runs over
    TURB_DEC_STEPS steps against turb4_N64 (64 ranks; both every step): the largest difference up to each step
    relative to the largest |value| of the column."""
    out = {}
    g = gen_ic.Grid(64, 64, 64)
    ref10 = os.path.join(root, 'dec_N64_np64_n10')
    for P in (1, 8):
        name = 'dec_N64_np%d' % P
        d = {}
        if cases.get(name + '_n10', {}).get('status') == 'exit 0' and \
                cases.get('dec_N64_np64_n10', {}).get('status') == 'exit 0':
            d['n10'] = {}
            for pl in TURB_PLANES:
                ra, rb = plane_records(os.path.join(root, name + '_n10'), pl), plane_records(ref10, pl)
                for (ta, sa, pa, ia), (tb, sb, pb, ib) in zip(ra, rb):
                    diff = plane_diff(plane_at(g, pa, ia, pl), plane_at(g, pb, ib, pl), 2*np.pi/64)
                    d['n10']['%s %d' % (pl, sb)] = dict(step=(sa, sb), **diff)
            d['ok10'] = (all(v['step'][0] == v['step'][1] for v in d['n10'].values())
                         and {v['step'][1] for v in d['n10'].values()} == {0, 10}
                         and max(v['max'] for v in d['n10'].values()) <= TOL_DEC10)
        ref = os.path.join(root, 'turb4_N64')
        if cases.get(name, {}).get('status') == 'exit 0' and os.path.exists(os.path.join(ref, 'conservation.dat')):
            a = np.concatenate(load_runs(os.path.join(root, name), 'conservation.dat'))
            b = np.concatenate(load_runs(ref, 'conservation.dat'))
            n = min(len(a), len(b))
            assert np.allclose(a[:n, 0], b[:n, 0], rtol=0, atol=1e-12), name
            cols = dict(E_K=7, E_M=8, E=10, rho_min=11, rho_max=12, zeta_max=13)
            d['cons'] = {}
            for ns in TURB_DEC_STEPS:
                j = min(n, ns + 1)   # row 0: t = 0
                d['cons'][str(ns)] = {q: float(np.max(np.abs(a[:j, c] - b[:j, c]))
                                               / max(float(np.max(np.abs(b[:j, c]))), 1e-300))
                                      for q, c in cols.items()}
        out[name] = d
    return out


def turb_figures(root, runs):
    """eps(t) with its shock-viscous part, E_tot(k) at t = tau0 and 3 tau0, and log10 zeta on the plane z0 at
    t = tau0 (largest N of each problem) as PNG files in root."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import gen_stage4 as gs
    probs = sorted({r['info']['problem'] for r in runs.values()})
    Ns = sorted({r['info']['N'] for r in runs.values()})
    shade = {N: plt.cm.Blues(0.35 + 0.65*i/max(1, len(Ns) - 1)) for i, N in enumerate(Ns)}
    order = sorted(runs.items(), key=lambda x: x[1]['info']['N'])
    fig, ax = plt.subplots(1, len(probs), figsize=(4.6*len(probs), 3.4), squeeze=False)
    for j, prob in enumerate(probs):
        a = ax[0, j]
        for name, r in order:
            if r['info']['problem'] != prob:
                continue
            tau0, N = r['info']['tau0'], r['info']['N']
            t = np.array(r['series']['t'])/tau0
            a.plot(t, r['series']['eps'], color=shade[N], lw=1.5, label='N = %d' % N)
            a.plot(t, r['series']['visc'], color=shade[N], lw=1.0, ls='--')
        a.set_title('%s: eps (solid), shock viscosity (dashed)' % gs.turb_label(prob), fontsize=9)
        a.set_xlabel('t/tau0')
        a.legend(frameon=False, fontsize=8)
        a.grid(alpha=0.3, lw=0.5)
    ax[0, 0].set_ylabel('eps = -dE/dt')
    fig.tight_layout()
    fig.savefig(os.path.join(root, 'fig_turb_eps.png'), dpi=110)
    plt.close(fig)
    # only the problems with a record at t = tau0 (the runs of M_s = 4 stopped before it)
    sprobs = [p for p in probs if any(r['info']['problem'] == p and r['spec'].get('1.00') for r in runs.values())]
    fig, ax = plt.subplots(len(sprobs), 2, figsize=(9, 3.4*len(sprobs)), squeeze=False, sharey='row')
    for i, prob in enumerate(sprobs):
        for j, key in enumerate(('1.00', '3.00')):
            for name, r in order:
                s = r['spec'].get(key)
                if r['info']['problem'] != prob or not s:
                    continue
                N = r['info']['N']
                k = np.arange(len(s['E_tot']))
                ax[i, j].loglog(k[1:], np.array(s['E_tot'][1:]), color=shade[N], lw=1.5, label='N = %d' % N)
                ax[i, j].axvline(N/3, color=shade[N], lw=0.6, ls=':')
            ax[i, j].set_title('E_tot(k), %s, t = %s tau0' % (gs.turb_label(prob), key[0]), fontsize=9)
            ax[i, j].set_xlabel('k')
            ax[i, j].legend(frameon=False, fontsize=8)
            ax[i, j].grid(alpha=0.3, lw=0.5, which='both')
    fig.tight_layout()
    fig.savefig(os.path.join(root, 'fig_turb_spec.png'), dpi=110)
    plt.close(fig)
    top = [max((x for x in order if x[1]['info']['problem'] == p), key=lambda x: x[1]['info']['N'])[0]
           for p in probs]
    fig, ax = plt.subplots(1, len(top), figsize=(4.4*len(top), 4), squeeze=False)
    for a, name in zip(ax[0], top):
        r = runs[name]
        N, tau0 = r['info']['N'], r['info']['tau0']
        recs = plane_records(os.path.join(root, name), 'z0')
        if not recs:
            continue
        tr, step, path, rec = min(recs, key=lambda x: abs(x[0] - tau0))
        f = plane_at(gen_ic.Grid(N, N, N), path, rec, 'z0')
        im = a.imshow(np.log10(np.maximum(f['zeta'], 1e-12)).T, origin='lower', cmap='Greys',
                      extent=(0, 2*np.pi, 0, 2*np.pi), vmin=np.log10(f['zeta'].max()) - 4)
        a.set_title('log10 zeta (z = 0), %s N = %d, t = %.2f tau0' % (r['info']['problem'], N, tr/tau0), fontsize=9)
        a.set_xlabel('x')
        a.set_ylabel('y')
        fig.colorbar(im, ax=a, shrink=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(root, 'fig_turb_zeta.png'), dpi=110)
    plt.close(fig)


def turb_print(name, status, r, superseded_by=None):
    """One run of turb_run (with its residual r['R']): stability, drifts, dissipation, dt limit, s_out, zeta."""
    st, sh, dl = r['stats'], r['share'], r['dt_lim']
    print('%-12s %-12s t %.4f/%.4f dM %.1e dP %.1e R %.1e dt/lim %.2f dt*lam %.2f rho [%.4f, %.4f] divB %.1e '
          'eps_peak %.4f at %.3f share visc %.3f svv %.3f svv_b %.3f ideal %.3f err_t %.1e'
          % (name, status, r['t_last'], r['t_end'],
             st['drift_M'], max(st['drift_Px'], st['drift_Py'], st['drift_Pz']), r['R'],
             st['max_dt_over_limit'], st['max_dt_lam_svv'], st['rho_min'], st['rho_max'], st['divb_max'],
             r['eps_peak']['eps'], r['eps_peak']['t'], sh['visc'], sh['svv'], sh['svv_b'], sh['ideal'],
             st['time_err']))
    if superseded_by:
        print('   superseded by %s (not judged)' % superseded_by)
    print('   dt limit: min %.3e at t = %.4f (adv %.3e visc %.3e lam_svv %.3e), max dt*lam %.2f'
          % (dl['min'], dl['t'], dl['adv'], dl['visc'], dl['lam_svv'], dl['max_dt_lam']))
    print('   s_out %s' % {k: '%.1e' % v for k, v in r['s_out'].items()})
    for pl in TURB_PLANES:
        print('   %s zeta frac %s a90 %s' % (pl, {k: round(z['frac'], 4) for k, z in r['zeta'][pl].items()},
                                            {k: (round(z['a90'], 5) if z.get('a90') else None)
                                             for k, z in r['zeta'][pl].items()}))


def turb_print_conv(conv):
    for name, c in conv.items():
        print('conv %s: Ns %s pairs %s' % (name, c['Ns'], c['pairs']))
        for k in ('state_t', 'state', 'eps_peak', 't_peak', 'share_visc', 'keff_ratio', 'a90_last', 'ok'):
            if k in c:
                print('   %s %s' % (k, c[k]))


def turb_residual(r):
    """Largest residual of the budget (R_zeta, R_svv_A, R_eta) of a turb_run result."""
    st = r['stats']
    return max(v for v in (st['R_zeta'], st['R_svv_A'], st['R_eta'], 0.0) if v is not None)


def turb(root):
    cases, runs = {}, {}
    for f in sorted(glob.glob(os.path.join(root, '*', 'case.json'))):
        case = os.path.dirname(f)
        name = os.path.basename(case)
        m = case_info(case)
        cases[name] = dict(status=status(case), role=m['role'], superseded_by=superseded(case))
        if name.startswith('turb') and os.path.exists(os.path.join(case, 'conservation.dat')):
            runs[name] = turb_run(case, m)
    checks = {}
    for name, r in runs.items():
        st = r['stats']
        r['R'] = res = turb_residual(r)
        if cases[name]['superseded_by']:
            continue
        checks['run ' + name] = dict(complete=r['complete'], status=cases[name]['status'], R=res,
                                     drift_M=st['drift_M'], ok=bool(r['complete'] and res <= TOL_RES))
    for name, c in cases.items():
        if name.startswith('turb') and name not in runs and not c['superseded_by']:
            checks['run ' + name] = dict(ok=False, status=c['status'])
    judged = {k: r for k, r in runs.items() if not cases[k]['superseded_by']}
    conv = turb_convergence(judged)
    dec = turb_decomp(root, cases)
    for name, d in dec.items():
        if 'ok10' in d:
            checks['decomp ' + name] = dict(ok=d['ok10'], max=max(v['max'] for v in d['n10'].values()))
    for name, c in conv.items():
        if 'ok' in c:
            checks['conv ' + name] = dict(ok=c['ok'])
    for name, r in sorted(runs.items(), key=lambda x: (x[1]['info']['problem'], x[1]['info']['N'])):
        turb_print(name, cases[name]['status'], r, cases[name]['superseded_by'])
    turb_print_conv(conv)
    for name, d in dec.items():
        print('decomp %s: n10 %s' % (name, {k: '%.1e' % v['max'] for k, v in d.get('n10', {}).items()}))
        if 'cons' in d:
            print('   cons %s' % {s: '%.1e' % max(v.values()) for s, v in d['cons'].items()})
    for k, c in checks.items():
        if not c['ok']:
            print('%-30s FAIL %s' % (k, c))
    if judged:
        turb_figures(root, judged)
    ok_all = bool(checks) and all(c['ok'] for c in checks.values())
    print('turb:', 'PASS' if ok_all else 'FAIL')
    return dict(cases=cases, runs=runs, checks=checks, conv=conv, decomp=dec, ok=ok_all)


# ------------------------------------------------------------ S4.5 parallel performance
SCALE_PARTS = ('step (solve)', 'div u, zeta, eta (halo)', 'get_svv_es', 'svv_power', 'stage_budget')
SCALE_MAIN = ('FFT', 'nonlinear terms')  # timers of proc0 in calliope.out.std [min]


def std_timers(case):
    """Timers of the .std of a run: the parts of MHD_COMP_ISOTH (max over ranks [s]) and FFT, nonlinear terms
    (proc0 [s]); None if the run has no timer block (not finished)."""
    fs = [f for d in run_dirs(case) for f in glob.glob(os.path.join(d, '*.std'))]
    if len(fs) != 1:
        return None
    txt = open(fs[0], errors='replace').read()
    out = {}
    for name in SCALE_PARTS:
        m = re.search(r'^\s+%s\s+(\S+)\s+(\S+)\s+\S+ %%' % re.escape(name), txt, re.M)
        if m:
            out[name] = float(m.group(2))
    for name in SCALE_MAIN:
        m = re.search(r'^\s+%s\s+(\S+) min' % re.escape(name), txt, re.M)
        if m:
            out[name] = 60*float(m.group(1))
    return out if len(out) == len(SCALE_PARTS) + len(SCALE_MAIN) else None


def scale(root):
    """Cost per step (T(n = 110) - T(n = 10))/100 of each part, strong-scaling efficiency
    T(P0) P0/(T(P) P) on the fixed grid and weak-scaling efficiency T(P0)/T(P) at fixed points per rank (with
    and without the factor log2(n)/log2(n0) of the FFT)."""
    import gen_stage4 as gs
    cases, t = {}, {}
    for f in sorted(glob.glob(os.path.join(root, '*', 'case.json'))):
        case = os.path.dirname(f)
        m = case_info(case)
        cases[os.path.basename(case)] = st = status(case)
        key = (m['problem'], m['N'], m['nproc'])
        if st == 'exit 0' and os.path.isdir(os.path.join(case, 'chunk_1')):
            # continued from restart/: the timers of calliope.out.std cover only the last chunk
            cases[os.path.basename(case)] = 'exit 0 after a restart'
        elif st == 'exit 0':
            t.setdefault(key, {})[m['nrun']] = std_timers(case)
    cost = {}
    for key, d in t.items():
        a, b = d.get(gs.SCALE_NSTEP[0]), d.get(gs.SCALE_NSTEP[1])
        if a and b:
            n = gs.SCALE_NSTEP[1] - gs.SCALE_NSTEP[0]
            cost[key] = {q: (b[q] - a[q])/n for q in b}
    out = dict(cases=cases, strong={}, weak={})
    for prob, tab in (('ot', gs.SCALE_2D), ('turb4', gs.SCALE_3D)):
        N, Ps = tab['strong']
        rows = [(P, cost[(prob, N, P)]) for P in Ps if (prob, N, P) in cost]
        if rows:
            P0, c0 = rows[0]
            out['strong'][prob] = [dict(P=P, eff=c0['step (solve)']*P0/(c['step (solve)']*P), **c) for P, c in rows]
        pts = sorted(set(tab['weak']) | {(N, P) for P in Ps if (N**(2 if prob == 'ot' else 3))//P ==
                                         (tab['weak'][0][0]**(2 if prob == 'ot' else 3))//tab['weak'][0][1]})
        rows = [(n, P, cost[(prob, n, P)]) for n, P in pts if (prob, n, P) in cost]
        if rows:
            n0, P0, c0 = rows[0]
            dim = 2 if prob == 'ot' else 3
            out['weak'][prob] = [dict(N=n, P=P, eff=c0['step (solve)']/c['step (solve)'],
                                      eff_log=c0['step (solve)']/c['step (solve)']*math.log2(n**dim)/math.log2(n0**dim),
                                      **c) for n, P, c in rows]
    for kind in ('strong', 'weak'):
        for prob, rows in out[kind].items():
            print('%s %s' % (kind, prob))
            for r in rows:
                print('   N %5d P %5d  step %.4e s  eff %.3f%s  ' % (r.get('N', gs.SCALE_2D['strong'][0] if prob == 'ot'
                                                                         else gs.SCALE_3D['strong'][0]), r['P'],
                                                                  r['step (solve)'], r['eff'],
                                                                  ' (log %.3f)' % r['eff_log'] if 'eff_log' in r else '')
                      + ' '.join('%s %.1f%%' % (q, 100*r[q]/r['step (solve)']) for q in SCALE_PARTS[1:] + SCALE_MAIN))
    bad = {k: v for k, v in cases.items() if v != 'exit 0'}
    print('scale: %d cases, not exit 0: %s' % (len(cases), bad or 'none'))
    if out['strong'] or out['weak']:
        scale_figure(root, out)
    return out


def scale_figure(root, out):
    """Cost per step against P (strong, weak) with the ideal lines, as fig_scale.png in root."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    style = dict(ot=('C0', 'o', 'OT 2D'), turb4=('C1', 's', 'turbulence 3D'))
    for j, kind in enumerate(('strong', 'weak')):
        for prob, rows in out[kind].items():
            c, mk, lab = style[prob]
            P = np.array([r['P'] for r in rows], float)
            T = np.array([r['step (solve)'] for r in rows])
            ax[j].loglog(P, T, color=c, marker=mk, ms=8, lw=2, label=lab)
            ideal = T[0]*P[0]/P if kind == 'strong' else np.full_like(P, T[0])
            ax[j].loglog(P, ideal, color=c, lw=1, ls='--')
        ax[j].set_xlabel('ranks P')
        ax[j].set_title('%s scaling (dashed: ideal)' % kind, fontsize=9)
        ax[j].grid(alpha=0.3, lw=0.5, which='both')
        ax[j].legend(frameon=False, fontsize=8)
    ax[0].set_ylabel('wall time per step [s]')
    fig.tight_layout()
    fig.savefig(os.path.join(root, 'fig_scale.png'), dpi=110)
    plt.close(fig)


BATCHES = {'readback': readback, 'regress': regress, 'alfven': alfven, 'cpdiag': cpdiag, 'planar': planar,
           'ot': ot, 'turb': turb, 'scale': scale}


if __name__ == '__main__':
    batch = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.join(paths.RUNS, 'stage4')
    root = os.path.join(base, batch)
    res = BATCHES[batch](root)
    with open(os.path.join(root, 'summary.json'), 'w') as f:
        json.dump(res, f, indent=1, default=str)
