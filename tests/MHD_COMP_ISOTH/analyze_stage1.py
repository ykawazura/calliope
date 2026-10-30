#!/usr/bin/env python3
"""Analysis of Stage 1 runs (1D isothermal Euler, B = 0) of MHD_COMP_ISOTH.

    python3 analyze_stage1.py [--plot] case_dir [case_dir ...]

Reads conservation.dat, budget.dat, cfl.dat and out1d/line_*.dat of each case
and prints one summary block per case (also saved as <case>/summary.json).
For init_type = 'riemann' the final line is compared with the exact solution
(riemann_exact.py). Only numpy is required (matplotlib for --plot).
"""
import sys
import os
import re
import glob
import io
import json
import math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import riemann_exact as rx


# ---------------------------------------------------------------- input files
def read_nml(path):
    """Flat dict {name: value-string} of every namelist entry in calliope.in."""
    d = {}
    for line in open(path):
        line = line.split('!')[0].strip()
        m = re.match(r'(\w+)\s*=\s*(.+)', line)
        if m:
            d[m.group(1).lower()] = m.group(2).strip()
    return d


def num(s):
    return float(s.replace('d', 'e').replace('D', 'e'))


def nums(s):
    return [num(x) for x in s.split()]


def strv(s):
    return s.strip().strip("'\"")


def load(path):
    """Numeric table; also accepts Fortran Ew.d output whose 3-digit exponent drops the 'E'
    (e.g. 1.797693134862315708145+308 written by binaries older than the es26.17e3 format)."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    txt = re.sub(r'(?<=[0-9.])([+-]\d{3})(?=\s|$)', r'e\1', open(path).read())
    a = np.loadtxt(io.StringIO(txt), comments='#', ndmin=2)
    return a if a.size else None


def read_lines(case):
    out = []
    for fn in sorted(glob.glob(os.path.join(case, 'out1d', 'line_*.dat'))):
        with open(fn) as f:
            hdr = f.readline()
        t = float(re.search(r't =\s*([-+0-9.eE]+)', hdr).group(1))
        step = int(re.search(r'step =\s*(\d+)', hdr).group(1))
        out.append((t, step, np.loadtxt(fn, comments='#')))
    return out


# ---------------------------------------------------------------- shocks
def unwrap(x, x0, L):
    return x0 + ((x - x0 + 0.5 * L) % L - 0.5 * L)


def features(t, L, ifaces):
    """Exact positions of all wave edges at time t: list of (x, interface, side, kind)."""
    f = []
    for n, (x0, w) in enumerate(ifaces):
        for side in ('L', 'R'):
            if w[side]['type'] == 'shock':
                f.append(((x0 + w[side]['speed'] * t) % L, n + 1, side, 'shock'))
            else:
                f.append(((x0 + w[side]['head'] * t) % L, n + 1, side, 'head'))
                f.append(((x0 + w[side]['tail'] * t) % L, n + 1, side, 'tail'))
    return f


def shock_metrics(x, rho, xs, rho_up, rho_dn, down_left, W, dl, L):
    """Shock at exact position xs with upstream rho_up, downstream rho_dn (downstream on the
    left if down_left). Window |x - xs| <= W must lie inside the two plateaus.
    Position from the integral of H = (rho - rho_up)/(rho_dn - rho_up) over the window
    (sub-cell accurate, exact for any monotone profile), thickness = jump / max slope,
    over/undershoot relative to the jump."""
    xw = unwrap(x, xs, L)
    sel = np.abs(xw - xs) <= W
    if sel.sum() < 3:
        return None
    xw, rw = xw[sel], rho[sel]
    o = np.argsort(xw)
    xw, rw = xw[o], rw[o]
    H = (rw - rho_up) / (rho_dn - rho_up)
    xa, xb = xw[0] - 0.5 * dl, xw[-1] + 0.5 * dl
    I = np.sum(H) * dl
    pos = xa + I if down_left else xb - I
    jump = abs(rho_dn - rho_up)
    hi, lo = max(rho_up, rho_dn), min(rho_up, rho_dn)
    return dict(x=float(pos), dx_cells=float(unwrap(pos, xs, L) - xs) / dl,
                thickness_cells=float(jump / (np.max(np.abs(np.diff(rw))) / dl) / dl),
                overshoot=float((rw.max() - hi) / jump), undershoot=float((lo - rw.min()) / jump),
                window_cells=float(W / dl))


def shock_window(xs, feats, me, L, dl, halfwin):
    d = [abs(unwrap(xf, xs, L) - xs) for (xf, n, side, kind) in feats if (n, side, kind) != me]
    return min(halfwin * dl, 0.5 * min(d))


def riemann_analysis(nml, lines, out, halfwin=24):
    dirc = strv(nml.get('dir', "'x'"))
    idx = 'xyz'.index(dirc)
    L = 2.0 * math.pi * num(nml['l' + dirc])  # grid.F90: lx is in units of 2 pi
    N = int(nml['nl' + dirc])
    dl = L / N
    cs = math.sqrt(num(nml.get('beta0', '1.0')) / 2.0)
    rho_L, rho_R = num(nml['rho_l']), num(nml['rho_r'])
    u_L, u_R = nums(nml['u_l'])[idx], nums(nml['u_r'])[idx]
    states = [((rho_L, u_L), (rho_R, u_R)), ((rho_R, u_R), (rho_L, u_L))]
    ifaces = [(0.25 * L, rx.waves(rho_L, u_L, rho_R, u_R, cs)),
              (0.75 * L, rx.waves(rho_R, u_R, rho_L, u_L, cs))]
    out['exact'] = dict(interface1=ifaces[0][1], interface2=ifaces[1][1])

    t, step, a = lines[-1]
    x, rho, u = a[:, 0], a[:, 1], a[:, 2 + idx]
    if t <= 0.0:
        out['final'] = None
        return x, rho, u, None, None
    rho_ex, u_ex, valid = rx.periodic_solution(x, t, rho_L, u_L, rho_R, u_R, L, cs)
    out['final'] = dict(t=t, step=step, valid=bool(valid),
                        L1_rho=float(np.mean(np.abs(rho - rho_ex)) * L),
                        L1_u=float(np.mean(np.abs(u - u_ex)) * L),
                        rel_L1_rho=float(np.mean(np.abs(rho - rho_ex)) / np.mean(np.abs(rho_ex))),
                        Linf_rho=float(np.max(np.abs(rho - rho_ex))),
                        TV_rho=float(np.sum(np.abs(np.diff(np.append(rho, rho[0]))))),
                        TV_rho_exact=float(np.sum(np.abs(np.diff(np.append(rho_ex, rho_ex[0]))))),
                        rho_min=float(rho.min()), rho_max=float(rho.max()),
                        rho_min_exact=float(rho_ex.min()))

    # star states: middle half of the region between the L and R waves of each interface
    stars = []
    feats = features(t, L, ifaces)
    for n, (x0, w) in enumerate(ifaces):
        eL = w['L']['speed'] if w['L']['type'] == 'shock' else w['L']['tail']
        eR = w['R']['speed'] if w['R']['type'] == 'shock' else w['R']['tail']
        xa, xb = x0 + eL * t + 0.25 * (eR - eL) * t, x0 + eR * t - 0.25 * (eR - eL) * t
        xu = unwrap(x, 0.5 * (xa + xb), L)
        sel = (xu >= xa) & (xu <= xb)
        st = dict(interface=n + 1, rho_s_exact=w['rho_s'], u_s_exact=w['u_s'], ncells=int(sel.sum()))
        if sel.sum() >= 2:
            st.update(rho_s=float(rho[sel].mean()), rho_s_rel_err=float(rho[sel].mean() / w['rho_s'] - 1.0),
                      rho_s_spread=float((rho[sel].max() - rho[sel].min()) / w['rho_s']),
                      u_s=float(u[sel].mean()), u_s_err=float(u[sel].mean() - w['u_s']))
        stars.append(st)
    out['stars'] = stars

    shocks = []
    for n, (x0, w) in enumerate(ifaces):
        for side, (rup, uup) in (('L', states[n][0]), ('R', states[n][1])):
            if w[side]['type'] != 'shock':
                continue
            S = w[side]['speed']
            down_left = (side == 'R')
            me = (n + 1, side, 'shock')
            xs = (x0 + S * t) % L
            W = shock_window(xs, features(t, L, ifaces), me, L, dl, halfwin)
            m = shock_metrics(x, rho, xs, rup, w['rho_s'], down_left, W, dl, L)
            if m is None:
                continue
            ts, ps = [], []
            for (tn, _, an) in lines:
                if tn <= 0.0:
                    continue
                xsn = (x0 + S * tn) % L
                Wn = shock_window(xsn, features(tn, L, ifaces), me, L, dl, halfwin)
                if Wn < 6 * dl:
                    continue
                mn = shock_metrics(an[:, 0], an[:, 1], xsn, rup, w['rho_s'], down_left, Wn, dl, L)
                if mn is None:
                    continue
                ts.append(tn); ps.append(unwrap(mn['x'], x0 + S * tn, L))
            if len(ts) >= 3:
                p = np.polyfit(ts, ps, 1)
                m['speed_fit'] = float(p[0])
                m['speed_rel_err'] = float(p[0] / S - 1.0)
                m['n_fit'] = len(ts)
            m.update(interface=n + 1, side=side, speed_exact=S, mach=w[side]['mach'],
                     rho2_over_rho1_exact=w['rho_s'] / rup)
            shocks.append(m)
    out['shocks'] = shocks
    return x, rho, u, rho_ex, u_ex


# ---------------------------------------------------------------- smooth waves
def wave1d_setup(nml):
    dirc = strv(nml.get('dir', "'x'"))
    idx = 'xyz'.index(dirc)
    L = 2.0 * math.pi * num(nml['l' + dirc])  # grid.F90: lx is in units of 2 pi
    k = 2.0 * math.pi * num(nml.get('mode', '1')) / L
    cs = math.sqrt(num(nml.get('beta0', '1.0')) / 2.0)
    return idx, L, k, cs, num(nml['amp']), num(nml.get('rho0', '1.0'))


def sound_analysis(nml, lines, out):
    """Linear sound wave s = ln rho0 + A cos(kx) cos(wt), u = c A sin(kx) sin(wt), w = c k.
    Z(t) = a_s(t) + i b_u(t)/c with a_s the cos coefficient of s and b_u the sin coefficient
    of u is A exp(i w t) exactly; amplitude and phase errors are |Z|/A - 1 and arg(Z e^{-iwt}/A)."""
    idx, L, k, cs, A, rho0 = wave1d_setup(nml)
    w = cs * k
    rows = []
    for (t, step, a) in lines:
        x, s, u = a[:, 0], np.log(a[:, 1]), a[:, 2 + idx]
        N = len(x)
        a_s = 2.0 / N * np.sum((s - math.log(rho0)) * np.cos(k * x))
        b_u = 2.0 / N * np.sum(u * np.sin(k * x))
        Z = complex(a_s, b_u / cs) * np.exp(-1j * w * t) / A
        # power outside the fundamental (nonlinear harmonics + noise), relative to A
        sk = np.fft.rfft(s - math.log(rho0)) * 2.0 / N
        m1 = int(round(k * L / (2.0 * math.pi)))
        others = np.sqrt(np.sum(np.abs(sk) ** 2) - np.abs(sk[m1]) ** 2) / A
        rows.append(dict(t=t, amp_err=float(abs(Z) - 1.0), phase_err=float(np.angle(Z)),
                         harmonics=float(others)))
    out['sound'] = dict(omega=w, rows=rows,
                        max_abs_amp_err=max(abs(r['amp_err']) for r in rows),
                        max_abs_phase_err=max(abs(r['phase_err']) for r in rows))


def simple_wave_exact(x, t, A, k, cs, rho0, L):
    """Right-going isothermal simple wave: J- = u - c ln(rho/rho0) = 0, u = U0(xi),
    x = xi + (U0(xi) + c) t, U0 = c A sin(k xi). Valid for t < 1/(c A k)."""
    U0 = lambda xi: cs * A * np.sin(k * xi)
    dU0 = lambda xi: cs * A * k * np.cos(k * xi)
    xi = np.asarray(x, dtype=float) - cs * t
    for _ in range(100):
        F = xi + (U0(xi) + cs) * t - x
        xi = xi - F / (1.0 + dU0(xi) * t)
    u = U0(xi)
    return np.log(rho0) + u / cs, u


def simple_wave_analysis(nml, lines, out):
    idx, L, k, cs, A, rho0 = wave1d_setup(nml)
    tb = 1.0 / (cs * A * k)
    rows = []
    for (t, step, a) in lines:
        x, s, u = a[:, 0], np.log(a[:, 1]), a[:, 2 + idx]
        s_ex, u_ex = simple_wave_exact(x, t, A, k, cs, rho0, L)
        rows.append(dict(t=t, t_over_tb=t / tb, Linf_s=float(np.max(np.abs(s - s_ex))),
                         Linf_u=float(np.max(np.abs(u - u_ex))),
                         zeta_max=float(np.max(a[:, 8]))))
    out['simple'] = dict(t_break=tb, rows=rows, final=rows[-1])


# ---------------------------------------------------------------- per case
def analyze(case, plot=False):
    out = dict(case=case)
    nml = read_nml(os.path.join(case, 'calliope.in'))
    st = os.path.join(case, 'exit_status')
    out['exit'] = open(st).read().strip() if os.path.exists(st) else None
    std = os.path.join(case, 'calliope.out.std')
    if os.path.exists(std):
        txt = open(std, errors='replace').read()
        if 'This run is stopped' in txt:  # time.F90 reset_dt, reset_method = 'stop' (exit code 0)
            out['exit'] = '%s, stopped by the CFL check' % out['exit']
        if 'Wall time has exceeded' in txt:  # terminate.F90, max_wall_time (exit code 0)
            out['exit'] = '%s, terminated by max_wall_time' % out['exit']

    cons = load(os.path.join(case, 'conservation.dat'))
    if cons is not None:
        M0 = cons[0, 2]
        Px, Py, Pz = cons[:, 4], cons[:, 5], cons[:, 6]
        P0 = np.array([Px[0], Py[0], Pz[0]])
        # momentum scale: sum |m| dV ~ M * max|u|; use sqrt(2 M E_kin,max)
        pscale = math.sqrt(2.0 * M0 * max(cons[:, 7].max(), 1e-300))
        out['cons'] = dict(n=len(cons), t_end=float(cons[-1, 0]),
                           max_abs_dM_over_M0=float(np.max(np.abs(cons[:, 3]))),
                           final_dM_over_M0=float(cons[-1, 3]),
                           max_abs_dP_over_scale=float(np.max(np.abs(np.stack([Px, Py, Pz], 1) - P0)) / pscale),
                           rho_min=float(cons[:, 11].min()), rho_max=float(cons[:, 12].max()),
                           zeta_max=float(cons[:, 13].max()),
                           E0=float(cons[0, 10]), E_end=float(cons[-1, 10]),
                           dmass_corr_tot=float(cons[-1, 17]))

    bud = load(os.path.join(case, 'budget.dat'))
    if bud is not None:
        dEdt = bud[:, 6]
        visc = -bud[:, 11]
        scale_E = max(np.max(np.abs(dEdt)), np.max(np.abs(visc)), 1e-300)
        dMs = max(np.max(bud[:, 5]), 1e-300)
        dt = bud[:, 1]
        out['budget'] = dict(
            n=len(bud),
            max_R_zeta_ratio=float(np.max(bud[:, 13])),
            max_abs_R_zeta_over_visc=float(np.max(np.abs(bud[:, 12]) / np.maximum(visc, 1e-300))) if np.any(visc > 0) else None,
            max_abs_dMdt_raw_over_abs=float(np.max(np.abs(bud[:, 3]) / np.maximum(bud[:, 5], 1e-300))),
            max_abs_dMdt_proj_over_abs=float(np.max(np.abs(bud[:, 4]) / np.maximum(bud[:, 5], 1e-300))),
            max_abs_dMdt_actual_over_abs=float(np.max(np.abs(bud[:, 2]) / np.maximum(bud[:, 5], 1e-300))),
            max_abs_err_t_over_scale=float(np.max(np.abs(bud[:, 8])) / scale_E),
            int_err_t=float(np.sum(bud[:, 8] * dt)),
            # time-integration error of E relative to |E(0)| (E is not positive definite in general)
            int_err_t_over_E0=float(np.sum(bud[:, 8] * dt) / max(abs(out['cons']['E0']), 1e-300)) if 'cons' in out else None,
            int_eps_ideal=float(np.sum(bud[:, 9] * dt)),
            int_eps_shock=float(np.sum(bud[:, 10] * dt)),
            int_visc=float(np.sum(visc * dt)),
            int_dEdt_actual=float(np.sum(dEdt * dt)),
            int_eps_nu=float(np.sum(bud[:, 16] * dt)),
            int_eps_svv=float(np.sum(bud[:, 20] * dt)) if bud.shape[1] > 20 else 0.0,
            dE_total=float(bud[-1, 19] - bud[0, 18]),
        )

    cfl = load(os.path.join(case, 'cfl.dat'))
    if cfl is not None:
        out['cfl'] = dict(dt_min=float(cfl[:, 1].min()), dt_max=float(cfl[:, 1].max()),
                          min_dt_adv=float(cfl[:, 2].min()), min_dt_visc=float(cfl[:, 3].min()),
                          max_dt_over_limit=float(np.max(cfl[:, 1] / np.minimum(cfl[:, 2], cfl[:, 3]))),
                          max_eps_svv=float(np.max(cfl[:, 9:12])) if cfl.shape[1] > 9 else 0.0)

    lines = read_lines(case)
    if lines and strv(nml.get('init_type', "'zero'")) == 'riemann':
        x, rho, u, rho_ex, u_ex = riemann_analysis(nml, lines, out)
        if plot and rho_ex is not None:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(3, 1, figsize=(7, 8), sharex=True)
            ax[0].plot(x, rho_ex, 'k-', lw=0.8, label='exact')
            ax[0].plot(x, rho, 'r.', ms=2, label='CALLIOPE')
            ax[0].set_ylabel(r'$\rho$'); ax[0].legend()
            ax[1].plot(x, u_ex, 'k-', lw=0.8); ax[1].plot(x, u, 'r.', ms=2)
            ax[1].set_ylabel(r'$u$')
            a = lines[-1][2]
            ax[2].plot(x, a[:, 8], 'b-', lw=0.8)
            ax[2].set_ylabel(r'$\zeta$'); ax[2].set_xlabel('x')
            ax[0].set_title('%s  t = %.4f' % (os.path.basename(case.rstrip('/')), lines[-1][0]))
            fig.tight_layout()
            fig.savefig(os.path.join(case, 'profile.png'), dpi=120)
            plt.close(fig)
    elif lines and strv(nml.get('init_type', "'zero'")) == 'sound_wave':
        sound_analysis(nml, lines, out)
    elif lines and strv(nml.get('init_type', "'zero'")) == 'simple_wave':
        simple_wave_analysis(nml, lines, out)
    elif lines:
        out['lines'] = dict(n=len(lines), t_last=lines[-1][0])

    with open(os.path.join(case, 'summary.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)
    return out


def brief(o):
    s = ['== %s  (exit: %s)' % (o['case'], o.get('exit'))]
    c = o.get('cons')
    if c:
        s.append('  cons: n=%d t_end=%.5g max|dM/M0|=%.2e max|dP|/scale=%.2e rho=[%.4g, %.4g] zeta_max=%.3g corr=%.2e'
                 % (c['n'], c['t_end'], c['max_abs_dM_over_M0'], c['max_abs_dP_over_scale'],
                    c['rho_min'], c['rho_max'], c['zeta_max'], c['dmass_corr_tot']))
    b = o.get('budget')
    if b:
        s.append('  budget: max R_zeta/visc=%s  |dMdt raw/proj/actual|/abs = %.1e/%.1e/%.1e  max|err_t|/scale=%.2e'
                 % ('%.2e' % b['max_R_zeta_ratio'], b['max_abs_dMdt_raw_over_abs'], b['max_abs_dMdt_proj_over_abs'],
                    b['max_abs_dMdt_actual_over_abs'], b['max_abs_err_t_over_scale']))
        s.append('          int: dE=%.4e = ideal %.3e + shock %.4e (visc %.4e) + svv %.4e + nu %.2e + err_t %.2e'
                 % (b['int_dEdt_actual'], b['int_eps_ideal'], b['int_eps_shock'], -b['int_visc'],
                    b['int_eps_svv'], b['int_eps_nu'], b['int_err_t']))
        if b.get('int_err_t_over_E0') is not None:
            s.append('          int err_t/|E0| = %.2e' % b['int_err_t_over_E0'])
    c = o.get('cfl')
    if c:
        s.append('  cfl: dt=[%.3e, %.3e] min dt_adv=%.3e min dt_visc=%.3e max dt/limit=%.3f max eps_svv=%.3e'
                 % (c['dt_min'], c['dt_max'], c['min_dt_adv'], c['min_dt_visc'], c['max_dt_over_limit'],
                    c['max_eps_svv']))
    w = o.get('sound')
    if w:
        s.append('  sound: max|amp err|=%.2e max|phase err|=%.2e rad; at t=%.3f harmonics/A=%.2e'
                 % (w['max_abs_amp_err'], w['max_abs_phase_err'], w['rows'][-1]['t'], w['rows'][-1]['harmonics']))
    w = o.get('simple')
    if w:
        r = w['final']
        s.append('  simple wave: t/t_b=%.3f Linf(s)=%.3e Linf(u)=%.3e zeta_max=%.3e'
                 % (r['t_over_tb'], r['Linf_s'], r['Linf_u'], r['zeta_max']))
    f = o.get('final')
    if f:
        s.append('  final t=%.4f valid=%s L1(rho)=%.3e rel=%.3e L1(u)=%.3e TV/TVex=%.4f rho=[%.4f, %.4f]'
                 % (f['t'], f['valid'], f['L1_rho'], f['rel_L1_rho'], f['L1_u'], f['TV_rho'] / f['TV_rho_exact'],
                    f['rho_min'], f['rho_max']))
        for st in o['stars']:
            if 'rho_s' in st:
                s.append('  star %d: rho*=%.6f (exact %.6f, err %+.2e, spread %.2e) u* err=%+.2e  [%d cells]'
                         % (st['interface'], st['rho_s'], st['rho_s_exact'], st['rho_s_rel_err'], st['rho_s_spread'],
                            st['u_s_err'], st['ncells']))
        for m in o['shocks']:
            s.append('  shock %d%s (M1=%.3f): dx=%+.2f cells thick=%.2f cells over=%.3f under=%.3f speed err=%s [win %.0f]'
                     % (m['interface'], m['side'], m['mach'], m['dx_cells'], m['thickness_cells'], m['overshoot'],
                        m['undershoot'], '%+.2e (n=%d)' % (m['speed_rel_err'], m['n_fit']) if 'speed_rel_err' in m else '-',
                        m['window_cells']))
    return '\n'.join(s)


if __name__ == '__main__':
    args = sys.argv[1:]
    plot = '--plot' in args
    for case in [a for a in args if a != '--plot']:
        print(brief(analyze(case, plot)))
