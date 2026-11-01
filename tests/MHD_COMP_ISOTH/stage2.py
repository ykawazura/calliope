#!/usr/bin/env python3
"""Stage 2: the formulations A, B, C, the form of the shock viscosity, the place of the dealiasing filter and the
form of the SVV, compared with the 1D prototype proto1d.py (REPORT_stage2.md).

    python3 stage2.py check                 identities of proto1d.budget and the ideal parts of A, B, C
    python3 stage2.py run [substring ...]   run the cases (6 processes) into runs/stage2/<case>.json and .npz
                                            (cases with a .json are skipped; -f reruns them)
    python3 stage2.py table                 Markdown tables runs/stage2/tables.md and figures fig/stage2_*.png
    python3 stage2.py list                  the case names

Every variant of a (problem, N) uses the same dt, SAFETY*DTMIN of gen_stage1.py (the ES-SVV default of C);
'lin_own' uses the dt of the linear SVV in gen_stage1.py. A run that fails (non-finite values, CFL stop,
dt lambda_svv > 2.51) is repeated with dt/2 and dt/4; all attempts are kept.
The semi-discrete rates are proto1d.budget() of the stage-1 state of each step; the time integrals are sums of
dt times that value (not the stage-weighted sums of the Fortran budget.dat).
"""
import sys
import os
import json
import math
import time
import glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from proto1d import Isoth1D
import analyze_stage1 as an
from gen_stage1 import RIEMANN, DTMIN, DT512_LINEAR, SAFETY, fixed_dt

OUT = os.path.normpath(os.path.join(HERE, '..', '..', 'runs', 'stage2'))
VARIANTS = dict(
    C={}, B=dict(formulation='B'), A=dict(formulation='A'),
    Cexp=dict(visc_form='expanded'), Aexp=dict(formulation='A', visc_form='expanded'),
    D0=dict(filter_m=False), D2=dict(filter_s=True), D3=dict(filter_m=False, deriv_mask='mask'),
    D4=dict(filter_rhs=True), D5=dict(pad=2),
    lin=dict(svv_form='linear'), lin_own=dict(svv_form='linear'))
SHOCK_VARIANTS = ['C', 'B', 'A', 'Cexp', 'Aexp', 'D0', 'D2', 'D3', 'D5']
NS = (512, 1024, 2048)


def cases():
    cs = []
    for name in RIEMANN:
        for N in NS:
            for v in SHOCK_VARIANTS:
                cs.append(dict(kind='riemann', prob=name, N=N, var=v))
    cs.append(dict(kind='riemann', prob='U2', N=1024, var='D4'))
    for name in ('U2', 'U4'):
        for N in (1024, 2048):
            for v in ('lin', 'lin_own'):
                cs.append(dict(kind='riemann', prob=name, N=N, var=v))
    for N in (32, 64):
        for tag in ('C0', 'C2'):
            for f in 'ABC':
                cs.append(dict(kind='sound_wave', N=N, tag=tag, var=f))
    for N in (32, 64, 128, 256, 512):
        for tag in ('C0_nosvv', 'C0_svv', 'C2_svv'):
            for f in 'ABC':
                cs.append(dict(kind='simple_wave', N=N, tag=tag, var=f))
    for c in cs:
        c['name'] = ('%s_N%d_%s' % (c['prob'], c['N'], c['var']) if c['kind'] == 'riemann'
                     else '%s_N%d_%s_%s' % (c['kind'].split('_')[0], c['N'], c['tag'], c['var']))
    return cs


def setup(c):
    """(proto, dt, nstep, nwrite, nml) of a case; nml in the form read by analyze_stage1."""
    N = c['N']
    nml = dict(init_type="'%s'" % c['kind'], dir="'x'", lx='1.0', nlx=str(N), beta0='2.0')
    if c['kind'] == 'riemann':
        rL, rR, uL, uR, te = RIEMANN[c['prob']]
        dtmax = DT512_LINEAR[c['prob']] * 512 / N if c['var'] == 'lin_own' else DTMIN[(c['prob'], N)]
        tm = fixed_dt(te, SAFETY * dtmax)
        p = Isoth1D(N, **VARIANTS[c['var']])
        p.riemann(rL, rR, uL, uR)
        nml.update(rho_l=repr(rL), rho_r=repr(rR), u_l='%r 0.0 0.0' % uL, u_r='%r 0.0 0.0' % uR)
        return p, tm['dt'], tm['nstep'], 20, nml
    # the smooth cases of gen_stage1.batch_smooth
    kw = dict(formulation=c['var'])
    if c['tag'].startswith('C0'):
        kw['C'] = 0.0
    if c['tag'].endswith('nosvv'):
        kw['svv_eps0'] = 0.0
    p = Isoth1D(N, **kw)
    if c['kind'] == 'sound_wave':
        amp, dt, nstep, nwrite = 1e-3, 2.0 * math.pi / 4000, 4000, 8
    else:
        amp, dt, nstep, nwrite = 0.5, 2.5e-4, 4000, 4
    p.wave1d(c['kind'], amp=amp, mode=1, rho0=1.0)
    nml.update(amp=repr(amp), mode='1', rho0='1.0')
    return p, dt, nstep, nwrite, nml


def profile(p):
    """Columns of out1d/line_*.dat read by analyze_stage1: 0 x, 1 rho, 2 u_x, 8 zeta (N points for pad > 1)."""
    st = p.state()
    k = p.pad
    a = np.zeros((p.N, 9))
    a[:, 0], a[:, 1], a[:, 2], a[:, 8] = st['x'][::k], st['rho'][::k], st['u'][::k], st['zeta'][::k]
    return a


class Stats:
    """Running maxima, minima and time integrals over the hist entries (stage-1 state of each step)."""
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
                self.d.update(M0=h['M'], P0=h['P'], E0=h['E'], n_E_up=0)
            h0, dt = self.h0, h['dt']
            self.mx('max_dM_over_M0', abs(h['M'] - h0['M']) / h0['M'])
            self.mx('max_dP', abs(h['P'] - h0['P']))
            self.mx('max_E_kin', h['E_kin'])
            self.mn('rho_min', h['rho_min'])
            self.mx('rho_max', h['rho_max'])
            self.mx('zeta_max', h['zeta_max'])
            self.mx('max_dt_over_limit', dt / min(h['dt_adv'], h['dt_visc']))
            self.mx('max_dt_lam_svv', dt * h['lam_svv'])
            if self.E_prev is not None and h['E'] > self.E_prev:
                self.d['n_E_up'] += 1
                self.mx('max_E_up_over_E0', (h['E'] - self.E_prev) / abs(h0['E']))
            self.E_prev = h['E']
            self.d.update(E_end=h['E'], M_end=h['M'], P_end=h['P'])
            for tag in ('', '_raw'):
                self.mx('max_dM' + tag, abs(h['dM' + tag]) / h['A_M'] if h['A_M'] > 0 else 0.0)
                self.mx('max_dP' + tag, abs(h['dP' + tag]) / h['A_P'] if h['A_P'] > 0 else 0.0)
                if h['visc'] > 0.0:
                    self.mx('max_R_zeta' + tag, abs(h['R_zeta' + tag]) / h['visc'])
                if h['D_svv'] != 0.0 and h['A_svv'] > 0.0:
                    self.mx('max_R_svv' + tag, abs(h['R_svv' + tag]) / h['A_svv'])
                for k in ('ideal_eps', 'visc_eps', 'svv_eps', 'dM', 'dP'):
                    self.add('int_' + k + tag, dt * h[k + tag])
            if h['D_svv'] != 0.0:
                self.mn('min_D_svv', h['D_svv'])
                if h['svv_eps'] > 0.0:
                    self.add('n_eps_svv_pos', 1)
            self.add('int_visc', dt * h['visc'])
            self.add('int_D_svv', dt * h['D_svv'])

    def result(self):
        d = dict(self.d)
        if 'M0' in d:
            d['max_dP_over_scale'] = d['max_dP'] / math.sqrt(2.0 * d['M0'] * max(d['max_E_kin'], 1e-300))
        return {k: (None if isinstance(v, float) and not math.isfinite(v) else v) for k, v in d.items()}


def rh_residuals(out, nml):
    """Rankine-Hugoniot speeds of each shock from the measured downstream (star) state and the initial upstream
    state: S_mass = [m]/[rho], S_mom = [rho u^2 + cs^2 rho]/[m]; relative difference from the fitted speed."""
    cs2 = an.num(nml['beta0']) / 2.0
    rL, rR = an.num(nml['rho_l']), an.num(nml['rho_r'])
    uL, uR = an.nums(nml['u_l'])[0], an.nums(nml['u_r'])[0]
    states = [((rL, uL), (rR, uR)), ((rR, uR), (rL, uL))]
    for sh in out.get('shocks', []):
        st = out['stars'][sh['interface'] - 1]
        if 'speed_fit' not in sh or 'rho_s' not in st:
            continue
        ru, uu = states[sh['interface'] - 1][0 if sh['side'] == 'L' else 1]
        rd, ud = st['rho_s'], st['u_s']
        S_mass = (rd * ud - ru * uu) / (rd - ru)
        S_mom = ((rd * ud * ud + cs2 * rd) - (ru * uu * uu + cs2 * ru)) / (rd * ud - ru * uu)
        sh.update(S_mass=S_mass, S_mom=S_mom, rh_mass=sh['speed_fit'] / S_mass - 1.0,
                  rh_mom=sh['speed_fit'] / S_mom - 1.0)


def run_case(c):
    t_wall = time.time()
    attempts = []
    for fac in (1, 2, 4):
        p, dt, nstep, nwrite, nml = setup(c)
        dt, nstep = dt / fac, nstep * fac
        p.track_budget = True
        stats = Stats()
        lines = [(0.0, 0, profile(p))]
        status, t0 = 'ok', time.time()
        with np.errstate(all='ignore'):
            for j in range(1, nwrite + 1):
                status = p.run(dt, int(round(j * nstep / nwrite)) - p.step)
                stats.update(p.hist)
                p.hist = []
                if status != 'ok':
                    break
                lines.append((p.t, p.step, profile(p)))
        attempts.append(dict(dt_factor=fac, dt=dt, nstep=nstep, status=status, t_stop=p.t, steps=p.step,
                             wall=time.time() - t0, stats=stats.result()))
        if status == 'ok':
            break
    out = dict(case=c, nml=nml, attempts=attempts, ok=status == 'ok', t_last_line=lines[-1][0])
    if len(lines) > 1:
        if c['kind'] == 'riemann':
            an.riemann_analysis(nml, lines, out)
            rh_residuals(out, nml)
        elif c['kind'] == 'sound_wave':
            an.sound_analysis(nml, lines, out)
        else:
            an.simple_wave_analysis(nml, lines, out)
    a = lines[-1][2]
    np.savez_compressed(os.path.join(OUT, c['name'] + '.npz'), t=lines[-1][0], x=a[:, 0], rho=a[:, 1], u=a[:, 2],
                        zeta=a[:, 8])
    with open(os.path.join(OUT, c['name'] + '.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)
    return c['name'], [(at['dt_factor'], at['status'], round(at['t_stop'], 4)) for at in attempts], time.time() - t_wall


def cost(c):
    p, dt, nstep, nwrite, nml = setup(c)
    return nstep * p.Ng


def run(args):
    import multiprocessing as mp
    force = '-f' in args
    pats = [a for a in args if a != '-f']
    os.makedirs(OUT, exist_ok=True)
    todo = [c for c in cases() if (not pats or any(s in c['name'] for s in pats))
            and (force or not os.path.exists(os.path.join(OUT, c['name'] + '.json')))]
    todo.sort(key=cost, reverse=True)
    print('%d cases' % len(todo), flush=True)
    with mp.Pool(6) as pool:
        for name, att, wall in pool.imap_unordered(run_case, todo):
            print('%-28s %s  %.0fs' % (name, att, wall), flush=True)


def check():
    """Identities of budget() at a U4 N = 1024 state (t = 0.127, ES-SVV, C) mapped to every variant, the
    generic code path of C against the Stage 1 path, and the ideal parts of A, B against C on a smooth state."""
    p = Isoth1D(1024)
    p.riemann(1, 1, 4, -4)
    p.run(SAFETY * DTMIN[('U4', 1024)], 8000, stop_on_cfl=False)
    p.prepare()
    rho = np.exp(p.bwd(p.sh))
    u = p.bwd(p.mh) / rho

    def clone(**kw):
        q = Isoth1D(1024, **kw)
        q.eps, q.lam_svv, q.prepared = p.eps, p.lam_svv, 0
        q.sh = p.sh.copy()
        q.mh = p.mh.copy() if q.formulation == 'C' else q.fwd(u) * q.Fm
        return q
    g = clone(deriv_mask='mask')
    a, b = p.rhs(p.sh, p.mh), g.rhs(g.sh, g.mh)
    print('U4 N=1024 t=%.4f: generic C path vs Stage 1 path: ds/dt %.1e, dm/dt %.1e (max relative, after the m filter)' % (
        p.t, np.max(np.abs(a[0] - b[0])) / np.max(np.abs(a[0])),
        np.max(np.abs(p.Fm * (a[1] - b[1]))) / np.max(np.abs(a[1]))))
    print('%-6s %-19s %-19s %-19s %-19s %s' % ('', '|dM|/A_M', '|dP|/A_P', '|R_zeta|/visc', '|R_svv|/A_svv',
                                              'eps_ideal (raw)'))
    for v in ('C', 'B', 'A', 'Cexp', 'Aexp', 'D0', 'D2', 'D3', 'D4', 'lin'):
        q = clone(**VARIANTS[v])
        B = q.budget()
        f = lambda k, s: '%.1e (%.1e)' % (abs(B[k]) / s, abs(B[k + '_raw']) / s) if s > 0 else '-'
        print('%-6s %-19s %-19s %-19s %-19s %+.3e (%+.3e)' % (
            v, f('dM', B['A_M']), f('dP', B['A_P']), f('R_zeta', B['visc']),
            f('R_svv', B['A_svv']) if B['D_svv'] else '-', B['ideal_eps'], B['ideal_eps_raw']))
    print('(after the stage filter; before it in parentheses)')
    for N in (32, 64, 128):
        res = {}
        for f in 'ABC':
            q = Isoth1D(N, formulation=f, C=0.0, svv_eps0=0.0)
            x = q.x
            q.set_su(0.3 * np.sin(x) + 0.1 * np.cos(2 * x), 0.5 * np.cos(x))
            P, fl = q.parts(q.sh, q.mh)
            a, b = P['ideal']
            rdot = fl['rho'] * q.bwd(a * q.Fs)
            b = q.bwd(b * q.Fm)
            res[f] = (rdot, b if f == 'C' else fl['rho'] * b + fl['u'] * rdot)
        e = lambda f, i: np.max(np.abs(res[f][i] - res['C'][i])) / np.max(np.abs(res['C'][i]))
        print('smooth state N=%3d, ideal (d rho/dt, dm/dt) relative to C: A %.1e %.1e, B %.1e %.1e' % (
            N, e('A', 0), e('A', 1), e('B', 0), e('B', 1)))
    # D3 against C (U4 N = 1024, to t = 0.05): every term of dm/dt carries M_d, so the k >= k_c modes of m have
    # no tendency in D3; D3 = C once the initial m is band-limited
    dt = SAFETY * DTMIN[('U4', 1024)]
    ref = None
    for tag, kw, pre in (('C', {}, False), ('D3', VARIANTS['D3'], False), ('D3, initial m filtered', VARIANTS['D3'], True)):
        q = Isoth1D(1024, **kw)
        q.riemann(1, 1, 4, -4)
        if pre:
            q.mh = q.mh * q.mask
        hi0 = np.max(np.abs(q.mh * (1.0 - q.mask)))
        with np.errstate(all='ignore'):
            st = q.run(dt, int(round(0.05 / dt)))
        line = 'U4 N=1024 %-22s %-4s t=%.4f  max|m^(k >= k_c)| %.2e -> %.2e' % (
            tag, st, q.t, hi0, np.max(np.abs(q.mh * (1.0 - q.mask))))
        if ref is None:
            ref = q
        elif st == 'ok':
            line += '  vs C: s^ %.1e, m^ %.1e (max relative)' % (
                np.max(np.abs(q.sh - ref.sh)) / np.max(np.abs(ref.sh)), np.max(np.abs(q.mh - ref.mh)) / np.max(np.abs(ref.mh)))
        print(line)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'run':
        run(sys.argv[2:])
    elif cmd == 'check':
        check()
    elif cmd == 'list':
        for c in cases():
            print(c['name'])
    elif cmd == 'table':
        import stage2_table
        stage2_table.main(OUT)
    else:
        sys.exit(__doc__)
