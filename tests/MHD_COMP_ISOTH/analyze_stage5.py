#!/usr/bin/env python3
"""Stage 5 analysis: the hole of rho of turb4 and turb4b on 64^3 (runs/stage4/turb) against 128^3
(runs/stage5/res, gen_stage5.py res), and the upper limit of M_s on 64^3 (runs/stage5/ms, gen_stage5.py ms).

    python3 tests/MHD_COMP_ISOTH/analyze_stage5.py res [runs] [--json f]   time series and hole shapes (on ps1)
    python3 tests/MHD_COMP_ISOTH/analyze_stage5.py plot f.json [out.png]    figure from the JSON
    python3 tests/MHD_COMP_ISOTH/analyze_stage5.py ms [runs] [--json f]    the same for M_s = 1, 2, 3, 4 on 64^3
                                                                            and the verdict of each case
    python3 tests/MHD_COMP_ISOTH/analyze_stage5.py msplot f.json [out.png]  figure from the JSON of ms
    python3 tests/MHD_COMP_ISOTH/analyze_stage5.py fcheck <case> <eig.json> [<case> <eig.json> ...]
                                                                            F_svv of cfl.dat against svv3d.py

res reads, for each case, conservation.dat and cfl.dat of all runs (written at every step) and the fields of
every out3d record and of restart/, and prints
- the times at which rho_min falls below RHO_CROSS, and F_svv rho_min, dt lambda_svv, cmax in the deep phase;
- the shape of the hole at each record with rho_min < SHAPE_RHO: position, neighbours along x, y, z (B0 along
  z), w10 (cells through the min with rho < 10 rho_min along each axis), the half width (distance at which
  ln rho is back to ln rho_min / 2, in cells and in cells of 64^3) and the counts of cells below 1e-2, 1e-3;
- the change of M, max div B and non-finite values.
ms adds, for each M_s, dt0 (the first dt), min dt/dt0, the min of rho_min and its t/tau0, max F_svv and the
verdict fixed before the runs: S = t_end (TURB_NTAU tau0) reached with dt >= dt0/MS_DT_FALL throughout, H =
rho_min < DEEP_RHO and dt < dt0/MS_DT_FALL; 'not S' otherwise once dt < dt0/MS_DT_FALL (dt does not grow
again with increase_dt off). It also prints when dt first fell below dt0/MS_DT_FALL and min(dt_adv, dt_visc)/dt0,
which can grow again.
"""
import os
import sys
import json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_ic  # noqa: E402
from analyze_stage1 import load_runs, read_nml  # noqa: E402
from analyze_stage4 import out3d_records, out3d_at  # noqa: E402
from gen_stage4 import TURB_NTAU  # noqa: E402

PAIRS = (('turb4', 'stage4/turb/turb4_N64_a', 'stage5/res/turb4_N128_a'),
         ('turb4b', 'stage4/turb/turb4b_N64', 'stage5/res/turb4b_N128'))
RHO_CROSS = (0.5, 0.3, 0.1, 1e-2, 1e-3, 1e-4)
SHAPE_RHO = 0.1
DEEP_RHO = 1e-3
NSERIES = 3000   # rows of each time series kept in the JSON
# column (0-based) of conservation.dat and cfl.dat
C_T, C_M, C_RMIN, C_DIVB = 0, 2, 11, 19
F_T, F_DT, F_DTADV, F_DTVISC, F_CX, F_LAM, F_F, F_DTLAM = 0, 1, 2, 3, 4, 12, 13, 14
AXES = (('x', 0), ('y', 2), ('z', 1))   # arrays are (x, z, y)
# M_s and case of the problems with rms|dB| = B0 on 64^3 (turb1, turb2b, turb3b, turb4b)
MS_CASES = ((1.0, 'stage4/turb/turb1_N64'), (2.0, 'stage5/ms/turb2b_N64'), (3.0, 'stage5/ms/turb3b_N64'),
            (4.0, 'stage4/turb/turb4b_N64'))
MS_DT_FALL = 10.0


def series(case):
    """Time series of the case (all runs): conservation.dat and cfl.dat rows at the same t."""
    cons = np.vstack(load_runs(case, 'conservation.dat'))
    cfl = np.vstack(load_runs(case, 'cfl.dat'))
    i = np.clip(np.searchsorted(cons[:, C_T], cfl[:, F_T]), 0, len(cons) - 1)
    on = np.abs(cons[i, C_T] - cfl[:, F_T]) <= 1e-12*max(1.0, cfl[-1, F_T])
    c, f = cons[i[on]], cfl[on]
    return dict(t=f[:, F_T], dt=f[:, F_DT], dt_adv=f[:, F_DTADV], dt_visc=f[:, F_DTVISC],
                cmax=np.max(f[:, F_CX:F_CX + 3], axis=1), lam=f[:, F_LAM], F=f[:, F_F], dtlam=f[:, F_DTLAM],
                rho_min=c[:, C_RMIN]), cons, cfl


def crossings(s):
    out = {}
    for th in RHO_CROSS:
        j = np.nonzero(s['rho_min'] < th)[0]
        out['%g' % th] = float(s['t'][j[0]]) if len(j) else None
    return out


def deep(s):
    sel = s['rho_min'] < DEEP_RHO
    if not np.any(sel):
        return None
    fr = s['F'][sel]*s['rho_min'][sel]
    return dict(n=int(np.sum(sel)), t0=float(s['t'][sel][0]), F_rho_min=[float(fr.min()), float(np.median(fr)), float(fr.max())],
                F_rho_min_last=float(fr[-1]), dtlam_max=float(s['dtlam'].max()), cmax_last=float(s['cmax'][-1]))


def health(cons, cfl):
    return dict(dM=float(np.max(np.abs(cons[:, C_M] - cons[0, C_M]))/cons[0, C_M]),
                divb_max=float(np.max(cons[:, C_DIVB])),
                finite=bool(np.all(np.isfinite(cons)) and np.all(np.isfinite(cfl))),
                t_last=float(cfl[-1, F_T] + cfl[-1, F_DT]), steps=int(len(cfl)), dt_last=float(cfl[-1, F_DT]))


def half_width(l, i0, target):
    """Distance (cells) from i0 on each side at which the periodic profile l is back to target (linear in
    between); None when it is not within half the period."""
    n = len(l)
    out = []
    for sgn in (1, -1):
        d = None
        for o in range(1, n//2 + 1):
            a, b = l[(i0 + sgn*(o - 1)) % n], l[(i0 + sgn*o) % n]
            if b >= target:
                d = o - 1 + (target - a)/(b - a)
                break
        out.append(d)
    return None if None in out else out


def shape(rho, N):
    """Shape of the hole at the min of rho (see the module docstring)."""
    i0 = np.unravel_index(np.argmin(rho), rho.shape)
    rmin = float(rho[i0])
    lr = np.log(rho)
    r = dict(rho_min=rmin, at={a: float(i0[ax])/N for a, ax in AXES},
             below={'%g' % th: int(np.sum(rho < th)) for th in (1e-2, 1e-3)}, axes={})
    for a, ax in AXES:
        sl = list(i0)
        sl[ax] = slice(None)
        p = rho[tuple(sl)]
        lp = lr[tuple(sl)]
        i = i0[ax]
        w10 = 1
        for sgn in (1, -1):
            o = 1
            while o < N and p[(i + sgn*o) % N] < 10*rmin:
                w10 += 1
                o += 1
        hw = half_width(lp, i, 0.5*np.log(rmin))
        r['axes'][a] = dict(nb=[float(p[(i - 1) % N]), float(p[(i + 1) % N])], w10=w10,
                            half=None if hw is None else [float(hw[0]), float(hw[1])],
                            width=None if hw is None else float(sum(hw)),
                            width64=None if hw is None else float(sum(hw))*64/N)
    return r


def shapes(case):
    """shape() at each out3d record and at restart/ with rho_min < SHAPE_RHO (the repeated record at the start
    of a continued run is left out)."""
    d = read_nml(os.path.join(case, 'calliope.in'))
    N = int(d['nlx'])
    g = gen_ic.Grid(N, N, N)
    out, seen = [], set()
    for t, path, rec in out3d_records(case):
        if round(t, 12) in seen:
            continue
        seen.add(round(t, 12))
        rho = g.ifft(out3d_at(g, path, rec)['rho'])
        if rho.min() < SHAPE_RHO:
            out.append(dict(t=t, src='%s:%d' % (os.path.relpath(path, case), rec), **shape(rho, N)))
    lk, mk, bk, tt = gen_ic.read_restart(os.path.join(case, 'restart'), g)
    rho = np.exp(g.ifft(lk))
    if round(float(tt), 12) not in seen and rho.min() < SHAPE_RHO:
        out.append(dict(t=float(tt), src='restart', **shape(rho, N)))
    return N, out


def thin(s, n=NSERIES):
    k = max(1, len(s['t'])//n)
    idx = np.unique(np.r_[np.arange(0, len(s['t']), k), len(s['t']) - 1])
    return {q: v[idx].tolist() for q, v in s.items()}


def fmt(x, f='%.3e'):
    return '-' if x is None else f % x


def summary(base, c):
    """(result, full time series) of the case base/c as printed by res; None when it has no output."""
    case = os.path.join(base, c)
    if not os.path.exists(os.path.join(case, 'cfl.dat')) and not os.path.isdir(os.path.join(case, 'chunk_1')):
        print('%s: no output' % case)
        return None
    s, cons, cfl = series(case)
    N, sh = shapes(case)
    r = dict(case=c, N=N, cross=crossings(s), deep=deep(s), health=health(cons, cfl), shapes=sh,
             series=thin(s))
    h = r['health']
    print('%s  N = %d  t = %.6f  %d steps  dt %.3e  dM %.1e  max div B %.1e  finite %s'
          % (c, N, h['t_last'], h['steps'], h['dt_last'], h['dM'], h['divb_max'], h['finite']))
    print('  t at rho_min <', ', '.join('%s: %s' % (k, fmt(v, '%.4f')) for k, v in r['cross'].items()))
    if r['deep']:
        dp = r['deep']
        print('  rho_min < %g from t = %.5f (%d steps): F rho_min min/median/max %.3f %.3f %.3f, last %.3f, '
              'max dt lam %.3f, cmax %.1f' % (DEEP_RHO, dp['t0'], dp['n'], *dp['F_rho_min'], dp['F_rho_min_last'],
                                              dp['dtlam_max'], dp['cmax_last']))
    for x in sh:
        ax = x['axes']
        print('  %-22s t = %.5f rho_min %.2e at (%.3f, %.3f, %.3f)  n(<1e-2) %d n(<1e-3) %d'
              % (x['src'], x['t'], x['rho_min'], x['at']['x'], x['at']['y'], x['at']['z'],
                 x['below']['0.01'], x['below']['0.001']))
        print('    ' + '  '.join('%s: nb %.1e %.1e w10 %d width %s (%s cells of 64^3)'
                                 % (a, *ax[a]['nb'], ax[a]['w10'], fmt(ax[a]['width'], '%.2f'),
                                    fmt(ax[a]['width64'], '%.2f')) for a, _ in AXES))
    return r, s


def res(base, jout=None):
    out = {}
    for prob, c64, c128 in PAIRS:
        out[prob] = {}
        for c in (c64, c128):
            x = summary(base, c)
            if x is not None:
                out[prob][str(x[0]['N'])] = x[0]
    if jout:
        with open(jout, 'w') as f:
            json.dump(out, f)
    return out


def verdict(reached, fell, deep_hole):
    """S, H or neither (module docstring). dt does not grow again, so a case with dt < dt0/MS_DT_FALL is not S."""
    if deep_hole and fell:
        return 'H'
    if not fell:
        return 'S' if reached else 'running'
    return 'not S: dt < dt0/%g, %s' % (MS_DT_FALL, 't_end' if reached else 'before t_end')


def ms(base, jout=None):
    """summary() of MS_CASES with the verdict (module docstring)."""
    out = {}
    for Ms, c in MS_CASES:
        x = summary(base, c)
        if x is None:
            continue
        r, s = x
        h = r['health']
        tau0 = np.pi/Ms
        dt0 = float(s['dt'][0])
        i = int(np.argmin(s['rho_min']))
        below = np.nonzero(s['dt'] < dt0/MS_DT_FALL)[0]
        fell = below.size > 0
        reached = h['t_last'] >= TURB_NTAU*tau0 - 0.5*h['dt_last']
        dtcfl = np.minimum(s['dt_adv'], s['dt_visc'])/dt0   # dt_cfl can grow again, dt cannot (increase_dt off)
        r.update(Ms=Ms, tau0=tau0, dt0=dt0, dt_min_rel=float(s['dt'].min()/dt0), rho_min=float(s['rho_min'][i]),
                 t_rho_min=float(s['t'][i]/tau0), F_max=float(s['F'].max()), t_last_tau0=h['t_last']/tau0,
                 fall=dict(t=float(s['t'][below[0]]/tau0), rho_min=float(s['rho_min'][below[0]]),
                           F=float(s['F'][below[0]])) if fell else None,
                 dtcfl_rel=dict(min=float(dtcfl.min()), last=float(dtcfl[-1])),
                 verdict=verdict(reached, fell, bool(s['rho_min'].min() < DEEP_RHO)))
        print('  M_s = %g: t = %.4f tau0, dt0 %.3e, min dt/dt0 %.2e, min rho_min %.3e at t = %.3f tau0, max F_svv %.3g'
              ' -> %s' % (Ms, r['t_last_tau0'], dt0, r['dt_min_rel'], r['rho_min'], r['t_rho_min'], r['F_max'],
                          r['verdict']))
        if fell:
            print('    dt < dt0/%g from t = %.4f tau0 (rho_min %.3e, F_svv %.3g); dt_cfl/dt0 min %.2e, last %.2e'
                  % (MS_DT_FALL, r['fall']['t'], r['fall']['rho_min'], r['fall']['F'], dtcfl.min(), dtcfl[-1]))
        out['%g' % Ms] = r
    if jout:
        with open(jout, 'w') as f:
            json.dump(out, f)
    return out


# ordinal blue ramp of M_s = 1, 2, 3, 4 (dataviz reference palette, steps 250, 400, 550, 700; validated --ordinal)
MS_COLORS = ('#86b6ef', '#3987e5', '#1c5cab', '#0d366b')


def msplot(jin, png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = json.load(open(jin))
    fig, axs = plt.subplots(3, 1, figsize=(7, 8), sharex=True)
    for k, (Ms, _) in enumerate(MS_CASES):
        r = d.get('%g' % Ms)
        if r is None:
            continue
        s = r['series']
        t = np.array(s['t'])/r['tau0']
        kw = dict(color=MS_COLORS[k], lw=1.5, label='$M_s = %g$ (%s)' % (Ms, r['verdict'].split(':')[0]))
        for a, y in zip(axs, (np.array(s['rho_min']), np.array(s['dt'])/r['dt0'], np.array(s['F']))):
            a.semilogy(t, y, **kw)
        axs[0].annotate('%g' % Ms, (t[-1], s['rho_min'][-1]), xytext=(3, 0), textcoords='offset points',
                        color='#333333', fontsize=9, va='center')
    axs[1].axhline(1/MS_DT_FALL, color='#888888', lw=0.8, ls='--')
    axs[0].axhline(DEEP_RHO, color='#888888', lw=0.8, ls='--')
    for a, lab in zip(axs, (r'$\rho_{\min}$', r'$\Delta t/\Delta t_0$', r'$F_{\rm svv}$')):
        a.set_ylabel(lab)
        a.grid(alpha=0.3, lw=0.5)
    axs[2].set_xlabel(r'$t/\tau_0$')
    axs[0].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(png, dpi=120)
    print(png)


def fcheck(pairs):
    """F_svv in cfl.dat at the restart time of each case against svv3d.py state <case> restart --eig (global/rho).
    The cfl.dat row nearest to the restart time is taken; its distance from the restart in dt is printed."""
    for case, ej in pairs:
        cfl = np.vstack(load_runs(case, 'cfl.dat'))
        e = json.load(open(ej))
        i = int(np.argmin(np.abs(cfl[:, F_T] - e['t'])))
        F = [v['F'] for v in e['variants'] if v['eps'] == 'global' and v['weight'] == 'rho'][0]
        print('%-30s restart t %.9f  cfl.dat t %.9f (%+.2f dt)  F cfl.dat %.6g  F svv3d %.6g  rel %+.1e'
              % (case, e['t'], cfl[i, F_T], (cfl[i, F_T] - e['t'])/cfl[i, F_DT], cfl[i, F_F], F, F/cfl[i, F_F] - 1))


def plot(jin, png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = json.load(open(jin))
    col = {'64': '#2a6fdb', '128': '#d1495b'}
    fig, axs = plt.subplots(3, 2, figsize=(9, 8), sharex='col')
    for j, prob in enumerate(('turb4', 'turb4b')):
        for N, r in sorted(d.get(prob, {}).items(), key=lambda x: int(x[0])):
            s = r['series']
            t = np.array(s['t'])
            kw = dict(color=col[N], lw=1.5, label='%s$^3$' % N)
            axs[0, j].semilogy(t, s['rho_min'], **kw)
            axs[1, j].semilogy(t, s['dt'], **kw)
            axs[2, j].semilogy(t, s['F'], **kw)
        axs[0, j].set_title(prob)
        axs[2, j].set_xlabel('t')
    for i, lab in enumerate((r'$\rho_{\min}$', r'$\Delta t$', r'$F_{\rm svv}$')):
        axs[i, 0].set_ylabel(lab)
    for a in axs.flat:
        a.grid(alpha=0.3, lw=0.5)
    axs[0, 0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(png, dpi=120)
    print(png)


if __name__ == '__main__':
    args = sys.argv[1:]
    jout = None
    if '--json' in args:
        k = args.index('--json')
        jout = args[k + 1]
        args = args[:k] + args[k + 2:]
    if args[0] == 'res':
        res(args[1] if len(args) > 1 else 'runs', jout)
    elif args[0] == 'plot':
        plot(args[1], args[2] if len(args) > 2 else os.path.splitext(args[1])[0] + '.png')
    elif args[0] == 'ms':
        ms(args[1] if len(args) > 1 else 'runs', jout)
    elif args[0] == 'msplot':
        msplot(args[1], args[2] if len(args) > 2 else os.path.splitext(args[1])[0] + '.png')
    elif args[0] == 'fcheck':
        fcheck(list(zip(args[1::2], args[2::2])))
