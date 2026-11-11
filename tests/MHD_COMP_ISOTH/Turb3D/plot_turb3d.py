#!/usr/bin/env python3
"""Figures of the 3D decaying isothermal MHD turbulence (README.md) in the style of the paper (../stages/figstyle.py).

    python3 tests/MHD_COMP_ISOTH/Turb3D/plot_turb3d.py fetch                                inputs from ps1
    python3 tests/MHD_COMP_ISOTH/Turb3D/plot_turb3d.py [plot] [name ...] [--runs runs]     -> Turb3D/fig/<name>.pdf

turb1 (M_s = 1) unless noted; t in tau0 = pi/M_s. The local runs/ already hold the summaries of Stages 6-8, 13, 15
and 16 (turb_eps, turb_spec, turb_rho_pdf, turb_K_spec, turb_K, turb_kappa, turb_coef, turb_nomass). fetch adds
the plane z0 of the 256^3 run at t = tau0 (turb_slice), the cfl.dat tables of the Stage 4 and 'track' runs
(turb_track) and the summaries of Stage 5 (turb_ms) under runs/. Every figure prints the numbers to check against
the REPORTs: Stage 4 (8.3), 5 (6, 12), 6 (7, 8), 7, 8 (1), 13 (1, 7), 15 (1, 6) and 16.
"""
import os
import sys
import json
import numpy as np
TEST = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(TEST, '..', 'stages'))
import figstyle as fs
import figdata as fd

plt = fs.setup()
TURB = {64: 'turb1_N64', 128: 'turb1_N128_a', 256: 'turb1_N256'}       # stage6q/turb/summary.json
COLOR = {N: fs.SERIES[i] for i, N in enumerate(sorted(TURB))}         # the colour follows N
TAIL_T = (0.5, 1.0)                                                     # t of the density PDFs
KRUNS = ('ES', 'H_r8', 'L_r1', 'L_r2')                                 # 128^3 runs of Stage 15 (B* = H_r8)
KCOLOR = dict(ES=fs.COLORS[4], H_r8=fs.COLORS[5], L_r1=fs.COLORS[6], L_r2=fs.COLORS[2])
COEF = (('e0.5_C2', 'k', '-'), ('e1.0_C2', fs.COLORS[5], '-'), ('e0.25_C2', fs.COLORS[0], '-'),
        ('e0.5_C4', fs.COLORS[4], '--'), ('e0.5_C1', fs.COLORS[1], '--'))  # Stage 8: the control first
COEF_WIN, COEF_DELTA = (0.6, 1.2), 0.1                                  # window of L and threshold of R (Stage 8)
SLICE = ('stage6q/turb/turb1_N256', 'z0', np.pi, 256)                   # turb_slice: case, plane, t, N
LNRHO_LIM, ZETA_LIM = 2.0, (1e-4, 4e-2)                                 # colour ranges of turb_slice
TRACK = {64: (('stage4/turb/turb1_N64', r'fixed $\Delta t$'), ('stage6q/track/turb1_N64_t', 'track')),
         128: (('stage4/turb/turb1_N128_a', 'decrement'), ('stage6q/track/turb1_N128_t', 'track'))}
TRACK_SAFETY, TRACK_RAISE = 0.8, 0.7


def load(runs, *path):
    with open(os.path.join(runs, *path)) as f:
        return json.load(f)


def scheme(lab):
    """'H_r8' -> 'H(8)'."""
    s, _, r = lab.partition('_r')
    return '%s(%s)' % (s, r) if r else s


def turb_eps(runs):
    """eps(t) = -dE/dt (solid) and its shock-viscous part (dashed), N = 64, 128, 256."""
    from matplotlib.lines import Line2D
    S = load(runs, 'stage6q', 'turb', 'summary.json')['runs']
    fig, ax = fs.panels(1, 1)
    a = ax[0, 0]
    for N, name in sorted(TURB.items()):
        r = S[name]
        t = np.array(r['series']['t']) / r['info']['tau0']
        a.plot(t, r['series']['eps'], color=COLOR[N], lw=fs.LW, label=r'$N = %d$' % N)
        a.plot(t, r['series']['visc'], color=COLOR[N], lw=fs.LW - 1, ls='--')
        p = r['eps_peak']
        print('N=%-3d eps peak %.3f at t = %.4f tau0, t_last %.4f tau0, share of the shock viscosity %.4f'
              % (N, p['eps'], p['t'] / r['info']['tau0'], t[-1], r['share']['visc']))
    a.set_xlim(0.0, 3.0)
    a.set_ylim(0.0, 42.0)
    a.set_yticks((10, 20, 30, 40))
    a.set_xlabel(r'$t/\tau_0$')
    a.set_ylabel(r'$\varepsilon = -\mathrm{d}E/\mathrm{d}t$')
    h, lab = a.get_legend_handles_labels()
    h += [Line2D([], [], color='k', lw=fs.LW), Line2D([], [], color='k', lw=fs.LW - 1, ls='--')]
    lab += [r'$\varepsilon$', 'shock viscosity']
    fs.legend(a, handles=h, labels=lab, loc='upper right')
    fs.save(fig, TEST, 'turb_eps')


def turb_spec(runs, keys=('1.00', '3.00')):
    """E(k) = E_sqrt(rho)u(k) + E_B(k) (shell sums) at t = 1 and 3, N = 64, 128, 256; dotted lines at k = N/3."""
    S = load(runs, 'stage6q', 'turb', 'summary.json')['runs']
    fig, ax = fs.panels(1, 2, sharey=True, wspace=0.08, top=0.8)
    for j, key in enumerate(keys):
        a = ax[0, j]
        E = {}
        for N, name in sorted(TURB.items()):
            s = S[name]['spec'][key]
            E[N] = np.array(s['E_tot'])
            k = np.arange(len(E[N]))
            a.loglog(k[1:], E[N][1:], color=COLOR[N], lw=fs.LW, label=r'$N = %d$' % N)
            a.axvline(N / 3.0, color=COLOR[N], lw=1.5, ls=':')
            print('t = %.4f tau0 N=%-3d sum E(k) %.5f' % (s['t'] / np.pi, N, E[N].sum()))
        k = np.arange(1, len(E[128]))
        bad = np.abs(np.log(E[128][k] / E[256][k])) > 0.1
        print('t = %s tau0: N = 128 and 256 within 10%% for k < %d (first k off)' % (key, k[bad][0]))
        a.set_title(r'$t = %d\,\tau_0$' % round(float(key)), fontsize=fs.FONTSIZE)
        a.set_xlabel(r'$k$')
        a.set_xlim(1, 200)
    ax[0, 0].set_ylim(3e-9, 20.0)                      # no label at the corner
    ax[0, 0].set_ylabel(r'$E(k)$')
    fs.legend(ax[0, 0], loc='lower left')
    fs.save(fig, TEST, 'turb_spec')


def turb_rho_pdf(runs):
    """Volume PDF of s = ln rho at t = 0.5 and 1, N = 64, 128, 256 (Stage 7); dashed: the Gaussian with the mean
    and the standard deviation of s of N = 256."""
    from analyze_stage7 import PDF_EDGES
    sc = 0.5 * (PDF_EDGES[1:] + PDF_EDGES[:-1])
    T = {N: load(runs, 'stage7', 'tail_N%d.json' % N) for N in TURB}
    fig, ax = fs.panels(1, 2, sharey=True, wspace=0.08, top=0.8)
    for j, tt in enumerate(TAIL_T):
        a = ax[0, j]
        for N in sorted(TURB):
            s = next(x for x in T[N]['snaps'] if x['t_req'] == tt)
            pdf = np.array(s['pdf'])
            a.semilogy(sc[pdf > 0], pdf[pdf > 0], color=COLOR[N], lw=fs.LW, label=r'$N = %d$' % N)
            print('t = %.4f tau0 N=%-3d sigma_s %.4f mean_s %.4f (-sigma^2/2 %.4f) rho_min %.5f F_V(rho < 0.1) %.2e '
                  'PDF outside [%g, %g]: %s, sum PDF ds %.6f'
                  % (s['t_tau'], N, s['sigma_s'], s['mean_s'], -0.5 * s['sigma_s'] ** 2, s['rho_min'],
                     s['F_V']['0.1'], PDF_EDGES[0], PDF_EDGES[-1], s['pdf_out'], pdf.sum() * 0.1))
        g = T[256]['snaps'][[x['t_req'] for x in T[256]['snaps']].index(tt)]
        x = np.linspace(-4.0, 2.0, 400)
        a.plot(x, np.exp(-0.5 * ((x - g['mean_s']) / g['sigma_s']) ** 2) / (np.sqrt(2 * np.pi) * g['sigma_s']),
               label='Gaussian', **fs.GUIDE)
        a.set_title(r'$t = %g\,\tau_0$' % tt, fontsize=fs.FONTSIZE)
        a.set_xlabel(r'$s = \ln\rho$')
        a.set_xlim(-4.0, 2.0)
        a.set_xticks((-3, -2, -1, 0, 1))                # none at the edges between panels
    ax[0, 0].set_ylim(1e-7, 5.0)
    ax[0, 0].set_ylabel(r'PDF of $s$')
    fs.legend(ax[0, 0], loc='upper left')
    fs.save(fig, TEST, 'turb_rho_pdf')


def stage15(runs):
    """(spec of the 128^3 runs {label: spec}, ES-256, judge.json['ms1'], analyze_stage15)."""
    import analyze_stage15 as a15
    sp = load(runs, 'stage15', 'spec.json')
    R = {v['info']['label']: v for v in sp.values()
         if v['info']['problem'] == 'turb1' and v['info']['N'] == 128 and v['info']['role'] == 'b'}
    return R, sp[a15.ES256], load(runs, 'stage15', 'judge.json')['ms1'], a15


def turb_K_spec(runs):
    """Window means of E_B(n) and E_m(n) (m = rho u) of the 128^3 runs ES, H(8) (B*), L(1), L(2) and of ES-256
    (grey), windows W1 = [0.75, 1.25] and W2 = [1.5, 2] (Stage 15); dotted: k_c = 128/3."""
    R, ref, J, a15 = stage15(runs)
    fig, ax = fs.panels(2, 2, sharex=True, sharey=True, wspace=0.08, hspace=0.12, top=0.8)
    for i, q in enumerate(a15.QJ):
        for j, W in enumerate(a15.WIN):
            a = ax[i, j]
            E = np.array(ref['win'][W]['mean'][q])
            a.loglog(np.arange(1, len(E)), E[1:], label=r'ES, $N = 256$', **fs.REF)
            for lab in KRUNS:
                E = np.array(R[lab]['win'][W]['mean'][q])
                a.loglog(np.arange(1, len(E)), E[1:], color=KCOLOR[lab], lw=fs.LW - 1, label=scheme(lab))
            a.axvline(128 / 3.0, color='k', lw=1.5, ls=':')
            if i == 0:
                a.set_title(r'$W_%d = [%g, %g]\,\tau_0$' % (j + 1, *a15.WIN[W]), fontsize=fs.FONTSIZE)
            if i == 1:
                a.set_xlabel(r'$n$')
        ax[i, 0].set_ylabel(r'$E_{%s}(n)$' % q[2:])
    ax[0, 0].set_xlim(1, 128)
    ax[0, 0].set_ylim(2e-8, 2.0)                       # no label at the corner
    fs.legend(ax[0, 0], loc='lower left')
    fs.save(fig, TEST, 'turb_K_spec')


def turb_K(runs):
    """ln T(n), T = E_128(n)/E_ES-256(n), of E_B and E_m in W1 and W2 (Stage 15, M_s = 1); grey band: |ln T| <=
    ln 1.25; circles: k_x of ES and B* = H(8); dotted: k_c = 128/3. The kappa = k_x/k_c are recomputed with
    analyze_stage15.pair and checked against judge.json."""
    R, ref, J, a15 = stage15(runs)
    kc = 128 / 3.0
    thr = np.log(a15.T_TOL)
    fig, ax = fs.panels(2, 2, sharex=True, sharey=True, wspace=0.08, hspace=0.12, top=0.8)
    for i, q in enumerate(a15.QJ):
        for j, W in enumerate(a15.WIN):
            a = ax[i, j]
            a.axhspan(-thr, thr, color=fs.GRAY, alpha=0.3, lw=0)
            key = '%s %s ES256' % (q, W)
            for lab in KRUNS:
                E, Er = np.array(R[lab]['win'][W]['mean'][q]), np.array(ref['win'][W]['mean'][q])
                n = np.arange(1, min(len(E), len(Er)))
                y = np.log(E[n] / Er[n])
                a.semilogx(n, y, color=KCOLOR[lab], lw=fs.LW - 1, label=scheme(lab))
                p = a15.pair(R[lab], ref, q, W, a15.CAP_K, a15.D_HI, kc)
                pj = J['pairs'][lab][key]
                print('%-5s %s kappa %.4f (judge %.4f) k_x %.2f censored %s D %.4f (judge %.4f)'
                      % (lab, key, p['kappa'], pj['kappa'], p['k'], p['cens'], p['D'], pj['D']))
                if lab in ('ES', J['B']):
                    a.plot(p['k'], np.sign(np.interp(p['k'], n, y)) * thr, 'o', color=KCOLOR[lab], ms=fs.MS - 3,
                           mfc='w', mew=fs.MEW, zorder=5)
            a.axvline(kc, color='k', lw=1.5, ls=':')
            if i == 0:
                a.set_title(r'$W_%d = [%g, %g]\,\tau_0$' % (j + 1, *a15.WIN[W]), fontsize=fs.FONTSIZE)
            if i == 1:
                a.set_xlabel(r'$n$')
            fs.tag(a, r'$E_{%s}$' % q[2:], x=0.05, y=0.93)
        ax[i, 0].set_ylabel(r'$\ln T(n)$')
    ax[0, 0].set_xlim(1, 64)
    ax[0, 0].set_ylim(-1.0, 1.0)
    ax[0, 0].set_yticks((-0.5, 0.0, 0.5))
    fs.legend(ax[1, 0], loc='lower left', ncol=2)
    print('B* = %s, verdict %r (reference ES-256 only), counts %s' % (J['B'], J['verdict'], J['counts']))
    fs.save(fig, TEST, 'turb_K')


def turb_kappa(runs):
    """kappa = k_x/k_c of each record of ES and B* = H(8) against the record of ES-256 at the same time (Stage 16,
    time-matched, post hoc); bars: the kappa of the window."""
    import analyze_stage15 as a15
    A = load(runs, 'stage16', 'judge_tm.json')['A']
    fig, ax = fs.panels(2, 1, sharex=True, hspace=0.12)
    for i, q in enumerate(a15.QJ):
        a = ax[i, 0]
        for lab in ('ES', A['B']):
            for j, W in enumerate(a15.WIN):
                p = A['pairs'][lab]['%s %s ES256' % (q, W)]
                a.plot(p['t'], p['kappas'], 'o-', color=KCOLOR[lab], lw=fs.LW - 1, ms=fs.MS - 5,
                       label=scheme(lab) if j == 0 else None)
                a.hlines(p['kappa'], *a15.WIN[W], color=KCOLOR[lab], lw=8, alpha=0.4)
                print('%-5s %s %s kappa %.4f censored %s D %.4f sigma %.4f records %d (t %.3f-%.3f)'
                      % (lab, q, W, p['kappa'], p['cens'], p['D'], p['sigma'], len(p['t']), p['t'][0], p['t'][-1]))
        a.set_ylim(0.0, 1.0)
        a.set_ylabel(r'$\kappa$')
        fs.tag(a, r'$E_{%s}$' % q[2:], x=0.04, y=0.93)
    ax[1, 0].set_xlim(0.7, 2.05)
    ax[1, 0].set_xlabel(r'$t/\tau_0$')
    fs.legend(ax[0, 0], loc='lower right')
    print('B* = %s, verdict %r, counts %s, kbar %s, cost %s' % (A['B'], A['verdict'], A['counts'], A['kbar'],
                                                                 A['cost']))
    fs.save(fig, TEST, 'turb_kappa')


def turb_coef(runs):
    """rho_min(t) of the 128^3 runs of Stage 8 (svv_eps0 = epsilon_0 and C_shock doubled and halved; black: the
    control) and ln(rho_min/rho_min of the control); grey: the window of L, bars: R = L - L_control (L: the time
    mean of ln rho_min over the window; e1.0_C2 and e0.5_C1 overlap at R = -0.08), dotted: +-delta of R."""
    from analyze_stage8 import mean_log
    D = load(runs, 'stage8', 'rhomin.json')
    rows = load(runs, 'stage8', 'summary.json')['rows']
    fig, ax = fs.panels(1, 2, wspace=0.3)
    c = D['e0.5_C2']['series']
    tc, lc = np.array(c['t_tau']), np.log(np.array(c['rho_min']))
    r0 = D['ref']['series']
    print('control against the reference stage6q/track/turb1_N128_t: max |ln rho_min diff| %.2e'
          % np.abs(np.interp(tc, r0['t_tau'], np.log(r0['rho_min'])) - lc).max())
    t1 = min([COEF_WIN[1]] + [D[n]['series']['t_tau'][-1] for n, _, _ in COEF])    # analyze_stage8.windows
    Lc = mean_log(tc, np.exp(lc), COEF_WIN[0], t1)
    for name, col, ls in COEF:
        e, C = name[1:].split('_C')
        lab = r'$\epsilon_0 = %s$, $C_\mathrm{shock} = %s$' % (e, C)
        s = D[name]['series']
        t, lr = np.array(s['t_tau']), np.log(np.array(s['rho_min']))
        ax[0, 0].semilogy(t, np.exp(lr), color=col, ls=ls, lw=fs.LW - (0 if name == 'e0.5_C2' else 1), label=lab)
        w = rows[name]
        if name != 'e0.5_C2':
            m = tc <= t[-1]
            ax[0, 1].plot(tc[m], np.interp(tc[m], t, lr) - lc[m], color=col, ls=ls, lw=fs.LW - 1)
            ax[0, 1].hlines(w['R'], *COEF_WIN, color=col, lw=8, alpha=0.6, zorder=4)
        L = mean_log(t, np.exp(lr), COEF_WIN[0], t1)
        print('%-9s m %.5f at t %.4f (rho_min.json %.5f at %.4f) L %.4f (summary %.4f) R %.4f (summary %.4f) '
              't_final %.4f' % (name, w['m'], w['t_m'], D[name]['m'], D[name]['t_m'], L, w['L'], L - Lc, w['R'], t[-1]))
    for a in ax[0]:
        a.axvspan(*COEF_WIN, color=fs.GRAY, alpha=0.25, lw=0)
        a.set_xlim(0.0, 1.2)
        a.set_xlabel(r'$t/\tau_0$')
    for s in (1, -1):
        ax[0, 1].axhline(s * COEF_DELTA, color='k', lw=1.5, ls=':')
    ax[0, 0].set_ylim(0.03, 1.2)
    fs.logticks(ax[0, 0].yaxis, (0.05, 0.1, 0.2, 0.5, 1))
    ax[0, 0].set_ylabel(r'$\rho_\mathrm{min}$')
    ax[0, 1].set_ylim(-0.4, 0.4)
    ax[0, 1].set_yticks((-0.2, 0.0, 0.2))
    ax[0, 1].set_ylabel(r'$\ln(\rho_\mathrm{min}/\rho_\mathrm{min,control})$')
    fs.legend(ax[0, 0], loc='lower left', fontsize=fs.LEGEND_FONTSIZE - 2)
    fs.save(fig, TEST, 'turb_coef')


def turb_nomass(runs):
    """svv_mass = F against the control (Stage 13, 128^3): rho_min(t), and dt against rho_min with a rho^1.5
    guide."""
    D = load(runs, 'stage13d', 'rhomin.json')
    fig, ax = fs.panels(1, 2, wspace=0.3)
    for k, col, lab in (('C', 'k', 'control'), ('N', fs.COLORS[5], r'\texttt{svv\_mass} = F')):
        s, e = D[k]['series'], D[k]['extras']['series']
        t, r = np.array(s['t_tau']), np.array(s['rho_min'])
        te, dt = np.array(e['t_tau']), np.array(e['dt'])
        ri = np.interp(te, t, r)
        ax[0, 0].semilogy(t, r, color=col, lw=fs.LW, label=lab)
        ax[0, 1].loglog(ri, dt, color=col, lw=fs.LW, label=lab)
        m = ri < 1e-2
        slope = np.polyfit(np.log(ri[m]), np.log(dt[m]), 1)[0] if m.sum() > 10 else np.nan
        print('%s t_final %.4f tau0 rho_min %.3e dt [%.3e, %.3e] steps %d; d ln dt/d ln rho_min (rho_min < 1e-2) %.4f'
              % (k, t[-1], r.min(), dt.min(), dt.max(), len(dt), slope))
    fs.slope(ax[0, 1], 1e-7, 1e-4, 3e-10, 1.5, label=r'$\rho_\mathrm{min}^{3/2}$')     # 100 times above the data
    ax[0, 0].set_xlim(0.0, 1.2)
    ax[0, 0].set_xlabel(r'$t/\tau_0$')
    ax[0, 0].set_ylabel(r'$\rho_\mathrm{min}$')
    ax[0, 1].set_xlabel(r'$\rho_\mathrm{min}$')
    ax[0, 1].set_ylabel(r'$\Delta t$')
    fs.legend(ax[0, 0], loc='lower left')
    fs.save(fig, TEST, 'turb_nomass')


def slice_file():
    return os.path.join('figdata', 'Turb3D', '%s_%s_tau1.npz' % (os.path.basename(SLICE[0]), SLICE[1]))


def cfl_file(case):
    return os.path.join('figdata', 'Turb3D', 'cfl_%s.npz' % os.path.basename(case))


def turb_slice(runs):
    """ln rho and the shock-viscosity coefficient zeta = C_shock Delta^2 S[max(0, -div u)] on the plane z = 0 of
    the 256^3 run at t = tau0; zeta below ZETA_LIM[0] (zeros included) in the lowest colour."""
    from matplotlib.colors import LogNorm
    with np.load(os.path.join(runs, slice_file())) as f:
        p = {k: f[k].astype(np.float64) for k in f.files}
    r = load(runs, 'stage6q', 'turb', 'summary.json')['runs'][os.path.basename(SLICE[0])]
    s, z = r['series'], r['zeta']['z0']['1.00']
    i = int(np.argmin(np.abs(np.array(s['t']) - p['t'])))
    rho, zeta = p['rho'], p['zeta']
    q = np.sort((rho * zeta * p['D'] ** 2).ravel())[::-1]
    a90 = (np.searchsorted(np.cumsum(q) / q.sum(), 0.9) + 1) / q.size
    print('t = %.6f tau0 (summary %.6f); rho [%.4f, %.4f] in the box [%.4f, %.4f] at t = %.4f tau0' % (
        p['t'] / np.pi, z['t'] / np.pi, rho.min(), rho.max(), s['rho_min'][i], s['rho_max'][i], s['t'][i] / np.pi))
    print('zeta_max %.6e (summary %.6e, box %.4e); area of zeta > 0.01 zeta_max %.5f (summary %.5f), a90 %.5f '
          '(summary %.5f); zeta = 0 on %.4f of the plane' % (zeta.max(), z['zeta_max'], s['zeta_max'][i],
                                                              (zeta > 0.01 * zeta.max()).mean(), z['frac'], a90,
                                                              z['a90'], (zeta == 0).mean()))
    fig, ax = fs.maps(1, 2, labels=('x', 'y'), gap=(2.6, 1.1))
    im = fs.image(ax[0, 0], np.log(rho), cmap='RdBu_r', vmin=-LNRHO_LIM, vmax=LNRHO_LIM)
    fs.colorbar(fig, ax[0, 0], im, r'$\ln\rho$')
    cmap = fs.parula().with_extremes(under=fs.parula()(0.0), bad=fs.parula()(0.0))
    cmap.colorbar_extend = 'min'
    im = fs.image(ax[0, 1], np.ma.masked_less_equal(zeta, 0.0), cmap=cmap, norm=LogNorm(*ZETA_LIM))
    fs.colorbar(fig, ax[0, 1], im, r'$\zeta$')
    fs.save(fig, TEST, 'turb_slice')


def turb_track(runs):
    """dt/dt_cfl (dt_cfl = min(dt_adv, dt_visc) of cfl.dat, every step) of turb1 at N = 64 and 128: the Stage 4
    runs (fixed dt and 'decrement') and the 'track' runs of Stage 6; dotted: track_safety and track_raise."""
    fig, ax = fs.panels(1, 2, sharey=True, wspace=0.08, top=0.8)
    for j, (N, cases) in enumerate(sorted(TRACK.items())):
        a = ax[0, j]
        for (case, lab), col in zip(cases, (fs.COLORS[4], fs.COLORS[5])):
            with np.load(os.path.join(runs, cfl_file(case))) as f:
                c = np.vstack(sorted((f['part%d' % i] for i in range(int(f['n']))), key=lambda x: x[0, 0]))
            t, r = c[:, 0] / np.pi, c[:, 1] / np.minimum(c[:, 2], c[:, 3])
            a.plot(t, r, color=col, lw=fs.LW - 1, label='%s (%s steps)' % (lab, '{:,}'.format(len(c))))
            print('N=%-3d %-26s steps %6d t_last %.4f tau0 max dt/dt_cfl %.4f mean dt/dt_cfl over [j, j+1) tau0 %s'
                  % (N, case, len(c), t[-1], r.max(), ' '.join('%.3f' % r[(t >= k) & (t < k + 1)].mean()
                                                               for k in range(3))))
        for y in (TRACK_SAFETY, TRACK_RAISE):
            a.axhline(y, color='k', lw=1.5, ls=':')
        a.set_title(r'$N = %d$' % N, fontsize=fs.FONTSIZE)
        a.set_xlabel(r'$t/\tau_0$')
        a.set_xlim(0.0, 3.0)
        fs.legend(a, loc='lower left')
    ax[0, 1].set_xticks((0.5, 1.0, 1.5, 2.0, 2.5, 3.0))         # no 0 under the 3 of the left panel
    ax[0, 0].set_ylim(0.0, 1.05)
    ax[0, 0].set_yticks((0.2, 0.4, 0.6, 0.8, 1.0))
    ax[0, 0].set_ylabel(r'$\Delta t/\Delta t_\mathrm{cfl}$')
    fs.save(fig, TEST, 'turb_track')


def turb_ms(runs):
    """M_s = 1-4 at 64^3 (Stage 5, 12; the seed and delta B of turb1, the amplitude of u scaled) and M_s = 4 at
    128^3 (dashed, Stage 5, 6): rho_min and dt/dt0 against t/tau0; dotted: rho_min = 1e-3 and dt0/10 (the
    criteria H and S); dt0: the first dt of each run. M_s = 1 runs on to 3 tau0 with dt = dt0."""
    M = load(runs, 'stage5', 'ana', 'ms.json')
    runs4 = [(k, v, '-', r'$M_s = %d$' % float(k)) for k, v in sorted(M.items())]
    runs4.append(('4', load(runs, 'stage5', 'ana', 'res.json')['turb4b']['128'], '--', r'$M_s = 4$, $N = 128$'))
    fig, ax = fs.panels(1, 2, wspace=0.3)
    for k, v, ls, lab in runs4:
        s, tau0 = v['series'], M[k]['tau0']
        t, dt = np.array(s['t']) / tau0, np.array(s['dt'])
        col = fs.SERIES[int(k) - 1]
        ax[0, 0].semilogy(t, s['rho_min'], color=col, ls=ls, lw=fs.LW, label=lab)
        ax[0, 1].semilogy(t, dt / dt[0], color=col, ls=ls, lw=fs.LW)
        h, i = v['health'], int(np.argmin(s['rho_min']))
        fall = np.nonzero(dt < dt[0] / 10)[0]
        print('M_s=%s %-24s dt0 %.4e t_last %.3f tau0 steps %d min dt/dt0 %.1e min rho_min %.3e at %.3f tau0; '
              'dt < dt0/10 from %s' % (k, v['case'], dt[0], h['t_last'] / tau0, h['steps'], dt.min() / dt[0],
                                      s['rho_min'][i], t[i], '%.3f tau0 (rho_min %.3e)' % (t[fall[0]],
                                      s['rho_min'][fall[0]]) if len(fall) else '-'))
        if ls == '-':
            print('    ms.json: min dt/dt0 %.1e rho_min %.3e at %.3f tau0 F_max %.3g fall %s verdict %s'
                  % (v['dt_min_rel'], v['rho_min'], v['t_rho_min'], v['F_max'], v['fall'], v['verdict']))
    ax[0, 0].axhline(1e-3, color='k', lw=1.5, ls=':')
    ax[0, 1].axhline(0.1, color='k', lw=1.5, ls=':')
    for a in ax[0]:
        a.set_xlim(0.0, 1.0)
        a.set_xlabel(r'$t/\tau_0$')
    ax[0, 0].set_ylim(5e-5, 2.0)
    ax[0, 0].set_ylabel(r'$\rho_\mathrm{min}$')
    ax[0, 1].set_ylim(5e-5, 2.0)
    ax[0, 1].set_ylabel(r'$\Delta t/\Delta t_0$')
    h, lab = ax[0, 0].get_legend_handles_labels()
    fs.legend(ax[0, 1], handles=h, labels=lab, loc='lower right', ncol=2)     # under the ends of the curves
    fs.save(fig, TEST, 'turb_ms')


FIGS = dict(turb_eps=turb_eps, turb_spec=turb_spec, turb_slice=turb_slice, turb_track=turb_track,
            turb_rho_pdf=turb_rho_pdf, turb_K_spec=turb_K_spec, turb_K=turb_K, turb_kappa=turb_kappa,
            turb_coef=turb_coef, turb_nomass=turb_nomass, turb_ms=turb_ms)
FETCH = ([fd.Extract(slice_file(), 'plane', 'runs/' + SLICE[0], SLICE[1], '%.12f' % SLICE[2],
                     '%d,%d,%d' % ((SLICE[3],) * 3), 'rho,zeta,D')]
         + [fd.Extract(cfl_file(c), 'table', 'runs/' + c, 'cfl.dat') for N in sorted(TRACK) for c, _ in TRACK[N]]
         + [fd.Copy('stage5/ana/res.json', 'stage5/ana/ms.json')])

if __name__ == '__main__':
    fd.main(TEST, FIGS, FETCH)
