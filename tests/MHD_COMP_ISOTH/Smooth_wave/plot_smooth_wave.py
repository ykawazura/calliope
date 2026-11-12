#!/usr/bin/env python3
"""Figures of the circularly polarized Alfven wave and the oblique planar problems (README.md) in the style of the
paper (../stages/figstyle.py).

    python3 tests/MHD_COMP_ISOTH/Smooth_wave/plot_smooth_wave.py fetch                   inputs from ps1
    python3 tests/MHD_COMP_ISOTH/Smooth_wave/plot_smooth_wave.py [plot] [name ...] [--runs runs]
                                                                                -> Smooth_wave/fig/<name>.pdf

The inputs are the summaries of Stage 4 (runs/stage4/{alfven,cpdiag,planar}/summary.json: analyze_stage4.py, and
alfven/cp_rates.json: the largest growth rate of the physical system on the grid, cp_stability.py), which fetch
copies from ps1 (0.5 MB), and the last line of the planar runs phy_<prob>_k12_N256_{2d,1d,1dh}, which fetch extracts
on ps1 with analyze_stage4.planar_snaps (figdata.py planar, 0.3 MB); the FV reference of C2 is in the local runs/
(runs/stage3/fvref). Every figure prints the numbers to check against REPORT_stage4 (sections 5.2-5.4, 6.3).
Growth rates are in units of w0 = |k0| v_A (B0 = rho0 = 1), so that a mode grows by exp(2 pi gamma) per period.
"""
import os
import sys
import json
import numpy as np
TEST = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(TEST, '..', 'stages'))
import figstyle as fs
import figdata as fd
import stage3 as s3
import gen_stage4 as gs

plt = fs.setup()
K0 = ((1, 0, 1), (1, 0, 2), (1, 1, 2))
KLAB = {k: r'$\bm{k}_0 = (%d, %d, %d)$' % k for k in K0}
KNAME = {(1, 0, 1): '2d_k101', (1, 0, 2): '2d_k102', (1, 1, 2): '3d_k112'}
KMARK = dict(zip(K0, 'osD'))
PLANAR = ('F1', 'C2', 'TDS_w1', 'RDpc')
# planar_profiles: variable (index of W), its label, the place of the panel tag (axes fraction)
PROFILE = dict(F1=(0, r'$\rho$', (0.62, 0.92)), C2=(4, r'$B_{t1}$', (0.35, 0.92)),
               TDS_w1=(0, r'$\rho$', (0.33, 0.6)), RDpc=(4, r'$B_{t1}$', (0.32, 0.6)))


def load(runs, batch, name='summary.json'):
    with open(os.path.join(runs, 'stage4', batch, name)) as f:
        return json.load(f)


def cases(A, k0, amp, setting, N):
    """{name: case} of the alfven runs of (k0, amplitude, setting, N), all n_p and rank counts."""
    pre = '%s_a%s_N%d_%s_p' % (KNAME[k0], '01' if amp == 0.1 else '1', N, setting)
    return {n: c for n, c in A['cases'].items() if n.startswith(pre)}


def series(c, key):
    m = sorted(int(j) for j in c['err'] if int(j) > 0)
    return np.array(m), np.array([c['err'][str(j)][key] for j in m])


def rk3_loss(n_p):
    """1 - |R|^n_p, the amplitude lost in one period by RK3 on an oscillation of n_p steps (R = 1 + z + z^2/2 + z^3/6,
    z = 2 pi i/n_p); (2 pi)^4/(24 n_p^3) to leading order, the phase error being smaller by 2 pi/n_p."""
    z = 2j * np.pi / np.asarray(n_p, float)
    return 1.0 - np.abs(1.0 + z + z ** 2 / 2.0 + z ** 3 / 6.0) ** n_p


def cp_error(runs):
    """(a) Amplitude 0.1, ideal, N = 16: max|B - B_ex|/A at t = m T (filled) against m, the same of the RK3 solution of
    the eigenmode (B_pred, dashed) and the difference from that solution (open: round-off). (b) Amplitude 1: the
    difference max|(u, B) - RK3 solution|/A against m for k0 = (1, 0, 2), N = 64 and (1, 1, 2), N = 32, ideal and
    default, with exp(2 pi gamma m) of the rate fitted to the ideal run (the truncation instability, section 5.3)."""
    A = load(runs, 'alfven')
    fig, ax = fs.panels(1, 2, wspace=0.32)
    a = ax[0, 0]
    for i, (k0, n_p) in enumerate(zip(K0, (400, 200, 200))):
        c = A['cases']['%s_a01_N16_ideal_p%d' % (KNAME[k0], n_p)]
        m, e = series(c, 'B')
        _, ep = series(c, 'B_pred')
        _, pr = series(c, 'pred')
        mk = dict(color=fs.SERIES[i], ms=fs.MS * 0.7, mew=fs.MEW)
        a.loglog(m, e, KMARK[k0], label=r'%s, $n_p = %d$' % (KLAB[k0], n_p), **mk)
        a.loglog(m, ep, **fs.EXACT)
        a.loglog(m, pr, KMARK[k0], mfc='none', **mk)
        print('A=0.1 %-28s err_B(1, 10) %.3e %.3e  B_pred %.3e %.3e  |err/B_pred - 1| <= %.1e  pred <= %.1e'
              % ('%s_p%d' % (KNAME[k0], n_p), e[0], e[-1], ep[0], ep[-1], np.max(np.abs(e / ep - 1)), pr.max()))
    a.plot([], [], label='RK3 of the eigenmode', **fs.EXACT)
    a.plot([], [], 'o', color=fs.GRAY, ms=fs.MS * 0.7, mew=fs.MEW, mfc='none', label='difference from it')
    a.set_xlim(0.8, 12)
    fs.logticks(a.xaxis, (1, 2, 5, 10))
    a.set_ylim(1e-16, 1e-1)
    a.set_yticks([1e-15, 1e-12, 1e-9, 1e-6, 1e-3])
    a.set_xlabel(r'$t/T$')
    a.set_ylabel(r'$\max|\bm{B} - \bm{B}_\mathrm{ex}|/A$')
    fs.tag(a, r'(a) $A = 0.1$, ideal, $N = 16$')
    fs.legend(a, loc='center left', bbox_to_anchor=(0.0, 0.40))
    a = ax[0, 1]
    for i, (k0, N, sfx) in enumerate((((1, 0, 2), 64, ''), ((1, 1, 2), 32, '_np8'))):
        for setting, mk, mfc in (('ideal', 'o-', None), ('default', 's-', 'white')):
            name = '%s_a1_N%d_%s_p400%s' % (KNAME[k0], N, setting, sfx)
            c = A['cases'][name]
            m, pr = series(c, 'pred')
            a.semilogy(m, pr, mk, color=fs.SERIES[i], lw=fs.LW, ms=fs.MS * 0.7, mew=fs.MEW, mfc=mfc,
                       label=r'%s, $N = %d$' % (KLAB[k0], N) if setting == 'ideal' else None)
            g = c['growth']
            print('A=1  %-30s %-8s pred(1) %.1e pred(end) %.1e at m = %d  gamma %s galerkin %.4f'
                  % (name, c['status'], pr[0], pr[-1], m[-1], '-' if g['gamma'] is None else '%.4f' % g['gamma'],
                     g['galerkin']))
            if setting == 'ideal':
                mm = np.array([m[-1] - 4, m[-1]])
                y = 8 * pr[-1] * np.exp(2 * np.pi * g['gamma'] * (mm - m[-1]))
                a.semilogy(mm, y, **fs.GUIDE)
                if i == 0:
                    a.text(mm[1] + 0.25, y[1], r'$\mathrm{e}^{2\pi\gamma t/T}$', fontsize=fs.LEGEND_FONTSIZE,
                           va='center')
    a.plot([], [], 'o-', color='k', lw=fs.LW, ms=fs.MS * 0.7, mew=fs.MEW, label='ideal')
    a.plot([], [], 's-', color='k', lw=fs.LW, ms=fs.MS * 0.7, mew=fs.MEW, mfc='white', label='default')
    a.set_xlim(0, 10.5)
    a.set_ylim(1e-16, 1e2)
    a.set_yticks([1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1])
    a.set_xlabel(r'$t/T$')
    a.set_ylabel(r'$\max|(\bm{u}, \bm{B}) - \mathrm{RK3}|/A$')
    fs.tag(a, r'(b) $A = 1$')
    fs.legend(a, loc='upper left', bbox_to_anchor=(0.0, 0.88))
    fs.save(fig, TEST, 'cp_error')


def cp_dt(runs):
    """max|B - B_ex|/A after one period against the steps per period n_p, default setting, amplitude 0.1 and 1, with
    the amplitude lost by RK3 (rk3_loss): the error is that of the time integration alone, independent of k0 and of
    the amplitude; in 3D it is 0.91 of it, about (1 - kx^2/k0^2)^1/2 = 0.913, since the error is the largest of the
    components. The ideal runs agree to four digits."""
    A = load(runs, 'alfven')
    fig, ax = fs.panels(1, 1)
    a = ax[0, 0]
    for i, k0 in enumerate(K0):
        for amp, mfc in ((0.1, None), (1, 'white')):
            for setting in ('default', 'ideal'):
                key = next(k for k in A['checks'] if k.startswith('dt3 %s_a%s_N16_%s_p'
                                                                  % (KNAME[k0], '01' if amp == 0.1 else '1', setting)))
                d = A['checks'][key]
                n_p = int(key.rsplit('_p', 1)[1]) * np.array([1, 2, 4])
                print('%-36s err(T) %s /rk3_loss %s ratios(T) %s ratios(10T) %s'
                      % (key[4:], ' '.join('%.3e' % v for v in d['err']['1']),
                         ' '.join('%.4f' % v for v in np.array(d['err']['1']) / rk3_loss(n_p)),
                         ' '.join('%.3f' % v for v in d['ratio']['1']), ' '.join('%.2f' % v for v in d['ratio']['10'])))
                if setting == 'default':
                    a.loglog(n_p, d['err']['1'], KMARK[k0], color=fs.SERIES[i], ms=fs.MS * (1.0 if amp == 1 else 0.45),
                             mew=fs.MEW, mfc=mfc, zorder=4 - amp)
    n_p = np.geomspace(170, 3800, 50)
    a.loglog(n_p, rk3_loss(n_p), **fs.EXACT)
    h = [a.plot([], [], KMARK[k0], color=fs.SERIES[i], ms=fs.MS * 0.6, mew=fs.MEW, label=KLAB[k0])[0]
         for i, k0 in enumerate(K0)]
    h += [a.plot([], [], 'o', color='k', ms=fs.MS * 0.45, mew=fs.MEW, label=r'$A = 0.1$')[0],
          a.plot([], [], 'o', color='k', ms=fs.MS * 1.0, mew=fs.MEW, mfc='white', label=r'$A = 1$')[0],
          a.plot([], [], label=r'RK3, $1 - |R|^{n_p}$', **fs.EXACT)[0]]
    a.set_xlim(150, 4500)
    fs.logticks(a.xaxis, (200, 400, 800, 1600, 3200))
    a.set_ylim(5e-10, 1e-4)
    a.set_xlabel(r'$n_p = T/\Delta t$')
    a.set_ylabel(r'$\max|\bm{B} - \bm{B}_\mathrm{ex}|/A$ at $t = T$')
    fs.legend(a, handles=h, loc='upper right', ncol=2, columnspacing=1.0)
    fs.save(fig, TEST, 'cp_dt')


def cp_growth(runs):
    """Growth rate gamma (10 periods, amplitude 1) against N for the three k0: ideal (circles) on the largest Galerkin
    rate (dashed: the truncation instability), default (squares; triangles on the bottom: no measurable growth), and
    the largest rate of the physical system on the grid with the bound 1.25 of it for default at N >= 32 (grey)."""
    A = load(runs, 'alfven')
    R = load(runs, 'alfven', 'cp_rates.json')
    fig, ax = fs.panels(1, 3, w=7.0, sharey=True, wspace=0.08, top=1.0)
    floor = 6e-3
    for j, k0 in enumerate(K0):
        a = ax[0, j]
        Ns = (16, 32, 64) if k0[1] == 0 else (16, 32)
        gal, gid = [], []
        for N in Ns:
            for setting in ('ideal', 'default'):
                cs = cases(A, k0, 1, setting, N)
                n_p = min(int(n.split('_p')[1].split('_')[0]) for n in cs)
                for n, c in sorted(cs.items()):
                    if '_p%d' % n_p not in n:
                        continue
                    g = c['growth']
                    print('%-36s %-8s gamma %s galerkin %.4f%s' % (
                        n, c['status'], '-' if g['gamma'] is None else '%.4f' % g['gamma'], g['galerkin'],
                        '  physical_max %.4f ok %s' % (g['physical_max']['rate'], g['ok'])
                        if 'physical_max' in g else ''))
                    if setting == 'ideal':
                        gal.append(g['galerkin'])
                        gid.append(g['gamma'])
                    elif g['gamma'] is None:
                        a.semilogy(N, floor, 'v', color=fs.SERIES[1], ms=fs.MS * 0.7, mew=fs.MEW, mfc='white')
                    else:
                        a.semilogy(N, g['gamma'], 's', color=fs.SERIES[1], ms=fs.MS * 0.8, mew=fs.MEW,
                                   mfc='white')
        a.semilogy(Ns, gal, **fs.EXACT)
        a.semilogy(Ns, gid, 'o', color=fs.SERIES[0], ms=fs.MS * 0.8, mew=fs.MEW, zorder=2)
        for key, r in R.items():
            if key.startswith(str(k0)):
                N = int(key.split('N=')[1].split()[0])
                x = N * np.array([0.85, 1.18])
                a.semilogy(x, [r['rate']] * 2, **fs.REF)
                a.semilogy(x, [1.25 * r['rate']] * 2, color=fs.GRAY, lw=1.5, ls=':')
                print('%s N=%d physical max on the grid %.4f (q %s), bound 1.25x = %.4f'
                      % (k0, N, r['rate'], r['q'], 1.25 * r['rate']))
        a.set_xscale('log')
        a.set_xlim(12, 85)
        fs.logticks(a.xaxis, (16, 32, 64))
        a.set_xlabel(r'$N$')
        fs.tag(a, KLAB[k0], x=0.06, y=0.95)
    a = ax[0, 0]
    a.set_ylim(4e-3, 1.5)
    fs.logticks(a.yaxis, (0.01, 0.03, 0.1, 0.3, 1))
    a.set_ylabel(r'$\gamma/\omega_0$')
    h = [a.plot([], [], 'o', color=fs.SERIES[0], ms=fs.MS * 0.8, mew=fs.MEW, label='ideal')[0],
         a.plot([], [], label='Galerkin, largest', **fs.EXACT)[0],
         a.plot([], [], 's', color=fs.SERIES[1], ms=fs.MS * 0.8, mew=fs.MEW, mfc='white', label='default')[0],
         a.plot([], [], 'v', color=fs.SERIES[1], ms=fs.MS * 0.7, mew=fs.MEW, mfc='white',
                label='default, no growth')[0],
         a.plot([], [], label='physical, largest', **fs.REF)[0],
         a.plot([], [], color=fs.GRAY, lw=1.5, ls=':', label=r'$1.25\times$ physical')[0]]
    fig.legend(handles=h, loc='upper center', ncol=6, frameon=False, fontsize=fs.LEGEND_FONTSIZE, handlelength=1.6,
               columnspacing=1.4, bbox_to_anchor=(0.5, 1.0))
    fs.save(fig, TEST, 'cp_growth')


def cp_diag(runs):
    """Amplitude 1, 10 periods: gamma over the Galerkin rate of the dominant family with the same dissipation, with
    the dissipation terms one at a time (cpdiag, section 5.4). Adding the zeta term changes nothing (noeta = svvonly,
    nosvv = etaonly: the wave has div u = 0). k0 = (1, 0, 1), N = 16 has no measurable growth with SVV."""
    D = load(runs, 'cpdiag')['runs']
    groups = (('2d_k102_a1_N16', '_p400', r'$\bm{k}_0 = (1, 0, 2)$, $N = 16$'),
              ('3d_k112_a1_N16', '_p400', r'$\bm{k}_0 = (1, 1, 2)$, $N = 16$'),
              ('3d_k112_a1_N32', '_p400_np8', r'$\bm{k}_0 = (1, 1, 2)$, $N = 32$'))
    cats = (('ideal',), ('svvonly', 'noeta'), ('etaonly', 'nosvv'), ('default',))
    fig, ax = fs.panels(1, 1, w=11.0)
    a = ax[0, 0]
    for i, (g, sfx, lab) in enumerate(groups):
        for x, cat in enumerate(cats):
            r = [D['%s_%s%s' % (g, s, sfx)]['growth'] for s in cat]
            print('%-16s %-16s gamma %s galerkin %s ratio %.2f' % (
                g, '='.join(cat), ' '.join('%.4f' % q['gamma'] for q in r), ' '.join('%.4f' % q['galerkin'] for q in r),
                r[0]['gamma'] / r[0]['galerkin']))
            a.plot(x + 0.12 * (i - 1), r[0]['gamma'] / r[0]['galerkin'], 'osD'[i], color=fs.SERIES[i], ms=fs.MS * 0.8,
                   mew=fs.MEW, label=lab if x == 0 else None)
    for g in ('2d_k101_a1_N16',):
        print('%s: gamma %s' % (g, {s: D['%s_%s_p800' % (g, s)]['growth']['gamma']
                                    for s in ('ideal', 'noeta', 'etaonly', 'default')}))
    a.axhline(1.0, **fs.GUIDE)
    a.set_xticks(range(4))
    a.set_xticklabels(('ideal', 'SVV', r'$\eta_\mathrm{shock}$', 'default'))
    a.tick_params(axis='x', which='minor', bottom=False, top=False)
    a.set_xlim(-0.5, 3.5)
    a.set_ylim(0, 8)
    a.set_ylabel(r'$\gamma/\gamma_\mathrm{Galerkin}$')
    fs.legend(a, loc='upper left')
    fs.save(fig, TEST, 'cp_diag')


def planar_file(prob, dim, k='12', N=256):
    return os.path.join('figdata', 'Smooth_wave', 'phy_%s_k%s_N%d_%s.npz' % (prob, k, N, dim))


def fv_C2(runs, t):
    """(x, W) of the finest FV reference of C2 at the time nearest to t (stage3.fv_reference)."""
    for n in (65536, 32768):
        fn = os.path.join(runs, 'stage3', 'fvref', 'C2_N%d.npz' % n)
        if os.path.exists(fn):
            with np.load(fn) as f:
                k = int(np.argmin(np.abs(f['t'] - t)))
                assert abs(f['t'][k] - t) < 1e-9, (fn, t)
                return s3.grid(n), f['W'][k]
    raise FileNotFoundError(os.path.join(runs, 'stage3', 'fvref', 'C2_N65536.npz'))


def planar_profiles(runs):
    """F1, C2, TDS_w1, RDpc at t_end: the 2D run with k = (1, 2), N = 256 along the normal (open circles), the 1D runs
    on N/2 points (the bandwidth of the 2D run, the reference of the criterion) and on N points, with the exact
    solution where it is known (dashed; F1: outside the fan) or the FV reference of C2 (grey). Prints the 6-variable
    L1 against summary.json."""
    S = load(runs, 'planar')['checks']
    fig, ax = fs.panels(2, 2, hspace=0.25)
    for p, a in zip(PLANAR, ax.flat):
        iv, lab, (tx, ty) = PROFILE[p]
        c3 = gs.s3_case(p, 256)
        row = S['phy phy_%s_k12_N256' % p]
        l1 = {}
        for dim in ('2d', '1dh', '1d'):
            with np.load(os.path.join(runs, planar_file(p, dim))) as f:
                t, W = float(f['t']), f['W']
            N = W.shape[1]
            x = s3.grid(N)
            c = gs.planar_bw_case(c3, (1, 2)) if dim == '1dh' else c3
            We, valid = s3.exact(c, x, t)
            if We is None:
                xf, Wf = fv_C2(runs, t)
                We, valid = Wf[:, ::len(xf) // N], np.ones(N, bool)
            l1[dim] = sum(np.sum(np.abs(W[q] - We[q])[valid]) for q in range(6)) * s3.L / N
            if dim == '2d':
                a.plot(x, W[iv], 'o', color=fs.SERIES[0], ms=fs.MS * 0.4, mew=fs.MEW * 0.5, mfc='none', zorder=2,
                       label=r'2D, $N = %d$' % N)
            else:
                a.plot(x, W[iv], color=fs.SERIES[1 if dim == '1dh' else 2], lw=fs.LW * (0.8 if dim == '1dh' else 0.6),
                       zorder=4 if dim == '1dh' else 3, label=r'1D, $N = %d$' % N)
        if p == 'C2':
            xf, Wf = fv_C2(runs, t)
            a.plot(xf[::16], Wf[iv, ::16], label=r'FV, $N = %d$' % len(xf), **fs.REF)
        else:
            x = np.linspace(0.0, s3.L, 4097)
            We, valid = s3.exact(c3, x, t)
            a.plot(x, np.where(valid, We[iv], np.nan), label='exact', **fs.EXACT)
        print('%-7s t=%.4f L1 2d %.6f (summary %.6f) 1dh %.6f (%.6f) 1d %.6f (%.6f)  2d/1dh %.3f (%.3f) 2d/1d %.2f'
              % (p, t, l1['2d'], row['L1_2d'], l1['1dh'], row['L1_c1dh'], l1['1d'], row['L1_c1d'],
                 l1['2d'] / l1['1dh'], row['r_2d_ref'], l1['2d'] / l1['1d']))
        a.set_xlim(0.0, s3.L)
        a.set_ylabel(lab)
        fs.tag(a, r'%s, $t = %g$' % (p.replace('_', r'\_'), round(t, 6)), x=tx, y=ty)
    for a in ax[-1]:
        a.set_xlabel(r'$x_n$')
    for j in range(2):
        fig.align_ylabels(ax[:, j])
    fs.legend(ax[0, 0], loc='upper right', bbox_to_anchor=(1.0, 0.8))
    fs.legend(ax[0, 1], handles=ax[0, 1].get_lines()[-1:], loc='upper center', bbox_to_anchor=(0.5, 0.75))
    fs.save(fig, TEST, 'planar_profiles')


FIGS = dict(cp_error=cp_error, cp_dt=cp_dt, cp_growth=cp_growth, cp_diag=cp_diag, planar_profiles=planar_profiles)
FETCH = [fd.Copy('stage4/alfven/summary.json', 'stage4/alfven/cp_rates.json', 'stage4/cpdiag/summary.json',
                 'stage4/planar/summary.json')]
FETCH += [fd.Extract(planar_file(p, d), 'planar', 'stage4/planar/phy_%s_k12_N256_%s' % (p, d))
          for p in PLANAR for d in ('2d', '1d', '1dh')]

if __name__ == '__main__':
    fd.main(TEST, FIGS, FETCH)
