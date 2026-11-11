#!/usr/bin/env python3
"""Figures of the 2D isothermal Orszag-Tang vortex (README.md) in the style of the paper (../stages/figstyle.py).

    python3 tests/MHD_COMP_ISOTH/OT2D/plot_ot2d.py fetch                                  planes from ps1
    python3 tests/MHD_COMP_ISOTH/OT2D/plot_ot2d.py [plot] [name ...] [--runs runs]       -> OT2D/fig/<name>.pdf

ot_eps and ot_spec read the local runs/stage14/ot/summary.json (analyze_stage4.py ot runs/stage14). ot_rho and ot_J
read the planes y0 of the C_eta = 0 runs at t = 3 and 6 (Mignone 2007, t = 1.5 and 3), which fetch extracts from the
outputs on ps1 into runs/figdata/OT2D/ (rho, B_x, B_z in float32, 20 MB). Every figure prints the numbers to check
against REPORT_stage4 (section 7) and REPORT_stage14 (section 6); ot_J checks max|J_y| and the range of rho against
the 'late' entries of summary.json, computed on ps1 from the same records in float64.
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
CASES = {128: 'ot_N128_eta0', 256: 'ot_N256_eta0', 512: 'ot_N512_eta0', 1024: 'ot_N1024_eta0_p256'}
COLOR = {N: fs.SERIES[i] for i, N in enumerate(sorted(CASES))}          # the colour follows N in every figure
PLANES = ((3, (256, 512, 1024)), (6, (128, 256, 512)))                 # (t, N) of the maps
J_LIM = 50.0                                                            # colour range of J_y (the sheets saturate)
RHO_LIM = (0.15, 6.5)                                                   # colour range of rho (all of ot_rho)


def summary(runs):
    with open(os.path.join(runs, 'stage14', 'ot', 'summary.json')) as f:
        return json.load(f)


def by_eta(S, eta):
    """{N: run} of the runs with C_eta = eta."""
    return {r['info']['N']: r for r in S['runs'].values() if r['info']['C_eta'] == eta}


def plane_file(N, t):
    return os.path.join('figdata', 'OT2D', '%s_t%d.npz' % (CASES[N], t))


def plane(runs, N, t):
    with np.load(os.path.join(runs, plane_file(N, t))) as f:
        return {k: f[k] for k in f.files}


def current(f):
    """J_y = d_z B_x - d_x B_z on the (x, z) plane, spectrally (analyze_stage4.ot_late)."""
    N = f['bx'].shape[0]
    k = np.fft.fftfreq(N, 1.0 / N)
    KX, KZ = np.meshgrid(k, k, indexing='ij')
    bx, bz = (np.fft.fft2(f[q].astype(np.float64)) for q in ('bx', 'bz'))
    return np.real(np.fft.ifft2(1j * KZ * bx - 1j * KX * bz))


def ot_eps(runs):
    """eps(t) = -dE/dt for C_eta = 0 (N = 128-1024; N = 1024 to t = pi) and C_eta = 2 (N = 128-512)."""
    S = summary(runs)
    fig, ax = fs.panels(1, 2, sharey=True, wspace=0.08, top=0.8)
    for j, eta in enumerate((0.0, 2.0)):
        a = ax[0, j]
        for N, r in sorted(by_eta(S, eta).items()):
            a.plot(r['series']['t'], r['series']['eps'], color=COLOR[N], lw=fs.LW, label=r'$N = %d$' % N)
            p = r['eps_peak']
            print('C_eta=%g N=%-4d eps peak %.4f at t = %.4f, t_last %.4f' % (eta, N, p['eps'], p['t'], r['t_last']))
        a.set_title(r'$C_\eta = %d$' % eta, fontsize=fs.FONTSIZE)
        a.set_xlabel(r'$t$')
        a.set_xlim(0.0, 2.0 * np.pi)
        a.set_xticks(range(7))
    ax[0, 0].set_ylim(0.0, 32.0)
    ax[0, 0].set_yticks((10, 20, 30))                   # no 0 on the y axis, under the 0 of t
    ax[0, 0].set_ylabel(r'$\varepsilon = -\mathrm{d}E/\mathrm{d}t$')
    fs.legend(ax[0, 0], loc='upper right')
    fs.save(fig, TEST, 'ot_eps')


def ot_spec(runs, keys=('1.0000', '3.0000')):
    """E(k) = E_sqrt(rho)u(k) + E_B(k) (shell sums) at t = 1 and 3, C_eta = 0, N = 128-1024; dotted lines at k = N/3
    and a k^-2 guide."""
    S = summary(runs)
    R = by_eta(S, 0.0)
    fig, ax = fs.panels(1, 2, sharey=True, wspace=0.08, top=0.8)
    for j, key in enumerate(keys):
        a = ax[0, j]
        for N, r in sorted(R.items()):
            E = np.array(r['spec'][key]['E_tot'])
            k = np.arange(len(E))
            a.loglog(k[1:], E[1:], color=COLOR[N], lw=fs.LW, label=r'$N = %d$' % N)
            a.axvline(N / 3.0, color=COLOR[N], lw=1.5, ls=':')
            print('t = %.4f N=%-4d sum E(k) %.5f' % (r['spec'][key]['t'], N, E.sum()))
        e1, e2 = (np.array(R[N]['spec'][key]['E_tot']) for N in (512, 1024))
        k = np.arange(1, len(e1))
        bad = np.abs(np.log(e1[k] / e2[k])) > 0.1
        print('t = %s: N = 512 and 1024 within 10%% for k < %d (first k off)' % (key, k[bad][0]))
        a.set_title(r'$t = %d$' % round(float(key)), fontsize=fs.FONTSIZE)
        a.set_xlabel(r'$k$')
        a.set_xlim(1, 600)
    fs.slope(ax[0, 1], 10, 200, 1.5e-2, -2, label=r'$k^{-2}$')      # above the spectra
    ax[0, 0].set_ylim(1e-10, 2.0)
    ax[0, 0].set_ylabel(r'$E(k)$')
    fs.legend(ax[0, 0], loc='lower left')
    fs.save(fig, TEST, 'ot_spec')


def check(S, N, t, f, J=None):
    """max|J_y| and the range of rho of a plane against the 'late' entry of summary.json."""
    key = next(k for k, v in S['late'].items() if v['N'] == N and abs(v['t'] - float(f['t'])) < 1e-9)
    v = S['late'][key]
    msg = 'N=%-4d t = %.4f rho [%.5f, %.5f] (late [%.5f, %.5f])' % (N, f['t'], f['rho'].min(), f['rho'].max(),
                                                                      v['rho_min'], v['rho_max'])
    if J is not None:
        msg += ' max|J_y| %.3f (late %.3f, rel. diff %.1e)' % (np.abs(J).max(), v['J_max'],
                                                               abs(np.abs(J).max() / v['J_max'] - 1.0))
    print(msg)


def ot_rho(runs, N=512):
    """rho on the (x, z) plane at t = 3 and 6 (Mignone 2007, t = 1.5 and 3), N = 512, C_eta = 0, on a common
    logarithmic colour scale (at t = 6, rho = 6.3 at one spot in the centre: 99 cells above 3 within 13 x 9 cells)."""
    from matplotlib.colors import LogNorm
    S = summary(runs)
    fig, ax = fs.maps(1, 2)
    norm = LogNorm(RHO_LIM[0], RHO_LIM[1])
    for j, (t, _) in enumerate(PLANES):
        f = plane(runs, N, t)
        im = fs.image(ax[0, j], f['rho'], cmap=fs.parula(), norm=norm)
        ax[0, j].set_title(r'$t = %d$' % t, fontsize=fs.FONTSIZE)
        check(S, N, t, f)
    cb = fs.colorbar(fig, ax[0, -1], im, r'$\rho$')
    cb.set_ticks((0.2, 0.5, 1, 2, 5))
    cb.set_ticklabels(('0.2', '0.5', '1', '2', '5'))
    fs.save(fig, TEST, 'ot_rho')


def ot_J(runs):
    """J_y = d_z B_x - d_x B_z on the (x, z) plane, C_eta = 0: t = 3 for N = 256, 512, 1024 (top) and t = 6 for
    N = 128, 256, 512 (bottom); the colours saturate at |J_y| = J_LIM."""
    S = summary(runs)
    fig, ax = fs.maps(len(PLANES), 3)
    for i, (t, Ns) in enumerate(PLANES):
        for j, N in enumerate(Ns):
            f = plane(runs, N, t)
            J = current(f)
            im = fs.image(ax[i, j], J, cmap='RdBu_r', vmin=-J_LIM, vmax=J_LIM)
            ax[i, j].set_title(r'$N = %d$, $t = %d$' % (N, t), fontsize=fs.LEGEND_FONTSIZE + 2)
            check(S, N, t, f, J)
    fs.colorbar(fig, ax[:, -1], im, r'$J_y$')
    fs.save(fig, TEST, 'ot_J')


FIGS = dict(ot_eps=ot_eps, ot_spec=ot_spec, ot_rho=ot_rho, ot_J=ot_J)
FETCH = [fd.Extract(plane_file(N, t), 'plane', 'runs/stage14/ot/' + CASES[N], 'y0', str(t), '%d,2,%d' % (N, N),
                    'rho,bx,bz') for t, Ns in PLANES for N in Ns]

if __name__ == '__main__':
    fd.main(TEST, FIGS, FETCH)
