#!/usr/bin/env python3
"""Figures of the 1D isothermal Euler tests (README.md) in the style of the paper (../stages/figstyle.py).

    python3 tests/MHD_COMP_ISOTH/Euler1D/plot_euler1d.py fetch                              inputs from ps1
    python3 tests/MHD_COMP_ISOTH/Euler1D/plot_euler1d.py [plot] [name ...] [--runs runs]   -> Euler1D/fig/<name>.pdf

The inputs are the judge table of Stage 1 (runs/stage1_es/judge.json, judge_stage1.py) and the last out1d line of
the shock cases at N = 512 and 2048, which fetch copies from ps1 (2.7 MB); the calliope.in of the cases are in the local
runs/. euler_profiles recomputes L1(rho) and the star states with analyze_stage1.riemann_analysis and prints them
next to judge.json; the other figures print the orders to check against REPORT_stage1 (sections 7.1, 7.3, 7.4).
"""
import os
import sys
import json
import numpy as np
TEST = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(TEST, '..', 'stages'))
import figstyle as fs
import figdata as fd
import analyze_stage1 as a1
import riemann_exact as rx

plt = fs.setup()
PROBS = ('Sod', 'U1', 'U2', 'U4')
LAST = dict(Sod='line_00020.dat', U1='line_00020.dat', U2='line_00021.dat', U4='line_00020.dat')   # t_end
RHO_TICKS = dict(Sod=(1, 2, 5, 10), U1=(0.5, 1, 2), U2=(0.2, 0.5, 1, 2, 5), U4=(0.02, 0.1, 1, 10))
SIMPLE = (('C0_nosvv', r'$C_\mathrm{shock} = 0$, no SVV'), ('C0_svv', r'$C_\mathrm{shock} = 0$'),
          ('C2_svv', r'$C_\mathrm{shock} = 2$'))


def judge(runs):
    with open(os.path.join(runs, 'stage1_es', 'judge.json')) as f:
        return json.load(f)


def shock_case(runs, prob, N=512):
    return os.path.join(runs, 'stage1_es', 'shock', '%s_N%d' % (prob, N))


def fit(N, y):
    """Order p of y ~ N^-p, least squares in log-log (judge_stage1.py)."""
    return -np.polyfit(np.log(N), np.log(y), 1)[0]


def euler_profiles(runs, Ns=(512, 2048)):
    """Sod, U1, U2, U4: rho and u at t_end, N = 512 and 2048, with the exact solution of the two interfaces (x = L/4
    and 3L/4); rho on a logarithmic axis."""
    J = {r['case']: r for r in judge(runs)['shock']['rows']}
    fig, ax = fs.panels(len(PROBS), 2, h=4.6, sharex=True, hspace=0.12, wspace=0.28)
    for r, prob in enumerate(PROBS):
        for i, N in enumerate(Ns):
            case = shock_case(runs, prob, N)
            lines = a1.read_lines(case)
            if not lines:
                raise FileNotFoundError(2, 'no out1d lines', os.path.join(case, 'out1d', LAST[prob]))
            nml = a1.read_nml(os.path.join(case, 'calliope.in'))
            out = {}
            x, rho, u, _, _ = a1.riemann_analysis(nml, lines, out)
            t = out['final']['t']
            ax[r, 0].semilogy(x, rho, color=fs.SERIES[i], lw=fs.LW, label=r'$N = %d$' % N)
            ax[r, 1].plot(x, u, color=fs.SERIES[i], lw=fs.LW, label=r'$N = %d$' % N)
            st = ', '.join('%+.2e' % s['rho_s_rel_err'] for s in out['stars'])
            print('%-3s N=%-4d t=%.6f step %-5d valid %s  L1(rho) %.5f (judge %.5f)  rho* err %s  rho [%.4f, %.4f]'
                  % (prob, N, t, lines[-1][1], out['final']['valid'], out['final']['L1_rho'],
                     J['%s_N%d' % (prob, N)]['L1'], st, rho.min(), rho.max()))
        cs = np.sqrt(a1.num(nml['beta0']) / 2.0)
        args = (a1.num(nml['rho_l']), a1.nums(nml['u_l'])[0], a1.num(nml['rho_r']), a1.nums(nml['u_r'])[0])
        xf = np.linspace(0.0, 2.0 * np.pi, 8001)
        re, ue, valid = rx.periodic_solution(xf, t, *args, 2.0 * np.pi, cs)
        assert valid, prob
        ax[r, 0].semilogy(xf, re, label='exact', **fs.EXACT)
        ax[r, 1].plot(xf, ue, label='exact', **fs.EXACT)
        fs.logticks(ax[r, 0].yaxis, RHO_TICKS[prob])
        fs.tag(ax[r, 0], r'%s, $t = %g$' % (prob, round(t, 6)), x=0.6, ha='center')
    for a in ax[-1]:
        a.set_xlabel(r'$x$')
        a.set_xlim(0.0, 2.0 * np.pi)
        a.set_xticks(np.pi * np.arange(5) / 2.0)
        a.set_xticklabels((r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'))
    for a in ax[:, 0]:
        a.set_ylabel(r'$\rho$')
    for a in ax[:, 1]:
        a.set_ylabel(r'$u$')
    for j in range(2):
        fig.align_ylabels(ax[:, j])
    fs.legend(ax[0, 1], loc='upper right')
    fs.save(fig, TEST, 'euler_profiles')


def euler_l1(runs):
    """L1(rho) against N = 256-4096 for the four problems with an N^-1 guide; the shock criteria apply from N = 512
    (U4: 1024)."""
    C = judge(runs)['shock']['convergence']
    fig, ax = fs.panels(1, 1)
    a = ax[0, 0]
    for i, prob in enumerate(PROBS):
        N, L1 = np.array(C[prob]['N']), np.array(C[prob]['L1'])
        a.loglog(N, L1, 'o-', color=fs.SERIES[i], lw=fs.LW, ms=fs.MS * 0.8, mew=fs.MEW, label=prob)
        print('%-3s L1 %s  order N>=512 %.3f (judge %.3f), N>=1024 %.3f' % (
            prob, ' '.join('%.3e' % v for v in L1), fit(N[N >= 512], L1[N >= 512]), C[prob]['fit_order_512'],
            fit(N[N >= 1024], L1[N >= 1024])))
    fs.slope(a, 512, 2048, 2e-2, -1, label=r'$N^{-1}$')
    a.set_xlim(180, 6500)
    fs.logticks(a.xaxis, (256, 512, 1024, 2048, 4096))
    a.set_ylim(3e-3, 1.5)
    a.set_xlabel(r'$N$')
    a.set_ylabel(r'$L^1$ error of $\rho$')
    fs.legend(a, loc='upper right', ncol=2, columnspacing=1.2)
    fs.save(fig, TEST, 'euler_l1')


def euler_simple(runs):
    """Simple wave before breaking (t = t_b/2): max|u - u_exact| against N = 32-512 for C_shock = 0 without and with
    SVV and C_shock = 2, with an N^-2 guide; C_shock = 0 reaches the floor 2e-11 of the time integration."""
    S = judge(runs)['smooth']['simple']
    fig, ax = fs.panels(1, 1)
    a = ax[0, 0]
    for i, (key, lab) in enumerate(SIMPLE):
        N, e = np.array(S[key]['N']), np.array(S[key]['Linf_u'])
        a.loglog(N, e, 'o-', color=fs.SERIES[i], lw=fs.LW, ms=fs.MS * 0.8, mew=fs.MEW, label=lab)
        print('%-8s Linf(u) %s  orders %s' % (key, ' '.join('%.2e' % v for v in e),
                                              ' '.join('%.2f' % v for v in S[key]['orders'])))
    fs.slope(a, 64, 512, 0.15, -2, label=r'$N^{-2}$')
    a.set_xlim(25, 900)
    fs.logticks(a.xaxis, (32, 64, 128, 256, 512))
    a.set_ylim(5e-12, 1.0)
    a.set_xlabel(r'$N$')
    a.set_ylabel(r'$\max|u - u_\mathrm{exact}|$')
    fs.legend(a, loc='center right')
    fs.save(fig, TEST, 'euler_simple')


def euler_dt(runs):
    """U2, N = 512, t = 0.2: the mass error max|M/M_0 - 1| and max|rho - rho_ref| (the smallest dt the reference)
    against dt (in units of 1e-4), with a dt^3 guide."""
    D = judge(runs)['dtconv']
    dt = np.array(D['dt'])
    x = dt * 1e4
    fig, ax = fs.panels(1, 1)
    a = ax[0, 0]
    a.loglog(x, D['dM'], 'o-', color=fs.SERIES[0], lw=fs.LW, ms=fs.MS * 0.8, mew=fs.MEW,
             label=r'$\max|M/M_0 - 1|$')
    a.loglog(x[:-1], D['d_rho'][:-1], 's-', color=fs.SERIES[1], lw=fs.LW, ms=fs.MS * 0.7, mew=fs.MEW,
             label=r'$\max|\rho - \rho_\mathrm{ref}|$')
    fs.slope(a, 0.31, 2.5, 1.2e-10 * 0.3, 3, label=r'$\Delta t^3$')
    print('dt %s\ndM %s orders %s\nd_rho %s orders %s' % (
        ' '.join('%.4e' % v for v in dt), ' '.join('%.3e' % v for v in D['dM']),
        ' '.join('%.2f' % v for v in D['order_M']), ' '.join('%.2e' % v for v in D['d_rho'][:-1]),
        ' '.join('%.2f' % v for v in D['order_rho'])))
    a.set_xlim(0.2, 5.0)
    a.set_ylim(1e-11, 1e-4)
    fs.logticks(a.xaxis, (0.25, 0.5, 1, 2, 4))
    a.set_xlabel(r'$10^4\,\Delta t$')
    a.set_ylabel('error')
    fs.legend(a, loc='upper left')
    fs.save(fig, TEST, 'euler_dt')


FIGS = dict(euler_profiles=euler_profiles, euler_l1=euler_l1, euler_simple=euler_simple, euler_dt=euler_dt)
FETCH = [fd.Copy('stage1_es/judge.json', *['stage1_es/shock/%s_N%d/out1d/%s' % (p, N, f) for N in (512, 2048)
                                            for p, f in LAST.items()])]

if __name__ == '__main__':
    fd.main(TEST, FIGS, FETCH)
