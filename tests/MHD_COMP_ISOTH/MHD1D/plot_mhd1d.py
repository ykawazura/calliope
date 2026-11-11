#!/usr/bin/env python3
"""Figures of the 1D MHD Riemann problems (README.md) in the style of the paper (../stages/figstyle.py).

    python3 tests/MHD_COMP_ISOTH/MHD1D/plot_mhd1d.py [plot] [name ...] [--runs runs]    -> MHD1D/fig/<name>.pdf

The inputs are in the local runs/ (no fetch): the final states of the Stage 14 runs (runs/stage14/mhd1d/suite_last,
eta_suite_last: analyze_stage3.py last), the FV solutions of Stage 3 (runs/stage3/fvref, fvwin) and the table of
Stage 16 (runs/stage16/conv_table.json: analyze_stage16.py table). Every figure prints the numbers to check against
REPORT_stage3 (sections 13, 15) and REPORT_stage16 (section 3).
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
import mhd_exact as mx
import stage3_aux as aux

plt = fs.setup()
VAR = {0: r'$\rho$', 1: r'$u_x$', 2: r'$u_y$', 4: r'$B_y$'}
TEND = r'$t = %g$'
LOWER = dict(x=0.04, y=0.06, va='bottom')         # panel tags in the lower left, below the profiles


def last(runs, name, eta=False):
    d = os.path.join(runs, 'stage14', 'mhd1d', 'eta_suite_last' if eta else 'suite_last')
    with np.load(os.path.join(d, name + ('_eta2' if eta else '') + '.npz')) as f:
        return float(f['t']), f['W']


def fvref(runs, prob):
    """(N, t, W (len t, 6, N)) of the finest FV reference of a compound problem (stage3.fv_reference)."""
    for n in (65536, 32768):
        fn = os.path.join(runs, 'stage3', 'fvref', '%s_N%d.npz' % (prob, n))
        if os.path.exists(fn):
            with np.load(fn) as f:
                return n, f['t'], f['W']
    raise FileNotFoundError(os.path.join(runs, 'stage3', 'fvref', '%s_N65536.npz' % prob))


def fvwin(runs, prob):
    with np.load(os.path.join(runs, 'stage3', 'fvwin', aux.FV_CONV[prob])) as f:
        return {k: f[k] for k in f.files}


def shock_x(prob):
    c = dict(fam='coplanar', prob=prob, N=256, C_eta=0.0)
    return s3.cop_interval(prob)[0] + mx.coplanar_case(prob)[0] * s3.t_end(c)


def l6(W, Wr):
    """6-variable L1 (per unit length) of W against Wr sampled at the points of W (stage3.analyze)."""
    return float(sum(np.mean(np.abs(W[i] - Wr[i, ::Wr.shape[1] // W.shape[1]])) for i in range(6)))


def labels(ax, rows, cols):
    for a in ax[-1, :]:
        a.set_xlabel(r'$x$')
    for j, k in enumerate(cols):
        for a in ax[:, j]:
            a.set_ylabel(VAR[k])


def mhd1d_coplanar(runs, Ns=(256, 1024)):
    """F1, F2, S1, S2, I1: rho and B_y near the shock at t_end, N = 256, 1024, with the exact solution, against the
    distance x - x_s from the shock of the exact solution."""
    probs = ('F1', 'F2', 'S1', 'S2', 'I1')
    fig, ax = fs.panels(len(probs), 2, h=4.6, sharex=True, hspace=0.12, wspace=0.28)
    lim = (-0.1, 0.06)
    for r, prob in enumerate(probs):
        xs = shock_x(prob)
        msg = []
        for i, N in enumerate(Ns):
            t, W = last(runs, '%s_N%d' % (prob, N))
            x = s3.grid(N)
            We, valid = s3.exact(dict(fam='coplanar', prob=prob, N=N, C_eta=0.0), x, t)
            msg.append('N=%d %.2e' % (N, np.sum(np.abs(W[0] - We[0])[valid]) / N))
            for j, k in enumerate((0, 4)):
                ax[r, j].plot(x - xs, W[k], color=fs.SERIES[i], lw=fs.LW, label=r'$N = %d$' % N)
        xf = xs + np.linspace(*lim, 4001)
        We, valid = s3.exact(dict(fam='coplanar', prob=prob, N=len(xf), C_eta=0.0), xf, t)
        for j, k in enumerate((0, 4)):
            ax[r, j].plot(xf - xs, np.where(valid, We[k], np.nan), label='exact', **fs.EXACT)
        fs.tag(ax[r, 0], r'%s, %s' % (prob, TEND % t), **LOWER)
        print('%-3s x_s=%.4f rho L1 vs exact (valid points): %s' % (prob, xs, ', '.join(msg)))
    ax[0, 0].set_xlim(*lim)
    labels(ax, probs, (0, 4))
    for a in ax[-1]:
        a.set_xlabel(r'$x - x_\mathrm{s}$')
    fs.legend(ax[0, 1], loc='lower left')
    fs.save(fig, TEST, 'mhd1d_coplanar')


def mhd1d_compound(runs, Ns=(512, 2048)):
    """C1, C3: rho and B_y at t = 0.1, N = 512, 2048, with the FV reference (x < 1/2: the problem of x = 1/4)."""
    fig, ax = fs.panels(2, 2, hspace=0.12, wspace=0.28)
    for r, prob in enumerate(('C1', 'C3')):
        Nf, tf, Wf = fvref(runs, prob)
        xf = s3.grid(Nf)
        for j, k in enumerate((0, 4)):
            ax[r, j].plot(xf, Wf[-1][k], label=r'FV, $N = %d$' % Nf, **fs.REF)
        msg = []
        for i, N in enumerate(Ns):
            t, W = last(runs, '%s_N%d' % (prob, N))
            assert abs(t - tf[-1]) < 1e-9, (prob, N, t)
            msg.append('N=%d %.2e' % (N, l6(W, Wf[-1])))
            for j, k in enumerate((0, 4)):
                ax[r, j].plot(s3.grid(N), W[k], color=fs.SERIES[i], lw=fs.LW, label=r'$N = %d$' % N)
        for a in ax[r]:
            a.set_xlim(0.0, 0.5)
        fs.tag(ax[r, 0], r'%s, %s' % (prob, TEND % tf[-1]), **LOWER)
        print('%s 6-variable L1 vs FV N=%d: %s' % (prob, Nf, ', '.join(msg)))
    for a in ax[0]:
        a.tick_params(labelbottom=False)
    labels(ax, 2, (0, 4))
    fs.legend(ax[0, 0], loc='upper right')
    fs.save(fig, TEST, 'mhd1d_compound')


def mhd1d_rd(runs, N=512):
    """RDh, RDp (w = 0, 8 cells) and RDpc: B_y and |B_t| after one period (t = 1), N = 512, with the exact solution
    (the initial data moved by one period)."""
    rows = (('RDh', ('RDh_w0', 'RDh_w8')), ('RDp', ('RDp_w0', 'RDp_w8')), ('RDpc', ('RDpc',)))
    fig, ax = fs.panels(3, 2, h=5.2, hspace=0.12, wspace=0.3)
    for r, (prob, names) in enumerate(rows):
        for i, nm in enumerate(names):
            t, W = last(runs, '%s_N%d' % (nm, N))
            c = dict(fam='rd', prob=prob[:3], N=N, C_eta=0.0, w=int(nm[-1]) if '_w' in nm else 0, cop=nm == 'RDpc')
            x = s3.grid(N)
            We = s3.exact(c, x, t)[0]
            bt, bte = np.hypot(W[4], W[5]), np.hypot(We[4], We[5])
            lab = r'$w = %s$' % nm.split('_w')[1] if '_w' in nm else r'$N = %d$' % N
            ax[r, 0].plot(x, W[4], color=fs.SERIES[i], lw=fs.LW, label=lab)
            ax[r, 1].plot(x, bt, color=fs.SERIES[i], lw=fs.LW, label=lab)
            ax[r, 0].plot(x, We[4], **fs.EXACT)
            ax[r, 1].plot(x, bte, **fs.EXACT)
            print('%-7s t=%.6f  B_y L1 %.2e  max||B_t|-1| %.2e' % (nm, t, np.mean(np.abs(W[4] - We[4])),
                                                                np.abs(bt - 1.0).max()))
        fs.tag(ax[r, 0], r'%s, %s' % (prob, TEND % round(t, 6)), **LOWER)
        fs.legend(ax[r, 1], loc='lower left')
        for a in ax[r]:
            a.set_xlim(0.0, 0.5)
        if r < 2:
            for a in ax[r]:
                a.tick_params(labelbottom=False)
    for a in ax[-1]:
        a.set_xlabel(r'$x$')
    for a in ax[:, 0]:
        a.set_ylabel(r'$B_y$')
    for a in ax[:, 1]:
        a.set_ylabel(r'$|\mathbf{B}_t|$')
    fs.save(fig, TEST, 'mhd1d_rd')


def mhd1d_reversal(runs, N=1024):
    """I1 (near the shock, t_end) and RDpc (t = 1): rho, u_y, B_y at N = 1024, C_eta = 0 and 2, with the converged
    FV solution and the ideal one (stage3_aux.fig)."""
    fig, ax = fs.panels(2, 3, w=7.2, hspace=0.3, wspace=0.42)
    cols = (0, 2, 4)
    for r, prob in enumerate(('I1', 'RDpc')):
        fv = fvwin(runs, prob)
        Wf = aux.last(fv['W'])
        xf, x = s3.grid(Wf.shape[1]), s3.grid(N)
        if prob == 'I1':
            c = dict(fam='coplanar', prob='I1', N=Wf.shape[1], C_eta=0.0)
            We, valid = s3.exact(c, xf, s3.t_end(c))
            We = np.where(valid, We, np.nan)
            xs = shock_x('I1')
            lim = (xs - 0.12, xs + 0.1)
        else:
            We, lim = fv['We'], (0.0, 0.5)
        for j, k in enumerate(cols):
            a = ax[r, j]
            a.plot(xf, Wf[k], label=r'FV, $N = %d$' % Wf.shape[1], **fs.REF)
            for i, eta in enumerate((False, True)):
                t, W = last(runs, '%s_N%d' % (prob, N), eta)
                a.plot(x, W[k], color=fs.SERIES[i], lw=fs.LW, label=r'$C_\eta = %d$' % (2 if eta else 0))
            a.plot(xf, We[k], label='ideal', **fs.EXACT)
            a.set_xlim(*lim)
        fs.tag(ax[r, 0], r'%s, %s' % (prob, TEND % round(t, 6)), **LOWER)
        for eta in (False, True):
            W = last(runs, '%s_N%d' % (prob, N), eta)[1]
            print('%-4s C_eta=%d rho L1 vs FV N=%d: %.2e' % (prob, 2 if eta else 0, Wf.shape[1], aux.l1(W, Wf)))
    labels(ax, 2, cols)
    for a in ax[0]:
        a.set_xlabel(r'$x$')
    fs.legend(ax[0, 2], loc='lower right')
    fs.save(fig, TEST, 'mhd1d_reversal')


def mhd1d_conv(runs):
    """I1 and RDpc: rho L1 against the converged FV solution (solid) and between N and 2N (dashed) for N = 256-8192,
    C_eta = 0 and 2 (Stage 16, section 3); the dotted line is the distance of the FV solution at N = 4096."""
    with open(os.path.join(runs, 'stage16', 'conv_table.json')) as f:
        T = json.load(f)
    fig, ax = fs.panels(1, 2, sharey=True, wspace=0.08, top=0.8)
    for i, prob in enumerate(('I1', 'RDpc')):
        a = ax[0, i]
        for j, (key, lab) in enumerate(((prob, r'$C_\eta = 0$'), (prob + '_eta2', r'$C_\eta = 2$'))):
            new = T[key]['new']
            N = sorted(int(n) for n in new)
            Nall = [n for n in (256, 512, 1024, 2048, 4096, 8192)]
            fv = [new[str(n)]['fv'] if str(n) in new else np.nan for n in Nall]
            sf = [new[str(n)].get('self', np.nan) if str(n) in new else np.nan for n in Nall]
            sf = [np.nan if v is None else v for v in sf]
            a.loglog(Nall, fv, 'o-', color=fs.SERIES[j], lw=fs.LW, ms=fs.MS * 0.8, mew=fs.MEW, label=lab)
            a.loglog(Nall, sf, 's--', color=fs.SERIES[j], lw=fs.LW, ms=fs.MS * 0.7, mew=fs.MEW, mfc='none')
            print('%-9s N %s\n          fv   %s\n          self %s' % (key, N, ' '.join('%.2e' % v for v in fv),
                                                                   ' '.join('%.2e' % v for v in sf)))
        floor = T[prob + '_fv'][aux.FV_SELF[prob][-1]]
        a.axhline(floor, color='k', lw=1.5, ls=':')
        print('%-9s FV N=4096 vs 16384: %.2e' % (prob, floor))
        a.set_title(prob, fontsize=fs.FONTSIZE)
        a.set_xlabel(r'$N$')
        a.set_xlim(180, 12000)
        fs.logticks(a.xaxis, (256, 512, 1024, 2048, 4096, 8192))
    fs.slope(ax[0, 0], 256, 1024, 6e-3, -1, label=r'$N^{-1}$')        # under the data (shared axes)
    ax[0, 0].set_ylim(6e-4, 1.5e-1)
    ax[0, 0].set_ylabel(r'$L^1$ distance of $\rho$')
    h = ax[0, 0].get_legend_handles_labels()[0]
    h += [plt.Line2D([], [], color='k', lw=fs.LW, marker='o', ms=fs.MS * 0.8, mew=fs.MEW),
          plt.Line2D([], [], color='k', lw=fs.LW, ls='--', marker='s', ms=fs.MS * 0.7, mew=fs.MEW, mfc='none')]
    fs.legend(ax[0, 0], handles=h, labels=[r'$C_\eta = 0$', r'$C_\eta = 2$', 'vs FV', r'$N$ vs $2N$'], ncol=2,
              loc='upper right', columnspacing=1.2)
    fs.save(fig, TEST, 'mhd1d_conv')


def mhd1d_c2(runs, Ns=(256, 1024, 4096)):
    """C2 (B_t reverses): rho and B_y at t = 0.1 for N = 256, 1024, 4096 with the FV reference, C_eta = 0 (top) and
    2 (bottom); x < 1/2 as in mhd1d_compound."""
    Nf, tf, Wf = fvref(runs, 'C2')
    xf = s3.grid(Nf)
    fig, ax = fs.panels(2, 2, hspace=0.12, wspace=0.28)
    for r, eta in enumerate((False, True)):
        for j, k in enumerate((0, 4)):
            ax[r, j].plot(xf, Wf[-1][k], label=r'FV, $N = %d$' % Nf, **fs.REF)
        msg = []
        for i, N in enumerate(Ns):
            t, W = last(runs, 'C2_N%d' % N, eta)
            msg.append('N=%d %.2e' % (N, l6(W, Wf[-1])))
            for j, k in enumerate((0, 4)):
                ax[r, j].plot(s3.grid(N), W[k], color=fs.SERIES[i], lw=fs.LW, label=r'$N = %d$' % N)
        for a in ax[r]:
            a.set_xlim(0.0, 0.5)
        fs.tag(ax[r, 0], r'C2, $C_\eta = %d$' % (2 if eta else 0), **LOWER)
        print('C2 C_eta=%d 6-variable L1 vs FV N=%d: %s' % (2 if eta else 0, Nf, ', '.join(msg)))
    for a in ax[0]:
        a.tick_params(labelbottom=False)
    labels(ax, 2, (0, 4))
    fs.legend(ax[0, 1], loc='lower left')
    fs.save(fig, TEST, 'mhd1d_c2')


FIGS = dict(mhd1d_coplanar=mhd1d_coplanar, mhd1d_compound=mhd1d_compound, mhd1d_rd=mhd1d_rd,
            mhd1d_reversal=mhd1d_reversal, mhd1d_conv=mhd1d_conv, mhd1d_c2=mhd1d_c2)

if __name__ == '__main__':
    fd.main(TEST, FIGS, [])
