#!/usr/bin/env python3
"""Style of the paper figures drawn by the plot_<test>.py scripts of the test folders, after plot_flow.py of the
gallope strong-scaling figures (~/Dropbox/Codes/scaling/gallope/flow/strong scaling/202610): LaTeX text in Times
(newtx, the preamble of diagnostics/MHD_COMP_ISOTH/latex_preamble.py), 30 pt, axes and ticks 2 pt wide, ticks
inside on all four sides with minor ticks, the Okabe-Ito colours, legends without frame, PDF output.

    import figstyle as fs
    plt = fs.setup()
    fig, ax = fs.panels(1, 2)
    ...
    fs.save(fig, test_dir, 'name')      -> <test_dir>/fig/name.pdf

The rendering needs latex, newtx and mathrsfs (TeX Live has them; ps1 does not).
"""
import os
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
DIAG = os.path.normpath(os.path.join(HERE, '..', '..', '..', 'diagnostics', 'MHD_COMP_ISOTH'))

LINEWIDTH = 2
FONTSIZE = 30
LEGEND_FONTSIZE = 22
LW = 3                                   # data lines
MS, MEW = 15, 3                          # markers (the errorbar of plot_flow.py)
GUIDE = dict(color='k', lw=1.5, ls='--')     # reference lines (slopes, theory)
# Okabe-Ito, in the order of plot_flow.py
COLORS = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#999999']
GRAY = '#999999'
# series of lines (resolutions, cases) in this fixed order: Okabe-Ito without the yellow, too light for lines
SERIES = [COLORS[i] for i in (0, 1, 2, 4, 5, 6)]
REF = dict(color=GRAY, lw=6, solid_capstyle='butt', zorder=1)     # reference solutions (FV), under the data
EXACT = dict(color='k', lw=1.5, ls='--', zorder=3)                 # exact or ideal solutions, over the data
PANEL = (9.5, 6.6)                       # one panel of the figures, inches (plot_flow.py: 9.5 wide)
MARGIN = dict(left=1.5, right=0.6, bottom=1.1, top=0.35)     # inches around the grid of panels


def setup():
    """matplotlib.pyplot with the rc settings of plot_flow.py (Agg backend)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from cycler import cycler
    sys.path.insert(0, DIAG)
    from latex_preamble import preamble
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('font', serif='Times')
    plt.rc('font', size=FONTSIZE)
    plt.rc('axes', linewidth=LINEWIDTH)
    plt.rc('axes', labelsize=FONTSIZE)
    plt.rc('legend', fontsize=LEGEND_FONTSIZE)
    for ax in ('xtick', 'ytick'):
        plt.rc(ax, labelsize=FONTSIZE)
        plt.rc(ax, direction='in')
        plt.rc(ax + '.major', width=LINEWIDTH)
        plt.rc(ax + '.major', size=17)
        plt.rc(ax + '.minor', width=LINEWIDTH)
        plt.rc(ax + '.minor', visible=True)
        plt.rc(ax + '.minor', size=8)
    plt.rc('xtick', top=True)
    plt.rc('ytick', right=True)
    plt.rcParams['axes.prop_cycle'] = cycler(color=COLORS)
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['text.latex.preamble'] = preamble
    plt.rcParams['savefig.dpi'] = 300            # rasterized maps
    return plt


def parula():
    sys.path.insert(0, DIAG)
    from parula import parula_map
    return parula_map


def panels(nrow, ncol, w=PANEL[0], h=PANEL[1], sharex=False, sharey=False, wspace=0.3, hspace=0.3, **margin):
    """(fig, ax[nrow, ncol]) on a gridspec with fixed margins (no tight_layout, as plot_flow.py). The figure is
    w*ncol by h*nrow inches; the margins are MARGIN (inches) unless given."""
    import matplotlib.pyplot as plt
    W, H = w * ncol, h * nrow
    m = dict(MARGIN)
    m.update(margin)
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(nrow, ncol, wspace=wspace, hspace=hspace, left=m['left'] / W, right=1.0 - m['right'] / W,
                          bottom=m['bottom'] / H, top=1.0 - m['top'] / H)
    ax = gs.subplots(sharex=sharex, sharey=sharey, squeeze=False)
    return fig, ax


def legend(ax, **kw):
    kw.setdefault('frameon', False)
    kw.setdefault('fontsize', LEGEND_FONTSIZE)
    kw.setdefault('handlelength', 1.6)
    return ax.legend(**kw)


def tag(ax, s, x=0.04, y=0.92, **kw):
    """Panel label in the upper left corner, e.g. tag(ax, r'(a) Sod')."""
    kw.setdefault('fontsize', LEGEND_FONTSIZE + 2)
    kw.setdefault('va', 'top')
    return ax.text(x, y, s, transform=ax.transAxes, **kw)


def slope(ax, x0, x1, y0, p, label=None, **kw):
    """Guide line y = y0 (x/x0)^p on [x0, x1] (log axes), labelled at its end."""
    import numpy as np
    k = dict(GUIDE)
    k.update(kw)
    x = np.array([x0, x1])
    ax.plot(x, y0 * (x / x0) ** p, **k)
    if label:
        ax.text(x1 * 1.08, y0 * (x1 / x0) ** p, label, fontsize=LEGEND_FONTSIZE, va='center')


def logticks(axis, values):
    """Ticks at values with plain labels (1, 0.5, 512) on a logarithmic axis (ax.xaxis or ax.yaxis), the minor
    ticks unlabelled; matplotlib otherwise labels one tick per decade only, or the minor ticks as 2x10^0."""
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
    axis.set_major_locator(FixedLocator(values))
    axis.set_major_formatter(FixedFormatter(['%g' % v for v in values]))
    axis.set_minor_formatter(NullFormatter())


def maps(nrow, ncol, labels=('x', 'z'), side=6.0, gap=(0.9, 1.1), margin=(1.5, 2.6, 1.1, 0.9)):
    """(fig, ax[nrow, ncol]) of square map axes, side inches, with gaps (between columns, rows) and margins (left,
    right, bottom, top) in inches; the coordinates labels (horizontal, vertical) run from 0 to 2 pi, tick labels on
    the outer axes only, one 0 at the corner and the vertical label upright."""
    import numpy as np
    import matplotlib.pyplot as plt
    lm, rm, bm, tm = margin
    W = lm + ncol * side + (ncol - 1) * gap[0] + rm
    H = bm + nrow * side + (nrow - 1) * gap[1] + tm
    fig = plt.figure(figsize=(W, H))
    ax = np.empty((nrow, ncol), dtype=object)
    for i in range(nrow):
        for j in range(ncol):
            a = fig.add_axes([(lm + j * (side + gap[0])) / W, (bm + (nrow - 1 - i) * (side + gap[1])) / H,
                              side / W, side / H])
            for s, lab in (('x', r'$0$'), ('y', '')):
                getattr(a, 'set_%sticks' % s)([0.0, np.pi, 2.0 * np.pi])
                getattr(a, 'set_%sticklabels' % s)([lab, r'$\pi$', r'$2\pi$'])
            a.minorticks_off()
            a.tick_params(labelbottom=(i == nrow - 1), labelleft=(j == 0))
            if i == nrow - 1:
                a.set_xlabel(r'$%s$' % labels[0])
            if j == 0:
                a.set_ylabel(r'$%s$' % labels[1], rotation=0, va='center', labelpad=16)
            ax[i, j] = a
    return fig, ax


def image(a, q, **kw):
    """q[i, j] (i along the horizontal coordinate) on a map axes of maps()."""
    import numpy as np
    return a.imshow(q.T, origin='lower', extent=(0.0, 2.0 * np.pi, 0.0, 2.0 * np.pi), interpolation='none', **kw)


def colorbar(fig, axs, im, label, pad=0.25, width=0.3):
    """Colour bar right of the axes axs (one or a column), the label above it."""
    import numpy as np
    W = fig.get_figwidth()
    p = [a.get_position() for a in np.atleast_1d(axs)]
    x1, y0, y1 = max(q.x1 for q in p), min(q.y0 for q in p), max(q.y1 for q in p)
    cax = fig.add_axes([x1 + pad / W, y0, width / W, y1 - y0])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.minorticks_off()
    cb.ax.tick_params(length=8)
    cb.ax.set_title(label, fontsize=FONTSIZE, pad=22)     # clear of the top tick label
    return cb


def save(fig, test_dir, name):
    """<test_dir>/fig/<name>.pdf, written to a temporary name first."""
    import matplotlib.pyplot as plt
    d = os.path.join(test_dir, 'fig')
    os.makedirs(d, exist_ok=True)
    out = os.path.join(d, name + '.pdf')
    tmp = os.path.join(d, '.%s.tmp.pdf' % name)
    fig.savefig(tmp)
    os.replace(tmp, out)
    plt.close(fig)
    print('wrote', os.path.relpath(out))
    return out
