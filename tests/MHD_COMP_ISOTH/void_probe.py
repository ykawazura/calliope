#!/usr/bin/env python3
"""State around the density minimum and the fastest wave of a 3D run, read from restart/ (or ic/).

    python3 tests/MHD_COMP_ISOTH/void_probe.py runs/stage4/turb/turb4_N64_a [restart|ic]

Prints rho, |u|, |B|, v_A, beta and p_th + p_B at the cell of min rho and at the cell of max_d(|u_d| + c_f,d),
the profile of rho along each axis through the min of rho, and the fraction of cells below rho thresholds."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_ic  # noqa: E402
from analyze_stage1 import read_nml, num  # noqa: E402


def probe(case, sub='restart'):
    d = read_nml(os.path.join(case, 'calliope.in'))
    g = gen_ic.Grid(int(d['nlx']), int(d['nly']), int(d['nlz']),
                    num(d.get('lx', '1.0')), num(d.get('ly', '1.0')), num(d.get('lz', '1.0')))
    cs2 = num(d.get('beta0', '1.0'))/2  # params.F90: cs2va2 = beta0/2
    lk, mk, bk, tt = gen_ic.read_restart(os.path.join(case, sub), g)
    rho, u, b = gen_ic.real_from_state(g, lk, mk, bk)
    u2 = sum(x**2 for x in u)
    b2 = sum(x**2 for x in b)
    va2 = b2/rho
    # max over d of |u_d| + c_f,d (advance.F90)
    cmax = np.max([np.abs(u[i]) + np.sqrt(0.5*(cs2 + va2 + np.sqrt(np.maximum(0, (cs2 + va2)**2 - 4*cs2*b[i]**2/rho))))
                   for i in range(3)], axis=0)
    pth, pb = rho*cs2, 0.5*b2
    print('%s/%s  t = %.6f  N = %s  cs^2 = %g' % (case, sub, tt, g.nl, cs2))
    print('  mean p_th %.4g  mean p_B %.4g  max p_th + p_B %.4g' % (pth.mean(), pb.mean(), (pth + pb).max()))
    for lab, idx in (('min rho', np.unravel_index(np.argmin(rho), rho.shape)),
                     ('max |u_d| + c_f,d', np.unravel_index(np.argmax(cmax), rho.shape))):
        print('  at %-18s %s: rho %.3e |u| %.3g |m| %.3g |B| %.3g v_A %.3g |u_d| + c_f,d %.3g beta %.3g p_th + p_B %.4g'
              % (lab, tuple(int(i) for i in idx), rho[idx], np.sqrt(u2[idx]), rho[idx]*np.sqrt(u2[idx]), np.sqrt(b2[idx]),
                 np.sqrt(va2[idx]), cmax[idx], 2*pth[idx]/b2[idx], pth[idx] + pb[idx]))
    i0 = np.unravel_index(np.argmin(rho), rho.shape)
    for ax, name in ((0, 'x'), (1, 'z'), (2, 'y')):
        for lab, f in (('rho', rho), ('|B|', np.sqrt(b2))):
            sl = list(i0)
            vals = []
            for o in range(-4, 5):
                sl[ax] = (i0[ax] + o) % rho.shape[ax]
                vals.append(f[tuple(sl)])
            print('  %-3s along %s through the min (offsets -4..4): %s' % (lab, name, ' '.join('%.2e' % v for v in vals)))
    print('  power fraction of s outside the 2/3 mask: %.2e' % (np.sum(np.abs(lk)**2*(1 - g.mask))/np.sum(np.abs(lk)**2)))
    for th in (1e-1, 1e-2, 1e-3):
        print('  fraction of cells with rho < %g: %.2e  (mass %.2e)' % (th, (rho < th).mean(), rho[rho < th].sum()/rho.sum()))


if __name__ == '__main__':
    probe(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'restart')
