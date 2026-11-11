#!/usr/bin/env python3
"""The generated inputs (profile1d.dat or ic/) of a case of the per-test folders (tests/MHD_COMP_ISOTH/<Test>), written
again and checked against the md5 of the ps1 runs in its generated.md5.

    python3 tests/MHD_COMP_ISOTH/stages/regen_inputs.py <case> ...

A case is a directory with calliope.in, case.json and generated.md5 (a copy of a case of the per-test folders,
e.g. under runs/). A profile1d.dat or ic/ that exists is not written again, only checked. The fields come from the
generators of the Stage scripts, chosen by case.json:
    'case' and 'dir'               1D MHD of Stage 3, 14, 16: gen_stage3.write_profile
    batch 'planar', dim '1d'       1D reference of a planar problem of Stage 4: as gen_stage4.planar_1d_case
    batch 'planar', dim '2d'       gen_stage4.planar_fields
    batch 'alfven'                 gen_stage4.cp_fields (the batches alfven and cpdiag of Stage 4)
    batch 'ot'                     gen_stage4.ot_fields
    batch 'turb'                   gen_stage4.problem_fields
ic/ is written by gen_ic.write_restart on the grid of calliope.in (analyze_stage4.grid_of), as gen_stage4.ic_case.
A case without generated.md5 (Euler1D: built-in initial conditions) has nothing to write. One line per case; the
exit status is 1 if a file is missing or its md5 differs from generated.md5.

The md5 match bitwise with the numpy of the ps1 runs (2.0.2 on x86_64). Another numpy or CPU changes the last bits
of the FFT and of elementary functions: with numpy 1.26.4 on macOS arm64, 53 of the 266 profile1d.dat and every ic/
differ from generated.md5 by round-off (README.md).
"""
import sys
import os
import json
import math
import hashlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_ic
import gen_stage3 as g3
import gen_stage4 as s4
import stage3 as s3
from analyze_stage4 import grid_of

IC_NAMES = sorted('ic/%s.dat' % n for n in gen_ic.NAMES + ('time',))


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def profile(case, info):
    """profile1d.dat of a 1D case."""
    if 'dir' in info:
        g3.write_profile(case, info['case'], info['dir'])
    elif info.get('batch') == 'planar' and info.get('dim') == '1d':
        # the table of gen_stage4.planar_1d_case
        c, kvec = info['case'], info['kvec']
        N, kn = c['N'], math.hypot(*kvec)
        W = s3.initial(c)
        u, b = g3.rotate(W, s3.bx_of(c), 'x')
        np.savetxt(os.path.join(case, 'profile1d.dat'), np.column_stack([np.arange(N)*(s3.L/kn/N), W[0], *u, *b]),
                   fmt='%.17e', header='%s compressed by 1/|k| = 1/%.6f along x: x, rho, u_x, u_y, u_z, b_x, b_y, b_z'
                                       % (c['name'], kn))
    else:
        raise SystemExit('%s: no generator of profile1d.dat for case.json %s' % (case, sorted(info)))


def fields(case, info, g):
    """Real lnrho, u, B of a 2D or 3D case."""
    batch = info.get('batch')
    if batch == 'ot':
        assert info['B0'] == s4.OT_B0, (case, info['B0'])
        return s4.ot_fields(g)
    if batch == 'turb':
        Ms, dB = s4.TURB_PROBS[info['problem']]
        got = (info['Ms'], info.get('dB', Ms), info['seed'], info['B0'])
        assert got == (Ms, dB, s4.TURB_SEED, 1.0), (case, got)
        return s4.problem_fields(info['problem'], g)[0]
    if batch == 'alfven':
        return s4.cp_fields(g, info['kvec'], info['amp'], tuple(np.array(v) for v in info['basis']))
    if batch == 'planar' and info.get('dim') == '2d':
        c = info['case']
        return s4.planar_fields(s3.initial(c), s3.bx_of(c), c['N'], info['kvec'])
    raise SystemExit('%s: no generator for batch %r' % (case, batch))


def restart(case, info):
    """ic/ of a 2D or 3D case."""
    g = grid_of(case)
    assert list(g.nl) == info['nl'], (case, g.nl, info['nl'])
    gen_ic.write_restart(os.path.join(case, 'ic'), g, *gen_ic.state_from_real(g, *fields(case, info, g)))


def regen(case):
    """'wrote', 'kept' or 'nothing to write', and the names that are missing or differ from generated.md5."""
    if not os.path.exists(os.path.join(case, 'generated.md5')):
        return 'nothing to write', []
    want = {name: h for h, name in (ln.split() for ln in open(os.path.join(case, 'generated.md5')))}
    if sorted(want) == ['profile1d.dat']:
        target, write = 'profile1d.dat', profile
    elif sorted(want) == IC_NAMES:
        target, write = 'ic', restart
    else:
        raise SystemExit('%s: generated.md5 lists %s' % (case, sorted(want)))
    done = os.path.exists(os.path.join(case, target))
    if not done:
        write(case, json.load(open(os.path.join(case, 'case.json'))))
    paths = {n: os.path.join(case, n) for n in sorted(want)}
    return 'kept' if done else 'wrote', [n for n, p in paths.items() if not os.path.exists(p) or md5(p) != want[n]]


def main(cases):
    nbad = 0
    for case in cases:
        what, bad = regen(case)
        print('%s: %s, %s' % (case, what, ('missing or different: %s' % ' '.join(bad)) if bad else 'md5 OK'))
        nbad += bool(bad)
    if nbad:
        sys.exit('%d of %d cases differ from generated.md5' % (nbad, len(cases)))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
