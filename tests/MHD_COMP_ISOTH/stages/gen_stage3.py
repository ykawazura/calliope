#!/usr/bin/env python3
"""Stage 3 (1D isothermal MHD) case generator for the Fortran code.

    python3 gen_stage3.py nstep                 (here) the fixed steps of the prototype runs -> stage3_nstep.json
                                                (the cases of stage3.cases(ETAS))
    python3 gen_stage3.py <batch> [base]        writes <base>/<batch>/<case>/{calliope.in, nproc, profile1d.dat,
                                                case.json} (default base: runs/stage3_f) and prints the CASES string
                                                for job.pbs-ps-fes-dev, by decreasing cost

Every case is a case of stage3.py: the initial data stage3.initial(c) go to profile1d.dat (init_type = 'profile1d',
one row x, rho, u_x, u_y, u_z, B_x, B_y, B_z per grid point, written with 17 digits), and the run uses the same
fixed dt = t_end/nstep as the prototype run that reached t_end (stage3_nstep.json, written from runs/stage3/*.json),
with reset_method = 'stop', so that a run that violates min(dt_adv, dt_visc) stops. The profiles along y or z (decomp)
are the x profile with the components rotated cyclically (normal, t1, t2) -> (y, z, x) or (z, x, y). c = 1
(beta0 = 2), L = 1 (lx = 1/(2 pi)), out1d at the steps of the NWRITE + 1 prototype snapshots.

Batches
    smoke   the prototype against the Fortran code: the eight cases at N = 256 of the first prototype run
    agree   the same at larger N: F1 (N = 1024), RDh_w4, TDA_w1pre and C3 (N = 512)
    decomp  F1, RDh_w4, C3 at N = 256 along x with 1, 2, 4 ranks, along y with 2 and along z with 4 ranks
    suite   every case of stage3.cases() (bulk viscosity only)
    eta_smoke, eta_decomp, eta_suite   the same with the shock resistivity C_eta = ETA (stage3 name suffix _eta<C>):
            smoke on I1, F1, C3, RDpc at N = 256, decomp on I1 and C3, suite on the cases of the eta scan
    conv    I1 and RDpc at N = 4096, 8192 with C_eta = 0 and ETA (no prototype runs: dt from the largest-N one)
"""
import sys
import os
import json
import math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths
from isoth_case import write_case
from gen_stage1 import box_1d
import stage3 as s3

HERE = os.path.dirname(os.path.abspath(__file__))
NSTEP_FILE = os.path.join(HERE, 'stage3_nstep.json')
# ranks of the longest runs (1 elsewhere)
NPROC = {2048: 2, 4096: 4, 8192: 8}
# C_eta of the prototype eta scan (stage3.py run -eta) and the one of the eta batches
ETAS = (0.0, 0.5, 2.0, 8.0)
ETA = 2.0
# problems and resolutions of the batch conv
CONV = ('I1', 'RDpc')
CONV_N = (4096, 8192)


def write_nstep():
    """{case: [nstep, dt]} of the attempt that reached t_end in runs/stage3/<case>.json."""
    d = {}
    for c in s3.cases(ETAS):
        fn = os.path.join(s3.OUT, c['name'] + '.json')
        if not os.path.exists(fn):
            continue
        r = json.load(open(fn))
        ok = [a for a in r['attempts'] if a['status'] == 'ok']
        if ok:
            d[c['name']] = [ok[0]['nstep'], ok[0]['dt']]
    with open(NSTEP_FILE, 'w') as f:
        json.dump(d, f, indent=0, sort_keys=True)
    print('%d cases -> %s' % (len(d), NSTEP_FILE))


def rotate(W, bx, dir):
    """(u, B) in x, y, z of the 1D profile W = (rho, u_n, u_t1, u_t2, B_t1, B_t2) along dir."""
    u = [W[1], W[2], W[3]]
    b = [bx + 0.0 * W[0], W[4], W[5]]
    i = 'xyz'.index(dir)
    perm = [(j - i) % 3 for j in range(3)]       # component j of x, y, z is the (j - i)th of (n, t1, t2)
    return [u[p] for p in perm], [b[p] for p in perm]


def write_profile(path, c, dir='x'):
    W = s3.initial(c)
    x = s3.grid(c['N'])
    u, b = rotate(W, s3.bx_of(c), dir)
    tab = np.column_stack([x, W[0], *u, *b])
    np.savetxt(os.path.join(path, 'profile1d.dat'), tab, fmt='%.17e',
               header='%s along %s: x, rho, u_x, u_y, u_z, b_x, b_y, b_z' % (c['name'], dir))


def case_dir(root, c, dir='x', nproc=None, tag='', nstep=None):
    """Write the Fortran case of the stage3 case c (nstep: in place of the one of the prototype run)."""
    te = s3.t_end(c)
    if nstep is None:
        nstep, dt = json.load(open(NSTEP_FILE))[c['name']]
    else:
        dt = te / nstep
    assert abs(nstep * dt - te) < 1e-12 * te, (c['name'], nstep, dt, te)
    nproc = NPROC.get(c['N'], 1) if nproc is None else nproc
    path = os.path.join(root, c['name'] + tag)
    write_case(path, nproc=nproc, header='Stage 3 case %s along %s (gen_stage3.py)' % (c['name'], dir),
               box=dict(box_1d(c['N'], dir, nproc), lx=s3.LX, ly=s3.LX, lz=s3.LX),
               time=dict(dt=dt, nstep=nstep),
               shock=dict(C_eta=c['C_eta']) if c['C_eta'] > 0.0 else {},
               # write_line writes the state at tt when (tt, tt + dt] contains a multiple of the interval: with
               # (t_end/NWRITE)(1 + 1/(2 nstep)) the j-th multiple is j dt/(2 NWRITE) after step j nstep/NWRITE,
               # so the lines are at the steps of the prototype snapshots (the last one at t_end, istage = 0)
               isoth_diag=dict(write_intvl_1D=te / s3.NWRITE * (1.0 + 0.5 / nstep)),
               ic=('profile1d', dict(dir="'%s'" % dir, table="'profile1d.dat'")))
    write_profile(path, c, dir)
    with open(os.path.join(path, 'case.json'), 'w') as f:
        json.dump(dict(case=c, dir=dir, nstep=nstep, dt=dt, t_end=te), f, indent=1)
    return path


def find(name):
    return next(c for c in s3.cases(ETAS) if c['name'] == name)


SMOKE = ('F1_N256', 'I1_N256', 'P1b_N256', 'C3_N256', 'RDh_w4_N256', 'RDpc_N256', 'TDS_w1_N256', 'TDA_w1pre_N256')


def batch_smoke(root):
    return [case_dir(root, find(n), nproc=1) for n in SMOKE]


def batch_agree(root):
    return [case_dir(root, find(n)) for n in ('F1_N1024', 'RDh_w4_N512', 'TDA_w1pre_N512', 'C3_N512')]


def batch_decomp(root):
    out = []
    for n in ('F1_N256', 'RDh_w4_N256', 'C3_N256'):
        for dir, np_ in (('x', 1), ('x', 2), ('x', 4), ('y', 2), ('z', 4)):
            out.append(case_dir(root, find(n), dir=dir, nproc=np_, tag='_%s_np%d' % (dir, np_)))
    return out


def batch_suite(root):
    nst = json.load(open(NSTEP_FILE))
    return [case_dir(root, c) for c in s3.cases() if c['name'] in nst]


def eta_name(n):
    return n + '_eta%g' % ETA


def batch_eta_smoke(root):
    return [case_dir(root, find(eta_name(n)), nproc=1) for n in ('I1_N256', 'F1_N256', 'C3_N256', 'RDpc_N256')]


def batch_eta_decomp(root):
    out = []
    for n in ('I1_N256', 'C3_N256'):
        for dir, np_ in (('x', 1), ('x', 2), ('x', 4), ('y', 2), ('z', 4)):
            out.append(case_dir(root, find(eta_name(n)), dir=dir, nproc=np_, tag='_%s_np%d' % (dir, np_)))
    return out


def batch_eta_suite(root):
    nst = json.load(open(NSTEP_FILE))
    return [case_dir(root, c) for c in s3.cases((ETA,)) if c['name'] in nst]


def batch_conv(root):
    """I1 and RDpc beyond the prototype resolutions, with C_eta = 0 and ETA: nstep of the largest-N prototype run
    of the case times N/N_ref (dt in proportion to the grid spacing)."""
    nst = json.load(open(NSTEP_FILE))
    out = []
    for p in CONV:
        for eta in (0.0, ETA):
            ref = max((c for c in s3.cases(ETAS) if c['name'].startswith(p + '_N') and c['C_eta'] == eta
                       and c['name'] in nst), key=lambda c: c['N'])
            for N in CONV_N:
                c = dict(ref, N=N, name=ref['name'].replace('_N%d' % ref['N'], '_N%d' % N))
                out.append(case_dir(root, c, nstep=nst[ref['name']][0] * N // ref['N']))
    return out


BATCHES = {'smoke': batch_smoke, 'agree': batch_agree, 'decomp': batch_decomp, 'suite': batch_suite,
           'eta_smoke': batch_eta_smoke, 'eta_decomp': batch_eta_decomp, 'eta_suite': batch_eta_suite,
           'conv': batch_conv}


def cost(case):
    """nstep * (grid points) * log2(grid points), for listing the cases by decreasing cost."""
    from analyze_stage1 import read_nml
    d = read_nml(os.path.join(case, 'calliope.in'))
    n = int(d['nlx']) * int(d['nly']) * int(d['nlz'])
    return int(d['nstep']) * n * math.log2(n)


if __name__ == '__main__':
    batch = sys.argv[1]
    if batch == 'nstep':
        write_nstep()
        sys.exit()
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.join(paths.RUNS, 'stage3_f')
    print(' '.join(sorted(BATCHES[batch](os.path.join(base, batch)), key=cost, reverse=True)))
