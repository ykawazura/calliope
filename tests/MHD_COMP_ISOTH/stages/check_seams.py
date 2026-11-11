"""Energy budget across the restart seams of a chained case (job.pbs-ps-fes-dev).

  python3 tests/MHD_COMP_ISOTH/stages/check_seams.py <case> ... [--ref <case>]

For each case, the budget rows of its runs (chunk_1, chunk_2, ..., case directory) are checked.

- Steps are counted by the rows of cfl.dat (one per step).  A binary before the fix of deallocate_advance
  (bin/calliope_isoth_s4 and older) leaves out the budget row of the last step of a run that check_terminate
  stopped (max_wall_time or the file "terminate"), since that row is written at the next step.
- At each seam the gap is the time between the end of the last row of a run and the start of the first
  row of the next one, and the jump is E_n of that first row minus E_np1 of that last row.
- residual = dE - sum dt (dE/dt_sum + err) - sum of the jumps, which is 0 to round-off when the rows of
  each run close.

With --ref, the rows of each case are compared with those of the reference case at the same t.
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_stage1 import load_runs

# budget.dat columns (0-based)
T, DT, RATE, ERR, EN, ENP1 = 0, 1, 7, 8, 18, 19


def check(case):
    bs = load_runs(case, 'budget.dat')
    nstep = sum(len(a) for a in load_runs(case, 'cfl.dat'))
    b = np.concatenate(bs)
    seams = []
    for p, n in zip(bs, bs[1:]):
        gap = n[0, T] - (p[-1, T] + p[-1, DT])
        seams.append((n[0, T], gap/p[-1, DT], n[0, EN] - p[-1, ENP1]))
    dE = b[-1, ENP1] - b[0, EN]
    integral = np.sum(b[:, DT]*(b[:, RATE] + b[:, ERR]))
    jumps = sum(s[2] for s in seams)
    res = dE - integral - jumps
    ok = len(b) == nstep and all(abs(s[1]) < 1e-6 for s in seams) and abs(res) <= 1e-12*np.max(np.abs(b[:, EN]))
    print('%s: %d runs, %d steps, %d budget rows (%d missing)' % (case, len(bs), nstep, len(b), nstep - len(b)))
    for t, g, j in seams:
        print('  seam t = %.6f: gap %.3f dt, E jump %.4e' % (t, g, j))
    print('  dE %.10e  sum dt (rate + err) %.10e  sum jumps %.4e  residual %.2e  %s'
          % (dE, integral, jumps, res, 'PASS' if ok else 'FAIL'))
    return b, ok


def compare(b, r):
    """max |b - r| over the rows of b at a t of r, by column, relative to max |r| of the column."""
    i = np.searchsorted(r[:, T], b[:, T])
    i = np.clip(i, 0, len(r) - 1)
    on = np.abs(r[i, T] - b[:, T]) <= 1e-12*max(1.0, np.max(np.abs(r[:, T])))
    d = np.abs(b[on] - r[i[on]])/np.maximum(np.max(np.abs(r), axis=0), 1e-300)
    return int(np.sum(on)), np.max(d, axis=0)


if __name__ == '__main__':
    args = sys.argv[1:]
    ref = None
    if '--ref' in args:
        k = args.index('--ref')
        ref = args[k + 1]
        args = args[:k] + args[k + 2:]
    allok = True
    rb = check(ref)[0] if ref else None
    for c in args:
        b, ok = check(c)
        allok &= ok
        if rb is not None:
            n, d = compare(b, rb)
            print('  vs %s: %d of %d rows at its t, max relative difference %.2e (E_np1 %.2e)'
                  % (ref, n, len(b), np.max(d), d[ENP1]))
    print('all PASS' if allok else 'some FAIL')
