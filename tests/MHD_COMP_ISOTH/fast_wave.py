#!/usr/bin/env python3
"""Frequency and damping of the fast eigenmode (kx = kz = 1) from calliope.modes.out.

    python3 fast_wave.py <run_dir> [...]

Fits the unwrapped phase of the mx coefficient at k = (1, 0, 1) linearly in t and prints
omega_fit, the theoretical omega = sqrt(2 + sqrt 2) (cs = vA = 1, k^2 = 2, k at 45 deg to B0) and
|mx(t_end)|/|mx(t_0)|. Used for the Stage 0 baseline and the new model.
"""
import sys
import os
import math
import numpy as np

OMEGA = math.sqrt(2.0 + math.sqrt(2.0))


def fit(run, field='mx', k=(1.0, 0.0, 1.0)):
    t, a = [], []
    with open(os.path.join(run, 'calliope.modes.out')) as f:
        for line in f:
            w = line.split()
            if w[1] != field or tuple(float(x) for x in w[2:5]) != k:
                continue
            t.append(float(w[0]))
            a.append(complex(float(w[5]), float(w[6])))
    t, a = np.array(t), np.array(a)
    ph = np.unwrap(np.angle(a))
    slope = np.polyfit(t, ph, 1)[0]
    return dict(n=len(t), t0=t[0], t1=t[-1], omega_fit=abs(slope), rel=abs(slope) / OMEGA - 1.0,
                amp_ratio=abs(a[-1]) / abs(a[0]))


if __name__ == '__main__':
    for run in sys.argv[1:]:
        r = fit(run)
        print('%s: n=%d t=[%.3f, %.3f] omega_fit=%.8f (theory %.8f, rel %+.2e) |a(t1)|/|a(t0)|=%.6f'
              % (run, r['n'], r['t0'], r['t1'], r['omega_fit'], OMEGA, r['rel'], r['amp_ratio']))
