#!/bin/bash
# Post-processing of runs/stage1_es on ps1 (no network).  Writes runs/stage1_es/post/*.txt and
# runs/stage1_es/post.tgz (judge table, summaries, these texts) so that only a few MB go back.
#   nohup bash tests/MHD_COMP_ISOTH/post_stage1_es.sh > runs/stage1_es/post.log 2>&1 < /dev/null &
set -u
cd ~/work/calliope_isoth
B=runs/stage1_es
P=$B/post
T=tests/MHD_COMP_ISOTH
mkdir -p $P

# 1. MHD_INCOMP regression against the develop baseline (job 477247)
{
  echo "exit_status: $(cat runs/regress/incomp_branch_es/exit_status)  exe: $(cat runs/regress/incomp_branch_es/exe 2>/dev/null)"
  for f in $(cd runs/regress/incomp_base && find . -type f ! -name exit_status ! -name exe ! -name nproc | sort); do
    # job.pbs-ps-fes-dev runs the input calliope_run.in, so the branch case writes calliope_run.* in place of
    # calliope.* (calliope.in and the redirected calliope.out.std keep their names)
    case $f in
      ./calliope.in|./calliope.out.std) g=$f ;;
      ./calliope.*) g=./calliope_run.${f#./calliope.} ;;
      *) g=$f ;;
    esac
    a=runs/regress/incomp_base/$f; b=runs/regress/incomp_branch_es/$g
    if [ ! -f $b ]; then echo "MISSING $f"; continue; fi
    if cmp -s $a $b; then echo "same  $f"; else echo "DIFF  $f"; diff $a $b 2>/dev/null | head -20; fi
  done
} > $P/regress.txt 2>&1

# 2. fast eigenmode: new runs and the Stage 0 baseline
{
  for d in $B/fast/ideal $B/fast/default $B/fast/default_linear $(find runs/stage0 -name calliope.modes.out -printf '%h\n' | sort); do
    echo "== $d ($(cat $d/exit_status 2>/dev/null))"
    python3 $T/fast_wave.py $d
  done
} > $P/fast.txt 2>&1

# 3. cross-check with the prototype (first 1000 steps)
for c in collide_U2_N512 collide_U2_N512_C1_nosvv_fs; do
  python3 $T/proto1d.py $B/smoke/$c 1000 > $P/proto_$c.txt 2>&1
done

# 4. power method, viscous dt and cost of every case (all chunks)
python3 - $B > $P/svv_stats.txt 2>&1 << 'EOF'
import sys, os, glob, re
import numpy as np
sys.path.insert(0, 'tests/MHD_COMP_ISOTH')
from analyze_stage1 import load
base = sys.argv[1]
cases = sorted(os.path.dirname(p) for p in glob.glob(os.path.join(base, '*', '*', 'calliope.in')))
print('%-44s %-12s %6s %9s %9s %9s %9s %6s %6s %9s %9s %8s %7s' % (
    'case', 'exit', 'chunks', 'steps', 'dt_adv', 'dt_visc', 'max dtlam', 'n_it0', 'n_chk', 'F_min', 'F_max',
    'wall[s]', 'stop'))
for c in cases:
    dirs = sorted(glob.glob(os.path.join(c, 'chunk_*')), key=lambda s: int(s.rsplit('_', 1)[1])) + [c]
    steps, wall, dtl, nit0, nchk, fmin, fmax, stop = 0, 0.0, 0.0, -1, 0, np.inf, 0.0, ''
    dta = dtv = np.nan
    for k, d in enumerate(dirs):
        std = os.path.join(d, 'calliope.out.std')
        if os.path.exists(std):
            s = open(std, errors='replace').read()
            m = re.search(r'# of steps advanced\s+(\d+)', s)
            steps += int(m.group(1)) if m else 0
            m = re.search(r'wall ([0-9.]+) s', s)
            wall += float(m.group(1)) if m else 0.0
            if 'ABORT' in s:
                stop = 'abort'
            elif 'This run is stopped' in s:
                stop = 'cfl'
        f = load(os.path.join(d, 'cfl.dat'))
        if f is None or f.ndim != 2 or f.shape[1] < 16:
            continue
        if k == 0:
            dta, dtv = f[0, 2], f[0, 3]
            nit0 = int(f[0, 15])
        dtl = max(dtl, np.nanmax(f[:, 14]))
        nchk = max(nchk, int(np.nanmax(f[1:, 15])) if len(f) > 1 else 0)
        sv = f[:, 13][f[:, 13] > 0]
        if sv.size:
            fmin, fmax = min(fmin, sv.min()), max(fmax, sv.max())
    ex = open(os.path.join(c, 'exit_status')).read().strip() if os.path.exists(os.path.join(c, 'exit_status')) else '-'
    print('%-44s %-12s %6d %9d %9.3e %9.3e %9.3e %6d %6d %9.3e %9.3e %8.1f %7s' % (
        os.path.relpath(c, base), ex, len(dirs), steps, dta, dtv, dtl, nit0, nchk,
        fmin if np.isfinite(fmin) else 0.0, fmax, wall, stop))
EOF

# 5. the SVV-off path against the linear-SVV build (common columns), when the same case exists there
{
  for c in scan/Sod_N1024_nosvv scan/U2_N1024_nosvv scan/U4_N1024_nosvv; do
    a=runs/stage1_svv/$c; b=$B/$c
    [ -d $a ] || { echo "== $c: no runs/stage1_svv counterpart"; continue; }
    echo "== $c"
    diff <(grep -v -i -E "max_wall_time|svv_form|cfl_svv|svv_pi" $a/calliope.in) \
         <(grep -v -i -E "max_wall_time|svv_form|cfl_svv|svv_pi" $b/calliope.in) | head -10
    python3 - $a $b << 'EOF'
import sys, os
import numpy as np
sys.path.insert(0, 'tests/MHD_COMP_ISOTH')
from analyze_stage1 import load
for name in ('conservation.dat', 'cfl.dat', 'budget.dat'):
    x, y = load(os.path.join(sys.argv[1], name)), load(os.path.join(sys.argv[2], name))
    if x is None or y is None:
        print(name, 'missing'); continue
    n, m = min(len(x), len(y)), min(x.shape[1], y.shape[1])
    d = np.abs(x[:n, :m] - y[:n, :m])
    print('%-17s rows %d/%d cols %d: max|diff| = %.3e, identical = %s' % (name, len(x), len(y), m, np.nanmax(d),
          bool(np.array_equal(x[:n, :m], y[:n, :m], equal_nan=True))))
EOF
  done
} > $P/nosvv.txt 2>&1

ls -la bin > $P/bin.txt
tar czf $B/post.tgz $B/judge.md $B/judge.json $B/judge.err $P $(find $B -name summary.json) \
  $(find $B -name exit_status) 2>/dev/null
echo "=== post done $(date)"
