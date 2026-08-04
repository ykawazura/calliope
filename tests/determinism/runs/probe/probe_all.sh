#!/bin/sh
#=====================================================================
#  Collect every comparison the three probe jobs make possible:
#
#    probe_*   shear on,  nonlinear on    (job.pbs-shearprobe)
#    order_*   shear off, nonlinear on    (job.pbs-orderprobe)
#    lin_*     shear on,  nonlinear off   (job.pbs-linorder)
#
#  For each series and each code, the two successive dt halvings give the
#  self-convergence order, which is the sharp diagnostic: it is a property
#  of one code alone and does not depend on the two codes agreeing.  The
#  code-to-code difference at each dt is reported alongside it, so that the
#  two can be compared directly.
#
#  Usage:  ./probe_all.sh [OUT]
#=====================================================================
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:-probe_all.out}
GRID="--nlx 128 --nly 128 --nlz 128 --lx 1.0 --ly 1.0 --lz 1.0"

: > "${OUT}"
for spec in "probe 1.5" "order 0.0" "lin 1.5"; do
  set -- ${spec}; pre=$1; q=$2
  for code in gallope calliope; do
    for pair in "dt1e3 dt5e4" "dt5e4 dt25e5"; do
      set -- ${pair}; a=$1; b=$2
      { echo "############ ${pre}: ${code} self-convergence, ${a} vs ${b}"
        python3 "${HERE}/compare.py" --a "${pre}_${a}_${code}" --a-code ${code} \
                                     --b "${pre}_${b}_${code}" --b-code ${code} \
                                     ${GRID} --q ${q}
        echo; } >> "${OUT}" 2>&1
    done
  done
  for tag in dt1e3 dt5e4 dt25e5; do
    { echo "############ ${pre}: gallope vs calliope at ${tag}"
      python3 "${HERE}/compare.py" --a "${pre}_${tag}_gallope" --a-code gallope \
                                   --b "${pre}_${tag}_calliope" --b-code calliope \
                                   ${GRID} --q ${q}
      echo; } >> "${OUT}" 2>&1
  done
done

echo "wrote ${OUT}"
grep -c '^############' "${OUT}"
grep -n 'MISMATCH\|Traceback\|SystemExit' "${OUT}" || echo "no gate failures"
