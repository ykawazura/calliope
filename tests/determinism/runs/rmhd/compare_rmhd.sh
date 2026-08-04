#!/bin/sh
#=====================================================================
#  Run compare.py over every pair produced by job.pbs-rmhd-gallope and
#  job.pbs-rmhd-calliope, and collect the output in one file.
#
#  Usage:  ./compare_rmhd.sh [OUT]      (run from the directory holding
#                                        run_*/ and probe_*/)
#
#  Three groups of comparisons:
#
#    step count   gallope 4 GPU vs calliope, at nstep = 1 .. 1000, and
#                 gallope 1 GPU vs gallope 4 GPU as the summation-order
#                 control.  The control must be far below the code-to-code
#                 difference, otherwise that difference is a decomposition
#                 artefact rather than a statement about the two codes.
#
#    self-conv.   each code against itself at dt and dt/2, at fixed final
#                 time t = 0.1.  This is the sharp diagnostic: it is a
#                 property of one code alone and so catches an error the
#                 two codes share, which code-to-code agreement cannot.
#
#    code-to-code at each dt of the probe, to show that the difference is
#                 flat in dt (roundoff) rather than proportional to it.
#=====================================================================
set -u
OUT=${1:-compare_rmhd.out}
HERE=$(cd "$(dirname "$0")" && pwd)
CMP="${HERE}/../compare.py"
NP=16
GRID="--model RMHD --nlx 128 --nly 128 --nlz 128 --lx 1.0 --ly 1.0 --lz 1.0"

: > "${OUT}"
for n in 1 10 100 1000; do
  g4="run_n${n}_g4"
  g1="run_n${n}_g1"
  cc="run_n${n}_c${NP}"

  if [ -d "${cc}/restart" ] && [ -d "${g4}/restart" ]; then
    { echo "############ nstep=${n}  gallope(4 GPU) vs calliope"
      python3 "${CMP}" --a "${g4}" --a-code gallope \
                       --b "${cc}" --b-code calliope ${GRID}
      echo; } >> "${OUT}" 2>&1
  else
    echo "############ nstep=${n}  SKIPPED (missing ${cc} or ${g4})" >> "${OUT}"
  fi

  if [ -d "${g1}/restart" ] && [ -d "${g4}/restart" ]; then
    { echo "############ nstep=${n}  gallope 1 GPU vs 4 GPU [control]"
      python3 "${CMP}" --a "${g1}" --a-code gallope \
                       --b "${g4}" --b-code gallope ${GRID}
      echo; } >> "${OUT}" 2>&1
  else
    echo "############ nstep=${n}  control SKIPPED" >> "${OUT}"
  fi
done

for code in gallope calliope; do
  for pair in "dt1e3 dt5e4" "dt5e4 dt25e5"; do
    set -- ${pair}; a=$1; b=$2
    { echo "############ ${code} self-convergence, ${a} vs ${b}"
      python3 "${CMP}" --a "probe_${a}_${code}" --a-code ${code} \
                       --b "probe_${b}_${code}" --b-code ${code} ${GRID}
      echo; } >> "${OUT}" 2>&1
  done
done

for tag in dt1e3 dt5e4 dt25e5; do
  { echo "############ gallope vs calliope at ${tag}"
    python3 "${CMP}" --a "probe_${tag}_gallope" --a-code gallope \
                     --b "probe_${tag}_calliope" --b-code calliope ${GRID}
    echo; } >> "${OUT}" 2>&1
done

echo "wrote ${OUT}"
grep -c '^############' "${OUT}"
grep -n 'SKIPPED\|MISMATCH\|Traceback\|SystemExit' "${OUT}" || echo "no skips or gate failures"
