#!/usr/bin/env bash
# Run the nlminb benchmark: controls, then every start of each case under both
# R versions, one fresh R process per fit, one fit at a time. The order of the
# two versions alternates from start to start (as in opal's
# dev/run-benchmark-ab.R) so drift in machine speed affects both equally.
#
# Usage: bash run-benchmark.sh [case ...]   (default: opaka io_yft)
# Env:   R45, R46 (R prefixes; default /opt/r452, /opt/r462),
#        OPAL_SOURCE (opal checkout at the installed commit; for the stored
#        trajectories), NLMINB_N_OPAKA, NLMINB_N_IO_YFT, NLMINB_JF
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
R45=${R45:-/opt/r452}
R46=${R46:-/opt/r462}
export NLMINB_BENCH_DIR=$HERE OPAL_SOURCE=${OPAL_SOURCE:?set OPAL_SOURCE}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LANG=C.UTF-8 LC_ALL=C.UTF-8
CASES=${*:-opaka io_yft}
mkdir -p "$HERE/results/raw"
cd /tmp   # keep the repository's renv .Rprofile out of the way

run() { "$1/bin/Rscript" --vanilla "$HERE/$2" "${@:3}"; }
rver() { "$1/bin/Rscript" --vanilla -e 'cat(gsub("[^0-9]", "", getRversion()))'; }
declare -A V=([$R45]=$(rver "$R45") [$R46]=$(rver "$R46"))

[ -f "$HERE/results/starts.rds" ] || run "$R45" 01-make-starts.R
for R in "$R45" "$R46"; do run "$R" 03-controls.R; done

for case in $CASES; do
  n=$("$R45/bin/Rscript" --vanilla -e \
    "cat(nrow(readRDS('$HERE/results/starts.rds')[['$case']]\$starts))")
  for id in $(seq 0 "$n"); do
    if (( id % 2 == 0 )); then order="$R45 $R46"; else order="$R46 $R45"; fi
    for R in $order; do
      [ -f "$HERE/results/raw/$case-$(printf %02d "$id")-${V[$R]}.rds" ] && continue
      run "$R" 02-run-fit.R "$case" "$id"
    done
  done
done
"$R45/bin/Rscript" --vanilla "$HERE/04-summarise.R"
