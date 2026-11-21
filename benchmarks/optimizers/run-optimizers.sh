#!/usr/bin/env bash
# Optimizer benchmark: every method on every start of both cases, in R 4.6.2,
# one R process per case and method, one at a time. Uses the models and
# jittered starts of benchmarks/nlminb-r46. nlminb_hessian is run for opaka
# only (one IO yellowfin Hessian takes about 18 s).
# Usage: bash run-optimizers.sh      Env: R46 (default /opt/r462; needs nloptr)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
R46=${R46:-/opt/r462}
export OPT_BENCH_DIR=$HERE NLMINB_BENCH_DIR=$(dirname "$HERE")/nlminb-r46
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LANG=C.UTF-8 LC_ALL=C.UTF-8
METHODS="nlminb2 nlminb1 nlminb1_newton nlminb_hessian nlminb2_reltol8 lbfgsb lbfgsb_guard nloptr_lbfgs nloptr_tnewton nlminb2_atomic"
cd /tmp   # keep the repository's renv .Rprofile out of the way
for case in opaka io_yft; do
  for m in $METHODS; do
    [ "$case" = io_yft ] && [ "$m" = nlminb_hessian ] && continue
    [ -f "$HERE/results/raw/$case-$m.rds" ] && continue
    "$R46/bin/Rscript" --vanilla "$HERE/run-method.R" "$case" "$m"
  done
done
"$R46/bin/Rscript" --vanilla "$HERE/profile-gradient.R"
"$R46/bin/Rscript" --vanilla "$HERE/summarise.R"
