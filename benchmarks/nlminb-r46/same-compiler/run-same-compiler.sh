#!/usr/bin/env bash
# Same-compiler test: build the PORT variants, then run the plain-R check and
# every start of both cases with each variant in R 4.6.2, one R process per
# case and variant, one at a time.
# Usage: bash run-same-compiler.sh      Env: R46 (default /opt/r462)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
R46=${R46:-/opt/r462}
export NLMINB_BENCH_DIR=$(dirname "$HERE")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LANG=C.UTF-8 LC_ALL=C.UTF-8
[ -f "$HERE/lib/port-E.so" ] || bash "$HERE/build-port-variants.sh"
cd /tmp   # keep the repository's renv .Rprofile out of the way
"$R46/bin/Rscript" --vanilla "$HERE/rosen-variants.R"
for case in opaka io_yft; do
  for v in A B C D E; do
    [ -f "$HERE/results/raw/$case-$v.rds" ] && continue
    "$R46/bin/Rscript" --vanilla "$HERE/run-variants.R" "$case" "$v"
  done
done
"$R46/bin/Rscript" --vanilla "$HERE/summarise-variants.R"
