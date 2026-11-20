#!/usr/bin/env bash
# Build nlminb's PORT code (R's src/library/stats/src/port.c + portsrc.f) as
# standalone shared libraries, to separate R 4.6.2's nlminb fix (PR#19086)
# from the compiler change between the conda-forge R builds.
#
#   A  portsrc.f from R 4.5.2 (old), GCC 10.4, -O2   (as R 4.5.2's stats.so)
#   B  old,                         GCC 15.3, -O2
#   C  portsrc.f from R 4.6.2 (fix), GCC 10.4, -O2
#   D  fix,                          GCC 15.3, -O2   (as R 4.6.2's stats.so)
#   E  old,                         GCC 10.4, Fortran at -O0
#
# Only portsrc.f (where the fix is) varies. port.c, the C glue, is R 4.6.2's
# in all variants (it differs from 4.5.2's only in using R_getVar() for a
# variable lookup) and is compiled with GCC 15.3 at -O2 in all of them: R 4.6's
# headers need a C compiler with C23 enum base types, which GCC 10 lacks.
# Flags are R's SAFE_FFLAGS / CFLAGS for portsrc.f / port.c in the
# conda-forge builds, without the path options.
#
# Usage: bash build-port-variants.sh
# Env:   GCC10, GCC15 (toolchain prefixes; default /opt/gcc10, /opt/gcc15),
#        R46 (R used to load the libraries; default /opt/r462)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
GCC10=${GCC10:-/opt/gcc10}
GCC15=${GCC15:-/opt/gcc15}
R46=${R46:-/opt/r462}
SRC=$HOME/rsrc

# The two conda-forge toolchains that built R 4.5.2's and R 4.6.2's stats.so
# (micromamba from setup-r452.sh / setup-r462.sh)
MM=${MM:-$HOME/r-setup-r462/mm/bin/micromamba}
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/r-setup-r462/mamba}
[ -x "$GCC10/bin/x86_64-conda-linux-gnu-gfortran" ] || "$MM" create -y -q -p "$GCC10" \
  -c conda-forge "gcc_linux-64=10.4" "gfortran_linux-64=10.4"
[ -x "$GCC15/bin/x86_64-conda-linux-gnu-gfortran" ] || "$MM" create -y -q -p "$GCC15" \
  -c conda-forge "gcc_linux-64=15.3" "gfortran_linux-64=15.3"
OUT=$HERE/lib
mkdir -p "$SRC" "$OUT"

# R sources (sparse checkouts of the two tags on the GitHub mirror)
for tag in 4-5-2 4-6-2; do
  d=$SRC/r${tag//-/}
  if [ ! -d "$d" ]; then
    git clone -q --depth 1 --filter=blob:none --sparse -b "tags/R-$tag" \
      https://github.com/wch/r-source "$d"
    git -C "$d" sparse-checkout set src/library/stats/src src/library/stats/R
  fi
done
S45=$SRC/r452/src/library/stats/src
S46=$SRC/r462/src/library/stats/src
cp "$SRC/r462/src/library/stats/R/nlminb.R" "$OUT/nlminb.R"
diff "$S45/portsrc.f" "$S46/portsrc.f" > "$OUT/portsrc.diff" || true

FLAGS="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -ffunction-sections -pipe"
SAFE="-msse2 -mfpmath=sse"
RINC=$("$R46/bin/R" RHOME)/include
RLIB=$("$R46/bin/R" RHOME)/lib

build() {  # name gcc-prefix portsrc-dir fortran-opt
  local name=$1 gcc=$2 fsrc=$3 fopt=$4 w=$OUT/build-$1
  mkdir -p "$w"
  cp "$S46/port.c" "$S46/port.h" "$S46/statsErr.h" "$w/"
  cp "$fsrc/portsrc.f" "$w/"
  (cd "$w" &&
   "$GCC15/bin/x86_64-conda-linux-gnu-gcc" -I"$RINC" -DNDEBUG $FLAGS -O2 -c port.c -o port.o &&
   "$gcc/bin/x86_64-conda-linux-gnu-gfortran" $FLAGS $fopt $SAFE -c portsrc.f -o portsrc.o &&
   "$GCC15/bin/x86_64-conda-linux-gnu-gcc" -shared -o "$OUT/port-$name.so" port.o portsrc.o \
     -L"$RLIB" -L"$R46/lib" -lR $("$R46/bin/R" CMD config BLAS_LIBS) -lgfortran -lm \
     -Wl,-rpath,"$RLIB" -Wl,-rpath,"$R46/lib")
  "$gcc/bin/x86_64-conda-linux-gnu-gfortran" --version | head -1 > "$OUT/port-$name.compiler"
  echo "$name: portsrc.f from $(basename "$(dirname "$(dirname "$(dirname "$(dirname "$fsrc")")")")"), $(cat "$OUT/port-$name.compiler"), $fopt"
}
build A "$GCC10" "$S45" -O2
build B "$GCC15" "$S45" -O2
build C "$GCC10" "$S46" -O2
build D "$GCC15" "$S46" -O2
build E "$GCC10" "$S45" -O0
