#!/usr/bin/env bash
# R 4.5 + opal (dev) + r4ss + Quarto for the SS3 -> opal IO yellowfin work.
# For containers where CRAN and GitHub HTTP are blocked but conda-forge and
# `git clone` from GitHub work. Installs into /opt/renv45; ~20-30 min.
set -euo pipefail

PREFIX=/opt/renv45
WORK=${WORK:-$HOME/r-setup}
mkdir -p "$WORK/src" && cd "$WORK"

# 1. micromamba (static binary from conda-forge)
curl -fsSL -o mm.tar.bz2 https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.3.3-0.tar.bz2
mkdir -p mm && tar -xjf mm.tar.bz2 -C mm
export MAMBA_ROOT_PREFIX="$WORK/mamba"

# 2. R and every dependency available on conda-forge, plus compilers and Quarto
mm/bin/micromamba create -y -p "$PREFIX" -c conda-forge \
  "r-base=4.5" r-data.table r-yaml r-tmb r-dplyr r-forecast r-ggplot2 r-ggridges \
  r-reshape2 r-rlang r-scales r-cli r-coda r-corpcor r-forcats r-furrr r-future r-gh \
  r-glue r-kableextra r-lifecycle r-purrr r-stringr r-tidyr r-viridis r-ellipse \
  r-microbenchmark r-mvtnorm r-matrix r-posterior r-sparseinv r-rcpp r-rcppparallel \
  r-rcppeigen r-bh r-loo r-this.path r-testthat r-remotes \
  r-circular r-sn r-statmod r-rapidjsonr r-stanheaders r-rstan r-slam r-clue r-skmeans \
  r-rmarkdown r-knitr r-patchwork quarto \
  gxx_linux-64 gcc_linux-64 gfortran_linux-64 make

export PATH="$PREFIX/bin:$PATH"

# 3. Source packages (CRAN-only ones via the github.com/cran mirrors)
cd "$WORK/src"
git clone -q --depth 1 https://github.com/cran/RcppParallel
git clone -q --depth 1 https://github.com/cran/movMF
git clone -q --depth 1 https://github.com/kaskr/RTMB
git clone -q --depth 1 https://github.com/janoleko/RTMBdist
git clone -q --depth 1 https://github.com/andrjohns/StanEstimators
git clone -q --depth 1 https://github.com/noaa-afsc/SparseNUTS
git clone -q --depth 1 https://github.com/r4ss/r4ss
git clone -q -b dev https://github.com/N-DucharmeBarth-NOAA/opal
# (optional) pin opal to the commit the vignettes were built with:
# git -C opal checkout c493461

# RcppParallel must be rebuilt from source: conda-forge's build lacks the
# TBB_INTERFACE_NEW flag that StanEstimators' Stan headers need.
# Its bundled-TBB build fails under a parallel make, so no MAKEFLAGS here.
env -u MAKEFLAGS R CMD INSTALL RcppParallel

# Order matters: each line depends on the ones above it.
R CMD INSTALL RTMB/RTMB
R CMD INSTALL movMF
R CMD INSTALL RTMBdist
MAKEFLAGS=-j4 R CMD INSTALL StanEstimators
R CMD INSTALL SparseNUTS
R CMD INSTALL opal
R CMD INSTALL r4ss

# 4. Check
cd /tmp   # avoid any project .Rprofile (renv) in the current directory
Rscript -e 'suppressMessages({library(opal); library(r4ss); library(yaml); library(data.table)});
  cat("R", as.character(getRversion()), "| opal", as.character(packageVersion("opal")),
      "| r4ss", as.character(packageVersion("r4ss")), "| opal_obj:", exists("opal_obj"),
      "| age-specific K:", "log_k_mult" %in% names(formals(opal::get_growth)), "\n")'
QUARTO_SHARE_PATH="$PREFIX/share/quarto" quarto --version

cat <<MSG

Done. In each new shell:
  export PATH=$PREFIX/bin:\$PATH LANG=C.UTF-8 LC_ALL=C.UTF-8
  export QUARTO_SHARE_PATH=$PREFIX/share/quarto QUARTO_DENO=$PREFIX/bin/deno \\
         QUARTO_ESBUILD=$PREFIX/bin/esbuild QUARTO_DART_SASS=$PREFIX/bin/sass \\
         QUARTO_PANDOC=$PREFIX/bin/pandoc QUARTO_TYPST=$PREFIX/bin/typst
Run R scripts from outside the repo root (or with --no-init-file), because the
repo's .Rprofile activates renv, which tries to reach CRAN.
MSG
