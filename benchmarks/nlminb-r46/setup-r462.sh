#!/usr/bin/env bash
# R 4.6.2 + opal for the nlminb benchmark, in a container where CRAN is
# blocked and conda-forge has R 4.6.2 but (as of 2026-11-20) no R packages
# built for it. R and the compilers come from conda-forge; every R package is
# built from source at the version installed in the R 4.5.2 library, using the
# version tags of the github.com/cran mirrors. GitHub-only packages (RTMB,
# RTMBdist, StanEstimators, SparseNUTS) and opal reuse the clones made by the
# R 4.5.2 install, so both libraries hold the same package versions.
#
# Usage: bash setup-r462.sh <R 4.5.2 prefix> <R 4.5.2 source dir>
#   e.g. bash setup-r462.sh /opt/r452 $HOME/r-setup-r452/src
set -euo pipefail
R45=$1
SRC45=$2
PREFIX=/opt/r462
WORK=$HOME/r-setup-r462
mkdir -p "$WORK/cran" "$WORK/repo/src/contrib"
cd "$WORK"
export MAMBA_ROOT_PREFIX=$WORK/mamba
[ -x "$PREFIX/bin/R" ] || mm/bin/micromamba create -y -p "$PREFIX" -c conda-forge \
  r-base=4.6.2 gxx_linux-64 gcc_linux-64 gfortran_linux-64 make cmake pkg-config \
  zlib libcurl openssl libxml2
export PATH=$PREFIX/bin:$PATH

# CRAN packages at the versions in the R 4.5.2 library (dependency closure of
# opal's Imports, RTMB, RTMBdist, StanEstimators, SparseNUTS; base packages
# excluded, recommended packages such as Matrix included).
"$R45/bin/Rscript" --vanilla -e '
ip <- installed.packages()
base <- rownames(ip)[ip[, "Priority"] %in% "base"]
roots <- c("opal", "RTMB", "RTMBdist", "StanEstimators", "SparseNUTS", "TMB", "Matrix", "RcppEigen")
d <- unlist(tools::package_dependencies(roots, db = ip, which = c("Depends", "Imports", "LinkingTo"),
                                        recursive = TRUE))
d <- setdiff(unique(d), c(base, roots[1:5], "movMF"))
writeLines(paste(d, ip[d, "Version"]), "deps.txt")'
# RcppParallel's bundled TBB fails under a parallel make: built on its own
# below, from the R 4.5.2 install's clone (as in shared/setup/install-r-opal.sh).
grep -v '^RcppParallel ' deps.txt > deps-repo.txt
while read -r pkg ver; do
  [ -d "cran/$pkg" ] || git clone -q --depth 1 -b "$ver" "https://github.com/cran/$pkg" "cran/$pkg"
  tar -czf "repo/src/contrib/${pkg}_${ver}.tar.gz" -C cran "$pkg"
done < deps-repo.txt
Rscript --vanilla -e 'tools::write_PACKAGES("repo/src/contrib", type = "source")'
Rscript --vanilla -e 'install.packages(file.path("repo/src/contrib",
  grep("^Rcpp_", dir("repo/src/contrib"), value = TRUE)), repos = NULL, type = "source")'
env -u MAKEFLAGS R CMD INSTALL "$SRC45/RcppParallel"
MAKEFLAGS=-j4 Rscript --vanilla -e '
pk <- read.dcf("repo/src/contrib/PACKAGES")[, "Package"]
install.packages(pk, repos = paste0("file://", normalizePath("repo")), type = "source",
                 Ncpus = 2)'
R CMD INSTALL "$SRC45/RTMB/RTMB"
R CMD INSTALL "$SRC45/movMF"
R CMD INSTALL "$SRC45/RTMBdist"
MAKEFLAGS=-j4 R CMD INSTALL "$SRC45/StanEstimators"
R CMD INSTALL "$SRC45/SparseNUTS"
R CMD INSTALL "$SRC45/opal"

cd /tmp
Rscript --vanilla -e 'suppressMessages(library(opal)); cat("R", as.character(getRversion()),
  "| opal", as.character(packageVersion("opal")), "| RTMB", as.character(packageVersion("RTMB")),
  "| TMB", as.character(packageVersion("TMB")), "\n")'
