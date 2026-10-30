# Notes for Claude Code

## R environment in cloud containers

R is not preinstalled. `shared/setup/install-r-opal.sh` installs R 4.5,
opal (dev), r4ss and Quarto into `/opt/renv45` through conda-forge and
`git clone` (CRAN and GitHub HTTP are blocked). It takes 20-30 minutes, so run
it in the background:

```bash
bash shared/setup/install-r-opal.sh > /tmp/install-r.log 2>&1
```

In each new shell afterwards:

```bash
export PATH=/opt/renv45/bin:$PATH LANG=C.UTF-8 LC_ALL=C.UTF-8
export QUARTO_SHARE_PATH=/opt/renv45/share/quarto QUARTO_DENO=/opt/renv45/bin/deno \
       QUARTO_ESBUILD=/opt/renv45/bin/esbuild QUARTO_DART_SASS=/opt/renv45/bin/sass \
       QUARTO_PANDOC=/opt/renv45/bin/pandoc QUARTO_TYPST=/opt/renv45/bin/typst
```

The repository's `.Rprofile` activates renv, which tries to reach CRAN. Run R
with `--no-init-file` or from outside the repository root. Scripts find the
project root with `this.path::this.proj()`, so they still work from anywhere.

## Layout

See README.md. Each stock has its own folder under `case-studies/`; code used
by more than one case study goes in `shared/` (`ss3-to-opal/` converter,
`extractors/{ss3,mfcl}/`, `setup/`). `simulations/` holds the opal simulation
studies and `tests/testthat/` the tests of the shared code, which
`.github/workflows/tests.yaml` runs on pull requests and pushes to `main`.

## SS3 -> opal (shared/ss3-to-opal, case-studies/io-yft)

* `shared/ss3-to-opal/`: the SS3-to-opal converter and fit helpers.
* `case-studies/io-yft/scripts/01`-`03`: fits that write the caches in
  `case-studies/io-yft/vignettes/data/`.
* `case-studies/io-yft/scripts/04-export-io-yft-for-opal.r`: exports fit C2 as
  package data (including SS3's `age_selex` parameters in `io_yft_reference`), a stored
  fit and data docs for the opal package, under `opal-export/` (laid out
  like the opal repository; see its README). The age-based fits, the
  vignette cache and the two vignettes there are made inside opal's
  `age-selex` branch by its `data-raw/` scripts and copied back;
  `opal-export/patches/` holds the same changes as `git am` patches.
  Files in `opal-export/` and `opal-patches/` are snapshots of opal and still
  name the pre-reorganisation paths (`code/rtmb/...`).
* `case-studies/io-yft/scripts/05-io-yft-ss3-to-opal-workflow.r`: end-to-end
  conversion, tests and fit C2 with the settings as an R list; companion to the step-by-step
  guide `vignettes/ss3-to-opal-workflow.qmd` (HTML + PDF), which shows its
  chunks via `knitr::read_chunk()` and reads the cache it writes.
* `.gitignore` ignores everything not explicitly whitelisted; add new file
  types there.

## SS3 -> opal (case-studies/epo-bet)

* Quarterly SS3 model (each SS3 "year" is a quarter, 17-204 = 1975 Q1 -
  2021 Q4; ages 0-40 quarters) with growth cessation (GrowthModel 8) and
  SD-linear growth variability (CV_Growth_Pattern 2), converted with
  `cfg$growth$model: cessation` (needs opal >= 0.1.2).
* `scripts/01-build-ss3-bridge.r`: re-runs 00-base and builds 01-no-regime
  (~1 h with Hessians; a step with `admodel.cov` is skipped unless
  `EPO_BET_RERUN=true`). 01 starts from 00's par file: from the control
  file's INIT values the recruitment deviations are 0 and the stock crashes.
* `scripts/02-convert-fit-epo-bet-opal.r`: conversion, O0 and the opal fits
  (cached; `EPO_BET_REFIT=true` refits), the catch-timing replays, and the
  cache for `vignettes/epo-bet-ss3-to-opal.qmd`. The longline and index
  generalised size compositions are opal's second composition bin scheme.
