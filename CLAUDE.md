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

## MFCL -> opal (case-studies/swpo-alb)

* 2024 south Pacific albacore diagnostic model (two regions; MFCL 2.2.7.0 in
  `executables/mfcl/2.2.7.0/`). MFCL fits take 15-20 min (~0.5 s/evaluation).
* `scripts/01-build-mfcl-bridge.r`: M0 (evaluate 09.par), M1 (no movement),
  M2 (region 1 only; par transplanted from M1 with `build_region_par()` in
  `scripts/mfcl-edit-helpers.r`), M3 (growth fixed, no age-at-length).
  Fits cached; `SWPO_ALB_REFIT=true` refits.
* The `.frq` dimension line's field 9 is MFCL's month_1 (10 here): the model
  year runs October-September and calendar months >= 10 belong to the next
  model year. `mfcl_to_opal()` handles this.
* data.table: the incident tables have columns `y` and `s`; do not use those
  names as loop variables inside `[...]`.
* MFCL splits each year into periods at the fishing incidents and, only in
  multi-region models, at the movement weeks (`check_for_missing_node()`,
  `rnaux2.cpp`), and shares the year's M equally among the periods. Dropping
  regions therefore changes the within-year timing of M (the M1 -> M2
  transplant difference).

## SS3 -> opal (case-studies/npo-swo)

* 2023 ISC North Pacific swordfish base case: SS3 3.30.20, two sexes (growth
  and M only), four seasons with annual ages 0-15, spawning month 5,
  settlement month 7, selectivity blocks on F1 (1994-97) and F2 (2004-09),
  S6 a recruitment index (SS3 units 33).
* `scripts/01-build-ss3-bridge.r`: 00 (re-run from the supplied par, whose
  checksum line is fixed to 999) to 07 (single sex, blocks as fleets,
  spawning month 1, one season, combined compositions, settlement and index
  month 1), and the same chain with sex-averaged biology (01s-07s). ~2 h with
  Hessians; `NPO_SWO_RERUN=true` reruns, `NPO_SWO_NOHESS=true` skips Hessians.
  Each step starts from the previous par file where the layout allows: from
  the control file (rec devs at 0) SS3 finds a worse minimum of the base case.
  SS3's q mirror (link 2) breaks the gradient; S2 1994-97 has its own q.
* `scripts/02-convert-fit-npo-swo-opal.r`: converts B7 (`config/npo-swo-b7.yaml`;
  `NPO_SWO_VARIANT=sexavg` converts B7s, `npo-swo-b7s.yaml`), replays catch
  timing, fits O0-O2b and the phase D models (four seasons, quarterly catch
  from 03/03s, mid-year weight-at-age); fits cached (`NPO_SWO_REFIT=true`
  refits; `NPO_SWO_SS3_DIR` converts a copy elsewhere). opal's growth is
  fitted to SS3's mid-year lengths (biology stays start-of-year). S6 is an
  index in numbers with age-0 selectivity (age-based double normal in O2).
* Committed SS3 files: inputs, `*.ss_new`, `ss3.par` and `ss_summary.sso`
  per step; the full `Report.sso` only for 07 and 07s (`.gitignore` rules at
  the end of the file). `tests/testthat/test-npo-swo-ss3-bridge.r` reads
  `ss_summary.sso`, not `Report.sso`; `test-npo-swo-opal-cache.r` checks the
  vignette caches against the bridge and the stored fits.
* Quarto: a line starting `@label)` (e.g. a wrapped `(see\n@sec-x)`) is a
  pandoc example-list item and hijacks the cross-reference; keep `@sec-x`
  off the start of a line when followed by `)`.
