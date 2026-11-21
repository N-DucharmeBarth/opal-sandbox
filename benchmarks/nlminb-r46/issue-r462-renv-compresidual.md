# Move the renv environment to R 4.6.2 and drop compResidual

Two changes, best done as two PRs in this order so each can be checked on its
own: (1) remove `compResidual`, which removes most of the custom CI
environment code; (2) move `renv.lock` from R 4.5.2 to R 4.6.2.

## Why

* R 4.6.2 fixes a long-standing bug in bounded `nlminb()` (R NEWS: "Fix long
  standing (Fortran compiler optimization related) bug in nlminb with
  constraints", PR#19086). In a jitter study of the Opakapaka quickstart and
  IO yellowfin AGE-SF fits (`opal-sandbox/benchmarks/nlminb-r46/`), R 4.6.2
  needed fewer gradient calls from 29 of 30 starts (median fit time
  Opakapaka 0.71 → 0.22 s, IO yellowfin 14.2 → 9.4 s). The model arithmetic
  is identical between the two versions; the difference is in the optimiser's
  path.
* `compResidual` is a `Suggests` dependency used in one place (a guarded
  chunk of `dev/vignettes/bet.qmd`), but it needs two CI workarounds: a
  source install outside renv (renv 1.2.x fails to stage its GitHub
  sub-directory) and TMB contrib OSA headers grafted into the TMB
  installation.

## Assessment (state of `dev` at `87c5494`)

### compResidual: safe to remove

* `DESCRIPTION`: `Suggests` only. Not in `NAMESPACE`.
* No use in `R/`, `tests/`, `vignettes/`, `data-raw/` or `inst/`.
* One use: `dev/vignettes/bet.qmd`, section "Composition OSA residuals"
  (`compResidual::resMulti()` row by row, inside
  `if (!requireNamespace("compResidual", ...))`). `dev/vignettes/` is not
  rendered by any workflow.
* opal has its own composition residuals: `opal_osa()` (RTMB
  `oneStepPredict()`), tested for every composition family in
  `test-osa-residuals.R`. They are a different method from `resMulti()`, so
  the BET plot would change.
* CI code that exists only for compResidual (`.github/scripts/setup-r-env.R`):
  the `composition` phase, the compResidual record removed from the lockfile
  in the `restore` phase, `graft_osa()` (its comment: "Required by
  compResidual"; nothing else in the repository includes
  `OSA_multivariate_dists` or `distr.hpp`), `set_locked_provenance()`, and
  `compResidual` in the `verify` list. `locked_github_downloader()` is also
  used for SparseNUTS (GitHub source) and stays.

### R 4.6.2: safe, with known follow-ups

Evidence:

* The R-CMD-check `release` job (R 4.6.x) passes on `dev` at `87c5494`
  (run 42081049289). It restores the R 4.5.2 lockfile on R 4.6 with one
  override: Matrix 1.7-6 instead of the locked 1.7-1, then rebuilds TMB and
  RTMB (`matrix` phase of `setup-r-env.R`). The exact patch version that
  "release" resolved to was not checked.
* In the sandbox, opal and its runtime dependencies (102 packages) were built
  from source under R 4.6.2. The models give the same objective and gradient
  under R 4.5.2 and R 4.6.2. Note: that library had newer versions of several
  packages than the lockfile (RTMB 2.0, TMB 1.9.23, Matrix 1.7-6, RTMBdist
  1.2.1, ...), not the locked ones.

Known issues to handle:

1. **Matrix.** Locked 1.7-1 can segfault on R 4.6 (comment in
   `setup-r-env.R`). The lockfile must move to Matrix ≥ 1.7-6, and TMB and
   RTMB must be built against it.
2. **Fits change.** Bounded `nlminb()` takes different paths under R 4.6.2.
   Opakapaka has two minima (objectives 1831.669 and 1832.412); from jittered
   starts each version reaches either, so anything that stores or checks a
   fit from a non-stored start can change. From the stored starts both
   versions reached the same end points in the study (Opakapaka 1832.412, IO
   yellowfin 15601.99).
3. **Stored fits and objects.** Fits in `inst/extdata/` and `data/` record
   `r_version = "4.5.2"`. Under R 4.6.2, `read_opal_fit(integrity =
   "portable")` warns ("serialized under a different R version; verifying by
   rebuild") and verifies by rebuilding; `integrity = "exact"` fails. They
   should be regenerated under R 4.6.2.
4. **Benchmark baselines and trajectories.** `tests/testthat/_benchmark/`
   (`baseline.rds`, `trajectories/*.rds`) and any self-test calibration were
   recorded under R 4.5.2; fit counts will drop (Opakapaka gradient calls per
   fit about 500 → 140 in the study). A/B comparisons run both sides on the
   same R and are unaffected, but stored baselines need refreshing.
5. **CI time.** The `release` job's environment step took 34 min against 4
   min for the `renv` job on the same run, consistent with packages being
   built from source on R 4.6 (the setup action notes P3M had no jammy
   binaries for R 4.6). Whether Posit Package Manager now has jammy binaries
   for R 4.6 could not be checked from the sandbox. Without binaries every
   cold-cache job pays this, and the devcontainer build too.
6. **RTMB version.** The lockfile has RTMB 1.9 / TMB 1.9.21. The gradient
   refactor issue (`issue-atomic-composition-products.md`) was measured with
   RTMB 2.0. If that change goes in, decide here whether to move RTMB to 2.0
   (and check `TapeConfig(list(matmul = ...))` on the version locked).
7. **Packages not exercised under R 4.6.2 in the sandbox.** About 155
   development packages (devtools, pkgdown, rmarkdown, quarto, profvis, ...)
   were not built there. The CI `release` job installs the full lockfile on
   R 4.6 and passes R CMD check, which covers installation but not every
   tool (pkgdown site, README render, benchmark and self-test workflows).

## PR 1: remove compResidual

1. `DESCRIPTION`: remove `compResidual` from `Suggests`.
2. `renv.lock`: remove the record (`renv::remove("compResidual")` then
   `renv::snapshot()`, or edit the JSON; check that no other record is
   dropped or changed).
3. `dev/vignettes/bet.qmd`: replace the "Composition OSA residuals" chunk with
   one based on `opal_osa()` (or delete the section). Say in the text that
   the residuals are now OSA residuals from RTMB's `oneStepPredict()`.
4. `.github/scripts/setup-r-env.R`:
   * remove the `composition` phase and `graft_osa()` (both calls);
   * `restore` phase: restore from `renv.lock` directly (no trimmed lockfile);
   * remove `set_locked_provenance()`;
   * `verify`: drop `compResidual` from `required` and the distr.hpp check;
   * `phases <- c("tmb", "restore", "matrix", "verify")`; the `tmb` phase can
     probably be folded into `restore` once nothing is grafted into TMB, but
     keep it if TMB must still be restored first.
5. `tests/testthat/test-development-environment.R`: remove the
   `set_locked_provenance` test; keep the `locked_github_downloader` tests.
6. `.github/actions/setup-opal-env/action.yml`: update the comment on
   `RENV_CONFIG_PAK_ENABLED`.
7. `.github/CONTRIBUTING.md`: remove the TMB OSA-header and compResidual
   workaround text.
8. `NEWS.md`: note the removal.
9. Bump the renv cache key if needed (`renv-v2` → `renv-v3`) so a cache with a
   grafted TMB is not reused.

Checks: all workflows green (R-CMD-check renv and release, test-coverage,
pkgdown, document, benchmark, selftest, copilot-setup-steps); a cold-cache
run of the setup action (no restore key hit) succeeds; `dev/vignettes/bet.qmd`
renders locally.

## PR 2: move renv to R 4.6.2

1. Install R 4.6.2 (rig), restore the current lockfile with Matrix 1.7-6, then
   `renv::snapshot()` so `R.Version` becomes `4.6.2` and Matrix 1.7-6 is
   recorded. Decide whether to update other packages at the same time; the
   smallest change is R + Matrix (+ the TMB/RTMB rebuild, which does not
   change their versions). Updating RTMB to 2.0 is a separate decision (item
   6 above).
2. `.github/scripts/setup-r-env.R`: remove the `matrix` phase (the lockfile
   now has a compatible Matrix), or keep it only for R ≥ 4.7 if wanted.
3. `.github/workflows/R-CMD-check.yaml`: keep `renv` and `release`; consider
   adding an `oldrel` (R 4.5.x) job while users are still on 4.5, since
   `DESCRIPTION` says `R (>= 3.5)` and R < 4.6.2 still has the old `nlminb()`.
4. Devcontainer: nothing to change; the Dockerfile reads `R.Version` from
   `renv.lock`. Rebuild and time it.
5. Regenerate under R 4.6.2:
   * stored fits and objects in `inst/extdata/` and `data/` (`data-raw/`
     scripts; `dev/generate-object-portability.R`); check each refit reaches
     the same objective as the current stored fit (within 1e-6) and say so in
     the PR, or explain any difference;
   * benchmark baseline and trajectories (`dev/update-benchmark-baseline.R`,
     `dev/update-benchmark-trajectories.R`);
   * self-test calibration if it is version-specific
     (`dev/run-selftest-calibration.R`);
   * rendered vignettes and README if they show fit results.
6. Tests that compare fits with stored values: check none depends on
   `nlminb()`'s path from a non-stored start (Opakapaka's two minima). Run the
   full suite on R 4.6.2 and on R 4.5.x.
7. `NEWS.md`: note the R version, why (the `nlminb()` fix), and that fits
   from R < 4.6.2 can differ in path and, for multimodal models, in end
   point.

Checks: all workflows green; stored-fit objectives unchanged (or explained);
benchmark A/B of `dev` vs this branch on the new R shows no regression; CI
environment time recorded before and after.

## Out of scope

* Changing opal's code for R 4.6 (none needed in the study).
* The atomic matrix product change (separate issue).
