# Tape the catch-at-age × PLA products in the composition likelihoods as atomic matrix products

## Summary

In `get_length_like()` and `get_weight_like()` (`R/likelihoods.R`), the
product of predicted catch-at-age and the probability-of-length-at-age matrix
(`catch_ya %*% t(pla)`, and `catch_ya %*% t(rebin_pla)` in the weight
likelihood's legacy path) is taped by RTMB as scalar operations: one multiply
and one add per element of an `[n_obs × n_age] × [n_age × n_len]` product. For
the IO yellowfin AGE-SF model (1,241 length rows, 28 ages, 39 length bins)
this product accounts for about 2.1 of the 2.9 million tape operations and
most of the gradient time.

Taping only this product as one atomic matrix product (RTMB's
`TapeConfig(list(matmul = "atomic"))`, switched on around the product and
restored after) leaves the model unchanged and, in a sandbox study
(`opal-sandbox/benchmarks/gradient-refactor/`, report
`gradient-refactor.html`):

* cut the IO yellowfin tape from 2.88 to 0.90 million operations and the
  gradient from about 20 ms to about 8 ms per call; Opakapaka from 0.19 to
  0.07 million operations and the gradient by about half; the Hessian
  (`obj$he`, used by `sdreport()` and Newton steps) from 10.0 to 3.4 s for IO
  yellowfin and from 193 to 103 ms for Opakapaka;
* gave the same objective to the last bit for IO yellowfin and within one
  rounding unit (1.9e-16 relative) for Opakapaka, gradients within 1e-13
  (scaled by 1 + the gradient's largest element), Hessians within 3e-16, and
  identical REPORT output;
* matched with every parameter estimated (map removed: growth, steepness,
  sigma_r, variance terms, initial F and fixed selectivity elements on the
  tape, so the PLA and weights are AD, not constants): 159 (Opakapaka) and 477
  (IO yellowfin) parameters;
* passed the opal test files listed under *Tests* with the same results as
  the unmodified package;
* reached the same end points in `nlminb` fits (21 Opakapaka starts, 5 IO
  yellowfin starts) in roughly half (Opakapaka) and 40% (IO yellowfin) of the
  time.

The same change applies to the weight-composition likelihood, which has the
same structure; the study tested it by supplying each model's length data as
weight compositions (identical objective) and found the same gains.

## Background: where the gradient time goes

Measured at the stored estimate of the IO yellowfin AGE-SF model, R 4.6.2,
RTMB 2.0, one thread; cumulative cuts of `get_length_like()` with the rest of
the model unchanged:

| Step inside `get_length_like()` | tape ops added | gr ms added (approx.) |
|---|---|---|
| everything outside the length likelihood | 451k (total) | 3.3 (total) |
| `catch_ya %*% t(pla)`, one product per fishery | +2.08M | +15 |
| per-row rebinning (`lf_rebin_s[[...]] %*% pred`) | +0 | ~0 |
| add `lf_addtocomp`, normalise | +0.17M | +1–2 |
| multinomial density | +0.18M | +0.5–1 |

`pla` is an `advector` on the tape even when every growth parameter is mapped
(`RTMB:::isConstant(pla)` is `TRUE`, but the product is still expanded into
scalar operations), so the product's cost does not depend on whether growth is
estimated.

Two related observations from the study, so they are not re-tried:

* `RTMB::TapeConfig(atomic = "enable")` for the whole tape gives almost no
  gain for IO yellowfin: the per-row rebinning products (1,241 small products,
  1,114 of them with a 39 × 39 identity matrix) then become atomic too, and
  their cost replaces the saving. It also changes `mvnorm` taping. Do not set
  it globally.
* Vectorising the per-fishery loop in `do_dynamics()` (matrix operations over
  the fisheries with catch in each step, with plain or atomic sums) changed
  neither the gradient time nor, materially, the tape size. Everything outside
  the composition likelihoods costs about 3 ms of the IO yellowfin gradient.
  Out of scope here.

## Implementation

### 1. Internal helper

Add to `R/utils.R` (or next to the likelihoods):

```r
# Matrix product taped as one atomic RTMB operation, whatever the current
# tape configuration; the previous matmul setting is restored afterwards.
# Outside an AD context (plain numeric inputs) this is an ordinary %*%.
# getOption("opal.atomic_matmul", TRUE) = FALSE gives the scalar product
# (for tests and benchmarks).
.opal_matmul <- function(a, b) {
  if (!isTRUE(getOption("opal.atomic_matmul", TRUE))) return(a %*% b)
  prev <- RTMB::TapeConfig(list(matmul = "NA"))$matmul
  RTMB::TapeConfig(list(matmul = "atomic"))
  on.exit(RTMB::TapeConfig(list(matmul = prev)), add = TRUE)
  a %*% b
}
```

Notes:

* `TapeConfig(list(...))` passes the list to RTMB's `set_tape_config()`;
  `"NA"` leaves a setting unchanged and the call returns the current
  configuration, so the first line reads the previous `matmul` value. Check
  the minimum RTMB version that supports this (tested with RTMB 2.0) and set
  `RTMB (>= x.y)` in `DESCRIPTION` accordingly, or fall back to `a %*% b`
  when the list form is unavailable.
* `on.exit` restores the setting if the product errors.
* RTMB's `%*%` method for `advector` already handles a vector left operand
  (one observation for a fishery): `if (is.null(dim(x)) && length(x) ==
  nrow(y)) x <- t(x)`. The IO yellowfin model has a fishery with a single
  length row and passed.

### 2. `get_length_like()`

Replace

```r
pred_yl  <- catch_ya %*% t(pla)
```

with

```r
pred_yl  <- .opal_matmul(catch_ya, t(pla))
```

Leave everything else as it is: the per-row rebinning, normalisation, the
density calls, the `REPORT` objects and, in particular, the `OBS()` calls.

### 3. `get_weight_like()`

Replace both per-fishery products:

```r
pred_yl <- catch_ya %*% t(pla)          # schemes
pred_yw <- catch_ya %*% t(rebin_pla)    # legacy path
```

with `.opal_matmul(catch_ya, t(pla))` and `.opal_matmul(catch_ya,
t(rebin_pla))`. `rebin_pla <- wf_rebin_matrix %*% pla` is computed once per
call and can stay as it is.

### 4. Things that must not change

* **`OBS()` must stay on the original names.** RTMB's `OBS()` registers an
  observation under the text of its argument (`deparse(substitute(x))`). If
  the length and weight likelihoods are later merged into a shared helper,
  `OBS()` has to be called on `lf_obs_flat`, `lf_obs_prop`, `lf_obs_ints`,
  `wf_obs_*` in the exported functions, not on renamed copies inside the
  helper. In the study, a shared helper that called `OBS(obs_flat)` passed
  every objective, gradient and Hessian check but broke `obj$simulate()`
  (simulated `lf_obs_flat` missing) and OSA residuals ("'observation.name'
  must be in data component"); the tests below catch this.
* AD safety: nothing in the change may depend on a parameter being fixed. The
  helper works for constant and AD operands alike; do not add shortcuts that
  test whether `pla` is constant.
* No global `TapeConfig()` changes, and no change to the per-row rebinning in
  this issue (see *Optional follow-ups*).

## Tests

### Run unchanged (must give the same results as `dev` before the change)

`test-likelihoods-length.R`, `test-likelihoods-schemes.R`,
`test-likelihoods-weight.R`, `test-model-weight-integration.R`,
`test-composition-schemes-reporting.R`, `test-osa-residuals.R`,
`test-model-gradients.R`, `test-opaka-quickstart.R`, `test-golden-opaka.R`,
`test-dynamics*.R`, `test-harvest-rate-consistency.R`,
`test-selectivity-integration.R`, `test-assessment-regressions.R`,
`test-initial-harvest.R`, `test-lambdas.R`, `test-sex-ratio.R`, then the full
suite. In the study these gave 158 tests / 1,051 expectations with 0
failures both before and after (one error in `test-golden-opaka.R`, "The saved
fit runtime identity does not match its payload", occurred with and without
the change and is an environment issue of the study's container).

### Add

1. **Atomic and scalar products agree** (new test in
   `test-likelihoods-schemes.R` or a new `test-atomic-matmul.R`): for the
   Opakapaka quickstart model (legacy bin pooling) and the IO yellowfin
   AGE-SF model (schemes), build the object with
   `options(opal.atomic_matmul = FALSE)` and with the default, and check at
   the stored estimate and at a few perturbed parameter vectors:
   * objective: `expect_equal(..., tolerance = 1e-12)`;
   * gradient: `max(abs(g1 - g0)) / (1 + max(abs(g0))) < 1e-10`;
   * `obj$report()` items identical in names and shapes, values within
     `1e-12` relative;
   * Hessian at the estimate within `1e-8` relative.
2. **Same with every parameter estimated**: as 1, but with `map = list()`
   (growth, steepness, sigma_r, composition and CPUE variance terms on the
   tape). The study's script is `test-ad-safety.R`.
3. **Weight compositions**: as 1 for a weight-composition model. The study
   supplied each model's length data as weight-composition data
   (`wf-data.R`: copy `lf_*` to `wf_*`, `wf_rebin_matrix = diag(n_len)`,
   `n_wt = n_len`, `lf_switch = 0`), which reproduces the length objective
   exactly and covers the scheme and legacy paths at realistic size. A small
   synthetic version (`helper-synthetic.R`) is enough for the test suite; one
   at realistic size belongs in the benchmark.
4. **Tape size**: the IO yellowfin AGE-SF tape (`InfoADFunObject(...)
   $opstack_size`) is below 1.2 million operations (2.88 million before), so a
   later change that reverts to scalar taping is caught.
5. **Simulation and OSA**: covered by the existing tests above; keep them in
   the required set.

## Benchmark and acceptance

Run `dev/run-benchmark-ab.R` (A = `dev`, B = this branch) and report fn, gr,
Hessian and fit times for Opakapaka, IO yellowfin and BET. Expected, from the
study (single machine, paired alternating rounds):

| Model | gr per call | Hessian per call | Fit time (nlminb, 2 passes) |
|---|---|---|---|
| Opakapaka | 1.2 → 0.5 ms | 193 → 103 ms | about half |
| IO yellowfin AGE-SF | 21.5 → 8.5 ms | 10.0 → 3.4 s | about 40% |

(Medians of 6 paired rounds; Hessian 2 rounds for IO yellowfin. Weight
composition versions of both models gave the same ratios.)

Fit counts may change slightly (in the study, ±10% gradient calls per start,
same end points), because the objective and gradient differ by rounding. If
the A/B fit-count gate fails on Opakapaka for this reason, say so in the PR
and update the baseline with `dev/update-benchmark-baseline.R` after review;
the trajectory replays use stored parameter vectors and are not affected.

Acceptance:

* all tests above pass;
* objective, gradient, REPORT and Hessian agreement as in *Add* 1–3;
* `nlminb` fits of Opakapaka and IO yellowfin from the stored starts reach
  the same objective (within 1e-6) as `dev`;
* benchmark shows no regression for any model and the gains above for IO
  yellowfin.

## Optional follow-ups (measured, not recommended now)

Measured on top of the atomic product; each gave at most about 10% more and
was within run-to-run noise in single measurements:

* rebinning a fishery's rows with one product per scheme instead of one per
  row;
* skipping the rebinning product when a scheme's matrix is the identity
  (1,114 of 1,241 IO yellowfin rows); decided from data, so AD-safe;
* normalising a whole block at once;
* taping the block rebinning as atomic too.

They passed the same checks. They are worth revisiting only for models where
the rebinning, normalisation or density dominate after this change.

## References

* Study: `opal-sandbox/benchmarks/gradient-refactor/` (scripts, results,
  report `gradient-refactor.qmd` / `.html`).
* Earlier profile: `opal-sandbox/benchmarks/optimizers/` (its statement that
  the dynamics account for most of the non-composition gradient cost came from
  cut tapes that keep quantities the full model drops, and is corrected by the
  study above).
