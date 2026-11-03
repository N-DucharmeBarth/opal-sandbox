# Case study: eastern Pacific Ocean bigeye tuna

An exploratory Stock Synthesis (SS3) assessment model of bigeye tuna
(*Thunnus obesus*) in the eastern Pacific Ocean, prepared by IATTC staff for
the 14th meeting of the IATTC Scientific Advisory Committee (SAC-14, 2023),
bridged to opal. The model is quarterly: each SS3 "year" is a quarter
(17-204) and ages are quarters (0-40). It has 22 fisheries (longline by area
in numbers and in weight, floating-object purse seine by area, a discard
fleet and unassociated purse seine) and one longline index, with length
compositions for the floating-object fisheries and generalised size
compositions for the longline fisheries and the index.

The baseline SS3 model is already a simplified model. The models in the
IATTC's risk analysis for EPO bigeye have two sexes and time-varying
selectivity; this exploratory model has one sex and time-invariant
selectivity.

The SS3 model differs from opal in one structural feature, an
initial-equilibrium regime offset, so the bridge has a single step (B1)
that removes it, and opal is converted from and compared with B1. Growth
needs no bridge step: SS3's growth cessation model (GrowthModel 8) and its
SD-linear growth variability (CV_Growth_Pattern 2) are both in opal (0.1.2).

Any model runs contained within this case study represent preliminary model
explorations and should not be used *as is* for the consideration of
management advice. The model as supplied is an exploratory model; the IATTC's
advice on the stock is based on its benchmark assessments.

## Contents

| Path | What it is |
|---|---|
| `model-files/ss3/00-base/` | The model as supplied (`BET-EPO.dat` and `BET-EPO.ctl`, renamed `data.ss` and `control.ss`; `ss3-supplied.par`), re-run with SS3 3.30.22.1. The supplied run's likelihood components (an SS3 3.30.22 beta build) are in `supplied-likelihood-3.30.22beta.csv` |
| `model-files/ss3/01-no-regime/` | B1: the initial-equilibrium regime offset removed |
| `model-files/ss3/02-base-init-f/` | B0-F: B0 with the initial F of LL-n-A2 and OBJ-A5 estimated (sensitivity) |
| `model-files/ss3/03-no-regime-init-f/` | B1-F: B1 with the same initial Fs estimated (sensitivity) |
| `config/epo-bet-b1.yaml` | Converter settings for B1 -> opal |
| `scripts/01-build-ss3-bridge.r` | Re-runs 00, builds and runs 01 and the initial-F sensitivities 02 and 03; writes the bridge cache |
| `scripts/02-convert-fit-epo-bet-opal.r` | Converts B1, tests the conversion, fits opal (O0-O2, and O1-F with the initial F estimated) |
| `vignettes/epo-bet-ss3-to-opal.qmd` | The write-up (rendered HTML and PDF); reads only `vignettes/data/` |

## The SS3 bridge

| Step | Change from the previous step |
|---|---|
| 00 B0 | as supplied, re-run with 3.30.22.1 (identical likelihood) |
| 01 B1 | initial-equilibrium regime offset (block year 16, estimated at -0.41 on ln R0) removed |

## Initial F

SS3 fixes the initial F of the two fleets with an initial equilibrium catch
(LL-n-A2 and OBJ-A5; the catch has lambda 0, so the Fs are not tied to it).
Estimating them:

| Run | Regime offset | Initial F LL-n-A2 / OBJ-A5 (per quarter) | NLL / objective | Depletion 1975 Q1 / 2021 Q4 |
|---|---|---|---|---|
| B0 | -0.41 | 0.0033 / 1.8e-7 (fixed) | 271.20 | 0.649 / 0.221 |
| B0-F | -0.41 | 0.0040 (SE 0.072) / ~0 | 271.20 | 0.647 / 0.221 |
| B1 | off | 0.0033 / 1.8e-7 (fixed) | 274.60 | 0.868 / 0.249 |
| B1-F | off | ~0 / 0.045 (SE 0.022) | 272.40 | 0.640 / 0.226 |
| opal O1 | none | 0.0033 / 1.8e-7 (fixed) | 6369.85 | 0.979 / 0.333 |
| opal O1-F | none | ~0 / 0.0034 | 6369.83 | 0.976 / 0.333 |

With the regime offset on the initial F is not estimable. Without it (B1-F)
OBJ-A5's initial F recovers 2.2 of the 3.4 likelihood units the offset is
worth and brings the trajectory back close to B0's. In opal (O1-F) the
initial F stays near 0 and changes almost nothing; both of its starts (O1's
values and B1-F's) reach the same optimum.

## Gaps in opal found here

* No within-step catch timing: opal removes each quarter's catch at its
  start, which at SS3's values leaves spawning biomass about 40% above SS3's
  by 2021 (a mid-quarter removal would be within 5% in 2021 and 8% at most;
  SS3's catch equation reproduces SS3 to 1e-5).
* One weight-at-age vector for all fisheries, at one time in the step; SS3
  uses each fleet's catch weight-at-age.
* No unpenalised offset on the initial equilibrium recruitment (SS3's regime
  offset). Removing it (B1) moves SSB by +35% in 1975 and +14% in 2021; it is
  the reason opal is converted from B1 rather than the supplied model.
* No combined age and length selectivity (SS3 age pattern 10 sets age 0 to
  zero on top of each fleet's length curve).
* opal's objective is larger than the sum of the weighted terms it reports
  whenever the compositions are in the likelihood (by about 2,300 units for
  O1 here, and about 270 for the NPO blue shark O1); the gap disappears with
  the compositions' lambda at 0, so part of the composition likelihood is not
  in the reported `lp_lf`.

## Time

The data file gives model quarters only. This case study takes quarter 17
as 1975 Q1, so quarter 204 is 2021 Q4 (quarter q is year
1975 + (q - 17) / 4). Calendar years in the vignette's figures rest on that
assumption.

## Order of running

```
scripts/01-build-ss3-bridge.r      # ~1 h per step with Hessians (EPO_BET_NOHESS=true skips them)
scripts/02-convert-fit-epo-bet-opal.r
quarto render vignettes/epo-bet-ss3-to-opal.qmd
```

The tests `tests/testthat/test-epo-bet-ss3-bridge.r` and
`tests/testthat/test-ss3-to-opal-epo-bet.r` check the bridge models
(convergence, the change made in each step, unchanged data, agreement with
the supplied run) and the conversion (reconciliation with SS3, fleet layout,
growth cessation and SD-at-age against SS3, recruitment alignment, the size
compositions as a second bin scheme, selectivity, and opal's dynamics at
SS3's values).

## Credit

Model files: IATTC staff (2023), exploratory stock assessment model of
bigeye tuna in the eastern Pacific Ocean, prepared for the 14th meeting of
the IATTC Scientific Advisory Committee (SAC-14).
