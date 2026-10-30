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
| `config/epo-bet-b1.yaml` | Converter settings for B1 -> opal |
| `scripts/01-build-ss3-bridge.r` | Re-runs 00, builds and runs 01; writes the bridge cache |
| `scripts/02-convert-fit-epo-bet-opal.r` | Converts B1, tests the conversion, fits opal (O0-O2) |
| `vignettes/epo-bet-ss3-to-opal.qmd` | The write-up (rendered HTML and PDF); reads only `vignettes/data/` |

## The SS3 bridge

| Step | Change from the previous step |
|---|---|
| 00 B0 | as supplied, re-run with 3.30.22.1 (identical likelihood) |
| 01 B1 | initial-equilibrium regime offset (block year 16, estimated at -0.41 on ln R0) removed |

## Gaps in opal found here

* No within-step catch timing: opal removes each quarter's catch at its
  start, which at SS3's values leaves spawning biomass about 40% above SS3's
  by 2021 (a mid-quarter removal would be within 5%).
* One weight-at-age vector for all fisheries, at one time in the step; SS3
  uses each fleet's catch weight-at-age.
* No unpenalised offset on the initial equilibrium recruitment (SS3's regime
  offset). Removing it (B1) moves SSB by +35% in 1975 and +14% in 2021; it is
  the reason opal is converted from B1 rather than the supplied model.
* No combined age and length selectivity (SS3 age pattern 10 sets age 0 to
  zero on top of each fleet's length curve).

## Time

The data file gives model quarters only. This case study takes quarter 17
as 1975 Q1, so quarter 204 is 2021 Q4 (quarter q is year
1975 + (q - 17) / 4). Calendar years in the vignette's figures rest on that
assumption.

## Order of running

```
scripts/01-build-ss3-bridge.r      # ~1 h with Hessians (EPO_BET_NOHESS=true skips them)
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
