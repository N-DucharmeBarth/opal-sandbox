# Case study: North Pacific blue shark

The 2022 ISC stock assessment of blue shark (*Prionace glauca*) in the North
Pacific Ocean (Stock Synthesis 3.30.19) bridged to opal. The assessment's stock
status comes from an ensemble of three SS3 models (ISC 2022, Annex 12); this
case study uses the most heavily weighted of them, S6_base (S6 here). S6 has
two sexes, which opal does not model, so the bridge runs through a sequence of
SS3 models that end in a single-sex SS3 model opal can represent (B2e), and
opal is converted from and compared with that model.

Any model runs contained within this case study represent preliminary model
explorations and should not be used *as is* for the consideration of
management advice. The ISC's stock status advice is based on the three-model
ensemble, of which S6 as supplied in `model-files/ss3/00-s6-base/` is one
member. The assessment report is
<https://isc.fra.go.jp/pdf/ISC22/ISC22_ANNEX12_Stock_Assessment_for_Blue_Shark.pdf>.

## Contents

| Path | What it is |
|---|---|
| `model-files/ss3/00-s6-base/` | S6 as supplied (inputs, `ss.par`), re-run with SS3 3.30.22.1 (`Report.sso`, `CompReport.sso`) |
| `model-files/ss3/01-single-sex/` ... `02e-no-f-ballpark/` | The SS3 bridge, one change per step (see below) |
| `config/npo-bsh-b2e.yaml` | Converter settings for B2e -> opal |
| `scripts/01-build-ss3-bridge.r` | Builds and runs the SS3 bridge from S6; writes the bridge cache |
| `scripts/02-profile-ss3-r0.r` | Likelihood profiles over ln(R0) for S6, B1 and B2e |
| `scripts/03-convert-fit-npo-bsh-opal.r` | Converts B2e, tests the conversion, fits opal (O0-O2) |
| `vignettes/npo-bsh-ss3-to-opal.qmd` | The write-up (rendered HTML and PDF); reads only `vignettes/data/` |

## The SS3 bridge

| Step | Change from the previous step |
|---|---|
| 00 S6 | as supplied, re-run with 3.30.22.1 (identical likelihood) |
| 01 B1 | single sex (`Nsexes = -1`, SSB x female fraction 0.5); female growth, M and weight-length; compositions summed over sex; sex-offset selectivity removed |
| 02a B2a | F18 cubic spline -> double normal |
| 02b B2b | F8's selectivity blocks -> F8 (1973-81) and a new fleet 32 (1982-92) |
| 02c B2c | F19's annual selectivity deviations removed |
| 02d B2d | initial-equilibrium regime offset removed |
| 02e B2e | F-ballpark penalty switched off |

## Order of running

```
scripts/01-build-ss3-bridge.r      # ~45 min with Hessians (NPO_BSH_NOHESS=true skips them)
scripts/02-profile-ss3-r0.r        # ~1.5 h on 4 cores
scripts/03-convert-fit-npo-bsh-opal.r
quarto render vignettes/npo-bsh-ss3-to-opal.qmd
```

The tests `tests/testthat/test-npo-bsh-ss3-bridge.r` and
`tests/testthat/test-ss3-to-opal-npo-bsh.r` check the bridge models
(convergence, Hessian, the structure and data of each step) and the
conversion (reconciliation with SS3, fleet layout, recruitment alignment,
selectivity, and opal's dynamics at SS3's values).

## Credit

Model files: ISC Shark Working Group (2022), *Stock assessment of blue shark
in the North Pacific Ocean through 2020*, ISC/22/ANNEX/12.
