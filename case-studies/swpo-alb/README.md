# Case study: south-west Pacific albacore (MFCL, region 1)

Conversion of the [2024 south Pacific albacore stock assessment](https://meetings.wcpfc.int/node/23119) (WCPFC-SC20-2024/SA-WP-02; diagnostic model from [ofp-sam-alb-2024-diagnostic](https://github.com/PacificCommunity/ofp-sam-alb-2024-diagnostic)) to opal. The diagnostic model is a two-region MULTIFAN-CL (MFCL) model (region 1: the WCPFC Convention Area south of the equator, including the overlap area; region 2: the EPO). opal is single region, so an MFCL bridge first reduces the model to region 1. The MFCL-to-opal comparisons document what opal would need to represent the MFCL model; opal itself is not changed.

Any model runs contained within this case study represent preliminary model explorations and should not be used '*as is*' for the consideration of management advice. Only models that have been presented and accepted to the Western and Central Pacific Fisheries Commission ([WCPFC](https://www.wcpfc.int/)) Scientific Committee should be used as the basis for management advice.

## The diagnostic model

* 2 regions, 20 fisheries (13 extraction and 2 index fisheries in region 1; 4 extraction and 1 index fishery in region 2); fixed (SEAPODYM) movement four times a year; recruitment split 0.819 / 0.181 between the regions.
* 12 annual age classes, one recruitment a year; the model year runs October-September (`.frq` field 9, MFCL's `month_1` = 10), with recruitment and ageing on 1 October. Catches are taken by quarter (longline, driftnet) or month (troll), 1954-2022.
* Catch-conditioned (Baranov) fishing mortality; catch in numbers except troll and driftnet (weight).
* Length-based cubic-spline selectivity (3-4 nodes), shared across 16 groups, seasonal (4 seasons) for fisheries 1-7.
* Lorenzen M-at-age; spawning potential from maturity-at-length; von Bertalanffy growth with offsets at ages 2-4, estimated with conditional age-at-length data; Beverton-Holt steepness 0.8.
* Length compositions (50 2-cm bins from 30 cm; robust normal likelihood, sample sizes capped at 1000 and divided by Francis weights, samples below 50 fish dropped); two annual WCPFC-CA CPUE indices with time-varying CVs.

## MFCL bridge

| Step | Folder | Change | 
|---|---|---|
| M0 | `00-diagnostic` | the supplied model (`09.par`), evaluated with MFCL 2.2.7.0: reproduces the supplied objective (-41824.603) |
| M1 | `01-no-movement` | movement matrices set to the identity; re-estimated (phase 9 settings) |
| M2 | `02-region-1` | region 1 only: EPO fisheries 14-17 and index 20 dropped, indices 18 and 19 renumbered 14 and 15; par built from M1's estimates; re-estimated |
| M3 | `03-fixed-growth` | growth fixed at M2's estimates and the conditional age-at-length data dropped; re-estimated |

Each MFCL fit takes 15-20 minutes. M2's par file is built from M1's by section: fishery- and region-indexed values are subset, region 1's recruitment is kept by shifting the total population scaling parameter by log(0.818971), and the sections shaped by regions, movement or catch records come from `mfclo64 ... -makepar` for the one-region files (all zero or structural in this model).

With no movement the two regions are closed populations, but in M1 they still share recruitment (total recruitment deviations split by fixed proportions) and the stock-recruitment relationship, so the EPO data still inform region 1's recruitment. M2 removes that link.

The transplant check (`02-region-1/check.par`, M2's structure at M1's estimates) reproduces M1's region 1 recruitment, selectivity, growth, M and predicted catches exactly, but region 1's adult biomass differs by 0.1-1.3% (numbers at the oldest ages by up to 3.7%), growing from a 0.02% per age difference in the initial age structure. M1's region 1 is unaffected by region 2's catch or abundance (tested), so the cause is in the one-region MFCL structure; it is not yet identified.

## opal

`shared/mfcl-to-opal/` reads MFCL `.par`, `.frq` and plot reports in base R and builds opal inputs (`mfcl_to_opal()`): four seasons within MFCL's model year, catch by fishery and season, the two indices (log-normal with MFCL's time-varying CVs), annual length compositions (MFCL's caps, minimum sample size and divisors), and biology at age from MFCL's report, so opal's spawning biomass is MFCL's adult biomass. `mfcl_match_recruitment()` sets opal's B0 from MFCL's Beverton-Holt curve and its initial and recruitment deviations so that opal starts from MFCL's numbers-at-age and recruits MFCL's recruitment.

`scripts/03-convert-fit-swpo-alb-opal.r` builds O0 (opal at MFCL's values), the replays and O1 (opal fitted with MFCL's selectivity). The replays are standalone simulations with MFCL's estimates that add, one at a time, the two dynamics opal lacks: catch removed at each fishing incident's time (Baranov, M and F together) rather than at the start of each quarter, and seasonal selectivity.

### Differences between MFCL and opal

Features of the MFCL model that opal does not have, and how they are handled here:

| MFCL | opal | Here |
|---|---|---|
| Multiple regions, movement | one region | MFCL bridge to region 1 |
| Conditional age-at-length likelihood | none | dropped in M3, growth fixed |
| VB growth with independent mean lengths at ages 2-4 | VB, Richards, cessation, age-specific K | age-specific K at ages 1-4 fitted to MFCL's mean lengths |
| SD of length-at-age log-linear in mean length | CV or SD linear in mean length | SD linear, matched at ages 1 and 12 |
| Cubic-spline selectivity, selectivity groups | logistic, double normal, double Richards; no mirrors | MFCL's selectivity-at-age as `sel_fa_external` (O0, O1) |
| Seasonal selectivity (fisheries 1-7) | selectivity by year | catch-weighted average over seasons; effect measured by replay |
| Catch-conditioned Baranov at each incident's time (several per quarter for troll) | Pope removal at the start of each season | effect measured by replay |
| Length compositions per fishing incident (robust normal) | annual compositions (multinomial) | incident samples summed by fishery and model year |
| Index abundance at the incident's time | at the start of the season | as opal |
| Recruitment from the previous year's spawning potential (lag 1) | from the same year's | O0 imposes MFCL's recruitment; O1 as opal |
| Last two recruitments fixed at the mean | free deviations | last two deviations fixed at 0 |

## Results so far

MFCL bridge, region 1 (`scripts/02-compare-mfcl-bridge.r`):

| Step | Objective | Max gradient | SB 1954 (kt) | SB 2022 (kt) | SB / SB(F=0), 2022 |
|---|---|---|---|---|---|
| M0 | -41824.60 | 5.5e-7 | 871 | 251 | 0.50 |
| M1 | -41813.17 | 3.3e-6 | 1716 | 375 | 0.57 |
| M2 | -36315.10 | 6.2e-7 | 3413 | 460 | 0.62 |
| M3 | -37149.05 | 7.4e-6 | 6592 | 859 | 0.75 |

(M0's values are region 1 of the two-region model. M1 and M3 stop at MFCL's 2000-evaluation limit with the line search stalled.) Each step roughly doubles the scale: the stock size in this catch-conditioned model is weakly determined, and each source of information removed (movement, the shared recruitment with the EPO, the conditional age-at-length data) lets it drift up.

opal (`scripts/03-convert-fit-swpo-alb-opal.r`), against M3:

* O0 reproduces MFCL's initial numbers-at-age (5e-7) and recruitment (2e-14); its spawning biomass is within 2.1% of MFCL's (+1.6% in 2022) and its 2022 depletion 0.763 against 0.751.
* Replays (spawning biomass relative to MFCL, largest annual difference): opal timing with averaged selectivity 2.1% (and within 0.06% of O0, which checks the replay), with seasonal selectivity 2.2%; MFCL's Baranov timing 0.55% (averaged) and 0.53% (seasonal). Catch timing explains almost all of the O0 difference; seasonal selectivity hardly affects spawning biomass. The remaining 0.5% is likely the weight-at-age used for catches in weight (largest in the mid-1990s, after the driftnet years), not confirmed.
* Growth in opal's form (age-specific K at ages 1-4) matches MFCL's mean lengths with an RMSE of 0.52 cm and SDs within 0.47 cm.
* O1 (MFCL's selectivity fixed; B0, recruitment deviations 1955-2020 and catchabilities estimated) runs to the upper bound of log(B0). A profile over B0 shows the CPUE (best near MFCL's scale) and the recruitment penalty against the length compositions, which keep improving towards an unexploited stock. With annual multinomial compositions and MFCL's Francis divisors, the length data do not bound the scale in opal.

### Why the length compositions favour an unexploited stock

* MFCL predicts each length sample at its date, from lengths that grow through the year (from M3's `length.fit`: age 1 is 47.2 cm in model month 2 and 54.6 cm in month 7, age 5 88.7 to 91.9 cm); opal predicts them from start-of-year length-at-age (45.5 and 88.0 cm). opal's predicted catches are too small, which the fit offsets with more large, old fish, i.e. less fishing. The longline fisheries drive it (fishery 3 alone gains 16 likelihood units at a stock 20 times larger; the troll fishery 11 prefers MFCL's scale).
* With opal's growth fitted to MFCL's mid-year length-at-age (inputs only), the B0 profile has its minimum at MFCL's scale and the CPUE fit improves (15.7 against 44.7 at MFCL's scale). O1-mid (O1 with this growth) converges at log(B0) 13.95 against MFCL's 14.04.
* Splitting fisheries 1-7 into one fishery per season (`mfcl_split_seasonal_fisheries()`, the seasonal selectivity exactly, 36 fisheries) changes little: with start-of-year growth the profile is still unbounded; with mid-year growth its minimum is just below MFCL's scale.
* A quarterly-age opal model (`mfcl_quarterly_ages()`: 48 quarterly ages, 276 quarterly steps, recruitment kept annual by fixing the other quarters' deviations at -30) gives each quarter MFCL's length-at-age, weight and selectivity and quarterly length compositions; at MFCL's values it reproduces MFCL's spawning biomass within 2.4%. Its B0 and steepness lose their usual meaning (one recruitment per four steps). Its B0 profile is in `scripts/04-length-comp-diagnostics.r` (to be run).

## Next steps

* Run `04-length-comp-diagnostics.r` (quarterly-age profile) and compare it with the annual mid-year variants.
* O2: opal's own length-based selectivity, with MFCL's groupings as shared parameters.
* The 0.1-1.3% residual of the M2 transplant check.
* The vignette `vignettes/swpo-alb-mfcl-to-opal.qmd` from the caches.

## Running

From any working directory (see `CLAUDE.md` for the R environment):

```bash
Rscript --no-init-file case-studies/swpo-alb/scripts/01-build-mfcl-bridge.r
Rscript --no-init-file case-studies/swpo-alb/scripts/02-compare-mfcl-bridge.r
Rscript --no-init-file case-studies/swpo-alb/scripts/03-convert-fit-swpo-alb-opal.r
```

The MFCL fits and opal fits are cached and only repeated with `SWPO_ALB_REFIT=true`. MFCL 2.2.7.0 (`executables/mfcl/2.2.7.0/mfclo64`, the build supplied with the diagnostic model) is required.

## Contents

| Path | What it is |
|---|---|
| `model-files/mfcl/00-diagnostic/` | the supplied diagnostic model inputs, `08.par`/`09.par` and the plot report |
| `model-files/mfcl/01-...03-.../` | the bridge steps: inputs, start and fitted par files and plot reports |
| `scripts/mfcl-edit-helpers.r` | `.frq`/`.ini`/`.par` edits for the bridge (one region, no movement) |
| `scripts/01-build-mfcl-bridge.r` | builds and fits M0-M3 |
| `scripts/02-compare-mfcl-bridge.r` | bridge comparison, cache `vignettes/data/swpo-alb-mfcl-bridge.rds` |
| `scripts/swpo-alb-opal-setup.r` | conversion, growth, O0 and replays, sourced by 03 and 04 |
| `scripts/03-convert-fit-swpo-alb-opal.r` | O1, O1-mid and the B0 profile; caches `vignettes/data/swpo-alb-opal-*.rds` |
| `scripts/04-length-comp-diagnostics.r` | length-composition diagnostics: B0 profiles with mid-year growth, seasonal fisheries and quarterly ages; cache `vignettes/data/swpo-alb-lf-diagnostics.rds` |
