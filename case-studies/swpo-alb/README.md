# Case study: south-west Pacific albacore (MFCL, region 1)

Conversion of the [2024 south Pacific albacore stock assessment](https://meetings.wcpfc.int/node/23119) (WCPFC-SC20-2024/SA-WP-02; diagnostic model from [ofp-sam-alb-2024-diagnostic](https://github.com/PacificCommunity/ofp-sam-alb-2024-diagnostic)), a two-region MULTIFAN-CL (MFCL) model, to opal. An MFCL bridge first reduces the model to region 1 (the southern WCPFC Convention Area); opal itself is not changed, and the comparisons document what opal would need to represent the MFCL model.

The vignette [`vignettes/swpo-alb-mfcl-to-opal.qmd`](vignettes/swpo-alb-mfcl-to-opal.qmd) (rendered [HTML](vignettes/swpo-alb-mfcl-to-opal.html) and [PDF](vignettes/swpo-alb-mfcl-to-opal.pdf)) has the details and results.

Any model runs contained within this case study represent preliminary model explorations and should not be used '*as is*' for the consideration of management advice. Only models that have been presented and accepted to the Western and Central Pacific Fisheries Commission ([WCPFC](https://www.wcpfc.int/)) Scientific Committee should be used as the basis for management advice.

## Summary

* **MFCL bridge.** M0 the diagnostic model (reproduced with MFCL 2.2.7.0), M1 no movement, M2 region 1 only (par transplanted from M1), M3 growth fixed without the conditional age-at-length data. Region 1's spawning biomass roughly doubles at each step (1954: 0.87, 1.7, 3.4, 6.6 Mt). MFCL adds periods at the movement dates only in multi-region models and shares M equally among a year's periods, so dropping region 2 also changes the within-year timing of M.
* **opal at MFCL's values (O0)** reproduces MFCL's initial numbers and recruitment exactly and its spawning biomass within 2.1%; replays put the gap on catch timing (start-of-quarter removal against MFCL's Baranov at each incident).
* **Fits.** With MFCL's selectivity fixed (O1), opal's scale is unbounded: opal predicts the length compositions from start-of-year length-at-age, MFCL at each sample's date with within-year growth. With opal's growth at MFCL's mid-year lengths (O1-mid) the fit is bounded, with B0 at MFCL's Beverton-Holt B0, but its stock is about half MFCL's in 2022: MFCL's initial population is scaled apart from its stock-recruitment relationship (1954 SB about five times its B0). Estimating the initial population, as MFCL does (O1i), brings the stock close to MFCL's (709 against 860 kt in 2022, depletion 0.72 against 0.75). With opal's own selectivity (O2, O2q; three double-normal parameters free) the curves still collapse to knife-edges.
* **Seasonal set-ups** through opal's inputs only: one fishery per season for MFCL's seasonal selectivity (`mfcl_split_seasonal_fisheries()`), and quarterly age classes and time steps with annual recruitment (`mfcl_quarterly_ages()`), which profiles to MFCL's B0 and fits the indices best.

## Running

From any working directory (see `CLAUDE.md` for the R environment):

```bash
Rscript --no-init-file case-studies/swpo-alb/scripts/01-build-mfcl-bridge.r
Rscript --no-init-file case-studies/swpo-alb/scripts/02-compare-mfcl-bridge.r
Rscript --no-init-file case-studies/swpo-alb/scripts/03-convert-fit-swpo-alb-opal.r
Rscript --no-init-file case-studies/swpo-alb/scripts/04-length-comp-diagnostics.r
quarto render case-studies/swpo-alb/vignettes/swpo-alb-mfcl-to-opal.qmd
```

MFCL fits (15-20 min each), opal fits and profiles are cached and only repeated with `SWPO_ALB_REFIT=true`. MFCL 2.2.7.0 (`executables/mfcl/2.2.7.0/mfclo64`, the build supplied with the diagnostic model) is required.

## Contents

| Path | What it is |
|---|---|
| `model-files/mfcl/00-diagnostic/` | the supplied diagnostic model inputs, `08.par`/`09.par` and the plot report |
| `model-files/mfcl/01-...03-.../` | the bridge steps: inputs, start and fitted par files and plot reports |
| `scripts/mfcl-edit-helpers.r` | `.frq`/`.ini`/`.par` edits for the bridge (one region, no movement) |
| `scripts/01-build-mfcl-bridge.r` | builds and fits M0-M3 |
| `scripts/02-compare-mfcl-bridge.r` | bridge comparison; cache `vignettes/data/swpo-alb-mfcl-bridge.rds` |
| `scripts/swpo-alb-opal-setup.r` | conversion, growth, O0 and replays, sourced by 03 and 04 |
| `scripts/03-convert-fit-swpo-alb-opal.r` | O1, O1-mid and the B0 profile; caches `vignettes/data/swpo-alb-opal-*.rds` |
| `scripts/04-length-comp-diagnostics.r` | B0 profiles with mid-year growth, seasonal fisheries and quarterly ages; cache `vignettes/data/swpo-alb-lf-diagnostics.rds` |
| `vignettes/` | the vignette (`.qmd`, `.html`, `.pdf`), `references.bib` and the caches it reads (`data/`) |

The MFCL-to-opal converter and readers are in [`shared/mfcl-to-opal/`](../../shared/mfcl-to-opal).

## Next steps

* O2: opal's own length-based selectivity (with the seasonal fisheries or quarterly ages, and MFCL's groupings as shared parameters).
* A quarterly-age fit with B0 estimated, and a way to report its unfished spawning biomass.
