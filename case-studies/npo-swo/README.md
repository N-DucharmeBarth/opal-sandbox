# Case study: North Pacific swordfish

The 2023 ISC stock assessment of swordfish (*Xiphias gladius*) in the North
Pacific Ocean (Stock Synthesis 3.30.20; ISC 2023, Annex 11) bridged to opal.
The base case has two sexes and four seasons (quarters) with annual ages.
opal has one sex, and the SS3-to-opal converter handles one season, so the
bridge runs through a sequence of SS3 models that ends in a single-sex,
annual SS3 model opal can represent (B7). opal is converted from and compared
with that model.

Any model runs contained within this case study represent preliminary model
explorations and should not be used *as is* for the consideration of
management advice. The ISC's stock status advice is based on the base case
as supplied in `model-files/ss3/00-base/`. The assessment report is
ISC/23/ANNEX/11 (July 2023).

## Contents

| Path | What it is |
|---|---|
| `model-files/ss3/00-base/` | The base case as supplied (`swo2023_v004.dat`, `swo2023_v007.ctl` renamed `data.ss`, `control.ss`; `forecast.ss` is the supplied run's `forecast.ss_new`; `ss-supplied.par`), re-run with SS3 3.30.22.1 from the supplied par file |
| `model-files/ss3/01-single-sex/` ... `07-index-month-1/` | The SS3 bridge, one change per step (see below) |
| `model-files/ss3/01s-sex-averaged/` | Sensitivity: B1 with sex-averaged M and growth |
| `config/npo-swo-b7.yaml` | Converter settings for B7 -> opal |
| `scripts/01-build-ss3-bridge.r` | Builds and runs the SS3 bridge from the base case; writes the bridge cache |
| `scripts/02-convert-fit-npo-swo-opal.r` | Converts B7, tests the conversion, replays SS3's dynamics with opal's catch timing, fits opal (O0-O2 and two diagnostics) |
| `vignettes/npo-swo-ss3-to-opal.qmd` | The write-up (rendered HTML and PDF); reads only `vignettes/data/` |

## The SS3 bridge

| Step | Change from the previous step |
|---|---|
| 00 B0 | as supplied, re-run with 3.30.22.1 from the supplied par file (identical likelihood, 289.2996) |
| 01 B1 | single sex (`Nsexes = -1`, SSB x female fraction 0.5); female growth, M and weight-length; the sex-combined compositions keep one copy of their vector |
| 01s B1s | (sensitivity, from B0) as B1 with M and growth averaged over the sexes, weighted by unfished numbers-at-age |
| 02 B2 | selectivity blocks become fleets: F1, F12 and S2 in 1994-97 (F1's block), F2, F15, F16 and F19 in 2004-09 (F2's block); S2 1994-97 has its own floated q |
| 03 B3 | spawning in month 1 (season 1) instead of month 5 |
| 04 B4 | one season: catches summed over the quarters |
| 05 B5 | one length composition per fleet and year (quarterly vectors combined, Nsamp summed); same estimates as B4 |
| 06 B6 | recruits settle in month 1 (age 0 at the start of the year) instead of month 7 |
| 07 B7 | indices observed in month 1; no change for the fitted indices (months 1 and 4 are both in SS3's first subseason) |

Notes on the bridge:

* The supplied par file's checksum line holds 6.98 instead of 999, which SS3
  3.30.22 rejects; `01-build-ss3-bridge.r` writes a corrected copy.
* Estimated from the control file (whose INIT values are the supplied
  estimates, but whose recruitment deviations start at zero), the base case
  converges to a worse minimum (332.18 against 289.30). Each step therefore
  starts from the previous step's par file where the parameter layout allows
  (B2-B7), or from the control file plus the previous step's recruitment
  deviations (`read_recdevs`; B1, B1s), and is re-run from its own par file
  until the maximum gradient is below 1e-3.
* SS3's q mirror (link 2) does not carry the gradient of the shared q: with
  S2 1994-97 mirroring S2's estimated q, B2 reproduces B1's likelihood at
  B1's estimates but cannot converge. S2 1994-97 has its own floated q.
* SS3 does not allow spawning after the start of the season in which recruits
  settle, so the one-season model needs spawning in month 1 (B3) before the
  seasons can be merged (B4).

## Gaps in opal found here

See the vignette's "Limitations and differences from SS3" for the full list.

* **Two sexes** (sex-specific growth and M here; the compositions are
  sex-combined, so no data are lost in the single-sex bridge).
* **Seasons with annual ages.** opal has `n_season`, but predicts
  compositions from the year's total catch-at-age, uses one length- and
  weight-at-age for the whole year, and spawns and recruits at the start of
  season 1. SS3 here spawns in month 5, settles recruits in month 7 and
  observes quarterly compositions with each quarter's lengths. The SS3-to-opal
  converter handles one season only.
* **Catch timing.** opal removes each season's catch at its start (Pope);
  SS3's hybrid F removes it through the season.
* **Time-varying selectivity** (blocks): reproduced exactly by splitting the
  fleets, including every fleet that mirrors a blocked fleet.
* **Recruitment index** (SS3 units 33): represented by an index in numbers
  with selectivity 1 at age 0. `sel_fa_external` replaces the selectivity of
  every fishery, so in a fit with opal's own curves the age-0 selectivity is
  an age-based double normal pinned to that shape.
* **Initial equilibrium.** opal's initial equilibrium removes the initial-F
  catch at the start of the year; SS3's applies the initial F through the
  year.

## Order of running

```
scripts/01-build-ss3-bridge.r              # ~1 h with Hessians (NPO_SWO_NOHESS=true skips them)
scripts/02-convert-fit-npo-swo-opal.r      # fits cached; NPO_SWO_REFIT=true refits
quarto render vignettes/npo-swo-ss3-to-opal.qmd
```

## Credit

Model files: ISC Billfish Working Group (2023), *Stock assessment report for
swordfish (Xiphias gladius) in the North Pacific through 2021*,
ISC/23/ANNEX/11.
