# Case study: Indian Ocean yellowfin tuna (simulated data)

Converts a Stock Synthesis (SS3) model of simulated Indian Ocean yellowfin
into an [opal](https://github.com/N-DucharmeBarth-NOAA/opal) model, tests the
conversion, fits opal and compares it with SS3 and with the simulated truth.
It is the reference case study for the SS3-to-opal converter in
[`shared/ss3-to-opal/`](../../shared/ss3-to-opal), and the source of the
`io_yft_*` example data in the opal package.

The SS3 model was developed by the IATTC team (M.N. Maunder, H. Xu and
C. Minte-Vera) for the Spatial Stock Assessment Simulation Experiment; see
`model-files/ss3/xx-io-yft-sim-iattc-sample-tree-4/README.md` for credit and
citation.

## Contents

| Path | What it is |
|---|---|
| `model-files/ss3/xx-io-yft-sim-iattc-sample-tree-4/` | Baseline SS3 model (v3.24) |
| `model-files/ss3/xx-io-yft-sim-iattc-sample-tree-4-{lensel,simple-k}/` | SS3 variants (length-based selectivity; simple growth) |
| `model-files/sim-truth/io-yft-iattc-sample-tree-4/` | Operating-model truth for the simulated replicate |
| `config/` | Converter settings (YAML) |
| `scripts/00`-`03` | Convert, fit and extract the truth; write the caches in `vignettes/data/` |
| `scripts/04`, `06` | Export the fits and data as opal package data under `opal-export/` |
| `scripts/05-io-yft-ss3-to-opal-workflow.r` | End-to-end conversion, tests and fit; companion to `vignettes/ss3-to-opal-workflow.qmd` |
| `vignettes/` | The step-by-step workflow guide and the comparison reports (`.qmd` with rendered PDF/HTML); they read only the caches |
| `opal-export/` | Files laid out as in the opal repository (data, stored fits, docs, articles) plus `git am` patches; see its README |
| `opal-patches/io-yft-ss3-matched/` | Patch series and rendered article for the SS3-matched opal article |

The files under `opal-export/` (other than its README) and the patches are
snapshots of what was added to opal, so they still name the paths used before
this repository was reorganised (`code/rtmb/...`).

## Order of running

```
scripts/00-prep-io-yft-sim-data-for-rtmb.r   # convert, validate, build the opal object
scripts/01-fit-io-yft-ss3-opal.r             # fits -> vignettes/data/io-yft-ss3-opal-cache.rds
scripts/02-extract-io-yft-sim-truth.r        # truth -> vignettes/data/io-yft-sim-truth-cache.rds
scripts/03-fit-io-yft-age-k.r                # age-specific K fits -> vignettes/data/io-yft-age-k-cache.rds
scripts/05-io-yft-ss3-to-opal-workflow.r     # workflow guide cache
quarto render vignettes/<name>.qmd
```

The selectivity and catch-timing simulations in
[`simulations/`](../../simulations) use this case study's exported opal data
as their operating-model inputs.
