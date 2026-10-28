# opal-sandbox

[![tests](https://github.com/N-DucharmeBarth/opal-sandbox/actions/workflows/tests.yaml/badge.svg)](https://github.com/N-DucharmeBarth/opal-sandbox/actions/workflows/tests.yaml)
[![codecov](https://codecov.io/gh/N-DucharmeBarth/opal-sandbox/graph/badge.svg)](https://codecov.io/gh/N-DucharmeBarth/opal-sandbox)

Case studies and tests for [opal](https://github.com/N-DucharmeBarth-NOAA/opal),
an RTMB-based, age-structured stock assessment model.

Each case study starts from an existing assessment, usually a Stock Synthesis
(SS3) or MULTIFAN-CL (MFCL) model. It is converted into an opal model and the
conversion is checked, then opal is fitted and compared with the original. The
work is written up as reproducible vignettes. The case studies are also used to
find where opal differs from SS3 and MFCL, test new opal features, and produce
the example data and articles that go into the opal package.

Any model runs contained within this repository represent preliminary model explorations and should not be used '*as is*' for the consideration of management advice. Only models that have been presented and accepted by the relevant Regional Fisheries Management Organisation's scientific body should be used as the basis for management advice.

## Repository layout

```
case-studies/          one folder per stock; nothing in it is used by another case study
  wcpo-bet/            WCPO bigeye tuna: MFCL single-region model rebuilt in SS3
  io-yft/              Indian Ocean yellowfin (simulated data): SS3 -> opal
shared/                code used by more than one case study
  ss3-to-opal/         SS3-to-opal converter, validation and fit helpers
  extractors/ss3/      SS3 output extractors and plotting helpers
  extractors/mfcl/     MFCL output extractors and .frq parser
  setup/               install script for R, opal, r4ss and Quarto in a cloud container
simulations/           simulation studies of opal behaviour (selectivity, catch timing)
tests/testthat/        unit tests for the shared extractors and converter
executables/           SS3 (3.30.22.1-3.30.24.1) and MFCL binaries
```

Each case study is laid out the same way, with only the folders it needs:

| Folder | Contents |
|---|---|
| `model-files/` | the source assessment models (`ss3/`, `mfcl/`) and any reference data |
| `config/` | converter settings |
| `scripts/` | numbered R scripts that convert, fit and write caches |
| `vignettes/` | Quarto write-ups (rendered PDF/HTML alongside) that read only the caches in `vignettes/data/` |
| `opal-export/`, `opal-patches/` | material prepared for the opal repository |
| `README.md` | background, contents and the order to run things |

## Case studies

| Case study | Source model | Status |
|---|---|---|
| [`wcpo-bet`](case-studies/wcpo-bet) | 2020 WCPFC bigeye tuna, simplified to one region in MFCL, then rebuilt in SS3 (from a southwest Pacific striped marlin base file) | SS3 bridge complete; data prepared for RTMB; no opal model yet |
| [`io-yft`](case-studies/io-yft) | IATTC team's SS3 model of simulated Indian Ocean yellowfin (Spatial Stock Assessment Simulation Experiment) | Converted to opal, tested and fitted; source of opal's `io_yft_*` data and articles |

### Adding a case study

1. Create `case-studies/<stock>/` with a `README.md`, and put the source model under `model-files/`.
2. Use the converter in `shared/ss3-to-opal/` and follow `case-studies/io-yft/vignettes/ss3-to-opal-workflow.qmd`, the step-by-step guide. Write the settings to `config/` and the scripts to `scripts/`.
3. Code that a second case study needs goes in `shared/`. Code specific to one stock stays in its own folder.
4. `.gitignore` ignores everything that is not explicitly whitelisted, so add any new file types there.

## Shared code

### SS3-to-opal converter (`shared/ss3-to-opal/`)

These functions read an SS3 model with r4ss and build opal's data, parameter and map lists from it. They also check the result against SS3 (`ss3-validate.R`), build and fit the opal object (`ss3-opal-obj.R`, `ss3-fit-tools.R`), and handle SS3 generalised size-frequency data (`ss3-sizefreq.R`). The guide `case-studies/io-yft/vignettes/ss3-to-opal-workflow.qmd` explains how to use them.

### Extractors (`shared/extractors/`)

The SS3 and MFCL extractors return outputs in a common format, so that models can be compared directly:

- `extract_ss3_biomass()` and `extract_mfcl_biomass()` return `model, year, ts, season, ssb, ssb_se, depletion, depletion_se`. With `quarterly = FALSE`, quarterly values are aggregated to annual (season = 1).
- `extract_ss3_length_comp()` / `extract_mfcl_length_comp()` and the matching weight-composition extractors return `id, Fleet, Fleet_name, Used, Kind, Sex, Bin, Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj`. They can harmonise bins with `rebin_composition()`.
- There are also CPUE and selectivity extractors and modular plotting functions.

## Running the code

R is managed with `renv` (`renv.lock`). Scripts find the project root with `this.path::this.proj()`, so they can be run from any working directory. In a cloud container with no R installed, `shared/setup/install-r-opal.sh` installs R 4.5, opal, r4ss and Quarto (see `CLAUDE.md`).

## Unit testing

The tests in `tests/testthat/` cover the shared code:

- The extractors are checked for the structure, types and value ranges of their output, and for compatibility between the SS3 and MFCL versions. These tests use the WCPO bigeye model files.
- The SS3-to-opal converter is checked on the IO yellowfin case study (`test-ss3-to-opal-io-yft.r`). The test converts the SS3 model and runs the structural checks and the reconciliation with SS3. It then builds the opal object and checks that replaying the converted inputs with SS3's catch timing reproduces SS3's spawning biomass to within 1%.

The SS3 length- and weight-composition tests read `CompReport.sso`, which is not in the repository. Those two test files are skipped until `case-studies/wcpo-bet/model-files/ss3/01-bet-base/CompReport.sso` is added.

The GitHub Actions workflow `.github/workflows/tests.yaml` runs the tests on every pull request and on every push to `main`, and caches the R package library between runs. It then measures the tests' coverage of `shared/` with `tests/coverage.r` and uploads it to [Codecov](https://codecov.io/gh/N-DucharmeBarth/opal-sandbox). Changes are listed in [`NEWS.md`](NEWS.md). To run them locally from the project root:

```r
testthat::test_dir("tests/testthat")
# or
source("tests/run-tests.r")
```

## License

The code contained in this repository is licensed under the GNU GENERAL PUBLIC LICENSE version 3 ([GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)).

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
