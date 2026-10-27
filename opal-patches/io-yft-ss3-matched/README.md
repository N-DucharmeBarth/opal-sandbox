# opal patches: SS3-matched Indian Ocean yellowfin article

New opal article `vignettes/io-yft-ss3-matched.qmd` ("Matching the Stock
Synthesis configuration: Indian Ocean yellowfin"). It starts from the
age-based fit A-SS3 of the age-selectivity article (called **AGE** here; no
opal fit in the article has "SS3" in its name) and adds SS3's settings:

| Fit | Priors | Index size compositions (SS3 fleet 17, lambda 0.51) |
|---|---|---|
| AGE (= A-SS3) | opal defaults | no |
| AGE-NP | off (`lambdas$prior = 0`, as SS3) | no |
| AGE-NP-SF | off | yes, on SS3's 30-160 cm bins (second bin scheme) |
| AGE-SF | opal defaults | yes |

`io-yft-ss3-matched.html` here is the rendered article (self-contained).

## Base

`N-DucharmeBarth-NOAA/opal` `dev` at `b65b165` (Merge pull request #151).

## Commits

1. `0001` Add the SS3-matched Indian Ocean yellowfin article: the article,
   `data-raw/io-yft-ss3-matched-setup.R`,
   `data-raw/generate-io-yft-matched-fits.R`,
   `data-raw/generate-io-yft-matched-s3-fit.R`, the extended
   `data-raw/generate-io-yft-vignette-cache.R`, `io_yft_reference$ss3`
   gaining `sizefreq` and `lambdas` (with docs in `R/data-io-yft.R` and
   `man/io_yft_reference.Rd`, edited by hand: roxygen2 was not available),
   `_pkgdown.yml` and `NEWS.md`.
2. `0002` The stored fits and cache (binary, `--binary`):
   `inst/extdata/io_yft_matched_s{1,2,3}_fit.rds` (AGE-NP, AGE-NP-SF,
   AGE-SF), `io_yft_matched_s1_free.rds`, `io_yft_matched_starts.rds` and the
   rebuilt `io_yft_vignette_cache.rds`. The file and cache names keep the
   working labels s1/s2/s3.

## Apply

```bash
cd /path/to/opal
git fetch origin && git checkout dev && git pull
git am /path/to/opal-sandbox/opal-patches/io-yft-ss3-matched/*.patch
```

## Verified

* The series applies with `git am` to a clean checkout of `b65b165` and
  gives a tree identical to the one it was exported from.
* `tests/testthat/test-lambdas.R` (the tests that use the IO yellowfin data)
  pass: 0 failures.
* Rebuilding the cache leaves the existing entries (C2, A, A-SS3, A_start,
  E_start) unchanged: same objectives, B0 and trajectories.
* At SS3's solution, opal's predictions of the index size compositions are
  within 0.010 of SS3's expected values; the SS3-style likelihood is 1064.0
  (SS3's own expected values: 1067.5).
* AGE-SF passes `opal_check()`. AGE-NP and AGE-NP-SF do not (AGE-NP:
  gradient 3.2e-3 > 1e-4; AGE-NP-SF: Hessian not positive definite); the
  article says so.
* The article renders (quarto 1.9) with no chunk errors. The full package
  check (`R CMD check`) and a pkgdown build were not run.

## Not done / notes

* The commit hash the article prints (`cache$commit`) is the local one; it
  changes when the patches are applied. Rerun
  `data-raw/generate-io-yft-vignette-cache.R` to refresh it.
* The stored C2 and A fits load with a model-checksum warning (made on the
  `age-selex` branch); their cached values are unchanged.
* SS3's inputs come from `code/rtmb/scripts/06-export-io-yft-sizefreq-for-opal.r`
  in this repository.
