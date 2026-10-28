# Indian Ocean yellowfin example for the opal package

Everything in this folder is laid out as it should be in the
[opal](https://github.com/N-DucharmeBarth-NOAA/opal) repository. It adds two
self-contained articles on simulated Indian Ocean yellowfin:

* "Comparing opal with Stock Synthesis: Indian Ocean yellowfin": opal fit C2
  (length-based selectivity) against the SS3 model and the simulated truth.
* "Length- or age-based selectivity: Indian Ocean yellowfin": the same model
  with age-based selectivity (`sel_basis_f = 2`, opal's `age-selex` branch),
  started from SS3's own `age_selex` parameters. Fit A uses C2's map; fit
  A-SS3 estimates every selectivity parameter SS3 estimated.

| File | Size | What it is |
|---|---|---|
| `data/io_yft_data.rda` | 135 kB | opal data list (catch, index, length compositions, biology, priors) |
| `data/io_yft_parameters.rda` | 3 kB | starting parameters at the SS3 solution, with seeded selectivity |
| `data/io_yft_map.rda` | <1 kB | fit C2's map (310 estimated parameters) |
| `data/io_yft_reference.rda` | 21 kB | SS3 and simulated-truth comparison tables, including SS3's `age_selex` patterns and parameters (`ss3$sel_pattern`, `ss3$sel_par`) |
| `inst/extdata/io_yft_c2_fit.rds` | 143 kB | fitted C2 `opal_obj` |
| `inst/extdata/io_yft_age_fit.rds` | 143 kB | fitted A `opal_obj` (age-based, C2's map) |
| `inst/extdata/io_yft_age_ss3_fit.rds` | 143 kB | fitted A-SS3 `opal_obj` (age-based, SS3's estimated parameters) |
| `inst/extdata/io_yft_vignette_cache.rds` | 63 kB | everything the two articles plot, precomputed from the fits |
| `R/data-io-yft.R` | | roxygen documentation for the four datasets |
| `man/io_yft_*.Rd` | | the Rd files roxygen generates from it |
| `data-raw/generate-io-yft-c2-fit.R` | | refits C2 from the bundled data (~2-3 min) |
| `data-raw/generate-io-yft-age-fits.R` | | fits A and A-SS3 (~5 min); needs `sel_basis_f` |
| `data-raw/io-yft-age-setup.R` | | age-based inputs, maps and bounds shared by the two scripts |
| `data-raw/generate-io-yft-vignette-cache.R` | | rebuilds the cache from the stored fits (~1 min) |
| `vignettes/io-yft-ss3-comparison.qmd` | | the comparison article |
| `vignettes/io-yft-age-selectivity.qmd` | | the age-selectivity article |
| `patches/*.patch` | | the same changes as two commits on opal's `age-selex` branch |

The articles use only exported opal functions, `ggplot2` (an Import) and
`knitr` (Suggested). They read the cache instead of rebuilding the RTMB
model, so each renders in about 15 s; the fitting code is shown with
`eval: false`.

The fits and the cache were made on opal's `age-selex` branch (31360db). The
fits' model checksum matches that branch; a fit made with an older opal
loads with a checksum warning.

## Adding it to opal (age-selex branch)

Either apply the patches, which also edit `_pkgdown.yml` and `NEWS.md`:

```bash
cd /path/to/opal
git checkout age-selex
git am /path/to/opal-export/patches/*.patch
```

or copy the folder by hand:

1. Copy the folder's contents into the opal repository root (leave this
   README and `patches/` behind):

   ```bash
   cp -r data inst R man data-raw vignettes /path/to/opal/
   ```

2. Regenerate the documentation with opal's roxygen version:

   ```r
   devtools::document()
   ```

3. Register the articles and the datasets in `_pkgdown.yml`. pkgdown fails
   on documented topics missing from the reference index, so both are needed:

   ```yaml
   articles:
     ...
     - title: Case studies
       navbar: ~
       contents:
         - io-yft-ss3-comparison
         - io-yft-age-selectivity
   ```

   and under `reference:` → `Bundled data`:

   ```yaml
         - io_yft_data
         - io_yft_map
         - io_yft_parameters
         - io_yft_reference
   ```

4. Add a line to `NEWS.md` (the patches add one under opal 0.1.1).

5. Check and preview:

   ```r
   devtools::check()
   pkgdown::build_article("io-yft-ss3-comparison")
   pkgdown::build_article("io-yft-age-selectivity")
   ```

Nothing in `DESCRIPTION` needs to change.

## Regenerating

* **The data and the reference tables** come from this repository:
  `case-studies/io-yft/scripts/04-export-io-yft-for-opal.r` (needs the SS3 files, r4ss
  and the caches written by scripts 01-03). It refits C2, stops unless the
  objective matches script 03's and `opal_check()` passes, and rewrites
  `data/` and `inst/extdata/io_yft_c2_fit.rds` here. The C2 fit it writes
  carries the checksum of the opal it ran with, so refit it inside opal
  afterwards (next point) if that was not `age-selex`.
* **The stored fits** are regenerated inside opal, from the package root:
  `data-raw/generate-io-yft-c2-fit.R` and
  `data-raw/generate-io-yft-age-fits.R`.
* **The cache** must be rebuilt whenever the fits or opal's model code
  change: `data-raw/generate-io-yft-vignette-cache.R`. It records the opal
  version and commit, which the articles print.

Last run on opal `age-selex` (31360db): C2 objective 14406.585, A 14658.102,
A-SS3 14625.286, all passing `opal_check()`.
