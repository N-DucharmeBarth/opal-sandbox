# Indian Ocean yellowfin example for the opal package

Everything in this folder is laid out as it should be in the
[opal](https://github.com/N-DucharmeBarth-NOAA/opal) repository. It adds a
self-contained article, "Comparing opal with Stock Synthesis: Indian Ocean
yellowfin", which compares opal fit C2 with the SS3 model and with the
simulated truth.

| File | Size | What it is |
|---|---|---|
| `data/io_yft_data.rda` | 135 kB | opal data list (catch, index, length compositions, biology, priors) |
| `data/io_yft_parameters.rda` | 3 kB | starting parameters at the SS3 solution, with seeded selectivity |
| `data/io_yft_map.rda` | <1 kB | fit C2's map (310 estimated parameters) |
| `data/io_yft_reference.rda` | 20 kB | SS3 and simulated-truth comparison tables |
| `inst/extdata/io_yft_c2_fit.rds` | 143 kB | the fitted C2 `opal_obj` (`opal_save()`) |
| `R/data-io-yft.R` | | roxygen documentation for the four datasets |
| `man/io_yft_*.Rd` | | the Rd files roxygen generates from it |
| `data-raw/generate-io-yft-c2-fit.R` | | refits C2 from the bundled data inside opal (~2-3 min) |
| `vignettes/io-yft-ss3-comparison.qmd` | | the article |

The vignette uses only exported opal functions, `ggplot2` (an Import) and
`knitr` (Suggested). It loads the stored fit rather than refitting, and
renders in about 50 s, most of it rebuilding the RTMB tape once for
`opal_report()`.

## Adding it to opal (dev branch)

1. Copy the folder's contents into the opal repository root:

   ```bash
   cp -r data inst R man data-raw vignettes /path/to/opal/
   ```

   (Leave this README behind.)

2. Regenerate the documentation with opal's roxygen version, which also
   refreshes the copied `man/` files:

   ```r
   devtools::document()
   ```

3. Register the article and the datasets in `_pkgdown.yml`. pkgdown fails
   on documented topics missing from the reference index, so both are needed:

   ```yaml
   articles:
     - title: Assessment workflow
       navbar: ~
       contents:
         - quickstart
         - opakapaka-fit-review
         - assessment-tools
         - projections
     - title: Case studies
       navbar: ~
       contents:
         - io-yft-ss3-comparison
   ```

   and under `reference:` → `Bundled data`:

   ```yaml
         - io_yft_data
         - io_yft_map
         - io_yft_parameters
         - io_yft_reference
   ```

4. Check and preview:

   ```r
   devtools::check()                 # data, Rd and vignette build
   pkgdown::build_article("io-yft-ss3-comparison")
   ```

5. Optionally add a line to `NEWS.md`, e.g. "New article comparing opal with
   Stock Synthesis for simulated Indian Ocean yellowfin, with bundled data
   `io_yft_data`, `io_yft_parameters`, `io_yft_map` and `io_yft_reference`."

Nothing in `DESCRIPTION` needs to change.

## Regenerating

* **The data and the reference tables** come from this repository:
  `code/rtmb/scripts/04-export-io-yft-for-opal.r` (needs the SS3 files, r4ss
  and the caches written by scripts 01-03). It refits C2, stops unless the
  objective matches script 03's and `opal_check()` passes, and rewrites this
  folder.
* **The stored fit alone** can be regenerated inside opal from the bundled
  data with `data-raw/generate-io-yft-c2-fit.R`, for example after a change
  to opal's model code. Run it from the package root.

Both were last run with opal 0.1.1 (dev, 6ecfef4) and give the same
objective, 14406.585.
