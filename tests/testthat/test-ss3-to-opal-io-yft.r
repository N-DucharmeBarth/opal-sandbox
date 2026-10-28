# Tests for the SS3-to-opal converter (shared/ss3-to-opal) on the IO
# yellowfin case study. Mirrors case-studies/io-yft/scripts/00-prep-io-yft-sim-data-for-rtmb.r
# without writing any files.

library(testthat)
library(data.table)
library(magrittr)

skip_if_not_installed("opal")
skip_if_not_installed("r4ss")
skip_if_not_installed("yaml")
library(opal)

# Define project paths
proj_dir = this.path::this.proj()
dir_helper_fns_rtmb = file.path(proj_dir, "shared", "ss3-to-opal")
dir_config = file.path(proj_dir, "case-studies", "io-yft", "config")

# Source helper functions
invisible(sapply(list.files(dir_helper_fns_rtmb, pattern = "\\.[Rr]$", full.names = TRUE), source))

context("SS3-to-opal converter: IO yellowfin")

# Convert once and reuse across tests
cfg = ss3_read_config(file.path(dir_config, "io-yft-sim-iattc-sample-tree-4.yaml"))
cfg$bio$source = "endgrowth"
converted = ss3_to_opal(cfg, dir_root = proj_dir)
lf_wide = dcast(converted$lf, fishery + year + month + ts ~ bin,
                value.var = "value", fill = 0) %>%
          setorder(fishery, ts) %>%
          as.data.frame()
converted$data = prep_lf_data(
    converted$data,
    lf_wide,
    lf_switch         = 1L,
    lf_keep_fisheries = as.integer(cfg$lf$keep_fisheries),
    lf_var_adjust     = converted$meta$lf_var_adjust,
    lf_cap            = cfg$lf$cap,
    lf_addtocomp      = as.numeric(cfg$lf$addtocomp)
) %>% ss3_make_portable()

test_that("converted inputs pass every structural check", {
    checks = validate_opal_inputs(converted, strict = FALSE)
    expect_gt(nrow(checks), 0)
    expect_true(all(checks$pass), info = paste(checks[pass == FALSE, check], collapse = "; "))
})

test_that("converted inputs reconcile with SS3", {
    recon = ss3_reconcile(converted, tol = 1e-6)
    expect_gt(nrow(recon), 0)
    # rows with pass = NA are informational notes (growth approximations)
    expect_true(all(recon$pass, na.rm = TRUE), info = paste(recon[pass %in% FALSE, quantity], collapse = "; "))
    expect_gt(sum(!is.na(recon$pass)), 0)
})

assessment = ss3_opal_obj(converted, cfg, map = NULL, bounds = NULL,
                          control = list(eval.max = 10000, iter.max = 10000),
                          build = TRUE)

# opal_build() stops if the initial objective or gradient is not finite
test_that("the opal object builds with a finite objective", {
    expect_s3_class(assessment, "opal_obj")
    expect_true(is.finite(assessment$build$objective))
})

test_that("replaying the converted inputs with SS3 timing reproduces SS3 spawning biomass", {
    expect_false(is.null(converted$selectivity$sel_fa_external))
    external = ss3_use_external_selectivity(assessment, converted$selectivity$sel_fa_external) %>%
        opal_build()
    bridge = ss3_bridge_table(converted, opal_report(external))
    expect_lt(max(abs(bridge$rel_replica_ss3)), 0.01)
    # opal's own dynamics differ through catch timing (start-of-season harvest),
    # documented in case-studies/io-yft/vignettes; guard against gross breakage only
    expect_lt(max(abs(bridge$rel_opal)), 0.25)
})
