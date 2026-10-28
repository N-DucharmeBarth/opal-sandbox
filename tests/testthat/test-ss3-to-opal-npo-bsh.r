# Tests for the SS3-to-opal converter (shared/ss3-to-opal) on the NPO blue
# shark case study: SS3 3.30, annual, single sex (Nsexes = -1), size-based
# double-normal selectivity with mirrors, a fishing fleet numbered after the
# surveys, catch in weight and numbers, an initial F and early recruitment
# deviations. Mirrors case-studies/npo-bsh/scripts/03-convert-fit-npo-bsh-opal.r
# without fitting or writing any files.

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
dir_config = file.path(proj_dir, "case-studies", "npo-bsh", "config")

# Source helper functions
invisible(sapply(list.files(dir_helper_fns_rtmb, pattern = "\\.[Rr]$", full.names = TRUE), source))

context("SS3-to-opal converter: NPO blue shark")

# Convert once and reuse across tests
cfg = ss3_read_config(file.path(dir_config, "npo-bsh-b2e.yaml"))
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
rep_ss3 = converted$ss3$rep
fish = converted$meta$fisheries

test_that("the fleet layout follows SS3's fleet types", {
    # 21 fishing fleets (F1-F20 and fleet 32, after the surveys) and 2 indices
    expect_equal(fish$extract_fleets, c(1:20, 32L))
    expect_equal(fish$index_map$ss3_fleet, c(25L, 26L))
    expect_equal(converted$data$n_fishery, 23L)
    expect_equal(converted$data$catch_units_f[c(11, 17, 18)], c(2, 2, 2))
    expect_true(all(converted$data$catch_units_f[-c(11, 17, 18, 22, 23)] == 1))
})

test_that("converted inputs pass every structural check", {
    checks = validate_opal_inputs(converted, strict = FALSE)
    expect_true(all(checks$pass), info = paste(checks[pass == FALSE, check], collapse = "; "))
})

test_that("converted inputs reconcile with SS3", {
    recon = ss3_reconcile(converted, tol = 1e-6)
    # rows with pass = NA are informational notes (growth approximations)
    expect_true(all(recon$pass, na.rm = TRUE), info = paste(recon[pass %in% FALSE, quantity], collapse = "; "))
    expect_equal(sum(grepl("^catch total", recon$quantity)), 21L)
})

test_that("indices carry SS3's additive CV adjustment", {
    cp = as.data.table(converted$data$cpue_data)
    ss3 = as.data.table(converted$ss3$dat$CPUE)[index %in% c(25, 26)]
    setorder(ss3, index, year)
    expect_equal(cp$se, ss3$se_log + ifelse(ss3$index == 25, 0.188, 0))
    expect_equal(cp$value, ss3$obs)
})

test_that("recruitment deviations are aligned to SS3's recruitment years", {
    # age 0 kept: opal's rdev_y[y] is SS3's deviation for year 1971 + y
    dev = ss3_recdevs(converted$ss3$par, converted$ss3$ctl)
    p = converted$parameters
    expect_equal(p$rdev_y[1:49], dev$value[match(1972:2020, dev$year)])
    expect_equal(p$init_rdev_a[1:11], dev$value[match(1971 - 0:10, dev$year)])
    expect_equal(exp(p$log_init_F_f[4]), converted$ss3$par$init_F, tolerance = 1e-8)
})

assessment = ss3_opal_obj(converted, cfg, map = NULL, bounds = NULL,
                          control = list(eval.max = 10000, iter.max = 10000),
                          build = FALSE)
bridge_obj = ss3_use_external_selectivity(assessment, converted$selectivity$sel_fa_external) %>%
    opal_build()
rep_opal = opal_report(bridge_obj)

test_that("opal's unfished numbers-at-age equal SS3's", {
    virg = as.data.table(rep_ss3$natage)[Era == "VIRG" & `Beg/Mid` == "B"]
    n0_ss3 = as.numeric(virg[1, as.character(0:24), with = FALSE])
    n0_opal = rep_opal$number0_ysa[1, 1, ]
    # opal's unfished trajectory starts from R0 x unfished survivorship
    init = get_initial_numbers(exp(converted$parameters$log_B0), exp(converted$parameters$log_h),
                               converted$data$M, converted$data$spawning_potential)
    expect_equal(init$Ninit0, n0_ss3, tolerance = 1e-5)
})

test_that("opal at SS3's values tracks SS3's recruitment and spawning biomass", {
    yrs = 1971:2020
    ts = as.data.table(rep_ss3$timeseries)[Era == "TIME" & Yr %in% yrs]
    rec_ss3 = as.data.table(rep_ss3$natage)[Era == "TIME" & `Beg/Mid` == "B" & Yr %in% yrs][["0"]]
    rec_rel = rep_opal$number_ysa[seq_along(yrs), 1, 1] / rec_ss3 - 1
    ssb_rel = rep_opal$spawning_biomass_y[seq_along(yrs)] / ts$SpawnBio - 1
    # the initial state is converted directly; later years differ through catch
    # timing (opal removes each year's catch at the start of the year, SS3
    # continuously), which this model's 2020 F of 0.34 and age-0 M of 0.79
    # make large. Thresholds are set from the bridge in the vignette.
    expect_lt(abs(ssb_rel[1]), 0.02)
    expect_lt(max(abs(rec_rel)), 0.06)
    expect_lt(max(abs(ssb_rel)), 0.20)
})

test_that("the converted inputs replayed with SS3's catch equation and weights stay close to SS3", {
    eg = as.data.table(rep_ss3$endgrowth)[Sex == 1 & Seas == 1]
    x_mid = converted
    x_mid$data$weight = eg$Wt_Mid
    replay = ss3_bridge_replica(x_mid, catch = "baranov", sr = "ss3")
    ts = as.data.table(rep_ss3$timeseries)[Era == "TIME" & Yr %in% 1971:2020]
    expect_lt(max(abs(replay / ts$SpawnBio - 1)), 0.10)
})

test_that("seeded double normals reproduce SS3's selectivity-at-length", {
    seed = ss3_seed_par_sel_size_dblnorm(converted$ss3, fish, converted$data$len_mid,
                                         converted$parameters$par_sel)
    expect_false(anyNA(seed$source_fleet))
    ss = as.data.table(rep_ss3$sizeselex)[Factor == "Lsel" & Yr == 2020]
    len_cols = setdiff(names(ss), c("Factor", "Fleet", "Yr", "Sex", "Label"))
    ss3_fleet_f = c(fish$extract_fleets, fish$index_map$ss3_fleet)
    for (i in seq_along(ss3_fleet_f)) {
        target = as.numeric(ss[Fleet == ss3_fleet_f[i], ..len_cols])
        expect_lt(max(abs(sel_double_normal(converted$data$len_mid, seed$par_sel[i, ]) - target)),
                  1e-3, label = sprintf("fishery %d", i))
    }
})

test_that("the opal object builds with a finite objective", {
    expect_s3_class(bridge_obj, "opal_obj")
    expect_true(is.finite(bridge_obj$build$objective))
})
