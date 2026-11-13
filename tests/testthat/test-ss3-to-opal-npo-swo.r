# Tests for the SS3-to-opal converter (shared/ss3-to-opal) on the NPO
# swordfish case study: SS3 3.30 bridge model B7 (annual, single sex), with
# fleets split at selectivity blocks and numbered after the surveys, size-based
# double-normal and logistic selectivity with mirrors, catch in weight and
# numbers, an initial F, early recruitment deviations and a recruitment index
# (SS3 units 33). Mirrors case-studies/npo-swo/scripts/02-convert-fit-npo-swo-opal.r
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
dir_config = file.path(proj_dir, "case-studies", "npo-swo", "config")

# Source helper functions
invisible(sapply(list.files(dir_helper_fns_rtmb, pattern = "\\.[Rr]$", full.names = TRUE), source))

context("SS3-to-opal converter: NPO swordfish")

# Convert once and reuse across tests
cfg = ss3_read_config(file.path(dir_config, "npo-swo-b7.yaml"))
converted = suppressWarnings(ss3_to_opal(cfg, dir_root = proj_dir))
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
yrs = 1975:2021

test_that("the fleet layout follows SS3's fleet types", {
    # 25 fishing fleets (F1-F19 and the six block fleets numbered after the
    # surveys) and 7 indices
    expect_equal(fish$extract_fleets, c(1:19, 28L, 29L, 31:34))
    expect_equal(fish$index_map$ss3_fleet, c(20:23, 25L, 26L, 30L))
    expect_equal(converted$data$n_fishery, 32L)
    expect_equal(converted$data$n_year, 47L)
    expect_equal(converted$data$n_age, 16L)
    # catch in numbers: Japanese and US longline fleets (and F1 1994-97)
    expect_equal(which(converted$data$catch_units_f[1:25] == 2), c(1, 3, 5, 6, 8, 9, 17, 18, 20))
})

test_that("converted inputs pass every structural check", {
    checks = validate_opal_inputs(converted, strict = FALSE)
    expect_true(all(checks$pass), info = paste(checks[pass == FALSE, check], collapse = "; "))
})

test_that("converted inputs reconcile with SS3", {
    recon = ss3_reconcile(converted, tol = 1e-6)
    expect_true(all(recon$pass, na.rm = TRUE), info = paste(recon[pass %in% FALSE, quantity], collapse = "; "))
    expect_equal(sum(grepl("^catch total", recon$quantity)), 25L)
})

test_that("the converter reads the recruitment index's units as biomass", {
    # SS3 units 33 (recruitment) has no opal equivalent; the fitting script
    # sets it to numbers with age-0 selectivity
    f_s6 = fish$index_map[ss3_fleet == 25L, opal_fleet]
    expect_equal(converted$data$catch_units_f[f_s6], 1)
    expect_true(all(converted$selectivity$sel_fa_external[f_s6, ] == 1))
})

test_that("indices and compositions are SS3's", {
    cp = as.data.table(converted$data$cpue_data)
    ss3 = as.data.table(converted$ss3$dat$CPUE)[index %in% cfg$index$keep_ss3_fleets]
    setorder(ss3, index, year)
    expect_equal(cp$value, ss3$obs)
    expect_equal(cp$se, ss3$se_log)
    # one composition per fleet and year, on SS3's 5 cm bins
    lc = as.data.table(converted$lf)
    expect_false(anyDuplicated(unique(lc[, .(fishery, ts)])) > 0)
    expect_equal(converted$data$len_bin_start, 10)
    expect_equal(converted$data$len_bin_width, 5)
})

test_that("recruitment deviations are aligned to SS3's recruitment years", {
    dev = ss3_recdevs(converted$ss3$par, converted$ss3$ctl)
    p = converted$parameters
    expect_equal(p$rdev_y[1:46], dev$value[match(1976:2021, dev$year)])
    expect_equal(p$init_rdev_a[1:11], dev$value[match(1975 - 0:10, dev$year)])
    expect_equal(exp(p$log_init_F_f[6]), as.numeric(converted$ss3$par$init_F), tolerance = 1e-8)
})

assessment = ss3_opal_obj(converted, cfg, map = NULL, bounds = NULL,
                          control = list(eval.max = 10000, iter.max = 10000),
                          build = FALSE)
bridge_obj = ss3_use_external_selectivity(assessment, converted$selectivity$sel_fa_external) %>%
    opal_build()
rep_opal = opal_report(bridge_obj)

test_that("opal's unfished numbers-at-age equal SS3's", {
    virg = as.data.table(rep_ss3$natage)[Era == "VIRG" & `Beg/Mid` == "B"]
    n0_ss3 = as.numeric(virg[1, as.character(0:15), with = FALSE])
    init = get_initial_numbers(exp(converted$parameters$log_B0), exp(converted$parameters$log_h),
                               converted$data$M, converted$data$spawning_potential)
    expect_equal(init$Ninit0, n0_ss3, tolerance = 1e-5)
})

test_that("opal at SS3's values tracks SS3's spawning biomass", {
    ts = as.data.table(rep_ss3$timeseries)[Era == "TIME" & Yr %in% yrs]
    ssb_rel = rep_opal$spawning_biomass_y[seq_along(yrs)] / ts$SpawnBio - 1
    # opal's initial equilibrium removes F6's initial catch at the start of
    # the year (SS3: through the year), and later years differ through catch
    # timing. Thresholds are set from the bridge in the vignette.
    expect_lt(abs(ssb_rel[1]), 0.06)
    expect_lt(max(abs(ssb_rel)), 0.30)
})

test_that("SS3's dynamics replayed with Baranov catch and SS3's catch weights reproduce SS3", {
    natage = as.data.table(rep_ss3$natage)[Era == "TIME" & `Beg/Mid` == "B" & Yr %in% yrs]
    ages = as.character(0:15)
    sel = converted$selectivity$sel_fa_external
    bodywt = as.data.table(rep_ss3$ageselex)[Factor == "bodywt" & Yr == 2021]
    ss3_fleet_f = c(fish$extract_fleets, fish$index_map$ss3_fleet)
    W = as.matrix(bodywt[match(ss3_fleet_f, Fleet), ..ages])
    d = converted$data
    num = d$catch_units_f != 1
    N = as.numeric(natage[Yr == yrs[1], ..ages])
    ssb = numeric(length(yrs))
    for (y in seq_along(yrs)) {
        N[1] = natage[Yr == yrs[y]][["0"]]
        ssb[y] = sum(N * d$spawning_potential)
        cy = d$catch_obs_ysf[y, 1, ]
        fl = which(cy > 0)
        cw = function(f) if (num[f]) 1 else W[f, ]
        F_f = numeric(length(cy))
        for (f in fl) F_f[f] = cy[f] / sum(N * sel[f, ] * cw(f))
        for (it in 1:50) {
            Z = d$M + colSums(F_f * sel)
            for (f in fl) F_f[f] = F_f[f] * cy[f] / sum(F_f[f] * sel[f, ] / Z * N * (1 - exp(-Z)) * cw(f))
        }
        N_end = N * exp(-(d$M + colSums(F_f * sel)))
        N = c(0, N_end[-16])
        N[16] = N[16] + N_end[16]
    }
    ts = as.data.table(rep_ss3$timeseries)[Era == "TIME" & Yr %in% yrs]
    expect_lt(max(abs(ssb / ts$SpawnBio - 1)), 1e-4)
})

test_that("seeded double normals reproduce SS3's selectivity-at-length", {
    seed = ss3_seed_par_sel_size_dblnorm(converted$ss3, fish, converted$data$len_mid,
                                         converted$parameters$par_sel)
    ss = as.data.table(rep_ss3$sizeselex)[Factor == "Lsel" & Yr == 2021]
    len_cols = setdiff(names(ss), c("Factor", "Fleet", "Yr", "Sex", "Label"))
    ss3_fleet_f = c(fish$extract_fleets, fish$index_map$ss3_fleet)
    dn = which(!is.na(seed$source_fleet))
    expect_true(all(c(1, 3, 6, 7, 8, 20) %in% dn))
    for (i in dn) {
        target = as.numeric(ss[Fleet == ss3_fleet_f[i], ..len_cols])
        expect_lt(max(abs(sel_double_normal(converted$data$len_mid, seed$par_sel[i, ]) - target)),
                  1e-3, label = sprintf("fishery %d", i))
    }
})

test_that("the opal object builds with a finite objective", {
    expect_s3_class(bridge_obj, "opal_obj")
    expect_true(is.finite(bridge_obj$build$objective))
})
