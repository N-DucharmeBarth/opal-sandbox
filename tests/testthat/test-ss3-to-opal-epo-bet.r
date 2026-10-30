# Tests for the SS3-to-opal converter (shared/ss3-to-opal) on the EPO bigeye
# tuna case study: SS3 3.30, quarterly (each SS3 "year" is a quarter), ages
# 0-40 quarters, growth cessation (GrowthModel 8) with SD-linear growth
# variability (CV_Growth_Pattern 2), size double normals with pattern-15
# mirrors, an age-based discard fleet, catch in weight and numbers, early and
# late recruitment deviations, and length and generalised size compositions as
# two bin schemes. Mirrors case-studies/epo-bet/scripts/02-convert-fit-epo-bet-opal.r
# without fitting or writing any files.

library(testthat)
library(data.table)
library(magrittr)

skip_if_not_installed("opal")
skip_if_not_installed("r4ss")
skip_if_not_installed("yaml")
skip_if_not("sd_model" %in% names(formals(opal::get_sd_at_age)),
            "needs opal >= 0.1.2 (SD-linear growth variability)")
library(opal)

# Define project paths
proj_dir = this.path::this.proj()
dir_helper_fns_rtmb = file.path(proj_dir, "shared", "ss3-to-opal")
dir_config = file.path(proj_dir, "case-studies", "epo-bet", "config")

# Source helper functions
invisible(sapply(list.files(dir_helper_fns_rtmb, pattern = "\\.[Rr]$", full.names = TRUE), source))

context("SS3-to-opal converter: EPO bigeye")

# Convert once and reuse across tests
cfg = ss3_read_config(file.path(dir_config, "epo-bet-b1.yaml"))
converted = suppressWarnings(ss3_to_opal(cfg, dir_root = proj_dir))
d = converted$data
rep_ss3 = converted$ss3$rep
fish = converted$meta$fisheries
gf = converted$meta$growth_fit

lf_wide = dcast(converted$lf, fishery + year + month + ts ~ bin,
                value.var = "value", fill = 0) %>%
          setorder(fishery, ts) %>%
          as.data.frame()
lf_wide = lf_wide[lf_wide$fishery %in% as.integer(cfg$lf$keep_fisheries), ]
sf = ss3_sizefreq_to_opal(converted$ss3, cfg, converted$meta$dims, fish)
sf_fisheries = as.integer(cfg$sizefreq$keep_fisheries)
var_adjust = converted$meta$lf_var_adjust
var_adjust[sf_fisheries] = sf$var_adjust[sf_fisheries]
top = d$len_bin_start + d$n_len * d$len_bin_width
converted$data = prep_lf_data(
    d, list(lf_wide, sf$wide),
    lf_switch         = 1L,
    lf_keep_fisheries = sort(c(as.integer(cfg$lf$keep_fisheries), sf_fisheries)),
    lf_var_adjust     = var_adjust,
    lf_addtocomp      = as.numeric(cfg$lf$addtocomp),
    lf_upper_edge     = c(top, sf$upper_edge)
)
converted$data$lambdas = list(prior = 0)
converted$data$sel_basis_f = replace(rep(1L, d$n_fishery), as.integer(cfg$selectivity$age_based), 2L)
converted$data = ss3_make_portable(converted$data)

test_that("the model is quarterly with ages 0-40", {
    expect_equal(d$n_year, 188L)
    expect_equal(d$n_season, 1L)
    expect_equal(d$n_age, 41L)
    expect_equal(d$min_age, 0L)
})

test_that("the fleet layout follows SS3's fleet types", {
    # 22 fishing fleets and the longline survey
    expect_equal(fish$extract_fleets, 1:22)
    expect_equal(fish$index_map$ss3_fleet, 23L)
    expect_equal(d$n_fishery, 23L)
    # LL-n fleets and the survey in numbers, the rest in weight
    expect_equal(d$catch_units_f, c(rep(2, 7), rep(1, 15), 2))
})

test_that("converted inputs pass every structural check", {
    checks = validate_opal_inputs(converted, strict = FALSE)
    expect_true(all(checks$pass), info = paste(checks[pass == FALSE, check], collapse = "; "))
})

test_that("converted inputs reconcile with SS3", {
    recon = ss3_reconcile(converted, tol = 1e-6)
    expect_true(all(recon$pass, na.rm = TRUE), info = paste(recon[pass %in% FALSE, quantity], collapse = "; "))
    expect_equal(sum(grepl("^catch total", recon$quantity)), 22L)
    expect_true(any(grepl("growth cessation", recon$quantity)))
})

test_that("SS3's growth cessation curve is carried exactly", {
    expect_equal(gf$source, "ss3_cessation")
    expect_equal(d$growth_model, "cessation")
    # SS3 evaluates the curve at real age minus AFIX (0); opal at positions 1-41
    expect_equal(d$growth_age_offset, -1L)
    expect_equal(c(d$A1, d$A2), c(1L, 41L))
    # t50 as SS3 computes it from Linf
    expect_equal(gf$t50, log(exp((gf$Linf - gf$L0) * gf$k_gc / gf$rmax) - 1) / gf$k_gc)
    # opal's objective function reproduces SS3's length-at-age below the plus group
    p = converted$parameters
    mu = get_growth(d$n_age, d$A1, d$A2, growth_model = "cessation",
                    L0 = exp(p$log_L0), log_rmax = p$log_rmax, log_k_gc = p$log_k_gc,
                    t50 = p$t50_gc, age_offset = d$growth_age_offset)
    eg = ss3_endgrowth_table(rep_ss3, converted$meta$dims, "Beg")
    expect_lt(max(abs(mu[-41] / eg$length[-41] - 1)), 5e-5)
    # SS3 averages the plus group's length over the ages it accumulates
    expect_lt(abs(mu[41] - eg$length[41]), 0.2)
})

test_that("SS3's SD-linear growth variability is carried", {
    expect_equal(d$growth_sd_model, "sd")
    p = converted$parameters
    mu = gf$laa_opal
    sd_a = get_sd_at_age(mu, mu[d$A1], mu[d$A2], sd_model = "sd",
                         log_SD1 = p$log_SD1, log_SD2 = p$log_SD2, A1 = d$A1, A2 = d$A2)
    eg = ss3_endgrowth_table(rep_ss3, converted$meta$dims, "Beg")
    expect_lt(max(abs(sd_a[-41] / eg$sd[-41] - 1)), 5e-5)
    # SS3 holds the oldest age at SD_old itself
    expect_lt(abs(sd_a[41] / eg$sd[41] - 1), 0.005)
})

test_that("the growth reference ages are passed to opal as positions", {
    # opal_model() evaluates growth at positions 1, ..., n_age, so with
    # min_age = 0 the first position is SS3 age 0
    expect_equal(ss3_age_position(c(0L, 40L), converted$meta$dims), c(1L, 41L))
})

test_that("recruitment deviations are aligned to SS3's recruitment quarters", {
    # age 0 kept: opal's rdev_y[y] is SS3's deviation for quarter 17 + y;
    # main deviations to 202, late deviations 203-204
    dev = ss3_recdevs(converted$ss3$par, c(converted$ss3$ctl, list(endyr = converted$ss3$dat$endyr)))
    p = converted$parameters
    expect_equal(p$rdev_y[1:187], dev$value[match(18:204, dev$year)])
    expect_equal(p$rdev_y[188], 0)
    expect_true(all(c(203L, 204L) %in% dev$year))
    expect_false(205L %in% dev$year)
    # initial age structure: quarter 17 (age 0) and the early deviations -7 to 16
    expect_equal(p$init_rdev_a[1:25], dev$value[match(17 - 0:24, dev$year)])
    expect_equal(exp(p$log_init_F_f[c(2, 19)]), converted$ss3$par$init_F, tolerance = 1e-8)
})

test_that("generalised size compositions become a second bin scheme", {
    expect_equal(sf$bins, seq(60, 190, by = 10))
    expect_equal(sort(unique(sf$wide$fishery)), sf_fisheries)
    # rows sum to SS3's input sample size
    sfo = as.data.table(converted$ss3$dat$sizefreq_data_list[[1]])[fleet %in% sf_fisheries]
    expect_equal(sum(as.matrix(sf$wide[, -(1:4)])), sum(sfo$Nsamp), tolerance = 1e-8)
    # SS3's mult_by_generalized_sizecomp, inverted
    expect_equal(sf$var_adjust[c(2, 3, 23)], 1 / c(0.8, 0.59316, 0.8))
    expect_equal(length(converted$data$lf_scheme_edges), 2L)
})

assessment = ss3_opal_obj(converted, cfg, control = list(eval.max = 10000, iter.max = 10000),
                          build = FALSE)

test_that("the growth cessation parameters are fixed at SS3's values", {
    for (nm in c("log_L0", "log_rmax", "log_k_gc", "t50_gc", "log_SD1", "log_SD2")) {
        expect_true(all(is.na(assessment$map[[nm]])), info = nm)
    }
})

test_that("opal's dynamics at SS3's values equal the replayed conversion", {
    x = assessment
    p = x$parameters
    q = as.data.table(rep_ss3$parameters)[grepl("^LnQ_base", Label)]$Value
    p$log_cpue_q = rep(0, length(q))
    x = ss3_use_external_selectivity(opal_update(x, parameters = p),
                                     converted$selectivity$sel_fa_external) %>% opal_build()
    r = opal_report(x)
    replay = ss3_bridge_replica(converted, catch = "start", sr = "opal")
    expect_lt(max(abs(r$spawning_biomass_y[seq_len(d$n_year)] / replay - 1)), 1e-6)
    # the unfished state is SS3's
    expect_equal(exp(converted$parameters$log_B0), ss3_derived_value(rep_ss3, "SSB_Virgin"))
})

test_that("seeded double normals reproduce SS3's selectivity-at-length", {
    seed = ss3_seed_par_sel_size_dblnorm(converted$ss3, fish, converted$data$len_mid,
                                         converted$parameters$par_sel)
    ss = as.data.table(rep_ss3$sizeselex)[Factor == "Lsel" & Yr == 204]
    len_cols = setdiff(names(ss), c("Factor", "Fleet", "Yr", "Sex", "Label"))
    len_cols = len_cols[match(converted$data$len_mid, as.numeric(len_cols))]
    for (f in setdiff(seq_len(23), 20L)) {
        target = as.numeric(ss[Fleet == seed$source_fleet[f], ..len_cols])
        expect_lt(max(abs(sel_double_normal(converted$data$len_mid, seed$par_sel[f, ]) - target)),
                  1e-3, label = sprintf("fishery %d", f))
    }
})
