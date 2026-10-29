# Tests for the NPO blue shark SS3 bridge models
# (case-studies/npo-bsh/model-files/ss3, built by
# case-studies/npo-bsh/scripts/01-build-ss3-bridge.r).
#
# Each step must have converged with an invertible Hessian and no parameter
# on a bound, must contain exactly the structural change it is meant to make,
# must keep the data unchanged, and must agree with the bridge summary cached
# for the vignette.

library(testthat)
library(data.table)

skip_if_not_installed("r4ss")

# Define project paths
proj_dir = this.path::this.proj()
dir_ss3 = file.path(proj_dir, "case-studies", "npo-bsh", "model-files", "ss3")
cache_file = file.path(proj_dir, "case-studies", "npo-bsh", "vignettes", "data",
                       "npo-bsh-ss3-bridge-cache.rds")

context("NPO blue shark SS3 bridge")

steps = c("00-s6-base", "01-single-sex", "02a-f18-dblnorm", "02b-f8-split",
          "02c-no-f19-devs", "02d-no-regime", "02e-no-f-ballpark",
          "03-start-of-year-obs")

# Report.sso header lines: final gradient and ln(det(Hessian))
report_header = function(step) {
    h = readLines(file.path(dir_ss3, step, "Report.sso"), n = 40L)
    list(
        version  = h[1],
        gradient = as.numeric(sub("^Convergence_Level: ([^ ]+).*", "\\1",
                                  grep("^Convergence_Level:", h, value = TRUE))),
        log_det  = suppressWarnings(as.numeric(sub("^Hessian: ([^ ]+).*", "\\1",
                                                   grep("^Hessian:", h, value = TRUE))))
    )
}
inputs = lapply(steps, function(s) suppressWarnings(
    r4ss::SS_read(file.path(dir_ss3, s), ss_new = TRUE, verbose = FALSE)))
names(inputs) = steps
reps = lapply(steps, function(s) suppressWarnings(
    r4ss::SS_output(file.path(dir_ss3, s), covar = FALSE, verbose = FALSE, printstats = FALSE)))
names(reps) = steps

test_that("every bridge step was run with SS3 3.30.22.1", {
    for (s in steps) expect_match(report_header(s)$version, "3.30.22.1", fixed = TRUE, info = s)
})

test_that("every bridge step converged with an invertible Hessian", {
    for (s in steps) {
        h = report_header(s)
        expect_lt(h$gradient, 1e-3)
        expect_true(is.finite(h$log_det), info = paste(s, "Hessian not inverted"))
    }
})

test_that("no estimated parameter is on a bound", {
    for (s in steps) {
        st = reps[[s]]$parameters$Status
        expect_false(any(grepl("LO|HI", st)), info = s)
    }
})

test_that("S6 reproduces the 2022 assessment's likelihood with SS3 3.30.22.1", {
    expect_equal(reps[["00-s6-base"]]$likelihoods_used["TOTAL", "values"], 32.2958, tolerance = 1e-4)
    dq = reps[["00-s6-base"]]$derived_quants
    expect_equal(dq$Value[dq$Label == "SSB_Virgin"], 243875, tolerance = 1e-4)
    expect_equal(dq$Value[dq$Label == "Bratio_2020"], 0.563918, tolerance = 1e-4)
})

test_that("B1 is single sex with female biology and no sex offsets", {
    dat = inputs[["01-single-sex"]]$dat
    ctl = inputs[["01-single-sex"]]$ctl
    expect_equal(dat$Nsexes, -1)
    expect_true(all(dat$lencomp$sex == 0))
    expect_false(any(grepl("^[fm][0-9]", names(dat$lencomp))))
    expect_equal(nrow(ctl$natM), 1L)
    expect_equal(as.numeric(ctl$natM[1, ]), as.numeric(inputs[["00-s6-base"]]$ctl$natM[1, ]))
    expect_false(any(grepl("_Mal_", rownames(ctl$MG_parms))))
    expect_true(all(ctl$size_selex_types$Male == 0))
    expect_false(any(grepl("PFemOff|PMalOff", rownames(ctl$size_selex_parms))))
    # female growth and weight-length are S6's
    mg0 = inputs[["00-s6-base"]]$ctl$MG_parms
    mg1 = ctl$MG_parms
    fem = grep("_Fem_GP_1", rownames(mg0), value = TRUE)
    expect_equal(mg1[fem, "INIT"][grepl("Wtlen|L_at|VonBert|Mat", fem)],
                 mg0[fem, "INIT"][grepl("Wtlen|L_at|VonBert|Mat", fem)])
})

test_that("B1 sums the sex x length compositions over sex", {
    lc0 = as.data.table(inputs[["00-s6-base"]]$dat$lencomp)
    lc1 = as.data.table(inputs[["01-single-sex"]]$dat$lencomp)
    f = grep("^f[0-9]", names(lc0), value = TRUE)
    # r4ss names single-sex bins l<bin>
    l = sub("^f", "l", f)
    expect_equal(sum(as.matrix(lc1[, ..l])),
                 sum(as.matrix(lc0[, ..f])) + sum(as.matrix(lc0[, sub("^f", "m", f), with = FALSE])))
    expect_equal(lc1$Nsamp, lc0$Nsamp)
})

test_that("B2a replaces F18's spline with a double normal", {
    st = inputs[["02a-f18-dblnorm"]]$ctl$size_selex_types
    expect_equal(st["F18_USA_Lonline_SH", "Pattern"], 24)
    expect_false(any(grepl("Spline|Spine", rownames(inputs[["02a-f18-dblnorm"]]$ctl$size_selex_parms))))
})

test_that("B2b splits F8 at its blocks without changing the data", {
    dat0 = inputs[["02a-f18-dblnorm"]]$dat
    dat  = inputs[["02b-f8-split"]]$dat
    ctl  = inputs[["02b-f8-split"]]$ctl
    expect_equal(dat$Nfleets, 32L)
    ct = dat$catch[dat$catch$year > 0 & dat$catch$catch > 0, ]
    expect_true(all(ct$year[ct$fleet == 8] <= 1981))
    expect_true(all(ct$year[ct$fleet == 32] >= 1982))
    expect_equal(sum(ct$catch[ct$fleet %in% c(8, 32)]),
                 sum(dat0$catch$catch[dat0$catch$fleet == 8 & dat0$catch$year > 0]))
    expect_equal(sum(dat$catch$catch[dat$catch$year > 0]), sum(dat0$catch$catch[dat0$catch$year > 0]))
    expect_equal(sum(dat$lencomp$Nsamp), sum(dat0$lencomp$Nsamp))
    expect_true(all(dat$lencomp$year[dat$lencomp$fleet == 32] >= 1982))
    i8 = grep("_F8_JPN_LG_MESH_EARLY\\(8\\)", rownames(ctl$size_selex_parms))
    expect_true(all(ctl$size_selex_parms[i8, "Block"] == 0))
    expect_equal(nrow(ctl$size_selex_parms[grep("\\(32\\)", rownames(ctl$size_selex_parms)), ]), 6L)
})

test_that("B2c removes F19's selectivity deviations", {
    sp = inputs[["02c-no-f19-devs"]]$ctl$size_selex_parms
    expect_true(all(sp$dev_link == 0))
    expect_null(inputs[["02c-no-f19-devs"]]$ctl$size_selex_parms_tv)
})

test_that("B2d removes the initial regime offset", {
    sr = inputs[["02d-no-regime"]]$ctl$SR_parms
    expect_equal(sr["SR_regime", "Block"], 0)
    expect_equal(sr["SR_regime", "INIT"], 0)
    expect_equal(reps[["02d-no-regime"]]$likelihoods_used["InitEQ_Regime", "values"], 0)
})

test_that("B2e switches the F-ballpark penalty off", {
    expect_lt(inputs[["02e-no-f-ballpark"]]$ctl$F_ballpark_year, 0)
    expect_false("F_Ballpark" %in% rownames(reps[["02e-no-f-ballpark"]]$likelihoods_used))
})

test_that("B3 moves every observation to month 1 and changes nothing else", {
    d2 = inputs[["02e-no-f-ballpark"]]$dat
    d3 = inputs[["03-start-of-year-obs"]]$dat
    expect_true(all(d3$CPUE$month == 1))
    expect_true(all(d3$lencomp$month == 1))
    expect_true(all(vapply(d3$sizefreq_data_list, function(sf) all(sf$month == 1), logical(1))))
    expect_true(all(d3$fleetinfo$surveytiming == 1))
    expect_equal(d2$fleetinfo$surveytiming[d2$fleetinfo$type == 1],
                 rep(-1, sum(d2$fleetinfo$type == 1)))
    # settlement (month 7) and the length-at-age are B2e's
    expect_equal(inputs[["03-start-of-year-obs"]]$ctl$recr_dist_pattern,
                 inputs[["02e-no-f-ballpark"]]$ctl$recr_dist_pattern)
    expect_equal(reps[["03-start-of-year-obs"]]$endgrowth$Len_Beg,
                 reps[["02e-no-f-ballpark"]]$endgrowth$Len_Beg)
    expect_equal(reps[["03-start-of-year-obs"]]$N_estimated_parameters,
                 reps[["02e-no-f-ballpark"]]$N_estimated_parameters)
})

test_that("every step keeps the catch and index data", {
    tot = function(s) {
        d = inputs[[s]]$dat
        c(catch = sum(d$catch$catch[d$catch$year > 0]), cpue = sum(d$CPUE$obs),
          n_cpue = nrow(d$CPUE), nsamp = sum(d$lencomp$Nsamp))
    }
    ref = tot("00-s6-base")
    for (s in steps[-1]) expect_equal(tot(s), ref, info = s)
})

test_that("the vignette's bridge cache matches the committed SS3 runs", {
    skip_if_not(file.exists(cache_file), "bridge cache not built")
    cache = readRDS(cache_file)$summary
    for (s in steps) {
        dq = reps[[s]]$derived_quants
        expect_equal(cache[step == s]$nll, reps[[s]]$likelihoods_used["TOTAL", "values"],
                     tolerance = 1e-6, info = s)
        expect_equal(cache[step == s]$depletion, dq$Value[dq$Label == "Bratio_2020"],
                     tolerance = 1e-6, info = s)
    }
})
