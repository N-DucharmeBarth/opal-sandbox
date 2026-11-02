# Tests for the EPO bigeye tuna SS3 bridge models
# (case-studies/epo-bet/model-files/ss3, built by
# case-studies/epo-bet/scripts/01-build-ss3-bridge.r).
#
# Each step must have converged with an invertible Hessian, must contain
# exactly the structural change it is meant to make, must keep the data
# unchanged, and must agree with the bridge summary cached for the vignette.
# The re-run of the supplied model must reproduce the supplied run.

library(testthat)
library(data.table)

skip_if_not_installed("r4ss")

# Define project paths
proj_dir = this.path::this.proj()
dir_ss3 = file.path(proj_dir, "case-studies", "epo-bet", "model-files", "ss3")
cache_file = file.path(proj_dir, "case-studies", "epo-bet", "vignettes", "data",
                       "epo-bet-ss3-bridge-cache.rds")

context("EPO bigeye SS3 bridge")

steps = c("00-base", "01-no-regime")

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
        # 01's final gradient is 0.02 (the model has 263 parameters, 211 of
        # them recruitment deviations); 00's is 0.0018
        expect_lt(h$gradient, 0.05)
        expect_true(is.finite(h$log_det), info = paste(s, "Hessian not inverted"))
    }
})

test_that("the re-run of the supplied model reproduces the supplied run", {
    # The supplied Report.sso came from an SS3 3.30.22 beta build
    supplied = fread(file.path(dir_ss3, "00-base", "supplied-likelihood-3.30.22beta.csv"))
    lk = reps[["00-base"]]$likelihoods_used
    for (comp in c("TOTAL", "Survey", "Length_comp", "SizeFreq", "Recruitment")) {
        expect_equal(lk[comp, "values"], supplied[component == comp, logL_lambda],
                     tolerance = 1e-5, info = comp)
    }
})

test_that("01 removes the initial-equilibrium regime offset and nothing else", {
    c0 = inputs[["00-base"]]$ctl
    c1 = inputs[["01-no-regime"]]$ctl
    expect_true(!is.null(c0$SR_parms_tv) && nrow(c0$SR_parms_tv) == 1L)
    expect_null(c1$SR_parms_tv)
    expect_equal(unname(unlist(c1$SR_parms["SR_regime", c("INIT", "Block", "Block_Fxn")])), c(0, 0, 0))
    expect_false("SR_regime_BLK1add_16" %in% rownames(reps[["01-no-regime"]]$parameters))
    # every other parameter block is structurally unchanged
    for (nm in c("MG_parms", "size_selex_parms", "age_selex_parms", "Q_parms")) {
        expect_equal(rownames(c1[[nm]]), rownames(c0[[nm]]), info = nm)
        expect_equal(c1[[nm]][, c("LO", "HI", "PHASE")], c0[[nm]][, c("LO", "HI", "PHASE")], info = nm)
    }
    # with no offset, the initial equilibrium is the unfished state's
    ts1 = as.data.table(reps[["01-no-regime"]]$timeseries)
    expect_equal(ts1[Era == "INIT", Recruit_0], ts1[Era == "VIRG", Recruit_0], tolerance = 1e-6)
})

test_that("the data are the same in every step", {
    d0 = inputs[["00-base"]]$dat
    d1 = inputs[["01-no-regime"]]$dat
    for (nm in c("catch", "CPUE", "lencomp", "agecomp")) {
        expect_equal(d1[[nm]], d0[[nm]], info = nm)
    }
    expect_equal(d1$sizefreq_data_list, d0$sizefreq_data_list)
})

test_that("the bridge cache matches the SS3 runs", {
    skip_if_not(file.exists(cache_file))
    cache = readRDS(cache_file)
    expect_equal(cache$steps, steps)
    for (i in seq_along(steps)) {
        expect_equal(cache$summary$nll[i], reps[[i]]$likelihoods_used["TOTAL", "values"],
                     tolerance = 1e-8, info = steps[i])
    }
})

# Initial-F sensitivities: 02 is 00 and 03 is 01 with the initial F of
# LL-n-A2 and OBJ-A5 estimated
init_f_steps = c("02-base-init-f", "03-no-regime-init-f")
init_f_parent = c("02-base-init-f" = "00-base", "03-no-regime-init-f" = "01-no-regime")
have_init_f = all(file.exists(file.path(dir_ss3, init_f_steps, "Report.sso")))

test_that("the initial-F sensitivities converged with an invertible Hessian", {
    skip_if_not(have_init_f)
    for (s in init_f_steps) {
        h = report_header(s)
        expect_match(h$version, "3.30.22.1", fixed = TRUE, info = s)
        expect_lt(h$gradient, 0.05)
        expect_true(is.finite(h$log_det), info = paste(s, "Hessian not inverted"))
    }
})

test_that("the initial-F sensitivities change only the initial F's phase", {
    skip_if_not(have_init_f)
    for (s in init_f_steps) {
        cs = suppressWarnings(r4ss::SS_read(file.path(dir_ss3, s), ss_new = TRUE, verbose = FALSE))
        cp = inputs[[init_f_parent[[s]]]]$ctl
        expect_true(all(cs$ctl$init_F$PHASE > 0), info = s)
        expect_true(all(cp$init_F$PHASE < 0), info = s)
        expect_equal(cs$ctl$init_F[, c("LO", "HI")], cp$init_F[, c("LO", "HI")], info = s)
        for (nm in c("MG_parms", "SR_parms", "size_selex_parms", "age_selex_parms", "Q_parms")) {
            expect_equal(cs$ctl[[nm]][, c("LO", "HI", "PHASE")], cp[[nm]][, c("LO", "HI", "PHASE")],
                         info = paste(s, nm))
        }
        expect_equal(is.null(cs$ctl$SR_parms_tv), is.null(cp$SR_parms_tv), info = s)
        for (nm in c("catch", "CPUE", "lencomp")) {
            expect_equal(cs$dat[[nm]], inputs[[init_f_parent[[s]]]]$dat[[nm]], info = paste(s, nm))
        }
        rep = suppressWarnings(r4ss::SS_output(file.path(dir_ss3, s), covar = FALSE,
                                               verbose = FALSE, printstats = FALSE))
        # estimating more parameters cannot make the fit worse than the parent's
        expect_lte(rep$likelihoods_used["TOTAL", "values"],
                   reps[[init_f_parent[[s]]]]$likelihoods_used["TOTAL", "values"] + 1e-3)
    }
})

test_that("the bridge cache holds the initial-F sensitivities", {
    skip_if_not(file.exists(cache_file) && have_init_f)
    cache = readRDS(cache_file)
    expect_equal(cache$init_f$steps, init_f_steps)
    expect_equal(cache$init_f$summary$step, c(steps, init_f_steps))
})
