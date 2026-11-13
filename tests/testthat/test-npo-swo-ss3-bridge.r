# Tests for the NPO swordfish SS3 bridge models
# (case-studies/npo-swo/model-files/ss3, built by
# case-studies/npo-swo/scripts/01-build-ss3-bridge.r).
#
# Each step must have converged with an invertible Hessian and no parameter
# on a bound, must contain exactly the structural change it is meant to make,
# must keep the data unchanged, and must agree with the bridge summary cached
# for the vignette. The tests read each step's inputs and *.ss_new files,
# ss3.par (final gradient) and ss_summary.sso (likelihood, parameters with
# their standard errors and position in their bounds, derived quantities);
# the full Report.sso is committed only for B7 and B7s.

library(testthat)
library(data.table)

skip_if_not_installed("r4ss")

# Define project paths
proj_dir = this.path::this.proj()
dir_ss3 = file.path(proj_dir, "case-studies", "npo-swo", "model-files", "ss3")
cache_file = file.path(proj_dir, "case-studies", "npo-swo", "vignettes", "data",
                       "npo-swo-ss3-bridge-cache.rds")

context("NPO swordfish SS3 bridge")

chain = c("02-block-split", "03-spawn-month-1", "04-annual", "05-annual-comps",
          "06-settle-month-1", "07-index-month-1")
chain_s = sub("^([0-9]+)", "\\1s", chain)
steps = c("00-base", "01-single-sex", chain, "01s-sex-averaged", chain_s)

inputs = lapply(steps, function(s) suppressWarnings(
    r4ss::SS_read(file.path(dir_ss3, s), ss_new = TRUE, verbose = FALSE)))
names(inputs) = steps
sums = lapply(steps, function(s) r4ss::SS_read_summary(file.path(dir_ss3, s, "ss_summary.sso")))
names(sums) = steps
dq = function(s, label) sums[[s]]$derived_quants[label, "Value"]
nll = function(s) sums[[s]]$likelihoods["TOTAL_LogL", 1]
active = function(s) {
    p = sums[[s]]$parameters
    p[p[["Active?"]] == "Act", ]
}
# the first line of ss3.par: "... Maximum gradient component = <value>"
gradient = function(s) {
    h = readLines(file.path(dir_ss3, s, "ss3.par"), n = 1L)
    as.numeric(sub(".*Maximum gradient component = ", "", h))
}

test_that("every bridge step was run with SS3 3.30.22.1", {
    for (s in steps) expect_match(sums[[s]]$header[1], "3.30.22.1", fixed = TRUE, info = s)
})

test_that("every bridge step converged with an invertible Hessian", {
    for (s in steps) {
        expect_lt(gradient(s), 2e-3)
        # standard errors come from the inverted Hessian
        se = active(s)$SE
        expect_true(all(is.finite(se) & se > 0), info = paste(s, "Hessian not inverted"))
    }
})

test_that("no estimated parameter is on a bound (female chain)", {
    # Range is the position within the bounds (0 = lower, 1 = upper); SS3
    # flags a parameter within 1% of a bound. In the sex-averaged chain F4's
    # logistic inflection lies within 1% of its 150 cm upper bound (148.8 cm;
    # the base case's is 148.3 cm).
    for (s in c("00-base", "01-single-sex", chain)) {
        r = active(s)$Range
        expect_true(all(r > 0.01 & r < 0.99), info = s)
    }
})

test_that("the base case reproduces the 2023 assessment with SS3 3.30.22.1", {
    expect_equal(nll("00-base"), 289.2996, tolerance = 1e-5)
    # ISC (2023), Table S2
    expect_equal(dq("00-base", "SSB_2021"), 35778, tolerance = 1e-4)
    expect_equal(dq("00-base", "SSB_MSY"), 16388, tolerance = 1e-4)
    expect_equal(dq("00-base", "annF_MSY"), 0.178, tolerance = 1e-2)
})

test_that("B1 is single sex with female biology", {
    dat = inputs[["01-single-sex"]]$dat
    ctl = inputs[["01-single-sex"]]$ctl
    expect_equal(dat$Nsexes, -1)
    expect_true(all(dat$lencomp$sex == 0))
    expect_false(any(grepl("^[fm][0-9]", names(dat$lencomp))))
    expect_equal(nrow(ctl$natM), 1L)
    expect_equal(as.numeric(ctl$natM[1, ]), as.numeric(inputs[["00-base"]]$ctl$natM[1, ]))
    expect_false(any(grepl("_Mal_", rownames(ctl$MG_parms))))
    mg0 = inputs[["00-base"]]$ctl$MG_parms
    fem = grep("_Fem_GP_1", rownames(mg0), value = TRUE)
    expect_equal(ctl$MG_parms[fem, "INIT"], mg0[fem, "INIT"])
})

test_that("B1 keeps one copy of each sex-combined composition", {
    lc0 = as.data.table(inputs[["00-base"]]$dat$lencomp)
    lc1 = as.data.table(inputs[["01-single-sex"]]$dat$lencomp)
    f = grep("^f[0-9]", names(lc0), value = TRUE)
    l = sub("^f", "l", f)
    expect_equal(unname(as.matrix(lc1[, ..l])), unname(as.matrix(lc0[, ..f])))
    expect_equal(lc1$Nsamp, lc0$Nsamp)
})

test_that("B1s averages M and growth over the sexes", {
    m = as.numeric(inputs[["01s-sex-averaged"]]$ctl$natM[1, ])
    m_f = as.numeric(inputs[["00-base"]]$ctl$natM[1, ])
    m_m = as.numeric(inputs[["00-base"]]$ctl$natM[2, ])
    expect_true(all(m >= pmin(m_f, m_m) - 1e-6 & m <= pmax(m_f, m_m) + 1e-6))
    mg = inputs[["01s-sex-averaged"]]$ctl$MG_parms
    l2 = mg["L_at_Amax_Fem_GP_1", "INIT"]
    expect_true(l2 > 206.4 && l2 < 226.3)
})

test_that("B2 replaces the selectivity blocks with fleets", {
    dat1 = inputs[["01-single-sex"]]$dat
    dat  = inputs[["02-block-split"]]$dat
    ctl  = inputs[["02-block-split"]]$ctl
    expect_equal(dat$Nfleets, 34L)
    expect_equal(ctl$N_Block_Designs, 0)
    expect_null(ctl$size_selex_parms_tv)
    expect_true(all(ctl$size_selex_parms$Block == 0))
    # catch moved to the block fleets only within the block years
    ct = dat$catch[dat$catch$year > 0 & dat$catch$catch > 0, ]
    expect_true(all(ct$year[ct$fleet %in% 28:30] %in% 1994:1997))
    expect_true(all(ct$year[ct$fleet %in% 31:34] %in% 2004:2009))
    expect_false(any(ct$fleet %in% c(1, 12) & ct$year %in% 1994:1997))
    expect_false(any(ct$fleet %in% c(2, 15, 16, 19) & ct$year %in% 2004:2009))
    expect_true(all(dat$CPUE$year[dat$CPUE$index == 30] %in% 1994:1997))
    # mirrors of the new fleets
    st = ctl$size_selex_types
    expect_equal(as.numeric(st[c(29, 30, 32, 33, 34), "Special"]), c(28, 28, 31, 31, 31))
    # the split fleets start from B1's estimated block replacement values
    # (B2's input control file; control.ss_new holds B2's estimates)
    sp1 = inputs[["01-single-sex"]]$ctl$size_selex_parms_tv
    sp2 = suppressWarnings(r4ss::SS_read(file.path(dir_ss3, "02-block-split"), ss_new = FALSE,
                                         verbose = FALSE))$ctl$size_selex_parms
    expect_equal(sp2[grep("\\(28\\)$", rownames(sp2)), "INIT"][1:4],
                 sp1[grep("^SizeSel_P_[1-4]_F1_", rownames(sp1)), "INIT"])
    # S2 1994-97 is a floated index of its own
    qo = ctl$Q_options
    expect_equal(as.numeric(qo[qo$fleet == 30, c("link", "float")]), c(1, 1))
})

test_that("B3-B7 make their one change each", {
    expect_equal(inputs[["02-block-split"]]$dat$spawn_month, 5)
    expect_equal(inputs[["03-spawn-month-1"]]$dat$spawn_month, 1)
    expect_equal(inputs[["03-spawn-month-1"]]$dat$nseas, 4)
    expect_equal(inputs[["04-annual"]]$dat$nseas, 1)
    expect_equal(inputs[["04-annual"]]$dat$months_per_seas, 12)
    expect_true(all(inputs[["04-annual"]]$dat$catch$seas == 1))
    expect_false(anyDuplicated(inputs[["05-annual-comps"]]$dat$lencomp[, c("year", "fleet")]) > 0)
    expect_equal(inputs[["05-annual-comps"]]$ctl$recr_dist_pattern$month, 7)
    expect_equal(inputs[["06-settle-month-1"]]$ctl$recr_dist_pattern$month, 1)
    expect_true(any(inputs[["06-settle-month-1"]]$dat$CPUE$month != 1))
    expect_true(all(inputs[["07-index-month-1"]]$dat$CPUE$month == 1))
})

test_that("combining the compositions (B5, B5s) leaves the estimates unchanged", {
    for (sfx in c("", "s")) {
        s4 = paste0("04", sfx, "-annual")
        s5 = paste0("05", sfx, "-annual-comps")
        expect_equal(active(s5)$Value, active(s4)$Value, tolerance = 1e-4, info = s5)
        expect_equal(dq(s5, "SSB_2021"), dq(s4, "SSB_2021"), tolerance = 1e-5, info = s5)
    }
})

test_that("the sex-averaged chain differs from the female chain only in M and growth", {
    for (i in seq_along(chain)) {
        a = inputs[[chain[i]]]
        b = inputs[[chain_s[i]]]
        expect_equal(a$dat$catch, b$dat$catch, info = chain_s[i])
        expect_equal(a$dat$CPUE, b$dat$CPUE, info = chain_s[i])
        expect_equal(a$dat$lencomp, b$dat$lencomp, info = chain_s[i])
        expect_equal(a$dat$spawn_month, b$dat$spawn_month, info = chain_s[i])
        expect_equal(a$ctl$recr_dist_pattern, b$ctl$recr_dist_pattern, info = chain_s[i])
        expect_equal(as.numeric(b$ctl$natM[1, ]), as.numeric(inputs[["01s-sex-averaged"]]$ctl$natM[1, ]))
        mg = c("L_at_Amin_Fem_GP_1", "L_at_Amax_Fem_GP_1", "VonBert_K_Fem_GP_1")
        expect_equal(b$ctl$MG_parms[mg, "INIT"], inputs[["01s-sex-averaged"]]$ctl$MG_parms[mg, "INIT"])
    }
})

test_that("every step keeps the catch, index and composition data", {
    tot = function(s) {
        d = inputs[[s]]$dat
        c(catch = sum(d$catch$catch[d$catch$year > 0]), cpue = sum(d$CPUE$obs),
          n_cpue = nrow(d$CPUE), nsamp = sum(d$lencomp$Nsamp[d$lencomp$year > 0]))
    }
    ref = tot("00-base")
    for (s in steps[-1]) expect_equal(tot(s), ref, tolerance = 1e-8, info = s)
})

test_that("the vignette's bridge cache matches the committed SS3 runs", {
    skip_if_not(file.exists(cache_file), "bridge cache not built")
    cache = readRDS(cache_file)$summary
    for (s in steps) {
        # ss_summary.sso prints six significant digits
        expect_equal(cache[step == s]$nll, nll(s), tolerance = 1e-5, info = s)
        expect_equal(cache[step == s]$depletion, dq(s, "Bratio_2021"), tolerance = 1e-5, info = s)
    }
})
