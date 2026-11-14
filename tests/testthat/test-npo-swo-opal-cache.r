# Tests for the cached opal results of the NPO swordfish case study
# (case-studies/npo-swo/vignettes/data, written by
# case-studies/npo-swo/scripts/02-convert-fit-npo-swo-opal.r for both
# variants). The vignette reads only these caches, so they must agree with the
# SS3 bridge and with the stored fits, and the final fits must have converged.

library(testthat)
library(data.table)

skip_if_not_installed("opal")
library(opal)

# Define project paths
proj_dir = this.path::this.proj()
dir_data = file.path(proj_dir, "case-studies", "npo-swo", "vignettes", "data")

context("NPO swordfish opal caches")

bridge = readRDS(file.path(dir_data, "npo-swo-ss3-bridge-cache.rds"))$summary
variants = list(
    female = list(cache = "npo-swo-opal-cache.rds", fits = "npo-swo-opal-fits.rds",
                  ss3 = "SS3 B7", step = "07-index-month-1"),
    sexavg = list(cache = "npo-swo-opal-cache-sexavg.rds", fits = "npo-swo-opal-fits-sexavg.rds",
                  ss3 = "SS3 B7s", step = "07s-index-month-1"))
fitted_models = c("O1", "O2b", "O1-noLF", "O1-noS6", "D")

for (v in names(variants)) {
    vv = variants[[v]]
    cache = readRDS(file.path(dir_data, vv$cache))
    fits = readRDS(file.path(dir_data, vv$fits))
    S = cache$summary

    test_that(sprintf("%s: the cache's SS3 row is the bridge's %s", v, vv$step), {
        b = bridge[step == vv$step]
        expect_equal(S[model == vv$ss3]$depletion_2021, b$depletion, tolerance = 1e-6)
        expect_equal(S[model == vv$ss3]$B0, b$ssb0, tolerance = 1e-6)
        expect_equal(S[model == vv$ss3]$ssb_2021, b$ssb_2021, tolerance = 1e-6)
    })

    test_that(sprintf("%s: the final fits converged with positive-definite Hessians", v), {
        for (m in fitted_models) {
            expect_true(S[model == m]$check_passes, info = paste(v, m))
            expect_true(S[model == m]$positive_hessian, info = paste(v, m))
            expect_lt(S[model == m]$max_gradient, 1e-4)
        }
        # O2 is kept to show why O2b is needed
        expect_false(S[model == "O2"]$positive_hessian)
    })

    test_that(sprintf("%s: the cache agrees with the stored fits", v), {
        for (m in intersect(fitted_models, names(fits))) {
            expect_equal(S[model == m]$objective, fits[[m]]$fit$opt$objective,
                         tolerance = 1e-8, info = paste(v, m))
            expect_equal(S[model == m]$B0, exp(fits[[m]]$fit$opt$par[["log_B0"]]),
                         tolerance = 1e-8, info = paste(v, m))
        }
    })

    test_that(sprintf("%s: phase D runs four seasons with the annual models' catch", v), {
        dD = fits$D$data
        d2 = fits$O2b$data
        expect_equal(dD$n_season, 4L)
        expect_equal(d2$n_season, 1L)
        expect_equal(apply(dD$catch_obs_ysf, c(1, 3), sum), d2$catch_obs_ysf[, 1, ], tolerance = 1e-10)
        # the indices are observed in season 1 of their year
        expect_true(all((dD$cpue_data$ts - 1L) %% 4L == 0L))
        expect_equal((dD$cpue_data$ts - 1L) %/% 4L + 1L, d2$cpue_data$ts)
        # mid-year weight-at-age in the quarterly models
        expect_true(all(dD$weight[-1] > d2$weight[-1]))
        # O2b's fixed selectivity parameters stay fixed in D
        fixed_o2b = is.na(fits$O2b$map$par_sel)
        expect_true(all(is.na(fits$D$map$par_sel)[fixed_o2b]))
    })

    test_that(sprintf("%s: the cached trajectories agree with the summary and the fits", v), {
        # opal_report() would rebuild the stored fits, but opal's
        # object-identity check rejects stored fits on other machines (it
        # passes where they were made); these checks need no rebuild
        tr = cache$trajectory
        for (m in c(vv$ss3, "O0", fitted_models)) {
            x = tr[model == m]
            expect_equal(nrow(x), 47L, info = paste(v, m))
            expect_equal(last(x$ssb), S[model == m]$ssb_2021, tolerance = 1e-8, info = paste(v, m))
            expect_equal(x$depletion, x$ssb / S[model == m]$B0, tolerance = 1e-8, info = paste(v, m))
        }
    })
}
