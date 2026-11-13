# Extracted from test-ss3-to-opal-npo-swo.r:193

# prequel ----------------------------------------------------------------------
library(testthat)
library(data.table)
library(magrittr)
skip_if_not_installed("opal")
skip_if_not_installed("r4ss")
skip_if_not_installed("yaml")
library(opal)
proj_dir = this.path::this.proj()
dir_helper_fns_rtmb = file.path(proj_dir, "shared", "ss3-to-opal")
dir_config = file.path(proj_dir, "case-studies", "npo-swo", "config")
invisible(sapply(list.files(dir_helper_fns_rtmb, pattern = "\\.[Rr]$", full.names = TRUE), source))
context("SS3-to-opal converter: NPO swordfish")
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
assessment = ss3_opal_obj(converted, cfg, map = NULL, bounds = NULL,
                          control = list(eval.max = 10000, iter.max = 10000),
                          build = FALSE)
bridge_obj = ss3_use_external_selectivity(assessment, converted$selectivity$sel_fa_external) %>%
    opal_build()
rep_opal = opal_report(bridge_obj)
baranov_replay_rel = function(converted) {
    rep_ss3 = converted$ss3$rep
    fish = converted$meta$fisheries
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
    ssb / ts$SpawnBio - 1
}
cfg_s = ss3_read_config(file.path(dir_config, "npo-swo-b7s.yaml"))
converted_s = suppressWarnings(ss3_to_opal(cfg_s, dir_root = proj_dir))

# test -------------------------------------------------------------------------
expect_equal(converted_s$meta$fisheries$extract_fleets, fish$extract_fleets)
expect_equal(converted_s$meta$fisheries$index_map, fish$index_map)
