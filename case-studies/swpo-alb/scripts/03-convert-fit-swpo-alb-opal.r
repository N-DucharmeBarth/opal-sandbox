# Nicholas Ducharme-Barth
# 2026/11/03
# Convert the final MFCL bridge model M3 (region 1 only, growth fixed, no
# conditional age-at-length; scripts/01-build-mfcl-bridge.r) to opal, test the
# conversion, replay the dynamics that opal cannot represent, and fit opal.
#
#   O0  opal at MFCL's values: MFCL's initial numbers-at-age and recruitment
#       (through init_rdev_a and rdev_y), biology at age, and MFCL's
#       selectivity-at-age averaged over seasons (weighted by catch) as
#       sel_fa_external; nothing estimated.
#   Replays: standalone simulations of the population with MFCL's estimates
#       (not opal code), crossing (a) catch removal at the start of each
#       quarter (opal) or in each of MFCL's fishing periods (the year's
#       incident dates, each with an equal share of the year's M, F and M
#       together: MFCL's Baranov), with (b) MFCL's
#       season-averaged selectivity (opal) or its seasonal selectivity. The
#       opal/average replay reproduces O0; the MFCL/seasonal replay
#       approaches MFCL. The others measure the two missing opal features.
#   O1  MFCL's (season-averaged) selectivity-at-age fixed; B0, the recruitment
#       deviations (1955-2020; the last two at the curve, as MFCL fixes them
#       at the mean) and the two index catchabilities estimated. Growth,
#       steepness (0.8), sigma_R, M and the biology fixed at MFCL's values.
#   O1-mid  O1 with opal's growth at MFCL's mid-year length-at-age.
#
# opal is quarterly here: n_season = 4 within MFCL's October-September model
# year (recruitment and ageing on 1 October), with annual age classes.
#
# OUTPUTS:
#   case-studies/swpo-alb/vignettes/data/swpo-alb-opal-cache.rds
#   case-studies/swpo-alb/vignettes/data/swpo-alb-opal-fits.rds  (fitted opal_obj's)
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Conversion, growth, O0 and the replays
    source(file.path(this.path::this.proj(), "case-studies", "swpo-alb", "scripts", "swpo-alb-opal-setup.r"))

#_____________________________________________________________________________
# O1: fit with MFCL's selectivity fixed
    est_map = function(p, rdev_years) {
        m = fixed_map(p)
        m$log_B0 = factor(1)
        m$log_cpue_q = factor(seq_along(p$log_cpue_q))
        i = rep(NA_integer_, length(p$rdev_y))
        i[rdev_years] = seq_along(rdev_years)
        m$rdev_y = factor(i)
        m
    }
    p1 = p0
    p1$init_rdev_a[] = 0        # MFCL: unfished equilibrium at the start
    p1$rdev_y[] = 0
    p1$log_B0 = p0$log_B0
    # rdev_y[y] sets recruitment at the start of year y + 1: 1955-2020 estimated
    rdev_years = which(d0$years[-1] <= 2020)
    o1 = opal_obj(data = d0, parameters = p1, map = est_map(p1, rdev_years),
                  metadata = list(stock = "SWPO albacore (MFCL M3, region 1)", source = "MFCL 2.2.7.0"))
    o1 = opal_build(o1)
    fits_file = file.path(dir_cache, "swpo-alb-opal-fits.rds")
    saved = if (!refit && file.exists(fits_file)) readRDS(fits_file) else list()
    fit_o1 = saved$O1 %||% opal_fit(o1, n_passes = 3L, check = FALSE,
                                    control = list(eval.max = 1e5, iter.max = 1e5))
    if (is.null(saved$O1)) saveRDS(list(O1 = fit_o1), fits_file, compress = "xz")
    rep1 = opal_report(fit_o1)

    # O1-mid: as O1, with opal's growth at MFCL's mid-year length-at-age (the
    # length compositions are otherwise predicted from start-of-year lengths;
    # see 04-length-comp-diagnostics.r)
    p1m = with_growth(p1, g_mid$par)
    o1m = opal_obj(data = d0, parameters = p1m, map = est_map(p1m, rdev_years),
                   metadata = list(stock = "SWPO albacore (MFCL M3, region 1)", source = "MFCL 2.2.7.0",
                                   growth = "MFCL mid-year length-at-age"))
    o1m = opal_build(o1m)
    fit_o1m = saved$`O1-mid` %||% opal_fit(o1m, n_passes = 3L, check = FALSE,
                                           control = list(eval.max = 1e5, iter.max = 1e5))
    if (is.null(saved$`O1-mid`)) saveRDS(list(O1 = fit_o1, `O1-mid` = fit_o1m), fits_file, compress = "xz")
    rep1m = opal_report(fit_o1m)
    opt1m = fit_o1m$fit$opt
    cat("O1-mid objective", opt1m$objective, " log_B0", opt1m$par[["log_B0"]], " (MFCL", p0$log_B0, ") max gradient",
        max(abs(opal_rtmb(fit_o1m, fresh = TRUE)$gr(opt1m$par))), "\n")
    opt1 = fit_o1$fit$opt
    cat("O1 objective", opt1$objective, " max gradient",
        max(abs(opal_rtmb(fit_o1, fresh = TRUE)$gr(opt1$par))), "\n")

#_____________________________________________________________________________
# O1 runs to the upper bound of log_B0. Profile: B0 fixed on a grid around
# MFCL's, recruitment deviations and catchabilities estimated
    lp_parts = function(r) c(cpue = sum(r$lp_cpue), lf = sum(r$lp_lf), rec = sum(r$lp_rec),
                             penalty = sum(r$lp_penalty))
    profile_file = file.path(dir_cache, "swpo-alb-opal-profile.rds")
    log_B0_mfcl = p0$log_B0
    grid = log_B0_mfcl + c(-0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5)
    profile = if (!refit && file.exists(profile_file)) readRDS(profile_file) else {
        start = p1
        out = rbindlist(lapply(grid, function(lb) {
            pp = start
            pp$log_B0 = lb
            mm = est_map(pp, rdev_years)
            mm$log_B0 = factor(NA)
            ob = MakeADFun(cmb(opal_model, d0), pp, map = mm, silent = TRUE)
            op = nlminb(ob$par, ob$fn, ob$gr, control = list(eval.max = 1e4, iter.max = 1e4))
            r = ob$report(ob$env$last.par.best)
            data.table(log_B0 = lb, objective = op$objective, max_gradient = max(abs(ob$gr(op$par))),
                       t(lp_parts(r)), depletion_end = r$dynamic_depletion_y[n_year],
                       sb_end = r$spawning_biomass_y[n_year])
        }))
        saveRDS(out, profile_file)
        out
    }
    profile[, delta := objective - min(objective)]
    print(profile[, .(log_B0 = round(log_B0, 2), delta = round(delta, 2), cpue = round(cpue, 2),
                      lf = round(lf, 2), rec = round(rec, 2), depletion_end = round(depletion_end, 3))])

#_____________________________________________________________________________
# Results
    sb0_mfcl = mf$adult_biomass0_y
    trajectory = rbind(
        data.table(model = "MFCL M3", year = d0$years, sb = mf$adult_biomass_y,
                   dynamic_depletion = mf$adult_biomass_y / sb0_mfcl,
                   recruits = mf$recruitment_y),
        rbindlist(lapply(list(O0 = rep0, O1 = rep1, `O1-mid` = rep1m), function(r) data.table(
            year = d0$years, sb = r$spawning_biomass_y[seq_len(n_year)],
            dynamic_depletion = r$dynamic_depletion_y[seq_len(n_year)],
            recruits = r$number_ysa[seq_len(n_year), 1, 1])), idcol = "model"))
    cpue_fit = rbind(
        cpue[, .(model = "MFCL M3", index, year, ts, obs = value, se,
                 pred = NA_real_)],
        cpue[, .(model = "O0", index, year, ts, obs = value, se, pred = rep0$cpue_pred)],
        cpue[, .(model = "O1", index, year, ts, obs = value, se, pred = rep1$cpue_pred)],
        cpue[, .(model = "O1-mid", index, year, ts, obs = value, se, pred = rep1m$cpue_pred)])
    # MFCL's predicted index (mean-centred log) on opal's scale
    mfcl_pred = unlist(lapply(cv$meta$index_fisheries, function(f) mf$pred_cpue[[f]]))
    cpue_fit[model == "MFCL M3", pred := exp(mfcl_pred + log(obs) - mf$obs_cpue[cv$meta$index_fisheries] %>% unlist())]

    fisheries = data.table(
        fishery = seq_len(d0$n_fishery), mfcl_diagnostic = c(1:13, 18, 19), name = fishery_names,
        units = ifelse(d0$catch_units_f == 1, "weight (t)", "numbers"),
        sel_seasons = vapply(seq_len(d0$n_fishery), function(f)
            length(unique(apply(round(mf$sel_fsa[f, , ], 10), 1, paste, collapse = " "))), integer(1)),
        lf_divisor = cv$meta$lf_var_adjust,
        lf_quarters = vapply(seq_len(d0$n_fishery), function(f) sum(cv$lf_season$fishery == f), integer(1)),
        lf_years = vapply(seq_len(d0$n_fishery), function(f) sum(cv$lf$fishery == f), integer(1)),
        catch_total = vapply(seq_len(d0$n_fishery), function(f) sum(d0$catch_obs_ysf[, , f]), numeric(1)))
    saveRDS(list(
        meta = cv$meta, fishery_names = fishery_names, fisheries = fisheries,
        n_lf_incident_samples = nrow(cv$lf_season), growth = growth_check,
        o0_check = o0_check, trajectory = trajectory, replays = replays, cpue = cpue_fit,
        selectivity = list(seasonal = mf$sel_fsa, average = d0$sel_fa_external,
                           weights = cv$meta$sel_weights_fs),
        profile = profile, log_B0_mfcl = log_B0_mfcl,
        o1 = list(objective = opt1$objective, par = opt1$par,
                  lp = rep1[intersect(c("lp_cpue", "lp_lf", "lp_rec", "lp_penalty", "nll"), names(rep1))],
                  sigma_r = sigma_r, rdev_years = d0$years[rdev_years + 1L]),
        o1_mid = list(objective = opt1m$objective, par = opt1m$par,
                      lp = rep1m[intersect(c("lp_cpue", "lp_lf", "lp_rec", "lp_penalty", "nll"), names(rep1m))],
                      growth = data.table(age = seq_len(n_age), mfcl_mid_year = mid_year,
                                          opal = get_growth(n_age, A1 = 1, A2 = n_age, L1 = exp(g_mid$par[1]),
                                                            L2 = exp(g_mid$par[2]), log_k = g_mid$par[3],
                                                            log_k_mult = g_mid$par[4:7], k_mult_ages = k_mult_ages)))
    ), file.path(dir_cache, "swpo-alb-opal-cache.rds"))
