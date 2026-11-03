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
#       quarter (opal) or at each fishing incident's time with M and F acting
#       together until the next incident (MFCL's Baranov), with (b) MFCL's
#       season-averaged selectivity (opal) or its seasonal selectivity. The
#       opal/average replay reproduces O0; the MFCL/seasonal replay
#       approaches MFCL. The others measure the two missing opal features.
#   O1  MFCL's (season-averaged) selectivity-at-age fixed; B0, the recruitment
#       deviations (1955-2020; the last two at the curve, as MFCL fixes them
#       at the mean) and the two index catchabilities estimated. Growth,
#       steepness (0.8), sigma_R, M and the biology fixed at MFCL's values.
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
# Packages and paths
    library(data.table)
    library(magrittr)
    library(opal)
    library(RTMB)

    proj_dir  = this.path::this.proj()
    dir_cs    = file.path(proj_dir, "case-studies", "swpo-alb")
    dir_m3    = file.path(dir_cs, "model-files", "mfcl", "03-fixed-growth")
    dir_cache = file.path(dir_cs, "vignettes", "data")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
    source(file.path(proj_dir, "shared", "extractors", "mfcl", "parse-frq.r"))
    for (f in list.files(file.path(proj_dir, "shared", "mfcl-to-opal"), "\\.R$", full.names = TRUE)) source(f)

    refit = identical(tolower(Sys.getenv("SWPO_ALB_REFIT")), "true")
    # MFCL's recruitment penalty is equivalent to a CV of 0.71 (assessment
    # report, section 5.2.2): sigma_R = sqrt(log(1 + 0.71^2))
    sigma_r = sqrt(log(1 + 0.71^2))

#_____________________________________________________________________________
# Convert
    cv = mfcl_to_opal(dir_m3, "alb.frq", "03.par", "plot-03.par.rep",
                      index_fisheries = 14:15, sigma_r = sigma_r)
    d = cv$data
    mf = cv$mfcl
    n_age = d$n_age
    n_year = d$n_year
    fishery_names = c("1 LL DWFN 1a", "2 LL DWFN 1b", "3 LL DWFN 1c", "4 LL DWFN 1d",
                      "5 LL PICT 1ab", "6 LL PICT 1ab", "7 LL AU/NZ 1ab", "8 LL DWFN 1cd",
                      "9 LL PICT 1cd", "10 LL AU/NZ 1cd", "11 TR 1c", "12 TR 1d",
                      "13 DN 1cd", "14 LL index 1abcd", "15 TR index 1ef")

    # Length compositions: one per fishery and model year, multinomial, with
    # MFCL's divisors (fish flag 49) on the capped sample sizes
    d = prep_lf_data(d, as.data.frame(cv$lf), lf_switch = 1L,
                     lf_var_adjust = cv$meta$lf_var_adjust)
    d$lambdas = list(prior = 0)   # MFCL has no parameter priors here

#_____________________________________________________________________________
# Growth in opal's form (only the length compositions use it: weight,
# maturity and M are supplied at age)
    # VB (Schnute) with age-specific K multipliers at ages 1-4 for MFCL's
    # independent mean lengths at ages 2-4; SD linear in mean length between
    # MFCL's SDs at ages 1 and 12
    k_mult_ages = 1:4
    growth_obj = function(x) {
        mu = get_growth(n_age, A1 = 1, A2 = n_age, L1 = exp(x[1]), L2 = exp(x[2]), log_k = x[3],
                        log_k_mult = x[4:7], k_mult_ages = k_mult_ages)
        sum((mu - mf$mu_a)^2)
    }
    g = optim(c(log(mf$mu_a[1]), log(mf$mu_a[n_age]), cv$parameters$log_k, rep(0, 4)),
              growth_obj, method = "BFGS", control = list(maxit = 5000, reltol = 1e-14))
    mu_opal = get_growth(n_age, A1 = 1, A2 = n_age, L1 = exp(g$par[1]), L2 = exp(g$par[2]),
                         log_k = g$par[3], log_k_mult = g$par[4:7], k_mult_ages = k_mult_ages)
    sd_opal = get_sd_at_age(mu_opal, mu_opal[1], mu_opal[n_age], sd_model = "sd",
                            log_SD1 = log(mf$sd_a[1]), log_SD2 = log(mf$sd_a[n_age]),
                            A1 = 1, A2 = n_age)
    growth_check = data.table(age = seq_len(n_age), mfcl_length = mf$mu_a, opal_length = mu_opal,
                              mfcl_sd = mf$sd_a, opal_sd = as.numeric(sd_opal))
    print(growth_check[, .(rmse_length = sqrt(mean((opal_length - mfcl_length)^2)),
                           max_sd_diff = max(abs(opal_sd - mfcl_sd)))])
    d$k_mult_ages = k_mult_ages
    d$growth_sd_model = "sd"

    p_base = cv$parameters
    p_base$log_L1 = g$par[1]
    p_base$log_L2 = g$par[2]
    p_base$log_k = g$par[3]
    p_base$log_k_mult = g$par[4:7]
    p_base$log_SD1 = log(mf$sd_a[1])
    p_base$log_SD2 = log(mf$sd_a[n_age])

#_____________________________________________________________________________
# O0: opal at MFCL's values
    d0 = d
    d0$sel_fa_external = cv$selectivity$sel_fa_external
    fixed_map = function(p) lapply(p, function(z) factor(rep(NA, length(z))))
    report_at = function(p, data = d0) {
        MakeADFun(cmb(opal_model, data), p, map = fixed_map(p), silent = TRUE)$report()
    }
    p0 = mfcl_match_recruitment(list(data = d0, parameters = p_base, mfcl = mf), report_at)
    rep0 = report_at(p0)
    o0_check = data.table(
        quantity = c("numbers at age, first year", "recruitment"),
        max_rel_diff = c(max(abs(rep0$number_ysa[1, 1, ] / mf$N_ya[1, ] - 1)),
                         max(abs(rep0$number_ysa[2:n_year, 1, 1] / mf$recruitment_y[-1] - 1))))
    print(o0_check)

    # The index catchabilities that match MFCL's (on opal's normalised scale,
    # q is where opal's prediction centres on the observations)
    cpue = as.data.table(d0$cpue_data)
    cpue[, pred := rep0$cpue_pred]
    q0 = cpue[, .(log_q = mean(log(value)) - mean(log(pred))), by = index]
    p0$log_cpue_q = q0$log_q
    rep0 = report_at(p0)

#_____________________________________________________________________________
# Replays (standalone simulations, not opal code)
    inc = copy(mf$incidents)[!(fishery %in% cv$meta$index_fisheries) & catch > 0,
                             .(year, month, fishery, catch)]
    m1 = cv$meta$month_1
    inc[, `:=`(model_year = year + (m1 > 1 & month >= m1), model_month = (month - m1) %% 12 + 1)]
    inc[, `:=`(t = (model_month - 0.5) / 12, season = as.integer(ceiling(model_month / 3)))]
    sp = d0$spawning_potential
    replay = function(timing = c("opal", "mfcl"), selectivity = c("average", "seasonal")) {
        timing = match.arg(timing)
        selectivity = match.arg(selectivity)
        M = d0$M
        w = d0$weight
        sel_of = function(f, s) if (selectivity == "seasonal") mf$sel_fsa[f, s, ] else d0$sel_fa_external[f, ]
        N = mf$N_ya[1, ]
        sb = numeric(n_year)
        for (iy in seq_len(n_year)) {
            sb[iy] = sum(N * sp)
            yi = inc[model_year == d0$years[iy]]
            if (timing == "opal") {
                for (ss in seq_len(d0$n_season)) {
                    ys = yi[season == ss, .(catch = sum(catch)), by = fishery]
                    hr = numeric(n_age)
                    for (k in seq_len(nrow(ys))) {
                        f = ys$fishery[k]
                        cw = if (d0$catch_units_f[f] == 1) w else 1
                        hr = hr + ys$catch[k] / (sum(N * sel_of(f, ss) * cw) + 1e-6) * sel_of(f, ss)
                    }
                    N = N * (1 - hr) * exp(-M / d0$n_season)
                }
            } else {
                tt = sort(unique(yi$t))
                t_end = c(tt[-1], 1)
                t_now = 0
                for (j in seq_along(tt)) {
                    N = N * exp(-M * (tt[j] - t_now))
                    ii = yi[t == tt[j]]
                    K = nrow(ii)
                    sel = matrix(0, K, n_age)
                    cw = matrix(1, K, n_age)
                    for (k in seq_len(K)) {
                        sel[k, ] = sel_of(ii$fishery[k], ii$season[k])
                        if (d0$catch_units_f[ii$fishery[k]] == 1) cw[k, ] = w
                    }
                    Mp = M * (t_end[j] - tt[j])
                    F = ii$catch / rowSums(sweep(sel * cw, 2, N, "*"))
                    for (it in 1:30) {
                        Z = colSums(F * sel) + Mp
                        F = F * ii$catch / rowSums(sweep(F * sel * cw, 2, N / Z * (1 - exp(-Z)), "*"))
                    }
                    N = N * exp(-colSums(F * sel) - Mp)
                    t_now = t_end[j]
                }
                N = N * exp(-M * (1 - t_now))
            }
            if (iy < n_year) {
                N_next = c(mf$recruitment_y[iy + 1], N[-n_age])
                N_next[n_age] = N_next[n_age] + N[n_age]
                N = N_next
            }
        }
        sb
    }
    replays = rbindlist(lapply(list(c("opal", "average"), c("opal", "seasonal"),
                                    c("mfcl", "average"), c("mfcl", "seasonal")), function(a)
        data.table(timing = a[1], selectivity = a[2], year = d0$years,
                   sb = replay(a[1], a[2]), mfcl = mf$adult_biomass_y,
                   opal_o0 = rep0$spawning_biomass_y[seq_len(n_year)])))
    replays[, rel_mfcl := sb / mfcl - 1]
    print(replays[, .(rel_end = rel_mfcl[.N], rel_max = max(abs(rel_mfcl)),
                      vs_o0_max = max(abs(sb / opal_o0 - 1))), by = .(timing, selectivity)])

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
        rbindlist(lapply(list(O0 = rep0, O1 = rep1), function(r) data.table(
            year = d0$years, sb = r$spawning_biomass_y[seq_len(n_year)],
            dynamic_depletion = r$dynamic_depletion_y[seq_len(n_year)],
            recruits = r$number_ysa[seq_len(n_year), 1, 1])), idcol = "model"))
    cpue_fit = rbind(
        cpue[, .(model = "MFCL M3", index, year, ts, obs = value, se,
                 pred = NA_real_)],
        cpue[, .(model = "O0", index, year, ts, obs = value, se, pred = rep0$cpue_pred)],
        cpue[, .(model = "O1", index, year, ts, obs = value, se, pred = rep1$cpue_pred)])
    # MFCL's predicted index (mean-centred log) on opal's scale
    mfcl_pred = unlist(lapply(cv$meta$index_fisheries, function(f) mf$pred_cpue[[f]]))
    cpue_fit[model == "MFCL M3", pred := exp(mfcl_pred + log(obs) - mf$obs_cpue[cv$meta$index_fisheries] %>% unlist())]

    saveRDS(list(
        meta = cv$meta, fishery_names = fishery_names, growth = growth_check,
        o0_check = o0_check, trajectory = trajectory, replays = replays, cpue = cpue_fit,
        selectivity = list(seasonal = mf$sel_fsa, average = d0$sel_fa_external,
                           weights = cv$meta$sel_weights_fs),
        profile = profile, log_B0_mfcl = log_B0_mfcl,
        o1 = list(objective = opt1$objective, par = opt1$par,
                  lp = rep1[intersect(c("lp_cpue", "lp_lf", "lp_rec", "lp_penalty", "nll"), names(rep1))],
                  sigma_r = sigma_r, rdev_years = d0$years[rdev_years + 1L])
    ), file.path(dir_cache, "swpo-alb-opal-cache.rds"))
