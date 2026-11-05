# Nicholas Ducharme-Barth
# 2026/11/03
# Why do the length compositions push opal (O1, MFCL's selectivity fixed)
# towards an unexploited stock? Profiles over B0 (recruitment deviations
# 1955-2020 and the index catchabilities estimated at each value) for:
#
#   annual       opal as in O1: four seasons, annual age classes, length-at-age
#                at the start of the year, annual length compositions
#   annual-mid   as annual, with opal's growth fitted to MFCL's length-at-age
#                at mid-year (MFCL predicts each sample at its date, with
#                lengths that grow through the year; from M3's length.fit)
#   split        as annual, fisheries 1-7 (seasonal selectivity in MFCL) split
#                into one fishery per season (mfcl_split_seasonal_fisheries())
#   split-mid    split with mid-year growth
#   quarterly    quarterly age classes and time steps (mfcl_quarterly_ages()):
#                recruitment kept annual, MFCL's quarterly length-at-age and
#                seasonal selectivity, quarterly length compositions. B0 does
#                not have its usual meaning there (see mfcl-quarterly.R): its
#                grid is centred where opal's unfished recruitment matches
#                MFCL's (about 4 x the annual B0) and reported on the annual
#                model's scale.
#
# plus the length likelihood by fishery at MFCL's scale and at a large one.
# The opal code is not changed: every variant is a different set of opal
# inputs.
#
# OUTPUT: case-studies/swpo-alb/vignettes/data/swpo-alb-lf-diagnostics.rds
# (cached; SWPO_ALB_REFIT=true recomputes)
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Set-up (conversion, growth, O0)
    source(file.path(this.path::this.proj(), "case-studies", "swpo-alb", "scripts", "swpo-alb-opal-setup.r"))
    out_file = file.path(dir_cache, "swpo-alb-lf-diagnostics.rds")
    if (!refit && file.exists(out_file)) {
        message("Using cached ", out_file)
        diag = readRDS(out_file)
        print(diag$summary)
        quit(save = "no")
    }

    rdev_years = which(d0$years[-1] <= 2020)
    est_map = function(p, est_rdev) {
        m = fixed_map(p)
        m$log_cpue_q = factor(seq_along(p$log_cpue_q))
        i = rep(NA_integer_, length(p$rdev_y))
        i[est_rdev] = seq_along(est_rdev)
        m$rdev_y = factor(i)
        m
    }
    profile_fit = function(data, start, grid, est_rdev, sb_steps = seq_len(n_year)) {
        rbindlist(lapply(grid, function(lb) {
            pp = start
            pp$log_B0 = lb
            ob = MakeADFun(cmb(opal_model, data), pp, map = est_map(pp, est_rdev), silent = TRUE)
            op = tryCatch(nlminb(ob$par, ob$fn, ob$gr, control = list(eval.max = 1e4, iter.max = 1e4)),
                          error = function(e) NULL)
            if (is.null(op)) return(data.table(log_B0 = lb))
            r = ob$report(ob$env$last.par.best)
            sb = r$spawning_biomass_y[sb_steps]
            sb0 = r$spawning_biomass0_y[sb_steps]
            data.table(log_B0 = lb, objective = op$objective, convergence = op$convergence,
                       max_gradient = max(abs(ob$gr(op$par))),
                       cpue = sum(r$lp_cpue), lf = sum(r$lp_lf), rec = sum(r$lp_rec),
                       depletion_end = sb[length(sb)] / sb0[length(sb0)], sb_end = sb[length(sb)],
                       lp_lf = list(r$lp_lf))
        }), fill = TRUE)
    }
    p_start = p0
    p_start$init_rdev_a[] = 0
    p_start$rdev_y[] = 0
    grid = p0$log_B0 + c(-0.5, -0.25, 0, 0.25, 0.5, 1, 2, 3, 5)

#_____________________________________________________________________________
# Annual profiles (start-of-year and mid-year growth; whole and split fisheries)
    prof = list()
    prof$annual = profile_fit(d0, p_start, grid, rdev_years)
    prof$`annual-mid` = profile_fit(d0, with_growth(p_start, g_mid$par), grid, rdev_years)

    cs = mfcl_split_seasonal_fisheries(cv)
    ds = prep_lf_data(cs$data, as.data.frame(cs$lf), lf_switch = 1L, lf_var_adjust = cs$meta$lf_var_adjust)
    ds$lambdas = d0$lambdas
    ds$k_mult_ages = d0$k_mult_ages
    ds$growth_sd_model = d0$growth_sd_model
    ds$sel_fa_external = cs$selectivity$sel_fa_external
    ps = p_start
    ps$log_lf_tau = cs$parameters$log_lf_tau
    ps$log_init_F_f = cs$parameters$log_init_F_f
    ps$par_sel = cs$parameters$par_sel
    prof$split = profile_fit(ds, ps, grid, rdev_years)
    prof$`split-mid` = profile_fit(ds, with_growth(ps, g_mid$par), grid, rdev_years)

#_____________________________________________________________________________
# Quarterly ages
    qm = mfcl_quarterly_ages(cv, laa_month)
    dq = prep_lf_data(qm$data, as.data.frame(cv$lf_season), lf_switch = 1L,
                      lf_var_adjust = cv$meta$lf_var_adjust)
    dq$lambdas = d0$lambdas
    report_q = function(p) MakeADFun(cmb(opal_model, dq), p, map = fixed_map(p), silent = TRUE)$report()
    pq0 = mfcl_match_recruitment_quarterly(qm, mf, report_q, cv = cv)
    rq0 = report_q(pq0)
    q1 = which(qm$step_season == 1L)
    o0q_check = data.table(
        quantity = c("numbers at age, first year", "recruitment", "spawning biomass (max)",
                     "spawning biomass (2022)"),
        rel_diff = c(max(abs(rq0$number_ysa[1, 1, qm$age_season == 1L] / mf$N_ya[1, ] - 1)),
                     max(abs(rq0$number_ysa[q1, 1, 1] / mf$recruitment_y - 1)),
                     max(abs(rq0$spawning_biomass_y[q1] / mf$adult_biomass_y - 1)),
                     rq0$spawning_biomass_y[q1][n_year] / mf$adult_biomass_y[n_year] - 1))
    print(o0q_check)
    pq_start = pq0
    pq_start$init_rdev_a[qm$age_season == 1L] = 0
    pq_start$rdev_y[qm$rdev_steps] = 0
    # the plus group starts at its unfished equilibrium (12+), not one cohort
    pq_start$init_rdev_a[max(which(qm$age_season == 1L))] = -log(1 - exp(-d0$M[n_age]))
    est_q = qm$rdev_steps[d0$years[qm$rec_year_of_step] <= 2020]
    centre_q = pq0$log_B0
    prof$quarterly = profile_fit(dq, pq_start, centre_q + (grid - p0$log_B0), est_q, sb_steps = q1)
    prof$quarterly[, log_B0 := log_B0 - (pq0$log_B0 - p0$log_B0)]   # on the annual model's B0 scale

#_____________________________________________________________________________
# Length likelihood by fishery: annual model at MFCL's scale and at +3
    lf_fish = rep(d0$lf_fishery_f, d0$lf_n_f)
    by_fishery = rbindlist(lapply(c("annual", "annual-mid"), function(v) {
        pr = prof[[v]]
        i0 = which.min(abs(pr$log_B0 - p0$log_B0))
        i3 = which.min(abs(pr$log_B0 - (p0$log_B0 + 3)))
        data.table(variant = v, fishery = sort(unique(lf_fish)),
                   name = fishery_names[sort(unique(lf_fish))],
                   n_obs = as.vector(table(lf_fish)),
                   at_mfcl = as.vector(tapply(pr$lp_lf[[i0]], lf_fish, sum)),
                   at_plus3 = as.vector(tapply(pr$lp_lf[[i3]], lf_fish, sum)))
    }))
    by_fishery[, change := at_plus3 - at_mfcl]

#_____________________________________________________________________________
# Summary and cache
    profiles = rbindlist(lapply(prof, function(x) x[, !"lp_lf"]), idcol = "variant", fill = TRUE)
    profiles[, dlog_B0 := log_B0 - p0$log_B0]
    profiles[, `:=`(d_objective = objective - min(objective, na.rm = TRUE),
                    d_lf = lf - min(lf, na.rm = TRUE)), by = variant]
    summary_dt = profiles[, .(best_dlog_B0 = dlog_B0[which.min(objective)],
                              rise_at_plus5 = d_objective[which.max(dlog_B0)],
                              depletion_at_mfcl_scale = depletion_end[which.min(abs(dlog_B0))]),
                          by = variant]
    print(summary_dt)
    print(by_fishery[order(variant, change)])
    saveRDS(list(laa_month = laa_month, mid_year = mid_year, length_fit = lfit,
                 profiles = profiles, by_fishery = by_fishery, summary = summary_dt,
                 o0q_check = o0q_check,
                 quarterly_growth = qm$growth_check,
                 fishery_map_split = cs$meta$fishery_map,
                 log_B0_mfcl = p0$log_B0), out_file)
