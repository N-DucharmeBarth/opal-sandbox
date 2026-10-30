# Nicholas Ducharme-Barth
# 2026/10/30
# Convert the SS3 bridge model B1 of the EPO bigeye tuna assessment to opal,
# test the conversion, and fit opal.
#
#   O0  bridge: opal at SS3's parameter values with SS3's selectivity-at-age;
#       nothing estimated
#   O1  SS3's selectivity-at-age held fixed; B0, recruitment deviations
#       (initial, main and late), the initial age structure and the survey
#       catchability estimated
#   O2  opal's own length-based double normals, started exactly at SS3's
#       curves; the selectivity parameters SS3 estimated are estimated, plus
#       everything O1 estimates. OBJ-Disc's age-based selectivity (SS3 ages
#       1-3) is an opal age-based double normal held at SS3's curve.
#
# Growth (SS3's growth cessation curve and SD-linear growth variability),
# steepness (1), sigma_R (0.6), M and the biology stay at SS3's values, as in
# SS3 (which estimates only SD_young of the growth parameters; it is held at
# SS3's estimate here). opal's default priors are switched off (SS3 has no
# parameter priors in the likelihood). The model is quarterly: each opal
# timestep is one SS3 "year" (a quarter). opal removes each quarter's catch at
# the start of the quarter, SS3 continuously (hybrid F).
#
# The compositions are the OBJ length compositions (2 cm bins) and the
# longline and survey generalised size compositions (10 cm bins), as two
# composition bin schemes.
#
# OUTPUTS:
#   case-studies/epo-bet/vignettes/data/epo-bet-opal-cache.rds
#   case-studies/epo-bet/vignettes/data/epo-bet-opal-fits.rds  (fitted opal_obj's)
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Packages
    library(data.table)
    library(magrittr)
    library(opal)
    library(r4ss)
    library(yaml)

#_____________________________________________________________________________
# Project paths
    proj_dir  = this.path::this.proj()
    dir_fns   = file.path(proj_dir, "shared", "ss3-to-opal")
    dir_cfg   = file.path(proj_dir, "case-studies", "epo-bet", "config")
    dir_cache = file.path(proj_dir, "case-studies", "epo-bet", "vignettes", "data")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
    invisible(sapply(list.files(dir_fns, pattern = "\\.[Rr]$", full.names = TRUE), source))

#_____________________________________________________________________________
# Convert and test the conversion
    cfg = ss3_read_config(file.path(dir_cfg, "epo-bet-b1.yaml"))
    converted = ss3_to_opal(cfg, dir_root = proj_dir)
    d = converted$data
    n_f = d$n_fishery

    # Length compositions (2 cm) and generalised size compositions (10 cm) as
    # two composition bin schemes. Each fishery has one or the other, so the
    # SS3 variance adjustments (factors 4 and 7) combine into one vector.
    lf_wide = dcast(converted$lf, fishery + year + month + ts ~ bin,
                    value.var = "value", fill = 0) %>%
              setorder(fishery, ts) %>%
              as.data.frame()
    lf_wide = lf_wide[lf_wide$fishery %in% as.integer(cfg$lf$keep_fisheries), ]
    sf = ss3_sizefreq_to_opal(converted$ss3, cfg, converted$meta$dims, converted$meta$fisheries)
    sf_fisheries = as.integer(cfg$sizefreq$keep_fisheries)
    stopifnot("A fishery has both length and size compositions" =
                  !length(intersect(sf_fisheries, cfg$lf$keep_fisheries)))
    var_adjust = converted$meta$lf_var_adjust
    var_adjust[sf_fisheries] = sf$var_adjust[sf_fisheries]
    top = d$len_bin_start + d$n_len * d$len_bin_width
    d = prep_lf_data(
        d, list(lf_wide, sf$wide),
        lf_switch         = 1L,                       # multinomial, as SS3
        lf_keep_fisheries = sort(c(as.integer(cfg$lf$keep_fisheries), sf_fisheries)),
        lf_var_adjust     = var_adjust,
        lf_cap            = cfg$lf$cap,
        lf_addtocomp      = as.numeric(cfg$lf$addtocomp),
        lf_upper_edge     = c(top, sf$upper_edge)
    )
    # SS3 has no parameter priors in the likelihood
    d$lambdas = list(prior = 0)
    # OBJ-Disc is age-based in opal (SS3 age pattern 11, ages 1-3)
    d$sel_basis_f = rep(1L, n_f)
    d$sel_basis_f[as.integer(cfg$selectivity$age_based)] = 2L
    converted$data = ss3_make_portable(d)

    checks = validate_opal_inputs(converted, strict = FALSE)
    recon  = ss3_reconcile(converted, tol = 1e-6)
    print(recon[, .(quantity, rel_diff, pass)])
    stopifnot("Structural validation failed" = all(checks$pass),
              "Reconciliation with SS3 failed" = all(recon$pass, na.rm = TRUE))

    assessment = ss3_opal_obj(converted, cfg,
                              control = list(eval.max = 10000, iter.max = 10000),
                              build = FALSE)

#_____________________________________________________________________________
# Selectivity
    rep_ss3 = converted$ss3$rep
    end_q = cfg$time$last_ts
    ss3_fleet_f = c(converted$meta$fisheries$extract_fleets,
                    converted$meta$fisheries$index_map$ss3_fleet)
    fleet_names = converted$meta$dims$.fleet_names[ss3_fleet_f]
    sel_ss3 = converted$selectivity$sel_fa_external

    # opal's double normals at SS3's curves (see ss3_seed_par_sel_size_dblnorm())
    seed = ss3_seed_par_sel_size_dblnorm(converted$ss3, converted$meta$fisheries,
                                         converted$data$len_mid,
                                         converted$parameters$par_sel)
    sizeselex = as.data.table(rep_ss3$sizeselex)[Factor == "Lsel" & Yr == end_q]
    len_cols = setdiff(names(sizeselex), c("Factor", "Fleet", "Yr", "Sex", "Label"))
    # sizeselex columns are labelled by population bin midpoints
    len_cols = len_cols[match(converted$data$len_mid, as.numeric(len_cols))]
    stopifnot("SS3 sizeselex does not cover opal's length bins" = !anyNA(len_cols))

    # OBJ-Disc: an age-based double normal fitted to SS3's selectivity-at-age
    # (1 at ages 1-3, 0 elsewhere) and held there
    age_a = converted$data$age_a
    for (f in as.integer(cfg$selectivity$age_based)) {
        target = sel_ss3[f, ]
        obj_fn = function(p) sum((sel_double_normal(age_a, p) - target)^2)
        # peak at the first selected age, top to the last, steep sides
        sel_ages = age_a[target > 0.5]
        s_a = stats::sd(age_a)
        start = c((min(sel_ages) - mean(age_a)) / s_a,
                  -log((0.99 * max(age_a) - min(sel_ages) - 1) / (max(sel_ages) - min(sel_ages) - 1 + 0.5) - 1),
                  -8, -8, -15, -15)
        o = stats::optim(start, obj_fn, method = "Nelder-Mead", control = list(maxit = 20000, reltol = 1e-14))
        seed$par_sel[f, ] = o$par
    }
    seed_check = data.table(
        fishery    = seq_along(ss3_fleet_f),
        ss3_fleet  = ss3_fleet_f,
        name       = fleet_names,
        source     = seed$source_fleet,
        max_abs_diff = vapply(seq_along(ss3_fleet_f), function(i) {
            if (d$sel_basis_f[i] == 2L)
                return(max(abs(sel_double_normal(age_a, seed$par_sel[i, ]) - sel_ss3[i, ])))
            max(abs(sel_double_normal(converted$data$len_mid, seed$par_sel[i, ]) -
                    as.numeric(sizeselex[Fleet == seed$source_fleet[i], ..len_cols])))
        }, numeric(1L)))
    print(seed_check)

    # Which double-normal parameters SS3 estimated, by opal fishery (same order:
    # peak, top, ascending width, descending width, start, end)
    pars_ss3 = as.data.table(rep_ss3$parameters)
    est_cols = lapply(seq_along(ss3_fleet_f), function(i) {
        sf_i = seed$source_fleet[i]
        if (is.na(sf_i) || i %in% as.integer(names(cfg$selectivity$mirror))) return(integer())
        if (!i %in% as.integer(cfg$selectivity$estimate)) return(integer())
        p = pars_ss3[grepl("^Size_DblN_", Label) & grepl(sprintf("\\(%d\\)$", sf_i), Label)]
        which(!is.na(p$Active_Cnt))
    })

#_____________________________________________________________________________
# Maps
    n_age = converted$data$n_age
    n_year = converted$data$n_year
    # Initial age structure: SS3 estimates the first quarter's recruitment (age
    # 0, the first main deviation) and the 24 early deviations (ages 1-24)
    early_start = converted$ss3$ctl$MainRdevYrFirst + converted$ss3$ctl$recdev_early_start
    n_init = cfg$time$first_ts - early_start + 1L
    map_init_rdev = factor(c(seq_len(n_init), rep(NA, n_age - n_init)))
    # Recruitment deviations: opal's rdev_y[y] sets the age-0 class at the
    # start of timestep y + 1, so the last one belongs to a quarter after the
    # model (205) and is fixed. SS3 estimates the main (to 202) and late
    # (203-204) deviations.
    map_rdev = factor(c(seq_len(n_year - 1L), NA))
    base_map = function(x) {
        m = x$map %||% get_map(x$parameters)
        m$init_rdev_a = map_init_rdev
        m$rdev_y = map_rdev
        m
    }

    # Selectivity map: SS3's estimated columns; mirrors share their source's
    sel_map = function(par_sel) {
        m = matrix(NA_integer_, nrow(par_sel), ncol(par_sel))
        lvl = 0L
        for (f in seq_along(est_cols)) {
            cols = est_cols[[f]]
            if (!length(cols)) next
            m[f, cols] = lvl + seq_along(cols)
            lvl = lvl + length(cols)
        }
        for (f in names(cfg$selectivity$mirror)) {
            m[as.integer(f), ] = m[as.integer(cfg$selectivity$mirror[[f]]), ]
        }
        factor(as.vector(m))
    }

    # Bounds for the estimated selectivity parameters: SS3's own bounds on
    # opal's scale
    sel_bounds = function() {
        x = converted$data$len_mid
        mu = mean(x)
        sd_x = stats::sd(x)
        to_opal = function(v, col) switch(col, (v - mu) / sd_x, v, v - 2 * log(sd_x), v - 2 * log(sd_x))
        out = list(lower = numeric(), upper = numeric())
        for (f in seq_along(est_cols)) {
            sf_i = seed$source_fleet[f]
            for (col in est_cols[[f]]) {
                p = pars_ss3[grepl("^Size_DblN_", Label) & grepl(sprintf("\\(%d\\)$", sf_i), Label)][col]
                out$lower = c(out$lower, to_opal(p$Min, col))
                out$upper = c(out$upper, to_opal(p$Max, col))
            }
        }
        out
    }

#_____________________________________________________________________________
# O0: bridge at SS3's values
    # SS3's survey catchability (floated) as opal's starting value
    q_ss3 = pars_ss3[grepl("^LnQ_base", Label)]$Value
    a0 = assessment
    p0 = a0$parameters
    p0$log_cpue_q = q_ss3
    a0 = opal_update(a0, parameters = p0)
    # SS3's ln(Q) (9.3) is outside opal's default bounds for log_cpue_q
    # (log(0.001) to log(10)), which opal_build() checks. Build a copy at q = 1
    # for opal's default bounds, widen the catchability's to SS3's +/-15, and
    # build with them. Run after every map change (opal resets bounds then).
    build_with_q_bounds = function(x, sel = NULL) {
        # bounds stored for an earlier map no longer match its parameters
        x = opal_update(x, bounds = NULL)
        p = x$parameters
        # opal's default bounds come from a copy whose catchability and
        # selectivity start inside them
        tmp = p
        tmp$log_cpue_q = rep(0, length(p$log_cpue_q))
        tmp$par_sel = p$par_sel * 0
        x1 = opal_build(opal_update(x, parameters = tmp))
        b = x1$bounds
        nm = names(opal_rtmb(x1, fresh = TRUE)$par)
        i = which(nm == "log_cpue_q")
        b$lower[i] = -15
        b$upper[i] = 15
        # the estimated selectivity parameters: SS3's own bounds (sel_bounds())
        if (!is.null(sel)) {
            i = which(nm == "par_sel")
            stopifnot(length(i) == length(sel$lower))
            b$lower[i] = sel$lower
            b$upper[i] = sel$upper
        }
        opal_update(x, bounds = list(lower = b$lower, upper = b$upper)) %>% opal_build()
    }
    o0 = ss3_use_external_selectivity(a0, sel_ss3) %>% build_with_q_bounds()
    bridge = ss3_bridge_table(converted, opal_report(o0))
    print(summary(bridge$rel_opal))

    # Catch timing and catch weight at SS3's values: the converted inputs
    # replayed with SS3's selectivity-at-age and recruitment, removing each
    # quarter's catch (a) at the start of the quarter (opal), (b) at mid-quarter
    # (Pope: half the quarter's M, the catch, the other half) or (c) through
    # the quarter (Baranov, as SS3), with catch in weight turned into numbers
    # at start-of-quarter weight, mid-quarter weight (opal's weight here) or
    # each fleet's SS3 catch weight-at-age (ageselex "bodywt").
    replay_timing = function(timing = c("start", "pope", "baranov"), W) {
        timing = match.arg(timing)
        dd = converted$data
        p  = converted$parameters
        sp = dd$spawning_potential
        M  = dd$M
        n  = dd$n_age
        S_a = exp(-M)
        rel = cumprod(c(1, S_a[-n]))
        rel[n] = rel[n] / (1 - S_a[n])
        R0  = exp(p$log_B0) / sum(sp * rel)
        sig = exp(p$log_sigma_r)
        num = dd$catch_units_f != 1
        cw = function(f) if (num[f]) 1 else W[f, ]
        N = R0 * rel * exp(p$init_rdev_a - dd$init_bias_adj_a * 0.5 * sig^2)
        ssb = numeric(dd$n_year)
        for (y in seq_len(dd$n_year)) {
            ssb[y] = sum(N * sp)
            cy = dd$catch_obs_ysf[y, 1, ]
            fl = which(cy > 0)
            if (timing == "baranov") {
                F_f = numeric(length(cy))
                for (f in fl) F_f[f] = cy[f] / sum(N * sel_ss3[f, ] * cw(f))
                for (it in 1:50) {
                    Z = M + colSums(F_f * sel_ss3)
                    for (f in fl) F_f[f] = F_f[f] * cy[f] /
                        sum(F_f[f] * sel_ss3[f, ] / Z * N * (1 - exp(-Z)) * cw(f))
                }
                N_end = N * exp(-(M + colSums(F_f * sel_ss3)))
            } else {
                Nc = if (timing == "pope") N * exp(-M / 2) else N
                hr = numeric(n)
                for (f in fl) hr = hr + cy[f] / sum(Nc * sel_ss3[f, ] * cw(f)) * sel_ss3[f, ]
                hr = pmin(hr, 0.999)
                N_end = if (timing == "pope") Nc * (1 - hr) * exp(-M / 2) else N * (1 - hr) * S_a
            }
            N_next = c(0, N_end[-n])
            N_next[n] = N_next[n] + N_end[n]
            N_next[1] = R0 * exp(p$rdev_y[y] - dd$bias_adj_y[y] * 0.5 * sig^2)
            N = N_next
        }
        ssb
    }
    eg_beg = ss3_endgrowth_table(rep_ss3, converted$meta$dims, "Beg")
    eg_mid = ss3_endgrowth_table(rep_ss3, converted$meta$dims, "Mid")
    bodywt = as.data.table(rep_ss3$ageselex)[Factor == "bodywt" & Yr == end_q]
    W_fleet = as.matrix(bodywt[match(ss3_fleet_f, Fleet), as.character(converted$meta$dims$.ss3_ages), with = FALSE])
    W_beg = matrix(eg_beg$weight, n_f, d$n_age, byrow = TRUE)
    W_mid = matrix(eg_mid$weight, n_f, d$n_age, byrow = TRUE)
    ssb_ss3_q = bridge$ss3
    timing = rbindlist(lapply(list(
        list("start", W_beg,   "start of quarter", "start-of-quarter"),
        list("start", W_mid,   "start of quarter", "mid-quarter (opal here)"),
        list("start", W_fleet, "start of quarter", "SS3 fleet catch weight"),
        list("pope",  W_fleet, "mid-quarter (Pope)", "SS3 fleet catch weight"),
        list("baranov", W_fleet, "through the quarter (Baranov, SS3)", "SS3 fleet catch weight")),
        function(a) data.table(ts = seq_len(n_year), removal = a[[3]], weight = a[[4]],
                               ssb = replay_timing(a[[1]], a[[2]]), ss3 = ssb_ss3_q)))
    timing[, rel := ssb / ss3 - 1]
    print(timing[, .(rel_end = rel[.N], rel_max = max(abs(rel))), by = .(removal, weight)])

#_____________________________________________________________________________
# Fits
    fit_opal = function(x) {
        t0  = proc.time()[["elapsed"]]
        fit = opal_fit(x, n_passes = 3L, check = FALSE,
                       control = list(eval.max = 1e5, iter.max = 1e5))
        obj = opal_rtmb(fit, fresh = TRUE)
        opt = fit$fit$opt
        # Newton steps to bring the gradient down; keep the nlminb optimum if a
        # step lands where the objective is undefined or outside the bounds
        pol = tryCatch(ss3_newton_polish(obj, opt$par), error = function(e) NULL)
        b = x$bounds
        if (is.null(pol) || !is.finite(pol$objective) || pol$objective > opt$objective + 1e-8 ||
            any(pol$par < b$lower | pol$par > b$upper)) {
            g = obj$gr(opt$par)
            pol = list(par = opt$par, objective = opt$objective,
                       max_gradient = max(abs(g)), steps = 0L)
        }
        opt$par = setNames(pol$par, names(opt$par))
        opt$objective = pol$objective
        fit = opal_attach_fit(fit, opt, check = TRUE,
                              check_args = list(gradient_tolerance = 1e-4))
        attr(fit, "seconds") = proc.time()[["elapsed"]] - t0
        attr(fit, "max_gradient") = pol$max_gradient
        fit
    }

    # The fits are saved as soon as they finish and reused on later runs
    # unless EPO_BET_REFIT=true
    fits_file = file.path(dir_cache, "epo-bet-opal-fits.rds")
    refit = !file.exists(fits_file) || identical(Sys.getenv("EPO_BET_REFIT"), "true")
    saved = if (refit) list() else readRDS(fits_file)

    # O1: SS3 selectivity fixed
    o1 = o0
    o1 = opal_update(o1, map = base_map(o1)) %>% build_with_q_bounds()
    fit_o1 = if (refit) fit_opal(o1) else saved$O1
    utils::str(fit_o1$validation$fit, max.level = 2, give.attr = FALSE)

    # O2: opal selectivity, started at SS3's curves and from O1's estimates
    # for everything else
    o2 = ss3_use_opal_selectivity(a0, par_sel = seed$par_sel)
    p2 = fit_o1$fit$parameters
    p2$par_sel = seed$par_sel
    o2 = opal_update(o2, parameters = p2)
    m2 = base_map(o2)
    m2$par_sel = sel_map(seed$par_sel)
    o2 = opal_update(o2, map = m2) %>% build_with_q_bounds(sel = sel_bounds())
    fit_o2 = if (refit) fit_opal(o2) else saved$O2
    if (refit) saveRDS(list(O1 = fit_o1, O2 = fit_o2), fits_file, compress = "xz")

#_____________________________________________________________________________
# Diagnostics: what moves O1 away from SS3
    # O1-noComp  O1 without the compositions (lambda 0), so only the survey
    #            and the recruitment penalty inform the estimates
    with_lambdas = function(x, lam) {
        dd = x$data
        dd$lambdas = c(list(prior = 0), lam)
        opal_update(x, data = dd)
    }
    o1_nocomp = with_lambdas(o1, list(lf = 0)) %>% opal_build()
    fit_o1_nocomp = fit_opal(o1_nocomp)

#_____________________________________________________________________________
# Results
    qtr_to_year = function(q) 1975 + (q - 17) / 4
    qtrs = cfg$time$first_ts:cfg$time$last_ts
    ts_ss3 = as.data.table(rep_ss3$timeseries)[Era == "TIME" & Yr %in% qtrs]
    ssb0_ss3 = ss3_derived_value(rep_ss3, c("SSB_Virgin", "SSB_Unfished"))

    models = list(O0 = o0, O1 = fit_o1, O2 = fit_o2, `O1-noComp` = fit_o1_nocomp)
    trajectory = rbind(
        data.table(model = "SS3 B1", quarter = qtrs, year = qtr_to_year(qtrs),
                   ssb = ts_ss3$SpawnBio, depletion = ts_ss3$SpawnBio / ssb0_ss3,
                   recruits = as.data.table(rep_ss3$natage)[Era == "TIME" & `Beg/Mid` == "B" &
                                                                Yr %in% qtrs][["0"]]),
        rbindlist(lapply(names(models), function(id) {
            r = opal_report(models[[id]])
            data.table(model = id, quarter = qtrs, year = qtr_to_year(qtrs),
                       ssb = r$spawning_biomass_y[seq_len(n_year)],
                       depletion = r$static_depletion_y[seq_len(n_year)],
                       recruits = r$number_ysa[seq_len(n_year), 1, 1])
        })))

    summary_dt = rbindlist(lapply(names(models), function(id) {
        x = models[[id]]
        r = opal_report(x)
        fitted = !is.null(x$fit$opt)
        data.table(
            model        = id,
            n_par        = if (fitted) length(x$fit$opt$par) else 0L,
            objective    = if (fitted) x$fit$opt$objective else as.numeric(opal_rtmb(x, fresh = TRUE)$fn()),
            max_gradient = if (fitted) attr(x, "max_gradient") else NA_real_,
            check_passes = if (fitted) isTRUE(x$validation$fit$passes) else NA,
            positive_hessian = if (fitted) isTRUE(x$validation$fit$metrics$positive_hessian) else NA,
            seconds      = if (fitted) attr(x, "seconds") else NA_real_,
            B0           = exp((if (fitted) x$fit$opt$par[["log_B0"]] else x$parameters$log_B0)),
            ssb_end      = r$spawning_biomass_y[n_year],
            depletion_end = r$static_depletion_y[n_year],
            lp_cpue      = sum(r$lp_cpue),
            lp_lf        = sum(r$lp_lf),
            lp_rec       = sum(r$lp_rec)
        )
    }))
    summary_dt = rbind(data.table(model = "SS3 B1", n_par = rep_ss3$N_estimated_parameters,
                                  B0 = ssb0_ss3, ssb_end = last(ts_ss3$SpawnBio),
                                  depletion_end = last(ts_ss3$SpawnBio) / ssb0_ss3),
                       summary_dt, fill = TRUE)
    print(summary_dt)

    # Index fits
    cpue = rbindlist(lapply(c("O0", "O1", "O2"), function(id) {
        dd = as.data.table(models[[id]]$data$cpue_data)
        dd[, .(model = id, quarter = year, index, obs = value, se,
               pred = opal_report(models[[id]])$cpue_pred)]
    }))
    cpue_ss3 = as.data.table(rep_ss3$cpue)[Fleet %in% cfg$index$keep_ss3_fleets,
                                           .(model = "SS3 B1", quarter = Yr,
                                             index = match(Fleet, cfg$index$keep_ss3_fleets),
                                             obs = Obs, se = SE, pred = Exp)]
    cpue = rbind(cpue_ss3, cpue)
    cpue[, year := qtr_to_year(quarter)]

    # Composition fits, aggregated over time by fishery: observed and
    # predicted proportions on each scheme's bins
    comp_agg = function(x, id) {
        dd = x$data
        r = opal_report(x)
        # opal's predictions are a list by fishery (rows in lf_row_fi order);
        # lf_scheme follows the same fishery-grouped order
        pos = 0L
        rbindlist(lapply(seq_along(dd$lf_row_fi), function(k) {
            rows = dd$lf_row_fi[[k]]
            s_k = dd$lf_scheme[pos + 1L]
            pos <<- pos + length(rows)
            edges = dd$lf_scheme_edges[[s_k]]
            obs = do.call(rbind, lapply(rows, function(i) dd$lf_obs_in[[i]]))
            pred = r$lf_pred[[k]]
            n = dd$lf_n[rows]
            o = colSums(obs / rowSums(obs) * n)
            p = colSums(pred / rowSums(pred) * n)
            data.table(model = id, fishery = as.integer(names(dd$lf_row_fi)[k]),
                       bin = utils::head(edges, -1L), obs = o / sum(o), pred = p / sum(p))
        }))
    }
    comps = tryCatch(rbindlist(lapply(c("O1", "O2"), function(id) comp_agg(models[[id]], id))),
                     error = function(e) { message("comp_agg: ", conditionMessage(e)); NULL })

    # Selectivity-at-age: SS3 against opal (O2)
    sel_o2  = opal_report(fit_o2)$sel_fya[, n_year, ]
    selectivity = rbind(
        data.table(model = "SS3 B1", fishery = rep(seq_len(nrow(sel_ss3)), times = ncol(sel_ss3)),
                   age = rep(age_a, each = nrow(sel_ss3)), sel = as.vector(sel_ss3)),
        data.table(model = "O2", fishery = rep(seq_len(nrow(sel_o2)), times = ncol(sel_o2)),
                   age = rep(age_a, each = nrow(sel_o2)), sel = as.vector(sel_o2)))
    # opal's starting age selectivity (its PLA) against SS3's, before fitting
    sel_o2_start = opal_report(o2)$sel_fya[, n_year, ]
    sel_start_diff = data.table(fishery = seq_len(nrow(sel_ss3)), name = fleet_names,
                                max_abs_diff = apply(abs(sel_o2_start - sel_ss3), 1, max),
                                age0_opal = sel_o2_start[, 1])

    gf = converted$meta$growth_fit
    cache = list(
        created     = format(Sys.time(), tz = "UTC", usetz = TRUE),
        session     = list(r = as.character(getRversion()),
                           opal = as.character(packageVersion("opal")),
                           r4ss = as.character(packageVersion("r4ss"))),
        fleets      = data.table(fishery = seq_along(ss3_fleet_f), ss3_fleet = ss3_fleet_f,
                                 name = fleet_names, units = converted$data$catch_units_f,
                                 comps = fifelse(seq_along(ss3_fleet_f) %in% cfg$lf$keep_fisheries, "length (2 cm)",
                                         fifelse(seq_along(ss3_fleet_f) %in% sf_fisheries, "size (10 cm)", "none")),
                                 var_adjust = var_adjust),
        checks      = checks,
        recon       = recon,
        growth_fit  = gf[c("L0", "Linf", "rmax", "k_gc", "t50", "age_offset", "cv_pattern",
                           "sd_model", "v1", "v2", "v2_at_oldest", "rmse", "max_abs_error",
                           "max_rel_error_below_plus", "sd_max_rel_error_below_plus",
                           "laa_ss3", "laa_opal", "sd_ss3", "sd_opal")],
        bridge      = bridge,
        timing      = timing,
        weights     = data.table(age = age_a, beg = eg_beg$weight, mid = eg_mid$weight),
        seed_check  = seed_check,
        sel_start_diff = sel_start_diff,
        summary     = summary_dt,
        trajectory  = trajectory,
        cpue        = cpue,
        comps       = comps,
        selectivity = selectivity
    )
    saveRDS(cache, file.path(dir_cache, "epo-bet-opal-cache.rds"), compress = "xz")
