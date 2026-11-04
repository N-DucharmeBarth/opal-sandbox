# Nicholas Ducharme-Barth
# 2026/11/04
# Fits of the south-west Pacific albacore opal models beyond O1 (MFCL bridge
# model M3; see 03 and 04). All estimate the initial population, as MFCL does
# (it scales its initial population apart from its stock-recruitment
# relationship): one deviation for the first year's recruits (age 1) and one
# shared by ages 2 and over (the level of the initial age structure; in the
# quarterly models the plus group keeps its equilibrium offset through
# init_bias_adj_a), unpenalised (lambda init_rec = 0).
#
#   O1i       annual ages, mid-year growth, MFCL's selectivity fixed (O1-mid
#             with the initial population estimated).
#   O1q       quarterly ages and time steps (mfcl_quarterly_ages(); recruitment
#             annual), MFCL's seasonal selectivity-at-age fixed; B0, the
#             first-quarter recruitment deviations (1955-2020) and the two
#             index catchabilities estimated.
#   O2        annual ages, four seasons, opal's growth at MFCL's mid-year
#             lengths (O1-mid), opal's own length-based selectivity: double
#             normal (logistic for the longline index, non-decreasing in MFCL),
#             fisheries in the same MFCL selectivity group sharing parameters.
#   O2-split  as O2 with fisheries 1-7 split into one fishery per season
#             (mfcl_split_seasonal_fisheries()); seasonal fisheries share
#             parameters within MFCL's groups and seasons.
#   O2q       O1q with opal's own length-based selectivity (one curve per MFCL
#             group, whose selectivity-at-age varies by quarter through growth).
#
# The double normals estimate three parameters: the peak and the ascending
# and descending widths (columns 1, 3 and 4; the peak bounded to +/-2 SD of
# the length bins, the widths to [-6, 3]). The plateau and final selectivity
# (columns 2 and 6) are held at their starting values and the initial
# selectivity (column 5) at -8 (near zero). The logistic estimates both its
# parameters. Starting curves (all of columns 1-4 and 6)
# are fitted to MFCL's selectivity-at-age through opal's own
# probability-of-length-at-age. Everything else (growth, steepness, sigma_R,
# M, biology) is held as in O1.
#
# OUTPUTS (cached; SWPO_ALB_REFIT=true refits):
#   case-studies/swpo-alb/vignettes/data/swpo-alb-opal-sel-fits.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Set-up (conversion, growth, O0, replays)
    source(file.path(this.path::this.proj(), "case-studies", "swpo-alb", "scripts", "swpo-alb-opal-setup.r"))
    out_file = file.path(dir_cache, "swpo-alb-opal-sel-fits.rds")
    saved = if (!refit && file.exists(out_file)) readRDS(out_file)$fits else list()
    if (is.null(saved)) saved = list()

    rdev_years = which(d0$years[-1] <= 2020)
    group_f = c(1, 2, 1, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)   # MFCL selectivity groups (fish flag 24)
    sel_type_f = ifelse(seq_len(d0$n_fishery) == 14L, 1L, 2L)       # longline index non-decreasing

    # Length-at-age, SD and probability of length-at-age from opal's growth inputs
    pla_of = function(data, p) {
        mu = get_growth(data$n_age, A1 = data$A1, A2 = data$A2, L1 = exp(p$log_L1), L2 = exp(p$log_L2),
                        log_k = p$log_k, log_k_mult = p$log_k_mult, k_mult_ages = data$k_mult_ages)
        sd = get_sd_at_age(mu, mu[data$A1], mu[data$A2], sd_model = "sd",
                           log_SD1 = p$log_SD1, log_SD2 = p$log_SD2, A1 = data$A1, A2 = data$A2)
        get_pla(data$len_lower, data$len_upper, mu, as.numeric(sd))
    }
    # Starting curve: fit to target selectivity-at-age (rows = targets sharing
    # the curve), from a grid of starts (peak across the length range, widths
    # from narrow to wide), then refined
    seed_curve = function(type, targets, pla, len_mid) {
        targets = matrix(targets, ncol = ncol(pla))
        full = function(x) if (type == 1L) c(x, 0, 0, 0, 0) else c(x[1:4], -8, x[5])
        obj = function(x) {
            sa = c(t(pla) %*% sel_length(len_mid, full(x), type))
            sum(sweep(targets, 2, sa)^2)
        }
        starts = if (type == 1L) as.matrix(expand.grid(a = seq(-1.5, 1.5, 0.5), b = c(-2, -1, 0))) else
            as.matrix(expand.grid(a = seq(-1.5, 1.5, 0.375), b = c(-4, -1), c = c(-4, -2.5, -1), d = c(-4, -2.5, -1, 1),
                                  f = c(-6, 0, 4)))
        val = apply(starts, 1, obj)
        best = order(val)[1:5]
        lo = if (type == 1L) c(-2, -7) else c(-2, -7, -6, -6, -8)
        hi = if (type == 1L) c(2, 7) else c(2, 7, 3, 3, 8)
        fits = lapply(best, function(i) optim(pmin(pmax(starts[i, ], lo), hi), obj, method = "L-BFGS-B",
                                              lower = lo, upper = hi, control = list(maxit = 5000, factr = 1e2)))
        o = fits[[which.min(vapply(fits, `[[`, 0, "value"))]]
        list(par = full(o$par), rmse = sqrt(o$value / length(targets)))
    }
    # Selectivity map: estimated columns per type, shared within curve groups
    sel_map = function(curve_f, type_f) {
        m = matrix(NA_integer_, length(curve_f), 6)
        k = 0L
        for (cg in unique(curve_f)) {
            f = which(curve_f == cg)
            cols = if (type_f[f[1]] == 1L) 1:2 else c(1, 3, 4)
            m[f, cols] = rep(k + seq_along(cols), each = length(f))
            k = k + length(cols)
        }
        out = factor(as.vector(m))
        attr(out, "col") = vapply(seq_len(k), function(j) which(m == j, arr.ind = TRUE)[1, 2], 1L)
        out
    }

    fixed_all = function(p) lapply(p, function(z) factor(rep(NA, length(z))))
    fit_model = function(data, start, map, label) {
        if (!is.null(saved[[label]])) return(saved[[label]])
        t0 = proc.time()[["elapsed"]]
        ob = MakeADFun(cmb(opal_model, data), start, map = map, silent = TRUE)
        nm = names(ob$par)
        lower = ifelse(nm == "log_B0", 0, ifelse(nm == "par_sel", -7, ifelse(nm == "log_cpue_q", -15, -10)))
        upper = ifelse(nm == "log_B0", 22, ifelse(nm == "par_sel", 7, ifelse(nm == "log_cpue_q", 15, 10)))
        sc = attr(map$par_sel, "col")
        if (!is.null(sc)) {
            i = which(nm == "par_sel")
            lower[i] = c(-2, -7, -6, -6, -10, -8)[sc]
            upper[i] = c(2, 7, 3, 3, 10, 8)[sc]
        }
        par = pmin(pmax(ob$par, lower + 1e-6), upper - 1e-6)
        for (pass in 1:3) {
            op = nlminb(par, ob$fn, ob$gr, lower = lower, upper = upper,
                        control = list(eval.max = 1e5, iter.max = 1e5, rel.tol = 1e-12))
            par = op$par
        }
        r = ob$report(op$par)
        g = ob$gr(op$par)
        list(label = label, objective = op$objective, convergence = op$convergence,
             max_gradient = max(abs(g)), n_par = length(op$par), par = op$par,
             parameters = ob$env$parList(op$par), at_bound = sum(op$par <= lower + 1e-4 | op$par >= upper - 1e-4),
             lp = c(cpue = sum(r$lp_cpue), lf = sum(r$lp_lf), rec = sum(r$lp_rec), penalty = sum(r$lp_penalty)),
             sb = r$spawning_biomass_y, sb0 = r$spawning_biomass0_y, recruits = r$number_ysa[, 1, 1],
             cpue_pred = r$cpue_pred, sel_fya = r$sel_fya[, 1, ], seconds = proc.time()[["elapsed"]] - t0)
    }
    fits = saved
    # Initial population: age 1 and a shared level for the older ages
    init_lambdas = function(lam) { lam$init_rec = 0; lam }

#_____________________________________________________________________________
# O1q: quarterly ages, MFCL's selectivity fixed
    qm = mfcl_quarterly_ages(cv, laa_month)
    dq = prep_lf_data(qm$data, as.data.frame(cv$lf_season), lf_switch = 1L,
                      lf_var_adjust = cv$meta$lf_var_adjust)
    dq$lambdas = init_lambdas(d0$lambdas)
    report_q = function(p) MakeADFun(cmb(opal_model, dq), p, map = fixed_all(p), silent = TRUE)$report()
    pq0 = mfcl_match_recruitment_quarterly(qm, mf, report_q, cv = cv)
    pq_start = pq0
    pq_start$init_rdev_a[qm$age_season == 1L] = 0
    # The plus group's class (quarterly age 45) holds the 12+ equilibrium, a
    # factor 1 / (1 - exp(-M)) above one cohort. opal's map cannot tie it to
    # the shared level with an offset, so the offset is supplied through
    # init_bias_adj_a (opal multiplies initial numbers at age a by
    # exp(init_rdev_a[a] - init_bias_adj_a[a] sigma_R^2 / 2)): a fixed input,
    # not a model change. The initial-recruitment penalty does not use it.
    q_plus = max(which(qm$age_season == 1L))
    plus_offset = -log(1 - exp(-d0$M[n_age]))
    dq$init_bias_adj_a[q_plus] = -2 * plus_offset / exp(pq0$log_sigma_r)^2
    pq_start$rdev_y[qm$rdev_steps] = 0
    est_q = qm$rdev_steps[d0$years[qm$rec_year_of_step] <= 2020]
    q_map = function(p, sel = NULL) {
        m = fixed_all(p)
        m$log_B0 = factor(1)
        m$log_cpue_q = factor(seq_along(p$log_cpue_q))
        i = rep(NA_integer_, length(p$rdev_y))
        i[est_q] = seq_along(est_q)
        m$rdev_y = factor(i)
        # initial numbers: quarterly ages of existing cohorts (first quarter of
        # each annual class); age 1 on its own, the older classes (the plus
        # group with its equilibrium offset) sharing one level
        e = which(qm$age_season == 1L)
        j = rep(NA_integer_, length(p$init_rdev_a))
        j[e] = c(1L, rep(2L, length(e) - 1L))
        m$init_rdev_a = factor(j)
        if (!is.null(sel)) m$par_sel = sel
        m
    }
    fits$O1q = fit_model(dq, pq_start, q_map(pq_start), "O1q")
    cat("O1q:", fits$O1q$objective, "log_B0", fits$O1q$par[["log_B0"]], "(quarterly scale; MFCL-matched",
        pq0$log_B0, ") max grad", fits$O1q$max_gradient, "\n")
    saveRDS(list(fits = fits), out_file)

#_____________________________________________________________________________
# O2: annual ages, mid-year growth, opal's selectivity
    annual_start = function(p) {
        p = with_growth(p, g_mid$par)
        p$init_rdev_a[] = 0
        p$rdev_y[] = 0
        p
    }
    a_map = function(p, sel) {
        m = fixed_all(p)
        m$log_B0 = factor(1)
        m$log_cpue_q = factor(seq_along(p$log_cpue_q))
        i = rep(NA_integer_, length(p$rdev_y))
        i[rdev_years] = seq_along(rdev_years)
        m$rdev_y = factor(i)
        m$init_rdev_a = factor(c(1L, rep(2L, length(p$init_rdev_a) - 1L)))
        if (!is.null(sel)) m$par_sel = sel
        m
    }

    # O1i: O1-mid (MFCL's selectivity) with the initial population estimated
    d1 = d0
    d1$lambdas = init_lambdas(d0$lambdas)
    p1i = annual_start(p0)
    fits$O1i = fit_model(d1, p1i, a_map(p1i, NULL), "O1i")
    cat("O1i:", fits$O1i$objective, "log_B0", fits$O1i$par[["log_B0"]], "max grad", fits$O1i$max_gradient, "\n")
    saveRDS(list(fits = fits), out_file)

    d2 = d0
    d2$lambdas = init_lambdas(d0$lambdas)
    d2$sel_fa_external = NULL
    d2$sel_type_f = sel_type_f
    d2$sel_basis_f = rep(1L, d2$n_fishery)
    p2 = annual_start(p0)
    pla_a = pla_of(d2, p2)
    seeds_a = lapply(unique(group_f), function(cg) {
        f = which(group_f == cg)
        seed_curve(sel_type_f[f[1]], cv$selectivity$sel_fa_external[f, , drop = FALSE], pla_a, d2$len_mid)
    })
    for (k in seq_along(seeds_a)) p2$par_sel[group_f == unique(group_f)[k], ] =
        matrix(seeds_a[[k]]$par, sum(group_f == unique(group_f)[k]), 6, byrow = TRUE)
    fits$O2 = fit_model(d2, p2, a_map(p2, sel_map(group_f, sel_type_f)), "O2")
    cat("O2:", fits$O2$objective, "log_B0", fits$O2$par[["log_B0"]], "max grad", fits$O2$max_gradient, "\n")
    saveRDS(list(fits = fits), out_file)

#_____________________________________________________________________________
# O2-split: one fishery per season for fisheries 1-7
    cs = mfcl_split_seasonal_fisheries(cv)
    fm = cs$meta$fishery_map
    ds = prep_lf_data(cs$data, as.data.frame(cs$lf), lf_switch = 1L, lf_var_adjust = cs$meta$lf_var_adjust)
    ds$lambdas = init_lambdas(d0$lambdas)
    ds$k_mult_ages = d0$k_mult_ages
    ds$growth_sd_model = d0$growth_sd_model
    ds$sel_type_f = sel_type_f[fm$mfcl]
    ds$sel_basis_f = rep(1L, ds$n_fishery)
    curve_s = paste(group_f[fm$mfcl], ifelse(is.na(fm$season), 0L, fm$season))
    ps = annual_start(p0)
    ps$log_lf_tau = cs$parameters$log_lf_tau
    ps$log_init_F_f = cs$parameters$log_init_F_f
    ps$par_sel = matrix(0, ds$n_fishery, 6)
    for (cg in unique(curve_s)) {
        f = which(curve_s == cg)
        sd_seed = seed_curve(ds$sel_type_f[f[1]], cs$selectivity$sel_fa_external[f, , drop = FALSE], pla_a, ds$len_mid)
        ps$par_sel[f, ] = matrix(sd_seed$par, length(f), 6, byrow = TRUE)
    }
    fits$`O2-split` = fit_model(ds, ps, a_map(ps, sel_map(curve_s, ds$sel_type_f)), "O2-split")
    cat("O2-split:", fits$`O2-split`$objective, "log_B0", fits$`O2-split`$par[["log_B0"]],
        "max grad", fits$`O2-split`$max_gradient, "\n")
    saveRDS(list(fits = fits), out_file)

#_____________________________________________________________________________
# O2q: quarterly ages, opal's selectivity
    dq2 = dq
    dq2$sel_fa_external = NULL
    dq2$sel_type_f = sel_type_f
    dq2$sel_basis_f = rep(1L, dq2$n_fishery)
    pq2 = fits$O1q$parameters
    pla_q = pla_of(dq2, pq2)
    for (cg in unique(group_f)) {
        f = which(group_f == cg)
        sq = seed_curve(sel_type_f[f[1]], qm$data$sel_fa_external[f, , drop = FALSE], pla_q, dq2$len_mid)
        pq2$par_sel[f, ] = matrix(sq$par, length(f), 6, byrow = TRUE)
    }
    fits$O2q = fit_model(dq2, pq2, q_map(pq2, sel_map(group_f, sel_type_f)), "O2q")
    cat("O2q:", fits$O2q$objective, "log_B0", fits$O2q$par[["log_B0"]], "max grad", fits$O2q$max_gradient, "\n")

#_____________________________________________________________________________
# Summary and cache
    q1 = which(qm$step_season == 1L)
    yrs = d0$years
    traj = rbindlist(lapply(names(fits), function(id) {
        x = fits[[id]]
        steps = if (id %in% c("O1q", "O2q")) q1 else seq_len(n_year)
        data.table(model = id, year = yrs, sb = x$sb[steps], dynamic_depletion = x$sb[steps] / x$sb0[steps],
                   recruits = x$recruits[steps])
    }))
    summary_dt = rbindlist(lapply(fits, function(x) data.table(
        model = x$label, n_par = x$n_par, objective = x$objective, convergence = x$convergence,
        max_gradient = x$max_gradient, at_bound = x$at_bound, cpue = x$lp[["cpue"]], lf = x$lp[["lf"]],
        rec = x$lp[["rec"]], log_B0 = x$par[["log_B0"]], minutes = x$seconds / 60)))
    summary_dt = merge(summary_dt, traj[year == max(year), .(model, sb_end = sb, depletion_end = dynamic_depletion)],
                       by = "model")
    print(summary_dt)
    saveRDS(c(list(fits = fits), list(summary = summary_dt, trajectory = traj, group_f = group_f,
                         sel_type_f = sel_type_f, split_map = fm, curve_split = curve_s,
                         len_mid = d2$len_mid, age_class_q = qm$age_class, age_season_q = qm$age_season,
                         mfcl_sel_fsa = mf$sel_fsa)), out_file)
