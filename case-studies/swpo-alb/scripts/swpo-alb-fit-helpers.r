# Nicholas Ducharme-Barth
# 2026/11/04
# Fit helpers for the south-west Pacific albacore opal fits with opal's own
# selectivity (sourced by 07-test-lf-likelihood.r; the same definitions as in
# 05-fit-opal-selectivity.r). Needs the objects of swpo-alb-opal-setup.r and
# `saved` (a list of cached fits, possibly empty).
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

    rdev_years = which(d0$years[-1] <= 2020)
    group_f = c(1, 2, 1, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)   # MFCL selectivity groups (fish flag 24)
    sel_type_f = ifelse(seq_len(d0$n_fishery) == 14L, 1L, 2L)       # longline index non-decreasing
    init_lambdas = function(lam) { lam$init_rec = 0; lam }

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

