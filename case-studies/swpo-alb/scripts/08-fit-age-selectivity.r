# Nicholas Ducharme-Barth
# 2026/11/04
# Age-based double normals for the south-west Pacific albacore opal model
# (MFCL bridge model M3). Like O2 (05-fit-opal-selectivity.r: annual ages,
# mid-year growth, initial population estimated, MFCL's selectivity groups
# sharing parameters), but each curve is evaluated on age (opal's
# sel_basis_f = 2) rather than on length through the length-at-age key, so
# growth affects only the predicted length compositions:
#
#   O3       age-based double normals (a logistic for the longline index):
#            the peak and the two widths estimated, the plateau and final
#            selectivity held at starting curves fitted to MFCL's
#            selectivity-at-age, the initial selectivity near zero.
#   O3-SD2   O3 with the SD of length at the oldest age (log_SD2) estimated.
#   O3-all   O3 with tail pooling, the Dirichlet-multinomial and log_SD2
#            estimated (as O2-all, 07-test-lf-likelihood.r).
#   O4       hybrid: age-based curves for the longline fisheries and index,
#            length-based for the troll and driftnet fisheries (11, 12, 13,
#            15), which catch young, fast-growing fish whose compositions an
#            age-based curve at whole ages cannot follow through the year.
#   O4-all   O4 with tail pooling, the Dirichlet-multinomial and log_SD2.
#
# All start from O1i's estimates (MFCL's selectivity).
#
# OUTPUT (cached; SWPO_ALB_REFIT=true refits):
#   case-studies/swpo-alb/vignettes/data/swpo-alb-age-sel-fits.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Set-up
    source(file.path(this.path::this.proj(), "case-studies", "swpo-alb", "scripts", "swpo-alb-opal-setup.r"))
    out_file = file.path(dir_cache, "swpo-alb-age-sel-fits.rds")
    saved = if (!refit && file.exists(out_file)) readRDS(out_file)$fits else list()
    if (is.null(saved)) saved = list()
    source(file.path(proj_dir, "case-studies", "swpo-alb", "scripts", "swpo-alb-fit-helpers.r"))
    sel_fits = readRDS(file.path(dir_cache, "swpo-alb-opal-sel-fits.rds"))$fits
    lf_tests = readRDS(file.path(dir_cache, "swpo-alb-lf-tests.rds"))
    stopifnot("Run 05 and 07 first" = !is.null(sel_fits$O1i) && !is.null(lf_tests$tail_bins))
    fits = saved
    ages = seq_len(n_age)
    tb = lf_tests$tail_bins

    length_based = c(11L, 12L, 13L, 15L)     # troll and driftnet (O4)
    make_data = function(pool = FALSE, dm = FALSE, basis = rep(2L, d0$n_fishery)) {
        dd = prep_lf_data(cv$data, as.data.frame(cv$lf), lf_switch = if (dm) 3L else 1L,
                          lf_var_adjust = cv$meta$lf_var_adjust,
                          lf_minbin = if (pool) tb$minbin else NULL, lf_maxbin = if (pool) tb$maxbin else NULL)
        dd$lambdas = init_lambdas(d0$lambdas)
        dd$k_mult_ages = d0$k_mult_ages
        dd$growth_sd_model = d0$growth_sd_model
        dd$sel_type_f = sel_type_f
        dd$sel_basis_f = basis                         # 2 = age-based, 1 = length-based
        dd
    }

    # Starting curves: fitted directly to MFCL's selectivity-at-age (an identity
    # "length-at-age key" on the ages)
    p_start = sel_fits$O1i$parameters
    seeds = lapply(unique(group_f), function(cg) {
        f = which(group_f == cg)
        seed_curve(sel_type_f[f[1]], cv$selectivity$sel_fa_external[f, , drop = FALSE], diag(n_age), ages)
    })
    p_start$par_sel = matrix(0, d0$n_fishery, 6)
    for (k in seq_along(seeds)) {
        f = which(group_f == unique(group_f)[k])
        p_start$par_sel[f, ] = matrix(seeds[[k]]$par, length(f), 6, byrow = TRUE)
    }
    seed_rmse = data.table(group = unique(group_f), fisheries = vapply(unique(group_f), function(cg)
        paste(which(group_f == cg), collapse = ","), ""), rmse = vapply(seeds, `[[`, 0, "rmse"))
    print(seed_rmse)

    # Hybrid starting curves: length-based ones fitted through opal's
    # length-at-age key (mid-year growth), as in 05
    basis_hybrid = ifelse(seq_len(d0$n_fishery) %in% length_based, 1L, 2L)
    p_hybrid = p_start
    pla_a = pla_of(make_data(), p_start)
    for (cg in unique(group_f[length_based])) {
        f = which(group_f == cg)
        sq = seed_curve(sel_type_f[f[1]], cv$selectivity$sel_fa_external[f, , drop = FALSE], pla_a, d0$len_mid)
        p_hybrid$par_sel[f, ] = matrix(sq$par, length(f), 6, byrow = TRUE)
    }

    variant = function(label, pool = FALSE, dm = FALSE, sd2 = FALSE, basis = rep(2L, d0$n_fishery),
                       start = p_start) {
        dd = make_data(pool, dm, basis)
        p = start
        if (dm) p$log_lf_tau = rep(log(50), dd$n_fishery)
        m = a_map(p, sel_map(group_f, sel_type_f))
        if (dm) {
            k = rep(NA_integer_, dd$n_fishery)
            k[dd$lf_fishery_f] = seq_along(dd$lf_fishery_f)
            m$log_lf_tau = factor(k)
        }
        if (sd2) m$log_SD2 = factor(1)
        fit = fit_model(dd, p, m, label)
        cat(sprintf("%-7s objective %.2f  log_B0 %.3f  max grad %.2g  at bound %d  SD2 %.2f\n", label,
                    fit$objective, fit$par[["log_B0"]], fit$max_gradient, fit$at_bound,
                    exp(fit$parameters$log_SD2)))
        fit
    }

#_____________________________________________________________________________
# Fits
    fits$O3 = variant("O3")
    saveRDS(list(fits = fits), out_file)
    fits$`O3-SD2` = variant("O3-SD2", sd2 = TRUE)
    saveRDS(list(fits = fits), out_file)
    fits$`O3-all` = variant("O3-all", pool = TRUE, dm = TRUE, sd2 = TRUE)
    saveRDS(list(fits = fits), out_file)
    fits$O4 = variant("O4", basis = basis_hybrid, start = p_hybrid)
    saveRDS(list(fits = fits), out_file)
    fits$`O4-all` = variant("O4-all", pool = TRUE, dm = TRUE, sd2 = TRUE, basis = basis_hybrid, start = p_hybrid)

#_____________________________________________________________________________
# Summary
    all_fits = c(sel_fits[c("O1i", "O2")], lf_tests$fits["O2-all"], fits)
    summary_dt = rbindlist(lapply(names(all_fits), function(id) {
        x = all_fits[[id]]
        ps = x$parameters$par_sel[!duplicated(group_f), , drop = FALSE]
        dn = sel_type_f[!duplicated(group_f)] == 2L
        data.table(model = id, n_par = x$n_par, objective = x$objective, max_gradient = x$max_gradient,
                   at_bound = x$at_bound,
                   widths_at_bound = if (id == "O1i") NA_integer_ else
                       sum(ps[dn, 3:4] <= -6 + 1e-3 | ps[dn, 3:4] >= 3 - 1e-3),
                   cpue = x$lp[["cpue"]], lf = x$lp[["lf"]], sd2 = exp(x$parameters$log_SD2),
                   sb_1954 = x$sb[1], sb_2022 = x$sb[n_year], depletion_2022 = x$sb[n_year] / x$sb0[n_year])
    }))
    print(summary_dt)
    traj = rbindlist(lapply(names(all_fits), function(id) data.table(
        model = id, year = d0$years, sb = all_fits[[id]]$sb[seq_len(n_year)],
        depletion = all_fits[[id]]$sb[seq_len(n_year)] / all_fits[[id]]$sb0[seq_len(n_year)])))
    sel_age = rbindlist(c(
        list(data.table(model = "MFCL", fishery = rep(seq_len(d0$n_fishery), each = n_age),
                        age = rep(ages, d0$n_fishery), sel = as.vector(t(cv$selectivity$sel_fa_external)))),
        lapply(names(fits), function(id) data.table(model = id, fishery = rep(seq_len(d0$n_fishery), each = n_age),
                                                   age = rep(ages, d0$n_fishery),
                                                   sel = as.vector(t(fits[[id]]$sel_fya))))))
    sel_age[, name := fishery_names[fishery]]
    saveRDS(list(fits = fits, summary = summary_dt, trajectory = traj, sel_age = sel_age, seed_rmse = seed_rmse),
            out_file)
