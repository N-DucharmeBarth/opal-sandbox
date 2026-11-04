# Nicholas Ducharme-Barth
# 2026/11/04
# Tests of the length-composition treatment behind O2's knife-edge
# selectivities (06-diagnose-selectivity.r: the longline compositions are
# narrower than the model's spread of length at the ages they catch). Each is
# O2 (annual ages, mid-year growth, opal's double normals with the peak and
# two widths estimated, initial population estimated; 05-fit-opal-
# selectivity.r) with one change, through opal's inputs only:
#
#   O2-pool  tail pooling: for each fishery the bins below the 0.5% and above
#            the 99.5% points of its pooled observed composition are summed
#            into the end bins (opal's lf_minbin / lf_maxbin), in the spirit of
#            MFCL's tail compression.
#   O2-DM    Dirichlet-multinomial likelihood (lf_switch = 3) with each
#            fishery's dispersion (log_lf_tau) estimated.
#   O2-SD2   the SD of length at the oldest age (log_SD2; growth is
#            otherwise fixed) estimated.
#   O2-all   all three.
#   O1i-all  O1i (MFCL's selectivity) with all three, for comparison.
#
# OUTPUT (cached; SWPO_ALB_REFIT=true refits):
#   case-studies/swpo-alb/vignettes/data/swpo-alb-lf-tests.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Set-up
    source(file.path(this.path::this.proj(), "case-studies", "swpo-alb", "scripts", "swpo-alb-opal-setup.r"))
    out_file = file.path(dir_cache, "swpo-alb-lf-tests.rds")
    saved = if (!refit && file.exists(out_file)) readRDS(out_file)$fits else list()
    if (is.null(saved)) saved = list()
    source(file.path(proj_dir, "case-studies", "swpo-alb", "scripts", "swpo-alb-fit-helpers.r"))
    sel_fits = readRDS(file.path(dir_cache, "swpo-alb-opal-sel-fits.rds"))$fits
    stopifnot("Run 05-fit-opal-selectivity.r first" = all(c("O1i", "O2") %in% names(sel_fits)))
    fits = saved

    # Tail pooling bins from each fishery's pooled observed composition
    bins = setdiff(names(cv$lf), c("fishery", "year", "month", "ts"))
    tail_bins = function(p = 0.005) {
        lo = rep(1L, d0$n_fishery)
        hi = rep(d0$n_len, d0$n_fishery)
        for (f in unique(cv$lf$fishery)) {
            x = colSums(as.matrix(cv$lf[fishery == f, ..bins]))
            cx = cumsum(x) / sum(x)
            lo[f] = max(1L, which(cx >= p)[1])
            hi[f] = min(d0$n_len, which(cx >= 1 - p)[1])
        }
        list(minbin = lo, maxbin = hi)
    }
    tb = tail_bins()
    print(data.table(fishery = seq_len(d0$n_fishery), from = d0$len_bin_start + (tb$minbin - 1) * d0$len_bin_width,
                     to = d0$len_bin_start + tb$maxbin * d0$len_bin_width))

    # Data and parameters for a variant
    make_data = function(sel = c("opal", "mfcl"), pool = FALSE, dm = FALSE) {
        sel = match.arg(sel)
        dd = cv$data
        dd = prep_lf_data(dd, as.data.frame(cv$lf), lf_switch = if (dm) 3L else 1L,
                          lf_var_adjust = cv$meta$lf_var_adjust,
                          lf_minbin = if (pool) tb$minbin else NULL, lf_maxbin = if (pool) tb$maxbin else NULL)
        dd$lambdas = init_lambdas(d0$lambdas)
        dd$k_mult_ages = d0$k_mult_ages
        dd$growth_sd_model = d0$growth_sd_model
        if (sel == "mfcl") {
            dd$sel_fa_external = cv$selectivity$sel_fa_external
        } else {
            dd$sel_type_f = sel_type_f
            dd$sel_basis_f = rep(1L, dd$n_fishery)
        }
        dd
    }
    make_map = function(p, data, sel = TRUE, dm = FALSE, sd2 = FALSE) {
        m = a_map(p, if (sel) sel_map(group_f, sel_type_f) else NULL)
        if (dm) {
            k = rep(NA_integer_, data$n_fishery)
            k[data$lf_fishery_f] = seq_along(data$lf_fishery_f)
            m$log_lf_tau = factor(k)
        }
        if (sd2) m$log_SD2 = factor(1)
        m
    }
    variant = function(label, start, sel = "opal", pool = FALSE, dm = FALSE, sd2 = FALSE) {
        dd = make_data(sel, pool, dm)
        p = start
        if (dm) p$log_lf_tau = rep(log(50), dd$n_fishery)
        fit = fit_model(dd, p, make_map(p, dd, sel == "opal", dm, sd2), label)
        cat(sprintf("%-8s objective %.2f  log_B0 %.3f  max grad %.2g  at bound %d  SD2 %.2f\n", label,
                    fit$objective, fit$par[["log_B0"]], fit$max_gradient, fit$at_bound,
                    exp(fit$parameters$log_SD2)))
        fit
    }

#_____________________________________________________________________________
# Fits (started from O2's or O1i's estimates)
    p_o2 = sel_fits$O2$parameters
    p_o1 = sel_fits$O1i$parameters
    fits$`O2-pool` = variant("O2-pool", p_o2, pool = TRUE)
    saveRDS(list(fits = fits), out_file)
    fits$`O2-DM` = variant("O2-DM", p_o2, dm = TRUE)
    saveRDS(list(fits = fits), out_file)
    fits$`O2-SD2` = variant("O2-SD2", p_o2, sd2 = TRUE)
    saveRDS(list(fits = fits), out_file)
    fits$`O2-all` = variant("O2-all", p_o2, pool = TRUE, dm = TRUE, sd2 = TRUE)
    saveRDS(list(fits = fits), out_file)
    fits$`O1i-all` = variant("O1i-all", p_o1, sel = "mfcl", pool = TRUE, dm = TRUE, sd2 = TRUE)

#_____________________________________________________________________________
# Summary
    all_fits = c(sel_fits[c("O1i", "O2")], fits)
    lo_col = c(-2, -7, -6, -6, -10, -8)
    hi_col = c(2, 7, 3, 3, 10, 8)
    summary_dt = rbindlist(lapply(names(all_fits), function(id) {
        x = all_fits[[id]]
        ps = x$parameters$par_sel[!duplicated(group_f), , drop = FALSE]
        dn = sel_type_f[!duplicated(group_f)] == 2L
        width_bound = sum(ps[dn, 3:4] <= -6 + 1e-3 | ps[dn, 3:4] >= 3 - 1e-3)
        data.table(model = id, n_par = x$n_par, objective = x$objective, max_gradient = x$max_gradient,
                   at_bound = x$at_bound, widths_at_bound = if (grepl("^O2", id)) width_bound else NA_integer_,
                   cpue = x$lp[["cpue"]], lf = x$lp[["lf"]], sd2 = exp(x$parameters$log_SD2),
                   sb_1954 = x$sb[1], sb_2022 = x$sb[n_year], depletion_2022 = x$sb[n_year] / x$sb0[n_year])
    }))
    print(summary_dt)
    traj = rbindlist(lapply(names(all_fits), function(id) data.table(
        model = id, year = d0$years, sb = all_fits[[id]]$sb[seq_len(n_year)],
        depletion = all_fits[[id]]$sb[seq_len(n_year)] / all_fits[[id]]$sb0[seq_len(n_year)])))
    saveRDS(list(fits = fits, summary = summary_dt, trajectory = traj, tail_bins = tb), out_file)
