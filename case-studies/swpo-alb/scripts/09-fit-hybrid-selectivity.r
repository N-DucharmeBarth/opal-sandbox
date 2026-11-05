# Nicholas Ducharme-Barth
# 2026/11/05
# The hybrid selectivity model (O4, 08-fit-age-selectivity.r: age-based
# double normals for the longline fisheries and index, length-based for the
# troll and driftnet fisheries 11, 12, 13 and 15) refined:
#
#   O4-fix    O4 with each selectivity parameter that O4 left at a bound fixed
#             there, and the rest refitted from O4's estimates.
#   O4q       the quarterly-age version of O4 (quarterly ages and steps,
#             recruitment annual, initial population estimated with the plus
#             group tied; as O1q, 05-fit-opal-selectivity.r). The age-based
#             curves run over the 48 quarterly ages; started from O1q's
#             estimates and curves fitted to MFCL's selectivity by quarterly
#             age.
#   O4q-fix   O4q with its selectivity parameters at bounds fixed there.
#
# OUTPUT (cached; SWPO_ALB_REFIT=true refits):
#   case-studies/swpo-alb/vignettes/data/swpo-alb-hybrid-fits.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Set-up
    source(file.path(this.path::this.proj(), "case-studies", "swpo-alb", "scripts", "swpo-alb-opal-setup.r"))
    out_file = file.path(dir_cache, "swpo-alb-hybrid-fits.rds")
    saved = if (!refit && file.exists(out_file)) readRDS(out_file)$fits else list()
    if (is.null(saved)) saved = list()
    source(file.path(proj_dir, "case-studies", "swpo-alb", "scripts", "swpo-alb-fit-helpers.r"))
    sel_fits = readRDS(file.path(dir_cache, "swpo-alb-opal-sel-fits.rds"))$fits
    age_fits = readRDS(file.path(dir_cache, "swpo-alb-age-sel-fits.rds"))$fits
    stopifnot("Run 05 and 08 first" = !is.null(sel_fits$O1q) && !is.null(age_fits$O4))
    fits = saved
    length_based = c(11L, 12L, 13L, 15L)
    basis_hybrid = ifelse(seq_len(d0$n_fishery) %in% length_based, 1L, 2L)
    lo_col = c(-2, -7, -6, -6, -10, -8)
    hi_col = c(2, 7, 3, 3, 10, 8)

    # Selectivity map with the parameters a fit left at a bound fixed there
    fix_bounded = function(sel, par_sel, tol = 1e-3) {
        m = matrix(as.integer(as.character(sel)), nrow(par_sel))
        col = attr(sel, "col")
        fixed = integer()
        for (j in sort(unique(m[!is.na(m)]))) {
            v = par_sel[which(m == j)[1]]
            if (v <= lo_col[col[j]] + tol || v >= hi_col[col[j]] - tol) {
                m[m == j] = NA
                fixed = c(fixed, j)
            }
        }
        keep = sort(unique(m[!is.na(m)]))
        m[] = match(m, keep)
        out = factor(as.vector(m))
        attr(out, "col") = col[keep]
        attr(out, "n_fixed") = length(fixed)
        out
    }
    report_summary = function(label, fit) {
        cat(sprintf("%-8s objective %.2f  log_B0 %.3f  max grad %.2g  at bound %d\n", label,
                    fit$objective, fit$par[["log_B0"]], fit$max_gradient, fit$at_bound))
    }

#_____________________________________________________________________________
# O4-fix: annual
    d4 = prep_lf_data(cv$data, as.data.frame(cv$lf), lf_switch = 1L, lf_var_adjust = cv$meta$lf_var_adjust)
    d4$lambdas = init_lambdas(d0$lambdas)
    d4$k_mult_ages = d0$k_mult_ages
    d4$growth_sd_model = d0$growth_sd_model
    d4$sel_type_f = sel_type_f
    d4$sel_basis_f = basis_hybrid
    p4 = age_fits$O4$parameters
    m4 = a_map(p4, fix_bounded(sel_map(group_f, sel_type_f), p4$par_sel))
    cat("O4-fix: selectivity parameters fixed at bounds:", attr(m4$par_sel, "n_fixed"), "\n")
    fits$`O4-fix` = fit_model(d4, p4, m4, "O4-fix")
    fits$`O4-fix`$n_sel_fixed = attr(m4$par_sel, "n_fixed")
    report_summary("O4-fix", fits$`O4-fix`)
    saveRDS(list(fits = fits), out_file)

#_____________________________________________________________________________
# O4q: quarterly ages (as O1q in 05)
    qm = mfcl_quarterly_ages(cv, laa_month)
    dq = prep_lf_data(qm$data, as.data.frame(cv$lf_season), lf_switch = 1L,
                      lf_var_adjust = cv$meta$lf_var_adjust)
    dq$lambdas = init_lambdas(d0$lambdas)
    report_q = function(p) MakeADFun(cmb(opal_model, dq), p, map = fixed_all(p), silent = TRUE)$report()
    pq0 = mfcl_match_recruitment_quarterly(qm, mf, report_q, cv = cv)
    q_plus = max(which(qm$age_season == 1L))
    dq$init_bias_adj_a[q_plus] = -2 * (-log(1 - exp(-d0$M[n_age]))) / exp(pq0$log_sigma_r)^2
    est_q = qm$rdev_steps[d0$years[qm$rec_year_of_step] <= 2020]
    q_map = function(p, sel = NULL) {
        m = fixed_all(p)
        m$log_B0 = factor(1)
        m$log_cpue_q = factor(seq_along(p$log_cpue_q))
        i = rep(NA_integer_, length(p$rdev_y))
        i[est_q] = seq_along(est_q)
        m$rdev_y = factor(i)
        e = which(qm$age_season == 1L)
        j = rep(NA_integer_, length(p$init_rdev_a))
        j[e] = c(1L, rep(2L, length(e) - 1L))
        m$init_rdev_a = factor(j)
        if (!is.null(sel)) m$par_sel = sel
        m
    }
    dq4 = dq
    dq4$sel_fa_external = NULL
    dq4$sel_type_f = sel_type_f
    dq4$sel_basis_f = basis_hybrid
    pq4 = sel_fits$O1q$parameters
    pla_q = pla_of(dq4, pq4)
    q_ages = seq_len(dq4$n_age)
    for (cg in unique(group_f)) {
        f = which(group_f == cg)
        target = qm$data$sel_fa_external[f, , drop = FALSE]
        sq = if (basis_hybrid[f[1]] == 2L) seed_curve(sel_type_f[f[1]], target, diag(dq4$n_age), q_ages) else
            seed_curve(sel_type_f[f[1]], target, pla_q, dq4$len_mid)
        pq4$par_sel[f, ] = matrix(sq$par, length(f), 6, byrow = TRUE)
    }
    fits$O4q = fit_model(dq4, pq4, q_map(pq4, sel_map(group_f, sel_type_f)), "O4q")
    report_summary("O4q", fits$O4q)
    saveRDS(list(fits = fits), out_file)

    pq4f = fits$O4q$parameters
    mq4f = q_map(pq4f, fix_bounded(sel_map(group_f, sel_type_f), pq4f$par_sel))
    cat("O4q-fix: selectivity parameters fixed at bounds:", attr(mq4f$par_sel, "n_fixed"), "\n")
    fits$`O4q-fix` = fit_model(dq4, pq4f, mq4f, "O4q-fix")
    fits$`O4q-fix`$n_sel_fixed = attr(mq4f$par_sel, "n_fixed")
    report_summary("O4q-fix", fits$`O4q-fix`)

#_____________________________________________________________________________
# Summary
    q1 = which(qm$step_season == 1L)
    all_fits = c(sel_fits[c("O1i", "O1q")], age_fits["O4"], fits)
    steps_of = function(id) if (grepl("q", id)) q1 else seq_len(n_year)
    summary_dt = rbindlist(lapply(names(all_fits), function(id) {
        x = all_fits[[id]]
        ps = x$parameters$par_sel[!duplicated(group_f), , drop = FALSE]
        dn = sel_type_f[!duplicated(group_f)] == 2L
        st = steps_of(id)
        data.table(model = id, n_par = x$n_par, objective = x$objective, max_gradient = x$max_gradient,
                   at_bound = x$at_bound, sel_fixed = if (is.null(x$n_sel_fixed)) NA_integer_ else x$n_sel_fixed,
                   widths_at_bound = if (id %in% c("O1i", "O1q")) NA_integer_ else
                       sum(ps[dn, 3:4] <= -6 + 1e-3 | ps[dn, 3:4] >= 3 - 1e-3),
                   cpue = x$lp[["cpue"]], lf = x$lp[["lf"]],
                   sb_1954 = x$sb[st][1], sb_2022 = x$sb[st][n_year],
                   depletion_2022 = x$sb[st][n_year] / x$sb0[st][n_year])
    }))
    print(summary_dt)
    traj = rbindlist(lapply(names(all_fits), function(id) {
        st = steps_of(id)
        data.table(model = id, year = d0$years, sb = all_fits[[id]]$sb[st],
                   depletion = all_fits[[id]]$sb[st] / all_fits[[id]]$sb0[st])
    }))
    saveRDS(list(fits = fits, summary = summary_dt, trajectory = traj), out_file)
