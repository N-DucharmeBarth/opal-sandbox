# An opal model with quarterly age classes for an MFCL model with annual age
# classes and sub-annual fishing (e.g. case-studies/swpo-alb).
#
# opal's seasons share one length-at-age and one selectivity per year, and
# its length compositions are annual. MFCL predicts each length sample at its
# incident's date, from lengths that grow through the year, with selectivity
# that can differ by season. Making opal's time step a quarter and its age
# classes quarters (n_season = 1, 4 x n_age ages, 4 x n_year steps) gives
# each quarter its own length-at-age, weight-at-age, selectivity and length
# compositions, through opal's inputs alone:
#   * recruitment only at the start of MFCL's model year: the recruitment
#     deviations of the other three quarters are fixed at -30 and the initial
#     numbers of the quarterly ages that do not exist (cohorts recruit once a
#     year) are set to ~0 through init_rdev_a;
#   * quarterly age q = 4 (a - 1) + s holds annual age class a in season s, so
#     seasonal selectivity is sel[f, q] = MFCL's selectivity of season s at
#     age a (exact, except the plus group, which is present in every season
#     and takes season n_season's selectivity);
#   * M per quarter is MFCL's M-at-age / 4; spawning potential at every
#     quarterly age of class a is MFCL's start-of-year value (MFCL's adult
#     biomass is at the start of the year);
#   * length-at-age is MFCL's own (from length.fit) at each quarter's fishing
#     month; opal's VB with age-specific K multipliers is fitted to it.
#
# What does not carry over: opal's stock-recruitment relationship treats every
# step as a recruitment, so with recruitment in one quarter of four B0 and
# steepness lose their usual meaning (unfished spawning biomass is about B0 /
# 4 and the curve is evaluated at a quarter of the usual relative stock size).
# Use dynamic depletion (against opal's unfished projection) and treat B0 as
# a scale parameter. Catch is still removed at the start of each quarter.

#' Read MFCL's length.fit: per sample, the fishery, year index, model month,
#' sample size and the mean length of each age's predicted component.
mfcl_read_length_fit = function(path) {
    x = mfcl_read_lines(path)
    tok = strsplit(trimws(x), "[[:space:]]+")
    nf = lengths(tok)
    nf[!nzchar(trimws(x))] = 0
    hdr_lines = grep("^# fishery [0-9]+$", x)
    fish_of = findInterval(seq_along(x), hdr_lines)
    n_bins = as.numeric(tok[[3]][1])
    first = as.numeric(tok[[3]][2])
    width = as.numeric(tok[[3]][3])
    n_age = as.numeric(tok[[6]][1])
    mid = first + width * (seq_len(n_bins) - 0.5)
    i_rec = which(nf == 3 & !startsWith(x, "#") & seq_along(x) > hdr_lines[1] &
                  nf[pmin(seq_along(x) + 3L, length(x))] == n_age)
    num = function(i) as.numeric(tok[[i]])
    out = lapply(i_rec, function(i) {
        L = vapply(seq_len(n_age), function(a) {
            v = num(i + 7L + a)
            if (sum(v) > 1e-8) sum(v * mid) / sum(v) else NA_real_
        }, numeric(1))
        c(fishery = fish_of[i], year_index = num(i)[1], model_month = num(i)[2],
          n = num(i + 2L), L)
    })
    out = do.call(rbind, out)
    colnames(out)[-(1:4)] = paste0("age_", seq_len(n_age))
    data.table::as.data.table(out)
}

#' MFCL's mean length-at-age by model month [12, n_age]: for each age, a line
#' in model month through the mean predicted lengths of length.fit's samples.
mfcl_length_at_age_by_month = function(lfit) {
    ages = grep("^age_", names(lfit), value = TRUE)
    vapply(ages, function(a) {
        z = lfit[is.finite(get(a)), .(L = mean(get(a))), by = model_month]
        stats::predict(stats::lm(L ~ model_month, data = z), data.frame(model_month = 1:12))
    }, numeric(12))
}

#' Turn a converted annual-age MFCL model (mfcl_to_opal()) into a
#' quarterly-age opal model.
#'
#' @param cv output of mfcl_to_opal() with n_season = 4.
#' @param laa_month MFCL length-at-age by model month (mfcl_length_at_age_by_month()).
#' @param fish_months model month of each season's fishing (MFCL's mid-quarter incidents).
#' @param n_k_mult quarterly ages given their own K multiplier in opal's growth.
#' @return list(data, parameters, rdev_steps, map_fixed, growth_check, steps)
mfcl_quarterly_ages = function(cv, laa_month, fish_months = c(2, 5, 8, 11), n_k_mult = 16L) {
    d = cv$data
    ns = d$n_season
    na = d$n_age
    nq = na * ns
    nt = d$n_year * ns
    a_of = (seq_len(nq) - 1L) %/% ns + 1L
    s_of = (seq_len(nq) - 1L) %% ns + 1L
    step_year = (seq_len(nt) - 1L) %/% ns + 1L
    step_season = (seq_len(nt) - 1L) %% ns + 1L

    catch = array(0, c(nt, 1, d$n_fishery))
    for (f in seq_len(d$n_fishery)) catch[, 1, f] = as.vector(t(d$catch_obs_ysf[, , f]))
    cpue = cv$data$cpue_data          # ts is already (y - 1) * n_season + s
    L_q = laa_month[cbind(fish_months[s_of], a_of)]
    w_a = d$weight
    w_q = w_a[a_of] * (L_q / laa_month[1, a_of])^d$lw_b     # weight at the fishing month

    sel_fsa = cv$mfcl$sel_fsa
    sel_q = t(vapply(seq_len(d$n_fishery), function(f) sel_fsa[cbind(f, s_of, a_of)], numeric(nq)))

    dq = list(
        age_a = seq_len(nq), n_age = nq, min_age = 1L, max_age = nq,
        n_season = 1L, n_fishery = d$n_fishery, n_year = nt,
        first_yr = 1L, last_yr = nt, years = seq_len(nt), first_yr_catch = 1L,
        len_bin_start = d$len_bin_start, len_bin_width = d$len_bin_width, n_len = d$n_len,
        catch_units_f = d$catch_units_f, catch_obs_ysf = catch,
        cpue_switch = 1L, n_index = d$n_index, cpue_data = cpue,
        lw_a = d$lw_a, lw_b = d$lw_b,
        weight = w_q, maturity = d$maturity[a_of], fecundity = w_q, sex_ratio = rep(1, nq),
        spawning_potential = d$spawning_potential[a_of], M = d$M[a_of] / ns,
        A1 = 1L, A2 = as.integer(nq), sel_type_f = d$sel_type_f,
        removal_switch_f = d$removal_switch_f,
        bias_adj_y = rep(0, nt), init_bias_adj_a = rep(0, nq),
        sel_fa_external = sel_q
    )

    # Growth: opal VB with K multipliers at quarterly ages 1..n_k_mult, fitted to L_q
    k_ages = seq_len(n_k_mult)
    gfun = function(x) get_growth(nq, A1 = 1, A2 = nq, L1 = exp(x[1]), L2 = exp(x[2]),
                                  log_k = x[3], log_k_mult = x[-(1:3)], k_mult_ages = k_ages)
    g = stats::optim(c(log(L_q[1]), log(L_q[nq]), log(0.3 / ns), rep(0, n_k_mult)),
                     function(x) sum((gfun(x) - L_q)^2), method = "BFGS",
                     control = list(maxit = 20000, reltol = 1e-14))
    mu_q = gfun(g$par)
    sd_a = cv$mfcl$sd_a
    dq$k_mult_ages = k_ages
    dq$growth_sd_model = "sd"

    p = cv$parameters
    pq = list(
        log_B0 = p$log_B0, log_h = p$log_h, log_sigma_r = p$log_sigma_r,
        log_cpue_q = p$log_cpue_q, cpue_creep = p$cpue_creep,
        log_cpue_tau = p$log_cpue_tau, log_cpue_omega = p$log_cpue_omega,
        log_lf_tau = p$log_lf_tau,
        rdev_y = rep(0, nt), init_rdev_a = rep(0, nq),
        log_init_F_f = p$log_init_F_f, par_sel = p$par_sel,
        log_L1 = g$par[1], log_L2 = g$par[2], log_k = g$par[3], log_k_mult = g$par[-(1:3)],
        log_CV1 = log(sd_a[1] / mu_q[1]), log_CV2 = log(sd_a[na] / mu_q[nq]),
        log_SD1 = log(sd_a[1]), log_SD2 = log(sd_a[na])
    )
    # rdev_y[t] sets recruitment at the start of step t + 1: recruitment steps
    # are the first quarter of each model year
    rec_step = which(step_season == 1L)
    rdev_steps = rec_step[-1] - 1L
    pq$rdev_y[-rdev_steps] = -30

    list(data = dq, parameters = pq, rdev_steps = rdev_steps,
         rec_year_of_step = step_year[rdev_steps + 1L],
         age_class = a_of, age_season = s_of, step_year = step_year, step_season = step_season,
         growth_check = data.table::data.table(q = seq_len(nq), age = a_of, season = s_of,
                                               mfcl_length = L_q, opal_length = mu_q))
}

#' Match a quarterly-age model to MFCL's initial numbers and recruitment
#' (the quarterly analogue of mfcl_match_recruitment()).
#'
#' B0 is set so that opal's unfished recruitment per step equals MFCL's
#' annual unfished recruitment (R0 = alpha - beta / SPR0 from MFCL's
#' Beverton-Holt curve and the annual biology in `cv`): with one recruitment a
#' year, unfished spawning biomass is then MFCL's B0, and opal's B0 is about
#' four times it.
mfcl_match_recruitment_quarterly = function(qm, mf, report_fn, cv = NULL, max_iter = 50L, tol = 1e-12) {
    d = qm$data
    p = qm$parameters
    nq = d$n_age
    S = exp(-d$M)
    rel = cumprod(c(1, S[-nq]))
    rel[nq] = rel[nq] / (1 - S[nq])
    if (!is.null(cv)) {
        da = cv$data
        Sa = exp(-da$M)
        rel_a = cumprod(c(1, Sa[-da$n_age]))
        rel_a[da$n_age] = rel_a[da$n_age] / (1 - Sa[da$n_age])
        R0_annual = mf$bh[["alpha"]] - mf$bh[["beta"]] / sum(da$spawning_potential * rel_a)
        p$log_B0 = log(R0_annual * sum(d$spawning_potential * rel))
        p$log_h = log(mf$bh[["steepness"]])
    }
    R0 = exp(p$log_B0) / sum(d$spawning_potential * rel)
    # first-year numbers: class a in quarterly age 4 (a - 1) + 1; the plus
    # group (class n_age) in the first quarter of its class; the rest ~0
    target = rep(R0 * rel * 1e-12, length.out = nq)
    first_q = which(qm$age_season == 1L)
    target[first_q] = mf$N_ya[1, qm$age_class[first_q]]
    p$init_rdev_a = log(target / (R0 * rel))
    S0 = exp(p$log_B0); h = exp(p$log_h)
    alpha = 4 * h * R0 / (5 * h - 1)
    beta = S0 * (1 - h) / (5 * h - 1)
    R_target = mf$recruitment_y[qm$rec_year_of_step]
    for (it in seq_len(max_iter)) {
        rp = report_fn(p)
        sb = rp$spawning_biomass_y[qm$rdev_steps + 1L]
        new = log(R_target / (alpha * sb / (beta + sb)))
        delta = max(abs(new - p$rdev_y[qm$rdev_steps]))
        p$rdev_y[qm$rdev_steps] = new
        if (delta < tol) break
    }
    attr(p, "iterations") = it
    p
}
