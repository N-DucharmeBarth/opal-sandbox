# Nicholas Ducharme-Barth
# 2026/10/14
# A small quarterly age-structured model with switchable conventions for
#   (a) how catch is removed within a time step, and
#   (b) when biology (length-at-age) and the index are evaluated,
# used both as the operating model (OM) and, through RTMB, as the
# estimation model (EM) in the timing simulation.
#
# catch_method:
#   "start"     opal: the step's catch is removed instantaneously at the
#               start, then natural mortality acts for the whole step.
#   "pope_mid"  Pope's approximation: M for half a step, the catch removed
#               instantaneously, then M for the second half.
#   "baranov"   SS3 (F_Method 2/3): F and M act together through the step;
#               F is solved from the observed catch by fixed-point iteration.
#
# bio_timing:
#   "start"     the age-length key (and so selectivity-at-age and the
#               predicted length compositions) uses length at the start of
#               the step, and the index samples abundance at the start.
#   "mid"       both use the middle of the step: length half-way to the
#               next age, abundance after half the step's mortality.
#
# Everything else is shared: Beverton-Holt recruitment at age 1 from SSB at
# the start of the next step (opal's timing), SSB at the start of the step,
# an unfished initial state, and opal's length-based selectivity curves.
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
#' Probability of length at age (rows ages, columns length bins)
#' The first and last bins absorb the tails.
pla_matrix = function(mu_a, sd_a, len_lower, len_upper) {
    "[<-" = RTMB::ADoverload("[<-")
    n_len = length(len_lower)
    out = matrix(0, length(mu_a), n_len)
    for (a in seq_along(mu_a)) {
        cdf = RTMB::pnorm((len_upper[-n_len] - mu_a[a]) / sd_a[a])
        out[a, ] = diff(c(0, cdf, 1))
    }
    out
}

#' Selectivity at length for one fishery, from opal's curve families
#' (1 = logistic, 2 = double normal).
sel_at_length = function(par, type, len_mid) {
    if (type == 1) opal::sel_logistic(len_mid, par) else opal::sel_double_normal(len_mid, par)
}

#' Selectivity at length and at age for every fishery
#' @param par_sel [n_fishery, 6] opal par_sel matrix.
#' @param pla [n_age, n_len] age-length key.
selectivity = function(par_sel, sel_type, pla, len_mid) {
    "[<-" = RTMB::ADoverload("[<-")
    n_f = nrow(par_sel); n_a = nrow(pla)
    sel_fl = matrix(0, n_f, length(len_mid))
    sel_fa = matrix(0, n_f, n_a)
    for (f in seq_len(n_f)) {
        sel_fl[f, ] = sel_at_length(par_sel[f, ], sel_type[f], len_mid)
        for (a in seq_len(n_a)) sel_fa[f, a] = sum(pla[a, ] * sel_fl[f, ])
    }
    list(sel_fl = sel_fl, sel_fa = sel_fa)
}

#' Length-at-age at the start and middle of a step
#' Mid-step length is half-way to the next age's start-of-step length; the
#' plus group continues at the last increment.
length_at_age = function(len_beg, cv = 0.1) {
    n = length(len_beg)
    nxt = c(len_beg[-1], len_beg[n] + 0.5 * (len_beg[n] - len_beg[n - 1]))
    len_mid = (len_beg + nxt) / 2
    list(beg = list(mu = len_beg, sd = cv * len_beg),
         mid = list(mu = len_mid, sd = cv * len_mid))
}

#' One step of catch removal
#' @return list(catch_fa, surv_a = survival over the step, penalty)
remove_catch = function(N_a, C_f, sel_fa, M_a, method, n_iter = 8L) {
    "[<-" = RTMB::ADoverload("[<-")
    n_f = length(C_f)
    eps = 1e-8
    penalty = 0
    cap = function(H) {           # keep total harvest fraction below 0.95
        p = opal::posfun(0.95 - H, eps = 0.01)
        list(H = 0.95 - p$new, penalty = sum(p$penalty))
    }
    n_a = length(N_a)
    catch_fa = matrix(0, n_f, n_a)
    # sum_f x_f * sel_fa[f, ] and sum_a sel_fa[f, a] * y_a, without matrix()
    # reshapes that drop RTMB's AD class
    by_age = function(x_f) { out = 0 * N_a; for (f in seq_len(n_f)) out = out + x_f[f] * sel_fa[f, ]; out }
    by_fleet = function(y_a) { out = 0 * C_f + 0 * y_a[1]; for (f in seq_len(n_f)) out[f] = sum(sel_fa[f, ] * y_a); out }
    if (method == "start" || method == "pope_mid") {
        Nx = if (method == "start") N_a else N_a * exp(-M_a / 2)
        u_f = C_f / (by_fleet(Nx) + eps)
        H = cap(by_age(u_f))
        penalty = H$penalty
        for (f in seq_len(n_f)) catch_fa[f, ] = u_f[f] * sel_fa[f, ] * Nx
        surv_a = if (method == "start") (1 - H$H) * exp(-M_a) else
            exp(-M_a / 2) * (1 - H$H) * exp(-M_a / 2)
    } else {
        # Start from Pope's approximation, then iterate
        # F_f = C_f / sum_a sel_fa N_a (1 - exp(-Z_a)) / Z_a.
        Nm = N_a * exp(-M_a / 2)
        u_f = C_f / (by_fleet(Nm) + eps)
        F_f = -log(1 - pmin_ad(u_f, 0.9))
        for (k in seq_len(n_iter)) {
            Z_a = M_a + by_age(F_f)
            F_f = C_f / (by_fleet(N_a * (1 - exp(-Z_a)) / Z_a) + eps)
        }
        Z_a = M_a + by_age(F_f)
        for (f in seq_len(n_f)) catch_fa[f, ] = F_f[f] * sel_fa[f, ] * N_a * (1 - exp(-Z_a)) / Z_a
        surv_a = exp(-Z_a)
        # A Baranov catch the stock cannot supply shows up as a huge F.
        penalty = sum(opal::posfun(3 - F_f, eps = 0.01)$penalty)
    }
    list(catch_fa = catch_fa, surv_a = surv_a, penalty = penalty)
}

# Smooth min(x, b) that RTMB can differentiate.
pmin_ad = function(x, b) (x + b - sqrt((x - b)^2 + 1e-6)) / 2

#' Population dynamics under a given catch convention
#' @param sel_fa selectivity at age for this run's bio_timing.
#' @return list of numbers at the start (N_ta) and middle (Nmid_ta) of each
#'   step, catch at age, SSB, penalty.
dynamics = function(B0, h, sigma_r, rdev_y, M_a, sp_a, sel_fa, catch_tf, method,
                    n_iter = 8L) {
    "[<-" = RTMB::ADoverload("[<-")
    n_t = nrow(catch_tf); n_a = length(M_a); n_f = ncol(catch_tf)
    rel = numeric(n_a); rel[1] = 1
    for (a in 2:n_a) rel[a] = rel[a - 1] * exp(-M_a[a - 1])
    rel[n_a] = rel[n_a] / (1 - exp(-M_a[n_a]))
    R0 = B0 / sum(sp_a * rel)
    alpha = 4 * h * R0 / (5 * h - 1)
    beta = B0 * (1 - h) / (5 * h - 1)

    N_ta = matrix(0, n_t + 1, n_a)
    Nmid_ta = matrix(0, n_t, n_a)
    catch_tfa = array(0, c(n_t, n_f, n_a))
    ssb_t = numeric(n_t + 1)
    N_ta[1, ] = R0 * rel
    ssb_t[1] = B0
    penalty = 0 * B0
    for (t in seq_len(n_t)) {
        N = N_ta[t, ]
        step = remove_catch(N, catch_tf[t, ], sel_fa, M_a, method, n_iter)
        penalty = penalty + step$penalty
        catch_tfa[t, , ] = step$catch_fa
        Nmid_ta[t, ] = N * sqrt(step$surv_a)
        surv = N * step$surv_a
        N_ta[t + 1, 2:n_a] = surv[1:(n_a - 1)]
        N_ta[t + 1, n_a] = N_ta[t + 1, n_a] + surv[n_a]
        ssb_t[t + 1] = sum(N_ta[t + 1, ] * sp_a)
        N_ta[t + 1, 1] = alpha * ssb_t[t + 1] / (beta + ssb_t[t + 1]) *
            exp(rdev_y[t] - 0.5 * sigma_r^2)
    }
    list(N_ta = N_ta, Nmid_ta = Nmid_ta, catch_tfa = catch_tfa, ssb_t = ssb_t,
         dep_t = ssb_t / B0, R0 = R0, penalty = penalty)
}

#' Expected observations: the index and the length compositions of the catch
expected_obs = function(dyn, sel, pla, idx_fleet, q, timing, comp_fleets) {
    "[<-" = RTMB::ADoverload("[<-")
    Nobs = if (timing == "start") dyn$N_ta[seq_len(nrow(dyn$Nmid_ta)), ] else dyn$Nmid_ta
    n_t = nrow(Nobs); n_a = ncol(Nobs); n_l = ncol(pla)
    index_t = numeric(n_t)
    for (t in seq_len(n_t)) index_t[t] = q * sum(Nobs[t, ] * sel$sel_fa[idx_fleet, ])
    # catch-at-length: sum_a C_fa * pla[a, l] * s_f(l) / sel_fa[a]
    comp = list()
    for (f in comp_fleets) {
        w_al = matrix(0, n_a, n_l)
        for (a in seq_len(n_a)) w_al[a, ] = pla[a, ] * sel$sel_fl[f, ] / (sel$sel_fa[f, a] + 1e-12)
        cl = matrix(0, n_t, n_l)
        for (t in seq_len(n_t)) {
            row = 0 * w_al[1, ]
            for (a in seq_len(n_a)) row = row + dyn$catch_tfa[t, f, a] * w_al[a, ]
            cl[t, ] = row / (sum(row) + 1e-12)
        }
        comp[[f]] = cl
    }
    list(index_t = index_t, comp = comp)
}

#' Build the full model for one convention and parameter set
run_model = function(par, inputs, method, timing, n_iter = 8L) {
    lw = inputs$laa[[if (timing == "start") "beg" else "mid"]]
    pla = pla_matrix(lw$mu, lw$sd, inputs$len_lower, inputs$len_upper)
    sel = selectivity(par$par_sel, inputs$sel_type, pla, inputs$len_mid)
    dyn = dynamics(exp(par$log_B0), inputs$h, inputs$sigma_r, par$rdev_y, inputs$M_a,
                   inputs$sp_a, sel$sel_fa, inputs$catch_tf, method, n_iter)
    obs = expected_obs(dyn, sel, pla, inputs$index_fleet, exp(par$log_q), timing,
                       inputs$comp_fleets)
    c(dyn, obs, list(sel = sel))
}
