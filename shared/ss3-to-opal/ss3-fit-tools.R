# Nicholas Ducharme-Barth
# 2026/10/13
# Fit a configured opal_obj from the SS3 converter and extract what the
# vignettes plot: trajectories, index and composition fits, selectivity,
# biology, optimiser diagnostics and timings.
#
# Used by case-studies/io-yft/scripts/03-fit-io-yft-age-k.r. (01-fit-io-yft-ss3-opal.r
# predates these helpers and keeps its own copies.)
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
#' Evaluate an expression and time it
#' @keywords internal
ss3_elapsed = function(expr) {
    t0 = proc.time()[["elapsed"]]
    value = force(expr)
    list(value = value, seconds = proc.time()[["elapsed"]] - t0)
}

#_____________________________________________________________________________
#' Newton refinement from an nlminb optimum
#'
#' nlminb stops on relative function change, which for objectives of order
#' 1e4 leaves gradients of order 1e-3. A few Newton steps with the RTMB
#' Hessian take them below \code{tol}.
#'
#' @param object an RTMB object from \code{opal_rtmb()}.
#' @param par starting parameter vector.
#' @return a \code{list} with \code{par}, \code{max_gradient}, \code{objective}
#'   and the number of \code{steps} taken.
ss3_newton_polish = function(object, par, max_steps = 5L, tol = 1e-8) {
    steps = 0L
    for (i in seq_len(max_steps)) {
        g = as.numeric(object$gr(par))
        if (max(abs(g)) < tol) break
        H = if (is.function(object$he)) object$he(par) else
            stats::optimHess(par, object$fn, object$gr)
        par = par - as.numeric(solve(H, g))
        steps = i
    }
    list(par = par, max_gradient = max(abs(as.numeric(object$gr(par)))),
         objective = object$fn(par), steps = steps)
}

#_____________________________________________________________________________
#' Observed and predicted length compositions in long format
#' @keywords internal
ss3_lf_long = function(data, lf_pred) {
    data.table::rbindlist(lapply(seq_along(data$lf_fishery_f), function(j) {
        obs  = data$lf_obs_data$obs[[j]]
        obs  = obs / rowSums(obs)
        data.table::data.table(
            fishery = data$lf_fishery_f[j],
            ts      = rep(data$lf_year_fi[[j]], times = ncol(obs)),
            n       = rep(data$lf_n_fi[[j]],    times = ncol(obs)),
            bin     = rep(data$len_lower, each = nrow(obs)),
            obs     = as.vector(obs),
            pred    = as.vector(lf_pred[[j]]))
    }))
}

#_____________________________________________________________________________
#' Fit an opal_obj from its starting values and extract vignette quantities
#'
#' Runs \code{opal_fit()} (two nlminb passes), polishes with Newton steps,
#' attaches the result with \code{opal_attach_fit()} and runs
#' \code{opal_check()} at a gradient tolerance of 1e-4.
#'
#' @param x a built \code{opal_obj} at its starting values.
#' @param id \code{character} label for this fit.
#' @param converted output of \code{ss3_to_opal()} (with LF prepared).
#' @param setup_seconds named \code{numeric} of setup timings to report.
#' @return a \code{list} of \code{data.table}s, each with an \code{id} column.
ss3_fit_extract = function(x, id, converted, setup_seconds = numeric()) {

    cat(sprintf("\n==== %s ====\n", id))
    d = converted$data
    n = d$n_year

    start_object = opal::opal_rtmb(x, fresh = TRUE)
    start_par    = start_object$par
    rep_start    = start_object$report(start_par)

    fit    = ss3_elapsed(opal::opal_fit(x, n_passes = 2L, check = FALSE))
    passes = fit$value$fit$diagnostics$optimisation$passes
    tape   = ss3_elapsed(opal::opal_rtmb(fit$value, fresh = TRUE))
    object = tape$value
    polish = ss3_elapsed(ss3_newton_polish(object, fit$value$fit$opt$par))

    opt = fit$value$fit$opt
    opt$par = stats::setNames(polish$value$par, names(opt$par))
    opt$objective = polish$value$objective
    check  = ss3_elapsed(opal::opal_attach_fit(fit$value, opt, check = TRUE,
                                               check_args = list(gradient_tolerance = 1e-4)))
    fitted = check$value
    rep_fit = object$report(opt$par)
    pl = object$env$parList(opt$par)
    val = fitted$validation$fit

    cat(sprintf("  objective %.3f -> %.3f, max|grad| %.2e, check %s; %s\n",
                start_object$fn(start_par), opt$objective, polish$value$max_gradient,
                val$passes, fitted$fit$estimability$message %||% ""))

    traj = function(r, state) data.table::data.table(
        id = id, state = state, ts = seq_len(n),
        ssb         = r$spawning_biomass_y[seq_len(n)],
        ssb0        = r$spawning_biomass0_y[seq_len(n)],
        static_dep  = r$static_depletion_y[seq_len(n)],
        dynamic_dep = r$dynamic_depletion_y[seq_len(n)],
        recruits    = r$number_ysa[seq_len(n), 1, 1])

    cpue = data.table::as.data.table(d$cpue_data)[, .(id = id, ts, fishery, index, obs = value, se)]
    cpue[, `:=`(pred_fit = rep_fit$cpue_pred, pred_start = rep_start$cpue_pred,
                sigma = rep_fit$cpue_sigma)]

    sel_age = function(r, state) data.table::data.table(
        id = id, state = state,
        fishery = rep(seq_len(d$n_fishery), times = d$n_age),
        age     = rep(d$age_a, each = d$n_fishery),
        sel     = as.vector(r$sel_fya[, 1, ]))

    growth = opal::get_growth(d$n_age, d$A1, d$A2, exp(pl$log_L1), exp(pl$log_L2), pl$log_k,
                              min_age = d$min_age, log_k_mult = pl$log_k_mult,
                              k_mult_ages = d$k_mult_ages,
                              plus_group_adjust = isTRUE(d$growth_plus_group_adjust))

    list(
        timing = data.table::data.table(
            id = id,
            stage = c(names(setup_seconds), "opal_fit: tape + nlminb (2 passes)",
                      "rebuild tape", "Newton polish", "attach fit + opal_check"),
            seconds = c(unname(setup_seconds), fit$seconds, tape$seconds,
                        polish$seconds, check$seconds)),
        optimiser = data.table::data.table(
            id = id, pass = seq_along(passes),
            objective   = vapply(passes, `[[`, numeric(1), "objective"),
            convergence = vapply(passes, `[[`, numeric(1), "convergence"),
            iterations  = vapply(passes, function(p) p$iterations, numeric(1)),
            fn_evals    = vapply(passes, function(p) p$evaluations[["function"]], numeric(1)),
            gr_evals    = vapply(passes, function(p) p$evaluations[["gradient"]], numeric(1))),
        summary = data.table::data.table(
            id = id, n_par = length(opt$par),
            n_par_sel = sum(grepl("^par_sel", names(opt$par))),
            objective_start = start_object$fn(start_par),
            objective_fit   = opt$objective,
            newton_steps    = polish$value$steps,
            max_gradient_nlminb = max(abs(as.numeric(object$gr(fit$value$fit$opt$par)))),
            max_gradient    = polish$value$max_gradient,
            check_passes    = isTRUE(val$passes),
            positive_hessian = isTRUE(val$metrics$positive_hessian),
            estimability    = fitted$fit$estimability$message %||% NA_character_,
            B0_start = exp(start_par[["log_B0"]]), B0_fit = exp(pl$log_B0),
            R0_start = rep_start$R0, R0_fit = rep_fit$R0,
            lp_cpue = sum(rep_fit$lp_cpue), lp_lf = sum(rep_fit$lp_lf),
            lp_rec = rep_fit$lp_rec, lp_prior = rep_fit$lp_prior,
            lp_penalty = rep_fit$lp_penalty),
        parameters = data.table::data.table(
            id = id, ts = seq_len(n),
            rdev_start = start_object$env$parList(start_par)$rdev_y,
            rdev_fit = pl$rdev_y),
        q = data.table::data.table(id = id, index = seq_along(pl$log_cpue_q), log_q = pl$log_cpue_q),
        par_sel = data.table::data.table(
            id = id,
            fishery = rep(seq_len(nrow(pl$par_sel)), times = ncol(pl$par_sel)),
            col     = rep(seq_len(ncol(pl$par_sel)), each = nrow(pl$par_sel)),
            start   = as.vector(start_object$env$parList(start_par)$par_sel),
            fit     = as.vector(pl$par_sel),
            estimated = !is.na(as.vector(fitted$map$par_sel %||%
                                             rep(NA, length(pl$par_sel))))),
        trajectory = rbind(traj(rep_start, "start"), traj(rep_fit, "fitted")),
        cpue = cpue,
        lf = if (is.null(rep_fit$lf_pred)) NULL else
            ss3_lf_long(d, rep_fit$lf_pred)[, id := id],
        sel_age = rbind(sel_age(rep_start, "start"), sel_age(rep_fit, "fitted")),
        biology = data.table::data.table(
            id = id, age = d$age_a, M = rep_fit$M_a,
            weight = rep_fit$weight_fya_mod[1, 1, ],
            spawning_potential = rep_fit$spawning_potential_a,
            length = growth),
        reconcile = ss3_reconcile(converted)[, id := id]
    )
}
