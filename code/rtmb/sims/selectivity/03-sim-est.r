# Nicholas Ducharme-Barth
# 2026/10/14
# Part 2: simulation-estimation for the selectivity basis.
#
# Two operating models, identical except for selectivity:
#   age      age-based double normal (SS3's curves)
#   length   length-based double normal (C2's curves)
# Catch and biology timing are opal's ("start", "start") in every model.
# Each replicate draws new recruitment deviations (shared by both OMs) and
# observation error, with the timing simulation's data: a longline index
# (CV 0.2, quarters 81-256) and multinomial length compositions (ESS 20) for
# each fishery in every quarter with catch.
#
# Three estimation models per dataset, differing only in selectivity:
#   age          age-based double normal; columns 1-4 (peak, plateau and
#                both widths) estimated for every fishery, as SS3 did
#   length_c2    length-based; C2's choices (1-4, longline 1 and 3)
#   length_all   length-based; 1-4 for every fishery plus the initial and
#                final selectivity (5-6) of the two dome-shaped fisheries
# All put opal's N(0, 2) prior on the estimated selectivity parameters.
# Unestimated columns stay at the start values: the truth when the basis
# matches, otherwise the least-squares approximation (C2's seeding step).
#
# Usage: Rscript 03-sim-est.r [n_rep] [n_cores]
# OUTPUT: code/rtmb/sims/selectivity/data/sim-est.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

    args    = commandArgs(trailingOnly = TRUE)
    n_rep   = if (length(args) >= 1) as.integer(args[1]) else 30L
    n_cores = if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 1L)

    proj_dir = this.path::this.proj()
    dir_sim  = file.path(proj_dir, "code", "rtmb", "sims", "selectivity")
    source(file.path(proj_dir, "code", "rtmb", "sims", "timing", "R", "timing-model.R"))
    inputs = readRDS(file.path(dir_sim, "data", "om-inputs.rds"))
    meth = inputs$method; tim = inputs$timing

    oms = c("age", "length")
    ems = list(age        = list(basis = "age",    free = list(1:4, 1:4, 1:4)),
               length_c2  = list(basis = "length", free = list(1:4, c(1, 3), 1:4)),
               length_all = list(basis = "length", free = list(1:6, 1:4, 1:6)))

    start_sel = function(om, basis)
        if (om == basis) inputs$sel$truth[[basis]] else inputs$sel$approx[[om]][[basis]]

#_____________________________________________________________________________
# Simulate one dataset
    simulate = function(rdev, om) {
        par = inputs$truth
        par$rdev_y = rdev
        par$par_sel = inputs$sel$truth[[om]]
        par$log_q = inputs$sel$log_q[[om]]
        m = run_model(par, inputs, meth, tim, sel_basis = om)
        sd = sqrt(log(1 + inputs$index_cv^2))
        it = inputs$index_t
        index = m$index_t[it] * exp(rnorm(length(it), -sd^2 / 2, sd))
        comps = lapply(inputs$comp_fleets, function(f) {
            keep = inputs$comp_t[inputs$catch_tf[inputs$comp_t, f] > 0]
            obs = t(vapply(keep, function(t) as.numeric(rmultinom(1, inputs$comp_ess, m$comp[[f]][t, ])),
                           numeric(length(inputs$len_mid))))
            list(t = keep, obs = obs)
        })
        hr = apply(m$catch_tfa, c(1, 3), sum) / m$N_ta[seq_len(inputs$n_t), ]
        list(index = index, index_sd = sd, comps = comps, om = om,
             truth = list(B0 = exp(par$log_B0), ssb = m$ssb_t, dep = m$dep_t, rdev = rdev,
                          log_q = par$log_q, sel_fa = m$sel$sel_fa),
             max_hrate = max(hr), penalty = as.numeric(m$penalty))
    }

#_____________________________________________________________________________
# Fit one estimation model to one dataset (starts at the true B0, q and
# recruitment deviations; selectivity as described above)
    fit_em = function(dat, em) {
        sel0 = start_sel(dat$om, em$basis)
        free = em$free
        par0 = list(log_B0 = inputs$truth$log_B0, rdev_y = dat$truth$rdev, log_q = dat$truth$log_q,
                    sel = unlist(lapply(seq_along(free), function(f) sel0[f, free[[f]]])))
        nll = function(p) {
            "[<-" = RTMB::ADoverload("[<-")
            par_sel = sel0
            k = 0
            for (f in seq_along(free)) {
                par_sel[f, free[[f]]] = p$sel[k + seq_along(free[[f]])]
                k = k + length(free[[f]])
            }
            m = run_model(list(log_B0 = p$log_B0, rdev_y = p$rdev_y, log_q = p$log_q,
                               par_sel = par_sel), inputs, meth, tim, sel_basis = em$basis)
            it = inputs$index_t
            out = -sum(RTMB::dnorm(log(dat$index), log(m$index_t[it]) - dat$index_sd^2 / 2,
                                   dat$index_sd, log = TRUE))
            for (f in inputs$comp_fleets) {
                cf = dat$comps[[f]]
                out = out - sum(cf$obs * log(m$comp[[f]][cf$t, ] + 1e-10))
            }
            out = out - sum(RTMB::dnorm(p$rdev_y, 0, inputs$sigma_r, log = TRUE)) -
                sum(RTMB::dnorm(p$sel, 0, 2, log = TRUE)) + m$penalty
            RTMB::REPORT(m$ssb_t)
            RTMB::REPORT(m$dep_t)
            RTMB::REPORT(m$sel$sel_fa)
            out
        }
        obj = RTMB::MakeADFun(nll, par0, silent = TRUE)
        ctrl = list(eval.max = 5000, iter.max = 5000)
        t0 = proc.time()[["elapsed"]]
        opt = tryCatch(nlminb(obj$par, obj$fn, obj$gr, control = ctrl), error = function(e) NULL)
        if (is.null(opt)) return(NULL)
        opt = nlminb(opt$par, obj$fn, obj$gr, control = ctrl)
        for (i in 1:3) {
            g = obj$gr(opt$par)
            if (max(abs(g)) < 1e-6) break
            step = tryCatch(solve(obj$he(opt$par), as.numeric(g)), error = function(e) NULL)
            if (is.null(step)) break
            trial = opt$par - step
            if (!is.finite(obj$fn(trial)) || obj$fn(trial) > opt$objective + 1e-8) break
            opt$par = trial
            opt$objective = obj$fn(trial)
        }
        rep = obj$report(opt$par)
        list(B0 = exp(opt$par[["log_B0"]]), ssb = rep[["m$ssb_t"]], dep = rep[["m$dep_t"]],
             sel_fa = rep[["m$sel$sel_fa"]], objective = opt$objective,
             convergence = opt$convergence, max_gradient = max(abs(obj$gr(opt$par))),
             seconds = proc.time()[["elapsed"]] - t0)
    }

#_____________________________________________________________________________
# Run
    datasets = lapply(seq_len(n_rep), function(r) {
        set.seed(2000 + r)
        for (try in 1:50) {
            rdev = rnorm(inputs$n_t, 0, inputs$sigma_r)
            sims = setNames(lapply(oms, function(om) simulate(rdev, om)), oms)
            if (all(vapply(sims, function(s) s$max_hrate < 0.5 && s$penalty < 1e-6, logical(1)))) break
        }
        list(sims = sims, redraws = try - 1)
    })
    tasks = expand.grid(em = names(ems), om = oms, rep = seq_len(n_rep), stringsAsFactors = FALSE)

    one_fit = function(i) {
        k = tasks[i, ]
        dat = datasets[[k$rep]]$sims[[k$om]]
        fit = tryCatch(fit_em(dat, ems[[k$em]]), error = function(err) NULL)
        tru = dat$truth
        n = inputs$n_t + 1
        ok = !is.null(fit)
        summ = data.frame(
            rep = k$rep, om = k$om, em = k$em, redraws = datasets[[k$rep]]$redraws,
            B0_true = tru$B0, dep_T_true = tru$dep[n], ssb_T_true = tru$ssb[n],
            B0 = if (ok) fit$B0 else NA, dep_T = if (ok) fit$dep[n] else NA,
            ssb_T = if (ok) fit$ssb[n] else NA,
            mare_dep = if (ok) mean(abs(fit$dep - tru$dep)) else NA,
            objective = if (ok) fit$objective else NA,
            convergence = if (ok) fit$convergence else NA,
            max_gradient = if (ok) fit$max_gradient else NA,
            seconds = if (ok) fit$seconds else NA)
        traj = if (ok) data.frame(rep = k$rep, om = k$om, em = k$em, ts = 0:inputs$n_t,
                                  ssb = fit$ssb, dep = fit$dep) else NULL
        sel = if (ok) data.frame(rep = k$rep, om = k$om, em = k$em,
                                 fishery = rep(inputs$fleet_names, times = inputs$n_age),
                                 age = rep(seq_len(inputs$n_age), each = 3),
                                 sel = as.vector(fit$sel_fa)) else NULL
        list(summary = summ, traj = traj, sel = sel)
    }

    t0 = proc.time()[["elapsed"]]
    res = parallel::mclapply(seq_len(nrow(tasks)), one_fit, mc.cores = n_cores,
                             mc.preschedule = FALSE)
    bad = vapply(res, inherits, logical(1), "try-error")
    if (any(bad)) message(sum(bad), " fit(s) failed: ", conditionMessage(attr(res[[which(bad)[1]]], "condition")))
    pick = function(what) do.call(rbind, lapply(res[!bad], `[[`, what))
    truth = do.call(rbind, lapply(seq_len(n_rep), function(r) do.call(rbind, lapply(oms, function(o) {
        tr = datasets[[r]]$sims[[o]]$truth
        data.frame(rep = r, om = o, ts = 0:inputs$n_t, ssb = tr$ssb, dep = tr$dep)
    }))))
    out = list(results = pick("summary"), traj = pick("traj"), sel = pick("sel"), truth = truth,
               n_rep = n_rep, created = format(Sys.time(), tz = "UTC", usetz = TRUE),
               minutes = (proc.time()[["elapsed"]] - t0) / 60)
    saveRDS(out, file.path(dir_sim, "data", "sim-est.rds"))

    r = out$results
    cat(sprintf("%d replicates x %d fits in %.1f min\n", n_rep, length(oms) * length(ems), out$minutes))
    r$rel_B0 = r$B0 / r$B0_true - 1
    r$err_dep = r$dep_T - r$dep_T_true
    print(aggregate(cbind(rel_B0, err_dep, mare_dep, max_gradient, seconds) ~ om + em,
                    data = r, FUN = median), digits = 3)
