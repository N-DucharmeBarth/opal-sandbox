# Nicholas Ducharme-Barth
# 2026/10/14
# Part 2: simulation-estimation.
#
# Two operating models generate data:
#   ss3_like   catch "baranov", biology "mid"   (SS3's conventions)
#   opal_like  catch "start",   biology "start" (opal's conventions)
# Each replicate draws new recruitment deviations (shared by both OMs, so the
# OMs differ only in their conventions) and new observation error: a
# longline index (lognormal, CV 0.2, quarters 81-256) and multinomial length
# compositions of each fishery's catch (ESS 20, quarters 81-256). Catch is
# known exactly.
#
# Every OM dataset is fitted by six estimation models: catch "start",
# "pope_mid" or "baranov" crossed with biology "start" or "mid". They
# estimate log B0, 256 recruitment deviations (N(0, sigma_R) penalty), log q
# and ten selectivity parameters; growth, M, steepness and sigma_R are
# fixed at their true values. The EM matching its OM is the control.
#
# Usage: Rscript 03-sim-est.r [n_rep] [n_cores]
# OUTPUT: code/rtmb/sims/timing/data/sim-est.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

    args    = commandArgs(trailingOnly = TRUE)
    n_rep   = if (length(args) >= 1) as.integer(args[1]) else 30L
    n_cores = if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 1L)

    proj_dir = this.path::this.proj()
    dir_sim  = file.path(proj_dir, "code", "rtmb", "sims", "timing")
    source(file.path(dir_sim, "R", "timing-model.R"))
    inputs = readRDS(file.path(dir_sim, "data", "om-inputs.rds"))

    oms = list(ss3_like = c(method = "baranov", timing = "mid"),
               opal_like = c(method = "start", timing = "start"))
    ems = expand.grid(method = c("start", "pope_mid", "baranov"),
                      timing = c("start", "mid"), stringsAsFactors = FALSE)

#_____________________________________________________________________________
# Simulate one dataset
    simulate = function(rdev, om) {
        par = inputs$truth
        par$rdev_y = rdev
        m = run_model(par, inputs, om[["method"]], om[["timing"]], n_iter = 30L)
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
        list(index = index, index_sd = sd, comps = comps,
             truth = list(B0 = exp(par$log_B0), ssb = m$ssb_t, dep = m$dep_t, rdev = rdev),
             max_hrate = max(hr), penalty = as.numeric(m$penalty))
    }

#_____________________________________________________________________________
# Fit one estimation model to one dataset
    fit_em = function(dat, em, n_iter = 8L) {
        # Every EM starts at the true values (with the replicate's true
        # recruitment deviations); mismatched EMs then move to their own
        # optimum. From zero deviations nlminb can stall far from it.
        tr = inputs$truth
        free = inputs$sel_free
        par0 = list(log_B0 = tr$log_B0, rdev_y = dat$truth$rdev, log_q = tr$log_q,
                    sel = unlist(lapply(seq_along(free), function(f) tr$par_sel[f, free[[f]]])))
        nll = function(p) {
            "[<-" = RTMB::ADoverload("[<-")
            par_sel = tr$par_sel
            k = 0
            for (f in seq_along(free)) {
                par_sel[f, free[[f]]] = p$sel[k + seq_along(free[[f]])]
                k = k + length(free[[f]])
            }
            m = run_model(list(log_B0 = p$log_B0, rdev_y = p$rdev_y, log_q = p$log_q,
                               par_sel = par_sel), inputs, em$method, em$timing, n_iter)
            it = inputs$index_t
            out = -sum(RTMB::dnorm(log(dat$index), log(m$index_t[it]) - dat$index_sd^2 / 2,
                             dat$index_sd, log = TRUE))
            for (f in inputs$comp_fleets) {
                cf = dat$comps[[f]]
                out = out - sum(cf$obs * log(m$comp[[f]][cf$t, ] + 1e-10))
            }
            out = out - sum(RTMB::dnorm(p$rdev_y, 0, inputs$sigma_r, log = TRUE)) + m$penalty
            RTMB::REPORT(m$ssb_t)
            RTMB::REPORT(m$dep_t)
            out
        }
        obj = RTMB::MakeADFun(nll, par0, silent = TRUE)
        ctrl = list(eval.max = 5000, iter.max = 5000)
        t0 = proc.time()[["elapsed"]]
        opt = tryCatch(nlminb(obj$par, obj$fn, obj$gr, control = ctrl), error = function(e) NULL)
        if (is.null(opt)) return(NULL)
        opt = nlminb(opt$par, obj$fn, obj$gr, control = ctrl)
        # Newton polish from the nlminb optimum
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
             objective = opt$objective, convergence = opt$convergence,
             max_gradient = max(abs(obj$gr(opt$par))),
             seconds = proc.time()[["elapsed"]] - t0)
    }

#_____________________________________________________________________________
# Run
    # Datasets first (cheap, serial), then every fit as its own parallel task
    datasets = lapply(seq_len(n_rep), function(r) {
        set.seed(1000 + r)
        # Redraw recruitment until neither OM exceeds a 50% quarterly harvest
        # fraction at any age (the observed catch is fixed).
        for (try in 1:50) {
            rdev = rnorm(inputs$n_t, 0, inputs$sigma_r)
            sims = lapply(oms, function(om) simulate(rdev, om))
            if (all(vapply(sims, function(s) s$max_hrate < 0.5 && s$penalty < 1e-6, logical(1)))) break
        }
        list(sims = sims, redraws = try - 1)
    })
    tasks = expand.grid(e = seq_len(nrow(ems)), om = names(oms), rep = seq_len(n_rep),
                        stringsAsFactors = FALSE)

    one_fit = function(i) {
        k = tasks[i, ]
        dat = datasets[[k$rep]]$sims[[k$om]]
        fit = tryCatch(fit_em(dat, ems[k$e, ]), error = function(err) NULL)
        tru = dat$truth
        n = inputs$n_t + 1
        ok = !is.null(fit)
        summ = data.frame(
            rep = k$rep, om = k$om, em_method = ems$method[k$e], em_timing = ems$timing[k$e],
            redraws = datasets[[k$rep]]$redraws,
            B0_true = tru$B0, dep_T_true = tru$dep[n], ssb_T_true = tru$ssb[n],
            B0 = if (ok) fit$B0 else NA, dep_T = if (ok) fit$dep[n] else NA,
            ssb_T = if (ok) fit$ssb[n] else NA,
            mare_dep = if (ok) mean(abs(fit$dep - tru$dep)) else NA,
            objective = if (ok) fit$objective else NA,
            convergence = if (ok) fit$convergence else NA,
            max_gradient = if (ok) fit$max_gradient else NA,
            seconds = if (ok) fit$seconds else NA)
        traj = if (ok) data.frame(rep = k$rep, om = k$om, em_method = ems$method[k$e],
                                  em_timing = ems$timing[k$e], ts = 0:inputs$n_t,
                                  ssb = fit$ssb, dep = fit$dep) else NULL
        list(summary = summ, traj = traj)
    }

    t0 = proc.time()[["elapsed"]]
    res = parallel::mclapply(seq_len(nrow(tasks)), one_fit, mc.cores = n_cores,
                             mc.preschedule = FALSE)
    bad = vapply(res, inherits, logical(1), "try-error")
    if (any(bad)) message(sum(bad), " fit(s) failed: ", conditionMessage(attr(res[[which(bad)[1]]], "condition")))
    traj = do.call(rbind, lapply(res[!bad], `[[`, "traj"))
    res = do.call(rbind, lapply(res[!bad], `[[`, "summary"))
    truth = do.call(rbind, lapply(seq_len(n_rep), function(r) do.call(rbind, lapply(names(oms), function(o) {
        tr = datasets[[r]]$sims[[o]]$truth
        data.frame(rep = r, om = o, ts = 0:inputs$n_t, ssb = tr$ssb, dep = tr$dep)
    }))))
    saveRDS(list(results = res, traj = traj, truth = truth, n_rep = n_rep, created = format(Sys.time(), tz = "UTC", usetz = TRUE),
                 minutes = (proc.time()[["elapsed"]] - t0) / 60),
            file.path(dir_sim, "data", "sim-est.rds"))
    cat(sprintf("%d replicates x %d fits in %.1f min\n", n_rep, 2 * nrow(ems),
                (proc.time()[["elapsed"]] - t0) / 60))
    res$rel_B0 = res$B0 / res$B0_true - 1
    res$err_dep = res$dep_T - res$dep_T_true
    print(aggregate(cbind(rel_B0, err_dep, mare_dep, seconds) ~ om + em_method + em_timing,
                    data = res, FUN = median), digits = 3)
