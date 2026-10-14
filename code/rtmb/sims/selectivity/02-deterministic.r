# Nicholas Ducharme-Barth
# 2026/10/14
# Part 1: the structural effect of the selectivity basis, with no estimation.
#
# Each true selectivity (age-based, length-based) is projected alongside its
# least-squares approximation in the other basis, with the same B0,
# recruitment deviations and catch (scaled 0.25-1.25), and opal's timing
# conventions. This is the error C2's seeding step makes before any fitting.
#
# OUTPUT: code/rtmb/sims/selectivity/data/deterministic.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

    proj_dir = this.path::this.proj()
    dir_sim  = file.path(proj_dir, "code", "rtmb", "sims", "selectivity")
    source(file.path(proj_dir, "code", "rtmb", "sims", "timing", "R", "timing-model.R"))
    inputs = readRDS(file.path(dir_sim, "data", "om-inputs.rds"))

    cases = list(
        list(truth = "age",    basis = "age",    par = inputs$sel$truth$age,           label = "age truth"),
        list(truth = "age",    basis = "length", par = inputs$sel$approx$age$length,   label = "length approximation"),
        list(truth = "length", basis = "length", par = inputs$sel$truth$length,        label = "length truth"),
        list(truth = "length", basis = "age",    par = inputs$sel$approx$length$age,   label = "age approximation"))
    mults = c(0.25, 0.5, 0.75, 1, 1.25)

    traj = do.call(rbind, lapply(cases, function(k) do.call(rbind, lapply(mults, function(mu) {
        inp = inputs
        inp$catch_tf = inputs$catch_tf * mu
        par = inputs$truth
        par$par_sel = k$par
        m = run_model(par, inp, inputs$method, inputs$timing, sel_basis = k$basis)
        hr = apply(m$catch_tfa, c(1, 3), sum) / m$N_ta[seq_len(inputs$n_t), ]
        data.frame(truth = k$truth, curve = k$label, basis = k$basis, mult = mu,
                   ts = 0:inputs$n_t, ssb = m$ssb_t, dep = m$dep_t, max_hrate = max(hr),
                   penalty = as.numeric(m$penalty))
    }))))

    lw  = inputs$laa$beg
    pla = pla_matrix(lw$mu, lw$sd, inputs$len_lower, inputs$len_upper)
    sel_age = do.call(rbind, lapply(cases, function(k) {
        s = selectivity(k$par, inputs$sel_type, pla, inputs$len_mid, k$basis)$sel_fa
        data.frame(truth = k$truth, curve = k$label, fishery = rep(inputs$fleet_names, times = inputs$n_age),
                   age = rep(seq_len(inputs$n_age), each = nrow(s)), sel = as.vector(s))
    }))

    saveRDS(list(traj = traj, sel_age = sel_age, mults = mults),
            file.path(dir_sim, "data", "deterministic.rds"))

    term = traj[traj$ts == inputs$n_t, c("truth", "curve", "mult", "dep", "max_hrate", "penalty")]
    print(term[order(term$truth, term$mult), ], digits = 3, row.names = FALSE)
