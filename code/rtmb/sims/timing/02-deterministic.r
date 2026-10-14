# Nicholas Ducharme-Barth
# 2026/10/14
# Part 1: the structural effect of each convention, with no estimation.
#
# The same parameters (B0, recruitment deviations, selectivity curves) and
# the same catches are projected under every combination of catch method and
# biology timing, at several multiples of the observed catch. Differences
# are therefore purely structural: how much SSB and depletion a convention
# implies for the same removals.
#
# OUTPUT: code/rtmb/sims/timing/data/deterministic.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

    proj_dir = this.path::this.proj()
    dir_sim  = file.path(proj_dir, "code", "rtmb", "sims", "timing")
    source(file.path(dir_sim, "R", "timing-model.R"))
    inputs = readRDS(file.path(dir_sim, "data", "om-inputs.rds"))

    methods = c("start", "pope_mid", "baranov")
    timings = c("start", "mid")
    mults   = c(0.25, 0.5, 0.75, 1, 1.25)

    grid = expand.grid(method = methods, timing = timings, mult = mults,
                       stringsAsFactors = FALSE)
    runs = lapply(seq_len(nrow(grid)), function(i) {
        g = grid[i, ]
        inp = inputs
        inp$catch_tf = inputs$catch_tf * g$mult
        m = run_model(inputs$truth, inp, g$method, g$timing, n_iter = 30L)
        hr = apply(m$catch_tfa, c(1, 3), sum) / m$N_ta[seq_len(inputs$n_t), ]
        pred_catch = apply(m$catch_tfa, c(1, 2), sum)
        data.frame(g, ts = 0:inputs$n_t, ssb = m$ssb_t, dep = m$dep_t,
                   max_hrate = max(hr),
                   catch_err = max(abs(pred_catch / pmax(inp$catch_tf, 1e-9) - 1)[inp$catch_tf > 0]),
                   penalty = as.numeric(m$penalty))
    })
    traj = do.call(rbind, runs)

    # selectivity-at-age under each biology timing
    sel_age = do.call(rbind, lapply(timings, function(tm) {
        lw = inputs$laa[[if (tm == "start") "beg" else "mid"]]
        pla = pla_matrix(lw$mu, lw$sd, inputs$len_lower, inputs$len_upper)
        s = selectivity(inputs$truth$par_sel, inputs$sel_type, pla, inputs$len_mid)$sel_fa
        data.frame(timing = tm, fishery = rep(inputs$fleet_names, times = inputs$n_age),
                   age = rep(seq_len(inputs$n_age), each = nrow(s)), sel = as.vector(s))
    }))

    saveRDS(list(traj = traj, sel_age = sel_age, mults = mults),
            file.path(dir_sim, "data", "deterministic.rds"))

    term = traj[traj$ts == inputs$n_t, c("method", "timing", "mult", "dep", "max_hrate",
                                         "catch_err", "penalty")]
    print(term[order(term$mult, term$timing, term$method), ], digits = 3, row.names = FALSE)
