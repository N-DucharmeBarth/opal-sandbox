# Nicholas Ducharme-Barth
# 2026/10/14
# Operating-model inputs for the selectivity simulation.
#
# Reuses the timing simulation's inputs (biology, three pooled fisheries,
# catches, B0) and model (../timing/R/timing-model.R), with the catch and
# biology timing fixed at opal's conventions ("start", "start") so only
# selectivity differs. Two true selectivities:
#
#   length  opal's length-based double normal with C2's parameters (as in
#           the timing simulation), reaching age through the age-length key.
#   age     the same double-normal family evaluated on age (as SS3's
#           age-based pattern 20), fitted by least squares to SS3's
#           selectivity-at-age for the same source fleets (11, 4, 8).
#
# Each is also approximated in the other basis by least squares on
# selectivity-at-age (C2's seeding step): these give Part 1's comparison and
# the starting values of the mismatched estimation models.
#
# OUTPUT: simulations/selectivity/data/om-inputs.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

    proj_dir = this.path::this.proj()
    dir_sim  = file.path(proj_dir, "simulations", "selectivity")
    dir_tim  = file.path(proj_dir, "simulations", "timing")
    dir.create(file.path(dir_sim, "data"), showWarnings = FALSE)
    source(file.path(dir_tim, "R", "timing-model.R"))

    inputs = readRDS(file.path(dir_tim, "data", "om-inputs.rds"))
    load(file.path(proj_dir, "case-studies", "io-yft", "opal-export", "data", "io_yft_reference.rda"))
    ss3_sel = io_yft_reference$ss3$sel_age
    sel_source = c(11, 4, 8)

    lw  = inputs$laa$beg
    pla = pla_matrix(lw$mu, lw$sd, inputs$len_lower, inputs$len_upper)
    sel_age_of = function(par_sel, basis)
        selectivity(par_sel, inputs$sel_type, pla, inputs$len_mid, basis)$sel_fa

    # Least-squares fit of one fishery's curve, in a given basis, to a target
    # selectivity-at-age. Parameters stay inside opal's default bounds.
    fit_curve = function(target, basis, start = c(0, 0, -3, -1, -5, 0)) {
        obj = function(p) {
            s = sel_age_of(matrix(p, 1), basis)[1, ]
            sum((s - target)^2)
        }
        best = NULL
        for (st in list(start, c(-0.5, -2, -4, -1, -7, 0), c(-1, 1, -3, -3, -7, 7))) {
            o = optim(st, obj, method = "L-BFGS-B", lower = -7, upper = 7,
                      control = list(maxit = 5000))
            if (is.null(best) || o$value < best$value) best = o
        }
        best
    }

    # True age-based curves, fitted to SS3's selectivity-at-age
    age_fit = lapply(sel_source, function(f) fit_curve(ss3_sel$sel[ss3_sel$fishery == f], "age"))
    par_age = do.call(rbind, lapply(age_fit, `[[`, "par"))

    par_len = inputs$truth$par_sel
    # Each truth approximated in the other basis
    len_of_age = lapply(seq_len(3), function(f) fit_curve(sel_age_of(par_age, "age")[f, ], "length"))
    age_of_len = lapply(seq_len(3), function(f) fit_curve(sel_age_of(par_len, "length")[f, ], "age"))

    inputs$method = "start"
    inputs$timing = "start"
    inputs$sel = list(
        truth = list(age = par_age, length = par_len),
        approx = list(age = list(length = do.call(rbind, lapply(len_of_age, `[[`, "par"))),
                      length = list(age = do.call(rbind, lapply(age_of_len, `[[`, "par")))),
        rmse = data.frame(
            fishery = inputs$fleet_names,
            age_truth_vs_ss3 = sqrt(vapply(age_fit, `[[`, numeric(1), "value") / inputs$n_age),
            length_approx_of_age = sqrt(vapply(len_of_age, `[[`, numeric(1), "value") / inputs$n_age),
            age_approx_of_length = sqrt(vapply(age_of_len, `[[`, numeric(1), "value") / inputs$n_age)),
        ss3_sel_age = ss3_sel[ss3_sel$fishery %in% sel_source, ])

    # q so that the index averages one, for each truth
    for (b in c("age", "length")) {
        par = inputs$truth; par$par_sel = inputs$sel$truth[[b]]; par$log_q = 0
        m = run_model(par, inputs, "start", "start", sel_basis = b)
        inputs$sel$log_q[[b]] = -log(mean(m$index_t[inputs$index_t]))
        cat(sprintf("%-6s truth: terminal depletion %.3f, max harvest fraction %.3f\n", b,
                    m$dep_t[inputs$n_t + 1],
                    max(apply(m$catch_tfa, c(1, 3), sum) / m$N_ta[1:inputs$n_t, ])))
    }
    print(inputs$sel$rmse, digits = 3)
    saveRDS(inputs, file.path(dir_sim, "data", "om-inputs.rds"))
