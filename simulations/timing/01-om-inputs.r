# Nicholas Ducharme-Barth
# 2026/10/14
# Operating-model inputs for the catch- and biology-timing simulation.
#
# Biology, length-at-age, catches and true selectivity come from the IO
# yellowfin conversion (case-studies/io-yft/opal-export/data and the fitted C2 opal
# object), simplified to three fisheries by gear so the estimation models fit
# in seconds:
#   1 purse seine   SS3 fleets 11-13 (ps_1, ps_2, ps_4); selectivity of fleet 11
#   2 longline      SS3 fleets 4-7   (ll_1-ll_4);       selectivity of fleet 4
#   3 other         all remaining fleets (gillnet, handline, bait boat, troll,
#                   other);                              selectivity of fleet 8
# The longline fishery also provides the index (its selectivity, one q).
#
# OUTPUT: simulations/timing/data/om-inputs.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

    proj_dir = this.path::this.proj()
    dir_sim  = file.path(proj_dir, "simulations", "timing")
    dir_exp  = file.path(proj_dir, "case-studies", "io-yft", "opal-export")
    dir.create(file.path(dir_sim, "data"), showWarnings = FALSE)
    source(file.path(dir_sim, "R", "timing-model.R"))

    for (f in list.files(file.path(dir_exp, "data"), full.names = TRUE)) load(f)
    d   = io_yft_data
    fit = opal::opal_read(file.path(dir_exp, "inst", "extdata", "io_yft_c2_fit.rds"))
    pl  = opal::opal_rtmb(fit, fresh = TRUE)$env$parList(fit$fit$opt$par)

    groups = list(purse_seine = 11:13, longline = 4:7,
                  other = c(1, 2, 3, 8, 9, 10, 14, 15, 16))
    sel_source = c(purse_seine = 11, longline = 4, other = 8)

    catch_tf = sapply(groups, function(g) rowSums(d$catch_obs_ysf[, 1, g, drop = FALSE]))

    inputs = list(
        n_t        = d$n_year,
        n_age      = d$n_age,
        M_a        = d$M,
        sp_a       = d$spawning_potential,
        laa        = length_at_age(io_yft_reference$ss3$length_at_age$length, cv = 0.1),
        len_lower  = d$len_lower,
        len_upper  = d$len_upper,
        len_mid    = d$len_mid,
        h          = exp(pl$log_h),
        sigma_r    = exp(pl$log_sigma_r),
        catch_tf   = unname(catch_tf),
        fleet_names = names(groups),
        groups     = groups,
        sel_type   = c(2L, 2L, 2L),
        index_fleet = 2L,
        comp_fleets = 1:3,
        # observations
        index_t    = 81:d$n_year,
        index_cv   = 0.2,
        comp_t     = seq_len(d$n_year),
        comp_ess   = 20,
        # truth
        truth = list(
            log_B0  = pl$log_B0,
            rdev_y  = pl$rdev_y,             # C2's estimates; replicates redraw
            par_sel = unname(pl$par_sel[sel_source, ]),
            log_q   = 0),
        # which selectivity columns the estimation models estimate
        # (C2's choices for the source fleets)
        sel_free = list(c(1, 2, 3, 4), c(1, 3), c(1, 2, 3, 4))
    )

    # q so that the index averages one under the SS3-like OM
    om = run_model(inputs$truth, inputs, "baranov", "mid", n_iter = 30L)
    inputs$truth$log_q = -log(mean(om$index_t[inputs$index_t]))

    saveRDS(inputs, file.path(dir_sim, "data", "om-inputs.rds"))
    cat(sprintf("B0 %.0f t, terminal depletion (SS3-like OM) %.3f, max harvest fraction %.3f\n",
                exp(inputs$truth$log_B0), om$dep_t[inputs$n_t + 1],
                max(apply(om$catch_tfa, c(1, 3), sum) / om$N_ta[1:inputs$n_t, ])))
