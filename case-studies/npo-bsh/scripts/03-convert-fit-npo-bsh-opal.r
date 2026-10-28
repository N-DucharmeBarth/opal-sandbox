# Nicholas Ducharme-Barth
# 2026/10/28
# Convert the single-sex SS3 bridge model B2e of the NPO blue shark
# assessment to opal, test the conversion, and fit opal.
#
#   O0  bridge: opal at SS3's parameter values with SS3's selectivity-at-age;
#       nothing estimated
#   O1  SS3's selectivity-at-age held fixed; B0, recruitment deviations,
#       initial age structure and the two catchabilities estimated (F4's
#       initial F stays at SS3's value; see the maps below)
#   O2  opal's own length-based double normals, started exactly at SS3's
#       curves; the selectivity parameters SS3 estimated are estimated, plus
#       everything O1 estimates
#
# Growth, steepness (0.613), sigma_R (0.4), M and the biology stay at SS3's
# values, as in SS3. opal's default priors are switched off (SS3's prior
# lambda is 0). All three use one season per year: opal removes each year's
# catch at the start of the year, SS3 continuously (hybrid F).
#
# OUTPUTS:
#   case-studies/npo-bsh/vignettes/data/npo-bsh-opal-cache.rds
#   case-studies/npo-bsh/vignettes/data/npo-bsh-opal-fits.rds  (fitted opal_obj's)
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Packages
    library(data.table)
    library(magrittr)
    library(opal)
    library(r4ss)
    library(yaml)

#_____________________________________________________________________________
# Project paths
    proj_dir  = this.path::this.proj()
    dir_fns   = file.path(proj_dir, "shared", "ss3-to-opal")
    dir_cfg   = file.path(proj_dir, "case-studies", "npo-bsh", "config")
    dir_cache = file.path(proj_dir, "case-studies", "npo-bsh", "vignettes", "data")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
    invisible(sapply(list.files(dir_fns, pattern = "\\.[Rr]$", full.names = TRUE), source))

#_____________________________________________________________________________
# Convert and test the conversion
    cfg = ss3_read_config(file.path(dir_cfg, "npo-bsh-b2e.yaml"))
    converted = ss3_to_opal(cfg, dir_root = proj_dir)

    lf_wide = dcast(converted$lf, fishery + year + month + ts ~ bin,
                    value.var = "value", fill = 0) %>%
              setorder(fishery, ts) %>%
              as.data.frame()
    converted$data = prep_lf_data(
        converted$data, lf_wide,
        lf_switch         = 1L,                       # multinomial, as SS3
        lf_keep_fisheries = as.integer(cfg$lf$keep_fisheries),
        lf_var_adjust     = converted$meta$lf_var_adjust,
        lf_cap            = cfg$lf$cap,
        lf_addtocomp      = as.numeric(cfg$lf$addtocomp)
    ) %>% ss3_make_portable()
    # SS3 has no parameter priors in the likelihood (prior lambda 0)
    converted$data$lambdas = list(prior = 0)

    checks = validate_opal_inputs(converted, strict = FALSE)
    recon  = ss3_reconcile(converted, tol = 1e-6)
    stopifnot("Structural validation failed" = all(checks$pass),
              "Reconciliation with SS3 failed" = all(recon$pass, na.rm = TRUE))

    assessment = ss3_opal_obj(converted, cfg,
                              control = list(eval.max = 10000, iter.max = 10000),
                              build = FALSE)

#_____________________________________________________________________________
# Selectivity
    # opal's double normals at SS3's curves (exact to < 1e-4 at the SS3
    # parameters; see ss3_seed_par_sel_size_dblnorm())
    seed = ss3_seed_par_sel_size_dblnorm(converted$ss3, converted$meta$fisheries,
                                         converted$data$len_mid,
                                         converted$parameters$par_sel)
    sizeselex = as.data.table(converted$ss3$rep$sizeselex)[Factor == "Lsel" & Yr == 2020]
    len_cols = setdiff(names(sizeselex), c("Factor", "Fleet", "Yr", "Sex", "Label"))
    ss3_fleet_f = c(converted$meta$fisheries$extract_fleets,
                    converted$meta$fisheries$index_map$ss3_fleet)
    seed_check = data.table(
        fishery    = seq_along(ss3_fleet_f),
        ss3_fleet  = ss3_fleet_f,
        name       = converted$meta$dims$.fleet_names[ss3_fleet_f],
        source     = seed$source_fleet,
        max_abs_diff = vapply(seq_along(ss3_fleet_f), function(i)
            max(abs(sel_double_normal(converted$data$len_mid, seed$par_sel[i, ]) -
                    as.numeric(sizeselex[Fleet == ss3_fleet_f[i], ..len_cols]))), numeric(1L)))

    # Which double-normal parameters SS3 estimated, by opal fishery (same order:
    # peak, top, ascending width, descending width, start, end)
    pars_ss3 = as.data.table(converted$ss3$rep$parameters)
    est_cols = lapply(seq_along(ss3_fleet_f), function(i) {
        sf = seed$source_fleet[i]
        if (is.na(sf) || i %in% as.integer(names(cfg$selectivity$mirror))) return(integer())
        if (!i %in% as.integer(cfg$selectivity$estimate)) return(integer())
        p = pars_ss3[grepl("^Size_DblN_", Label) & grepl(sprintf("\\(%d\\)$", sf), Label)]
        which(!is.na(p$Active_Cnt))
    })

#_____________________________________________________________________________
# Maps
    # Initial age structure: SS3 estimates the 1971 recruitment (age 0) and
    # the ten early deviations (ages 1-10)
    n_init = sum(converted$parameters$init_rdev_a != 0 | seq_len(converted$data$n_age) <= 11L)
    map_init_rdev = factor(c(seq_len(n_init), rep(NA, converted$data$n_age - n_init)))
    # F4's initial F stays at SS3's value: SS3 ties it to the initial-
    # equilibrium catch (40,000 t, SE 0.05) through a likelihood term opal does
    # not have, and left free in opal it doubles (0.11 -> 0.22)
    base_map = function(x) {
        m = x$map %||% get_map(x$parameters)
        m$init_rdev_a = map_init_rdev
        m
    }

    # Selectivity map: SS3's estimated columns; mirrors share their source's
    sel_map = function(par_sel) {
        m = matrix(NA_integer_, nrow(par_sel), ncol(par_sel))
        lvl = 0L
        for (f in seq_along(est_cols)) {
            cols = est_cols[[f]]
            if (!length(cols)) next
            m[f, cols] = lvl + seq_along(cols)
            lvl = lvl + length(cols)
        }
        for (f in names(cfg$selectivity$mirror)) {
            m[as.integer(f), ] = m[as.integer(cfg$selectivity$mirror[[f]]), ]
        }
        factor(as.vector(m))
    }

    # Bounds for the estimated selectivity parameters: SS3's own bounds on
    # opal's scale (opal's default of +/-7 is narrower than SS3's +/-15 on
    # the widths once they are rescaled by the spread of the length bins)
    sel_bounds = function() {
        x = converted$data$len_mid
        mu = mean(x)
        sd_x = stats::sd(x)
        to_opal = function(v, col) switch(col, (v - mu) / sd_x, v, v - 2 * log(sd_x), v - 2 * log(sd_x))
        out = list(lower = numeric(), upper = numeric())
        for (f in seq_along(est_cols)) {
            sf = seed$source_fleet[f]
            for (col in est_cols[[f]]) {
                p = pars_ss3[grepl("^Size_DblN_", Label) & grepl(sprintf("\\(%d\\)$", sf), Label)][col]
                out$lower = c(out$lower, to_opal(p$Min, col))
                out$upper = c(out$upper, to_opal(p$Max, col))
            }
        }
        out
    }

#_____________________________________________________________________________
# O0: bridge at SS3's values
    o0 = ss3_use_external_selectivity(assessment, converted$selectivity$sel_fa_external) %>%
        opal_build()
    bridge = ss3_bridge_table(converted, opal_report(o0))

#_____________________________________________________________________________
# Fits
    fit_opal = function(x) {
        t0  = proc.time()[["elapsed"]]
        fit = opal_fit(x, n_passes = 3L, check = FALSE,
                       control = list(eval.max = 1e5, iter.max = 1e5))
        obj = opal_rtmb(fit, fresh = TRUE)
        opt = fit$fit$opt
        # Newton steps to bring the gradient down; keep the nlminb optimum if a
        # step lands where the objective is undefined
        pol = tryCatch(ss3_newton_polish(obj, opt$par), error = function(e) NULL)
        if (is.null(pol) || !is.finite(pol$objective) || pol$objective > opt$objective + 1e-8) {
            g = obj$gr(opt$par)
            pol = list(par = opt$par, objective = opt$objective,
                       max_gradient = max(abs(g)), steps = 0L)
        }
        opt$par = setNames(pol$par, names(opt$par))
        opt$objective = pol$objective
        fit = opal_attach_fit(fit, opt, check = TRUE,
                              check_args = list(gradient_tolerance = 1e-4))
        attr(fit, "seconds") = proc.time()[["elapsed"]] - t0
        attr(fit, "max_gradient") = pol$max_gradient
        fit
    }

    # The fits take tens of minutes; they are saved as soon as they finish and
    # reused on later runs unless NPO_BSH_REFIT=true
    fits_file = file.path(dir_cache, "npo-bsh-opal-fits.rds")
    refit = !file.exists(fits_file) || identical(Sys.getenv("NPO_BSH_REFIT"), "true")
    saved = if (refit) list() else readRDS(fits_file)

    # O1: SS3 selectivity fixed
    o1 = o0
    m1 = base_map(o1)
    o1 = opal_update(o1, map = m1) %>% opal_build()
    fit_o1 = if (refit) fit_opal(o1) else saved$O1

    # O2: opal length-based selectivity, started at SS3's curves
    # O2 starts from O1's estimates for everything but selectivity
    o2 = ss3_use_opal_selectivity(assessment, par_sel = seed$par_sel)
    p2 = fit_o1$fit$parameters
    p2$par_sel = seed$par_sel
    o2 = opal_update(o2, parameters = p2)
    m2 = base_map(o2)
    m2$par_sel = sel_map(seed$par_sel)
    o2 = opal_update(o2, map = m2) %>% opal_build()
    b = o2$bounds
    i_sel = which(names(opal_rtmb(o2, fresh = TRUE)$par) == "par_sel")
    sb = sel_bounds()
    stopifnot(length(i_sel) == length(sb$lower))
    b$lower[i_sel] = sb$lower
    b$upper[i_sel] = sb$upper
    o2 = opal_update(o2, bounds = list(lower = b$lower, upper = b$upper)) %>% opal_build()
    fit_o2 = if (refit) fit_opal(o2) else saved$O2
    if (refit) saveRDS(list(O1 = fit_o1, O2 = fit_o2), fits_file, compress = "xz")

#_____________________________________________________________________________
# Results
    rep_ss3 = converted$ss3$rep
    n = converted$data$n_year
    yrs = cfg$time$first_ts:cfg$time$last_ts
    ts_ss3 = as.data.table(rep_ss3$timeseries)[Era == "TIME" & Yr %in% yrs]
    ssb0_ss3 = ss3_derived_value(rep_ss3, c("SSB_Virgin", "SSB_Unfished"))

    models = list(O0 = o0, O1 = fit_o1, O2 = fit_o2)
    trajectory = rbind(
        data.table(model = "SS3 B2e", year = yrs, ssb = ts_ss3$SpawnBio,
                   depletion = ts_ss3$SpawnBio / ssb0_ss3,
                   recruits = as.data.table(rep_ss3$natage)[Era == "TIME" & `Beg/Mid` == "B" &
                                                                Yr %in% yrs][["0"]]),
        rbindlist(lapply(names(models), function(id) {
            r = opal_report(models[[id]])
            data.table(model = id, year = yrs, ssb = r$spawning_biomass_y[seq_len(n)],
                       depletion = r$static_depletion_y[seq_len(n)],
                       recruits = r$number_ysa[seq_len(n), 1, 1])
        })))

    summary_dt = rbindlist(lapply(names(models), function(id) {
        x = models[[id]]
        r = opal_report(x)
        fitted = !is.null(x$fit$opt)
        data.table(
            model        = id,
            n_par        = if (fitted) length(x$fit$opt$par) else 0L,
            objective    = if (fitted) x$fit$opt$objective else as.numeric(opal_rtmb(x, fresh = TRUE)$fn()),
            max_gradient = if (fitted) attr(x, "max_gradient") else NA_real_,
            check_passes = if (fitted) isTRUE(x$validation$fit$passes) else NA,
            positive_hessian = if (fitted) isTRUE(x$validation$fit$metrics$positive_hessian) else NA,
            seconds      = if (fitted) attr(x, "seconds") else NA_real_,
            B0           = exp((if (fitted) x$fit$opt$par[["log_B0"]] else x$parameters$log_B0)),
            ssb_2020     = r$spawning_biomass_y[n],
            depletion_2020 = r$static_depletion_y[n],
            lp_cpue      = sum(r$lp_cpue),
            lp_lf        = sum(r$lp_lf),
            lp_rec       = sum(r$lp_rec)
        )
    }))
    summary_dt = rbind(data.table(model = "SS3 B2e", n_par = rep_ss3$N_estimated_parameters,
                                  B0 = ssb0_ss3, ssb_2020 = last(ts_ss3$SpawnBio),
                                  depletion_2020 = last(ts_ss3$SpawnBio) / ssb0_ss3),
                       summary_dt, fill = TRUE)
    print(summary_dt)

    # Index fits
    cpue = rbindlist(lapply(c("O1", "O2"), function(id) {
        d = as.data.table(models[[id]]$data$cpue_data)
        d[, .(model = id, year, index, obs = value, se, pred = opal_report(models[[id]])$cpue_pred)]
    }))
    cpue_ss3 = as.data.table(rep_ss3$cpue)[Fleet %in% cfg$index$keep_ss3_fleets,
                                           .(model = "SS3 B2e", year = Yr,
                                             index = match(Fleet, cfg$index$keep_ss3_fleets),
                                             obs = Obs, se = SE, pred = Exp)]
    cpue = rbind(cpue_ss3, cpue)

    # Selectivity-at-age: SS3 against opal (O2)
    sel_ss3 = converted$selectivity$sel_fa_external
    sel_o2  = opal_report(fit_o2)$sel_fya[, n, ]
    selectivity = rbind(
        data.table(model = "SS3 B2e", fishery = rep(seq_len(nrow(sel_ss3)), times = ncol(sel_ss3)),
                   age = rep(converted$data$age_a, each = nrow(sel_ss3)), sel = as.vector(sel_ss3)),
        data.table(model = "O2", fishery = rep(seq_len(nrow(sel_o2)), times = ncol(sel_o2)),
                   age = rep(converted$data$age_a, each = nrow(sel_o2)), sel = as.vector(sel_o2)))
    # opal's starting age selectivity (its PLA) against SS3's, before fitting
    sel_o2_start = opal_report(o2)$sel_fya[, n, ]
    sel_start_diff = data.table(fishery = seq_len(nrow(sel_ss3)),
                                name = converted$meta$dims$.fleet_names[ss3_fleet_f],
                                max_abs_diff = apply(abs(sel_o2_start - sel_ss3), 1, max))

    cache = list(
        created     = format(Sys.time(), tz = "UTC", usetz = TRUE),
        session     = list(r = as.character(getRversion()),
                           opal = as.character(packageVersion("opal")),
                           r4ss = as.character(packageVersion("r4ss"))),
        fleets      = data.table(fishery = seq_along(ss3_fleet_f), ss3_fleet = ss3_fleet_f,
                                 name = converted$meta$dims$.fleet_names[ss3_fleet_f],
                                 units = converted$data$catch_units_f),
        checks      = checks,
        recon       = recon,
        growth_fit  = converted$meta$growth_fit[c("k", "CV1", "CV2", "rmse", "max_abs_error",
                                                  "laa_ss3", "laa_opal", "sd_ss3", "sd_opal")],
        bridge      = bridge,
        seed_check  = seed_check,
        sel_start_diff = sel_start_diff,
        summary     = summary_dt,
        trajectory  = trajectory,
        cpue        = cpue,
        selectivity = selectivity
    )
    saveRDS(cache, file.path(dir_cache, "npo-bsh-opal-cache.rds"), compress = "xz")
