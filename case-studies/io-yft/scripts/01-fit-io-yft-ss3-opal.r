# Nicholas Ducharme-Barth
# 2026/10/13
# Fit the converted IO yellowfin model in opal on both biology routes and
# cache everything the vignette needs, alongside the matching SS3 quantities.
#
# Both routes are fitted in the ASPM-R configuration: SS3's selectivity-at-age
# is injected as sel_fa_external and held fixed, and opal estimates log_B0, the
# 256 recruitment deviations and the four quarterly catchabilities. Every fit
# starts from the converted SS3 optimum (SS3's unfished SSB, recruitment
# deviations and bias ramp), so the timings below are for refining a good
# solution, not for finding one from scratch.
#
# Two diagnostic variants of the endgrowth route are fitted as well, to
# attribute the difference from SS3: one with the length-composition
# likelihood switched off, and one with the growth parameters (and so the PLA
# that maps catch-at-age to predicted length compositions) estimated.
#
# Two further endgrowth variants test the growth explanation directly: opal's
# PLA is built from SS3's own length-at-age (start-of-timestep Len_Beg/SD_Beg,
# and mid-timestep Len_Mid/SD_Mid) instead of the Schnute curve. opal has no
# option for this, so these use a LOCAL TEST PATCH to opal_model(), applied in
# this R session only by ss3_patch_opal_external_laa() (ss3-opal-obj.R). They
# run last; the patch is inert for any model without the external vectors.
#
# The vignette (case-studies/io-yft/vignettes/ss3-to-opal-io-yft.qmd) only reads the
# cache written here, so it rebuilds in seconds. Rerun this script when the
# converter, the config or the SS3 files change.
#
# OUTPUT:
#   case-studies/io-yft/vignettes/data/io-yft-ss3-opal-cache.rds
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
    proj_dir            = this.path::this.proj()
    dir_helper_fns_rtmb = file.path(proj_dir, "shared", "ss3-to-opal")
    dir_config          = file.path(proj_dir, "case-studies", "io-yft", "config")
    dir_cache           = file.path(proj_dir, "case-studies", "io-yft", "vignettes", "data")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)

    invisible(sapply(
        list.files(dir_helper_fns_rtmb, pattern = "\\.[Rr]$", full.names = TRUE),
        source
    ))

    cfg_path = file.path(dir_config, "io-yft-sim-iattc-sample-tree-4.yaml")

#_____________________________________________________________________________
# Helpers
    elapsed = function(expr) {
        t0 = proc.time()[["elapsed"]]
        value = force(expr)
        list(value = value, seconds = proc.time()[["elapsed"]] - t0)
    }

    # Newton refinement from the nlminb optimum. nlminb stops on relative
    # function change, which for an objective near 1.7e4 leaves gradients of
    # order 1e-3; two Newton steps take them below 1e-8.
    newton_polish = function(object, par, max_steps = 5L, tol = 1e-8) {
        for (i in seq_len(max_steps)) {
            g = as.numeric(object$gr(par))
            if (max(abs(g)) < tol) break
            H = if (is.function(object$he)) object$he(par) else
                stats::optimHess(par, object$fn, object$gr)
            par = par - as.numeric(solve(H, g))
        }
        list(par = par, max_gradient = max(abs(as.numeric(object$gr(par)))),
             objective = object$fn(par), steps = i - 1L)
    }

    lf_long = function(data, lf_pred) {
        rbindlist(lapply(seq_along(data$lf_fishery_f), function(j) {
            f    = data$lf_fishery_f[j]
            obs  = data$lf_obs_data$obs[[j]]
            obs  = obs / rowSums(obs)
            pred = lf_pred[[j]]
            ts   = data$lf_year_fi[[j]]
            n    = data$lf_n_fi[[j]]
            data.table(fishery = f,
                       ts      = rep(ts, times = ncol(obs)),
                       n       = rep(n,  times = ncol(obs)),
                       bin     = rep(data$len_lower, each = nrow(obs)),
                       obs     = as.vector(obs),
                       pred    = as.vector(pred))
        }))
    }

#_____________________________________________________________________________
# Fit one biology route
    fit_route = function(route, variant = c("base", "no_lf", "free_growth",
                                            "laa_beg", "laa_mid")) {

        variant = match.arg(variant)
        id = switch(variant, base = route, no_lf = paste0(route, "_no_lf"),
                    free_growth = paste0(route, "_free_growth"),
                    laa_beg = paste0(route, "_laa_beg"),
                    laa_mid = paste0(route, "_laa_mid"))
        cat(sprintf("\n==== %s ====\n", id))
        cfg = ss3_read_config(cfg_path)
        cfg$bio$source = route

        conv = elapsed(ss3_to_opal(cfg, dir_root = proj_dir))
        converted = conv$value

        lf_wide = dcast(converted$lf, fishery + year + month + ts ~ bin,
                        value.var = "value", fill = 0) %>%
                  setorder(fishery, ts) %>%
                  as.data.frame()
        converted$data = prep_lf_data(
            converted$data, lf_wide,
            lf_switch         = 1L,
            lf_keep_fisheries = as.integer(cfg$lf$keep_fisheries),
            lf_var_adjust     = converted$meta$lf_var_adjust,
            lf_cap            = cfg$lf$cap,
            lf_addtocomp      = as.numeric(cfg$lf$addtocomp)
        ) %>% ss3_make_portable()

        build = elapsed({
            a = ss3_opal_obj(converted, cfg,
                             control = list(eval.max = 10000, iter.max = 10000),
                             build = TRUE)
            ss3_use_external_selectivity(a, converted$selectivity$sel_fa_external) %>%
                opal_build()
        })
        bridge = build$value
        if (variant == "no_lf") {
            d0 = bridge$data
            d0$lf_switch = 0L
            bridge = opal_build(opal_update(bridge, data = d0))
        }
        if (variant == "free_growth") {
            map = bridge$map
            for (p in c("log_L1", "log_L2", "log_k", "log_CV1", "log_CV2")) map[[p]] = NULL
            bridge = opal_build(opal_update(bridge, map = map))
        }
        laa = NULL
        if (variant %in% c("laa_beg", "laa_mid")) {
            ss3_patch_opal_external_laa()
            laa = ss3_endgrowth_table(converted$ss3$rep, converted$meta$dims,
                                      basis = if (variant == "laa_beg") "Beg" else "Mid")
            bridge = opal_build(ss3_use_external_length_at_age(bridge, laa$length, laa$sd))
        }
        start_object = opal_rtmb(bridge, fresh = TRUE)
        start_par    = start_object$par
        rep_start    = start_object$report(start_par)

        # opal_fit() rebuilds the RTMB tape, then runs the nlminb passes.
        fit = elapsed(opal_fit(bridge, n_passes = 2L, check = FALSE))
        passes = fit$value$fit$diagnostics$optimisation$passes

        tape = elapsed(opal_rtmb(fit$value, fresh = TRUE))
        object = tape$value
        polish = elapsed(newton_polish(object, fit$value$fit$opt$par))

        opt = fit$value$fit$opt
        opt$par = stats::setNames(polish$value$par, names(opt$par))
        opt$objective = polish$value$objective
        check = elapsed(opal_attach_fit(fit$value, opt, check = TRUE,
                                        check_args = list(gradient_tolerance = 1e-4)))
        fitted = check$value
        rep_fit = object$report(opt$par)
        pl = object$env$parList(opt$par)

        cat(sprintf("  objective %.3f -> %.3f, max|grad| %.2e, check %s\n",
                    start_object$fn(start_par), opt$objective,
                    polish$value$max_gradient, fitted$validation$fit$passes))

        d = converted$data
        n = d$n_year
        traj = function(r, what) data.table(
            route       = route,
            state       = what,
            ts          = seq_len(n),
            ssb         = r$spawning_biomass_y[seq_len(n)],
            ssb0        = r$spawning_biomass0_y[seq_len(n)],
            static_dep  = r$static_depletion_y[seq_len(n)],
            dynamic_dep = r$dynamic_depletion_y[seq_len(n)],
            recruits    = r$number_ysa[seq_len(n), 1, 1]
        )

        route = id
        cpue = as.data.table(d$cpue_data)[, .(route, ts, fishery, index, obs = value, se)]
        cpue[, pred_fit   := rep_fit$cpue_pred]
        cpue[, pred_start := rep_start$cpue_pred]
        cpue[, sigma      := rep_fit$cpue_sigma]

        # Predicted compositions are reported only when the LF likelihood is on.
        lf = if (is.null(rep_fit$lf_pred)) NULL else lf_long(d, rep_fit$lf_pred)[, route := route]

        mu_a = with(pl, exp(log_L1) + (exp(log_L2) - exp(log_L1)) *
                        (1 - exp(-exp(log_k) * (d$age_a - d$A1))) /
                        (1 - exp(-exp(log_k) * (d$A2 - d$A1))))
        if (!is.null(laa)) mu_a = as.numeric(laa$length)

        list(
            route      = route,
            timing     = data.table(
                route = route,
                stage = c("convert SS3 files", "build opal objects (two tapes)",
                          "opal_fit: tape + nlminb (2 passes)", "rebuild tape",
                          "Newton polish", "attach fit + opal_check"),
                seconds = c(conv$seconds, build$seconds, fit$seconds,
                            tape$seconds, polish$seconds, check$seconds)),
            optimiser  = data.table(
                route = route,
                pass  = seq_along(passes),
                objective   = vapply(passes, `[[`, numeric(1), "objective"),
                convergence = vapply(passes, `[[`, numeric(1), "convergence"),
                iterations  = vapply(passes, function(p) p$iterations, numeric(1)),
                fn_evals    = vapply(passes, function(p) p$evaluations[["function"]], numeric(1)),
                gr_evals    = vapply(passes, function(p) p$evaluations[["gradient"]], numeric(1))),
            summary    = data.table(
                route              = route,
                n_par              = length(opt$par),
                objective_start    = start_object$fn(start_par),
                objective_fit      = opt$objective,
                newton_steps       = polish$value$steps,
                max_gradient_nlminb = max(abs(as.numeric(object$gr(fit$value$fit$opt$par)))),
                max_gradient       = polish$value$max_gradient,
                check_passes       = isTRUE(fitted$validation$fit$passes),
                estimability       = fitted$fit$estimability$message %||% NA_character_,
                B0_start           = exp(start_par[["log_B0"]]),
                B0_fit             = exp(pl$log_B0),
                R0_start           = rep_start$R0,
                R0_fit             = rep_fit$R0,
                lp_cpue            = sum(rep_fit$lp_cpue),
                lp_lf              = sum(rep_fit$lp_lf),
                lp_rec             = rep_fit$lp_rec,
                lp_penalty         = rep_fit$lp_penalty),
            parameters = data.table(route = route, ts = seq_len(n),
                                    rdev_start = start_object$env$parList(start_par)$rdev_y,
                                    rdev_fit   = pl$rdev_y),
            q          = data.table(route = route, index = seq_along(pl$log_cpue_q),
                                    log_q = pl$log_cpue_q),
            trajectory = rbind(traj(rep_start, "SS3 values"), traj(rep_fit, "fitted")),
            cpue       = cpue,
            lf         = lf,
            biology    = data.table(route = route, age = d$age_a,
                                    M = rep_fit$M_a,
                                    weight = rep_fit$weight_fya_mod[1, 1, ],
                                    spawning_potential = rep_fit$spawning_potential_a,
                                    length_opal = mu_a),
            reconcile  = ss3_reconcile(converted)[, route := route],
            converted_meta = list(growth_fit = converted$meta$growth_fit,
                                  var_adjust = converted$meta$lf_var_adjust,
                                  fleet_names = converted$meta$dims$.fleet_names),
            ss3        = converted$ss3
        )
    }

#_____________________________________________________________________________
# Run both routes
    routes = list(
        endgrowth             = fit_route("endgrowth"),
        control               = fit_route("control"),
        endgrowth_no_lf       = fit_route("endgrowth", "no_lf"),
        endgrowth_free_growth = fit_route("endgrowth", "free_growth"),
        # Local opal patch from here on; see the header.
        endgrowth_laa_beg     = fit_route("endgrowth", "laa_beg"),
        endgrowth_laa_mid     = fit_route("endgrowth", "laa_mid")
    )

#_____________________________________________________________________________
# SS3 quantities for comparison
    ss3 = routes$endgrowth$ss3
    rep = ss3$rep
    cfg = ss3_read_config(cfg_path)
    dims = ss3_dims_to_opal(ss3, cfg)

    ts_ss3 = as.data.table(rep$timeseries)[Era == "TIME"]
    ssb0_ss3 = ss3_derived_value(rep, c("SSB_Virgin", "SSB_Unfished"))
    ss3_traj = ts_ss3[, .(ts = Yr, ssb = SpawnBio, static_dep = SpawnBio / ssb0_ss3,
                          recruits_age0 = Recruit_0)]

    # SS3's CPUE fit, carved into the same quarterly series and on the same
    # mean-one scale as the opal indices so the two can share an axis.
    ss3_cpue = as.data.table(rep$cpue)[, .(ts = Yr, ss3_fleet = Fleet, obs = Obs, exp = Exp, se = SE)]
    ss3_cpue[, fishery := ss3_apply_index_split(ss3_fleet, ts, cfg)]
    ss3_cpue[, `:=`(exp = exp / mean(obs), obs = obs / mean(obs)), by = fishery]

    # SS3's expected length compositions are in the expected-value section of
    # data.ss_new (CompReport.sso was not written for this run).
    dir_ss3 = file.path(proj_dir, cfg$ss3_dir)
    dat_exp = suppressWarnings(r4ss::SS_readdat_3.24(
        file.path(dir_ss3, cfg$files$dat), section = 2, verbose = FALSE))
    lc = as.data.table(dat_exp$lencomp)
    bin_cols = grep("^l[0-9]", names(lc), value = TRUE)
    lc_long = melt(lc[, c("Yr", "FltSvy", bin_cols), with = FALSE],
                   id.vars = c("Yr", "FltSvy"), variable.name = "bin_col", value.name = "exp")
    lc_long[, bin := dims$.len_bins[match(as.numeric(sub("^l", "", bin_col)), dims$.len_bins_ss3)]]
    lc_long[, exp := exp / sum(exp), by = .(Yr, FltSvy)]
    ss3_lf = lc_long[, .(ts = Yr, fishery = FltSvy, bin, ss3_pred = exp)]

    # Run time of the SS3 model, from the Report.sso header.
    hdr = readLines(file.path(dir_ss3, "Report.sso"), n = 6)
    stamp = function(tag) {
        x = sub(paste0("^", tag, ":\\s*"), "", grep(paste0("^", tag), hdr, value = TRUE))
        as.POSIXct(x, format = "%a %b %d %H:%M:%S %Y", tz = "UTC")
    }
    ss3_seconds = as.numeric(difftime(stamp("EndTime"), stamp("StartTime"), units = "secs"))

    lik = rep$likelihoods_used

    ss3_out = list(
        trajectory   = ss3_traj,
        ssb0         = ssb0_ss3,
        r0           = ss3_derived_value(rep, c("Recr_Virgin", "Recr_Unfished")),
        cpue         = ss3_cpue,
        lf           = ss3_lf,
        endgrowth    = as.data.table(rep$endgrowth)[, .(age = Age, length = Len_Beg,
                                                        weight = Wt_Beg, M = M,
                                                        mat_fec = `Mat*Fecund`)],
        recdevs      = as.data.table(ss3$par$recdev1),
        likelihoods  = data.table(component = rownames(lik), value = lik$values),
        run_seconds  = ss3_seconds,
        n_par        = as.numeric(sub(".*Number of parameters = ([0-9]+).*", "\\1",
                                      readLines(file.path(dir_ss3, cfg$files$par), n = 1))),
        f_method     = ss3$ctl$F_Method,
        version      = ss3$dat$SSversion %||% "3.24AB"
    )

#_____________________________________________________________________________
# Cache
    stack = function(what) rbindlist(lapply(routes, `[[`, what), use.names = TRUE, fill = TRUE)

    cache = list(
        created     = format(Sys.time(), tz = "UTC", usetz = TRUE),
        session     = list(r = as.character(getRversion()),
                           opal = as.character(packageVersion("opal")),
                           r4ss = as.character(packageVersion("r4ss")),
                           cpu = tryCatch(trimws(sub(".*: ", "", grep("model name",
                               readLines("/proc/cpuinfo"), value = TRUE)[1])),
                               error = function(e) Sys.info()[["machine"]]),
                           cores = parallel::detectCores()),
        config_md5  = unname(tools::md5sum(cfg_path)),
        fleet_names = routes$endgrowth$converted_meta$fleet_names,
        growth_fit  = lapply(routes, function(r) r$converted_meta$growth_fit),
        timing      = stack("timing"),
        optimiser   = stack("optimiser"),
        summary     = stack("summary"),
        parameters  = stack("parameters"),
        q           = stack("q"),
        trajectory  = stack("trajectory"),
        cpue        = stack("cpue"),
        lf          = stack("lf"),
        biology     = stack("biology"),
        reconcile   = stack("reconcile"),
        ss3         = ss3_out
    )

    saveRDS(cache, file.path(dir_cache, "io-yft-ss3-opal-cache.rds"), compress = "xz")
    cat(sprintf("\nCache written: %s (%.0f kB)\n",
                file.path(dir_cache, "io-yft-ss3-opal-cache.rds"),
                file.size(file.path(dir_cache, "io-yft-ss3-opal-cache.rds")) / 1024))
    print(cache$summary)
    print(cache$timing)
