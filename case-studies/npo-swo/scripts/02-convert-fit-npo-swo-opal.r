# Nicholas Ducharme-Barth
# 2026/11/13
# Convert the single-sex, annual SS3 bridge model B7 of the 2023 ISC North
# Pacific swordfish assessment to opal, test the conversion, replay SS3's
# dynamics with opal's catch timing, and fit opal.
#
#   O0  bridge: opal at SS3's parameter values with SS3's selectivity-at-age;
#       nothing estimated
#   O1  SS3's selectivity-at-age held fixed; B0, recruitment deviations,
#       initial age structure, F6's initial F and the seven catchabilities
#       estimated (as SS3, which floats q and estimates the rest)
#   O2  opal's own selectivity curves (length-based double normals and
#       logistics), started exactly at SS3's curves; the selectivity
#       parameters SS3 estimated are estimated, plus everything O1 estimates
#   O2b O2 with the selectivity parameters that reach a bound in O2 fixed
#       there
#
# The recruitment index S6 (SS3 units 33: the year's age-0 recruitment) has no
# opal equivalent. It is an index in numbers whose selectivity is 1 at age 0
# and 0 at older ages: a row of SS3's selectivity-at-age in O0 and O1, and in
# O2 an age-based (sel_basis_f = 2) double normal fixed at that shape.
#
# Growth, steepness (0.9), sigma_R (0.42), M and the biology stay at SS3's
# values, as in SS3. opal's default priors are switched off (SS3's prior
# lambda is 0). The length compositions carry SS3's lambda of 0.5.
#
# Diagnostics: O1-noLF (O1 without the compositions) and O1-noS6 (O1 without
# the recruitment index; the report's sensitivity 15).
#
# Phase D: the same opal models with four seasons (quarters) per year and
# each fleet's catch by quarter (from the quarterly SS3 model B3 of the same
# chain, whose catch is B7's before the quarters were summed). opal then
# removes a quarter of the year's M and that quarter's catch in each season.
# Everything else is B7's: spawning and recruitment at the start of season 1,
# the indices observed in season 1, compositions predicted from the year's
# catch-at-age, one weight-at-age for the year (SS3's mid-year weights, see
# below).
#   D0       O0 with quarterly catch (SS3's values)
#   D0-even  O0 with each year's catch split evenly over the quarters
#   D        O2b (the final opal configuration) with quarterly catch,
#            refitted
#
# Two variants (NPO_SWO_VARIANT): "female" (default; B7, female biology) and
# "sexavg" (B7s, sex-averaged biology). The fits are cached; NPO_SWO_REFIT=true
# refits.
#
# OUTPUTS (<sfx> is "" for female, "-sexavg" for sexavg):
#   case-studies/npo-swo/vignettes/data/npo-swo-opal-cache<sfx>.rds
#   case-studies/npo-swo/vignettes/data/npo-swo-opal-fits<sfx>.rds  (fitted opal_obj's)
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
    dir_cfg   = file.path(proj_dir, "case-studies", "npo-swo", "config")
    dir_cache = file.path(proj_dir, "case-studies", "npo-swo", "vignettes", "data")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
    invisible(sapply(list.files(dir_fns, pattern = "\\.[Rr]$", full.names = TRUE), source))

#_____________________________________________________________________________
# Variant
    variant = Sys.getenv("NPO_SWO_VARIANT", "female")
    stopifnot(variant %in% c("female", "sexavg"))
    sfx  = if (variant == "female") "" else "-sexavg"
    ssfx = if (variant == "female") "" else "s"       # SS3 step suffix
    ss3_label = paste0("SS3 B7", ssfx)
    quarterly_step = paste0("03", ssfx, "-spawn-month-1")

#_____________________________________________________________________________
# Convert
    cfg = ss3_read_config(file.path(dir_cfg, paste0("npo-swo-b7", ssfx, ".yaml")))
    if (nzchar(Sys.getenv("NPO_SWO_SS3_DIR"))) {
        # development only: convert a copy of the SS3 model held elsewhere
        cfg$absolute_paths = TRUE
        cfg$ss3_dir = Sys.getenv("NPO_SWO_SS3_DIR")
    }
    converted = ss3_to_opal(cfg, dir_root = proj_dir)

    # opal's growth curve: fitted to SS3's mid-year length-at-age. Growth only
    # feeds opal's PLA here (the biology is SS3's start-of-year endgrowth
    # table), which maps length-based selectivity to age and predicts the
    # compositions; SS3 does both with mid-year lengths (sample timing -1: the
    # fishery's catch over the year). The start-of-year lengths also cannot be
    # fitted: SS3 grows age-0 fish linearly from 10 cm at settlement to L1 at
    # age 0.5, and the Schnute fit misses them by up to 42 cm (RMSE 11.7 cm);
    # at mid-year the RMSE is 0.45 cm.
    cfg_mid = cfg
    cfg_mid$bio$endgrowth$length_basis = "Mid"
    converted_mid = ss3_to_opal(cfg_mid, dir_root = proj_dir)
    growth_beg = converted$meta$growth_fit
    for (nm in c("log_L1", "log_L2", "log_k", "log_CV1", "log_CV2")) {
        converted$parameters[[nm]] = converted_mid$parameters[[nm]]
    }
    converted$meta$growth_fit = converted_mid$meta$growth_fit
    d = converted$data
    rep_ss3 = converted$ss3$rep
    fish = converted$meta$fisheries
    ss3_fleet_f = c(fish$extract_fleets, fish$index_map$ss3_fleet)
    n_f = d$n_fishery
    n_age = d$n_age
    n_year = d$n_year
    yrs = cfg$time$first_ts:cfg$time$last_ts

    # S6, the recruitment index: numbers at age 0. The converter reads SS3's
    # units code 33 as biomass (it knows 0 = numbers and 1 = biomass) and
    # SS3's selectivity-at-age for S6 (size pattern 0) as 1 at every age.
    f_s6 = fish$index_map[ss3_fleet == as.integer(cfg$index$recruitment_fleet), opal_fleet]
    stopifnot(length(f_s6) == 1L)
    d$catch_units_f[f_s6] = 2
    d$cpue_data$units[d$cpue_data$fishery == f_s6] = 2
    sel_ss3 = converted$selectivity$sel_fa_external
    sel_ss3[f_s6, ] = c(1, rep(0, n_age - 1L))
    converted$selectivity$sel_fa_external = sel_ss3

    # Length compositions (multinomial, as SS3) and SS3's lambdas: 0.5 on the
    # compositions, 0 on parameter priors
    lf_wide = dcast(converted$lf, fishery + year + month + ts ~ bin,
                    value.var = "value", fill = 0) %>%
              setorder(fishery, ts) %>%
              as.data.frame()
    d = prep_lf_data(
        d, lf_wide,
        lf_switch         = 1L,
        lf_keep_fisheries = as.integer(cfg$lf$keep_fisheries),
        lf_var_adjust     = converted$meta$lf_var_adjust,
        lf_cap            = cfg$lf$cap,
        lf_addtocomp      = as.numeric(cfg$lf$addtocomp)
    ) %>% ss3_make_portable()
    lf_lambda = ifelse(seq_len(n_f) %in% cfg$lf$keep_fisheries, 0.5, 1)
    d$lambdas = list(prior = 0, lf = lf_lambda)
    converted$data = d

#_____________________________________________________________________________
# Test the conversion
    checks = validate_opal_inputs(converted, strict = FALSE)
    recon  = ss3_reconcile(converted, tol = 1e-6)
    stopifnot("Structural validation failed" = all(checks$pass),
              "Reconciliation with SS3 failed" = all(recon$pass, na.rm = TRUE))

    # SS3's own selectivity-at-age for the split fleets and the mirrors
    # (the converter takes each fleet's row of ageselex in 2021)
    asel = as.data.table(rep_ss3$ageselex)[Factor == "Asel2" & Yr == 2021]

    assessment = ss3_opal_obj(converted, cfg,
                              control = list(eval.max = 10000, iter.max = 10000),
                              build = FALSE)

#_____________________________________________________________________________
# Bounds: opal's defaults, with catchability widened (the indices are not
# normalised) and the estimated selectivity parameters at SS3's bounds
    build_with_bounds = function(x, sel = NULL) {
        x = opal_update(x, bounds = NULL)
        p = x$parameters
        tmp = p
        tmp$log_cpue_q = rep(0, length(p$log_cpue_q))
        tmp$par_sel = p$par_sel * 0
        x1 = opal_build(opal_update(x, parameters = tmp))
        b = x1$bounds
        nm = names(opal_rtmb(x1, fresh = TRUE)$par)
        i = which(nm == "log_cpue_q")
        b$lower[i] = -15
        b$upper[i] = 15
        if (!is.null(sel)) {
            i = which(nm == "par_sel")
            stopifnot(length(i) == length(sel$lower))
            b$lower[i] = sel$lower
            b$upper[i] = sel$upper
        }
        opal_update(x, bounds = list(lower = b$lower, upper = b$upper)) %>% opal_build()
    }

#_____________________________________________________________________________
# O0: bridge at SS3's values
    o0 = ss3_use_external_selectivity(assessment, sel_ss3) %>% build_with_bounds()
    rep_o0 = opal_report(o0)
    bridge = ss3_bridge_table(converted, rep_o0)
    print(bridge[, .(ts, ss3, opal, rel_opal)])

    # Initial numbers-at-age: SS3 (start of 1975) against opal (O0). SS3 solves
    # the initial equilibrium with F6's initial F applied through the year
    # (Baranov); opal removes the equilibrium catch at the start of the year.
    natage = as.data.table(rep_ss3$natage)
    age_cols = as.character(converted$meta$dims$.ss3_ages)
    n_ss3 = function(yr) as.numeric(natage[Era == "TIME" & Yr == yr & `Beg/Mid` == "B", ..age_cols][1])
    init_numbers = data.table(age = d$age_a, ss3 = n_ss3(yrs[1]), opal = rep_o0$number_ysa[1, 1, ])
    init_numbers[, rel := opal / ss3 - 1]

#_____________________________________________________________________________
# Catch timing and catch weight at SS3's values. SS3's numbers-at-age at the
# start of 1975 and SS3's recruitment each year are replayed with SS3's
# selectivity-at-age, removing each year's catch (a) at the start of the year
# (opal), (b) at mid-year (Pope: half the year's M, the catch, the other
# half) or (c) through the year (Baranov, as SS3's hybrid F), with catch in
# weight turned into numbers at start-of-year weight (opal here), mid-year
# weight or each fleet's SS3 catch weight-at-age (ageselex "bodywt").
    recruits_ss3 = vapply(yrs, function(y) n_ss3(y)[1], numeric(1))
    replay_timing = function(timing = c("start", "pope", "baranov"), W) {
        timing = match.arg(timing)
        sp = d$spawning_potential
        M  = d$M
        S_a = exp(-M)
        num = d$catch_units_f != 1
        cw = function(f) if (num[f]) 1 else W[f, ]
        N = n_ss3(yrs[1])
        ssb = numeric(n_year)
        for (y in seq_len(n_year)) {
            N[1] = recruits_ss3[y]
            ssb[y] = sum(N * sp)
            cy = d$catch_obs_ysf[y, 1, ]
            fl = which(cy > 0)
            if (timing == "baranov") {
                F_f = numeric(length(cy))
                for (f in fl) F_f[f] = cy[f] / sum(N * sel_ss3[f, ] * cw(f))
                for (it in 1:50) {
                    Z = M + colSums(F_f * sel_ss3)
                    for (f in fl) F_f[f] = F_f[f] * cy[f] /
                        sum(F_f[f] * sel_ss3[f, ] / Z * N * (1 - exp(-Z)) * cw(f))
                }
                N_end = N * exp(-(M + colSums(F_f * sel_ss3)))
            } else {
                Nc = if (timing == "pope") N * exp(-M / 2) else N
                hr = numeric(n_age)
                for (f in fl) hr = hr + cy[f] / sum(Nc * sel_ss3[f, ] * cw(f)) * sel_ss3[f, ]
                hr = pmin(hr, 0.999)
                N_end = if (timing == "pope") Nc * (1 - hr) * exp(-M / 2) else N * (1 - hr) * S_a
            }
            N_next = c(0, N_end[-n_age])
            N_next[n_age] = N_next[n_age] + N_end[n_age]
            N = N_next
        }
        ssb
    }
    eg_beg = ss3_endgrowth_table(rep_ss3, converted$meta$dims, "Beg")
    eg_mid = ss3_endgrowth_table(rep_ss3, converted$meta$dims, "Mid")
    bodywt = as.data.table(rep_ss3$ageselex)[Factor == "bodywt" & Yr == 2021]
    W_fleet = as.matrix(bodywt[match(ss3_fleet_f, Fleet), ..age_cols])
    W_fleet[is.na(W_fleet)] = 0
    W_beg = matrix(eg_beg$weight, n_f, n_age, byrow = TRUE)
    W_mid = matrix(eg_mid$weight, n_f, n_age, byrow = TRUE)
    ssb_ss3 = as.numeric(bridge$ss3)
    timing = rbindlist(lapply(list(
        list("start", W_beg,   "start of year (opal)", "start of year (opal)"),
        list("start", W_fleet, "start of year (opal)", "SS3 fleet catch weight"),
        list("pope",  W_fleet, "mid-year (Pope)", "SS3 fleet catch weight"),
        list("baranov", W_mid, "through the year (Baranov, SS3)", "mid-year"),
        list("baranov", W_fleet, "through the year (Baranov, SS3)", "SS3 fleet catch weight")),
        function(a) data.table(year = yrs, removal = a[[3]], weight = a[[4]],
                               ssb = replay_timing(a[[1]], a[[2]]), ss3 = ssb_ss3)))
    timing[, rel := ssb / ss3 - 1]
    print(timing[, .(rel_end = rel[.N], rel_max = max(abs(rel))), by = .(removal, weight)])

#_____________________________________________________________________________
# Selectivity: opal's curves at SS3's
    # Double normals (SS3 pattern 24): opal's double normal at SS3's curve
    seed = ss3_seed_par_sel_size_dblnorm(converted$ss3, fish, d$len_mid,
                                         converted$parameters$par_sel)
    par_sel0 = seed$par_sel
    source_fleet = seed$source_fleet
    # Logistics (SS3 pattern 1): the same curve, re-scaled
    types = as.data.frame(converted$ss3$ctl$size_selex_types)
    pars_ss3 = as.data.table(rep_ss3$parameters)
    resolve = function(sf) {
        while (types$Pattern[sf] %in% c(5L, 15L)) sf = as.integer(types$Special[sf])
        sf
    }
    for (i in seq_len(n_f)) {
        sf = resolve(ss3_fleet_f[i])
        if (types$Pattern[sf] != 1L) next
        p = pars_ss3[grepl("^Size_(inflection|95%width)_", Label) &
                         grepl(sprintf("\\(%d\\)$", sf), Label)]$Value
        par_sel0[i, 1:2] = opal::convert_ss3_selex_to_rtmb(matrix(c(p, 0, 0, 0, 0), 1L), 1L, d$len_mid)[1, 1:2]
        source_fleet[i] = sf
    }
    # S6: age-based double normal with selectivity 1 at age 0 and < 1e-4
    # from age 1 (peak at age -1; the plateau ends one age later, at age 0;
    # steep limbs; selectivity at the last age 0). opal's plateau is at least
    # one length bin (here one age) wide, so the peak sits at age -1.
    ages = d$age_a
    par_sel0[f_s6, ] = c((-1 - mean(ages)) / sd(ages), -20, -15, -15, -20, -20)
    sel_basis_f = rep(1L, n_f)
    sel_basis_f[f_s6] = 2L
    sel_s6 = sel_double_normal(ages, par_sel0[f_s6, ])
    stopifnot(abs(sel_s6[1] - 1) < 1e-4, max(sel_s6[-1]) < 1e-4, min(sel_s6) >= 0)

    sel_curve = function(i, par) {
        x = if (sel_basis_f[i] == 2L) ages else d$len_mid
        if (d$sel_type_f[i] == 1L) sel_logistic(x, par) else sel_double_normal(x, par)
    }
    sizeselex = as.data.table(rep_ss3$sizeselex)[Factor == "Lsel" & Yr == 2021]
    len_cols = setdiff(names(sizeselex), c("Factor", "Fleet", "Yr", "Sex", "Label"))
    seed_check = data.table(
        fishery    = seq_len(n_f),
        ss3_fleet  = ss3_fleet_f,
        name       = converted$meta$dims$.fleet_names[ss3_fleet_f],
        source     = source_fleet,
        max_abs_diff = vapply(seq_len(n_f), function(i) {
            if (i == f_s6) return(max(abs(sel_s6 - sel_ss3[f_s6, ])))
            max(abs(sel_curve(i, par_sel0[i, ]) -
                    as.numeric(sizeselex[Fleet == ss3_fleet_f[i], ..len_cols])))
        }, numeric(1L)))
    print(seed_check)

    # Which parameters SS3 estimated, by opal fishery (double normal: peak,
    # top, ascending width, descending width, start, end; logistic:
    # inflection, width)
    est_cols = lapply(seq_len(n_f), function(i) {
        sf = source_fleet[i]
        if (is.na(sf) || i == f_s6 || i %in% as.integer(names(cfg$selectivity$mirror))) return(integer())
        if (!i %in% as.integer(cfg$selectivity$estimate)) return(integer())
        p = pars_ss3[grepl("^Size_(DblN_|inflection|95%width)", Label) & grepl(sprintf("\\(%d\\)$", sf), Label)]
        which(!is.na(p$Active_Cnt))
    })
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
    # SS3's bounds on opal's scale (double normal: peak, top, widths; logistic:
    # inflection, width)
    sel_bounds = function() {
        x = d$len_mid
        mu = mean(x)
        sd_x = stats::sd(x)
        to_opal = function(v, col, type) {
            if (type == 1L) return(switch(col, (v - mu) / sd_x, log(pmax(v, 1e-3) / sd_x)))
            switch(col, (v - mu) / sd_x, v, v - 2 * log(sd_x), v - 2 * log(sd_x))
        }
        out = list(lower = numeric(), upper = numeric())
        for (f in seq_along(est_cols)) {
            sf = source_fleet[f]
            for (col in est_cols[[f]]) {
                p = pars_ss3[grepl("^Size_(DblN_|inflection|95%width)", Label) &
                                 grepl(sprintf("\\(%d\\)$", sf), Label)][col]
                out$lower = c(out$lower, to_opal(p$Min, col, d$sel_type_f[f]))
                out$upper = c(out$upper, to_opal(p$Max, col, d$sel_type_f[f]))
            }
        }
        out
    }

#_____________________________________________________________________________
# Maps
    # Initial age structure: SS3 estimates the deviations of the cohorts born
    # 1965-1975 (ages 0-10 at the start of 1975); opal's last recruitment
    # deviation (2022) has no data and is fixed at 0. F6's initial F is
    # estimated, as in SS3 (whose initial-equilibrium catch has lambda 0).
    n_init = sum(converted$parameters$init_rdev_a != 0 | seq_len(n_age) <= 11L)
    # fisheries with an SS3 initial F (F6 here)
    f_init = match(ss3_init_F(converted$ss3)$fleet, fish$extract_fleets)
    stopifnot(length(f_init) >= 1L, !anyNA(f_init))
    base_map = function(x) {
        m = x$map %||% get_map(x$parameters)
        m$init_rdev_a = factor(c(seq_len(n_init), rep(NA, n_age - n_init)))
        m$rdev_y = factor(c(seq_len(n_year - 1L), NA))
        m$log_init_F_f = factor(ifelse(seq_len(n_f) %in% f_init, seq_along(f_init)[match(seq_len(n_f), f_init)], NA))
        m
    }

#_____________________________________________________________________________
# Fits
    fit_opal = function(x) {
        t0  = proc.time()[["elapsed"]]
        fit = opal_fit(x, n_passes = 3L, check = FALSE,
                       control = list(eval.max = 1e5, iter.max = 1e5))
        obj = opal_rtmb(fit, fresh = TRUE)
        opt = fit$fit$opt
        # Newton steps to bring the gradient down; keep the nlminb optimum if a
        # step lands where the objective is undefined or outside the bounds
        pol = tryCatch(ss3_newton_polish(obj, opt$par), error = function(e) NULL)
        b = x$bounds
        if (is.null(pol) || !is.finite(pol$objective) || pol$objective > opt$objective + 1e-8 ||
            any(pol$par < b$lower | pol$par > b$upper)) {
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

    fits_file = file.path(dir_cache, paste0("npo-swo-opal-fits", sfx, ".rds"))
    refit = !file.exists(fits_file) || identical(Sys.getenv("NPO_SWO_REFIT"), "true")
    saved = if (refit) list() else readRDS(fits_file)
    get_fit = function(id, x) {
        if (!refit && !is.null(saved[[id]])) return(saved[[id]])
        message("fitting ", id)
        fit = fit_opal(x)
        saved[[id]] <<- fit
        saveRDS(saved, fits_file, compress = "xz")
        fit
    }
    with_lambdas = function(x, lam) {
        dd = x$data
        dd$lambdas = modifyList(dd$lambdas, lam)
        opal_update(x, data = dd)
    }

    # O1: SS3 selectivity fixed
    o1 = opal_update(o0, map = base_map(o0)) %>% build_with_bounds()
    fit_o1 = get_fit("O1", o1)

    # O2: opal selectivity, started at SS3's curves and O1's other estimates
    o2 = ss3_use_opal_selectivity(assessment, par_sel = par_sel0)
    dd = o2$data
    dd$sel_basis_f = sel_basis_f
    p2 = fit_o1$fit$parameters
    p2$par_sel = par_sel0
    o2 = opal_update(o2, data = dd, parameters = p2)
    m2 = base_map(o2)
    m2$par_sel = sel_map(par_sel0)
    o2 = opal_update(o2, map = m2) %>% build_with_bounds(sel = sel_bounds())
    o2_start = o2
    fit_o2 = get_fit("O2", o2)

    # O2b: O2 refitted with the selectivity parameters it cannot determine
    # fixed at their O2 values: those that end on a bound (in O2 F2's
    # inflection reaches SS3's upper bound of 150 cm, and the descending
    # widths of F3, F6 and F8 SS3's lower bounds, where the likelihood is flat
    # because their plateaus reach the largest length), then, repeatedly,
    # the one that dominates a direction in which the Hessian is flat
    # (eigenvalue below 1e-6)
    lvl_matrix = matrix(as.integer(as.character(m2$par_sel)), n_f)
    # positions (among the estimated selectivity parameters) that lie within
    # 1e-3 of a bound on opal's scale
    sel_on_bound = function(fit, tol = 1e-3) {
        opt = fit$fit$opt
        i = which(names(opt$par) == "par_sel")
        b = fit$bounds
        which(abs(opt$par[i] - b$lower[i]) < tol | abs(opt$par[i] - b$upper[i]) < tol)
    }
    sel_lvl_at_bound = function(fit) sel_on_bound(fit)
    sel_lvl_flat = function(fit, free_lvl) {
        obj = opal_rtmb(fit, fresh = TRUE)
        opt = fit$fit$opt
        H = stats::optimHess(opt$par, obj$fn, obj$gr)
        e = eigen(H, symmetric = TRUE)
        i = which(names(opt$par) == "par_sel")
        out = integer()
        for (k in which(e$values < 1e-6)) {
            v = e$vectors[, k]
            j = which.max(abs(v))
            if (j %in% i && abs(v[j]) > 0.9) out = c(out, free_lvl[match(j, i)])
        }
        unique(out)
    }
    refit_fixing = function(fixed_lvl) {
        lv = as.integer(as.character(m2$par_sel))
        lv[lv %in% fixed_lvl] = NA
        free_lvl = sort(unique(lv[!is.na(lv)]))
        m = m2
        m$par_sel = factor(match(lv, free_lvl))
        sb = sel_bounds()
        x = opal_update(o2_start, parameters = fit_o2$fit$parameters, map = m) %>%
            build_with_bounds(sel = list(lower = sb$lower[free_lvl], upper = sb$upper[free_lvl]))
        list(x = x, free_lvl = free_lvl)
    }
    fixed_lvl = sel_lvl_at_bound(fit_o2)
    if (refit || is.null(saved$O2b)) {
        for (pass in 1:4) {
            r = refit_fixing(fixed_lvl)
            fit_o2b = fit_opal(r$x)
            flat = sel_lvl_flat(fit_o2b, r$free_lvl)
            if (!length(flat)) break
            fixed_lvl = c(fixed_lvl, flat)
        }
        attr(fit_o2b, "fixed_levels") = fixed_lvl
        saved$O2b = fit_o2b
        saveRDS(saved, fits_file, compress = "xz")
    } else {
        fit_o2b = saved$O2b
        fixed_lvl = attr(fit_o2b, "fixed_levels")
    }
    sel_fixed = rbindlist(lapply(fixed_lvl, function(k) {
        rc = which(lvl_matrix == k, arr.ind = TRUE)[1, ]
        data.table(fishery = rc[[1]], column = rc[[2]],
                   name = converted$meta$dims$.fleet_names[ss3_fleet_f[rc[[1]]]],
                   value = fit_o2$fit$parameters$par_sel[rc[[1]], rc[[2]]],
                   reason = if (k %in% sel_lvl_at_bound(fit_o2)) "bound" else "flat")
    }))
    print(sel_fixed)

    # Diagnostics
    o1_nolf = with_lambdas(o1, list(lf = 0)) %>% build_with_bounds()
    fit_o1_nolf = get_fit("O1-noLF", o1_nolf)
    lam_cpue = rep(1, d$n_index)
    lam_cpue[fish$index_map[opal_fleet == f_s6, which = TRUE]] = 0
    o1_nos6 = with_lambdas(o1, list(cpue = lam_cpue))
    m_nos6 = o1_nos6$map
    q_lvl = seq_len(d$n_index)
    q_lvl[lam_cpue == 0] = NA
    m_nos6$log_cpue_q = factor(q_lvl)
    o1_nos6 = opal_update(o1_nos6, map = m_nos6) %>% build_with_bounds()
    fit_o1_nos6 = get_fit("O1-noS6", o1_nos6)


#_____________________________________________________________________________
# Phase D: four seasons per year
    # each fleet's catch by quarter, from the quarterly SS3 model of this chain
    dat_q = SS_read(file.path(proj_dir, "case-studies", "npo-swo", "model-files", "ss3",
                              quarterly_step), ss_new = TRUE, verbose = FALSE)$dat
    cq = as.data.table(dat_q$catch)[year %in% yrs]
    cq[, f := match(fleet, fish$extract_fleets)]
    cq = cq[!is.na(f) & catch > 0]
    catch_q = array(0, dim = c(n_year, 4L, n_f))
    for (r in seq_len(nrow(cq))) catch_q[cq$year[r] - yrs[1] + 1L, cq$seas[r], cq$f[r]] = cq$catch[r]
    # quarterly shares applied to B7's annual catch, so the annual totals are
    # B7's exactly (SS3's echoed data file rounds catches to six digits)
    tot_q = apply(catch_q, c(1, 3), sum)
    stopifnot(max(abs(tot_q - d$catch_obs_ysf[, 1, ]) / pmax(d$catch_obs_ysf[, 1, ], 1)) < 1e-5)
    for (s in 1:4) catch_q[, s, ] = ifelse(tot_q > 0, catch_q[, s, ] / tot_q, 0) * d$catch_obs_ysf[, 1, ]
    catch_even = array(rep(d$catch_obs_ysf[, 1, ] / 4, each = 1), dim = c(n_year, n_f))
    catch_even = aperm(array(rep(catch_even, 4), dim = c(n_year, n_f, 4L)), c(1, 3, 2))
    # opal has one weight-at-age for the whole year. The annual models use
    # SS3's start-of-year weights, consistent with removing the catch at the
    # start of the year. With the catch taken through the year those weights
    # turn each tonne into far too many fish (age 0 weighs 0.03 kg at the
    # start of the year and 9 kg at mid-year), and at SS3's values the
    # quarterly population collapses; the quarterly models use SS3's
    # mid-year weight-at-age instead (the spawning output is unchanged).
    quarterly = function(x, catch) {
        dd = x$data
        dd$n_season = 4L
        dd$catch_obs_ysf = catch
        dd$weight = eg_mid$weight
        # indices in season 1 of their year (as B7: month 1); the compositions
        # stay indexed by year (opal predicts them from the year's catch)
        cp = dd$cpue_data
        cp$ts = (cp$ts - 1L) * 4L + 1L
        dd$cpue_data = cp
        opal_update(x, data = dd)
    }
    # SS3's dynamics replayed with each quarter's catch removed at the start of
    # the quarter (as quarterly opal), by catch weight-at-age
    replay_quarterly = function(catch, W) {
        sp = d$spawning_potential
        S4 = exp(-d$M / 4)
        num = d$catch_units_f != 1
        N = n_ss3(yrs[1])
        ssb = numeric(n_year)
        for (y in seq_len(n_year)) {
            N[1] = recruits_ss3[y]
            ssb[y] = sum(N * sp)
            for (s in 1:4) {
                cy = catch[y, s, ]
                hr = numeric(n_age)
                for (f in which(cy > 0)) hr = hr + cy[f] / sum(N * sel_ss3[f, ] * (if (num[f]) 1 else W[f, ])) * sel_ss3[f, ]
                N = N * (1 - pmin(hr, 0.999)) * S4
            }
            N_next = c(0, N[-n_age])
            N_next[n_age] = N_next[n_age] + N[n_age]
            N = N_next
        }
        ssb
    }
    timing_q = rbindlist(lapply(list(
        list(catch_q, W_beg, "by quarter", "start of year"),
        list(catch_q, W_mid, "by quarter", "mid-year"),
        list(catch_q, W_fleet, "by quarter", "SS3 fleet catch weight"),
        list(catch_even, W_beg, "even", "start of year"),
        list(catch_even, W_mid, "even", "mid-year")),
        function(a) data.table(year = yrs, catch = a[[3]], weight = a[[4]],
                               ssb = replay_quarterly(a[[1]], a[[2]]), ss3 = ssb_ss3)))
    timing_q[, rel := ssb / ss3 - 1]
    print(timing_q[, .(rel_end = rel[.N], rel_max = max(abs(rel))), by = .(catch, weight)])

    d0 = quarterly(o0, catch_q) %>% build_with_bounds()
    d0_even = quarterly(o0, catch_even) %>% build_with_bounds()
    # D: O2b's configuration (map, bounds) and estimates as the start
    o2b_obj = refit_fixing(fixed_lvl)$x
    fl_free = refit_fixing(fixed_lvl)$free_lvl
    sb = sel_bounds()
    # With the catch taken through the year, M acts on fish before they are
    # caught, so the same catch costs more fish. For the female variant
    # O2b's estimates describe a stock too small for that (the quarterly
    # population collapses after 2010), so D starts from O2b's estimates with
    # B0 raised in steps of 0.1 on the log scale until the objective and its
    # gradient are finite (one step of 0 for the sex-averaged variant).
    p_d = fit_o2b$fit$parameters
    d_obj = NULL
    for (k in 0:30) {
        p_try = p_d
        p_try$log_B0 = p_d$log_B0 + 0.1 * k
        # only a non-finite start is caught; any other error stops the script
        d_obj = tryCatch(quarterly(opal_update(o2b_obj, parameters = p_try), catch_q) %>%
                             build_with_bounds(sel = list(lower = sb$lower[fl_free], upper = sb$upper[fl_free])),
                         error = function(e) {
                             if (grepl("not finite", conditionMessage(e))) NULL else stop(e)
                         })
        if (!is.null(d_obj)) break
    }
    stopifnot(!is.null(d_obj))
    d_start_log_B0_step = 0.1 * k
    fit_d = get_fit("D", d_obj)
    # D refitted, as O2b, with any further selectivity parameters it cannot
    # determine fixed at its estimates: those on a bound, then those that
    # dominate a flat direction of the Hessian (the sex-averaged D has a
    # descending width just inside its lower bound with a slightly negative
    # Hessian eigenvalue)
    d_fixed_lvl = fixed_lvl
    if (!isTRUE(fit_d$validation$fit$passes) && is.null(attr(fit_d, "fixed_levels"))) {
        for (pass in 1:4) {
            free_now = refit_fixing(d_fixed_lvl)$free_lvl
            at_b = free_now[sel_on_bound(fit_d)]
            new_lvl = unique(c(at_b, sel_lvl_flat(fit_d, free_now)))
            if (!length(new_lvl)) break
            d_fixed_lvl = c(d_fixed_lvl, new_lvl)
            r = refit_fixing(d_fixed_lvl)
            x = quarterly(opal_update(r$x, parameters = fit_d$fit$parameters), catch_q) %>%
                build_with_bounds(sel = list(lower = sb$lower[r$free_lvl], upper = sb$upper[r$free_lvl]))
            message("refitting D with selectivity levels ", paste(new_lvl, collapse = ", "), " fixed")
            fit_d = fit_opal(x)
            if (isTRUE(fit_d$validation$fit$passes)) break
        }
        attr(fit_d, "fixed_levels") = d_fixed_lvl
        saved$D = fit_d
        saveRDS(saved, fits_file, compress = "xz")
    } else if (!is.null(attr(fit_d, "fixed_levels"))) {
        d_fixed_lvl = attr(fit_d, "fixed_levels")
    }
    d_sel_fixed = rbindlist(lapply(setdiff(d_fixed_lvl, fixed_lvl), function(k) {
        rc = which(lvl_matrix == k, arr.ind = TRUE)[1, ]
        data.table(fishery = rc[[1]], column = rc[[2]],
                   name = converted$meta$dims$.fleet_names[ss3_fleet_f[rc[[1]]]])
    }))
    catch_season = data.table(season = 1:4, share = apply(catch_q, 2, sum) / sum(catch_q))

#_____________________________________________________________________________
# Results
    ts_ss3 = as.data.table(rep_ss3$timeseries)[Era == "TIME" & Yr %in% yrs]
    ssb0_ss3 = ss3_derived_value(rep_ss3, c("SSB_Virgin", "SSB_Unfished"))

    models = list(O0 = o0, O1 = fit_o1, O2 = fit_o2, O2b = fit_o2b, `O1-noLF` = fit_o1_nolf,
                  `O1-noS6` = fit_o1_nos6, D0 = d0, `D0-even` = d0_even, D = fit_d)
    trajectory = rbind(
        data.table(model = ss3_label, year = yrs, ssb = ts_ss3$SpawnBio,
                   depletion = ts_ss3$SpawnBio / ssb0_ss3, recruits = recruits_ss3),
        rbindlist(lapply(names(models), function(id) {
            r = opal_report(models[[id]])
            data.table(model = id, year = yrs, ssb = r$spawning_biomass_y[seq_len(n_year)],
                       depletion = r$static_depletion_y[seq_len(n_year)],
                       recruits = r$number_ysa[seq_len(n_year), 1, 1])
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
            init_F6      = r$init_F_f[f_init[1]],
            ssb_2021     = r$spawning_biomass_y[n_year],
            depletion_2021 = r$static_depletion_y[n_year],
            lp_cpue      = sum(r$lp_cpue),
            lp_lf        = sum(r$lp_lf),
            lp_rec       = sum(r$lp_rec)
        )
    }))
    init_f_ss3 = pars_ss3[grepl("^InitF_", Label), Value]
    summary_dt = rbind(data.table(model = ss3_label, n_par = rep_ss3$N_estimated_parameters,
                                  B0 = ssb0_ss3, init_F6 = init_f_ss3,
                                  ssb_2021 = last(ts_ss3$SpawnBio),
                                  depletion_2021 = last(ts_ss3$SpawnBio) / ssb0_ss3),
                       summary_dt, fill = TRUE)
    print(summary_dt)

    # Index fits
    keep_idx = cfg$index$keep_ss3_fleets
    cpue = rbindlist(lapply(c("O1", "O2b", "O1-noS6", "D"), function(id) {
        dd = as.data.table(models[[id]]$data$cpue_data)
        dd[, .(model = id, year, index, obs = value, se, pred = opal_report(models[[id]])$cpue_pred)]
    }))
    cpue_ss3 = as.data.table(rep_ss3$cpue)[Fleet %in% keep_idx,
                                           .(model = ss3_label, year = Yr,
                                             index = match(Fleet, keep_idx),
                                             obs = Obs, se = SE, pred = Exp)]
    cpue = rbind(cpue_ss3, cpue)
    index_names = data.table(index = seq_along(keep_idx), ss3_fleet = keep_idx,
                             name = converted$meta$dims$.fleet_names[keep_idx])

    # Selectivity-at-age: SS3 against opal (O2b) and opal's start (O2 before fitting)
    sel_o2  = opal_report(fit_o2b)$sel_fya[, n_year, ]
    sel_o2_start = opal_report(o2_start)$sel_fya[, n_year, ]
    long_sel = function(m, id) data.table(model = id, fishery = rep(seq_len(nrow(m)), times = ncol(m)),
                                          age = rep(d$age_a, each = nrow(m)), sel = as.vector(m))
    selectivity = rbind(long_sel(sel_ss3, ss3_label), long_sel(sel_o2, "O2b"))
    sel_start_diff = data.table(fishery = seq_len(n_f),
                                name = converted$meta$dims$.fleet_names[ss3_fleet_f],
                                max_abs_diff = apply(abs(sel_o2_start - sel_ss3), 1, max))

    # Length-at-age: SS3 (start of year, mid-year) and opal's Schnute fit
    laa = data.table(age = d$age_a, ss3_beg = eg_beg$length, ss3_mid = eg_mid$length,
                     opal_beg_fit = growth_beg$laa_opal,
                     opal = converted$meta$growth_fit$laa_opal)

    cache = list(
        created     = format(Sys.time(), tz = "UTC", usetz = TRUE),
        session     = list(r = as.character(getRversion()),
                           opal = as.character(packageVersion("opal")),
                           r4ss = as.character(packageVersion("r4ss"))),
        fleets      = data.table(fishery = seq_len(n_f), ss3_fleet = ss3_fleet_f,
                                 name = converted$meta$dims$.fleet_names[ss3_fleet_f],
                                 units = d$catch_units_f,
                                 lf = seq_len(n_f) %in% cfg$lf$keep_fisheries),
        index_names = index_names,
        checks      = checks,
        recon       = recon,
        growth_fit  = converted$meta$growth_fit[c("k", "CV1", "CV2", "rmse", "max_abs_error",
                                                  "laa_ss3", "laa_opal", "sd_ss3", "sd_opal")],
        growth_fit_beg = growth_beg[c("rmse", "max_abs_error")],
        laa         = laa,
        bridge      = bridge,
        init_numbers = init_numbers,
        timing      = timing,
        seed_check  = seed_check,
        sel_fixed   = sel_fixed,
        sel_start_diff = sel_start_diff,
        summary     = summary_dt,
        trajectory  = trajectory,
        cpue        = cpue,
        selectivity = selectivity,
        catch_season = catch_season,
        timing_q    = timing_q,
        d_sel_fixed = d_sel_fixed,
        d_start_log_B0_step = d_start_log_B0_step,
        variant     = variant
    )
    saveRDS(cache, file.path(dir_cache, paste0("npo-swo-opal-cache", sfx, ".rds")), compress = "xz")
