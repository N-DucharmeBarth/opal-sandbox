# Nicholas Ducharme-Barth
# 2026/09/18
# Validate converted opal inputs, and reconcile them against the source SS3
# model. The reconciliation is the part worth showing in the vignette: it is
# the evidence that the conversion preserved the data, separately from any
# question about whether the two models produce the same dynamics.
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
#' Structural validation of an opal data / parameter pair
#'
#' Checks internal dimensional consistency only. It does not compare against
#' SS3; use \code{ss3_reconcile} for that.
#'
#' @param x output of \code{ss3_to_opal}.
#' @param strict \code{logical}; error on failure rather than returning the
#'   check table.
#' @return a \code{data.table} of checks, invisibly when all pass.
validate_opal_inputs = function(x, strict = TRUE) {

    d  = x$data
    p  = x$parameters
    lf = x$lf

    chk = list()
    add = function(name, pass, detail = "") {
        chk[[length(chk) + 1L]] <<- data.table::data.table(
            check = name, pass = isTRUE(pass), detail = detail)
    }

    # opal_obj() refuses any data element holding a closure, environment or
    # external pointer. A data.table's .internal.selfref is the usual culprit,
    # so check here where the failure can name the element.
    nonportable = ss3_nonportable_elements(d)
    add("data list is portable for opal_obj()", length(nonportable) == 0,
        paste(nonportable, collapse = ", "))

    add("n_age matches age_a", length(d$age_a) == d$n_age,
        sprintf("length(age_a)=%d, n_age=%d", length(d$age_a), d$n_age))
    add("years length matches n_year", length(d$years) == d$n_year)
    add("M length is n_age or n_len",
        length(d$M) %in% c(d$n_age, d$n_len),
        sprintf("length(M)=%d", length(d$M)))
    add("maturity length is n_age or n_len",
        length(d$maturity) %in% c(d$n_age, d$n_len),
        sprintf("length(maturity)=%d", length(d$maturity)))
    if (is.null(d$spawning_potential) || !is.null(d$fecundity)) {
        add("fecundity length is n_age or n_len",
            length(d$fecundity) %in% c(d$n_age, d$n_len),
            sprintf("length(fecundity)=%d", length(d$fecundity)))
    }
    add("maturity is a probability",
        all(d$maturity >= 0 & d$maturity <= 1, na.rm = TRUE))
    add("M is positive and finite", all(is.finite(d$M) & d$M > 0))

    # Age-basis biology, supplied only by the endgrowth route.
    if (!is.null(d$weight)) {
        add("weight length is n_age or n_len",
            length(d$weight) %in% c(d$n_age, d$n_len),
            sprintf("length(weight)=%d", length(d$weight)))
        add("weight is positive, finite and increasing with age",
            all(is.finite(d$weight) & d$weight > 0) &&
                (length(d$weight) != d$n_age || all(diff(d$weight) >= -1e-8)))
    }
    if (!is.null(d$spawning_potential)) {
        add("spawning_potential length is n_age or n_len",
            length(d$spawning_potential) %in% c(d$n_age, d$n_len),
            sprintf("length(spawning_potential)=%d", length(d$spawning_potential)))
        add("spawning_potential is non-negative and finite",
            all(is.finite(d$spawning_potential) & d$spawning_potential >= 0))
        add("spawning_potential is non-zero somewhere",
            any(d$spawning_potential > 0))
    }
    if (identical(x$meta$bio_source, "endgrowth")) {
        add("endgrowth route supplies M on an age basis",
            length(d$M) == d$n_age)
        add("endgrowth route supplies maturity on an age basis",
            length(d$maturity) == d$n_age)
    }

    add("catch array dimensions",
        identical(dim(d$catch_obs_ysf), c(d$n_year, d$n_season, d$n_fishery)),
        paste(dim(d$catch_obs_ysf), collapse = " x "))
    add("catch is non-negative and finite",
        all(is.finite(d$catch_obs_ysf) & d$catch_obs_ysf >= 0))

    add("catch_units_f length", length(d$catch_units_f) == d$n_fishery)
    add("catch_units_f codes are 1 or 2", all(d$catch_units_f %in% c(1, 2)))
    add("sel_type_f length", length(d$sel_type_f) == d$n_fishery)
    add("sel_type_f uses opal codes (1-3), not SS3 patterns",
        all(d$sel_type_f %in% 1:3),
        paste(sort(unique(d$sel_type_f)), collapse = ", "))

    cp = d$cpue_data
    add("cpue_data has required columns",
        all(c("year", "month", "ts", "fishery", "index", "metric",
              "units", "value", "se") %in% names(cp)))
    add("cpue index runs 1:n_index",
        identical(sort(unique(as.integer(cp$index))), seq_len(d$n_index)),
        sprintf("indices present: %s", paste(sort(unique(cp$index)), collapse = ", ")))
    add("cpue ts within model window",
        all(cp$ts >= 1 & cp$ts <= d$n_year))
    add("cpue values positive and finite", all(is.finite(cp$value) & cp$value > 0))
    add("cpue se positive and finite", all(is.finite(cp$se) & cp$se > 0))

    if (!is.null(lf)) {
        expected_bins = seq(d$len_bin_start, by = d$len_bin_width, length.out = d$n_len)
        add("lf bins match the model length structure",
            all(sort(unique(lf$bin)) %in% expected_bins))
        add("lf ts within model window", all(lf$ts >= 1 & lf$ts <= d$n_year))
        add("lf fisheries within 1:n_fishery",
            all(lf$fishery >= 1 & lf$fishery <= d$n_fishery))
        add("lf values non-negative and finite",
            all(is.finite(lf$value) & lf$value >= 0))
    }

    add("rdev_y length matches n_year", length(p$rdev_y) == d$n_year)
    add("init_rdev_a length matches n_age", length(p$init_rdev_a) == d$n_age)
    add("log_cpue_q length matches n_index", length(p$log_cpue_q) == d$n_index)
    add("par_sel rows match n_fishery", nrow(p$par_sel) == d$n_fishery)
    add("growth parameters are finite",
        all(is.finite(c(p$log_L1, p$log_L2, p$log_k, p$log_CV1, p$log_CV2))))
    add("L2 exceeds L1", exp(p$log_L2) > exp(p$log_L1))

    sel = x$selectivity$sel_fa_external
    if (!is.null(sel)) {
        add("sel_fa_external dimensions",
            identical(dim(sel), c(as.integer(d$n_fishery), as.integer(d$n_age))),
            paste(dim(sel), collapse = " x "))
        add("sel_fa_external within [0, 1]",
            all(sel >= 0 & sel <= 1 + 1e-8, na.rm = TRUE))
    }

    out = data.table::rbindlist(chk)
    if (strict && any(!out$pass)) {
        stop(paste0("Validation failed:\n",
                    paste(sprintf("  - %s (%s)", out$check[!out$pass],
                                  out$detail[!out$pass]), collapse = "\n")),
             call. = FALSE)
    }
    out
}

#_____________________________________________________________________________
#' Reconcile converted inputs against the source SS3 model
#'
#' Compares totals that should survive conversion exactly, and reports the
#' quantities that can only be approximated. Anything with a relative
#' difference above its tolerance is flagged.
#'
#' Unfished spawning potential per recruit is held to \code{spr0_tol} rather
#' than \code{tol}, for two reasons. Report.sso prints six significant figures,
#' which alone limits the comparison to about 5e-6. And the two models close
#' the plus group differently: SS3 3.24 iterates the equilibrium out to
#' 3 x Nages and closes the tail with the annual rather than the seasonal Z,
#' while opal uses the closed-form geometric tail at the per-timestep M. For
#' this model that is a 4.7e-5 relative difference. The errors this check
#' exists to catch are far larger: a misplaced female fraction is a factor of
#' two, M in the wrong units roughly a factor of four.
#'
#' @param x output of \code{ss3_to_opal}.
#' @param tol relative tolerance for the exact comparisons.
#' @param spr0_tol relative tolerance for the SPR0 comparison.
#' @return a \code{data.table} with one row per comparison.
ss3_reconcile = function(x, tol = 1e-6, spr0_tol = 1e-4) {

    d    = x$data
    ss3  = x$ss3
    dims = x$meta$dims
    fish = x$meta$fisheries

    rows = list()
    add = function(quantity, ss3_value, opal_value, note = "", row_tol = tol) {
        rel = if (is.finite(ss3_value) && abs(ss3_value) > 0)
            abs(opal_value - ss3_value) / abs(ss3_value) else NA_real_
        rows[[length(rows) + 1L]] <<- data.table::data.table(
            quantity   = quantity,
            ss3        = ss3_value,
            opal       = opal_value,
            rel_diff   = rel,
            tol        = row_tol,
            pass       = if (is.na(row_tol)) NA else is.na(rel) | rel <= row_tol,
            note       = note
        )
    }

    # ---- Catch: total by extraction fleet over the model window ------------
    catch_ss3 = data.table::as.data.table(ss3$dat$catch)
    yr_col = intersect(c("year", "Yr", "yr"), names(catch_ss3))[1]
    se_col = intersect(c("seas", "Seas", "season"), names(catch_ss3))[1]
    in_win = catch_ss3[[yr_col]] >= dims$.first_ts & catch_ss3[[yr_col]] <= dims$.last_ts
    if (all(c("fleet", "catch") %in% names(catch_ss3))) {
        # SS 3.30 long format
        for (f in seq_len(fish$n_extract)) {
            sf = fish$extract_fleets[f]
            add(sprintf("catch total, fleet %d (%s)", sf, dims$.fleet_names[sf]),
                sum(as.numeric(catch_ss3$catch[in_win & catch_ss3$fleet == sf])),
                sum(d$catch_obs_ysf[, , f]))
        }
    } else {
        fleet_cols = setdiff(names(catch_ss3), c(yr_col, se_col))
        for (f in seq_along(fleet_cols)) {
            add(sprintf("catch total, fleet %d (%s)", f, fleet_cols[f]),
                sum(as.numeric(catch_ss3[[fleet_cols[f]]][in_win])),
                sum(d$catch_obs_ysf[, , f]))
        }
    }

    # ---- CPUE: observation counts, and values before normalisation ---------
    cpue_ss3 = data.table::as.data.table(ss3$dat$CPUE)
    data.table::setnames(cpue_ss3, tolower(names(cpue_ss3)))
    fl_col = intersect(c("index", "fleet"), names(cpue_ss3))[1]
    yr_col = intersect(c("year", "yr"), names(cpue_ss3))[1]
    n_ss3 = sum(cpue_ss3[[fl_col]] %in% fish$index_map$ss3_fleet &
                    cpue_ss3[[yr_col]] >= dims$.first_ts &
                    cpue_ss3[[yr_col]] <= dims$.last_ts)
    add("CPUE observation count (kept SS3 index fleets)", n_ss3, nrow(d$cpue_data),
        "index splitting preserves observation count")

    # ---- Length compositions: observation-years and total sample size ------
    if (!is.null(x$lf) && !is.null(ss3$dat$lencomp)) {
        lc = data.table::as.data.table(ss3$dat$lencomp)
        yr_col = intersect(c("Yr", "year"), names(lc))[1]
        fl_col = intersect(c("FltSvy", "fleet"), names(lc))[1]
        n_col  = intersect(c("Nsamp", "nsamp"), names(lc))[1]
        keep = lc[[fl_col]] > 0 & lc[[yr_col]] >= dims$.first_ts &
            lc[[yr_col]] <= dims$.last_ts
        add("LF total input sample size",
            sum(as.numeric(lc[[n_col]][keep])),
            sum(x$lf$value),
            "opal pseudo-counts should sum to SS3 Nsamp")
        add("LF observation rows",
            sum(keep),
            x$lf[, data.table::uniqueN(paste(ts, fishery))])
    }

    # ---- Unfished spawning potential per recruit ---------------------------
    # This single number tests three things at once: the units of M, the
    # female-fraction convention in spawning potential, and the plus-group
    # treatment. It is the most informative check in the table.
    spr0_opal = opal_spr0(d)
    spr0_ss3  = ss3_spr0(ss3$rep, dims)
    if (is.finite(spr0_opal) && is.finite(spr0_ss3)) {
        add("unfished spawning potential per recruit", spr0_ss3, spr0_opal,
            paste0("tests M units, female fraction and plus group together; ",
                   "per recruit at the first opal age; tolerance allows for ",
                   "SS3's truncated plus-group equilibrium"),
            row_tol = spr0_tol)
    }

    # ---- Weight-at-age -----------------------------------------------------
    if (!is.null(d$weight) && length(d$weight) == d$n_age && !is.null(ss3$rep)) {
        wb = x$meta$config$bio$endgrowth$weight_basis %||%
            x$meta$config$bio$endgrowth$length_basis %||% "Beg"
        eg = tryCatch(ss3_endgrowth_table(ss3$rep, dims, basis = wb), error = function(e) NULL)
        if (!is.null(eg)) {
            add(sprintf("total weight summed over ages (Wt_%s)", wb),
                sum(eg$weight), sum(d$weight),
                "age-basis weight should be identical to SS3")
        }
    }

    # ---- Biology: quantities that can only be approximated -----------------
    gf = x$meta$growth_fit
    if (identical(gf$source, "ss3_age_k") && is.finite(gf$max_rel_error %||% NA_real_)) {
        # SS3's own growth form in opal: an exact comparison, limited only by
        # Report.sso printing Len_Beg to six significant figures.
        laa_tol = 5e-5
        rows[[length(rows) + 1L]] = data.table::data.table(
            quantity = "mean length-at-age, SS3 age-specific K (max relative error)",
            ss3      = NA_real_,
            opal     = NA_real_,
            rel_diff = gf$max_rel_error,
            tol      = laa_tol,
            pass     = gf$max_rel_error <= laa_tol,
            note     = sprintf("opal get_growth() vs SS3 Len_Beg; max |error| %.1e cm",
                               gf$max_abs_error)
        )
    } else if (identical(gf$source, "ss3_cessation") &&
               is.finite(gf$max_rel_error_below_plus %||% NA_real_)) {
        # SS3's growth cessation curve and SD-linear (or CV-linear) growth
        # variability in opal: exact below the plus group, to Report.sso print
        # precision. SS3 averages the plus group's length over the ages it
        # accumulates and holds its SD at SD_old; opal does neither.
        laa_tol = 5e-5
        n = length(gf$laa_opal)
        rows[[length(rows) + 1L]] = data.table::data.table(
            quantity = "mean length-at-age, SS3 growth cessation (max relative error below the plus group)",
            ss3 = NA_real_, opal = NA_real_,
            rel_diff = gf$max_rel_error_below_plus, tol = laa_tol,
            pass = gf$max_rel_error_below_plus <= laa_tol,
            note = "opal get_growth(growth_model = \"cessation\") vs SS3 Len_Beg")
        rows[[length(rows) + 1L]] = data.table::data.table(
            quantity = sprintf("SD of length-at-age, SS3 CV_Growth_Pattern %d (max relative error below the plus group)",
                               gf$cv_pattern),
            ss3 = NA_real_, opal = NA_real_,
            rel_diff = gf$sd_max_rel_error_below_plus, tol = laa_tol,
            pass = gf$sd_max_rel_error_below_plus <= laa_tol,
            note = sprintf("opal get_sd_at_age(sd_model = \"%s\") vs SS3 SD_Beg", gf$sd_model))
        add("mean length at the plus-group age", gf$laa_ss3[n], gf$laa_opal[n],
            "SS3 averages the plus group's length over the ages it accumulates; opal does not",
            row_tol = NA_real_)
        add("SD of length at the plus-group age", gf$sd_ss3[n], gf$sd_opal[n],
            "SS3 holds the oldest age at SD_old; opal's upper reference is SS3's line at that age's length",
            row_tol = NA_real_)
    } else if (!is.null(gf) && is.finite(gf$rmse %||% NA_real_)) {
        # Whether this residual matters depends on what growth is still used
        # for. On the endgrowth route it affects only the PLA for selectivity;
        # on the control route it propagates into spawning biomass.
        consequence = if (identical(x$meta$bio_source, "endgrowth"))
            "affects selectivity PLA only" else
            "propagates into weight-, maturity- and spawning-biomass-at-age"
        rows[[length(rows) + 1L]] = data.table::data.table(
            quantity = "mean length-at-age, Schnute approximation",
            ss3      = NA_real_,
            opal     = NA_real_,
            rel_diff = gf$rmse / mean(gf$laa_ss3),
            pass     = NA,
            note     = sprintf("RMSE %.3f cm, max |error| %.3f cm; %s",
                               gf$rmse, gf$max_abs_error, consequence)
        )
    }
    if (!is.null(gf) && !identical(gf$source, "ss3_cessation") &&
        is.finite(gf$sd_rmse %||% NA_real_)) {
        rows[[length(rows) + 1L]] = data.table::data.table(
            quantity = "SD of length-at-age, CV approximation",
            ss3 = NA_real_, opal = NA_real_,
            rel_diff = gf$sd_rmse / mean(gf$sd_ss3),
            pass = NA,
            note = sprintf("RMSE %.3f cm; controls how far age selectivity smears across length bins",
                           gf$sd_rmse)
        )
    }

    out = data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
    out[]
}

#_____________________________________________________________________________
#' Print a compact conversion report
#'
#' @param x output of \code{ss3_to_opal}.
#' @param checks output of \code{validate_opal_inputs}.
#' @param recon output of \code{ss3_reconcile}.
#' @param assessment optional \code{opal_obj} built from the conversion.
ss3_conversion_report = function(x, checks = NULL, recon = NULL,
                                 assessment = NULL) {
    d = x$data
    m = x$meta

    cat("\n=== SS3 -> opal conversion ===\n")
    cat(sprintf("  source      : %s\n", m$provenance$dir_ss3))
    cat(sprintf("  SS3 version : %s (r4ss %s)\n",
                m$provenance$ss3_version, m$provenance$r4ss_version))
    cat(sprintf("  config      : %s\n", m$config$config_path %||% "<none>"))
    cat(sprintf("  timesteps   : %d (SS3 %d-%d, %.2f per calendar year)\n",
                d$n_year, m$dims$.first_ts, m$dims$.last_ts, m$dims$.ts_per_year))
    cat(sprintf("  ages        : opal %d-%d = SS3 %d-%d\n",
                min(d$age_a), max(d$age_a),
                min(m$dims$.ss3_ages), max(m$dims$.ss3_ages)))
    cat(sprintf("  fisheries   : %d (%d extraction + %d index)\n",
                d$n_fishery, m$fisheries$n_extract, m$fisheries$n_index))
    cat(sprintf("  length bins : %d, %g cm wide from %g cm\n",
                d$n_len, d$len_bin_width, d$len_bin_start))
    cat(sprintf("  biology     : %s\n",
                if (identical(m$bio_source, "endgrowth"))
                    "endgrowth (age), PLA not used for biology" else
                    "control file (length), converted through the PLA"))
    if (!is.null(d$spawning_potential)) {
        cat("                spawning_potential supplied directly from SS3 Mat*Fecund\n")
    }
    cat(sprintf("  growth      : L1=%.2f L2=%.2f k=%.4f CV1=%.3f CV2=%.3f (%s)\n",
                exp(x$parameters$log_L1), exp(x$parameters$log_L2),
                exp(x$parameters$log_k), exp(x$parameters$log_CV1),
                exp(x$parameters$log_CV2), m$growth_fit$source))
    if (is.finite(m$growth_fit$rmse %||% NA_real_)) {
        cat(sprintf("                %s vs SS3 length-at-age: RMSE %.2g cm, max %.2g cm\n",
                    if (identical(m$growth_fit$source, "ss3_age_k"))
                        "SS3 age-specific K growth" else "Schnute fit",
                    m$growth_fit$rmse, m$growth_fit$max_abs_error))
    }
    if (length(m$lf_zero_lambda_fleets)) {
        cat(sprintf("  NOTE        : SS3 gives zero LF lambda to fleet(s) %s\n",
                    paste(m$lf_zero_lambda_fleets, collapse = ", ")))
    }

    if (!is.null(checks)) {
        n_fail = sum(!checks$pass)
        cat(sprintf("\n  validation  : %d/%d checks passed\n",
                    sum(checks$pass), nrow(checks)))
        if (n_fail) print(checks[pass == FALSE])
    }
    if (!is.null(recon)) {
        bad = recon[pass %in% FALSE]
        cat(sprintf("  reconcile   : %d/%d comparisons within tolerance\n",
                    sum(recon$pass, na.rm = TRUE), sum(!is.na(recon$pass))))
        if (nrow(bad)) print(bad)
        notes = recon[!is.na(note) & nzchar(note) & is.na(pass)]
        if (nrow(notes)) print(notes[, .(quantity, note)])
    }
    if (!is.null(assessment)) {
        cat("\n")
        print(summary(assessment))
    }
    cat("\n")
    invisible(x)
}

#_____________________________________________________________________________
#' Replay the converted dynamics outside opal, under either model's timing
#'
#' The bridging check compares opal's spawning biomass with SS3's. Two
#' structural differences remain even when every input is converted exactly:
#'
#'   catch   opal removes each timestep's catch instantaneously at the start of
#'           the step, before natural mortality; SS3 (F_Method 2 or 3) applies
#'           F and M together through the Baranov equation.
#'   SR      opal's first age class at t + 1 is recruited from spawning biomass
#'           at t + 1. With SS3 age 0 dropped, the matching SS3 fish were
#'           recruited from spawning biomass at t.
#'
#' This replays the population on the converted inputs (age-basis biology and
#' sel_fa_external) with each difference switched on or off. With opal timing
#' it reproduces opal_report() to numerical precision; with SS3 timing it
#' should track SS3's SSB closely. A large residual under SS3 timing points at
#' the conversion; a gap that closes when the timing is switched does not.
#'
#' @param x output of \code{ss3_to_opal} (endgrowth route).
#' @param catch \code{"start"} (opal) or \code{"baranov"} (SS3).
#' @param sr \code{"opal"} or \code{"ss3"} spawner-recruit timing.
#' @return \code{numeric} spawning biomass at the start of each timestep.
ss3_bridge_replica = function(x, catch = c("start", "baranov"), sr = c("opal", "ss3")) {
    catch = match.arg(catch)
    sr    = match.arg(sr)
    d   = x$data
    p   = x$parameters
    sel = x$selectivity$sel_fa_external
    sp  = d$spawning_potential
    if (is.null(sel) || is.null(sp) || length(d$M) != d$n_age) {
        stop("ss3_bridge_replica() needs age-basis biology and sel_fa_external.", call. = FALSE)
    }

    n_age = d$n_age
    n_yr  = d$n_year
    M     = d$M
    S     = exp(-M)
    rel   = cumprod(c(1, S[-n_age]))
    rel[n_age] = rel[n_age] / (1 - S[n_age])
    B0    = exp(p$log_B0)
    h     = exp(p$log_h)
    sig   = exp(p$log_sigma_r)
    R0    = B0 / sum(sp * rel)
    alpha = 4 * h * R0 / (5 * h - 1)
    beta  = B0 * (1 - h) / (5 * h - 1)
    badj  = d$bias_adj_y %||% rep(1, n_yr)

    fleets = which(apply(d$catch_obs_ysf[, 1, , drop = FALSE], 3, sum) > 0)
    num    = d$catch_units_f != 1
    wt     = d$weight %||% stop("ss3_bridge_replica() needs age-basis weight.", call. = FALSE)

    # Initial age structure: equilibrium under the initial F (as opal's
    # get_initial_numbers(), one season), then the initial deviations
    N = R0 * rel
    init_F = exp(p$log_init_F_f %||% rep(-Inf, nrow(sel)))
    if (any(init_F > 1e-6)) {
        S_init = S * (1 - colSums((1 - exp(-init_F)) * sel))
        rel_f = cumprod(c(1, S_init[-n_age]))
        rel_f[n_age] = rel_f[n_age] / (1 - S_init[n_age])
        R_eq = alpha - beta / sum(sp * rel_f)
        N = R_eq * rel_f
    }
    init_dev = p$init_rdev_a %||% rep(0, n_age)
    init_badj = d$init_bias_adj_a %||% rep(0, n_age)
    N = N * exp(init_dev - init_badj * 0.5 * sig^2)
    # SS3 recruits age 0 from the same year's spawning biomass, which is opal's
    # timing when age 0 is the first opal age
    age0 = isTRUE(x$meta$dims$.ss3_ages[1] == 0L)

    ssb = numeric(n_yr)
    for (y in seq_len(n_yr)) {
        ssb[y] = sum(N * sp)
        cy = d$catch_obs_ysf[y, 1, ]
        active = fleets[cy[fleets] > 0]
        if (catch == "start") {
            hr = numeric(n_age)
            for (f in active) {
                avail = N * sel[f, ] * (if (num[f]) 1 else wt)
                hr = hr + cy[f] / sum(avail) * sel[f, ]
            }
            N_end = N * (1 - hr) * S
        } else {
            # Solve each fleet's F so the Baranov catch matches the observation.
            F_f = numeric(nrow(sel))
            for (f in active) F_f[f] = cy[f] / sum(N * sel[f, ] * (if (num[f]) 1 else wt))
            for (it in 1:50) {
                Z = M + colSums(F_f * sel)
                for (f in active) {
                    pred = sum(F_f[f] * sel[f, ] / Z * N * (1 - exp(-Z)) * (if (num[f]) 1 else wt))
                    F_f[f] = F_f[f] * cy[f] / pred
                }
            }
            N_end = N * exp(-(M + colSums(F_f * sel)))
        }
        N_next = c(0, N_end[-n_age])
        N_next[n_age] = N_next[n_age] + N_end[n_age]
        sb_R = if (sr == "ss3" && !age0) ssb[y] else sum(N_next * sp)
        N_next[1] = alpha * sb_R / (beta + sb_R) *
            exp(p$rdev_y[y] - badj[y] * 0.5 * sig^2)
        N = N_next
    }
    ssb
}

#' Bridging table: SS3 spawning biomass against opal and the replayed dynamics
#'
#' @param x output of \code{ss3_to_opal}.
#' @param report output of \code{opal_report()} on the SS3-selectivity object.
#' @return a \code{data.table} with one row per timestep.
ss3_bridge_table = function(x, report) {
    ts_ss3 = x$ss3$rep$timeseries
    ts_ss3 = ts_ss3[ts_ss3$Era == "TIME" &
                        ts_ss3$Yr >= x$meta$dims$.first_ts &
                        ts_ss3$Yr <= x$meta$dims$.last_ts, ]
    n = x$data$n_year
    # The replay needs age-basis biology, so it is NA on the control route,
    # where spawning potential only exists inside opal (through the PLA).
    replay = function(...) tryCatch(ss3_bridge_replica(x, ...),
                                    error = function(e) rep(NA_real_, n))
    out = data.table::data.table(
        ts              = seq_len(n),
        ss3             = as.numeric(ts_ss3$SpawnBio),
        opal            = as.numeric(report$spawning_biomass_y[seq_len(n)]),
        replica_opal    = replay(catch = "start",   sr = "opal"),
        replica_ss3     = replay(catch = "baranov", sr = "ss3")
    )
    out[, rel_opal := opal / ss3 - 1]
    out[, rel_replica_ss3 := replica_ss3 / ss3 - 1]
    out[]
}
