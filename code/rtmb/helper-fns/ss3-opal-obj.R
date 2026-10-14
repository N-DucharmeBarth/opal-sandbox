# Nicholas Ducharme-Barth
# 2026/10/13
# Bridge the SS3 converter output into opal's staged S3 assessment object.
#
# opal_obj() is the single carrier for data, parameters, map, bounds, controls,
# fitted results and provenance. Two consequences for a converter:
#
#   1. The object refuses non-portable state. A data.table carries an
#      externalptr in its .internal.selfref attribute, so anything that reaches
#      the data list must be coerced to a plain data.frame first.
#   2. Scientific fields must not be assigned into directly. validate_opal_obj()
#      compares a stored identity hash and errors with "Model state was modified
#      directly; use opal_update()". Switching the selectivity route therefore
#      goes through opal_update(), not `x$data$sel_fa_external <- ...`.
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
#' Strip non-portable state from a nested list
#'
#' Recursively coerces data.tables and tibbles to plain data.frames and drops
#' attributes that opal classifies as non-portable. Matrices, arrays and atomic
#' vectors pass through untouched.
#'
#' @param x any R object, typically the opal data list.
#' @return \code{x} with portable replacements substituted.
ss3_make_portable = function(x) {
    if (inherits(x, "data.table") || inherits(x, "tbl_df") || inherits(x, "spec_tbl_df")) {
        x = as.data.frame(x, stringsAsFactors = FALSE)
        attr(x, ".internal.selfref") = NULL
        attr(x, "spec")              = NULL
        attr(x, "problems")          = NULL
        rownames(x) = NULL
        return(x)
    }
    if (is.list(x) && !is.data.frame(x)) {
        for (i in seq_along(x)) x[i] = list(ss3_make_portable(x[[i]]))
        return(x)
    }
    x
}

#_____________________________________________________________________________
#' Report any element of a list that opal would reject as non-portable
#'
#' Mirrors opal's internal test so a conversion failure names the offending
#' element instead of failing inside the constructor.
#'
#' @param x a named \code{list}.
#' @return \code{character} vector of element names that are not portable.
ss3_nonportable_elements = function(x) {
    bad = function(v) {
        if (typeof(v) %in% c("closure", "environment", "externalptr", "weakref")) return(TRUE)
        if (is.list(v) && any(vapply(v, bad, logical(1L)))) return(TRUE)
        a = attributes(v)
        if (!is.null(a)) return(any(vapply(a, bad, logical(1L))))
        FALSE
    }
    names(x)[vapply(x, bad, logical(1L))]
}

#_____________________________________________________________________________
#' Build a configured opal_obj from converted SS3 inputs
#'
#' Produces an unfitted, built object: parameters, map and bounds are resolved
#' and the initial objective and gradient are checked for finiteness, but no
#' optimisation happens. Fitting belongs in the vignette, not in data prep.
#'
#' The SS3 provenance, converter config checksum and growth-approximation
#' diagnostic are carried in \code{$provenance$metadata}, so a saved object can
#' be traced back to the exact SS3 files it came from.
#'
#' @param converted output of \code{ss3_to_opal}.
#' @param cfg converter configuration list.
#' @param map optional RTMB map. \code{NULL} lets \code{opal_build} resolve the
#'   opal default via \code{get_map()}.
#' @param bounds optional bounds. \code{NULL} defers to \code{get_bounds()}.
#' @param control optional optimiser controls stored on the object.
#' @param build \code{logical}; run \code{opal_build} so the configuration is
#'   checked now rather than at first fit.
#' @return an \code{opal_obj}.
ss3_opal_obj = function(converted, cfg, map = NULL, bounds = NULL,
                        control = NULL, build = TRUE) {

    if (!requireNamespace("opal", quietly = TRUE)) {
        stop("Package 'opal' is required.", call. = FALSE)
    }

    data = ss3_make_portable(converted$data)

    offenders = ss3_nonportable_elements(data)
    if (length(offenders)) {
        stop(sprintf(
            "These data elements are not portable and opal_obj() will reject them: %s",
            paste(offenders, collapse = ", ")), call. = FALSE)
    }

    parameters = converted$parameters
    data$priors = opal::get_priors(parameters = parameters, data = data)

    m = converted$meta
    metadata = c(
        list(
            stock              = cfg$model,
            source             = "SS3",
            ss3_dir            = m$provenance$dir_ss3,
            ss3_version        = m$provenance$ss3_version,
            ss3_file_md5       = m$provenance$md5,
            r4ss_version       = m$provenance$r4ss_version,
            converter_config   = m$config$config_path,
            converter_config_md5 = m$config$config_md5,
            converted_at       = m$converted_at,
            biology_route      = m$bio_source,
            growth_source      = m$growth_fit$source,
            growth_rmse_cm     = m$growth_fit$rmse,
            growth_sd_rmse_cm  = m$growth_fit$sd_rmse,
            selectivity_route  = "opal_pla"
        ),
        as.list(cfg$metadata %||% list())
    )

    x = opal::opal_obj(
        data       = data,
        parameters = parameters,
        map        = map,
        bounds     = bounds,
        control    = control,
        metadata   = metadata
    )

    if (isTRUE(build)) x = opal::opal_build(x)
    x
}

#_____________________________________________________________________________
#' Switch an assessment onto SS3's selectivity-at-age
#'
#' Injects \code{sel_fa_external}, bypassing opal's PLA selectivity entirely.
#' This is the bridging configuration: with every other parameter held at its
#' SS3 value, the resulting dynamics should reproduce SS3 directly, which
#' separates "the conversion is wrong" from "the two models differ".
#'
#' Selectivity parameters are fixed in the map, since \code{par_sel} has no
#' effect while an external matrix is supplied.
#'
#' @param x an \code{opal_obj}.
#' @param sel_fa \code{[n_fishery, n_age]} matrix, from
#'   \code{converted$selectivity$sel_fa_external}.
#' @return an updated \code{opal_obj}.
ss3_use_external_selectivity = function(x, sel_fa) {
    stopifnot("sel_fa must be a matrix" = is.matrix(sel_fa))
    stopifnot("sel_fa rows must match n_fishery" = nrow(sel_fa) == x$data$n_fishery)
    stopifnot("sel_fa columns must match n_age"  = ncol(sel_fa) == x$data$n_age)

    data = x$data
    data$sel_fa_external = sel_fa

    map = x$map
    if (is.null(map)) map = opal::get_map(x$parameters)
    map$par_sel = factor(rep(NA, length(x$parameters$par_sel)))

    opal::opal_update(x,
                      data     = data,
                      map      = map,
                      metadata = list(selectivity_route = "ss3_external"))
}

#_____________________________________________________________________________
#' Report what the growth parameters are still doing
#'
#' Growth enters opal in two places: the PLA that converts length-basis biology
#' to age, and the PLA that maps length-based selectivity onto ages. The
#' endgrowth biology route removes the first; \code{sel_fa_external} removes
#' the second. When both are gone, the growth parameters affect nothing in the
#' likelihood and should be fixed - leaving them free would make the Hessian
#' singular and fail the estimability check.
#'
#' @param x an \code{opal_obj}.
#' @return a one-row \code{data.frame} describing the two PLA uses.
ss3_growth_role = function(x) {
    biology_uses_pla = is.null(x$data$weight) ||
        length(x$data$maturity) != x$data$n_age
    sel_uses_pla = is.null(x$data$sel_fa_external)
    data.frame(
        biology_uses_pla   = biology_uses_pla,
        selectivity_uses_pla = sel_uses_pla,
        growth_is_inert    = !biology_uses_pla && !sel_uses_pla,
        stringsAsFactors   = FALSE
    )
}

#' Fix the growth parameters in the map
#'
#' The opal default map already fixes `log_L1`, `log_L2`, `log_k`, `log_CV1`
#' and `log_CV2`. Use this only when a custom map has freed them and growth has
#' since become inert.
#'
#' @param x an \code{opal_obj}.
#' @return an updated \code{opal_obj}.
ss3_fix_growth = function(x) {
    map = x$map
    if (is.null(map)) map = opal::get_map(x$parameters)
    for (p in c("log_L1", "log_L2", "log_k", "log_CV1", "log_CV2")) {
        if (!is.null(x$parameters[[p]])) {
            map[[p]] = factor(rep(NA, length(x$parameters[[p]])))
        }
    }
    opal::opal_update(x, map = map)
}

#_____________________________________________________________________________
#' Switch an assessment back onto opal's own length-based selectivity
#'
#' Removes \code{sel_fa_external} and frees the selectivity parameters named in
#' \code{free}. Starting values should first be seeded from
#' \code{converted$selectivity$sel_target_fl}, which is SS3's selectivity
#' projected onto the length bins.
#'
#' @param x an \code{opal_obj}.
#' @param par_sel optional replacement \code{par_sel} starting values.
#' @param free \code{integer} vector of fisheries whose selectivity is estimated.
#' @param n_par_f optional \code{integer} vector giving how many parameters to
#'   free per fishery; defaults to all columns of \code{par_sel}.
#' @return an updated \code{opal_obj}.
ss3_use_opal_selectivity = function(x, par_sel = NULL, free = integer(),
                                    n_par_f = NULL) {

    data = x$data
    data$sel_fa_external = NULL

    parameters = x$parameters
    if (!is.null(par_sel)) parameters$par_sel = as.matrix(par_sel)

    map_sel = matrix(NA_integer_, nrow(parameters$par_sel), ncol(parameters$par_sel))
    if (length(free)) {
        if (is.null(n_par_f)) n_par_f = rep(ncol(map_sel), nrow(map_sel))
        counter = 0L
        for (f in free) {
            cols = seq_len(n_par_f[f])
            map_sel[f, cols] = counter + seq_along(cols)
            counter = counter + length(cols)
        }
    }

    map = x$map
    if (is.null(map)) map = opal::get_map(parameters)
    map$par_sel = factor(as.vector(map_sel))

    opal::opal_update(x,
                      data       = data,
                      parameters = parameters,
                      map        = map,
                      metadata   = list(selectivity_route = "opal_pla"))
}

#_____________________________________________________________________________
#' Seed par_sel by fitting an opal selectivity curve to SS3's implied
#' selectivity-at-length
#'
#' SS3 parameterises selectivity on age here; opal parameterises on length.
#' There is no algebraic conversion, so the starting values are obtained by
#' least squares against \code{sel_target_fl}.
#'
#' The curve functions are taken from opal itself rather than reimplemented, so
#' the seeded parameters are on exactly the scale opal will interpret. If the
#' call signature differs from the assumption below, this errors loudly instead
#' of silently seeding nonsense.
#'
#' @param sel_target_fl \code{[n_fishery, n_len]} target matrix from
#'   \code{ss3_sel_to_opal}.
#' @param len_mid \code{numeric} length bin midpoints.
#' @param sel_type_f \code{integer} opal selectivity codes per fishery.
#' @param n_par \code{integer} number of par_sel columns.
#' @param start optional \code{[n_fishery, n_par]} matrix of starting values.
#' @param lower,upper box constraints on every parameter. The defaults are
#'   opal's own \code{get_bounds()} limits for \code{par_sel}, so seeded values
#'   are valid starting values for an opal fit.
#' @return a \code{list} with \code{par_sel}, \code{n_par_f} and a per-fishery
#'   \code{rmse} vector.
ss3_seed_par_sel = function(sel_target_fl, len_mid, sel_type_f, n_par = 6L,
                            start = NULL, lower = -7, upper = 7) {

    sel_fn = list(
        `1` = list(fn = opal::sel_logistic,      n = 2L),
        `2` = list(fn = opal::sel_double_normal, n = 6L),
        `3` = list(fn = opal::sel_double_richards, n = 6L)
    )

    n_fishery = nrow(sel_target_fl)
    par_sel   = if (is.null(start)) matrix(0, n_fishery, n_par) else as.matrix(start)
    n_par_f   = integer(n_fishery)
    rmse      = rep(NA_real_, n_fishery)

    for (f in seq_len(n_fishery)) {
        spec = sel_fn[[as.character(sel_type_f[f])]]
        if (is.null(spec)) {
            stop(sprintf("Unknown opal selectivity code %s for fishery %d.",
                         sel_type_f[f], f), call. = FALSE)
        }
        n_par_f[f] = spec$n
        target = sel_target_fl[f, ]
        if (all(!is.finite(target)) || max(target, na.rm = TRUE) <= 0) next

        objective = function(p) {
            pred = tryCatch(as.numeric(spec$fn(len_mid, p)),
                            error = function(e) stop(sprintf(
                                "Could not evaluate the opal selectivity function as sel_fn(len_mid, pars): %s",
                                conditionMessage(e)), call. = FALSE))
            if (length(pred) != length(target)) {
                stop("The opal selectivity function returned the wrong length.", call. = FALSE)
            }
            sum((pred - target)^2)
        }

        # Starting from zeros puts the curve at mean(len_mid) with a 0.5 floor
        # and ceiling, far from a narrow dome, and BFGS stalls there. Start
        # instead from the target's own shape, and keep the best of a few.
        starts = c(list(par_sel[f, seq_len(spec$n)]),
                   ss3_sel_starts(target, len_mid, sel_type_f[f]))
        best = NULL
        for (s0 in starts) {
            s0 = pmin(pmax(s0, lower + 1e-6), upper - 1e-6)
            o = tryCatch(stats::optim(par = s0, fn = objective, method = "L-BFGS-B",
                                      lower = lower, upper = upper,
                                      control = list(maxit = 5000, factr = 1e3)),
                         error = function(e) NULL)
            if (!is.null(o) && is.finite(o$value) &&
                (is.null(best) || o$value < best$value)) best = o
        }
        if (is.null(best)) next
        par_sel[f, seq_len(spec$n)] = best$par
        rmse[f] = sqrt(best$value / length(target))
    }

    list(par_sel = par_sel, n_par_f = n_par_f, rmse = rmse)
}

#' Data-driven starting values for ss3_seed_par_sel()
#'
#' Reads the peak, the half-maximum crossings and the end values off the
#' target curve and maps them onto opal's parameter scales (positions in
#' units of sd(len_mid) about mean(len_mid); widths on the log scale).
#' @keywords internal
ss3_sel_starts = function(target, len_mid, sel_type) {
    mu   = mean(len_mid)
    s    = stats::sd(len_mid)
    tgt  = pmax(target / max(target), 0)
    i_pk = which.max(tgt)
    x_pk = len_mid[i_pk]
    logit = function(p) stats::qlogis(min(max(p, 1e-4), 1 - 1e-4))

    # Half-maximum crossings either side of the peak (fall back to a bin width).
    lo = suppressWarnings(max(which(tgt[seq_len(i_pk)] < 0.5)))
    hi = suppressWarnings(min(which(tgt[i_pk:length(tgt)] < 0.5)) + i_pk - 1L)
    bw = len_mid[2] - len_mid[1]
    w_up = if (is.finite(lo)) max(x_pk - len_mid[lo], bw) else bw
    w_dn = if (is.finite(hi)) max(len_mid[hi] - x_pk, bw) else bw

    if (sel_type == 1L) {
        x50 = if (is.finite(lo)) len_mid[lo] else x_pk
        return(list(c((x50 - mu) / s, log(max(2 * w_up, bw) / s)),
                    c((x50 - mu) / s, log(max(w_up, bw) / s))))
    }
    # Double normal: widths enter as exp(c) * sd^2 = w^2 / log(2) at half max.
    c0 = log(w_up^2 / log(2) / s^2)
    d0 = log(w_dn^2 / log(2) / s^2)
    e0 = logit(tgt[1])
    f0 = logit(tgt[length(tgt)])
    out = list(c((x_pk - mu) / s, -4, c0, d0, e0, f0),
               c((x_pk - mu) / s, -1, c0, d0, e0, f0),
               c((x_pk - mu) / s - 0.25, -4, c0, d0, -8, f0))
    if (sel_type == 3L) {
        # Double Richards: l50_1, log width_1, log shape_1, log gap, log width_2, log shape_2.
        x50 = if (is.finite(lo)) len_mid[lo] else x_pk
        out = list(c((x50 - mu) / s, log(w_up / s), 0,
                     log(max(w_up + w_dn, bw) / s), log(w_dn / s), 0))
    }
    out
}

#_____________________________________________________________________________
#' LOCAL TEST PATCH: let opal_model() take length-at-age as data
#'
#' opal always builds its probability-of-length-at-age matrix (PLA) from the
#' Schnute growth curve. For a model like this one, where SS3's growth
#' (GrowthModel 3, age-specific K) cannot be represented by a single k, that
#' curve misplaces the predicted length compositions even when biology and
#' selectivity are taken from SS3 exactly.
#'
#' This inserts one statement into opal's namespace copy of opal_model(),
#' directly after it computes mean and SD of length-at-age:
#'
#'   if (exists("laa_mu_external", inherits = FALSE)) {
#'       mu_a <- laa_mu_external; sd_a <- laa_sd_external
#'   }
#'
#' It is inert unless the data list carries both vectors, so other fits in the
#' same session are unaffected. It lasts only for the R session, and the
#' installed opal is not modified. Objects built while it is active record a
#' different model signature from stock opal, which flags them as patched.
#' This is a test of the growth explanation, not a supported opal feature; the
#' supported version would be an opal data option analogous to sel_fa_external.
#'
#' @return \code{TRUE}, invisibly. Safe to call more than once.
ss3_patch_opal_external_laa = function() {
    f = utils::getFromNamespace("opal_model", "opal")
    b = as.list(body(f))
    if (any(vapply(b, function(e) any(grepl("laa_mu_external", deparse(e))), logical(1L)))) {
        return(invisible(TRUE))
    }
    i = which(vapply(b, function(e) any(grepl("get_sd_at_age", deparse(e))), logical(1L)))
    if (length(i) != 1L) {
        stop("Could not find the sd_a <- get_sd_at_age(...) statement in opal_model(); ",
             "the patch does not match this opal version.", call. = FALSE)
    }
    ins = quote(if (exists("laa_mu_external", inherits = FALSE)) {
        mu_a <- laa_mu_external
        sd_a <- laa_sd_external
    })
    body(f) = as.call(append(b, list(ins), after = i))
    utils::assignInNamespace("opal_model", f, "opal")
    invisible(TRUE)
}

#' Supply SS3's length-at-age to a patched opal model
#'
#' Requires \code{ss3_patch_opal_external_laa()} in the session. The growth
#' parameters then affect nothing and stay fixed (they are in the default map).
#'
#' @param x an \code{opal_obj}.
#' @param mu_a,sd_a \code{numeric} mean and SD of length at each opal age,
#'   e.g. from \code{ss3_endgrowth_table(rep, dims, basis = "Beg")}.
#' @return an updated \code{opal_obj}.
ss3_use_external_length_at_age = function(x, mu_a, sd_a) {
    stopifnot("mu_a and sd_a must have one value per opal age" =
                  length(mu_a) == x$data$n_age && length(sd_a) == x$data$n_age)
    data = x$data
    data$laa_mu_external = as.numeric(mu_a)
    data$laa_sd_external = as.numeric(sd_a)
    opal::opal_update(x, data = data,
                      metadata = list(length_at_age = "ss3_external (local opal patch)"))
}

#_____________________________________________________________________________
#' Switch an assessment onto estimated opal selectivity, SS3-style
#'
#' Sets up the opal-native selectivity configuration described by
#' \code{cfg$selectivity}: starting values from \code{ss3_seed_par_sel()},
#' the curve columns estimated per opal family (\code{free_par}), the
#' fisheries estimated (\code{estimate}), and mirrored fisheries sharing their
#' source's parameters (\code{mirror}). The defaults in the IO yellowfin config
#' estimate the opal equivalents of exactly the selectivity parameters SS3
#' estimated.
#'
#' @param x an \code{opal_obj}.
#' @param par_sel \code{[n_fishery, n_par]} starting values, typically
#'   \code{ss3_seed_par_sel(...)$par_sel}.
#' @param cfg converter configuration list.
#' @param fix optional named list, \code{fishery = par_sel columns}, of
#'   parameters to hold at their \code{par_sel} values even though
#'   \code{free_par} would estimate them. Mirrors of a fishery follow it.
#' @return an updated \code{opal_obj} with \code{sel_fa_external} removed.
ss3_use_estimated_selectivity = function(x, par_sel, cfg, fix = NULL) {

    sel      = cfg$selectivity
    type_f   = x$data$sel_type_f
    n_f      = nrow(par_sel)
    estimate = as.integer(sel$estimate %||% integer())
    mirror   = sel$mirror %||% list()
    free_par = sel$free_par %||% list()

    # Mirrored fisheries start from, and stay equal to, their source.
    for (f in names(mirror)) {
        src = as.integer(mirror[[f]])
        if (type_f[as.integer(f)] != type_f[src]) {
            stop(sprintf("Fishery %s mirrors %d but has a different curve family.", f, src),
                 call. = FALSE)
        }
        par_sel[as.integer(f), ] = par_sel[src, ]
    }

    map_sel = matrix(NA_integer_, n_f, ncol(par_sel))
    level = 0L
    for (f in estimate) {
        cols = as.integer(free_par[[as.character(type_f[f])]] %||% integer())
        cols = setdiff(cols, as.integer(fix[[as.character(f)]] %||% integer()))
        if (!length(cols)) next
        map_sel[f, cols] = level + seq_along(cols)
        level = level + length(cols)
    }
    for (f in names(mirror)) map_sel[as.integer(f), ] = map_sel[as.integer(mirror[[f]]), ]

    data = x$data
    data$sel_fa_external = NULL
    parameters = x$parameters
    parameters$par_sel = par_sel
    map = x$map
    if (is.null(map)) map = opal::get_map(x$parameters)
    map$par_sel = factor(as.vector(map_sel))

    opal::opal_update(x, data = data, parameters = parameters, map = map,
                      metadata = list(selectivity_route = "opal_estimated"))
}
