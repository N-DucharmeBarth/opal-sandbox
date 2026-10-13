# Nicholas Ducharme-Barth
# 2026/10/13
# Two routes for getting biology out of SS3 and into opal.
#
# WHY TWO ROUTES
#
# opal always computes mean length-at-age from the Schnute parameters, then
# builds a probability-of-length-at-age matrix (PLA) from it. Length-basis
# vectors in the data list (weight, maturity, fecundity) are converted to age
# through that PLA. So with the control-file route, any error in the Schnute
# approximation to SS3's growth propagates into weight-at-age, maturity-at-age
# and therefore spawning biomass.
#
# This SS3 model uses GrowthModel 3 (age-specific K, nine multipliers over ages
# 2-10). A single-k Schnute curve cannot reproduce that, so the approximation
# error is real and it lands on the quantity of interest.
#
# The endgrowth route sidesteps it. opal's resolve_bio_vector() returns a
# length-n_age vector unchanged, so supplying weight, maturity, spawning
# potential and M on an age basis from SS3's own endgrowth table means the PLA
# is never used for biology. Growth parameters then affect only selectivity,
# through sel_length() and the PLA. Supply sel_fa_external as well and growth
# becomes entirely inert - at which point the growth parameters should be
# mapped off, since nothing in the likelihood depends on them.
#
# WHAT THIS COSTS
#
# Age-basis biology is fixed data, not a function of estimated growth. You can
# no longer estimate growth internally, and you cannot propagate growth
# uncertainty. For a bridging exercise against a known operating model that is
# the right trade; for an assessment where growth is being estimated it is not.
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
#' Extract a tidy age-basis table from SS3's endgrowth output
#'
#' endgrowth can carry several rows per age (sex, season, growth morph,
#' platoon). Rows are filtered to one sex and one season; if duplicate ages
#' remain, they are averaged and a warning is issued. For a platoon model this
#' unweighted average is wrong - the correct aggregation weights platoons by
#' abundance, which endgrowth does not carry - so the warning should be acted
#' on rather than ignored.
#'
#' @param rep output of \code{r4ss::SS_output}.
#' @param dims output of \code{ss3_dims_to_opal}.
#' @param basis \code{"Beg"} or \code{"Mid"}; which within-timestep columns to
#'   read.
#' @return a \code{data.table} with one row per SS3 age in \code{dims$.ss3_ages}.
ss3_endgrowth_table = function(rep, dims, basis = c("Beg", "Mid")) {

    basis = match.arg(basis)
    eg = rep$endgrowth
    if (is.null(eg)) {
        stop("Report.sso contains no endgrowth table.", call. = FALSE)
    }
    eg = data.table::as.data.table(eg)

    pick_col = function(candidates, what, required = TRUE) {
        hit = intersect(candidates, names(eg))
        if (!length(hit)) {
            if (!required) return(NA_character_)
            stop(sprintf("endgrowth has no %s column. Tried: %s.\nAvailable: %s",
                         what, paste(candidates, collapse = ", "),
                         paste(names(eg), collapse = ", ")), call. = FALSE)
        }
        hit[1]
    }

    age_col = pick_col(c("Age_Beg", "int_Age", "Age", "Real_Age"), "age")
    len_col = pick_col(c(paste0("Len_", basis), "Len_Beg", "Len_Mid"), "length-at-age")
    sd_col  = pick_col(c(paste0("SD_",  basis), "SD_Beg",  "SD_Mid"),  "SD of length-at-age")
    wt_col  = pick_col(c(paste0("Wt_",  basis), "Wt_Beg",  "Wt_Mid"),  "weight-at-age")
    mat_col = pick_col(c("Mat", "Mat_F", "Age_Mat"), "maturity-at-age")
    m_col   = pick_col(c("M", "M_nums", "Nat_Mort"), "natural mortality")
    mf_col  = pick_col(c("Mat*Fecund", "Mat*Fec", "Mat_Fecund"),
                       "maturity x fecundity", required = FALSE)
    fec_col = pick_col(c("Fecund", "Fec", "Fecundity"), "fecundity", required = FALSE)

    # One sex and one season. Keep females where sexes are distinguished,
    # because SS3's spawning output is female-based.
    sex_col = intersect(c("Sex", "Gender"), names(eg))[1]
    if (!is.na(sex_col)) eg = eg[get(sex_col) == min(get(sex_col))]
    seas_col = intersect(c("Seas", "Season"), names(eg))[1]
    if (!is.na(seas_col)) eg = eg[get(seas_col) == min(get(seas_col))]

    cols = c(age = age_col, length = len_col, sd = sd_col, weight = wt_col,
             maturity = mat_col, M = m_col)
    if (!is.na(mf_col))  cols = c(cols, mat_fec = mf_col)
    if (!is.na(fec_col)) cols = c(cols, fecundity = fec_col)

    out = eg[, cols, with = FALSE]
    data.table::setnames(out, names(cols))
    out[, age := as.numeric(age)]

    if (anyDuplicated(out$age)) {
        warning(paste0("endgrowth has multiple rows per age after filtering sex and ",
                       "season (growth morphs or platoons). Averaging them without ",
                       "abundance weights, which is only correct for a single-platoon ",
                       "model. Check N_Morphs_Within_GrowthPattern."), call. = FALSE)
        out = out[, lapply(.SD, mean, na.rm = TRUE), by = age]
    }

    i = match(dims$.ss3_ages, out$age)
    if (anyNA(i)) {
        stop(sprintf("endgrowth does not cover SS3 ages %s.",
                     paste(dims$.ss3_ages[is.na(i)], collapse = ", ")), call. = FALSE)
    }
    out[i]
}

#_____________________________________________________________________________
#' Fit opal's Schnute growth and CV to SS3's length-at-age
#'
#' L1 and L2 are pinned at SS3's mean length at the reference ages, then k is
#' fitted by least squares to the full length-at-age vector and, optionally,
#' CV1 and CV2 are fitted to SS3's SD-at-age. Fitting the CVs matters when the
#' PLA is still used for selectivity, since it controls how much the age-based
#' SS3 selectivity smears across length bins.
#'
#' @param ages opal age vector.
#' @param laa SS3 mean length at those ages.
#' @param sd_a SS3 SD of length at those ages, or \code{NULL}.
#' @param A1,A2 Schnute reference ages.
#' @param L1,L2 length at the reference ages.
#' @param fit_cv \code{logical}; fit CV1 and CV2 to \code{sd_a}.
#' @param cv_fallback CV used for both ends when \code{sd_a} is unavailable.
#' @return a \code{list} with \code{k}, \code{CV1}, \code{CV2}, fit residuals
#'   and the SS3 / opal vectors for plotting.
fit_schnute_growth = function(ages, laa, sd_a = NULL, A1, A2, L1, L2,
                              fit_cv = TRUE, cv_fallback = 0.1) {

    mu_of = function(k) {
        L1 + (L2 - L1) * (1 - exp(-k * (ages - A1))) / (1 - exp(-k * (A2 - A1)))
    }
    o  = stats::optimize(function(lk) sum((mu_of(exp(lk)) - laa)^2),
                         interval = c(log(1e-4), log(5)))
    k  = exp(o$minimum)
    mu = mu_of(k)

    CV1 = cv_fallback
    CV2 = cv_fallback
    sd_opal = NULL
    sd_rmse = NA_real_

    if (!is.null(sd_a) && isTRUE(fit_cv) && all(is.finite(sd_a)) && all(sd_a > 0)) {
        sd_of = function(p) mu * (p[1] + (mu - L1) / (L2 - L1) * (p[2] - p[1]))
        start = c(sd_a[1] / mu[1], sd_a[length(sd_a)] / mu[length(mu)])
        oo = stats::optim(start, function(p) sum((sd_of(p) - sd_a)^2),
                          method = "L-BFGS-B",
                          lower = c(1e-4, 1e-4), upper = c(2, 2))
        CV1 = oo$par[1]
        CV2 = oo$par[2]
        sd_opal = sd_of(oo$par)
        sd_rmse = sqrt(mean((sd_opal - sd_a)^2))
    }

    list(k             = k,
         CV1           = CV1,
         CV2           = CV2,
         rmse          = sqrt(mean((mu - laa)^2)),
         max_abs_error = max(abs(mu - laa)),
         sd_rmse       = sd_rmse,
         source        = "fitted_to_endgrowth",
         laa_ss3       = laa,
         laa_opal      = mu,
         sd_ss3        = sd_a,
         sd_opal       = sd_opal)
}

#_____________________________________________________________________________
#' Biology from SS3's endgrowth table (age basis)
#'
#' @inheritParams ss3_bio_to_opal
#' @return a \code{list} of biology elements.
ss3_bio_from_endgrowth = function(ss3, cfg, dims) {

    opts  = cfg$bio$endgrowth %||% list()
    basis = opts$length_basis %||% "Beg"
    eg    = ss3_endgrowth_table(ss3$rep, dims, basis = basis)
    n_age = dims$n_age

    mg   = ss3_pick(ss3$ctl, c("MG_parms", "MG_parms_ctl"), what = "MG_parms")
    lw_a = unname(ss3_parm_value(mg, "^Wtlen_1", n_expected = 1L))
    lw_b = unname(ss3_parm_value(mg, "^Wtlen_2", n_expected = 1L))
    frac_female = ss3_spawning_female_fraction(ss3, cfg)

    A1 = as.integer(cfg$growth$A1 %||% dims$min_age)
    A2 = as.integer(cfg$growth$A2 %||% (dims$min_age + n_age - 1L))
    L1 = eg$length[which(dims$age_a == A1)]
    L2 = eg$length[which(dims$age_a == A2)]

    growth_fit = fit_schnute_growth(
        ages   = dims$age_a,
        laa    = eg$length,
        sd_a   = eg$sd,
        A1 = A1, A2 = A2, L1 = L1, L2 = L2,
        fit_cv = isTRUE(cfg$growth$fit_cv %||% TRUE)
    )

    # Weight-at-age, straight from SS3, in SS3's own weight units.
    weight = if (isTRUE(opts$weight %||% TRUE)) as.numeric(eg$weight) else NULL

    # Maturity-at-age.
    maturity = as.numeric(eg$maturity)

    # Spawning potential. opal uses spawning_potential directly and ignores
    # sex_ratio when it is supplied, so any female fraction is folded in here
    # (see ss3_spawning_female_fraction(); it is 1 for a single-sex SS3 model).
    # Fecundity-at-age. opal_model() resolves fecundity whether or not
    # spawning_potential is supplied, so it is always carried; it only drives
    # spawning potential when Mat*Fecund is absent, in which case opal forms
    # spawning_potential = sex_ratio * maturity * fecundity itself.
    fecundity = if (!is.null(eg$fecundity)) as.numeric(eg$fecundity) else
        as.numeric(eg$weight)
    spawning_potential = NULL
    if (isTRUE(opts$spawning_potential %||% TRUE) && !is.null(eg$mat_fec)) {
        spawning_potential = as.numeric(eg$mat_fec) * frac_female
    }

    # endgrowth reports M as an annual rate; SS3 applies it over the season
    # duration, so it is rescaled to the model timestep like the control value.
    M = if (isTRUE(opts$m %||% TRUE)) as.numeric(eg$M) else
        as.numeric(ss3_pick(ss3$ctl, c("natM", "natM_vector", "M_at_age"),
                            what = "natM"))[dims$.ss3_ages + 1L]
    M = M / ss3_m_rate_divisor(cfg, dims)

    list(
        source             = "endgrowth",
        M                  = M,
        lw_a               = lw_a,
        lw_b               = lw_b,
        weight             = weight,
        maturity           = maturity,
        fecundity          = fecundity,
        spawning_potential = spawning_potential,
        sex_ratio          = rep(frac_female, n_age),
        A1                 = A1,
        A2                 = A2,
        L1                 = L1,
        L2                 = L2,
        k                  = growth_fit$k,
        CV1                = growth_fit$CV1,
        CV2                = growth_fit$CV2,
        growth_fit         = growth_fit
    )
}

#_____________________________________________________________________________
#' Biology from the SS3 control file (length basis, conventional translation)
#'
#' Weight and maturity are supplied on a length basis and converted to age by
#' opal's own PLA, so the Schnute approximation propagates into spawning
#' biomass. Use \code{cfg$bio$source = "endgrowth"} when that matters.
#'
#' @inheritParams ss3_bio_to_opal
#' @return a \code{list} of biology elements.
ss3_bio_from_control = function(ss3, cfg, dims) {

    ss3_ages = dims$.ss3_ages
    n_age    = dims$n_age

    m_ctl = as.numeric(ss3_pick(ss3$ctl, c("natM", "natM_vector", "M_at_age"),
                                what = "natM"))
    M = m_ctl[ss3_ages + 1L] / ss3_m_rate_divisor(cfg, dims)
    stopifnot("M has the wrong length" = length(M) == n_age)

    mg   = ss3_pick(ss3$ctl, c("MG_parms", "MG_parms_ctl"), what = "MG_parms")
    lw_a = unname(ss3_parm_value(mg, "^Wtlen_1", n_expected = 1L))
    lw_b = unname(ss3_parm_value(mg, "^Wtlen_2", n_expected = 1L))

    A1 = as.integer(cfg$growth$A1 %||% dims$min_age)
    A2 = as.integer(cfg$growth$A2 %||% (dims$min_age + n_age - 1L))

    L_at_Amin = unname(ss3_parm_value(mg, "^L_at_Amin", n_expected = 1L))
    L_at_Amax = unname(ss3_parm_value(mg, "^L_at_Amax", n_expected = 1L))
    # VonBert_K is annual and, like M, applied over the season duration; opal's
    # Schnute k is per timestep (age unit), so it takes the same divisor. Any
    # age-specific K multipliers (GrowthModel 3) are not represented.
    k         = unname(ss3_parm_value(mg, "^VonBert_K", n_expected = 1L)) /
        ss3_m_rate_divisor(cfg, dims)
    age_L1    = as.numeric(ss3_pick(ss3$ctl, c("Growth_Age_for_L1"), default = A1))
    age_L2    = as.numeric(ss3_pick(ss3$ctl, c("Growth_Age_for_L2"), default = A2))

    Linf = if (age_L2 >= 999) L_at_Amax else
        L_at_Amin + (L_at_Amax - L_at_Amin) / (1 - exp(-k * (age_L2 - age_L1)))
    L1 = L_at_Amin + (Linf - L_at_Amin) * (1 - exp(-k * (A1 - age_L1)))
    L2 = L_at_Amin + (Linf - L_at_Amin) * (1 - exp(-k * (A2 - age_L1)))

    CV1 = unname(ss3_parm_value(mg, "^CV_young", n_expected = 1L))
    CV2 = unname(ss3_parm_value(mg, "^CV_old",   n_expected = 1L))

    growth_fit = list(k = k, CV1 = CV1, CV2 = CV2,
                      rmse = NA_real_, max_abs_error = NA_real_,
                      sd_rmse = NA_real_, source = "control_file_vb",
                      laa_ss3 = NULL, laa_opal = NULL,
                      sd_ss3 = NULL, sd_opal = NULL)

    # If Report.sso is available, report how far the control-file curve sits
    # from SS3's realised length-at-age. This is the cost of this route.
    if (!is.null(ss3$rep)) {
        eg = tryCatch(ss3_endgrowth_table(ss3$rep, dims), error = function(e) NULL)
        if (!is.null(eg)) {
            mu = L1 + (L2 - L1) * (1 - exp(-k * (dims$age_a - A1))) /
                (1 - exp(-k * (A2 - A1)))
            growth_fit$laa_ss3       = eg$length
            growth_fit$laa_opal      = mu
            growth_fit$rmse          = sqrt(mean((mu - eg$length)^2))
            growth_fit$max_abs_error = max(abs(mu - eg$length))
        }
    }

    len_mid = dims$.len_bins + dims$len_bin_width / 2

    mat_opt = as.integer(ss3_pick(ss3$ctl, c("maturity_option", "Maturity_option"),
                                  default = 1L))
    if (mat_opt %in% c(3L, 4L)) {
        mat_vec  = as.numeric(ss3_pick(ss3$ctl, c("Age_Maturity", "Age_Maturity_Fem",
                                                  "maturity_age"),
                                       what = "age maturity vector"))
        maturity = mat_vec[ss3_ages + 1L]
    } else {
        mat50    = unname(ss3_parm_value(mg, "^Mat50", n_expected = 1L))
        mslope   = unname(ss3_parm_value(mg, "^Mat_slope", n_expected = 1L))
        maturity = 1 / (1 + exp(mslope * (len_mid - mat50)))
    }

    fec_a = tryCatch(unname(ss3_parm_value(mg, "Eggs.*inter", n_expected = 1L)),
                     error = function(e) 1)
    fec_b = tryCatch(unname(ss3_parm_value(mg, "Eggs.*slope", n_expected = 1L)),
                     error = function(e) 0)
    wt_at_len = lw_a * len_mid^lw_b
    fecundity = wt_at_len * (fec_a + fec_b * wt_at_len)

    frac_female = ss3_spawning_female_fraction(ss3, cfg)

    list(
        source             = "control",
        M                  = M,
        lw_a               = lw_a,
        lw_b               = lw_b,
        weight             = NULL,
        maturity           = maturity,
        fecundity          = fecundity,
        spawning_potential = NULL,
        sex_ratio          = rep(frac_female, n_age),
        A1 = A1, A2 = A2, L1 = L1, L2 = L2,
        k = k, CV1 = CV1, CV2 = CV2,
        growth_fit = growth_fit
    )
}

#_____________________________________________________________________________
#' Unfished spawning potential per recruit from a converted data list
#'
#' Equilibrium survivorship under M alone, with a plus group at the oldest age,
#' times spawning potential at age. Used to reconcile against SS3's reported
#' unfished spawning biomass per recruit, which tests the female-fraction
#' convention and the units of M in one number.
#'
#' @param data an opal data list.
#' @return a single \code{numeric}.
opal_spr0 = function(data) {

    n_age = data$n_age
    M     = data$M
    if (length(M) != n_age) {
        stop("opal_spr0() needs M on an age basis.", call. = FALSE)
    }

    sp = data$spawning_potential
    if (is.null(sp)) {
        if (length(data$maturity) != n_age || length(data$fecundity) != n_age) {
            # Length-basis biology cannot be resolved without the PLA, which
            # only exists inside opal_model(). Skip rather than guess.
            return(NA_real_)
        }
        sr = data$sex_ratio %||% rep(1, n_age)
        sp = sr * data$maturity * data$fecundity
    }

    n = numeric(n_age)
    n[1] = 1
    for (a in 2:n_age) n[a] = n[a - 1] * exp(-M[a - 1])
    n[n_age] = n[n_age] / (1 - exp(-M[n_age]))

    sum(n * sp)
}

#' Divisor that turns SS3's annual rates into per-timestep rates
#'
#' SS3 reads M as an annual rate and applies it over each season's duration
#' (months_per_seas / 12). opal's timestep is one SS3 year of nseas seasons, so
#' the divisor is 12 / (nseas * months_per_seas): 4 for this model's
#' three-month "years". This holds for the control natM vector and for
#' endgrowth's M column alike; it was confirmed here from the virgin
#' numbers-at-age, where survival from age a to a + 1 is exp(-M_a / 4).
#' \code{cfg$bio$m_rate_divisor} overrides the derived value.
#' @keywords internal
ss3_m_rate_divisor = function(cfg, dims) {
    as.numeric(cfg$bio$m_rate_divisor %||% dims$.ts_per_year)
}

#' Female fraction to fold into spawning potential
#'
#' In a single-sex SS3 model every fish carries sex 1 and spawning biomass is
#' sum(N x Mat*Fecund) with no female fraction applied; fracfemale is not used.
#' Folding it in would halve SSB (checked against SSB_Virgin for this model).
#' With two sexes the numbers are split by sex upstream, so opal, which tracks
#' one sex, needs fracfemale. \code{cfg$bio$endgrowth$apply_frac_female}
#' overrides the choice when set.
#' @keywords internal
ss3_spawning_female_fraction = function(ss3, cfg) {
    apply = cfg$bio$endgrowth$apply_frac_female %||%
        (as.integer(ss3$dat$Nsexes %||% 1L) == 2L)
    if (!isTRUE(apply)) return(1)
    as.numeric(ss3_pick(ss3$ctl, c("fracfemale", "fracfem"), default = 1))
}

#' One value from SS3's derived quantities, by label
#' @keywords internal
ss3_derived_value = function(rep, labels) {
    dq = if (is.null(rep)) NULL else rep$derived_quants
    if (is.null(dq)) return(NA_real_)
    lab = intersect(c("Label", "LABEL"), names(dq))[1]
    val = intersect(c("Value", "value"), names(dq))[1]
    if (is.na(lab) || is.na(val)) return(NA_real_)
    hit = which(dq[[lab]] %in% labels)
    if (!length(hit)) return(NA_real_)
    as.numeric(dq[[val]][hit[1]])
}

#' SS3's unfished spawning biomass per recruit at the first opal age
#'
#' opal's SPR0 counts survivorship from its first age class with N = 1, so the
#' comparable SS3 number is SSB0 divided by the virgin numbers at the first
#' SS3 age carried, not by SS3's age-0 R0. The two differ by survival from
#' age 0 to the first carried age; when that age is 0 they coincide.
#'
#' @param rep output of \code{r4ss::SS_output}.
#' @param dims output of \code{ss3_dims_to_opal}; \code{NULL} means age 0.
#' @return a single \code{numeric}.
#' @keywords internal
ss3_spr0 = function(rep, dims = NULL) {
    if (is.null(rep)) return(NA_real_)
    ssb = ss3_derived_value(rep, c("SSB_Unfished", "SPB_Unfished", "SSB_Virgin", "SPB_Virgin"))
    if (!is.finite(ssb)) return(NA_real_)

    first_age = if (is.null(dims)) 0L else min(dims$.ss3_ages)
    na = rep$natage
    if (!is.null(na) && as.character(first_age) %in% names(na)) {
        bm  = intersect(c("Beg/Mid", "BegMid"), names(na))[1]
        sel = na$Era == "VIRG"
        if (!is.na(bm)) sel = sel & na[[bm]] == "B"
        n_first = sum(as.numeric(na[sel, as.character(first_age)]))
        if (is.finite(n_first) && n_first > 0) return(ssb / n_first)
    }
    if (first_age != 0L) return(NA_real_)
    recr = ss3_derived_value(rep, c("Recr_Unfished", "Recr_Virgin"))
    if (!is.finite(recr) || recr <= 0) return(NA_real_)
    ssb / recr
}

#' SS3's recruitment bias-adjustment ramp, aligned to the opal window
#'
#' SS3 recruits are R_BH x exp(dev - 0.5 b_y sigma_R^2) with b_y from the bias
#' ramp; opal's bias_adj_y plays the same role. Taken from Report.sso's
#' recruitment table, so it is exactly what SS3 used.
#'
#' @return \code{numeric} vector of length \code{n_year}, or \code{NULL}.
#' @keywords internal
ss3_bias_adj = function(rep, dims) {
    rc = if (is.null(rep)) NULL else rep$recruit
    if (is.null(rc) || !all(c("Yr", "biasadjuster") %in% names(rc))) return(NULL)
    hit = match(dims$.first_ts:dims$.last_ts, rc$Yr)
    if (anyNA(hit)) return(NULL)
    as.numeric(rc$biasadjuster[hit])
}
