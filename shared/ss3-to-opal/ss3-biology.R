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

    # Weight-at-age, straight from SS3, in SS3's own weight units. opal uses
    # it to turn catch in weight into numbers. cfg$bio$endgrowth$weight_basis
    # ("Beg" or "Mid", default length_basis) picks the timing within the
    # timestep: "Mid" is closer to the weight of SS3's catch, which is taken
    # through the timestep, where opal removes it at the start.
    weight_basis = opts$weight_basis %||% basis
    eg_wt = if (identical(weight_basis, basis)) eg else
        ss3_endgrowth_table(ss3$rep, dims, basis = weight_basis)
    weight = if (isTRUE(opts$weight %||% TRUE)) as.numeric(eg_wt$weight) else NULL

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
        weight_basis       = weight_basis,
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
    # Nsexes = 2 splits the numbers by sex upstream; Nsexes = -1 (SS3 3.30) is
    # a single sex whose spawning biomass SS3 multiplies by the female fraction
    apply = cfg$bio$endgrowth$apply_frac_female %||%
        (as.integer(ss3$dat$Nsexes %||% 1L) %in% c(2L, -1L))
    if (!isTRUE(apply)) return(1)
    ff = ss3_pick(ss3$ctl, c("fracfemale", "fracfem"), default = NULL)
    if (is.null(ff)) {
        mg = ss3_pick(ss3$ctl, c("MG_parms"), default = NULL)
        if (!is.null(mg) && any(grepl("^FracFemale", rownames(mg)))) {
            ff = ss3_parm_value(mg, "^FracFemale", n_expected = 1L)
        }
    }
    as.numeric(ff %||% 1)
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
    # aligned like rdev_y (see ss3_recdev_shift()); a year past the end of
    # SS3's table (the first forecast year) takes no bias adjustment
    yrs = dims$.first_ts:dims$.last_ts + ss3_recdev_shift(dims)
    hit = match(yrs, rc$Yr)
    if (anyNA(hit[-length(hit)])) return(NULL)
    out = as.numeric(rc$biasadjuster[hit])
    out[is.na(out)] = 0
    out
}

#' SS3's bias adjustment for the initial age structure (cfg$recdev$initial)
#' @keywords internal
ss3_init_bias_adj = function(rep, dims) {
    rc = if (is.null(rep)) NULL else rep$recruit
    if (is.null(rc) || !all(c("Yr", "biasadjuster") %in% names(rc))) return(NULL)
    hit = match(dims$.first_ts - dims$.ss3_ages, rc$Yr)
    out = as.numeric(rc$biasadjuster[hit])
    out[is.na(out)] = 0
    out
}

#_____________________________________________________________________________
#' Replace opal's Schnute growth with SS3's age-specific K (GrowthModel 3)
#'
#' Requires opal >= 0.1.1, where get_growth() implements Stock Synthesis
#' GrowthModel 3: cumulative K multipliers from given ages, Linf derived from
#' the oldest-age K, and SS3's plus-group length approximation (Methot and
#' Wetzel 2013). Everything is read from the SS3 files, nothing is fitted:
#'
#'   L1, L2      L_at_Amin / L_at_Amax, at Growth_Age_for_L1 / _L2
#'   log_k       log(VonBert_K / divisor); SS3's K is annual and applied over
#'               the season duration, like M (see ss3_m_rate_divisor())
#'   log_k_mult  log of the Age_K multipliers (par file values when read)
#'   k_mult_ages the control file's K-multiplier ages
#'   CV1, CV2    CV_young / CV_old
#'
#' Applies to either biology route. With it, opal's length-at-age equals SS3's
#' Len_Beg to Report.sso print precision, so the PLA (and through it the
#' predicted length compositions, and on the control route weight and
#' spawning potential) no longer carries a growth approximation.
#'
#' @param bio output of \code{ss3_bio_from_endgrowth()} or
#'   \code{ss3_bio_from_control()}.
#' @inheritParams ss3_bio_to_opal
#' @return \code{bio} with growth elements replaced and \code{k_mult_ages},
#'   \code{log_k_mult} and \code{plus_group_adjust} added.
ss3_apply_age_k_growth = function(bio, ss3, cfg, dims) {

    if (utils::packageVersion("opal") < "0.1.1") {
        stop("cfg$growth$model = 'ss3_age_k' needs opal >= 0.1.1 (age-specific K in get_growth()).",
             call. = FALSE)
    }
    ctl = ss3$ctl
    if (!isTRUE(as.integer(ss3_pick(ctl, c("GrowthModel"), default = NA)) == 3L)) {
        stop("cfg$growth$model = 'ss3_age_k' requires an SS3 model with GrowthModel 3.",
             call. = FALSE)
    }

    # Parameter values: the par file holds the run's values; the control file
    # only the INIT values the run started from.
    mg = ss3_pick(ss3$par %||% list(), c("MG_parms"), default = NULL) %||%
        ss3_pick(ctl, c("MG_parms"), what = "MG_parms")
    mg = as.data.frame(mg)

    age_L1 = as.integer(ss3_pick(ctl, c("Growth_Age_for_L1"), what = "Growth_Age_for_L1"))
    age_L2 = as.integer(ss3_pick(ctl, c("Growth_Age_for_L2"), what = "Growth_Age_for_L2"))
    if (age_L2 >= 999L) {
        stop("Growth_Age_for_L2 = 999 (L2 as Linf) is not supported by opal's get_growth().",
             call. = FALSE)
    }
    # opal age index = min_age + position; convert SS3 ages to opal ages.
    to_opal_age = function(a) dims$age_a[match(a, dims$.ss3_ages)]
    A1 = to_opal_age(age_L1)
    A2 = to_opal_age(age_L2)
    if (anyNA(c(A1, A2))) {
        stop("SS3 growth reference ages fall outside cfg$age$keep_ss3_ages.", call. = FALSE)
    }

    k_ages_ss3 = as.integer(ss3_pick(ctl, c("Age_K_points"), what = "Age_K_points"))
    k_mult     = unname(ss3_parm_value(mg, "^Age_K", n_expected = length(k_ages_ss3)))
    k_ages     = to_opal_age(k_ages_ss3)
    if (anyNA(k_ages)) {
        stop("SS3 K-multiplier ages fall outside cfg$age$keep_ss3_ages.", call. = FALSE)
    }

    L1  = unname(ss3_parm_value(mg, "^L_at_Amin", n_expected = 1L))
    L2  = unname(ss3_parm_value(mg, "^L_at_Amax", n_expected = 1L))
    k   = unname(ss3_parm_value(mg, "^VonBert_K", n_expected = 1L)) / ss3_m_rate_divisor(cfg, dims)
    CV1 = unname(ss3_parm_value(mg, "^CV_young", n_expected = 1L))
    CV2 = unname(ss3_parm_value(mg, "^CV_old",   n_expected = 1L))
    plus_adjust = isTRUE(cfg$growth$plus_group_adjust %||% TRUE)

    mu = opal::get_growth(dims$n_age, A1, A2, L1, L2, log(k), min_age = dims$min_age,
                          log_k_mult = log(k_mult), k_mult_ages = k_ages,
                          plus_group_adjust = plus_adjust)

    growth_fit = list(k = k, CV1 = CV1, CV2 = CV2,
                      rmse = NA_real_, max_abs_error = NA_real_, sd_rmse = NA_real_,
                      source = "ss3_age_k", laa_ss3 = NULL, laa_opal = mu,
                      sd_ss3 = NULL, sd_opal = NULL,
                      k_mult = k_mult, k_mult_ages = k_ages)
    if (!is.null(ss3$rep)) {
        eg = tryCatch(ss3_endgrowth_table(ss3$rep, dims, basis = "Beg"), error = function(e) NULL)
        if (!is.null(eg)) {
            growth_fit$laa_ss3       = eg$length
            growth_fit$rmse          = sqrt(mean((mu - eg$length)^2))
            growth_fit$max_abs_error = max(abs(mu - eg$length))
            growth_fit$max_rel_error = max(abs(mu / eg$length - 1))
        }
    }

    bio$A1 = A1
    bio$A2 = A2
    bio$L1 = L1
    bio$L2 = L2
    bio$k  = k
    bio$CV1 = CV1
    bio$CV2 = CV2
    bio$log_k_mult = log(k_mult)
    bio$k_mult_ages = k_ages
    bio$plus_group_adjust = plus_adjust
    bio$growth_fit = growth_fit
    bio
}

#_____________________________________________________________________________
#' SS3 growth cessation (GrowthModel 8), carried exactly into opal
#'
#' opal's \code{growth_model = "cessation"} is Eq. 2 of Maunder et al.
#' (2018), the curve SS3 GrowthModel 8 uses. The SS3 parameters map onto it
#' directly:
#'
#'   L_at_Amin (SS3 age AFIX) -> L0
#'   rmax (labelled VonBert_K in the control file) -> rmax
#'   Cessation               -> k_gc
#'   Linf (Growth_Age_for_L2 = 999)
#'                           -> t50 = log(exp((Linf - L0) k / rmax) - 1) / k,
#'                              as SS3 computes it (SS_biofxn.tpl)
#'
#' SS3 evaluates the curve at real age minus AFIX, so opal's age offset is the
#' first kept SS3 age minus 1 minus AFIX (opal evaluates growth at positions
#' 1, ..., n_age). With it, opal's length-at-age equals SS3's Len_Beg to
#' Report.sso print precision, except SS3's plus-group length (SS3 averages
#' lengths over the ages accumulated in the plus group; opal does not).
#'
#' The SD of length-at-age follows SS3's CV_Growth_Pattern:
#'   0 (CV linear in length): opal growth_sd_model "cv"
#'   2 (SD linear in length): opal growth_sd_model "sd" (opal >= 0.1.2)
#' With Growth_Age_for_L2 = 999 SS3's slope runs from L0 to Linf, which is not
#' a modelled length. opal's upper reference length is the mean length at the
#' oldest age, so the value there is set to SS3's linear form evaluated at that
#' length. opal then matches SS3 exactly at every age below the oldest, while
#' SS3 holds the oldest age at CV_old / SD_old itself (a difference of a few
#' tenths of a percent).
#'
#' @param bio output of \code{ss3_bio_from_endgrowth()} or
#'   \code{ss3_bio_from_control()}.
#' @inheritParams ss3_bio_to_opal
#' @return \code{bio} with the growth elements replaced and
#'   \code{growth_model}, \code{growth_age_offset}, \code{L0}, \code{rmax},
#'   \code{k_gc}, \code{t50} and \code{growth_sd_model} added (plus
#'   \code{SD1}, \code{SD2} under "sd").
ss3_apply_cessation_growth = function(bio, ss3, cfg, dims) {

    if (!"sd_model" %in% names(formals(opal::get_sd_at_age))) {
        stop("cfg$growth$model = 'cessation' needs opal >= 0.1.2 (growth cessation and SD-linear growth variability).",
             call. = FALSE)
    }
    ctl = ss3$ctl
    if (!isTRUE(as.integer(ss3_pick(ctl, c("GrowthModel"), default = NA)) == 8L)) {
        stop("cfg$growth$model = 'cessation' requires an SS3 model with GrowthModel 8.",
             call. = FALSE)
    }
    if (!isTRUE(all.equal(ss3_m_rate_divisor(cfg, dims), 1))) {
        stop("cfg$growth$model = 'cessation' supports models whose SS3 year is the model timestep only.",
             call. = FALSE)
    }

    # Parameter values: the par file holds the run's estimates (ESTIM); the
    # control file only the INIT values the run started from.
    mg_par = ss3_pick(ss3$par %||% list(), c("MG_parms"), default = NULL)
    mg = as.data.frame(mg_par %||% ss3_pick(ctl, c("MG_parms"), what = "MG_parms"))
    col = if (!is.null(mg_par) && "ESTIM" %in% names(mg)) "ESTIM" else "INIT"
    val = function(pattern) unname(ss3_parm_value(mg, pattern, column = col, n_expected = 1L))

    afix  = as.integer(ss3_pick(ctl, c("Growth_Age_for_L1"), what = "Growth_Age_for_L1"))
    afix2 = as.integer(ss3_pick(ctl, c("Growth_Age_for_L2"), what = "Growth_Age_for_L2"))
    if (afix2 < 999L) {
        stop("GrowthModel 8 with Growth_Age_for_L2 < 999 (L at Amax rather than Linf) is not supported.",
             call. = FALSE)
    }
    if (dims$.ss3_ages[1] < afix) {
        stop("SS3 ages below Growth_Age_for_L1 grow linearly in SS3; keep_ss3_ages must start at or above it.",
             call. = FALSE)
    }

    L0   = val("^L_at_Age|^L_at_Amin")
    Linf = val("^Linf")
    rmax = val("^rmax|^VonBert_K")
    k_gc = val("^Cessation")
    t50  = log(exp((Linf - L0) * k_gc / rmax) - 1) / k_gc
    age_offset = dims$.ss3_ages[1] - 1L - afix

    # opal positions of the reference ages: SS3 age AFIX and the oldest age
    n_age = dims$n_age
    A1_pos = match(afix, dims$.ss3_ages)
    A2_pos = n_age
    mu = opal::get_growth(n_age, A1_pos, A2_pos, growth_model = "cessation",
                          L0 = L0, log_rmax = log(rmax), log_k_gc = log(k_gc),
                          t50 = t50, age_offset = age_offset)

    # Growth variability at the reference lengths
    cv_pattern = as.integer(ss3_pick(ctl, c("CV_Growth_Pattern"), what = "CV_Growth_Pattern"))
    v1 = val("^CV_young|^SD_young")
    v2 = val("^CV_old|^SD_old")
    # SS3's linear form (in CV or SD) at the oldest age's mean length
    v2_at = v1 + (mu[A2_pos] - L0) / (Linf - L0) * (v2 - v1)
    if (identical(cv_pattern, 2L)) {
        sd_model = "sd"
        sd_a = opal::get_sd_at_age(mu, mu[A1_pos], mu[A2_pos], sd_model = "sd",
                                   log_SD1 = log(v1), log_SD2 = log(v2_at),
                                   A1 = A1_pos, A2 = A2_pos)
    } else if (identical(cv_pattern, 0L)) {
        sd_model = "cv"
        sd_a = opal::get_sd_at_age(mu, mu[A1_pos], mu[A2_pos], log(v1), log(v2_at))
    } else {
        stop(sprintf("CV_Growth_Pattern %d is not supported (opal has patterns 0 and 2).", cv_pattern),
             call. = FALSE)
    }

    growth_fit = list(k = NA_real_, CV1 = NA_real_, CV2 = NA_real_,
                      rmse = NA_real_, max_abs_error = NA_real_, sd_rmse = NA_real_,
                      source = "ss3_cessation", laa_ss3 = NULL, laa_opal = mu,
                      sd_ss3 = NULL, sd_opal = sd_a,
                      L0 = L0, Linf = Linf, rmax = rmax, k_gc = k_gc, t50 = t50,
                      age_offset = age_offset, cv_pattern = cv_pattern,
                      sd_model = sd_model, v1 = v1, v2 = v2, v2_at_oldest = v2_at)
    if (!is.null(ss3$rep)) {
        eg = tryCatch(ss3_endgrowth_table(ss3$rep, dims, basis = "Beg"), error = function(e) NULL)
        if (!is.null(eg)) {
            growth_fit$laa_ss3       = eg$length
            growth_fit$sd_ss3        = eg$sd
            growth_fit$rmse          = sqrt(mean((mu - eg$length)^2))
            growth_fit$max_abs_error = max(abs(mu - eg$length))
            growth_fit$max_rel_error = max(abs(mu / eg$length - 1))
            growth_fit$sd_rmse       = sqrt(mean((sd_a - eg$sd)^2))
            growth_fit$sd_max_rel_error = max(abs(sd_a / eg$sd - 1))
            # Below the plus group, where opal's curve and SS3's agree exactly
            below = seq_len(n_age - 1L)
            growth_fit$max_rel_error_below_plus    = max(abs(mu[below] / eg$length[below] - 1))
            growth_fit$sd_max_rel_error_below_plus = max(abs(sd_a[below] / eg$sd[below] - 1))
        }
    }

    bio$A1 = dims$age_a[A1_pos]
    bio$A2 = dims$age_a[A2_pos]
    bio$L1 = mu[A1_pos]
    bio$L2 = mu[A2_pos]
    bio$growth_model = "cessation"
    bio$growth_age_offset = age_offset
    bio$L0 = L0
    bio$rmax = rmax
    bio$k_gc = k_gc
    bio$t50 = t50
    bio$growth_sd_model = sd_model
    if (identical(sd_model, "sd")) {
        bio$SD1 = v1
        bio$SD2 = v2_at
    } else {
        bio$CV1 = v1
        bio$CV2 = v2_at
    }
    bio$growth_fit = growth_fit
    bio
}
