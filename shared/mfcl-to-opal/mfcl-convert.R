# MULTIFAN-CL (MFCL) to opal converter for single-region MFCL models with
# annual age classes and sub-annual fishing incidents (e.g. the south-west
# Pacific albacore bridge model, case-studies/swpo-alb).
#
# opal is single region, with n_season seasons per year, catch removed at the
# start of each season (Pope), selectivity by year (not season) and annual
# length compositions. The converter maps the MFCL model onto that structure
# and records each approximation in `meta$approximations`:
#   * seasons: MFCL's model year starts in calendar month month_1 (field 9
#     of the .frq dimension line; October for albacore), where recruitment
#     and ageing happen. Fishing incidents are assigned to model years and
#     months as MFCL does (calendar months >= month_1 belong to the next model
#     year) and to opal seasons by model month (season = ceiling(model month /
#     (12 / n_season))); incidents in the same fishery and season are summed;
#   * catch: MFCL's catch-conditioned Baranov removals at each incident's time
#     become opal start-of-season Pope removals;
#   * selectivity: MFCL's seasonal selectivity-at-age (fish flag 74) is
#     averaged over seasons, weighted by each season's catch (index fisheries:
#     their survey seasons), for opal's sel_fa_external;
#   * length compositions: MFCL's per-incident samples (robust normal
#     likelihood, sample sizes capped then divided by fish flag 49) become one
#     composition per fishery and year (opal's multinomial), each sample
#     weighted by its capped sample size;
#   * CPUE: MFCL's survey index likelihood (log-normal, variance lambda_i x
#     sigma^2, with lambda_i the .frq penalty column over its mean and sigma
#     fish flag 92 / 100) becomes opal's log-normal with se_i =
#     sigma x sqrt(lambda_i).
#
# Biology (mean weight, maturity and natural mortality at age) is taken from
# the MFCL plot report at age, so opal's spawning biomass is MFCL's adult
# biomass (numbers x weight x maturity at the start of the year).

#' Read and convert an MFCL model to opal inputs.
#'
#' @param dir folder with the MFCL files.
#' @param frq,par,rep file names (in `dir`) of the .frq, fitted .par and its
#'   plot report.
#' @param index_fisheries MFCL fisheries that are CPUE indices (no catch is
#'   removed for them).
#' @param n_season opal seasons per year.
#' @param lf_cap MFCL's cap on length sample sizes before the divisor.
#' @param sigma_r opal recruitment SD (MFCL's penalty on recruitment
#'   deviations from the stock-recruitment relationship).
#' @return list(data, parameters, lf, selectivity, mfcl, meta).
mfcl_to_opal = function(dir, frq, par, rep, index_fisheries,
                        n_season = 4L, lf_cap = 1000, sigma_r = 0.6) {
    frq_path = file.path(dir, frq)
    p = mfcl_read_par(file.path(dir, par))
    r = mfcl_read_rep(file.path(dir, rep))
    fr = parse_frq(frq_path)
    ce = data.table::as.data.table(cateffpen(fr))
    ff = mfcl_fish_flags(p)
    pf = mfcl_flags(p, "# The parest_flags")

    n_reg = mfcl_rep(r, "# Number of regions")[1, 1]
    if (n_reg != 1) stop("mfcl_to_opal() needs a single-region MFCL model, not ", n_reg, " regions")
    if (mfcl_rep(r, "# Number of recruitments per year")[1, 1] != 1)
        stop("mfcl_to_opal() needs annual recruitment (age classes in years)")

    first_yr = mfcl_rep(r, "# Year 1")[1, 1]
    n_year = mfcl_rep(r, "# Number of time periods")[1, 1]
    years = first_yr + seq_len(n_year) - 1L
    n_age = mfcl_rep(r, "# Number of age classes")[1, 1]
    n_fish = mfcl_rep(r, "# Number of fisheries")[1, 1]
    lfr = lf_range(fr)
    month_1 = mfcl_frq_month1(frq_path)
    month_season = function(mm) as.integer(ceiling(mm / (12 / n_season)))

    # Catch (index fisheries remove nothing)
    data_flags = mfcl_frq_data_flags(frq_path)
    catch_units_f = ifelse(data_flags[1, ] == 1, 1L, 2L)        # opal: 1 weight, 2 numbers
    ce[, `:=`(model_year = year + (month_1 > 1 & month >= month_1),
              model_month = (month - month_1) %% 12 + 1)]
    ce[, `:=`(y = model_year - first_yr + 1L, s = month_season(model_month))]
    if (any(ce$y < 1 | ce$y > n_year)) stop("Fishing incidents outside the model years")
    catch = array(0, c(n_year, n_season, n_fish))
    ext = ce[!(fishery %in% index_fisheries) & catch > 0, .(catch = sum(catch)), by = .(y, s, fishery)]
    catch[cbind(ext$y, ext$s, ext$fishery)] = ext$catch

    # CPUE
    sigma_f = ff[, 92] / 100
    cp = ce[fishery %in% index_fisheries & effort > 0]
    cp[, lambda := penalty / mean(penalty), by = fishery]
    cpue = cp[, .(year = model_year, month = model_month, ts = (y - 1L) * n_season + s, fishery,
                  index = match(fishery, index_fisheries), metric = "cpue",
                  units = ifelse(data_flags[1, fishery] == 1, 1L, 2L),
                  value = catch / effort, se = sigma_f[fishery] * sqrt(lambda))]
    data.table::setorder(cpue, index, ts)
    if (anyDuplicated(cpue[, .(index, ts)])) stop("Two CPUE observations share an index and season")

    # Biology at age (plot report)
    mu_a = mfcl_rep(r, "# Mean lengths at age")[1, ]
    sd_a = mfcl_rep(r, "# SD of length at age")[1, ]
    w_a = mfcl_rep(r, "# Mean weights at age")[1, ] / 1000     # tonnes
    m_a = mfcl_rep(r, "# Natural mortality at age")[1, ]
    mat_a = mfcl_num(p[[mfcl_key(p, "# percent maturity")]]$body)
    lw = mfcl_num(mfcl_read_lines(sub("\\.frq$", ".ini", frq_path))[
        which(startsWith(mfcl_read_lines(sub("\\.frq$", ".ini", frq_path)), "# Length-weight")) + 1L])
    vb = mfcl_matrix(p[[mfcl_key(p, "# The von Bertalanffy parameters")]]$body)

    # Selectivity at age by fishery and season, and its catch-weighted average
    sel_rows = mfcl_rep(r, "# Selectivity by age class")
    n_sel_season = pmax(1, ff[, 74])
    first_row = cumsum(c(1, n_sel_season))[seq_len(n_fish)]
    sel_fsa = array(NA_real_, c(n_fish, n_season, n_age))
    for (f in seq_len(n_fish)) for (s in seq_len(n_season)) {
        k = if (n_sel_season[f] == 1) 1 else as.integer(ceiling(s / (n_season / n_sel_season[f])))
        sel_fsa[f, s, ] = sel_rows[first_row[f] + k - 1L, ]   # by model-month season
    }
    w_fs = matrix(0, n_fish, n_season)
    for (f in seq_len(n_fish)) {
        if (f %in% index_fisheries) {
            ss = unique(cpue[fishery == f, (ts - 1L) %% n_season + 1L])
            w_fs[f, ss] = 1
        } else {
            w_fs[f, ] = colSums(catch[, , f])
        }
    }
    if (any(rowSums(w_fs) == 0)) w_fs[rowSums(w_fs) == 0, ] = 1
    sel_fa = t(vapply(seq_len(n_fish), function(f) colSums(w_fs[f, ] / sum(w_fs[f, ]) * sel_fsa[f, , ]), numeric(n_age)))

    # Length compositions: per fishery-year, samples of >= the MFCL minimum
    # (parest flag 312), each rescaled to its capped sample size, then summed
    ln = data.table::as.data.table(lnfrq(fr))
    bins = setdiff(names(ln), c("year", "month", "week", "fishery"))
    ln[, n := rowSums(.SD), .SDcols = bins]
    lf_min_n = if (pf[311] > 0) pf[312] else 0
    ln = ln[n >= max(lf_min_n, 1e-12)]
    ln[, scale := pmin(n, lf_cap) / n]
    for (b in bins) data.table::set(ln, j = b, value = ln[[b]] * ln$scale)
    ln[, y := year + (month_1 > 1 & month >= month_1) - first_yr + 1L]
    ln[, year := y + first_yr - 1L]
    lf_wide = ln[, lapply(.SD, sum), by = .(fishery, year, y), .SDcols = bins]
    lf_wide[, `:=`(month = 1L, ts = y)]
    lf_wide[, y := NULL]
    data.table::setcolorder(lf_wide, c("fishery", "year", "month", "ts", bins))
    data.table::setorder(lf_wide, fishery, ts)
    lf_var_adjust = ff[, 49]
    lf_var_adjust[lf_var_adjust == 0] = 1

    data = list(
        age_a = seq_len(n_age), n_age = n_age, min_age = 1L, max_age = n_age,
        n_season = n_season, n_fishery = n_fish, n_year = n_year,
        first_yr = first_yr, last_yr = max(years), years = years, first_yr_catch = first_yr,
        len_bin_start = lfr[["LFFirst"]], len_bin_width = lfr[["LFWidth"]], n_len = lfr[["LFIntervals"]],
        catch_units_f = catch_units_f, catch_obs_ysf = catch,
        cpue_switch = 1L, n_index = length(index_fisheries),
        cpue_data = as.data.frame(cpue),
        lw_a = lw[1] / 1000, lw_b = lw[2],
        weight = w_a, maturity = mat_a, fecundity = w_a, sex_ratio = rep(1, n_age),
        spawning_potential = mat_a * w_a, M = m_a,
        A1 = 1L, A2 = as.integer(n_age),
        sel_type_f = rep(2L, n_fish),
        removal_switch_f = rep(0L, n_fish),
        bias_adj_y = rep(0, n_year), init_bias_adj_a = rep(0, n_age)
    )

    parameters = list(
        log_B0 = log(1), log_h = log(0.8), log_sigma_r = log(sigma_r),
        log_cpue_q = rep(0, length(index_fisheries)), cpue_creep = rep(0, length(index_fisheries)),
        log_cpue_tau = rep(-Inf, length(index_fisheries)), log_cpue_omega = rep(0, length(index_fisheries)),
        log_lf_tau = rep(log(0.1), n_fish),
        rdev_y = rep(0, n_year), init_rdev_a = rep(0, n_age),
        log_init_F_f = rep(log(1e-8), n_fish),
        par_sel = matrix(0, n_fish, 6),
        log_L1 = log(mu_a[1]), log_L2 = log(mu_a[n_age]), log_k = log(vb[3, 1]),
        log_CV1 = log(sd_a[1] / mu_a[1]), log_CV2 = log(sd_a[n_age] / mu_a[n_age])
    )

    mfcl = list(
        N_ya = mfcl_rep(r, "# Population Number by age", region = 1),
        recruitment_y = mfcl_rep(r, "# Recruitment")[, 1],
        adult_biomass_y = mfcl_rep(r, "# Adult biomass")[, 1],
        total_biomass_y = mfcl_rep(r, "# Total biomass")[, 1],
        adult_biomass0_y = tryCatch(mfcl_rep(r, "# Adult biomass in absence of fishing")[, 1],
                                    error = function(e) NULL),
        mu_a = mu_a, sd_a = sd_a, sel_fsa = sel_fsa,
        realization_time = mfcl_rep_rows(r, "# Time of each realization by fishery"),
        obs_catch = mfcl_rep_rows(r, "# Observed catch by fishery"),
        pred_catch = mfcl_rep_rows(r, "# Predicted catch by fishery"),
        obs_cpue = mfcl_rep_rows(r, "# Observed CPUE by fishery"),
        pred_cpue = mfcl_rep_rows(r, "# Predicted CPUE by fishery"),
        bh = mfcl_bh(r), incidents = ce, objective = mfcl_num(p[[mfcl_key(p, "# Objective function value")]]$body)
    )

    meta = list(
        dir = dir, files = c(frq = frq, par = par, rep = rep),
        md5 = tools::md5sum(file.path(dir, c(frq, par, rep))),
        index_fisheries = index_fisheries, month_1 = month_1, lf_var_adjust = lf_var_adjust,
        lf_min_n = lf_min_n, lf_cap = lf_cap, sel_weights_fs = w_fs,
        approximations = c("seasons by .frq month", "start-of-season Pope removals",
                           "catch-weighted season-average selectivity",
                           "annual length compositions, multinomial",
                           "CPUE at start of season")
    )
    list(data = data, parameters = parameters, lf = lf_wide,
         selectivity = list(sel_fa_external = sel_fa), mfcl = mfcl, meta = meta)
}

#' First calendar month of MFCL's model year (field 9 of the .frq dimension
#' line, MFCL's dataswitch(9); 0 or absent means January).
mfcl_frq_month1 = function(frq_path) {
    x = mfcl_read_lines(frq_path)
    v = mfcl_num(x[which(!startsWith(x, "#"))[1]])
    m = if (length(v) >= 9) v[9] else 0
    if (m < 1 || m > 12) 1L else as.integer(m)
}

#' The five data-flag rows of an .frq file (row 1: 0 catch in numbers, 1 weight).
mfcl_frq_data_flags = function(frq_path) {
    x = mfcl_read_lines(frq_path)
    i = which(startsWith(x, "# Data flags")) + 1L
    do.call(rbind, lapply(x[i:(i + 4L)], mfcl_num))
}

#' MFCL's Beverton-Holt parameters from the plot report header line.
mfcl_bh = function(rep) {
    h = names(rep)[startsWith(names(rep), "# alpha =")]
    if (!length(h)) return(NULL)
    v = regmatches(h, gregexpr("[0-9.]+e[+-][0-9]+", h))[[1]]
    stats::setNames(as.numeric(v[1:4]), c("alpha", "beta", "steepness", "variance"))
}

#' Set opal's B0, initial deviations and recruitment deviations so that opal
#' starts from MFCL's initial numbers-at-age and recruits MFCL's recruitment
#' every year (with bias adjustment off, deviations multiply the BH curve).
#'
#' B0 is MFCL's unfished spawning biomass implied by its Beverton-Holt curve:
#' R0 = alpha - beta / SPR0 and B0 = R0 x SPR0, with opal's annual SPR0. The
#' recruitment deviations depend on opal's spawning biomass, so they are found
#' by fixed-point iteration on opal's own report.
#'
#' @param report_fn function(parameters) returning an opal report.
mfcl_match_recruitment = function(converted, report_fn, max_iter = 50L, tol = 1e-12) {
    d = converted$data
    p = converted$parameters
    mf = converted$mfcl
    n = d$n_age
    S = exp(-d$M)
    rel = cumprod(c(1, S[-n]))
    rel[n] = rel[n] / (1 - S[n])
    spr0 = sum(d$spawning_potential * rel)
    R0 = mf$bh[["alpha"]] - mf$bh[["beta"]] / spr0
    p$log_B0 = log(R0 * spr0)
    p$log_h = log(mf$bh[["steepness"]])
    # opal's unfished equilibrium numbers are R0 x rel; MFCL's first-year numbers
    p$init_rdev_a = log(mf$N_ya[1, ] / (R0 * rel))
    B0 = R0 * spr0
    h = mf$bh[["steepness"]]
    alpha = 4 * h * R0 / (5 * h - 1)
    beta = B0 * (1 - h) / (5 * h - 1)
    R_target = mf$recruitment_y[-1]            # recruits at the start of years 2..n
    for (it in seq_len(max_iter)) {
        rp = report_fn(p)
        sb = rp$spawning_biomass_y[2:d$n_year]
        new = log(R_target / (alpha * sb / (beta + sb)))
        delta = max(abs(new - p$rdev_y[seq_len(d$n_year - 1L)]))
        p$rdev_y[seq_len(d$n_year - 1L)] = new
        if (delta < tol) break
    }
    attr(p, "iterations") = it
    p
}
