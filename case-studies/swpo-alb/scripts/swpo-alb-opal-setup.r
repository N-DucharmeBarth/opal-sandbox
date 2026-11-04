# Nicholas Ducharme-Barth
# 2026/11/03
# Shared set-up for the opal scripts of the south-west Pacific albacore case
# study (sourced by 03-convert-fit-swpo-alb-opal.r and
# 04-length-comp-diagnostics.r): converts MFCL bridge model M3 to opal, fits
# opal's growth form to MFCL's, builds O0 (opal at MFCL's values) and runs
# the catch-timing and seasonal-selectivity replays. Defines the objects the
# two scripts use (cv, d0, p0, p_base, rep0, replays, fixed_map(),
# report_at(), ...).
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Packages and paths
    library(data.table)
    library(magrittr)
    library(opal)
    library(RTMB)

    proj_dir  = this.path::this.proj()
    dir_cs    = file.path(proj_dir, "case-studies", "swpo-alb")
    dir_m3    = file.path(dir_cs, "model-files", "mfcl", "03-fixed-growth")
    dir_cache = file.path(dir_cs, "vignettes", "data")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
    source(file.path(proj_dir, "shared", "extractors", "mfcl", "parse-frq.r"))
    for (f in list.files(file.path(proj_dir, "shared", "mfcl-to-opal"), "\\.R$", full.names = TRUE)) source(f)

    refit = identical(tolower(Sys.getenv("SWPO_ALB_REFIT")), "true")
    # MFCL's recruitment penalty is equivalent to a CV of 0.71 (assessment
    # report, section 5.2.2): sigma_R = sqrt(log(1 + 0.71^2))
    sigma_r = sqrt(log(1 + 0.71^2))

#_____________________________________________________________________________
# Convert
    cv = mfcl_to_opal(dir_m3, "alb.frq", "03.par", "plot-03.par.rep",
                      index_fisheries = 14:15, sigma_r = sigma_r)
    d = cv$data
    mf = cv$mfcl
    n_age = d$n_age
    n_year = d$n_year
    fishery_names = c("1 LL DWFN 1a", "2 LL DWFN 1b", "3 LL DWFN 1c", "4 LL DWFN 1d",
                      "5 LL PICT 1ab", "6 LL PICT 1ab", "7 LL AU/NZ 1ab", "8 LL DWFN 1cd",
                      "9 LL PICT 1cd", "10 LL AU/NZ 1cd", "11 TR 1c", "12 TR 1d",
                      "13 DN 1cd", "14 LL index 1abcd", "15 TR index 1ef")

    # Length compositions: one per fishery and model year, multinomial, with
    # MFCL's divisors (fish flag 49) on the capped sample sizes
    d = prep_lf_data(d, as.data.frame(cv$lf), lf_switch = 1L,
                     lf_var_adjust = cv$meta$lf_var_adjust)
    d$lambdas = list(prior = 0)   # MFCL has no parameter priors here

#_____________________________________________________________________________
# Growth in opal's form (only the length compositions use it: weight,
# maturity and M are supplied at age)
    # VB (Schnute) with age-specific K multipliers at ages 1-4 for MFCL's
    # independent mean lengths at ages 2-4; SD linear in mean length between
    # MFCL's SDs at ages 1 and 12
    k_mult_ages = 1:4
    growth_obj = function(x) {
        mu = get_growth(n_age, A1 = 1, A2 = n_age, L1 = exp(x[1]), L2 = exp(x[2]), log_k = x[3],
                        log_k_mult = x[4:7], k_mult_ages = k_mult_ages)
        sum((mu - mf$mu_a)^2)
    }
    g = optim(c(log(mf$mu_a[1]), log(mf$mu_a[n_age]), cv$parameters$log_k, rep(0, 4)),
              growth_obj, method = "BFGS", control = list(maxit = 5000, reltol = 1e-14))
    mu_opal = get_growth(n_age, A1 = 1, A2 = n_age, L1 = exp(g$par[1]), L2 = exp(g$par[2]),
                         log_k = g$par[3], log_k_mult = g$par[4:7], k_mult_ages = k_mult_ages)
    sd_opal = get_sd_at_age(mu_opal, mu_opal[1], mu_opal[n_age], sd_model = "sd",
                            log_SD1 = log(mf$sd_a[1]), log_SD2 = log(mf$sd_a[n_age]),
                            A1 = 1, A2 = n_age)
    growth_check = data.table(age = seq_len(n_age), mfcl_length = mf$mu_a, opal_length = mu_opal,
                              mfcl_sd = mf$sd_a, opal_sd = as.numeric(sd_opal))
    print(growth_check[, .(rmse_length = sqrt(mean((opal_length - mfcl_length)^2)),
                           max_sd_diff = max(abs(opal_sd - mfcl_sd)))])
    d$k_mult_ages = k_mult_ages
    d$growth_sd_model = "sd"

    p_base = cv$parameters
    p_base$log_L1 = g$par[1]
    p_base$log_L2 = g$par[2]
    p_base$log_k = g$par[3]
    p_base$log_k_mult = g$par[4:7]
    p_base$log_SD1 = log(mf$sd_a[1])
    p_base$log_SD2 = log(mf$sd_a[n_age])

#_____________________________________________________________________________
# O0: opal at MFCL's values
    d0 = d
    d0$sel_fa_external = cv$selectivity$sel_fa_external
    fixed_map = function(p) lapply(p, function(z) factor(rep(NA, length(z))))
    report_at = function(p, data = d0) {
        MakeADFun(cmb(opal_model, data), p, map = fixed_map(p), silent = TRUE)$report()
    }
    p0 = mfcl_match_recruitment(list(data = d0, parameters = p_base, mfcl = mf), report_at)
    rep0 = report_at(p0)
    o0_check = data.table(
        quantity = c("numbers at age, first year", "recruitment"),
        max_rel_diff = c(max(abs(rep0$number_ysa[1, 1, ] / mf$N_ya[1, ] - 1)),
                         max(abs(rep0$number_ysa[2:n_year, 1, 1] / mf$recruitment_y[-1] - 1))))
    print(o0_check)

    # The index catchabilities that match MFCL's (on opal's normalised scale,
    # q is where opal's prediction centres on the observations)
    cpue = as.data.table(d0$cpue_data)
    cpue[, pred := rep0$cpue_pred]
    q0 = cpue[, .(log_q = mean(log(value)) - mean(log(pred))), by = index]
    p0$log_cpue_q = q0$log_q
    rep0 = report_at(p0)

#_____________________________________________________________________________
# Replays (standalone simulations, not opal code)
    inc = copy(mf$incidents)[!(fishery %in% cv$meta$index_fisheries) & catch > 0,
                             .(year, month, fishery, catch)]
    m1 = cv$meta$month_1
    inc[, `:=`(model_year = year + (m1 > 1 & month >= m1), model_month = (month - m1) %% 12 + 1)]
    inc[, `:=`(t = (model_month - 0.5) / 12, season = as.integer(ceiling(model_month / 3)))]
    sp = d0$spawning_potential
    replay = function(timing = c("opal", "mfcl"), selectivity = c("average", "seasonal")) {
        timing = match.arg(timing)
        selectivity = match.arg(selectivity)
        M = d0$M
        w = d0$weight
        sel_of = function(f, s) if (selectivity == "seasonal") mf$sel_fsa[f, s, ] else d0$sel_fa_external[f, ]
        N = mf$N_ya[1, ]
        sb = numeric(n_year)
        for (iy in seq_len(n_year)) {
            sb[iy] = sum(N * sp)
            yi = inc[model_year == d0$years[iy]]
            if (timing == "opal") {
                for (ss in seq_len(d0$n_season)) {
                    ys = yi[season == ss, .(catch = sum(catch)), by = fishery]
                    hr = numeric(n_age)
                    for (k in seq_len(nrow(ys))) {
                        f = ys$fishery[k]
                        cw = if (d0$catch_units_f[f] == 1) w else 1
                        hr = hr + ys$catch[k] / (sum(N * sel_of(f, ss) * cw) + 1e-6) * sel_of(f, ss)
                    }
                    N = N * (1 - hr) * exp(-M / d0$n_season)
                }
            } else {
                # MFCL splits each year's M equally among the year's fishing
                # periods (rnaux2.cpp: fraction = -log(n periods)), with F and
                # M acting together in each period (Baranov)
                tt = sort(unique(yi$t))
                Mp = M / length(tt)
                for (j in seq_along(tt)) {
                    ii = yi[t == tt[j]]
                    K = nrow(ii)
                    sel = matrix(0, K, n_age)
                    cw = matrix(1, K, n_age)
                    for (k in seq_len(K)) {
                        sel[k, ] = sel_of(ii$fishery[k], ii$season[k])
                        if (d0$catch_units_f[ii$fishery[k]] == 1) cw[k, ] = w
                    }
                    F = ii$catch / rowSums(sweep(sel * cw, 2, N, "*"))
                    for (it in 1:30) {
                        Z = colSums(F * sel) + Mp
                        F = F * ii$catch / rowSums(sweep(F * sel * cw, 2, N / Z * (1 - exp(-Z)), "*"))
                    }
                    N = N * exp(-colSums(F * sel) - Mp)
                }
            }
            if (iy < n_year) {
                N_next = c(mf$recruitment_y[iy + 1], N[-n_age])
                N_next[n_age] = N_next[n_age] + N[n_age]
                N = N_next
            }
        }
        sb
    }
    replays = rbindlist(lapply(list(c("opal", "average"), c("opal", "seasonal"),
                                    c("mfcl", "average"), c("mfcl", "seasonal")), function(a)
        data.table(timing = a[1], selectivity = a[2], year = d0$years,
                   sb = replay(a[1], a[2]), mfcl = mf$adult_biomass_y,
                   opal_o0 = rep0$spawning_biomass_y[seq_len(n_year)])))
    replays[, rel_mfcl := sb / mfcl - 1]
    print(replays[, .(rel_end = rel_mfcl[.N], rel_max = max(abs(rel_mfcl)),
                      vs_o0_max = max(abs(sb / opal_o0 - 1))), by = .(timing, selectivity)])


#_____________________________________________________________________________
# MFCL's length-at-age through the year (M3 length.fit): MFCL predicts each
# length sample at its date, from lengths that grow through the year; opal's
# growth fitted to the mid-year lengths (g_mid) is used by O1-mid and the
# diagnostics
    lfit = mfcl_read_length_fit(file.path(dir_m3, "length.fit"))
    laa_month = mfcl_length_at_age_by_month(lfit)
    mid_year = (laa_month[6, ] + laa_month[7, ]) / 2
    g_mid = optim(g$par, function(x) sum((get_growth(n_age, A1 = 1, A2 = n_age, L1 = exp(x[1]),
                  L2 = exp(x[2]), log_k = x[3], log_k_mult = x[4:7], k_mult_ages = k_mult_ages) - mid_year)^2),
                  method = "BFGS", control = list(maxit = 5000, reltol = 1e-14))
    with_growth = function(p, x) {
        p$log_L1 = x[1]; p$log_L2 = x[2]; p$log_k = x[3]; p$log_k_mult = x[4:7]
        p
    }
