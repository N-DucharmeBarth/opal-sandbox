# Nicholas Ducharme-Barth
# 2026/11/04
# Why do opal's double-normal selectivities collapse to knife-edges (O2)?
# Compares O2 (opal's selectivity) with O1i (MFCL's selectivity; both with
# mid-year growth and the initial population estimated; 05-fit-opal-
# selectivity.r) at their estimates, fishery by fishery:
#
#   * the length-composition negative log-likelihood by fishery;
#   * observed and predicted mean length of each composition (fishery-year);
#   * selectivity-at-age: MFCL's (season-averaged) and O2's;
#   * O2's selectivity-at-length curves.
#
# OUTPUT: case-studies/swpo-alb/vignettes/data/swpo-alb-sel-diagnostics.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Set-up and the two fits
    source(file.path(this.path::this.proj(), "case-studies", "swpo-alb", "scripts", "swpo-alb-opal-setup.r"))
    fits = readRDS(file.path(dir_cache, "swpo-alb-opal-sel-fits.rds"))$fits
    stopifnot("Run 05-fit-opal-selectivity.r first" = all(c("O1i", "O2") %in% names(fits)))

    lam = d0$lambdas
    lam$init_rec = 0
    d1 = d0
    d1$lambdas = lam
    d2 = d0
    d2$lambdas = lam
    d2$sel_fa_external = NULL
    d2$sel_type_f = ifelse(seq_len(d0$n_fishery) == 14L, 1L, 2L)
    d2$sel_basis_f = rep(1L, d0$n_fishery)
    fixed_all = function(p) lapply(p, function(z) factor(rep(NA, length(z))))
    report_of = function(data, p) MakeADFun(cmb(opal_model, data), p, map = fixed_all(p), silent = TRUE)$report()
    r1 = report_of(d1, fits$O1i$parameters)
    r2 = report_of(d2, fits$O2$parameters)
    stopifnot(isTRUE(all.equal(sum(r1$lp_lf), fits$O1i$lp[["lf"]])),
              isTRUE(all.equal(sum(r2$lp_lf), fits$O2$lp[["lf"]])))

#_____________________________________________________________________________
# Compositions: observed and predicted, sample by sample
    # flat (fishery-grouped) order of opal's length likelihood
    lf_fish = rep(d0$lf_fishery_f, d0$lf_n_f)
    lf_year = unlist(d0$lf_year_fi)
    lf_n = unlist(d0$lf_n_fi)
    nb = d0$n_len
    obs = matrix(d0$lf_obs_prop, ncol = nb, byrow = TRUE)
    pred1 = do.call(rbind, r1$lf_pred_row)
    pred2 = do.call(rbind, r2$lf_pred_row)
    mid = d0$len_mid
    samples = data.table(fishery = lf_fish, year = d0$years[lf_year], n = lf_n,
                         obs_mean = c(obs %*% mid), O1i_mean = c(pred1 %*% mid), O2_mean = c(pred2 %*% mid),
                         O1i_nll = r1$lp_lf, O2_nll = r2$lp_lf)
    samples[, name := fishery_names[fishery]]
    by_fishery = samples[, .(n_obs = .N, O1i_nll = sum(O1i_nll), O2_nll = sum(O2_nll),
                             obs_mean = weighted.mean(obs_mean, n), O1i_mean = weighted.mean(O1i_mean, n),
                             O2_mean = weighted.mean(O2_mean, n),
                             O1i_rmse = sqrt(weighted.mean((O1i_mean - obs_mean)^2, n)),
                             O2_rmse = sqrt(weighted.mean((O2_mean - obs_mean)^2, n))),
                         by = .(fishery, name)]
    by_fishery[, gain := O1i_nll - O2_nll]
    setorder(by_fishery, -gain)
    print(by_fishery, digits = 4)
    # composition shapes pooled over years (n-weighted) per fishery
    pooled = rbindlist(lapply(sort(unique(lf_fish)), function(f) {
        i = which(lf_fish == f)
        w = lf_n[i] / sum(lf_n[i])
        data.table(fishery = f, name = fishery_names[f], length = mid,
                   observed = colSums(obs[i, , drop = FALSE] * w),
                   O1i = colSums(pred1[i, , drop = FALSE] * w),
                   O2 = colSums(pred2[i, , drop = FALSE] * w))
    }))

#_____________________________________________________________________________
# Selectivity
    sel_age = rbind(
        data.table(source = "MFCL", fishery = rep(seq_len(d0$n_fishery), each = n_age),
                   age = rep(seq_len(n_age), d0$n_fishery), sel = as.vector(t(cv$selectivity$sel_fa_external))),
        data.table(source = "O2", fishery = rep(seq_len(d0$n_fishery), each = n_age),
                   age = rep(seq_len(n_age), d0$n_fishery), sel = as.vector(t(r2$sel_fya[, 1, ]))))
    sel_age[, name := fishery_names[fishery]]
    ps = fits$O2$parameters$par_sel
    sel_len = rbindlist(lapply(seq_len(d0$n_fishery), function(f) data.table(
        fishery = f, name = fishery_names[f], length = mid,
        sel = sel_length(mid, ps[f, ], d2$sel_type_f[f]))))
    par_sel = data.table(fishery = seq_len(d0$n_fishery), name = fishery_names, type = d2$sel_type_f,
                         peak = ps[, 1], plateau = ps[, 2], asc_width = ps[, 3], desc_width = ps[, 4],
                         init = ps[, 5], final = ps[, 6])

    saveRDS(list(samples = samples, by_fishery = by_fishery, pooled = pooled, sel_age = sel_age,
                 sel_len = sel_len, par_sel = par_sel,
                 totals = c(O1i_lf = sum(r1$lp_lf), O2_lf = sum(r2$lp_lf),
                            O1i_cpue = sum(r1$lp_cpue), O2_cpue = sum(r2$lp_cpue))),
            file.path(dir_cache, "swpo-alb-sel-diagnostics.rds"))
