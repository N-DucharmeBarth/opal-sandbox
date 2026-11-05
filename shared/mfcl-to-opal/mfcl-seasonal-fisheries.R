# Split MFCL fisheries with seasonal selectivity into one opal fishery per
# season (for an annual-age opal model with n_season seasons).
#
# opal's selectivity varies by year, not season. A fishery whose MFCL
# selectivity differs by season becomes n_season fisheries, each with its
# season's catch only, its season's selectivity and its season's length
# compositions (annual within each seasonal fishery). The CPUE index
# fisheries keep one fishery each. Fishery order: the split fisheries (season
# within fishery), then the others in their MFCL order.

#' @param cv output of mfcl_to_opal().
#' @param split MFCL fisheries to split (default: those with more than one
#'   selectivity season, fish flag 74).
#' @return a copy of cv with data, parameters, lf, selectivity and meta
#'   rewritten, and meta$fishery_map (opal fishery, MFCL fishery, season).
mfcl_split_seasonal_fisheries = function(cv, split = NULL) {
    d = cv$data
    ns = d$n_season
    n_f = d$n_fishery
    sel_fsa = cv$mfcl$sel_fsa
    if (is.null(split)) {
        split = which(vapply(seq_len(n_f), function(f) {
            m = sel_fsa[f, , , drop = FALSE]
            any(apply(m, 2, function(v) max(abs(v - m[1, 1, ]))) > 0)
        }, logical(1)))
    }
    fmap = rbind(data.frame(mfcl = rep(split, each = ns), season = rep(seq_len(ns), length(split))),
                 data.frame(mfcl = setdiff(seq_len(n_f), split), season = NA_integer_))
    fmap$opal = seq_len(nrow(fmap))
    n_new = nrow(fmap)

    catch = array(0, c(d$n_year, ns, n_new))
    sel = matrix(0, n_new, d$n_age)
    for (i in seq_len(n_new)) {
        f = fmap$mfcl[i]
        if (is.na(fmap$season[i])) {
            catch[, , i] = d$catch_obs_ysf[, , f]
            sel[i, ] = cv$selectivity$sel_fa_external[f, ]
        } else {
            catch[, fmap$season[i], i] = d$catch_obs_ysf[, fmap$season[i], f]
            sel[i, ] = sel_fsa[f, fmap$season[i], ]
        }
    }
    out = cv
    dn = d
    dn$n_fishery = n_new
    dn$catch_obs_ysf = catch
    dn$catch_units_f = d$catch_units_f[fmap$mfcl]
    dn$sel_type_f = d$sel_type_f[fmap$mfcl]
    dn$removal_switch_f = d$removal_switch_f[fmap$mfcl]
    cp = d$cpue_data
    rest = fmap[is.na(fmap$season), ]
    if (any(cp$fishery %in% split)) stop("A CPUE index fishery has seasonal selectivity")
    cp$fishery = rest$opal[match(cp$fishery, rest$mfcl)]
    dn$cpue_data = cp
    out$data = dn

    # length compositions: seasonal ones for the split fisheries (annual within
    # each), the annual ones for the rest
    lfs = data.table::copy(cv$lf_season)
    lfs[, season := (ts - 1L) %% ns + 1L]
    lfs[, y := (ts - 1L) %/% ns + 1L]
    lf_split = lfs[fishery %in% split]
    lf_split[, fishery := fmap$opal[match(paste(fishery, season), paste(fmap$mfcl, fmap$season))]]
    lf_split[, `:=`(ts = y, month = 1L)]
    lf_split[, c("season", "y") := NULL]
    lf_rest = data.table::copy(cv$lf)[!(fishery %in% split)]
    lf_rest[, fishery := rest$opal[match(fishery, rest$mfcl)]]
    lf_new = rbind(lf_split, lf_rest)
    data.table::setorder(lf_new, fishery, ts)
    out$lf = lf_new

    p = cv$parameters
    p$log_lf_tau = rep(p$log_lf_tau[1], n_new)
    p$log_init_F_f = rep(p$log_init_F_f[1], n_new)
    p$par_sel = matrix(0, n_new, ncol(p$par_sel))
    out$parameters = p
    out$selectivity = list(sel_fa_external = sel)
    out$meta$lf_var_adjust = cv$meta$lf_var_adjust[fmap$mfcl]
    out$meta$fishery_map = fmap
    out
}
