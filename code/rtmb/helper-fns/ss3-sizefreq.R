#' SS3 generalised size compositions and likelihood lambdas, as reference
#' tables for opal
#'
#' SS3's generalised size composition data (sizefreq) sit on their own bins,
#' which opal can now take as a second composition bin scheme
#' (prep_lf_data() with a list of wide data frames). These helpers read the
#' observations, SS3's expected values and the control file's lambda
#' changes, so the opal package can rebuild an SS3-matched model without
#' SS3 or r4ss.
#'
#' Copyright (c) 2026 Nicholas Ducharme-Barth
#' GPL-3.0

#_____________________________________________________________________________
#' Generalised size compositions of one SS3 method
#'
#' @param ss3 list from \code{ss3_read_model()}.
#' @param dir_ss3 SS3 model directory; \code{dat_file} is read again with
#'   \code{section = 2} for SS3's expected values (data.ss_new only).
#' @param cfg converter configuration list (for the index split, which
#'   assigns each observation to its quarterly opal index fishery).
#' @param method SS3 sizefreq method.
#' @return list: \code{method}, \code{fleet} (SS3), \code{units}
#'   ("numbers" or "biomass"), \code{scale} (SS3 code, 3 = cm),
#'   \code{mincomp}, \code{bins} (lower edges, sign removed),
#'   \code{first_bin_negative}, \code{lambda} (control file, 1 if not
#'   changed) and \code{data}, a long data.frame (\code{ts}, \code{fishery},
#'   \code{nsamp}, \code{bin}, \code{obs}, \code{expected}) with observed and
#'   expected proportions.
ss3_sizefreq_reference = function(ss3, dir_ss3, cfg, dat_file = "data.ss_new",
                                  method = 1L) {
    dat = ss3$dat
    if (is.null(dat$N_sizefreq_methods) || dat$N_sizefreq_methods < method) {
        stop(sprintf("SS3 data have no sizefreq method %d.", method), call. = FALSE)
    }
    exp_dat = suppressWarnings(r4ss::SS_readdat(file.path(dir_ss3, dat_file),
                                                version = ss3$provenance$ss3_version %||% "3.24",
                                                section = 2, verbose = FALSE))
    bins = as.numeric(dat$sizefreq_bins_list[[method]])
    obs  = as.data.frame(dat$sizefreq_data_list[[method]])
    expd = as.data.frame(exp_dat$sizefreq_data_list[[method]])
    stopifnot("SS3 observed and expected sizefreq rows differ" =
                  nrow(obs) == nrow(expd) &&
                  all(obs$Yr == expd$Yr) && all(obs$FltSvy == expd$FltSvy))
    nb = length(bins)
    comp_cols = (ncol(obs) - nb + 1L):ncol(obs)
    fleet = unique(as.integer(obs$FltSvy))
    if (length(fleet) != 1L) stop("Sizefreq method spans several fleets.", call. = FALSE)

    ts = as.integer(obs$Yr) - as.integer(cfg$time$first_ts %||% 1L) + 1L
    fishery = ss3_apply_index_split(rep(fleet, nrow(obs)), ts, cfg)
    long = data.table::rbindlist(lapply(seq_len(nrow(obs)), function(i) {
        o = as.numeric(obs[i, comp_cols])
        data.table::data.table(ts = ts[i], fishery = fishery[i],
                               nsamp = as.numeric(obs$Nsamp[i]),
                               bin = abs(bins),
                               obs = o / sum(o),
                               expected = as.numeric(expd[i, comp_cols]))
    }))

    lam = as.data.frame(ss3$ctl$lambdas)
    hit = lam$like_comp == 6 & lam[["fleet/survey"]] == fleet & lam$sizefreq_method == method
    units = as.integer(dat$units_per_method[method])

    list(method = as.integer(method),
         fleet  = fleet,
         units  = if (units == 2L) "numbers" else "biomass",
         scale  = as.integer(dat$scale_per_method[method]),
         mincomp = as.numeric(dat$mincomp_per_method[method]),
         bins   = abs(bins),
         first_bin_negative = bins[1] < 0,
         lambda = if (any(hit)) as.numeric(lam$value[hit][1]) else 1,
         data   = as.data.frame(long))
}

#_____________________________________________________________________________
#' Lambda changes in the SS3 control file
#'
#' @return data.frame: \code{like_comp} (SS3 code), \code{component},
#'   \code{fleet}, \code{phase}, \code{value}, \code{sizefreq_method}.
ss3_lambda_reference = function(ss3) {
    lam = as.data.frame(ss3$ctl$lambdas)
    codes = c("survey", "discard", "mean body weight", "length", "age",
              "size frequency", "size at age", "catch", "initial equilibrium catch",
              "recruitment deviations", "parameter priors", "parameter deviations",
              "crash penalty", "morph composition", "tag composition",
              "tag negative binomial", "F ballpark")
    data.frame(like_comp = as.integer(lam$like_comp),
               component = codes[as.integer(lam$like_comp)],
               fleet = as.integer(lam[["fleet/survey"]]),
               phase = as.integer(lam$phase),
               value = as.numeric(lam$value),
               sizefreq_method = as.integer(lam$sizefreq_method),
               row.names = NULL)
}
