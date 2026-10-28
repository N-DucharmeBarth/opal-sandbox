# Nicholas Ducharme-Barth
# 2026/09/18
# Convert an SS3 model (read by ss3_read_model) into opal data, parameter and
# length-frequency objects.
#
# Design rules:
#   1. Derived quantities come from Report.sso where SS3 reports them
#      (M-at-age, length-at-age, selectivity-at-age). The control file is used
#      only for things Report.sso does not expose. This removes ambiguity about
#      SS3's internal transformations.
#   2. Every choice that is not derivable from the SS3 files is a config entry.
#      Nothing model-specific is hard-coded in these functions.
#   3. Where opal cannot represent an SS3 structure exactly (age-specific K,
#      age-based selectivity), the function returns both the approximation and
#      a diagnostic of how good it is.
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
#' Model dimensions
#'
#' Note on time: this SS3 model runs 256 "years" of 3 months each (nseas = 1,
#' months_per_seas = 3), i.e. the SS3 year IS the quarter. opal therefore takes
#' n_year = number of quarters and n_season = 1. `ts_per_year` records the
#' implied real-time scaling for anyone converting rates later.
#'
#' Note on ages: opal internal ages are 1:n_age. We carry SS3 ages
#' cfg$age$keep_ss3_ages[1]:cfg$age$keep_ss3_ages[2] so that, with the default
#' of keeping ages 1-28, the opal age index equals the SS3 age. This is the
#' self-documenting choice; set keep_ss3_ages to c(0, Nages) and min_age to 0
#' if you want the full SS3 age range instead.
#'
#' @param ss3 output of \code{ss3_read_model}.
#' @param cfg converter configuration list.
#' @return a \code{list} of dimension elements plus bookkeeping.
ss3_dims_to_opal = function(ss3, cfg) {

    dat = ss3$dat

    styr_ss3  = as.integer(ss3_pick(dat, c("styr"),  what = "styr"))
    endyr_ss3 = as.integer(ss3_pick(dat, c("endyr"), what = "endyr"))
    nseas     = as.integer(ss3_pick(dat, c("nseas", "N_seas"), default = 1L))
    mps       = as.integer(ss3_pick(dat, c("months_per_seas", "month_per_seas"), default = 12L))
    nages_ss3 = as.integer(ss3_pick(dat, c("Nages", "Nage", "N_ages"), what = "Nages"))

    first_ts = as.integer(cfg$time$first_ts %||% styr_ss3)
    last_ts  = as.integer(cfg$time$last_ts  %||% endyr_ss3)
    stopifnot("cfg$time window falls outside the SS3 model period" =
                  first_ts >= styr_ss3 && last_ts <= endyr_ss3)

    keep_ages = as.integer(cfg$age$keep_ss3_ages %||% c(1L, nages_ss3))
    stopifnot("cfg$age$keep_ss3_ages must be length 2" = length(keep_ages) == 2L)
    ss3_ages  = keep_ages[1]:keep_ages[2]
    n_age     = length(ss3_ages)
    min_age   = as.integer(cfg$age$min_age %||% keep_ages[1])
    age_a     = seq(min_age, by = 1L, length.out = n_age)

    # Data length bins, straight from the SS3 data file. SS3's last data bin is
    # an open plus bin and is often placed off the regular grid (here 195, 198).
    # opal assumes a constant width, so cfg$lf$regularise_last_bin moves the
    # last lower edge onto the grid; the SS3 -> opal bin map is kept so the
    # composition converter can relabel that column.
    len_bins_ss3 = sort(unique(as.numeric(ss3_pick(dat, c("lbin_vector", "lbin_vector_pop"),
                                                   what = "lbin_vector"))))
    len_bins = len_bins_ss3
    widths   = diff(len_bins)
    n_lb     = length(len_bins)
    if (isTRUE(cfg$lf$regularise_last_bin) && n_lb > 2L &&
        length(unique(round(widths[-(n_lb - 1L)], 8))) == 1L) {
        len_bins[n_lb] = len_bins[n_lb - 1L] + widths[1]
        widths = diff(len_bins)
    }
    stopifnot("SS3 length bins are not evenly spaced; opal assumes a constant width (see cfg$lf$regularise_last_bin)" =
                  length(unique(round(widths, 8))) == 1L)

    fleet_names = as.character(ss3_pick(dat, c("fleetnames", "fleetname"),
                                        default = character()))

    list(
        # opal dimensions
        n_year        = last_ts - first_ts + 1L,
        n_season      = 1L,
        n_age         = n_age,
        min_age       = min_age,
        age_a         = age_a,
        years         = seq_len(last_ts - first_ts + 1L),
        first_yr      = 1L,
        last_yr       = last_ts - first_ts + 1L,
        n_len         = length(len_bins),
        len_bin_start = len_bins[1],
        len_bin_width = widths[1],
        # bookkeeping (not passed to opal_model, used by the converters)
        .ss3_ages     = ss3_ages,
        .first_ts     = first_ts,
        .last_ts      = last_ts,
        .len_bins     = len_bins,
        .len_bins_ss3 = len_bins_ss3,
        .ts_per_year  = 12 / (nseas * mps),
        .fleet_names  = fleet_names,
        .n_ss3_fleet  = length(fleet_names)
    )
}

#_____________________________________________________________________________
#' Fishery definitions
#'
#' opal fisheries = SS3 extraction fleets, plus one opal fishery per opal index.
#' Index splitting (e.g. one SS3 CPUE series carved into four quarterly series)
#' is declared in cfg$index$split and applied consistently here, in the CPUE
#' converter and in the composition converter.
#'
#' @param ss3 output of \code{ss3_read_model}.
#' @param cfg converter configuration list.
#' @param dims output of \code{ss3_dims_to_opal}.
#' @return a \code{list} describing the opal fishery layout.
ss3_fisheries_to_opal = function(ss3, cfg, dims) {

    dat = ss3$dat

    # Extraction fleets become opal fisheries 1..n_extract, in SS3 order. In
    # SS3 3.24 they are fleets 1..Nfleet; in 3.30 fleets of any type can be in
    # any order, so they are the fleets of type 1 (catch) or 2 (bycatch).
    fi = ss3_pick(dat, c("fleetinfo"), default = NULL)
    if (is.data.frame(fi) && "type" %in% names(fi)) {
        extract_fleets = which(as.integer(fi$type) %in% c(1L, 2L))
        units_extract  = as.numeric(fi$units[extract_fleets])
    } else {
        extract_fleets = seq_len(as.integer(ss3_pick(dat, c("Nfleet", "Nfleets"), what = "Nfleet")))
        # Units of catch: SS3 codes 1 = biomass, 2 = numbers, as in opal.
        units_extract = as.numeric(ss3_pick(dat, c("units_of_catch"),
                                            default = rep(1, length(extract_fleets))))
        units_extract = units_extract[extract_fleets]
    }
    n_extract = length(extract_fleets)

    # Index fleets, after any splitting.
    split_spec = cfg$index$split %||% list()
    idx_map = data.table::rbindlist(lapply(split_spec, function(s) {
        data.table::data.table(
            ss3_fleet  = as.integer(s$fleet),
            opal_fleet = as.integer(s$into),
            slot       = seq_along(as.integer(s$into))
        )
    }), use.names = TRUE, fill = TRUE)

    keep_idx = as.integer(cfg$index$keep_ss3_fleets %||% integer())
    if (length(keep_idx)) {
        unsplit = setdiff(keep_idx, unique(idx_map$ss3_fleet))
        if (length(unsplit)) {
            idx_map = data.table::rbindlist(list(idx_map, data.table::data.table(
                ss3_fleet  = unsplit,
                opal_fleet = n_extract + seq_along(unsplit),
                slot       = 1L
            )), use.names = TRUE, fill = TRUE)
        }
    }
    data.table::setorder(idx_map, opal_fleet)

    n_index  = nrow(idx_map)
    n_fishery = n_extract + n_index

    # opal index units, from the SS3 CPUE setup (SS3: 0 = numbers, 1 = biomass).
    cpue_info = ss3_pick(dat, c("CPUEinfo", "CPUE_info"), default = NULL)
    units_index = rep(1, n_index)
    if (!is.null(cpue_info)) {
        ci = as.data.frame(cpue_info)
        fleet_col = intersect(c("Fleet", "fleet", "index"), names(ci))[1]
        unit_col  = intersect(c("Units", "units"), names(ci))[1]
        if (!is.na(fleet_col) && !is.na(unit_col)) {
            ss3_units = stats::setNames(as.numeric(ci[[unit_col]]), as.integer(ci[[fleet_col]]))
            u = ss3_units[as.character(idx_map$ss3_fleet)]
            # SS3 0 = numbers -> opal 2; SS3 1 = biomass -> opal 1.
            units_index = ifelse(is.na(u), 1, ifelse(u == 0, 2, 1))
        }
    }

    list(
        n_extract     = n_extract,
        extract_fleets = extract_fleets,
        n_index       = n_index,
        n_fishery     = n_fishery,
        catch_units_f = unname(c(units_extract, units_index)),
        index_map     = idx_map
    )
}

#_____________________________________________________________________________
#' opal fishery for each observation of an SS3 fleet
#'
#' Extraction fleets map to their position among the extraction fleets; index
#' fleets to their opal fishery, through any split. NA for fleets opal does
#' not carry.
#'
#' @param ss3_fleet \code{integer} vector of SS3 fleet ids.
#' @param ts \code{integer} vector of opal timesteps (1-based).
#' @param cfg converter configuration list.
#' @param fish output of \code{ss3_fisheries_to_opal}.
#' @return \code{integer} vector of opal fishery ids.
ss3_opal_fishery = function(ss3_fleet, ts, cfg, fish) {
    ss3_fleet = as.integer(ss3_fleet)
    out = match(ss3_fleet, fish$extract_fleets)
    split_fleets = vapply(cfg$index$split %||% list(), function(s) as.integer(s$fleet), integer(1L))
    is_split = ss3_fleet %in% split_fleets
    out[is_split] = ss3_apply_index_split(ss3_fleet[is_split], ts[is_split], cfg)
    im = fish$index_map
    unsplit = im[!im$ss3_fleet %in% split_fleets, ]
    hit = match(ss3_fleet, unsplit$ss3_fleet)
    out[is.na(out) & !is.na(hit)] = unsplit$opal_fleet[hit[is.na(out) & !is.na(hit)]]
    out
}

#_____________________________________________________________________________
#' Assign an opal fishery to each observation of a split SS3 index fleet
#'
#' @param ss3_fleet \code{integer} vector of SS3 fleet ids.
#' @param ts \code{integer} vector of opal timesteps (1-based).
#' @param cfg converter configuration list.
#' @return \code{integer} vector of opal fishery ids.
ss3_apply_index_split = function(ss3_fleet, ts, cfg) {
    out = as.integer(ss3_fleet)
    for (s in (cfg$index$split %||% list())) {
        target = which(out == as.integer(s$fleet))
        if (!length(target)) next
        into   = as.integer(s$into)
        method = s$by %||% "ts_modulo"
        if (!identical(method, "ts_modulo")) {
            stop(sprintf("Unsupported index split method '%s'.", method), call. = FALSE)
        }
        offset = as.integer(s$offset %||% 0L)
        m      = as.integer(s$modulo %||% length(into))
        out[target] = into[((ts[target] - offset) %% m) + 1L]
    }
    out
}

#_____________________________________________________________________________
#' Catch
#'
#' @return a \code{[n_year, n_season, n_fishery]} array.
ss3_catch_to_opal = function(ss3, cfg, dims, fish) {

    catch = ss3_pick(ss3$dat, c("catch"), what = "catch")
    catch = data.table::as.data.table(catch)

    if (all(c("year", "seas", "fleet", "catch") %in% tolower(names(catch)))) {
        # SS 3.30 long format.
        data.table::setnames(catch, tolower(names(catch)))
        long = catch[, .(year = as.integer(year), fleet = as.integer(fleet),
                         value = as.numeric(catch))]
    } else {
        # SS 3.24 wide format: one column per fleet, then year and seas.
        yr_col = intersect(c("year", "Yr", "yr"), names(catch))[1]
        se_col = intersect(c("seas", "Seas", "season"), names(catch))[1]
        fleet_cols = setdiff(names(catch), c(yr_col, se_col))
        stopifnot("Catch fleet columns do not match Nfleet" =
                      length(fleet_cols) == fish$n_extract)
        long = data.table::melt(
            catch[, c(yr_col, fleet_cols), with = FALSE],
            id.vars = yr_col, variable.name = "fleet_name", value.name = "value"
        )
        data.table::setnames(long, yr_col, "year")
        long[, fleet := match(as.character(fleet_name), fleet_cols)]
        long[, year := as.integer(year)]
        long[, value := as.numeric(value)]
    }

    long = long[year >= dims$.first_ts & year <= dims$.last_ts]
    long[, ts := year - dims$.first_ts + 1L]
    long[, fleet := match(fleet, fish$extract_fleets)]
    long = long[!is.na(fleet)]

    catch_obs_ysf = array(0, dim = c(dims$n_year, dims$n_season, fish$n_fishery))
    for (r in seq_len(nrow(long))) {
        catch_obs_ysf[long$ts[r], 1L, long$fleet[r]] = long$value[r]
    }
    catch_obs_ysf
}

#_____________________________________________________________________________
#' CPUE indices
#'
#' @return a \code{tibble}-like \code{data.table} in opal cpue_data format.
ss3_index_to_opal = function(ss3, cfg, dims, fish) {

    cpue = data.table::as.data.table(ss3_pick(ss3$dat, c("CPUE"), what = "CPUE"))
    data.table::setnames(cpue, tolower(names(cpue)))

    yr_col  = intersect(c("year", "yr"), names(cpue))[1]
    fl_col  = intersect(c("index", "fleet"), names(cpue))[1]
    ob_col  = intersect(c("obs"), names(cpue))[1]
    se_col  = intersect(c("se_log", "se", "err"), names(cpue))[1]
    stopifnot("Unrecognised CPUE columns" = !any(is.na(c(yr_col, fl_col, ob_col, se_col))))

    cpue = cpue[, .(year = as.integer(get(yr_col)),
                    ss3_fleet = as.integer(get(fl_col)),
                    obs = as.numeric(get(ob_col)),
                    se  = as.numeric(get(se_col)))]

    # SS3 encodes "fit this observation but do not include it" with a negative
    # fleet id. Drop those unless the config says otherwise.
    if (!isTRUE(cfg$index$keep_negative_fleets)) cpue = cpue[ss3_fleet > 0]

    cpue = cpue[year >= dims$.first_ts & year <= dims$.last_ts]
    cpue[, ts := year - dims$.first_ts + 1L]
    cpue[, fishery := ss3_opal_fishery(ss3_fleet, ts, cfg, fish)]
    cpue = cpue[fishery %in% fish$index_map$opal_fleet]

    # Index number is the position within the opal index set, not the fishery id.
    cpue[, index := match(fishery, fish$index_map$opal_fleet)]
    cpue[, units := fish$catch_units_f[fishery]]

    # SS3 adds the control file's add_to_survey_CV (variance adjustment factor
    # 1) to each observation's log-scale SE. opal combines its own tau in
    # quadrature, so the SS3 adjustment is added to the input SE here.
    if (isTRUE(cfg$index$add_survey_var_adjust)) {
        cpue[, se := se + ss3_survey_se_add(ss3, ss3_fleet)]
    }

    if (identical(cfg$index$normalise %||% "mean1", "mean1")) {
        cpue[, obs := obs / mean(obs), by = fishery]
    }

    month_labels = as.numeric(cfg$index$month_labels %||% 1)
    cpue[, month := month_labels[((index - 1L) %% length(month_labels)) + 1L]]

    cpue[, metric := "cpue"]
    data.table::setorder(cpue, index, ts)
    cpue[, .(year, month, ts, fishery, index, metric, units, value = obs, se)]
}

#_____________________________________________________________________________
#' Length compositions
#'
#' Returns both the long format opal expects for \code{prep_lf_data} and the
#' per-fishery variance adjustments implied by SS3's \code{mult_by_lencomp_N}
#' and composition lambdas. Fleets with a zero lambda in SS3 are reported so
#' the caller can decide whether to drop them.
#'
#' @return a \code{list} with \code{lf} (long data.table), \code{var_adjust}
#'   and \code{zero_lambda_fleets}.
ss3_lencomp_to_opal = function(ss3, cfg, dims, fish) {

    lc = ss3_pick(ss3$dat, c("lencomp"), default = NULL)
    if (is.null(lc) || !nrow(lc)) {
        return(list(lf = NULL, var_adjust = rep(1, fish$n_fishery),
                    zero_lambda_fleets = integer()))
    }
    lc = data.table::as.data.table(lc)

    yr_col = intersect(c("Yr", "year", "yr"), names(lc))[1]
    fl_col = intersect(c("FltSvy", "fleet", "FltSvy_"), names(lc))[1]
    se_col = intersect(c("Seas", "seas", "month"), names(lc))[1]
    n_col  = intersect(c("Nsamp", "nsamp"), names(lc))[1]
    meta   = c(yr_col, se_col, fl_col, n_col,
               intersect(c("Gender", "Sex", "Part"), names(lc)))

    bin_cols = setdiff(names(lc), meta)
    # r4ss prefixes bin columns with a letter ("l" for 3.24 single-sex, "f" in
    # other versions); strip any leading letters.
    bin_vals = suppressWarnings(as.numeric(gsub("^[A-Za-z]+", "", bin_cols)))
    keep     = !is.na(bin_vals) & bin_vals %in% dims$.len_bins_ss3
    bin_cols = bin_cols[keep]
    # Relabel onto the opal grid (only differs when the last bin was regularised).
    bin_vals = dims$.len_bins[match(bin_vals[keep], dims$.len_bins_ss3)]
    stopifnot("No length bin columns matched the SS3 lbin_vector" = length(bin_cols) > 0)

    wide = lc[, c(meta, bin_cols), with = FALSE]
    data.table::setnames(wide, c(yr_col, fl_col, n_col), c("year", "ss3_fleet", "Nsamp"))
    wide[, year := as.integer(year)]
    wide[, ss3_fleet := as.integer(ss3_fleet)]
    if (!isTRUE(cfg$lf$keep_negative_fleets)) wide = wide[ss3_fleet > 0]
    wide = wide[year >= dims$.first_ts & year <= dims$.last_ts]
    wide[, ts := year - dims$.first_ts + 1L]
    wide[, fishery := ss3_opal_fishery(ss3_fleet, ts, cfg, fish)]
    wide = wide[!is.na(fishery)]

    long = data.table::melt(
        wide[, c("year", "ts", "fishery", "Nsamp", bin_cols), with = FALSE],
        id.vars       = c("year", "ts", "fishery", "Nsamp"),
        variable.name = "bin_col",
        value.name    = "prop"
    )
    long[, bin := bin_vals[match(as.character(bin_col), bin_cols)]]

    # SS3 accepts compositions on any scale (this data file carries raw counts)
    # and normalises each row internally. opal wants pseudo-counts that sum to
    # the input sample size, so normalise first: value = prop / sum(prop) * Nsamp.
    long[, prop := as.numeric(prop)]
    long[, row_total := sum(prop), by = .(year, ts, fishery)]
    long[, value := ifelse(row_total > 0, prop / row_total, 0) * as.numeric(Nsamp)]
    long[, month := 1]
    long[, week  := 1]

    lf = long[, .(year, month, ts, fishery, bin, value, week)]
    data.table::setorder(lf, fishery, ts, bin)

    # Variance adjustment: SS3 multiplies input N; opal divides lf_n by
    # lf_var_adjust. So opal var_adjust = 1 / SS3 multiplier.
    mult = ss3_lencomp_n_mult(ss3, fish)
    var_adjust = rep(1, fish$n_fishery)
    ok = !is.na(mult) & mult > 0
    var_adjust[ok] = 1 / mult[ok]

    list(lf                 = lf,
         var_adjust         = var_adjust,
         zero_lambda_fleets = ss3_zero_lambda_fleets(ss3, like_comp = 4L))
}

#' SS3 mult_by_lencomp_N, expanded to the opal fishery vector
#'
#' Extraction fleets map one-to-one; split index fisheries inherit the value
#' of the SS3 fleet they were carved from.
#' @keywords internal
ss3_lencomp_n_mult = function(ss3, fish) {
    n_fishery = fish$n_fishery
    ss3_fleet = c(fish$extract_fleets, fish$index_map$ss3_fleet)
    va = ss3_pick(ss3$ctl, c("Variance_adjustments", "Variance_adjustment_list",
                             "var_adjust", "Variance_adjustment"), default = NULL)
    out = rep(NA_real_, n_fishery)
    if (is.null(va)) return(out)
    if (!isTRUE(as.integer(ss3$ctl$DoVar_adjust %||% 1L) == 1L)) return(out)
    if (is.data.frame(va) && all(c("factor", "fleet", "value") %in% names(va))) {
        # 3.30 layout: one row per (factor, fleet); factor 4 is mult_by_lencomp_N
        v = va[va$factor == 4, ]
        out = as.numeric(v$value[match(ss3_fleet, v$fleet)])
    } else if (is.matrix(va) || is.data.frame(va)) {
        m = as.matrix(va)
        # 3.24 layout: 6 rows x Nfleet columns. Find the mult_by_lencomp_N row
        # by name (r4ss keeps the control-file comment as the row name), falling
        # back to row 4.
        row = grep("mult_by_lencomp_N", rownames(m))[1]
        if (is.na(row) && nrow(m) >= 4L) row = 4L
        if (!is.na(row)) {
            v = as.numeric(m[row, ])
            out = v[ss3_fleet]
        }
    }
    out
}

#' SS3 add_to_survey_CV (variance adjustment factor 1) for each SS3 fleet id
#' @keywords internal
ss3_survey_se_add = function(ss3, ss3_fleet) {
    va = ss3_pick(ss3$ctl, c("Variance_adjustment_list", "Variance_adjustments"), default = NULL)
    out = rep(0, length(ss3_fleet))
    if (is.null(va) || !isTRUE(as.integer(ss3$ctl$DoVar_adjust %||% 1L) == 1L)) return(out)
    if (is.data.frame(va) && all(c("factor", "fleet", "value") %in% names(va))) {
        v = va[va$factor == 1, ]
        hit = match(ss3_fleet, v$fleet)
    } else {
        m = as.matrix(va)
        row = grep("add_to_survey_CV", rownames(m))[1]
        if (is.na(row)) row = 1L
        v = list(value = as.numeric(m[row, ]))
        hit = ss3_fleet
    }
    out[!is.na(hit)] = as.numeric(v$value[hit[!is.na(hit)]])
    out
}

#' SS3 fleets with a zero composition lambda
#' @keywords internal
ss3_zero_lambda_fleets = function(ss3, like_comp = 4L) {
    lam = ss3_pick(ss3$ctl, c("lambdas", "lambda_changes", "Lambdas"), default = NULL)
    if (is.null(lam)) return(integer())
    lam = as.data.frame(lam)
    lc_col = intersect(c("like_comp", "Like_comp", "likecomp"), names(lam))[1]
    fl_col = intersect(c("fleet/survey", "fleet", "fleet.survey", "Fleet"), names(lam))[1]
    v_col  = intersect(c("value", "Value", "lambda"), names(lam))[1]
    if (any(is.na(c(lc_col, fl_col, v_col)))) return(integer())
    as.integer(lam[[fl_col]][lam[[lc_col]] == like_comp & lam[[v_col]] == 0])
}

#_____________________________________________________________________________
#' Biology
#'
#' Dispatches on \code{cfg$bio$source}:
#'
#'   "endgrowth" : take length, SD, weight, maturity, spawning potential and M
#'                 at age directly from SS3's endgrowth table. Biology is then
#'                 exact and growth parameters affect only the PLA used for
#'                 selectivity. Combined with sel_fa_external, growth becomes
#'                 entirely inert.
#'   "control"   : derive biology from the control file, as a conventional
#'                 SS3 translation. Weight and maturity are then length-basis
#'                 vectors that opal converts through its own PLA, so the
#'                 Schnute approximation to SS3's age-specific K propagates
#'                 into spawning biomass.
#'
#' Implementations are in ss3-biology.R.
#'
#' @param ss3 output of \code{ss3_read_model}.
#' @param cfg converter configuration list.
#' @param dims output of \code{ss3_dims_to_opal}.
#' @return a \code{list} of biology elements plus a \code{growth_fit}
#'   diagnostic and a \code{source} label.
ss3_bio_to_opal = function(ss3, cfg, dims) {

    source = cfg$bio$source %||% "endgrowth"
    if (!source %in% c("endgrowth", "control")) {
        stop(sprintf("cfg$bio$source must be 'endgrowth' or 'control', not '%s'.",
                     source), call. = FALSE)
    }

    if (identical(source, "endgrowth")) {
        if (is.null(ss3$rep)) {
            stop(paste0("cfg$bio$source is 'endgrowth' but Report.sso was not read. ",
                        "Run the SS3 model to produce Report.sso, or set ",
                        "cfg$bio$source to 'control'."), call. = FALSE)
        }
        bio = ss3_bio_from_endgrowth(ss3, cfg, dims)
    } else {
        bio = ss3_bio_from_control(ss3, cfg, dims)
    }

    # Growth form: opal's single-k Schnute curve (fitted or from the control
    # file, per route), or SS3's own age-specific K read from the SS3 files.
    growth_model = cfg$growth$model %||% "schnute"
    if (!growth_model %in% c("schnute", "ss3_age_k")) {
        stop(sprintf("cfg$growth$model must be 'schnute' or 'ss3_age_k', not '%s'.",
                     growth_model), call. = FALSE)
    }
    if (identical(growth_model, "ss3_age_k")) bio = ss3_apply_age_k_growth(bio, ss3, cfg, dims)
    bio
}

#_____________________________________________________________________________
#' Selectivity
#'
#' This SS3 model uses age-based selectivity (size_selex pattern 0 for every
#' fleet, age_selex pattern 20/12). opal parameterises selectivity on length
#' and maps to age through the PLA, so the two routes are:
#'
#'   external : sel_fa_external, the SS3 selectivity-at-age matrix injected
#'              directly. Exact, and the right object for the bridging check.
#'   estimated: opal's own length-based curves, with starting values obtained
#'              by fitting the opal curve family to SS3's implied
#'              selectivity-at-length.
#'
#' Both are returned. The vignette can run the external version first to show
#' the dynamics match, then free the length-based parameters against the
#' composition data.
#'
#' @return a \code{list} with \code{sel_fa_external}, \code{sel_type_f},
#'   \code{par_sel} and \code{sel_target_fl}.
ss3_sel_to_opal = function(ss3, cfg, dims, fish) {

    n_fishery = fish$n_fishery
    n_age     = dims$n_age

    sel_fa_external = NULL
    sel_factor      = NULL
    if (!is.null(ss3$rep)) {
        sel_fa_external = ss3_age_selex_matrix(ss3$rep, dims, fish, cfg)
        if (!is.null(sel_fa_external)) {
            sel_factor = attr(sel_fa_external, "factor")
            attr(sel_fa_external, "factor") = NULL
        }
    }

    # opal selectivity type codes: 1 = logistic, 2 = double-normal,
    # 3 = double Richards. These are NOT the SS3 pattern codes.
    sel_type_f = as.integer(cfg$selectivity$opal_type_f %||% rep(1L, n_fishery))
    stopifnot("cfg$selectivity$opal_type_f must have one entry per opal fishery" =
                  length(sel_type_f) == n_fishery)

    n_sel_par = as.integer(cfg$selectivity$n_par %||% 6L)
    par_sel   = matrix(0, nrow = n_fishery, ncol = n_sel_par)

    # Implied selectivity-at-length, used as the fitting target for starting
    # values. sel_l = sel_a projected through the SS3 age-length key.
    sel_target_fl = NULL
    if (!is.null(sel_fa_external) && !is.null(ss3$rep)) {
        sel_target_fl = ss3_sel_at_age_to_length(sel_fa_external, ss3$rep, dims)
    }

    list(sel_fa_external = sel_fa_external,
         sel_factor      = sel_factor,
         sel_type_f      = sel_type_f,
         par_sel         = par_sel,
         sel_target_fl   = sel_target_fl)
}

#' SS3 selectivity-at-age matrix [n_fishery, n_age] for one year
#' @keywords internal
ss3_age_selex_matrix = function(rep, dims, fish, cfg) {
    as_dt = rep$ageselex
    if (is.null(as_dt)) return(NULL)
    as_dt = data.table::as.data.table(as_dt)

    factor_col = intersect(c("Factor", "factor"), names(as_dt))[1]
    fleet_col  = intersect(c("Fleet", "fleet"), names(as_dt))[1]
    yr_col     = intersect(c("Yr", "year"), names(as_dt))[1]
    if (any(is.na(c(factor_col, fleet_col, yr_col)))) return(NULL)

    # SS3 3.24 reports Asel2 (age selectivity combined with size selectivity
    # through the ALK) for fishing fleets only; surveys carry only Asel. Fall
    # back to Asel per fleet, not just when the whole factor is absent.
    want_factor = cfg$selectivity$report_factor %||% "Asel2"
    factors     = unique(c(want_factor, "Asel"))
    ref_yr      = cfg$selectivity$reference_year

    age_cols = as.character(dims$.ss3_ages)
    missing_cols = setdiff(age_cols, names(as_dt))
    if (length(missing_cols)) {
        stop(sprintf("ageselex is missing age columns: %s",
                     paste(missing_cols, collapse = ", ")), call. = FALSE)
    }

    out = matrix(NA_real_, nrow = fish$n_fishery, ncol = dims$n_age)
    used = rep(NA_character_, fish$n_fishery)
    for (f in seq_len(fish$n_fishery)) {
        ss3_f = if (f <= fish$n_extract) fish$extract_fleets[f] else
            fish$index_map$ss3_fleet[match(f, fish$index_map$opal_fleet)]
        for (fac in factors) {
            rows = as_dt[get(factor_col) == fac & get(fleet_col) == ss3_f]
            if (!nrow(rows)) next
            yr = as.integer(ref_yr %||% max(rows[[yr_col]]))
            row = rows[get(yr_col) == yr]
            if (!nrow(row)) next
            out[f, ] = as.numeric(row[1L, age_cols, with = FALSE])
            used[f]  = fac
            break
        }
    }
    if (anyNA(used)) {
        stop(sprintf("No ageselex rows (%s) for opal fisheries %s.",
                     paste(factors, collapse = "/"),
                     paste(which(is.na(used)), collapse = ", ")), call. = FALSE)
    }
    attr(out, "factor") = used
    out
}

#' Project selectivity-at-age onto the length bins using SS3's ALK
#' @keywords internal
ss3_sel_at_age_to_length = function(sel_fa, rep, dims) {
    eg = rep$endgrowth
    if (is.null(eg)) return(NULL)
    age_col = intersect(c("Age_Beg", "int_Age", "Age"), names(eg))[1]
    mu_col  = intersect(c("Len_Beg", "Len_Mid"), names(eg))[1]
    sd_col  = intersect(c("SD_Beg", "SD_Mid", "sd_Beg"), names(eg))[1]
    if (any(is.na(c(age_col, mu_col, sd_col)))) return(NULL)

    eg = eg[!duplicated(eg[[age_col]]), ]
    i  = match(dims$.ss3_ages, eg[[age_col]])
    mu = as.numeric(eg[[mu_col]][i])
    sd = as.numeric(eg[[sd_col]][i])
    if (anyNA(mu) || anyNA(sd)) return(NULL)

    lower = dims$.len_bins
    upper = lower + dims$len_bin_width
    # pla[l, a] = P(length in bin l | age a)
    pla = sapply(seq_along(mu), function(a) {
        stats::pnorm(upper, mu[a], sd[a]) - stats::pnorm(lower, mu[a], sd[a])
    })
    pla = sweep(pla, 2, pmax(colSums(pla), 1e-12), "/")

    # Weight each age's contribution to a bin by P(age | length).
    num = pla %*% t(sel_fa)                  # [n_len, n_fishery]
    den = pmax(rowSums(pla), 1e-12)
    t(num / den)                             # [n_fishery, n_len]
}

#_____________________________________________________________________________
#' Estimated parameters
#'
#' @return a \code{list} of opal parameter starting values.
ss3_pars_to_opal = function(ss3, cfg, dims, fish, bio) {

    # SS3's par file holds the estimates; the control file only the INIT values
    # the run started from. Prefer the par file where it is available.
    sr = ss3_pick(ss3$par %||% list(), c("SR_parms"), default = NULL) %||%
        ss3_pick(ss3$ctl, c("SR_parms", "SRparm"), what = "SR_parms")
    sr = as.data.frame(sr)
    ln_R0   = unname(ss3_parm_value(sr, "SR_LN\\(R0\\)|SR_LN_R0", n_expected = 1L))
    steep   = unname(ss3_parm_value(sr, "SR_BH_steep|steep",      n_expected = 1L))
    sigma_r = unname(ss3_parm_value(sr, "SR_sigmaR|sigmaR",       n_expected = 1L))

    # opal's scale parameter is unfished spawning biomass B0, not R0: opal
    # derives R0 = B0 / SPR0 internally. Take SS3's own unfished SSB where
    # Report.sso is available; otherwise form it as R0 x SPR0 from the
    # converted age-basis biology.
    ssb0 = ss3_derived_value(ss3$rep, c("SSB_Unfished", "SSB_Virgin",
                                        "SPB_Unfished", "SPB_Virgin"))
    if (!is.finite(ssb0)) {
        spr0 = opal_spr0(list(n_age = dims$n_age, M = bio$M,
                              spawning_potential = bio$spawning_potential,
                              maturity = bio$maturity, fecundity = bio$fecundity,
                              sex_ratio = bio$sex_ratio))
        if (!is.finite(spr0)) {
            stop(paste0("Cannot set log_B0: Report.sso was not read and the biology is ",
                        "on a length basis, so SPR0 cannot be formed outside opal."),
                 call. = FALSE)
        }
        ssb0 = exp(ln_R0) * spr0
    }
    log_B0 = log(ssb0)

    # Recruitment deviations, aligned to the opal timestep window. opal's
    # rdev_y[y] sets the first age class at the start of timestep y + 1. When
    # SS3 age 0 is kept that class is SS3's recruitment of y + 1; when age 0
    # is dropped it is SS3's age 1, recruited (and given its deviation) in y.
    # See ss3_recdev_shift().
    shift  = ss3_recdev_shift(dims)
    rdev_y = rep(0, dims$n_year)
    init_rdev_a = rep(0, dims$n_age)
    dev = if (is.null(ss3$par)) NULL else ss3_recdevs(ss3$par, ss3$ctl)
    if (!is.null(dev)) {
        hit = match(dims$.first_ts:dims$.last_ts + shift, dev$year)
        rdev_y = ifelse(is.na(hit), 0, dev$value[hit])
        # The initial age structure: SS3 age s at the start of the first
        # timestep was recruited in first_ts - s, with that year's deviation
        # (the early deviations, and for age 0 the first main deviation).
        # Only applied with cfg$recdev$initial: true.
        if (isTRUE(cfg$recdev$initial)) {
            hit = match(dims$.first_ts - dims$.ss3_ages, dev$year)
            init_rdev_a = ifelse(is.na(hit), 0, dev$value[hit])
        }
    }

    # Initial-equilibrium fishing mortality (SS3 init_F, annual and by fleet)
    log_init_F_f = rep(log(1e-8), fish$n_fishery)
    init_F = ss3_init_F(ss3)
    if (length(init_F)) {
        f = match(init_F$fleet, fish$extract_fleets)
        ok = !is.na(f)
        log_init_F_f[f[ok]] = log(pmax(init_F$value[ok], 1e-8))
    }

    n_index = fish$n_index

    pars = list(
        log_B0         = log_B0,
        log_h          = log(steep),
        log_sigma_r    = log(sigma_r),
        log_cpue_q     = rep(0, n_index),          # q = 1 on the natural scale
        cpue_creep     = rep(0, n_index),
        log_cpue_tau   = rep(-Inf, n_index),
        log_cpue_omega = rep(0, n_index),
        log_lf_tau     = rep(log(0.1), fish$n_fishery),
        rdev_y         = rdev_y,
        init_rdev_a    = init_rdev_a,
        log_init_F_f   = log_init_F_f,
        par_sel        = matrix(0, nrow = fish$n_fishery,
                                ncol = as.integer(cfg$selectivity$n_par %||% 6L)),
        log_L1         = log(bio$L1),
        log_L2         = log(bio$L2),
        log_k          = log(bio$k),
        log_CV1        = log(bio$CV1),
        log_CV2        = log(bio$CV2)
    )
    # Age-specific K multipliers, cfg$growth$model = "ss3_age_k" only.
    # Assigning NULL adds nothing, so the Schnute path is unchanged.
    pars$log_k_mult = bio$log_k_mult
    pars
}

#' Year offset between opal's rdev_y index and SS3's recruitment year
#'
#' 1 when SS3 age 0 is the first opal age (opal's first age class at t + 1 is
#' SS3's recruitment of t + 1), 0 when age 0 is dropped (that class is SS3's
#' age 1, recruited in t).
#' @keywords internal
ss3_recdev_shift = function(dims) as.integer(dims$.ss3_ages[1] == 0L)

#' Initial-equilibrium F by fleet from the par file (or the control file)
#' @keywords internal
ss3_init_F = function(ss3) {
    ctl_if = ss3_pick(ss3$ctl, c("init_F"), default = NULL)
    if (is.null(ctl_if) || !NROW(ctl_if)) return(NULL)
    # 3.30 names the rows InitF_seas_<s>_flt_<f><name>; 3.24 models (whose
    # conversions predate this) are left at opal's default of no initial F
    is_330 = grepl("_flt_[0-9]+", rownames(ctl_if))
    if (!any(is_330)) return(NULL)
    fleet = as.integer(sub(".*_flt_([0-9]+).*", "\\1", rownames(ctl_if)))
    value = as.numeric(ctl_if[, "INIT"])
    par_if = ss3_pick(ss3$par %||% list(), c("init_F"), default = NULL)
    if (!is.null(par_if) && length(par_if) == length(value)) value = as.numeric(par_if)
    data.frame(fleet = fleet, value = value)
}

#' Recruitment deviations from the SS3 par file, with their years
#' @keywords internal
ss3_recdevs = function(par, ctl) {
    early = ss3_pick(par, c("recdev_early", "recdev_early_parm"), default = NULL)
    main  = ss3_pick(par, c("recdev1", "recdev2", "recdev"),      default = NULL)
    if (is.null(main)) return(NULL)

    to_vec = function(x) {
        # SS_readpar_3.24() returns a two-column matrix (year, recdev).
        if (is.matrix(x) && !is.null(colnames(x))) x = as.data.frame(x)
        if (is.data.frame(x)) {
            v_col = intersect(c("recdev", "value", "V1"), names(x))[1]
            y_col = intersect(c("year", "Yr"), names(x))[1]
            if (!is.na(y_col)) return(list(year = as.integer(x[[y_col]]),
                                           value = as.numeric(x[[v_col]])))
            return(list(year = NULL, value = as.numeric(x[[v_col]])))
        }
        list(year = NULL, value = as.numeric(x))
    }

    m = to_vec(main)
    if (!is.null(m$year)) {
        out = data.frame(year = m$year, value = m$value)
    } else {
        first_main = as.integer(ss3_pick(ctl, c("MainRdevYrFirst", "first_yr_main_recr_devs",
                                                "MainRdevYrFirst_"), default = NA))
        if (is.na(first_main)) return(NULL)
        out = data.frame(year = seq(first_main, length.out = length(m$value)),
                         value = m$value)
    }

    if (!is.null(early)) {
        e = to_vec(early)
        if (!is.null(e$year)) {
            out = rbind(data.frame(year = e$year, value = e$value), out)
        } else {
            first_early = as.integer(ss3_pick(ctl, c("recdev_early_start",
                                                     "first_yr_early_recr_devs"),
                                              default = NA))
            if (!is.na(first_early)) {
                out = rbind(data.frame(year = seq(first_early, length.out = length(e$value)),
                                       value = e$value), out)
            }
        }
    }
    out[!duplicated(out$year), ]
}

#_____________________________________________________________________________
#' Orchestrator: SS3 directory plus config -> opal inputs
#'
#' @param cfg converter configuration, from \code{ss3_read_config}.
#' @param dir_root \code{character} repository root used to resolve relative
#'   paths in the config.
#' @return a \code{list} with \code{data}, \code{parameters}, \code{lf},
#'   \code{selectivity}, \code{meta} and \code{ss3}.
ss3_to_opal = function(cfg, dir_root = getwd()) {

    dir_ss3 = if (isTRUE(cfg$absolute_paths)) cfg$ss3_dir else
        file.path(dir_root, cfg$ss3_dir)

    ss3 = ss3_read_model(dir_ss3      = dir_ss3,
                         dat_file     = cfg$files$dat %||% "data.ss_new",
                         ctl_file     = cfg$files$ctl %||% "control.ss_new",
                         par_file     = cfg$files$par %||% "ss.par",
                         version      = as.character(cfg$ss3_version %||% "3.24"),
                         read_report  = isTRUE(cfg$files$read_report %||% TRUE))

    dims  = ss3_dims_to_opal(ss3, cfg)
    fish  = ss3_fisheries_to_opal(ss3, cfg, dims)
    catch = ss3_catch_to_opal(ss3, cfg, dims, fish)
    cpue  = ss3_index_to_opal(ss3, cfg, dims, fish)
    lfout = ss3_lencomp_to_opal(ss3, cfg, dims, fish)
    bio   = ss3_bio_to_opal(ss3, cfg, dims)
    sel   = ss3_sel_to_opal(ss3, cfg, dims, fish)
    pars  = ss3_pars_to_opal(ss3, cfg, dims, fish, bio)

    data = list(
        age_a          = dims$age_a,
        n_age          = dims$n_age,
        min_age        = dims$min_age,
        n_season       = dims$n_season,
        n_fishery      = fish$n_fishery,
        n_year         = dims$n_year,
        first_yr       = dims$first_yr,
        last_yr        = dims$last_yr,
        years          = dims$years,

        len_bin_start  = dims$len_bin_start,
        len_bin_width  = dims$len_bin_width,
        n_len          = dims$n_len,

        first_yr_catch = 1L,
        catch_units_f  = fish$catch_units_f,
        catch_obs_ysf  = catch,

        cpue_switch    = 1L,
        n_index        = fish$n_index,
        # Plain data.frame, not data.table: opal_obj() rejects non-portable
        # state, and a data.table carries an externalptr .internal.selfref
        # attribute that trips .opal_contains_nonportable().
        cpue_data      = as.data.frame(cpue),

        lw_a           = bio$lw_a,
        lw_b           = bio$lw_b,
        maturity       = bio$maturity,
        fecundity      = bio$fecundity,
        sex_ratio      = bio$sex_ratio,
        M              = bio$M,

        # Supplied only by the endgrowth route. opal_model() checks for these
        # with exists(), so a NULL element is not the same as an absent one:
        # the NULLs are stripped below.
        weight             = bio$weight,
        spawning_potential = bio$spawning_potential,

        A1             = bio$A1,
        A2             = bio$A2,

        # Age-specific K growth (cfg$growth$model = "ss3_age_k"); NULL and
        # stripped below otherwise.
        k_mult_ages              = bio$k_mult_ages,
        growth_plus_group_adjust = bio$plus_group_adjust,

        sel_type_f     = sel$sel_type_f,

        # SS3's recruitment bias-adjustment ramp (Report.sso only). Without it
        # opal applies the full lognormal correction in every timestep.
        bias_adj_y     = ss3_bias_adj(ss3$rep, dims),
        init_bias_adj_a = if (isTRUE(cfg$recdev$initial)) ss3_init_bias_adj(ss3$rep, dims) else NULL
    )
    data = data[!vapply(data, is.null, logical(1L))]

    meta = list(
        config       = cfg,
        provenance   = ss3$provenance,
        dims         = dims,
        fisheries    = fish,
        bio_source   = bio$source,
        growth_fit   = bio$growth_fit,
        lf_var_adjust        = lfout$var_adjust,
        lf_zero_lambda_fleets = lfout$zero_lambda_fleets,
        converted_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )

    # The data list is handed to opal_obj(), which refuses any element holding
    # a closure, environment or external pointer. Sweep it once here so a
    # stray data.table does not surface as an opaque construction failure.
    data = ss3_make_portable(data)

    list(data        = data,
         parameters  = pars,
         lf          = lfout$lf,
         selectivity = sel,
         meta        = meta,
         ss3         = ss3)
}

#' Null-coalescing operator
#' @keywords internal
`%||%` = function(x, y) if (is.null(x)) y else x
