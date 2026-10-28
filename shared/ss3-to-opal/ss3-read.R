# Nicholas Ducharme-Barth
# 2026/09/18
# Read SS3 model files via r4ss and capture provenance.
#
# These functions are deliberately thin. They do not reshape anything; they
# just get the SS3 objects into memory in a predictable container and record
# enough metadata that a downstream opal object can be traced back to the
# exact files it came from.
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
#' Defensive element accessor for r4ss lists
#'
#' r4ss element names drift between SS3 versions and between r4ss releases.
#' Rather than indexing a single hard-coded name, look up the first match from
#' a candidate vector and fail loudly with the names that were actually there.
#'
#' @param x a \code{list} returned by an r4ss reader.
#' @param candidates \code{character} vector of candidate element names, in
#'   order of preference.
#' @param default value returned when no candidate matches, which may itself
#'   be \code{NULL}. If not supplied, the function errors instead.
#' @param what \code{character} label used in the error message.
#' @return the matched element, or \code{default}.
ss3_pick = function(x, candidates, default, what = "element") {
    hit = intersect(candidates, names(x))
    if (length(hit) == 0) {
        if (missing(default)) {
            stop(sprintf(
                "Could not find %s. Tried: %s.\nAvailable names: %s",
                what,
                paste(candidates, collapse = ", "),
                paste(names(x), collapse = ", ")
            ), call. = FALSE)
        }
        return(default)
    }
    x[[hit[1]]]
}

#_____________________________________________________________________________
#' Pull an INIT value out of an r4ss parameter data.frame by row name
#'
#' r4ss returns MG_parms / SR_parms / selex parms as data.frames whose row
#' names are the SS3 parameter labels. Match on a regex so that partial or
#' pattern-based lookups work across models.
#'
#' @param parms a \code{data.frame} of SS3 parameter lines.
#' @param pattern \code{character} regex matched against \code{rownames(parms)}.
#' @param column \code{character} column to extract (default \code{"INIT"}).
#' @param n_expected \code{integer} expected number of matches, or \code{NA} to
#'   accept any number.
#' @return a named \code{numeric} vector of the matched values.
ss3_parm_value = function(parms, pattern, column = "INIT", n_expected = NA_integer_) {
    stopifnot(is.data.frame(parms))
    col = intersect(c(column, tolower(column), toupper(column)), names(parms))
    if (length(col) == 0) {
        stop(sprintf("Column '%s' not present. Available: %s",
                     column, paste(names(parms), collapse = ", ")), call. = FALSE)
    }
    idx = grep(pattern, rownames(parms))
    if (length(idx) == 0) {
        stop(sprintf("No parameter row matched '%s'.\nAvailable rows: %s",
                     pattern, paste(rownames(parms), collapse = ", ")), call. = FALSE)
    }
    if (!is.na(n_expected) && length(idx) != n_expected) {
        stop(sprintf("Pattern '%s' matched %d rows; expected %d (%s).",
                     pattern, length(idx), n_expected,
                     paste(rownames(parms)[idx], collapse = ", ")), call. = FALSE)
    }
    stats::setNames(as.numeric(parms[idx, col[1]]), rownames(parms)[idx])
}

#_____________________________________________________________________________
#' Read an SS3 model directory into a single container
#'
#' Reads the data file, control file, parameter file and, when a Report.sso is
#' present, the full SS_output. Derived quantities (M-at-age, length-at-age,
#' selectivity-at-age) are preferred from Report.sso downstream because they
#' resolve SS3 internal transformations that the control file only implies.
#'
#' @param dir_ss3 \code{character} path to the SS3 model directory.
#' @param dat_file \code{character} data file name.
#' @param ctl_file \code{character} control file name.
#' @param par_file \code{character} parameter file name, or \code{NA} to skip.
#' @param version \code{character} SS3 version string passed to r4ss.
#' @param read_report \code{logical}; read Report.sso via \code{SS_output}.
#' @return a \code{list} with elements \code{dat}, \code{ctl}, \code{par},
#'   \code{rep} and \code{provenance}.
ss3_read_model = function(dir_ss3,
                          dat_file    = "data.ss_new",
                          ctl_file    = "control.ss_new",
                          par_file    = "ss.par",
                          version     = "3.24",
                          read_report = TRUE) {

    if (!requireNamespace("r4ss", quietly = TRUE)) {
        stop("Package 'r4ss' is required. renv::install('r4ss/r4ss').", call. = FALSE)
    }
    if (!dir.exists(dir_ss3)) {
        stop(sprintf("SS3 directory not found: %s", dir_ss3), call. = FALSE)
    }

    path_dat = file.path(dir_ss3, dat_file)
    path_ctl = file.path(dir_ss3, ctl_file)

    dat = r4ss::SS_readdat(file = path_dat, version = version, verbose = FALSE)

    # For 3.24 the model dimensions are passed explicitly rather than via
    # use_datlist = TRUE: SS_readctl_3.24() only sets Do_AgeKey when the data
    # file has ageing-error definitions, and otherwise fails on if (NULL).
    if (as.numeric(version) >= 3.3) {
        ctl = r4ss::SS_readctl(file = path_ctl, version = version,
                               use_datlist = TRUE, datlist = dat, verbose = FALSE)
    } else {
        n_survey = as.integer(dat$Nsurveys %||% 0L)
        ctl = r4ss::SS_readctl(file         = path_ctl,
                               version      = version,
                               use_datlist  = FALSE,
                               nseas        = dat$nseas,
                               N_areas      = dat$N_areas,
                               Nages        = dat$Nages,
                               Nsexes       = dat$Nsexes,
                               Npopbins     = dat$N_lbinspop %||% NA,
                               Nfleet       = dat$Nfleet,
                               Nsurveys     = n_survey,
                               Do_AgeKey    = if (isTRUE(dat$N_ageerror_definitions > 0))
                                   as.integer(any(dat$ageerror[seq(2, nrow(dat$ageerror), by = 2), 1] < 0))
                                   else 0L,
                               N_tag_groups = dat$N_tag_groups %||% 0L,
                               N_CPUE_obs   = vapply(seq_len(dat$Nfleet + n_survey), function(i)
                                   sum(dat$CPUE[, "index"] == i), integer(1L)),
                               verbose      = FALSE)
        ctl$fleetnames = dat$fleetnames
    }

    par = NULL
    if (!is.na(par_file) && file.exists(file.path(dir_ss3, par_file))) {
        # SS_readpar_3.24 needs the dat/ctl lists to know the parameter layout.
        # Pass the lists already read: given paths, it re-reads the control
        # file with use_datlist = TRUE and hits the Do_AgeKey failure above.
        par = r4ss::SS_readpar_3.24(parfile   = file.path(dir_ss3, par_file),
                                    datsource = dat,
                                    ctlsource = ctl,
                                    verbose   = FALSE)
    }

    rep = NULL
    if (isTRUE(read_report) && file.exists(file.path(dir_ss3, "Report.sso"))) {
        rep = r4ss::SS_output(dir         = dir_ss3,
                              covar       = FALSE,
                              forecast    = FALSE,
                              verbose     = FALSE,
                              printstats  = FALSE,
                              hidewarn    = TRUE)
    }

    files = list.files(dir_ss3, full.names = TRUE)
    files = files[!dir.exists(files)]

    provenance = list(
        dir_ss3      = normalizePath(dir_ss3, winslash = "/", mustWork = FALSE),
        dat_file     = dat_file,
        ctl_file     = ctl_file,
        par_file     = par_file,
        ss3_version  = version,
        r4ss_version = as.character(utils::packageVersion("r4ss")),
        r_version    = as.character(getRversion()),
        read_at      = format(Sys.time(), tz = "UTC", usetz = TRUE),
        md5          = tools::md5sum(files)
    )

    list(dat = dat, ctl = ctl, par = par, rep = rep, provenance = provenance)
}

#_____________________________________________________________________________
#' Read a converter configuration file
#'
#' All decisions that are not derivable from the SS3 files themselves live in
#' a YAML sidecar so the conversion is reproducible and reviewable.
#'
#' @param path \code{character} path to the YAML configuration.
#' @return a \code{list}.
ss3_read_config = function(path) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
        stop("Package 'yaml' is required to read converter configs.", call. = FALSE)
    }
    if (!file.exists(path)) {
        stop(sprintf("Config not found: %s", path), call. = FALSE)
    }
    cfg = yaml::read_yaml(path)
    cfg$config_path = normalizePath(path, winslash = "/", mustWork = FALSE)
    cfg$config_md5  = unname(tools::md5sum(path))
    cfg
}
