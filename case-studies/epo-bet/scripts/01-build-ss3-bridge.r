# Nicholas Ducharme-Barth
# 2026/10/30
# Build and run the SS3 bridge from the IATTC SAC-14 exploratory assessment
# model of bigeye tuna in the eastern Pacific Ocean to an SS3 model whose
# structure opal can represent.
#
# Each step changes one thing, starting from the previous step's estimates
# (control.ss_new), so its effect can be measured on its own:
#   00-base        the model as supplied (run with an SS3 3.30.22 beta build),
#                  re-run with 3.30.22.1
#   01-no-regime   initial-equilibrium regime offset (block year 16) removed
# 01 is the SS3 model the opal conversion starts from.
#
# Growth needs no bridge step: SS3 GrowthModel 8 (growth cessation) and
# CV_Growth_Pattern 2 (SD linear in length) are both available in opal.
#
# The model is quarterly: each SS3 "year" is a quarter (17-204), and ages are
# quarters (0-40).
#
# OUTPUTS:
#   case-studies/epo-bet/model-files/ss3/<step>/   SS3 inputs and outputs
#   case-studies/epo-bet/vignettes/data/epo-bet-ss3-bridge-cache.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Packages
    library(data.table)
    library(r4ss)

#_____________________________________________________________________________
# Project paths
    proj_dir  = this.path::this.proj()
    dir_ss3   = file.path(proj_dir, "case-studies", "epo-bet", "model-files", "ss3")
    dir_cache = file.path(proj_dir, "case-studies", "epo-bet", "vignettes", "data")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
    # EPO_BET_NOHESS=true skips the Hessian (for quick iteration only; the
    # committed runs and the convergence tests use the full run)
    extras = if (identical(Sys.getenv("EPO_BET_NOHESS"), "true")) "-nohess" else ""
    exe = file.path(proj_dir, "executables", "stock-synthesis", "3.30.22.1",
                    if (.Platform$OS.type == "windows") "ss3_win.exe" else "ss3_linux")
    if (.Platform$OS.type != "windows") Sys.chmod(exe, "0755")

#_____________________________________________________________________________
# Helpers

    # Write an SS3 input list to dir and estimate the model there. The Hessian
    # is inverted so the bridge records whether each step is positive definite.
    # SS3 starts from par_lines (the previous step's par file, edited to the
    # new parameter layout) when given, otherwise from the control file's INIT
    # values. The control file does not hold the recruitment deviations, and
    # starting them at 0 crashes this stock under its catches.
    run_step = function(inputs, step, par_lines = NULL) {
        dir = file.path(dir_ss3, step)
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        inputs$start$init_values_src = if (is.null(par_lines)) 0 else 1
        SS_write(inputs, dir = dir, overwrite = TRUE)
        if (!is.null(par_lines)) writeLines(par_lines, file.path(dir, "ss3.par"))
        status = r4ss::run(dir, exe = exe, extras = extras, skipfinished = FALSE, verbose = FALSE)
        if (!file.exists(file.path(dir, "Report.sso"))) {
            stop(sprintf("SS3 did not finish for %s (%s).", step, status), call. = FALSE)
        }
        invisible(dir)
    }

    # Read a finished step back as the starting point for the next one
    read_step = function(step) SS_read(file.path(dir_ss3, step), ss_new = TRUE)

    # A step's par file without the parameter whose header line is `label`
    # (each par-file parameter is a "# label:" line followed by its value)
    par_without = function(step, label) {
        x = readLines(file.path(dir_ss3, step, "ss3.par"))
        i = which(x == paste0("# ", label, ":"))
        stopifnot("parameter not found in the par file" = length(i) == 1L)
        x[-c(i, i + 1L)]
    }

#_____________________________________________________________________________
# 00: the model as supplied. Its starter reads ss3.par (init_values_src = 1).
# The supplied data and control files (BET-EPO.dat, BET-EPO.ctl) are renamed
# data.ss and control.ss. The supplied run's likelihood components are kept in
# supplied-likelihood-3.30.22beta.csv.
# Re-run only while Report.sso still comes from the supplied beta run.
    dir_b0 = file.path(dir_ss3, "00-base")
    if (!grepl("3.30.22.1", readLines(file.path(dir_b0, "Report.sso"), n = 1L), fixed = TRUE)) {
        r4ss::run(dir_b0, exe = exe, extras = extras, skipfinished = FALSE, verbose = FALSE)
    }

#_____________________________________________________________________________
# 01: no regime offset on the initial equilibrium recruitment. SS3 estimates a
# block offset (-0.41) on ln(R0) for block year 16, the year before the model
# starts, so the initial age structure comes from a lower equilibrium
# recruitment than the stock-recruitment curve. opal has no such offset.
    inputs = read_step("00-base")
    ctl = inputs$ctl
    ctl$SR_parms["SR_regime", c("INIT", "Block", "Block_Fxn")] = 0
    ctl$SR_parms_tv = NULL
    inputs$ctl = ctl
    # start from 00's estimates, without the offset (SR_parm 6). Skipped when
    # the step has already been run (EPO_BET_RERUN=true runs it again).
    done = function(step) file.exists(file.path(dir_ss3, step, "admodel.cov")) &&
        !identical(Sys.getenv("EPO_BET_RERUN"), "true")
    if (!done("01-no-regime")) {
        run_step(inputs, "01-no-regime", par_lines = par_without("00-base", "SR_parm[6]"))
    }

#_____________________________________________________________________________
# Summarise the bridge
    steps  = c("00-base", "01-no-regime")
    labels = c("B0 SAC-14 model", "B1 no regime offset")
    # Model quarter -> calendar time: quarter 17 = 1975 Q1 (see README)
    qtr_to_year = function(q) 1975 + (q - 17) / 4
    end_q = 204L
    reps = lapply(steps, function(s)
        SS_output(file.path(dir_ss3, s), covar = TRUE, verbose = FALSE, printstats = FALSE))
    names(reps) = steps

    dq = function(rep, label) {
        x = rep$derived_quants
        v = x$Value[match(label, x$Label)]
        if (length(v) == 0L || is.na(v)) NA_real_ else v
    }
    lk_val = function(lk, comp) if (comp %in% rownames(lk)) lk[comp, "values"] else NA_real_
    summary_dt = rbindlist(lapply(seq_along(steps), function(i) {
        rep = reps[[i]]
        lk = rep$likelihoods_used
        data.table(
            step          = steps[i],
            label         = labels[i],
            n_par         = rep$N_estimated_parameters,
            nll           = lk_val(lk, "TOTAL"),
            nll_survey    = lk_val(lk, "Survey"),
            nll_length    = lk_val(lk, "Length_comp"),
            nll_sizefreq  = lk_val(lk, "SizeFreq"),
            nll_recruit   = lk_val(lk, "Recruitment"),
            max_gradient  = rep$maximum_gradient_component,
            hessian_ok    = !is.null(rep$CoVar) && isTRUE(rep$log_det_hessian > -Inf),
            ln_r0         = rep$parameters["SR_LN(R0)", "Value"],
            regime        = if ("SR_regime_BLK1add_16" %in% rownames(rep$parameters))
                                rep$parameters["SR_regime_BLK1add_16", "Value"] else 0,
            regime_se     = if ("SR_regime_BLK1add_16" %in% rownames(rep$parameters))
                                rep$parameters["SR_regime_BLK1add_16", "Parm_StDev"] else NA_real_,
            # mean of the early recruitment deviations (the initial age structure)
            early_dev_mean = mean(rep$parameters$Value[grepl("^Early_InitAge_", rep$parameters$Label)]),
            ssb0          = dq(rep, "SSB_Virgin"),
            ssb_end       = dq(rep, paste0("SSB_", end_q)),
            depletion     = dq(rep, paste0("Bratio_", end_q)),
            # F_<yr> is F / F_MSY (starter F_report_basis 2; F = sum of full Fs)
            f_end         = dq(rep, paste0("F_", end_q)),
            ssb_msy       = dq(rep, "SSB_MSY"),
            n_on_bounds   = sum(grepl("LO|HI", rep$parameters$Status))
        )
    }))
    print(summary_dt)

    ts_dt = rbindlist(lapply(seq_along(steps), function(i) {
        ts = as.data.table(reps[[i]]$timeseries)[Era %in% c("VIRG", "INIT", "TIME")]
        data.table(step = steps[i], label = labels[i], era = ts$Era, quarter = ts$Yr,
                   year = qtr_to_year(ts$Yr), ssb = ts$SpawnBio, recruits = ts$Recruit_0)
    }))
    ts_dt[, depletion := ssb / summary_dt$ssb0[match(step, summary_dt$step)]]

    supplied = fread(file.path(dir_b0, "supplied-likelihood-3.30.22beta.csv"))

    saveRDS(list(summary = summary_dt, timeseries = ts_dt, steps = steps, labels = labels,
                 supplied_likelihood = supplied,
                 ss3_version = "3.30.22.1", made_at = format(Sys.time(), tz = "UTC", usetz = TRUE)),
            file.path(dir_cache, "epo-bet-ss3-bridge-cache.rds"))
