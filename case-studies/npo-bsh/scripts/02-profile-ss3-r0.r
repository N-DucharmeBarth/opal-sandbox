# Nicholas Ducharme-Barth
# 2026/10/28
# Likelihood profiles over ln(R0) for the NPO blue shark SS3 bridge models.
#
# The single-sex bridge models' depletion moved from 0.28 to 0.54 when the
# F-ballpark penalty (worth < 1 likelihood unit) was switched off, so the
# stock's scale is only weakly informed. Profiling ln(R0) shows how flat the
# likelihood is, which data components pull where, and that each step's
# estimate is the profile's minimum rather than a local optimum.
#
# Each profile point is a copy of the step with SR_LN(R0) fixed and every
# other parameter re-estimated from the step's estimates (no Hessian). The
# runs are made in a temporary directory; only the summary is kept.
#
# OUTPUTS:
#   case-studies/npo-bsh/vignettes/data/npo-bsh-ss3-r0-profile-cache.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Packages
    library(data.table)
    library(r4ss)
    library(parallel)

#_____________________________________________________________________________
# Project paths
    proj_dir  = this.path::this.proj()
    dir_ss3   = file.path(proj_dir, "case-studies", "npo-bsh", "model-files", "ss3")
    dir_cache = file.path(proj_dir, "case-studies", "npo-bsh", "vignettes", "data")
    dir_runs  = file.path(tempdir(), "npo-bsh-r0-profile")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
    exe = file.path(proj_dir, "executables", "stock-synthesis", "3.30.22.1",
                    if (.Platform$OS.type == "windows") "ss3_win.exe" else "ss3_linux")

#_____________________________________________________________________________
# Profile settings
    steps   = c("00-s6-base", "01-single-sex", "02e-no-f-ballpark")
    ln_r0   = seq(9.2, 10.8, by = 0.1)
    n_cores = max(1L, min(4L, parallel::detectCores()))

#_____________________________________________________________________________
# One profile point
    profile_point = function(step, value) {
        dir = file.path(dir_runs, sprintf("%s_%.2f", step, value))
        inputs = SS_read(file.path(dir_ss3, step), ss_new = TRUE)
        inputs$start$init_values_src = 0
        i = grep("SR_LN", rownames(inputs$ctl$SR_parms))
        inputs$ctl$SR_parms[i, c("INIT", "PHASE")] = c(value, -1)
        SS_write(inputs, dir = dir, overwrite = TRUE)
        r4ss::run(dir, exe = exe, extras = "-nohess", skipfinished = FALSE, verbose = FALSE)
        if (!file.exists(file.path(dir, "Report.sso"))) {
            return(data.table(step = step, ln_r0 = value, converged = FALSE))
        }
        rep = SS_output(dir, covar = FALSE, verbose = FALSE, printstats = FALSE)
        lk = rep$likelihoods_used
        dq = rep$derived_quants
        val = function(label) dq$Value[match(label, dq$Label)]
        comp = function(name) if (name %in% rownames(lk)) lk[name, "values"] else 0
        data.table(
            step         = step,
            ln_r0        = value,
            converged    = rep$maximum_gradient_component < 1e-3,
            max_gradient = rep$maximum_gradient_component,
            nll          = comp("TOTAL"),
            survey       = comp("Survey"),
            length_comp  = comp("Length_comp"),
            size_freq    = comp("SizeFreq"),
            recruitment  = comp("Recruitment"),
            f_ballpark   = comp("F_Ballpark"),
            ssb0         = val("SSB_Virgin"),
            depletion    = val("Bratio_2020"),
            f_2020       = val("F_2020")
        )
    }

    grid = CJ(step = steps, ln_r0 = ln_r0)
    res = mclapply(seq_len(nrow(grid)), function(i)
        tryCatch(profile_point(grid$step[i], grid$ln_r0[i]),
                 error = function(e) data.table(step = grid$step[i], ln_r0 = grid$ln_r0[i],
                                                converged = FALSE)),
        mc.cores = n_cores)
    profile_dt = rbindlist(res, fill = TRUE)
    setorder(profile_dt, step, ln_r0)
    profile_dt[, delta_nll := nll - min(nll[converged], na.rm = TRUE), by = step]
    print(profile_dt[, .(step, ln_r0, converged, delta_nll, depletion, survey, length_comp)])

    # Each step's own estimate, for reference
    estimates = rbindlist(lapply(steps, function(s) {
        rep = SS_output(file.path(dir_ss3, s), covar = FALSE, verbose = FALSE, printstats = FALSE)
        data.table(step = s, ln_r0 = rep$parameters[grep("SR_LN", rep$parameters$Label), "Value"],
                   nll = rep$likelihoods_used["TOTAL", "values"],
                   depletion = rep$derived_quants$Value[match("Bratio_2020", rep$derived_quants$Label)])
    }))

    saveRDS(list(profile = profile_dt, estimates = estimates, ln_r0 = ln_r0,
                 made_at = format(Sys.time(), tz = "UTC", usetz = TRUE)),
            file.path(dir_cache, "npo-bsh-ss3-r0-profile-cache.rds"))
