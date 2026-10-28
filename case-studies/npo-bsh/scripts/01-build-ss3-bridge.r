# Nicholas Ducharme-Barth
# 2026/10/28
# Build and run the SS3 bridge from the 2022 ISC North Pacific blue shark
# assessment (ensemble model S6_base, two sexes) to a single-sex SS3 model whose
# structure opal can represent.
#
# Each step changes one thing, starting from the previous step's estimates
# (control.ss_new), so its effect can be measured on its own:
#   00-s6-base          S6 as supplied (SS3 3.30.19), re-run with 3.30.22.1
#   01-single-sex       Nsexes = -1 (SSB x female fraction 0.5), female
#                       growth, M and weight-length; sex x length compositions
#                       summed over sex; sex-offset selectivity parameters
#                       removed (each curve becomes the male base curve)
#   02a-f18-dblnorm     F18 cubic spline (pattern 27) -> double normal (24)
#   02b-f8-split        F8's two selectivity blocks (1973-81, 1982-92) -> F8
#                       and a new fleet 32 (catch and compositions 1982-92)
#   02c-no-f19-devs     F19 annual deviations on the peak removed
#   02d-no-regime       initial-equilibrium regime offset (1970) removed
#   02e-no-f-ballpark   F-ballpark penalty switched off
# 02e is the SS3 model the opal conversion starts from.
#
# OUTPUTS:
#   case-studies/npo-bsh/model-files/ss3/<step>/   SS3 inputs and outputs
#   case-studies/npo-bsh/vignettes/data/npo-bsh-ss3-bridge-cache.rds
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
    dir_ss3   = file.path(proj_dir, "case-studies", "npo-bsh", "model-files", "ss3")
    dir_cache = file.path(proj_dir, "case-studies", "npo-bsh", "vignettes", "data")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
    # NPO_BSH_NOHESS=true skips the Hessian (for quick iteration only; the
    # committed runs and the convergence tests use the full run)
    extras = if (identical(Sys.getenv("NPO_BSH_NOHESS"), "true")) "-nohess" else ""
    exe = file.path(proj_dir, "executables", "stock-synthesis", "3.30.22.1",
                    if (.Platform$OS.type == "windows") "ss3_win.exe" else "ss3_linux")
    if (.Platform$OS.type != "windows") Sys.chmod(exe, "0755")

#_____________________________________________________________________________
# Helpers

    # Write an SS3 input list to dir and estimate the model there. The Hessian
    # is inverted so the bridge records whether each step is positive definite.
    run_step = function(inputs, step) {
        dir = file.path(dir_ss3, step)
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        # start from the control file's INIT values (the previous step's
        # estimates), not a par file whose layout may have changed
        inputs$start$init_values_src = 0
        SS_write(inputs, dir = dir, overwrite = TRUE)
        status = r4ss::run(dir, exe = exe, extras = extras, skipfinished = FALSE, verbose = FALSE)
        if (!file.exists(file.path(dir, "Report.sso"))) {
            stop(sprintf("SS3 did not finish for %s (%s).", step, status), call. = FALSE)
        }
        invisible(dir)
    }

    # Read a finished step back as the starting point for the next one
    read_step = function(step) SS_read(file.path(dir_ss3, step), ss_new = TRUE)

    # Rows of a parameter table whose names match a pattern
    drop_rows = function(x, pattern) x[!grepl(pattern, rownames(x)), , drop = FALSE]

#_____________________________________________________________________________
# 00: S6 as supplied. Its starter reads the par file (init_values_src = 1);
# SS3 3.30.22 and later look for ss3.par.
# Re-run only while Report.sso still comes from the supplied 3.30.19 run.
    dir_b0 = file.path(dir_ss3, "00-s6-base")
    if (!grepl("3.30.22.1", readLines(file.path(dir_b0, "Report.sso"), n = 1L), fixed = TRUE)) {
        file.copy(file.path(dir_b0, "ss.par"), file.path(dir_b0, "ss3.par"), overwrite = TRUE)
        r4ss::run(dir_b0, exe = exe, extras = extras, skipfinished = FALSE, verbose = FALSE)
    }

#_____________________________________________________________________________
# 01: single sex
    inputs = read_step("00-s6-base")
    dat = inputs$dat
    ctl = inputs$ctl

    # Nsexes = -1: one sex, spawning biomass multiplied by the female fraction
    dat$Nsexes = -1L
    ctl$Nsexes = -1L

    # Sex x length compositions (sex code 3: females then males) summed over
    # sex into the female columns and recoded as combined sex
    lc = as.data.table(dat$lencomp)
    f_cols = grep("^f[0-9]", names(lc), value = TRUE)
    m_cols = sub("^f", "m", f_cols)
    for (j in seq_along(f_cols)) {
        set(lc, j = f_cols[j], value = lc[[f_cols[j]]] + lc[[m_cols[j]]])
    }
    lc[, (m_cols) := NULL]
    lc[, sex := 0L]
    dat$lencomp = as.data.frame(lc)

    # The two-sex data file carries female and male vectors for the
    # size-frequency observation too (sex 0, so the male half is zero)
    dat$sizefreq_data_list = lapply(dat$sizefreq_data_list, function(sf) {
        sf = as.data.table(sf)
        f_cols = grep("^f[0-9]", names(sf), value = TRUE)
        m_cols = sub("^f", "m", f_cols)
        for (j in seq_along(f_cols)) {
            set(sf, j = f_cols[j], value = sf[[f_cols[j]]] + sf[[m_cols[j]]])
        }
        sf[, (m_cols) := NULL]
        sf[, sex := 0L]
        as.data.frame(sf)
    })

    # Female M-at-age and biology; drop the male growth and weight parameters
    ctl$natM = ctl$natM[1, , drop = FALSE]
    ctl$MG_parms = drop_rows(ctl$MG_parms, "_Mal_GP_1")

    # Selectivity: no sex-specific offsets. In S6 the double-normal base
    # parameters describe males and females are offsets from them (Male
    # option 4); F18's spline has a male dogleg (option 2).
    ctl$size_selex_types$Male = 0L
    ctl$size_selex_parms = drop_rows(ctl$size_selex_parms, "PFemOff|PMalOff")
    ctl$size_selex_parms_tv = drop_rows(ctl$size_selex_parms_tv, "PFemOff|PMalOff")

    inputs$dat = dat
    inputs$ctl = ctl
    run_step(inputs, "01-single-sex")

#_____________________________________________________________________________
# 02a: F18 cubic spline -> double normal, started from F17 (USA longline,
# deep sets), the other USA longline fleet
    inputs = read_step("01-single-sex")
    ctl = inputs$ctl
    sp = ctl$size_selex_parms
    f17 = sp[grep("_F17_", rownames(sp)), ]
    f18 = f17
    rownames(f18) = sub("F17_USA_Lonline_DP\\(17\\)", "F18_USA_Lonline_SH(18)", rownames(f17))
    i18 = grep("_F18_", rownames(sp))
    sp = rbind(sp[seq_len(min(i18) - 1L), ], f18, sp[(max(i18) + 1L):nrow(sp), ])
    ctl$size_selex_parms = sp
    ctl$size_selex_types["F18_USA_Lonline_SH", c("Pattern", "Special")] = c(24L, 0L)
    inputs$ctl = ctl
    run_step(inputs, "02a-f18-dblnorm")

#_____________________________________________________________________________
# 02b: split F8 at its selectivity blocks. F8 has catch only in 1973-1992,
# all inside the two blocks, so F8 keeps 1973-1981 (its first block's peak)
# and a new fishing fleet 32 takes 1982-1992 (the second block's peak).
# Appending the fleet keeps every other fleet number unchanged.
    inputs = read_step("02a-f18-dblnorm")
    dat = inputs$dat
    ctl = inputs$ctl
    new_name = "F8_JPN_LG_MESH_EARLY_1982"
    new_f = dat$Nfleets + 1L
    tv = ctl$size_selex_parms_tv
    peak_73 = tv[grep("F8_.*BLK2repl_1973", rownames(tv)), "INIT"]
    peak_82 = tv[grep("F8_.*BLK2repl_1982", rownames(tv)), "INIT"]
    stopifnot(length(peak_73) == 1L, length(peak_82) == 1L)

    # data: fleet tables
    dat$Nfleets = new_f
    dat$fleetnames = c(dat$fleetnames, new_name)
    dat$fleetinfo = rbind(dat$fleetinfo, dat$fleetinfo[8, ])
    dat$fleetinfo$fleetname[new_f] = new_name
    rownames(dat$fleetinfo)[new_f] = new_name
    dat$CPUEinfo = rbind(dat$CPUEinfo, dat$CPUEinfo[8, ])
    dat$CPUEinfo$fleet[new_f] = new_f
    rownames(dat$CPUEinfo)[new_f] = new_name
    dat$len_info = rbind(dat$len_info, dat$len_info[8, ])
    rownames(dat$len_info)[new_f] = new_name
    # data: observations from 1982 on move to the new fleet
    dat$catch$fleet[dat$catch$fleet == 8 & dat$catch$year >= 1982] = new_f
    # the initial-equilibrium catch row (year -999) is zero for F8; give the
    # new fleet one as well, as SS3 expects one per fishing fleet
    eq = dat$catch[dat$catch$fleet == 8 & dat$catch$year == -999, ]
    if (nrow(eq)) { eq$fleet = new_f; dat$catch = rbind(dat$catch, eq) }
    dat$catch = dat$catch[order(dat$catch$fleet, dat$catch$year), ]
    dat$lencomp$fleet[dat$lencomp$fleet == 8 & dat$lencomp$year >= 1982] = new_f

    # control: selectivity types, parameters, weighting
    ctl$Nfleets = new_f
    ctl$fleetnames = dat$fleetnames
    ctl$size_selex_types = rbind(ctl$size_selex_types, ctl$size_selex_types[8, ])
    rownames(ctl$size_selex_types)[new_f] = new_name
    ctl$age_selex_types = rbind(ctl$age_selex_types, ctl$age_selex_types[8, ])
    rownames(ctl$age_selex_types)[new_f] = new_name
    sp = ctl$size_selex_parms
    i8 = grep("_F8_", rownames(sp))
    sp[i8, c("Block", "Block_Fxn")] = 0
    p8 = sp[i8, ]
    sp[i8[1], "INIT"] = peak_73
    p8[1, "INIT"] = peak_82
    # the block replacements were estimated (phase 4); estimate both peaks
    sp[i8[1], "PHASE"] = p8[1, "PHASE"] = abs(tv[grep("F8_.*BLK2repl_1973", rownames(tv)), "PHASE"])
    rownames(p8) = sub("F8_JPN_LG_MESH_EARLY\\(8\\)", paste0(new_name, "(", new_f, ")"), rownames(p8))
    ctl$size_selex_parms = rbind(sp, p8)
    ctl$size_selex_parms_tv = drop_rows(tv, "F8_")
    if (!nrow(ctl$size_selex_parms_tv)) ctl$size_selex_parms_tv = NULL
    va = ctl$Variance_adjustment_list
    va = rbind(va, transform(va[va$fleet == 8, ], fleet = new_f))
    ctl$Variance_adjustment_list = va[order(va$factor, va$fleet), ]
    lam = ctl$lambdas
    lam = rbind(lam, transform(lam[lam$fleet == 8, ], fleet = new_f))
    ctl$lambdas = lam[order(lam$like_comp, lam$fleet), ]
    ctl$N_lambdas = nrow(ctl$lambdas)

    inputs$dat = dat
    inputs$ctl = ctl
    run_step(inputs, "02b-f8-split")

#_____________________________________________________________________________
# 02c: no annual deviations on F19's peak
    inputs = read_step("02b-f8-split")
    ctl = inputs$ctl
    i19 = grep("SizeSel_P_1_F19_", rownames(ctl$size_selex_parms))
    ctl$size_selex_parms[i19, c("dev_link", "dev_minyr", "dev_maxyr")] = 0
    ctl$size_selex_parms_tv = drop_rows(ctl$size_selex_parms_tv, "F19_")
    if (!nrow(ctl$size_selex_parms_tv)) ctl$size_selex_parms_tv = NULL
    inputs$ctl = ctl
    run_step(inputs, "02c-no-f19-devs")

#_____________________________________________________________________________
# 02d: no regime offset on the initial equilibrium recruitment
    inputs = read_step("02c-no-f19-devs")
    ctl = inputs$ctl
    ctl$SR_parms["SR_regime", c("INIT", "Block", "Block_Fxn")] = 0
    ctl$SR_parms_tv = NULL
    inputs$ctl = ctl
    run_step(inputs, "02d-no-regime")

#_____________________________________________________________________________
# 02e: no F-ballpark penalty (a negative year switches it off)
    inputs = read_step("02d-no-regime")
    inputs$ctl$F_ballpark_year = -abs(inputs$ctl$F_ballpark_year)
    run_step(inputs, "02e-no-f-ballpark")

#_____________________________________________________________________________
# Summarise the bridge
    steps = c("00-s6-base", "01-single-sex", "02a-f18-dblnorm", "02b-f8-split",
              "02c-no-f19-devs", "02d-no-regime", "02e-no-f-ballpark")
    labels = c("S6", "B1 single sex", "B2a F18 double normal", "B2b F8 split",
               "B2c no F19 devs", "B2d no regime", "B2e no F ballpark")
    reps = lapply(steps, function(s)
        SS_output(file.path(dir_ss3, s), covar = TRUE, verbose = FALSE, printstats = FALSE))
    names(reps) = steps

    dq = function(rep, label) {
        x = rep$derived_quants
        v = x$Value[match(label, x$Label)]
        if (length(v) == 0L || is.na(v)) NA_real_ else v
    }
    summary_dt = rbindlist(lapply(seq_along(steps), function(i) {
        rep = reps[[i]]
        lk = rep$likelihoods_used
        data.table(
            step          = steps[i],
            label         = labels[i],
            n_par         = rep$N_estimated_parameters,
            nll           = lk["TOTAL", "values"],
            nll_survey    = lk["Survey", "values"],
            nll_length    = lk["Length_comp", "values"],
            max_gradient  = rep$maximum_gradient_component,
            hessian_ok    = !is.null(rep$CoVar) && isTRUE(rep$log_det_hessian > -Inf),
            ssb0          = dq(rep, "SSB_Virgin"),
            ssb_2020      = dq(rep, "SSB_2020"),
            depletion     = dq(rep, "Bratio_2020"),
            # F_2020 is F / F_MSY (starter F_report_basis 2; F = sum of full Fs)
            f_2020        = dq(rep, "F_2020"),
            f_msy         = dq(rep, "annF_MSY"),
            ssb_msy       = dq(rep, "SSB_MSY"),
            n_on_bounds   = sum(grepl("LO|HI", rep$parameters$Status))
        )
    }))
    print(summary_dt)

    ts_dt = rbindlist(lapply(seq_along(steps), function(i) {
        ts = as.data.table(reps[[i]]$timeseries)[Era == "TIME"]
        data.table(step = steps[i], label = labels[i], year = ts$Yr,
                   ssb = ts$SpawnBio, recruits = ts$Recruit_0)
    }))
    ts_dt[, depletion := ssb / summary_dt$ssb0[match(step, summary_dt$step)]]

    saveRDS(list(summary = summary_dt, timeseries = ts_dt, steps = steps, labels = labels,
                 ss3_version = "3.30.22.1", made_at = format(Sys.time(), tz = "UTC", usetz = TRUE)),
            file.path(dir_cache, "npo-bsh-ss3-bridge-cache.rds"))
