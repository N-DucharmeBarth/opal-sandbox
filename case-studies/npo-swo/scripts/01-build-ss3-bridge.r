# Nicholas Ducharme-Barth
# 2026/11/13
# Build and run the SS3 bridge from the 2023 ISC North Pacific swordfish
# base case (SS3 3.30.20, two sexes, four seasons) to a single-sex, annual SS3
# model whose structure opal can represent.
#
# Each step changes one thing, starting from the previous step's estimates
# (control.ss_new, plus the previous step's recruitment deviations read
# through the control file's read_recdevs option: the control file alone
# starts the deviations at zero, and from there SS3 finds a worse minimum of
# the base case, 332.18 against 289.30):
#   00-base             as supplied (SS3 3.30.20), re-run with 3.30.22.1 from
#                       the supplied par file
#   01-single-sex       Nsexes = -1 (SSB x female fraction 0.5), female
#                       growth, M and weight-length; the sex-combined
#                       compositions keep one copy of their vector
#   02-block-split      selectivity blocks replaced by fleets: 1994-97 for F1
#                       and the fleets that mirror it (F12, S2), 2004-09 for F2
#                       and the fleets that mirror it (F15, F16, F19). S2's
#                       1994-97 index has its own floated q (see 02 below)
#   03-spawn-month-1    spawning at the start of the year (season 1) instead
#                       of month 5 (season 2): SS3 does not allow spawning
#                       after the start of the season in which recruits
#                       settle, so a one-season model with settlement in
#                       month 7 needs spawning in month 1
#   04-annual           one season: catches summed over quarters; months of
#                       the observations and settlement (7) kept
#   05-annual-comps     the quarterly compositions of each fleet and year
#                       combined (Nsamp summed, proportions weighted by Nsamp)
#   06-settle-month-1   recruits settle at the start of the year (age 0)
#   07-index-month-1    indices observed at the start of the year
# A sensitivity branches from 00:
#   01s-sex-averaged    as 01, with M and growth averaged over the sexes
#                       (weighted by unfished numbers-at-age)
#
# Runs are skipped when the step's Report.sso exists, unless
# NPO_SWO_RERUN=true; NPO_SWO_NOHESS=true skips the Hessian.
#
# OUTPUTS:
#   case-studies/npo-swo/model-files/ss3/<step>/   SS3 inputs and outputs
#   case-studies/npo-swo/vignettes/data/npo-swo-ss3-bridge-cache.rds
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
    dir_ss3   = file.path(proj_dir, "case-studies", "npo-swo", "model-files", "ss3")
    dir_cache = file.path(proj_dir, "case-studies", "npo-swo", "vignettes", "data")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
    extras = if (identical(Sys.getenv("NPO_SWO_NOHESS"), "true")) "-nohess" else ""
    rerun  = identical(Sys.getenv("NPO_SWO_RERUN"), "true")
    exe = file.path(proj_dir, "executables", "stock-synthesis", "3.30.22.1",
                    if (.Platform$OS.type == "windows") "ss3_win.exe" else "ss3_linux")
    if (.Platform$OS.type != "windows") Sys.chmod(exe, "0755")

#_____________________________________________________________________________
# Helpers

    read_rep = function(step) {
        SS_output(file.path(dir_ss3, step), covar = FALSE, verbose = FALSE,
                  printstats = FALSE, forecast = FALSE)
    }

    # Previous step's recruitment deviations (early and main) as the control
    # file's read_recdevs table: starting values for the next estimation
    seed_recdevs = function(ctl, step) {
        rep = read_rep(step)
        p = rep$parameters
        # Early_InitAge_a is the deviation of the cohort born a years before
        # the first model year; Early_RecrDev_y and Main_RecrDev_y are by year
        i_age = grep("^Early_InitAge_", p$Label)
        i_yr  = grep("^(Early|Main)_RecrDev_", p$Label)
        yrs = c(rep$startyr - as.integer(sub(".*_", "", p$Label[i_age])),
                as.integer(sub(".*_", "", p$Label[i_yr])))
        dev = data.frame(Year = yrs, recdev = p$Value[c(i_age, i_yr)])
        ctl$N_Read_recdevs = nrow(dev)
        ctl$recdev_input = dev[order(dev$Year), ]
        ctl
    }

    # SS3 par files: a named list of value lines (names are the "# ..." lines)
    read_par = function(file) {
        x = readLines(file)
        keys = x[seq(2L, length(x), by = 2L)]
        structure(as.list(x[seq(3L, length(x), by = 2L)]), names = keys)
    }
    write_par = function(par, file) {
        writeLines(c("# start values", as.vector(rbind(names(par), unlist(par)))), file)
    }
    max_gradient = function(dir) {
        as.numeric(sub(".*Maximum gradient component = ", "", readLines(file.path(dir, "ss3.par"), n = 1L)))
    }
    ss3 = function(dir, args) {
        r4ss::run(dir, exe = exe, extras = args, skipfinished = FALSE, verbose = FALSE)
        if (!file.exists(file.path(dir, "Report.sso"))) {
            stop(sprintf("SS3 did not finish in %s.", dir), call. = FALSE)
        }
    }

    # Write an SS3 input list to dir and estimate the model there.
    # Start values: `par` (a par list with the new step's layout) when given,
    # else the control file's INIT values (the previous step's estimates)
    # plus the previous step's recruitment deviations. The control file's
    # INIT values are rounded to six digits, which is enough to leave SS3's
    # minimiser short of the optimum (`Step size too small`), so the step is
    # re-run from its own par file until the maximum gradient is below
    # 1e-3 (at most `n_refine` times), and the Hessian is computed last.
    run_step = function(inputs, step, from, par = NULL, n_refine = 4L) {
        dir = file.path(dir_ss3, step)
        if (!rerun && file.exists(file.path(dir, "Report.sso"))) return(invisible(dir))
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        inputs$par = NULL
        inputs$ctl$N_Read_recdevs = 0L
        inputs$ctl$recdev_input = NULL
        if (is.null(par)) {
            inputs$start$init_values_src = 0
            seeded = inputs
            seeded$ctl = seed_recdevs(inputs$ctl, from)
            SS_write(seeded, dir = dir, overwrite = TRUE)
        } else {
            inputs$start$init_values_src = 1
            SS_write(inputs, dir = dir, overwrite = TRUE)
            write_par(par, file.path(dir, "ss3.par"))
        }
        ss3(dir, "-nohess")
        # from here on, start from the step's own par file
        inputs$start$init_values_src = 1
        SS_write(inputs, dir = dir, overwrite = TRUE)
        k = 0L
        while (max_gradient(dir) > 1e-3 && k < n_refine) {
            ss3(dir, "-nohess")
            k = k + 1L
        }
        if (!identical(extras, "-nohess")) ss3(dir, extras)
        invisible(dir)
    }

    read_step = function(step) SS_read(file.path(dir_ss3, step), ss_new = TRUE)

    drop_rows = function(x, pattern) x[!grepl(pattern, rownames(x)), , drop = FALSE]

    # Female (first) half of each composition vector: in the two-sex model the
    # sex-combined (sex 0) observations carry the same vector twice
    single_sex_comps = function(lc) {
        lc = as.data.table(lc)
        f_cols = grep("^f[0-9]", names(lc), value = TRUE)
        m_cols = sub("^f", "m", f_cols)
        stopifnot(all(lc$sex == 0L),
                  isTRUE(all.equal(as.matrix(lc[, ..f_cols]), as.matrix(lc[, ..m_cols]),
                                   check.attributes = FALSE)))
        lc[, (m_cols) := NULL]
        as.data.frame(lc)
    }

    # Von Bertalanffy with L1 at age A1 and L2 at age A2 (SS3 GrowthModel 1)
    vb_l1l2 = function(age, L1, L2, K, A1 = 0.5, A2 = 15) {
        L2 + (L1 - L2) * exp(-K * (age - A1))
    }

#_____________________________________________________________________________
# 00: as supplied, re-run from the supplied par file. Its checksum line holds
# 6.98 instead of 999, which SS3 3.30.22 rejects; the starter reads the par
# (init_values_src = 1).
    dir_b0 = file.path(dir_ss3, "00-base")
    if (rerun || !file.exists(file.path(dir_b0, "Report.sso"))) {
        par = readLines(file.path(dir_b0, "ss-supplied.par"))
        par[length(par)] = "999"
        writeLines(par, file.path(dir_b0, "ss3.par"))
        r4ss::run(dir_b0, exe = exe, extras = extras, skipfinished = FALSE, verbose = FALSE)
    }

#_____________________________________________________________________________
# 01: single sex
    make_single_sex = function(inputs) {
        inputs$dat$Nsexes = -1L
        inputs$ctl$Nsexes = -1L
        inputs$dat$lencomp = single_sex_comps(inputs$dat$lencomp)
        inputs$ctl$natM = inputs$ctl$natM[1, , drop = FALSE]
        inputs$ctl$MG_parms = drop_rows(inputs$ctl$MG_parms, "_Mal_GP_1")
        inputs
    }
    inputs = make_single_sex(read_step("00-base"))
    run_step(inputs, "01-single-sex", from = "00-base")

#_____________________________________________________________________________
# 01s: single sex with sex-averaged biology (sensitivity). M-at-age is the
# rate that takes the unfished numbers of both sexes from one age to the
# next; L1, L2 and K are fitted to the mean of the two sexes' length-at-age
# weighted by their unfished numbers. Weight-length is the same for both
# sexes; maturity and fecundity stay female.
    inputs = read_step("00-base")
    ctl = inputs$ctl
    m_f = unlist(ctl$natM[1, ])
    m_m = unlist(ctl$natM[2, ])
    n_f = exp(-c(0, cumsum(m_f[-length(m_f)])))
    n_m = exp(-c(0, cumsum(m_m[-length(m_m)])))
    n_all = n_f + n_m
    m_avg = -log(n_all[-1] / n_all[-length(n_all)])
    # plus group: the two sexes' M weighted by their numbers at age 15
    a_max = length(m_f)
    m_avg = c(m_avg, (n_f[a_max] * m_f[a_max] + n_m[a_max] * m_m[a_max]) / n_all[a_max])
    mg = ctl$MG_parms
    ages = seq(0.5, 15, by = 0.25)
    len_f = vb_l1l2(ages, mg["L_at_Amin_Fem_GP_1", "INIT"], mg["L_at_Amax_Fem_GP_1", "INIT"],
                    mg["VonBert_K_Fem_GP_1", "INIT"])
    len_m = vb_l1l2(ages, mg["L_at_Amin_Mal_GP_1", "INIT"], mg["L_at_Amax_Mal_GP_1", "INIT"],
                    mg["VonBert_K_Mal_GP_1", "INIT"])
    w_f = approx(0:15, n_f, ages, rule = 2)$y
    w_m = approx(0:15, n_m, ages, rule = 2)$y
    len_avg = (w_f * len_f + w_m * len_m) / (w_f + w_m)
    vb_fit = optim(c(80, 220, 0.25), function(p) sum((len_avg - vb_l1l2(ages, p[1], p[2], p[3]))^2))
    inputs = make_single_sex(inputs)
    inputs$ctl$natM[1, ] = round(m_avg, 5)
    inputs$ctl$MG_parms["L_at_Amin_Fem_GP_1", c("INIT", "PRIOR")] = round(vb_fit$par[1], 3)
    inputs$ctl$MG_parms["L_at_Amax_Fem_GP_1", c("INIT", "PRIOR")] = round(vb_fit$par[2], 3)
    inputs$ctl$MG_parms["VonBert_K_Fem_GP_1", c("INIT", "PRIOR")] = round(vb_fit$par[3], 5)
    run_step(inputs, "01s-sex-averaged", from = "00-base")

#_____________________________________________________________________________
# 02: selectivity blocks -> fleets. Fleets are appended (numbers 28-34) so
# every existing fleet keeps its number, and each mirror points to a
# lower-numbered fleet, as SS3 requires:
#   28 F1 1994-97 (block 2's double normal)   29 F12 1994-97 (mirrors 28)
#   30 S2 1994-97 (mirrors 28; q mirrors S2)  31 F2 2004-09 (block 1's logistic)
#   32 F15 2004-09, 33 F16 2004-09, 34 F19 2004-09 (all mirror 31)
# F13 and F14 (mirrors of F2 and F1) have no catch in the block years. S5
# (mirror of F2, lambda 0) is not split. S2 1994-97 gets its own floated q:
# sharing S2's q through SS3's q mirror (link 2, with S2's q estimated)
# reproduces B1's likelihood exactly at B1's estimates, but the gradient of
# the shared q is wrong there (4.3 at B1's optimum) and the minimiser stops
# short ("Step size too small", gradient 5-60); with its own q the step
# converges, and its likelihood is lower than B1's by the extra q's gain.
    inputs = read_step("01-single-sex")
    dat = inputs$dat
    ctl = inputs$ctl
    tv = ctl$size_selex_parms_tv

    splits = data.table(
        src   = c(1L, 12L, 21L, 2L, 15L, 16L, 19L),
        new   = 28:34,
        first = c(1994L, 1994L, 1994L, 2004L, 2004L, 2004L, 2004L),
        last  = c(1997L, 1997L, 1997L, 2009L, 2009L, 2009L, 2009L),
        mirror = c(NA, 28L, 28L, NA, 31L, 31L, 31L))
    splits[, name := paste0(dat$fleetnames[src], "_", first, "_", substr(last, 3, 4))]

    add_row = function(x, src, name) {
        r = x[src, , drop = FALSE]
        rownames(r) = name
        rbind(x, r)
    }
    for (k in seq_len(nrow(splits))) {
        s = splits[k]
        dat$fleetinfo = add_row(dat$fleetinfo, s$src, s$name)
        dat$fleetinfo$fleetname[s$new] = s$name
        dat$CPUEinfo = add_row(dat$CPUEinfo, s$src, s$name)
        dat$CPUEinfo$fleet[s$new] = s$new
        dat$len_info = add_row(dat$len_info, s$src, s$name)
        ctl$size_selex_types = add_row(ctl$size_selex_types, s$src, s$name)
        ctl$age_selex_types = add_row(ctl$age_selex_types, s$src, s$name)
        if (!is.na(s$mirror)) ctl$size_selex_types[s$new, c("Pattern", "Special")] = c(15L, s$mirror)
        in_block = function(year) year >= s$first & year <= s$last
        dat$catch$fleet[dat$catch$fleet == s$src & in_block(dat$catch$year)] = s$new
        dat$lencomp$fleet[dat$lencomp$fleet == s$src & in_block(dat$lencomp$year)] = s$new
        dat$CPUE$index[dat$CPUE$index == s$src & in_block(dat$CPUE$year)] = s$new
        # one initial-equilibrium catch row per season for each fishing fleet
        eq = dat$catch[dat$catch$fleet == s$src & dat$catch$year == -999, ]
        if (nrow(eq)) { eq$fleet = s$new; eq$catch = 0; dat$catch = rbind(dat$catch, eq) }
    }
    dat$Nfleets = ctl$Nfleets = max(splits$new)
    dat$fleetnames = ctl$fleetnames = dat$fleetinfo$fleetname
    dat$catch = dat$catch[order(dat$catch$fleet, dat$catch$year, dat$catch$seas), ]
    dat$lencomp = dat$lencomp[order(dat$lencomp$fleet, dat$lencomp$year, dat$lencomp$month), ]
    dat$CPUE = dat$CPUE[order(dat$CPUE$index, dat$CPUE$year), ]

    # selectivity parameters: blocks removed; F1 and F2 1994-97 / 2004-09 take
    # the block replacement values (appended in fleet order: 28 then 31)
    sp = ctl$size_selex_parms
    sp[, c("Block", "Block_Fxn")] = 0
    new_sel = function(src_tag, src_fleet, new_fleet, name, blk) {
        p = sp[grep(paste0("_", src_tag, "_"), rownames(sp)), , drop = FALSE]
        rep_vals = tv[grep(paste0(src_tag, "_.*", blk), rownames(tv)), "INIT"]
        p$INIT[seq_along(rep_vals)] = rep_vals
        rownames(p) = sub(sprintf("\\(%d\\)$", src_fleet), sprintf("_%s(%d)", blk, new_fleet), rownames(p))
        p
    }
    ctl$size_selex_parms = rbind(sp,
        new_sel("F1", 1L, 28L, splits$name[1], "BLK2repl_1994"),
        new_sel("F2", 2L, 31L, splits$name[4], "BLK1repl_2004"))
    ctl$size_selex_parms_tv = NULL
    ctl$N_Block_Designs = 0L
    ctl$blocks_per_pattern = NULL
    ctl$Block_Design = NULL

    # catchability: S2 1994-97 is its own index with its own floated q
    qo = ctl$Q_options
    q_new = qo["S2_JPN_WCNPO_OSDWCOLL_late_Area1", , drop = FALSE]
    q_new$fleet = 30L
    rownames(q_new) = splits$name[3]
    ctl$Q_options = rbind(qo, q_new)
    qp = ctl$Q_parms
    q_par_new = qp[grep("S2_JPN", rownames(qp)), , drop = FALSE]
    rownames(q_par_new) = sprintf("LnQ_base_%s(30)", splits$name[3])
    ctl$Q_parms = rbind(qp, q_par_new)

    # length composition lambdas for the new F1 and F2 fleets
    lam = ctl$lambdas
    lam = rbind(lam, transform(lam[lam$like_comp == 4 & lam$fleet == 1, ], fleet = 28L),
                     transform(lam[lam$like_comp == 4 & lam$fleet == 2, ], fleet = 31L))
    ctl$lambdas = lam[order(lam$like_comp, lam$fleet), ]
    ctl$N_lambdas = nrow(ctl$lambdas)

    # Start values: B1's par in B2's layout. Q_parm[9] (S2 1994-97)
    # follows the eight indices' q; the block replacement parameters
    # (selparm 35-40 in B1) become F1 1994-97 (35-40: block values, then F1's
    # two fixed end logits) and F2 2004-09 (41-42).
    p1 = read_par(file.path(dir_ss3, "01-single-sex", "ss3.par"))
    sel = function(i) p1[[sprintf("# selparm[%d]:", i)]]
    p2 = p1[seq_len(match("# Q_parm[8]:", names(p1)))]
    p2[["# Q_parm[9]:"]] = p1[["# Q_parm[2]:"]]
    for (i in 1:42) {
        p2[[sprintf("# selparm[%d]:", i)]] = sel(if (i <= 38) i else c(5, 6, 39, 40)[i - 38])
    }
    p2 = c(p2, p1[setdiff(names(p1), names(p2))])
    stopifnot(identical(tail(names(p2), 1), "# checksum999:"))

    inputs$dat = dat
    inputs$ctl = ctl
    run_step(inputs, "02-block-split", from = "01-single-sex", par = p2)

#_____________________________________________________________________________
# 03: spawning in month 1 (season 1)
    inputs = read_step("02-block-split")
    inputs$dat$spawn_month = 1
    run_step(inputs, "03-spawn-month-1", from = "02-block-split",
             par = read_par(file.path(dir_ss3, "02-block-split", "ss3.par")))

#_____________________________________________________________________________
# 04: one season. Catches summed over the quarters (the initial-equilibrium
# catch as well); observation months and settlement month (7) kept.
    inputs = read_step("03-spawn-month-1")
    dat = inputs$dat
    ctl = inputs$ctl
    dat$nseas = ctl$nseas = 1L
    dat$months_per_seas = 12
    ct = as.data.table(dat$catch)
    ct = ct[, .(seas = 1L, catch = sum(catch), catch_se = catch_se[1]), by = .(year, fleet)]
    setcolorder(ct, c("year", "seas", "fleet", "catch", "catch_se"))
    dat$catch = as.data.frame(ct[order(fleet, year)])
    # Start values: 03's estimates, except the initial F and ln(R0): from
    # 03's initial F (1.39) the annual model's initial equilibrium is
    # depleted to a few percent of SSB0 and the stock cannot support the
    # catch; a small grid of starts (initial F 0.1-0.6, ln(R0) +0 or +0.5)
    # found the best minimum from initial F 0.6 and ln(R0) + 0.5.
    p4 = read_par(file.path(dir_ss3, "03-spawn-month-1", "ss3.par"))
    p4[["# init_F[1]:"]] = "0.6"
    p4[["# SR_parm[1]:"]] = as.character(as.numeric(p4[["# SR_parm[1]:"]]) + 0.5)
    inputs$dat = dat
    inputs$ctl = ctl
    run_step(inputs, "04-annual", from = "03-spawn-month-1", par = p4)

#_____________________________________________________________________________
# 05: one composition per fleet and year. In a one-season model every
# observation of a fleet in a year has the same expectation, so combining
# them (Nsamp summed, proportions weighted by Nsamp) leaves the multinomial
# likelihood unchanged apart from add_to_comp.
    inputs = read_step("04-annual")
    lc = as.data.table(inputs$dat$lencomp)
    bins = grep("^[fl][0-9]", names(lc), value = TRUE)
    lc[, tot := rowSums(.SD), .SDcols = bins]
    lc[, (bins) := lapply(.SD, function(x) x / tot * Nsamp), .SDcols = bins]
    agg = lc[, c(list(month = 7L, sex = sex[1], part = part[1], Nsamp = sum(Nsamp)),
                 lapply(.SD, sum)), by = .(year, fleet), .SDcols = bins]
    setcolorder(agg, c("year", "month", "fleet", "sex", "part", "Nsamp", bins))
    inputs$dat$lencomp = as.data.frame(agg[order(fleet, year)])
    run_step(inputs, "05-annual-comps", from = "04-annual",
             par = read_par(file.path(dir_ss3, "04-annual", "ss3.par")))

#_____________________________________________________________________________
# 06-07: timing aligned with opal, one change per step
    inputs = read_step("05-annual-comps")
    inputs$ctl$recr_dist_pattern$month = 1
    run_step(inputs, "06-settle-month-1", from = "05-annual-comps",
             par = read_par(file.path(dir_ss3, "05-annual-comps", "ss3.par")))

    inputs = read_step("06-settle-month-1")
    inputs$dat$CPUE$month = 1
    run_step(inputs, "07-index-month-1", from = "06-settle-month-1",
             par = read_par(file.path(dir_ss3, "06-settle-month-1", "ss3.par")))

#_____________________________________________________________________________
# Summarise the bridge
    steps = c("00-base", "01-single-sex", "01s-sex-averaged", "02-block-split",
              "03-spawn-month-1", "04-annual", "05-annual-comps",
              "06-settle-month-1", "07-index-month-1")
    labels = c("B0 base case", "B1 single sex", "B1s sex-averaged biology",
               "B2 blocks as fleets", "B3 spawning month 1", "B4 annual",
               "B5 annual compositions", "B6 settlement month 1",
               "B7 index month 1")
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
            nll_recruit   = lk["Recruitment", "values"],
            max_gradient  = rep$maximum_gradient_component,
            hessian_ok    = !is.null(rep$CoVar) && isTRUE(rep$log_det_hessian > -Inf),
            ln_r0         = rep$parameters["SR_LN(R0)", "Value"],
            init_f        = rep$parameters$Value[grep("^InitF_", rep$parameters$Label)],
            ssb0          = dq(rep, "SSB_Virgin"),
            ssb_2021      = dq(rep, "SSB_2021"),
            depletion     = dq(rep, "Bratio_2021"),
            # F_2021: average F over ages 1-10 (starter F_report_units 4)
            f_2021        = dq(rep, "F_2021"),
            f_msy         = dq(rep, "annF_MSY"),
            ssb_msy       = dq(rep, "SSB_MSY"),
            n_on_bounds   = sum(grepl("LO|HI", rep$parameters$Status))
        )
    }))
    print(summary_dt)

    # Annual series: SSB at spawning, total age-0 recruitment over the seasons
    ts_dt = rbindlist(lapply(seq_along(steps), function(i) {
        ts = as.data.table(reps[[i]]$timeseries)[Era == "TIME"]
        ts[, .(step = steps[i], label = labels[i],
               ssb = sum(SpawnBio, na.rm = TRUE), recruits = sum(Recruit_0)), by = .(year = Yr)]
    }))
    ts_dt[, depletion := ssb / summary_dt$ssb0[match(step, summary_dt$step)]]

    saveRDS(list(summary = summary_dt, timeseries = ts_dt, steps = steps, labels = labels,
                 ss3_version = "3.30.22.1", made_at = format(Sys.time(), tz = "UTC", usetz = TRUE)),
            file.path(dir_cache, "npo-swo-ss3-bridge-cache.rds"))
