# Nicholas Ducharme-Barth
# 2026/11/03
# Compare the MFCL bridge steps M0-M3 (scripts/01-build-mfcl-bridge.r) for
# model region 1: likelihood components, spawning (adult) biomass, depletion,
# recruitment and growth, plus the M2 transplant check (M2's par at M1's
# estimates against M1's region 1).
#
# Each step's fitted par is evaluated once more without estimation to print
# its likelihood components (MFCL writes them to standard output only).
#
# OUTPUT: case-studies/swpo-alb/vignettes/data/swpo-alb-mfcl-bridge.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Paths
    library(data.table)
    proj_dir  = this.path::this.proj()
    dir_cs    = file.path(proj_dir, "case-studies", "swpo-alb")
    dir_mfcl  = file.path(dir_cs, "model-files", "mfcl")
    dir_cache = file.path(dir_cs, "vignettes", "data")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
    mfcl_exe  = file.path(proj_dir, "executables", "mfcl", "2.2.7.0", "mfclo64")
    source(file.path(proj_dir, "shared", "mfcl-to-opal", "mfcl-read.R"))

    steps = data.table(
        step  = c("M0", "M1", "M2", "M3"),
        dir   = c("00-diagnostic", "01-no-movement", "02-region-1", "03-fixed-growth"),
        par   = c("09.par", "01.par", "02.par", "03.par"),
        label = c("Diagnostic (2 regions, movement)", "No movement",
                  "Region 1 only", "Growth fixed, no age-at-length")
    )

#_____________________________________________________________________________
# Likelihood components from an evaluation of each fitted par
    eval_components = function(dir, par) {
        old = setwd(file.path(dir_mfcl, dir))
        on.exit(setwd(old))
        log = "compare-eval.log"
        system2(mfcl_exe, c("alb.frq", par, "compare-eval.par", "-switch", "1", "1", "1", "0"),
                stdout = log, stderr = log)
        x = readLines(log)
        pick = function(pattern) {
            v = grep(pattern, x, value = TRUE)
            if (!length(v)) return(NA_real_)
            as.numeric(sub(".*[=a ] *([-0-9.e+]+) *$", "\\1", trimws(tail(v, 1))))
        }
        data.table(length_comp = pick("^Length frequency data"),
                   age_length  = pick("^Total age length data"),
                   cpue        = pick("Survey index term"),
                   penalties   = pick("After call_penalties pen"),
                   total       = pick("^Total func"))
    }
    like = rbindlist(lapply(seq_len(nrow(steps)), function(i)
        cbind(steps[i, .(step)], eval_components(steps$dir[i], steps$par[i]))))
    par_stats = rbindlist(lapply(seq_len(nrow(steps)), function(i) {
        p = mfcl_read_par(file.path(dir_mfcl, steps$dir[i], steps$par[i]))
        data.table(step = steps$step[i],
                   objective = mfcl_num(p[[mfcl_key(p, "# Objective function value")]]$body),
                   max_gradient = mfcl_num(p[[mfcl_key(p, "# Maximum magnitude gradient value")]]$body),
                   n_par = mfcl_num(p[[mfcl_key(p, "# The number of parameters")]]$body))
    }))
    like = merge(steps[, .(step, label)], merge(par_stats, like, by = "step"), by = "step")

#_____________________________________________________________________________
# Region 1 time series
    rep_of = function(i) {
        f = sub("\\.par$", ".par.rep", paste0("plot-", steps$par[i]))
        mfcl_read_rep(file.path(dir_mfcl, steps$dir[i], f))
    }
    reps = lapply(seq_len(nrow(steps)), rep_of)
    names(reps) = steps$step
    yrs = mfcl_rep(reps$M0, "# Year 1")[1, 1] + seq_len(mfcl_rep(reps$M0, "# Number of time periods")[1, 1]) - 1L
    ts = rbindlist(lapply(names(reps), function(s) {
        r = reps[[s]]
        data.table(step = s, year = yrs,
                   sb = mfcl_rep(r, "# Adult biomass")[, 1],
                   sb_f0 = mfcl_rep(r, "# Adult biomass in absence of fishing")[, 1],
                   total_biomass = mfcl_rep(r, "# Total biomass")[, 1],
                   recruitment = mfcl_rep(r, "# Recruitment")[, 1])
    }))
    ts[, depletion := sb / sb_f0]
    growth = rbindlist(lapply(names(reps), function(s) data.table(
        step = s, age = seq_len(mfcl_rep(reps[[s]], "# Number of age classes")[1, 1]),
        length = mfcl_rep(reps[[s]], "# Mean lengths at age")[1, ],
        sd = mfcl_rep(reps[[s]], "# SD of length at age")[1, ],
        M = mfcl_rep(reps[[s]], "# Natural mortality at age")[1, ])))
    sel = rbindlist(lapply(names(reps), function(s) {
        m = mfcl_rep(reps[[s]], "# Selectivity by age class")
        data.table(step = s, row = rep(seq_len(nrow(m)), ncol(m)),
                   age = rep(seq_len(ncol(m)), each = nrow(m)), sel = as.vector(m))
    }))

#_____________________________________________________________________________
# Transplant check: M2's par at M1's estimates (02-region-1/check.par) vs M1
    r_m1 = mfcl_read_rep(file.path(dir_mfcl, "01-no-movement", "plot-01.par.rep"))
    r_ck = mfcl_read_rep(file.path(dir_mfcl, "02-region-1", "plot-check-eval.par.rep"))
    keep = c(1:13, 18, 19)
    pc1 = mfcl_rep_rows(r_m1, "# Predicted catch by fishery")[keep]
    pc2 = mfcl_rep_rows(r_ck, "# Predicted catch by fishery")
    check = data.table(
        quantity = c("recruitment", "numbers at age (1954)", "numbers at age (all years)",
                     "adult biomass", "predicted catch (extraction fisheries)",
                     "selectivity at age", "mean length at age", "natural mortality at age"),
        max_rel_diff = c(
            max(abs(mfcl_rep(r_ck, "# Recruitment")[, 1] / mfcl_rep(r_m1, "# Recruitment")[, 1] - 1)),
            max(abs(mfcl_rep(r_ck, "# Population Number by age", region = 1)[1, ] /
                    mfcl_rep(r_m1, "# Population Number by age", region = 1)[1, ] - 1)),
            max(abs(mfcl_rep(r_ck, "# Population Number by age", region = 1) /
                    mfcl_rep(r_m1, "# Population Number by age", region = 1) - 1)),
            max(abs(mfcl_rep(r_ck, "# Adult biomass")[, 1] / mfcl_rep(r_m1, "# Adult biomass")[, 1] - 1)),
            max(mapply(function(a, b) max(abs(b / a - 1)), pc1[1:13], pc2[1:13])),
            {
                s1 = mfcl_rep(r_m1, "# Selectivity by age class")
                n_rows = c(rep(4, 7), rep(1, 6), 2, rep(1, 6))
                rows = split(seq_len(sum(n_rows)), rep(1:20, n_rows))
                max(abs(mfcl_rep(r_ck, "# Selectivity by age class") - s1[unlist(rows[keep]), ]))
            },
            max(abs(mfcl_rep(r_ck, "# Mean lengths at age") / mfcl_rep(r_m1, "# Mean lengths at age") - 1)),
            max(abs(mfcl_rep(r_ck, "# Natural mortality at age") / mfcl_rep(r_m1, "# Natural mortality at age") - 1)))
    )
    print(like)
    print(check)
    print(ts[year %in% c(min(year), 1990, max(year)), .(step, year, sb = round(sb), depletion = round(depletion, 3))])

    saveRDS(list(steps = steps, likelihood = like, timeseries = ts, growth = growth,
                 selectivity = sel, transplant_check = check),
            file.path(dir_cache, "swpo-alb-mfcl-bridge.rds"))
