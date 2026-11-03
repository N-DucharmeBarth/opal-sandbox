# Nicholas Ducharme-Barth
# 2026/11/03
# Build the MULTIFAN-CL (MFCL) bridge from the 2024 south Pacific albacore
# diagnostic model (two regions) to a single-region model of model region 1
# (the WCPFC Convention Area south of the equator, including the overlap area
# with the EPO).
#
#   M0  00-diagnostic: the supplied diagnostic model (09.par), evaluated with
#       MFCL 2.2.7.0 without estimation; checks the executable reproduces the
#       supplied objective function value.
#   M1  01-no-movement: M0 with the (SEAPODYM-derived, fixed) movement
#       matrices set to the identity, so the two regions are closed
#       populations; warm-started from 09.par and re-estimated with the
#       diagnostic model's phase 9 settings.
#   M2  02-region-1: region 1 only. The EPO fisheries (14-17) and EPO index
#       (20) are dropped, the WCPFC-CA indices 18 and 19 become 14 and 15. The
#       par is built from M1's estimates (region 1's recruitment kept by
#       shifting the total population scaling parameter by log(0.818971), the
#       region 1 recruitment share); growth is still estimated and the
#       conditional age-at-length data kept. Re-estimated with the phase 9
#       settings.
#   M3  03-fixed-growth: M2 with growth fixed at M2's estimates and the
#       conditional age-at-length data dropped (opal has no conditional
#       age-at-length likelihood). Re-estimated with the phase 9 settings.
#
# Each fit takes 15-20 minutes. The fitted par files are kept and the fits are
# not repeated unless SWPO_ALB_REFIT=true.
#
# OUTPUTS (model-files/mfcl/):
#   00-diagnostic/plot-m0-eval.par.rep
#   01-no-movement/start.par, 01.par, plot-01.par.rep
#   02-region-1/alb.frq, alb.ini, template.par, check.par, start.par, 02.par, plot-02.par.rep
#   03-fixed-growth/start.par, 03.par, plot-03.par.rep
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Paths and settings
    proj_dir = this.path::this.proj()
    dir_cs   = file.path(proj_dir, "case-studies", "swpo-alb")
    dir_mfcl = file.path(dir_cs, "model-files", "mfcl")
    mfcl_exe = file.path(proj_dir, "executables", "mfcl", "2.2.7.0", "mfclo64")
    source(file.path(dir_cs, "scripts", "mfcl-edit-helpers.r"))

    refit = identical(tolower(Sys.getenv("SWPO_ALB_REFIT")), "true")
    keep_fish = c(1:13, 18, 19)   # region 1 fisheries, in their new order
    rec_share_r1 = 0.818971       # region 1 recruitment share (alb.ini, from SEAPODYM)

    dir_m0 = file.path(dir_mfcl, "00-diagnostic")
    dir_m1 = file.path(dir_mfcl, "01-no-movement")
    dir_m2 = file.path(dir_mfcl, "02-region-1")
    dir_m3 = file.path(dir_mfcl, "03-fixed-growth")

# Run MFCL in `dir`; `args` after the frq file. Output goes to `log`.
    run_mfcl = function(dir, args, log) {
        old = setwd(dir)
        on.exit(setwd(old))
        status = system2(mfcl_exe, args, stdout = log, stderr = log)
        invisible(status)
    }
    objective = function(par) {
        secs = mfcl_read_par(par)
        mfcl_num(secs[[mfcl_key(secs, "# Objective function value")]]$body)
    }

# Phase 9 of the diagnostic doitall.sh (without the Hessian).
    phase9 = function(dir, start, out) {
        run_mfcl(dir, c("alb.frq", start, out, "-switch", "2", "1", "1", "200", "2", "116", "30"),
                 file.path(dir, "run-a.log"))
        run_mfcl(dir, c("alb.frq", out, out, "-switch", "2", "1", "1", "2000", "2", "116", "300"),
                 file.path(dir, "run-b.log"))
    }

#_____________________________________________________________________________
# M0: evaluate the supplied diagnostic model
    run_mfcl(dir_m0, c("alb.frq", "09.par", "m0-eval.par", "-switch", "1", "1", "1", "0"),
             file.path(dir_m0, "eval.log"))
    stopifnot("MFCL 2.2.7.0 does not reproduce the diagnostic objective" =
                  isTRUE(all.equal(objective(file.path(dir_m0, "m0-eval.par")),
                                   objective(file.path(dir_m0, "09.par")), tolerance = 1e-10)))

#_____________________________________________________________________________
# M1: no movement
    dir.create(dir_m1, showWarnings = FALSE)
    invisible(file.copy(file.path(dir_m0, c("alb.frq", "alb.ini", "alb.age_length", "mfcl.cfg")),
                        dir_m1, overwrite = TRUE))
    mfcl_write_par(par_identity_movement(mfcl_read_par(file.path(dir_m0, "09.par"))),
              file.path(dir_m1, "start.par"))
    if (refit || !file.exists(file.path(dir_m1, "01.par"))) {
        phase9(dir_m1, "start.par", "01.par")
    }

#_____________________________________________________________________________
# M2: region 1 only
    dir.create(dir_m2, showWarnings = FALSE)
    invisible(file.copy(file.path(dir_m0, c("mfcl.cfg", "alb.age_length")), dir_m2, overwrite = TRUE))
    # The age-at-length samples are all from fisheries 3-11, whose numbers do
    # not change.
    frq_map = subset_frq_region(file.path(dir_m0, "alb.frq"), file.path(dir_m2, "alb.frq"),
                                keep_fish = keep_fish, region = 1L)
    subset_ini_region(file.path(dir_m0, "alb.ini"), file.path(dir_m2, "alb.ini"), region = 1L)
    unlink(file.path(dir_m2, "template.par"))
    run_mfcl(dir_m2, c("alb.frq", "alb.ini", "template.par", "-makepar"),
             file.path(dir_m2, "makepar.log"))
    stopifnot("makepar failed" = file.exists(file.path(dir_m2, "template.par")))

    p2 = build_region_par(mfcl_read_par(file.path(dir_m1, "01.par")),
                          mfcl_read_par(file.path(dir_m2, "template.par")),
                          keep_fish = keep_fish, region = 1L, rec_share = rec_share_r1)
    # Transplant check at M1's estimates: with no movement, region 1's
    # dynamics should be unchanged (see 02-compare-mfcl-bridge.r).
    mfcl_write_par(p2, file.path(dir_m2, "check.par"))
    run_mfcl(dir_m2, c("alb.frq", "check.par", "check-eval.par", "-switch", "1", "1", "1", "0"),
             file.path(dir_m2, "check-eval.log"))
    file.copy(file.path(dir_m2, "check.par"), file.path(dir_m2, "start.par"), overwrite = TRUE)
    if (refit || !file.exists(file.path(dir_m2, "02.par"))) {
        phase9(dir_m2, "start.par", "02.par")
    }

#_____________________________________________________________________________
# M3: growth fixed at M2's estimates, no conditional age-at-length
    dir.create(dir_m3, showWarnings = FALSE)
    invisible(file.copy(file.path(dir_m2, c("alb.frq", "alb.ini", "mfcl.cfg")), dir_m3, overwrite = TRUE))
    # Growth fixed (K, L1, L2, the age 2-4 offsets, generic and length-dependent
    # SD, Richards) and conditional age-at-length off. The second copy of the
    # parest flags is MFCL's history and is left as is.
    no_growth = setNames(as.list(rep(0, 8)), c(12, 13, 14, 15, 16, 184, 227, 240))
    p3 = set_flags(mfcl_read_par(file.path(dir_m2, "02.par")), "# The parest_flags", no_growth)
    mfcl_write_par(p3, file.path(dir_m3, "start.par"))
    if (refit || !file.exists(file.path(dir_m3, "03.par"))) {
        phase9(dir_m3, "start.par", "03.par")
    }
