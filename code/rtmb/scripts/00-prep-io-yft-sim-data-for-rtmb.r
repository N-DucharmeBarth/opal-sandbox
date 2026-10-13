# Nicholas Ducharme-Barth
# 2026/10/13
# Turn the IO yellowfin simulated SS3 model (IATTC sample_tree/4) into a
# configured opal assessment object.
#
# The conversion logic lives in code/rtmb/helper-fns/ss3-*.R. Everything that
# is not derivable from the SS3 files lives in the YAML config. This script is
# just the driver: read config, convert, validate, reconcile, build, save.
#
# The deliverable is an unfitted opal_obj written with opal_save(). No
# optimisation happens here - fitting belongs in the vignette, where the
# selectivity route and the parameter map are the thing under discussion.
#
# OUTPUTS (into cfg$out_dir):
#   <prefix>.rds              configured opal_obj (opal_read() to resume)
#   <prefix>_external.rds     same, on SS3 selectivity-at-age (bridging check)
#   <prefix>_lf.rda           long-format length compositions, for plotting
#   <prefix>_selectivity.rda  sel_fa_external + implied length-sel target
#   <prefix>_meta.rda         config, provenance, checksums, diagnostics
#   <prefix>_checks.csv       structural validation table
#   <prefix>_reconcile.csv    SS3 reconciliation table
#   <prefix>_bridge.csv       SS3 vs opal spawning biomass, with replayed dynamics
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Packages
    library(data.table)
    library(magrittr)
    library(opal)
    library(r4ss)
    library(yaml)

#_____________________________________________________________________________
# Project paths
    proj_dir            = this.path::this.proj()
    dir_helper_fns_rtmb = file.path(proj_dir, "code", "rtmb", "helper-fns")
    dir_config          = file.path(proj_dir, "code", "rtmb", "config")

#_____________________________________________________________________________
# Source converter helpers
    invisible(sapply(
        list.files(dir_helper_fns_rtmb, pattern = "\\.[Rr]$", full.names = TRUE),
        source
    ))

#_____________________________________________________________________________
# WORKFLOW OPTIONS
#
# bio_source  "endgrowth" takes length, SD, weight, maturity, Mat*Fecund and M
#             at age straight from SS3's endgrowth table. Biology is then exact
#             and independent of opal's Schnute growth approximation, which
#             matters here because SS3 uses age-specific K that a single-k
#             curve cannot reproduce. Growth is left doing only one job: the
#             PLA that maps length-based selectivity onto ages.
#
#             "control" is the conventional translation. Weight and maturity
#             are length-basis vectors converted through opal's own PLA, so
#             growth error propagates into spawning biomass. Use it when growth
#             is to be estimated, since age-basis biology is fixed data and
#             cannot carry growth uncertainty.
#
# Set to NULL to use whatever the YAML config specifies.
    bio_source = "endgrowth"

#_____________________________________________________________________________
# Configuration
    cfg = ss3_read_config(file.path(dir_config,
                                    "io-yft-sim-iattc-sample-tree-4.yaml"))

    if (!is.null(bio_source)) cfg$bio$source = bio_source

    dir_out = file.path(proj_dir, cfg$out_dir)
    dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    prefix  = paste0(cfg$prefix, if (identical(cfg$bio$source, "control")) "_ctlbio" else "")

    cat(sprintf("\nBiology route: %s\n", cfg$bio$source))

#_____________________________________________________________________________
# Convert
    converted = ss3_to_opal(cfg, dir_root = proj_dir)

#_____________________________________________________________________________
# Length compositions
# prep_lf_data() attaches lf_obs_in, lf_n, lf_year, lf_fishery and the length
# bin boundaries. Variance adjustments come from SS3's mult_by_lencomp_N; the
# fleets to keep come from the config, which is where SS3's zero composition
# lambdas are acted on. The wide table is passed as a plain data.frame because
# everything downstream of here ends up inside an opal_obj.
    if (!is.null(converted$lf)) {
        lf_wide = dcast(converted$lf, fishery + year + month + ts ~ bin,
                        value.var = "value", fill = 0) %>%
                  setorder(fishery, ts) %>%
                  as.data.frame()

        converted$data = prep_lf_data(
            converted$data,
            lf_wide,
            lf_switch         = 1L,
            lf_keep_fisheries = as.integer(cfg$lf$keep_fisheries),
            lf_var_adjust     = converted$meta$lf_var_adjust,
            lf_cap            = cfg$lf$cap,
            lf_addtocomp      = as.numeric(cfg$lf$addtocomp)
        ) %>% ss3_make_portable()
    }

#_____________________________________________________________________________
# Validate and reconcile
# Both run on the plain lists, before object construction, so a failure points
# at the conversion rather than at opal_obj()'s internals.
    checks = validate_opal_inputs(converted, strict = FALSE)
    recon  = ss3_reconcile(converted, tol = 1e-6)

    if (any(!checks$pass)) {
        print(checks[pass == FALSE])
        stop("Structural validation failed; see the table above.", call. = FALSE)
    }

#_____________________________________________________________________________
# Build the opal assessment object
# ss3_opal_obj() resolves priors, carries the SS3 provenance and converter
# checksums into $provenance$metadata, and runs opal_build() so the initial
# objective and gradient are checked now rather than at first fit.
    assessment = ss3_opal_obj(converted, cfg,
                              map     = NULL,                # get_map() default
                              bounds  = NULL,                # get_bounds() default
                              control = list(eval.max = 10000, iter.max = 10000),
                              build   = TRUE)

#_____________________________________________________________________________
# Bridging configuration: SS3 selectivity-at-age injected directly
# This bypasses opal's PLA selectivity, so with every other parameter held at
# its SS3 value the dynamics should reproduce SS3. That separates "the
# conversion is wrong" from "the two models genuinely differ".
#
# Note the opal_update() route: assigning into assessment$data would break the
# object's identity hash and fail validate_opal_obj().
    stopifnot("SS3 age selectivity unavailable - is Report.sso present?" =
                  !is.null(converted$selectivity$sel_fa_external))

    assessment_external = ss3_use_external_selectivity(
        assessment, converted$selectivity$sel_fa_external) %>%
        opal_build()

# The residual that remains is structural: opal removes catch at the start of
# each timestep, SS3 through the Baranov equation, and opal's first age class
# is recruited from the following timestep's spawning biomass. The replayed
# dynamics switch both to SS3's convention on the same converted inputs, so
# a small replica_ss3 residual is the evidence that the conversion is right.
    bridge = ss3_bridge_table(converted, opal_report(assessment_external))
    cat("\nBridging check, spawning biomass relative to SS3:\n")
    print(bridge[, .(series = c("opal", "replay, SS3 timing"),
                     max_abs_rel  = c(max(abs(rel_opal)), max(abs(rel_replica_ss3))),
                     mean_abs_rel = c(mean(abs(rel_opal)), mean(abs(rel_replica_ss3))))])

#_____________________________________________________________________________
# Starting values for the estimated-selectivity configuration
# SS3's selectivity is on age and opal's is on length, so par_sel is seeded by
# least squares against SS3's selectivity projected onto the length bins. The
# per-fishery RMSE says how well each opal curve family can represent it; a
# poor fit here is a signal about the curve family, not about the data.
    seed = ss3_seed_par_sel(
        sel_target_fl = converted$selectivity$sel_target_fl,
        len_mid       = converted$data$len_mid,
        sel_type_f    = converted$data$sel_type_f,
        n_par         = as.integer(cfg$selectivity$n_par)
    )

    cat("\npar_sel seeding RMSE against SS3 implied selectivity-at-length:\n")
    print(data.table(fishery = seq_along(seed$rmse),
                     sel_type = converted$data$sel_type_f,
                     rmse = round(seed$rmse, 4)))

#_____________________________________________________________________________
# Report
    ss3_conversion_report(converted,
                          checks     = checks,
                          recon      = recon,
                          assessment = assessment)

#_____________________________________________________________________________
# What the growth parameters are still doing
# Growth enters opal through two PLA uses: converting length-basis biology to
# age, and mapping length-based selectivity onto ages. The endgrowth route
# removes the first; sel_fa_external removes the second. When both are gone the
# growth parameters affect nothing, and must stay fixed in the map or the
# Hessian is singular. The opal default map already fixes them.
    cat("\nGrowth parameter role:\n")
    print(rbind(
        cbind(configuration = "estimated selectivity", ss3_growth_role(assessment)),
        cbind(configuration = "SS3 external selectivity", ss3_growth_role(assessment_external))
    ), row.names = FALSE)

#_____________________________________________________________________________
# Save
    opal_save(assessment,
              file.path(dir_out, paste0(prefix, ".rds")),
              overwrite = TRUE)
    opal_save(assessment_external,
              file.path(dir_out, paste0(prefix, "_external.rds")),
              overwrite = TRUE)

    assign(paste0(prefix, "_lf"),          converted$lf)
    assign(paste0(prefix, "_selectivity"),
           c(converted$selectivity, list(par_sel_seed = seed)))
    assign(paste0(prefix, "_meta"),        converted$meta)

    save(list = paste0(prefix, "_lf"),
         file = file.path(dir_out, paste0(prefix, "_lf.rda")))
    save(list = paste0(prefix, "_selectivity"),
         file = file.path(dir_out, paste0(prefix, "_selectivity.rda")))
    save(list = paste0(prefix, "_meta"),
         file = file.path(dir_out, paste0(prefix, "_meta.rda")))

    fwrite(checks, file.path(dir_out, paste0(prefix, "_checks.csv")))
    fwrite(recon,  file.path(dir_out, paste0(prefix, "_reconcile.csv")))
    fwrite(bridge, file.path(dir_out, paste0(prefix, "_bridge.csv")))

#_____________________________________________________________________________
# Downstream, in the vignette:
#
#   assessment = opal_read(file.path(dir_out, "ioyft.rds"), rebuild = TRUE)
#
#   # 1. bridging check: SS3 selectivity, nothing estimated
#   bridge = opal_read(file.path(dir_out, "ioyft_external.rds"))
#   report = opal_report(bridge)            # compare spawning_biomass_y to SS3
#
#   # 2. estimate scale and recruitment only
#   aspm = opal_update(assessment, map = <map with rdev_y, log_B0 freed>)
#   aspm = opal_fit(aspm)
#   summary(aspm); aspm$validation$fit
#
#   # 3. free selectivity against the length compositions
#   lf_model = ss3_use_opal_selectivity(assessment,
#                  par_sel = seed$par_sel,
#                  free    = cfg$lf$keep_fisheries,
#                  n_par_f = seed$n_par_f)
#   lf_model = opal_fit(lf_model)
