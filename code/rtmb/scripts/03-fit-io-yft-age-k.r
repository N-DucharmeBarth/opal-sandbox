# Nicholas Ducharme-Barth
# 2026/10/13
# Fit the IO yellowfin SS3 conversion in opal with SS3's age-specific K growth
# (opal >= 0.1.1) and cache the results for the age-specific K vignette.
#
# Three fits, all with growth fixed at SS3's GrowthModel 3 values read from
# the SS3 files (cfg$growth$model = "ss3_age_k"), starting from the converted
# SS3 optimum:
#
#   A1  endgrowth biology, SS3 selectivity-at-age fixed (sel_fa_external)
#   A2  control biology,   SS3 selectivity-at-age fixed
#   C   endgrowth biology, opal's own length-based selectivity estimated.
#       The opal equivalents of exactly the selectivity parameters SS3
#       estimated are freed (cfg$selectivity: estimate / mirror / free_par),
#       starting from curves seeded against SS3's selectivity.
#
# A1 and A2 estimate log_B0, the recruitment deviations and the quarterly
# catchabilities (opal's default map); C adds the 58 selectivity parameters.
#
# The SS3 reference quantities (index fit, expected compositions,
# trajectories) and the earlier Schnute-growth fit are read from the cache of
# 01-fit-io-yft-ss3-opal.r rather than recomputed.
#
# OUTPUT:
#   code/rtmb/vignettes/data/io-yft-age-k-cache.rds
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

    stopifnot("opal >= 0.1.1 is required for age-specific K growth" =
                  utils::packageVersion("opal") >= "0.1.1")

#_____________________________________________________________________________
# Project paths
    proj_dir            = this.path::this.proj()
    dir_helper_fns_rtmb = file.path(proj_dir, "code", "rtmb", "helper-fns")
    dir_config          = file.path(proj_dir, "code", "rtmb", "config")
    dir_cache           = file.path(proj_dir, "code", "rtmb", "vignettes", "data")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)

    invisible(sapply(
        list.files(dir_helper_fns_rtmb, pattern = "\\.[Rr]$", full.names = TRUE),
        source
    ))

    cfg_path = file.path(dir_config, "io-yft-sim-iattc-sample-tree-4.yaml")

#_____________________________________________________________________________
# Convert, with SS3 growth, and prepare the length compositions
    convert = function(bio_source) {
        cfg = ss3_read_config(cfg_path)
        cfg$bio$source     = bio_source
        cfg$growth$model   = "ss3_age_k"
        converted = ss3_to_opal(cfg, dir_root = proj_dir)
        lf_wide = dcast(converted$lf, fishery + year + month + ts ~ bin,
                        value.var = "value", fill = 0) %>%
                  setorder(fishery, ts) %>%
                  as.data.frame()
        converted$data = prep_lf_data(
            converted$data, lf_wide,
            lf_switch         = 1L,
            lf_keep_fisheries = as.integer(cfg$lf$keep_fisheries),
            lf_var_adjust     = converted$meta$lf_var_adjust,
            lf_cap            = cfg$lf$cap,
            lf_addtocomp      = as.numeric(cfg$lf$addtocomp)
        ) %>% ss3_make_portable()
        validate_opal_inputs(converted, strict = TRUE)
        list(cfg = cfg, converted = converted)
    }

    build = function(cv) {
        ss3_opal_obj(cv$converted, cv$cfg,
                     control = list(eval.max = 10000, iter.max = 10000),
                     build = TRUE)
    }

#_____________________________________________________________________________
# A1 and A2: SS3 selectivity-at-age fixed
    fits = list()
    conv = list()
    for (route in c("endgrowth", "control")) {
        t_conv  = ss3_elapsed(convert(route))
        cv      = t_conv$value
        t_build = ss3_elapsed(build(cv) %>%
                                  ss3_use_external_selectivity(cv$converted$selectivity$sel_fa_external) %>%
                                  opal_build())
        id = c(endgrowth = "A1", control = "A2")[[route]]
        fits[[id]] = ss3_fit_extract(
            t_build$value, id, cv$converted,
            setup_seconds = c("convert SS3 files" = t_conv$seconds,
                              "build opal objects" = t_build$seconds))
        conv[[route]] = cv
    }

#_____________________________________________________________________________
# C: opal's own selectivity, estimated (endgrowth biology)
    cv = conv$endgrowth
    t_seed = ss3_elapsed(ss3_seed_par_sel(
        sel_target_fl = cv$converted$selectivity$sel_target_fl,
        len_mid       = cv$converted$data$len_mid,
        sel_type_f    = cv$converted$data$sel_type_f,
        n_par         = as.integer(cv$cfg$selectivity$n_par)))
    t_build = ss3_elapsed(build(cv) %>%
                              ss3_use_estimated_selectivity(t_seed$value$par_sel, cv$cfg) %>%
                              opal_build())
    fits[["C"]] = ss3_fit_extract(
        t_build$value, "C", cv$converted,
        setup_seconds = c("seed selectivity" = t_seed$seconds,
                          "build opal objects" = t_build$seconds))

#_____________________________________________________________________________
# Reference quantities
    stack = function(what) rbindlist(lapply(fits, `[[`, what), use.names = TRUE, fill = TRUE)

    prev = readRDS(file.path(dir_cache, "io-yft-ss3-opal-cache.rds"))
    ss3  = prev$ss3
    ss3$sel_age = data.table(
        fishery = rep(seq_len(nrow(cv$converted$selectivity$sel_fa_external)),
                      times = ncol(cv$converted$selectivity$sel_fa_external)),
        age     = rep(cv$converted$data$age_a,
                      each = nrow(cv$converted$selectivity$sel_fa_external)),
        sel     = as.vector(cv$converted$selectivity$sel_fa_external))
    ss3$sel_age_factor = cv$converted$selectivity$sel_factor
    ss3$age_selex_parms = {
        a = cv$converted$ss3$ctl$age_selex_parms
        data.table(parameter = rownames(a), value = a$INIT, phase = a$PHASE)
    }

    growth_spec = list(
        k_mult      = cv$converted$meta$growth_fit$k_mult,
        k_mult_ages = cv$converted$meta$growth_fit$k_mult_ages,
        k_per_step  = cv$converted$meta$growth_fit$k,
        L1 = exp(cv$converted$parameters$log_L1),
        L2 = exp(cv$converted$parameters$log_L2),
        CV1 = exp(cv$converted$parameters$log_CV1),
        CV2 = exp(cv$converted$parameters$log_CV2),
        laa_opal = cv$converted$meta$growth_fit$laa_opal,
        laa_ss3  = cv$converted$meta$growth_fit$laa_ss3)

    cache = list(
        created     = format(Sys.time(), tz = "UTC", usetz = TRUE),
        session     = list(r = as.character(getRversion()),
                           opal = as.character(packageVersion("opal")),
                           opal_sha = tryCatch(utils::packageDescription("opal")$RemoteSha %||%
                                                   NA_character_, error = function(e) NA_character_),
                           r4ss = as.character(packageVersion("r4ss")),
                           cpu = tryCatch(trimws(sub(".*: ", "", grep("model name",
                               readLines("/proc/cpuinfo"), value = TRUE)[1])),
                               error = function(e) Sys.info()[["machine"]]),
                           cores = parallel::detectCores()),
        config_md5  = unname(tools::md5sum(cfg_path)),
        fleet_names = cv$converted$meta$dims$.fleet_names,
        growth      = growth_spec,
        selectivity_config = cv$cfg$selectivity,
        seed_rmse   = t_seed$value$rmse,
        timing      = stack("timing"),
        optimiser   = stack("optimiser"),
        summary     = stack("summary"),
        parameters  = stack("parameters"),
        q           = stack("q"),
        par_sel     = stack("par_sel"),
        trajectory  = stack("trajectory"),
        cpue        = stack("cpue"),
        lf          = stack("lf"),
        sel_age     = stack("sel_age"),
        biology     = stack("biology"),
        reconcile   = stack("reconcile"),
        # earlier fits with opal's Schnute growth, for reference lines
        schnute     = list(trajectory = prev$trajectory[route %in% c("endgrowth", "control") &
                                                            state == "fitted"],
                           summary    = prev$summary[route %in% c("endgrowth", "control")]),
        ss3         = ss3
    )

    out = file.path(dir_cache, "io-yft-age-k-cache.rds")
    saveRDS(cache, out, compress = "xz")
    cat(sprintf("\nCache written: %s (%.0f kB)\n", out, file.size(out) / 1024))
    print(cache$summary[, .(id, n_par, n_par_sel, objective_start, objective_fit, max_gradient,
                            check_passes, B0_fit, R0_fit)])
    print(dcast(cache$timing, stage ~ id, value.var = "seconds"))
