# Nicholas Ducharme-Barth
# 2026/10/14
# Export the IO yellowfin fit C2 (SS3 age-specific K growth, opal's own
# length-based selectivity, prior-dominated selectivity parameters fixed) as
# self-contained inputs for an opal package vignette.
#
# Everything opal needs is written under code/rtmb/opal-export/, laid out as
# the opal repository is, so the folder can be copied across directly:
#
#   data/io_yft_data.rda           opal data list (catch, index, LF, biology)
#   data/io_yft_parameters.rda     starting parameters: the converted SS3
#                                  optimum with seeded selectivity curves
#   data/io_yft_map.rda            fit C2's map: log_B0, 256 rec devs, 4 q and
#                                  49 selectivity parameters estimated
#   data/io_yft_reference.rda      SS3 and simulated-truth comparison tables,
#                                  with SS3's age_selex patterns and parameters
#   inst/extdata/io_yft_c2_fit.rds the fitted C2 opal_obj (opal_save())
#
# The vignette, the data documentation and the opal-side fit script in the
# same folder use only these objects and exported opal functions. SS3, r4ss
# and the converter in code/rtmb/helper-fns stay in this repository.
#
# Needs the caches from 01-fit-io-yft-ss3-opal.r, 02-extract-io-yft-sim-truth.r
# and 03-fit-io-yft-age-k.r (the SS3 reference quantities, the truth and the
# earlier fits are read from them, not recomputed). The C2 fit here must
# reproduce 03's C2 objective; the script stops if it does not.
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
    dir_out             = file.path(proj_dir, "code", "rtmb", "opal-export")
    dir.create(file.path(dir_out, "data"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(dir_out, "inst", "extdata"), recursive = TRUE, showWarnings = FALSE)

    invisible(sapply(
        list.files(dir_helper_fns_rtmb, pattern = "\\.[Rr]$", full.names = TRUE),
        source
    ))

    cfg_path = file.path(dir_config, "io-yft-sim-iattc-sample-tree-4.yaml")

    # plain data.frame without data.table's self-reference, for package data
    plain = function(x) {
        x = as.data.frame(x, stringsAsFactors = FALSE)
        attr(x, ".internal.selfref") = NULL
        rownames(x) = NULL
        x
    }

#_____________________________________________________________________________
# Convert (endgrowth biology, SS3 age-specific K growth) exactly as in 03
    cfg = ss3_read_config(cfg_path)
    cfg$bio$source   = "endgrowth"
    cfg$growth$model = "ss3_age_k"
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

#_____________________________________________________________________________
# C2: opal selectivity seeded against SS3's, the nine prior-dominated
# parameters held at their seeded values
    seed = ss3_seed_par_sel(
        sel_target_fl = converted$selectivity$sel_target_fl,
        len_mid       = converted$data$len_mid,
        sel_type_f    = converted$data$sel_type_f,
        n_par         = as.integer(cfg$selectivity$n_par))
    x = ss3_opal_obj(converted, cfg, control = list(eval.max = 10000, iter.max = 10000),
                     build = FALSE) %>%
        ss3_use_estimated_selectivity(seed$par_sel, cfg,
                                      fix = cfg$selectivity$fix_prior_dominated)

#_____________________________________________________________________________
# The bundled inputs: what opal_obj() receives in the vignette
    io_yft_data       = x$data
    io_yft_parameters = x$parameters
    io_yft_map        = x$map

    # The rebuilt object must start where the converter's did.
    x_bundled = opal_obj(io_yft_data, io_yft_parameters, io_yft_map,
                         control = list(eval.max = 10000, iter.max = 10000)) %>%
                opal_build()
    x = opal_build(x)
    o1 = opal_rtmb(x, fresh = TRUE)
    o2 = opal_rtmb(x_bundled, fresh = TRUE)
    stopifnot("Bundled inputs do not reproduce the converter's starting objective" =
                  isTRUE(all.equal(o1$fn(o1$par), o2$fn(o2$par), tolerance = 1e-10)),
              "Bundled inputs do not estimate the same parameters" =
                  identical(names(o1$par), names(o2$par)))
    n_sel = sum(grepl("^par_sel", names(o2$par)))
    cat(sprintf("Bundled inputs: %d estimated parameters (%d selectivity), objective at start %.3f\n",
                length(o2$par), n_sel, o2$fn(o2$par)))
    rm(o1, x)

#_____________________________________________________________________________
# Fit C2 from the bundled inputs: two nlminb passes, Newton polish, check.
# data-raw/generate-io-yft-c2-fit.R in the export repeats this inside opal.
    t0  = proc.time()[["elapsed"]]
    fit = opal_fit(x_bundled, n_passes = 2L, check = FALSE)
    obj = opal_rtmb(fit, fresh = TRUE)
    pol = ss3_newton_polish(obj, fit$fit$opt$par)
    opt = fit$fit$opt
    opt$par = stats::setNames(pol$par, names(opt$par))
    opt$objective = pol$objective
    fit = opal_attach_fit(fit, opt, check = TRUE,
                          check_args = list(gradient_tolerance = 1e-4))
    fit = opal_update(fit, metadata = list(
        stock  = "Indian Ocean yellowfin (simulated, replicate 4)",
        model  = "C2: SS3 age-specific K growth, opal length-based selectivity",
        source = "SS3 model by M.N. Maunder, H. Xu and C. Minte-Vera (IATTC), converted to opal",
        converter = "https://github.com/N-DucharmeBarth/wcpfc-bet-single-region (code/rtmb)"))
    fit_seconds = proc.time()[["elapsed"]] - t0

    prev = readRDS(file.path(dir_cache, "io-yft-age-k-cache.rds"))
    c2_prev = prev$summary[id == "C2", objective_fit]
    cat(sprintf("C2 objective %.4f (03 cache: %.4f), max |grad| %.2e, check %s, %.0f s\n",
                opt$objective, c2_prev, pol$max_gradient,
                isTRUE(fit$validation$fit$passes), fit_seconds))
    stopifnot("C2 does not reproduce the 03 cache's objective" =
                  abs(opt$objective - c2_prev) < 1e-3,
              "C2 fails opal_check()" = isTRUE(fit$validation$fit$passes))

    fit_path = file.path(dir_out, "inst", "extdata", "io_yft_c2_fit.rds")
    opal_save(fit, fit_path, compress = "xz", overwrite = TRUE)
    stopifnot(isTRUE(all.equal(opal_read(fit_path, rebuild = TRUE)$fit$opt$objective,
                               opt$objective)))

#_____________________________________________________________________________
# Comparison tables: SS3, the simulated truth and the earlier opal fits
    truth = readRDS(file.path(dir_cache, "io-yft-sim-truth-cache.rds"))
    sch   = readRDS(file.path(dir_cache, "io-yft-ss3-opal-cache.rds"))
    ss3   = prev$ss3
    d     = io_yft_data

    # SS3's expected length compositions aggregated per fishery, weighted by
    # the observations' sample sizes, as the vignette aggregates opal's.
    n_obs = rbindlist(lapply(seq_along(d$lf_fishery_f), function(j)
        data.table(fishery = d$lf_fishery_f[j], ts = d$lf_year_fi[[j]], n = d$lf_n_fi[[j]])))
    ss3_lf = merge(ss3$lf, n_obs, by = c("fishery", "ts"))
    ss3_lf = ss3_lf[, .(expected = sum(ss3_pred * n)), by = .(fishery, bin)]
    ss3_lf[, expected := expected / sum(expected), by = fishery]
    setorder(ss3_lf, fishery, bin)

    end_ts = function(dt, v) dt[ts == max(ts)][[v]]
    history = data.table(
        model = c("opal, Schnute growth, SS3 selectivity",
                  "opal, SS3 age-specific K growth, SS3 selectivity",
                  "opal, SS3 age-specific K growth, opal selectivity (C2)"),
        growth = c("schnute", "ss3_age_k", "ss3_age_k"),
        selectivity = c("SS3 (fixed)", "SS3 (fixed)", "opal (estimated)"),
        B0 = c(sch$summary[route == "endgrowth", B0_fit],
               prev$summary[id == "A1", B0_fit],
               prev$summary[id == "C2", B0_fit]),
        terminal_depletion = c(
            end_ts(sch$trajectory[route == "endgrowth" & state == "fitted"], "static_dep"),
            end_ts(prev$trajectory[id == "A1" & state == "fitted"], "static_dep"),
            end_ts(prev$trajectory[id == "C2" & state == "fitted"], "static_dep")))

    # SS3's age_selex patterns and parameters (pattern 20 double normal, 12
    # logistic, 15 mirror), so opal's age-based selectivity (sel_basis_f = 2)
    # can start from them through convert_ss3_selex_to_rtmb(basis = "age").
    # control.ss_new carries the estimates as INIT.
    ctl = converted$ss3$ctl
    sel_types = as.data.frame(ctl$age_selex_types)
    sel_pattern = data.table(fleet   = seq_len(nrow(sel_types)),
                             name    = prev$fleet_names[seq_len(nrow(sel_types))],
                             pattern = as.integer(sel_types[[1]]),
                             special = as.integer(sel_types[[4]]))
    sp = as.data.frame(ctl$age_selex_parms)
    sel_par = data.table(label = rownames(sp), value = sp$INIT, lo = sp$LO, hi = sp$HI,
                         phase = as.integer(sp$PHASE))
    sel_par[, fleet := as.integer(sub("^AgeSel_([0-9]+)P_.*", "\\1", label))]
    sel_par[, par   := as.integer(sub("^AgeSel_[0-9]+P_([0-9]+)_.*", "\\1", label))]
    sel_par = merge(sel_par, sel_pattern[, .(fleet, pattern)], by = "fleet")
    setcolorder(sel_par, c("fleet", "pattern", "par", "value", "phase", "lo", "hi", "label"))
    setorder(sel_par, fleet, par)
    stopifnot("Unexpected SS3 age selectivity parameters" =
                  sel_par[pattern == 20L, .N] == 6L * sel_pattern[pattern == 20L, .N],
              sel_par[pattern == 12L, .N] == 2L * sel_pattern[pattern == 12L, .N])

    # SS3's generalised size compositions for the index (fleet 17), with
    # SS3's expected values and its lambda, for an SS3-matched opal model
    # (the opal article io-yft-ss3-matched). Also written on their own by
    # 06-export-io-yft-sizefreq-for-opal.r.
    sizefreq = ss3_sizefreq_reference(converted$ss3, file.path(proj_dir, cfg$ss3_dir), cfg)
    sizefreq$data = plain(sizefreq$data)

    io_yft_reference = list(
        fleet_names = c(prev$fleet_names[1:16], paste0("llcpue_q", 1:4)),
        ss3 = list(
            version     = ss3$version,
            ssb0        = ss3$ssb0,
            trajectory  = plain(ss3$trajectory[, .(ts, ssb, static_dep)]),
            recdevs     = plain(ss3$recdevs[, .(ts = year, recdev)]),
            cpue        = plain(ss3$cpue[, .(ts, fishery, obs, expected = exp, se)]),
            lf          = plain(ss3_lf),
            sel_age     = plain(ss3$sel_age[, .(fishery, age, sel)]),
            sel_pattern = plain(sel_pattern),
            sel_par     = plain(sel_par),
            sizefreq    = sizefreq,
            lambdas     = plain(ss3_lambda_reference(converted$ss3)),
            length_at_age = plain(data.table(age = io_yft_data$age_a,
                                             length = prev$growth$laa_ss3))),
        truth = list(
            b0         = truth$b0,
            trajectory = plain(truth$total[, .(ts, ssb, static_dep)])),
        schnute_length_at_age = plain(data.table(
            age = io_yft_data$age_a,
            length = sch$biology[route == "endgrowth", length_opal])),
        history = plain(history),
        source = list(
            repository = "https://github.com/N-DucharmeBarth/wcpfc-bet-single-region",
            script     = "code/rtmb/scripts/04-export-io-yft-for-opal.r",
            commit     = tryCatch(system2("git", c("-C", proj_dir, "rev-parse", "--short", "HEAD"),
                                          stdout = TRUE), error = function(e) NA_character_),
            created    = format(Sys.time(), tz = "UTC", usetz = TRUE),
            opal       = as.character(packageVersion("opal")),
            r4ss       = as.character(packageVersion("r4ss")))
    )

#_____________________________________________________________________________
# Write
    for (nm in c("io_yft_data", "io_yft_parameters", "io_yft_map", "io_yft_reference")) {
        save(list = nm, file = file.path(dir_out, "data", paste0(nm, ".rda")),
             compress = "xz", version = 3L)
    }
    sizes = file.size(c(list.files(file.path(dir_out, "data"), full.names = TRUE), fit_path))
    cat(sprintf("\nExport written to %s (%.0f kB in total)\n", dir_out, sum(sizes) / 1024))
    print(data.frame(file = basename(c(list.files(file.path(dir_out, "data")), fit_path)),
                     kB = round(sizes / 1024)))
