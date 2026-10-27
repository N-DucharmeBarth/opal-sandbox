# Nicholas Ducharme-Barth
# 2026/10/27
# Add SS3's generalised size compositions for the IO yellowfin index (fleet
# 17, sizefreq method 1) and the control file's lambda changes to the
# exported io_yft_reference, for opal's SS3-matched article
# (vignettes/io-yft-ss3-matched.qmd):
#
#   io_yft_reference$ss3$sizefreq  bins, units, lambda (0.51) and a long
#                                  table of observed and SS3-expected
#                                  proportions, one row per observation x bin
#   io_yft_reference$ss3$lambdas   SS3 lambda changes (priors 0, sizefreq 0.51)
#
# 04-export-io-yft-for-opal.r builds the same two elements when it rebuilds
# the whole export; this script adds them to the existing
# opal-export/data/io_yft_reference.rda without refitting C2.
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Packages
    library(data.table)
    library(magrittr)

#_____________________________________________________________________________
# Project paths
    proj_dir            = this.path::this.proj()
    dir_helper_fns_rtmb = file.path(proj_dir, "code", "rtmb", "helper-fns")
    dir_config          = file.path(proj_dir, "code", "rtmb", "config")
    dir_out             = file.path(proj_dir, "code", "rtmb", "opal-export")
    ref_path            = file.path(dir_out, "data", "io_yft_reference.rda")

    invisible(sapply(
        list.files(dir_helper_fns_rtmb, pattern = "\\.[Rr]$", full.names = TRUE),
        source
    ))

    plain = function(x) {
        x = as.data.frame(x, stringsAsFactors = FALSE)
        attr(x, ".internal.selfref") = NULL
        rownames(x) = NULL
        x
    }

#_____________________________________________________________________________
# Read SS3
    cfg = ss3_read_config(file.path(dir_config, "io-yft-sim-iattc-sample-tree-4.yaml"))
    dir_ss3 = file.path(proj_dir, cfg$ss3_dir)
    ss3 = suppressWarnings(ss3_read_model(dir_ss3,
                                          dat_file = cfg$files$dat,
                                          ctl_file = cfg$files$ctl,
                                          par_file = cfg$files$par,
                                          version  = cfg$ss3_version,
                                          read_report = FALSE))

    sizefreq = ss3_sizefreq_reference(ss3, dir_ss3, cfg, dat_file = cfg$files$dat)
    sizefreq$data = plain(sizefreq$data)
    lambdas = plain(ss3_lambda_reference(ss3))

    # Checks against the SS3 files: 127 observations of 14 bins, all from the
    # survey, every row a proportion, and SS3's lambda of 0.51
    sf = as.data.table(sizefreq$data)
    stopifnot(sf[, uniqueN(ts)] == 127L, sf[, .N, by = ts][, all(N == 14L)],
              sf[, abs(sum(obs) - 1) < 1e-6 && abs(sum(expected) - 1) < 1e-4, by = ts][, all(V1)],
              sizefreq$fleet == 17L, all(sf$fishery %in% 17:20),
              isTRUE(all.equal(sizefreq$lambda, 0.51)))

#_____________________________________________________________________________
# Add to the exported reference
    load(ref_path)
    io_yft_reference$ss3$sizefreq = sizefreq
    io_yft_reference$ss3$lambdas  = lambdas
    io_yft_reference$source$sizefreq_script = "code/rtmb/scripts/06-export-io-yft-sizefreq-for-opal.r"
    save(io_yft_reference, file = ref_path, compress = "xz", version = 3L)

    cat(sprintf("Sizefreq method %d: fleet %d, %d observations, bins %s, lambda %.2f\n",
                sizefreq$method, sizefreq$fleet, sf[, uniqueN(ts)],
                paste(range(sizefreq$bins), collapse = "-"), sizefreq$lambda))
    print(lambdas)
