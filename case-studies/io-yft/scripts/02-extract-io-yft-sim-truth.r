# Nicholas Ducharme-Barth
# 2026/10/13
# Extract the simulated truth for the IO yellowfin operating model
# (replicate 4) and cache it for the SS3-to-opal vignette.
#
# The truth is a spatially explicit SPM population model on 159 5x5 degree
# cells; its outputs here are reported by the 4 IOTC assessment regions. The
# SS3 and opal models are single-region, so quantities are summed over those
# regions. Years 1-256 of the
# operating model are quarters and align one-to-one with SS3 years 1-256.
# See case-studies/io-yft/model-files/sim-truth/io-yft-iattc-sample-tree-4/README.md for the
# definitions of SSB and B0.
#
# OUTPUT:
#   case-studies/io-yft/vignettes/data/io-yft-sim-truth-cache.rds
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

#_____________________________________________________________________________
# Packages
    library(data.table)

#_____________________________________________________________________________
# Project paths
    proj_dir  = this.path::this.proj()
    dir_truth = file.path(proj_dir, "case-studies", "io-yft", "model-files", "sim-truth", "io-yft-iattc-sample-tree-4")
    dir_cache = file.path(proj_dir, "case-studies", "io-yft", "vignettes", "data")
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)

#_____________________________________________________________________________
# Helpers
    read_by_region = function(file, value) {
        # The first column is R's row names, which fread also labels V1.
        x = fread(file.path(dir_truth, file))[, -1]
        stopifnot(all(grepl("^V[0-9]+$", names(x))))
        x = x[, which(colSums(!is.na(x)) > 0), with = FALSE]  # drop empty regions
        x[, ts := .I]
        long = melt(x, id.vars = "ts", variable.name = "region", value.name = value)
        long[, region := as.integer(sub("^V", "", region))]
        long[]
    }

    # SPM report lines are "<label>: v1 v2 ..."; return the numeric vectors.
    spm_lines = function(lines, label) {
        hit = grep(paste0("^", label), lines, value = TRUE)
        lapply(hit, function(l) as.numeric(strsplit(trimws(sub("^[^:]*:", "", l)), "\\s+")[[1]]))
    }

#_____________________________________________________________________________
# By region
    ssb   = read_by_region("ssb_yr_reg_ESS5_run4.csv", "ssb")
    bio   = read_by_region("biomass_yr_reg_ESS5_run4.csv", "biomass")
    vuln  = read_by_region("vuln_bio_reg_ESS5_run4.csv", "vuln_biomass")
    recr  = read_by_region("recruitment_yr_reg_ESS5_run4.csv", "recruitment")
    by_region = Reduce(function(a, b) merge(a, b, by = c("ts", "region")),
                       list(ssb, bio, vuln, recr))

#_____________________________________________________________________________
# B0 from the SPM SSB report: SSB at the end of initialisation phase 2, the
# phase the recruitment process references (b0 phase2). No fishing occurs in
# the initialisation time steps.
    ssb_report = readLines(file.path(dir_truth, "SSB_all.txt"))
    init       = spm_lines(ssb_report, "initialisation_phase")
    values     = spm_lines(ssb_report, "values")[[1]]
    stopifnot("SSB_all.txt should hold two initialisation phases" = length(init) == 2L)
    b0 = utils::tail(init[[2]], 1)

#_____________________________________________________________________________
# Totals over regions
    total = by_region[, .(ssb = sum(ssb), biomass = sum(biomass),
                          vuln_biomass = sum(vuln_biomass),
                          recruitment = sum(recruitment)), by = ts]
    setorder(total, ts)
    total[, static_dep := ssb / b0]

    # The region sums must reproduce the operating model's own SSB report,
    # which is printed to six significant figures.
    check = max(abs(total$ssb / values - 1))
    stopifnot("Region SSB does not sum to the SPM SSB report" = check < 1e-5)

#_____________________________________________________________________________
# Cache
    truth = list(
        created    = format(Sys.time(), tz = "UTC", usetz = TRUE),
        source     = "IO yellowfin operating model (SPM), replicate 4; Spatial Stock Assessment Simulation Experiment",
        url        = "https://aaronmberger-nwfsc.github.io/Spatial-Assessment-Modeling-Workshop/",
        files_md5  = tools::md5sum(list.files(dir_truth, full.names = TRUE)) |>
                         stats::setNames(list.files(dir_truth)),
        b0         = b0,
        b0_note    = "SSB at the end of initialisation phase 2 (no fishing); tonnes",
        n_region   = uniqueN(by_region$region),
        total      = total,
        by_region  = by_region,
        init_ssb   = init,
        report_check_max_rel = check
    )
    saveRDS(truth, file.path(dir_cache, "io-yft-sim-truth-cache.rds"))

    cat(sprintf("Truth cache written: B0 = %s t, terminal SSB = %s t, terminal SSB/B0 = %.3f\n",
                format(b0, big.mark = ","), format(round(total[.N, ssb]), big.mark = ","),
                total[.N, static_dep]))
