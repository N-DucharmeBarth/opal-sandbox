# Nicholas Ducharme-Barth
# 2026/10/15
# SS3 -> opal, end to end, for the IO yellowfin SS3 model: read the SS3
# files, convert, test the conversion, then set up, fit and save the final
# opal model (fit C2 of 03-fit-io-yft-age-k.r).
#
# This is the companion script to the guide
# case-studies/io-yft/vignettes/ss3-to-opal-workflow.qmd, which shows its code chunk by
# chunk (the "## ---- label ----" lines) and reads the cache it writes. It is
# meant to be copied and edited for another SS3 model: the conversion settings
# are an ordinary R list below, and every step is a plain function call.
#
# Needs opal >= 0.1.1 (SS3 age-specific K growth), r4ss and data.table, and
# the converter functions in shared/ss3-to-opal. Runs in about 3-4 minutes,
# almost all of it the fit.
#
# OUTPUTS:
#   case-studies/io-yft/vignettes/data/io-yft-workflow-cache.rds  tables for the guide
#   case-studies/io-yft/vignettes/data/io-yft-workflow-fit.rds    the fitted opal_obj
#
# Copyright (c) 2026 Nicholas Ducharme-Barth
# GPL-3.0

## ---- packages ----
library(data.table)
library(magrittr)
library(opal)
library(r4ss)

stopifnot("opal >= 0.1.1 is required for SS3 age-specific K growth" =
              utils::packageVersion("opal") >= "0.1.1")

## ---- paths ----
proj_dir = this.path::this.proj()
dir_fns  = file.path(proj_dir, "shared", "ss3-to-opal")
dir_out  = file.path(proj_dir, "case-studies", "io-yft", "vignettes", "data")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# The converter: ss3-read.R, ss3-convert.R, ss3-biology.R, ss3-validate.R,
# ss3-opal-obj.R and ss3-fit-tools.R.
for (f in list.files(dir_fns, pattern = "\\.[Rr]$", full.names = TRUE)) source(f)

## ---- settings ----
# Every choice that cannot be read from the SS3 files. The converter reads
# nothing else, so this list is the full record of the conversion decisions.
cfg = list(
    model       = "io-yft-sim-iattc-sample-tree-4",
    # SS3 model folder, relative to dir_root in ss3_to_opal() (or set
    # absolute_paths = TRUE and give a full path).
    ss3_dir     = "case-studies/io-yft/model-files/ss3/xx-io-yft-sim-iattc-sample-tree-4",
    ss3_version = "3.24",
    files = list(
        dat         = "data.ss_new",
        ctl         = "control.ss_new",
        par         = "ss3.par",
        read_report = TRUE        # Report.sso: endgrowth, ageselex, SSB0, bias ramp
    ),

    # Time: 256 SS3 "years" that are quarters (nseas = 1, months_per_seas = 3).
    # opal runs the same 256 steps with n_season = 1.
    time = list(first_ts = 1, last_ts = 256),

    # Ages: drop SS3 age 0 so that the opal age index equals the SS3 age.
    age = list(keep_ss3_ages = c(1, 28), min_age = 1),

    # Growth: SS3 GrowthModel 3 (age-specific K) reproduced exactly in opal.
    # "schnute" would instead fit opal's single-k curve to SS3's length-at-age.
    growth = list(
        model             = "ss3_age_k",
        plus_group_adjust = TRUE,   # SS3's plus-group mean length
        A1 = 1, A2 = 28, fit_cv = TRUE   # used only by model = "schnute"
    ),

    # Biology: M, weight, maturity and Mat*Fecund at age straight from SS3's
    # endgrowth table ("control" would rebuild them from the control file).
    bio = list(
        source = "endgrowth",
        endgrowth = list(
            length_basis       = "Beg",
            weight             = TRUE,
            maturity           = TRUE,
            spawning_potential = TRUE,
            apply_frac_female  = FALSE,   # single-sex model: SSB has no female fraction
            m                  = TRUE
        ),
        m_rate_divisor = NULL   # NULL: 12 / (nseas * months_per_seas) = 4 here
    ),

    # Index: SS3 fleet 17 (longline CPUE) split into four quarterly series,
    # opal fisheries 17-20, each with its own catchability.
    index = list(
        keep_ss3_fleets = 17,
        split = list(list(fleet = 17, into = 17:20, by = "ts_modulo",
                          modulo = 4, offset = 81)),
        normalise            = "mean1",
        month_labels         = c(2, 5, 8, 11),
        keep_negative_fleets = FALSE
    ),

    # Length compositions: fleets with length data in SS3; the last SS3 bin
    # (195-198 cm) is regularised onto the 5 cm grid.
    lf = list(
        keep_negative_fleets = FALSE,
        regularise_last_bin  = TRUE,
        keep_fisheries       = c(1:14, 16),
        cap                  = NULL,
        addtocomp            = 1e-10
    ),

    # Selectivity. opal codes (not SS3 patterns): 1 logistic, 2 double-normal,
    # 3 double Richards. Fisheries 1-16 double-normal, index 17-20 logistic.
    selectivity = list(
        report_factor  = "Asel2",   # Report.sso ageselex rows to read
        reference_year = NULL,      # NULL: terminal year
        n_par          = 6,
        opal_type_f    = c(rep(2, 16), rep(1, 4)),
        # Estimate the opal equivalents of the parameters SS3 estimated.
        estimate = c(1:14, 17),
        mirror   = list(`15` = 14, `16` = 14, `18` = 17, `19` = 17, `20` = 17),
        free_par = list(`1` = 1:2, `2` = 1:4),
        # Parameters the length data do not inform (they settled on the prior
        # mean when estimated), held at their seeded values instead.
        fix_prior_dominated = list(`4` = c(2, 4), `5` = c(2, 4), `7` = c(2, 4),
                                   `9` = 4, `13` = c(2, 4))
    ),

    # SS3's generalised size compositions are not converted.
    sizefreq = list(include = FALSE),

    metadata = list(species = "Yellowfin tuna", region = "Indian Ocean",
                    data_type = "simulated")
)
# Provenance stored on the opal object.
cfg$config_path = "case-studies/io-yft/scripts/05-io-yft-ss3-to-opal-workflow.r"
cfg$config_md5  = unname(tools::md5sum(file.path(proj_dir, cfg$config_path)))

## ---- read-ss3 ----
# ss3_to_opal() reads the files itself; reading them here as well is for
# looking at the model before converting it.
ss3 = ss3_read_model(dir_ss3     = file.path(proj_dir, cfg$ss3_dir),
                     dat_file    = cfg$files$dat,
                     ctl_file    = cfg$files$ctl,
                     par_file    = cfg$files$par,
                     version     = cfg$ss3_version,
                     read_report = TRUE)
names(ss3)   # dat, ctl, par, rep (r4ss::SS_output) and provenance

## ---- model-features ----
# The SS3 features the converter depends on. Run this on your own model and
# compare it with the yellowfin column in the guide.
ss3_model_features = function(ss3) {
    dat = ss3$dat; ctl = ss3$ctl; rep = ss3$rep
    get = function(x, nm, default = NA) {
        hit = intersect(nm, names(x))
        if (length(hit)) x[[hit[1]]] else default
    }
    n_rows = function(x) if (is.null(x) || !length(x)) 0L else NROW(x)
    patterns = function(x) {
        if (is.null(x) || !NROW(x)) return("none")
        p = as.data.frame(x)[[1]]
        paste(sort(unique(p)), collapse = ", ")
    }
    yn = function(x) if (isTRUE(x)) "yes" else "no"
    data.table::data.table(
        feature = c("SS3 version", "Model period (SS3 years)", "Seasons per year",
                    "Months per season", "Sexes", "Areas", "Ages (max)",
                    "Fishing fleets", "Survey / index fleets",
                    "Catch units (1 biomass, 2 numbers)",
                    "Index observations", "Length composition rows",
                    "Age composition rows", "Generalised size comp methods",
                    "Tag release groups", "Growth model (1 vB, 2 Richards, 3 age-K)",
                    "Growth patterns", "Platoons per growth pattern",
                    "Length selectivity patterns", "Age selectivity patterns",
                    "Time-block designs", "F method", "Stock-recruit function",
                    "Report.sso available"),
        value = c(
            as.character(get(rep, "SS_version", get(dat, "ReadVersion", NA))),
            paste(get(dat, "styr"), get(dat, "endyr"), sep = "-"),
            get(dat, "nseas"),
            paste(unique(get(dat, c("months_per_seas", "month_per_seas"))), collapse = ", "),
            get(dat, "Nsexes"),
            get(dat, "N_areas"),
            get(dat, "Nages"),
            get(dat, c("Nfleet", "Nfleets")),
            get(dat, "Nsurveys", 0),
            paste(unique(get(dat, "units_of_catch", NA)), collapse = ", "),
            n_rows(get(dat, "CPUE", NULL)),
            n_rows(get(dat, "lencomp", NULL)),
            n_rows(get(dat, "agecomp", NULL)),
            get(dat, "N_sizefreq_methods", 0),
            get(dat, "N_tag_groups", 0),
            get(ctl, "GrowthModel"),
            get(ctl, "N_GP"),
            get(ctl, "N_platoon", 1),
            patterns(get(ctl, "size_selex_types", NULL)),
            patterns(get(ctl, "age_selex_types", NULL)),
            get(ctl, "N_Block_Designs", 0),
            get(ctl, "F_Method"),
            get(ctl, "SR_function"),
            yn(!is.null(rep))
        )
    )
}
features = ss3_model_features(ss3)
features

## ---- convert ----
converted = ss3_to_opal(cfg, dir_root = proj_dir)
names(converted)   # data, parameters, lf, selectivity, meta, ss3

# Length compositions go into the data list through opal's prep_lf_data().
# SS3's mult_by_lencomp_N becomes opal's lf_var_adjust (its reciprocal).
lf_wide = dcast(converted$lf, fishery + year + month + ts ~ bin,
                value.var = "value", fill = 0) %>%
          setorder(fishery, ts) %>%
          as.data.frame()                 # opal objects cannot hold data.tables
converted$data = prep_lf_data(
    converted$data, lf_wide,
    lf_switch         = 1L,
    lf_keep_fisheries = as.integer(cfg$lf$keep_fisheries),
    lf_var_adjust     = converted$meta$lf_var_adjust,
    lf_cap            = cfg$lf$cap,
    lf_addtocomp      = cfg$lf$addtocomp
) %>% ss3_make_portable()

## ---- check-structure ----
# Test 1: is the opal input internally consistent?
checks = validate_opal_inputs(converted, strict = FALSE)
checks[, .(n_checks = .N, passed = sum(pass))]
checks[pass == FALSE]

## ---- check-reconcile ----
# Test 2: did the data and biology survive the conversion?
recon = ss3_reconcile(converted, tol = 1e-6)
recon[, .(quantity, ss3, opal, rel_diff, pass)]

## ---- build ----
# The opal assessment object. ss3_opal_obj() adds opal's default priors and
# the SS3 provenance, and opal_build() checks that the objective and gradient
# are finite at the starting values (the converted SS3 optimum).
assessment = ss3_opal_obj(converted, cfg,
                          control = list(eval.max = 10000, iter.max = 10000),
                          build   = TRUE)
summary(assessment)

## ---- check-bridge ----
# Test 3: with SS3's own selectivity-at-age and parameters, does opal
# reproduce SS3's spawning biomass? Nothing is estimated here.
bridge_obj = ss3_use_external_selectivity(
    assessment, converted$selectivity$sel_fa_external) %>%
    opal_build()
bridge = ss3_bridge_table(converted, opal_report(bridge_obj))
# ss3: SS3 SSB. opal: opal at SS3's values. replica_ss3: the same converted
# inputs replayed with SS3's catch and recruitment timing.
bridge[, .(series   = c("opal", "replay with SS3 timing"),
           max_abs  = c(max(abs(rel_opal)), max(abs(rel_replica_ss3))),
           mean_abs = c(mean(abs(rel_opal)), mean(abs(rel_replica_ss3))))]

## ---- start-objective ----
# Test 4: the starting objective and gradient. A non-finite value or a very
# large gradient points at a parameter or data element on the wrong scale.
start_check = function(x) {
    obj = opal_rtmb(x, fresh = TRUE)
    data.table(n_par = length(obj$par), objective = obj$fn(obj$par),
               max_abs_gradient = max(abs(obj$gr(obj$par))))
}
start_check(bridge_obj)

## ---- seed-selectivity ----
# SS3's selectivity is on age, opal's on length. Starting values come from
# fitting each opal curve to SS3's selectivity projected onto length bins.
seed = ss3_seed_par_sel(
    sel_target_fl = converted$selectivity$sel_target_fl,
    len_mid       = converted$data$len_mid,
    sel_type_f    = converted$data$sel_type_f,
    n_par         = cfg$selectivity$n_par)
data.table(fishery = seq_along(seed$rmse), rmse = round(seed$rmse, 4))

## ---- final-model ----
# Fit C2: opal's own length-based selectivity, estimated where SS3 estimated
# it, mirrored where SS3 mirrored it, with the prior-dominated parameters
# fixed. log_B0, the recruitment deviations and the four catchabilities are
# estimated as in opal's default map; growth stays at SS3's values.
final = assessment %>%
    ss3_use_estimated_selectivity(seed$par_sel, cfg,
                                  fix = cfg$selectivity$fix_prior_dominated) %>%
    opal_build()
start_check(final)

# What is estimated: map entries that are not NA.
estimated = rbindlist(lapply(names(final$parameters), function(p) {
    m = final$map[[p]]
    n = if (is.null(m)) length(final$parameters[[p]]) else length(unique(na.omit(m)))
    data.table(parameter = p, n_estimated = n)
}))[n_estimated > 0]
estimated

## ---- fit ----
# Two nlminb passes, then Newton steps to bring the gradient down, then
# opal_check() on the polished optimum.
t0  = proc.time()[["elapsed"]]
fit = opal_fit(final, n_passes = 2L, check = FALSE)
obj = opal_rtmb(fit, fresh = TRUE)
pol = ss3_newton_polish(obj, fit$fit$opt$par)
opt = fit$fit$opt
opt$par       = setNames(pol$par, names(opt$par))
opt$objective = pol$objective
fit = opal_attach_fit(fit, opt, check = TRUE,
                      check_args = list(gradient_tolerance = 1e-4))
fit_seconds = proc.time()[["elapsed"]] - t0

convergence = data.table(
    objective_start = start_check(final)$objective,
    objective_fit   = opt$objective,
    max_gradient    = pol$max_gradient,
    newton_steps    = pol$steps,
    check_passes    = isTRUE(fit$validation$fit$passes),
    positive_hessian = isTRUE(fit$validation$fit$metrics$positive_hessian),
    seconds         = fit_seconds)
convergence

## ---- save ----
fit_path = file.path(dir_out, "io-yft-workflow-fit.rds")
opal_save(fit, fit_path, compress = "xz", overwrite = TRUE)
# Later: fit = opal_read(fit_path, rebuild = TRUE)

## ---- results ----
# Quantities to compare with SS3, from opal_report() and Report.sso.
rep_opal = opal_report(fit)
rep_ss3  = converted$ss3$rep
d        = fit$data
n        = d$n_year

ts_ss3 = as.data.table(rep_ss3$timeseries)[Era == "TIME" & Yr >= cfg$time$first_ts &
                                               Yr <= cfg$time$last_ts]
ssb0_ss3 = ss3_derived_value(rep_ss3, c("SSB_Virgin", "SSB_Unfished"))
trajectory = rbind(
    data.table(model = "SS3", ts = seq_len(n), ssb = ts_ss3$SpawnBio,
               depletion = ts_ss3$SpawnBio / ssb0_ss3),
    data.table(model = "opal", ts = seq_len(n),
               ssb = rep_opal$spawning_biomass_y[seq_len(n)],
               depletion = rep_opal$static_depletion_y[seq_len(n)]))

# Index: opal's input is each quarterly series scaled to mean 1, so SS3's
# expected values are scaled by the same factor.
cpue = as.data.table(d$cpue_data)[, .(ts, fishery, index, obs = value, se)]
cpue[, opal := rep_opal$cpue_pred]
cpue_ss3 = as.data.table(rep_ss3$cpue)[Fleet == cfg$index$keep_ss3_fleets, .(ts = Yr, Obs, Exp)]
cpue = merge(cpue, cpue_ss3, by = "ts")
cpue[, ss3 := Exp * obs / Obs][, c("Obs", "Exp") := NULL]

# Length compositions aggregated over time, weighted by sample size.
lf_opal = ss3_lf_long(d, rep_opal$lf_pred)
lf_agg  = lf_opal[, .(obs = sum(obs * n), opal = sum(pred * n)), by = .(fishery, bin)]
lf_agg[, `:=`(obs = obs / sum(obs), opal = opal / sum(opal)), by = fishery]
# SS3's expected compositions: this run wrote no CompReport.sso (so r4ss has
# no lendbase), but data.ss_new carries them in its second (expected-value)
# section. The reader is SS3 3.24-specific.
dims    = converted$meta$dims
dat_exp = suppressWarnings(r4ss::SS_readdat_3.24(
    file.path(proj_dir, cfg$ss3_dir, cfg$files$dat), section = 2, verbose = FALSE))
lc = as.data.table(dat_exp$lencomp)[FltSvy %in% cfg$lf$keep_fisheries]
bin_cols = grep("^l[0-9]", names(lc), value = TRUE)
lb = melt(lc[, c("Yr", "FltSvy", bin_cols), with = FALSE], id.vars = c("Yr", "FltSvy"),
          variable.name = "bin_col", value.name = "Exp")
lb[, bin := dims$.len_bins[match(as.numeric(sub("^l", "", bin_col)), dims$.len_bins_ss3)]]
lb[, Exp := Exp / sum(Exp), by = .(Yr, FltSvy)]
lb = merge(lb, unique(lf_opal[, .(fishery, ts, n)]), by.x = c("FltSvy", "Yr"),
           by.y = c("fishery", "ts"))
lf_ss3 = lb[, .(ss3 = sum(Exp * n)), by = .(fishery = FltSvy, bin)]
lf_ss3[, ss3 := ss3 / sum(ss3), by = fishery]
lf_agg = merge(lf_agg, lf_ss3, by = c("fishery", "bin"), all.x = TRUE)

# Selectivity-at-age: SS3's against opal's estimated (length-based) curves.
sel_ss3 = converted$selectivity$sel_fa_external
selectivity = rbind(
    data.table(model = "SS3", fishery = rep(seq_len(nrow(sel_ss3)), times = ncol(sel_ss3)),
               age = rep(d$age_a, each = nrow(sel_ss3)), sel = as.vector(sel_ss3)),
    data.table(model = "opal", fishery = rep(seq_len(d$n_fishery), times = d$n_age),
               age = rep(d$age_a, each = d$n_fishery),
               sel = as.vector(rep_opal$sel_fya[, 1, ])))

summary_table = data.table(
    quantity = c("B0 (unfished SSB)", "Terminal SSB", "Terminal SSB / B0"),
    SS3  = c(ssb0_ss3, last(ts_ss3$SpawnBio), last(ts_ss3$SpawnBio) / ssb0_ss3),
    opal = c(exp(opt$par[["log_B0"]]), last(rep_opal$spawning_biomass_y[seq_len(n)]),
             last(rep_opal$static_depletion_y[seq_len(n)])))

## ---- cache ----
cache = list(
    created     = format(Sys.time(), tz = "UTC", usetz = TRUE),
    session     = list(r = as.character(getRversion()),
                       opal = as.character(packageVersion("opal")),
                       r4ss = as.character(packageVersion("r4ss"))),
    fleet_names = converted$meta$dims$.fleet_names,
    features    = features,
    checks      = checks,
    recon       = recon,
    bridge      = bridge,
    start       = rbind(cbind(model = "bridge (SS3 selectivity)", start_check(bridge_obj)),
                        cbind(model = "final (C2)", start_check(final))),
    seed_rmse   = data.table(fishery = seq_along(seed$rmse),
                             sel_type = converted$data$sel_type_f, rmse = seed$rmse),
    estimated   = estimated,
    convergence = convergence,
    trajectory  = trajectory,
    cpue        = cpue,
    lf          = lf_agg,
    selectivity = selectivity,
    summary     = summary_table
)
saveRDS(cache, file.path(dir_out, "io-yft-workflow-cache.rds"), compress = "xz")

## ---- reproduce-c2 ----
# This must be the same model as fit C2 in 03-fit-io-yft-age-k.r.
c2 = readRDS(file.path(dir_out, "io-yft-age-k-cache.rds"))$summary[id == "C2"]
cat(sprintf("Objective %.4f (03 cache C2: %.4f); %d parameters (C2: %d)\n",
            opt$objective, c2$objective_fit, length(opt$par), c2$n_par))
stopifnot("Fit does not reproduce C2" = abs(opt$objective - c2$objective_fit) < 1e-3,
          "Fit fails opal_check()" = convergence$check_passes)
