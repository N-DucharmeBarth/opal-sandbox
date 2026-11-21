# Cost of each step inside get_length_like() (multinomial, lf_switch = 1, the
# setting of both models). A profiling copy of the function stops after step
# k and returns the sum of everything computed so far, so its tape holds steps
# 1..k: 1 catch-at-age x PLA product per fishery; 2 rebinning (schemes) or bin
# pooling (legacy); 3 add constant and normalise; 4 multinomial density.
# The rest of the model is unchanged. Writes results/profile-lf.csv.
source(file.path(Sys.getenv("REFACTOR_DIR"), "refactor-helpers.R"))

lf_profile_fn <- function(stage) {
  force(stage)
  function(lf_obs_flat, lf_obs_ints, lf_obs_prop, catch_pred_fya, pla,
           lf_n_f, lf_fishery_f, lf_year_fi, lf_n_fi, lf_minbin, lf_maxbin,
           removal_switch_f, lf_switch, n_len, n_lf, log_lf_tau,
           lf_addtocomp = 1e-08, lf_n_int_fi = NULL, lf_rebin_s = NULL,
           lf_scheme = NULL, lf_nbins_i = NULL) {
    "[<-" <- ADoverload("[<-"); "c" <- ADoverload("c")
    stopifnot(lf_switch == 1)
    use_schemes <- !is.null(lf_rebin_s)
    lf_obs_flat <- OBS(lf_obs_flat)
    acc <- 0
    lp <- numeric(n_lf)
    idx <- 0L; obs_offset <- 0L
    for (j in seq_len(length(lf_n_f))) {
      f <- lf_fishery_f[j]
      pred_yl <- catch_pred_fya[f, lf_year_fi[[j]], ] %*% t(pla)
      if (stage == 1L) { acc <- acc + sum(pred_yl); next }
      if (!use_schemes) { bmin <- lf_minbin[f]; bmax <- lf_maxbin[f]; nbins <- bmax - bmin + 1L }
      for (i in seq_len(lf_n_f[j])) {
        idx <- idx + 1L
        pred <- pred_yl[i, ]
        if (use_schemes) {
          nbins <- lf_nbins_i[idx]
          pred <- c(lf_rebin_s[[lf_scheme[idx]]] %*% pred)
        } else {
          if (bmin > 1) pred[bmin] <- sum(pred[1:bmin])
          if (bmax < n_len) pred[bmax] <- sum(pred[bmax:n_len])
          pred <- pred[bmin:bmax]
        }
        if (stage == 2L) { acc <- acc + sum(pred); obs_offset <- obs_offset + nbins; next }
        pred <- pred + lf_addtocomp
        pred <- pred / sum(pred)
        if (stage == 3L) { acc <- acc + sum(log(pred)); obs_offset <- obs_offset + nbins; next }
        obs_i <- lf_obs_flat[(obs_offset + 1):(obs_offset + nbins)]
        lp[idx] <- -opal:::.opal_dmultinom(x = obs_i, prob = pred, log = TRUE)
        obs_offset <- obs_offset + nbins
      }
    }
    if (stage < 4L) lp[1] <- acc
    lp
  }
}

rows <- list()
for (case in c("opaka", "io_yft")) {
  x <- stored(case)
  p <- x$fit$opt$par
  n <- if (case == "opaka") 200L else 20L
  for (tape in c("default", "atomic")) {
    RTMB::TapeConfig(atomic = if (tape == "atomic") "enable" else "disable")
    for (st in c(0:4, 99)) {
      # stage 0: no length likelihood at all (n_lf set to 0 skips it)
      if (st == 0) { swap_in(); d <- x$data; d$n_lf <- 0L; obj <- build(x, d) }
      else if (st == 99) { swap_in(); obj <- build(x) }
      else { swap_in(list(get_length_like = lf_profile_fn(st))); obj <- build(x) }
      r <- time_calls(obj, p, n)
      lab <- c("none", "product", "rebin/pool", "normalise", "density (copy)")[st + 1]
      if (st == 99) lab <- "original function"
      rows[[length(rows) + 1L]] <- data.frame(case = case, tape = tape, stage = st,
                                              step = lab, t(r))
      cat(sprintf("%-6s %-7s %-18s ops %8.0f fn %7.3f gr %7.3f\n", case, tape, lab,
                  r[["tape"]], r[["fn_ms"]], r[["gr_ms"]]))
    }
  }
}
swap_in()
out <- file.path(Sys.getenv("REFACTOR_DIR"), "results")
dir.create(out, showWarnings = FALSE)
write.csv(do.call(rbind, rows), file.path(out, "profile-lf.csv"), row.names = FALSE)
