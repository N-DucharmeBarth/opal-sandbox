# Helpers for the R 4.5.2 vs R 4.6.2 nlminb benchmark. Sourced by every
# script in this folder. Needs only an installed opal (no renv).

bench_dir <- function() {
  d <- Sys.getenv("NLMINB_BENCH_DIR")
  if (!nzchar(d)) stop("Set NLMINB_BENCH_DIR to benchmarks/nlminb-r46")
  normalizePath(d, mustWork = TRUE)
}

r_label <- function() paste0("R ", getRversion())

# The two models, rebuilt from the fits stored in the opal package.
#   opaka   the quickstart fit (inst/extdata/opaka_quickstart_fit.rds)
#   io_yft  AGE-SF, the io-yft-ss3-comparison article's model
#           (inst/extdata/io_yft_age_sf_fit.rds)
# Returns the RTMB object, the bounds used for the stored fit, its estimate
# and the nlminb control list.
bench_model <- function(case) {
  suppressPackageStartupMessages(library(opal))
  file <- switch(case,
    opaka = system.file("extdata", "opaka_quickstart_fit.rds", package = "opal"),
    io_yft = system.file("extdata", "io_yft_age_sf_fit.rds", package = "opal"),
    stop("Unknown case: ", case))
  # read_opal_fit() refused the Opakapaka file in this build ("runtime identity
  # does not match its payload"), so the RTMB object is rebuilt from the stored
  # payload and checked against the stored objective.
  x <- readRDS(file)
  if (identical(x$format, "opal_obj")) x <- x$object
  obj <- RTMB::MakeADFun(cmb(opal_model, x$data), x$parameters, map = x$map,
                         silent = TRUE)
  mle <- x$fit$opt$par
  stopifnot(identical(names(obj$par), names(mle)), length(x$random) == 0L,
            abs(obj$fn(mle) - x$fit$opt$objective) < 1e-8)
  ctrl <- x$control
  if (is.null(ctrl)) ctrl <- list(eval.max = 10000L, iter.max = 10000L)
  list(obj = obj, lower = unname(x$bounds$lower), upper = unname(x$bounds$upper),
       mle = mle, mle_objective = x$fit$opt$objective, start0 = obj$par,
       control = ctrl)
}

# Jitter as described for SS3 3.30: a normal with Pr(lower) = 0.001 and
# Pr(upper) = 0.999; the estimate's cumulative probability K is moved to
# J ~ U(K - jf, K + jf) and mapped back. J is clamped to [0.0011, 0.9989] so
# the start stays inside the bounds. Parameters without two finite bounds are
# left at the estimate.
jitter_start <- function(mle, lower, upper, jf = 0.1) {
  ok <- is.finite(lower) & is.finite(upper)
  mu <- (lower + upper) / 2
  sd <- (upper - lower) / (2 * stats::qnorm(0.999))
  out <- mle
  k <- stats::pnorm(mle[ok], mu[ok], sd[ok])
  j <- stats::runif(sum(ok), k - jf, k + jf)
  j <- pmin(pmax(j, 0.0011), 0.9989)
  out[ok] <- stats::qnorm(j, mu[ok], sd[ok])
  out
}

# Two bounded nlminb passes, as in opal_fit() and benchmark_fit(), counting
# objective and gradient calls (wrapper adapted from opal's
# tests/testthat/helper-benchmark.R). Only the nlminb calls are timed.
counted_fit <- function(obj, start, lower, upper, control) {
  counts <- c(fn_evals = 0L, gr_evals = 0L)
  fn <- function(p) { counts[["fn_evals"]] <<- counts[["fn_evals"]] + 1L; obj$fn(p) }
  gr <- function(p) { counts[["gr_evals"]] <<- counts[["gr_evals"]] + 1L; obj$gr(p) }
  invisible(obj$fn(start))
  counts[] <- 0L
  passes <- vector("list", 2L)
  p <- start
  err <- NA_character_
  elapsed <- system.time(tryCatch({
    for (i in 1:2) {
      passes[[i]] <- stats::nlminb(p, fn, gr, lower = lower, upper = upper,
                                   control = control)
      p <- passes[[i]]$par
    }
  }, error = function(e) err <<- conditionMessage(e)))[["elapsed"]]
  if (!is.na(err)) {
    return(list(row = data.frame(
      elapsed_s = elapsed, fn_evals = counts[["fn_evals"]], gr_evals = counts[["gr_evals"]],
      iterations = NA_integer_, iter_pass1 = NA_integer_, iter_pass2 = NA_integer_,
      conv_pass1 = NA_integer_, conv_pass2 = NA_integer_, message = paste("error:", err),
      objective = NA_real_, max_abs_grad = NA_real_, n_at_bound = NA_integer_,
      stringsAsFactors = FALSE), par = rep(NA_real_, length(start))))
  }
  opt <- passes[[2]]
  g <- as.numeric(obj$gr(opt$par))
  list(
    row = data.frame(
      elapsed_s = elapsed,
      fn_evals = counts[["fn_evals"]], gr_evals = counts[["gr_evals"]],
      iterations = passes[[1]]$iterations + passes[[2]]$iterations,
      iter_pass1 = passes[[1]]$iterations, iter_pass2 = passes[[2]]$iterations,
      conv_pass1 = passes[[1]]$convergence, conv_pass2 = passes[[2]]$convergence,
      message = opt$message, objective = opt$objective,
      max_abs_grad = max(abs(g)),
      n_at_bound = sum(opt$par <= lower | opt$par >= upper),
      stringsAsFactors = FALSE),
    par = opt$par)
}
