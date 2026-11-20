# Draw the jittered starting values once, so both R versions start from the
# same points. Usage: Rscript 01-make-starts.R
# Env: NLMINB_BENCH_DIR, NLMINB_N_OPAKA (20), NLMINB_N_IO_YFT (10), NLMINB_JF (0.1)
source(file.path(Sys.getenv("NLMINB_BENCH_DIR"), "bench-helpers.R"))
n <- c(opaka = as.integer(Sys.getenv("NLMINB_N_OPAKA", "20")),
       io_yft = as.integer(Sys.getenv("NLMINB_N_IO_YFT", "10")))
jf <- as.numeric(Sys.getenv("NLMINB_JF", "0.1"))

starts <- list()
for (case in names(n)) {
  m <- bench_model(case)
  # Starts with a non-finite objective or gradient stop nlminb at once; they
  # are redrawn and counted.
  set.seed(20261120)
  s <- matrix(NA_real_, 0, length(m$mle))
  rejected <- 0L
  while (nrow(s) < n[[case]]) {
    p <- jitter_start(m$mle, m$lower, m$upper, jf)
    if (is.finite(m$obj$fn(p)) && all(is.finite(m$obj$gr(p)))) {
      s <- rbind(s, p)
    } else {
      rejected <- rejected + 1L
    }
  }
  colnames(s) <- names(m$mle)
  rownames(s) <- NULL
  starts[[case]] <- list(
    starts = s, jf = jf, rejected = rejected,
    mle = m$mle, mle_objective = m$mle_objective,
    n_par = length(m$mle),
    n_jittered = sum(is.finite(m$lower) & is.finite(m$upper)),
    made_with = r_label())
  cat(sprintf("%s: %d starts (%d non-finite draws rejected), %d of %d parameters jittered\n",
              case, n[[case]], rejected, starts[[case]]$n_jittered, starts[[case]]$n_par))
}
saveRDS(starts, file.path(bench_dir(), "results", "starts.rds"))
