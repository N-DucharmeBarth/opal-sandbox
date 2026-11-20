# Fit every start of one case with one PORT variant, in R 4.6.2.
# Usage: Rscript run-variants.R <case> <variant>
# Writes same-compiler/results/raw/<case>-<variant>.rds. Starts and fit
# settings are those of the main benchmark (results/starts.rds; start 0 is the
# stored model's initial values). Counts do not vary between runs; elapsed
# time is the median of NLMINB_REPS repeats (default 3 for opaka, 1 for io_yft).
dir <- Sys.getenv("NLMINB_BENCH_DIR")
source(file.path(dir, "bench-helpers.R"))
source(file.path(dir, "same-compiler", "port-variants.R"))
args <- commandArgs(trailingOnly = TRUE)
case <- args[[1]]
variant <- args[[2]]
reps <- as.integer(Sys.getenv("NLMINB_REPS", if (case == "opaka") "3" else "1"))
nl <- port_nlminb(variant)
m <- bench_model(case)
starts <- readRDS(file.path(dir, "results", "starts.rds"))[[case]]$starts
rows <- list()
pars <- list()
for (id in 0:nrow(starts)) {
  start <- if (id == 0L) m$start0 else starts[id, ]
  start <- stats::setNames(as.numeric(start), names(m$mle))
  fits <- lapply(seq_len(reps), function(r)
    counted_fit(m$obj, start, m$lower, m$upper, m$control, nlminb = nl))
  row <- fits[[1]]$row
  row$elapsed_s <- stats::median(vapply(fits, function(f) f$row$elapsed_s, 0))
  rows[[length(rows) + 1L]] <- cbind(data.frame(case = case, variant = variant,
                                                start_id = id,
                                                mle_objective = m$mle_objective), row)
  pars[[length(pars) + 1L]] <- fits[[1]]$par
  cat(sprintf("%s %s start %d: %.1f s, %d gr, obj %.6f, conv %s\n", variant, case, id,
              row$elapsed_s, row$gr_evals, row$objective, row$conv_pass2))
}
out <- file.path(dir, "same-compiler", "results", "raw")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(rows = do.call(rbind, rows), pars = pars),
        file.path(out, sprintf("%s-%s.rds", case, variant)))
