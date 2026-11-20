# One fit in a fresh R process. Usage: Rscript 02-run-fit.R <case> <start_id>
# start_id 0 is the model's bundled initial values (opal_obj parameters),
# 1..N are the jittered starts in results/starts.rds.
source(file.path(Sys.getenv("NLMINB_BENCH_DIR"), "bench-helpers.R"))
args <- commandArgs(trailingOnly = TRUE)
case <- args[[1]]
id <- as.integer(args[[2]])
m <- bench_model(case)
start <- if (id == 0L) {
  m$obj$env$par[m$obj$env$lfixed()]
} else {
  readRDS(file.path(bench_dir(), "results", "starts.rds"))[[case]]$starts[id, ]
}
start <- stats::setNames(as.numeric(start), names(m$mle))
f0 <- m$obj$fn(start)
g0 <- as.numeric(m$obj$gr(start))
fit <- counted_fit(m$obj, start, m$lower, m$upper, m$control)
row <- cbind(data.frame(case = case, start_id = id, r_version = r_label(),
                        start_objective = f0, mle_objective = m$mle_objective),
             fit$row)
out <- file.path(bench_dir(), "results", "raw",
                 sprintf("%s-%02d-%s.rds", case, id, gsub("[^0-9]", "", getRversion())))
saveRDS(list(row = row, par = fit$par, start_fn = f0, start_gr = g0), out)
cat(sprintf("%s %s start %d: %.1f s, %d fn, %d gr, obj %.6f (MLE %.6f), conv %d\n",
            r_label(), case, id, row$elapsed_s, row$fn_evals, row$gr_evals,
            row$objective, row$mle_objective, row$conv_pass2))
