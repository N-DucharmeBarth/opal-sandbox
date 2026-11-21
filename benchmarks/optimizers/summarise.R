# Collect results/raw/*.rds into
#   results/fits.csv     one row per case, method and start, with the final
#                        objective relative to the best any method reached
#                        from that start (d_best) and to the stored estimate
#   results/summary.csv  per case and method, jittered starts only
source(file.path(Sys.getenv("OPT_BENCH_DIR"), "opt-methods.R"))
res <- file.path(Sys.getenv("OPT_BENCH_DIR"), "results")
raw <- lapply(dir(file.path(res, "raw"), full.names = TRUE), readRDS)
fits <- do.call(rbind, lapply(raw, `[[`, "rows"))
key <- paste(fits$case, fits$start_id)
best <- tapply(fits$objective, key, min, na.rm = TRUE)
fits$d_best <- fits$objective - best[key]
fits$d_mle <- fits$objective - fits$mle_objective
fits$reached_best <- !is.na(fits$d_best) & fits$d_best < 1e-4
fits$label <- method_labels[fits$method]
fits <- fits[order(fits$case, fits$method, fits$start_id), ]
write.csv(fits, file.path(res, "fits.csv"), row.names = FALSE)

jit <- fits[fits$start_id > 0, ]
summ <- do.call(rbind, lapply(split(jit, list(jit$case, jit$method), drop = TRUE), function(d) {
  ok <- !grepl("^error", d$message)
  data.frame(case = d$case[1], method = d$method[1], label = d$label[1], starts = nrow(d),
             errors = sum(!ok),
             reported_converged = sum(d$converged %in% TRUE),
             reached_best = sum(d$reached_best),
             reached_mle = sum(!is.na(d$d_mle) & d$d_mle < 1e-3),
             median_time_s = stats::median(d$elapsed_s[ok]),
             total_time_s = sum(d$elapsed_s),
             median_gr = stats::median(d$gr_evals[ok]),
             median_fn = stats::median(d$fn_evals[ok]),
             total_he = sum(d$he_evals),
             median_max_grad = stats::median(d$max_abs_grad[ok]))
}))
summ <- summ[order(summ$case, summ$total_time_s), ]
write.csv(summ, file.path(res, "summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE, digits = 3)
