# Fits with the original functions and with a candidate set, from the same
# starts, with opal's default nlminb settings (2 passes). Compares end points,
# evaluation counts and time. The original and candidate fits of each start
# run back to back, alternating which goes first.
# Usage: Rscript test-fits.R <case> <set> [start ids]   (set: lf, lf_dyn_plain)
# Writes results/fits-<case>-<set>.csv.
source(file.path(Sys.getenv("REFACTOR_DIR"), "refactor-helpers.R"))
source(file.path(Sys.getenv("REFACTOR_DIR"), "cand-comp.R"))
source(file.path(Sys.getenv("REFACTOR_DIR"), "cand-dyn.R"))
args <- commandArgs(trailingOnly = TRUE)
case <- args[[1]]; set <- args[[2]]
lf_best <- Sys.getenv("LF_BEST", "atomic_product")
fns <- switch(set,
  lf = list(get_length_like = make_lf(candidates[[lf_best]])),
  lf_dyn_plain = list(get_length_like = make_lf(candidates[[lf_best]]), do_dynamics = do_dynamics_vec))
x <- stored(case)
st <- starts_for(case)
ids <- if (length(args) > 2) as.integer(strsplit(args[[3]], ",")[[1]]) else 0:nrow(st)
swap_in(); obj0 <- build(x)
swap_in(fns); obj1 <- build(x)
swap_in()
lower <- unname(x$bounds$lower); upper <- unname(x$bounds$upper)
ctrl <- if (is.null(x$control)) list(eval.max = 10000L, iter.max = 10000L) else x$control
rows <- list()
for (k in seq_along(ids)) {
  id <- ids[k]
  start <- stats::setNames(as.numeric(if (id == 0L) obj0$par else st[id, ]), names(obj0$par))
  ord <- if (k %% 2 == 1) c("original", set) else c(set, "original")
  res <- list()
  for (v in ord) res[[v]] <- counted_fit(if (v == "original") obj0 else obj1, start, lower, upper, ctrl)
  for (v in c("original", set)) rows[[length(rows) + 1L]] <- cbind(case = case, start_id = id, version = v, res[[v]]$row)
  d <- max(abs(res[[set]]$par - res$original$par))
  cat(sprintf("%s start %2d: original %.2f s %4d gr obj %.6f | %s %.2f s %4d gr obj %.6f | max |par diff| %.2e\n",
              case, id, res$original$row$elapsed_s, res$original$row$gr_evals, res$original$row$objective,
              set, res[[set]]$row$elapsed_s, res[[set]]$row$gr_evals, res[[set]]$row$objective, d))
  rows[[length(rows)]]$max_par_diff <- d
  rows[[length(rows) - 1L]]$max_par_diff <- d
}
write.csv(do.call(rbind, rows), file.path(Sys.getenv("REFACTOR_DIR"), "results",
                                          sprintf("fits-%s-%s.csv", case, set)), row.names = FALSE)
