# Paired timing of the original and candidate functions: fn, gr and Hessian
# per call at the stored estimate, in alternating rounds (original first in
# odd rounds), as in opal's benchmark harness. Length-composition models as
# stored, and the weight-composition versions of test-wf.R.
# Usage: Rscript timing.R [rounds]   Writes results/timing.csv.
source(file.path(Sys.getenv("REFACTOR_DIR"), "refactor-helpers.R"))
source(file.path(Sys.getenv("REFACTOR_DIR"), "cand-comp.R"))
rounds <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (is.na(rounds)) rounds <- 6L
best <- Sys.getenv("LF_BEST", "atomic_product")
as_wf <- local({ e <- new.env(); sys.source(file.path(Sys.getenv("REFACTOR_DIR"), "wf-data.R"), e); e$as_wf })
rows <- list()
for (case in c("opaka", "io_yft")) for (kind in c("lf", "wf")) {
  x <- stored(case)
  d <- if (kind == "wf") as_wf(x$data) else x$data
  p <- x$fit$opt$par
  swap_in(); obj0 <- build(x, d)
  swap_in(if (kind == "lf") list(get_length_like = make_lf(candidates[[best]]))
          else list(get_weight_like = make_wf(candidates[[best]])))
  obj1 <- build(x, d)
  swap_in()
  n <- if (case == "opaka") 100L else 10L
  he_rounds <- if (case == "opaka") rounds else 2L
  for (r in seq_len(rounds)) {
    ord <- if (r %% 2 == 1) c("original", "candidate") else c("candidate", "original")
    for (v in ord) {
      tt <- time_calls(if (v == "original") obj0 else obj1, p, n, blocks = 3L, he = r <= he_rounds)
      rows[[length(rows) + 1L]] <- data.frame(case = case, kind = kind, round = r, version = v,
                                              tape = tt[["tape"]], fn_ms = tt[["fn_ms"]], gr_ms = tt[["gr_ms"]],
                                              he_ms = if ("he_ms" %in% names(tt)) tt[["he_ms"]] else NA)
    }
  }
  s <- do.call(rbind, rows); s <- s[s$case == case & s$kind == kind, ]
  m <- aggregate(cbind(fn_ms, gr_ms, he_ms) ~ version, s, median, na.action = na.pass)
  cat(case, kind, "\n"); print(m, row.names = FALSE)
}
write.csv(do.call(rbind, rows), file.path(Sys.getenv("REFACTOR_DIR"), "results", "timing.csv"), row.names = FALSE)
