# Test the dynamics candidate (plain and atomic sums) against the original,
# with the original and with the recommended length-composition candidate.
# Same checks as test-candidates.R. Writes results/candidates-dyn-<case>.csv.
source(file.path(Sys.getenv("REFACTOR_DIR"), "refactor-helpers.R"))
source(file.path(Sys.getenv("REFACTOR_DIR"), "cand-comp.R"))
source(file.path(Sys.getenv("REFACTOR_DIR"), "cand-dyn.R"))
args <- commandArgs(trailingOnly = TRUE)
cases <- if (length(args)) args else c("opaka", "io_yft")
lf_best <- Sys.getenv("LF_BEST", "atomic_product_block_skip")
flat <- function(r) { v <- unlist(r, use.names = FALSE); v[is.finite(v)] }
rel <- function(a, b) max(abs(a - b)) / max(abs(b), 1e-300)
gscaled <- function(a, b) max(abs(a - b)) / (1 + max(abs(b)))
variants <- list(
  dyn_plain = list(fns = list(do_dynamics = do_dynamics_vec), mode = "plain"),
  dyn_atomic = list(fns = list(do_dynamics = do_dynamics_vec), mode = "atomic"),
  lf_best = list(fns = list(get_length_like = make_lf(candidates[[lf_best]])), mode = "plain"),
  lf_best_dyn_plain = list(fns = list(get_length_like = make_lf(candidates[[lf_best]]),
                                      do_dynamics = do_dynamics_vec), mode = "plain"),
  lf_best_dyn_atomic = list(fns = list(get_length_like = make_lf(candidates[[lf_best]]),
                                       do_dynamics = do_dynamics_vec), mode = "atomic"))
for (case in cases) {
  x <- stored(case)
  P <- rbind(x$fit$opt$par, build(x)$par, starts_for(case))
  n <- if (case == "opaka") 200L else 20L
  swap_in(); obj0 <- build(x)
  f0 <- apply(P, 1, obj0$fn); g0 <- t(apply(P, 1, function(p) as.numeric(obj0$gr(p))))
  r0 <- obj0$report(P[1, ]); h0 <- obj0$he(P[1, ])
  t0 <- time_calls(obj0, P[1, ], n)
  rows <- list(data.frame(case = case, candidate = "original", fn_identical = TRUE, gr_identical = TRUE,
                          fn_max_rel = 0, gr_max_abs = 0, gr_max_scaled = 0, report_same_shape = TRUE,
                          report_max_rel = 0, he_max_rel = 0,
                          he_min_eig = min(eigen(h0, TRUE, only.values = TRUE)$values), t(t0)))
  cat(sprintf("%s original: ops %.0f fn %.3f gr %.3f\n", case, t0[["tape"]], t0[["fn_ms"]], t0[["gr_ms"]]))
  for (nm in names(variants)) {
    options(opal_dyn_matmul = variants[[nm]]$mode)
    swap_in(variants[[nm]]$fns); obj <- build(x)
    f1 <- apply(P, 1, obj$fn); g1 <- t(apply(P, 1, function(p) as.numeric(obj$gr(p))))
    r1 <- obj$report(P[1, ]); h1 <- obj$he(P[1, ])
    shape <- identical(names(r1), names(r0)) &&
      identical(rapply(r1, function(z) c(length(z), dim(z)), how = "list"),
                rapply(r0, function(z) c(length(z), dim(z)), how = "list"))
    tt <- time_calls(obj, P[1, ], n)
    row <- data.frame(case = case, candidate = nm, fn_identical = identical(f1, f0),
                      gr_identical = identical(g1, g0), fn_max_rel = max(abs(f1 - f0) / abs(f0)),
                      gr_max_abs = max(abs(g1 - g0)),
                      gr_max_scaled = max(sapply(seq_len(nrow(P)), function(i) gscaled(g1[i, ], g0[i, ]))),
                      report_same_shape = shape, report_max_rel = if (shape) rel(flat(r1), flat(r0)) else NA,
                      he_max_rel = rel(h1, h0), he_min_eig = min(eigen(h1, TRUE, only.values = TRUE)$values), t(tt))
    rows[[length(rows) + 1L]] <- row
    with(row, cat(sprintf("%s %-20s fn= %s gr= %s fn %.1e gr %.1e rep %.1e he %.1e | ops %.0f fn %.3f gr %.3f\n",
                          case, nm, fn_identical, gr_identical, fn_max_rel, gr_max_scaled,
                          report_max_rel, he_max_rel, tape, fn_ms, gr_ms)))
  }
  swap_in(); options(opal_dyn_matmul = NULL)
  write.csv(do.call(rbind, rows), file.path(Sys.getenv("REFACTOR_DIR"), "results",
                                            sprintf("candidates-dyn-%s.csv", case)), row.names = FALSE)
}
