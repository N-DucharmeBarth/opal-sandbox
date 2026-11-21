# Weight compositions: neither model has any, so each model's length data are
# supplied as weight-composition data (same observations, bins and schemes;
# for Opakapaka's legacy path wf_rebin_matrix is the identity), with the
# length likelihood switched off. The weight-composition version of each
# model should reproduce the length-composition objective. The candidates of
# cand-comp.R are then tested in get_weight_like() as in test-candidates.R.
# Writes results/candidates-wf-<case>.csv.
source(file.path(Sys.getenv("REFACTOR_DIR"), "refactor-helpers.R"))
source(file.path(Sys.getenv("REFACTOR_DIR"), "cand-comp.R"))
args <- commandArgs(trailingOnly = TRUE)
cases <- if (length(args)) args else c("opaka", "io_yft")
only <- Sys.getenv("CANDIDATES", "refactor,atomic_product,atomic_product_block_skip,atomic_all_block_skip")
cand <- candidates[strsplit(only, ",")[[1]]]
flat <- function(r) { v <- unlist(r, use.names = FALSE); v[is.finite(v)] }
rel <- function(a, b) max(abs(a - b)) / max(abs(b), 1e-300)
gscaled <- function(a, b) max(abs(a - b)) / (1 + max(abs(b)))

source(file.path(Sys.getenv("REFACTOR_DIR"), "wf-data.R"))

for (case in cases) {
  x <- stored(case)
  dw <- as_wf(x$data)
  P <- rbind(x$fit$opt$par, build(x)$par, starts_for(case))
  n <- if (case == "opaka") 200L else 20L
  swap_in()
  obj_lf <- build(x)
  obj0 <- build(x, dw)
  cat(sprintf("%s: length-comp objective %.10f, as weight comps %.10f (difference %.2e)\n", case,
              obj_lf$fn(P[1, ]), obj0$fn(P[1, ]), obj0$fn(P[1, ]) - obj_lf$fn(P[1, ])))
  f0 <- apply(P, 1, obj0$fn); g0 <- t(apply(P, 1, function(p) as.numeric(obj0$gr(p))))
  r0 <- obj0$report(P[1, ]); h0 <- obj0$he(P[1, ])
  t0 <- time_calls(obj0, P[1, ], n)
  rows <- list(data.frame(case = case, candidate = "original", wf_equals_lf = obj0$fn(P[1, ]) == obj_lf$fn(P[1, ]),
                          fn_identical = TRUE, gr_identical = TRUE, fn_max_rel = 0, gr_max_abs = 0,
                          gr_max_scaled = 0, report_same_shape = TRUE, report_max_rel = 0, he_max_rel = 0,
                          he_min_eig = min(eigen(h0, TRUE, only.values = TRUE)$values), t(t0)))
  cat(sprintf("%s original (wf): ops %.0f fn %.3f gr %.3f\n", case, t0[["tape"]], t0[["fn_ms"]], t0[["gr_ms"]]))
  for (nm in names(cand)) {
    swap_in(list(get_weight_like = make_wf(cand[[nm]])))
    obj <- build(x, dw)
    f1 <- apply(P, 1, obj$fn); g1 <- t(apply(P, 1, function(p) as.numeric(obj$gr(p))))
    r1 <- obj$report(P[1, ]); h1 <- obj$he(P[1, ])
    shape <- identical(names(r1), names(r0)) &&
      identical(rapply(r1, function(z) c(length(z), dim(z)), how = "list"),
                rapply(r0, function(z) c(length(z), dim(z)), how = "list"))
    tt <- time_calls(obj, P[1, ], n)
    row <- data.frame(case = case, candidate = nm, wf_equals_lf = NA, fn_identical = identical(f1, f0),
                      gr_identical = identical(g1, g0), fn_max_rel = max(abs(f1 - f0) / abs(f0)),
                      gr_max_abs = max(abs(g1 - g0)),
                      gr_max_scaled = max(sapply(seq_len(nrow(P)), function(i) gscaled(g1[i, ], g0[i, ]))),
                      report_same_shape = shape, report_max_rel = if (shape) rel(flat(r1), flat(r0)) else NA,
                      he_max_rel = rel(h1, h0), he_min_eig = min(eigen(h1, TRUE, only.values = TRUE)$values), t(tt))
    rows[[length(rows) + 1L]] <- row
    with(row, cat(sprintf("%s wf %-28s fn= %s gr= %s fn %.1e gr %.1e rep %.1e he %.1e | ops %.0f fn %.3f gr %.3f\n",
                          case, nm, fn_identical, gr_identical, fn_max_rel, gr_max_scaled,
                          report_max_rel, he_max_rel, tape, fn_ms, gr_ms)))
  }
  swap_in()
  write.csv(do.call(rbind, rows), file.path(Sys.getenv("REFACTOR_DIR"), "results",
                                            sprintf("candidates-wf-%s.csv", case)), row.names = FALSE)
}
