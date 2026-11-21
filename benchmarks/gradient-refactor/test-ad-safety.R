# AD safety: any fixed parameter could be estimated, so each candidate is also
# compared with the original with the map removed (every parameter, including
# growth, steepness, sigma_r, variance terms, initial F and the fixed
# selectivity elements, on the tape as an estimated parameter). With growth
# estimated the PLA and the weights are AD, not constants.
# Objective and gradient are compared at the stored estimate (fixed values at
# their stored values), at start 0 and at the jittered starts (jitter applied
# to the estimated elements only). Writes results/ad-safety.csv.
source(file.path(Sys.getenv("REFACTOR_DIR"), "refactor-helpers.R"))
source(file.path(Sys.getenv("REFACTOR_DIR"), "cand-comp.R"))
source(file.path(Sys.getenv("REFACTOR_DIR"), "cand-dyn.R"))
args <- commandArgs(trailingOnly = TRUE)
cases <- if (length(args)) args else c("opaka", "io_yft")
sets <- list(
  lf_atomic_product = list(fns = function() list(get_length_like = make_lf(candidates$atomic_product)), mode = "plain"),
  lf_best = list(fns = function() list(get_length_like = make_lf(candidates$atomic_product_block_skip)), mode = "plain"),
  lf_atomic_all = list(fns = function() list(get_length_like = make_lf(candidates$atomic_all_block_skip)), mode = "plain"),
  dyn_plain = list(fns = function() list(do_dynamics = do_dynamics_vec), mode = "plain"),
  dyn_atomic = list(fns = function() list(do_dynamics = do_dynamics_vec), mode = "atomic"))
gscaled <- function(a, b) max(abs(a - b)) / (1 + max(abs(b)))
rows <- list()
for (case in cases) {
  x <- stored(case)
  swap_in()
  obj_m <- build(x)                       # mapped, as stored
  xf <- x; xf$map <- list()               # everything estimated
  # full parameter vectors: stored estimate / start 0 / jittered starts placed
  # into the full parameter list (fixed elements at their stored values)
  full_at <- function(p) {
    pl <- obj_m$env$parList(p)
    unlist(pl[names(xf$parameters)])
  }
  P <- rbind(full_at(x$fit$opt$par), full_at(obj_m$par),
             t(apply(starts_for(case), 1, full_at)))
  obj0 <- build(xf)
  stopifnot(length(obj0$par) == ncol(P))
  f0 <- apply(P, 1, obj0$fn); g0 <- t(apply(P, 1, function(p) as.numeric(obj0$gr(p))))
  cat(sprintf("%s: %d parameters estimated (stored model: %d); objective at estimate %.8f (stored %.8f)\n",
              case, ncol(P), length(obj_m$par), f0[1], x$fit$opt$objective))
  for (nm in names(sets)) {
    options(opal_dyn_matmul = sets[[nm]]$mode)
    swap_in(sets[[nm]]$fns()); obj <- build(xf)
    f1 <- apply(P, 1, obj$fn); g1 <- t(apply(P, 1, function(p) as.numeric(obj$gr(p))))
    rows[[length(rows) + 1L]] <- data.frame(
      case = case, candidate = nm, n_par = ncol(P), fn_identical = identical(f1, f0),
      fn_max_rel = max(abs(f1 - f0) / abs(f0)), gr_identical = identical(g1, g0),
      gr_max_abs = max(abs(g1 - g0)),
      gr_max_scaled = max(sapply(seq_len(nrow(P)), function(i) gscaled(g1[i, ], g0[i, ]))),
      growth_gr_nonzero = any(abs(g0[, grepl("^log_(L1|L2|k|CV)", colnames(P))]) > 0))
    with(rows[[length(rows)]], cat(sprintf("  %-18s fn= %s (%.1e)  gr= %s (scaled %.1e, abs %.1e)\n",
                                           nm, fn_identical, fn_max_rel, gr_identical, gr_max_scaled, gr_max_abs)))
  }
  swap_in(); options(opal_dyn_matmul = NULL)
}
write.csv(do.call(rbind, rows), file.path(Sys.getenv("REFACTOR_DIR"), "results", "ad-safety.csv"),
          row.names = FALSE)
