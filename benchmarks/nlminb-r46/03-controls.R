# Checks that do not depend on the jitter fits, run once per R version.
# Usage: Rscript 03-controls.R    (env NLMINB_BENCH_DIR, OPAL_SOURCE = opal checkout)
#  1. Replay of opal's stored benchmark trajectories
#     (tests/testthat/_benchmark/trajectories/{opaka,io_yft}.rds): fn/gr values
#     are saved so the two R versions can be compared bit for bit.
#  2. Cost of one fn and one gr call at the estimate of each benchmark model.
#  3. nlminb on a pure-R test function, with and without bounds, so the
#     optimiser can be timed without opal or RTMB.
source(file.path(Sys.getenv("NLMINB_BENCH_DIR"), "bench-helpers.R"))
suppressPackageStartupMessages(library(opal))
src <- Sys.getenv("OPAL_SOURCE")
out <- list(r_version = r_label())

# 1. Trajectory replay (inputs as in opal's benchmark_inputs())
traj_inputs <- function(model) {
  if (model == "opaka") return(opal::opaka_quickstart_inputs())
  env <- new.env()
  utils::data(list = c("io_yft_data", "io_yft_parameters", "io_yft_map"),
              package = "opal", envir = env)
  list(data = env$io_yft_data, parameters = env$io_yft_parameters, map = env$io_yft_map)
}
out$replay <- list()
for (model in c("opaka", "io_yft")) {
  tr <- readRDS(file.path(src, "tests/testthat/_benchmark/trajectories", paste0(model, ".rds")))
  inputs <- traj_inputs(model)
  obj <- suppressWarnings(RTMB::MakeADFun(cmb(opal_model, inputs$data), inputs$parameters,
                                          map = inputs$map, silent = TRUE))
  kinds <- vapply(tr$events, `[[`, "", "kind")
  idx <- sort(unlist(lapply(c("fn", "gr"), function(k) {
    i <- which(kinds == k)
    i[unique(round(seq(1, length(i), length.out = min(100L, length(i)))))]
  })))
  values <- lapply(tr$events[idx], function(e) obj[[e$kind]](e$par))
  same_as_stored <- mapply(function(v, e) identical(v, e$value), values, tr$events[idx])
  out$replay[[model]] <- list(index = idx, kind = kinds[idx], values = values,
                              same_as_stored = same_as_stored)
}

# 2. Cost per call at the estimate
per_call <- function(f, n, blocks = 5L) {
  stats::median(vapply(seq_len(blocks), function(b)
    system.time(for (i in seq_len(n)) f())[["elapsed"]] / n, numeric(1)))
}
out$per_call <- do.call(rbind, lapply(c("opaka", "io_yft"), function(case) {
  m <- bench_model(case)
  p <- m$mle
  invisible(m$obj$fn(p)); invisible(m$obj$gr(p))
  n_fn <- if (case == "opaka") 200L else 20L
  data.frame(case = case, r_version = r_label(),
             fn_s = per_call(function() m$obj$fn(p), n_fn),
             gr_s = per_call(function() m$obj$gr(p), n_fn))
}))

# 3. Pure-R test problem: extended Rosenbrock in 200 dimensions; the bounded
#    version has bounds [-1.5, 0.8], so about half the solution is on a bound.
rosen <- function(x) {
  i <- seq(1, length(x), by = 2)
  sum(100 * (x[i + 1] - x[i]^2)^2 + (1 - x[i])^2)
}
rosen_gr <- function(x) {
  i <- seq(1, length(x), by = 2)
  g <- numeric(length(x))
  g[i] <- -400 * x[i] * (x[i + 1] - x[i]^2) - 2 * (1 - x[i])
  g[i + 1] <- 200 * (x[i + 1] - x[i]^2)
  g
}
x0 <- rep(c(-1.2, 1), 100)
out$pure_r <- do.call(rbind, lapply(c("unbounded", "bounded"), function(type) {
  lo <- if (type == "bounded") -1.5 else -Inf
  up <- if (type == "bounded") 0.8 else Inf
  opt <- stats::nlminb(x0, rosen, rosen_gr, lower = lo, upper = up,
                       control = list(eval.max = 10000, iter.max = 10000))
  reps <- 20L
  t <- system.time(for (r in seq_len(reps))
    stats::nlminb(x0, rosen, rosen_gr, lower = lo, upper = up,
                  control = list(eval.max = 10000, iter.max = 10000)))[["elapsed"]] / reps
  data.frame(problem = type, r_version = r_label(), seconds = t,
             iterations = opt$iterations, fn_evals = opt$evaluations[["function"]],
             gr_evals = opt$evaluations[["gradient"]], objective = opt$objective,
             convergence = opt$convergence)
}))
out$pure_r_par <- NULL

saveRDS(out, file.path(bench_dir(), "results", "raw",
                       sprintf("controls-%s.rds", gsub("[^0-9]", "", getRversion()))))
print(out$per_call); print(out$pure_r)
cat("Replay identical to stored values:",
    vapply(out$replay, function(r) mean(r$same_as_stored), 0), "\n")
