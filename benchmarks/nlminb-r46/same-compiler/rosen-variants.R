# The plain-R test problem of 03-controls.R (200-parameter extended
# Rosenbrock, unbounded and bounded in [-1.5, 0.8]) with stats::nlminb and
# each PORT variant, in R 4.6.2. Writes same-compiler/results/rosen.csv.
dir <- Sys.getenv("NLMINB_BENCH_DIR")
source(file.path(dir, "same-compiler", "port-variants.R"))
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
fit <- function(nl, lo, up) nl(x0, rosen, rosen_gr, lower = lo, upper = up,
                               control = list(eval.max = 10000, iter.max = 10000))
opt <- c(list(stats = stats::nlminb),
         sapply(c("A", "B", "C", "D", "E"), port_nlminb, simplify = FALSE))
rows <- list()
for (type in c("unbounded", "bounded")) {
  lo <- if (type == "bounded") -1.5 else -Inf
  up <- if (type == "bounded") 0.8 else Inf
  ref <- fit(stats::nlminb, lo, up)
  for (v in names(opt)) {
    o <- fit(opt[[v]], lo, up)
    rows[[length(rows) + 1L]] <- data.frame(
      problem = type, variant = v, iterations = o$iterations,
      fn_evals = o$evaluations[["function"]], gr_evals = o$evaluations[["gradient"]],
      objective = o$objective, convergence = o$convergence,
      identical_to_stats = identical(o[c("par", "objective", "iterations", "evaluations")],
                                     ref[c("par", "objective", "iterations", "evaluations")]))
  }
}
rows <- do.call(rbind, rows)
dir.create(file.path(dir, "same-compiler", "results"), showWarnings = FALSE)
write.csv(rows, file.path(dir, "same-compiler", "results", "rosen.csv"), row.names = FALSE)
print(rows)
