# Fit every start of one case with one optimiser, in a fresh R process.
# Usage: Rscript run-method.R <case> <method>
# Env:   OPT_BENCH_DIR (benchmarks/optimizers), NLMINB_BENCH_DIR
#        (benchmarks/nlminb-r46: model helpers and the jittered starts),
#        OPT_REPS (default 3 for opaka, 1 for io_yft)
# Writes results/raw/<case>-<method>.rds. Errors are recorded, not fatal.
source(file.path(Sys.getenv("NLMINB_BENCH_DIR"), "bench-helpers.R"))
source(file.path(Sys.getenv("OPT_BENCH_DIR"), "opt-methods.R"))
args <- commandArgs(trailingOnly = TRUE)
case <- args[[1]]
method <- args[[2]]
reps <- as.integer(Sys.getenv("OPT_REPS", if (case == "opaka") "3" else "1"))
if (method == "nlminb2_atomic") RTMB::TapeConfig(atomic = "enable")
m <- bench_model(case)
starts <- readRDS(file.path(Sys.getenv("NLMINB_BENCH_DIR"), "results", "starts.rds"))[[case]]$starts
f <- counted(m$obj)
run_one <- function(start) {
  f$reset()
  err <- NA_character_
  res <- NULL
  t <- system.time(res <- tryCatch(methods[[method]](start, f, m$lower, m$upper),
                                   error = function(e) { err <<- conditionMessage(e); NULL }))[["elapsed"]]
  list(res = res, err = err, t = t, n = f$counts())
}
rows <- list(); pars <- list()
for (id in 0:nrow(starts)) {
  start <- stats::setNames(as.numeric(if (id == 0L) m$start0 else starts[id, ]), names(m$mle))
  runs <- lapply(seq_len(reps), function(r) run_one(start))
  r1 <- runs[[1]]
  ok <- !is.null(r1$res)
  g <- if (ok) as.numeric(m$obj$gr(r1$res$par)) else NA_real_
  rows[[length(rows) + 1L]] <- data.frame(
    case = case, method = method, start_id = id,
    elapsed_s = stats::median(vapply(runs, `[[`, 0, "t")),
    fn_evals = r1$n[["fn"]], gr_evals = r1$n[["gr"]], he_evals = r1$n[["he"]],
    objective = if (ok) r1$res$objective else NA_real_,
    objective_check = if (ok) m$obj$fn(r1$res$par) else NA_real_,
    converged = if (ok) r1$res$converged else NA,
    message = if (ok) as.character(r1$res$message) else paste("error:", r1$err),
    max_abs_grad = if (ok) max(abs(g)) else NA_real_,
    n_at_bound = if (ok) sum(r1$res$par <= m$lower | r1$res$par >= m$upper) else NA_integer_,
    mle_objective = m$mle_objective, stringsAsFactors = FALSE)
  pars[[length(pars) + 1L]] <- if (ok) r1$res$par else NULL
  with(rows[[length(rows)]], cat(sprintf("%s %s start %d: %.2f s, %d fn, %d gr, %d he, obj %.4f, %s\n",
    method, case, id, elapsed_s, fn_evals, gr_evals, he_evals, objective, message)))
}
out <- file.path(Sys.getenv("OPT_BENCH_DIR"), "results", "raw")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(rows = do.call(rbind, rows), pars = pars,
             r_version = r_label(), nloptr = tryCatch(as.character(packageVersion("nloptr")),
                                                      error = function(e) NA_character_)),
        file.path(out, sprintf("%s-%s.rds", case, method)))
