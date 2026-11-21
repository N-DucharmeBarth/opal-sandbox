# Shared helpers for the gradient refactor study. Needs NLMINB_BENCH_DIR
# (benchmarks/nlminb-r46) for bench_model() and the jittered starts.
source(file.path(Sys.getenv("NLMINB_BENCH_DIR"), "bench-helpers.R"))
suppressPackageStartupMessages(library(opal))

# Original opal functions, kept so a variant can be swapped in and out
orig <- list(get_length_like = opal:::get_length_like,
             get_weight_like = opal:::get_weight_like,
             do_dynamics = opal:::do_dynamics,
             opal_model = opal:::opal_model)

# Put a set of functions into the opal namespace ("orig" restores all)
swap_in <- function(fns = list()) {
  for (nm in names(orig)) {
    f <- if (!is.null(fns[[nm]])) fns[[nm]] else orig[[nm]]
    # functions defined at top level get the namespace as their environment;
    # closures (e.g. from a factory) get a child of the namespace instead
    if (identical(environment(f), globalenv())) environment(f) <- asNamespace("opal")
    else if (identical(parent.env(environment(f)), globalenv()))
      parent.env(environment(f)) <- asNamespace("opal")
    utils::assignInNamespace(nm, f, ns = "opal")
  }
  invisible(NULL)
}

stored <- function(case) {
  file <- switch(case,
    opaka = system.file("extdata", "opaka_quickstart_fit.rds", package = "opal"),
    io_yft = system.file("extdata", "io_yft_age_sf_fit.rds", package = "opal"))
  x <- readRDS(file)
  if (identical(x$format, "opal_obj")) x <- x$object
  x
}

# Build the RTMB object from data/parameters/map with whatever opal_model and
# helpers are currently in the namespace. opal_model is looked up at call time
# so a swapped-in model function is used.
build <- function(x, data = x$data) {
  RTMB::MakeADFun(cmb(get("opal_model", asNamespace("opal")), data), x$parameters,
                  map = x$map, silent = TRUE)
}

tape_ops <- function(obj)
  .Call(get("InfoADFunObject", asNamespace("RTMB")), obj$env$ADFun$ptr)$opstack_size

# fn and gr time per call (ms): medians of `blocks` blocks of n calls, with the
# parameter vector alternated by +-1e-10 so no call reuses a cached result
time_calls <- function(obj, p, n, blocks = 5L, he = FALSE) {
  d <- -1
  q <- function() { d <<- -d; p + d * 1e-10 }
  invisible(obj$fn(q())); invisible(obj$gr(q()))
  blk <- function(f, n) stats::median(replicate(blocks,
    system.time(for (i in seq_len(n)) f(q()))[["elapsed"]] / n))
  out <- c(tape = tape_ops(obj), fn_ms = 1000 * blk(obj$fn, n), gr_ms = 1000 * blk(obj$gr, n))
  if (he) out <- c(out, he_ms = 1000 * blk(obj$he, 1L))
  out
}

starts_for <- function(case)
  readRDS(file.path(Sys.getenv("NLMINB_BENCH_DIR"), "results", "starts.rds"))[[case]]$starts
