# Where does a gradient call spend its time? opal_model() is cut after each
# stage; the cut returns the sum of every AD quantity computed so far, so the
# tape holds all stages up to the cut. The cost of a stage is the difference
# between consecutive cuts. Run with the default tape and with
# TapeConfig(atomic = "enable"), for both cases, at the stored estimate.
# Writes results/profile.csv. Usage: Rscript profile-gradient.R
source(file.path(Sys.getenv("NLMINB_BENCH_DIR"), "bench-helpers.R"))
suppressPackageStartupMessages(library(opal))

b <- body(opal::opal_model)
txt <- vapply(as.list(b), function(s) paste(deparse(s, width.cutoff = 500L), collapse = " "), "")
at <- function(pattern) {
  i <- grep(pattern, txt)
  if (length(i) < 1L) stop("Cut point not found: ", pattern)
  i[1L]
}
# Quantities summed at each cut (only those that exist are used)
cuts <- list(
  biology = list(after = at("^if \\(exists\\(\"spawning_potential\""),
                 vars = c("pla", "weight_fya_mod", "M_a", "spawning_potential_a")),
  selectivity = list(after = at("^if \\(exists\\(\"sel_fa_external\""), vars = "sel_fya"),
  initial_numbers = list(after = at("^init <- get_initial_numbers"), vars = "init"),
  dynamics = list(after = at("^dyn <- do_dynamics"), vars = "dyn"),
  priors = list(after = at("^if \\(exists\\(\"priors\""),
                vars = c("lp_rec", "lp_init_rec", "lp_prior")),
  cpue_likelihood = list(after = at("^if \\(cpue_switch > 0\\)"), vars = "lp_cpue"),
  length_comp_likelihood = list(after = at("^if \\(lf_switch > 0 && n_lf > 0\\)"), vars = "lp_lf"),
  weight_comp_likelihood = list(after = at("^if \\(wf_switch > 0 && n_wf > 0\\)"), vars = "lp_wf"))
for (k in seq_along(cuts)[-1]) cuts[[k]]$vars <- c(cuts[[k - 1L]]$vars, cuts[[k]]$vars)

truncated_model <- function(cut) {
  fn <- opal::opal_model
  # unlist() drops RTMB's advector class, so each object is summed separately
  ret <- substitute({
    .add <- function(z) {
      if (is.list(z)) return(Reduce(`+`, lapply(z, .add), 0))
      if (is.numeric(z) || inherits(z, "advector")) sum(z) else 0
    }
    return(Reduce(`+`, lapply(mget(VARS, ifnotfound = list(NULL)), .add), 0))
  }, list(VARS = cut$vars))
  body(fn) <- as.call(c(as.list(b[seq_len(cut$after)]), as.list(ret)[-1]))
  fn
}

measure <- function(obj, p, n) {
  info <- .Call(get("InfoADFunObject", asNamespace("RTMB")), obj$env$ADFun$ptr)
  d <- -1
  q <- function() { d <<- -d; p + d * 1e-10 }
  invisible(obj$fn(q())); invisible(obj$gr(q()))
  blk <- function(f) stats::median(replicate(5, system.time(for (i in seq_len(n)) f(q()))[["elapsed"]] / n))
  c(tape = info$opstack_size, fn_ms = 1000 * blk(obj$fn), gr_ms = 1000 * blk(obj$gr))
}

rows <- list()
for (case in c("opaka", "io_yft")) {
  file <- switch(case,
    opaka = system.file("extdata", "opaka_quickstart_fit.rds", package = "opal"),
    io_yft = system.file("extdata", "io_yft_age_sf_fit.rds", package = "opal"))
  x <- readRDS(file)
  if (identical(x$format, "opal_obj")) x <- x$object
  p <- x$fit$opt$par
  n <- if (case == "opaka") 200L else 20L
  for (tape in c("default", "atomic")) {
    RTMB::TapeConfig(atomic = if (tape == "atomic") "enable" else "disable")
    for (nm in c(names(cuts), "full")) {
      model <- if (nm == "full") opal::opal_model else truncated_model(cuts[[nm]])
      obj <- RTMB::MakeADFun(cmb(model, x$data), x$parameters, map = x$map, silent = TRUE)
      r <- measure(obj, p, n)
      rows[[length(rows) + 1L]] <- data.frame(case = case, tape = tape, cut = nm,
                                              t(r), check.names = FALSE)
      cat(sprintf("%-6s %-7s %-24s tape %8.0f  fn %7.3f ms  gr %7.3f ms\n",
                  case, tape, nm, r[["tape"]], r[["fn_ms"]], r[["gr_ms"]]))
    }
  }
}
prof <- do.call(rbind, rows)
prof$order <- match(prof$cut, c(names(cuts), "full"))
out <- file.path(Sys.getenv("OPT_BENCH_DIR"), "results")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write.csv(prof, file.path(out, "profile.csv"), row.names = FALSE)
