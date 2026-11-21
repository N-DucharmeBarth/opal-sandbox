# Run opal's own tests of the affected code with the original functions and
# with a candidate set swapped into the opal namespace, and compare results.
# Usage: Rscript run-opal-tests.R <opal source dir> <set>   (set: original,
#        lf, lf_dyn_plain, lf_dyn_atomic). Writes results/opal-tests-<set>.csv.
# Needs testthat (the R 4.5.2 library here has it; the candidates are R code
# and do not depend on the R version).
args <- commandArgs(trailingOnly = TRUE)
opal_src <- args[[1]]
set <- args[[2]]
here <- Sys.getenv("REFACTOR_DIR")
suppressPackageStartupMessages({ library(opal); library(testthat) })
orig <- list(get_length_like = opal:::get_length_like, get_weight_like = opal:::get_weight_like,
             do_dynamics = opal:::do_dynamics)
source(file.path(here, "cand-comp.R"))
source(file.path(here, "cand-dyn.R"))
best <- Sys.getenv("LF_BEST", "atomic_product_block_skip")
fns <- switch(set,
  original = list(),
  lf = list(get_length_like = make_lf(candidates[[best]]), get_weight_like = make_wf(candidates[[best]])),
  lf_dyn_plain = list(get_length_like = make_lf(candidates[[best]]),
                      get_weight_like = make_wf(candidates[[best]]), do_dynamics = do_dynamics_vec),
  lf_dyn_atomic = list(get_length_like = make_lf(candidates[[best]]),
                       get_weight_like = make_wf(candidates[[best]]), do_dynamics = do_dynamics_vec),
  stop("unknown set"))
options(opal_dyn_matmul = if (set == "lf_dyn_atomic") "atomic" else "plain")
ns <- asNamespace("opal")
pkg <- as.environment("package:opal")
for (nm in names(orig)) {
  f <- if (!is.null(fns[[nm]])) fns[[nm]] else orig[[nm]]
  utils::assignInNamespace(nm, f, ns = "opal")
  if (exists(nm, envir = pkg, inherits = FALSE)) {
    unlockBinding(nm, pkg); assign(nm, f, envir = pkg); lockBinding(nm, pkg)
  }
}
files <- c("likelihoods-length", "likelihoods-schemes", "likelihoods-weight",
           "model-weight-integration", "composition-schemes-reporting", "osa-residuals",
           "dynamics", "dynamic-depletion", "dynamics-equilibrium", "harvest-rate-consistency",
           "model-gradients", "golden-opaka", "opaka-quickstart", "selectivity-integration",
           "assessment-regressions", "initial-harvest", "lambdas", "sex-ratio")
res <- testthat::test_dir(file.path(opal_src, "tests", "testthat"),
                          filter = paste0("^(", paste(files, collapse = "|"), ")$"),
                          package = "opal", load_package = "none", stop_on_failure = FALSE,
                          reporter = testthat::SilentReporter$new())
df <- as.data.frame(res)
out <- df[, c("file", "test", "nb", "failed", "skipped", "error", "warning")]
cat(sprintf("%s: %d tests, %d expectations, %d failed, %d errors, %d skipped\n", set, nrow(out),
            sum(out$nb), sum(out$failed), sum(out$error), sum(out$skipped)))
if (any(out$failed > 0 | out$error)) print(out[out$failed > 0 | out$error, c("file", "test")])
write.csv(out, file.path(here, "results", sprintf("opal-tests-%s.csv", set)), row.names = FALSE)
