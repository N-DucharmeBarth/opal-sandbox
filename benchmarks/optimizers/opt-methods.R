# Optimisers compared in the optimizer benchmark. Each method takes the RTMB
# object, start, bounds and a list of counted fn/gr/he wrappers, and returns
# list(par, objective, converged, message). Calls are counted by the wrappers,
# so every method is charged for the same work in the same way.

ctrl_nlminb <- list(eval.max = 10000L, iter.max = 10000L)

nlminb_passes <- function(start, f, lower, upper, n_passes, control = ctrl_nlminb,
                          hessian = NULL) {
  p <- start
  for (i in seq_len(n_passes)) {
    opt <- stats::nlminb(p, f$fn, f$gr, hessian, lower = lower, upper = upper,
                         control = control)
    p <- opt$par
  }
  list(par = opt$par, objective = opt$objective, converged = opt$convergence == 0L,
       message = opt$message)
}

# Up to `max_steps` Newton steps from `opt`, each kept only if it stays inside
# the bounds and lowers the objective (as io_yft_fit_polished(guarded = TRUE)).
newton_finish <- function(opt, f, lower, upper, max_steps = 3L, tol = 1e-8) {
  p <- opt$par
  fx <- opt$objective
  g <- f$gr(p)
  for (i in seq_len(max_steps)) {
    if (max(abs(g)) < tol) break
    step <- tryCatch(solve(f$he(p), as.numeric(g)), error = function(e) NULL)
    if (is.null(step)) break
    p_new <- p - step
    if (any(p_new < lower | p_new > upper)) break
    f_new <- f$fn(p_new)
    if (!is.finite(f_new) || f_new > fx) break
    p <- p_new; fx <- f_new; g <- f$gr(p)
  }
  opt$par <- p
  opt$objective <- fx
  opt
}

nloptr_run <- function(start, f, lower, upper, algorithm) {
  o <- nloptr::nloptr(start, eval_f = f$fn, eval_grad_f = function(p) as.numeric(f$gr(p)),
                      lb = lower, ub = upper,
                      opts = list(algorithm = algorithm, ftol_rel = 1e-10, xtol_rel = 0,
                                  maxeval = 20000L))
  list(par = o$solution, objective = o$objective, converged = o$status > 0,
       message = o$message)
}

methods <- list(
  nlminb2 = function(start, f, lower, upper) nlminb_passes(start, f, lower, upper, 2L),
  nlminb1 = function(start, f, lower, upper) nlminb_passes(start, f, lower, upper, 1L),
  nlminb1_newton = function(start, f, lower, upper)
    newton_finish(nlminb_passes(start, f, lower, upper, 1L), f, lower, upper),
  nlminb_hessian = function(start, f, lower, upper)
    nlminb_passes(start, f, lower, upper, 1L, hessian = f$he),
  nlminb2_reltol8 = function(start, f, lower, upper)
    nlminb_passes(start, f, lower, upper, 2L, control = c(ctrl_nlminb, rel.tol = 1e-8)),
  lbfgsb = function(start, f, lower, upper) {
    o <- stats::optim(start, f$fn, function(p) as.numeric(f$gr(p)), method = "L-BFGS-B",
                      lower = lower, upper = upper, control = list(maxit = 10000L))
    list(par = o$par, objective = o$value, converged = o$convergence == 0L,
         message = o$message)
  },
  # L-BFGS-B stops at the first non-finite objective; the usual workaround
  # returns a large finite value (1e10) and a zero gradient element instead.
  lbfgsb_guard = function(start, f, lower, upper) {
    fn <- function(p) { v <- f$fn(p); if (is.finite(v)) v else 1e10 }
    gr <- function(p) { g <- as.numeric(f$gr(p)); g[!is.finite(g)] <- 0; g }
    o <- stats::optim(start, fn, gr, method = "L-BFGS-B",
                      lower = lower, upper = upper, control = list(maxit = 10000L))
    list(par = o$par, objective = f$fn(o$par), converged = o$convergence == 0L,
         message = o$message)
  },
  nloptr_lbfgs = function(start, f, lower, upper)
    nloptr_run(start, f, lower, upper, "NLOPT_LD_LBFGS"),
  nloptr_tnewton = function(start, f, lower, upper)
    nloptr_run(start, f, lower, upper, "NLOPT_LD_TNEWTON_PRECOND_RESTART"),
  # Same as nlminb2; the runner builds the model with TapeConfig(atomic = "enable")
  nlminb2_atomic = function(start, f, lower, upper) nlminb_passes(start, f, lower, upper, 2L)
)

method_labels <- c(
  nlminb2 = "nlminb, 2 passes (opal default)",
  nlminb1 = "nlminb, 1 pass",
  nlminb1_newton = "nlminb, 1 pass + up to 3 Newton steps",
  nlminb_hessian = "nlminb with analytic Hessian",
  nlminb2_reltol8 = "nlminb, 2 passes, rel.tol = 1e-8",
  lbfgsb = "optim L-BFGS-B",
  lbfgsb_guard = "optim L-BFGS-B, non-finite values guarded",
  nloptr_lbfgs = "nloptr LBFGS",
  nloptr_tnewton = "nloptr truncated Newton",
  nlminb2_atomic = "nlminb, 2 passes, atomic tape")

# Counted fn/gr/he wrappers around an RTMB object
counted <- function(obj) {
  n <- c(fn = 0L, gr = 0L, he = 0L)
  list(fn = function(p) { n[["fn"]] <<- n[["fn"]] + 1L; obj$fn(p) },
       gr = function(p) { n[["gr"]] <<- n[["gr"]] + 1L; obj$gr(p) },
       he = function(p) { n[["he"]] <<- n[["he"]] + 1L; obj$he(p) },
       counts = function() n,
       reset = function() n[] <<- 0L)
}
