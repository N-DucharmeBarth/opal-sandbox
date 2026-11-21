# Supply a model's length-composition data as weight-composition data (same
# observations, bins and schemes; wf_rebin_matrix is the identity for the
# legacy path) and switch the length likelihood off. Used by test-wf.R and
# timing.R.
as_wf <- function(d) {
  w <- d
  cp <- function(to, from) if (!is.null(d[[from]])) w[[to]] <<- d[[from]]
  for (s in c("obs_flat", "obs_ints", "obs_prop", "n_f", "fishery_f", "year", "season",
              "fishery", "n", "minbin", "maxbin", "addtocomp", "n_int", "n_int_fi",
              "year_fi", "n_fi", "rebin_s", "scheme", "nbins_i", "scheme_edges", "var_adj"))
    cp(paste0("wf_", s), paste0("lf_", s))
  w$wf_switch <- d$lf_switch
  w$n_wf <- d$n_lf
  w$n_wt <- d$n_len
  w$wt_bin_start <- d$len_bin_start
  w$wt_bin_width <- d$len_bin_width
  w$wf_rebin_matrix <- diag(d$n_len)
  if (is.null(w$log_wf_tau)) w$log_wf_tau <- rep(0, d$n_fishery)
  if (!is.null(d$lambdas$lf)) w$lambdas$wf <- d$lambdas$lf
  w$lf_switch <- 0L
  w
}
