# Candidate versions of get_length_like() and get_weight_like(). Each change
# is a switch so candidates can be tested alone and combined:
#   atomic_product  catch-at-age x PLA (and, for weight comps, the rebinning)
#                   taped as one atomic matrix product: matmul is switched to
#                   "atomic" around that product only and restored after
#   block_rebin     rebin a fishery's rows with one product per scheme block
#                   instead of one product per row
#   skip_identity   no rebinning product when a scheme's matrix is the identity
#                   (decided from data before taping)
#   atomic_rebin    with block_rebin, tape the block rebinning as atomic too
#   vec_normalise   add the constant and normalise a whole block at once
# The density is still evaluated row by row on the same OBS-marked vector, so
# OSA residuals and simulation see the same structure as before. REPORT items
# (lf_pred, lf_pred_row / wf_pred, wf_pred_row) keep their shapes.

# Matrix product taped with matmul = mode ("atomic" or "plain"), restoring the
# previous setting afterwards
mm <- function(a, b, mode = NULL) {
  if (is.null(mode)) return(a %*% b)
  prev <- RTMB::TapeConfig(list(matmul = "NA"))$matmul
  RTMB::TapeConfig(list(matmul = mode))
  on.exit(RTMB::TapeConfig(list(matmul = prev)), add = TRUE)
  a %*% b
}

is_identity <- function(m) is.matrix(m) && nrow(m) == ncol(m) &&
  identical(unname(m), diag(nrow(m)))

# Predictions for one fishery block. catch_ya [n x n_age] (AD), pla
# [n_len x n_age]. Returns a list of row vectors (pred, normalised) and the
# per-row bin counts. `rebin_rows` is a list of rebinning matrices, one per row
# (NULL for legacy pooling, which then uses bmin/bmax).
comp_block <- function(catch_ya, pla, rebin_rows, scheme_rows, n_len, bmin, bmax,
                       addtocomp, o) {
  "[<-" <- ADoverload("[<-"); "c" <- ADoverload("c")
  n <- nrow(catch_ya)
  pred_yl <- mm(catch_ya, t(pla), if (o$atomic_product) "atomic" else NULL)
  rows <- vector("list", n)
  if (!is.null(rebin_rows) && o$block_rebin) {
    for (sc in unique(scheme_rows)) {
      ii <- which(scheme_rows == sc)
      R <- rebin_rows[[ii[1]]]
      blk <- pred_yl[ii, , drop = FALSE]
      if (!(o$skip_identity && is_identity(R)))
        blk <- mm(blk, t(R), if (o$atomic_rebin) "atomic" else NULL)
      if (o$vec_normalise) {
        blk <- blk + addtocomp
        # row sums as a product with a ones vector (rowSums() drops the AD class)
        rs <- mm(blk, matrix(1, ncol(blk), 1L), "plain")
        blk <- blk / c(rs)
      }
      for (k in seq_along(ii)) rows[[ii[k]]] <- blk[k, ]
    }
    if (!o$vec_normalise) rows <- lapply(rows, function(p) { p <- p + addtocomp; p / sum(p) })
    return(rows)
  }
  for (i in seq_len(n)) {
    pred <- pred_yl[i, ]
    if (!is.null(rebin_rows)) {
      R <- rebin_rows[[i]]
      if (!(o$skip_identity && is_identity(R))) pred <- c(R %*% pred)
    } else {
      if (bmin > 1) pred[bmin] <- sum(pred[1:bmin])
      if (bmax < n_len) pred[bmax] <- sum(pred[bmax:n_len])
      pred <- pred[bmin:bmax]
    }
    pred <- pred + addtocomp
    rows[[i]] <- pred / sum(pred)
  }
  rows
}

# Shared likelihood loop for length (kind = "lf") and weight (kind = "wf")
comp_like <- function(kind, obs_flat, obs_ints, obs_prop, catch_pred_fya, pla,
                      n_f_, fishery_f, year_fi, n_fi, minbin, maxbin,
                      removal_switch_f, switch_, n_bins_model, n_obs, log_tau,
                      addtocomp, n_int_fi, rebin_s, scheme, nbins_i,
                      wf_rebin_matrix = NULL, input_ints, o) {
  "[<-" <- ADoverload("[<-"); "c" <- ADoverload("c")
  # obs_* arrive already OBS-marked by the caller: RTMB's OBS() registers an
  # observation under the text of its argument (deparse(substitute(x))), so it
  # must be called on lf_obs_flat / wf_obs_flat etc. in the exported function,
  # not on a renamed copy here. input_ints is the unmarked integer vector.
  n_f <- length(n_f_)
  pred_f <- vector("list", n_f)
  pred_row <- vector("list", n_obs)
  use_schemes <- !is.null(rebin_s) && !is.null(scheme) && !is.null(nbins_i)
  lp <- numeric(n_obs)
  # Weight comps without schemes: the original folds the length-to-weight
  # rebinning into the PLA once (rebin_pla = R %*% pla) and pools bins per row
  pla_use <- pla
  n_model <- n_bins_model
  if (kind == "wf" && !use_schemes) pla_use <- wf_rebin_matrix %*% pla
  idx <- 0L
  obs_offset <- 0L
  for (j in seq_len(n_f)) {
    f <- fishery_f[j]
    ys <- year_fi[[j]]
    rows_j <- idx + seq_len(n_f_[j])
    catch_ya <- catch_pred_fya[f, ys, ]
    if (is.null(dim(catch_ya))) dim(catch_ya) <- c(1L, length(catch_ya))
    if (use_schemes) {
      single_scheme <- length(unique(scheme[rows_j])) == 1L
      preds <- comp_block(catch_ya, pla_use, rebin_s[scheme[rows_j]], scheme[rows_j],
                          n_model, NULL, NULL, addtocomp, o)
      pred_f[[j]] <- if (single_scheme) matrix(0, n_f_[j], nbins_i[rows_j[1L]]) else vector("list", n_f_[j])
    } else {
      bmin <- minbin[f]; bmax <- maxbin[f]
      preds <- comp_block(catch_ya, pla_use, NULL, NULL, n_model, bmin, bmax, addtocomp, o)
      pred_f[[j]] <- matrix(0, n_f_[j], bmax - bmin + 1L)
    }
    for (i in seq_len(n_f_[j])) {
      idx <- idx + 1L
      pred <- preds[[i]]
      nbins <- if (use_schemes) nbins_i[idx] else maxbin[f] - minbin[f] + 1L
      pred_row[[idx]] <- pred
      if (use_schemes && !single_scheme) pred_f[[j]][[i]] <- pred else pred_f[[j]][i, ] <- pred
      n_i <- n_fi[[j]][i]
      n_int_i <- if (is.null(n_int_fi)) as.integer(floor(n_i + 0.5)) else n_int_fi[[j]][i]
      include <- removal_switch_f[f] == 0 & if (switch_ == 3) n_int_i > 0 else n_i > 0
      if (include) {
        sl <- (obs_offset + 1):(obs_offset + nbins)
        if (switch_ == 1) lp[idx] <- -.opal_dmultinom(x = obs_flat[sl], prob = pred, log = TRUE)
        if (switch_ == 2)
          lp[idx] <- -RTMBdist::ddirichlet(x = obs_prop[sl], alpha = pred * n_i * exp(log_tau[f]), log = TRUE)
        if (switch_ == 3) {
          if (sum(input_ints[sl]) != n_int_i)
            stop("Dirichlet-multinomial counts must sum to the integer sample size.", call. = FALSE)
          lp[idx] <- -.opal_ddirmult(x = obs_ints[sl], size = n_int_i, alpha = pred * exp(log_tau[f]), log = TRUE)
        }
      }
      obs_offset <- obs_offset + nbins
    }
  }
  list(lp = lp, pred = pred_f, pred_row = pred_row)
}

make_lf <- function(o) {
  force(o)
  function(lf_obs_flat, lf_obs_ints, lf_obs_prop, catch_pred_fya, pla,
           lf_n_f, lf_fishery_f, lf_year_fi, lf_n_fi, lf_minbin, lf_maxbin,
           removal_switch_f, lf_switch, n_len, n_lf, log_lf_tau,
           lf_addtocomp = 1e-08, lf_n_int_fi = NULL, lf_rebin_s = NULL,
           lf_scheme = NULL, lf_nbins_i = NULL) {
    lf_input_ints <- lf_obs_ints
    if (lf_switch == 1) lf_obs_flat <- OBS(lf_obs_flat)
    if (lf_switch == 2) lf_obs_prop <- OBS(lf_obs_prop)
    if (lf_switch == 3) lf_obs_ints <- OBS(lf_obs_ints)
    r <- comp_like("lf", lf_obs_flat, lf_obs_ints, lf_obs_prop, catch_pred_fya, pla,
                   lf_n_f, lf_fishery_f, lf_year_fi, lf_n_fi, lf_minbin, lf_maxbin,
                   removal_switch_f, lf_switch, n_len, n_lf, log_lf_tau, lf_addtocomp,
                   lf_n_int_fi, lf_rebin_s, lf_scheme, lf_nbins_i,
                   input_ints = lf_input_ints, o = o)
    lf_pred <- r$pred; lf_pred_row <- r$pred_row
    REPORT(lf_pred); REPORT(lf_pred_row)
    r$lp
  }
}

make_wf <- function(o) {
  force(o)
  function(wf_obs_flat, wf_obs_ints, wf_obs_prop, catch_pred_fya, pla, wf_rebin_matrix,
           wf_n_f, wf_fishery_f, wf_year_fi, wf_n_fi, wf_minbin, wf_maxbin,
           removal_switch_f, wf_switch, n_wt, n_wf, log_wf_tau,
           wf_addtocomp = 1e-08, wf_n_int_fi = NULL, wf_rebin_s = NULL,
           wf_scheme = NULL, wf_nbins_i = NULL) {
    wf_input_ints <- wf_obs_ints
    if (wf_switch == 1) wf_obs_flat <- OBS(wf_obs_flat)
    if (wf_switch == 2) wf_obs_prop <- OBS(wf_obs_prop)
    if (wf_switch == 3) wf_obs_ints <- OBS(wf_obs_ints)
    r <- comp_like("wf", wf_obs_flat, wf_obs_ints, wf_obs_prop, catch_pred_fya, pla,
                   wf_n_f, wf_fishery_f, wf_year_fi, wf_n_fi, wf_minbin, wf_maxbin,
                   removal_switch_f, wf_switch, n_wt, n_wf, log_wf_tau, wf_addtocomp,
                   wf_n_int_fi, wf_rebin_s, wf_scheme, wf_nbins_i,
                   wf_rebin_matrix = wf_rebin_matrix, input_ints = wf_input_ints, o = o)
    wf_pred <- r$pred; wf_pred_row <- r$pred_row
    REPORT(wf_pred); REPORT(wf_pred_row)
    r$lp
  }
}

opts <- function(...) {
  o <- list(atomic_product = FALSE, block_rebin = FALSE, skip_identity = FALSE,
            atomic_rebin = FALSE, vec_normalise = FALSE)
  a <- list(...)
  o[names(a)] <- a
  o
}

# Candidates tested (name -> options). "refactor" reproduces the original
# structure through the shared code with every switch off.
candidates <- list(
  refactor = opts(),
  atomic_product = opts(atomic_product = TRUE),
  block_rebin = opts(block_rebin = TRUE),
  skip_identity = opts(skip_identity = TRUE),
  atomic_product_block = opts(atomic_product = TRUE, block_rebin = TRUE),
  atomic_product_skip = opts(atomic_product = TRUE, skip_identity = TRUE),
  atomic_product_block_skip = opts(atomic_product = TRUE, block_rebin = TRUE, skip_identity = TRUE),
  atomic_all_block_skip = opts(atomic_product = TRUE, block_rebin = TRUE, skip_identity = TRUE,
                               atomic_rebin = TRUE),
  atomic_product_block_skip_vec = opts(atomic_product = TRUE, block_rebin = TRUE,
                                       skip_identity = TRUE, vec_normalise = TRUE))

# Helpers and factories evaluate in an environment whose parent is the opal
# namespace, so the candidates see opal's internals (.opal_dmultinom, ...)
cand_env <- new.env(parent = asNamespace("opal"))
for (nm in c("mm", "is_identity", "comp_block", "comp_like", "make_lf", "make_wf")) {
  fn <- get(nm)
  environment(fn) <- cand_env
  assign(nm, fn, envir = cand_env)
  assign(nm, fn)
}
