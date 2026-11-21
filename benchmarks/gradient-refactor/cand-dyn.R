# Candidate do_dynamics(): the per-fishery loops inside each year and season
# are replaced by matrix operations over the fisheries with catch in that
# step. The sums over ages (sum(N * sel [* w]) and the predicted catch) are
# products with a ones/N column; .dyn_mm tapes them plainly (scalar ops, same
# order of summation as sum()) or as one atomic product (DYN_MATMUL=atomic).
# Generated from opal R/dynamics.R by a text substitution of the two loops.
# n copies of vector v as the rows of an [n x length(v)] matrix. Indexing keeps
# RTMB's AD class (matrix() on an advector does not).
.rows_of <- function(v, n) {
  m <- v[rep(seq_along(v), each = n)]
  dim(m) <- c(n, length(v))
  m
}

.dyn_mm <- function(a, b) {
  mode <- getOption("opal_dyn_matmul", "plain")
  prev <- RTMB::TapeConfig(list(matmul = "NA"))$matmul
  RTMB::TapeConfig(list(matmul = mode))
  on.exit(RTMB::TapeConfig(list(matmul = prev)), add = TRUE)
  a %*% b
}

do_dynamics_vec <- function(data, parameters,
                        B0, R0, alpha, beta, h = 0.95, sigma_r = 0.6,
                        M_a, spawning_potential_a, weight_fya,
                        init_number_a, init_number0_a, sel_fya, bias_adj_y = NULL) {
  
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  getAll(data, parameters, warn = FALSE)
  if (!all(catch_units_f %in% c(1, 2))) {
    stop("`catch_units_f` must contain only 1 (weight) or 2 (numbers).", call. = FALSE)
  }
  if (is.null(bias_adj_y)) bias_adj_y <- rep(1.0, n_year)
  fy <- first_yr_catch - first_yr + 1
  n_age1 <- n_age - 1
  S_a <- exp(-M_a / n_season)
  number_ysa <- array(0, dim = c(n_year + 1, n_season, n_age))
  number_ysa[1, 1,] <- init_number_a
  spawning_biomass_y <- numeric(n_year + 1)
  spawning_biomass_y[1] <- sum(number_ysa[1, 1,] * spawning_potential_a)
  number0_ysa <- array(0, dim = c(n_year + 1, n_season, n_age))
  number0_ysa[1, 1,] <- init_number0_a
  spawning_biomass0_y <- numeric(n_year + 1)
  spawning_biomass0_y[1] <- sum(number0_ysa[1, 1,] * spawning_potential_a)
  hrate_ysa  <- array(0, dim = c(n_year + 1, n_season, n_age))
  hrate_ysfa  <- array(0, dim = c(n_year + 1, n_season, n_fishery, n_age))
  catch_pred_fya <- array(0, dim = c(n_fishery, n_year, n_age))
  catch_pred_ysf <- array(0, dim = c(n_year, n_season, n_fishery))
  lp_penalty <- 0
  eps_denom <- 1e-6
  F_f <- numeric(n_fishery)
  h_rate_fa <- array(0, dim = c(n_fishery, n_age))
  f_weight <- which(catch_units_f == 1)
  f_numbers <- which(catch_units_f != 1)
  
  for (y in seq_len(n_year)) {
    for (s in seq_len(n_season)) {
      if (y >= fy) {
        F_f[] <- 0
        h_rate_fa[] <- 0
        N_ys <- number_ysa[y, s,]
        # all fisheries with catch in this step at once; weight-unit fisheries
        # first, then number-unit ones, as in the original loops
        for (grp in list(f_weight, f_numbers)) {
          act <- grp[catch_obs_ysf[y, s, grp] > 0]
          if (!length(act)) next
          S <- sel_fya[act, y, , drop = FALSE]
          dim(S) <- c(length(act), n_age)
          # sel_N = N * sel (then * weight), summed over ages, as in the original
          SN <- .rows_of(N_ys, length(act)) * S
          if (identical(grp, f_weight)) {
            Wt <- weight_fya[act, y, , drop = FALSE]; dim(Wt) <- c(length(act), n_age)
            SN <- SN * Wt
          }
          Nsum <- c(.dyn_mm(SN, matrix(1, n_age, 1L))) + eps_denom
          F_f[act] <- catch_obs_ysf[y, s, act] / Nsum
          h_rate_fa[act, ] <- F_f[act] * S
        }

        sum_F <- sum(F_f)
        tmp <- posfun(x = 1 - sum_F, eps = 0.001)
        lp_penalty <- lp_penalty + tmp$penalty
        hrate_ysfa[y, s,,] <- h_rate_fa
        hrate_ysa[y, s,] <- colSums(h_rate_fa)

        for (grp in list(f_weight, f_numbers)) {
          act <- grp[catch_obs_ysf[y, s, grp] > 0]
          if (!length(act)) next
          H <- h_rate_fa[act, , drop = FALSE]
          dim(H) <- c(length(act), n_age)
          C <- H * .rows_of(N_ys, length(act))
          catch_pred_fya[act, y, ] <- catch_pred_fya[act, y, ] + C
          CW <- if (identical(grp, f_weight)) {
            Wt <- weight_fya[act, y, , drop = FALSE]; dim(Wt) <- c(length(act), n_age); C * Wt
          } else C
          catch_pred_ysf[y, s, act] <- c(.dyn_mm(CW, matrix(1, n_age, 1L)))
        }
      }
      if (s < n_season) {
        number_ysa[y, s + 1,] <- number_ysa[y, s,] * (1 - hrate_ysa[y, s,]) * S_a
        number0_ysa[y, s + 1,] <- number0_ysa[y, s,] * S_a
      }
    }

    number_ysa[y + 1, 1, 2:n_age] <- number_ysa[y, n_season, 1:n_age1] * (1 - hrate_ysa[y, n_season, 1:n_age1]) * S_a[1:n_age1]
    number_ysa[y + 1, 1, n_age] <- number_ysa[y + 1, 1, n_age] + (number_ysa[y, n_season, n_age] * (1 - hrate_ysa[y, n_season, n_age]) * S_a[n_age])
    spawning_biomass_y[y + 1] <- sum(number_ysa[y + 1, 1,] * spawning_potential_a)
    number0_ysa[y + 1, 1, 2:n_age] <- number0_ysa[y, n_season, 1:n_age1] * S_a[1:n_age1]
    number0_ysa[y + 1, 1, n_age] <- number0_ysa[y + 1, 1, n_age] + (number0_ysa[y, n_season, n_age] * S_a[n_age])
    spawning_biomass0_y[y + 1] <- sum(number0_ysa[y + 1, 1,] * spawning_potential_a)

    number_ysa[y + 1, 1, 1] <- get_recruitment(sbio = spawning_biomass_y[y + 1], rdev = rdev_y[y], B0 = B0, alpha = alpha, beta = beta, sigma_r = sigma_r, bias_adj = bias_adj_y[y])
    number0_ysa[y + 1, 1, 1] <- get_recruitment(sbio = spawning_biomass0_y[y + 1], rdev = rdev_y[y], B0 = B0, alpha = alpha, beta = beta, sigma_r = sigma_r, bias_adj = bias_adj_y[y])
  }
  static_depletion_y <- spawning_biomass_y / B0
  dynamic_depletion_y <- spawning_biomass_y / spawning_biomass0_y
  
  REPORT(catch_pred_ysf)
  REPORT(catch_pred_fya)
  REPORT(hrate_ysa)
  REPORT(hrate_ysfa)
  REPORT(number0_ysa)
  REPORT(spawning_biomass_y)
  REPORT(spawning_biomass0_y)
  REPORT(static_depletion_y)
  REPORT(dynamic_depletion_y)
  ADREPORT(spawning_biomass_y)
  ADREPORT(spawning_biomass0_y)
  ADREPORT(static_depletion_y)
  ADREPORT(dynamic_depletion_y)
  
  return(list(number_ysa = number_ysa, number0_ysa = number0_ysa, lp_penalty = lp_penalty,
              catch_pred_fya = catch_pred_fya,
              hrate_ysfa = hrate_ysfa, hrate_ysa = hrate_ysa,
              catch_pred_ysf = catch_pred_ysf,
              spawning_biomass_y = spawning_biomass_y,
              spawning_biomass0_y = spawning_biomass0_y,
              static_depletion_y = static_depletion_y,
              dynamic_depletion_y = dynamic_depletion_y))
}

dyn_env <- new.env(parent = asNamespace("opal"))
for (nm in c(".rows_of", ".dyn_mm", "do_dynamics_vec")) {
  fn <- get(nm); environment(fn) <- dyn_env; assign(nm, fn, envir = dyn_env); assign(nm, fn)
}
