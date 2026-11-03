# Helpers to edit MULTIFAN-CL (MFCL) input files for the south-west Pacific
# albacore bridge: collapse the two-region model to model region 1 (.frq,
# .ini, .par), set movement to zero, and run MFCL.
#
# The .par editing works section by section: a par file is a sequence of
# "# header" lines, each followed by its value lines. Sections are matched by
# header text and occurrence, so the order and the unedited sections are kept
# exactly as MFCL wrote them. The generic par readers are in
# shared/mfcl-to-opal/mfcl-read.R.

source(file.path(this.path::this.proj(), "shared", "mfcl-to-opal", "mfcl-read.R"))

# Index of the line after the first header that starts with `pattern`.
after_header = function(lines, pattern, from = 1L) {
  i = which(startsWith(lines, pattern))
  i = i[i >= from]
  if (!length(i)) stop("Header not found: ", pattern)
  i[1] + 1L
}

# .frq -------------------------------------------------------------------------

#' Subset an MFCL .frq file to one region and a set of fisheries.
#'
#' @param keep_fish fisheries to keep, in their new order (renumbered 1..n).
#' @param region the region to keep; all kept fisheries must be in it.
subset_frq_region = function(frq_in, frq_out, keep_fish, region = 1L) {
  x = mfcl_read_lines(frq_in)
  # Dimension line: the first non-comment line.
  i_dim = which(!startsWith(x, "#"))[1]
  dims = mfcl_num(x[i_dim])
  n_reg = dims[1]
  n_fish = dims[2]
  dims[1] = 1
  dims[2] = length(keep_fish)
  x[i_dim] = paste0("\t", paste(mfcl_fmt(dims), collapse = "\t"), "\t")

  i = after_header(x, "# Relative Region Size")
  x[i] = mfcl_fmt(mfcl_num(x[i])[region])

  i = after_header(x, "# Region in which each fishery is located")
  fish_reg = mfcl_num(x[i])
  stopifnot(length(fish_reg) == n_fish, all(fish_reg[keep_fish] == region))
  x[i] = paste(rep(1, length(keep_fish)), collapse = " ")

  # Incidence matrix: n_reg - 1 rows (plus a blank line); none for one region.
  i = after_header(x, "# Incidence matrix")
  j = i
  while (!startsWith(x[j], "#")) j = j + 1L
  drop_inc = seq(i, j - 1L)

  i = after_header(x, "# Data flags")
  for (k in 0:4) {
    v = mfcl_num(x[i + k])
    stopifnot(length(v) == n_fish)
    x[i + k] = paste(mfcl_fmt(v[keep_fish]), collapse = " ")
  }

  i = after_header(x, "# Season-region flags")
  j = i
  while (!startsWith(x[j], "#")) {
    v = mfcl_num(x[j])
    if (length(v) == n_reg) x[j] = mfcl_fmt(v[region])
    j = j + 1L
  }

  i_ds = after_header(x, "# Datasets")
  i_rec = after_header(x, "# age_nage") + 1L
  recs = x[i_rec:length(x)]
  recs = recs[nzchar(trimws(recs))]
  fish = vapply(strsplit(trimws(recs), "[[:space:]]+"), function(v) as.integer(v[4]), 1L)
  keep = fish %in% keep_fish
  new_id = match(fish[keep], keep_fish)
  recs = recs[keep]
  recs = vapply(seq_along(recs), function(k) {
    v = strsplit(trimws(recs[k]), "[[:space:]]+")[[1]]
    v[4] = as.character(new_id[k])
    paste(v, collapse = " ")
  }, "")
  ds = mfcl_num(x[i_ds])
  ds[1] = length(recs)
  x[i_ds] = paste(mfcl_fmt(ds), collapse = " ")

  out = c(x[setdiff(seq_len(i_rec - 1L), drop_inc)], recs)
  writeLines(out, frq_out)
  invisible(list(n_records = length(recs), old_fishery = fish[keep], new_fishery = new_id))
}

# .ini -------------------------------------------------------------------------

#' Reduce an MFCL .ini file to one region (no tags in this model).
subset_ini_region = function(ini_in, ini_out, region = 1L) {
  x = mfcl_read_lines(ini_in)
  i = after_header(x, "# region flags")
  while (!startsWith(x[i], "#")) {
    x[i] = mfcl_fmt(mfcl_num(x[i])[region])
    i = i + 1L
  }
  # One region: the diffusion matrix keeps one (unused) column per movement.
  i = after_header(x, "# diffusion coffs")
  while (!startsWith(x[i], "#")) {
    x[i] = "0"
    i = i + 1L
  }
  i = after_header(x, "# recruitment distribution by region")
  x[i] = "1"
  writeLines(x, ini_out)
}

# .par -------------------------------------------------------------------------

matrix_body = function(m) {
  apply(m, 1, function(r) paste0(" ", paste(mfcl_fmt(r), collapse = " ")))
}

#' Movement fixed at zero: every movement matrix set to the identity.
par_identity_movement = function(secs) {
  k = names(secs)[startsWith(names(secs), "# Movement period")]
  for (kk in k) {
    n = length(secs[[kk]]$body)
    m = diag(n)
    secs[[kk]]$body = matrix_body(m)
  }
  secs
}

flag_body = function(v) c(paste0(" ", paste(mfcl_fmt(v), collapse = " ")), "")

set_flags = function(secs, pattern, values, occ = 1L) {
  k = mfcl_key(secs, pattern, occ)
  v = mfcl_num(secs[[k]]$body[1])
  v[as.integer(names(values))] = unlist(values)
  secs[[k]]$body = flag_body(v)
  secs
}

# Keep the columns `keep` of every non-blank line, preserving blank lines.
subset_cols_body = function(body, keep) {
  vapply(body, function(b) {
    if (!nzchar(trimws(b))) return(b)
    paste0(" ", paste(mfcl_fmt(mfcl_num(b)[keep]), collapse = " "))
  }, "", USE.NAMES = FALSE)
}

#' Build a one-region par from a fitted two-region par.
#'
#' `template` is the par written by `mfclo64 <frq> <ini> 00.par -makepar` for
#' the one-region .frq/.ini: it supplies every section shaped by the number of
#' regions, movement periods or catch records (all structural or zero in this
#' model, which the function checks). Estimates and flags are taken from `src`:
#' fish flags, selectivity, catchability and other fishery-indexed values are
#' subset to `keep_fish`; region-indexed values to `region`. The total
#' population scaling parameter is shifted by log(rec_share) so that region
#' `region`'s recruitment is unchanged.
build_region_par = function(src, template, keep_fish, region = 1L, rec_share) {
  n_fish_src = nrow(mfcl_fish_flags(src))
  out = template
  copy = c("# The parest_flags", "# age flags", "# percent maturity",
            "# relative recruitment", "# natural mortality coefficient",
            "# age-class related parameters", "# maturity at length",
            "# The von Bertalanffy parameters", "# Extra par for Richards",
            "# Seasonal growth parameters", "# Cohort specific growth deviations",
            "# Variance parameters", "# log_length variance", "# log_weight variance",
            "# length_tot_exp", "# year_flags", "# season_flags",
            "# First year in model", "# Do multsp Historical_flags",
            "# implicit total populations scaling", "# rec init pop level difference",
            "# recruitment times", "# new orthogonal coefficients")
  # Named sections whose values sit in the following unnamed "#" section.
  value_next = c("# relative recruitment", "# natural mortality coefficient",
                  "# average catchability coefficients",
                  "# initial trend in catchability coefficients", "# q0_miss",
                  "# extra fishery parameters", "# age-class related parameters")
  for (k in names(out)) {
    base = sub(" \\[[0-9]+\\]$", "", k)
    hit = copy[startsWith(base, copy)]
    if (length(hit)) {
      occ = as.integer(sub(".*\\[([0-9]+)\\]$", "\\1", k))
      ks = mfcl_key(src, hit[1], occ)
      out[[k]]$body = src[[ks]]$body
      if (hit[1] %in% value_next) out[[mfcl_next_key(out, k)]]$body = src[[mfcl_next_key(src, ks)]]$body
    }
  }
  # Fishery-column sections.
  for (h in c("# average catchability coefficients",
              "# initial trend in catchability coefficients", "# q0_miss",
              "# extra fishery parameters")) {
    ks = mfcl_next_key(src, mfcl_key(src, h))
    out[[mfcl_next_key(out, mfcl_key(out, h))]]$body = subset_cols_body(src[[ks]]$body, keep_fish)
  }
  # Fish flags: rows, then renumber the selectivity groups (flag 24).
  for (occ in 1:2) {
    ff = mfcl_fish_flags(src, occ)[keep_fish, , drop = FALSE]
    ff[, 24] = match(ff[, 24], unique(ff[, 24]))
    out[[mfcl_key(out, "# fish flags", occ)]]$body <-
      apply(ff, 1, function(r) paste(sprintf("%5s", mfcl_fmt(r)), collapse = ""))
  }
  # Selectivity: one row per fishery-season (fish flag 74 seasons, 0 = 1).
  ff = mfcl_fish_flags(src)
  n_rows = pmax(1, ff[, 74])
  rows = split(seq_len(sum(n_rows)), rep(seq_len(n_fish_src), n_rows))
  ks = mfcl_key(src, "# fishery selectivity")
  out[[mfcl_key(out, "# fishery selectivity")]]$body = src[[ks]]$body[unlist(rows[keep_fish])]
  # Region 1 relative initial population (a blank line, then one row per region).
  ks = mfcl_key(src, "# relative initial population")
  b = src[[ks]]$body
  out[[mfcl_key(out, "# relative initial population")]]$body = c(b[1], b[1 + region])
  # Total population scale.
  ks = mfcl_key(src, "# total populations scaling parameter")
  out[[mfcl_key(out, "# total populations scaling parameter")]]$body <-
    mfcl_fmt(mfcl_num(src[[ks]]$body) + log(rec_share))
  # Check: template-supplied sections in src hold only zeros or structure.
  used_src = c(copy, value_next, "# fish flags", "# fishery selectivity",
                "# relative initial population", "# total populations scaling")
  structural = c("# The number of age classes", "# movement", "# diff_coffs",
                  "# xdiff_coffs", "# zdiff_coffs", "# Movement", "# region parameters",
                  "#Recruitment standard", "# Objective", "# The number of parameters",
                  "# Maximum magnitude", "# The parest_flags", "# Average fish mort",
                  "# The diffusion coefficients")
  for (k in names(src)) {
    if (any(startsWith(k, c(used_src, structural)))) next
    if (k %in% unlist(lapply(value_next, function(h) mfcl_next_key(src, mfcl_key(src, h))))) next
    v = unlist(lapply(src[[k]]$body, mfcl_num))
    if (length(v) && any(v != 0, na.rm = TRUE) && !grepl("^# [0-9.e+-]+ ", k))
      warning("Non-zero values in ", k, " are replaced by the template's")
  }
  out
}
