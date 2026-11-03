# Readers for MULTIFAN-CL (MFCL) par and plot-report files, in base R (no
# FLR4MFCL). Both files are sequences of "# header" lines, each followed by
# lines of numbers, so both are read as named sections. Section names are the
# trimmed header plus its occurrence, e.g. "# fish flags [1]" and
# "# fish flags [2]"; unnamed "#" headers are "# [1]", "# [2]", ...
#
# Used by the MFCL-to-opal converter (mfcl-convert.R) and by case-study
# scripts that edit MFCL inputs (e.g. case-studies/swpo-alb).

#' Read a text file, dropping Windows line endings.
mfcl_read_lines = function(path) {
    sub("\r$", "", readLines(path, warn = FALSE))
}

#' Numbers on one line (empty for a blank line).
mfcl_num = function(x) {
    x = trimws(x)
    if (!nzchar(x)) return(numeric(0))
    as.numeric(strsplit(x, "[[:space:]]+")[[1]])
}

#' Format numbers as MFCL writes them: integers plain, others with 14 digits.
mfcl_fmt = function(x) {
    ifelse(x == round(x) & abs(x) < 1e9, format(x, scientific = FALSE, trim = TRUE),
           formatC(x, digits = 14, format = "e"))
}

#' Split a header/body file into named sections.
mfcl_sections = function(path) {
    x = mfcl_read_lines(path)
    is_h = startsWith(x, "#")
    if (!is_h[1]) stop("First line of ", path, " is not a header")
    grp = cumsum(is_h)
    secs = lapply(split(seq_along(x), grp), function(ix) list(header = x[ix[1]], body = x[ix[-1]]))
    key = trimws(vapply(secs, `[[`, "", "header"))
    occ = stats::ave(seq_along(key), key, FUN = seq_along)
    names(secs) = paste0(key, " [", occ, "]")
    secs
}

#' Read a par file into sections (see mfcl_sections()).
mfcl_read_par = function(path) mfcl_sections(path)

#' Write sections back to a par file.
mfcl_write_par = function(secs, path) {
    out = unlist(lapply(secs, function(s) c(s$header, s$body)), use.names = FALSE)
    writeLines(out, path)
}

#' Key of the section whose name starts with `pattern`, at occurrence `occ`.
mfcl_key = function(secs, pattern, occ = 1L) {
    exact = paste0(trimws(pattern), " [", occ, "]")
    if (exact %in% names(secs)) return(exact)
    k = names(secs)[startsWith(names(secs), pattern)]
    k = k[endsWith(k, paste0("[", occ, "]"))]
    if (length(k) != 1L) stop("Section not unique: ", pattern, " [", occ, "]")
    k
}

#' Key of the section after `key` (the values of named par headers sit in the
#' unnamed "#" section that follows them).
mfcl_next_key = function(secs, key) names(secs)[match(key, names(secs)) + 1L]

#' Non-blank body lines as a numeric matrix (rows of equal length).
mfcl_matrix = function(body) {
    rows = lapply(body[nzchar(trimws(body))], mfcl_num)
    if (!length(rows)) return(matrix(numeric(0), 0, 0))
    do.call(rbind, rows)
}

#' Fish flags as a matrix [fishery, flag] (occurrence 1 is the current set).
mfcl_fish_flags = function(par, occ = 1L) mfcl_matrix(par[[mfcl_key(par, "# fish flags", occ)]]$body)

#' A flag vector ("# The parest_flags" or "# age flags").
mfcl_flags = function(par, pattern, occ = 1L) mfcl_num(par[[mfcl_key(par, pattern, occ)]]$body[1])

#' Read a plot report (plot-*.par.rep) into sections.
mfcl_read_rep = function(path) mfcl_sections(path)

#' A plot-report block as a matrix.
#'
#' @param header start of the section header, e.g. "# Mean lengths at age".
#' @param region for blocks split by "# Region r" sub-headers, the region.
#' @param species for blocks split by "# Species s" sub-headers.
mfcl_rep = function(rep, header, region = NULL, species = NULL, occ = 1L) {
    k = mfcl_key(rep, header, occ)
    sub_hdr = if (!is.null(region)) paste("# Region", region) else
        if (!is.null(species)) paste("# Species", species) else NULL
    if (!is.null(sub_hdr)) {
        i = match(k, names(rep))
        j = i + 1L
        while (j <= length(rep) && !startsWith(names(rep)[j], sub_hdr)) {
            if (j > i + 50L) stop("Sub-header ", sub_hdr, " not found after ", header)
            j = j + 1L
        }
        return(mfcl_matrix(rep[[j]]$body))
    }
    mfcl_matrix(rep[[k]]$body)
}

#' Ragged plot-report block (rows of different lengths) as a list.
mfcl_rep_rows = function(rep, header, occ = 1L) {
    b = rep[[mfcl_key(rep, header, occ)]]$body
    lapply(b[nzchar(trimws(b))], mfcl_num)
}
