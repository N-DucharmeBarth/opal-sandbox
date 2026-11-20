# Collect same-compiler/results/raw/*.rds into
#   same-compiler/results/variants.csv  one row per case, start and variant
#   same-compiler/results/match.csv     per case and variant: starts whose fn
#       and gr counts and final objective are identical to the R 4.5.2 and the
#       R 4.6.2 stats::nlminb fits of the main benchmark (results/raw)
dir <- Sys.getenv("NLMINB_BENCH_DIR")
res <- file.path(dir, "same-compiler", "results")
raw <- lapply(dir(file.path(res, "raw"), full.names = TRUE), readRDS)
v <- do.call(rbind, lapply(raw, `[[`, "rows"))
v <- v[order(v$case, v$start_id, v$variant), ]
write.csv(v, file.path(res, "variants.csv"), row.names = FALSE)

# Full-precision rows of the main benchmark (fits.csv rounds to 15 digits)
main <- do.call(rbind, lapply(dir(file.path(dir, "results", "raw"),
                                  pattern = "^(opaka|io_yft)-", full.names = TRUE),
                              function(f) readRDS(f)$row))
key <- function(d) paste(d$case, d$start_id)
same_as <- function(d, rv) {
  m <- main[main$r_version == rv, ]
  m <- m[match(key(d), key(m)), ]
  d$fn_evals == m$fn_evals & d$gr_evals == m$gr_evals &
    (d$objective == m$objective | (is.na(d$objective) & is.na(m$objective)))
}
match_tab <- do.call(rbind, lapply(split(v, list(v$case, v$variant), drop = TRUE), function(d)
  data.frame(case = d$case[1], variant = d$variant[1], starts = nrow(d),
             same_as_r452 = sum(same_as(d, "R 4.5.2"), na.rm = TRUE),
             same_as_r462 = sum(same_as(d, "R 4.6.2"), na.rm = TRUE))))
match_tab <- match_tab[order(match_tab$case, match_tab$variant), ]
write.csv(match_tab, file.path(res, "match.csv"), row.names = FALSE)
print(match_tab, row.names = FALSE)
