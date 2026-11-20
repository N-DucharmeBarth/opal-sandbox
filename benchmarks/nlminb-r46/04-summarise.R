# Collect results/raw/*.rds into CSV files read by nlminb-r46.qmd.
#   results/fits.csv      one row per fit (case, start, R version)
#   results/paired.csv    one row per case and start, both R versions side by side
#   results/controls-{per-call,pure-r,replay}.csv  per-call cost of fn and gr,
#                         nlminb on a pure-R function, trajectory replay checks
source(file.path(Sys.getenv("NLMINB_BENCH_DIR"), "bench-helpers.R"))
raw <- file.path(bench_dir(), "results", "raw")
res <- file.path(bench_dir(), "results")

files <- dir(raw, pattern = "^(opaka|io_yft)-.*\\.rds$", full.names = TRUE)
fits <- lapply(files, readRDS)
rows <- do.call(rbind, lapply(fits, `[[`, "row"))
rows <- rows[order(rows$case, rows$start_id, rows$r_version), ]
write.csv(rows, file.path(res, "fits.csv"), row.names = FALSE)

key <- sub("-[0-9]+\\.rds$", "", basename(files))
ver <- sub(".*-([0-9]+)\\.rds$", "\\1", basename(files))
paired <- do.call(rbind, lapply(split(seq_along(files), key), function(i) {
  if (length(i) != 2L) return(NULL)
  i <- i[order(ver[i])]
  a <- fits[[i[1]]]; b <- fits[[i[2]]]
  data.frame(
    case = a$row$case, start_id = a$row$start_id,
    same_start_fn = identical(a$start_fn, b$start_fn),
    same_start_gr = identical(a$start_gr, b$start_gr),
    elapsed_452 = a$row$elapsed_s, elapsed_462 = b$row$elapsed_s,
    time_ratio = b$row$elapsed_s / a$row$elapsed_s,
    fn_evals_452 = a$row$fn_evals, fn_evals_462 = b$row$fn_evals,
    gr_evals_452 = a$row$gr_evals, gr_evals_462 = b$row$gr_evals,
    iter_452 = a$row$iterations, iter_462 = b$row$iterations,
    conv_452 = a$row$conv_pass2, conv_462 = b$row$conv_pass2,
    obj_452 = a$row$objective, obj_462 = b$row$objective,
    dobj_mle_452 = a$row$objective - a$row$mle_objective,
    dobj_mle_462 = b$row$objective - b$row$mle_objective,
    maxgrad_452 = a$row$max_abs_grad, maxgrad_462 = b$row$max_abs_grad,
    max_abs_par_diff = max(abs(a$par - b$par)))
}))
paired <- paired[order(paired$case, paired$start_id), ]
write.csv(paired, file.path(res, "paired.csv"), row.names = FALSE)

ctl <- lapply(dir(raw, pattern = "^controls-", full.names = TRUE), readRDS)
names(ctl) <- vapply(ctl, `[[`, "", "r_version")
per_call <- do.call(rbind, lapply(ctl, `[[`, "per_call"))
pure_r <- do.call(rbind, lapply(ctl, `[[`, "pure_r"))
replay <- do.call(rbind, lapply(names(ctl[[1]]$replay), function(m) {
  r <- lapply(ctl, function(x) x$replay[[m]])
  data.frame(model = m, events = length(r[[1]]$values),
             identical_between_versions = mean(mapply(identical, r[[1]]$values, r[[2]]$values)),
             identical_to_stored_452 = mean(r[[1]]$same_as_stored),
             identical_to_stored_462 = mean(r[[2]]$same_as_stored))
}))
write.csv(per_call, file.path(res, "controls-per-call.csv"), row.names = FALSE)
write.csv(pure_r, file.path(res, "controls-pure-r.csv"), row.names = FALSE)
write.csv(replay, file.path(res, "controls-replay.csv"), row.names = FALSE)
print(paired); print(per_call); print(pure_r); print(replay)
