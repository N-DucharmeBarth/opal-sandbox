# nlminb() bound to one of the PORT builds made by build-port-variants.sh.
# R 4.6.2's nlminb.R is sourced into an environment (child of the stats
# namespace) in which C_port_ivset and C_port_nlminb point to the variant's
# shared library, so only the compiled PORT code differs from stats::nlminb().
port_nlminb <- function(variant, lib = file.path(Sys.getenv("NLMINB_BENCH_DIR"),
                                                 "same-compiler", "lib")) {
  so <- file.path(lib, paste0("port-", variant, ".so"))
  dll <- dyn.load(so, local = TRUE, now = TRUE)
  env <- new.env(parent = asNamespace("stats"))
  env$C_port_ivset <- getNativeSymbolInfo("port_ivset", dll)
  env$C_port_nlminb <- getNativeSymbolInfo("port_nlminb", dll)
  sys.source(file.path(lib, "nlminb.R"), envir = env)
  env$nlminb
}
