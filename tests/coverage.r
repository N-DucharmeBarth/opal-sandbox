# Coverage of the shared code (shared/) by the tests in tests/testthat
# Writes cobertura.xml in the project root for Codecov.
# Usage: Rscript tests/coverage.r (from anywhere; renv must not be activated)

# Define project paths; source files are listed relative to the project root
# so the report's paths match the repository's for Codecov
proj_dir = this.path::this.proj()
setwd(proj_dir)
source_files = list.files("shared", pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)

# covr sources the shared files into an environment, instruments them, then
# sources the "test files" there. Several test files skip() at top level,
# which only works under testthat, so the one test file covr sees runs
# test_dir() inside that environment; the tests' own source() calls go to the
# global environment and do not mask the instrumented functions.
driver = tempfile(fileext = ".R")
writeLines(sprintf(
    'testthat::test_dir("%s", env = environment(), reporter = "summary", stop_on_failure = FALSE)',
    file.path(proj_dir, "tests", "testthat")), driver)

cov = covr::file_coverage(source_files, driver)
print(cov)
covr::to_cobertura(cov, filename = "cobertura.xml")
