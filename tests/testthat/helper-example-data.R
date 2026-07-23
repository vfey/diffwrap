# Shared fixtures for the test suite.
#
# The example data are read from inst/extdata rather than from data/, so that the tests
# work directly after installation without data-raw/make_example_data.R having been run,
# and so that the file-reading code path is exercised at the same time.
#
# Ground truth of the simulated data (see data-raw/make_example_data.R):
#   400 genes, 8 samples, groups "control" (4) and "treated" (4)
#   subjects P1-P4, each contributing one control and one treated sample
#   genes ENSG00000000001 - ENSG00000000060 are truly differentially expressed
#   5 htseq-count summary rows (__no_feature etc.) that filtering must remove
###############################################################################

N_TRUE_DE   <- 60L
N_GENES     <- 400L
N_SPECIAL   <- 5L
N_SAMPLES   <- 8L

ex_counts_file <- function() {
  f <- system.file("extdata", "example_counts.tsv", package = "diffwrap")
  if (!nzchar(f)) skip("example count file not found")
  f
}

ex_samp_info_file <- function() {
  f <- system.file("extdata", "example_samp_info.tsv", package = "diffwrap")
  if (!nzchar(f)) skip("example sample sheet not found")
  f
}

# raw sample sheet, as a user would supply it
ex_samp_info_raw <- function() {
  read.delim(ex_samp_info_file(), stringsAsFactors = FALSE)
}

# sample sheet standardised to the package conventions
ex_samp_info <- function(pairs = FALSE) {
  diff_expr_get_samp_info(ex_samp_info_raw(), samples = "SampleName", groups = "Group")
}

ex_groups <- function(si = ex_samp_info()) {
  stats::relevel(si$Groups, ref = "control")
}

# filtered counts + normalised DGEList, the usual starting point for modelling
ex_dge <- function(strict = TRUE) {
  si     <- ex_samp_info()
  counts <- diff_expr_read_counts(ex_counts_file(), si)
  counts <- diff_expr_filter_counts(counts, si, strict = strict)
  d      <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts, group = ex_groups(si)))
  list(samp.info = si, counts = counts, d = d, groups = ex_groups(si))
}

# is a gene id one of the truly differentially expressed ones?
is_true_de <- function(ids) {
  n <- suppressWarnings(as.integer(sub("^ENSG0*", "", ids)))
  !is.na(n) & n <= N_TRUE_DE
}

# quieten the package logger for the duration of a test
quiet_log <- function() {
  dw_log_start(FALSE)
  withr::defer(dw_log_end(), envir = parent.frame())
}
