# Build the packaged example objects from the shipped plain-text files.
#
# The tab-separated files under inst/extdata are the SINGLE SOURCE OF TRUTH for the
# example data. They are shipped as-is so that examples and tests can also exercise the
# file-reading code path of diff_expr_read_counts(). The objects in data/ are derived
# from them by this script, so the two can never drift apart.
#
# The data are SIMULATED, not real measurements. They were generated from a negative
# binomial (gamma-Poisson) model with:
#   400 genes x 8 samples, two groups ("control", "treated") of four samples each
#   four subjects (P1-P4), each contributing one control and one treated sample, so
#     that the paired/blocked analysis modes have a real subject effect to remove
#   60 truly differentially expressed genes (30 up, 30 down), |log2FC| between 1.2 and 3
#   biological dispersion 0.15 (BCV ~0.39)
#   library sizes of roughly 2-3 million reads, so that the default 'strict' filter
#     (> 5 counts per million in at least half the samples) sits at ~13 counts and
#     actually discriminates: it removes 28 low-expression genes
#   five htseq-count style summary rows (__no_feature etc.) that
#     diff_expr_filter_counts() is expected to strip
#
# Re-run with:  source("data-raw/make_example_data.R")

counts_file    <- system.file("extdata", "example_counts.tsv",    package = "diffwrap")
samp_info_file <- system.file("extdata", "example_samp_info.tsv", package = "diffwrap")
if (!nzchar(counts_file)) {          # not installed yet: fall back to the source tree
  counts_file    <- "inst/extdata/example_counts.tsv"
  samp_info_file <- "inst/extdata/example_samp_info.tsv"
}

diffwrap_counts <- as.matrix(read.delim(counts_file, row.names = 1, check.names = FALSE))
storage.mode(diffwrap_counts) <- "integer"

diffwrap_samp_info <- read.delim(samp_info_file, stringsAsFactors = FALSE)

stopifnot(
  ncol(diffwrap_counts) == nrow(diffwrap_samp_info),
  identical(colnames(diffwrap_counts), diffwrap_samp_info$SampleName),
  !anyNA(diffwrap_counts),
  all(diffwrap_counts >= 0)
)

usethis::use_data(diffwrap_counts, diffwrap_samp_info, overwrite = TRUE)
