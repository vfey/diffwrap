test_that("a counts matrix file is read with samples in sample-sheet order", {
  quiet_log()
  si <- ex_samp_info()
  counts <- diff_expr_read_counts(ex_counts_file(), si)

  expect_true(is.matrix(counts))
  expect_equal(ncol(counts), N_SAMPLES)
  expect_equal(nrow(counts), N_GENES + N_SPECIAL)
  # columns must follow the sample sheet, not the file
  expect_identical(colnames(counts), as.character(si$SampleNames))
  expect_false(anyNA(counts))
  expect_true(all(counts >= 0))
})

test_that("an in-memory count matrix is accepted as well", {
  quiet_log()
  si  <- ex_samp_info()
  m   <- as.matrix(read.delim(ex_counts_file(), row.names = 1, check.names = FALSE))
  out <- diff_expr_read_counts(m, si)

  expect_true(is.matrix(out))
  expect_identical(colnames(out), as.character(si$SampleNames))
  expect_equal(nrow(out), N_GENES + N_SPECIAL)
})

test_that("every input mode returns the same type, so the class does not depend on the input", {
  # regression guard: read.delim()/read.table() return a data.frame while readDGE()/
  # getCounts() return a matrix, so the return type used to differ per input mode
  quiet_log()
  si <- ex_samp_info()

  from_file   <- diff_expr_read_counts(ex_counts_file(), si)
  from_matrix <- diff_expr_read_counts(
    as.matrix(read.delim(ex_counts_file(), row.names = 1, check.names = FALSE)), si)

  expect_true(is.matrix(from_file))
  expect_true(is.matrix(from_matrix))
  expect_identical(class(from_file), class(from_matrix))
  expect_true(is.numeric(from_file))
  expect_equal(from_file, from_matrix)
})

test_that("non-numeric count columns are rejected with a clear message", {
  quiet_log()
  si  <- ex_samp_info()
  bad <- read.delim(ex_counts_file(), row.names = 1, check.names = FALSE)
  # NOTE: merely as.character()-ing a column is not enough, because writing it out and
  # reading it back parses the numeric-looking strings straight back to integers. The
  # values have to be genuinely non-numeric for the column to arrive as character.
  bad[[1]] <- paste0("x", bad[[1]])
  f <- file.path(tempdir(), "bad_counts.tsv")
  write.table(bad, f, sep = "\t", quote = FALSE, col.names = NA)
  on.exit(unlink(f), add = TRUE)

  # must stop at the source rather than silently coercing the whole matrix to character
  expect_error(diff_expr_read_counts(f, si), "Counts must be numeric|Non-numeric column")
})

test_that("a well-formed count file really does arrive numeric", {
  # the companion to the test above: confirms the guard is not firing spuriously
  quiet_log()
  counts <- diff_expr_read_counts(ex_counts_file(), ex_samp_info())
  expect_true(is.numeric(counts))
  expect_false(is.character(counts))
})

test_that("unknown sample names are rejected rather than silently ignored", {
  quiet_log()
  si <- ex_samp_info()
  si$SampleNames <- factor(paste0("not_a_sample_", seq_len(nrow(si))))
  expect_error(diff_expr_read_counts(ex_counts_file(), si), "Sample names not found")
})

test_that("a subset of samples can be imported", {
  quiet_log()
  si     <- ex_samp_info()[1:4, ]
  counts <- diff_expr_read_counts(ex_counts_file(), si)
  expect_equal(ncol(counts), 4L)
  expect_identical(colnames(counts), as.character(si$SampleNames))
})

test_that("the miRSEQ path reads a CAP-miRSeq-style summary file", {
  # regression guard: this path previously referenced an undefined 'expression.raw' object
  quiet_log()
  si <- ex_samp_info()
  # minimal CAP-miRSeq summary: a Mature.miRNA column plus one count column per sample
  mir <- data.frame(Mature.miRNA = paste0("hsa-miR-", 1:5),
                    check.names = FALSE, stringsAsFactors = FALSE)
  for (s in as.character(si$SampleNames)) mir[[s]] <- sample(0:100, 5)
  f <- file.path(tempdir(), "mirseq_summary.tsv")
  write.table(mir, f, sep = "\t", quote = FALSE, row.names = FALSE)
  on.exit(unlink(f), add = TRUE)

  counts <- diff_expr_read_counts(f, si, miRSEQ = TRUE)
  expect_true(is.matrix(counts))
  expect_equal(nrow(counts), 5L)
  expect_equal(ncol(counts), nrow(si))
  expect_setequal(rownames(counts), paste0("hsa-miR-", 1:5))
})

test_that("filtering removes the htseq-count summary rows", {
  quiet_log()
  si     <- ex_samp_info()
  counts <- diff_expr_read_counts(ex_counts_file(), si)
  expect_equal(length(grep("^__", rownames(counts))), N_SPECIAL)   # present before

  filt <- diff_expr_filter_counts(counts, si, strict = TRUE)
  expect_equal(length(grep("^__", rownames(filt))), 0L)            # gone after
})

test_that("strict filtering removes low-expression genes but keeps the true positives", {
  quiet_log()
  si     <- ex_samp_info()
  counts <- diff_expr_read_counts(ex_counts_file(), si)
  filt   <- diff_expr_filter_counts(counts, si, strict = TRUE)

  expect_lt(nrow(filt), nrow(counts))
  expect_gt(nrow(filt), N_TRUE_DE)
  # every truly differentially expressed gene must survive filtering,
  # otherwise the downstream tests cannot detect them
  true_de <- sprintf("ENSG%011d", seq_len(N_TRUE_DE))
  expect_true(all(true_de %in% rownames(filt)))
})

test_that("non-strict filtering is more permissive than strict filtering", {
  quiet_log()
  si     <- ex_samp_info()
  counts <- diff_expr_read_counts(ex_counts_file(), si)

  strict <- diff_expr_filter_counts(counts, si, strict = TRUE)
  loose  <- diff_expr_filter_counts(counts, si, strict = FALSE)
  expect_gte(nrow(loose), nrow(strict))
})

test_that("min.samp controls the non-strict filter", {
  quiet_log()
  si     <- ex_samp_info()
  counts <- diff_expr_read_counts(ex_counts_file(), si)

  lenient <- diff_expr_filter_counts(counts, si, strict = FALSE, min.samp = 1)
  harsh   <- diff_expr_filter_counts(counts, si, strict = FALSE, min.samp = N_SAMPLES)
  expect_gte(nrow(lenient), nrow(harsh))
})
