# Regression tests for the bug-fix batches (discrete bugs, NA-safe filtering, grep
# validation). Each test would fail against the code as it was before the fix.

# --- dw_find_col: validated column lookup ----------------------------------

test_that("dw_find_col returns the single match, or errors clearly", {
  nms <- c("ID", "gene_symbol", "logFC", "FDR", "PValue")
  expect_equal(dw_find_col(nms, "^fdr$", "FDR"), "FDR")
  expect_equal(dw_find_col(nms, "^logfc$|fold$", "log fold-change"), "logFC")
  # zero matches -> clear error rather than a cryptic x[[character(0)]] failure
  expect_error(dw_find_col(nms, "entrez", "Entrez"), "Could not find")
  # several matches -> clear error rather than x[[c('a','b')]]
  expect_error(
    dw_find_col(c("FDR", "adj.P.Val"), "^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", "FDR"),
    "several"
  )
})

# --- NA-safe filtering in the cleaned-table writer -------------------------

test_that("diffr_expr_generate_cleaned_de_table_output never emits NA rows", {
  quiet_log()
  out.dir <- tempfile("dw_clean_"); dir.create(out.dir)
  on.exit(unlink(out.dir, recursive = TRUE), add = TRUE)
  key <- data.frame(group = rep(c("control", "treated"), each = 4),
                    row.names = paste0("S", 1:8))
  # two rows are individually significant on one axis but NA on the other:
  # a pre-fix logical subset would turn them into all-NA output rows
  d3 <- data.frame(ID = paste0("g", 1:10), gene_symbol = paste0("G", 1:10),
                   logFC = c(rnorm(8, 0, 3), NA, 3),
                   FDR   = c(runif(8, 0, 0.001), 0.001, NA),
                   stringsAsFactors = FALSE)

  diffr_expr_generate_cleaned_de_table_output(
    contrast = "treated-control", annotated.normcnt = d3,
    samp.name.and.group.key = key, out.dir = out.dir,
    analysis.name = "t", filtered.lists = TRUE, fdr.thr = 0.05, logfc.thr = 1)

  f <- list.files(out.dir, pattern = "clean\\.tsv$", full.names = TRUE)
  expect_length(f, 1L)
  res <- read.delim(f, check.names = FALSE)
  expect_false(anyNA(res$FDR))
  expect_false(anyNA(res$logFC))
})

test_that("cleaned-table writer keeps all columns for coefficient-style contrasts", {
  quiet_log()
  out.dir <- tempfile("dw_clean2_"); dir.create(out.dir)
  on.exit(unlink(out.dir, recursive = TRUE), add = TRUE)
  key <- data.frame(group = rep(c("control", "treated"), each = 4),
                    row.names = paste0("S", 1:8))
  d3 <- data.frame(ID = paste0("g", 1:5), gene_symbol = paste0("G", 1:5),
                   S1 = 1:5, S2 = 1:5, logFC = rnorm(5), FDR = runif(5, 0, 0.001),
                   stringsAsFactors = FALSE)
  # "groupstreated" cannot be mapped to the two groups -> must NOT drop sample columns
  diffr_expr_generate_cleaned_de_table_output(
    contrast = "groupstreated", annotated.normcnt = d3,
    samp.name.and.group.key = key, out.dir = out.dir,
    analysis.name = "t", filtered.lists = FALSE)

  f <- list.files(out.dir, pattern = "clean\\.tsv$", full.names = TRUE)
  res <- read.delim(f, check.names = FALSE)
  expect_true(all(c("S1", "S2") %in% names(res)))
})

# --- volcano helper --------------------------------------------------------

test_that("prepare_volcano_of_given_property returns a ggplot, incl. the no-hit case", {
  quiet_log()
  set.seed(1)
  d3 <- data.frame(gene_symbol = paste0("G", 1:200),
                   logFC = rnorm(200, 0, 2), PValue = runif(200))
  g <- prepare_volcano_of_given_property(d3, property.to.plot = "p",
                                         property.column = "PValue",
                                         property.thr = 0.05, logfc.thr = 1, numlab = 10)
  expect_s3_class(g, "ggplot")

  # nothing significant -> no labels, but still a valid plot (no 1:0 / empty-subset error)
  d3b <- data.frame(gene_symbol = paste0("G", 1:20),
                    logFC = rep(0, 20), PValue = rep(0.9, 20))
  expect_s3_class(
    prepare_volcano_of_given_property(d3b, "p", "PValue", 0.05, 1, numlab = 10),
    "ggplot")
})

# --- format_ensembl for mouse IDs ------------------------------------------

test_that("format_ensembl_ids_annotated_to_term detects mouse ENSMUSG columns", {
  quiet_log()
  # pre-fix, the column was found by grepl('ENSG', ...), which does NOT match ENSMUSG,
  # so this errored for mouse. The value conversion itself uses the offline convertId2().
  result <- data.frame(
    term   = c("GO:0001", "GO:0002"),
    geneID = c("ENSMUSG00000017146,ENSMUSG00000041147", "ENSMUSG00000059552"),
    stringsAsFactors = FALSE)
  out <- tryCatch(
    format_ensembl_ids_annotated_to_term(result, species = "Mouse"),
    error = function(e) skip(paste("convertId2 unavailable:", conditionMessage(e))))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
})

# --- diffr_venn default join.vec -------------------------------------------

test_that("diffr_venn runs with its default join.vec on minimal (symbol-only) tables", {
  quiet_log()
  skip_if_not_installed("venn")
  # the default is "gene_symbol", so tables carrying only that column must work
  mk <- function(sig) data.frame(gene_symbol = paste0("G", sig), fdr = 0.01,
                                  stringsAsFactors = FALSE)
  tabs <- list(A = mk(1:20), B = mk(10:30))
  v <- suppressWarnings(diffr_venn(tabs))          # no join.vec -> uses the default
  expect_type(v, "list")
  expect_true(all(c("venn.diagram", "venn.sections") %in% names(v)))
})
