# Pure helpers with a well-defined return value: real assertions, no plotting device.

test_that("diff_expr_PCA returns a prcomp on the top-variance genes", {
  quiet_log()
  ex  <- ex_dge()
  lcpm <- edgeR::cpm(ex$d, log = TRUE)

  pca <- diff_expr_PCA(lcpm, n = 100)
  expect_s3_class(pca, "prcomp")
  # n = 100 restricts to the 100 most variable genes -> 100 rotation rows
  expect_equal(nrow(pca$rotation), 100L)
  expect_equal(nrow(pca$x), ncol(lcpm))          # one point per sample

  # n = NULL uses all genes
  pca_all <- diff_expr_PCA(lcpm, n = NULL)
  expect_equal(nrow(pca_all$rotation), nrow(lcpm))
})

test_that("quantile_breaks returns increasing, unique, in-range break points", {
  quiet_log()
  set.seed(1)
  xs <- rnorm(500)
  br <- quantile_breaks(xs, n = 20)

  expect_type(br, "double")
  expect_lte(length(br), 20L)
  expect_false(any(duplicated(br)))
  expect_identical(br, sort(br))
  expect_gte(min(br), min(xs))
  expect_lte(max(br), max(xs))
})

test_that("quantile_breaks collapses duplicate quantiles", {
  quiet_log()
  # a spiky vector produces many identical quantiles, which must be de-duplicated
  xs <- c(rep(0, 90), 1:10)
  br <- quantile_breaks(xs, n = 20)
  expect_false(any(duplicated(br)))
  expect_lt(length(br), 20L)
})

test_that("reorderFactors reorders levels without touching the values", {
  quiet_log()
  df <- data.frame(g = factor(c("b", "a", "c", "a")), stringsAsFactors = FALSE)
  out <- reorderFactors(df, column = "g", desired_level_order = c("c", "b", "a"))

  expect_identical(levels(out$g), c("c", "b", "a"))
  # the observed values themselves are unchanged, only the level order differs
  expect_identical(as.character(out$g), as.character(df$g))
})

test_that("make_pheatmap_anno_color builds one colour vector per annotation column", {
  quiet_log()
  clin <- data.frame(Group   = c("control", "control", "treated", "treated"),
                     Subject = c("P1", "P2", "P1", "P2"),
                     stringsAsFactors = FALSE)
  cols <- make_pheatmap_anno_color(clin)

  expect_type(cols, "list")
  expect_named(cols, c("Group", "Subject"))
  expect_length(cols$Group,   2L)   # two groups
  expect_length(cols$Subject, 2L)   # two subjects
  expect_named(cols$Group, c("control", "treated"))
  expect_true(all(grepl("^#", cols$Group)))   # colours are hex strings
})

test_that("diff_expr_pseudo_counts removes the blocking effect and keeps the matrix shape", {
  quiet_log()
  ex <- ex_dge()
  ex$samp.info$Subject <- ex_samp_info_raw()$Subject[order(ex_samp_info_raw()$SampleName)]

  design <- diff_expr_make_design(samp.info = ex$samp.info, groups = ex$groups,
                                  pairs = "Subject")

  pc <- diff_expr_pseudo_counts(d = ex$d, design = design, pairs = "pairs",
                                disp = "tagwise.dispersion", do.cpm = TRUE)

  expect_true(is.matrix(pc))
  expect_equal(dim(pc), dim(ex$d$counts))
  expect_false(anyNA(pc))
})
