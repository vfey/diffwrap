# Smoke tests for the plotting functions: they assert that a plot is produced
# (the right class, or simply no error) without inspecting its visual content, which
# would be brittle against ggplot2 changes. Anything drawn to a device is captured in
# a throw-away PDF so that no Rplots.pdf is left behind.

with_null_pdf <- function(code) {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  force(code)
}

test_that("diff_expr_ggplot_mds returns a ggplot", {
  quiet_log()
  ex   <- ex_dge()
  lcpm <- edgeR::cpm(ex$d, log = TRUE)
  g <- diff_expr_ggplot_mds(counts = lcpm,
                            samp.name = colnames(lcpm),
                            groups = ex$groups)
  expect_s3_class(g, "ggplot")
})

test_that("the PCA ggplot helpers return ggplot objects", {
  quiet_log()
  ex   <- ex_dge()
  lcpm <- edgeR::cpm(ex$d, log = TRUE)
  pca  <- diff_expr_PCA(lcpm, n = 100)

  gb <- diff_expr_PCA_ggbiplot(PCA = pca, groups = ex$groups)
  expect_s3_class(gb, "ggplot")

  # samp.name must be NULL, NA, or a named vector matching rownames(PCA$x); pass NULL
  # so the row names are used, which is the common case
  gp <- diff_expr_PCA_ggplot(PCA = pca, samp.name = NULL,
                             groups = ex$groups, do.plot = FALSE)
  expect_s3_class(gp, "ggplot")
})

test_that("diff_expr_dendro_plot draws without error", {
  quiet_log()
  ex   <- ex_dge()
  lcpm <- edgeR::cpm(ex$d, log = TRUE)
  expect_no_error(with_null_pdf(
    diff_expr_dendro_plot(counts = lcpm, groups = ex$groups)
  ))
})

test_that("diff_expr_mds_plot draws without error", {
  quiet_log()
  ex <- ex_dge()
  expect_no_error(with_null_pdf(
    diff_expr_mds_plot(d = ex$d, groups = ex$groups, do.pdf = FALSE,
                       out.dir = tempdir())
  ))
})

test_that("diff_expr_pval_hist_plot draws a histogram from a p-value column", {
  quiet_log()
  d3 <- data.frame(ID = paste0("g", 1:100),
                   PValue = runif(100),
                   logFC  = rnorm(100))
  h <- with_null_pdf(diff_expr_pval_hist_plot(d3))
  expect_s3_class(h, "histogram")
})
