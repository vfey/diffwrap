test_that("the edgeR GLM path returns the documented objects", {
  quiet_log()
  ex     <- ex_dge()
  design <- diff_expr_make_design(samp.info = ex$samp.info, groups = ex$groups)

  fit.l <- diff_expr_fit(counts = ex$counts, d = ex$d, design = design,
                         do.voom = FALSE, quasi.likelihood = TRUE)

  expect_named(fit.l, c("d", "d2", "fit"))
  expect_s4_class_or_s3 <- function(x, cls) expect_true(methods::is(x, cls))
  expect_s4_class_or_s3(fit.l$d,  "DGEList")
  expect_s4_class_or_s3(fit.l$d2, "DGEList")
  expect_s4_class_or_s3(fit.l$fit, "DGEGLM")
  expect_true(!is.null(fit.l$d2$tagwise.dispersion))
})

test_that("quasi-likelihood and likelihood-ratio fits differ in class", {
  quiet_log()
  ex     <- ex_dge()
  design <- diff_expr_make_design(samp.info = ex$samp.info, groups = ex$groups)

  ql  <- diff_expr_fit(ex$counts, ex$d, design, do.voom = FALSE, quasi.likelihood = TRUE)
  lrt <- diff_expr_fit(ex$counts, ex$d, design, do.voom = FALSE, quasi.likelihood = FALSE)

  expect_true(methods::is(ql$fit,  "DGEGLM"))
  expect_true(methods::is(lrt$fit, "DGEGLM"))
  # glmQLFit carries quasi-likelihood specific slots that glmFit does not
  expect_true(!is.null(ql$fit$df.prior))
  expect_null(lrt$fit$df.prior)
})

test_that("the voom path returns the documented objects", {
  quiet_log()
  ex     <- ex_dge()
  design <- diff_expr_make_design(samp.info = ex$samp.info, groups = ex$groups)
  contr  <- diff_expr_make_contrasts(design = design, groups = ex$groups)

  fit.l <- diff_expr_fit(counts = ex$counts, d = ex$d, design = design,
                         do.voom = TRUE, voom.fun = "voomLmFit",
                         norm.method = "tmm", contrasts = contr)

  expect_named(fit.l, c("v", "fit", "fit2"))
  expect_true(methods::is(fit.l$fit,  "MArrayLM"))
  expect_true(methods::is(fit.l$fit2, "MArrayLM"))
  expect_true(!is.null(fit.l$fit2$p.value))   # eBayes has been applied
  expect_equal(nrow(fit.l$v$E), nrow(ex$counts))
})

test_that("norm.method must be supplied explicitly when calling diff_expr_fit() directly", {
  # NOTE: the default is norm.method = c("quantile", "tmm"); the function tests it with
  # if (norm.method == "tmm"), which is a length-2 condition and an error in R >= 4.2.
  # diffExpr() always match.arg()s before calling, so this only bites direct callers.
  quiet_log()
  ex     <- ex_dge()
  design <- diff_expr_make_design(samp.info = ex$samp.info, groups = ex$groups)

  expect_error(
    diff_expr_fit(ex$counts, ex$d, design, do.voom = TRUE, voom.fun = "voomLmFit"),
    regexp = "condition has length|the condition has length > 1"
  )
})

test_that("both engines recover the truly differentially expressed genes", {
  quiet_log()
  ex     <- ex_dge()
  design <- diff_expr_make_design(samp.info = ex$samp.info, groups = ex$groups)
  contr  <- diff_expr_make_contrasts(design = design, groups = ex$groups)

  # edgeR
  glm  <- diff_expr_fit(ex$counts, ex$d, design, do.voom = FALSE, quasi.likelihood = TRUE)
  tt_g <- edgeR::topTags(edgeR::glmQLFTest(glm$fit, contrast = contr[, 1]), n = Inf)
  top_g <- rownames(as.data.frame(tt_g))[seq_len(N_TRUE_DE)]

  # limma/voom
  vm   <- diff_expr_fit(ex$counts, ex$d, design, do.voom = TRUE, voom.fun = "voomLmFit",
                        norm.method = "tmm", contrasts = contr)
  tt_v <- limma::topTable(vm$fit2, coef = 1, number = Inf, sort.by = "P")
  top_v <- rownames(tt_v)[seq_len(N_TRUE_DE)]

  # with 60 true positives out of ~372 genes, a working pipeline should put most of
  # them in the top 60; the threshold is deliberately loose to avoid a flaky test
  expect_gt(sum(is_true_de(top_g)), 40)
  expect_gt(sum(is_true_de(top_v)), 40)

  # the two engines should broadly agree
  expect_gt(length(intersect(top_g, top_v)), 35)
})

test_that("the paired fit removes the subject effect", {
  quiet_log()
  ex <- ex_dge()
  ex$samp.info$Subject <- ex_samp_info_raw()$Subject[order(ex_samp_info_raw()$SampleName)]

  unpaired <- diff_expr_make_design(samp.info = ex$samp.info, groups = ex$groups)
  paired   <- diff_expr_make_design(samp.info = ex$samp.info, groups = ex$groups,
                                    pairs = "Subject")

  fit_u <- diff_expr_fit(ex$counts, ex$d, unpaired, do.voom = FALSE, quasi.likelihood = TRUE)
  fit_p <- diff_expr_fit(ex$counts, ex$d, paired,   do.voom = FALSE, quasi.likelihood = TRUE)

  # accounting for the subject offset should reduce residual dispersion
  expect_lt(median(fit_p$d2$tagwise.dispersion), median(fit_u$d2$tagwise.dispersion))
})
