test_that("the sample sheet is standardised to the package conventions", {
  quiet_log()
  si <- diff_expr_get_samp_info(ex_samp_info_raw(), samples = "SampleName", groups = "Group")

  expect_s3_class(si, "data.frame")
  expect_true(all(c("SampleNames", "Groups") %in% names(si)))
  expect_true(is.factor(si$SampleNames))
  expect_true(is.factor(si$Groups))
  expect_setequal(levels(si$Groups), c("control", "treated"))
  # rows are ordered by sample name
  expect_identical(as.character(si$SampleNames), sort(as.character(si$SampleNames)))
})

test_that("an ellipse grouping column is picked up when requested", {
  quiet_log()
  si <- diff_expr_get_samp_info(ex_samp_info_raw(), samples = "SampleName",
                                groups = "Group", ellipse.mapping.groups = "Subject")
  expect_true("Ellipse" %in% names(si))
  expect_true(is.factor(si$Ellipse))
  expect_setequal(levels(si$Ellipse), c("P1", "P2", "P3", "P4"))
})

test_that("numeric groups are made syntactically safe", {
  quiet_log()
  raw <- ex_samp_info_raw()
  raw$Group <- rep(c(1, 2), each = 4)
  si <- diff_expr_get_samp_info(raw, samples = "SampleName", groups = "Group")
  expect_true(all(grepl("^group_", levels(si$Groups))))
})

# --- design matrices -------------------------------------------------------

test_that("an unpaired design is a means model without intercept", {
  quiet_log()
  si <- ex_samp_info(); grp <- ex_groups(si)
  design <- diff_expr_make_design(samp.info = si, groups = grp)

  expect_true(is.matrix(design))
  expect_equal(nrow(design), N_SAMPLES)
  expect_false("(Intercept)" %in% colnames(design))
  expect_setequal(colnames(design), c("control", "treated"))
  # each row belongs to exactly one group
  expect_true(all(rowSums(design) == 1))
})

test_that("a paired design is an additive model with intercept", {
  quiet_log()
  si <- ex_samp_info(); grp <- ex_groups(si)
  si$Subject <- ex_samp_info_raw()$Subject[order(ex_samp_info_raw()$SampleName)]

  design <- diff_expr_make_design(samp.info = si, groups = grp, pairs = "Subject")

  expect_true("(Intercept)" %in% colnames(design))
  expect_true(any(grepl("^pairs", colnames(design))))   # subject covariates
  expect_true(any(grepl("^groups", colnames(design))))  # group effect vs control
  expect_equal(nrow(design), N_SAMPLES)
})

test_that("block = TRUE falls back to the means model even when pairs is given", {
  quiet_log()
  si <- ex_samp_info(); grp <- ex_groups(si)
  si$Subject <- ex_samp_info_raw()$Subject[order(ex_samp_info_raw()$SampleName)]

  design <- diff_expr_make_design(samp.info = si, groups = grp,
                                  pairs = "Subject", block = TRUE)
  expect_false("(Intercept)" %in% colnames(design))
  expect_setequal(colnames(design), c("control", "treated"))
})

test_that("the design matches what dw_resolve_mode() promises", {
  quiet_log()
  si <- ex_samp_info(); grp <- ex_groups(si)
  si$Subject <- ex_samp_info_raw()$Subject[order(ex_samp_info_raw()$SampleName)]

  for (spec in list(list(pairs = NULL,      block = FALSE, want = "means"),
                    list(pairs = "Subject", block = FALSE, want = "additive"),
                    list(pairs = "Subject", block = TRUE,  want = "means"))) {
    mode   <- dw_resolve_mode(spec$pairs, spec$block, do.voom = TRUE)
    design <- diff_expr_make_design(samp.info = si, groups = grp,
                                    pairs = spec$pairs, block = spec$block)
    has_intercept <- "(Intercept)" %in% colnames(design)
    expect_equal(mode$design, spec$want)
    expect_equal(has_intercept, identical(spec$want, "additive"),
                 info = paste("mode:", mode$role))
  }
})

# --- contrasts -------------------------------------------------------------

test_that("an unpaired two-group comparison yields a single contrast", {
  quiet_log()
  si <- ex_samp_info(); grp <- ex_groups(si)
  design    <- diff_expr_make_design(samp.info = si, groups = grp)
  contrasts <- diff_expr_make_contrasts(design = design, groups = grp)

  expect_true(is.matrix(contrasts))
  expect_equal(ncol(contrasts), 1L)
  expect_equal(colnames(contrasts), "treated-control")
  # a contrast column sums to zero
  expect_equal(unname(colSums(contrasts)), 0, tolerance = 1e-12)
})

test_that("a two-group paired design needs no contrast matrix", {
  quiet_log()
  si <- ex_samp_info(); grp <- ex_groups(si)
  si$Subject <- ex_samp_info_raw()$Subject[order(ex_samp_info_raw()$SampleName)]
  design <- diff_expr_make_design(samp.info = si, groups = grp, pairs = "Subject")

  # the group effect is inherent to the intercept design, so it is extracted by
  # coefficient name and no explicit contrast matrix is required
  contrasts <- diff_expr_make_contrasts(design = design, groups = grp, pairs = "Subject")
  expect_null(contrasts)
})

test_that("explicitly requested contrasts are honoured", {
  quiet_log()
  si <- ex_samp_info(); grp <- ex_groups(si)
  design    <- diff_expr_make_design(samp.info = si, groups = grp)
  contrasts <- diff_expr_make_contrasts(design = design, groups = grp,
                                        contrasts = "treated-control")
  expect_equal(colnames(contrasts), "treated-control")
})
