test_that("the five intended scenarios resolve correctly", {
  expect_equal(dw_resolve_mode(NULL,   FALSE, FALSE)[c("role","engine","design")],
               list(role = "none",  engine = "edger", design = "means"))
  expect_equal(dw_resolve_mode(NULL,   FALSE, TRUE )[c("role","engine","design")],
               list(role = "none",  engine = "limma", design = "means"))
  expect_equal(dw_resolve_mode("subj", FALSE, FALSE)[c("role","engine","design")],
               list(role = "fixed", engine = "edger", design = "additive"))
  expect_equal(dw_resolve_mode("subj", FALSE, TRUE )[c("role","engine","design")],
               list(role = "fixed", engine = "limma", design = "additive"))
  expect_equal(dw_resolve_mode("subj", TRUE,  TRUE )[c("role","engine","correlation")],
               list(role = "block", engine = "limma", correlation = TRUE))
})

test_that("block forces the limma engine even when do.voom = FALSE", {
  m <- dw_resolve_mode("subj", block = TRUE, do.voom = FALSE)
  expect_equal(m$engine, "limma")
  expect_true(m$do.voom)
  expect_match(m$notes, "voom enforced")
})

test_that("block without pairs is an error, not a silent no-op", {
  expect_error(dw_resolve_mode(NULL, block = TRUE, do.voom = TRUE), "requires a 'pairs' column")
})
