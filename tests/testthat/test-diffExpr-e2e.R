# End-to-end runs of the main wrapper. These exercise the whole pipeline into a
# temporary directory, including the logging wiring and the required-out.dir behaviour.
# Enrichment is switched off so the run needs no annotation packages; biomart is off so
# it needs no network (gene symbols come from the offline convertId2()).

run_dir <- function() {
  d <- file.path(tempdir(), paste0("dw_e2e_", as.integer(runif(1, 1, 1e7))))
  dir.create(d, showWarnings = FALSE)
  d
}

test_that("out.dir is required and its omission is an error, not a working-directory write", {
  quiet_log()
  expect_error(
    suppressMessages(diffExpr(expr.dat  = ex_counts_file(),
                              samp.info = ex_samp_info_raw(),
                              samples   = "SampleName",
                              groups    = "Group",
                              control   = "control",
                              analysis.name = "demo",
                              enr.do = FALSE)),
    "out.dir"
  )
})

test_that("an unpaired run produces the documented outputs and a log file", {
  skip_on_cran()
  quiet_log()
  out.dir <- run_dir()
  on.exit(unlink(out.dir, recursive = TRUE), add = TRUE)

  res <- suppressMessages(
    diffExpr(expr.dat  = ex_counts_file(),
             samp.info = ex_samp_info_raw(),
             samples   = "SampleName",
             groups    = "Group",
             control   = "control",
             analysis.name = "demo",
             out.dir   = out.dir,
             enr.do = FALSE)
  )

  # return value
  expect_type(res, "list")
  expect_true("contrasts" %in% names(res))
  expect_true("treated-control" %in% names(res$contrasts))
  de <- res$contrasts[["treated-control"]]
  expect_s3_class(de, "data.frame")
  expect_true(any(grepl("gene_symbol", names(de))))

  # side effects: a per-contrast directory with a differential expression table
  contr_dir <- file.path(out.dir, "treated-control")
  expect_true(dir.exists(contr_dir))
  tsvs <- list.files(contr_dir, pattern = "differential_expression.*\\.tsv$")
  expect_gt(length(tsvs), 0L)

  # the run log is written below out.dir and is not empty
  logs <- list.files(out.dir, pattern = "\\.log$")
  expect_gt(length(logs), 0L)
  expect_gt(length(readLines(file.path(out.dir, logs[1]))), 0L)

  # the truly differentially expressed genes should dominate the top of the table
  pcol <- grep("^p\\.?val", tolower(names(de)))[1]
  de   <- de[order(de[[pcol]]), ]
  expect_gt(sum(is_true_de(de$ID[seq_len(50)])), 30)
})

test_that("verbose = FALSE is silent on the console but still writes the full log", {
  skip_on_cran()
  quiet_log()
  out.dir <- run_dir()
  on.exit(unlink(out.dir, recursive = TRUE), add = TRUE)

  # capture only our logger's channel: the major-step banners must not appear on the
  # console when verbose = FALSE. Unrelated warnings from edgeR/limma are ignored.
  msgs <- testthat::capture_messages(suppressWarnings(
    diffExpr(expr.dat  = ex_counts_file(),
             samp.info = ex_samp_info_raw(),
             samples   = "SampleName",
             groups    = "Group",
             control   = "control",
             analysis.name = "quiet",
             out.dir   = out.dir,
             verbose   = FALSE,
             enr.do = FALSE)
  ))
  expect_false(any(grepl("PREPROCESSING|STARTUP CHECKS|LINEAR MODELLING", msgs)))

  logs <- list.files(out.dir, pattern = "\\.log$", full.names = TRUE)
  expect_gt(length(logs), 0L)
  expect_gt(length(readLines(logs[1])), 10L)   # detail is logged even when quiet
})

test_that("a custom log.file location is honoured", {
  skip_on_cran()
  quiet_log()
  out.dir <- run_dir()
  on.exit(unlink(out.dir, recursive = TRUE), add = TRUE)
  lf <- file.path(out.dir, "custom_run.log")

  suppressMessages(
    diffExpr(expr.dat  = ex_counts_file(),
             samp.info = ex_samp_info_raw(),
             samples   = "SampleName",
             groups    = "Group",
             control   = "control",
             analysis.name = "demo",
             out.dir   = out.dir,
             log.file  = lf,
             enr.do = FALSE)
  )
  expect_true(file.exists(lf))
})

test_that("a paired run resolves to the fixed-effect mode and completes", {
  skip_on_cran()
  quiet_log()
  out.dir <- run_dir()
  on.exit(unlink(out.dir, recursive = TRUE), add = TRUE)

  res <- suppressMessages(
    diffExpr(expr.dat  = ex_counts_file(),
             samp.info = ex_samp_info_raw(),
             samples   = "SampleName",
             groups    = "Group",
             pairs     = "Subject",
             control   = "control",
             analysis.name = "paired",
             out.dir   = out.dir,
             enr.do = FALSE)
  )
  expect_type(res, "list")
  expect_true("contrasts" %in% names(res))
})
