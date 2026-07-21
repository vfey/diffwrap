tmp_log <- function() file.path(tempdir(), paste0("dw_test_", as.integer(runif(1, 1, 1e6)), ".log"))

test_that("verbose is validated", {
  expect_error(dw_log_start("yes"), "must be TRUE, FALSE")
  expect_silent(dw_log_start(TRUE))
  expect_silent(dw_log_start(FALSE))
  expect_silent(dw_log_start("all"))
  dw_log_end()
})

test_that("verbose = TRUE shows steps on the console but not detail", {
  dw_log_start(TRUE)
  on.exit(dw_log_end(), add = TRUE)
  expect_message(dw_step("major\n"), "major")
  expect_silent(dw_log("detail\n"))
})

test_that("verbose = FALSE silences the console entirely", {
  dw_log_start(FALSE)
  on.exit(dw_log_end(), add = TRUE)
  expect_silent(dw_step("major\n"))
  expect_silent(dw_log("detail\n"))
})

test_that("verbose = 'all' mirrors detail to the console", {
  dw_log_start("all")
  on.exit(dw_log_end(), add = TRUE)
  expect_message(dw_step("major\n"), "major")
  expect_message(dw_log("detail\n"), "detail")
})

test_that("both tiers always reach the log file regardless of verbosity", {
  f <- tmp_log()
  dw_log_start(FALSE)
  dw_log_file(f)
  dw_step("major\n")
  dw_log("detail\n")
  dw_log_end()

  lines <- readLines(f)
  expect_true(any(grepl("^# diffwrap log", lines)))
  expect_true("major"  %in% lines)
  expect_true("detail" %in% lines)
})

test_that("lines logged before the file is opened are buffered and flushed into it", {
  f <- tmp_log()
  dw_log_start(FALSE)
  dw_log("early line\n")     # logged while no file exists yet
  dw_log_file(f)             # phase 2: buffer must be flushed here
  dw_log("late line\n")
  dw_log_end()

  lines <- readLines(f)
  expect_true("early line" %in% lines)
  expect_true("late line"  %in% lines)
  expect_lt(which(lines == "early line"), which(lines == "late line"))
})

test_that("call-site formatting is preserved exactly", {
  f <- tmp_log()
  dw_log_start(FALSE)
  dw_log_file(f)
  dw_log("  Using", sQuote("voom"), "with", 3, "groups\n")
  dw_log_end()

  expect_true(paste("  Using", sQuote("voom"), "with", 3, "groups") %in% readLines(f))
})

test_that("dw_log_obj captures printed objects", {
  f <- tmp_log()
  dw_log_start(FALSE)
  dw_log_file(f)
  dw_log_obj(factor(c("a", "b")))
  dw_log_end()

  expect_true(any(grepl("Levels", readLines(f))))
})

test_that("dw_log_end closes the connection and can be called twice", {
  f <- tmp_log()
  dw_log_start(FALSE)
  dw_log_file(f)
  dw_log_end()
  expect_silent(dw_log_end())
  # writing after close must not error - it falls back to the in-memory buffer
  expect_silent(dw_log("after close\n"))
  dw_log_end()
})
