# Internal logging: one always-on log file + a tiered, suppressible console.
#
#   dw_step()    major step  -> console (verbose TRUE/"all") + log file
#   dw_log()     detail      -> log file (console only when verbose = "all")
#   dw_log_obj() object dump -> as dw_log(), via print()
#
# Console output goes to stdout via cat(), gated by 'verbose' (so verbose = FALSE silences it).
# stdout is used rather than message()/stderr because dependencies loaded during a run (e.g.
# VennDiagram/futile.logger) can disturb the stderr stream.
# Both helpers capture cat()/print() internally, so call sites keep their exact original formatting.
#
# Author: vidal
###############################################################################

.dw_state <- new.env(parent = emptyenv())
.dw_state$con     <- NULL          # open file connection, or NULL
.dw_state$verbose <- TRUE          # TRUE | FALSE | "all"
.dw_state$buffer  <- character()   # lines logged before the file was opened

# --- lifecycle ----------------------------------------------------------------

# Phase 1: called at the very top of diffExpr(), before out.dir exists.
dw_log_start <- function(verbose = TRUE) {
  if (!(isTRUE(verbose) || isFALSE(verbose) || identical(verbose, "all"))) {
    stop("'verbose' must be TRUE, FALSE or \"all\".", call. = FALSE)
  }
  dw_log_end()
  .dw_state$verbose <- verbose
  .dw_state$buffer  <- character()
  invisible(NULL)
}

# Phase 2: called once out.dir exists; opens the file and flushes the buffer.
dw_log_file <- function(path, append = FALSE) {
  con <- file(path, open = if (append) "at" else "wt")
  .dw_state$con <- con
  writeLines(sprintf("# diffwrap log | %s",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S")), con)
  if (length(.dw_state$buffer)) {
    writeLines(.dw_state$buffer, con)
    .dw_state$buffer <- character()
  }
  flush(con)
  invisible(path)
}

dw_log_end <- function() {
  if (!is.null(.dw_state$con)) {
    try(close(.dw_state$con), silent = TRUE)
    .dw_state$con <- NULL
  }
  .dw_state$buffer <- character()
  invisible(NULL)
}

# --- emitters -----------------------------------------------------------------

.dw_capture <- function(...) utils::capture.output(cat(...))

.dw_emit <- function(txt, console) {
  if (!length(txt)) return(invisible(NULL))
  if (!is.null(.dw_state$con)) {
    writeLines(txt, .dw_state$con)
    flush(.dw_state$con)          # survive an interrupted run
  } else {
    .dw_state$buffer <- c(.dw_state$buffer, txt)
  }
  # Console channel is stdout via cat(), NOT message()/stderr: some dependencies loaded in a run
  # (notably VennDiagram/futile.logger) disturb the stderr/message stream whereas stdout is
  # diffwrap's own channel and always shows. Still gated by 'verbose'.
  if (console) cat(paste(txt, collapse = "\n"), "\n", sep = "")
  invisible(NULL)
}

# detail: log file only (console only in "all" mode)
dw_log <- function(...) {
  .dw_emit(.dw_capture(...), console = identical(.dw_state$verbose, "all"))
}

# major step: console + log file
dw_step <- function(...) {
  .dw_emit(.dw_capture(...),
           console = isTRUE(.dw_state$verbose) || identical(.dw_state$verbose, "all"))
}

# object dump (replaces bare print() calls used for diagnostics)
dw_log_obj <- function(x, ...) {
  .dw_emit(utils::capture.output(print(x, ...)),
           console = identical(.dw_state$verbose, "all"))
}
