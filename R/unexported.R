# Unexported helpers
#
# Author: Vidal
################################################################################

#' Close a graphics device only if it is actually open
#'
#' @description Defensive wrapper around \code{dev.off()}. Plotting helpers occasionally close the
#'   device they were handed, and a bare \code{dev.off()} with nothing but the null device left
#'   fails with \dQuote{cannot shut down device 1 (the null device)}. This closes the requested
#'   device when it is still in \code{dev.list()} and is a silent no-op otherwise, so cleanup code
#'   and \code{on.exit()} handlers cannot error.
#' @param which (\code{integer}). Device number to close. Defaults to the current device.
#' @return Invisibly, \code{TRUE} if a device was closed and \code{FALSE} otherwise.
#' @keywords internal
dw_dev_off <- function(which = NULL) {
  open <- grDevices::dev.list()
  if (is.null(open)) return(invisible(FALSE))
  if (is.null(which)) which <- grDevices::dev.cur()
  if (!which %in% open) return(invisible(FALSE))   # never the null device: dev.list() excludes it
  grDevices::dev.off(which)
  invisible(TRUE)
}

#' Convert a numeric value of seconds to human-readable duration
#'
#' @param secs (\code{numeric}). run time in seconds
#' @return Character string with run time information.
#' @keywords internal
fmt_dur <- function(secs) {
  if (!is.finite(secs) || secs < 0) return(NA_character_)
  if (secs < 60) return(sprintf("%.1f sec", secs))
  s  <- round(secs)                 # whole seconds first -> no per-field rounding carry
  h  <- s %/% 3600
  m  <- (s %% 3600) %/% 60
  ss <- s %% 60
  if (h > 0) sprintf("%d h %02d min %02d sec", h, m, ss)
  else       sprintf("%d min %02d sec", m, ss)
}
