# Unexported helpers
#
# Author: Vidal
################################################################################

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
