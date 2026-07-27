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
  if (secs < 60) sprintf("%.1f sec", secs)
  else if (secs < 3600) sprintf("%d min %02.0f sec", secs %/% 60, secs %% 60)
  else sprintf("%d h %02d min %02.0f sec", secs %/% 3600, (secs %% 3600) %/% 60, secs %% 60)
}
