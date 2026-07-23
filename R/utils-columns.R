# Internal helper for locating a single expected column by name pattern.
#
# Throughout the package, statistic columns (p-value, FDR, log fold-change, gene symbol,
# etc.) are located by grepping the (lower-cased) column names. A bare
#   nms[grep(pattern, tolower(nms))]
# fails cryptically when the pattern matches no column (returns character(0), so a later
# x[[col]] errors) or several columns (returns length > 1, so x[[col]] errors). This helper
# validates the match and reports a clear message instead.
#
# Author: vidal
###############################################################################

dw_find_col <- function(nms, pattern, what = "required") {
  hits <- nms[grep(pattern, tolower(nms))]
  if (length(hits) == 0L) {
    stop("Could not find the ", what, " column (search pattern: ", sQuote(pattern), ").",
         call. = FALSE)
  }
  if (length(hits) > 1L) {
    stop("Found several candidate ", what, " columns (",
         paste(sQuote(hits), collapse = ", "), "); expected exactly one.", call. = FALSE)
  }
  hits
}
