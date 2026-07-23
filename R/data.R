# Documentation for the packaged example data
#
# Author: vidal
###############################################################################

#' Simulated RNA-Seq read counts for the package examples
#'
#' @description A small simulated count matrix used throughout the examples, tests and the
#'   vignette. The data are \strong{not} real measurements: they were generated from a negative
#'   binomial (gamma-Poisson) model so that the package can be demonstrated and tested without
#'   any external data or network access.
#'
#' @details The experiment consists of eight samples in two groups of four
#'   (\code{"control"} and \code{"treated"}). Each of the four subjects (\code{P1}-\code{P4})
#'   contributes one control and one treated sample, so the data support paired and blocked
#'   analyses as well as the simple unpaired comparison. A subject-specific offset is built in,
#'   giving the paired designs a real effect to remove.
#'
#'   Sixty of the 400 genes are truly differentially expressed (30 up-, 30 down-regulated) with
#'   absolute log2 fold changes between 1.2 and 3. The biological dispersion is 0.15
#'   (biological coefficient of variation of about 0.39) and library sizes are roughly 2-3
#'   million reads.
#'
#'   Forty genes are given deliberately low expression so that they straddle the default
#'   \code{strict} filtering threshold of \code{diff_expr_filter_counts()} (more than 5 counts
#'   per million in at least half of the samples, which corresponds to about 13 counts here);
#'   28 of them are removed by that filter.
#'
#'   The matrix additionally carries the five \command{htseq-count} summary rows
#'   (\code{__no_feature}, \code{__ambiguous}, \code{__too_low_aQual}, \code{__not_aligned} and
#'   \code{__alignment_not_unique}) that \code{diff_expr_filter_counts()} is expected to strip.
#'   They are retained on purpose so that the filtering step can be demonstrated.
#'
#'   Gene identifiers are well-formed but fictitious Ensembl gene IDs; they do not correspond
#'   to real genes and will not return annotation from 'biomart'.
#'
#' @format An \code{integer} matrix with 405 rows (400 genes plus 5 summary rows) and 8 columns
#'   (samples \code{S01}-\code{S08}). Row names are gene identifiers, column names sample names.
#' @source Simulated by \file{data-raw/make_example_data.R} from
#'   \file{inst/extdata/example_counts.tsv}.
#' @seealso \code{\link{diffwrap_samp_info}} for the matching sample sheet.
#' @examples
#' data(diffwrap_counts)
#' dim(diffwrap_counts)
#' head(diffwrap_counts[, 1:4])
#' # the htseq-count summary rows that get filtered out:
#' diffwrap_counts[grep("^__", rownames(diffwrap_counts)), 1:4]
"diffwrap_counts"

#' Sample sheet accompanying the simulated example counts
#'
#' @description The sample information table matching \code{\link{diffwrap_counts}}, in the
#'   layout expected by \code{diffExpr()} and \code{diff_expr_get_samp_info()}.
#'
#' @details The \code{Subject} column pairs each control sample with a treated sample from the
#'   same subject and is the column to pass to the \code{pairs} argument. \code{PlotName}
#'   supplies prettier labels for plotting and is the column to pass to
#'   \code{sample.plot.names}.
#'
#' @format A \code{data.frame} with 8 rows and 4 columns:
#' \describe{
#'   \item{SampleName}{\code{character}. Sample identifier, matching the column names of
#'     \code{\link{diffwrap_counts}}.}
#'   \item{Group}{\code{character}. Experimental group, \code{"control"} or \code{"treated"}.}
#'   \item{Subject}{\code{character}. Subject identifier (\code{P1}-\code{P4}); pass to
#'     \code{pairs} for paired or blocked designs.}
#'   \item{PlotName}{\code{character}. Human-readable sample label for plots.}
#' }
#' @source Simulated by \file{data-raw/make_example_data.R} from
#'   \file{inst/extdata/example_samp_info.tsv}.
#' @seealso \code{\link{diffwrap_counts}} for the matching count matrix.
#' @examples
#' data(diffwrap_samp_info)
#' diffwrap_samp_info
#' table(diffwrap_samp_info$Group, diffwrap_samp_info$Subject)
"diffwrap_samp_info"
