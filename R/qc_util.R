# TODO: Add comment
#
# Author: vidal
###############################################################################


#' Function to do PCA using `stats::prcomp`
#' @param counts Counts matrix.
#' @param n Number of rows to be selected from the sorted variance matrix (by default, the top 500 rows are selected from the matrix sorted in decreasing order).
#' @param scale. A logical value passed to \code{prcomp}; should the variables be scaled to have unit variance?
#' @return An object of class \code{prcomp} as returned by \code{\link[stats]{prcomp}}, holding the
#'   principal component decomposition of the (optionally variance-filtered) count matrix.
#' @examples
#' si <- diff_expr_get_samp_info(diffwrap_samp_info, "SampleName", "Group")
#' counts <- diff_expr_filter_counts(diff_expr_read_counts(diffwrap_counts, si), si)
#' pca <- diff_expr_PCA(edgeR::cpm(counts, log = TRUE), n = 100)
#' pca$sdev[1:3]
#' @export
diff_expr_PCA <-
		function(counts, n=500, scale.=FALSE)
{
	select <- seq_len(nrow(counts))
	if (!is.null(n)) {
		dw_log("    Getting row-wise variances...\n")
		# per-row (per-gene) sample variance, base-R equivalent of genefilter::rowVars()
		Pvars <- rowSums((counts - rowMeans(counts))^2) / (ncol(counts) - 1L)
		dw_log("    Selecting", n, "rows with highest variances...\n")
		select <- order(Pvars, decreasing = TRUE)[seq_len(min(n, length(Pvars)))]
	}
	dw_log("    Calculating principal components...\n")
	PCA <- prcomp(t(counts[select, ]), scale.=scale.)
	return(PCA)
}

#' Function to calculate pseudo counts representing batch-corrected normalised but untransformed values
#' @param d Passed to \code{edgeR} functions: matrix of counts or a DGEList object.
#' @param design numeric design matrix
#' @param pairs \code{character}. Name of column with identifiers specifying paired samples for paired or other block designs, or batch effects.
#'   Defaults to "pairs" as this is the prefix added by 'diff_expr_make_design()'. When used outside the package's scope the user
#'   must supply the correct prefix.
#' @param disp \code{character}. one of "tagwise.dispersion", "trended.dispersion", "bin.dispersion"
#' @param do.cpm \code{logical}. should the pseudo counts be transformed to CPMs?
#' @return A \code{matrix} of pseudo counts in which the effect of the blocking variable has been
#'   removed; on the log2 counts per million scale if \code{do.cpm=TRUE} and on the count scale
#'   otherwise.
#' @examples
#' \donttest{
#' si <- diff_expr_get_samp_info(diffwrap_samp_info, "SampleName", "Group")
#' si$Subject <- diffwrap_samp_info$Subject[order(diffwrap_samp_info$SampleName)]
#' counts <- diff_expr_filter_counts(diff_expr_read_counts(diffwrap_counts, si), si)
#' groups <- stats::relevel(si$Groups, ref = "control")
#' d <- edgeR::calcNormFactors(edgeR::DGEList(counts, group = groups))
#' design <- diff_expr_make_design(si, groups, pairs = "Subject")
#' pc <- diff_expr_pseudo_counts(d = d, design = design, pairs = "pairs")
#' dim(pc)
#' }
#' @export
diff_expr_pseudo_counts <-
		function(d, design, pairs="pairs", disp="tagwise.dispersion", do.cpm=TRUE)
{
	dw_log("   --> Using pairs column identifier:", dQuote(pairs), "\n")
	dw_log("    Estimating dispersion...\n")
	disp.mat <- edgeR::estimateDisp(d, design)
	dw_log("    Fitting generalised linear model...\n")
	fit0 <- edgeR::glmFit(d, design, dispersion=disp.mat[[disp]])
	old.fitted <- fit0$fitted.values
	batch.coefs <- grep(pairs, colnames(design) )
	new.coefs <- fit0$unshrunk.coefficients
	dw_log("    Setting coefficients for blocking variable to 0...\n")
	new.coefs[, batch.coefs] <- 0
	dw_log("    Refitting coefficients...\n")
	new.fitted <- exp(new.coefs %*% t(design) + as.vector(fit0$offset))
	dw_log("    Getting pseudo-counts...\n")
	if (is(d, "DGEList")) {
	  dc <- d$counts
	} else if (is.data.frame(d)) {
	  dc <- try(as.matrix(d))
	  if (is(dc, "try-error")) stop("An input data frame must only consist of numeric values so it can be coerced to a (numeric) matrix.")
	} else if (is.matrix(d)) {
	  dc <- d
	}
	pseudo.counts <- edgeR::q2qnbinom(dc, old.fitted, new.fitted, dispersion=disp.mat[[disp]])
	if (any(pseudo.counts < 0)) {
	  dw_log("    #! Negative pseudo-counts detected...\n")
	  lpc <- length(which(pseudo.counts<0))/length(pseudo.counts)*100
	  if (lpc < 1) {
	    dw_log(paste0("     --> Less than 1% negative pseudo-counts (", lpc, "). Correcting...\n"))
	    dw_log("         Setting negative pseudo-counts to 0.1...\n")
	    pseudo.counts[which(pseudo.counts<0)] <- 0.1
	  } else {
	    dw_log("     --> More than 1% negative pseudo-counts. Skipping...\n")
	    stop("Too many negative pseudo-counts!")
	  }
	}
	if (do.cpm) {
		dw_log("    Getting CPMs...\n")
		pseudo.counts <- edgeR::cpm(pseudo.counts, log=TRUE, prior.count=3)
	}
	return(pseudo.counts)
}


