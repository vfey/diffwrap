# TODO: Add comment
#
# Author: vidal
###############################################################################


#' Function to do PCA using `stats::prcomp`
#' @param counts Counts matrix.
#' @param n Number of rows to be selected from the sorted variance matrix (by default, the top 500 rows are selected from the matrix sorted in decreasing order).
#' @param scale. A logical value passed to \code{prcomp}; should the variables be scaled to have unit variance?
#' @export
diff_expr_PCA <-
		function(counts, n=500, scale.=FALSE)
{
	select <- 1:nrow(counts)
	if (!is.null(n)) {
		cat("    Getting row-wise variances...")
		Pvars <- genefilter::rowVars(counts)
		cat("done\n    Selecting", n, "rows with highest variances...")
		select <- order(Pvars, decreasing = TRUE)[seq_len(min(n, length(Pvars)))]
		cat("done\n")
	}
	cat("    Calculating principal components...")
	PCA <- prcomp(t(counts[select, ]), scale.=scale.)
	cat("done\n")
	return(PCA)
}

#' Function to calculate pseudo counts representing batch-corrected normalised but untransformed values
#' @param d Passed to \code{edgeR} functions: matrix of counts, or a DGEList object, or a SummarizedExperiment object.
#' @param design numeric design matrix
#' @param pairs Factor of identifiers specifying paired samples for paired or other block designs, or batch effects.
#' @param disp Character; one of "tagwise.dispersion", "trended.dispersion", "bin.dispersion"
#' @param do.cpm Logical; should the pseudo counts be transformed to CPMs?
#' @export
diff_expr_pseudo_counts <-
		function(d, design, pairs, disp="tagwise.dispersion", do.cpm=TRUE)
{
	cat("   --> Using pairs:", pairs, "\n")
	cat("    Estimating dispersion...")
	disp.mat <- edgeR::estimateDisp(d, design)
	cat("done\n    Fitting generalised linear model...")
	fit0 <- edgeR::glmFit(d, design, dispersion=disp.mat[[disp]])
	old.fitted <- fit0$fitted.values
	batch.coefs <- grep(pairs, colnames(design))
	new.coefs <- fit0$unshrunk.coefficients
	cat("done\n    Set coefficients for blocking variable to 0...")
	new.coefs[, batch.coefs] <- 0
	cat("done\n    Refit coefficients...")
	new.fitted <- exp(new.coefs %*% t(design) + as.vector(fit0$offset))
	cat("done\n    Getting pseudo-counts...")
	pseudo.counts <- edgeR::q2qnbinom(d$counts, old.fitted, new.fitted, dispersion=disp.mat[[disp]])
	cat("done\n")
	if (any(pseudo.counts < 0)) {
	  cat("    #! Negative pseudo-counts detected...\n")
	  lpc <- length(which(pseudo.counts<0))/length(pseudo.counts)*100
	  if (lpc < 1) {
	    cat(paste0("     --> Less than 1% negative pseudo-counts (", lpc, "). Correcting...\n"))
	    cat("         Setting negative pseudo-counts to 0.1...\n")
	    pseudo.counts[which(pseudo.counts<0)] <- 0.1
	    cat("done\n")
	  } else {
	    cat("     --> More than 1% negative pseudo-counts. Skipping...\n")
	    stop("Too many negative pseudo-counts!")
	  }
	}
	if (do.cpm) {
		cat("    Getting CPMs...")
		pseudo.counts <- edgeR::cpm(pseudo.counts, log=TRUE, prior.count=3)
		cat("done\n")
	}
	return(pseudo.counts)
}


