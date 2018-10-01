# TODO: Add comment
# 
# Author: vidal
###############################################################################


#' Function to do PCA using `stats::prcomp`
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
#' @export
diff_expr_pseudo_counts <-
		function(design, d, pairs, disp="tagwise.dispersion", do.cpm=TRUE)
{
	cat("    Estimating dispersion...")
	disp.mat <- estimateDisp(d, design)
	cat("done\n    Fitting generalised linear model...")
	fit0 <- glmFit(d, design, dispersion=disp.mat[[disp]])
	old.fitted <- fit0$fitted.values
	batch.coefs <- grep(pairs, colnames(design))
	new.coefs <- fit0$unshrunk.coefficients
	cat("done\n    Set coefficients for blocking variable to 0...")
	new.coefs[, batch.coefs] <- 0
	cat("done\n    Refit coefficients...")
	new.fitted <- exp(new.coefs %*% t(design) + as.vector(fit0$offset))
	cat("done\n    Getting pseudo-counts...")
	pseudo.counts <- q2qnbinom(d$counts, old.fitted, new.fitted, dispersion=disp.mat[[disp]])
	cat("done\n")
	if (do.cpm) {
		cat("    Getting CPMs...")
		pseudo.counts <- cpm(pseudo.counts, log=TRUE, prior.count=3)
		cat("done\n")
	}
	return(pseudo.counts)
}


