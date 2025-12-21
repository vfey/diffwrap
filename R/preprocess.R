# Preprocessing functions preparing necessary objects such as the design matrix or the contrast matrix
#
# Author: vidal
###############################################################################

utils::globalVariables("pairs.col")

#' Function to standardize samp.info sample information data frame
#' @param samp.info \code{data.frame}. samp.info object containing information of the project's sample sheet.
#' @param samples \code{character}. Name of the column in 'samp.info' containing sample names. If 'samp.info' is not supplied
#'     vector of sample names.
#' @param groups \code{character}. Name of the column in 'samp.info' containing grouping information. If 'samp.info' is not supplied
#'     vector of groups.
#' @param ellipse.mapping.groups Optional additional grouping for ellipse drawing. Overrides sample groups.
#'   Use for selective highlighting of user-defined sample groups.
#' @export
diff_expr_get_samp_info <-
		function(samp.info, samples, groups, ellipse.mapping.groups = NULL)
{
	if (!missing(samp.info)) {
		cat("Using user-provided sample information...\n")
		cat("  Renaming sample name column...\n")
		names(samp.info)[names(samp.info) == samples] <- "SampleNames"
		cat("  Renaming groups column...\n")
		names(samp.info)[names(samp.info) == groups] <- "Groups"
		if (!is.null(ellipse.mapping.groups)) {
			cat("  Renaming ellipse mapping column...\n")
			names(samp.info)[names(samp.info) == ellipse.mapping.groups] <- "Ellipse"
		}
	} else {
		samp.info <- data.frame(SampleNames=samples, Groups=groups, Ellipse={ if (is.null(ellipse.mapping.groups)) { groups } else { ellipse.mapping.groups } })
	}
	if (is.numeric(samp.info$Groups)) {
		samp.info$Groups <- paste0("group_", samp.info$Groups)
	}
	if (!is.null(ellipse.mapping.groups) && is.numeric(samp.info$Ellipse)) {
		samp.info$Ellipse <- paste0("ellipse_group_", samp.info$Ellipse)
	}
	cat("  Factorizing columns...\n")
	samp.info$SampleNames <- edgeR::dropEmptyLevels(as.factor(samp.info$SampleNames))
	samp.info$Groups <- edgeR::dropEmptyLevels(as.factor(samp.info$Groups))
	if (!is.null(ellipse.mapping.groups)) {
		samp.info$Ellipse <- edgeR::dropEmptyLevels(as.factor(samp.info$Ellipse))
	}
	cat("  Reordering by sample names...\n")
	samp.info <- samp.info[order(samp.info$SampleNames), ]
	return(samp.info)
}

#' Function to create design matrix
#' @param samp.info \code{data.frame}. samp.info object containing information of the project's sample sheet.
#' @param groups \code{character}. Name of the column in 'samp.info' containing grouping information. If 'samp.info' is not supplied
#'     vector of groups.
#' @param pairs \code{character}. Name of the column in 'samp.info' containing paired sample information.
#' @param block \code{logical}. Are the samples not independent? See Details section.
#' @param use_weights \code{logical}. Are sample-specific quality weights used? See Details section. (Placeholder for future versions)
#' @details
#' The 'block' argument is used to specify whether the comparisons are to be made within AND between subjects or in the case of
#' technical replicates, i.e., if the samples are not independent, in other words, correlated.
#' If sample-specific quality weights are to be estimated by means of 'voomWithQualityWeights()' or 'voomLmFit()' and 'sample.weights' set
#' to TRUE, 'use_weights' will be TRUE, enforcing a design matrix containing an 'intercept' column, i.e., where the columns reflect contrasts.
#' The choice of the design matrix type impacts the estimated weights due to the effect of the intercept on the residual variance per sample
#' in more complex designs, e.g., involving blocking factors, interactions or continuous covariates. The recommendation by the limma authors
#' is to use the default design matrix, i.e., with intercept. With simple designs, the type of design matrix is not relevant.
#' NOTE: This argument is not yet functional but a mere place-holder for future versions allowing for readily implemented more complex designs.
#' @seealso [model.matrix()]
#' @export
diff_expr_make_design <-
		function(samp.info, groups, pairs=NULL, block=FALSE, use_weights=FALSE)
{
	if (block || is.null(pairs)) {
		cat("Creating simple design matrix...\n")
		design <- model.matrix(~0+groups)
		colnames(design) <- gsub("groups", "", colnames(design))
	} else if (!is.null(pairs)) {
		cat("Creating design matrix for paired samples. Using column", sQuote(pairs), "as 'pairs' variable...\n")
		pairs.col <- pairs
		pairs <- samp.info[[pairs]]
		if (!is.numeric(pairs)) {
			pairs <- edgeR::dropEmptyLevels(as.factor(samp.info[[pairs.col]]))
			cat(" ", levels(pairs), "\n")
			design <- model.matrix(~pairs+groups)
		} else {
			design <- model.matrix(~groups+pairs)
		}
	}
	return(design)
}

#' Function to make contrast matrix
#' @param design Numeric design matrix.
#' @param groups \code{character}. Vector of group names.
#' @param pairs \code{character}. Name of the column in 'samp.info' containing paired sample information.
#' @param block \code{logical}. Are the samples not independent? See Details section.
#' @param contrasts Character vector specifying group name pairs to be compared in the format expected by
#'   \code{makeContrasts()}, i.e., "group2-group1".
#' @param user_design \code{logical}. Has the user provided a design matrix? This impacts the choosing of column names used for automatic contrasting.
#' @details
#' The 'block' argument is used to specify whether the comparisons are to be made within AND between subjects or in the case of
#' technical replicates, i.e., if the samples are not independent, in other words, correlated.
#'   @seealso [makeContrasts()]
#' @export
diff_expr_make_contrasts <-
		function(design, groups, pairs=NULL, block=FALSE, contrasts=NULL, user_design=FALSE)
{
	if (is.null(contrasts)) {
	  cat("Comparing all groups vs. all...\n")
		if (block || is.null(pairs)) {
			cat("  Using all levels...\n")
			n <- colnames(design)
		} else if (!is.null(pairs) || user_design) {
		  cat("  Extracting 'groups' levels...\n")
		  n <- grep("^groups.+", colnames(design), value=TRUE)
		  ic <- length(grep("Intercept", colnames(design)))
		  if (ic) {
		    cat("Intercept detected;", length(n), "comparison(s) found (inherent to the design)")
		    if (length(n)>1L) cat("  ", length(n)*(length(n)-1)/2, "comparisons to be made additionally")
		  } else {
		    cat("Non-intercept design;", length(n)*(length(n)-1)/2, "comparison(s) to be made")
		  }
		}
		if (length(n)==1L || length(grep("__", n))) {
			cat("  ! Skipping contrast matrix creation !\n")
			return(NULL)
		}
		contrasts <- unlist(lapply(n[1:(length(n)-1)], function(x) {
							contr <- paste(n[n!=x], x, sep="-")
							cyc <- which(n==x)
							contr[cyc:length(contr)]
						}))
		cat("  ", contrasts, "\n")
	} else {
	  cat("  Comparing selected groups", contrasts, "\n")
	  if (user_design) {
	    browser()

	  }
	  if (!is.null(pairs)) {
	    # the pairs argument being non-NULL means the design matrix has an intercept, but only in the context
	    # of this package; if a custom design matrix is provided that may not be the case
	    cat("  Extracting existing contrasts of interest from design matrix...\n")
	    cat("  - Complementing 'existing groups' vector...")
	    # Column names in the design matrix starting with "groups" represent comparisons to the baseline which is the grand mean
	    # of the first level of the groups factor and, in the context of this package, always the control, e.g., "normal",
	    # and represented by the "Intercept". In order to grep all existing comparisons, i.e., those inherent to the design,
	    # the name of the control, so, e.g., "groupsnormal", needs to be added in place of the Intercept.
	    n <- grep("^groups.+", colnames(design), value=TRUE)
	    n <- c(paste0("groups", levels(groups)[1]), n)
	    cat("done\n")
	    cat("  - Correcting 'contrasts' names...")
	    # The contrasts are provided based on the group names, so the word "groups" precedes column names in the design matrix.
	    # Hence, that word needs to be added to the provided contrasts.
	    cont <- unlist(strsplit(contrasts, "-"))
	    cont <- paste0("groups", cont)
	    cont <- unlist(lapply(seq(1, length(cont), 2), function(x) {
	      paste(cont[x], cont[x+1], sep="-")
	    }))
	    cat("done\n")
	    cat("  - Getting existing names...")
	    contrasts <- unlist(plyr::llply(cont, function(x) {
	      y <- unlist(strsplit(x, "-"))
	      if (n[1] != y[2]) {
	        cat("  -->", x, "is not inherent to the design matrix and will be extracted in addition to the 'group vs control' contrasts.\n")
	        return(x)
	      } else {
	        cat("  -->", x, "is inherent to the design matrix and will be extracted by default.\n")
	      }
	    }))
	    cat("done\n")
	}
	}
	cat("  Creating contrast matrix...\n")
	contrasts <- suppressWarnings(makeContrasts(contrasts=contrasts, levels=design))
	return(contrasts)
}


