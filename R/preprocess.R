# Preprocessing functions preparing necessary objects such as the design matrix or the contrast matrix
# 
# Author: vidal
###############################################################################


#' Function to standardize `samp.info` sample information data frame
#' @export
diff_expr_get_samp_info <-
		function(samp.info, samples, groups, ellipse.mapping.groups=NULL)
{
	if (!is.null(samp.info)) {
		cat("Using user-provided sample information...\n")
		cat("  Renaming sample name column...\n")
		names(samp.info)[names(samp.info) == samples] <- "SampleNames"
		cat("  Renaming groups column...\n")
		names(samp.info)[names(samp.info) == groups] <- "Groups"
		if (!is.null(ellipse.mapping.groups)) {
			cat("  Renaming ellipse mapping column...\n")
			names(samp.info)[names(samp.info) == ellipse.mapping.groups] <- "Ellipse"
		}
		cat("  Factorizing columns... \n")
		samp.info$SampleNames <- edgeR::dropEmptyLevels(as.factor(samp.info$SampleNames))
		samp.info$Groups <- edgeR::dropEmptyLevels(as.factor(samp.info$Groups))
		if (!is.null(ellipse.mapping.groups)) {
			samp.info$Ellipse <- edgeR::dropEmptyLevels(as.factor(samp.info$Ellipse))
		}
	} else {
		samp.info <- data.frame(SampleNames=samples, Groups=groups, Ellipse={ if (is.null(ellipse.mapping.groups)) { groups } else { ellipse.mapping.groups } })
	}
	if (is.numeric(samp.info$Groups)) {
		samp.info$Groups <- paste0("group_", samp.info$Groups)
	}
	if (is.numeric(samp.info$Ellipse)) {
		samp.info$Ellipse <- paste0("ellipse_group_", samp.info$Ellipse)
	}
	cat("  Reordering by sample names...\n")
	samp.info <- samp.info[order(samp.info$SampleNames), ]
	return(samp.info)
}

#' Function to create design matrix
#' @export
diff_expr_make_design <-
		function(samp.info, groups, pairs=NULL, block=FALSE)
{
	if (block || is.null(pairs)) {
		cat("Creating simple design matrix...\n")
		design <- model.matrix(~0+groups)
		colnames(design) <- gsub("groups", "", colnames(design))
	} else if (!is.null(pairs)) {
		cat("Creating design matrix for paired samples. Using column", sQuote(pairs), "as 'pairs' variable...\n")
		pairs <- samp.info[[pairs]]
		if (!is.numeric(pairs)) {
			pairs <- dropEmptyLevels(as.factor(samp.info[[pairs]]))
			cat(" ", levels(pairs), "\n")
			design <- model.matrix(~pairs+groups)
		} else {
			design <- model.matrix(~groups+pairs)
		}
	}
	return(design)
}

#' Function to make contrast matrix
#' @export
diff_expr_make_contrasts <-
		function(design, pairs=NULL, block=FALSE, contrasts=NULL)
{
	if (is.null(contrasts)) {
		cat("Comparing all groups vs. all...\n")
		if (block || is.null(pairs)) {
			cat("  Using all levels...\n")
			n <- colnames(design)
		} else if (!is.null(pairs)) {
			cat("  Extracting 'groups' levels...\n")
			n <- grep("^groups.+", colnames(design), value=TRUE)
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
	}
	cat("  Creating contrast matrix...\n")
	contrasts <- suppressWarnings(makeContrasts(contrasts=contrasts, levels=design))
	return(contrasts)
}


