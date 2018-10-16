# Functions for performing differential expression analysis
# 
# Author: vidal
###############################################################################


#' Function to compute linear model fit and optionally apply 'voom' beforehand
#' @export
diff_expr_fit <-
		function(counts, d, design, do.voom=TRUE, voom.fun=voom, norm.method=c("quantile", "tmm"), bayes.trend=FALSE, bayes.robust=FALSE, pairs=NULL,
				pairs_col=NULL, block=FALSE, contrasts=NULL, disp="tagwise.dispersion")
{
	if (do.voom) {
		if (norm.method=="tmm") {
			cts <- d
			norm.method="none"
		} else {
			cts <- counts
		}
		v <- voom.fun(cts, design=design, plot=FALSE, normalize.method=norm.method)
		
		# For paired samples or block designs:
		if (block && !is.null(pairs)) {
			if (!is.factor(pairs)) {
				stop("For block designs 'pairs' needs to be a factor.")
			}
			cat("  Samples are paired. Using column", sQuote(pairs_col), "as block variable...\n")
			cat("    ", levels(pairs), "\n")
			cat("  Calculating consensus correlation...\n")
			corfit <- duplicateCorrelation(v$E , design, block=pairs)
			cat("  Fitting linear model using consensus correlation...\n")
			fit <- lmFit(v, design, block=pairs, correlation=corfit$consensus)
		} else {
			cat("  Fitting linear model...\n")
			fit <- lmFit(v, design)
		}
		if (!is.null(contrasts)) {
			if (length(grep("^Intercept$", rownames(contrasts)))) {
				rownames(contrasts)[grep("^Intercept$", rownames(contrasts))] <- "(Intercept)"
			}
			cat("  Fitting contrasts...\n")
			fit2 <- contrasts.fit(fit, contrasts)
			cat("  eBayes fit...\n")
			fit2 <- eBayes(fit2, trend = bayes.trend, robust = bayes.robust)
		} else {
			cat("  eBayes fit...\n")
			fit2 <- eBayes(fit, trend = bayes.trend, robust = bayes.robust)
		}
		return(list(v=v, fit=fit, fit2=fit2))
	} else {
		## Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR)-adjusted likelihood
		cat("  Estimating dispersion values...\n")
		d2 <- estimateGLMTrendedDisp(d, design)
		d2 <- estimateGLMTagwiseDisp(d2, design)
		# Given the design matrix and dispersion estimates, fit a GLM to each feature
		cat("  Fitting linear model...\n")
		fit <- glmFit(d2, design, dispersion=d2[[disp]])
		return(list(d=d, d2=d2, fit=fit))
	}
}



#' Helper function to generate an output table with only most relevant columns
##' @param samp.name.and.group.key \code{data.frame}. Sample names (matching with annotated.normcnt cols) as rownames, groups in 1st column
#' @export
diffr_expr_generate_cleaned_de_table_output <-
  function(contrast, annotated.normcnt, out.dir=".", samp.name.and.group.key, analysis.name=NULL, filtered.lists = TRUE,  fdr.thr=0.05, logfc.thr=1 )
{
    cat("Generating a cleaned version of expression table...\n")
    DE.out <- paste(analysis.name, contrast, "differential_expression_clean.tsv", sep="_")
  
    group1 <- unlist(strsplit(contrast, split="-",fixed = TRUE))[1]
    group2 <- unlist(strsplit(contrast, split="-",fixed = TRUE))[2]
   
    contrast.samples <- rownames(samp.name.and.group.key)[samp.name.and.group.key$group == group1 | samp.name.and.group.key$group == group2]
    not.contrast.samples <- rownames(samp.name.and.group.key)[!(samp.name.and.group.key$group == group1 | samp.name.and.group.key$group == group2)]
    
    annotated.normcnt <- annotated.normcnt[,!colnames(annotated.normcnt) %in% not.contrast.samples]
    
    cleaned.colnames <- gsub("ReadsPerGene.out", "", colnames(annotated.normcnt))
    #print(cleaned.colnames)
    colnames(annotated.normcnt) <- cleaned.colnames 
    cat("   Samples outside the contrast were removed...\n")
    
    if (filtered.lists) {
      cat("     Filtering the DE list to be saved by fdr <", fdr.thr, "and logfc > ", logfc.thr , "...\n")
      fdr.col <- names(annotated.normcnt)[grep("^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", tolower(names(annotated.normcnt)))]
      fc.col <- names(annotated.normcnt)[grep("^logfc$|fold$", tolower(names(annotated.normcnt)))]
      
      #print(fdr.col)
      #print(fc.col)
      
      nrow.before.filtering <- nrow(annotated.normcnt)
      
      annotated.normcnt <- annotated.normcnt[(annotated.normcnt[[fdr.col]] < fdr.thr &
                                 abs(annotated.normcnt[[fc.col]]) > logfc.thr),]
      
      nrow.removed <- nrow.before.filtering - nrow(annotated.normcnt)
      cat("       ", nrow.removed, "unsignifcant genes were filtered out...\n")
      cat("        The dimensions of the table to be saved:", dim(annotated.normcnt), "...\n")
      
      
      DE.out <- paste("filtered", DE.out, sep="_")
    }
    
    cat("   Saving cleanded DE-list to", DE.out, "...\n")
    write.table(annotated.normcnt, file.path(out.dir, DE.out), sep="\t", quote=FALSE, row.names=FALSE)

}
        
#' Function to extract contrasts and generate top tables and plots
#' @export
diff_expr_extract_contrasts <-
		function(contrasts=NULL, fit, fit2=NULL, normcnt, out.l, do.voom=TRUE, out.dir=".",
				analysis.name=NULL, biomart=FALSE, biom.data.set="hsapiens_gene_ensembl", biom.mart=c("ensembl", "snp", "funcgen", "vega", "pride", "plants"),
				host="www.ensembl.org", biom.filter="ensembl_gene_id", biom.attributes=c("ensembl_gene_id","hgnc_symbol","description"), sym.col="hgnc_symbol",
				rm.dups=FALSE, p.thr=0.05, fdr.thr=0.05, logfc.thr=1, numlab=15, point.lab=TRUE, font.size=5, plots=TRUE, lists=TRUE, filtered.lists = TRUE)
{
	if (!length(grep("contrasts", names(out.l)))) {
		out.l$contrasts <- list()
	}
	if (is.null(contrasts)) {
		cat("  Extracting 'groups' levels from design...\n")
		cn <- grep("^groups.+", colnames(fit$design), value=TRUE)
	} else if (is.matrix(contrasts) && !is.null(colnames(contrasts))) {
		# if do.voom=TRUE contrasts have been calculated already in the fit function and only need to be extracted
		# in the case of paired samples, contrasts are inherent to the design and only need to be extracted
		# otherwise, we need to feed the contrast matrix to glmLRT to instruct it to calculate the desired contrasts
		# see the different behaviour of glmLRT below for do.voom=FALSE
		cn <- colnames(contrasts)
	} else {
		# contrasts is a character vector
		# this will only be called for the non-contrasted columns of the design matrix for paired samples
		# which means contrasts are inherent to the design and only need to be extracted
		cn <- contrasts
	}

	out.l$MAplots <- list()
	out.l$volcanoPlots <- list()
	for (contr in cn) {
		cat("Calculating differential expression for", contr, "\n")
		if (do.voom) {
			cat("Generating output table of differentially expressed features...\n")
			tt <- topTable(fit2, coef=contr, number=nrow(fit2), sort.by="P")
		} else {
			# Perform a likelihood ratio test, specifying the contrast of interest.
			# Needs to be done per-contrast as otherwise it will perform the LR test of all contrasts together while we want this done individually.
			cat("Performing likelihood ratio test on this contrast...\n")
			if (is.matrix(contrasts)) {
				# calculating the contrast
				de <- glmLRT(fit, contrast = contrasts[, contr])
			} else {
				# only extracting the contrast (for paired samples)
				de <- glmLRT(fit, coef=contr)
			}
			# Use the topTags function to present a tabular summary of the differential expression statistics (note that topTags
			# operates on the output of exactTest or glmLRT, but only the latter is shown here)
			cat("Generating tablular summary of differential expression statistics ('top table')\n")
			tt <- topTags(de, n = nrow(de))
#			} else {
#				stop("Need contrast matrix as input.")
#			}
		}
		if (is.null(contrasts)) {
			fitcnt <- fit$coefficients[, contr]
			addinfo <- cbind(normcnt, fitcnt)
			if (do.voom) {
				names(addinfo)[names(addinfo) %in% "fitcnt"] <- "Fitted Coefficients (log2)"
			} else {
				names(addinfo)[names(addinfo) %in% "fitcnt"] <- "Fitted Coefficients (ln)"
			}
		} else if (is.matrix(contrasts)) {
			fitcnt <- fit$coefficients[, rownames(contrasts)[contrasts[, contr]!=0]]
			addinfo <- cbind(normcnt, fitcnt)
		} else {
			addinfo <- normcnt
		}
		d3 <- merge(addinfo, tt, by="row.names", sort=FALSE)
		cat("Renaming ID column...\n")
		names(d3)[1] <- "ID"
		if (biomart) {
			d3 <- diff_expr_biomart(d3, biom.data.set, biom.mart, host, biom.filter, biom.attributes, sym.col, rm.dups)
			id.col <- names(d3)[names(d3) %in% biom.filter]
		} else {
			syms <- convertId2(as.character(d3$ID))
			syms <- data.frame(ID=names(syms), gene_symbol=as.character(syms), stringsAsFactors=FALSE)
			d3 <- merge(syms, d3, by="ID", all.y=TRUE, all.x=FALSE, sort=TRUE)
			if (any(d3$gene_symbol=="") || any(is.na(d3$gene_symbol))) {
				cat("  Replacing", length(which(d3$gene_symbol=="" | is.na(d3$gene_symbol))), "missing Gene Symbols by Ensembl IDs...\n")
				d3$gene_symbol[d3$gene_symbol=="" | is.na(d3$gene_symbol)] <- as.character(d3$ID[d3$gene_symbol=="" | is.na(d3$gene_symbol)])
			}
			id.col <- "ID"
		}
		cat("Setting unique row names...\n")
		rownames(d3) <- make.unique(as.character(d3[[id.col]]))
		out.l$contrasts[[contr]] <- d3
		if (lists) {
			DE.out <- paste(analysis.name, contr, "differential_expression.tsv", sep="_")
			cat("Saving list to", DE.out, "...\n")
			pv.col <- names(d3)[grep("^p\\.{0,1}val[e-u]{0,2}$", tolower(names(d3)))]
			cat("The result table is ordered in increasing order by column ", pv.col, "...\n")
			d3 = d3[order(d3[[pv.col]]),]
			write.table(d3, file.path(out.dir, DE.out), sep="\t", quote=FALSE, row.names=FALSE)
		
			diffr_expr_generate_cleaned_de_table_output(contrast=contr, annotated.normcnt=d3, samp.name.and.group.key = fit$samples, out.dir,
			                                            analysis.name, filtered.lists = TRUE, fdr.thr=fdr.thr, logfc.thr=logfc.thr)
			}
		
		# Create a graphical summary, such as an M (log-fold change) versus A (log-average expression) plot, here showing the
		# genes selected as differentially expressed (with a 5% false discovery rate)
		if (plots) {
			cat("Plotting...\n")
			pdf(file.path(out.dir, paste(analysis.name, contr, "_plots.pdf", sep="_")), width=11, height=8.5)
			par(mar = c(6,6,5,3))
			cat(" MA-plot...\n")
			out.l$MAplots[[contr]] <- diff_expr_ma_plot(d3, contr, id.col, sym.col, p.thr, fdr.thr, logfc.thr, numlab, out.dir, analysis.name, point.lab, biom.attributes, font.size, lists)
			
			## Volcano plot
			cat(" Volcano plot...\n")
			#The following three lines would probably be wise to do using regular expression
			volcano.name = gsub(".", " ", contr, fixed=TRUE)
			volcano.name = gsub("_", " ", volcano.name, fixed=TRUE)
			volcano.name = gsub("-", " vs. ", volcano.name, fixed=TRUE)
			out.l$volcanoPlots[[contr]] <- diff_expr_volcano_plot(d3, id.col, sym.col="gene_symbol", main=volcano.name, p.thr=p.thr, fdr.thr=fdr.thr, logfc.thr=logfc.thr, numlab=numlab, point.lab=point.lab)
			
#png(paste(out.dir,"/",analysis.name,".Pvalue_distribution.png",sep=""),width=1280,height=960,res=150)
			## Histogram of P-value distribution
			cat(" Dendrogram plot...\n")
			diff_expr_pval_hist_plot(d3)
			dev.off()
			
			cat("done\n")
		}
	}
	return(out.l)
}

#' Function to retrieve additional information from biomart
#' @export
diff_expr_biomart <-
		function(d3, biom.data.set="hsapiens_gene_ensembl", biom.mart=c("ensembl", "snp", "funcgen", "vega", "pride", "plants"),
				host="www.ensembl.org", biom.filter="ensembl_gene_id", biom.attributes=c("ensembl_gene_id","hgnc_symbol","description"),
				sym.col="hgnc_symbol", rm.dups=FALSE)
{
	gene.lab <- convert.bm(d3, "ID", biom.data.set, biom.mart, host, biom.filter, biom.attributes, sym.col, rm.dups)
	names(gene.lab)[names(gene.lab)==sym.col] <- "gene_symbol"
	cat("  Extended annotation:\n")
	biom.attributes[biom.attributes==sym.col] <- "gene_symbol"
	if (length(d3$ID)>8) {
		print(gene.lab[1:8, biom.attributes])
		cat("_truncated_ (", length(d3$ID), "features)\n")
	} else {
		print(gene.lab[, biom.attributes])
	}
	return(gene.lab)
}


