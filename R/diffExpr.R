#' Main wrapper for executing the entire pipeline from reading in expression data such as count
#' files to producing text files and graphs
#' @description \command{diffExpr} is a conveniency wrapper performing all steps automatically.
#'    Most sub-functions are exported and can be called by the user, as well, if desired.
#'    These functions may be applicable to different kinds of data/input, rely, however,
#'    on the conventions set for this package.
#' @param expr.dat \code{character} or \code{list}. String or vector or list of input file paths, or matrix of count values
#' @param samp.info \code{data.frame}. samp.info object containing information of the project's sample sheet
#' @param control \code{character}. Name of the control group
#' @param design \code{matrix}. design matrix
#' @param samples \code{character}. Name of the column in 'samp.info' containing sample names. If 'samp.info' is not supplied
#'     vector of sample names.
#' @param sample.plot.names \code{character}. Optional name of a column with "nice" sample names for plotting.
#' @param groups \code{character}. Name of the column in 'samp.info' containing grouping information. If 'samp.info' is not supplied
#'     vector of groups.
#' @param pairs \code{character}. Name of the column in 'samp.info' containing paired sample information.
#' @param block \code{logical}. Are the the samples not independent? See Details section.
#' @param contrasts \code{character}. Vector of contrasts to be made. If not provided, all possible contrasts will be made.
#'   This specifies group name pairs to be compared in the format expected by \code{makeContrasts()}, i.e., "group2-group1".
#' @param bayes.trend \code{logical}. Should an intensity-trend be allowed for the prior variance? Passed to 'limma::eBayes'.
#' @param bayes.robust \code{logical}. Should the estimation of df.prior and var.prior be robustified against outlier sample variances? Passed to 'limma::eBayes'.
#' @param quasi.likelihood Logical; should quasi-likelihood methods be used? See \emph{Details} section.
#'   Defaults to NA, which will determine the method based on the number of replicate samples.
#'   If more than 4 replicates are present, the likelihood ratio test is used, otherwise the quasi-likelihood methods.
#'   If \code{TRUE}, then the quasi-likelihood methods are used, if \code{FALSE}, then the likelihood ratio test is used.
#' @param sym.col \code{character}. Name of the column in the query result with gene symbols
#' @param ellipse \code{logical}. Should an ellipse be plotted around samples belonging to the same sample group? Defaults to \code{TRUE}.
#' @param ellipse.mapping.groups \code{character} The name of the column in 'samp.info' with group names for ellipse drawing. If \code{NULL} (default)
#'     will use the \code{groups} column. If 'samp.info' is not supplied vector of groups.
#' @param label.samples \code{logical}. Should points in appropriate QC plots be labelled. So far, applies only to PCA ggplot. Defaults to \code{TRUE}.
#' @param geom.point.size \code{numeric}. Size of points in appropriate QC plots. So far, applies only to PCA ggplot. Defaults to 2.
#' @param label.font.size \code{numeric}. Font size used for point labels in appropriate QC plots. So far, applies to PCA ggplot and M-A plots. Defaults to 5.
#' @param plot.ellipse.legend \code{logical}. Should a legend be addded for ellipses in PCA plots? NA, the default, includes
#'     if any aesthetics are mapped. FALSE never includes, and TRUE always includes. It can also be a named logical vector to finely select
#'     the aesthetics to display.
#' @param do.enrichment \code{logical}. Whether or not to call enrichment wrapper. Defaults to \code{TRUE}.
#' @param enrichment.plot \code{logical}. Whether or not to draw a network plot for the enrichment results. Defaults to \code{FALSE}.
#' @param enrichment.plot.fdr.thr \code{numeric}. FDR threshold used in the enrichment plot. This may be useful to tweak to produce a more informative plot.
#' Defaults to 0.05.
#' @param enrichment.plot.logfc.thr \code{numeric}. FC threshold on the log2-scale used in the enrichment plot. Defaults to 1.
#' @param enrichment.plot.num.terms \code{integer}. Number of terms shown in the plot. Defaults to 5.
#' @param enrichment.methods \code{character}. One or more of the following: c("clusterProfilerGO", "clusterProfilerKEGG", "gProfileR", "topGO"). By default, uses them all.
#' @param dry.run \code{logical}. If \code{TRUE}, the function will not create any output files or directories.
#' @param ... \code{ANY}. Additional arguments passed to functions.
#' @param out.dir \code{character}. Path to the output directory. If not provided, a new directory will be created in the project directory.
#' @param project.dir \code{character}. Path to the project directory. Defaults to the current working directory.
#' @param analysis.name \code{character}. Name of the analysis. If not provided, a default name will be generated.
#' @param biomart \code{logical}. Should the biomart be used for gene annotation? Defaults to \code{FALSE}.
#' @param biom.data.set \code{character}. The biomart dataset to be used. Defaults to "hsapiens_gene_ensembl".
#' @param biom.mart \code{character}. The biomart to be used.
#' @param host \code{character}. The host to be used for the biomart. Defaults to "www.ensembl.org".
#' @param biom.filter \code{character}. The biomart filter to be used. Defaults to "ensembl_gene_id".
#' @param biom.attributes \code{character}. The biomart attributes to be used. Defaults to c("ensembl_gene_id", "hgnc_symbol", "description", "entrezgene_id").
#' @param biom.force.ensg \code{logical}. Should Ensembl Gene IDs be checked for (and stripped of) ensemble version numbers? Defaults to \code{FALSE}.
#' @param biom.cache \code{character}. Path name giving the location of the cache \command{getBM()} uses if \code{use.cache=TRUE}. Defaults to the value in the \emph{BIOMART_CACHE} environment variable.
#' @param use.cache (\code{logical}). Should \command{getBM()} use the cache? Defaults to \code{TRUE} as in the \command{getBM()} function and is passed on to that.
#' @param rm.dups \code{logical}. Should duplicates be removed from the output of the biomart request? Defaults to \code{FALSE}.
#' @param p.thr \code{numeric}. Threshold for p-values. Defaults to 0.05.
#' @param fdr.thr \code{numeric}. Threshold for FDR values. Defaults to 0.05.
#' @param logfc.thr \code{numeric}. Threshold for fold-change values on the log2-scale. Defaults to 1.
#' @param numlab \code{numeric}. Maximum number of point labels to be shown in the plot. This overrides/limits
#'   values calculated by any thresholds. Defaults to 15.
#' @param point.lab \code{logical}. Should point labels be shown in the plot? Defaults to \code{TRUE}.
#' @param heatmap.topn \code{numeric}. Number of top values to be plotted. Defaults to 100.
#' @param heatmap.split.expr \code{logical}. Should the top up- and top down-regulated genes be displayed at equal numbers (50/50),
#' if they meet the significance threshold (regardless of the actual significance)? Defaults to \code{FALSE}.
#' @param min.samp \code{integer}. Number of samples in which a feature needs to be covered by at least one read per million.
#'   Defaults to the size of the smallest group of replicates. See \emph{details}.
#' @param strict \code{logical}. For miRNA analysis: only keep a miRNA if there are > 5 reads per million in at least half of the samples?
#' @param disp \code{character}. The dispersion method to be used. Defaults to "gene".
#' @param do.voom \code{logical}. Should the voom function be used? Defaults to \code{FALSE}.
#' @param voom.fun \code{character}. The voom function to be used. Should be one of 'limma::voom', 'limma::voomWithQualityWeights' or
#' 'edgeR::voomLmFit'. Defaults to \code{edgeR::voomLmFit}. See 'Details'.
#' @param use_weights \code{logical}. Should sample-specific quality weights be estimated?
#' @param norm.method \code{character}. The normalisation method to be used. Defaults to "tmm".
#' @param n \code{integer}. Passed to \code{plotMDS()} (\code{top}): number of top genes used to calculate pairwise distances. Defaults to 500.
#' @param gene.selection \code{character}. passed to \code{plotMDS()} specifying the mode to select genes for comparisons. Defaults to "common".
#' @param circle \code{logical}. Draw a correlation circle around points representing correlating samples? Only applies when prcomp was called with scale = TRUE and when var.scale = 1. Defaults to \code{TRUE}.
#' @param varname.size \code{numeric}. Size of the text for variable names. Defaults to 0.
#' @param var.axes \code{logical}. Draw arrows for the variables? Defaults to \code{FALSE}.
#' @param PC \code{numeric}. Which principal components to plot. Defaults to 1:3.
#' @param type \code{character}. Which type of plot to produce. Needs to be one of "both", "uncorrected", "pseudo-corrected"
#'   describing which values should be plotted. "uncorrected" will plot the input counts matrix while
#'   "pseudo-corrected" will plot pseudo counts for blocked designs (e.g., paired samples or batch factors).
#'   Defaults to "both".
#' @param plots \code{logical}. Should plots be produced? Defaults to \code{TRUE}.
#' @param lists \code{logical}. Should lists of differentially expressed genes be produced? Defaults to \code{TRUE}.
#' @param filtered.lists \code{logical}. Should lists of differentially expressed genes be produced for filtered data? Defaults to \code{TRUE}.
#'
#' @details For experimental designs involving comparisons within as well as between subjects inter-subject needs to be computed.
#'     In this case, the column specified in the 'pairs' argument must assign the subjects to the treatment/tissue/etc groups.
#'     For example, if we have two treatments the effects of which are to be observed in each two tissues, this design would apply.
#'     The 'pairs' factor is passed to the functions 'duplicateCorrelation()' and 'lmFit()'.
#'     The 'block' argument is used to specify whether the comparisons are to be made within AND between subjects or in the case of
#'     technical replicates, i.e., if the samples are not independent, in other words, correlated. That correlation is addressed by
#'     means of the 'duplicateCorrelation()' function in the limma package.
#'     If 'block' is set to \code{TRUE}, the (selected) 'voom' function is enforced.
#'     The (selected) 'voom' function is also enforced if paired samples are detected.
#'     As of version 0.4, the 'edgeR::voomLmFit()' function is incorporated, which replaces 'voom()', 'lmFit()' and
#'     'voomWithQualityWeights()'. voomLmFit()' ensures unbiased estimation of the residual variances and automates the estimation
#'     of sample weights and intrablock correlations.
#'     In edgeR, it is recommended to remove features without at least 1 read per million in n of the
#'     samples, where n is the size of the smallest group of replicates (determined from the 'groups' vector).
#'     The 'min.samp' argument is used to specify the number of samples in which a feature needs to be covered by at least one read per million.
#'     Quasi-likelihood pipeline:
#'     While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL
#'     F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It
#'     provides more robust and reliable error rate control when the number of replicates is small.
#'     The QL dispersion estimation and hypothesis testing is done by using the functions
#'     glmQLFit() and glmQLFTest().
#' @return A list of all relevant objects generated in the course of the workflow:
#' \itemize{
#'   \item if `voom=TRUE`, voomed counts, otherwise DGEList object with TMM-normalisation factors
#'   \item if `voom=TRUE`, lmFit object, otherwise DGEList object with estimated dispersion
#'   \item if `voom=TRUE`, eBayes output of contrasted groups, otherwise glmFit object
#'   \item annotated topTable/topTags output
#' }
#' However, the function is first and foremost called for its side effects of generating results tables and plots.
#'
#' @export
diffExpr <-
  function(expr.dat,
           samp.info,
           control,
           design = NULL,
           samples = NULL,
           sample.plot.names = NULL,
           groups = NULL,
           pairs = NULL,
           block = FALSE,
           contrasts = NULL,
           out.dir = NULL,
           project.dir = ".",
           analysis.name = NULL,
           biomart = FALSE,
           biom.data.set = "hsapiens_gene_ensembl",
           biom.mart = c("ensembl", "snp", "funcgen", "vega", "pride", "plants"),
           host = "https://www.ensembl.org",
           biom.filter = "ensembl_gene_id",
           biom.attributes = c("ensembl_gene_id", "hgnc_symbol", "description", "entrezgene_id"),
           biom.force.ensg = FALSE,
           biom.cache = NULL,
           use.cache = FALSE,
           sym.col = "hgnc_symbol",
           rm.dups = FALSE,
           p.thr = 0.05,
           fdr.thr = 0.05,
           logfc.thr = 1,
           numlab = 20,
           point.lab = TRUE,
           heatmap.topn = 100,
           heatmap.split.expr = FALSE,
           min.samp = NULL,
           strict = TRUE,
           disp = c("gene", "trend", "common"),
           do.voom = FALSE,
           voom.fun = "voomLmFit",
           use_weights = FALSE,
           norm.method = c("tmm", "quantile"),
           bayes.trend = FALSE,
           bayes.robust = FALSE,
           quasi.likelihood = NA,
           n = 500,
           gene.selection = "common",
           ellipse = TRUE,
           ellipse.mapping.groups = NULL,
           label.samples = TRUE,
           geom.point.size = 2,
           label.font.size = 5,
           plot.ellipse.legend = NA,
           circle = TRUE,
           varname.size = 0,
           var.axes = FALSE,
           PC = c(1, 2, 3),
           type = c("both", "uncorrected", "pseudo-corrected"),
           plots = TRUE,
           lists = TRUE,
           filtered.lists = TRUE,
           do.enrichment = TRUE,
           enrichment.plot = FALSE,
           enrichment.plot.fdr.thr = fdr.thr,
           enrichment.plot.logfc.thr = logfc.thr,
           enrichment.plot.num.terms = 5,
           enrichment.methods = c("clusterProfilerGO", "clusterProfilerKEGG", "gProfileR", "topGO"),
           ...,
           dry.run=FALSE)
  {
    ## initial checks
    cat("@ -- STARTUP CHECKS --\n\n")
    # test if needed packages are installed
    cat("  Checking if needed packages are installed...")
    if (do.enrichment && biom.data.set == "hsapiens_gene_ensembl" && !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop(
        paste("Package", sQuote("org.Hs.eg.db"), "must be installed to perform enrichment analysis.
              Please set 'do.enrichment' to FALSE to disable it but still run the rest of the pipeline."),
        call. = FALSE
      )
    }
    if (do.enrichment && biom.data.set == "mmusculus_gene_ensembl" && !requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop(
        paste("Package", sQuote("org.Mm.eg.db"), "must be installed to perform enrichment analysis.
              Please set 'do.enrichment' to FALSE to disable it but still run the rest of the pipeline."),
        call. = FALSE
      )
    }
    if (use.cache && !requireNamespace("rappdirs", quietly = TRUE)) {
      stop(
        paste("Package", sQuote("rappdirs"), "must be installed to use the biomart cache."),
        call. = FALSE
      )
    }
    cat("done\n")
    if (use.cache && is.null(biom.cache)) {
      biom.cache <- rappdirs::user_cache_dir("biomaRt")
    }
    cat("  Checking for necessary user input...")
    if (missing(expr.dat) || !is.character(unlist(expr.dat))) {
      stop("Need input file with raw expression values (read counts)!")
    }
    if (missing(samp.info)) {
      if (is.null(samples)) {
        stop("Need sample names!")
      }
      if (length(samples) == 1) {
        stop("Need more than 1 sample name!")
      }
      if (is.null(groups)) {
        stop("Need groups!")
      }
      if (length(groups) == 1) {
        stop("Need more than 1 groups!")
      }
      if (!is.null(ellipse.mapping.groups) && length(ellipse.mapping.groups) != length(groups)) {
        stop("Ellipse mapping groups factor needs to be of same length as groups factor")
      }
    } else {
      if (is.null(samples) || length(samples) > 1) {
        stop("Need name of colum containing sample names!")
      }
      if (is.null(groups) || length(groups) > 1) {
        stop("Need name of colum containing groupings!")
      }
      if (!is.data.frame(samp.info)) {
        stop("'samp.info' needs to be a data.frame object containing at least two columns: 'SampleNames' and 'Groups' (the actual names can be provided\nby the 'samples' and 'groups' arguments, repsectively.)")
      }
    }
    if (length(grep("Weights", voom.fun))) {
      cat ("  'voomWithQualityWeights' is set as voom function; using settings for sample weights...\n")
      use_weights <- TRUE
    }

    #Checking whether the enrichment analyses can be run
    kegg.enrichment.will.be.performed = sum(grepl("KEGG", enrichment.methods)) > 0
    #kegg pathway enrichments needs entrez-IDs while other approaches are fine with ensemble IDs
    if (do.enrichment & kegg.enrichment.will.be.performed) {
      if (!biomart) {
        stop("Cannot run KEGG-pathway enrichment without entrez IDs. Set biomart to TRUE or remove KEGG enrichment approach from enrichment.methods")
      }
    }
    if (kegg.enrichment.will.be.performed) {
      #entrez ids has to be retrieved (if not already included in biom.attributes)
      if (!length(grep("entrezgene_id", biom.attributes))) {
        biom.attributes <- c(biom.attributes, "entrezgene_id")
      }
    }
    cat("done\n\n")

    cat("@ -- PREPARING RUN ENVIRONMENT --\n\n")
    # setting standard output folder
    if (is.null(out.dir)) {
      cat("No output folder set. Using default\n")
      out.dir <- file.path(normalizePath(project.dir), "differential_expression")
    }
    # checking for content in existing output folder
    if (file.exists(out.dir)) {
      if (length(dir(out.dir))) {
        cat("  Output directory exists and is not empty! Found:\n")
        found <- dir(out.dir, include.dirs = TRUE)
        print(found)
        out.dir <- paste(normalizePath(out.dir), "new",
                         gsub("-", "", unlist(strsplit(as.character(Sys.time()), " "))[1]),
                         gsub(":", "", unlist(strsplit(as.character(Sys.time()), " "))[2]),
                         sep = "_")
        cat("  Creating new output folder", sQuote(out.dir), "...")
      } else {
        cat("  Output folder exists and is empty.\n")
      }
    }
    # creating output folder if not existing
    if (!file.exists(out.dir)) {
      cat("  Output directory does not exist. Creating...\n")
      if (!dry.run) {
        dir.create(out.dir)
      } else {
        cat(" ~~DRY RUN~~\n")
      }
    }
    cat("  Saving output to", out.dir, "\n")

    if (block) {
      do.voom <- TRUE
      cat("*** Using 'blocked' design (i.e., samples are measured _in blocks_). Enforcing 'voom'! ***\n")
    }

    ## switch off plotting device on exit in case a plot fails
    ## the resulting file will be empty but not broken
    on.exit({plyr::l_ply(dev.list(), dev.off)})

    ## enforce voom for paired samples
    #	if (!is.null(pairs)) {
    #		cat("Samples are paired! Enforcing 'voom'\n")
    #		do.voom <- TRUE
    #	}
    #
    ## standardize samp.info
    ### needs to be a data.frame with (at least) two columns: 'SampleNames' and 'Groups'
    ### if a data.frame is provided by the 'samp.info' argument the columns are renamed to meet that convention
    #### save Groups name and ellipse mapping groups name first
    ellipse.grp.nam <- ellipse.mapping.groups
    grp.nam <- groups
    #### reformat samp.info data frame to meet down-stream conditions
    samp.info <- diff_expr_get_samp_info(samp.info, samples, groups, ellipse.mapping.groups)

    ## make control group first levels to ensure it becomes the '(Intercept)'
    cat("  Setting control group...\n")
    if (is.numeric(samp.info$Groups)) {
      control <- paste0("group_", control)
    }
    groups <- relevel(samp.info$Groups, ref=control)

    #### extract vector of ellipse group names from samp.info
    ellipse.mapping.groups <- samp.info$Ellipse

    # setting quasi-likelihood method
    ## if NA will use the number of replicates to determine whether to use QL or LRT
    if (is.na(quasi.likelihood)) {
      min.samp <- min(table(samp.info$Groups))
      if (min.samp > 4) {
        quasi.likelihood <- FALSE
      } else {
        quasi.likelihood <- TRUE
      }
    }

    if (is.null(analysis.name)) {
      analysis.name <- paste0(paste(levels(groups)[1:2], collapse="_"), "_")
      cat("Using default settings for the 'analysis name':", sQuote(analysis.name), "\n")
      warning("A unique and descriptive name for the analysis should always be provided!")
    }

    cat("\n@ -- PREPROCESSING --\n\n")
    ## read counts
    counts <- diff_expr_read_counts(expr.dat, samp.info)

    ## Filter weakly expressed and non-informative (e.g., non-aligned) features
    cat(" Filtering counts...\n")
    counts <- diff_expr_filter_counts(counts, samp.info, strict, min.samp)

    ## Create a DGEList object (edgeR’s container for RNA-seq count data)
    d <- edgeR::DGEList(counts = counts, group = groups)
    ## Estimate normalization factors
    d <- edgeR::calcNormFactors(d)
    ## Inspect the relationships between samples using a multidimensional scaling (MDS) plot
    if (plots) {
      cat(" MDS plot...\n")
      if (!is.null(sample.plot.names)) {
        spn <- samp.info[, c("SampleNames", sample.plot.names)]
        spn <- spn[match(colnames(d), spn$SampleNames), ]
      } else {
        sample.plot.names <- "SamplePlotNames"
        spn <- samp.info["SampleNames"]
        spn$SamplePlotNames <- spn$SampleNames
        spn <- spn[match(colnames(d), spn$SampleNames), ]
      }
      if (identical(as.character(spn$SampleNames), colnames(d))) {
        sample.plot.names <- as.character(spn[, sample.plot.names])
        if (length(sample.plot.names)>4) {
          spr <- paste0(paste(sQuote(sample.plot.names)[1:4], collapse=", "), "(, truncated...)\n")
        } else {
          spr <- paste0(paste(sQuote(sample.plot.names), collapse=", "), "(, truncated...)\n")
        }
        cat("  Using pretty labels for MDS plot:", spr)
      } else {
        cat("  NOTE: Sanity check for pretty names failed. Falling back to column names...\n")
        sample.plot.names <- as.character(colnames(d)) ## check!
      }
      names(sample.plot.names) <- colnames(d)

      diff_expr_mds_plot(d, groups=groups, n=n, sample.plot.names=sample.plot.names, analysis.name=analysis.name, do.pdf=TRUE, out.dir=out.dir)
    }

    # create design matrix
    if (is.null(design)) {
      user_design <- FALSE
      # create a design matix based on the choices for pairs or block
      ## if block is TRUE or there are no paired samples, a simple design matrix is created without an intercept
      ### the block argument here is logical and determines whether or not the samples are independent
      ## if there are paired samples, a design matrix with covariates is created which will have an intercept;
      ### this is used for any (block) design involving batch effects or similar
      ### NOTE that the term 'block' in this context refers to the batch or pair factor and must be provided in the pairs column/argument
      ### internally it is used with duplicateCorrelation()
      design <- diff_expr_make_design(samp.info = samp.info,
                                      groups = groups,
                                      pairs = pairs,
                                      block = block,
                                      use_weights = use_weights)
    } else {
      grp.col <- grep(paste0("^", grp.nam), colnames(design))
      if (length(grp.col)) {
        colnames(design)[grp.col] <- sub(paste0("^", grp.nam), "groups", colnames(design)[grp.col])
      }
      user_design <- TRUE
    }
    if (length(grep(":", colnames(design)))) {
      cat("  Checking design column names...\n")
      colnames(design) <- gsub(":", "__", colnames(design))
    }

    ## set contrasts
    contrasts <- diff_expr_make_contrasts(design = design,
                                          groups = groups,
                                          pairs = pairs,
                                          block = block,
                                          contrasts = contrasts)
                                          # user_design = user_design)

    ## set dispersion method
    disp <- match.arg(disp)
    disp <- switch(disp, gene="tagwise.dispersion", trend="trended.dispersion", common="bin.dispersion")

    ## extract 'pairs' factor from 'samp.info' and keep column name in 'pairs_col' variable
    if (!is.null(pairs)) {
      pairs_col <- pairs
      pairs <- samp.info[[pairs]]
    }

    cat("@ -- LINEAR MODELLING --\n\n")
    ## voom
    if (do.voom) {
      cat("Running 'voom'...\n")
      norm.method <- match.arg(norm.method)
      cat("  Running", sQuote(voom.fun), "with", norm.method, "normalisation...\n")
      ## compute linear model fit and optionally apply voom beforehand
      #### NOTE: voom generates log2-cpms

      fit.l <- diff_expr_fit(counts = counts,
                             d = d,
                             design = design,
                             do.voom=TRUE,
                             voom.fun = voom.fun,
                             norm.method = norm.method,
                             quasi.likelihood = quasi.likelihood,
                             bayes.trend = bayes.trend,
                             bayes.robust = bayes.robust,
                             pairs = pairs,
                             pairs_col = pairs_col,
                             block = block,
                             contrasts = contrasts,
                             use_weights = use_weights)

      out.l <- list(v=fit.l$v, fit=fit.l$fit, fit2=fit.l$fit2)
      if (!block && !is.null(pairs)) {
        cat("    Paired samples: additional eBayes on first fit object...\n")
        fit3 <- limma::eBayes(fit = fit.l$fit, trend = bayes.trend, robust = bayes.robust)
        out.l$fit3 <- fit.l$fit3
      }
      normcnt <- fit.l$v$E
    } else {
      fit.l <- diff_expr_fit(counts = counts,
                             d = d,
                             design = design,
                             do.voom = FALSE,
                             quasi.likelihood = quasi.likelihood,
                             bayes.trend = bayes.trend,
                             bayes.robust = bayes.robust,
                             disp = disp)
      out.l <- list(d=fit.l$d, d2=fit.l$d2, fit=fit.l$fit)
      # normalized expression
      # Get the depth-adjusted reads per million
      cat("Generating output table of differentially expressed features...\n")
      normcnt <- edgeR::cpm(fit.l$d2, normalized.lib.sizes=TRUE, log=TRUE, prior.count=3)
    }

    type <- match.arg(type)
    if (type=="both") {
      cat("Plotting both uncorrected and pseudo-corrected if applicable...\n")
    }
    if (plots) {
      if (type %in% c("both", "uncorrected")) {
        ## quality control plots for uncorrected data
        cat("Quality control plots...\n")
        cat("  Uncorrected, normalised data...\n")
        if (do.voom) {
          type.plot <- "uncorrected, normalised voom output"
        } else {
          type.plot <- "uncorrected, normalised DGEList"
        }
        out.l <- diff_expr_QC_plots(counts=normcnt,
                                    samp.info=samp.info,
                                    control=control,
                                    out.l=out.l,
                                    grp.nam=grp.nam,
                                    PC=PC,
                                    sample.plot.names=sample.plot.names,
                                    ellipse=ellipse,
                                    ellipse.mapping.groups=ellipse.mapping.groups,
                                    ellipse.grp.nam=ellipse.grp.nam,
                                    label.samples=label.samples,
                                    geom.point.size=geom.point.size,
                                    label.font.size = label.font.size,
                                    plot.ellipse.legend=plot.ellipse.legend,
                                    circle=circle,
                                    varname.size=varname.size,
                                    var.axes=var.axes,
                                    pairs=pairs,
                                    pairs.name=pairs_col,
                                    gene.selection=gene.selection,
                                    n=n,
                                    type=type.plot,
                                    analysis.name=analysis.name,
                                    out.dir=out.dir)
      }

      if (!block && !is.null(pairs) && (type %in% c("both", "pseudo-corrected"))) {
        ## quality control plots for pseudo-corrected data
        cat("  Corrected, normalised pseudo-counts...\n")
        ## get pseudo counts for blocked designs (e.g., paired samples or batch factors)
        cat("   Calculating pseudo-counts...\n")
        pseudo.counts <- try(diff_expr_pseudo_counts(design=design, d=d, pairs="pairs", disp=disp, do.cpm=TRUE))
        if (is(pseudo.counts, "try-error")) {
          cat(" !NOTE! - Pseudo-count calculation failed. Omitting from output...\n")
        } else {
          cat("   done\n")
          out.l <- diff_expr_QC_plots(counts=pseudo.counts,
                                      samp.info=samp.info,
                                      control=control,
                                      out.l=out.l,
                                      grp.nam=grp.nam,
                                      PC=PC,
                                      sample.plot.names=sample.plot.names,
                                      ellipse=ellipse,
                                      ellipse.mapping.groups=ellipse.mapping.groups,
                                      ellipse.grp.nam=ellipse.grp.nam,
                                      label.samples=label.samples,
                                      geom.point.size=geom.point.size,
                                      label.font.size = label.font.size,
                                      plot.ellipse.legend=plot.ellipse.legend,
                                      circle=circle,
                                      varname.size=varname.size,
                                      var.axes=var.axes,
                                      pairs=pairs,
                                      pairs.name=pairs_col,
                                      gene.selection=gene.selection,
                                      n=n,
                                      type="pseudo-corrected, normalised DGEList",
                                      analysis.name=analysis.name,
                                      out.dir=out.dir)
        }
      }
    }


    cat("@ -- EXTRACTING CONTRASTS --\n\n")
    if (do.voom) {
      if (!block && !is.null(pairs)) {
        cat("  ...for paired samples comparisons (voom)...\n")
        cont <- grep("^groups.+", colnames(fit3$coefficients), value=TRUE)
        cat("  ", cont, "\n")
        out.l <- diff_expr_extract_contrasts(contrasts = cont,
                                             fit = fit.l$fit,
                                             fit2 = fit3,
                                             normcnt = normcnt,
                                             out.l = out.l,
                                             do.voom=TRUE,
                                             quasi.likelihood = quasi.likelihood,
                                             out.dir = out.dir,
                                             analysis.name = analysis.name,
                                             biomart = biomart,
                                             biom.data.set = biom.data.set,
                                             biom.mart = biom.mart,
                                             host = host,
                                             biom.filter = biom.filter,
                                             biom.attributes = biom.attributes,
                                             biom.force.ensg = biom.force.ensg,
                                             biom.cache = biom.cache,
                                             use.cache = use.cache,
                                             sym.col = sym.col,
                                             rm.dups = rm.dups,
                                             p.thr = p.thr,
                                             fdr.thr = fdr.thr,
                                             logfc.thr = logfc.thr,
                                             numlab = numlab,
                                             point.lab = point.lab,
                                             heatmap.topn = heatmap.topn,
                                             heatmap.split.expr = heatmap.split.expr,
                                             heatmap.main = paste0("limma::", voom.fun),
                                             font.size = label.font.size,
                                             plots = plots,
                                             lists = lists,
                                             filtered.lists = filtered.lists,
                                             samp.info = samp.info,
                                             samples = "SampleNames",
                                             groups = groups,
                                             sample.plot.names = sample.plot.names)
      }

      if (!is.null(contrasts)) {
        cat("  ...for all (remaining) comparisons (voom)...\n")
        out.l <- diff_expr_extract_contrasts(contrasts = contrasts,
                                             fit = fit.l$fit,
                                             fit2 = fit.l$fit2,
                                             normcnt = normcnt,
                                             out.l = out.l,
                                             do.voom=TRUE,
                                             quasi.likelihood = quasi.likelihood,
                                             out.dir = out.dir,
                                             analysis.name = analysis.name,
                                             biomart = biomart,
                                             biom.data.set = biom.data.set,
                                             biom.mart = biom.mart,
                                             host = host,
                                             biom.filter = biom.filter,
                                             biom.attributes = biom.attributes,
                                             biom.force.ensg = biom.force.ensg,
                                             biom.cache = biom.cache,
                                             use.cache = use.cache,
                                             sym.col = sym.col,
                                             rm.dups = rm.dups,
                                             p.thr = p.thr,
                                             fdr.thr = fdr.thr,
                                             logfc.thr = logfc.thr,
                                             numlab = numlab,
                                             point.lab = point.lab,
                                             heatmap.topn = heatmap.topn,
                                             heatmap.split.expr = heatmap.split.expr,
                                             heatmap.main = paste0("limma::", voom.fun),
                                             font.size = label.font.size,
                                             plots = plots,
                                             lists = lists,
                                             filtered.lists,
                                             samp.info = samp.info,
                                             samples = "SampleNames",
                                             groups = groups,
                                             sample.plot.names = sample.plot.names)
      }
    } else {
      if (!is.null(pairs)) {
      cat("  ...for paired samples comparisons (GLM)...\n")
        cont <- grep("^groups.+", colnames(fit.l$fit$coefficients), value=TRUE)
        cat("  ", cont, "\n")
        out.l <- diff_expr_extract_contrasts(contrasts = cont,
                                             fit = fit.l$fit,
                                             fit2 = NULL,
                                             normcnt = normcnt,
                                             out.l = out.l,
                                             do.voom=FALSE,
                                             quasi.likelihood = quasi.likelihood,
                                             out.dir = out.dir,
                                             analysis.name = analysis.name,
                                             biomart = biomart,
                                             biom.data.set = biom.data.set,
                                             biom.mart = biom.mart,
                                             host = host,
                                             biom.filter = biom.filter,
                                             biom.attributes = biom.attributes,
                                             biom.force.ensg = biom.force.ensg,
                                             biom.cache = biom.cache,
                                             use.cache = use.cache,
                                             sym.col = sym.col,
                                             rm.dups = rm.dups,
                                             p.thr = p.thr,
                                             fdr.thr = fdr.thr,
                                             logfc.thr = logfc.thr,
                                             numlab = numlab,
                                             point.lab = point.lab,
                                             heatmap.topn = heatmap.topn,
                                             heatmap.split.expr = heatmap.split.expr,
                                             heatmap.main = "edgeR GLM",
                                             font.size = label.font.size,
                                             plots = plots,
                                             lists = lists,
                                             filtered.lists = filtered.lists,
                                             samp.info = samp.info,
                                             samples = "SampleNames",
                                             groups = groups,
                                             sample.plot.names = sample.plot.names)
      }

      if (!is.null(contrasts)) {
        cat("  ...for all (remaining) comparisons (GLM)...\n")
        out.l <- diff_expr_extract_contrasts(contrasts = contrasts,
                                             fit = fit.l$fit,
                                             fit2 = NULL,
                                             normcnt = normcnt,
                                             out.l = out.l,
                                             do.voom=FALSE,
                                             quasi.likelihood = quasi.likelihood,
                                             out.dir = out.dir,
                                             analysis.name = analysis.name,
                                             biomart = biomart,
                                             biom.data.set = biom.data.set,
                                             biom.mart = biom.mart,
                                             host = host,
                                             biom.filter = biom.filter,
                                             biom.attributes = biom.attributes,
                                             biom.force.ensg = biom.force.ensg,
                                             biom.cache = biom.cache,
                                             use.cache = use.cache,
                                             sym.col = sym.col,
                                             rm.dups = rm.dups,
                                             p.thr = p.thr,
                                             fdr.thr = fdr.thr,
                                             logfc.thr = logfc.thr,
                                             numlab = numlab,
                                             point.lab = point.lab,
                                             heatmap.topn = heatmap.topn,
                                             heatmap.split.expr = heatmap.split.expr,
                                             heatmap.main = "edgeR GLM",
                                             font.size = label.font.size,
                                             plots = plots,
                                             lists = lists,
                                             filtered.lists = filtered.lists,
                                             samp.info = samp.info,
                                             samples = "SampleNames",
                                             groups = groups,
                                             sample.plot.names = sample.plot.names)
      }
    }

    if (do.enrichment) {
      cat("@ -- ENRICHMENT ANALYSIS --\n\n")
      enrich_spec <- switch(biom.data.set,
                            hsapiens_gene_ensembl = "human",
                            mmusculus_gene_ensembl = "mouse")
      out.l$enrichment <- runEnrichmentAnalyses(diffr.wrapper.output = out.l,
                                                analysis.name = analysis.name,
                                                use.background.from.diffr.output = TRUE,
                                                out.dir = out.dir,
                                                use.pval.in.DE.filtering.if.no.sign.fdrs = FALSE,
                                                species = enrich_spec,
                                                enrichment.methods = enrichment.methods,
                                                do.plot = enrichment.plot,
                                                plot.fdr.thr=enrichment.plot.fdr.thr,
                                                plot.logfc.thr=enrichment.plot.logfc.thr,
                                                plot.num.terms=enrichment.plot.num.terms
      )
    }
    on.exit()
    return(out.l)
  }
