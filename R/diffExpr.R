#' Main wrapper for executing the entire pipeline from reading in expression data such as count
#' files to producing text files and graphs
#' @description \command{diffExpr} is a convenience wrapper performing all steps automatically.
#'    Most sub-functions are exported and can be called by the user, as well, if desired.
#'    These functions may be applicable to different kinds of data/input, rely, however,
#'    on the conventions set for this package.
#' @section Argument groups:
#' The arguments are grouped by name prefix:
#' \itemize{
#'   \item no prefix: core experiment definition, plus the primary differential-expression thresholds
#'     (\code{p.thr}, \code{fdr.thr}, \code{logfc.thr}, \code{numlab}, \code{point.lab});
#'   \item \code{filter.*}: count pre-filtering;
#'   \item \code{fit.*}: model fitting (voom / edgeR options);
#'   \item \code{biom.*}: 'biomart' gene annotation;
#'   \item \code{qc.*}: quality-control plots (MDS / PCA);
#'   \item \code{hm.*}: heatmaps, including their own significance thresholds and colour palette;
#'   \item \code{enr.*}: enrichment analysis;
#'   \item \code{out.*}, \code{verbose}, \code{log.file}, \code{dry.run}: output and run control.
#' }
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
#' @param block \code{logical}. Are the samples not independent? See Details section.
#' @param contrasts \code{character}. Vector of contrasts to be made. If not provided, all possible contrasts will be made.
#'   This specifies group name pairs to be compared in the format expected by \code{makeContrasts()}, i.e., "group2-group1".
#' @param out.dir \code{character}. Path to the output directory. This argument is required: all
#'   result tables, plots and the run log are written below it. There is deliberately no default,
#'   so that no files are ever created in the working directory unintentionally. Use, e.g.,
#'   \code{out.dir = tempdir()} to try the pipeline out without keeping the results.
#' @param analysis.name \code{character}. Name of the analysis. If not provided, a default name will be generated.
#' @param filter.strict \code{logical}. For miRNA analysis: only keep a miRNA if there are > 5 reads per million in at least half of the samples?
#' @param filter.min.samp \code{integer}. Number of samples in which a feature needs to be covered by at least one read per million.
#'   Defaults to the size of the smallest group of replicates. See \emph{details}.
#' @param fit.voom \code{logical}. Should the voom function be used? Defaults to \code{FALSE}.
#' @param fit.voom.fun \code{character}. The voom function to be used. Should be one of 'limma::voom', 'limma::voomWithQualityWeights' or
#' 'edgeR::voomLmFit'. Defaults to \code{edgeR::voomLmFit}. See 'Details'.
#' @param fit.use.weights \code{logical}. Should sample-specific quality weights be estimated?
#' @param fit.norm.method \code{character}. The normalisation method to be used. Defaults to "tmm".
#' @param fit.disp \code{character}. The dispersion method to be used. Defaults to "gene".
#' @param fit.bayes.trend \code{logical}. Should an intensity-trend be allowed for the prior variance? Passed to 'limma::eBayes'.
#' @param fit.bayes.robust \code{logical}. Should the estimation of df.prior and var.prior be robustified against outlier sample variances? Passed to 'limma::eBayes'.
#' @param fit.quasi.likelihood Logical; should quasi-likelihood methods be used? See \emph{Details} section.
#'   Defaults to NA, which will determine the method based on the number of replicate samples.
#'   If more than 4 replicates are present, the likelihood ratio test is used, otherwise the quasi-likelihood methods.
#'   If \code{TRUE}, then the quasi-likelihood methods are used, if \code{FALSE}, then the likelihood ratio test is used.
#' @param p.thr \code{numeric}. Threshold for p-values. Defaults to 0.05.
#' @param fdr.thr \code{numeric}. Threshold for FDR values. Defaults to 0.05.
#' @param logfc.thr \code{numeric}. Threshold for fold-change values on the log2-scale. Defaults to 1.
#' @param numlab \code{numeric}. Maximum number of point labels to be shown in the plot. This overrides/limits
#'   values calculated by any thresholds. Defaults to 25.
#' @param point.lab \code{logical}. Should point labels be shown in the plot? Defaults to \code{TRUE}.
#' @param de.plot.base.size \code{numeric}. Overall text/point scale (the `ggplot2` \code{base_size}) for the
#'   M-A and volcano plots, defaulting to 16. Raising it enlarges legend, axes, titles, point labels and
#'   points together. The default (16, up from the `ggplot2` default of ~11) keeps text readable on the
#'   default per-contrast plot PDF, which is opened at 15x15 inches to hold the heatmaps.
#' @param biom.use \code{logical}. Should the biomart be used for gene annotation? Defaults to \code{FALSE}.
#' @param biom.data.set \code{character}. The biomart dataset to be used. Defaults to "hsapiens_gene_ensembl".
#' @param biom.mart \code{character}. The biomart to be used.
#' @param biom.host \code{character}. The host to be used for the biomart. Defaults to "www.ensembl.org".
#' @param biom.filter \code{character}. The biomart filter to be used. Defaults to "ensembl_gene_id".
#' @param biom.attributes \code{character}. The biomart attributes to be used. Defaults to c("ensembl_gene_id", "hgnc_symbol", "description", "entrezgene_id").
#' @param biom.force.ensg \code{logical}. Should Ensembl Gene IDs be checked for (and stripped of) Ensembl version numbers? Defaults to \code{FALSE}.
#' @param biom.cache \code{character}. Path name giving the location of the cache \command{getBM()} uses if \code{use.cache=TRUE}. Defaults to the value in the \emph{BIOMART_CACHE} environment variable.
#' @param biom.use.cache (\code{logical}). Should \command{getBM()} use the cache? Defaults to \code{TRUE} as in the \command{getBM()} function and is passed on to that.
#' @param biom.sym.col \code{character}. Name of the column in the query result with gene symbols
#' @param biom.rm.dups \code{logical}. Should duplicates be removed from the output of the biomart request? Defaults to \code{FALSE}.
#' @param qc.top.n \code{integer}. Passed to \code{plotMDS()} (\code{top}): number of top genes used to calculate pairwise distances. Defaults to 500.
#' @param qc.gene.selection \code{character}. passed to \code{plotMDS()} specifying the mode to select genes for comparisons. Defaults to "common".
#' @param qc.pc \code{numeric}. Which principal components to plot. Defaults to 1:3.
#' @param qc.type \code{character}. Which type of plot to produce. Needs to be one of "both", "uncorrected", "pseudo-corrected"
#'   describing which values should be plotted. "uncorrected" will plot the input counts matrix while
#'   "pseudo-corrected" will plot pseudo counts for blocked designs (e.g., paired samples or batch factors).
#'   Defaults to "both".
#' @param qc.ellipse \code{logical}. Should an ellipse be plotted around samples belonging to the same sample group? Defaults to \code{TRUE}.
#' @param qc.ellipse.groups \code{character} The name of the column in 'samp.info' with group names for ellipse drawing. If \code{NULL} (default)
#'     will use the \code{groups} column. If 'samp.info' is not supplied vector of groups.
#' @param qc.ellipse.legend \code{logical}. Should a legend be added for ellipses in PCA plots? NA, the default, includes
#'     if any aesthetics are mapped. FALSE never includes, and TRUE always includes. It can also be a named logical vector to finely select
#'     the aesthetics to display.
#' @param qc.label.samples \code{logical}. Should points in appropriate QC plots be labelled. So far, applies only to PCA ggplot. Defaults to \code{TRUE}.
#' @param qc.point.size \code{numeric}. Size of points in appropriate QC plots. So far, applies only to PCA ggplot. Defaults to 2.
#' @param qc.label.size \code{numeric}. Font size used for point labels in appropriate QC plots. So far, applies to PCA ggplot and M-A plots. Defaults to 5.
#' @param qc.circle \code{logical}. Draw a correlation circle around points representing correlating samples? Only applies when prcomp was called with scale = TRUE and when var.scale = 1. Defaults to \code{TRUE}.
#' @param qc.varname.size \code{numeric}. Size of the text for variable names. Defaults to 0.
#' @param qc.var.axes \code{logical}. Draw arrows for the variables? Defaults to \code{FALSE}.
#' @param hm.topn \code{numeric}. Number of top values to be plotted. Defaults to 100.
#' @param hm.p.thr,hm.fdr.thr,hm.logfc.thr \code{numeric}. Significance and fold-change thresholds used
#'   specifically for selecting genes shown in the heatmaps, kept separate from \option{p.thr},
#'   \option{fdr.thr} and \option{logfc.thr} (which control the plots and tables) because heatmaps
#'   usually read best with a stricter gene set. Default to 0.05, 0.05 and 1.
#' @param hm.split.expr \code{logical}. Should the top up- and top down-regulated genes be displayed at equal numbers (50/50),
#' if they meet the significance threshold (regardless of the actual significance)? Defaults to \code{FALSE}.
#' @param hm.pal.blind string determining the RColorBrewer color blind palette (default = "PuOr");
#' other option can be visualized with the following command: brewer.pal.info[brewer.pal.info$colorblind,]
#' @param hm.pal.n desired length of the number of different colours in 'color.blind.pal'. Will also be used
#' as length of the numeric vector of probabilities in 'quantile_breaks()' (see ?quantile); defaults to 11
#' @param hm.pal.extremes character vector of length 2 giving the two extremes of a user-defined colour palette
#' varying from the first hue to the second via white.
#' @param hm.pal.length integer setting the desired length of the colour palette to be used in the heatmap
#' @param hm.anno.color list of named character vectors giving the colours used in the heatmap annotation bars. See 'annotation_colors'
#' in [pheatmap()]. Automatically generated if NULL (default).
#' @param hm.anno.name character string used as the column annotation legend title. If 'anno.color' is not NULL and of length 1 the slot name will be used if existing.
#' @param enr.do \code{logical}. Whether or not to call enrichment wrapper. Defaults to \code{TRUE}.
#' @param enr.methods \code{character}. One or more of the following: c("clusterProfilerGO", "clusterProfilerKEGG", "gProfileR", "topGO"). By default, uses them all.
#' @param enr.plot \code{logical}. Whether or not to draw a network plot for the enrichment results. Defaults to \code{FALSE}.
#' @param enr.plot.fdr.thr \code{numeric}. FDR threshold used in the enrichment plot. This may be useful to tweak to produce a more informative plot.
#' Defaults to 0.05.
#' @param enr.plot.logfc.thr \code{numeric}. FC threshold on the log2-scale used in the enrichment plot. Defaults to 1.
#' @param enr.plot.num.terms \code{integer}. Number of terms shown in the plot. Defaults to 5.
#' @param out.plots \code{logical}. Should plots be produced? Defaults to \code{TRUE}.
#' @param out.tables \code{logical}. Should lists of differentially expressed genes be produced? Defaults to \code{TRUE}.
#' @param out.filtered.tables \code{logical}. Should lists of differentially expressed genes be produced for filtered data? Defaults to \code{TRUE}.
#'
#' @param verbose \code{logical} or \code{character}. Controls console output. \code{TRUE} (default) prints
#'   major workflow steps, \code{FALSE} prints nothing, and \code{"all"} mirrors the full detail of the log
#'   file to the console. Console output is written to \code{stdout} and is controlled solely by this argument.
#'   Irrespective of this setting, the complete log is always written to \option{log.file}.
#' @param log.file \code{character}. Path to the run log file. If \code{NULL} (default), a log file named
#'   \emph{<analysis.name>_diffwrap.log} is created in \option{out.dir}. No log file is written if
#'   \code{dry.run=TRUE}.
#' @param dry.run \code{logical}. If \code{TRUE}, the function will not create any output files or directories.
#' @param ... \code{ANY}. Additional arguments passed to functions.
#' @details For experimental designs involving comparisons within as well as between subjects inter-subject needs to be computed.
#'     In this case, the column specified in the 'pairs' argument must assign the subjects to the treatment/tissue/etc groups.
#'     For example, if we have two treatments the effects of which are to be observed in each two tissues, this design would apply.
#'     The 'pairs' factor is passed to the functions 'duplicateCorrelation()' and 'lmFit()'.
#'     The 'block' argument is used to specify whether the comparisons are to be made within AND between subjects or in the case of
#'     technical replicates, i.e., if the samples are not independent, in other words, correlated. That correlation is addressed by
#'     means of the 'duplicateCorrelation()' function in the limma package.
#'     If 'block' is set to \code{TRUE}, the (selected) 'voom' function is enforced.
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
#' @examples
#' \donttest{
#' out.dir <- file.path(tempdir(), "diffwrap_demo")
#' dir.create(out.dir, showWarnings = FALSE)
#' res <- diffExpr(expr.dat = diffwrap_counts,
#'                 samp.info = diffwrap_samp_info,
#'                 samples = "SampleName", groups = "Group",
#'                 control = "control", analysis.name = "demo",
#'                 out.dir = out.dir, enr.do = FALSE)
#' names(res$contrasts)
#' }
#' @export
diffExpr <-
  function(## --- core experiment definition ---
           expr.dat,
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
           analysis.name = NULL,
           ## --- count pre-filtering (filter.*) ---
           filter.strict = TRUE,
           filter.min.samp = NULL,
           ## --- model fitting (fit.*) ---
           fit.voom = FALSE,
           fit.voom.fun = "voomLmFit",
           fit.use.weights = FALSE,
           fit.norm.method = c("tmm", "quantile"),
           fit.disp = c("gene", "trend", "common"),
           fit.bayes.trend = FALSE,
           fit.bayes.robust = FALSE,
           fit.quasi.likelihood = NA,
           ## --- differential-expression thresholds, tables and MA/volcano labels ---
           p.thr = 0.05,
           fdr.thr = 0.05,
           logfc.thr = 1,
           numlab = 25,
           point.lab = TRUE,
           de.plot.base.size = 16,
           ## --- biomart annotation (biom.*) ---
           biom.use = FALSE,
           biom.data.set = "hsapiens_gene_ensembl",
           biom.mart = "ensembl",
           biom.host = "https://www.ensembl.org",
           biom.filter = "ensembl_gene_id",
           biom.attributes = c("ensembl_gene_id", "hgnc_symbol", "description", "entrezgene_id"),
           biom.force.ensg = FALSE,
           biom.cache = NULL,
           biom.use.cache = FALSE,
           biom.sym.col = "hgnc_symbol",
           biom.rm.dups = FALSE,
           ## --- QC plots: MDS / PCA (qc.*) ---
           qc.top.n = 500,
           qc.gene.selection = "common",
           qc.pc = c(1, 2, 3),
           qc.type = c("both", "uncorrected", "pseudo-corrected"),
           qc.ellipse = TRUE,
           qc.ellipse.groups = NULL,
           qc.ellipse.legend = NA,
           qc.label.samples = TRUE,
           qc.point.size = 2,
           qc.label.size = 5,
           qc.circle = TRUE,
           qc.varname.size = 0,
           qc.var.axes = FALSE,
           ## --- heatmaps (hm.*) ---
           hm.topn = 100,
           hm.p.thr = 0.05,
           hm.fdr.thr = 0.05,
           hm.logfc.thr = 1,
           hm.split.expr = FALSE,
           hm.pal.blind = "PuOr",
           hm.pal.n = 11,
           hm.pal.extremes = c("#3182BD", "#E6550D"),
           hm.pal.length = NULL,
           hm.anno.color = NULL,
           hm.anno.name = "Sample Class",
           ## --- enrichment (enr.*) ---
           enr.do = TRUE,
           enr.methods = c("clusterProfilerGO", "clusterProfilerKEGG", "gProfileR", "topGO"),
           enr.plot = FALSE,
           enr.plot.fdr.thr = fdr.thr,
           enr.plot.logfc.thr = logfc.thr,
           enr.plot.num.terms = 5,
           ## --- output and run control ---
           out.plots = TRUE,
           out.tables = TRUE,
           out.filtered.tables = TRUE,
           verbose = TRUE,
           log.file = NULL,
           dry.run = FALSE,
           ...)
  {
    ## ----------------------------------------------------------------------------------------------
    ## Map the categorised public argument names to the internal names used throughout the body.
    ## The public interface (above) is grouped by prefix (filter./fit./biom./qc./hm./enr./out.);
    ## the internal logic keeps its original names, so only this block changes when arguments are
    ## renamed. Editing anything below this block should use the internal (left-hand) names.
    ## ----------------------------------------------------------------------------------------------
    strict                    <- filter.strict
    min.samp                  <- filter.min.samp
    do.voom                   <- fit.voom
    voom.fun                  <- fit.voom.fun
    use_weights               <- fit.use.weights
    norm.method               <- fit.norm.method
    disp                      <- fit.disp
    bayes.trend               <- fit.bayes.trend
    bayes.robust              <- fit.bayes.robust
    quasi.likelihood          <- fit.quasi.likelihood
    biomart                   <- biom.use
    host                      <- biom.host
    use.cache                 <- biom.use.cache
    sym.col                   <- biom.sym.col
    rm.dups                   <- biom.rm.dups
    n                         <- qc.top.n
    gene.selection            <- qc.gene.selection
    PC                        <- qc.pc
    type                      <- qc.type
    ellipse                   <- qc.ellipse
    ellipse.mapping.groups    <- qc.ellipse.groups
    plot.ellipse.legend       <- qc.ellipse.legend
    label.samples             <- qc.label.samples
    geom.point.size           <- qc.point.size
    label.font.size           <- qc.label.size
    circle                    <- qc.circle
    varname.size              <- qc.varname.size
    var.axes                  <- qc.var.axes
    heatmap.topn              <- hm.topn
    heatmap.split.expr        <- hm.split.expr
    color.blind.pal           <- hm.pal.blind
    n.pal.cols                <- hm.pal.n
    color.extremes            <- hm.pal.extremes
    palette.length            <- hm.pal.length
    anno.color                <- hm.anno.color
    anno.name                 <- hm.anno.name
    do.enrichment             <- enr.do
    enrichment.methods        <- enr.methods
    enrichment.plot           <- enr.plot
    enrichment.plot.fdr.thr   <- enr.plot.fdr.thr
    enrichment.plot.logfc.thr <- enr.plot.logfc.thr
    enrichment.plot.num.terms <- enr.plot.num.terms
    plots                     <- out.plots
    lists                     <- out.tables
    filtered.lists            <- out.filtered.tables

    ## start logging
    ## the log file itself is opened further below, once 'out.dir' is known; until then
    ## log lines are buffered in memory so that nothing from the startup checks is lost
    dw_log_start(verbose)
    on.exit(dw_log_end(), add = TRUE)
    run_start <- Sys.time()   # start timing the run phase

    ## initial checks
    dw_step("@ -- STARTUP CHECKS --\n\n")
    # test if needed packages are installed
    dw_log("  Checking if needed packages are installed...", "\n")
    if (do.enrichment && biom.data.set == "hsapiens_gene_ensembl" && !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop(
        paste("Package", sQuote("org.Hs.eg.db"), "must be installed to perform enrichment analysis.
              Please set 'enr.do' to FALSE to disable it but still run the rest of the pipeline."),
        call. = FALSE
      )
    }
    if (do.enrichment && biom.data.set == "mmusculus_gene_ensembl" && !requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop(
        paste("Package", sQuote("org.Mm.eg.db"), "must be installed to perform enrichment analysis.
              Please set 'enr.do' to FALSE to disable it but still run the rest of the pipeline."),
        call. = FALSE
      )
    }
    # enrichment currently only supports human and mouse; fail here rather than after the
    # whole analysis, when 'species' would resolve to NULL inside runEnrichmentAnalyses()
    if (do.enrichment && !biom.data.set %in% c("hsapiens_gene_ensembl", "mmusculus_gene_ensembl")) {
      stop("Enrichment analysis currently supports only 'hsapiens_gene_ensembl' or ",
           "'mmusculus_gene_ensembl' (via 'biom.data.set'). Set 'enr.do = FALSE' ",
           "to run the rest of the pipeline for other datasets.", call. = FALSE)
    }
    # The enrichment engines live in 'Suggests'. If enrichment is requested, check up front
    # that the packages needed by the selected methods are installed, so the run fails here
    # rather than after all the modelling and plotting is done.
    if (do.enrichment) {
      enr.needed <- c(
        if (any(grepl("clusterProfiler", enrichment.methods))) "clusterProfiler",
        if (any(grepl("gProfileR", enrichment.methods)))       "gprofiler2",
        if (any(grepl("topGO", enrichment.methods)))           "topGO",
        "WriteXLS"
      )
      enr.miss <- enr.needed[!vapply(enr.needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(enr.miss)) {
        stop("Package(s) ", paste(sQuote(enr.miss), collapse = ", "),
             " must be installed to perform enrichment analysis.\n",
             "  Please install them, or set 'enr.do = FALSE' to run the rest of the pipeline.",
             call. = FALSE)
      }
    }
    if (use.cache && !requireNamespace("rappdirs", quietly = TRUE)) {
      stop(
        paste("Package", sQuote("rappdirs"), "must be installed to use the biomart cache."),
        call. = FALSE
      )
    }
    if (use.cache && is.null(biom.cache)) {
      biom.cache <- rappdirs::user_cache_dir("biomaRt")
    }
    dw_log("  Checking for necessary user input...", "\n")
    ## 'expr.dat' may be file path(s) (character) OR an in-memory count matrix/data.frame;
    ## diff_expr_read_counts() handles all of these, so the guard must accept them too.
    if (missing(expr.dat) ||
        !(is.character(unlist(expr.dat)) || is.matrix(expr.dat) || is.data.frame(expr.dat))) {
      stop("Need input for 'expr.dat': file path(s) to read counts, or a count matrix/data.frame!")
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
      dw_log("  'voomWithQualityWeights' is set as voom function; using settings for sample weights...\n")
      use_weights <- TRUE
    }

    #Checking whether the enrichment analyses can be run
    ## Probe BioMart up front so the decision holds for the whole run: if it is unreachable, fall back
    ## to offline convertId2() gene symbols (biomart off) and drop KEGG, which needs Entrez IDs.
    if (biomart) {
      dw_log("  Checking BioMart availability...\n")
      if (!dw_biomart_available(biom.data.set)) {
        warning("BioMart is not reachable (primary host or mirrors). Gene annotation will use ",
                "offline convertId2() symbols instead; KEGG enrichment (which needs Entrez IDs) ",
                "will be skipped.", call. = FALSE)
        dw_step("  !NOTE! BioMart unavailable -> offline symbols via convertId2(); KEGG skipped.\n")
        biomart <- FALSE
      }
    }
    #kegg pathway enrichment needs entrez-IDs (from biomart) while other approaches are fine with Ensembl IDs
    kegg.enrichment.will.be.performed <- sum(grepl("KEGG", enrichment.methods)) > 0
    if (do.enrichment && kegg.enrichment.will.be.performed && !biomart) {
      dw_log("  Removing KEGG from enrichment methods (needs BioMart/Entrez IDs).\n")
      enrichment.methods <- enrichment.methods[!grepl("KEGG", enrichment.methods)]
      kegg.enrichment.will.be.performed <- FALSE
    }
    if (kegg.enrichment.will.be.performed) {
      #entrez ids has to be retrieved (if not already included in biom.attributes)
      if (!length(grep("entrezgene_id", biom.attributes))) {
        biom.attributes <- c(biom.attributes, "entrezgene_id")
      }
    }
    dw_log("\n")

    dw_step("@ -- PREPARING RUN ENVIRONMENT --\n\n")
    # setting output folder
    ## NOTE: there is deliberately no default output location. Writing results into the
    ## working directory (or anywhere else in the user's home filespace) without the user
    ## explicitly asking for it is not permitted by CRAN policy, so 'out.dir' is required.
    if (is.null(out.dir)) {
      stop("'out.dir' is required: please supply the directory the results should be written to.\n",
           "  No default is used so that no files are created in the working directory unintentionally.\n",
           "  To try the pipeline without keeping the output, use e.g. out.dir = tempdir().",
           call. = FALSE)
    }
    if (!is.character(out.dir) || length(out.dir) != 1L) {
      stop("'out.dir' must be a single path.", call. = FALSE)
    }
    # checking for content in existing output folder
    if (file.exists(out.dir)) {
      if (length(dir(out.dir))) {
        dw_log("  Output directory exists and is not empty! Found:\n")
        found <- dir(out.dir, include.dirs = TRUE)
        dw_log_obj(found)
        out.dir <- paste(normalizePath(out.dir), "new",
                         gsub("-", "", unlist(strsplit(as.character(Sys.time()), " "))[1]),
                         gsub(":", "", unlist(strsplit(as.character(Sys.time()), " "))[2]),
                         sep = "_")
        dw_log("  Creating new output folder", sQuote(out.dir), "...", "\n")
      } else {
        dw_log("  Output folder exists and is empty.\n")
      }
    }
    # creating output folder if not existing
    if (!file.exists(out.dir)) {
      dw_log("  Output directory does not exist. Creating...\n")
      if (!dry.run) {
        dir.create(out.dir)
      } else {
        dw_log(" ~~DRY RUN~~\n")
      }
    }
    dw_log("  Saving output to", out.dir, "\n")

    ## Resolve the analysis mode once, up front. This is the ONLY place where the
    ## (pairs, block, do.voom) flags are interpreted and validated; everything
    ## downstream follows the resolved spec. 'pairs' here is still the column name
    ## (or NULL), i.e. before it is turned into a factor, as the resolver expects.
    mode <- dw_resolve_mode(pairs = pairs, block = block, do.voom = do.voom)
    do.voom <- mode$do.voom                                # honour block -> voom enforcement
    dw_step(paste(format(mode), collapse = "\n"), "\n")    # log role/engine/design/contrasts/notes

    ## switch off plotting device on exit in case a plot fails
    ## the resulting file will be empty but not broken
    ## NOTE: 'add=TRUE' is essential here, otherwise this would replace the logging exit handler
    on.exit({plyr::l_ply(dev.list(), dev.off)}, add = TRUE)

    ## standardize samp.info
    ### needs to be a data.frame with (at least) two columns: 'SampleNames' and 'Groups'
    ### if a data.frame is provided by the 'samp.info' argument the columns are renamed to meet that convention
    #### save Groups name and ellipse mapping groups name first
    ellipse.grp.nam <- ellipse.mapping.groups
    grp.nam <- groups
    #### reformat samp.info data frame to meet down-stream conditions
    samp.info <- diff_expr_get_samp_info(samp.info, samples, groups, ellipse.mapping.groups)

    ## make control group first levels to ensure it becomes the '(Intercept)'
    dw_log("  Setting control group...\n")
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
      dw_log("Using default settings for the 'analysis name':", sQuote(analysis.name), "\n")
      warning("A unique and descriptive name for the analysis should always be provided!")
    }

    ## open the log file now that both 'out.dir' and 'analysis.name' are known
    ## everything logged so far is buffered and gets flushed into the file here
    if (!dry.run) {
      if (is.null(log.file)) {
        log.file <- file.path(out.dir, paste(analysis.name, "diffwrap.log", sep = "_"))
      }
      dw_log_file(log.file)
      dw_step("  Writing run log to", sQuote(log.file), "\n")
    }

    dw_step("\n@ -- PREPROCESSING --\n\n")
    ## read counts
    counts <- diff_expr_read_counts(expr.dat, samp.info)

    ## Filter weakly expressed and non-informative (e.g., non-aligned) features
    dw_log(" Filtering counts...\n")
    counts <- diff_expr_filter_counts(counts, samp.info, strict, min.samp)

    ## Create a DGEList object (edgeR’s container for RNA-seq count data)
    d <- edgeR::DGEList(counts = counts, group = groups)
    ## Estimate normalization factors
    d <- edgeR::calcNormFactors(d)
    ## Inspect the relationships between samples using a multidimensional scaling (MDS) plot
    if (plots) {
      dw_log(" MDS plot...\n")
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
        dw_log("  Using pretty labels for MDS plot:", spr, "\n")
      } else {
        dw_log("  NOTE: Sanity check for pretty names failed. Falling back to column names...\n")
        sample.plot.names <- as.character(colnames(d)) ## check!
      }
      names(sample.plot.names) <- colnames(d)

      diff_expr_mds_plot(d, groups=groups, n=n, sample.plot.names=sample.plot.names, analysis.name=analysis.name, do.pdf=TRUE, out.dir=out.dir)
    }

    # create design matrix
    if (is.null(design)) {
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
    }
    if (length(grep(":", colnames(design)))) {
      dw_log("  Checking design column names...\n")
      colnames(design) <- gsub(":", "__", colnames(design))
    }

    ## set contrasts
    contrasts <- diff_expr_make_contrasts(design = design,
                                          groups = groups,
                                          pairs = pairs,
                                          block = block,
                                          contrasts = contrasts)

    ## set dispersion method
    disp <- match.arg(disp, c("gene", "trend", "common"))
    disp <- switch(disp, gene="tagwise.dispersion", trend="trended.dispersion", common="bin.dispersion")

    ## extract 'pairs' factor from 'samp.info' and keep column name in 'pairs_col' variable
    if (!is.null(pairs)) {
      pairs_col <- pairs
      pairs <- samp.info[[pairs]]
    }

    dw_step("@ -- LINEAR MODELLING --\n\n")
    ## voom
    if (do.voom) {
      dw_log("Running 'voom'...\n")
      norm.method <- match.arg(norm.method, c("tmm", "quantile"))
      dw_log("  Running", sQuote(voom.fun), "with", norm.method, "normalisation...\n")
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
        dw_log("    Paired samples: additional eBayes on first fit object...\n")
        fit3 <- limma::eBayes(fit = fit.l$fit, trend = bayes.trend, robust = bayes.robust)
        out.l$fit3 <- fit3
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
      dw_log("Generating output table of differentially expressed features...\n")
      normcnt <- edgeR::cpm(fit.l$d2, normalized.lib.sizes=TRUE, log=TRUE, prior.count=3)
    }

    type <- match.arg(type, c("both", "uncorrected", "pseudo-corrected"))
    if (type=="both") {
      dw_log("Plotting both uncorrected and pseudo-corrected if applicable...\n")
    }
    if (plots) {
      if (type %in% c("both", "uncorrected")) {
        ## quality control plots for uncorrected data
        dw_log("Quality control plots...\n")
        dw_log("  Uncorrected, normalised data...\n")
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
        dw_log("  Corrected, normalised pseudo-counts...\n")
        ## get pseudo counts for blocked designs (e.g., paired samples or batch factors)
        dw_log("   Calculating pseudo-counts...\n")
        pseudo.counts <- try(diff_expr_pseudo_counts(design=design, d=d, pairs="pairs", disp=disp, do.cpm=TRUE))
        if (is(pseudo.counts, "try-error")) {
          dw_log(" !NOTE! - Pseudo-count calculation failed. Omitting from output...\n")
        } else {
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


    dw_step("@ -- EXTRACTING CONTRASTS --\n")
    if (do.voom) {
      if (!block && !is.null(pairs)) {
        dw_log("  ...for paired samples comparisons (voom)...\n")
        cont <- grep("^groups.+", colnames(fit3$coefficients), value=TRUE)
        dw_log("  ", cont, "\n")
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
                                             hm.p.thr = hm.p.thr,
                                             hm.fdr.thr = hm.fdr.thr,
                                             hm.logfc.thr = hm.logfc.thr,
                                             heatmap.split.expr = heatmap.split.expr,
                                             color.blind.pal = color.blind.pal,
                                             n.pal.cols = n.pal.cols,
                                             color.extremes = color.extremes,
                                             palette.length = palette.length,
                                             anno.color = anno.color,
                                             anno.name = anno.name,
                                             heatmap.main = paste0("limma::", voom.fun),
                                             font.size = label.font.size,
                                             de.plot.base.size = de.plot.base.size,
                                             plots = plots,
                                             lists = lists,
                                             filtered.lists = filtered.lists,
                                             samp.info = samp.info,
                                             samples = "SampleNames",
                                             groups = groups,
                                             sample.plot.names = sample.plot.names)
      }

      if (!is.null(contrasts)) {
        dw_log("  ...for all (remaining) comparisons (voom)...\n")
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
                                             hm.p.thr = hm.p.thr,
                                             hm.fdr.thr = hm.fdr.thr,
                                             hm.logfc.thr = hm.logfc.thr,
                                             heatmap.split.expr = heatmap.split.expr,
                                             color.blind.pal = color.blind.pal,
                                             n.pal.cols = n.pal.cols,
                                             color.extremes = color.extremes,
                                             palette.length = palette.length,
                                             anno.color = anno.color,
                                             anno.name = anno.name,
                                             heatmap.main = paste0("limma::", voom.fun),
                                             font.size = label.font.size,
                                             de.plot.base.size = de.plot.base.size,
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
      dw_log("  ...for paired samples comparisons (GLM)...\n")
        cont <- grep("^groups.+", colnames(fit.l$fit$coefficients), value=TRUE)
        dw_log("  ", cont, "\n")
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
                                             hm.p.thr = hm.p.thr,
                                             hm.fdr.thr = hm.fdr.thr,
                                             hm.logfc.thr = hm.logfc.thr,
                                             heatmap.split.expr = heatmap.split.expr,
                                             color.blind.pal = color.blind.pal,
                                             n.pal.cols = n.pal.cols,
                                             color.extremes = color.extremes,
                                             palette.length = palette.length,
                                             anno.color = anno.color,
                                             anno.name = anno.name,
                                             heatmap.main = "edgeR GLM",
                                             font.size = label.font.size,
                                             de.plot.base.size = de.plot.base.size,
                                             plots = plots,
                                             lists = lists,
                                             filtered.lists = filtered.lists,
                                             samp.info = samp.info,
                                             samples = "SampleNames",
                                             groups = groups,
                                             sample.plot.names = sample.plot.names)
      }

      if (!is.null(contrasts)) {
        dw_log("  ...for all (remaining) comparisons (GLM)...\n")
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
                                             hm.p.thr = hm.p.thr,
                                             hm.fdr.thr = hm.fdr.thr,
                                             hm.logfc.thr = hm.logfc.thr,
                                             heatmap.split.expr = heatmap.split.expr,
                                             color.blind.pal = color.blind.pal,
                                             n.pal.cols = n.pal.cols,
                                             color.extremes = color.extremes,
                                             palette.length = palette.length,
                                             anno.color = anno.color,
                                             anno.name = anno.name,
                                             heatmap.main = "edgeR GLM",
                                             font.size = label.font.size,
                                             de.plot.base.size = de.plot.base.size,
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
      dw_step("@ -- ENRICHMENT ANALYSIS --\n\n")
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
    run_secs <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
    dw_step(sprintf("@ >> RUN COMPLETED IN %s << @\n", fmt_dur(run_secs)))

    ## the run completed normally, so the device-cleanup handler is no longer needed.
    ## NOTE: do NOT call bare on.exit() here - that would also cancel dw_log_end()
    ## and leave the log file connection open.
    on.exit(dw_log_end(), add = FALSE)
    return(out.l)
  }
