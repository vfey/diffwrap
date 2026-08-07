# Utilities for performing various enrichment analyses and visualisations.
# Major function that wraps other enrichment functions is runEnrichmentAnalyses
#
# Author: Meeri Pekkarinen
#*********************************************************************************

#' Helper function for enrichment visualisations: prepare plot legend y coordinates and labels
#' @param scale.minimum Numeric; minimum y coordinates for text labels
#' @param scale.maximum Numeric; maximum y coordinates for text labels
#' @param int.values.for.ticks Not implemented.
#' @details
#' The coordinates are used in a network graph and minimum and maximum scale values correspond to
#'   minimum and maximum fold-changes, by default.
#'
#' @return A named \code{list} with the elements \code{scale.y.coordinates} (the vertical positions of
#'   the legend tick labels) and \code{scale.labels} (the corresponding label values).
#' @examples
#' prepare_scale_for_legend(scale.minimum = -3, scale.maximum = 4)
#' @export
prepare_scale_for_legend = function(scale.minimum, scale.maximum, int.values.for.ticks=NULL){

  scale.minimum = round(scale.minimum, digits = 0)
  scale.maximum = round(scale.maximum, digits = 0)

  int.values.for.ticks = NULL
  if (scale.maximum <= 3) {
    scale.range = seq(scale.minimum,scale.maximum,1)
  }
  else if ((scale.maximum %% 2) == 0) { #even
    scale.range = seq(scale.minimum,scale.maximum,2)
  }
  else { #odd
    scale.range = seq(scale.minimum,scale.maximum, (scale.maximum / 2))
  }

  legend.params = list()
  legend.params$scale.y.coordinates = seq(0.8,1,l = length(scale.range))
  legend.params$scale.labels = scale.range
  #legend.params$scale.labels = c(round(scale.minimum, digits = 2), 0 , round(scale.maximum,
  #                                                                           digits = 2))
  return(legend.params)
}


#' Function for making network visualisation based on enrichment result, DE gene table and thresholding.
#' Can take the input tables either as data frames or Excel files
#' @param enrichment.result data.frame; output of the enrichment tool. In general, a data frame with all
#'   columns needed for the plot.
#' @param DE.result data.frame; fold-change table with corresponding statistics, e.g., the output of \code{topTable()} for
#'   a particular contrast.
#' @param plot.filename character; name of the output image file. The file extension can be omitted and will be added internally.
#' @param show.terms integer; Number of enrichment terms to be plotted, i.e., the first \code{Integer} rows of \code{enrichment.result}.
#' @param logfc.thr \code{numeric}. Fold-change threshold on the log2-scale.
#' @param fdr.thr \code{numeric}. Threshold for adjusted p-values (applied in both filtering of DE-genes and in enrichment results). Default 0.05
#' @param pdf.width,pdf.height Numeric; width and height of the PDF graphics region in inches.
#' @param legend.cex.main,legend.cex.text Numeric; text sizes of legend title and text, given as character expansion
#'   (magnification) relative to the default.
#' @return No return value. Called for its side effect of writing a network visualisation of the
#'   enrichment result to the file given by \option{plot.filename}.
#' @examples
#' \donttest{
#' # builds a network plot from an enrichment result and writes it to 'plot.filename'
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   enr <- data.frame(Description = c("Pathway A", "Pathway B"),
#'                     genes = c("G1,G2,G3", "G3,G4,G5"))
#'   de <- data.frame(gene_symbol = paste0("G", 1:5),
#'                    logFC = c(2.5, -1.8, 1.2, -2.1, 0.9),
#'                    FDR = c(0.001, 0.002, 0.01, 0.003, 0.02))
#'   plot_enrichment_network(enr, de,
#'                           plot.filename = file.path(tempdir(), "network.pdf"))
#' }
#' }
#' @export
plot_enrichment_network <- function(enrichment.result, DE.result, plot.filename,
                                    show.terms = 5,
                                    logfc.thr = NULL,
                                    fdr.thr = NULL,
                                    pdf.width = 11,
                                    pdf.height = 8,
                                    legend.cex.main = 0.8,
                                    legend.cex.text = 0.7) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package ", sQuote("igraph"),
         " must be installed to draw the enrichment network plot.", call. = FALSE)
  }

  if (typeof(enrichment.result) == "character") {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package ", sQuote("readxl"),
           " must be installed to read enrichment results from an Excel file.", call. = FALSE)
    }
    # Read in enrichment result table
    enrichment.table <- data.frame(readxl::read_excel(enrichment.result))
    DE.table <- data.frame(readxl::read_excel(DE.result))
  } else{
    enrichment.table = enrichment.result
    DE.table = DE.result
  }


  geneSymbolColumn = dw_find_col(names(DE.table), "ymbol", "gene symbol")
  logFoldChangeColumn = dw_find_col(names(DE.table), "foldch|logfc", "log fold-change")
  adjPvalColumn = dw_find_col(names(DE.table), "^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", "FDR")
  termColumn = dw_find_col(names(enrichment.table), "description|^term_name", "term description")
  DEGcolumn = dw_find_col(names(enrichment.table), "^degs|^genes|^intersection$", "gene set")
  dw_log("         Detected following columns in data inputted for network visualisation: \n")
  dw_log("         DE (log) fold changes: ", logFoldChangeColumn,  "\n")
  dw_log("         FDR: ",adjPvalColumn, "\n")
  dw_log("         Genes associated with the term: ", DEGcolumn, "\n")
  dw_log("         Term descriptions used in nodes: ", termColumn, "\n", "\n")
  dw_log("         Filtering the genes (fdr < ", fdr.thr, "and abs. logFC >=", logfc.thr, ") and", "by the existence of symbolic names...\n")

  # Take only n enriched terms (never more than are available, to avoid NA rows)
  enrichment.table <- enrichment.table[seq_len(min(show.terms, nrow(enrichment.table))), , drop = FALSE]
  # Remove rows without HGNC symbol from DE table
  DE.table <- DE.table[!DE.table[[geneSymbolColumn]] == "",]
  DE.table <- DE.table[!is.na(DE.table[[geneSymbolColumn]]),]

  #Filter also by unique gene symbols (SHOULD THIS BE DONE OR NOT?)
  DE.table = DE.table[!duplicated(DE.table[[geneSymbolColumn]]),]


  # Filter DE.table by log2FC
  if (!is.null(logfc.thr)) {
    DE.table <- DE.table[which(abs(DE.table[[logFoldChangeColumn]]) >= logfc.thr), ]

  }

  # Filter DE.table by p-value
  if (!is.null(fdr.thr)) {
    DE.table <- DE.table[which(DE.table[[adjPvalColumn]] < fdr.thr), ]
  }

  # Outlier detection:
  all.foldchanges = DE.table[[logFoldChangeColumn]]

  # ***NOTE*** because genes are already potentially filtered based on logfc, the fc-histogram might have a "gap" in the middle and is not actually a histogram!!
  #print(hist(all.foldchanges, 50))
  # ---> PROBLEM: How to define outliers??? (pure boxplot outlier method i.e. 'boxplot(foldchanges, range = 1.5, plot=FALSE)$out'
  # "thinks" we have two separate distributions and calculates outliers for both of them
  # - thus, it removes also values that are members of "the other" distribution)

  #Current solution: separating negatives and positives and taking the relevant IQR outliers separately
  negatives = all.foldchanges[all.foldchanges < 0]
  positives = all.foldchanges[all.foldchanges > 0]
  lower_limit = quantile(negatives,1/4) - 3*IQR(negatives)  #originally quantile(negatives,1/4) - 1.5*IQR(negatives), but this filters quite a lot
  upper_limit = quantile(positives,3/4) + 3*IQR(positives)

  outliers = all.foldchanges[all.foldchanges > upper_limit | all.foldchanges < lower_limit]
  dw_log("         ", length(outliers), "  outliers (r-boxplot method) found... \n")

  if (length(outliers) > 0) {
    to_be_saturated = all.foldchanges %in% outliers

    negSaturated = to_be_saturated & all.foldchanges < 0
    #print(sum(negSaturated))
    all.foldchanges[negSaturated] <- lower_limit

    posSaturated = to_be_saturated & all.foldchanges > 0
    #print(sum(posSaturated))
    all.foldchanges[posSaturated] <- upper_limit

    #replacing the fold change column with saturated values
    DE.table[[logFoldChangeColumn]] <- all.foldchanges
  }



  # Parse the pathway table to create edges
  edges <- c()
  genes = c()
  for (i in seq_along(enrichment.table[[termColumn]])) {
    dw_log("\n")
    x <- strsplit(as.character(enrichment.table[i, DEGcolumn]), split = ",")[[1]]
    dw_log("         Genes associated with the term '", enrichment.table[i, termColumn], "': ", length(x), "\n")
    x <- x[x %in% DE.table[[geneSymbolColumn]]] # Remove genes that are not in DE table
    dw_log("         Genes associated with the term after removing those not found in filtered DE-table: ", length(x), "\n")
    y <- rep(as.character(enrichment.table[i, termColumn]), times =  length(x))
    edge <- as.vector(rbind(y,x))
    edges <- append(edges, edge)
    genes <- append(genes, x)
  }
  # Generate graph object
  g <- igraph::make_undirected_graph(edges=as.vector(edges))
  genes = unique(genes)

  #print(g)

  # Acquire vertices names
  vertices <- igraph::V(g)$name
  #print(vertices)

  # Assign categories to vertices
  categories <- c()
  for (i in seq_along(vertices)) {
    if (vertices[i] %in% enrichment.table[[termColumn]]) {
      categories <- append(categories, "pathway")
    } else {
      categories <- append(categories, "gene")
    }
  }
  igraph::V(g)$categories <- categories


  # Based on foldchange, create color palette for the genes
  # Set a separate color for pathways
  palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))
  #genes = unique(genes)

  #foldchanges <- DE.table[[logFoldChangeColumn]] #original plot
  gene_indices = which(DE.table[[geneSymbolColumn]] %in% genes) #upgrade: taking fold changes only for genes associated with the terms to be plotted
  foldchanges <- DE.table[gene_indices,logFoldChangeColumn]


  # Making symmetric range for fold changes
  if (max(foldchanges, na.rm = TRUE) > abs(min(foldchanges, na.rm = TRUE))) {
    foldchanges <- append(foldchanges, -max(foldchanges, na.rm = TRUE))
  } else {
    foldchanges <- append(foldchanges, -min(foldchanges, na.rm = TRUE))
  }
  #print(foldchanges)
  centered_pal <- palette(length(foldchanges) + 1)[as.numeric(cut(foldchanges,breaks = length(foldchanges)))]

  DE.table = DE.table[gene_indices,]
  DE.table = DE.table[order(DE.table[[logFoldChangeColumn]], decreasing = TRUE),]
  DE.table$color <- centered_pal[-length(centered_pal)]

  # Assign graphical properties to their respective vertices
  shapes <- c()
  colors <- c()
  sizes <- c()
  vert.label.sizes <- c()
  vertex.label.dists <- c()
  for (i in seq_along(igraph::V(g)$name)) {
    if (igraph::V(g)$categories[i] == "pathway") {
      #print(vertices[i])
      shapes <- append(shapes, "circle")
      colors <- append(colors, "white")
      sizes <- append(sizes, 7)
      vert.label.sizes <- append(vert.label.sizes, 0.55)
      vertex.label.dists = append(vertex.label.dists, 0)
      igraph::V(g)$name[i] = gsub('(.{1,15})(\\s|$)', '\\1\n', igraph::V(g)$name[i])
      #print(gsub('(.{1,15})(\\s|$)', '\\1\n', igraph::V(g)$name[i]))
      #print(igraph::V(g)$name[i])
    } else {
      shapes <- append(shapes, "circle")
      colors <- append(colors, DE.table[DE.table[[geneSymbolColumn]] == igraph::V(g)$name[i],"color"])
      sizes <- append(sizes, 7)
      vert.label.sizes <- append(vert.label.sizes, 0.6)
      vertex.label.dists = append(vertex.label.dists, -0.2)

    }
  }
  #print(length(colors))
  igraph::V(g)$shapes <- shapes
  igraph::V(g)$colors <- colors
  igraph::V(g)$sizes <- sizes
  igraph::V(g)$vert.label.sizes <- vert.label.sizes
  igraph::V(g)$vertex.label.dists <- vertex.label.dists


  # Produce the plot
  #tiff(paste0(plot.filename, ".tiff"), width = 10, height = 10, units = 'cm', res = 300)
  pdf(paste0(plot.filename, ".pdf"), width = pdf.width, height = pdf.height)
  # the device and the layout()/par() state set on it are owned by this function; the immediate
  # on.exit() closes it even if the plotting below fails, so the caller's device is never touched
  net.dev <- grDevices::dev.cur()
  on.exit(dw_dev_off(net.dev), add = TRUE)
  layout(matrix(1:2, ncol = 2), widths = c(0.85*pdf.width,0.15*pdf.width),heights = c(1,0.25))
  plot(g, vertex.label.color = "black", vertex.size = igraph::V(g)$sizes, vertex.frame.color = "white",
       vertex.color = igraph::V(g)$colors, vertex.label.cex = igraph::V(g)$vert.label.sizes, vertex.shape = igraph::V(g)$shapes,
       vertex.label.dist = igraph::V(g)$vertex.label.dists,
       layout = igraph::layout.graphopt(g, spring.length = 1.3, spring.constant = 1.3), vertex.label.family = "sans", vertex.label.font = 2)
  # Add the legend
  legend_image <- as.raster(matrix(rev(palette(100)), ncol = 1))

  # Add legend title
  plot(c(0,0.5), c(0,1), type = 'n', axes = FALSE, xlab = '', ylab = '', main = "Log2 foldchange", cex.main = legend.cex.main)

  # Add legend labels
  legend.element = prepare_scale_for_legend(min(foldchanges), max(foldchanges))
  text(x = 0.3, y = legend.element$scale.y.coordinates, cex = legend.cex.text, labels = legend.element$scale.labels)
  # Add legend image
  rasterImage(legend_image, 0.1, 0.8, 0.2,1)

  # the device is closed by the on.exit() handler above; return nothing so the value is not
  # auto-printed when the function is called at top level (e.g. in examples)
  invisible(NULL)
}



#' Helper function for formatting the gene ID column of enrichment data frame.
#' @param result Data.frame; a data frame with Ensembl Gene IDs in one (only one) column
#' @param species Character of length one; name of the species the IDs refer to. Only "human" and "mouse" are supported.
#' @param which.split Character; separator used in the ID column of \code{result}. See \emph{details}.
#' @details
#' Enrichment tools report all genes annotated to a certain term in one string, separated by comma or similar. The function
#'   splits each row into individual IDs before converting and later re-collapses converted IDs.
#'
#' @return The input enrichment result with the Ensembl gene identifiers in the gene column replaced by
#'   the corresponding gene symbols.
format_ensembl_ids_annotated_to_term <- function(result, species, which.split = ",")
{
  # locate the column holding Ensembl gene IDs; works for human (ENSG...), mouse (ENSMUSG...)
  # and other species, unlike a bare grep for the literal "ENSG"
  is.ens <- vapply(result,
                   function(col) any(grepl("ENS[A-Z]*G[0-9]", as.character(col))),
                   logical(1))
  gene.col <- colnames(result)[is.ens]
  if (length(gene.col) != 1L) {
    stop("Could not unambiguously identify the Ensembl gene ID column in the enrichment result ",
         "(found ", length(gene.col), " candidate column(s)).", call. = FALSE)
  }

  for (i in seq_along(result[[gene.col]])) {
    genes <- unlist(strsplit(result[[gene.col]][i], split = which.split))
    #print(head(genes))
    syms <- convertid::convertId2(genes, species)
    syms[is.na(syms)] <- genes[is.na(syms)]

    new_names <- paste(syms,collapse = ",") #writing with comma-separator
    #print(head(new_names))
    result[[gene.col]][i] <- new_names
  }
  return(result)
}


#' Wrapper for executing various enrichment analyses
#' @description \command{runEnrichmentAnalyses} enables the auto-run of some over-representation analysis (ORA) and
#' gene set enrichment analysis (GSEA) functions for the output of diffExpr main wrapper.
#' Functions in various R-packages (clusterProfiler, topGO, gProfileR)
#' are integrated. Currently supports human or mouse!
#' @param diffr.wrapper.output \code{list}.  Nested list system produced by diffrExpr-wrapper.
#' Has to contain element "contrasts" that contains contrast-specific expression tables
#' @param analysis.name \code{character}. Descriptive character-tag used in output file names
#' @param use.background.from.diffr.output \code{logical}. Whether the (ORA) analyses are run
#' with experiment-specific background obtained from pre-filtered expression matrix or with the default background of the functions (genome)
#' @param use.pval.in.DE.filtering.if.no.sign.fdrs \code{logical} Sometimes no DE genes with significant adjusted p-value is found.
#' In such cases, should uncorrected p-values be used in order to get at least some results
#' @param out.dir \code{character}. Root directory for the resulting subdirectories. Must contain subfolders for contrasts.
#'   Required; no default is used so that nothing is written to the working directory unintentionally.
#' @param species \code{character}. Currently valid options are "human" or "mouse".
#' @param p.thr \code{numeric}. Threshold for un-adjusted p-values (applied in both filtering of DE-genes and in enrichment results,
#'  when relevant (i.e. no significant fdr-entries are found)). Default 0.05
#' @param fdr.thr \code{numeric}. Threshold for adjusted p-values (applied in both filtering of DE-genes and in enrichment results). Default 0.05
#' @param logfc.thr \code{numeric}. Threshold for the log2-fold-change. Defaults to 1.
#' @param do.plot \code{logical}. Whether or not to draw a network plot for the enrichment results. Defaults to \code{FALSE}.
#' @param plot.fdr.thr \code{numeric}. FDR threshold used in the enrichment plot. This may be useful to tweak to produce a more informative plot.
#' Defaults to 0.05.
#' @param plot.logfc.thr \code{numeric}. FC threshold on the log2-scale used in the enrichment plot. Defaults to 1.
#' @param plot.num.terms \code{integer}. Number of terms shown in the plot. Defaults to 5.
#' @param enrichment.methods \code{character}. Enrichment methods to be run. One or more of the following: c("clusterProfilerGO", "clusterProfilerKEGG", "gProfileR", "topGO")
#' @param clusterProfilerGO.params \code{list}. Method-specific parameters for clusterProfilerGO. One or more of the following (default values shown
#' and used for all such elements not given in the call):
#' analysis.approach="ORA", do.similarity.filtering=F,min.gene.set.size=10,max.gene.set.size=1000, ontology="BP", min.overlap=2,p.adjust.method="BH".
#' Analysis approach can be "ORA" or "KEGG". If do.similarity.filtering is set to TRUE, clusterProfiler::simplify() is run.
#' @param clusterProfilerKEGG.params \code{list}.
#' @param gProfileR.params \code{list}.
#' @param topGO.params \code{list}.
#'
#' @return A list of all relevant objects generated in the course of the enrichment analyses
###############################################################################

#' @examples
#' \donttest{
#' # needs annotation packages; run on real data with mappable gene IDs
#' if (requireNamespace("clusterProfiler", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
#'     requireNamespace("AnnotationDbi", quietly = TRUE)) {
#'   out.dir <- file.path(tempdir(), "diffwrap_demo")
#'   dir.create(out.dir, showWarnings = FALSE)
#'   res <- diffExpr(diffwrap_counts, diffwrap_samp_info, samples = "SampleName",
#'                   groups = "Group", control = "control", analysis.name = "demo",
#'                   out.dir = out.dir, enr.do = FALSE)
#'   enr <- runEnrichmentAnalyses(res, analysis.name = "demo", out.dir = out.dir,
#'                                species = "human", enrichment.methods = "clusterProfilerGO")
#' }
#' }
#' @export
runEnrichmentAnalyses <- function(diffr.wrapper.output, analysis.name="enrichment",
                                  use.background.from.diffr.output=TRUE, use.pval.in.DE.filtering.if.no.sign.fdrs=FALSE,
                                  out.dir,
                                  species="human",
                                  p.thr=0.05, fdr.thr=0.05, logfc.thr=1,
                                  do.plot=FALSE, plot.fdr.thr=fdr.thr, plot.logfc.thr=logfc.thr, plot.num.terms=5,
                                  enrichment.methods=c("clusterProfilerGO", "clusterProfilerKEGG", "gProfileR", "topGO"),
                                  clusterProfilerGO.params=list(analysis.approach = "ORA",
                                                                  do.similarity.filtering = FALSE,
                                                                  min.gene.set.size = 10,
                                                                  max.gene.set.size = 1000,
                                                                  ontology = "BP",
                                                                  min.overlap = 2,
                                                                  p.adjust.method = "BH"
                                                                  ),
                                  clusterProfilerKEGG.params = list(analysis.approach = "ORA",
                                                                    min.gene.set.size = 10,
                                                                    max.gene.set.size = 1000,
                                                                    ontology = "BP",
                                                                    min.overlap = 2,
                                                                    p.adjust.method = "BH"
                                                                    ),
                                  gProfileR.params = list(data.sources = "GO:BP", show.only.significant = TRUE,
                                                          measure_underrepresentation = FALSE,
                                                          evidence_codes = TRUE,
                                                          domain_scope = "annotated",
                                                          highlight = TRUE),
                                  topGO.params = list(ontologies.used = c("BP"), org = "hsapiens")
                                  )
{

  # only human and mouse are supported; guard here so a NULL/other value fails with a clear
  # message rather than deep inside clusterProfiler (or at the 'species == "human"' test below)
  if (is.null(species) || length(species) != 1L || !species %in% c("human", "mouse")) {
    stop("Enrichment currently supports species 'human' or 'mouse'; got ",
         if (is.null(species)) "NULL" else sQuote(species), ".", call. = FALSE)
  }

  # test if data packages are installed
  if (species == "human" && !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop(
      "Package \"org.Hs.eg.db\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (species == "mouse" && !requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    stop(
      "Package \"org.Mm.eg.db\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # the enrichment engines and the Excel writer live in 'Suggests'; check the ones the
  # selected 'enrichment.methods' actually need, so a missing package fails early and clearly
  needed <- c(
    if (any(grepl("clusterProfiler", enrichment.methods))) "clusterProfiler",
    if (any(grepl("gProfileR", enrichment.methods)))       "gprofiler2",
    if (any(grepl("topGO", enrichment.methods)))           "topGO",
    "WriteXLS"
  )
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop("The following package(s) must be installed for the requested enrichment analysis: ",
         paste(sQuote(miss), collapse = ", "), ".", call. = FALSE)
  }

  ## TODO: invent a smarter way to do this...
  enrich.resource.terms <- data.frame("Organism" = c("human", "mouse"),
                                     "clusterProfilerGO" = c("org.Hs.eg.db", "org.Mm.eg.db"),
                                     "clusterProfilerKEGG" = c("hsa","mmu"),
                                     "gProfileR" = c("hsapiens", "mmusculus"),
                                     "topGO" = c("org.Hs.eg.db", "org.Mm.eg.db"),
                                     "species4conversion" = c("Human", "Mouse"))
  rownames(enrich.resource.terms) <- c("human", "mouse")

  enrichment_out.l <- list()
  enrichment_out.l$clusterProfiler_GO <- list()
  enrichment_out.l$clusterProfiler_KEGG <- list()
  enrichment_out.l$gProfileR <- list()
  enrichment_out.l$topGO <- list()


  contrast.names <- names(diffr.wrapper.output$contrasts)
  dw_log("Performing enrichment analyses for ", length(contrast.names), " comparisons: \n")
  dw_log("  ", contrast.names, " \n")

  dat <- diffr.wrapper.output$contrasts[[1]] ## extracting appropriate colnames using the first contrast
  pv.col <- names(dat)[grep("^p\\.{0,1}val[e-u]{0,2}$", tolower(names(dat)))]
  fdr.col <- dw_find_col(names(dat), "^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", "FDR")
  fc.col <- dw_find_col(names(dat), "^logfc$|fold$", "log fold-change")
  entrez.col <- names(dat)[grep("entrez", tolower(names(dat)))]


  background.genes <- ""
  if (use.background.from.diffr.output) {
    # drop the htseq-/STAR-style summary rows if present (keep everything that is NEITHER)
    background.genes <-
      rownames(dat)[!(rownames(dat) %in% c("N_multimapping", "N_noFeature"))]
    dw_log("The genes of full expression table (", length(background.genes), ") is used as the background (expect in GSEA-analyses, where no background is used)... \n")
  }

  for (contrast in contrast.names) {

    contr.out.dir  <- dir(out.dir, pattern = paste0("^", contrast, "$"), full.names = TRUE)
    #Creating sub-folders
    dw_log("  Creating per-enrichment sub-folders...", "\n")
    lapply(enrichment.methods, function(method) dir.create(file.path(contr.out.dir, method), showWarnings = FALSE))

    dw_log("***", contrast, "*** \n")
    de_table <- diffr.wrapper.output$contrasts[[contrast]]

    dw_log("   Filtering the DE genes (fdr < ", fdr.thr, "and abs. logFC >=", logfc.thr, ")...\n")
    filtered_de_table <- de_table[which(de_table[[fdr.col]] < fdr.thr & abs(de_table[[fc.col]]) >= logfc.thr), ]

    #if there are no entries with significant fdr, then filter by p.value
    if (!(nrow(filtered_de_table) > 0) & use.pval.in.DE.filtering.if.no.sign.fdrs) {
      dw_log("      No entries with significant fdr. P-values used instead...\n")
      filtered_de_table <- de_table[which(de_table[[pv.col]] < p.thr & abs(de_table[[fc.col]]) >= logfc.thr), ]
    }

    input.genes <- rownames(filtered_de_table)
    if (length(input.genes) == 0) {
      dw_log("  ", length(input.genes), " genes considered significant for the analysis in contrast", contrast,". Skipping the contrast...\n")
      next
    }

    dw_log("  ", length(input.genes), " genes used as input (expect in GSEA-analyses, where all measured genes is used) \n")

    for (method in enrichment.methods) {

     if (method == "clusterProfilerGO") {

       ## Replacing missing parameters
       default.params <- list(analysis.approach = "ORA", do.similarity.filtering = FALSE, min.gene.set.size = 10,
                             max.gene.set.size = 1000, ontology = "BP", min.overlap = 2, p.adjust.method = "BH")
       missing.params <- default.params[names(default.params)[!(names(default.params) %in% names(clusterProfilerGO.params))]]
       clusterProfilerGO.params <- c(clusterProfilerGO.params, missing.params)

       dw_log("   Performing GO BP enrichment with clusterProfiler... \n")
       dw_log("   ", "\n")
       dw_log_obj(unlist(clusterProfilerGO.params))
       org.db <- as.character(enrich.resource.terms[species, method])

       if (clusterProfilerGO.params$analysis.approach == "ORA") {

         ordered.query <- FALSE
         genes <- input.genes

       }
       else {

         ordered.query <- TRUE

         # For GSEA, an ordered gene list of ALL genes must be prepared:
         #### feature 1: numeric vector, named by the (unique) Ensembl gene IDs
         geneList <- de_table[[fc.col]]
         names(geneList) <- as.character(rownames(de_table))
         ## feature 2: drop genes without a fold change, then order decreasingly
         geneList <- geneList[!is.na(geneList)]
         geneList <- sort(geneList, decreasing = TRUE)
         genes <- geneList

       }
       dw_log("  Running tests...\n")
       result <- run_clusterProfiler_GO(input_genes = genes,
                                        background_genes = background.genes,
                                        ordered_query = ordered.query,
                                        OrgDb = org.db,
                                        pvalueCutoff = fdr.thr,
                                        ontology = clusterProfilerGO.params$ontology,
                                        min_set_size = clusterProfilerGO.params$min.gene.set.size,
                                        max_set_size = clusterProfilerGO.params$max.gene.set.size,
                                        min_overlap = clusterProfilerGO.params$min.overlap,
                                        pAdjustMethod = clusterProfilerGO.params$p.adjust.method,
                                        similarity_filtering = clusterProfilerGO.params$do.similarity.filtering)

       #Saving the table, if relevant
       method.dir <- dir(contr.out.dir, pattern = paste0("^",method), full.names = TRUE)
       if (length(result) > 0 && is.data.frame(result) && nrow(result) > 0) {

         filename <- paste0(analysis.name, ".", contrast,".", method, ".",
                            clusterProfilerGO.params$analysis.approach, ".",
                            clusterProfilerGO.params$ontology,".xls" )
         full.filename <- file.path(method.dir, filename)
         dw_log("      Saving ", method, " result table into ",  full.filename, "...\n")
         WriteXLS::WriteXLS(result, ExcelFileName = full.filename, SheetNames = NULL, BoldHeaderRow = TRUE)

         #Making the graph visualisation
         if (do.plot) {
           graph.name = gsub(".xls", ".network", full.filename, fixed = TRUE)
           dw_log("      Saving ", method, " network with the name ", graph.name, ".pdf", "...\n")
           plot_enrichment_network(enrichment.result = result, DE.result = de_table,
                                   plot.filename = graph.name, show.terms = plot.num.terms,
                                   logfc.thr = plot.logfc.thr, fdr.thr = plot.fdr.thr)
         }
       } else {
         result <- "No significant enrichment found"
       }

       enrichment_out.l$clusterProfiler_GO[[contrast]] <- result

       }

      if (method == "clusterProfilerKEGG") {

        ## Replacing missing parameters
        method.dir <- dir(contr.out.dir, pattern = paste0("^",method), full.names = TRUE)

        default.params <- list(analysis.approach = "ORA", min.gene.set.size = 10, max.gene.set.size = 1000,
                              ontology = "BP", min.overlap = 2, p.adjust.method = "BH")
        missing.params <- default.params[names(default.params)[!(names(default.params) %in% names(clusterProfilerKEGG.params))]]
        clusterProfilerKEGG.params <- c(clusterProfilerKEGG.params, missing.params)


        dw_log("   Performing KEGG enrichment with clusterProfiler... \n")
        dw_log("   ", "\n")
        dw_log_obj(unlist(clusterProfilerKEGG.params))
        org <- as.character(enrich.resource.terms[species, method])

        if (length(entrez.col) == 0) {
          dw_log("      No entrez IDs found  from the data. skipping KEGG-enrichment...\n")
          next   # skip only KEGG for this contrast, not the remaining enrichment methods
        }

        background.gene.entrez <- ""
        if (clusterProfilerKEGG.params$analysis.approach == "ORA") {

          if (length(background.genes)  > 1) {
            background.gene.entrez <- as.character(de_table[[entrez.col]][!is.na(de_table[[entrez.col]])])
            dw_log("      ",length(background.gene.entrez), "/", length(background.genes), " of the all measured genes having corresponding entrez id used as background...\n")

          }
          else{
            dw_log("      Using default background in KEGG-enrichment...\n")
          }

          input.gene.entrez <- as.character(filtered_de_table[[entrez.col]][!is.na(filtered_de_table[[entrez.col]])])
          dw_log("      ",length(input.gene.entrez), "/", length(input.genes), " of the input genes having corresponding entrez id used for the analysis...\n")


          ordered.query <- FALSE
          genes <- as.character(input.gene.entrez)

        }
        else {

          ordered.query <- TRUE

          # For GSEA, an ordered gene list of ALL genes must be prepared, named by Entrez ID.
          # Several Ensembl IDs can map to the same Entrez ID, so deduplicate by Entrez (keeping
          # the most extreme fold change) to avoid duplicate names, which gseKEGG() rejects.
          ent <- as.character(de_table[[entrez.col]])
          fc  <- de_table[[fc.col]]
          keep <- !is.na(ent) & !is.na(fc)
          ent <- ent[keep]; fc <- fc[keep]
          ord <- order(abs(fc), decreasing = TRUE)   # most extreme fold change first
          ent <- ent[ord]; fc <- fc[ord]
          dup <- duplicated(ent)                      # keep the first (largest |logFC|) per Entrez
          geneList <- fc[!dup]
          names(geneList) <- ent[!dup]
          geneList <- sort(geneList, decreasing = TRUE)

          genes <- geneList

        }

        dw_log("  Running tests...\n")
        result <- run_clusterProfiler_KEGG(input_genes = genes,
                                              background_genes = background.gene.entrez,
                                              ordered_query = ordered.query,
                                              organism = org,
                                              pvalueCutoff = fdr.thr,
                                              min_set_size = clusterProfilerKEGG.params$min.gene.set.size,
                                              max_set_size = clusterProfilerKEGG.params$max.gene.set.size,
                                              min_overlap = clusterProfilerKEGG.params$min.overlap,
                                              pAdjustMethod = clusterProfilerKEGG.params$p.adjust.method)

        #Saving the table, if relevant
        method.dir <- dir(contr.out.dir, pattern = paste0("^",method), full.names = TRUE)
        if (length(result) > 0 && is.data.frame(result) && nrow(result) > 0) {


          gene.col <- colnames(result)[grepl("^DEGs", colnames(result))] #TODO: invent more robust approach to this?
          #First, fetching corresponding Ensembl IDs for entrez-IDs
          for (i in seq_along(result[[gene.col]])) {
            #print(result[[gene.col]])
            genes <- unlist(strsplit(result[[gene.col]][i], split = ','))
            #print(head(genes))
            ensembls <- rownames(dat)[dat[[entrez.col]] %in% genes]
            new_names <- paste(ensembls,collapse = ",") #writing with comma-separator
            #print(head(new_names))
            result[[gene.col]][i] <- new_names
          }
          #Then, converting the Ensembl IDs into symbolic names
          spec.name <- as.character(enrich.resource.terms[species, "species4conversion"])
          result <- format_ensembl_ids_annotated_to_term(result, spec.name)

          filename <- paste0(analysis.name, ".", contrast,".", method, ".",
                             clusterProfilerKEGG.params$analysis.approach, ".",
                             clusterProfilerKEGG.params$ontology,".xls" )
          full.filename <- file.path(method.dir, filename)
          dw_log("      Saving ", method, " result table into ",  full.filename, "...\n")
          WriteXLS::WriteXLS(result, ExcelFileName = full.filename, SheetNames = NULL, BoldHeaderRow = TRUE)

          #Making the graph visualisation
          if (do.plot) {
            graph.name = gsub(".xls", ".network", full.filename, fixed = TRUE)
            dw_log("      Saving ", method, " network with the name ", graph.name, ".pdf", "...\n")
            plot_enrichment_network(enrichment.result = result, DE.result = de_table,
                                    plot.filename = graph.name, show.terms = plot.num.terms,
                                    logfc.thr = plot.logfc.thr, fdr.thr = plot.fdr.thr)
          }
        } else {
          result <- "No significant enrichment found"
        }

        enrichment_out.l$clusterProfiler_KEGG[[contrast]] <- result
      }

      if (method == "gProfileR") {

        dw_log("   Performing GO BP enrichment with gProfileR... \n")
        org <- as.character(enrich.resource.terms[species, method])
        dw_log(paste0("Organism: ",org), "\n")
        dw_log(paste0("Data sources: ",gProfileR.params$data.sources), "\n")

        dw_log("  Running tests...\n")
        results <- run_gprofiler(input.genes, background.genes,
                                organism = org,
                                data_sources = gProfileR.params$data.sources,
                                show_only_significant = gProfileR.params$show.only.significant,
                                measure_underrepresentation = gProfileR.params$measure_underrepresentation,
                                evidence_codes = gProfileR.params$evidence_codes,
                                domain_scope = gProfileR.params$domain_scope,
                                highlight = gProfileR.params$highlight)

        # Formatting and saving the table, if relevant
        method.dir <- dir(contr.out.dir, pattern = paste0("^",method), full.names = TRUE) # detecting output directory
        # getting 'result' data frame from gProfiler results list
        result <- results$result
        if (length(result) > 0 && is.data.frame(result) && nrow(result) > 0) {

          spec.name <- as.character(enrich.resource.terms[species, "species4conversion"])
          result <- format_ensembl_ids_annotated_to_term(result, spec.name)

          filename <- paste0(analysis.name, ".", contrast,".", method, ".", gProfileR.params$data.sources,".xls" )
          full.filename <- file.path(method.dir, filename)
          dw_log("      Saving ", method, " result table into ",  full.filename, "...\n")
          WriteXLS::WriteXLS(result, ExcelFileName = full.filename, SheetNames = NULL, BoldHeaderRow = TRUE)

          #Making the graph visualisation
          if (do.plot) {
            graph.name = gsub(".xls", ".network", full.filename, fixed = TRUE)
            dw_log("      Saving ", method, " network with the name ", graph.name, ".pdf", "...\n")
            plot_enrichment_network(enrichment.result = result, DE.result = de_table,
                                    plot.filename = graph.name, show.terms = plot.num.terms,
                                    logfc.thr = plot.logfc.thr, fdr.thr = plot.fdr.thr)
          }
        }
        else {
          result <- "No significant enrichment found"
        }

        enrichment_out.l$gProfileR[[contrast]] = result

      }

      if (method == "topGO") {

        method.dir <- dir(contr.out.dir, pattern = paste0("^",method), full.names = TRUE)

        dw_log("Performing topGO.... \n")
        dw_log(paste0("Ontologies used: ", topGO.params$ontologies.used), "\n")
        dw_log(paste0("Organism: ", topGO.params$org), "\n")

        org.db <- as.character(enrich.resource.terms[species, method])
        dw_log("  Running tests...\n")
        result <- run.topGO(background = background.genes, foreground = input.genes,ontologies =  topGO.params$ontologies.used, organism = org.db)

        #Formatting and saving the table, if relevant
        if (length(result) > 0 && is.data.frame(result) && nrow(result) > 0) {

          spec.name <- as.character(enrich.resource.terms[species, "species4conversion"])
          result <- format_ensembl_ids_annotated_to_term(result, spec.name)

          filename <- paste0(analysis.name, ".", contrast,".", method, ".", topGO.params$ontologies.used,".xls" )
          full.filename <- file.path(method.dir, filename)
          dw_log("      Saving ", method, " result table into ",  full.filename, "...\n")
          WriteXLS::WriteXLS(result, ExcelFileName = full.filename, SheetNames = NULL, BoldHeaderRow = TRUE)

          #Making the graph visualisation
          if (do.plot) {
            graph.name = gsub(".xls", ".network", full.filename, fixed = TRUE)
            dw_log("      Saving ", method, " network with the name ", graph.name, ".pdf", "...\n")
            plot_enrichment_network(enrichment.result = result, DE.result = de_table,
                                    plot.filename = graph.name, show.terms = plot.num.terms,
                                    logfc.thr = plot.logfc.thr, fdr.thr = plot.fdr.thr)
          }
        }
        else {
          result <- "No significant enrichment found"
        }

        enrichment_out.l$topGO[[contrast]] <- result
      }
    }
  }
  return(enrichment_out.l)
}



