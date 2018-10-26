# Utilities for for performing various enrichment analyses and visualisations.
# Major function that wraps other enrichment functions is runEnrichmentAnalyses
# 
# Author: Meeri Pekkarinen
#*********************************************************************************
library(WriteXLS)
library(medseqr)
library(igraph)
library(RColorBrewer)
library(readxl)




#Helper function for enrichment visualisations:
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


#Function for making network visualisation based on enrichment result, DE gene table and thresholding. 
#Can take the input tables either as dataframes or excel-files
plot_enrichment_network = function(enrichment.result, DE.result, plot.filename, show.terms = 5, logfc.thr = NULL, fdr.thr = NULL, 
                                   pdf.width = 20, pdf.heigth = 15, vertex.cex.label = 1.3, legend.cex.main = 2.0, legend.cex.text = 1.5) {
  
  if (typeof(enrichment.result) == "character") {
    # Read in enrichment result table
    enrichment.table <- data.frame(read_excel(enrichment.result))
    DE.table <- data.frame(read_excel(DE.result)) 
  } else{
    enrichment.table = enrichment.result
    DE.table = DE.result
  }
  
  
  geneSymbolColumn = names(DE.table)[grep("ymbol", tolower(names(DE.table)))]
  logFoldChangeColumn = names(DE.table)[grep("foldch|logfc", tolower(names(DE.table)))]
  adjPvalColumn = names(DE.table)[grep("^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", tolower(names(DE.table)))]
  termColumn = names(enrichment.table)[grep("description", tolower(names(enrichment.table)))]
  DEGcolumn =   names(enrichment.table)[grep("^degs|^genes", tolower(names(enrichment.table)))]
  cat("         Detected following columns in data inputted for network visualisation: \n")
  cat("         DE (log) fold changes: ", logFoldChangeColumn,  "\n")
  cat("         FDR: ",adjPvalColumn, "\n")
  cat("         Genes associated with the term: ", DEGcolumn, "\n")
  cat("         Term descriptions used in nodes: ", termColumn, "\n", "\n")
  cat("         Filtering the the genes (fdr < ", fdr.thr, "and abs. logFC >=", logfc.thr, ") and", "\n",
                "by the existence of symbolic names...\n")
  
  # Take only n enriched terms
  enrichment.table <- enrichment.table[1:show.terms,]
  # Remove rows without HGNC symbol from DE table
  DE.table <- DE.table[!DE.table[[geneSymbolColumn]] == "",]
  DE.table <- DE.table[!is.na(DE.table[[geneSymbolColumn]]),]
  
  
  # Filter DE.table by log2FC
  if (!is.null(logfc.thr)) {
    DE.table <- DE.table[abs(DE.table[[logFoldChangeColumn]]) > logfc.thr, ]
    
  }
  
  # Filter DE.table by p-value
  if (!is.null(fdr.thr)) {
    DE.table <- DE.table[DE.table[[adjPvalColumn]] < fdr.thr, ]
  }
  
  # Parse the pathway table to create edges
  edges <- c()
  genes = c()
  for (i in 1:length(enrichment.table[[termColumn]])) {
    cat("\n")
    x <- strsplit(as.character(enrichment.table[i, DEGcolumn]), split = ",")[[1]]
    cat("         Genes associated with the term '", enrichment.table[i, termColumn], "': ", length(x), "\n")
    x <- x[x %in% DE.table[[geneSymbolColumn]]] # Remove genes that are not in DE table
    cat("         Genes associated with the term after removing those not found in filtered DE-table: ", length(x), "\n")
    y <- rep(as.character(enrichment.table[i, termColumn]), times =  length(x))
    edge <- as.vector(rbind(y,x))
    edges <- append(edges, edge)
    genes <- append(genes, x)
  }
  # Generate graph object
  g <- graph(edges, directed = F)
  #print(g)
  
  # Acquire vertices names
  vertices <- V(g)$name
  #print(vertices)
  
  # Assign categories to vertices
  categories <- c()
  shapes <- c()
  for (i in 1:length(vertices)) {
    if (vertices[i] %in% enrichment.table[[termColumn]]) {
      categories <- append(categories, "pathway")
      shapes <- append(shapes, "circle")
    } else {
      categories <- append(categories, "gene")
      shapes <- append(shapes, "circle")
    }
  }
  V(g)$categories <- categories
  V(g)$shapes <- shapes
  
  # Based on foldchange, create color palette for the genes
  # Set a separate color for pathways
  palette <- colorRampPalette(rev(brewer.pal(11,"RdBu")))
  #genes = unique(genes)
  
  #foldchanges <- DE.table[[logFoldChangeColumn]] #origiginal plot
  gene_indices = which(DE.table[[geneSymbolColumn]] %in% genes) #upgrade: taking fold changes only for genes associated with the terms to be plotted 
  foldchanges <- DE.table[gene_indices,logFoldChangeColumn]
  
  
  # ***NOTE*** because genes are filtered based on logfc, the fc-histogram has "gap" in the middle and is not actually a histogram!!
  #print(hist(foldchanges, 50)) 
  # ---> PROBLEM: How to define outliers??? (pure boxplot outlier method i.e. 'boxplot(foldchanges, range = 1.5, plot=FALSE)$out' 
  # "thinks" we have two separate distributions and calculates outliers for both of them
  # - thus, it removes also values that are members of "the other" distribution)
  
  #Current solution: separating negatives and positives and taking the relevant IQR outliers separately
  negatives = foldchanges[foldchanges < 0]
  positives = foldchanges[foldchanges > 0]
  lower_limit = quantile(negatives,1/4) - 3*IQR(negatives)  #originally quantile(negatives,1/4) - 1.5*IQR(negatives), but this filters quite a lot
  upper_limit = quantile(positives,3/4) + 3*IQR(positives)
  
  outliers = foldchanges[foldchanges > upper_limit | foldchanges < lower_limit]
  cat("         ", length(outliers), "  outliers (r-boxplot method) found... \n")
  if (length(outliers) > 0) {
    to_be_removed = which(foldchanges %in% outliers)
    
    gene_indices = gene_indices[-to_be_removed]
    foldchanges = foldchanges[-to_be_removed]
  }
  
  
  # Making symmetric range for fold changes
  if (max(foldchanges, na.rm = T) > abs(min(foldchanges, na.rm = T))) {
    foldchanges <- append(foldchanges, -max(foldchanges, na.rm = T))
  } else {
    foldchanges <- append(foldchanges, -min(foldchanges, na.rm = T))
  }
  #print(foldchanges)
  centered_pal <- palette(length(foldchanges) + 1)[as.numeric(cut(foldchanges,breaks = length(foldchanges)))]
  
  DE.table = DE.table[gene_indices,]
  DE.table$color <- centered_pal[-length(centered_pal)]
  
  # Assign colors to their respective vertices
  colors <- c()
  for (i in 1:length(V(g)$name)) {
    if (V(g)$categories[i] == "pathway") {
      colors <- append(colors, "white")
    } else {
      colors <- append(colors, DE.table[DE.table[[geneSymbolColumn]] == V(g)$name[i],"color"])
    }
  }
  V(g)$colors <- colors
  
  # Produce the plot
  #tiff(paste0(plot.filename, ".tiff"), width = 10, height = 10, units = 'cm', res = 300)
  pdf(paste0(plot.filename, ".pdf"), width = pdf.width, height = pdf.heigth)
  layout(matrix(1:2, ncol = 2), widths = c(4,1),heights = c(1,1))  
  plot(g, vertex.label.color = "black", vertex.size = 10, vertex.frame.color = "white",
       vertex.color = V(g)$colors, vertex.label.cex = vertex.cex.label, vertex.shape = V(g)$shapes,
       layout = layout_nicely, vertex.label.font = 2, vertex.label.family = "sans")
  
  # Add the legend
  legend_image <- as.raster(matrix(rev(palette(100)), ncol = 1))
  
  # Add legend title
  plot(c(0,0.5), c(0,1), type = 'n', axes = F, xlab = '', ylab = '', main = "Log. foldchange", cex.main = legend.cex.main)
  
  # Add legend labels
  legend.element = prepare_scale_for_legend(min(foldchanges), max(foldchanges))
  text(x = 0.3, y = legend.element$scale.y.coordinates, cex = legend.cex.text, labels = legend.element$scale.labels)
  # Add legend image                                                                                                                              
  rasterImage(legend_image, 0.1, 0.8, 0.2,1)
  
  dev.off()
}



#Helper function for formatting the gene column of DAVID and gProfileR enrichment dataframe. 
#Can be used even when biomart is not run within the full pipeline. Requires medseqr!
format_ensembl_ids_annotated_to_term <- function(result, organism.term, which.split = ",") 
{
  gene.col <- colnames(result)[grepl("ENS", result)]
  
  for (i in 1:length(result[[gene.col]])) {
    genes <- unlist(strsplit(result[[gene.col]][i], split = which.split))
    #print(head(genes))
    syms <- convertId2(genes, organism.term)
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
#' Functions in various R-packages (clusterProfiler, topGO, gProfileR and RDAVIDWebService) 
#' are integrated. Currently supports human or mouse!
#' @param diffr.wrapper.output \code{list}.  Nested list system produced by diffrExpr-wrapper. 
#' Has to contain element "contrasts" that contains contrast-specific expression tables
#' @param analysis.name \code{character}. Descriptive character-tag used in output file names
#' @param use.background.from.diffr.output \code{logical}. Whether the (ORA) analyses are run 
#' with experiment-specific background obtained from pre-filtered expression matrix or with the default background of the functions (genome)
#' @param use.pval.in.DE.filtering.if.no.sign.fdrs \code{logical} Sometimes no DE genes with significant adjusted p-value is found. 
#' In such cases, should uncorrected p-values be used in order to get at least some results
#' @param out.dir \code{character}. Root directory for the resulting subdirectories 
#' @param species \code{character}. Currently valid options are "human" or "mouse".
#' @param p.thr \code{numeric}. Threshold for un-adjusted p-values (applied in both filtering of DE-genes and in enrichment results,
#'  when relevant (i.e. no significant fdr-entries are found)). Default 0.05 
#' @param fdr.thr \code{numeric}. Threshold for adjusted p-values (applied in both filtering of DE-genes and in enrichment results). Default 0.05
#' @param logfc.thr \code{numeric}. 
#' @param enrichment.methods \code{character}. Enrichment methods to be run. One or more of the following: c("clusterProfilerGO", "clusterProfilerKEGG","DAVID", "gProfileR", "topGO")
#' @param clusterProfilerGO.params \code{list}. Method-specific parameters for clusterProfilerGO. One or more of the following (default values shown 
#' and used for all such elements not given in the call):
#' analysis.approach="ORA", do.similarity.filtering=F,min.gene.set.size=10,max.gene.set.size=1000, ontology="BP", min.overlap=2,p.adjust.method="BH". 
#' Analysis approach can be "ORA" or "KEGG". If do.similarity.filtering is set to TRUE, clusterProfiler::simplify() is run.
#' @param clusterProfilerKEGG.params \code{list}.
#' @param david.params \code{list}.
#' @param gProfileR.params \code{list}.
#' @param topGO.params \code{list}.
#'
#' @details DAVID approach requires a registered email-address, correct java-version and an url configured with the settings. 
#' url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"
#' @return A list of all relevant objects generated in the course of the enrichment analyses
###############################################################################


runEnrichmentAnalyses <- function(diffr.wrapper.output, analysis.name="enrichment", 
                                  use.background.from.diffr.output=TRUE, use.pval.in.DE.filtering.if.no.sign.fdrs=FALSE,
                                  out.dir=NULL,
                                  species="human",
                                  p.thr=0.05, fdr.thr=0.05, logfc.thr=1, 
                                  enrichment.methods=c("clusterProfilerGO", "clusterProfilerKEGG","DAVID", "gProfileR", "topGO"), 
                                  clusterProfilerGO.params=list(analysis.approach = "ORA",
                                                                  do.similarity.filtering = F,
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
                                  david.params = list(email.address = "", 
                                                      url = "",
                                                      time.out.value = 60000,
                                                      annotation.category = "GOTERM_BP_FAT",
                                                      max.gene.set.size = 1000),
                                  gProfiler.params = list(data.sources = "GO:BP", show.only.significant = TRUE),
                                  topGO.params = list(ontologies.used = c("BP"), org = "org.Hs.eg.db")
                                  ) 
{
  
  ## TODO: invent a smarter way to do this...
  enrich.resource.terms <- data.frame("Organism" = c("human", "mouse"), 
                                     "clusterProfilerGO" = c("org.Hs.eg.db", "org.Mm.eg,db"),
                                     "clusterProfilerKEGG" = c("hsa","mmu"),
                                     "gProfileR" = c("hsapiens", "mmusculus"),
                                     "topGO" = c("org.Hs.eg.db", "org.Mm.eg,db"),
                                     "medseqr" = c("Human", "Mouse"))
  rownames(enrich.resource.terms) <- c("human", "mouse")
  
  enrichment_out.l <- list()
  enrichment_out.l$clusterProfiler_GO <- list()
  enrichment_out.l$clusterProfiler_KEGG <- list()
  enrichment_out.l$DAVID <- list()
  enrichment_out.l$gProfileR <- list()
  enrichment_out.l$topGO <- list()
  
  #Creating subfolders
  lapply(enrichment.methods, function(method) dir.create(file.path(out.dir, method), showWarnings = F))
  
  contrast.names <- names(diffr.wrapper.output$contrasts)
  cat("Performing enrichment analyses for ", length(contrast.names), " comparisons: \n")
  cat("  ", contrast.names, " \n")
  
  dat <- diffr.wrapper.output$contrasts[[1]] ## extracting appropriate colnames using the first contrast
  pv.col <- names(dat)[grep("^p\\.{0,1}val[e-u]{0,2}$", tolower(names(dat)))]
  fdr.col <- names(dat)[grep("^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", tolower(names(dat)))]
  fc.col <- names(dat)[grep("^logfc$|fold$", tolower(names(dat)))]
  entrez.col <- names(dat)[grep("entrez", tolower(names(dat)))]
  
  
  background.genes <- ""
  if (use.background.from.diffr.output) {
    # background.genes <- rownames(dat) ## THIS should be used after filtering is corrected. Until that, two extra elements need to be skipped:
    background.genes <- 
      rownames(dat)[rownames(dat) != "N_multimapping" | rownames(dat) != "N_noFeature"]
    cat("The genes of full expression table (", length(background.genes), ") is used as the background (expect in GSEA-analyses, where no background is used)... \n") 
  }
 
  for (contrast in contrast.names) {

    cat("***", contrast, "*** \n")
    de_table <- diffr.wrapper.output$contrasts[[contrast]]
    
    cat("   Filtering the DE genes (fdr < ", fdr.thr, "and abs. logFC >=", logfc.thr, ")...\n")
    filtered_de_table <- de_table[de_table[[fdr.col]] < fdr.thr & abs(de_table[[fc.col]]) >= logfc.thr,]
    
    #if there are no entries with significant fdr, then filter by p.value
    if (!(nrow(filtered_de_table) > 0) & use.pval.in.DE.filtering.if.no.sign.fdrs) {
      cat("      No entries with significant fdr. P-values used instead...\n")
      filtered_de_table <- de_table[de_table[[pv.col]] < p.thr & abs(de_table[[fc.col]]) >= logfc.thr,]
    }
    
    input.genes <- rownames(filtered_de_table)
    if (length(input.genes) == 0) {
      cat("  ", length(input.genes), " genes considered significant for the analysis in contrast", contrast,". Skipping the contrast...\n")
      next
    }
    
    cat("  ", length(input.genes), " genes used as input (expect in GSEA-analyses, where all measured genes is used) \n")
    
    for (method in enrichment.methods) { 
   
     if (method == "clusterProfilerGO") {
       
       ## Replacing missing parameters
       default.params <- list(analysis.approach = "ORA", do.similarity.filtering = F, min.gene.set.size = 10, 
                             max.gene.set.size = 1000, ontology = "BP", min.overlap = 2, p.adjust.method = "BH")
       missing.params <- default.params[names(default.params)[!(names(default.params) %in% names(clusterProfilerGO.params))]]
       clusterProfilerGO.params <- c(clusterProfilerGO.params, missing.params)
       
       cat("   Performing GO BP enrichment with clusterProfiler... \n")
       cat("   ")
       print(unlist(clusterProfilerGO.params))
       org.db <- as.character(enrich.resource.terms[species, method])
      
       if (clusterProfilerGO.params$analysis.approach == "ORA") {
         
         ordered.query <- FALSE
         genes <- input.genes
         
       }
       else {
         
         ordered.query <- TRUE
         
         # For GSEA,an ordered gene list of ALL genes must be prepared:
         #### feature 1: numeric vector
         geneList <- de_table[[fc.col]]
         ## feature 2: named vector
         names(geneList) <- as.character(rownames(de_table))
         ## feature 3: decreasing order
         geneList <- sort(geneList, decreasing = TRUE)
         genes <- geneList
         
       }
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
       method.dir <- dir(out.dir, pattern = paste0("^",method), full.names = TRUE)
       if (is.data.frame(result) ) {
         
         filename <- paste0(analysis.name, ".", contrast,".", method, ".", 
                            clusterProfilerGO.params$analysis.approach, ".", 
                            clusterProfilerGO.params$ontology,".xls" )
         full.filename <- file.path(method.dir, filename)
         cat("      Saving ", method, " result table into ",  full.filename, "...\n")
         WriteXLS(result, ExcelFileName = full.filename, SheetNames = NULL, BoldHeaderRow = T)
         
         #Making the graph visualisation
         graph.name = gsub(".xls", ".network", full.filename, fixed = TRUE)
         cat("      Saving ", method, " network with the name ", graph.name, ".pdf", "...\n")
         plot_enrichment_network(enrichment.result = result, DE.result = de_table, 
                                 plot.filename = graph.name, show.terms = 5, logfc.thr = 1, fdr.thr = 0.05)
        }
       
       enrichment_out.l$clusterProfiler_GO[[contrast]] <- result
       
       } 
      
      if (method == "clusterProfilerKEGG") {
      
        ## Replacing missing parameters
        method.dir <- dir(out.dir, pattern = paste0("^",method), full.names = TRUE)
      
        default.params <- list(analysis.approach = "ORA", min.gene.set.size = 10, max.gene.set.size = 1000, 
                              ontology = "BP", min.overlap = 2, p.adjust.method = "BH")
        missing.params <- default.params[names(default.params)[!(names(default.params) %in% names(clusterProfilerKEGG.params))]]
        clusterProfilerKEGG.params <- c(clusterProfilerKEGG.params, missing.params)
        
    
        cat("   Performing KEGG enrichment with clusterProfiler... \n")
        cat("   ")
        print(unlist(clusterProfilerKEGG.params))
        org <- as.character(enrich.resource.terms[species, method])
      
        if (length(entrez.col) == 0) {
          cat("      No entrez IDs found  from the data. skipping KEGG-enrichment...\n")
          break 
        }
        
        background.gene.entrez <- ""
        if (clusterProfilerKEGG.params$analysis.approach == "ORA") {
          
          if (length(background.genes)  > 1) {
            background.gene.entrez <- as.character(de_table[[entrez.col]][!is.na(de_table[[entrez.col]])])
            cat("      ",length(background.gene.entrez), "/", length(background.genes), " of the all measured genes having corresponding entrez id used as background...\n")
            
          }
          else{
            cat("      Using default background in KEGG-enrichment...\n")
          }
          
          input.gene.entrez <- as.character(filtered_de_table[[entrez.col]][!is.na(filtered_de_table[[entrez.col]])])
          cat("      ",length(input.gene.entrez), "/", length(input.genes), " of the input genes having corresponding entrez id used for the analysis...\n")
          
          
          ordered.query <- FALSE
          genes <- as.character(input.gene.entrez)
          
        }
        else {
          
          ordered.query <- TRUE
          
          #For GSEA,an ordered gene list of ALL genes must be prepared:
          #### feature 1: numeric vector
          geneList <- de_table[[fc.col]][!is.na(de_table[[entrez.col]])]
          
          ## feature 2: named vector
          names(geneList) <- as.character(de_table[[entrez.col]][!is.na(de_table[[entrez.col]])])
        
          ## feature 3: decreasing order
          geneList <- sort(geneList, decreasing = TRUE)
         
          genes <- geneList
      
        }
        
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
        method.dir <- dir(out.dir, pattern = paste0("^",method), full.names = TRUE)
        if (is.data.frame(result) & dim(result)[1] > 0) {
          
         
          gene.col <- colnames(result)[grepl("^DEGs", colnames(result))] #TODO: invent more robust approach to this?
          #First, fetching corresponding ensembles for entrez-IDs
          for (i in 1:length(result[[gene.col]])) {
            #print(result[[gene.col]])
            genes <- unlist(strsplit(result[[gene.col]][i], split = ','))
            #print(head(genes))
            ensembls <- rownames(dat)[dat[[entrez.col]] %in% genes]
            new_names <- paste(ensembls,collapse = ",") #writing with comma-separator
            #print(head(new_names))
            result[[gene.col]][i] <- new_names
          }
          #Then, converting the ensembles into symbolic names
          organism.term <- as.character(enrich.resource.terms[species, "medseqr"])
          result <- format_ensembl_ids_annotated_to_term(result, organism.term)
          
          filename <- paste0(analysis.name, ".", contrast,".", method, ".", 
                             clusterProfilerKEGG.params$analysis.approach, ".", 
                             clusterProfilerKEGG.params$ontology,".xls" )
          full.filename <- file.path(method.dir, filename)
          cat("      Saving ", method, " result table into ",  full.filename, "...\n")
          WriteXLS(result, ExcelFileName = full.filename, SheetNames = NULL, BoldHeaderRow = T)
          
          #Making the graph visualisation
          graph.name = gsub(".xls", ".network", full.filename, fixed = TRUE)
          cat("      Saving ", method, " network with the name ", graph.name, ".pdf", "...\n")
          plot_enrichment_network(enrichment.result = result, DE.result = de_table, 
                                  plot.filename = graph.name, show.terms = 5, logfc.thr = 1, fdr.thr = 0.05)
        } 
        
        enrichment_out.l$clusterProfiler_KEGG[[contrast]] <- result  
      }
      
      if (method == "DAVID") {
        
        default.params <- list(email.address = "", url = "",time.out.value = 60000, 
                              annotation.category = "GOTERM_BP_FAT", max.gene.set.size = 1000)
        missing.params = default.params[names(default.params)[!(names(default.params) %in% names(david.params))]]
        david.params = c(david.params, missing.params)
        
        cat("   Performing DAVID... \n")
        if (david.params$email.address != "") {
          
          result <- doDavidEnrichmentAnalysis(background.ensembl.ids = background.genes,
                                   foreground.ensembl.ids = input.genes, 
                                   email.address = david.params$email.address,
                                   url.address = david.params$url,
                                   time.out.value = david.params$time.out.value, 
                                   annotation.category = david.params$annotation.category, 
                                   pval.thr = p.thr, 
                                   max.gene.set.size = david.params$max.gene.set.size)
          
          
        
        }
        else{
          cat("      No e-mail address. No connection into DAVID...\n")
          result <- "No e-mail address. DAVID could not be performed."
        }
       
        # Formattig and saving the table, if relevant
        method.dir <- dir(out.dir, pattern = paste0("^",method), full.names = TRUE) # detecting output directory
        if (is.data.frame(result)) {
          
          organism.term <- as.character(enrich.resource.terms[species, "medseqr"])
          result <- format_ensembl_ids_annotated_to_term(result, organism.term)
          
          filename <- paste0(analysis.name, ".", contrast,".", method, ".", david.params$annotation.category,".xls" )
          full.filename <- file.path(method.dir, filename)
          cat("      Saving ", method, " result table into ",  full.filename, "...\n")
          WriteXLS(result, ExcelFileName = full.filename, SheetNames = NULL, BoldHeaderRow = T)
          
          #Making the graph visualisation
          graph.name = gsub(".xls", ".network", full.filename, fixed = TRUE)
          cat("      Saving ", method, " network with the name ", graph.name, ".pdf", "...\n")
          plot_enrichment_network(enrichment.result = result, DE.result = de_table, 
                                  plot.filename = graph.name, show.terms = 5, logfc.thr = 1, fdr.thr = 0.05)
        }
        
        enrichment_out.l$DAVID[[contrast]] = result
      }

      if (method == "gProfileR") {
        
        cat("   Performing GO BP enrichment with gProfileR... \n")
        org <- as.character(enrich.resource.terms[species, method])
        print(paste0("Organism: ",org))
        print(paste0("Data sources: ",gProfiler.params$data.sources))
        
        result <- run_gprofiler(input.genes, background.genes, 
                                organism = org,
                                data_sources = gProfiler.params$data.sources,
                                show_only_significant = gProfiler.params$show.only.significant
                                ) 
        
        # Formattig and saving the table, if relevant
        method.dir <- dir(out.dir, pattern = paste0("^",method), full.names = TRUE) # detecting output directory
        if (is.data.frame(result) & dim(result)[1] > 0) {
          
          organism.term <- as.character(enrich.resource.terms[species, "medseqr"])
          result <- format_ensembl_ids_annotated_to_term(result, organism.term)
          
          filename <- paste0(analysis.name, ".", contrast,".", method, ".", gProfiler.params$data.sources,".xls" )
          full.filename <- file.path(method.dir, filename)
          cat("      Saving ", method, " result table into ",  full.filename, "...\n")
          WriteXLS(result, ExcelFileName = full.filename, SheetNames = NULL, BoldHeaderRow = T)
          
          #Making the graph visualisation
          graph.name = gsub(".xls", ".network", full.filename, fixed = TRUE)
          cat("      Saving ", method, " network with the name ", graph.name, ".pdf", "...\n")
          plot_enrichment_network(enrichment.result = result, DE.result = de_table, 
                                  plot.filename = graph.name, show.terms = 5, logfc.thr = 1, fdr.thr = 0.05)
        }
        else {
          result <- "No significant enrichment found"
        }
    
        enrichment_out.l$gProfileR[[contrast]] = result
  
      }
      
      if (method == "topGO") {
        
        method.dir <- dir(out.dir, pattern = paste0("^",method), full.names = TRUE)
        
        cat("Performing topGO.... \n")
        print(paste0("Ontologies used: ",topGO.params$ontologies.used))
        print(paste0("Organism: ",topGO.params$org))
        
        org.db <- as.character(enrich.resource.terms[species, method])
        result <- run.topGO(background = background.genes, foreground = input.genes,ontologies =  topGO.params$ontologies.used, organism = org.db)
      
        #Formatting and saving the table, if relevant
        if (is.data.frame(result) & dim(result)[1] > 0) {
          
          organism.term <- as.character(enrich.resource.terms[species, "medseqr"])
          result <- format_ensembl_ids_annotated_to_term(result, organism.term)
          
          filename <- paste0(analysis.name, ".", contrast,".", method, ".", topGO.params$ontologies.used,".xls" )
          full.filename <- file.path(method.dir, filename)
          cat("      Saving ", method, " result table into ",  full.filename, "...\n")
          WriteXLS(result, ExcelFileName = full.filename, SheetNames = NULL, BoldHeaderRow = T)
          
          #Making the graph visualisation
          graph.name = gsub(".xls", ".network", full.filename, fixed = TRUE)
          cat("      Saving ", method, " network with the name ", graph.name, ".pdf", "...\n")
          plot_enrichment_network(enrichment.result = result, DE.result = de_table, 
                                  plot.filename = graph.name, show.terms = 5, logfc.thr = 1, fdr.thr = 0.05)
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



