# Wrapper function for performing various enrichment analyses
# 
# Author: Meeri Pekkarinen

#' Wrapper for executing various enrichment analyses
#' @description \command{runEnrichmentAnalyses} enables the auto-run of some over-representation analysis (ORA) and 
#' gene set enrichment analysis (GSEA) for the output of diffExpr main wrapper. 
#' Functions in various R-packages (clusterProfiler, topGO, gProfileR and RDAVIDWebService) 
#' are integrated. Currently supports human or mouse!
#' @param diffr.wrapper.output \code{list}.  Nested list system produced by diffrExpr-wrapper. 
#' Has to contain element "contrasts" that contains contrast-specific expression tables
#' @param analysis.name \code{character}. Descriptive character-tag used in output file names
#' @param use.background.from.diffr.output \code{logical}. Whether the (ORA) analyses are run 
#' with experiment-specific background obtained from pre-filtered expression matrix or with the default background of the functions (genome)
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
#' @details DAVID approach requires a registered email-address, correct java-version and an url configured with the settings
#' @return A list of all relevant objects generated in the course of the enrichment analyses
###############################################################################
runEnrichmentAnalyses <- function(diffr.wrapper.output, analysis.name="", 
                                  use.background.from.diffr.output=TRUE, 
                                  out.dir=NULL,
                                  species="human",
                                  p.thr=0.05, fdr.thr=0.05, logfc.thr=1, 
                                  enrichment.methods = c("clusterProfilerGO", "clusterProfilerKEGG","DAVID", "gProfileR", "topGO"), 
                                  clusterProfilerGO.params = list(analysis.approach="ORA",
                                                                  do.similarity.filtering=F,
                                                                  min.gene.set.size=10,
                                                                  max.gene.set.size=1000, 
                                                                  ontology="BP", 
                                                                  min.overlap=2,
                                                                  p.adjust.method="BH"
                                                                  ),
                                  clusterProfilerKEGG.params = list(analysis.approach="ORA",
                                                                    min.gene.set.size=10,
                                                                    max.gene.set.size=1000, 
                                                                    ontology="BP", 
                                                                    min.overlap=2,
                                                                    p.adjust.method="BH"
                                                                    ),
                                  david.params = list(email.address="", 
                                                      url="",
                                                      time.out.value = 60000,
                                                      annotation.category = "GOTERM_BP_FAT",
                                                      max.gene.set.size = 1000),
                                  gProfiler.params = list(data.sources="GO:BP", show.only.significant = TRUE),
                                  topGO.params = list(ontologies.used = c("BP"), org = "org.Hs.eg.db")
                                  ) 
{
  
  ## TODO: invent a smarter way to do this...
  enrich.resource.terms = data.frame("Organism" = c("human", "mouse"), 
                                     "clusterProfilerGO" = c("org.Hs.eg.db", "org.Mm.eg,db"),
                                     "clusterProfilerKEGG" = c("hsa","mmu"),
                                     "gProfileR" = c("hsapiens", "mmusculus"),
                                     "topGO" = c("org.Hs.eg.db", "org.Mm.eg,db"))
  rownames(enrich.resource.terms) = c("human", "mouse")
  
  enrichment_out.l <- list()
  enrichment_out.l$clusterProfiler_GO <- list()
  enrichment_out.l$clusterProfiler_KEGG <- list()
  enrichment_out.l$DAVID <- list()
  enrichment_out.l$gProfileR <- list()
  enrichment_out.l$topGO <- list()
  
  #Creating subfolders
  lapply(enrichment.methods, function(x) dir.create(file.path(out.dir, method), showWarnings = F))
  
  contrast.names = names(diffr.wrapper.output$contrasts)
  cat("Performing enrichment analyses for ", length(contrast.names), " comparisons: \n")
  cat("  ", contrast.names, " \n")
  
  dat <- diffr.wrapper.output$contrasts[[1]] ## extracting appropriate colnames using the first contrast
  pv.col <- names(dat)[grep("^p\\.{0,1}val[e-u]{0,2}$", tolower(names(dat)))]
  fdr.col <- names(dat)[grep("^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", tolower(names(dat)))]
  fc.col <- names(dat)[grep("^logfc$|fold$", tolower(names(dat)))]
  entrez.col <- names(dat)[grep("entrez", tolower(names(dat)))]
  
  
  background.genes = ""
  if(use.background.from.diffr.output) {
    # background.genes <- rownames(dat) ## THIS should be used after filtering is corrected. Until that, two extra elements need to be skipped:
    background.genes <- 
      rownames(dat)[rownames(dat) != "N_multimapping" | rownames(dat) != "N_noFeature"]
    cat("The genes of full expression table (", length(background.genes), ") is used as the background (expect in GSEA-analyses, where no background is used)... \n") 
  }
 
  for (contrast in contrast.names) {

    cat("***", contrast, "*** \n")
    de_table = diffr.wrapper.output$contrasts[[contrast]]
    
    cat("   Filtering the DE genes (fdr < ", fdr.thr, "and abs. logFC >=", logfc.thr, ")...\n")
    filtered_de_table = de_table[de_table[[fdr.col]] < fdr.thr & abs(de_table[[fc.col]]) >= logfc.thr,]
    
    #if there are no entries with significant fdr, then filter by p.value
    if (!(nrow(filtered_de_table) > 0)) {
      cat("      No entries with significant fdr. P-values used instead...\n")
      filtered_de_table = de_table[de_table[[pv.col]] < p.thr & abs(de_table[[fc.col]]) >= logfc.thr,]
    }
    
    input.genes = rownames(filtered_de_table)
    cat("  ", length(input.genes), " genes used as input (expect in GSEA-analyses, where all measured genes is used) \n")
    
    for (method in enrichment.methods) { 
   
     if(method == "clusterProfilerGO"){
       
       ## Initial preparations: setting result directory and replacing missing parameters
       method.dir = dir(out.dir, pattern = paste0("^",method), full.names = TRUE)
       default.params = list(analysis.approach="ORA", do.similarity.filtering=F, min.gene.set.size=10, 
                             max.gene.set.size=1000, ontology="BP", min.overlap=2, p.adjust.method="BH")
       missing.params = default.params[names(default.params)[!(names(default.params) %in% names(clusterProfilerGO.params))]]
       clusterProfilerGO.params = c(clusterProfilerGO.params, missing.params)
       
       cat("   Performing GO BP enrichment with clusterProfiler... \n")
       cat("   ")
       print(unlist(clusterProfilerGO.params))
       org.db = as.character(enrich.resource.terms[species, method])
      
       if(clusterProfilerGO.params$analysis.approach == "ORA"){
         
         ordered.query=FALSE
         genes = input.genes
         fname.no.suffix = file.path(method.dir, paste(analysis.name, contrast, method, "ORA", sep = "_"))
         
       }
       else {
         
         ordered.query = TRUE
         
         #For GSEA,an ordered gene list of ALL genes must be prepared:
         #### feature 1: numeric vector
         geneList <- de_table[[fc.col]]
         ## feature 2: named vector
         names(geneList) <- as.character(rownames(de_table))
         ## feature 3: decreasing order
         geneList <- sort(geneList, decreasing = TRUE)
         genes = geneList
         fname.no.suffix = file.path(method.dir, paste(analysis.name, contrast, method, "GSEA", sep = "_"))
       }
       enrichment_out.l$clusterProfiler_GO[[contrast]] = run_clusterProfiler_GO(input_genes = genes,
                                                                                 background_genes = background.genes,
                                                                                 file_name = fname.no.suffix,
                                                                                 ordered_query = ordered.query, 
                                                                                 OrgDb = org.db,
                                                                                 pvalueCutoff = fdr.thr,
                                                                                 ontology = clusterProfilerGO.params$ontology,
                                                                                 min_set_size = clusterProfilerGO.params$min.gene.set.size,
                                                                                 max_set_size = clusterProfilerGO.params$max.gene.set.size,
                                                                                 min_overlap = clusterProfilerGO.params$min.overlap,
                                                                                 pAdjustMethod = clusterProfilerGO.params$p.adjust.method,
                                                                                 similarity_filtering = clusterProfilerGO.params$do.similarity.filtering)

      
     } 
      
      if(method == "clusterProfilerKEGG") {
      
        ## Initial preparations: setting result directory and replacing missing parameters
        method.dir = dir(out.dir, pattern = paste0("^",method), full.names = TRUE)
        default.params = list(analysis.approach="ORA", min.gene.set.size=10, max.gene.set.size=1000, 
                              ontology="BP", min.overlap=2, p.adjust.method="BH")
        missing.params = default.params[names(default.params)[!(names(default.params) %in% names(clusterProfilerKEGG.params))]]
        clusterProfilerKEGG.params = c(clusterProfilerKEGG.params, missing.params)
        
    
        cat("   Performing KEGG enrichment with clusterProfiler... \n")
        cat("   ")
        print(unlist(clusterProfilerKEGG.params))
        org = as.character(enrich.resource.terms[species, method])
      
        if(length(entrez.col) == 0){
          cat("      No entrez IDs found  from the data. skipping KEGG-enrichment...\n")
          break 
        }
        
        
        background.gene.entrez = ""
        if(clusterProfilerKEGG.params$analysis.approach == "ORA"){
          
          if(length(background.genes)  > 1){
            background.gene.entrez = as.character(de_table[[entrez.col]][!is.na(de_table[[entrez.col]])])
            cat("      ",length(background.gene.entrez), "/", length(background.genes), " of the all measured genes having corresponding entrez id used as background...\n")
            
          }
          else{
            cat("      Using default background in KEGG-enrichment...\n")
          }
          
          input.gene.entrez=as.character(filtered_de_table[[entrez.col]][!is.na(filtered_de_table[[entrez.col]])])
          cat("      ",length(input.gene.entrez), "/", length(input.genes), " of the input genes having corresponding entrez id used for the analysis...\n")
          
          
          ordered.query=FALSE
          genes = as.character(input.gene.entrez)
          fname.no.suffix = file.path(method.dir, paste(analysis.name, contrast, method, "ORA", sep = "_"))
          
        }
        else {
          
          ordered.query = TRUE
          
          #For GSEA,an ordered gene list of ALL genes must be prepared:
          #### feature 1: numeric vector
          geneList <- de_table[[fc.col]][!is.na(de_table[[entrez.col]])]
          
          ## feature 2: named vector
          names(geneList) <- as.character(de_table[[entrez.col]][!is.na(de_table[[entrez.col]])])
        
          ## feature 3: decreasing order
          geneList <- sort(geneList, decreasing = TRUE)
         
          genes = geneList
          fname.no.suffix = file.path(method.dir, paste(analysis.name, contrast, method, "GSEA", sep = "_"))
        }
        
        enrichment_out.l$clusterProfiler_KEGG[[contrast]] = run_clusterProfiler_KEGG(input_genes=genes,
                                              background_genes = background.gene.entrez,
                                              file_name = NULL,
                                              ordered_query = ordered.query,
                                              organism = org,
                                              pvalueCutoff = fdr.thr,
                                              min_set_size = clusterProfilerKEGG.params$min.gene.set.size,
                                              max_set_size = clusterProfilerKEGG.params$max.gene.set.size,
                                              min_overlap = clusterProfilerKEGG.params$min.overlap,
                                              pAdjustMethod = clusterProfilerKEGG.params$p.adjust.method)
         
      }
      
      if(method == "DAVID") {
        
        method.dir = dir(out.dir, pattern = paste0("^",method), full.names = TRUE)
        default.params = list(email.address="", url="",time.out.value = 60000, 
                              annotation.category = "GOTERM_BP_FAT", max.gene.set.size = 1000)
        missing.params = default.params[names(default.params)[!(names(default.params) %in% names(david.params))]]
        david.params = c(david.params, missing.params)
        
        cat("   Performing DAVID... \n")
        if(david.params$email.address != "") {
          
          enrichment_out.l$DAVID[[contrast]] <- doDavidEnrichmentAnalysis(background.ensembl.ids = background.genes,
                                   foreground.ensembl.ids = input.genes, 
                                   email.address = david.params$email.address,
                                   url.address = david.params$url,
                                   time.out.value = david.params$time.out.value, 
                                   annotation.category = david.params$annotation.category, 
                                   pval.thr = p.thr, 
                                   max.gene.set.size = david.params$max.gene.set.size, 
                                   save.result.table = TRUE, out.dir=method.dir,analysis.name,
                                   contrast.name=contrast)
        
        }
        else{
          cat("      No e-mail address. No connection into DAVID...\n")
          enrichment_out.l$DAVID[[contrast]]="No e-mail address. DAVID could not be performed."
        }
      }

      if(method == "gProfileR") {
        
        method.dir = dir(out.dir, pattern = paste0("^",method), full.names = TRUE)
        
        cat("   Performing GO BP enrichment with gProfileR... \n")
        org = as.character(enrich.resource.terms[species, method])
        print(paste0("Organism: ",org))
        print(paste0("Data sources: ",gProfiler.params$data.sources))
        fname.no.suffix = file.path(method.dir, paste(analysis.name, contrast, method, sep = "_"))
        enrichment_out.l$gProfileR[[contrast]] <- run_gprofiler(input.genes, background.genes, 
                                                                file_name = fname.no.suffix,
                                                                organism = org,
                                                                data_sources = gProfiler.params$data.sources,
                                                                show_only_significant = gProfiler.params$show.only.significant
                                                                ) 
                                                        
      }
      
      if (method == "topGO") {
        
        method.dir = dir(out.dir, pattern = paste0("^",method), full.names = TRUE)
        
        cat("Performing topGO.... \n")
        print(paste0("Ontologies used: ",topGO.params$ontologies.used))
        print(paste0("Organism: ",topGO.params$org))
        
        org.db = as.character(enrich.resource.terms[species, method])
        enrichment_out.l$topGO[[contrast]] <- run.topGO(background=background.genes, foreground = input.genes,ontologies =  topGO.params$ontologies.used, organism = org.db)
      }
    }

  }
  return(enrichment_out.l)
}
