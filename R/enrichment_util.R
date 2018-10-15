# Wrapper function for performing various enrichment analyses
# 
# Author: Meeri Pekkarinen
###############################################################################
runEnrichmentAnalyses <- function(diffr_wrapper.output, analysis.name="", 
                                  use.background.from.diffr.output=TRUE, out.dir=NULL,
                                  species="human",
                                  p.thr=0.05, fdr.thr=0.05, logfc.thr=1, 
                                  enrichment.methods = c("clusterProfilerGO", "DAVID", "gProfileR", "topGO"), 
                                  david.email.address="", david.url="", max.gene.set.size = 1000, #params for David
                                  organism.db = "org.Hs.eg.db", do.similarity.filtering=FALSE, # params for clusterProfiler, organism.db used in topGO too
                                  ontologies.used = c("BP")) #params for topGO
{
  
  ## TODO: invent a smarter way to do this...
  enrich.resource.terms = data.frame("Organism" = c("human", "mouse"), 
                                     "clusterProfilerGO" = c("org.Hs.eg.db", "org.Mm.eg,db"),
                                     "clusterProfilerKEGG" = c("hsa","mmu"),
                                     "gProfileR" = c("hsapiens", "mmusculus"))
  rownames(enrich.resource.terms) = c("human", "mouse")
  
  enrichment_out.l <- list()
  enrichment_out.l$clusterProfiler_ORA <- list()
  enrichment_out.l$clusterProfiler_GSEA <- list()
  enrichment_out.l$DAVID <- list()
  enrichment_out.l$gProfileR <- list()
  enrichment_out.l$topGO <- list()
  
  contrast.names = names(diffr_wrapper.output$contrasts)
  cat("Performing enrichment analyses for ", length(contrast.names), " comparisons: \n")
  cat("  ", contrast.names, " \n")
  
  dat <- diffr_wrapper.output$contrasts[[1]] #extracting appropriate colnames by examining the first contrast
  pv.col <- names(dat)[grep("^p\\.{0,1}val[e-u]{0,2}$", tolower(names(dat)))]
  fdr.col <- names(dat)[grep("^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", tolower(names(dat)))]
  fc.col <- names(dat)[grep("^logfc$|fold$", tolower(names(dat)))]
  
  background.genes = ""
  if(use.background.from.diffr.output) {
    # background.genes <- rownames(dat) ## THIS should be used after filtering is corrected. Until that, two extra elements need to be skipped:
    background.genes <- 
      rownames(dat)[rownames(dat) != "N_multimapping" | rownames(dat) != "N_noFeature"]
    cat("The genes of full expression table (", length(background.genes), ") is used as the background... \n") 
  }
 
  for (contrast in contrast.names) {

    cat("***", contrast, "*** \n")
    de_table = diffr_wrapper.output$contrasts[[contrast]]
    
    filtered_de_table = de_table[de_table[[fdr.col]] < fdr.thr & de_table[[fc.col]] >= logfc.thr,]
    
    #if there are no entries with significant fdr, then filter by p.value
    if (!(nrow(filtered_de_table) > 0)) {
      filtered_de_table = de_table[de_table[[pv.col]] < p.thr ,]
    }
    
    input.genes = rownames(filtered_de_table)
    cat("  ", length(input.genes), " genes used as input \n")
    
    for (method in enrichment.methods) { 
   
     if(method == "clusterProfilerGO"){
       cat("   Performing GO BP enrichment with clusterProfiler... \n")
       org.db = as.character(enrich.resource.terms[species, method])
       fname.no.suffix = file.path(out.dir, paste(analysis.name, contrast, method, "ORA", sep = "_"))

       enrichment_out.l$clusterProfiler_ORA[[contrast]] = run_clusterProfiler_GO(input_genes = input.genes,
                                                                                 background_genes = background.genes,
                                                                                 file_name = fname.no.suffix,
                                                                                 ordered_query = FALSE,
                                                                                 ontology = ontologies.used,
                                                                                 OrgDb = org.db,
                                                                                 pvalueCutoff = p.thr,
                                                                                 max_set_size = max.gene.set.size,
                                                                                 pAdjustMethod = "BH",
                                                                                 similarity_filtering = do.similarity.filtering)

       #For GSEA,an ordered gene list must be prepared:

       #### feature 1: numeric vector
       geneList <- filtered_de_table[[fc.col]]
       ## feature 2: named vector
       names(geneList) <- as.character(rownames(filtered_de_table))
       ## feature 3: decreasing order
       geneList <- sort(geneList, decreasing = TRUE)

       
       fname.no.suffix = file.path(out.dir, paste(analysis.name, contrast, method, "GSEA", sep = "_"))
       enrichment_out.l$clusterProfiler_GSEA[[contrast]] <- run_clusterProfiler_GO(input_genes = geneList,
                                                                                  background_genes = background.genes,
                                                                                  file_name = fname.no.suffix,
                                                                                  ordered_query = TRUE,
                                                                                  ontology = ontologies.used,
                                                                                  OrgDb = org.db,
                                                                                  pvalueCutoff = p.thr,
                                                                                  max_set_size = max.gene.set.size,
                                                                                  pAdjustMethod = "BH",
                                                                                  similarity_filtering = do.similarity.filtering)
      } 
    
      if(method == "DAVID") {
        cat("   Performing DAVID... \n")
        if(david.email.address != "") {
          
          enrichment_out.l$DAVID[[contrast]] <- doDavidEnrichmentAnalysis(background.ensembl.ids = background.genes,
                                   foreground.ensembl.ids = input.genes, email.address=david.email.address,
                                   url.address = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/",
                                   time.out.value = 60000, annotation.category = "GOTERM_BP_FAT", pval.thr = p.thr, max.gene.set.size = max.gene.set.size, save.result.table = TRUE, out.dir=out.dir,analysis.name,
                                   contrast.name=contrast)
        
        }
        else{
          cat("      No e-mail address. No connection into DAVID...\n")
          enrichment_out.l$DAVID[[contrast]]="No e-mail address. DAVID could not be performed."
        }
      }

      if(method == "gProfileR") {
        cat("   Performing GO BP enrichment with gProfileR... \n")
        org = as.character(enrich.resource.terms[species, method])
        fname.no.suffix = file.path(out.dir, paste(analysis.name, contrast, method, sep = "_"))
        enrichment_out.l$gProfileR[[contrast]] <- run_gprofiler(input.genes, background.genes, 
                                                                file_name = fname.no.suffix, 
                                                                data_sources = "GO:BP", 
                                                                organism = org) 
                                                        
      }
      
      if (method == "topGO") {
   
        enrichment_out.l$topGO[[contrast]] <- run.topG0(background.genes, input.genes,ontologies = ontologies.used, organism = organism.db)
      }
    }

  }
  return(enrichment_out.l)
}
