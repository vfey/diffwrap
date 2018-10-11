runEnrichmentAnalyses <- function(diffr_wrapper.output, analysis.name="", use.background.from.diffr.output=TRUE, 
                                  enrichment.methods=c("clusterProfilerGO", "DAVID", "gProfileR"),
                                  p.thr=0.05, fdr.thr=0.05, logfc.thr=1, out.dir=NULL, david.email.address="", david.url="", max.gene.set.size=1000,
                                  do.similarity.filtering=FALSE)
{
  
  species="human" ## TODO: HOW to obtain this...
  
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
  
  #for (contrast in contrast.names) {
  contrast=contrast.names[1]
    cat("***", contrast, "*** \n")
    
    de_table = diffr_wrapper.output$contrasts[[contrast]]
    filtered_de_table = de_table[de_table[[fdr.col]] < fdr.thr & de_table[[fc.col]] >= logfc.thr,]
    input.genes = rownames(filtered_de_table)
    cat("  ", length(input.genes), " genes used as input \n")
    
    for (method in enrichment.methods) { 
   
     if(method == "clusterProfilerGO"){
       
       org.db = as.character(enrich.resource.terms[species, method])
       enrichment_out.l$clusterProfiler_ORA[[contrast]] = run_clusterProfiler_GO(input_genes = input.genes,
                                                                                 background_genes = background.genes,
                                                                                 file_name = NULL,
                                                                                 ordered_query = FALSE,
                                                                                 id_type = "ENSEMBL",
                                                                                 ontology = "BP",
                                                                                 OrgDb = org.db,
                                                                                 pvalueCutoff = p.thr,
                                                                                 min_set_size = 10,
                                                                                 max_set_size = max.gene.set.size,
                                                                                 min_overlap = 2,
                                                                                 pAdjustMethod = "BH",
                                                                                 similarity_filtering = do.similarity.filtering)

       #For GSEA,an ordered gene list must be prepared:
       
       #### feature 1: numeric vector 
       geneList <- filtered_de_table[[fc.col]]
       ## feature 2: named vector
       names(geneList) <- as.character(rownames(filtered_de_table))
       ## feature 3: decreasing order
       geneList <- sort(geneList, decreasing = TRUE)
       
       enrichment_out.l$clusterProfiler_GSEA[[contrast]] <- run_clusterProfiler_GO(input_genes = geneList, 
                                                                                  background_genes = background.genes, 
                                                                                  file_name = NULL, 
                                                                                  ordered_query = TRUE, 
                                                                                  id_type = "ENSEMBL", 
                                                                                  ontology = "BP", 
                                                                                  OrgDb = org.db, 
                                                                                  pvalueCutoff = p.thr, 
                                                                                  min_set_size = 10, 
                                                                                  max_set_size = 1000, 
                                                                                  min_overlap = 2, 
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
          
        }
      }
      if(method == "gProfileR") {
        enrichment_out.l$gProfileR[[contrast]] <- run_gprofiler(input.genes, background.genes, data_sources = "BP", file_name = paste(out.dir,"GOBP_enrichment", sep="/"))
      }
    }
    
  #}
  return(enrichment_out.l)
}
  
