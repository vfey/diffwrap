runEnrichmentAnalyses <- function(diffr_wrapper.output, analysis.name="", use.background.from.diffr.output=TRUE, 
                                  enrichment.methods=c("clusterProfiler", "DAVID"), organism.db.for.clusterProfiler=org.Hs.eg.db,
                                  p.thr=0.05, fdr.thr=0.05, logfc.thr=1, out.dir=NULL, david.email.address="", david.url="", max.gene.set.size=1000,
                                  do.similarity.filtering=FALSE)
{
  
  enrichment_out.l <- list()
  enrichment_out.l$clusterProfiler_ORA <- list()
  enrichment_out.l$clusterProfiler_GSEA <- list()
  enrichment_out.l$DAVID <- list()
  
  contrast.names = names(diffr_wrapper.output$contrasts)
  cat("Performing enrichment analyses for ", length(contrast.names), " comparisons: \n")
  cat("  ", contrast.names, " \n")
  
  dat <- diffr_wrapper.output$contrasts[[1]] #extracting appropriate colnames by examining the first contrast
  pv.col <- names(dat)[grep("^p\\.{0,1}val[e-u]{0,2}$", tolower(names(dat)))]
  fdr.col <- names(dat)[grep("^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", tolower(names(dat)))]
  fc.col <- names(dat)[grep("^logfc$|fold$", tolower(names(dat)))]
  
  background_genes = ""
  if(use.background.from.diffr.output) {
    # background.genes <- rownames(diffr_wrapper.output$d2$counts) ## THIS should be used after filtering is corrected. Until that, two extra elements need to be skipped:
    background.genes <- 
      rownames(diffr_wrapper.output$d2$counts)[rownames(diffr_wrapper.output$d2$counts) != "N_multimapping" |
                                                 rownames(diffr_wrapper.output$d2$counts) != "N_noFeature"]
    cat("The genes of full expression table (", length(background.genes), ") is used as the background... \n") 
  }
  
  #for (contrast in contrast.names) {
  contrast=contrast.names[2]
    cat("***", contrast, "*** \n")
    
    de_table = diffr_wrapper.output$contrasts[[contrast]]
    filtered_de_table = de_table[de_table[[fdr.col]] < fdr.thr & de_table[[fc.col]] >= logfc.thr,]
    input.genes = rownames(filtered_de_table)
    cat("  ", length(input.genes), " genes used as input \n")
    
    for (method in enrichment.methods) { 
   
     if(method == "clusterProfiler"){
       enrichment_out.l$clusterProfiler_ORA[[contrast]] = run_clusterProfiler_GO(input_genes = input.genes,
                                                                                 background_genes = background.genes,
                                                                                 file_name = NULL,
                                                                                 ordered_query = FALSE,
                                                                                 id_type = "ENSEMBL",
                                                                                 ontology = "BP",
                                                                                 OrgDb = organism.db.for.clusterProfiler,
                                                                                 pvalueCutoff = p.thr,
                                                                                 min_set_size = 10,
                                                                                 max_set_size = max.gene.set.size,
                                                                                 min_overlap = 2,
                                                                                 pAdjustMethod = "BH",
                                                                                 similarity_filtering = do.similarity.filtering)

       #For GSEA,an ordered gene list must be prepared:
       
       #### feature 1: numeric vector 
       geneList = filtered_de_table[[fc.col]]
       ## feature 2: named vector
       names(geneList) = as.character(rownames(filtered_de_table))
       ## feature 3: decreasing order
       geneList = sort(geneList, decreasing = TRUE)
       
       enrichment_out.l$clusterProfiler_GSEA[[contrast]] = run_clusterProfiler_GO(input_genes = geneList, 
                                                                                  background_genes = background.genes, 
                                                                                  file_name = NULL, 
                                                                                  ordered_query = TRUE, 
                                                                                  id_type = "ENSEMBL", 
                                                                                  ontology = "BP", 
                                                                                  OrgDb = organism.db.for.clusterProfiler, 
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
          
          enrichment_out.l$DAVID[[contrast]] = doDavidEnrichmentAnalysis(background.ensembl.ids = background.genes,
                                   foreground.ensembl.ids = input.genes, email.address=david.email.address,
                                   url.address = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/",
                                   time.out.value = 60000, annotation.category = "GOTERM_BP_FAT", pval.thr = p.thr, max.gene.set.size = max.gene.set.size, save.result.table = TRUE, out.dir=out.dir,analysis.name,
                                   contrast.name=contrast)
        
        }
        else{
          
        }
      }
      
    }
    
  #}
  return(enrichment_out.l)
}
  
