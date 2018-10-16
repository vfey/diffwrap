# Load packages
library(clusterProfiler)
library(WriteXLS)

# Function for running over-representation analysis (ORA) or gene set enrichment analysis (GSEA) of GO terms using clusterProfiler.

# Input ORA:  input_genes = differentially expressed genes as a vector of Ensembl IDs
#             background_genes = vector of Ensembl IDs of all studied genes (optional)

# Input GSEA: input_genes = a named vector of fold changes of ALL genes ranked in the order of decreasing fold change, 
#                           with Ensembl IDs as names

# Example: run_clusterProfiler_GO(DEG.results$ID, DEG.results.full$ID, file_name = "GOBP_enrichment").


#' Runs clusterProfiler GO enrichment function for a DEG list or for a ranked gene list.
#' @param input_genes A character vector of gene IDs (ORA) or a named, 
#'                    ordered vector of fold changes of ALL genes with gene IDs as names (GSEA).
#' @param background_genes A character vector of background gene IDs. If not specified, by default uses all human genes
#'                         annotated to term domain.
#' @param file_name A character string used as a file name.
#' @param ordered_query If set to TRUE, computes GSEA-style p-values for an input ranked gene list.
#' @param ontology A vector of ontology types to use. One of the following:
#'                 GO:BP = GO Biological Process
#'                 GO:MF = GO Molecular Function
#'                 GO:CC = GO Cellular Component
#' @param OrgDb Organism annotation package, by default uses "org.Hs.eg.db".
#' @param id_type By default uses "ENSEMBL". Can be any ID type that is supported by the corresponding OrgDb.
#' @param pvalueCutoff Adjusted p-value cut-off.
#' @param min_set_size Minimun size of the functional category, uses 10 by default.
#' @param max_set_size Maximum size of the functional categowy, uses 1000 by default.
#' @param min_overlap Minimum size of the overlap (intersection) between query and functional category, 
#'                    smaller intersections are excluded. By default uses 2.
#' @param pAdjustMethod The algorithm used for multiple testing correction, one of
#'                      "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". By default uses "BH".
#' @param similarity_filtering Similarity filtering method, either FALSE (default) or TRUE (uses simplify function).
#' @return A table listing statistically significant enrichment results according to threshold set in pvalueCutoff 
#'        The table is also saved in xlsx format with user-specified name.
#' @export



run_clusterProfiler_GO <- function(input_genes, 
                                   background_genes = "", 
                                   file_name = NULL, 
                                   ordered_query = FALSE, 
                                   id_type = "ENSEMBL", 
                                   ontology = "BP", 
                                   OrgDb = org.Hs.eg.db, 
                                   pvalueCutoff = 0.05, 
                                   min_set_size = 10, 
                                   max_set_size = 1000, 
                                   min_overlap = 2, 
                                   pAdjustMethod = "BH", 
                                   similarity_filtering = FALSE) {
  
  
  if (ordered_query) {
    print("Running gene set enrichment analysis...")
    set.seed(123)
    GSE.results <- clusterProfiler::gseGO(geneList = input_genes,
                                          OrgDb = OrgDb,
                                          ont = ontology,
                                          keyType = id_type,
                                          nPerm = 1000,
                                          minGSSize = min_set_size,
                                          maxGSSize = max_set_size,
                                          pvalueCutoff = pvalueCutoff,
                                          pAdjustMethod = "BH",
                                          seed = T)

    if(dim(GSE.results)[1] >= 1) {
    
      GSE.results.df <- data.frame(GSE.results)
      GSE.results.df$Count <- sapply(GSE.results.df$core_enrichment, function(x) length(unlist(strsplit(x, split = "/"))))
      # Filter by minimun overlap and extract only interesting columns
      GSE.results.df <- GSE.results.df[GSE.results.df$Count >= min_overlap, c(1:3,5:7,11:12)] 
      colnames(GSE.results.df) <- c("Term ID", "Term description", "Gene Set Size", "Normalized Enrichment Score", "P-Value", 
                                  "Adjusted P-Value", "DEGs Contributing to Enrichment", "No of DEGs Contributing to Enrichment")
    
      # Write to Excel file
      cat("         Saving result into ", paste0(file_name, ".xlsx"), "...\n")
      WriteXLS(GSE.results.df, ExcelFileName = paste0(file_name, ".xlsx"), SheetNames = NULL, BoldHeaderRow = T)
    
      return(GSE.results.df)
      
    } else {
      return("No enrichments found")
      
    }
  }
  
  if (!ordered_query) {
    
    print("Running over-representation analysis...")
    ORA.results <- clusterProfiler::enrichGO(gene = input_genes,
                                            universe = background_genes,
                                            OrgDb = OrgDb,
                                            keyType = id_type,
                                            ont = ontology,
                                            pAdjustMethod = pAdjustMethod,
                                            pvalueCutoff  = pvalueCutoff,
                                            qvalueCutoff = 0.2,
                                            minGSSize = min_set_size,
                                            maxGSSize = max_set_size,
                                            readable = TRUE)
    if(dim(ORA.results)[1] >= 1) {
    
      if (similarity_filtering) {
        print("Simplifying results...")
        ORA.results <- clusterProfiler::simplify(ORA.results)
      }

      ORA.results.df <- data.frame(ORA.results)
      # Filter by minimun overlap and extract only interesting columns
      ORA.results.df <- ORA.results.df[ORA.results.df$Count >= min_overlap, c(1:6,8:9)]
      colnames(ORA.results.df) <- c("Term ID", "Term description", "Gene Ratio", "Background Ratio", "P-Value", 
                                  "Adjusted P-Value", "DEGs Annotated to Term", "No of DEGs Annotated to Term")
      
      
      # Write to Excel file  
      cat("         Saving result into ", paste0(file_name, ".xlsx"), "...\n")
      WriteXLS(ORA.results.df, ExcelFileName = paste0(file_name, ".xlsx"), SheetNames = NULL, BoldHeaderRow = T)
      return(ORA.results.df)
      }
      else {
        return("No enrichments found")
      }
  }
}



  
    
    

