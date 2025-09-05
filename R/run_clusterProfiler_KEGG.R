# Function for running over-representation analysis (ORA) or gene set enrichment analysis (GSEA) of
# KEGG terms using clusterProfiler.

# Input ORA:  input_genes = differentially expressed genes as a vector of Entrez IDs
#             background_genes = vector of Entrez IDs of all studied genes (optional)

# Input GSEA: input_genes = a named vector of fold changes of ALL genes ranked in the order of decreasing fold change,
#                           with Entrez IDs as names

# Example: run_clusterProfiler_KEGG(DEG.results$ID, DEG.results.full$ID, file_name = "KEGG_enrichment").


#' Runs clusterProfiler KEGG enrichment function for a DEG list or for a ranked gene list.
#' @param input_genes A character vector of Entrez gene IDs (ORA) or a named,
#'                    ordered vector of fold changes of of ALL genes with Entrez IDs as names (GSEA).
#' @param background_genes A character vector of background gene IDs. If not specified, by default uses all human genes
#'                         annotated to term domain.
#' @param file_name A character string used as a file name.
#' @param ordered_query If set to TRUE, runs gene set enrichment analysis for an input ranked gene list.
#' @param organism "hsa" by default. Supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'.
#' @param id_type By default uses "kegg", which is Entrez ID for eukaryotes and Locus ID for prokaryotes. Other options:
#'                'ncbi-geneid’, ‘ncbi-proteinid’ or ‘uniprot’.
#' @param pvalueCutoff Adjusted p-value cut-off.
#' @param min_set_size Minimun size of the functional category, uses 10 by default.
#' @param max_set_size Maximum size of the functional category, uses 1000 by default.
#' @param min_overlap Minimum size of the overlap (intersection) between query and functional category,
#'                    smaller intersections are excluded. By default uses 2.
#' @param pAdjustMethod The algorithm used for multiple testing correction, one of
#'                      "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". By default uses "BH".
#' @return A table listing statistically significant enrichment results according to threshold set in pvalueCutoff
#'         The table is also saved in xlsx format with user-specified name.
#' @export


run_clusterProfiler_KEGG <- function(input_genes,
                                     background_genes = "",
                                     file_name = NULL,
                                     ordered_query = FALSE,
                                     id_type = "kegg",
                                     organism = "hsa",
                                     pvalueCutoff = 0.05,
                                     min_set_size = 10,
                                     max_set_size = 1000,
                                     min_overlap = 2,
                                     pAdjustMethod = "BH") {


  if (ordered_query) {
    print("Running gene set enrichment analysis...")
    set.seed(123)
    GSE.results <- try(clusterProfiler::gseKEGG(geneList = input_genes,
                                            organism = organism,
                                            keyType = id_type,
                                            nPerm = 1000,
                                            minGSSize = min_set_size,
                                            maxGSSize = max_set_size,
                                            pvalueCutoff = pvalueCutoff,
                                            pAdjustMethod = "BH",
                                            seed = T))
    if (is(GSE.results, "try-error")) {
      cat(" ~~ NOTE: KEGG GSE analysis failed ~~\n")
      return(NULL)
    }

    GSE.results.df <- data.frame(GSE.results)
    GSE.results.df$Count <- sapply(GSE.results.df$core_enrichment, function(x) length(unlist(strsplit(x, split = "/"))))
    # Filter by minimun overlap and extract only interesting columns
    GSE.results.df <- GSE.results.df[GSE.results.df$Count >= min_overlap, c(1:3,5:7,11:12)]
    colnames(GSE.results.df) <- c("Term ID", "Term description", "Gene Set Size", "Normalized Enrichment Score", "P-Value",
                                  "Adjusted P-Value", "DEGs Contributing to Enrichment", "No of DEGs Contributing to Enrichment")

    # switching default gene separator of clusterProfiler into comma
    GSE.results.df$`DEGs Contributing to Enrichment` = gsub('/', ',', GSE.results.df$`DEGs Contributing to Enrichment`)

    # Write to Excel file (in comment since wrapper does this)
    #WriteXLS(GSE.results.df, ExcelFileName = paste0(file_name, ".xlsx"), SheetNames = NULL, BoldHeaderRow = T)

    return(GSE.results.df)
  }

  if (!ordered_query) {

    print("Running over-representation analysis...")
    ORA.results <- try(clusterProfiler::enrichKEGG(gene = input_genes,
                                              universe = background_genes,
                                              organism = organism,
                                              keyType = id_type,
                                              pAdjustMethod = pAdjustMethod,
                                              pvalueCutoff  = pvalueCutoff,
                                              qvalueCutoff = 0.2,
                                              minGSSize = min_set_size,
                                              maxGSSize = max_set_size))
    if (is(ORA.results, "try-error")) {
      cat(" ~~ NOTE: KEGG ORA failed ~~\n")
      return(NULL)
    }

    ORA.results.df <- data.frame(ORA.results)
    # Filter by minimun overlap and extract only interesting columns
    sel.cols <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "geneID", "Count")
    names(sel.cols) <- c("Term ID", "Term description", "Gene Ratio", "Background Ratio", "P-Value",
                         "Adjusted P-Value", "DEGs Annotated to Term", "No of DEGs Annotated to Term")
    sel.cols <- sel.cols[sel.cols %in% names(ORA.results.df)]
    ORA.results.df <- ORA.results.df[ORA.results.df$Count >= min_overlap, sel.cols]
    colnames(ORA.results.df) <- names(sel.cols)

    # switching default gene separator of clusterProfiler into comma:
    ORA.results.df$`DEGs Annotated to Term` = gsub('/', ',', ORA.results.df$`DEGs Annotated to Term`)


    # Write to Excel file (in comment since wrapper does this)
    #WriteXLS(ORA.results.df, ExcelFileName = paste0(file_name, ".xlsx"), SheetNames = NULL, BoldHeaderRow = T)

    return(ORA.results.df)
  }
}
