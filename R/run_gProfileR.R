# Load packages
library(gProfileR)
library(WriteXLS)

# Function for running over-representation analysis (ORA) or gene set enrichment analysis (GSEA) using gprofiler R interface.

# Input ORA:  input_genes = vector of differentially expressed genes
#             background_genes = vector of all measured genes (optional)

# Input GSEA: input_genes = vector of genes ranked in the order of decreasing importance

# Input IDs can be Ensembl IDs or any other IDs that have been linked to genes
# in the Ensembl database, also mixed list of IDs

# Example: run_gprofiler(DEG.results.sig$ID, DEG.results.full$ID, dataSources = "GO:BP", fileName = "GOBP_enrichment").


#' Runs gprofiler function for a DEG list or for a ranked gene list.
#' @param input_genes A character vector of gene IDs, e.g. Ensembl or HGNC. Can be any ID type that has been linked
#' to genes in the Ensembl database, and also a mixed vector of IDs.
#' @param background_genes A character vector of background gene IDs. If not specified, by default uses all human genes
#' annotated to term domain.
#' @param file_name A character string used as a file name.
#' @param ordered_query If set to TRUE, computes GSEA-style p-values for an input ranked gene list.
#' @param show_only_significant Shows only significant results based on the set padjCutoff, default TRUE.
#' @param data_sources A vector of data sources to use. One or more of the following:
#'     GO:BP = GO Biological Process
#'     GO:MF = GO Molecular Function
#'     GO:CC = GO Cellular Component
#'     KEGG = Kyoto Encyclopedia of Genes and Genomes Pathway Database
#'     REAC = Reactome Pathway Database
#'     TF = TRANSFAC Database, putative transcription factor binding sites
#'     MI = miRTarBase, microRNA-Target Interactions
#'     CORUM = Database of manually annotated protein complexes in human, mouse and rat
#'     HPA = Human Protein Atlas, protein expression in normal tissues
#'     HP = Human Phenotype Ontology, human disease gene annotations
#'     OMIM = Online Mendelian Inheritance in Man, an online catalog of human genes and genetic disorders
#'     By default, uses all data sources.
#' @param organism By default uses "hsapiens".
#' @param exclude_iea If TRUE, excludes electronic GO annotations (with evidence code IEA).
#' @param max_p_value Adjusted p-value cut-off. Shows only terms with p-value under this cut-off if showOnlySignificant = TRUE.
#' @param min_set_size Minimun size of the functional category, uses 10 by default.
#' @param max_set_size Maximum size of the functional categowy, uses 1000 by default.
#' @param min_overlap Minimum size of the overlap (intersection) between query and functional category,
#'     smaller intersections are excluded.
#' @param correction_method The algorithm used for multiple testing correction,
#'     one of "gSCS", "fdr", "bonferroni". By default uses "fdr".
#' @param sort_by_structure If TRUE, hiearchical sorting of the results is enabled.
#' @param hier_filtering Hierarchical filtering strength, one of "none", "moderate", "strong", by default "none".
#'     The hierarchical filtering option is only available if the option of hierarchical sorting is activated.
#' @return A table listing statistically significant enrichment results according to threshold set in padjCutoff.
#'     The table is also saved in xlsx format with user-specified name.
#' @export

run_gprofiler <- function(input_genes,
                          background_genes = "",
                          file_name = NULL,
                          ordered_query = FALSE,
                          show_only_significant = F,
                          data_sources = NULL,
                          organism = "hsapiens",
                          exclude_iea = FALSE,
                          max_p_value = 0.5,
                          sort_by_structure = FALSE,
                          min_set_size = 10,
                          max_set_size = 1000,
                          min_overlap = 2,
                          correction_method = "fdr",
                          hier_filtering = "none") {

    results <- gprofiler(query = as.character(input_genes),
                     organism = organism,
                     sort_by_structure = sort_by_structure,
                     ordered_query = ordered_query,
                     significant = show_only_significant,
                     exclude_iea = exclude_iea,
                     max_p_value = max_p_value,
                     min_set_size = min_set_size,
                     max_set_size = max_set_size,
                     min_isect_size = min_overlap,
                     correction_method = correction_method,
                     hier_filtering = hier_filtering,
                     custom_bg = background_genes,
                     src_filter = data_sources)

    cols.res = colnames(results)
    cols.interest = c("term.id", "domain","term.name","term.size","query.size","overlap.size","p.value","intersection") 
    ids.list = list()
    for(col.intr in cols.interest) {
      ids.list[[col.intr]] = grep(col.intr, col.res)
    }
    ids.interest = as.vector(unlist(ids.list))
    
    results <- results[,ids.interest] # Extract only interesting columns
    colnames(results) <- c("Term ID", "Term Domain", "Term Description", "Term Size", "Query Size",
                       "No of DEGs annotated to Term", "Adjusted P-Value", "DEGs Annotated to Term")
    
    cat("         Saving result into ", paste0(file_name, ".xlsx"), "...\n")
    WriteXLS(results, ExcelFileName = paste0(file_name, ".xlsx"), SheetNames = NULL, BoldHeaderRow = T)
    return(results)
}
