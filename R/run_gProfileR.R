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
#' @param organism By default uses "hsapiens".
#' @param background_genes A character vector of background gene IDs. If not specified, by default uses all human genes
#' annotated to term domain.
#' @param file_name A character string used as a file name without file extension. If not NULL output will be saved
#' to an Excel file. Currently not in use.
#' @param ordered_query If set to TRUE, computes GSEA-style p-values for an input ranked gene list.
#' @param multi_query In case of multiple gene lists, returns comparison table of these lists.
#' If enabled, the result data frame has columns named 'p_values', 'gconvert_sizes', 'intersection_sizes' with
#' vectors showing values in the order of input queries.
#' @param show_only_significant Shows only significant results based on the set padjCutoff, default TRUE.
#' @param evidence_codes If set to TRUE, includes evidence codes to the results.
#' Note that this can decrease performance and make the query slower. In addition, a column
#' 'intersection' is created that contains the gene id-s that intersect between the query and term.
#' This parameter does not work if 'multi_query' is set to TRUE.
#' @param exclude_iea If TRUE, excludes electronic GO annotations (with evidence code IEA).
#' @param measure_underrepresentation If TRUE, measures under-representation.
#' @param max_p_value Adjusted p-value cut-off. Shows only terms with p-value under this cut-off if showOnlySignificant = TRUE.
#' @param correction_method The algorithm used for multiple testing correction,
#' one of "gSCS", "fdr", "bonferroni". By default uses "fdr".
#' @param domain_scope How to define statistical domain, one of "annotated", "known", "custom" or "custom_annotated".
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
#' @param highlight If set to TRUE, returns a TRUE-FALSE column called 'highlighted' to indicate driver terms in GO.
#' @return A table listing statistically significant enrichment results according to threshold set in padjCutoff.
#'     The table is also saved in xlsx format with user-specified name.
#' @export

run_gprofiler <- function(input_genes,
                          organism = "hsapiens",
                          background_genes = "",
                          file_name = NULL,
                          ordered_query = FALSE,
                          multi_query = FALSE,
                          show_only_significant = FALSE,
                          evidence_codes = FALSE,
                          exclude_iea = FALSE,
                          measure_underrepresentation = FALSE,
                          max_p_value = 0.5,
                          correction_method = "fdr",
                          domain_scope = "annotated",
                          data_sources = NULL,
                          highlight = FALSE)
{
  message("Set parameters:")
  message(paste0(">> ordered_query: ",ordered_query))
  message(paste0(">> multi_query: ",multi_query))
  message(paste0(">> show_only_significant: ",show_only_significant))
  message(paste0(">> evidence_codes: ",evidence_codes))
  message(paste0(">> exclude_iea: ",exclude_iea))
  message(paste0(">> measure_underrepresentation: ",measure_underrepresentation))
  message(paste0(">> max_p_value: ",max_p_value))
  message(paste0(">> correction_method: ", correction_method))
  message(paste0(">> domain_scope: ", domain_scope))
  message(paste0(">> data_sources: ",data_sources))
  message(paste0(">> highlight: ", highlight))

  results <- gprofiler2::gost(query = as.character(input_genes),
                              organism = organism,
                              ordered_query = ordered_query,
                              multi_query = multi_query,
                              significant = show_only_significant,
                              exclude_iea = exclude_iea,
                              measure_underrepresentation = measure_underrepresentation,
                              evcodes = evidence_codes,
                              user_threshold = max_p_value,
                              correction_method = correction_method,
                              domain_scope = domain_scope,
                              custom_bg = background_genes,
                              sources = data_sources,
                              highlight = highlight)

  if (!is.null(file_name)) {
    cat("         Saving result into ", paste0(file_name, ".xlsx"), "...\n")
    WriteXLS::WriteXLS(results, ExcelFileName = paste0(file_name, ".xlsx"), SheetNames = NULL, BoldHeaderRow = T)
  }
  return(results)
}
