# Function: run.topGO()
#
# Author: Bogdan Iancu - Genevia Technologies Oy
#
# Arguments:
#         background = path to the background set of genes
#         foreground = path to the foreground set of genes
#         ontologies = character string specifying the ontology of interest (BP,MF,CC), default = "BP"
#         organism = the organism database used, eg. for Human: "org.Hs.eg.db"
#

# Output: table with enriched terms and p-values from Fisher's exact test, using different algorithms: elim, classic, weight01.
#
#
# Details: columns in the final table produced by topGo (their description is fetched from topGo documentation):
#                                    - Annotated : number of genes in org.Hs.eg.db which are annotated with the GO-term.
#                                    - Significant : number of genes belonging to your input which are annotated with the GO-term.
#                                    - Expected : show an estimate of the number of genes a node of size Annotated would have if the significant genes were to be randomly selected from the gene universe.
#                                    - pvalues : pvalue obtained after the test

#           column: p.adj.weight01 - represents the adjusted p-value for the weight01 algorithm
#           requires packages: topGO, readxl, org.Hs.eg.db

utils::globalVariables("genesInTerm")

#' Function to run GO term enrichment analysis using the 'topGO' package.
#' @param background path to the background set of genes
#' @param foreground path to the foreground set of genes
#' @param ontologies character string specifying the ontology of interest (BP,MF,CC), default = "BP"
#' @param organism the organism database used, eg. for Human: "org.Hs.eg.db"
#' @param ID_type character; the type of gene ID used in the input data, default = "ENSEMBL"
#' @param pAdjustMethod character; the method used for p-value adjustment, default = "BH"
#' @details
#' The columns in the final table produced by topGO (their description is fetched from topGo documentation):
#'   - Annotated : number of genes in org.Hs.eg.db which are annotated with the GO-term.
#'   - Significant : number of genes belonging to your input which are annotated with the GO-term.
#'   - Expected : show an estimate of the number of genes a node of size Annotated would have if the significant genes were to be randomly selected from the gene universe.
#'   - pvalues : pvalue obtained after the test
#'
#' Column 'p.adj.weight01' represents the adjusted p-value for the weight01 algorithm.
#'
#' @author Bogdan Iancu - Genevia Technologies Oy
#' @return A data frame with enriched terms and p-values from Fisher's exact test, using different algorithms: elim, classic, weight01.
#' @export
run.topGO <- function(background, foreground, ontologies = c("BP"), organism, ID_type = "ENSEMBL", pAdjustMethod = "BH") {

  # test if data packages are installed
  if (!requireNamespace(organism, quietly = TRUE)) {
    stop(
      paste("Package", sQuote(organism), "must be installed to use this function."),
      call. = FALSE
    )
  }

  # genes.full = data.frame(readxl::read_excel(background, skip = 1))
  # gene.full.names = genes.full$Ensembl.ID
  #
  # genes.of.interest = data.frame(readxl::read_excel(foreground, skip = 1))
  # gene.of.interest.names = genes.of.interest$Ensembl.ID
  #
  print(paste0("pAdjustMethod: ", pAdjustMethod))
  print(paste0("ID_type: ",ID_type))

  gene.full.names = background
  gene.of.interest.names = foreground
  gene.list <- factor(as.integer(gene.full.names %in% gene.of.interest.names))
  names(gene.list) = gene.full.names


  table.go = as.list(ontologies)
  names(table.go) = ontologies

  result.topGO.elim <- NULL
  result.topGO.classic <- NULL
  result.topGO.weight01 <- NULL
  for (i in 1:length(table.go)) {
    cat(" @", sQuote(table.go[[i]]), "\n")
    ## prepare data
    cat("   Generating new object of class", sQuote("topGOdata"), "...\n")
    cat("   > Splitting GOTERM envrionment...\n")
    topGO::groupGOTerms()
    GOdata <- methods::new("topGOdata", ontology = ontologies[i], allGenes = gene.list, nodeSize = 10,
                   annot = topGO::annFUN.org, mapping = organism, ID = ID_type )
    cat("done\n")

    ## run tests
    cat("    topGO.elim...\n")
    result.topGO.elim <- topGO::runTest(GOdata, algorithm = "elim", statistic = "Fisher", cutOff = 0.05)
    cat("    topGO.classic...\n")
    result.topGO.classic <- topGO::runTest(GOdata, algorithm = "classic", statistic = "Fisher", cutOff = 0.05)
    cat("    topGO.weight01...\n")
    result.topGO.weight01 <- topGO::runTest(GOdata, algorithm = "weight01", statistic = "Fisher", cutOff = 0.05)
    #resultTopGO.elim

    ## look at results
    cat("    weight01 summary...\n")
    weight01.summary <- summary(attributes(result.topGO.weight01)$score <= 0.05)
    numsignif <- as.integer(weight01.summary[[3]])

    cat("    topGO GenTable...\n")
    table.go[[i]] <- topGO::GenTable( GOdata, Fisher.elim = result.topGO.elim,
                               Fisher.classic = result.topGO.classic, Fisher.weight01 = result.topGO.weight01,
                               orderBy = "Fisher.weight01", topNodes = numsignif)

  }

  # create the file with all the statistics from GO analysis
  topGO.results <- as.data.frame(rbind.fill(table.go))


  # list containing genes annotated to significant GO terms
  annotated.genes <- lapply(topGO.results$GO.ID, function(x) as.character(unlist(genesInTerm(object = GOdata, whichGO = x))))


  ## Selecting only the genes that are among foreground set:
  significant.genes <- lapply(annotated.genes, function(x) intersect(x, gene.of.interest.names))


  #performing BH correction on the weight01 p-values
  p.adj.weight01 <- round(p.adjust(topGO.results$Fisher.weight01, method = pAdjustMethod), digits = 5)

  #bind new p.adj.wieght01 col
  topGO.results$p.adj.weight01 <- p.adj.weight01

  #Binding the DEGs (ensembl ids) associated into the term into result table and switching the gene separator
  topGO.results$Genes <- significant.genes
  topGO.results$Genes <- unlist(lapply(topGO.results$Genes, function(x) paste(x, collapse = ","))) #To make the column compatible for wrapper formatting
  names(topGO.results)[2] = "Description" # To make the result more analogous with other enrichment approaches
  return(topGO.results)
}
