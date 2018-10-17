library(topGO)
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

# Output: table with enriched terms and p-values from Fisher's exact test, using differeng algorithms: elim, classic, weight01. 
#         
#
# Details: columns in the final table produced by topGo (their description is fetched from topGo documentation): 
#                                    - Annotated : number of genes in org.Hs.eg.db which are annotated with the GO-term.
#                                    - Significant : number of genes belonging to your input which are annotated with the GO-term.
#                                    - Expected : show an estimate of the number of genes a node of size Annotated would have if the significant genes were to be randomly selected from the gene universe.
#                                    - pvalues : pvalue obtained after the test

#           column: p.adj.weight01 - represents the adjusted p-value for the weight01 algorithm    
#           requires packages: topGO, readxl, org.Hs.eg.db
run.topGO <- function(background, foreground, ontologies = c("BP"), organism, ID_type = "ENSEMBL", pAdjustMethod = "BH") {
  
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
    
    ## prepare data
    GOdata <- new( "topGOdata", ontology = ontologies[i], allGenes = gene.list, nodeSize = 10,
                   annot = annFUN.org , mapping = organism, ID = ID_type )
    
    ## run tests
    result.topGO.elim <- runTest(GOdata, algorithm = "elim", statistic = "Fisher", cutOff = 0.05)
    result.topGO.classic <- runTest(GOdata, algorithm = "classic", statistic = "Fisher", cutOff = 0.05)
    result.topGO.weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "Fisher", cutOff = 0.05)
    #resultTopGO.elim
    
    ## look at results
    weight01.summary <- summary(attributes(result.topGO.weight01)$score <= 0.05)
    numsignif <- as.integer(weight01.summary[[3]])
    
    table.go[[i]] <- GenTable( GOdata, Fisher.elim = result.topGO.elim,
                               Fisher.classic = result.topGO.classic, Fisher.weight01 = result.topGO.weight01,
                               orderBy = "Fisher.weight01", topNodes = numsignif)
    
  }
  
  # create the file with all the statistics from GO analysis
  topGO.results <- rbind.fill(table.go)
  
  #performing BH correction on the weight01 p-values
  p.adj.weight01 = round(p.adjust(topGO.results$Fisher.weight01, method = pAdjustMethod), digits = 5)
  
  #bind new p.adj.wieght01 col
  topGO.results$p.adj.weight01 = p.adj.weight01
  
  return(topGO.results)
} 