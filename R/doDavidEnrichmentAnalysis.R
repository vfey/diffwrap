library("RDAVIDWebService")

doDavidEnrichmentAnalysis = function(background.ensembl.ids,
                                     foreground.ensembl.ids,
                                     email.address,
                                     url.address = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/",
                                     time.out.value = 60000,
                                     annotation.category = "GOTERM_BP_FAT",
                                     pval.thr = 0.05,
                                     max.gene.set.size = '',
                                     save.result.table = FALSE,
                                     out.dir=NULL,
                                     analysis.name,
                                     contrast.name){
  
  
  #Setting up the david interface (probably wiser to do outside a function:
  david <- DAVIDWebService$new(email= email.address, url = url.address)
  
  
  setTimeOut(david, time.out.value)
  
  #Setting up gene lists /ensembl_ids):
  FG <- addList(david, foreground.ensembl.ids, idType="ENSEMBL_GENE_ID", listName="isClass", listType="Gene")
  BG <- addList(david, background.ensembl.ids, idType="ENSEMBL_GENE_ID", listName="all", listType="Background")
  
  cat("The species is ", getSpecieNames(david), "\n" )
  
  cat("Proportion ", FG$inDavid, " of the foreground genes in DAVID \n")
  cat("Proportion ", BG$inDavid, " of the background genes in DAVID \n")
  
  setAnnotationCategories(david, annotation.category)
  FuncAnnotChart <- getFunctionalAnnotationChart(david, threshold=pval.thr, count=3L) #PValue < 0.05, at least 3 genes
  
  #Filtering and ordering the result based on adjusted p-value
  filtered_FuncAnnotChart = FuncAnnotChart[FuncAnnotChart$Benjamini < pval.thr,]
  filtered_FuncAnnotChart = filtered_FuncAnnotChart[order(filtered_FuncAnnotChart$Benjamini),]
  
  #filtering out lines with Large Gene sets:
  if(max.gene.set.size != ''){
    filtered_FuncAnnotChart = filtered_FuncAnnotChart[filtered_FuncAnnotChart$Count <= max.gene.set.size,]
  }
  
  ##Requires medseqr!! Replacing ENSMBL IDS with symbolic names if found (otherwise ENSEMBL ID is kept)
  for( i in 1:dim(filtered_FuncAnnotChart)[1]){
    genes = unlist(strsplit(filtered_FuncAnnotChart$Genes[i], split = ", "))
    syms = convertId2(genes, "Human")
    syms[is.na(syms)] = genes[is.na(syms)]
    new_names = paste(syms,collapse="/") #writing with /-separator (easier to read from file later)
    filtered_FuncAnnotChart$Genes[i] = new_names
  }
  
  
  #splitting the Term-column of original output into two columns
  GO_terms = sapply(strsplit(as.character(filtered_FuncAnnotChart$Term),'~'), "[", 1)
  descriptions = sapply(strsplit(as.character(filtered_FuncAnnotChart$Term),'~'), "[", 2)
  
  result.table = data.frame(Term = GO_terms,
                            Description = descriptions)
  result.table = cbind(result.table, filtered_FuncAnnotChart[, colnames(filtered_FuncAnnotChart) != "Term"])
  
  result.table = result.table[, c(1:2, 4, 8:10, 6, 13, 7)] #Taking only most interesting columns
  
  if(save.result.table){
    filename=paste0(contrast.name,".DAVID.", annotation.category,".enrichment_table.xls" )
    cat("Saving DAVID result table into ", file.path(out.dir, paste0(analysis.name, ".", filename)), "...\n")
    write.table(result.table,
                file = file.path(out.dir, paste0(analysis.name, ".", filename)), row.names = FALSE,
                col.names = TRUE)
  }
  return(filtered_FuncAnnotChart) # or would the result.table be better?
}

