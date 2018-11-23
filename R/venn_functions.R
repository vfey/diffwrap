library(data.table)
library(purrr)
library(dplyr)
library(venn)
library(VennDiagram)
# compar.table1 <- as.data.frame(fread("./diffr_test_results_new_20181109_184403/test_run_CRPC-Advanced.PC_differential_expression.tsv", stringsAsFactors = TRUE))
# compar.table2 <- as.data.frame(fread("./diffr_test_results_new_20181109_184403/test_run_Advanced.PC-Localized.PC_differential_expression.tsv", stringsAsFactors = TRUE))
# compar.table3 <- as.data.frame(fread("./diffr_test_results_new_20181109_184403/test_run_Localized.PC-BPH_differential_expression.tsv", stringsAsFactors = TRUE))  

compar.table1 <- as.data.frame(fread("./data/finvector_test_data/DEGresults_M1_vs_M0_SIGN.csv", stringsAsFactors = TRUE))
compar.table2 <- as.data.frame(fread("./data/finvector_test_data/DEGresults_Mreg_Kiel_vs_M0_SIGN.csv", stringsAsFactors = TRUE))
compar.table3 <- as.data.frame(fread("./data/finvector_test_data/DEGresults_Mreg_Kiel_vs_M1_SIGN.csv", stringsAsFactors = TRUE))
compar.table4 <- as.data.frame(fread("./data/finvector_test_data/DEGresults_Mreg_Kiel_vs_Mreg_ML_SIGN.csv", stringsAsFactors = TRUE))

compar.table1 <- compar.table1 %>%  dplyr::select("Ensembl ID","HGNC symbol","Gene description","Gene biotype","Average expression","Log2 foldchange","P-value","Adjusted p-value")
compar.table2 <- compar.table2 %>%  dplyr::select("Ensembl ID","HGNC symbol","Gene description","Gene biotype","Average expression","Log2 foldchange","P-value","Adjusted p-value")
compar.table3 <- compar.table3 %>%  dplyr::select("Ensembl ID","HGNC symbol","Gene description","Gene biotype","Average expression","Log2 foldchange","P-value","Adjusted p-value")
compar.table4 <- compar.table4 %>%  dplyr::select("Ensembl ID","HGNC symbol","Gene description","Gene biotype","Average expression","Log2 foldchange","P-value","Adjusted p-value")

list.comp.tables2 <- list("A" = compar.table1, "B" = compar.table2)
list.comp.tables3 <- list("A" = compar.table1, "B" = compar.table2, "C" = compar.table3)
list.comp.tables4 <- list("A" = compar.table1, "B" = compar.table2, "C" = compar.table3, "D" =compar.table4)

join_vec <- c("Ensembl ID","HGNC symbol","Gene description","Gene biotype")
list.comp.tables4 <- list("A" = compar.table1, "B" = compar.table2, "C" = compar.table3, "D" =compar.table4)

diffr_venn <- function(list.comp.tables, join_vec) {
  
  #extract DEGs list from the contrast tables list 
  DEGs.list <- lapply(list.comp.tables4, function(x) {
    #extract hgnc colname
    hgnc.col <- names(x)[grep("symbol|hgnc", tolower(names(x)))]
    #extract DEGs from table x based on hgnc.col
    DEGs <- as.character(x[[hgnc.col]])
    #remove all entry with empty string
    DEGs <- DEGs[DEGs!=""]
    return(DEGs)
  })
  
  venn.diag <- venn.diagram(DEGs.list, 
                            col = "transparent", fill = c("darkblue", "deepskyblue", "darkturquoise", "darkorchid4"),  height = 8000, width = 8000, print.mode = c("raw", "percent"),
                            cat.cex = 1.5, cex = 2.5, imagetype = "png", filename = NULL)
  venn.sets.lists = venn(DEGs.list)
  venn.sets.intersections = attr(venn.sets.lists, "intersections")
  
  venn.tt <- venn.sets.lists %>% dplyr::select(-counts)
  
  #go over all data frames in list x and, for all columns that are not found in join_vec, change them to be identifiable by adding the list name to the column
  #this is done so that after the joins in the intersection tables are 
  create.individual.ids = function(x,y) {
    colnames(x)[which(!(colnames(x) %in% join_vec))] = paste0(colnames(x)[which(!(colnames(x) %in% join_vec))],"-",y)
    return(x)
  }
  
  list.venn.tables <- list()
  list.comp.tables = map2(list.comp.tables,names(list.comp.tables),create.individual.ids)
  for(j in 1:nrow(venn.tt)) {
    #if all entries in the venn truth table are 0 thek skip
    if (sum(venn.tt[j,]) == 0) {
      next
    }
    #if all entries are 1 then join them all based on the join_vec
    else if (all(venn.tt[j,] == 1)) {
      print(rownames(venn.tt)[j])
      inters.all <- list.comp.tables %>% reduce(inner_join, by = join_vec)
      inters.all <- inters.all[inters.all$`HGNC symbol` != "",]
      list.venn.tables[[rownames(venn.tt)[j]]] <- inters.all
      
      print(paste0(rownames(venn.tt)[j],"-",nrow(inters.all)))
    }
    #otherwise calculate the intersection for the ones that have index 1 and the union for those who have 0 and extract those entries which are in the intersection and not union
    #this gives the 
    else {
      join.part <- list.comp.tables[which(venn.tt[j,] == 1)]
      inters <- join.part %>% reduce(dplyr::inner_join, by = join_vec)
      anti.part <- list.comp.tables[which(venn.tt[j,] == 0)]
      outsect <- anti.part %>% reduce(dplyr::full_join, by = join_vec)
      res.join <- inters %>% dplyr::anti_join(outsect, by = join_vec)
      res.join <- res.join[res.join$`HGNC symbol` != "",]
      list.venn.tables[[rownames(venn.tt)[j]]] <- res.join
      print(paste0("clams, clams everywhere - ",rownames(venn.tt)[j],"--",nrow(res.join)))
    }
    
  }
  
  return(venn.diag)
}
