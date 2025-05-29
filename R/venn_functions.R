# Function: diffr_venn
#
# Author: Bogdan Iancu - Genevia Technologies Oy
#
# Arguments:
#         list.comp.tables  = list of DE tables, preferably a list of data.frames
#         join.vec = vector to perform the join operation on; corresponds to column names in the DE tables
#

# Output: returns the Venn diagram of DE tables and the Venn intersections tables
#
#

utils::globalVariables("counts")

#' Function to produce a Venn diagram of differentially expressed gene tables
#' @param list.comp.tables list of DE tables, preferably a list of data.frames
#' @param join.vec vector to perform the join operation on; corresponds to column names in the DE tables.
#'   Defaults to \code{c("ensembl_gene_id","gene_symbol","description","entrezgene_id")}.
#' @param .log Logical; should logging be done?
#' @details
#' The actual plot is produced in the main plotting function by means of 'grid::grid.draw()' using the plot object as input.
#'
#' @author Bogdan Iancu - Genevia Technologies Oy
#' @return A list containing the Venn diagram grid object and the intersected (or joined) input tables.
#' @export
diffr_venn <- function(list.comp.tables, join.vec, .log = FALSE) {

  # test if logging packages are installed
  if (.log && !requireNamespace("futile.logger", quietly = TRUE)) {
    stop(
      paste("Package", sQuote("futile.logger"), "must be installed to use the logging functionality."),
      call. = FALSE
    )
  }

  #extract DEGs list from the contrast tables list
  DEGs.list <- lapply(list.comp.tables, function(x) {
    #extract hgnc colname
    hgnc.col <- names(x)[grep("symbol|hgnc", tolower(names(x)))]
    #extract DEGs from table x based on hgnc.col
    DEGs <- as.character(x[[hgnc.col]])
    #remove all entry with empty string
    DEGs <- DEGs[DEGs!=""]
    #choose only unique DEGs
    DEGs <- unique(DEGs)
    return(DEGs)
  })

  fill.color = c("darkblue", "deepskyblue", "darkturquoise", "darkorchid4")
  if (.log) {
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  }
  venn.diag <- VennDiagram::venn.diagram(DEGs.list,
                            col = "transparent", fill = fill.color[1:length(DEGs.list)],  height = 8000, width = 8000, print.mode = c("raw", "percent"),
                            cat.cex = 1.2, cex = 2.2, imagetype = "png", filename = NULL)
  venn.sets.lists = venn::venn(DEGs.list)
  venn.sets.intersections = attr(venn.sets.lists, "intersections")

  venn.tt <- venn.sets.lists %>% dplyr::select(-counts)

  #go over all data frames in list x and, for all columns that are not found in join.vec, change them to be identifiable by adding the list name to the column
  #this is done so that after the joins in the intersection tables are
  create.individual.ids = function(x,y) {
    colnames(x)[which(!(colnames(x) %in% join.vec))] = paste0(colnames(x)[which(!(colnames(x) %in% join.vec))],"-",y)
    return(x)
  }

  list.venn.tables <- list()
  list.comp.tables = purrr::map2(list.comp.tables,names(list.comp.tables),create.individual.ids)
  for(j in 1:nrow(venn.tt)) {
    #if all entries in the venn truth table are 0 then skip
    if (sum(venn.tt[j,]) == 0) {
      next
    }
    #if all entries are 1 then join them all based on the join.vec
    else if (all(venn.tt[j,] == 1)) {
      print(rownames(venn.tt)[j])
      inters.all <- list.comp.tables %>% purrr::reduce(dplyr::inner_join, by = join.vec)
      hgnc.col <- names(inters.all)[grep("symbol|hgnc", tolower(names(inters.all)))]
      inters.all <- inters.all[inters.all[hgnc.col] != "",]
      list.venn.tables[[rownames(venn.tt)[j]]] <- inters.all

      print(paste0(rownames(venn.tt)[j],"-",nrow(inters.all[unique(inters.all[[hgnc.col]]),])))
    }
    #otherwise calculate the intersection for the ones that have index 1 and the union for those who have 0 and extract those entries which are in the intersection and not in the union
    #this gives the intersection table for the Venn sections seen in the Venn diagram
    else {
      join.part <- list.comp.tables[which(venn.tt[j,] == 1)]
      inters <- join.part %>% purrr::reduce(dplyr::inner_join, by = join.vec)
      anti.part <- list.comp.tables[which(venn.tt[j,] == 0)]
      outsect <- anti.part %>% purrr::reduce(dplyr::full_join, by = join.vec)
      res.join <- inters %>% dplyr::anti_join(outsect, by = join.vec)
      hgnc.col <- names(res.join)[grep("symbol|hgnc", tolower(names(res.join)))]
      res.join <- res.join[res.join[[hgnc.col]] != "",]
      list.venn.tables[[rownames(venn.tt)[j]]] <- res.join
      print(paste0("Venn sections summary - ",rownames(venn.tt)[j],"--",nrow(res.join[unique(res.join[[hgnc.col]]),])))
    }

  }
  res.venn <- list("venn.diagram" = venn.diag, "venn.sections" = list.venn.tables)
  return(res.venn)
}
