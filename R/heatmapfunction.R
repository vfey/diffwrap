#' Function to define quantile breaks to be used for changing the palette of the heatmap.
#' @note Function by: Kamil Slowikowski, https://github.com/slowkow/slowkow.com/blob/master/_rmd/2017-02-16-heatmap-tutorial.R,
#' @param xs numeric vector to calculate quantiles for
#' @param n desired length of the numeric vector of probabilities (see ?quantile)
#'
#' @seealso [quantile()]
quantile_breaks <- function(xs, n = 20) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#' Function which reorders the levels of a column of a data frame specified as a factor
#' @details The desired order is a vector containing the levels of the factor in the desired order.
#' @note Written by https://stackoverflow.com/users/1701600/boern.
#' @param df data frame to be processed
#' @param column name of the column to be reordered
#' @param desired_level_order vector of factor levels in the desired order
reorderFactors <- function(df, column = "my_column_name",
                           desired_level_order) {

  x = df[[column]]
  lvls_src = levels(x)

  idxs_target <- vector(mode = "numeric", length = 0)
  for (target in desired_level_order) {
    idxs_target <- c(idxs_target, which(lvls_src == target))
  }

  x_new <- factor(x,levels(x)[idxs_target])

  df[[column]] <- x_new

  return(df)
}

#' Function to create the annotation colour list used in the heatmap
#' @param clinical.mat matrix with clinical annotation values in (clinical category, samples) format
make_pheatmap_anno_color = function(clinical.mat) {
  #create the annotation colour list
  anno.vars = list()
  for (j in 1:length(colnames(clinical.mat))) {
    #find unique annotation terms and exclude NAs
    unique.anno.terms = unique(clinical.mat[,j])
    anno.vars[[j]] = unique.anno.terms[!is.na(unique.anno.terms)]
  }

  #annotation palette of 11 colours for the categories in the annotation
  anno.palette  = c('#FFC125','#8B0A50','#8A2BE2','#FFA07A','#8B4500','#CD5555','#FF8C00','#8B7500','#CD6090','#8B0000','#53868B')

  anno.color = list()
  #choose number of colors from the annotation palette to be equal with the annotation.vars.count
  for (i in 1:length(colnames(clinical.mat))) {
    anno.color[[names(clinical.mat)[i]]] = anno.palette[1:length(anno.vars[[i]])]
    names(anno.color[[i]]) = levels(as.factor(clinical.mat[,i]))
  }
  return(anno.color)
}


# Function: correlogram_pheatmap
#
# Author: Bogdan Iancu - Genevia Technologies Oy
#
# Arguments:
#         expr.mat  = differential expression matrix in  (genes, samples) format
#

# Output: pheatmap object that represents a correlogram
#
#' Function to create a correlogram pheatmap, i.e., a plot to check randomness in the data set.
#' @param expr.mat differential expression matrix in  (genes, samples) format
#' @param clinical.mat matrix with clinical annotation values in (clinical category, samples) format
#' @param scale.fl character indicating if values should be centered and scaled in either
#' the row direction or the column direction, or none (values %in%
#' ("row","column","none"), default = none)
#' @param legend.fl logical to determine if legend should be drawn or not (default = TRUE)
#' @param row.clust boolean values determining if rows should be clustered
#' @param col.clust boolean values determining if cols should be clustered
#' and corresponding heatmap outputed, assuming clinical.mat is provided
#' the colors of the heatmap, otherwise min-max breaks are used by default
#' @param signif.stars.fl boolean determining whether significance stars of p-values are shown in the correlogram (default = FALSE)
#' @param cell.size double determining the widh and height of the cell and the row/col font size (default = 8)
#' @param font.size double determining the font size (default = 10)
#' @param color.blind.pal string determining the RColorBrewer color blind palette (default = "PuOr");
#' other option can be visualized with the following command: brewer.pal.info[brewer.pal.info$colorblind,]
#' @param color.extremes character vector of length 2 giving the two extremes of a user-defined colour palette
#' varying from the first hue to the second via white.
#' @param main.correl Character; main title of the Correlogram
#' @param sample.correl Logical;
#' @details
#' Calculates Pearson's rank correlation coefficients for all possible pairs of the input matrix and plots the resulting
#' correlation matrix, either with significance stars and coloured according to correlation or without stars and coloured
#' according to clinical annotation.
#'
#' @author Bogdan Iancu - Genevia Technologies Oy
#' @return Returns a pheatmap plot object used in the diffr_pheatmap() function.
#'
correlogram_pheatmap = function(expr.mat, clinical.mat, scale.fl = "none", legend.fl = TRUE,
                                row.clust = TRUE, col.clust = TRUE,
                                signif.stars.fl = FALSE, cell.size = 8,
                                font.size = 11, color.blind.pal = "PuOr",
                                color.extremes = c("#3182BD", "#E6550D"),
                                main.correl = "Correlogram", sample.correl = FALSE) {
  #calculate the correlation matrix using Hmisc package

  corr_hmap = Hmisc::rcorr(expr.mat)
  cor.mat = as.matrix(corr_hmap$r)
  pv.mat = as.matrix(corr_hmap$P)
  mat.breaks.correl = NULL
  #correlogram scale and breaks
  if(!sample.correl) {
    mat.breaks.correl  = seq(-1, 1, by = 0.05)
  } else {
    mat.breaks.correl = seq(0, 1, by = 0.05)
  }

  if (!is.null(color.blind.pal)) {
    colour.correl = colorRampPalette(rev(brewer.pal(n = 11, name = color.blind.pal)))(length(mat.breaks.correl) - 1)
  } else {
    colour.correl = colorRampPalette(c(color.extremes[1], "white", color.extremes[2]))(length(mat.breaks.correl) - 1)
  }

  #define new clustering distance
  dissimilarity <- 1 - cor.mat
  dist.clust <- as.dist(dissimilarity)

  if (scale.fl != "none") {
    mat.breaks.correl = NA
  }
  sign.stars = NULL
  p.cor = NULL
  if (signif.stars.fl) {
    # define significance levels
    sign.stars <- ifelse(pv.mat < .001, "***", ifelse(pv.mat < .01, "** ", ifelse(pv.mat < .05, "* ", " ")))
    sign.stars <- ifelse(is.na(sign.stars),"",sign.stars)
    #build the correlogram using pheatmap
    p.cor = pheatmap::pheatmap(cor.mat,
                               show_colnames = T, show_rownames = T,
                               cluster_rows = T, cluster_cols = T,
                               legend = legend.fl, border_color = "white",
                               scale = scale.fl, color = colour.correl,
                               clustering_distance_rows = dist.clust, clustering_distance_cols = dist.clust,
                               display_numbers = sign.stars, number_color = "#FFFF00", breaks = mat.breaks.correl,
                               cellwidth = cell.size, cellheight = cell.size, fontsize = font.size,
                               fontsize_row = cell.size, fontsize_col = cell.size, silent = TRUE, main = main.correl)

  }
  else {
    #print("mat.breaks.correl")
    #print(mat.breaks.correl)
    anno.color = make_pheatmap_anno_color(clinical.mat = clinical.mat)
    browser()

    if(sample.correl) {

      p.cor = pheatmap::pheatmap(cor.mat,
                                 show_colnames = T, show_rownames = T,
                                 cluster_rows = T, cluster_cols = T,
                                 legend = legend.fl, border_color = "white",
                                 scale = scale.fl, color = colour.correl,
                                 clustering_distance_rows = dist.clust, clustering_distance_cols = dist.clust,
                                 display_numbers = FALSE,
                                 cellwidth = cell.size, cellheight = cell.size, fontsize = font.size,
                                 fontsize_row = cell.size, fontsize_col = cell.size, silent = TRUE, main = main.correl,
                                 annotation_color = anno.color,
                                 annotation = clinical.mat, annotation_row = clinical.mat)

    } else{
      p.cor = pheatmap::pheatmap(cor.mat,
                                 show_colnames = T, show_rownames = T,
                                 cluster_rows = T, cluster_cols = T,
                                 legend = legend.fl, border_color = "white",
                                 scale = scale.fl, color = colour.correl,
                                 clustering_distance_rows = dist.clust, clustering_distance_cols = dist.clust,
                                 display_numbers = FALSE, breaks = mat.breaks.correl,
                                 cellwidth = cell.size, cellheight = cell.size, fontsize = font.size,
                                 fontsize_row = cell.size, fontsize_col = cell.size, silent = TRUE, main = main.correl,
                                 annotation_color = anno.color)
    }


    #
    # p = pheatmap::pheatmap(expr.mat,
    #                        show_colnames = T, show_rownames = T,
    #                        cluster_rows = row.clust, cluster_cols = col.clust,
    #                        legend = legend.fl, border_color = "white",
    #                        scale = scale.fl, color = colour,annotation_color = anno.color,
    #                        annotation = clinical.mat, breaks = breaks.hm,
    #                        cellwidth = cell.size, cellheight = cell.size, fontsize = font.size,
    #                        fontsize_row = cell.size, fontsize_col = cell.size, silent = TRUE, main = "Heatmap of genes over samples")
  }

  return(p.cor)
}




# Function: diffr_pheatmap()
#
# Author: Bogdan Iancu - Genevia Technologies Oy
#
# Arguments:
#         expr.mat = expression matrix in regular  (genes, samples) format
#         clinical.mat = clinical matrix, by default a column with sample types
#         scale.fl = character indicating if values should be centered and scaled in either
#                   the row direction or the column direction, or none (values %in%
#                   ("row","column","none"), default = none)
#          legend.fl = logical to determine if legend should be drawn or not (default = TRUE)
#          row.clust = boolean values determining if rows should be clustered
#          col.clust = boolean values determining if cols should be clustered
#          biserial.fl = boolean values determining if biserial correlation is calculated
#                       and corresponding heatmap outputed, assuming clinical.mat is provided
#          quantile.breaks.fl = boolean values determining if quantile breaks are used to change
#                               the colors of the heatmap, otherwise min-max breaks are used by default
#          signig.stars.fl = boolean determining whether significance stars of p-values are shown in the correlogram (default = FALSE)
#          cell.size = double determining the widh and height of the cell and the row/col font size (default = 8)
#          font.size = double determining the font size (default = 10)
#          color.blind.pal = string determining the RColorBrewer color blind palette (default = "PuOr");
#                           other option can be visualized with the following command: brewer.pal.info[brewer.pal.info$colorblind,]


# Output: list of pheatmap objects that should be preferably handled with lapply
#
# Details: calls function quantile_breaks
#          calls function reorderFactors written by https://stackoverflow.com/users/1701600/boern,
#                                        which reorders the levels of a column of a data frame specified as a factor,
#                                        the desired order is a vector containing the levels of the factor
#          requires R-packages: 'pheatmap','RColorBrewer','Hmisc','ltm'
#
#
#
#' Function to create a heatmap from differential gene expression values
#' @param expr.mat differential expression matrix in  (genes, samples) format
#' @param clinical.mat matrix with clinical annotation values in (clinical category, samples) format
#' @param scale.fl character indicating if values should be centered and scaled in either
#' the row direction or the column direction, or none (values %in%
#' ("row","column","none"), default = none)
#' @param legend.fl logical to determine if legend should be drawn or not (default = TRUE)
#' @param row.clust boolean values determining if rows should be clustered
#' @param col.clust boolean values determining if cols should be clustered
#' @param biserial.fl = boolean values determining if biserial correlation is calculated
#' and corresponding heatmap outputed, assuming clinical.mat is provided
#' @param quantile.breaks.fl boolean values determining if quantile breaks are used to change
#' the colors of the heatmap, otherwise min-max breaks are used by default
#' @param signif.stars.fl boolean determining whether significance stars of p-values are shown in the correlogram (default = FALSE)
#' @param cell.size double determining the widh and height of the cell and the row/col font size (default = 8)
#' @param font.size double determining the font size (default = 10)
#' @param color.blind.pal string determining the RColorBrewer color blind palette (default = "PuOr");
#' @param color.extremes character vector of length 2 giving the two extremes of a user-defined colour palette
#' varying from the first hue to the second via white.
#' @param main \code{character}. Main plot title.
#' @param add.main \code{character}. Optional text added in parenthesis to the plot title.
#' @param filt.info \code{character}. Optional information on the filtering strategy added in parenthesis to the plot title.
#' @details
#' The plot is produced as a side effect of function 'pheatmap_plots()' by printing the individual
#' plot objects generated by 'diffr_pheatmap()'.
#'
#' @author Bogdan Iancu - Genevia Technologies Oy
#' @return Returns a list of pheatmap plot objects used in the pheatmap_plots() function.
#' @export
diffr_pheatmap = function(expr.mat, clinical.mat,
                          scale.fl = "none", legend.fl = TRUE,
                          row.clust = TRUE, col.clust = TRUE,
                          biserial.fl = FALSE, quantile.breaks.fl = FALSE,
                          signif.stars.fl = FALSE, cell.size =8,
                          font.size = 10, color.blind.pal = "PuOr",
                          color.extremes = c("#3182BD", "#E6550D"),
                          main = NULL, add.main = NULL, filt.info = NULL) {

  # make additional information for main title
  if (!is.null(add.main) && is.null(filt.info)) {
    main.plus <- paste0("(", add.main, ")")
  } else if (!is.null(filt.info) && is.null(add.main)) {
    main.plus <- paste0("(", filt.info, ")")
  } else if (!is.null(add.main) && !is.null(filt.info)) {
    main.plus <- paste0("(", filt.info, ", ", add.main, ")")
  } else {
    main.plus <- NULL
  }

  #initially the breaks parameter is set to NA
  breaks.hm = NA

  #calculate min-max breaks
  palette.length <- 100
  if (!is.null(color.blind.pal)) {
    colour = colorRampPalette(rev(brewer.pal(n = 11, name = color.blind.pal)))( palette.length )
  } else {
    colour = colorRampPalette(c(color.extremes[1], "white", color.extremes[2]))( palette.length )
  }
  # # length(breaks) == length(paletteLength) + 1
  # # use floor and ceiling to deal with even/odd length pallettelengths
  mat.breaks <- c(seq(min(expr.mat), 0, length.out=ceiling(palette.length/2) + 1),
                  seq(max(expr.mat)/palette.length, max(expr.mat), length.out=floor(palette.length/2)))

  mat.breaks <- seq(min(expr.mat), max(expr.mat), by = 0.05)

  #colour = colorRampPalette(rev(brewer.pal(n = 11, name = color.blind.pal)))(length(mat.breaks) - 1)

  #change the colours of the matrix based on quantile breaks, the default is above
  if (quantile.breaks.fl) {
    mat.breaks = quantile_breaks(as.matrix(expr.mat), n = 11)
    if (!is.null(color.blind.pal)) {
      colour = colorRampPalette(rev(brewer.pal(n = 11, name = color.blind.pal)))( length(mat.breaks) - 1 )
    } else {
      colour = colorRampPalette(c(color.extremes[1], "white", color.extremes[2]))( length(mat.breaks) - 1 )
    }
  }
  #cat("mat.breaks")
  #print(mat.breaks)
  #cat("colour")
  #print(colour)

  #for the non-scaled case, if "use.breaks" is mentioned, then use breaks to create the heatmap, by default it only use breaks for colours and allows heatmap to do the breaks automatically
  if (scale.fl == "none") {
    breaks.hm = mat.breaks
  }

  p = NULL

  if (missing(clinical.mat)) {

    p = pheatmap::pheatmap(expr.mat,
                           show_colnames = T, show_rownames = T,
                           cluster_rows = row.clust, cluster_cols = col.clust,
                           legend = legend.fl, border_color = "white",
                           scale = scale.fl, color = colour, breaks = breaks.hm,
                           cellwidth = cell.size, cellheight = cell.size, fontsize = font.size,
                           fontsize_row = cell.size, fontsize_col = cell.size, silent = TRUE,
                           main = main)

  }
  else {
    anno.color = make_pheatmap_anno_color(clinical.mat = clinical.mat)

    # #create the annotation colour list
    # anno.vars = list()
    # for (j in 1:length(colnames(clinical.mat))) {
    #   #find unique annotation terms and exclude NAs
    #   unique.anno.terms = unique(clinical.mat[,j])
    #   anno.vars[[j]] = unique.anno.terms[!is.na(unique.anno.terms)]
    # }
    #
    # #annotation palette of 11 colours for the categories in the annotation
    # anno.palette  = c('#FFC125','#8B0A50','#8A2BE2','#FFA07A','#8B4500','#CD5555','#FF8C00','#8B7500','#CD6090','#8B0000','#53868B')
    #
    # anno.color = list()
    # #choose number of colors from the annotation palette to be equal with the annotation.vars.count
    # for (i in 1:length(colnames(clinical.mat))) {
    #   anno.color[[names(clinical.mat)[i]]] = anno.palette[1:length(anno.vars[[i]])]
    #   names(anno.color[[i]]) = levels(as.factor(clinical.mat[,i]))
    # }

    #if biserial flag is TRUE it calculates the biserial correlation based on the first column of the annotation file
    if (!biserial.fl) {

      p = pheatmap::pheatmap(expr.mat,
                             show_colnames = T, show_rownames = T,
                             cluster_rows = row.clust, cluster_cols = col.clust,
                             legend = legend.fl, border_color = "white",
                             scale = scale.fl, color = colour, annotation_color = anno.color,
                             annotation = clinical.mat, breaks = breaks.hm,
                             cellwidth = cell.size, cellheight = cell.size, fontsize = font.size,
                             fontsize_row = cell.size, fontsize_col = cell.size, silent = TRUE,
                             main = ifelse(is.null(main), paste("Heatmap of genes over samples", main.plus), main))

    }
    else {

      p = pheatmap::pheatmap(expr.mat,
                             show_colnames = T, show_rownames = T,
                             cluster_rows = row.clust, cluster_cols = col.clust,
                             legend = legend.fl, border_color = "white",
                             scale = scale.fl, color = colour, annotation_color = anno.color,
                             annotation = clinical.mat, breaks = breaks.hm,
                             cellwidth = cell.size, cellheight = cell.size, fontsize = font.size,
                             fontsize_row = cell.size, fontsize_col = cell.size, silent = TRUE,
                             main = ifelse(is.null(main), paste("Heatmap of genes over samples", main.plus), main))


      #get the annotation that will be transformed into a categorial variable
      biserial.vars = unique(clinical.mat[,1])
      #transform annotation into categorical variable
      biserial.vec = ifelse(!grepl(biserial.vars[1],clinical.mat[,1]),1,0)
      #calculate biserial correl for every gene
      biserial.cor.vec = apply(expr.mat, 1, function(x) {ltm::biserial.cor(as.numeric(x), biserial.vec,level = 2)})
      #build the annotation for the correlations
      anno.correl = ifelse(biserial.cor.vec > 0.75, "0.75<r<=1",
                           ifelse(biserial.cor.vec > 0.5 & biserial.cor.vec <= 0.75, "0.5<r<=0.75",
                                  ifelse(biserial.cor.vec > 0.25 & biserial.cor.vec <= 0.5 , "0.25<r<=0.5",
                                         ifelse(biserial.cor.vec <= 0.25 & biserial.cor.vec > 0, "0<r<=0.25",
                                                ifelse(biserial.cor.vec <= 0 & biserial.cor.vec > -0.25, "-0.25<r<=0",
                                                       ifelse(-0.5 < biserial.cor.vec & biserial.cor.vec <= -0.25, "-0.5<r<=-0.25",
                                                              ifelse(-0.75 < biserial.cor.vec & biserial.cor.vec <= -0.5,"-0.75<r<=-0.5","-1<r<=-0.75")))))))
      #build the correlation annotation
      anno.correl = as.data.frame(anno.correl)
      colnames(anno.correl)[1] = "Correlation"
      anno.correl$Correlation = as.factor(anno.correl$Correlation)
      anno.correl = reorderFactors(anno.correl,"Correlation", c("-1<r<=-0.75","-0.75<r<=-0.5","-0.5<r<=-0.25","-0.25<r<=0", "0<r<=0.25","0.25<r<=0.5","0.5<r<=0.75"))


      #build the heatmap with biserial correlation annotation
      p.biserial = pheatmap::pheatmap(expr.mat,
                                      show_colnames = T, show_rownames = T,
                                      cluster_rows = row.clust, cluster_cols = col.clust,
                                      legend = legend.fl, border_color = "white",
                                      scale = scale.fl, color = colour, breaks = breaks.hm, annotation_color = anno.color,
                                      annotation = clinical.mat, annotation_row = anno.correl,
                                      cellwidth = cell.size, cellheight = cell.size, fontsize = font.size,
                                      fontsize_row = cell.size, fontsize_col = cell.size, silent = TRUE,
                                      main = main)
    }

  }

  t_data_heat_map = t(expr.mat)

  #create the heatmap list and populate it
  heatmap.list = list()
  heatmap.list[["regular"]] = p
  heatmap.list[["correlogram"]] = correlogram_pheatmap(t_data_heat_map, clinical.mat,
                                                       scale.fl = scale.fl,
                                                       signif.stars.fl = signif.stars.fl,
                                                       color.blind.pal = color.blind.pal,
                                                       color.extremes = color.extremes,
                                                       main.correl = ifelse(is.null(main),
                                                                            paste("Gene Correlogram", main.plus), main))
  cell.new = cell.size + 2
  font.new = font.size + 2
  heatmap.list[["correlogram.small"]] = correlogram_pheatmap(t_data_heat_map, clinical.mat,
                                                             scale.fl = scale.fl,
                                                             signif.stars.fl = signif.stars.fl,
                                                             color.blind.pal = color.blind.pal,
                                                             color.extremes = color.extremes,
                                                             main.correl = ifelse(is.null(main),
                                                                                  paste("Gene Correlogram", main.plus), main),  cell.size =  cell.new, font.size = font.new)

  heatmap.list[["correlogram.sample"]] = correlogram_pheatmap(as.matrix(expr.mat), clinical.mat,
                                                              scale.fl = scale.fl,
                                                              signif.stars.fl = signif.stars.fl,
                                                              color.blind.pal = color.blind.pal,
                                                              color.extremes = color.extremes,
                                                              main.correl = ifelse(is.null(main),
                                                                                   paste("Sample Correlogram", main.plus), main), sample.correl = TRUE, cell.size =  cell.new, font.size = font.new)


  if (biserial.fl) {
    heatmap.list[["biserial.cor"]] = p.biserial
  }
  return(heatmap.list)
}

