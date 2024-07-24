
  # Function: heatmap_plots
  #
  # Author: Bogdan Iancu - Genevia Technologies Oy
  #
  # Arguments:
  #         d3 = differentia expression matrix in  (genes, samples) format
  #

  # Output: list of pheatmap objects that should be preferably handled with lapply
  #
  # Details: calls function diffr_pheatmap
  #

#' Function to generate heatmap of gene expression values
#' @details
#' The actual heatmap is produced as a side effect by printing the individual plot components.
#' @param d3 \code{data.frame}. Data frame containing all necessary columns to generate a heatmap with gene labels (at least p-values, FDR values, log-ratios and gene symbols or other IDs)
#' @param id \code{character}. Name of the gene ID column. Can be the same as `sym.col` but usually refers to an additional column with, e.g., Ensembl Gene IDs.
#' @param sym.col \code{character}. Name of column with gene symbols, e.g., HGNC Symbols.
#' @param samp.info \code{data.frame}. samp.info object containing information of the project's sample sheet
#' @param samples \code{character}. Name of the column in 'samp.info' containing sample names. If 'samp.info' is not supplied
#'     vector of sample names.
#' @param groups Factor of sample groups for colouring and legend.
#' @param sample.plot.names \code{character}. Optional name of a column with "nice" sample names for plotting.
#' @param main \code{character}. Main plot title. (Will be complemented with additional information, e.g., 'FDR' when labelling according to and FDR threshold.)
#' @param p.thr \code{numeric}. Plotted values with a P-Value below this threshold will be labelled in the P-Value plot.
#' @param fdr.thr \code{numeric}. Plotted values with a FDR below this threshold will be labelled in the FDR plot.
#' @param logfc.thr \code{numeric}. Plotted (`abs`olute) values above this threshold will have bigger dots.
#' @return Returns a list of pheatmap plot objects
#' @note This function is not exported as it is tailor-made for the diffwrap workflow (and called in the main wrapper function).
#'     The 'diffr_pheatmap' function used here is the main work horse and exported.
pheatmap_plots <-
    function(d3, id, sym.col="gene_symbol", samp.info = samp.info, samples, groups, sample.plot.names, main=NULL, p.thr=0.05, fdr.thr=0.05, logfc.thr=1)
    {
      # elim row with duplicate sym.col entries, set rownames to gene_symbol names

      d3 = as.data.table(d3)
      d3 = as.data.frame(unique(d3, by = sym.col))
      rownames(d3) = as.character(d3[,sym.col])

      # find columns with P-Values or FDR values
      pv.col = names(d3)[grep("^p\\.{0,1}val[e-u]{0,2}$", tolower(names(d3)))]
      fdr.col = names(d3)[grep("^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", tolower(names(d3)))]

      g.l = list()
      #cat("Heatmap plots...\n")

      #make two data frames for significant p-values and significant adj.p.values
      dat.sign.pv = d3[d3[[pv.col]] < 0.05,]
      dat.sign.fdr = d3[d3[[fdr.col]] < 0.05,]

      #print(dat.sign.pv)
      samp.info = as.data.frame(samp.info)
      #select only the columns representing the samples
      samp.names = as.character(samp.info[, samples])
      #print(samp.names)

      dat.sign.pv = dat.sign.pv[samp.names]
      dat.sign.fdr = dat.sign.fdr[samp.names]

      #get the heatmap column annotation
      samp.anno = as.data.frame(groups)
      rownames(samp.anno) = samp.names
      colnames(samp.anno)[1] = "Sample Class"

      #if pretty names for the plots are available change the names to the pretty names
      if (!is.null(sample.plot.names)) {
        colnames(dat.sign.pv) = sample.plot.names
        colnames(dat.sign.fdr) = sample.plot.names
        rownames(samp.anno) = sample.plot.names
      }
      #create lists of pheatmaps for significant p-values and adjusted p-values respectively
      pv_hm_list = list()
      fdr_hm_list = list()

      if (nrow(dat.sign.fdr) > 0) {
        if (nrow(dat.sign.fdr) > 100) {
          dat.sign.fdr = as.data.frame(dat.sign.fdr[1:100,])
          #make 2nd option smaller pheatmap for genes
          dat.sign.fdr.small = as.data.frame(dat.sign.fdr[1:50,])
        }
        fdr_hm_list[["row"]] = diffr_pheatmap(dat.sign.fdr, clinical.mat = samp.anno, scale.fl = "row")
        fdr_hm_list[["none"]] = diffr_pheatmap(dat.sign.fdr, clinical.mat = samp.anno, scale.fl = "none")
        fdr_hm_list[["rowsmall"]] = diffr_pheatmap(dat.sign.fdr.small, clinical.mat = samp.anno)
        fdr_hm_list[["nonesmall"]] = diffr_pheatmap(dat.sign.fdr.small, clinical.mat = samp.anno)

      } else {
        print("There are 0 entries with significant adjusted P-values in the differential expression data frame")
        print("Checking P-value entries...")

        if (nrow(dat.sign.pv) > 0) {
          if (nrow(dat.sign.pv) > 100) {
            dat.sign.pv = as.data.frame(dat.sign.pv[1:100,])
            dat.sign.pv.small = as.data.frame(dat.sign.pv[1:50,])
          }
          pv_hm_list[["row"]] = diffr_pheatmap(dat.sign.pv, clinical.mat = samp.anno, scale.fl = "row")
          pv_hm_list[["none"]] = diffr_pheatmap(dat.sign.pv, clinical.mat = samp.anno, scale.fl = "none")
          pv_hm_list[["rowsmall"]] = diffr_pheatmap(dat.sign.pv.small, clinical.mat = samp.anno, scale.fl = "row")
          pv_hm_list[["nonesmall"]] = diffr_pheatmap(dat.sign.pv.small, clinical.mat = samp.anno, scale.fl = "none")
          } else {
          print("There are 0 entries with significant P-values in the differential expression data frame")
        }

        #if there are no entries with significant adjusted p-values, but there are with significant p-values, then save in the g.l list the row-scaled regular and non-scaled correlograms
        if (length(pv_hm_list) != 0) {
          gl.pv = list()
          gl.pv$regular = pv_hm_list$row$regular
          gl.pv$gene.correlogram = pv_hm_list$none$correlogram
          gl.pv$gene.correlogram.small = pv_hm_list$nonesmall$correlogram.small
          gl.pv$samp.correlogram = pv_hm_list$none$correlogram.sample
          grid::grid.newpage()
          print(gl.pv$regular )
          grid::grid.newpage()
          print(gl.pv$gene.correlogram)
          grid::grid.newpage()
          print(gl.pv$gene.correlogram.small)
          grid::grid.newpage()
          print(gl.pv$samp.correlogram)
          g.l[["pval"]] = gl.pv
        }

      }
      #if there are any entries with significant adjusted p-values in the heatmap plots, then save in the g.l list the row-scaled regular and non-scaled correlograms
      if (length(fdr_hm_list) != 0) {
        gl.fdr = list()
        gl.fdr$regular = fdr_hm_list$row$regular
        gl.fdr$gene.correlogram = fdr_hm_list$none$correlogram
        gl.fdr$gene.correlogram.small = fdr_hm_list$nonesmall$correlogram.small
        gl.fdr$samp.correlogram = fdr_hm_list$none$correlogram.sample
        grid::grid.newpage()
        print(gl.fdr$regular )
        grid::grid.newpage()
        print(gl.fdr$gene.correlogram)
        grid::grid.newpage()
        print(gl.fdr$gene.correlogram.small)
        grid::grid.newpage()
        print(gl.fdr$samp.correlogram)
        g.l[["fdr"]] = gl.fdr
      }

      return(g.l)

    }
