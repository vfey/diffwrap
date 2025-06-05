
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
#' Need to be in the same order as sample column names!
#' @param main \code{character}. Main plot title. (Will be complemented with additional information, e.g., 'FDR' when
#' labelling according to and FDR threshold.)
#' @param p.thr \code{numeric}. Plotted values with a P-Value below this threshold will be labelled in the P-Value plot.
#' @param fdr.thr \code{numeric}. Plotted values with a FDR below this threshold will be labelled in the FDR plot.
#' @param logfc.thr \code{numeric}. FC threshold on the log2-scale. Only used if 'split.expr' is TRUE. Values above this
#' threshold will be retained. Defaults to 1.
#' @param topn \code{numeric}. Number of top values to be plotted. Defaults to 100.
#' @param split.expr \code{logical}. Should the top up- and top down-regulated genes be displayed at equal numbers (50/50),
#' if they meet the significance threshold (regardless of the actual significance)? Defaults to \code{FALSE}.
#' @return Returns a list of pheatmap plot objects
#' @note This function is not exported as it is tailor-made for the diffwrap workflow (and called in the main wrapper function).
#'     The 'diffr_pheatmap' function used here is the main work horse and exported.
pheatmap_plots <-
    function(d3, id, sym.col="gene_symbol", samp.info = samp.info, samples, groups, sample.plot.names, main=NULL, p.thr=0.05, fdr.thr=0.05,
             logfc.thr = 1, topn = 100, split.expr = FALSE)
    {
      # elim row with duplicate sym.col entries, set rownames to gene_symbol names

      d3 = as.data.table(d3)
      d3 = as.data.frame(unique(d3, by = sym.col))
      rownames(d3) = as.character(d3[,sym.col])

      # find columns with P-Values or FDR values
      pv.col = names(d3)[grep("^p\\.{0,1}val[e-u]{0,2}$", tolower(names(d3)))]
      fdr.col = names(d3)[grep("^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", tolower(names(d3)))]
      fc.col <- names(d3)[grep("^logfc$|fold$", tolower(names(d3)))]

      g.l = list()
      #cat("Heatmap plots...\n")

      #make two data frames for significant p-values and significant adj.p.values
      dat.sign.pv = d3[d3[[pv.col]] < p.thr,]
      dat.sign.fdr = d3[d3[[fdr.col]] < fdr.thr,]

      if (split.expr) {
        cat("  Splitting expression values into 50% up- and 50% down-regulated...")
        dat.sign.pv.up <- d3[d3[[pv.col]] < p.thr & d3[[fc.col]] >= logfc.thr, ]
        dat.sign.pv.up <- dat.sign.pv.up[order(dat.sign.pv.up[[pv.col]], rev(dat.sign.pv.up[[fc.col]])), ]
        dat.sign.pv.up <- dat.sign.pv.up[1:min(nrow(dat.sign.pv.up), topn %/% 2), ]
        dat.sign.pv.down <- d3[d3[[pv.col]] < p.thr & d3[[fc.col]] <= -1*logfc.thr, ]
        dat.sign.pv.down <- dat.sign.pv.down[order(dat.sign.pv.down[[pv.col]], rev(dat.sign.pv.down[[fc.col]])), ]
        dat.sign.pv.down <- dat.sign.pv.down[1:min(nrow(dat.sign.pv.down), topn-(topn %/% 2)), ]
        dat.sign.pv <- rbind(dat.sign.pv.up, dat.sign.pv.down)
        dat.sign.fdr.up <- d3[d3[[fdr.col]] < p.thr & d3[[fc.col]] >= logfc.thr, ]
        dat.sign.fdr.up <- dat.sign.fdr.up[order(dat.sign.fdr.up[[fdr.col]], rev(dat.sign.fdr.up[[fc.col]])), ]
        dat.sign.fdr.up <- dat.sign.fdr.up[1:min(nrow(dat.sign.fdr.up), topn %/% 2), ]
        dat.sign.fdr.down <- d3[d3[[fdr.col]] < p.thr & d3[[fc.col]] <= -1*logfc.thr, ]
        dat.sign.fdr.down <- dat.sign.fdr.down[order(dat.sign.fdr.down[[fdr.col]], rev(dat.sign.fdr.down[[fc.col]])), ]
        dat.sign.fdr.down <- dat.sign.fdr.down[1:min(nrow(dat.sign.fdr.down), topn-(topn %/% 2)), ]
        dat.sign.fdr <- rbind(dat.sign.fdr.up, dat.sign.fdr.down)
      }
      cat("done\n")

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
        cat("  Creating pheatmap objects for FDR-filtered genes...\n")
        if (nrow(dat.sign.fdr) > topn) {
          cat("    Labelling", topn, "genes...\n")
          dat.sign.fdr = as.data.frame(dat.sign.fdr[1:topn,])
        }
        #make 2nd option smaller pheatmap for genes
        cat("    Creating 2nd, smaller data frame labelling only 50 genes...")
        dat.sign.fdr.small = as.data.frame(dat.sign.fdr[1:50,])
        cat("done\n")
        browser()
        fdr_hm_list[["row"]] = diffr_pheatmap(dat.sign.fdr, clinical.mat = samp.anno, scale.fl = "row", sign.val = "FDR")
        fdr_hm_list[["none"]] = diffr_pheatmap(dat.sign.fdr, clinical.mat = samp.anno, scale.fl = "none", sign.val = "FDR")
        fdr_hm_list[["rowsmall"]] = diffr_pheatmap(dat.sign.fdr.small, clinical.mat = samp.anno, sign.val = "FDR")
        fdr_hm_list[["nonesmall"]] = diffr_pheatmap(dat.sign.fdr.small, clinical.mat = samp.anno, sign.val = "FDR")
      cat(" Pheatmap objects for FDR done\n")
      } else {
        print("There are 0 entries with significant adjusted P-values in the differential expression data frame")
        print("Checking P-value entries...")

        cat("Creating pheatmap objects for P-value-filtered genes...")
        if (nrow(dat.sign.pv) > 0) {
          if (nrow(dat.sign.pv) > topn) {
            dat.sign.pv = as.data.frame(dat.sign.pv[1:topn,])
          }
          dat.sign.pv.small = as.data.frame(dat.sign.pv[1:50,])
          pv_hm_list[["row"]] = diffr_pheatmap(dat.sign.pv, clinical.mat = samp.anno, scale.fl = "row", sign.val = "P-value")
          pv_hm_list[["none"]] = diffr_pheatmap(dat.sign.pv, clinical.mat = samp.anno, scale.fl = "none", sign.val = "P-value")
          pv_hm_list[["rowsmall"]] = diffr_pheatmap(dat.sign.pv.small, clinical.mat = samp.anno, scale.fl = "row", sign.val = "P-value")
          pv_hm_list[["nonesmall"]] = diffr_pheatmap(dat.sign.pv.small, clinical.mat = samp.anno, scale.fl = "none", sign.val = "P-value")
          cat("done\n")
        } else {
          print("There are 0 entries with significant P-values in the differential expression data frame")
        }

        #if there are no entries with significant adjusted p-values, but there are with significant p-values, then save in the g.l list the row-scaled regular and non-scaled correlograms
        if (length(pv_hm_list) != 0) {
          cat("  Plotting P-value-filtered genes...")
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
          cat("done\n")
        }

      }
      #if there are any entries with significant adjusted p-values in the heatmap plots, then save in the g.l list the row-scaled regular and non-scaled correlograms
      if (length(fdr_hm_list) != 0) {
        cat("  Plotting FDR-filtered genes...")
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
        cat("done\n")
      }

      return(g.l)

    }
