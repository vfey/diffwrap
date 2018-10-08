
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
pheatmap_plots <-
  function(d3, id, sym.col="gene_symbol", samp.info = samp.info, samples, groups, main=NULL, p.thr=0.05, fdr.thr=0.05, logfc.thr=1)
  {
    #print(samples)
    #print(groups)
    #print(d3)
    
    #set rownames to gene_symbol names
    rownames(d3) = as.character(d3[,sym.col])
    
    # find columns with P-Values or FDR values
    pv.col = names(d3)[grep("^p\\.{0,1}val[e-u]{0,2}$", tolower(names(d3)))]
    fdr.col = names(d3)[grep("^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", tolower(names(d3)))]
    
    g.l = list()
    cat("Heatmap plots...\n")
    #cat(pv.col)
    #cat(fdr.col)
    #make two data frames for significant p-values and significant adj.p.values
    dat.sign.pv = d3[d3[[pv.col]] < 0.05,]
    dat.sign.fdr = d3[d3[[fdr.col]] < 0.05,]
    #print(dat.sign.pv)
    #print(samp.info)
    
    #select only the columns representing the samples
    samp.names = as.character(samp.info[,"SampleNames"])
    #print(samp.names)    
    
    dat.sign.pv = dat.sign.pv[samp.names]
    dat.sign.fdr = dat.sign.fdr[samp.names]
    #print(groups)
    
    #get the heatmap column annotation
    samp.anno = as.data.frame(groups)
    rownames(samp.anno) = samp.names
    colnames(samp.anno)[1] = "Sample.Class"
    
    #create lists of pheatmaps for significant p-values and adjusted p-values respectively
    pv_hm_list = list()
    fdr_hm_list = list()
    
    if (nrow(dat.sign.fdr) > 0) {
      if (nrow(dat.sign.fdr) > 100) {
        dat.sign.fdr = as.data.frame(dat.sign.fdr[1:100,])
      }
      fdr_hm_list[["row"]] = diffr_pheatmap(dat.sign.fdr, clinical.mat = samp.anno, scale.fl = "row")
      fdr_hm_list[["none"]] = diffr_pheatmap(dat.sign.fdr, clinical.mat = samp.anno, scale.fl = "none")
    } else {
      print("There are 0 entries with significant adjusted P-values in the differential expression data frame")
      print("Checking P-value entries...")
      
      if (nrow(dat.sign.pv) > 0) {
        if (nrow(dat.sign.pv) > 100) {
          dat.sign.pv = as.data.frame(dat.sign.pv[1:100,])
        }
        pv_hm_list[["row"]] = diffr_pheatmap(dat.sign.pv, clinical.mat = samp.anno, scale.fl = "row")
        pv_hm_list[["none"]] = diffr_pheatmap(dat.sign.pv, clinical.mat = samp.anno, scale.fl = "none")
      } else {
        print("There are 0 entries with significant P-values in the differential expression data frame")
      }
      
      #if there are no entries with significant adjusted p-values, but there are with significant p-values, then save in the g.l list the row-scaled regular and non-scaled correlograms
      if (length(pv_hm_list) != 0) {
        gl.pv = list()
        gl.pv$regular = pv_hm_list$row$regular
        gl.pv$gene.correlogram = pv_hm_list$none$correlogram
        gl.pv$samp.correlogram = pv_hm_list$none$correlogram.sample
        grid.newpage()
        print(gl.pv$regular )
        grid.newpage()
        print(gl.pv$gene.correlogram)
        grid.newpage()
        print(gl.pv$samp.correlogram)
        g.l[["pval"]] = gl.pv
      }
      
    }
    #if there are any entries with significant adjusted p-values in the heatmap plots, then save in the g.l list the row-scaled regular and non-scaled correlograms
    if (length(fdr_hm_list) != 0) {
      gl.fdr = list()
      gl.fdr$regular = fdr_hm_list$row$regular
      gl.fdr$gene.correlogram = fdr_hm_list$none$correlogram
      gl.fdr$samp.correlogram = fdr_hm_list$none$correlogram.sample
      grid.newpage()
      print(gl.fdr$regular )
      grid.newpage()
      print(gl.fdr$gene.correlogram)
      grid.newpage()
      print(gl.fdr$samp.correlogram)
      g.l[["fdr"]] = gl.fdr
    }
    
    return(g.l)
    
  }
