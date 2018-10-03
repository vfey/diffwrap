
# Function: heatmap_plots
#
# Author: Bogdan Iancu - Genevia Technologies Oy
#
# Arguments:
#         d3 = differentia expression matrix in regular  (genes, samples) format 
#         

# Output: list of pheatmap objects that should be preferably handled with lapply
#
# Details: calls function diffr_pheatmap 
#
pheatmap_plots <-
  function(d3, id, sym.col="gene_symbol", samp.info = samp.info, samples, groups, main=NULL, p.thr=0.05, fdr.thr=0.05, logfc.thr=1)
  {
    
    # initial checks
    if ((missing(id) || !id %in% names(d3)) && (!sym.col %in% names(d3))) {
      stop("Need at least one ID column, i.e., one of 'id' or 'sym.col'.")
    }
    ## if ID column exists and symbol column is missing and different from ID column, create symbol column from ID column
    if (!missing(id) && id %in% names(d3) && !sym.col %in% names(d3) && sym.col != id) {
      d3[[sym.col]] <- id
    }
    
    #data for testing reasons
    samp.info = read.table("sampInfo_MAAKE_final.tsv",sep="\t", header=T, stringsAsFactor=F)
    d3 = read.table("MAAKE_contrast_Ctrl-OAB_dry_differential_expression.tsv", sep = "\t", header = TRUE)
    colnames(d3)[4:33] = gsub("X", "", colnames(d3)[4:33])
    sym.col="gene_symbol"
    samp.info$SampleName = gsub("-", ".", samp.info$SampleName)
    samples = "SampleName"
    groups = "Treatment"
    
    #set rownames to gene_symbol names
    rownames(d3) = as.character(d3[,sym.col])
   
     # find columns with P-Values or FDR values
    pv.col = names(d3)[grep("^p\\.{0,1}val[e-u]{0,2}$", tolower(names(d3)))]
    fdr.col = names(d3)[grep("^fdr$|^adj*\\.{0,1}p\\.{0,1}val[e-u]{0,2}$", tolower(names(d3)))]
    
    g.l = list()
    cat("Heatmap plots...\n")
    
    #make two data frames for signigicant p-values and significant adj.p.values
    dat.sign.pv = d3[d3[[pv.col]] < 0.05,]
    dat.sign.fdr = d3[d3[[fdr.col]] < 0.05,]
    
    #select only the columns representing the samples 
    dat.sign.pv = dat.sign.pv[samp.info[,samples]]
     
    dat.sign.fdr = dat.sign.fdr[samp.info[,samples]]
    
    #get the heatmap column annotation
    samp.anno = as.data.frame(samp.info[,groups])
    rownames(samp.anno) = samp.info[,samples]
    colnames(samp.anno)[1] = groups
    
    #create lists of pheatmaps for significant p-values and adjusted p-values respectively
    pv_hm_list = list()
    fdr_hm_list = list()
    
    if ((nrow(dat.sign.fdr) > 0) && (nrow(dat.sign.fdr) <= 100)) {
      fdr_hm_list[["row"]] = diffr_pheatmap(dat.sign.fdr, clinical.mat = samp.anno, scale.fl = "row")
      fdr_hm_list[["none"]] = diffr_pheatmap(dat.sign.fdr, clinical.mat = samp.anno, scale.fl = "none")
      
    } else if (nrow(dat.sign.fdr) > 100) {
      dat.sign.fdr = dat.sign.fdr[1:100,]
      fdr_hm_list[["row"]] = diffr_pheatmap(dat.sign.fdr, clinical.mat = samp.anno, scale.fl = "row")
      fdr_hm_list[["none"]] = diffr_pheatmap(dat.sign.fdr, clinical.mat = samp.anno, scale.fl = "none")
      
    } else {
      print("There are 0 entries with significant adjusted P-values in the differential expression data frame")
      print("Checking P-value entries...")
      
      
      if ((nrow(dat.sign.pv) > 0) && (nrow(dat.sign.pv) <= 100)) {
        pv_hm_list[["row"]] = diffr_pheatmap(dat.sign.pv, clinical.mat = samp.anno, scale.fl = "row")
        pv_hm_list[["none"]] = diffr_pheatmap(dat.sign.pv, clinical.mat = samp.anno, scale.fl = "none")
      } else if (nrow(dat.sign.pv) > 100) {
        dat.sign.pv = as.data.frame(dat.sign.pv[1:100,])
        pv_hm_list[["row"]] = diffr_pheatmap(dat.sign.pv, clinical.mat = samp.anno, scale.fl = "row")
        pv_hm_list[["none"]] = diffr_pheatmap(dat.sign.pv, clinical.mat = samp.anno, scale.fl = "none")
      } else{
        print("There are 0 entries with significant P-values in the differential expression data frame")
      }
      
      #if there are any noentries with significant adjusted p-values, but there are with significant p-values, then save in the g.l list the row-scaled regular and non-scaled correlograms
      if (length(pv_hm_list) != 0) {
        gl.pv = list()
        gl.pv$regular = pv_hm_list$row$regular
        gl.pv$gene.correlogram = pv_hm_list$none$correlogram
        g.l[["pval"]] = gl.pv
      }
      
      
    }
    
    #if there are any entries with significant adjusted p-values in the heatmap plots, then save in the g.l list the row-scaled regular and non-scaled correlograms
    if (length(fdr_hm_list) != 0) {
      gl.fdr = list()
      gl.fdr$regular = fdr_hm_list$row$regular
      gl.fdr$gene.correlogram = fdr_hm_list$none$correlogram
      g.l[["fdr"]] = gl.fdr
    }
    
    
    return(g.l)
    
  }


