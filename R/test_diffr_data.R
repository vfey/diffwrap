test_data_location = "/data/diffr_test_data"

# NOTE: Set this to your diffr git directory and checkout branch bogdan
setwd("~/Gitlab/diffr/")

# Load libraries
library(readxl)
library(readR)
library(medseqr)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)
library(scales)
library(RColorBrewer)
library(scatterplot3d)
library(dendextend)
library(org.Hs.eg.db)

# Load preprocessing and dE scripts
source("./R/counts_read_and_filter.R")
source("./R/de_plots.R")
source("./R/de_util.R")
source("./R/diffExpr.R")
source("./R/preprocess.R")
source("./R/qc_plots.R")
source("./R/qc_util.R")
source("./R/heatmapfunction.R")
source("./R/heatmap_plots.R")

# Load enrichment analysis scripts
# source("./R/run_clusterProfiler_GO.R")
# source("./R/run_clusterProfiler_KEGG.R")
# source("./R/doDavidEnrichmentAnalysis.R")
# source("./R/run_gProfileR.R")
# source("./R/enrichment_util.R")
# source("./R/run_topGo.R")

# Set output directory
out.dir="/data/diffr_test_results"
dir.create(out.dir, showWarnings = F)

# Read in and preprocess files
samp.info = read_excel(file.path(test_data_location, "Kim_dataset_sample_sheet_for_Diffr.xlsx"))
samp.info$SampleName <- sub("(.+)\\.tab$", "\\1", basename(samp.info$countFile))
samp.info$Progression_step = gsub(' ', '.', samp.info$Progression_step)
samp.info$countFile = gsub("/Users/meeripekkarinen/Example_RNAseq_project/STAR_files/",
                           paste0(test_data_location, '/'),
                           samp.info$countFile)
f = factor(samp.info$Progression_step, levels =c("BPH", "Localized.PC", "Advanced.PC", "CRPC"))
design1 <- model.matrix(~0+f)
colnames(design1)=c("BPH", "Localized.PC", "Advanced.PC", "CRPC")
#at this point samp.info is tibble, make it data.frame
samp.info = as.data.frame(samp.info)
out <- diffExpr(expr.file = samp.info$countFile, 
                control = "BPH", 
                design=design1,
                contrasts=c("Localized.PC-BPH", "Advanced.PC-Localized.PC", "CRPC-Advanced.PC"),
                samp.info = samp.info, 
                do.voom = TRUE,
                strict=FALSE,
                samples = "SampleName", 
                sample.plot.names = "Sample.name", 
                groups = "Progression_step", 
                block = FALSE, 
                out.dir = out.dir, 
                analysis.name = "test_run", 
                biomart = TRUE, norm.method = "tmm", 
                biom.data.set = "human", 
                biom.mart = "ensembl", 
                # host="http://plants.ensembl.org", 
                biom.attributes=c("ensembl_gene_id","hgnc_symbol","description", "entrezgene"), 
                # biom.filter = "ensembl_gene_id", 
                #sym.col="hgnc_symbol", 
                p.thr=0.05, 
                fdr.thr=0.05, 
                logfc.thr=1,
                numlab=15, 
                n=1000)

