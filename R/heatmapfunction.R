required_packages = c('ggplot2','pheatmap','RColorBrewer','kableExtra','viridis','Hmisc','ltm')
for(p in required_packages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

c

data_heatmap = read.table("./heatmap_data.tsv",
                          stringsAsFactors = F,
                          header = T)
data_heatmap = as.data.frame(data_heatmap)

sample_class =  read.table("./sample_class.tsv",
                           stringsAsFactors = F,
                           header = T)
sample_class = as.data.frame(sample_class)


my_pheatmap = diffr_pheatmap(expr.mat = data_heatmap, clinical.mat = sample_class)



#calculate the correlation matrix using Hmisc package
t_data_heat_map = t(data_heatmap)
corr_hmap = Hmisc::rcorr(t_data_heat_map) #round(cor(data_heatmap), 2)
cor.mat = as.matrix(corr_hmap$r)
pv.mat = as.matrix(corr_hmap$P)

cor.mat.breaks <- seq(min(cor.mat), max(cor.mat),by = 0.05)
colour = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(cor.mat.breaks) + 1)

#define new clustering distance
dissimilarity <- 1 - cor.mat
dist.clust <- as.dist(dissimilarity)


# define significance levels
sign.stars <- ifelse(pv.mat < .001, "***", ifelse(pv.mat < .01, "** ", ifelse(pv.mat < .05, "* ", " ")))
sign.stars <- ifelse(is.na(sign.stars),"",sign.stars)


pheatmap::pheatmap(as.data.frame(corr_hmap$r),
                   show_colnames = T, show_rownames = T,
                   cluster_rows = T, cluster_cols = T,
                   legend = T, border_color = "white",
                   scale = "none", color = colour,
                   clustering_distance_rows = dist.clust, clustering_distance_cols = dist.clust,
                   display_numbers = sign.stars, number_color = "#FFFF00",
                   treeheight_row = 100, treheight_col = 100 )
