library(pheatmap)

a <- matrix(rnorm(90), 9, 10)
txt <- matrix(sample(letters, replace = TRUE, 90), 9, 10)
pheatmap(a, cluster_rows = F, cluster_cols = F, display_numbers = txt)


corstarsl <- function(x){
  require(Hmisc)
  x <- as.matrix(x)
  R <- rcorr(x)$r
  p <- rcorr(x)$P

  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))

  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]

  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")

  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew)

  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew)
}

##Create table _insert your dataframe below
New_table<-corstarsl(yourdataframe)
#=====================================================================================================================
data_heatmap = read.table("./heatmap_data.tsv",
                          stringsAsFactors = F,
                          header = T)
data_heatmap = as.data.frame(data_heatmap)

sample_class =  read.table("./sample_class.tsv",
                           stringsAsFactors = F,
                           header = T)
sample_class = as.data.frame(sample_class)

t_data_heat_map = t(data_heatmap)

corr_hmap = Hmisc::rcorr(t_data_heat_map) #round(cor(data_heatmap), 2)
cor.mat = as.matrix(corr_hmap$r)
pv.mat = as.matrix(corr_hmap$P)

## define notions for significance levels; spacing is important.
sign.stars <- ifelse(pv.mat < .001, "***", ifelse(pv.mat < .01, "** ", ifelse(pv.mat < .05, "* ", " ")))
sign.stars <- ifelse(is.na(sign.stars),"",sign.stars)


my_pheatmap = diffr_pheatmap(data_heatmap,sample_class)

