
library(data.table)
library(reshape2)
library(readxl)
options(stringsAsFactors = FALSE)

#load the reference
reference <- read_xlsx("/Applications/Results/Article/6 main scRNAseq studies in human IBD/Smillie-2019-Intra--and-inter-cellular-rewiring-/Supplementary Table/Table S2.xlsx", sheet = 1, col_names = T)
reference <- reference[,1:2]

#grab the marker genes for each cell tyeps
markergenes <- dcast(setDT(reference), ident~rowid(ident, prefix="gene"), value.var="gene")
markergenes <- as.data.frame(t(markergenes))
colnames(markergenes) <- markergenes[1,]
markergenes <- markergenes[-1,]
dim(markergenes)
rownames(markergenes) <- paste0("gene", 1:202)
rm(reference)
markergenes[is.na(markergenes)] <- "+++"


#Load input_clusters (markergenes of your clusters)
input <- fread("/Applications/Results/Vedo2_batch1_biopsy_dataset/subclusters/epithelial/celltyping/markers_epi_log.csv",data.table = T)
input <- input[,-c(1:6)]
clusters <- dcast(setDT(input), cluster~rowid(cluster, prefix="gene"), value.var="gene")
clusters <- as.data.frame(t(clusters))
colnames(clusters) <- paste0("cluster", 0:13)
clusters <- clusters[-1,]
dim(clusters)
rownames(clusters) <- paste0("gene", 1:993)
rm(input)
clusters[is.na(clusters)] <- "---"


#matching genes
output=as.data.frame(matrix(nrow = ncol(clusters),ncol = ncol(markergenes)))
for(i in 1:ncol(clusters)){
  tmp.cluster=colnames(clusters)[i]
  tmp.cluster.data=clusters[,tmp.cluster]
  
  mm=c()
  for(j in 1:ncol(markergenes)){
    tmp.markergenes=colnames(markergenes)[j]
    tmp.markergenes.data=markergenes[,tmp.markergenes]
    aa=length(which(tmp.cluster.data %in% tmp.markergenes.data))
    mm=append(mm,aa)
  }
  output[i,]=mm
}
colnames(output)=colnames(markergenes)
rownames(output)=colnames(clusters)

write.csv(output, "/Applications/Results/Vedo2_batch1_biopsy_dataset/subclusters/epithelial/celltyping/output/Smillie_matching_markergenes.csv")

