###########################################################################################################################
# Author: WTC
# Date: July 2018
# Function: Sample 1000x 251 cells from 10x cytotoxic cell dataset, calculate DE with CD dataset and determine risk genes in DE genes
###########################################################################################################################
#
# Libraries
#
###########################################################################################################################
library(Seurat)
library(Matrix)
#library(Matrix.utils)
#library(ggplot2)
#library(pryr)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################

# download Cytotoxic T cell gene/cell matrix (raw) from https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/cytotoxic_t
#CTL_10x<-Read10X(data.dir = "~/Downloads/matrices_mex 3/hg19/")
#CTL_10x<-CreateSeuratObject(raw.data = CTL_10x, min.cells = 3, project = "healthy")
# select subset of 251 cells
#load("~/Desktop/Single_cell/CTL_10x_CD/CTL_10x.Rda")

setwd("~/randomised_permutation_10xsubset_cd/")
load("~/Desktop/Single_cell/CTL_10x_CD/CTL_10x.Rda")
load("~/Desktop/Single_cell/CTL_10x_CD/CTL_CD.Rda")
risk_use<-read.csv("~/Desktop/Single_cell/final_data_paper/DE/risk_genes_in_cd_and_10x.csv")
genes_use<-read.csv("~/Desktop/Single_cell/final_data_paper/DE/genes_in_cd_and_10x.csv")
mkdir("~/nonrandomised_permutation_10xsubset_cd/")
setwd("~/nonrandomised_permutation_10xsubset_cd/")
a<-"number_DEgenes"
b<-"number_DEriskgenes"
number_DE_genes<-data.frame(a)
number_DE_riskgenes<-data.frame(b)
for (i in 1:1000){
cat("iteration", i, "\n")
samp_col_idxs  <- sample(ncol(CTL_10x@data), 251)
samp_col_names <- colnames(CTL_10x@data) [samp_col_idxs]
CTL_10x_subset<-CTL_10x@data[,samp_col_names]
dim(CTL_10x_subset)
CTL_10x_subset<-CreateSeuratObject(raw.data = CTL_10x_subset, min.cells = 3, project = "subset")
mito.genes <- grep(pattern = "^MT-", x = rownames(x = CTL_10x_subset@data), value = TRUE)
percent.mito <- Matrix::colSums(CTL_10x_subset@raw.data[mito.genes, ])/Matrix::colSums(CTL_10x_subset@raw.data)
CTL_10x_subset <- AddMetaData(object = CTL_10x_subset, metadata = percent.mito, col.name = "percent.mito")

# for ngene: 200-2500, percent mito: <5%
CTL_10x_subset <- FilterCells(object = CTL_10x_subset, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

# log normalize
CTL_10x_subset <- NormalizeData(object = CTL_10x_subset, normalization.method = "LogNormalize", scale.factor = 10000)
CTL_10x_subset <- ScaleData(object = CTL_10x_subset, vars.to.regress = c("nUMI", "percent.mito"))
CTL_10x_subset@meta.data$dataset <- "subset10x"

save(CTL_10x_subset, file=paste0("CTL_10x_subset",i,".Rda"))

## assign name to CD dataset
CTL_CD@meta.data$dataset <- "cd"

## Combine the cdulated and ctrlulated cells into a single object
CTL_subset.combined <- RunCCA(CTL_10x_subset, CTL_CD, genes.use = intersect(rownames(CTL_10x_subset@data), rownames(CTL_CD@data)), num.cc = 30, scale.data=T)

###############################################################
# Process the data after combining the two datasets
###############################################################
CTL_subset.combined<-SetAllIdent(CTL_subset.combined, "dataset")
subset10x_cells<-WhichCells(CTL_subset.combined, "subset10x")
cd_cells<-WhichCells(CTL_subset.combined, "cd")

# calculate DE
# DE with genes present in the two datasets, in >1%, MAST
CTL_subset.combined_DE_markers_1perc = FindAllMarkers(CTL_subset.combined, only.pos = T, test.use = "MAST", genes.use = genes_use$gene)
write.table(CTL_subset.combined_DE_markers_1perc, file = paste0("CTL_subset",i,".blood.10x._DEmarkers_MAST.txt"))

# and the significant ones (p<0.05)
CTL_subset.combined_DE_markers_1perc<-CTL_subset.combined_DE_markers_1perc[CTL_subset.combined_DE_markers_1perc$p_val_adj <0.05,]
number_DE_genes[i,a]<-nrow(CTL_subset.combined_DE_markers_1perc)

# merge all DE results with risk genes
risk_all_filtered_subset<-merge(risk_use, CTL_subset.combined_DE_markers_1perc, by="gene", all=F)
write.csv(risk_all_filtered_subset, paste0("riskgenes_subset",i,"_vs_cd_filtered.csv"))
number_DE_riskgenes[i,b]<-nrow(risk_all_filtered_subset)}

write.csv(number_DE_genes, "number_DE_genes_cd_vs_10xsubsets.csv")
write.csv(number_DE_riskgenes, "number_DE_riskgenes_cd_vs_10xsubsets.csv")
