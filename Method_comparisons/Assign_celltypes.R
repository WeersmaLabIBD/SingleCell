## scPred on Sanger datasets

library(Seurat)
library(scPred)
library(tidyverse)
library(SingleCellExperiment)

HC_Sanger1<-readRDS("~/dataset_sct.rds")

load("~/scp_smillie_750_30pcs_categories.Rd")
load("~/res_smillie_750_30pcs_categories.Rd")

dataset2<-as.SingleCellExperiment(HC_Sanger1, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_Sanger1@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$seurat_clusters, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "1"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "2"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "4"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "5"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "6"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "7"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "8"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "9"]<-"fibroblast"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "10"]<-"Bcell"
colnames(CellsMeta)[12]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1, 12,13,14)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_Sanger1<-AddMetaData(HC_Sanger1, CellsMeta)
DimPlot(HC_Sanger1, group.by="Predicted_all")

saveRDS(HC_Sanger1, "~/dataset_sct.rds")



######

HC_Sanger2<-readRDS("~/dataset_sct.rds")


dataset2<-as.SingleCellExperiment(HC_Sanger2, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_Sanger2@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$seurat_clusters, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "1"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "2"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "4"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "5"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "6"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "7"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "8"]<-"epithelium"
colnames(CellsMeta)[12]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1,12,13,14)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_Sanger2<-AddMetaData(HC_Sanger2, CellsMeta)
DimPlot(HC_Sanger2, group.by="Predicted_all") # save at size of 5 inch x 10 inch

saveRDS(HC_Sanger2, "~/dataset_sct.rds")

######


HC_Sanger3<-readRDS("~/dataset_sct.rds")


dataset2<-as.SingleCellExperiment(HC_Sanger3, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_Sanger3@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$seurat_clusters, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "1"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "2"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "4"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "5"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "6"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "7"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "8"]<-"epithelium"
colnames(CellsMeta)[12]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1,12,13,14)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_Sanger3<-AddMetaData(HC_Sanger3, CellsMeta)
DimPlot(HC_Sanger3, group.by="Predicted_all") # save at size of 5 inch x 10 inch

saveRDS(HC_Sanger3, "~/dataset_sct.rds")

#######
HC_Sanger4<-readRDS("~/dataset_sct.rds")


dataset2<-as.SingleCellExperiment(HC_Sanger4, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_Sanger4@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$seurat_clusters, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "1"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "2"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "4"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "5"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "6"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "7"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "8"]<-"Tcell"
colnames(CellsMeta)[12]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1,12,13,14)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_Sanger4<-AddMetaData(HC_Sanger4, CellsMeta)
DimPlot(HC_Sanger4, group.by="Predicted_all") # save at size of 5 inch x 10 inch
saveRDS(HC_Sanger4, "~/dataset_sct.rds")

#######

HC_Sanger5<-readRDS("~/dataset_sct.rds")

dataset2<-as.SingleCellExperiment(HC_Sanger5, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_Sanger5@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$seurat_clusters, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "1"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "2"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "4"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "5"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "6"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "7"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "8"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "8"]<-"myeloid"

colnames(CellsMeta)[12]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1,12,13,14)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_Sanger5<-AddMetaData(HC_Sanger5, CellsMeta)
DimPlot(HC_Sanger5, group.by="Predicted_all") # save at size of 5 inch x 10 inch
saveRDS(HC_Sanger5, "~/dataset_sct.rds")

#######
HC_Sanger6<-readRDS("~/dataset_sct.rds")
dataset2<-as.SingleCellExperiment(HC_Sanger6, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_Sanger6@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$seurat_clusters, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "1"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "2"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "3"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "4"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "5"]<-"epithelium"


colnames(CellsMeta)[12]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1,12,13,14)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_Sanger6<-AddMetaData(HC_Sanger6, CellsMeta)
DimPlot(HC_Sanger6, group.by="Predicted_all") # save at size of 5 inch x 10 inch
saveRDS(HC_Sanger6, "~/dataset_sct.rds")

######
HC_Sanger7<-readRDS("~/dataset_sct.rds")

dataset2<-as.SingleCellExperiment(HC_Sanger7, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_Sanger7@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$seurat_clusters, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "1"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "2"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "4"]<-"fibroblast"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "5"]<-"myeloid"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "6"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "7"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "8"]<-"Tcell"

colnames(CellsMeta)[12]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1,12,13,14)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_Sanger7<-AddMetaData(HC_Sanger7, CellsMeta)
DimPlot(HC_Sanger7, group.by="Predicted_all") # save at size of 5 inch x 10 inch
saveRDS(HC_Sanger7, "~/dataset_sct.rds")

#######
HC_Sanger8<-readRDS("~/dataset_sct.rds")
dataset2<-as.SingleCellExperiment(HC_Sanger8, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_Sanger8@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$SCT_snn_res.0.1, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "1"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "2"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "4"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "5"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "6"]<-"epithelium"

colnames(CellsMeta)[18]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1,16,17,18)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_Sanger8<-AddMetaData(HC_Sanger8, CellsMeta)
DimPlot(HC_Sanger8, group.by="Predicted_all") # save at size of 5 inch x 10 inch
saveRDS(HC_Sanger8, "~/Ddataset_sct.rds")

#######

HC_Sanger9<-readRDS("~/dataset_sct.rds")
dataset2<-as.SingleCellExperiment(HC_Sanger9, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_Sanger9@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$SCT_snn_res.0.3, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.3 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.3 == "1"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.3 == "2"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.3 == "3"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.3 == "4"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.3 == "5"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.3 == "6"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.3 == "7"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.3 == "8"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.3 == "9"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.3 == "10"]<-"myeloid"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.3 == "11"]<-"Tcell"

colnames(CellsMeta)[13]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1,13,14,15)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_Sanger9<-AddMetaData(HC_Sanger9, CellsMeta)
DimPlot(HC_Sanger9, group.by="Predicted_all") # save at size of 5 inch x 10 inch
saveRDS(HC_Sanger9, "~/dataset_sct.rds")


## check CD3E expression on sanger files
DimPlot(HC_Sanger1, group.by="Predicted_all")
FeaturePlot(HC_Sanger1, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 

DimPlot(HC_Sanger2, group.by="Predicted_all")
FeaturePlot(HC_Sanger2, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 

DimPlot(HC_Sanger3, group.by="Predicted_all")
FeaturePlot(HC_Sanger3, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 

DimPlot(HC_Sanger4, group.by="Predicted_all")
FeaturePlot(HC_Sanger4, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 

DimPlot(HC_Sanger5, group.by="Predicted_all")
FeaturePlot(HC_Sanger5, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 

DimPlot(HC_Sanger6, group.by="Predicted_all")
FeaturePlot(HC_Sanger6, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 

DimPlot(HC_Sanger7, group.by="Predicted_all")
FeaturePlot(HC_Sanger7, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 

DimPlot(HC_Sanger8, group.by="Predicted_all")
FeaturePlot(HC_Sanger8, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 

DimPlot(HC_Sanger9, group.by="Predicted_all")
FeaturePlot(HC_Sanger9, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 


x<-data.frame(table(HC_Sanger1@meta.data$Predicted_all))
x$HC_Sanger2<-table(HC_Sanger2@meta.data$Predicted_all)
x$HC_Sanger3<-table(HC_Sanger3@meta.data$Predicted_all)
x$HC_Sanger4<-table(HC_Sanger4@meta.data$Predicted_all)
x$HC_Sanger5<-table(HC_Sanger5@meta.data$Predicted_all)
x$HC_Sanger6<-table(HC_Sanger6@meta.data$Predicted_all)
x[1,10]<-94
x[2,10]<-40
x[3,10]<-4447
x[4,10]<-149
x[5,10]<-132
x[6,10]<-93
x$HC_Sanger8<-table(HC_Sanger8@meta.data$Predicted_all)
x$HC_Sanger9<-table(HC_Sanger9@meta.data$Predicted_all)
write.csv(x, "~/cellcounts_Snager_data.csv")

#####
..<-readRDS("~/Desktop/PSC/Day2/sequencing/lane1/")
HC_NI_3129<-readRDS("~/HC-NI-3129_sct.rds")
HC_NI_3002<-readRDS("~/HC-NI-3002_sct.rds")
..<-readRDS("~/")
..<-readRDS("~/")

load("~/scp_smillie_750_30pcs_categories.Rd")
load("~/res_smillie_750_30pcs_categories.Rd")

dataset2<-as.SingleCellExperiment(HC_NI_3129, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_NI_3129@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$SCT_snn_res.0.5, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "1"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "2"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "4"]<-"fibroblast"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "5"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "6"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "7"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "8"]<-"endothelial"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "9"]<-"myeloid"

colnames(CellsMeta)[22]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1, 22,23,24)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_NI_3129<-AddMetaData(HC_NI_3129, CellsMeta)
DimPlot(HC_NI_3129, group.by="Predicted_all")

saveRDS(HC_NI_3129, "~/HC-NI-3129_sct.rds")

####
HC_NI_3002<-readRDS("~/Desktop/PSC/Day4/sequencing/lane1/HC-NI-3002_sct.rds")
dataset2<-as.SingleCellExperiment(HC_NI_3002, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_NI_3002@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$SCT_snn_res.1, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1 == "1"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1 == "2"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1 == "4"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1 == "5"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1 == "6"]<-"myeloid"

colnames(CellsMeta)[26]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1, 26,27,28)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_NI_3002<-AddMetaData(HC_NI_3002, CellsMeta)
DimPlot(HC_NI_3002, group.by="Predicted_all")

saveRDS(HC_NI_3002, "~/HC-NI-3002_sct.rds")


#####
####
HC_NI_3030<-readRDS("~/HC-NI-3030_sct.rds")
HC_NI_3030<-FindClusters(HC_NI_3030, resolution = 0.7)
dataset2<-as.SingleCellExperiment(HC_NI_3030, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_NI_3030@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$SCT_snn_res.0.5, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "0"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "1"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "2"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "4"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "5"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "6"]<-"fibroblast"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "7"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "8"]<-"endothelial"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "9"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.5 == "10"]<-"myeloid"


colnames(CellsMeta)[20]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1, 20,21,22)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_NI_3030<-AddMetaData(HC_NI_3030, CellsMeta)
DimPlot(HC_NI_3030, group.by="Predicted_all")

saveRDS(HC_NI_3030, "~/HC-NI-3030_sct.rds")

#####
####
HC_NI_3083<-readRDS("~/HC-NI-3083_sct.rds")
HC_NI_3083<-FindClusters(HC_NI_3083, resolution=0.8)
dataset2<-as.SingleCellExperiment(HC_NI_3083, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_NI_3083@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$SCT_snn_res.0.8, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "1"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "2"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "3"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "4"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "5"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "6"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "7"]<-"fibroblast"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "8"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "9"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "10"]<-"myeloid"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "11"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.8 == "12"]<-"endothelial"


colnames(CellsMeta)[24]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype


for (i in 1:nrow(CellsMeta)){ 
  if (CellsMeta[i,24] == "unassigned"){
    CellsMeta[i,23]<-CellsMeta[i,21] 
  } else {
    CellsMeta[i,23]<-CellsMeta[i,24]
  }
}



CellsMeta<-CellsMeta[c(1, 21,23,24)]
rownames(CellsMeta)<-CellsMeta$NAME

nr3083<-AddMetaData(nr3083, CellsMeta)
DimPlot(nr3083, group.by="Predicted_all")
table(nr3083@meta.data$SCT_snn_res.0.8, nr3083@meta.data$Predicted_all)
saveRDS(nr3083, "~/HC-NI-3083_sct.rds")

######
####
HC_NI_3049<-readRDS("~/HC-NI-3049_sct.rds")
HC_NI_3049<-FindClusters(HC_NI_3049, resolution = 0.7)
dataset2<-as.SingleCellExperiment(HC_NI_3049, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_NI_3049@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$SCT_snn_res.0.7, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.7 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.7 == "1"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.7 == "2"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.7 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.7 == "4"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.7 == "5"]<-"fibroblast"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.7 == "6"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.7 == "7"]<-"endothelial"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.7 == "8"]<-"myeloid"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.7 == "9"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.7 == "10"]<-"epithelium"


colnames(CellsMeta)[19]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1, 19,20,21)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_NI_3049<-AddMetaData(HC_NI_3049, CellsMeta)
DimPlot(HC_NI_3049, group.by="Predicted_all")

saveRDS(HC_NI_3049, "~/HC-NI-3049_sct.rds")


#######
####
HC_NI_3037<-readRDS("~/Desktop/CITEseq/D-3037-HC_sct.rds")
HC_NI_3037<-FindClusters(HC_NI_3037, resolution = 1.5)
dataset2<-as.SingleCellExperiment(HC_NI_3037, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_NI_3037@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$SCT_snn_res.1.5, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "0"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "1"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "2"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "4"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "5"]<-"endothelial"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "6"]<-"Bcell"

colnames(CellsMeta)[21]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1, 21,22,23)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_NI_3037<-AddMetaData(HC_NI_3037, CellsMeta)
DimPlot(HC_NI_3037, group.by="Predicted_all")

saveRDS(HC_NI_3037, "~/D-3037-HC_sct.rds")

######
HC_NI_3034<-readRDS("~/Desktop/CITEseq/C-3034-HC_sct.rds")
HC_NI_3034<-FindClusters(HC_NI_3034, resolution = 1.5)
dataset2<-as.SingleCellExperiment(HC_NI_3034, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_NI_3034@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$SCT_snn_res.1.5, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "0"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "1"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "2"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "4"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "5"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "6"]<-"epithelium"

colnames(CellsMeta)[23]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1, 23,24,25)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_NI_3034<-AddMetaData(HC_NI_3034, CellsMeta)
DimPlot(HC_NI_3034, group.by="Predicted_all")

saveRDS(HC_NI_3034, "~/C-3034-HC_sct.rds")



######
HC_NI_3296<-readRDS("~/HC-NI-3296_sct.rds")
HC_NI_3296<-FindClusters(HC_NI_3296, resolution = 1.5)
dataset2<-as.SingleCellExperiment(HC_NI_3296, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-HC_NI_3296@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(9,10)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$SCT_snn_res.1.5, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "0"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "1"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "2"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "4"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "5"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "6"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "7"]<-"fibroblast"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "8"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "9"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "10"]<-"myeloid"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "11"]<-"endothelial"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.1.5 == "12"]<-"Bcell"


colnames(CellsMeta)[21]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1, 21,22,23)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_NI_3296<-AddMetaData(HC_NI_3296, CellsMeta)
DimPlot(HC_NI_3296, group.by="Predicted_all")

saveRDS(HC_NI_3296, "~/HC-NI-3296_sct.rds")



######
DimPlot(HC_NI_3002, group.by="Predicted_all")
FeaturePlot(HC_NI_3002, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 
DimPlot(HC_NI_3129, group.by="Predicted_all")
FeaturePlot(HC_NI_3129, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 
DimPlot(HC_NI_3030, group.by="Predicted_all")
FeaturePlot(HC_NI_3030, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 
DimPlot(HC_NI_3083, group.by="Predicted_all")
FeaturePlot(HC_NI_3083, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 
DimPlot(HC_NI_3049, group.by="Predicted_all")
FeaturePlot(HC_NI_3049, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 
DimPlot(HC_NI_3034, group.by="Predicted_all")
FeaturePlot(HC_NI_3034, c("CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 
DimPlot(HC_NI_3037, group.by="Predicted_all")
FeaturePlot(HC_NI_3037, c("SELL", "CD19", "CD3E", "CD8A", "MS4A1", "EPCAM", "THY1", "MADCAM1"), pt.size = 0.01, max.cutoff = "q99") 


smillie_meta<-read.table("~/all.meta2.txt", header=T, sep="\t")
table(smillie_meta$Cluster)



#######
# T cell predcitions - integrate all data first and then extract T cells, prepare data and predict in bulk




load("~/res_smillie_allT_50pcs_clusters.Rd")
load("~/scp_smillie_allT_50pcs_clusters.Rd")


dataset2<-as.SingleCellExperiment(T_Sanger1, assay="RNA")
dataset_metadata <- as.data.frame(colData(dataset2))
dataset_counts <- counts(dataset2)
set.seed(1234)
dataset_counts_cpm  <- apply(dataset_counts, 2, function(x) (x/sum(x))*1000000)
prediction <- scPredict(scp, newData = dataset_counts_cpm, threshold = 0.7)

x<-getPredictions(prediction)
x$NAME<-row.names(x)
CellsMeta<-T_Sanger1@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(7,8)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
table(CellsMeta$SCT_snn_res.0.1, CellsMeta$predClass)
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "0"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "1"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "2"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "3"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "4"]<-"epithelium"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "5"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "6"]<-"Bcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "7"]<-"Tcell"
CellsMeta$Predicted_clusters[CellsMeta$SCT_snn_res.0.1 == "8"]<-"epithelium"
colnames(CellsMeta)[12]<-"Predicted_celltype"
CellsMeta$Predicted_all<-CellsMeta$Predicted_celltype
CellsMeta$Predicted_all[CellsMeta$Predicted_all == "unassigned"]<-CellsMeta$Predicted_clusters

CellsMeta<-CellsMeta[c(1,12,13,14)]
rownames(CellsMeta)<-CellsMeta$NAME

HC_Sanger2<-AddMetaData(HC_Sanger2, CellsMeta)
DimPlot(HC_Sanger2, group.by="Predicted_all") # save at size of 5 inch x 10 inch

saveRDS(HC_Sanger2, "~/dataset_sct.rds")

sanger_celltypes<-read.csv("~/cellcounts_Snager_data.csv", sep=";", dec=",")
sanger_50_50_celltypes<-sanger_celltypes[sanger_celltypes$sample == "HC_Sanger2",]
x<-sanger_celltypes[sanger_celltypes$sample == "HC_Sanger6",]
y<-sanger_celltypes[sanger_celltypes$sample == "HC_Sanger9",]
sanger_50_50_celltypes<-rbind(sanger_50_50_celltypes, x)
sanger_50_50_celltypes<-rbind(sanger_50_50_celltypes, y)
sanger_celltypes$mix_or_epi<-"epithelial_only"
sanger_celltypes$mix_or_epi[sanger_celltypes$sample == "HC_Sanger2"]<-"mix_epi_immune"
sanger_celltypes$mix_or_epi[sanger_celltypes$sample == "HC_Sanger6"]<-"mix_epi_immune"
sanger_celltypes$mix_or_epi[sanger_celltypes$sample == "HC_Sanger9"]<-"mix_epi_immune"
sanger_celltypes$proportion_total[sanger_celltypes$mix_or_epi == "mix_epi_immune"]<-(sanger_celltypes$proportion)/3
sanger_celltypes$proportion_total[sanger_celltypes$mix_or_epi == "epithelial_only"]<-(sanger_celltypes$proportion)/6
write.csv(sanger_celltypes, "~/cellcounts_Snager_data.csv")

ggplot(sanger_celltypes, aes(mix_or_epi,proportion_total, fill=celltype))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","deepskyblue2", "hotpink2")) 

ggplot(sanger_50_50_celltypes, aes(group,proportion, fill=celltype))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","deepskyblue2", "hotpink2")) 


## integrate all datasets using SCT integration (Seurat)





# make cell categories of smillie cell types

                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Best4+ Enterocytes"]<-"epithelium"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Enterocytes"]<-"epithelium"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Enterocyte Progenitors"]<-"epithelium"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Immature Enterocytes 1"]<-"epithelium"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Immature Enterocytes 2"]<-"epithelium"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "TA 1"]<-"epithelium"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "TA 2"]<-"epithelium"
                      
                      # epi general
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "M cells"]<-"epithelium"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Cycling TA"]<-"epithelium"
                      
                      # epi secretory
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Goblet"]<-"epithelium"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Enteroendocrine"]<-"epithelium"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Immature Goblet"]<-"epithelium"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Secretory TA"]<-"epithelium"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Stem"]<-"epithelium"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Tuft"]<-"epithelium"
                      
                      
                      #Tcell
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD4+ Activated Fos-lo"]<-"Tcell"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD4+ Activated Fos-hi"]<-"Tcell"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD4+ Memory"]<-"Tcell"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD4+ PD1+"]<-"Tcell"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD8+ IELs"]<-"Tcell"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD8+ IL17+"]<-"Tcell"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD8+ LP"]<-"Tcell"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "MT-hi"]<-"Tcell"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Tregs"]<-"Tcell"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Cycling T"]<-"Tcell"
                      
                      
                      # ILS?NK
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "ILCs"]<-"Lymphoid"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "NKs"]<-"Lymphoid"
                      
                      
                      
                      
                      #B
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Cycling B"]<-"Bcell"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "GC"]<-"Bcell"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Plasma"]<-"Bcell"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Follicular"]<-"Bcell"
                      
                      #Myeloid
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD69+ Mast"]<-"myeloid"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD69- Mast"]<-"myeloid"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Inflammatory Monocytes"]<-"myeloid"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Macrophages"]<-"myeloid"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Cycling Monocytes"]<-"myeloid"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "DC1"]<-"myeloid"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "DC2"]<-"myeloid"
                      
                      #mesenchymal
                      
                      
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Glia"]<-"Glia"
                      
                      
                      
                      # fibroblasts
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Inflammatory Fibroblasts"]<-"fibroblast"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "RSPO3+"]<-"fibroblast"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT2B+ Fos-hi"]<-"fibroblast"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT2B+ Fos-lo 1"]<-"fibroblast"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT2B+ Fos-lo 2"]<-"fibroblast"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT5B+ 1"]<-"fibroblast"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT5B+ 2"]<-"fibroblast"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT2B+ Fos-lo 1"]<-"fibroblast"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Myofibroblasts"]<-"fibroblast"
                      
                      #endothelial
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Post-capillary Venules"]<-"endothelial"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Pericytes"]<-"endothelial"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Endothelial"]<-"endothelial"
                      methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Microvascular"]<-"endothelial"
                      
                      
                      
               
