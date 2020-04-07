---
title: "Celltype prediction using SingleR"
author: W. Uniken Venema
date: April 7 2020
---

### Sample prep

install libraries
```
library(Seurat)
library(devtools)
```
download and load data from Smillie et al. (Cell,2019) as a reference
```
epi_fresh_split_smillie<-Read10X("/epi_smillie_2019/", unique.features = TRUE,  gene.column=1)
imm_fresh_split_smillie<-Read10X("/imm_smillie_2019/", unique.features = TRUE,  gene.column=1)
fib_fresh_split_smillie<-Read10X("/fibro_smillie_2019/", unique.features = TRUE,  gene.column=1)
```
create Seurat objects
```
epi_fresh_split <- CreateSeuratObject(counts = epi_fresh_split_smillie, project = "Smillie", min.cells = 3, min.features = 200)
imm_fresh_split <- CreateSeuratObject(counts = imm_fresh_split_smillie, project = "Smillie", min.cells = 3, min.features = 200)
fib_fresh_split <- CreateSeuratObject(counts = fib_fresh_split_smillie, project = "Smillie", min.cells = 3, min.features = 200)
```
filter datasets for <60% mitochondrial reads per cell, run PCA and UMAP and find clusters
```
epi_fresh_split[["percent.mt"]] <- PercentageFeatureSet(epi_fresh_split, pattern = "^MT-")
epi_fresh_split<-subset(epi_fresh_split, subset=percent.mt < 60)

epi_fresh_split<-SCTransform(epi_fresh_split,vars.to.regress = "percent.mt", verbose = F)
epi_fresh_split <- RunPCA(epi_fresh_split, verbose = FALSE)
epi_fresh_split <- RunUMAP(epi_fresh_split, dims = 1:30, verbose = FALSE)
epi_fresh_split <- FindNeighbors(epi_fresh_split, dims = 1:30, verbose = FALSE)
epi_fresh_split <- FindClusters(epi_fresh_split, resolution = 0.2, verbose = FALSE)

imm_fresh_split[["percent.mt"]] <- PercentageFeatureSet(imm_fresh_split, pattern = "^MT-")
imm_fresh_split<-subset(imm_fresh_split, subset=percent.mt < 60)

imm_fresh_split<-SCTransform(imm_fresh_split,vars.to.regress = "percent.mt", verbose = F)
imm_fresh_split <- RunPCA(imm_fresh_split, verbose = FALSE)
imm_fresh_split <- RunUMAP(imm_fresh_split, dims = 1:30, verbose = FALSE)
imm_fresh_split <- FindNeighbors(imm_fresh_split, dims = 1:30, verbose = FALSE)
imm_fresh_split <- FindClusters(imm_fresh_split, resolution = 0.2, verbose = FALSE)

fib_fresh_split[["percent.mt"]] <- PercentageFeatureSet(fib_fresh_split, pattern = "^MT-")
fib_fresh_split<-subset(fib_fresh_split, subset=percent.mt < 60)

fib_fresh_split<-SCTransform(fib_fresh_split,vars.to.regress = "percent.mt", verbose = F)
fib_fresh_split <- RunPCA(fib_fresh_split, verbose = FALSE)
fib_fresh_split <- RunUMAP(fib_fresh_split, dims = 1:30, verbose = FALSE)
fib_fresh_split <- FindNeighbors(fib_fresh_split, dims = 1:30, verbose = FALSE)
fib_fresh_split <- FindClusters(fib_fresh_split, resolution = 0.2, verbose = FALSE)
```
run RPCA integration of 3 datasets
```
smillie_integration_list<-list(fib_fresh_split, imm_fresh_split, epi_fresh_split)
features <- SelectIntegrationFeatures(object.list = smillie_integration_list, nfeatures = 3000)
smillie_integration_list <- PrepSCTIntegration(object.list = smillie_integration_list, anchor.features = features)
smillie_integration_list <- lapply(X = smillie_integration_list, FUN = RunPCA, verbose = FALSE, features = features)
anchors <- FindIntegrationAnchors(object.list = smillie_integration_list, anchor.features = features, normalization.method = "SCT", reduction = "rpca")
rm(smillie_integration_list)
smillie.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
rm(anchors)
smillie.integrated <- RunPCA(smillie.integrated, verbose = FALSE)
smillie.integrated <- RunUMAP(smillie.integrated, dims = 1:30)
```
add metadata from Smillie et al. to integrated dataset
```
CellsMeta<-data.frame(smillie.integrated@meta.data)
CellsMeta$NAME<-row.names(CellsMeta)
metadata<-read.table(metadata.txt)
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
row.names(CellsMeta)=NULL
row.names(x)<-NULL
CellsMeta<-keeping.order(CellsMeta, merge, y=metadata, by = "NAME", all=T)
row.names(CellsMeta)<-CellsMeta$NAME
CellsMeta$NAME<-NULL
smillie.integrated <- AddMetaData(smillie.integrated, CellsMeta)
```
subset this integrated dataset to max 300 cells per celltype of only healthy samples to enable celltype prediction analysis on PC
```
Idents(smillie.integrated)<-smillie.integrated[["Health"]]
smillie_subset_healthy<-subset(smillie.integrated, subset = Health == "Healthy")
Idents(smillie_subset_healthy)<-smillie.integrated[["Cluster"]]
smillie_subset_healthy<-subset(smillie.integrated, downsample = 300)
```

### Cell type prediction
prepare each dataset using the same preparation and filtering settings as described here-above
Then, load the necessary libraries
```
library(scRNAseq)
library(SingleR)
library(scater)
```
add cell type 'categories' and 'subcategories' (based on paper Smillie et al) to metadata to enable courser cell type classification
```
merge_meta<-smillie_subset_healthy@meta.data
merge_meta$category[merge_meta$Cluster == "Best4+ Enterocytes"]<-"epithelium"
merge_meta$category[merge_meta$Cluster == "Enterocytes"]<-"epithelium"
merge_meta$category[merge_meta$Cluster == "Enterocyte Progenitors"]<-"epithelium"
merge_meta$category[merge_meta$Cluster == "Immature Enterocytes 1"]<-"epithelium"
merge_meta$category[merge_meta$Cluster == "Immature Enterocytes 2"]<-"epithelium"
merge_meta$category[merge_meta$Cluster == "TA 1"]<-"epithelium"
merge_meta$category[merge_meta$Cluster == "TA 2"]<-"epithelium"

# epi general
merge_meta$category[merge_meta$Cluster == "M cells"]<-"epithelium"
merge_meta$category[merge_meta$Cluster == "Cycling TA"]<-"epithelium"

# epi secretory
merge_meta$category[merge_meta$Cluster == "Goblet"]<-"epithelium"
merge_meta$category[merge_meta$Cluster == "Enteroendocrine"]<-"epithelium"
merge_meta$category[merge_meta$Cluster == "Immature Goblet"]<-"epithelium"
merge_meta$category[merge_meta$Cluster == "Secretory TA"]<-"epithelium"
merge_meta$category[merge_meta$Cluster == "Stem"]<-"epithelium"
merge_meta$category[merge_meta$Cluster == "Tuft"]<-"epithelium"


#Tcell
merge_meta$category[merge_meta$Cluster == "CD4+ Activated Fos-lo"]<-"Tcell"
merge_meta$category[merge_meta$Cluster == "CD4+ Activated Fos-hi"]<-"Tcell"
merge_meta$category[merge_meta$Cluster == "CD4+ Memory"]<-"Tcell"
merge_meta$category[merge_meta$Cluster == "CD4+ PD1+"]<-"Tcell"
merge_meta$category[merge_meta$Cluster == "CD8+ IELs"]<-"Tcell"
merge_meta$category[merge_meta$Cluster == "CD8+ IL17+"]<-"Tcell"
merge_meta$category[merge_meta$Cluster == "CD8+ LP"]<-"Tcell"
merge_meta$category[merge_meta$Cluster == "MT-hi"]<-"Tcell"
merge_meta$category[merge_meta$Cluster == "Tregs"]<-"Tcell"
merge_meta$category[merge_meta$Cluster == "Cycling T"]<-"Tcell"


# ILS and NK
merge_meta$category[merge_meta$Cluster == "ILCs"]<-"Lymphoid"
merge_meta$category[merge_meta$Cluster == "NKs"]<-"Lymphoid"

#B
merge_meta$category[merge_meta$Cluster == "Cycling B"]<-"Bcell"
merge_meta$category[merge_meta$Cluster == "GC"]<-"Bcell"
merge_meta$category[merge_meta$Cluster == "Plasma"]<-"Bcell"
merge_meta$category[merge_meta$Cluster == "Follicular"]<-"Bcell"

#Myeloid
merge_meta$category[merge_meta$Cluster == "CD69+ Mast"]<-"myeloid"
merge_meta$category[merge_meta$Cluster == "CD69- Mast"]<-"myeloid"
merge_meta$category[merge_meta$Cluster == "Inflammatory Monocytes"]<-"myeloid"
merge_meta$category[merge_meta$Cluster == "Macrophages"]<-"myeloid"
merge_meta$category[merge_meta$Cluster == "Cycling Monocytes"]<-"myeloid"
merge_meta$category[merge_meta$Cluster == "DC1"]<-"myeloid"
merge_meta$category[merge_meta$Cluster == "DC2"]<-"myeloid"

#Glia
merge_meta$category[merge_meta$Cluster == "Glia"]<-"Glia"

# fibroblasts
merge_meta$category[merge_meta$Cluster == "Inflammatory Fibroblasts"]<-"fibroblast"
merge_meta$category[merge_meta$Cluster == "RSPO3+"]<-"fibroblast"
merge_meta$category[merge_meta$Cluster == "WNT2B+ Fos-hi"]<-"fibroblast"
merge_meta$category[merge_meta$Cluster == "WNT2B+ Fos-lo 1"]<-"fibroblast"
merge_meta$category[merge_meta$Cluster == "WNT2B+ Fos-lo 2"]<-"fibroblast"
merge_meta$category[merge_meta$Cluster == "WNT5B+ 1"]<-"fibroblast"
merge_meta$category[merge_meta$Cluster == "WNT5B+ 2"]<-"fibroblast"
merge_meta$category[merge_meta$Cluster == "WNT2B+ Fos-lo 1"]<-"fibroblast"
merge_meta$category[merge_meta$Cluster == "Myofibroblasts"]<-"fibroblast"

#endothelial
merge_meta$category[merge_meta$Cluster == "Post-capillary Venules"]<-"endothelial"
merge_meta$category[merge_meta$Cluster == "Pericytes"]<-"endothelial"
merge_meta$category[merge_meta$Cluster == "Endothelial"]<-"endothelial"
merge_meta$category[merge_meta$Cluster == "Microvascular"]<-"endothelial"

table(merge_meta$category)

merge_meta$subcategory<-merge_meta$category
merge_meta$subcategory[merge_meta$Cluster == "Best4+ Enterocytes"]<-"epithelium_A"
merge_meta$subcategory[merge_meta$Cluster == "Enterocytes"]<-"epithelium_A"
merge_meta$subcategory[merge_meta$Cluster == "Enterocyte Progenitors"]<-"epithelium_A"
merge_meta$subcategory[merge_meta$Cluster == "Immature Enterocytes 1"]<-"epithelium_A"
merge_meta$subcategory[merge_meta$Cluster == "Immature Enterocytes 2"]<-"epithelium_A"
merge_meta$subcategory[merge_meta$Cluster == "TA 1"]<-"epithelium_A"
merge_meta$subcategory[merge_meta$Cluster == "TA 2"]<-"epithelium_A"

# epi general
merge_meta$subcategory[merge_meta$Cluster == "M cells"]<-"epithelium"
merge_meta$subcategory[merge_meta$Cluster == "Cycling TA"]<-"epithelium"

# epi secretory
merge_meta$subcategory[merge_meta$Cluster == "Goblet"]<-"epithelium_S"
merge_meta$subcategory[merge_meta$Cluster == "Enteroendocrine"]<-"epithelium_S"
merge_meta$subcategory[merge_meta$Cluster == "Immature Goblet"]<-"epithelium_S"
merge_meta$subcategory[merge_meta$Cluster == "Secretory TA"]<-"epithelium_S"
merge_meta$subcategory[merge_meta$Cluster == "Stem"]<-"epithelium_S"
merge_meta$subcategory[merge_meta$Cluster == "Tuft"]<-"epithelium_S"

merge_meta<-merge_meta[c(14,15)]
smillie_subset_healthy <- AddMetaData(smillie_subset_healthy, merge_meta)
```
load the Seuratfile of the dataset you want to predict cell types for
```
dataset_seurat<-load("dataset.rds")
```
convert datasets to SingleCellExperiment format and lognormalize counts
``
smillie_subset<-as.SingleCellExperiment(smillie_subset_healthy, assay="RNA")
smillie_subset<-logNormCounts(smillie_subset)
dataset<-as.SingleCellExperiment(dataset_seurat, assay="RNA")
dataset<-logNormCounts(dataset)
```
and predict cell types
```
table(smillie_subset$subcategory) # cell types as reference
pred_celltypes <- SingleR(test=dataset, ref=smillie_subset, labels=smillie_subset$Cluster, de.method="wilcox")
pred_cat <- SingleR(test=dataset, ref=smillie_subset, labels=smillie_subset$category, de.method="wilcox")
pred_subcat <- SingleR(test=dataset, ref=smillie_subset, labels=smillie_subset$subcategory, de.method="wilcox")
```
add cell type predictions to metadata
```
CellsMeta<-dataset_seurat@meta.data

x<-data.frame(pred_subcat)
x$NAME<-row.names(x)
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
x<-x[,c(11:16)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
row.names(CellsMeta)<-CellsMeta$NAME
dim(CellsMeta)
CellsMeta<-CellsMeta[,c(23:27)]
head(CellsMeta)
colnames(CellsMeta)[1]<-"first_labels_subcategory"
colnames(CellsMeta)[2]<-"first_tuning_scores_subcategory"
colnames(CellsMeta)[3]<-"second_tuning_scores_subcategory"
colnames(CellsMeta)[4]<-"labels_subcategory"
colnames(CellsMeta)[5]<-"pruned_labels_subcategory"

x<-data.frame(pred_cat)
x$NAME<-row.names(x)
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
x<-x[,c(11:16)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
row.names(CellsMeta)<-CellsMeta$NAME
dim(CellsMeta)
head(CellsMeta)
colnames(CellsMeta)[6]<-"first_labels_category"
colnames(CellsMeta)[7]<-"first_tuning_scores_category"
colnames(CellsMeta)[8]<-"second_tuning_scores_category"
colnames(CellsMeta)[9]<-"labels_category"
colnames(CellsMeta)[10]<-"pruned_labels_category"

x<-data.frame(pred_celltypes)
x$NAME<-row.names(x)
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
x<-x[,c(11:16)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
row.names(CellsMeta)<-CellsMeta$NAME
dim(CellsMeta)
CellsMeta<-CellsMeta[,c(23:27)]
head(CellsMeta)
colnames(CellsMeta)[11]<-"first_labels_celltypes"
colnames(CellsMeta)[12]<-"first_tuning_celltypes"
colnames(CellsMeta)[13]<-"second_tuning_scores_celltypes"
colnames(CellsMeta)[14]<-"labels_celltypes"
colnames(CellsMeta)[15]<-"pruned_labels_celltypes"

dataset_seurat <- AddMetaData(dataset_seurat, CellsMeta)
```
visualize predicted Smillie celltypes on UMAP of dataset
```
DimPlot(dataset, group.by="pruned_labels_subcategory", label=T,pt.size = 0.2 )
```

