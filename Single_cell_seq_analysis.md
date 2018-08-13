**Seurat single cell seq analysis**
=====
Author: WTC .   
Date: 20180813 
  
Load libraries
----
```
library(fpc)  # Density based clustering dbscan
library(gplots)  # Colorpanel
library(scatterplot3d)  # 3D plotting
library(monocle)
library(tsne)  # Non-linear ordination
library(pheatmap)
library(MASS)
library(cluster)
library(mclust)
library(flexmix)
library(lattice)
library(amap)
library(RColorBrewer)
library(locfit)
library(Seurat)
library(vioplot)
library(dplyr) # Dataframe manipulation
library(Matrix) # Sparse matrices
library(useful) # Corner function
library(scater) # Single Cell QC
library(org.Hs.eg.db) # Gene name manipulation
library(ggplot2)
library(tidyr)
library(devtools)
library(biomaRt)
library(MAST) # for DE analyses
```
**Create new seurat (v2.0) object of raw data file + QC**
---
load data and metadata
```
data<-read.csv("~/..csv")
metadata<-read.csv("~/..csv")
```
erase collagenase digestion-induced genes from dataset
```
collagenasegenes<-read.csv("~/..csv")
genes<-as.data.frame(collagenasegenes$gene)
datafile$gene<-row.names(datafile)
row.names(datafile)<-NULL
no_collagen<-subset(data, !(data$gene %in% genes$`collagenasegenes$gene`))
row.names(no_collagen)<-no_collagen$gene
no_collagen$gene<-NULL
dim(no_collagen)
```
extract interesting cell subsets
```
bloodcells<-metadata[metadata$tissue == "Blood",]
ielcells<-metadata[metadata$tissue == "IEL",]
lplcells<-metadata[metadata$tissue == "LPL",]
mucosacells<-rbind(lplcells, ielcells)
ielcells<-as.vector(ielcells$sample)
lplcells<-as.vector(lplcells$sample)

pt1cells<-metadata[metadata$patient == "pt1",]
pt1cells<-as.vector(pt1cells$sample)
pt2cells<-metadata[metadata$patient == "pt2",]
pt2cells<-as.vector(pt2cells$sample)
pt3cells<-metadata[metadata$patient == "pt3",]
pt3cells<-as.vector(pt3cells$sample)

cd8poscells<-facs_meta[metadata$X640.CD8ab.APC.Height == 1,]
cd8negcells<-facs_meta[metadata$X640.CD8ab.APC.Height == 0,]
cd8posbloodcells<-cd8poscells[which(cd8poscells$tissue.x == "Blood"),]
cd8posmucosacells<-merge(mucosacells, cd8poscells, by="sample", all=F)
cd8negbloodcells<-cd8negcells[cd8negcells$tissue.x == "Blood",]
cd8negmucosacells<-merge(mucosacells, cd8negcells, by="sample", all=F)
cd8negmucosacells<-as.vector(cd8negmucosacells$sample)
cd8poscells<-as.vector(cd8poscells$sample)
cd8negcells<-as.vector(cd8negcells$sample)
cd8posbloodcells<-as.vector(cd8posbloodcells$sample)
cd8negbloodcells<-as.vector(cd8negbloodcells$sample)
cd8posmucosacells<-as.vector(cd8posmucosacells$sample)
mucosacells<-as.vector(mucosacells$sample)
bloodcells<-as.vector(bloodcells$sample)
```
create separate dataframes for these subsets, using them as a raw.data input in the CreateSeuratObject function
```
```
**Create seuratobject and filter for genes expressed in =>3 cells**
```
seuratfile <- CreateSeuratObject(raw.data = data, min.cells = 3, project = "CD")
```
**Calculate % mitochondrial genes (high percentage of total expressed genes can mark cell breakdown), number of UMI and number of genes**
```
mito.genes <- grep(pattern = "^MT\\.", x = rownames(x = seuratfile@data), value = TRUE)
percent.mito <- Matrix::colSums(seuratfile@raw.data[mito.genes, ])/Matrix::colSums(seuratfile@raw.data)
datafile <- AddMetaData(object = datafile, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = datafile, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = datafile, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = datafile, gene1 = "nUMI", gene2 = "nGene")
```
**Filter for cells with =>200 AND =<2500 genes per cell, and <5% mitochondrial genes**
```
seuratfile <- FilterCells(object = seuratfile, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
```
**Log normalize**
```
seuratfile <- NormalizeData(object = seuratfile, normalization.method = "LogNormalize", scale.factor = 10000)
```
**Add metadata**
---
```
metadata<-read.csv("..")
```
add a function that will keep the order of your cells when you extract the metadatafile to add information 
```
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
```
extract file that is a copy from @meta.data
```
CellsMeta = datafile@meta.data
head(CellsMeta)
CellsMeta$cellen=row.names(CellsMeta)
row.names(CellsMeta)=NULL
dim(CellsMeta)
```
merge the required additional metadata with the extracted meta.data file, using keeping.order
```
CellsMeta<-keeping.order(CellsMeta, merge, y=metadata, by = "cellen", all=FALSE)
```
make sure this file has the colnames of your CellsMeta file as rownames
```
row.names(CellsMeta)<-CellsMeta$cellen
dim(CellsMeta)
CellsMeta$cellen<-NULL
head(CellsMeta)
```
add the new metadatafile to your Seurat file
```
datafile <- AddMetaData(datafile, CellsMeta)
```
**Add patient**
```
CellsMeta = seuratfile@meta.data
head(CellsMeta)
```
if the patient number is in our cellnames, we extract the number using:
```
patient<-data.frame(substr(colnames(seuratfile@data), 1, 1))
row.names(patient)=colnames(seuratfile@data)
CellsMeta["patient"] <- patient[,1]
head(CellsMeta)
CellsMetaTrim <- subset(CellsMeta, select = c("patient"))
head(CellsMetaTrim)
file_per_patient <- AddMetaData(seuratfile, CellsMetaTrim)
head(seuratfile@meta.data)
```
**If necessary: change rownames from ensembl gene names to gene symbols**
-----
```
CellsMeta = seuratfile@data
head(CellsMeta)
```
create df with genes
```
gene<-data.frame(sapply(strsplit(rownames(seuratfile@data), split='-', fixed=TRUE), function(x) (x[2])))
```
make genes unique rownames of seurat file
```
rownames(gene) = make.names(gene[,1], unique=TRUE)
rownames(CellsMeta) <- rownames(gene)
head(CellsMeta)
seuratfile@data = CellsMeta
head(seurat@data)
```
**Scale, perform PCAs and JackStraw analysis to define significance of PCAs**
  
scale, regressing for nUMI and perc.mit, because these cause unwanted variation in you expression data
```
seuratfile <- ScaleData(object = seuratfile, vars.to.regress = c("nUMI", "percent.mito"))
seuratfile<-FindVariableGenes(seuratfile, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 2, y.cutoff = 2)
seuratfile <- RunPCA(object = seuratfile, pc.genes = row.names(seuratfile@data), do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
seuratfile <- ProjectPCA(object = seuratfile, do.print = T, pcs.store = 40)
seuratfile <- JackStraw(object = seuratfile, num.replicate = 100)
```
visualize results
```
PrintPCA(object = seuratfile, pcs.print = 1:7, genes.print = 5, use.full = FALSE)
VizPCA(object = seuratfile, pcs.use = 1:2)
PCAPlot(object = seuratfile, dim.1 = 1, dim.2 = 2)
JackStrawPlot(object = seuratfile, PCs = 1:15)
```
in the kink in the plot, the expected significance cutoff of the  PCAs lies:
```
PCElbowPlot(object = seuratfile)
```
**Find clusters**
  
using relevant PCAs, adapt resolution to amount of clusters that is wishful/makes sense
```
seuratfile <- FindClusters(object = seuratfile, reduction.type = "pca", dims.use = 1:18, resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = seuratfile)
```
explore gene expression in PCAs
```
PCHeatmap(object = seuratfile, pc.use = 6, cells.use = 500, do.balanced = TRUE,label.columns = FALSE)
PCHeatmap(object = seuratfile, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
```
**Run TSNE and visualize**
```
seuratfile <- RunTSNE(object = seuratfile, dims.use = 1:18, do.fast = TRUE)
TSNEPlot(object = seuratfile)
```
by 'patient', or other variable from meta.data
```
TSNEPlot(object =seuratfile, group.by="patient")
```
**Define cell types**
  
Cluster with epitopic CD8+ cells and epitopic CD8- cells separately, and cluster with tissue of origin. Assign cell types to clusters based on cell type marker genes known from literature. Make consensus clustering of epitopic marker clustering and tissue of origin clustering per cell.
  
Add consensus cell type names to metadata (see above)
  
**DE analysis with genes >1% expressed in dataset, using MAST**
  
set 'only.pos=F' to obtain both up- and downregulated genes
```
seuratfile<-SetAllIdent(seuratfile, "eight_cell_types")
allcells_DE_markers = FindAllMarkers(seuratfile, min.pct = 0.01, only.pos = T, test.use = "MAST")
```
subset CDriskgenes_DE for significant results only
```
allcells_DE_markers<-allcells_DE_markers[allcells_DE_markers$p_val_adj < 0.05,]
```
**CD risk gene analysis**
  
load genes of interest file
```
CDriskgenes<-read.csv("~/...txt")
```
obtain DE CDriskgenes
```
CDriskgenes_DE<-merge(CDriskgenes, allcells_DE_markers, by="gene", all=FALSE)
```
**IBD Drugtarget analysis**
  
load genes of interest file
```
DRUGTARGETGENES<-read.csv("~/...txt")
```
match DE list and DRUGTARGETGENES
```
DRUGTARGETGENES_DE<-merge(DRUGTARGETGENES, allcells_DE_markers, by="gene", all=FALSE)
```
