**Seurat single cell seq analysis**
=====
author: WTC
date: 20171020

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
library(scde) # Differential Expression
library(biomaRt)
library(org.Hs.eg.db) # Gene name manipulation
library(ggplot2)
library(tidyr)
library(devtools)
library(biomaRt)
```

**create new seurat (v2.0) object of raw data file + QC**
---
**filter for genes expressed in =>3 cells, with =>200 genes per cell**

```
datafile <- CreateSeuratObject(raw.data = x, min.cells = 3, min.genes = 200, project = "datafile")
```
**discover mitochondrial genes (high percentage of total expressed genes can mark cell breakdown)**
```
mito.genes <- grep(pattern = "^MT\\.", x = rownames(x = datafile@data), value = TRUE)
percent.mito <- colSums(datafile@raw.data[mito.genes, ])/colSums(datafile@raw.data)
datafile <- AddMetaData(object = datafile, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = datafile, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = datafile, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = datafile, gene1 = "nUMI", gene2 = "nGene")
```
**filter for ngene: 200-2500, percent mito: <5%**
```
datafile <- FilterCells(object = datafile, subset.names = "nGene", low.thresholds = 200, high.thresholds = 2500)
datafile <- FilterCells(object = datafile, subset.names = "percent.mito", low.thresholds = -Inf, high.thresholds =  0.05)
```
**log normalize the data**
```
datafile <- NormalizeData(object = datafile, normalization.method = "LogNormalize", scale.factor = 10000)
```

**add metadata**
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
**add location**
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

**add patient**
```
CellsMeta = datafile@meta.data
head(CellsMeta)
```
the patient number is in our cellnames, we extract the number using:
```
patient<-data.frame(substr(colnames(datafile@data), 1, 1))
```
```
row.names(patient)=colnames(datafile@data)
CellsMeta["patient"] <- patient[,1]
head(CellsMeta)
CellsMetaTrim <- subset(CellsMeta, select = c("patient"))
head(CellsMetaTrim)
file_per_patient <- AddMetaData(allpts_perplexity500, CellsMetaTrim)
head(file_per_patient@meta.data)
```

**change rownames from ensembl gene names to gene symbols**
-----
```
CellsMeta = datafile@data
head(CellsMeta)
```
create df with genes
```
gene<-data.frame(sapply(strsplit(rownames(datafile@data), split='-', fixed=TRUE), function(x) (x[2])))
```
make genes unique rownames of seurat file
```
rownames(gene) = make.names(gene[,1], unique=TRUE)
rownames(CellsMeta) <- rownames(gene)
head(CellsMeta)
datafile@data = CellsMeta
head(datafile@data)
```

**scale, perform PCAs and JackStraw analysis to define significance of PCAs**
---
scale, regressing for nUMI and percent.mit, because these cause unwanted variation in you expression data
```
datafile <- ScaleData(object = datafile, vars.to.regress = c("nUMI", "percent.mito"))
```
find the variable genes to perform a PCA on
```
datafile<-FindVariableGenes(datafile, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 2, y.cutoff = 2)
datafile <- RunPCA(object = datafile, pc.genes = row.names(datafile@data), do.print = TRUE, pcs.print = 1:5, genes.print = 5)
```
visualize results
```
PrintPCA(object = datafile, pcs.print = 1:7, genes.print = 5, use.full = FALSE)
VizPCA(object = datafile, pcs.use = 1:2)
PCAPlot(object = datafile, dim.1 = 1, dim.2 = 2)
```
perform a JackStraw with 5 replications to calculate the significance of these PCAs (takes a short while)
```
datafile = JackStraw(datafile, num.replicate = 5, do.print = T)
```
project your PCA 
```
datafile <- ProjectPCA(object = datafile, do.print = T)
```
explore gene expression in PCAs
```
PCHeatmap(object = datafile, pc.use = 6, cells.use = 500, do.balanced = TRUE,label.columns = FALSE)
PCHeatmap(object = datafile, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
```
run PCA on significantly expressed genes using the (here) 20 significant PCAs
```
datafile_siggenes<-PCASigGenes(datafile, 1:20, pval.cut=1e-5, max.per.pc = 250)
datafile_sig<-RunPCA(datafile, pc.genes = datafile_siggenes, do.print=T)
```
run jackstraw using 100 replicates, can take a long time. preferably  performed on cluster
```
datafile_sig <- JackStraw(object = datafile_sig, num.replicate = 100, do.print = FALSE)
```
visualize significance of PCAs
```
JackStrawPlot(object = datafile, PCs = 1:15)
```
in the kink in the plot, the expected significance cutoff of the  PCAs lies:
```
PCElbowPlot(object = datafile)
```

**Find clusters**
---
using relevant PCAs, adapt resolution to amount of clusters that is wishful/makes sense
```
datafile_sig <- FindClusters(object = datafile_sig, reduction.type = "pca", dims.use = 1:18, resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = datafile)
```

**run TSNE and visualize**
---
```
datafile_sig <- RunTSNE(object = datafile_sig, dims.use = 1:18, do.fast = TRUE)
TSNEPlot(object = datafile_sig)
```
by 'patient'
```
TSNEPlot(object = datafile_sig, group.by="patient")
```

**subset cluster cells**
---

**Subset cluster**
for example cluster 0
```
cluster0<-SubsetData(datafile, ident=0)
```
Run/Plot TSNE
```
cluster0.1=RunTSNE(cluster0, dims.use=1:30)
TSNEPlot(cluster0.1)
table(cluster0.1@meta.data$res.1, cluster0.1@meta.data$patient)
table(cluster0.1@meta.data$res.1, cluster0.1@meta.data$tissue)
```
```
dev.copy(pdf,paste("TSNEPlot.pdf"))
dev.off()
```

define expression markers
```
markers_cluster0.1<-FindAllMarkers(cluster0.1)
markers_cluster0.1 %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10
write.csv(markers_cluster0.1, "markers_cluster0.1.csv")
```


**After SCDE analysis: extract CD risk genes from DE genes table**
---

**open DE .txt file**
```
DE_mucosa_CD4_CD8<-read.table("..")
```
**create column 'Gene'**
```
DE_mucosa_CD4_CD8$Gene<-sapply(strsplit(row.names(DE_mucosa_CD4_CD8), split='-', fixed=TRUE), function(x) (x[2]))
```
load genes of interest file
```
GOI<-read.table("..", sep=";", header=TRUE)
```
match DE list and GOI
```
GOI_pos<-merge(GOI, DE_mucosa_CD4_CD8, by="Gene", all=FALSE)
```
subset DE_mucosa_CD4_CD8 with only genes above DE value for two-sided testing (<-1.96 and >1.96)
```
GOI_DE_CD4pos<-subset(GOI_pos, cZ>1.96)
write.csv(GOI_DE_CD4pos, "..")
GOI_DE_CD8pos<-subset(GOI_pos, cZ<(-1.96))
write.csv(GOI_DE_CD8pos, "..")
```
load genes of interest file
```
DRUGTARGETGENES<-read.table("..", sep=";", header=TRUE)
```
match DE list and DRUGTARGETGENES
```
DRUGTARGETGENES_pos<-merge(DRUGTARGETGENES, DE_blood_CD4_CD8, by="Gene", all=FALSE)
```
save matched list
```
write.csv(DRUGTARGETGENES_pos, "..")
```

subset DE_mucosa_CD4_CD8 with only genes above DE value for two-sided testing (<-1.96 and >1.96)
```
DRUGTARGETGENES_CD4pos<-subset(DRUGTARGETGENES_pos, cZ>1.95)
write.csv(DRUGTARGETGENES_CD4pos, "..")
DRUGTARGETGENES_CD8pos<-subset(DRUGTARGETGENES_pos, cZ<(-1.95))
write.csv(DRUGTARGETGENES_CD8pos, "..")
```
