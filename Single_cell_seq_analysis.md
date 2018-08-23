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
library(reactome.db)
library(clusterProfiler)
library(ReactomePA)
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
**Add metadata, for example on FACS surface markers, patient, gender, treatment, other**
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
**Clustering and cell types**
---
**Scale, perform PCAs and JackStraw analysis to define significance of PCAs**
  
scale, regressing for nUMI and perc.mit, because these cause unwanted variation in you expression data
```
seuratfile <- ScaleData(object = seuratfile, vars.to.regress = c("nUMI", "percent.mito")) # add "patient" or "gender" here for regression for patient or gender
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
  
Add consensus cell type names to metadata, and attach to seurat file (see above)
  
    
**DE analysis using MAST**
---
with genes >1% expressed in dataset; set 'only.pos=F' to obtain both up- and downregulated genes. Use SetAllIdent() to change DE parameters to for example "tissue" 
```
seuratfile<-SetAllIdent(seuratfile, "eight_cell_types")
allcells_DE_markers = FindAllMarkers(seuratfile, min.pct = 0.01, only.pos = T, test.use = "MAST")
```
subset allcells_DE_markers for significant (p<0.05) results only
```
allcells_DE_markers<-allcells_DE_markers[allcells_DE_markers$p_val_adj < 0.05,]
```
make seuratfile regressed for patient or gender, see above, and calculate overlapping genes between datasets that are regressed for resp. patient and gender, and that are not regressed for these factors
```
seuratfile_patientregressed<-SetAllIdent(seuratfile_patientregressed, "eight_cell_types")
allcells_patientregr_DE_markers = FindAllMarkers(seuratfile_patientregressed, min.pct = 0.01, only.pos = T, test.use = "MAST")
row.names(allcells_patientregr_DE_markers)<-NULL
allcells_patientregr_DE_markers<-allcells_patientregr_DE_markers[allcells_patientregr_DE_markers$p_val_adj < 0.05,]
dim(allcells_patientregr_DE_markers)
colnames(allcells_patientregr_DE_markers)[1]<-"p_val_patientregr"
colnames(allcells_patientregr_DE_markers)[2]<-"Av_logFC_patientregr"
colnames(allcells_patientregr_DE_markers)[3]<-"pct.1_patientregr"
colnames(allcells_patientregr_DE_markers)[4]<-"pct.2_patientregr"
colnames(allcells_patientregr_DE_markers)[5]<-"p_val_adj_patientregr"
colnames(allcells_patientregr_DE_markers)[6]<-"cluster_patientregr"
allcells_patientregr_DE_markers$gene_cluster<-paste0(allcells_patientregr_DE_markers$gene, sep="_", allcells_patientregr_DE_markers$cluster_patientregr)
allcells_DE_markers$gene_cluster<-paste0(allcells_DE_markers$gene, sep="_", allcells_DE_markers$cluster)
overlap_patient_w_o_regress<-merge(allcells_patientregr_DE_markers, allcells_DE_markers, by="gene_cluster", all=F)
percentage_overlapping_DE_genes<-nrow(overlap_patient_w_o_regress)/nrow(allcells_DE_markers)
```
**CD risk gene analysis**
  
load genes of interest file, as in supplementary data
```
CDriskgenes<-read.csv("~/...txt")
```
obtain DE CDriskgenes
```
CDriskgenes_DE<-merge(CDriskgenes, allcells_DE_markers, by="gene", all=FALSE)
```
**IBD Drugtarget analysis**
  
load genes of interest file, as extracted from opentargets.org - IBD drugtargets
```
DRUGTARGETGENES<-read.csv("~/...txt")
```
match DE list and DRUGTARGETGENES
```
DRUGTARGETGENES_DE<-merge(DRUGTARGETGENES, allcells_DE_markers, by="gene", all=FALSE)
```
**Pathway analyses using ReactomeDB**
---
for example for CTL_blood genes
```
CTL_Blood_genes<-subset(allcells_DE_markers, (allcells_DE_markers$cluster == "CTL_Blood")) 
CTL_Blood_genes<-CTL_Blood_genes$gene
```
convert gene symbols to entrez-id for Reactome
```
CTL_Blood_genes = bitr(CTL_Blood_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
```
take only Entrez-id column in list
```
CTL_Blood_entrez <- CTL_Blood_genes$ENTREZID
```
do pathway analysis
```
CTL_Blood_pathways <- enrichPathway(gene=CTL_Blood_entrez,pvalueCutoff=0.05, readable=T)
```
barplot pathways
```
barplot(CTL_Blood_pathways, showCategory=15)
```
dotplot enrichment
```
dotplot(CTL_Blood_pathways, showCategory=15)
```
**Comparison with previously published DE results**
---
create a df "genes_24799394", consisting of the DE genes between blood and IEL+LPL from PMID 24799394 in R, with column name "gene" and merge with DE genes between blood and mucosa (supplemental figure 3)
```
overlap_24799394<-merge(genes_24799394, allcells_DE_markers, by="gene", all=F)
number_overlapping_genes_previously_published<-nrow(overlap_24799394)
```
**Comparison with healthy CTL dataset**
---
download Cytotoxic T cell gene/cell matrix (raw) from https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/cytotoxic_t
```
CTL_10x<-Read10X(data.dir = "~/Downloads/matrices_mex 3/hg19/")
CTL_10x<-CreateSeuratObject(raw.data = CTL_10x, min.cells = 3, project = "healthy")
```
load 179 risk genes and all genes present in 1% of both healthy and CD CTL cells
```
risk_use<-read.csv("...csv")
genes_use<-read.csv("..csv")
```
**Select 1000x a subset of 251 cells, combine datasets using CCA and perform DE analyses on combined datasets**
```
```
subsetting
```
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
```
filter for ngene: 200-2500, percent mito: <5%
```
CTL_10x_subset <- FilterCells(object = CTL_10x_subset, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

```
log normalize
```
CTL_10x_subset <- NormalizeData(object = CTL_10x_subset, normalization.method = "LogNormalize", scale.factor = 10000)
CTL_10x_subset <- ScaleData(object = CTL_10x_subset, vars.to.regress = c("nUMI", "percent.mito"))
CTL_10x_subset@meta.data$dataset <- "subset10x"

save(CTL_10x_subset, file=paste0("CTL_10x_subset",i,".Rda"))
```
assign name to CD dataset
```
CTL_CD@meta.data$dataset <- "cd"
```
Combine the CD and subset10x cells into a single object using CCA
```
CTL_subset.combined <- RunCCA(CTL_10x_subset, CTL_CD, genes.use = intersect(rownames(CTL_10x_subset@data), rownames(CTL_CD@data)), num.cc = 30, scale.data=T)
```
calculate DE with genes present in the two datasets, in >1%, MAST, between CD  and subset10x cells
```
CTL_subset.combined<-SetAllIdent(CTL_subset.combined, "dataset")
CTL_subset.combined_DE_markers_1perc = FindAllMarkers(CTL_subset.combined, only.pos = T, test.use = "MAST", genes.use = genes_use$gene)
write.table(CTL_subset.combined_DE_markers_1perc, file = paste0("CTL_subset",i,".blood.10x._DEmarkers_MAST.txt"))
```
select the significant ones (p<0.05)
```
CTL_subset.combined_DE_markers_1perc<-CTL_subset.combined_DE_markers_1perc[CTL_subset.combined_DE_markers_1perc$p_val_adj <0.05,]
number_DE_genes[i,a]<-nrow(CTL_subset.combined_DE_markers_1perc)
```
determine risk genes
```
risk_all_filtered_subset<-merge(risk_use, CTL_subset.combined_DE_markers_1perc, by="gene", all=F)
write.csv(risk_all_filtered_subset, paste0("riskgenes_subset",i,"_vs_cd_filtered.csv"))
number_DE_riskgenes[i,b]<-nrow(risk_all_filtered_subset)}
```
write output file with numbers DE genes and numbers risk genes
```
write.csv(number_DE_genes, "number_DE_genes_cd_vs_10xsubsets.csv")
write.csv(number_DE_riskgenes, "number_DE_riskgenes_cd_vs_10xsubsets_newriskgenes.csv")
```
extract and count which genes are upregulated in CD 
```
for (n in 1:1000){
  cat("iteration", n, "\n")
  CTL_subset.combined_DE_markers_1perc<-read.table(paste0("~/nonrandomised_permutation_10xsubset_cd/CTL_subset", n, ".blood.10x._DEmarkers_MAST.txt"))
  # and the significant ones (p<0.05)
  CTL_subset.combined_DE_markers_1perc<-CTL_subset.combined_DE_markers_1perc[CTL_subset.combined_DE_markers_1perc$p_val_adj <0.05,]
  # and the CD upregulated ones
  CTL_subset.combined_DE_markers_1perc<-CTL_subset.combined_DE_markers_1perc[CTL_subset.combined_DE_markers_1perc$cluster == "cd",]
  number_DE_genes[n,a]<-nrow(CTL_subset.combined_DE_markers_1perc)
  # merge all DE results with risk genes
  risk_all_filtered_subset<-merge(risk_use, CTL_subset.combined_DE_markers_1perc, by="gene", all=F)
  write.csv(risk_all_filtered_subset, paste0("riskgenes_cd_subset",n ,"_vs_cd_filtered.csv"))
  number_DE_riskgenes[n,b]<-nrow(risk_all_filtered_subset)}
mean(number_DE_riskgenes$number_DEriskgenes)
mean(number_DE_genes$number_DEgenes)
```
merge all riskgenes found in 1000 permutations
```
riskgenes<- read.csv("riskgenes_cd_subset1_vs_cd_filtered.csv", row.names=1)
for(n in 2:1000){
  cat("iteration", n, "\n")
  x<-read.csv(paste0("riskgenes_cd_subset", n, "_vs_cd_filtered.csv"), row.names=1)
  riskgenes<-rbind(x, riskgenes)
}
```
write output file with counts of DE risk genes
```
write.csv(table(riskgenes$gene), "riskgene_counts_csv.csv")
```
merge all DE genes found in 1000 permutations
```
setwd("/..")
DEgenes<- read.table("CTL_subset1.blood.10x._DEmarkers_MAST.txt", row.names=1)
for(n in 2:1000){
  cat("iteration", n, "\n")
  x<-read.table(paste0("CTL_subset", n, ".blood.10x._DEmarkers_MAST.txt"), row.names=1)
  DEgenes<-rbind(x, DEgenes)
}
```
and the significant ones (p<0.05)
```
DEgenes<-DEgenes[DEgenes$p_val_adj <0.05,]
```
and the CD upregulated ones
```
DEgenes_CD<-DEgenes[DEgenes$cluster == "cd",]
```
write output file with counts of DE genes
```
write.csv(table(DEgenes_CD$gene), "gene_counts_cd.csv")
```
calculate means for subsequent permutation analyses
```
mean(number_DE_riskgenes$number_DEriskgenes)
mean(number_DE_genes$number_DEgenes)
```
same for subset10x
```
```
**Permutations for numer of risk genes found in DE genes**
```
genes_cd_10x<-genes_use # all genes present in at least 1% of CTL and 1% of CD samples
a<-"number"
results_permutation<-data.frame(a)
for(i in 1:100000){
  samp_row_idxs  <- sample(nrow(genes_cd_10x), mean(number_DE_genes$number_DEgenes))
  samp_gene_names <- as.data.frame(genes_cd_10x[samp_row_idxs,])
  colnames(samp_gene_names)[1]<-"gene"
  risk_genes_sample<- merge(samp_gene_names, riskgenes, by="gene", all=F)
  results_permutation[i,a]<-nrow(risk_genes_sample)}
```
and test whether it is higher than observed (mean(number_DE_riskgenes$number_DEriskgenes))
```
permutation_higher<-results_permutation[results_permutation$number > mean(number_DE_riskgenes$number_DEriskgenes),]
pval_10x_versus_cd=(nrow(permutation_higher))/100000 # pval= (# of obs > our mean value) / # of bootstraps
```
**Comparison with healthy naive CD4+ T cell dataset**
---
download Naive CD4-positive T cell gene/cell matrix (raw) from https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/naive_t
create new seurat (v2.0) object of raw data files
```
seuratfile <- CreateSeuratObject(raw.data = naive_cd4T_10x, min.cells = 3, project = "x")
```
filter and lognormalize
```
mito.genes <- grep(pattern = "^MT\\.", x = rownames(x = seuratfile@data), value = TRUE)
percent.mito <- Matrix::colSums(seuratfile@raw.data[mito.genes, ])/Matrix::colSums(seuratfile@raw.data)
seuratfile <- AddMetaData(object = seuratfile, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = seuratfile, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = seuratfile, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seuratfile, gene1 = "nUMI", gene2 = "nGene")
seuratfile <- FilterCells(object = seuratfile, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
seuratfile <- NormalizeData(object = seuratfile, normalization.method = "LogNormalize", scale.factor = 10000)
```
change synonymous gene name 'MAL' to 'CD8A'
```
seuratfile@meta.data$ident<-"Naive CD4T"
seuratfile<-SetAllIdent(seuratfile, "ident")
x<-seuratfile
MAL<-"MAL"
CD8A<-"CD8A"
row.names(x@data) <- plyr::mapvalues(x = row.names(x@data), from = MAL, to = CD8A)
```
plot expression of CD8A (also known as MAL), CD8B, CD4
```
VlnPlot(x, c("CD8A","CD8B", "CD4"))
```
plot expression of CD8A (also known as MAL), CD8B, CD4 on epitopic CD8-positive and -negative sets of CD patients. Seuratfile for epitopic CD8-positive of -negative cells only can be made following the script above 'Create new seurat (v2.0) object of raw data file + QC'
```
VlnPlot(seuratfile_cd8poscells, c("CD8A","CD8B", "CD4"))
VlnPlot(seuratfile_cd8negcells, c("CD8A","CD8B", "CD4"))
```
**Plotting**
---
tSNEL PBL, IEL, LPL
```
allcells_meta<-SetAllIdent(seuratfile_allcells_patientregr, "tissue")
TSNEPlot(allcells_meta, pt.size=2, colors.use =c("brightred", "brightblue", "darkgreen"))
```
tSNE mucosacells
```
seuratfile_mucosacells<-SetAllIdent(seuratfile_mucosacells, "eight_cell_types")
current.cluster.ids <- c("Cytotoxic_Blood", "Cytotoxic_mucosa", "Quiescent_Blood", "REG1A_REG1B_mucosa","Th17_mucosa", "Treg/EMC_Blood", "Treg/Quiescent_Blood", "Treg/Quiescent_mucosa")
new.cluster.ids <- c("CTL blood", "CTL mucosa", "Quiescent blood","REG1A/1B mucosa", "Th17 mucosa", "Effector/Treg blood", "Treg/Quiescent blood", "Treg/Quiescent mucosa")
seuratfile_mucosacells@ident <- plyr::mapvalues(x = seuratfile_mucosacells@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(seuratfile_mucosacells, pt.shape="tissue", pt.size=3, colors.use = c("orange", "magenta","green2", "purple"))
```
tSNE bloodcells, regressed for 'patient'
```
seuratfile_bloodcells<-SetAllIdent(seuratfile_bloodcells, "eight_cell_types")
current.cluster.ids <- c("Cytotoxic_Blood", "Cytotoxic_mucosa", "Quiescent_Blood", "REG1A_REG1B_mucosa","Th17_mucosa", "Treg/EMC_Blood", "Treg/Quiescent_Blood", "Treg/Quiescent_mucosa")
new.cluster.ids <- c("CTL blood", "CTL mucosa", "Quiescent blood","REG1A/1B mucosa", "Th17 mucosa", "Effector/Treg blood", "Treg/Quiescent blood", "Treg/Quiescent mucosa")
seuratfile_bloodcells@ident <- plyr::mapvalues(x = seuratfile_bloodcells@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(seuratfile_bloodcells, pt.size=3, colors.use = c("red", "yellow2", "cyan2", "blue" ) )
```
tSNE celltypes in all tissues
```
allcells_meta<-SetAllIdent(seuratfile_allcells_patientregr, "eight_cell_types")
current.cluster.ids <- c("Cytotoxic_Blood", "Cytotoxic_mucosa", "Quiescent_Blood", "REG1A_REG1B_mucosa","Th17_mucosa", "Treg/EMC_Blood", "Treg/Quiescent_Blood", "Treg/Quiescent_mucosa")
new.cluster.ids <- c("CTL blood", "CTL mucosa", "Quiescent blood","REG1A/1B mucosa", "Th17 mucosa", "Effector/Treg blood", "Treg/Quiescent blood", "Treg/Quiescent mucosa")
x<-allcells_meta
x<-SetAllIdent(x, "eight_cell_types")
Dutch<-c("m", "v")
English<-c("Male", "Female")
x@ident <- plyr::mapvalues(x = x@ident, from = current.cluster.ids, to = new.cluster.ids)
x@meta.data$gender <- plyr::mapvalues(x = x@meta.data$gender, from = Dutch, to = English)
TSNEPlot(x, pt.size=2, colors.use =c("red", "orange", "yellow2", "magenta", "green2", "cyan2", "blue", "purple"), pt.shape="gender")
```
make plots with various cell type markers
```
LIST<-c("CD62L", "CD8ab", "CD45RO")
for(i in LIST){
  FeaturePlot(seuratfile_bloodcells, c(i), cols.use = c( "lightgrey", "royalblue"),pt.size = 1, max.cutoff = "q75",no.legend = T, no.axes = T)
  dev.copy(pdf,width=3, height=2.5,paste0("Bloodcells_", i, ".pdf"))
  dev.off()
}

for(i in LIST){
  FeaturePlot(seuratfile_mucosacells, c(i), cols.use = c( "lightgrey", "royalblue"),pt.size = 1, max.cutoff = "q75",no.legend = T, no.axes = T)
  dev.copy(pdf,width=3, height=2.5,paste0("Mucosacells_", i, ".pdf"))
  dev.off()
}

LIST<-c("CD160", "XCL2", "XCL1", "ITGA1", "FASLG", "CCL5", "IL17A", "IL22" ,"IRF4", "NR4A2", "CREM", "TNFRSF4","FOXP3",  "NOG","CCR7", "SELL","IL6ST", "REG1A", "REG1B", "EPCAM","CD3G","EGF", "LYZ" )
for(i in LIST){
  FeaturePlot(seuratfile_mucosacells, c(i), cols.use = c( "lightgrey", "red"),pt.size = 1, max.cutoff = "q75",no.legend = F, no.axes = T)
  dev.copy(pdf,width=3, height=2.5,paste0("Mucosacells_", i, ".pdf"))
  dev.off()
}

LIST<-c("CCR7","CD160", "GNLY", "GZMH", "GZMB", "CCL4", "CCL5", "TBX21", "ITGAE", "IL32", "PRDM1", "CMTM6", "LTB", "MTA2","HNRNPH1", "TNFRSF4", "IL6ST", "LEF1", "TCF7",  "NOG")
for(i in LIST){
  FeaturePlot(seuratfile_bloodcells, c(i), cols.use = c( "lightgrey", "red"),pt.size = 1, max.cutoff = "q75",no.legend = F, no.axes = T)
  dev.copy(pdf,width=3, height=2.5,paste0("Bloodcells_", i, ".pdf"))
  dev.off()
}
