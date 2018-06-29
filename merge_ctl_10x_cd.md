**Compare cytotoxic cell dataset from 10X, and compare with 'own' CTLs**

# Author: WTC
# Date: 20180629

# Description:
using this script in R, one can merge two datasets from different sources (i.e. a subset of our CD dataset and a publicly available 10X genomics dataset) and explore gene expression differences between those

# Libraries
#
###########################################################################################################################
library(Seurat)
library(Matrix)
#library(Matrix.utils)
library(ggplot2)
#library(pryr)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################

# Name: prepare.lane
# Function: Process all the data from a lane and filter out the doublets found with the DeAnonymizer
# Input:
#   Name 	            Type          Description
#   input.loc       	character     The directory to the lane data
#   deAnonymizer.loc    character     Location to the deAnonymizer output for that lane
#   lane 		  	    numeric       The lane number
#   cdulated         	character     Either "all" if all samples are cdulated, "none" if no samples are cdulated, or a
#									  vector with the sample IDs that are cdulated
#   ctrlulated      	character     Either "all" if all samples are ctrlulated, "none" if no samples are ctrlulated, or a
#									  vector with the sample IDs that are ctrlulated
#
# Output:
#   A list with the data matrix, a vector with the cdulation information, the person associated to each cell, all mitochondrial 
#   genes and the lane number

# download Cytotoxic T cell gene/cell matrix (raw) from https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/cytotoxic_t
CTL_10x<-Read10X(data.dir = "~/Downloads/matrices_mex/hg19/")
CTL_10x<-CreateSeuratObject(raw.data = CTL_10x, min.cells = 3, project = "healthy")


# extract blood Cytotoxic T cells from CD Seurat file
#allcells_meta<-SetAllIdent(allcells_meta, "eight_cell_types")
#CTL_CD<-WhichCells(allcells_meta, "Cytotoxic_Blood")
#CTL_CD<-SubsetData(allcells_meta, CTL_CD)

# filter mito genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = CTL_10x@data), value = TRUE)
percent.mito <- Matrix::colSums(CTL_10x@raw.data[mito.genes, ])/Matrix::colSums(CTL_10x@raw.data)
CTL_10x <- AddMetaData(object = CTL_10x, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = CTL_10x, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# for ngene: 200-2500, percent mito: <5%
CTL_10x <- FilterCells(object = CTL_10x, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

# log normalize
CTL_10x <- NormalizeData(object = CTL_10x, normalization.method = "LogNormalize", scale.factor = 10000)

CTL_10x <- ScaleData(object = CTL_10x, vars.to.regress = c("nUMI", "percent.mito"))

## same for CD
mito.genes <- grep(pattern = "^MT\\.", x = rownames(x = CTL_CD@data), value = TRUE)
percent.mito <- Matrix::colSums(CTL_CD@raw.data[mito.genes, ])/Matrix::colSums(CTL_CD@raw.data)
CTL_CD <- AddMetaData(object = CTL_CD, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = CTL_CD, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# for ngene: 200-2500, percent mito: <5%
CTL_CD <- FilterCells(object = CTL_CD, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

# log normalize
CTL_CD <- NormalizeData(object = CTL_CD, normalization.method = "LogNormalize", scale.factor = 10000)

CTL_CD <- ScaleData(object = CTL_CD, vars.to.regress = c("nUMI", "percent.mito"))



CTL_10x@meta.data$dataset <- "ctrl"
CTL_CD@meta.data$dataset <- "cd"

save(CTL_10x, file="~/Desktop/Single_cell/CTL_10x_CD/CTL_10x.Rda")
save(CTL_CD, file="~/Desktop/Single_cell/CTL_10x_CD/CTL_CD.Rda")


## Combine the cdulated and ctrlulated cells into a single object
ctl.combined <- RunCCA(CTL_10x, CTL_CD, genes.use = intersect(rownames(CTL_10x@data), rownames(CTL_CD@data)), num.cc = 30, scale.data=T)

###############################################################
# Process the data after combining the two datasets
###############################################################

ctl.combined <- RunPCA(object = ctl.combined, pc.genes = ctl.combined@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
ctl.combined = RunTSNE(ctl.combined, dims.use = 1:15, do.fast = T)

## The DimPlot can use either "pca", "tsne" or "cca" for it's reduction.use option and you can add any meta data to the group.by for the colouring
DimPlot(object = ctl.combined, reduction.use = "cca", group.by = "dataset", pt.size = 0.5, do.return = TRUE)
VlnPlot(object = ctl.combined, features.plot = "CC1", group.by = "dataset", do.return = TRUE)

## CCA equivalent of the PCElbowPlot
MetageneBicorPlot(ctl.combined, grouping.var = "dataset", dims.eval = 1:30, display.progress = FALSE)

## Use the CCA reduction to better overlap the two original datasetes
ctl.combined <- AlignSubspace(ctl.combined, reduction.type = "cca", grouping.var = "dataset", dims.align = 1:15)

ctl.combined <- RunTSNE(ctl.combined, reduction.use = "cca.aligned", dims.use = 1:15, do.fast = T)
ctl.combined <- FindClusters(ctl.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:15)

setwd("~/Desktop/Single_cell/CTL_10x_CD/")
save(ctl.combined, file="~/Desktop/Single_cell/CTL_10x_CD/ctl_combined_aligned.Rda")
###############################################################
# Making several tSNE plots, coloured by different meta data
###############################################################
TSNEPlot(ctl.combined, do.return = T, pt.size = 0.5, group.by = "dataset")
TSNEPlot(ctl.combined, do.label = T, do.return = T, pt.size = 0.5)


###############################################################
# Make feature plots for some marker genes
###############################################################
genes.to.use <- c("MS4A1", "CD3E", "LYZ", "GNLY", "CD8A", "PF4")
for (gene in genes.to.use){
	cols <- colorRampPalette(c("lightgrey", "blue"))(100)
	ens.id <- gene
	expression <- ctl.combined@data[ens.id,]
	bins <- .bincode(expression, seq(min(expression)-0.00001,max(expression)+0.00001, length.out=100))
	plot(ctl.combined@dr$tsne@cell.embeddings[,"tSNE_1"], ctl.combined@dr$tsne@cell.embeddings[,"tSNE_2"], col=alpha(cols[bins],0.3), xlab="tSNE_1", ylab="tSNE_2", main=paste0("Expression of ", gene), pch=20)
}


# Name: get.violin.data
# Function: Get the data to make a multi-violin plot
# Input:
#   Name 	            Type          Description
#   seurat 		       	Seurat 	      The seurat object to take the data from
#   genes 		       	Character     The genes for which to make the violin plots
#   HGNC.names 			Character     The HGNC IDs for the genes
#
# Output:
#   The data to plot in the violin plots

get.violin.data <- function(seurat, genes, HGNC.names) {
  output = data.frame(gene = character(0), value= numeric(0), ident = character(0))
  for (i in 1:length(genes)) {
    data.use = data.frame(FetchData(seurat,genes[i]))
    data.use = t(data.use)
    data.melt=data.frame(gene=rep(genes[i], length(seurat@ident)), gene.HGNC=rep(HGNC.names[i], length(seurat@ident)))
    data.melt$value=as.numeric(data.use[1,1:length(seurat@ident)])
    data.melt$id=names(data.use)[1:length(seurat@ident)]
    data.melt$ident=seurat@ident
    noise = rnorm(length(data.melt$value))/100000
    data.melt$value=as.numeric(as.character(data.melt$value))+noise
    output = rbind(output, data.melt)
  }
  return(output)
}

# Name: gg.color.hue
# Function: Create ggplot colors
# Input:
#   Name 	            Type          Description
#   n  				  	numeric 	  The number of colours to make
#
# Output:
#   The ggplot colors

gg.color.hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

###############################################################
# Making the multi-violin plots
###############################################################
ctl.combined<-SetAllIdent(ctl.combined, "dataset")
identities <- ctl.combined@ident
colors.x <- gg.color.hue(length(unique(identities)))
plot.data <- get.violin.data(ctl.combined, c("CD14", "LYZ", "S100A9", "LYN", "ITGAX", "IFITM3",  "CST3", "GZMB", "PRF1", "CD3D", "NKG7", "KLRC1", "CD79A", "SELL", "CCR7", "S100A4", "CD27", "CD8A", "CD8B", "MS4A1"), c("CD14", "LYZ", "S100A9", "LYN", "ITGAX", "IFITM3", "CST3", "GZMB", "PRF1", "CD3D", "NKG7", "KLRC1", "CD79A", "SELL", "CCR7", "S100A4", "CD27", "CD8A", "CD8B", "MS4A1"))

ggplot(plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) + 
  ylab("") + xlab("") +
  coord_flip() + 
  facet_wrap(~ gene.HGNC,scales = "free_x", ncol = length(levels(plot.data$gene.HGNC))) + 
  theme(strip.text.x = element_text(size=18, angle=90),
        axis.text.y = element_text(size=24),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
  scale_x_discrete(limits = rev(levels(plot.data$ident)), position = "top") +
  scale_fill_manual(values = colors.x) #+ theme_set(theme_gray(base_size = 28))

###############################################################
# Compare the new Seurat object with the old Seurat object
###############################################################

# visualize
# add metadata own data
info<-allcells_meta@meta.data
info<-info[,c(7:18,35,38:41)]
info$cellen<-row.names(info)
row.names(info)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
## add location
#extract file that is a copy from @meta.data
CellsMeta = ctl.combined@meta.data
head(CellsMeta)
CellsMeta$cellen=row.names(CellsMeta)
row.names(CellsMeta)=NULL
dim(CellsMeta)
info<-keeping.order(CellsMeta, merge, y=info, by = "cellen", all=FALSE)

CellsMeta = ctl.combined@meta.data
head(CellsMeta)
CellsMeta$cellen=row.names(CellsMeta)
row.names(CellsMeta)=NULL
dim(CellsMeta)
CellsMeta<-keeping.order(CellsMeta, merge, y=info, by = "cellen", all=T)
row.names(CellsMeta)<-CellsMeta$cellen
dim(CellsMeta)
CellsMeta$cellen<-NULL
head(CellsMeta)
ctl.combined <- AddMetaData(ctl.combined, CellsMeta)

ctl.combined<-SetAllIdent(ctl.combined, "dataset")
cd_cells<-WhichCells(ctl.combined, "cd")
healthy_cells<-WhichCells(ctl.combined, "ctrl")
FeaturePlot(ctl.combined, "NKG7", cells.use = cd_cells, pt.size=3)
FeaturePlot(ctl.combined, "NKG7", cells.use = healthy_cells, pt.size=3)

# calculate DE
# DE with genes >0.01%, MAST
ctl.combined_DE_markers = FindAllMarkers(ctl.combined, min.pct = 0.001, only.pos = T, test.use = "MAST")
write.table(ctl.combined_DE_markers, file = "ctl.blood.10x.cd_DEmarkers_0.1percMAST.txt")

# and the significant ones (p<0.05)
ctl.combined_DE_markers<-ctl.combined_DE_markers[ctl.combined_DE_markers$p_val_adj <0.05,]

# top 10 of both sets
library(dplyr)
ctl.combined_DE_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top10 <- ctl.combined_DE_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# histogram of top 10
DoHeatmap(object = ctl.combined, genes.use = top10$gene, slim.col.label = TRUE, remove.key = FALSE)
