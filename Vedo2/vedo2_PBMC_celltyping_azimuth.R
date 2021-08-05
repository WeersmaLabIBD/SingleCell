# be sure to on the cluster, do these
# ml R
# ml HDF5

# All libraries required
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

# Load the reference Seurat object
reference <- LoadH5Seurat("/groups/umcg-weersma/tmp01/projects/vedo2_emilia/ongoing/cell-type-classification/azimuth/reference/pbmc_multimodal.h5se
urat")

# Read the PBMC object
PBMC <- readRDS("/groups/umcg-weersma/tmp01/Frank/batch1_integration/resoulution_1.2/vedo2_batch1_PBMC_integrated_noribo_withmito.rds")

# Apply the same normalization strategy as used in the reference
PBMC <- SCTransform(PBMC)
# Do principal component analysis
PBMC <- RunPCA(PBMC)

# Find transfer anchors between the reference and the query, the query is PBMC dataset
anchors <- FindTransferAnchors(
  reference = reference,
  query = PBMC,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,
  recompute.residuals = FALSE
)

# Map the query to the reference, this will project PBMC dataset onto the reference, and copy the labels with a certain amount of certainty (the score)
PBMC <- MapQuery(
  anchorset = anchors,
  query = PBMC,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    celltype.l3 = "celltype.l3",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

# Save the PBMC celltyping object
saveRDS(PBMC, file = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/vedo2_batch1_PBMC_integrated_Azimuth_celltyping.rds")



# Normalization
DefaultAssay(data) <- "RNA"
data <- NormalizeData(data, verbose = TRUE)



# The object now has the UMAP that was projected in 'ref.umap', and your own UMAP in 'umap', as dimensional reductions
DimPlot(data, reduction = "umap", group.by = "seurat_clusters", label = T) # this will show the clusters and unprojected UMAP
ggsave("unprojected_UMAP.pdf", width = 10, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/DimPlot/unprojected_UMAP/")

DimPlot(data, reduction = 'ref.umap', group.by = 'predicted.celltype.l1', label = TRUE) # this will show the predicted l1 cell type and the projected UMAP
ggsave("predicted_l1_cell_type.pdf", width = 10, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/DimPlot")

DimPlot(data, reduction = 'ref.umap', group.by = 'predicted.celltype.l2', label = FALSE) # this will show the predicted l1 cell type and the projected UMAP
ggsave("predicted_l2_cell_type.pdf", width = 10, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/DimPlot")


# You can do this to see if you find any 'confetti' that might be doublets, contamination, or an unknown cell type
P1 <- DimPlot(data, reduction = 'umap', group.by = 'predicted.celltype.l1', label = TRUE, label.box = T, cols = brewer, repel = T) # here you can check the azimuth prediction against your own UMAP and clustering
ggsave("predicted_clusters_l1.pdf", width = 10, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/DimPlot/predicted_celltype")

P2 <- DimPlot(data, reduction = 'umap', group.by = 'predicted.celltype.l2', label = TRUE, label.box = T, cols = brewer, repel = T)
ggsave("predicted_clusters_l2.pdf", width = 12, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/DimPlot/predicted_celltype")


## Check the prediction score

# Copy the l1 column in the metadata to a new column
data@meta.data$predicted.celltype.l1.wreject <- as.character(data@meta.data$predicted.celltype.l1)
# For all rows in metadata, where l1 score < 0.8, change the l1.reject column to NA
data@meta.data[data@meta.data$predicted.celltype.l1.score < 0.8, 'predicted.celltype.l1.wreject'] <- NA
# Plot the l1 cell type without changes
P3 <- DimPlot(data, reduction = 'umap', group.by = 'predicted.celltype.l1', label = TRUE, label.box = T, repel = T)
# Plot the l1 cell type, with the rejections (<0.8 turned into NA)
P4 <- DimPlot(data, reduction = 'umap', group.by = 'predicted.celltype.l1.wreject')
P3 + P4
# Save the plot
ggsave('predicted_clusters_l1+wreject.pdf', width = 20, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/DimPlot/check_azimuth_score/")

# Copy the l2 column in the metadata to a new column
data@meta.data$predicted.celltype.l2.wreject <- as.character(data@meta.data$predicted.celltype.l2)
# For all rows in metadata, where l1 score < 0.8, change the l1.reject column to NA
data@meta.data[data@meta.data$predicted.celltype.l2.score < 0.8, 'predicted.celltype.l2.wreject'] <- NA
# Plot the l2 cell type without changes
P5 <- DimPlot(data, reduction = 'umap', group.by = 'predicted.celltype.l2', label = TRUE, label.box = T, repel = T)
# Plot the l2 cell type, with the rejections (<0.8 turned into NA)
P6 <- DimPlot(data, reduction = 'umap', group.by = 'predicted.celltype.l2.wreject')
P5 + P6
# Save the plot
ggsave('predicted_clusters_l2_wreject.pdf', width = 20, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/DimPlot/check_azimuth_score/")


# Plot score FeaturePlot
FeaturePlot(data, reduction = "umap", features = "predicted.celltype.l1.score")
ggsave('l1_score.pdf', width = 10, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/FeaturePlot/l1_score/")

FeaturePlot(data, reduction = "umap", features = "predicted.celltype.l2.score")
ggsave('l2_score.pdf', width = 10, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/FeaturePlot/l1_score/")

## There are a lot of cells with a bad prediction score in l2 prediction, find the reason.

# 1. Do they have a high mitochondrial percent per cells?
P1 <- FeaturePlot(data, reduction = "umap", features = c('percent.mt'))
P2 <- DimPlot(data, reduction = 'umap', group.by = 'predicted.celltype.l2.wreject')
cowplot::plot_grid(plotlist = list(P2, P1), ncol = 2, nrow = 1)
ggsave("check_mt.pdf", width = 20, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/DimPlot/check_mt")
# Check the plot, the answer is no, because there is the same mitochondrial percent between the cells with better score and that with bad score.

# 2. Do they express marker genes of multiple cells types?
CD4exp <- FeaturePlot(data, reduction = 'umap', features = c('CD4'))
P2 <- DimPlot(data, reduction = 'umap', group.by = 'predicted.celltype.l2.wreject')
P2+CD4exp
ggsave("CD4.pdf", width = 20, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/DimPlot/check_doublets")
MS4A1exp <- FeaturePlot(data, reduction = 'umap', features = c('MS4A1'))
P2+ MS4A1exp
ggsave("MS4A1.pdf", width = 20, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/DimPlot/check_doublets")
CD8Aexp <- FeaturePlot(data, reduction = 'umap', features = c('CD8A'))
P2+ CD8Aexp
ggsave("CD8A.pdf", width = 20, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/DimPlot/check_doublets")


# 3. Do they have low expression?
P3 <- FeaturePlot(data, reduction = "umap", features = c('nFeature_RNA'))
P2+P3
ggsave("check_low_expression.pdf", width = 20, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/DimPlot/check_low_expression")
# Check the plot, the answer is no, because there is the same gene expression between the cells with better score and that with bad score.

# 4. Is the more coarse prediction better?
# Yes, the l1 prediction is better than l2 prediction

# 5. Do they express genes the better scored cells don't?

### CD4+ T cells
seurat_object_cd4 <- subset(data, predicted.celltype.l2 %in% c('CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 'CD4 TCT', 'CD4 TEM'))
seurat_object_cd4@meta.data$bad_score <- NA
seurat_object_cd4@meta.data$bad_score <- seurat_object_cd4@meta.data$predicted.celltype.l2.score < 0.8
Idents(seurat_object_cd4) <- 'bad_score'
head(seurat_object_cd4@active.ident)
cd4t_bad_vs_good <- FindMarkers(seurat_object_cd4, ident.1 = T, ident.2 = F, test.use = 'MAST')
head(cd4t_bad_vs_good)
table(seurat_object_cd4@meta.data$predicted.celltype.l2)
table(seurat_object_cd4@meta.data[seurat_object_cd4@meta.data$bad_score == T, ]$predicted.celltype.l2)
table(seurat_object_cd4@meta.data[seurat_object_cd4@meta.data$bad_score == F, ]$predicted.celltype.l2)

mean(seurat_object_cd4@meta.data[seurat_object_cd4@meta.data$predicted.celltype.l2 == 'CD4 TEM' & seurat_object_cd4@meta.data$bad_score == F, ]$predicted.celltype.l2.score)
mean(seurat_object_cd4@meta.data[seurat_object_cd4@meta.data$predicted.celltype.l2 == 'CD4 TEM' & seurat_object_cd4@meta.data$bad_score == T, ]$predicted.celltype.l2.score)

rownames(cd4t_bad_vs_good[cd4t_bad_vs_good$p_val_adj < 0.05, ])

rm(seurat_object_cd4, cd4t_bad_vs_good)

### CD8+ T cells
seurat_object_cd8 <- subset(data, predicted.celltype.l2 %in% c('CD8 Naive', 'CD8 TCM', 'CD8 TEM'))
seurat_object_cd8@meta.data$bad_score <- NA
seurat_object_cd8@meta.data$bad_score <- seurat_object_cd8@meta.data$predicted.celltype.l2.score < 0.8
Idents(seurat_object_cd8) <- 'bad_score'
head(seurat_object_cd8@active.ident)
cd8t_bad_vs_good <- FindMarkers(seurat_object_cd8, ident.1 = T, ident.2 = F, test.use = 'MAST')
head(cd8t_bad_vs_good)
table(seurat_object_cd8@meta.data$predicted.celltype.l2)
table(seurat_object_cd8@meta.data[seurat_object_cd8@meta.data$bad_score == T, ]$predicted.celltype.l2)
table(seurat_object_cd8@meta.data[seurat_object_cd8@meta.data$bad_score == F, ]$predicted.celltype.l2)

mean(seurat_object_cd8@meta.data[seurat_object_cd8@meta.data$predicted.celltype.l2 == 'CD8 Naive' & seurat_object_cd8@meta.data$bad_score == F, ]$predicted.celltype.l2.score)
mean(seurat_object_cd8@meta.data[seurat_object_cd8@meta.data$predicted.celltype.l2 == 'CD8 Naive' & seurat_object_cd8@meta.data$bad_score == T, ]$predicted.celltype.l2.score)

rownames(cd8t_bad_vs_good[cd8t_bad_vs_good$p_val_adj < 0.05, ])

rm(seurat_object_cd8, cd8t_bad_vs_good)








# Next are some methods to check specific markers for the Azimuth predictions
# Plot celltype markers l2 predictions
plot_celltype_markers_l2 <- function(seurat_object, celltype_marker_genes=c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34"), assay = "RNA", slot="data", plot_dir = "./", cluster_name='predicted.celltype.l2', reduction='umap'){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot, reduction = reduction)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      #p$data$clusters <- seurat_object@meta.data[[cluster_name]]
      #LabelClusters(plot = p, id = "clusters")
      ggsave(paste(plot_dir,gene,".png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}


# B intermediate
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('MS4A1', 'TNFRSF13B', 'IGHM', 'IGHD', 'AIM2', 'CD79A', 'LINC01857', 'RALGPS2', 'BANK1', 'CD79B'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/B intermediate/")

# B memory
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('MS4A1', 'COCH', 'AIM2', 'BANK1', 'SSPN', 'CD79A', 'TEX9', 'RALGPS2', 'TNFRSF13C', 'LINC01781'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/B memory/")

# B naive
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('IGHM', 'IGHD', 'CD79A', 'IL4R', 'MS4A1', 'CXCR4', 'BTG1', 'TCL1A', 'CD79B', 'YBX3'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/B naive/")

# Plasmablast
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('IGHA2', 'MZB1', 'TNFRSF17', 'DERL3', 'TXNDC5', 'TNFRSF13B', 'POU2AF1', 'CPNE5', 'HRASLS2', 'NT5DC2'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/Plasmablast/")

# CD4 CTL
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('GZMH', 'CD4', 'FGFBP2', 'ITGB1', 'GZMA', 'CST7', 'GNLY', 'B2M', 'IL32', 'NKG7'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/CD4 CTL/")

# CD4 Naive
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('TCF7', 'CD4', 'CCR7', 'IL7R', 'FHIT', 'LEF1', 'MAL', 'NOSIP', 'LDHB', 'PIK3IP1'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/CD4 Naive/")

# CD4 Proliferating
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('MKI67', 'TOP2A', 'PCLAF', 'CENPF', 'TYMS', 'NUSAP1', 'ASPM', 'PTTG1', 'TPX2', 'RRM2'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/CD4 Proliferating/")

# CD4 TCM
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('IL7R', 'TMSB10', 'CD4', 'ITGB1', 'LTB', 'TRAC', 'AQP3', 'LDHB', 'IL32', 'MAL'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/CD4 TCM/")

# CD4 TEM
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('IL7R', 'CCL5', 'FYB1', 'GZMK', 'IL32', 'GZMA', 'KLRB1', 'TRAC', 'LTB', 'AQP3'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/CD4 TEM/")

# Treg
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('RTKN2', 'FOXP3', 'AC133644.2', 'CD4', 'IL2RA', 'TIGIT', 'CTLA4', 'FCRL3', 'LAIR2', 'IKZF2'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/Treg/")

#CD8 Naive
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('CD8B', 'S100B', 'CCR7', 'RGS10', 'NOSIP', 'LINC02446', 'LEF1', 'CRTAM', 'CD8A', 'OXNAD1'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/CD8 Naive/")

# CD8 Proliferating
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('MKI67', 'CD8B', 'TYMS', 'TRAC', 'PCLAF', 'CD3D', 'CLSPN', 'CD3G', 'TK1', 'RRM2'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/CD8 Proliferating/")

# CD8 TCM
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('CD8B', 'ANXA1', 'CD8A', 'KRT1', 'LINC02446', 'YBX3', 'IL7R', 'TRAC', 'NELL2', 'LDHB'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/CD8 TCM/")

# CD8 TEM
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('CCL5', 'GZMH', 'CD8A', 'TRAC', 'KLRD1', 'NKG7', 'GZMK', 'CST7', 'CD8B', 'TRGC2'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/CD8 TEM/")

# ASDC
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('PPP1R14A', 'LILRA4', 'AXL', 'IL3RA', 'SCT', 'SCN9A', 'LGMN', 'DNASE1L3', 'CLEC4C', 'GAS6'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/ASDC/")

# cDC1
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('CLEC9A', 'DNASE1L3', 'C1orf54', 'IDO1', 'CLNK', 'CADM1', 'FLT3', 'ENPP1', 'XCR1', 'NDRG2'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/cDC1/")

# cDC2
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('FCER1A', 'HLA-DQA1', 'CLEC10A', 'CD1C', 'ENHO', 'PLD4', 'GSN', 'SLC38A1', 'NDRG2', 'AFF3'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/cDC2/")

# pDC
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('ITM2C', 'PLD4', 'SERPINF1', 'LILRA4', 'IL3RA', 'TPM2', 'MZB1', 'SPIB', 'IRF4', 'SMPD3'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/pDC/")

# CD14 Mono
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('S100A9', 'CTSS', 'S100A8', 'LYZ', 'VCAN', 'S100A12', 'IL1B', 'CD14', 'G0S2', 'FCN1'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/CD14 Mono/")

# CD16 Mono
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('CDKN1C', 'FCGR3A', 'PTPRC', 'LST1', 'IER5', 'MS4A7', 'RHOC', 'IFITM3', 'AIF1', 'HES4'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/CD16 Mono/")

# NK
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('GNLY', 'TYROBP', 'NKG7', 'FCER1G', 'GZMB', 'TRDC', 'PRF1', 'FGFBP2', 'SPON2', 'KLRF1'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/NK/")

# NK Proliferating
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('MKI67', 'KLRF1', 'TYMS', 'TRDC', 'TOP2A', 'FCER1G', 'PCLAF', 'CD247', 'CLSPN', 'ASPM'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/NK Proliferating/")

# NK_CD56bright
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('XCL2', 'FCER1G', 'SPINK2', 'TRDC', 'KLRC1', 'XCL1', 'SPTSSB', 'PPP1R9A', 'NCAM1', 'TNFRSF11A'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/NK_CD56bright/")

# Eryth
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('HBD', 'HBM', 'AHSP', 'ALAS2', 'CA1', 'SLC4A1', 'IFIT1B', 'TRIM58', 'SELENBP1', 'TMCC2'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/Eryth/")

# HSPC
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('SPINK2', 'PRSS57', 'CYTL1', 'EGFL7', 'GATA2', 'CD34', 'SMIM24', 'AVP', 'MYB', 'LAPTM4B'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/HSPC/")

# ILC
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('KIT', 'TRDC', 'TTLL10', 'LINC01229', 'SOX4', 'KLRB1', 'TNFRSF18', 'TNFRSF4', 'IL1R1', 'HPGDS'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/ILC/")

# Platelet
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('PPBP', 'PF4', 'NRGN', 'GNG11', 'CAVIN2', 'TUBB1', 'CLU', 'HIST1H2AC', 'RGS18', 'GP9'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/Platelet/")

# dnT
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('PTPN3', 'MIR4422HG', 'NUCB2', 'CAV1', 'DTHD1', 'GZMA', 'MYB', 'FXYD2', 'GZMK', 'AC004585.1'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/specific_markers_plots/dnT/")

# gdT
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('TRDC', 'TRGC1', 'TRGC2', 'KLRC1', 'NKG7', 'TRDV2', 'CD7', 'TRGV9', 'KLRD1', 'KLRG1'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/gdT/")

#MAIT
plot_celltype_markers_l2(seurat_object = data, celltype_marker_genes=c('KLRB1', 'NKG7', 'GZMK', 'IL7R', 'SLC4A10', 'GZMA', 'CXCR6', 'PRSS35', 'RBM24', 'NCR3'), plot_dir = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/markers_plots/MAIT/")



# create violin plots l2 predictions
plot_celltype_violins_l2 <- function(seurat_object, celltype_marker_genes=c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34"), assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # plot per gene
  for(gene in celltype_marker_genes){
    # only plot if the gene is present
    if(gene %in% rownames(seurat_object)){
      # plot the violins and save
      VlnPlot(seurat_object, features = c(gene), group.by = "predicted.celltype.l2", assay=assay, slot=slot)
      ggsave(paste(plot_dir,gene,"_violin.png", sep=""), dpi=600, width=10, height=10)
      VlnPlot(seurat_object, features = c(gene), group.by = "predicted.celltype.l2", assay=assay, slot=slot, pt.size = 0)
      ggsave(paste(plot_dir,gene,"_violin_nodot.png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
  # we don't want to error because of missing genes
  genes_to_plot <- intersect(rownames(seurat_object), celltype_marker_genes)
  # plot per cluster
  for(ident in unique(as.character(seurat_object@meta.data$predicted.celltype.l2))){
    # plot the violins and save
    Idents(seurat_object) <- 'predicted.celltype.l2'
    vln <- VlnPlot(seurat_object, features = genes_to_plot, idents = c(ident), assay=assay, slot=slot)
    ggsave(paste(plot_dir,ident,"_violin.png", sep=""), dpi=600, width=20, height=20, plot=vln)
  }
}
  
plot_celltype_violins_l2(seurat_object = data, plot_dir ="/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/VlnPlot/")
  
  
# run findmarkers method, to check the difference between two clusters or two cell types etc.
run_findmarkers <- function(seurat_object, output_loc_final, ident_to_set, ident.1, ident.2=NULL, assay='RNA', min.pct=0.1, logfc.threshold=0.25, latent.vars=NULL, features=NULL){
  # set the assay we wish to analyse
  DefaultAssay(seurat_object) <- assay
  # a hard column name is required for the model.matrix step, so let's add that
  seurat_object <- AddMetaData(seurat_object, seurat_object@meta.data[ident_to_set], 'condition')
  Idents(seurat_object) <- seurat_object@meta.data$condition
  # create an output directory
  #output_loc_final <- paste(output_loc, ident.1, ident.2, '.tsv', sep = '')
  # perform MAST if possible
  result <- NULL
  # we need to 'try' here, as MAST ungraciously exits if it cannot perform
  try({
    result <- FindMarkers(object = seurat_object, ident.1 = ident.1, ident.2 = ident.2, test.use = 'MAST', min.pct = min.pct, assay = assay, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  })
  # write the actual table if possible
  if(is.null(result)){
    print(paste('nothing to compare due to threshold for', output_loc_final))
  } else{
    write.table(result, output_loc_final, sep = '\t')
  }
}

run_findmarkers(data, output_loc_final = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/marker_genes/CD4TCM_vs_CD14Mono.tsv", ident_to_set = "predicted.celltype.l2", ident.1 = "CD4 TCM", ident.2 = "CD14 Mono")





create_confusion_matrix <- function(assignment_table, truth_column, prediction_column, truth_column_label=NULL, prediction_column_label=NULL){
  # init the table
  confusion_table <- NULL
  # check each truth
  for(truth in unique(assignment_table[[truth_column]])){
    # get these truths
    truth_rows <- assignment_table[assignment_table[[truth_column]] == truth, ]
    # check now many have this truth
    this_truth_number <- nrow(truth_rows)
    # check what was predicted for these truths
    #for(prediction in unique(truth_rows[[prediction_column]])){
    for(prediction in unique(assignment_table[[prediction_column]])){
      # check the number of this prediction
      this_prediction_number <- nrow(truth_rows[truth_rows[[prediction_column]] == prediction, ])
      # init variable
      fraction <- NULL
      # we can only calculate a fraction if the result is not zero
      if(this_prediction_number > 0){
        # calculate the fraction
        fraction <- this_prediction_number / this_truth_number
      }
      # otherwise we just set it to zero
      else{
        fraction <- 0
      }
      # turn into row
      this_row <- data.frame(truth=c(truth), prediction=c(prediction), freq=c(fraction), stringsAsFactors = F)
      # add this entry to the dataframe
      if(is.null(confusion_table)){
        confusion_table <- this_row
      }
      else{
        confusion_table <- rbind(confusion_table, this_row)
      }
    }
  }
  # round the frequency off to a sensible cutoff
  confusion_table$freq <- round(confusion_table$freq, digits=2)
  # turn into plot
  p <- ggplot(data=confusion_table, aes(x=truth, y=prediction, fill=freq)) + geom_tile() + scale_fill_gradient(low='red', high='blue') + geom_text(aes(label=freq))
  # some options
  if(!is.null(truth_column_label)){
    p <- p + xlab(truth_column_label)
  }
  if(!is.null(prediction_column_label)){
    p <- p + ylab(prediction_column_label)
  }
  return(p)
}

create_confusion_matrix(data@meta.data, 'seurat_clusters', 'predicted.celltype.l1', prediction_column_label = 'l2 predict', truth_column_label = 'seurat cluster')

ggsave('matrix1', width = 10, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/matrix/")

create_confusion_matrix(data@meta.data, 'seurat_clusters', 'predicted.celltype.l2')
ggsave('matrix2', width = 10, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/matrix/")

create_confusion_matrix(data@meta.data, 'predicted.celltype.l2', 'predicted.celltype.l1')
ggsave('matrix3', width = 10, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/matrix/")








