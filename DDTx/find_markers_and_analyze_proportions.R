## script for marker finding and proportion analysis DDTx project
# May 16 2023
# WTC

# open filtered dataset
library(Seurat)
data<-readRDS("/source/ddtx_merged_demultiplexed_clustered_compartment_azi_elmentaiteadultileum_below60pctmito_withoutEpithelialR.rds")

#Vlnplots to check donor/recipient cellorigin
library(ggplot2)
VlnPlot(data, "CCR6", assay="RNA", group.by="compartment_final", split.by="donor_recipient") # CCR6 higher in donor immune cells (as expected)
ggsave("CCR6_donor_recipient.png", width = 20, height = 5, dpi = 600)

VlnPlot(data, "SELL", assay="RNA", group.by="compartment_final", split.by="donor_recipient") # CD62L higher in recipient immune cells (as expected)
ggsave("SELL_donor_recipient.png", width = 20, height = 5, dpi = 600)

#find general markers for donor and recipient, all cells and datapoints taken together
Idents(data)<-"donor_recipient"
markers_data <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25)
write.csv(markers_data, "markers_donor_recipient_mito_epi_filtered_data.csv")

#find general markers per cell type, all participants and datapoints taken together
Idents(data)<-"predicted.celltype.elmentaiteadultileum"
markers_data <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25)
write.csv(markers_data, "markers_allcelltypes_mito_epi_filtered_data.csv")

#downsampling to equal cell numbers per sample for cell proportion analysis
library(Seurat)
Idents(data)<-"sample"
data_1949<-subset(data, downsample=1949)
saveRDS(data_1949, "DDTX_sample_1949cells.rds")

#generate CD4 T subsets

## find subsets T cells
Idents(data)<-"predicted.celltype.elmentaiteadultileum"
cd4T<-subset(data,idents="Activated CD4 T")
DimPlot(cd4T)
ggsave("cd4T_before_reclustering.png", width = 5, height = 5, dpi = 600)

DefaultAssay(object = cd4T) <- "RNA"
cd4T <- NormalizeData(cd4T, normalization.method = "LogNormalize", scale.factor = 10000)
cd4T <- NormalizeData(cd4T)
cd4T <- FindVariableFeatures(cd4T, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cd4T), 10)
all.genes <- rownames(cd4T)
cd4T <- ScaleData(cd4T, features = all.genes)
cd4T <- RunPCA(cd4T, features = VariableFeatures(object = cd4T))
ElbowPlot(cd4T)
ggsave("Elbowplot_cd4T_after_reclustering.png", width = 5, height = 5, dpi = 600)

cd4T <- FindNeighbors(cd4T, dims = 1:10)
cd4T <- FindClusters(cd4T, resolution = 0.4)
head(Idents(cd4T), 5)
cd4T <- RunUMAP(cd4T, dims = 1:10)
DimPlot(cd4T, reduction = "umap")
ggsave("cd4T_after_reclustering.png", width = 5, height = 5, dpi = 600)
markers_cd4T <- FindAllMarkers(cd4T, only.pos = TRUE, min.pct = 0.25)
write.csv(markers_cd4T, "markers_cd4T_mito_epi_filtered_data.csv")
FeaturePlot(cd4T, c("IL17A"))
ggsave("cd4T_after_reclustering_il17.png", width = 5, height = 5, dpi = 600)

## markers for different timepoints epithelial cells in patient 3
Idents(data)<-"compartment_final"
epi<-subset(data, ident="Epithelial")
Idents(epi)<-"patient"
epi_pt3<-subset(epi, ident="UMCGDDtx00005")
Idents(epi_pt3)<-"Timepoint_days"
markers_epi_pt3<-FindAllMarkers(epi_pt3, only.pos=T) 
write.csv(markers_epi_pt3, "/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/seurat_preprocess_samples/DE/markers_epi_pt3_timepoints.csv")

epi_pt2<-subset(epi, ident="UMCGDDtx00004")
Idents(epi_pt2)<-"Timepoint_days"
markers_epi_pt2<-FindAllMarkers(epi_pt2, only.pos=T) 
write.csv(markers_epi_pt2, "/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/seurat_preprocess_samples/DE/markers_epi_pt2_timepoints.csv")

epi_pt1<-subset(epi, ident="UMCGDDtx00003")
Idents(epi_pt1)<-"Timepoint_days"
markers_epi_pt1<-FindAllMarkers(epi_pt1, only.pos=T) 
write.csv(markers_epi_pt1, "/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/seurat_preprocess_samples/DE/markers_epi_pt1_timepoints.csv")

## find endothelial markers
Idents(data)<-"major_celltype"
endo<-subset(data, ident="Endothelial")
Idents(endo)<-"patient"
endo_pt3<-subset(endo, ident="UMCGDDtx00005")
Idents(endo_pt3)<-"Timepoint_days"
markers_endo_pt3<-FindAllMarkers(endo_pt3, only.pos=T) 
write.csv(markers_endo_pt3, "/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/seurat_preprocess_samples/DE/markers_endo_pt3_timepoints.csv")

endo_pt2<-subset(endo, ident="UMCGDDtx00004")
Idents(endo_pt2)<-"Timepoint_days"
markers_endo_pt2<-FindAllMarkers(endo_pt2, only.pos=T) 
write.csv(markers_endo_pt2, "/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/seurat_preprocess_samples/DE/markers_endo_pt2_timepoints.csv")

endo_pt1<-subset(endo, ident="UMCGDDtx00003")
Idents(endo_pt1)<-"Timepoint_days"
markers_endo_pt1<-FindAllMarkers(endo_pt1, only.pos=T) 
write.csv(markers_endo_pt1, "/groups/umcg-weersma/tmp01/projects/ddtx/ongoing/seurat_preprocess_samples/DE/markers_endo_pt1_timepoints.csv")





# check markers on gutcellatlas


