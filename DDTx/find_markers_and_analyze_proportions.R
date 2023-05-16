## script for marker finding and proportion analysis DDTx project
# May 16 2023
# WTC

# open filtered dataset
library(Seurat)
data<-readRDS("ddtx_merged_demultiplexed_clustered_compartment_azi_elmentaiteadultileum_below60pctmito_withoutEpithelialR.rds")

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
