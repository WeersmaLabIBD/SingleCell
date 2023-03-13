library(Seurat)

#open Seuratfile
data<-readRDS("/source/ddtx_merged_demultiplexed_clustered_compartment_azi_elmentaiteadultileum")

# filter 60% mito
data<-subset(data, subset=percent.mito < 60)

# exclude T8. patient 3 (wrongly sampled, only recipient cells)
data@meta.data$pt_project_lane<-paste(data@meta.data$lane,data@meta.data$donor_final, sep="_")
data<-subset(data, subset=pt_project_lane!="220504_lane04_UMCGDDtx00005r")

#save
saveRDS(data, "/source/ddtx_merged_demultiplexed_clustered_compartment_azi_elmentaiteadultileum_below60pctmito.rds")
