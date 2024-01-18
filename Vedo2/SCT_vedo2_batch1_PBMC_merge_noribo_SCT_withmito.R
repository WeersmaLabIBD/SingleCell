library(Seurat)
library(ggplot2)



# Load the 4 lanes datasets
lane18_1<-readRDS("/groups/umcg-weersma/tmp01/Frank/batch1_integration/preprocessed_lanes_new/200618_lane1_sct.rds")
lane18_2<-readRDS("/groups/umcg-weersma/tmp01/Frank/batch1_integration/preprocessed_lanes_new/200618_lane2_sct.rds")  
lane18_3<-readRDS("/groups/umcg-weersma/tmp01/Frank/batch1_integration/preprocessed_lanes_new/200618_lane3_sct.rds")
lane18_4<-readRDS("/groups/umcg-weersma/tmp01/Frank/batch1_integration/preprocessed_lanes_new/200618_lane4_sct.rds")

# Merge the dataset
dataset_combined <- merge(lane18_1, y = c(lane18_2, lane18_3, lane18_4), merge.data = TRUE)

# Remove doublets
dataset_combined <- subset(dataset_combined, subset = Final_HTO_status == 'Singlet')

# Normalize after merging (variable features are different for separate objects)
dataset_combined <- SCTransform(dataset_combined,vars.to.regress = c("percent.mt"))

# Dimensional reduction
dataset_combined <- RunPCA(dataset_combined)
dataset_combined <- FindNeighbors(dataset_combined, dims = 1:30)
dataset_combined <- FindClusters(dataset_combined, resolution = 1.2)
dataset_combined <- RunUMAP(dataset_combined, dims = 1:30)

# Check the batch effects
P1 <- DimPlot(megerd_data, group.by = "seurat_clusters", label = TRUE)
P2 <- DimPlot(megerd_data, group.by = "lane", label = FALSE)
P1+P2

# Save the plot
ggsave("megerd_batch_effects.pdf", width = 20, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_merge/plot")


saveRDS(dataset_combined, file = "/groups/umcg-weersma/tmp01/Frank/batch1_merge/final_output/SCT_vedo2_batch1_PBMC_merge_noribo_SCT_withmito.rds")



