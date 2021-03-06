library(Seurat)
data<-readRDS("PSC_202002_integrated_v2_noribo.rds")


data@meta.data$cell_cat<-"other"

data@meta.data$cell_cat[data@meta.data$celltypes == "Immature_enterocyte"]<-"epithelial"
data@meta.data$cell_cat[data@meta.data$celltypes == "Cycling_TA"]<-"epithelial"
data@meta.data$cell_cat[data@meta.data$celltypes == "Absorptive_enterocyte"]<-"epithelial"
data@meta.data$cell_cat[data@meta.data$celltypes == "Ribo_TA"]<-"epithelial"
data@meta.data$cell_cat[data@meta.data$celltypes == "Stem"]<-"epithelial"

data@meta.data$cell_cat[data@meta.data$celltypes == "Immature_goblet"]<-"epithelial"
data@meta.data$cell_cat[data@meta.data$celltypes == "BEST4_enterocyte"]<-"epithelial"
data@meta.data$cell_cat[data@meta.data$celltypes == "DUOX2_enterocyte"]<-"epithelial"
data@meta.data$cell_cat[data@meta.data$celltypes == "PLCG2_TA"]<-"epithelial"
data@meta.data$cell_cat[data@meta.data$celltypes == "Absorptive_TA"]<-"epithelial"
data@meta.data$cell_cat[data@meta.data$celltypes == "Tuft"]<-"epithelial"
data@meta.data$cell_cat[data@meta.data$celltypes == "REG_TA"]<-"epithelial"
data@meta.data$cell_cat[data@meta.data$celltypes == "Enteroendocrine"]<-"epithelial"

data@meta.data$cell_cat[data@meta.data$celltypes == "Goblet"]<-"epithelial"

Idents(data)<-"cell_cat"
DefaultAssay(data)<-"RNA"
markers_epi<-FindMarkers(data, subset.ident = "epithelial", group.by="state", test.use = "MAST", ident.1 = "PSC-NI", ident.2 = "HC-NI")
write.csv(markers_epi, "Results/DE/epithelium_pscNI_vs_HC.csv")

markers_epi<-FindMarkers(data, subset.ident = "epithelial", group.by="state", test.use = "MAST", ident.1 = "PSC-NI", ident.2 = "HC-NI")
write.csv(markers_epi, "Results/DE/pithelium_pscNI_vs_HC.csv")

markers_epi<-FindMarkers(data, subset.ident = "epithelial", group.by="state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "HC-NI")
write.csv(markers_epi, "Results/DE/epithelium_pscI_vs_HC.csv")

markers_epi<-FindMarkers(data, subset.ident = "epithelial", group.by="state", test.use = "MAST", ident.1 = "UC-NI", ident.2 = "HC-NI")
write.csv(markers_epi, "Results/DE/epithelium_ucNI_vs_HC.csv")

markers_epi<-FindMarkers(data, subset.ident = "epithelial", group.by="state", test.use = "MAST", ident.1 = "UC-I", ident.2 = "HC-NI")
write.csv(markers_epi, "Results/DE/epithelium_ucI_vs_HC.csv")

markers_epi<-FindMarkers(data, subset.ident = "epithelial", group.by="state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "PSC-NI")
write.csv(markers_epi, "Results/DE/epithelium_pscI_vs_pscNI.csv")

markers_epi<-FindMarkers(data, subset.ident = "epithelial", group.by="state", test.use = "MAST", ident.1 = "UC-I", ident.2 = "UC-NI")
write.csv(markers_epi, "Results/DE/epithelium_ucI_vs_ucNI.csv")
