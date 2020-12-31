## comparison rectum vs sigmoid
library(Seurat)
library(ggplot2)

options(future.globals.maxSize = 150000 * 1024^2)

## load immune, epithelial and mesenchyme set from Smillie et al (2019)
load("imm_smillie_seuratfile.Rd")
load("epi_smillie_seuratfile.Rd")
load("fibro_smillie_seuratfile.Rd")

UpdateSeuratObject(imm_smillie_meta)
UpdateSeuratObject(fibro_smillie_meta)
UpdateSeuratObject(epi_smillie_meta)

imm_smillie_meta@meta.data$dataset<-"imm_smillie_meta"
epi_smillie_meta@meta.data$dataset<-"epi_smillie_meta"
fibro_smillie_meta@meta.data$dataset<-"fibro_smillie_meta"

imm_smillie_meta@meta.data$method<-"splitcollagenase"
epi_smillie_meta@meta.data$method<-"splitcollagenase"
fibro_smillie_meta@meta.data$method<-"splitcollagenase"

imm_smillie_meta[["percent.ribo"]] <- PercentageFeatureSet(imm_smillie_meta, pattern = "^RPL|^RPS")
epi_smillie_meta[["percent.ribo"]] <- PercentageFeatureSet(epi_smillie_meta, pattern = "^RPL|^RPS")
fibro_smillie_meta[["percent.ribo"]] <- PercentageFeatureSet(fibro_smillie_meta, pattern = "^RPL|^RPS")

Idents(imm_smillie_meta)<-"Subject"
imm_rectum_and_sigmoid= subset(imm_smillie_meta, ident= c("N12","N26","N44","N661","N50","N106","N49", "N52", "N111", "N58", "N24"))

Idents(fibro_smillie_meta)<-"Subject"
fibro_rectum_and_sigmoid= subset(fibro_smillie_meta, ident= c("N12","N26","N44","N661","N50","N106","N49", "N52", "N111", "N58", "N24"))

Idents(epi_smillie_meta)<-"Subject"
epi_rectum_and_sigmoid= subset(epi_smillie_meta, ident= c("N12","N26","N44","N661","N50","N106","N49", "N52", "N111", "N58", "N24"))

# select sigmoid UC samples 
sigmoid<-c("N12","N26","N44","N661","N50")
# select rectum UC samples
rectum<-c("N106","N49", "N52", "N111", "N58", "N24")

imm_rectum_and_sigmoid@meta.data$location<-imm_rectum_and_sigmoid@meta.data$Subject
imm_rectum_and_sigmoid@meta.data$location<-droplevels(imm_rectum_and_sigmoid@meta.data$location)
imm_rectum_and_sigmoid@meta.data$location<-as.character(imm_rectum_and_sigmoid@meta.data$location)


for(a in rectum){
for (i in 1:nrow(imm_rectum_and_sigmoid@meta.data)){ 
  if (imm_rectum_and_sigmoid@meta.data[i,17] == a){
    imm_rectum_and_sigmoid@meta.data[i,17]<-"rectum"
  } else {
    imm_rectum_and_sigmoid@meta.data[i,17]<-imm_rectum_and_sigmoid@meta.data[i,17]
  }
}}

for(b in sigmoid){
for (i in 1:nrow(imm_rectum_and_sigmoid@meta.data)){ 
  if (imm_rectum_and_sigmoid@meta.data[i,17] == b){
    imm_rectum_and_sigmoid@meta.data[i,17]<-"sigmoid"
  } else {
    imm_rectum_and_sigmoid@meta.data[i,17]<-imm_rectum_and_sigmoid@meta.data[i,17]
  }
}}

fibro_rectum_and_sigmoid@meta.data$location<-fibro_rectum_and_sigmoid@meta.data$Subject
fibro_rectum_and_sigmoid@meta.data$location<-droplevels(fibro_rectum_and_sigmoid@meta.data$location)
fibro_rectum_and_sigmoid@meta.data$location<-as.character(fibro_rectum_and_sigmoid@meta.data$location)


for(a in rectum){
for (i in 1:nrow(fibro_rectum_and_sigmoid@meta.data)){ 
  if (fibro_rectum_and_sigmoid@meta.data[i,18] == a){
    fibro_rectum_and_sigmoid@meta.data[i,18]<-"rectum"
  } else {
    fibro_rectum_and_sigmoid@meta.data[i,18]<-fibro_rectum_and_sigmoid@meta.data[i,18]
  }
}}

for(b in sigmoid){
for (i in 1:nrow(fibro_rectum_and_sigmoid@meta.data)){ 
  if (fibro_rectum_and_sigmoid@meta.data[i,18] == b){
    fibro_rectum_and_sigmoid@meta.data[i,18]<-"sigmoid"
  } else {
    fibro_rectum_and_sigmoid@meta.data[i,18]<-fibro_rectum_and_sigmoid@meta.data[i,18]
  }
}}

table(fibro_rectum_and_sigmoid@meta.data$location)

epi_rectum_and_sigmoid@meta.data$location<-epi_rectum_and_sigmoid@meta.data$Subject
epi_rectum_and_sigmoid@meta.data$location<-droplevels(epi_rectum_and_sigmoid@meta.data$location)
epi_rectum_and_sigmoid@meta.data$location<-as.character(epi_rectum_and_sigmoid@meta.data$location)


for(a in rectum){
for (i in 1:nrow(epi_rectum_and_sigmoid@meta.data)){ 
  if (epi_rectum_and_sigmoid@meta.data[i,16] == a){
    epi_rectum_and_sigmoid@meta.data[i,16]<-"rectum"
  } else {
    epi_rectum_and_sigmoid@meta.data[i,16]<-epi_rectum_and_sigmoid@meta.data[i,16]
  }
}}

for(b in sigmoid){
for (i in 1:nrow(epi_rectum_and_sigmoid@meta.data)){ 
  if (epi_rectum_and_sigmoid@meta.data[i,16] == b){
    epi_rectum_and_sigmoid@meta.data[i,16]<-"sigmoid"
  } else {
    epi_rectum_and_sigmoid@meta.data[i,16]<-epi_rectum_and_sigmoid@meta.data[i,16]
  }
}}

table(epi_rectum_and_sigmoid@meta.data$location)

alldata_202002_integration_list <- list(epi_rectum_and_sigmoid, fibro_rectum_and_sigmoid, imm_rectum_and_sigmoid)

saveRDS(alldata_202002_integration_list, "integration_smillie_rectum_sigmoid.rds")

alldata_202002_integration_list<-readRDS("integration_smillie_rectum_sigmoid.rds")
for (i in 1:length(alldata_202002_integration_list)){
alldata_202002_integration_list[[i]] <- SCTransform(alldata_202002_integration_list[[i]],vars.to.regress = c("percent.mito", "percent.ribo"), verbose = FALSE)
}

#print("loaded")
features <- SelectIntegrationFeatures(object.list = alldata_202002_integration_list, nfeatures = 3000)
#print("features_selected")
alldata_202002_integration_list <- PrepSCTIntegration(object.list = alldata_202002_integration_list, anchor.features = features)
#print("integration_prepped")
saveRDS(alldata_202002_integration_list, "integration_smillie_rectum_sigmoid.rds")

anchors <- FindIntegrationAnchors(object.list = alldata_202002_integration_list, anchor.features = features, normalization.method = "SCT")
#print("achors_found")
rm(alldata_202002_integration_list)
alldata.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
#print("integrated")
rm(anchors)
alldata.integrated <- RunPCA(alldata.integrated, verbose = FALSE)
#print("pca2_ran")
alldata.integrated <- RunUMAP(alldata.integrated, dims = 1:30)
#print("umap_ran")
alldata.integrated <- FindNeighbors(alldata.integrated, dims = 1:30)
alldata.integrated <- FindClusters(alldata.integrated, resolution = 0.5)

DefaultAssay(alldata.integrated) <- "RNA"
alldata.integrated <- NormalizeData(alldata.integrated, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

saveRDS(alldata.integrated, "rectum_and_sigmoid_integrated.rds")

methodspaper<-readRDS("rectum_and_sigmoid_integrated.rds")

# group cell types Smillie
methodspaper@meta.data$Predicted_all<-methodspaper@meta.data$Cluster
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Best4+ Enterocytes"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Enterocytes"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Enterocyte Progenitors"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Immature Enterocytes 1"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Immature Enterocytes 2"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "TA 1"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "TA 2"]<-"epithelium"

# epi general
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "M cells"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Cycling TA"]<-"epithelium"

# epi secretory
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Goblet"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Enteroendocrine"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Immature Goblet"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Secretory TA"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Stem"]<-"epithelium"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Tuft"]<-"epithelium"


#Tcell
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD4+ Activated Fos-lo"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD4+ Activated Fos-hi"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD4+ Memory"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD4+ PD1+"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD8+ IELs"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD8+ IL17+"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD8+ LP"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "MT-hi"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Tregs"]<-"Tcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Cycling T"]<-"Tcell"


# ILS?NK
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "ILCs"]<-"Lymphoid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "NKs"]<-"Lymphoid"

#B
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Cycling B"]<-"Bcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "GC"]<-"Bcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Plasma"]<-"Bcell"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Follicular"]<-"Bcell"

#Myeloid
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD69+ Mast"]<-"myeloid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "CD69- Mast"]<-"myeloid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Inflammatory Monocytes"]<-"myeloid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Macrophages"]<-"myeloid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Cycling Monocytes"]<-"myeloid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "DC1"]<-"myeloid"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "DC2"]<-"myeloid"

#mesenchymal


methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Glia"]<-"Glia"



# fibroblasts
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Inflammatory Fibroblasts"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "RSPO3+"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT2B+ Fos-hi"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT2B+ Fos-lo 1"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT2B+ Fos-lo 2"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT5B+ 1"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT5B+ 2"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "WNT2B+ Fos-lo 1"]<-"fibroblast"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Myofibroblasts"]<-"fibroblast"

#endothelial
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Post-capillary Venules"]<-"endothelial"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Pericytes"]<-"endothelial"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Endothelial"]<-"endothelial"
methodspaper@meta.data$Predicted_all[methodspaper@meta.data$Cluster == "Microvascular"]<-"endothelial"

metadata<-methodspaper@meta.data


write.csv(metadata, "metadata_rectum_sigmoid.csv")

metadata_rectum_sigmoid<-read.csv("metadata_rectum_sigmoid.csv")
metadata_rectum<-metadata_rectum_sigmoid[metadata_rectum_sigmoid$location == "rectum",]
metadata_sigmoid<-metadata_rectum_sigmoid[metadata_rectum_sigmoid$location == "sigmoid",]
metadata_I<-metadata_rectum_sigmoid[metadata_rectum_sigmoid$Health == "Inflamed",]
metadata_NI<-metadata_rectum_sigmoid[metadata_rectum_sigmoid$Health == "Non-inflamed",]

metadata_rectum$sample<-paste0(metadata_rectum$Subject, metadata_rectum$Health)
metadata_sigmoid$sample<-paste0(metadata_sigmoid$Subject, metadata_sigmoid$Health)

proportion_celltypes_rectum<-data.frame(prop.table(table(metadata_rectum$sample, metadata_rectum$Health, metadata_rectum$Predicted_all), 1))
proportion_celltypes_sigmoid<-data.frame(prop.table(table(metadata_sigmoid$sample, metadata_sigmoid$Health, metadata_sigmoid$Predicted_all), 1))
proportion_celltypes_I<-data.frame(prop.table(table(metadata_I$Subject, metadata_I$location, metadata_I$Predicted_all), 1))
proportion_celltypes_NI<-data.frame(prop.table(table(metadata_NI$Subject, metadata_NI$location, metadata_NI$Predicted_all), 1))

proportion_celltypes_rectum<-proportion_celltypes_rectum[proportion_celltypes_rectum$Freq > 0,]
proportion_celltypes_sigmoid<-proportion_celltypes_sigmoid[proportion_celltypes_sigmoid$Freq > 0,]
proportion_celltypes_I<-proportion_celltypes_I[proportion_celltypes_I$Freq > 0,]
proportion_celltypes_NI<-proportion_celltypes_NI[proportion_celltypes_NI$Freq > 0,]

library(dplyr)
#Create ccomparison barplots for  NI rectum and sigmoid
ggplot(proportion_celltypes_NI, aes(Var2,Freq, fill=Var3))+
  geom_bar(position="fill", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))

library(ggpubr)
#create violingplots comparing NI rectum and sigmoid
ggplot (proportion_celltypes_NI, aes (Var2,Freq,fill=Var3)) + geom_violin() + geom_boxplot(width=0.06,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var3, nrow=2) + theme_bw() +scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))   + ylab ("Frequency") + xlab ("Protocol")


