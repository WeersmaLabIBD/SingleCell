library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyr)
library(MAST)
library(readxl)
library(openxlsx)
library(EnhancedVolcano)

data<- readRDS("Data/PSC_processed_oct.rds")
DefaultAssay(data) = "RNA"

# Paired test PSC
Idents(data) = "sample"
data_PSC= subset(data, ident = c("3019","3069","3086","3147","3191","3317","3325"))
data_UC= subset(data, ident = c("3006","3059","3076","3085","3107","3125","3249"))

#IgG selected genes
Idents(data_PSC) = "celltypes"
Idents(data_UC) = "celltypes"

# absorptive enterocytes
PSC_absentero_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_absentero_Pvalues) = "Gene"
PSC_absentero_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_absentero_logFC) = "Gene"
IDs_absentero_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Absorptive_enterocyte"])
for (i in 1:20){
  per = sample(IDs_absentero_PSC, size =399, replace=F) #size=75% of smallest group
  PSCINI_absentero_try <- FindMarkers(data_PSC[,per], subset.ident = "Absorptive_enterocyte", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_absentero_try) <- paste(i, colnames(PSCINI_absentero_try), sep = "_") 
  PSCINI_absentero_try$Gene <- rownames(PSCINI_absentero_try)
  PSC_absentero_Pvalues <- merge(PSC_absentero_Pvalues, PSCINI_absentero_try[,c(5,6)], by = "Gene")
  PSC_absentero_logFC <- merge(PSC_absentero_logFC, PSCINI_absentero_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_absentero_Pvalues, file="Results/DE_2/PSC_absentero_Pvalues.csv")
write.csv(PSC_absentero_logFC, file="Results/DE_2/PSC_absentero_logFC.csv")

UC_absentero_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_absentero_Pvalues) = "Gene"
UC_absentero_logFC <-data.frame(rownames(data_UC))
colnames(UC_absentero_logFC) = "Gene"
IDs_absentero_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Absorptive_enterocyte"])
for (i in 1:20){
  per = sample(IDs_absentero_UC, size =399, replace=F) #size=75% of smallest group
  UCINI_absentero_try <- FindMarkers(data_UC[,per], subset.ident = "Absorptive_enterocyte", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_absentero_try) <- paste(i, colnames(UCINI_absentero_try), sep = "_") 
  UCINI_absentero_try$Gene <- rownames(UCINI_absentero_try)
  UC_absentero_Pvalues <- merge(UC_absentero_Pvalues, UCINI_absentero_try[,c(5,6)], by = "Gene")
  UC_absentero_logFC <- merge(UC_absentero_logFC, UCINI_absentero_try[,c(2,6)], by = "Gene")
}
write.csv(UC_absentero_Pvalues, file="Results/DE_2/UC_absentero_Pvalues.csv")
write.csv(UC_absentero_logFC, file="Results/DE_2/UC_absentero_logFC.csv")

# immature enterocytes
PSC_immentero_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_immentero_Pvalues) = "Gene"
PSC_immentero_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_immentero_logFC) = "Gene"
IDs_immentero_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Immature_enterocyte"])
for (i in 1:20){
  per = sample(IDs_immentero_PSC, size =645, replace=F) #size=75% of smallest group
  PSCINI_immentero_try <- FindMarkers(data_PSC[,per], subset.ident = "Immature_enterocyte", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_immentero_try) <- paste(i, colnames(PSCINI_immentero_try), sep = "_") 
  PSCINI_immentero_try$Gene <- rownames(PSCINI_immentero_try)
  PSC_immentero_Pvalues <- merge(PSC_immentero_Pvalues, PSCINI_immentero_try[,c(5,6)], by = "Gene")
  PSC_immentero_logFC <- merge(PSC_immentero_logFC, PSCINI_immentero_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_immentero_Pvalues, file="Results/DE_2/PSC_immentero_Pvalues.csv")
write.csv(PSC_immentero_logFC, file="Results/DE_2/PSC_immentero_logFC.csv")

UC_immentero_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_immentero_Pvalues) = "Gene"
UC_immentero_logFC <-data.frame(rownames(data_UC))
colnames(UC_immentero_logFC) = "Gene"
IDs_immentero_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Immature_enterocyte"])
for (i in 1:20){
  per = sample(IDs_immentero_UC, size =645, replace=F) #size=75% of smallest group
  UCINI_immentero_try <- FindMarkers(data_UC[,per], subset.ident = "Immature_enterocyte", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_immentero_try) <- paste(i, colnames(UCINI_immentero_try), sep = "_") 
  UCINI_immentero_try$Gene <- rownames(UCINI_immentero_try)
  UC_immentero_Pvalues <- merge(UC_immentero_Pvalues, UCINI_immentero_try[,c(5,6)], by = "Gene")
  UC_immentero_logFC <- merge(UC_immentero_logFC, UCINI_immentero_try[,c(2,6)], by = "Gene")
}
write.csv(UC_immentero_Pvalues, file="Results/DE_2/UC_immentero_Pvalues.csv")
write.csv(UC_immentero_logFC, file="Results/DE_2/UC_immentero_logFC.csv")

# BEST4 enterocytes 
PSC_BEST4_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_BEST4_Pvalues) = "Gene"
PSC_BEST4_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_BEST4_logFC) = "Gene"
IDs_BEST4_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "BEST4_enterocyte"])
for (i in 1:20){
  per = sample(IDs_BEST4_PSC, size =296, replace=F) #size=75% of smallest group
  PSCINI_BEST4_try <- FindMarkers(data_PSC[,per], subset.ident = "BEST4_enterocyte", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_BEST4_try) <- paste(i, colnames(PSCINI_BEST4_try), sep = "_") 
  PSCINI_BEST4_try$Gene <- rownames(PSCINI_BEST4_try)
  PSC_BEST4_Pvalues <- merge(PSC_BEST4_Pvalues, PSCINI_BEST4_try[,c(5,6)], by = "Gene")
  PSC_BEST4_logFC <- merge(PSC_BEST4_logFC, PSCINI_BEST4_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_BEST4_Pvalues, file="Results/DE_2/PSC_BEST4_Pvalues.csv")
write.csv(PSC_BEST4_logFC, file="Results/DE_2/PSC_BEST4_logFC.csv")

UC_BEST4_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_BEST4_Pvalues) = "Gene"
UC_BEST4_logFC <-data.frame(rownames(data_UC))
colnames(UC_BEST4_logFC) = "Gene"
IDs_BEST4_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "BEST4_enterocyte"])
for (i in 1:20){
  per = sample(IDs_BEST4_UC, size =296, replace=F) #size=75% of smallest group
  UCINI_BEST4_try <- FindMarkers(data_UC[,per], subset.ident = "BEST4_enterocyte", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_BEST4_try) <- paste(i, colnames(UCINI_BEST4_try), sep = "_") 
  UCINI_BEST4_try$Gene <- rownames(UCINI_BEST4_try)
  UC_BEST4_Pvalues <- merge(UC_BEST4_Pvalues, UCINI_BEST4_try[,c(5,6)], by = "Gene")
  UC_BEST4_logFC <- merge(UC_BEST4_logFC, UCINI_BEST4_try[,c(2,6)], by = "Gene")
}
write.csv(UC_BEST4_Pvalues, file="Results/DE_2/UC_BEST4_Pvalues.csv")
write.csv(UC_BEST4_logFC, file="Results/DE_2/UC_BEST4_logFC.csv")

# DUOX2 enterocytes 
PSC_DUOX2_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_DUOX2_Pvalues) = "Gene"
PSC_DUOX2_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_DUOX2_logFC) = "Gene"
IDs_DUOX2_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "DUOX2_enterocyte"])
for (i in 1:20){
  per = sample(IDs_DUOX2_PSC, size =320, replace=F) #size=75% of smallest group
  PSCINI_DUOX2_try <- FindMarkers(data_PSC[,per], subset.ident = "DUOX2_enterocyte", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_DUOX2_try) <- paste(i, colnames(PSCINI_DUOX2_try), sep = "_") 
  PSCINI_DUOX2_try$Gene <- rownames(PSCINI_DUOX2_try)
  PSC_DUOX2_Pvalues <- merge(PSC_DUOX2_Pvalues, PSCINI_DUOX2_try[,c(5,6)], by = "Gene")
  PSC_DUOX2_logFC <- merge(PSC_DUOX2_logFC, PSCINI_DUOX2_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_DUOX2_Pvalues, file="Results/DE_2/PSC_DUOX2_Pvalues.csv")
write.csv(PSC_DUOX2_logFC, file="Results/DE_2/PSC_DUOX2_logFC.csv")

UC_DUOX2_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_DUOX2_Pvalues) = "Gene"
UC_DUOX2_logFC <-data.frame(rownames(data_UC))
colnames(UC_DUOX2_logFC) = "Gene"
IDs_DUOX2_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "DUOX2_enterocyte"])
for (i in 1:20){
  per = sample(IDs_DUOX2_UC, size =320, replace=F) #size=75% of smallest group
  UCINI_DUOX2_try <- FindMarkers(data_UC[,per], subset.ident = "DUOX2_enterocyte", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_DUOX2_try) <- paste(i, colnames(UCINI_DUOX2_try), sep = "_") 
  UCINI_DUOX2_try$Gene <- rownames(UCINI_DUOX2_try)
  UC_DUOX2_Pvalues <- merge(UC_DUOX2_Pvalues, UCINI_DUOX2_try[,c(5,6)], by = "Gene")
  UC_DUOX2_logFC <- merge(UC_DUOX2_logFC, UCINI_DUOX2_try[,c(2,6)], by = "Gene")
}
write.csv(UC_DUOX2_Pvalues, file="Results/DE_2/UC_DUOX2_Pvalues.csv")
write.csv(UC_DUOX2_logFC, file="Results/DE_2/UC_DUOX2_logFC.csv")

# Absorptive TA
PSC_absTA_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_absTA_Pvalues) = "Gene"
PSC_absTA_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_absTA_logFC) = "Gene"
IDs_absTA_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Absorptive_TA"])
for (i in 1:20){
  per = sample(IDs_absTA_PSC, size =200, replace=F) #size=75% of smallest group
  PSCINI_absTA_try <- FindMarkers(data_PSC[,per], subset.ident = "Absorptive_TA", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_absTA_try) <- paste(i, colnames(PSCINI_absTA_try), sep = "_") 
  PSCINI_absTA_try$Gene <- rownames(PSCINI_absTA_try)
  PSC_absTA_Pvalues <- merge(PSC_absTA_Pvalues, PSCINI_absTA_try[,c(5,6)], by = "Gene")
  PSC_absTA_logFC <- merge(PSC_absTA_logFC, PSCINI_absTA_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_absTA_Pvalues, file="Results/DE_2/PSC_absTA_Pvalues.csv")
write.csv(PSC_absTA_logFC, file="Results/DE_2/PSC_absTA_logFC.csv")

UC_absTA_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_absTA_Pvalues) = "Gene"
UC_absTA_logFC <-data.frame(rownames(data_UC))
colnames(UC_absTA_logFC) = "Gene"
IDs_absTA_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Absorptive_TA"])
for (i in 1:20){
  per = sample(IDs_absTA_UC, size =200, replace=F) #size=75% of smallest group
  UCINI_absTA_try <- FindMarkers(data_UC[,per], subset.ident = "Absorptive_TA", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_absTA_try) <- paste(i, colnames(UCINI_absTA_try), sep = "_") 
  UCINI_absTA_try$Gene <- rownames(UCINI_absTA_try)
  UC_absTA_Pvalues <- merge(UC_absTA_Pvalues, UCINI_absTA_try[,c(5,6)], by = "Gene")
  UC_absTA_logFC <- merge(UC_absTA_logFC, UCINI_absTA_try[,c(2,6)], by = "Gene")
}
write.csv(UC_absTA_Pvalues, file="Results/DE_2/UC_absTA_Pvalues.csv")
write.csv(UC_absTA_logFC, file="Results/DE_2/UC_absTA_logFC.csv")

# PLCG2 TA
PSC_PLCG2_TA_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_PLCG2_TA_Pvalues) = "Gene"
PSC_PLCG2_TA_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_PLCG2_TA_logFC) = "Gene"
IDs_PLCG2_TA_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "PLCG2_TA"])
for (i in 1:20){
  per = sample(IDs_PLCG2_TA_PSC, size =248, replace=F) #size=75% of smallest group
  PSCINI_PLCG2_TA_try <- FindMarkers(data_PSC[,per], subset.ident = "PLCG2_TA", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_PLCG2_TA_try) <- paste(i, colnames(PSCINI_PLCG2_TA_try), sep = "_") 
  PSCINI_PLCG2_TA_try$Gene <- rownames(PSCINI_PLCG2_TA_try)
  PSC_PLCG2_TA_Pvalues <- merge(PSC_PLCG2_TA_Pvalues, PSCINI_PLCG2_TA_try[,c(5,6)], by = "Gene")
  PSC_PLCG2_TA_logFC <- merge(PSC_PLCG2_TA_logFC, PSCINI_PLCG2_TA_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_PLCG2_TA_Pvalues, file="Results/DE_2/PSC_PLCG2_TA_Pvalues.csv")
write.csv(PSC_PLCG2_TA_logFC, file="Results/DE_2/PSC_PLCG2_TA_logFC.csv")

UC_PLCG2_TA_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_PLCG2_TA_Pvalues) = "Gene"
UC_PLCG2_TA_logFC <-data.frame(rownames(data_UC))
colnames(UC_PLCG2_TA_logFC) = "Gene"
IDs_PLCG2_TA_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "PLCG2_TA"])
for (i in 1:20){
  per = sample(IDs_PLCG2_TA_UC, size =248, replace=F) #size=75% of smallest group
  UCINI_PLCG2_TA_try <- FindMarkers(data_UC[,per], subset.ident = "PLCG2_TA", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_PLCG2_TA_try) <- paste(i, colnames(UCINI_PLCG2_TA_try), sep = "_") 
  UCINI_PLCG2_TA_try$Gene <- rownames(UCINI_PLCG2_TA_try)
  UC_PLCG2_TA_Pvalues <- merge(UC_PLCG2_TA_Pvalues, UCINI_PLCG2_TA_try[,c(5,6)], by = "Gene")
  UC_PLCG2_TA_logFC <- merge(UC_PLCG2_TA_logFC, UCINI_PLCG2_TA_try[,c(2,6)], by = "Gene")
}
write.csv(UC_PLCG2_TA_Pvalues, file="Results/DE_2/UC_PLCG2_TA_Pvalues.csv")
write.csv(UC_PLCG2_TA_logFC, file="Results/DE_2/UC_PLCG2_TA_logFC.csv")

# Ribo TA
PSC_Ribo_TA_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_Ribo_TA_Pvalues) = "Gene"
PSC_Ribo_TA_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_Ribo_TA_logFC) = "Gene"
IDs_Ribo_TA_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Ribo_TA"])
for (i in 1:20){
  per = sample(IDs_Ribo_TA_PSC, size =578, replace=F) #size=75% of smallest group
  PSCINI_Ribo_TA_try <- FindMarkers(data_PSC[,per], subset.ident = "Ribo_TA", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_Ribo_TA_try) <- paste(i, colnames(PSCINI_Ribo_TA_try), sep = "_") 
  PSCINI_Ribo_TA_try$Gene <- rownames(PSCINI_Ribo_TA_try)
  PSC_Ribo_TA_Pvalues <- merge(PSC_Ribo_TA_Pvalues, PSCINI_Ribo_TA_try[,c(5,6)], by = "Gene")
  PSC_Ribo_TA_logFC <- merge(PSC_Ribo_TA_logFC, PSCINI_Ribo_TA_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_Ribo_TA_Pvalues, file="Results/DE_2/PSC_Ribo_TA_Pvalues.csv")
write.csv(PSC_Ribo_TA_logFC, file="Results/DE_2/PSC_Ribo_TA_logFC.csv")

UC_Ribo_TA_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_Ribo_TA_Pvalues) = "Gene"
UC_Ribo_TA_logFC <-data.frame(rownames(data_UC))
colnames(UC_Ribo_TA_logFC) = "Gene"
IDs_Ribo_TA_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Ribo_TA"])
for (i in 1:20){
  per = sample(IDs_Ribo_TA_UC, size =578, replace=F) #size=75% of smallest group
  UCINI_Ribo_TA_try <- FindMarkers(data_UC[,per], subset.ident = "Ribo_TA", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_Ribo_TA_try) <- paste(i, colnames(UCINI_Ribo_TA_try), sep = "_") 
  UCINI_Ribo_TA_try$Gene <- rownames(UCINI_Ribo_TA_try)
  UC_Ribo_TA_Pvalues <- merge(UC_Ribo_TA_Pvalues, UCINI_Ribo_TA_try[,c(5,6)], by = "Gene")
  UC_Ribo_TA_logFC <- merge(UC_Ribo_TA_logFC, UCINI_Ribo_TA_try[,c(2,6)], by = "Gene")
}
write.csv(UC_Ribo_TA_Pvalues, file="Results/DE_2/UC_Ribo_TA_Pvalues.csv")
write.csv(UC_Ribo_TA_logFC, file="Results/DE_2/UC_Ribo_TA_logFC.csv")

# REG TA
PSC_REG_TA_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_REG_TA_Pvalues) = "Gene"
PSC_REG_TA_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_REG_TA_logFC) = "Gene"
IDs_REG_TA_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "REG_TA"])
for (i in 1:20){
  per = sample(IDs_REG_TA_PSC, size =224, replace=F) #size=75% of smallest group
  PSCINI_REG_TA_try <- FindMarkers(data_PSC[,per], subset.ident = "REG_TA", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_REG_TA_try) <- paste(i, colnames(PSCINI_REG_TA_try), sep = "_") 
  PSCINI_REG_TA_try$Gene <- rownames(PSCINI_REG_TA_try)
  PSC_REG_TA_Pvalues <- merge(PSC_REG_TA_Pvalues, PSCINI_REG_TA_try[,c(5,6)], by = "Gene")
  PSC_REG_TA_logFC <- merge(PSC_REG_TA_logFC, PSCINI_REG_TA_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_REG_TA_Pvalues, file="Results/DE_2/PSC_REG_TA_Pvalues.csv")
write.csv(PSC_REG_TA_logFC, file="Results/DE_2/PSC_REG_TA_logFC.csv")

UC_REG_TA_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_REG_TA_Pvalues) = "Gene"
UC_REG_TA_logFC <-data.frame(rownames(data_UC))
colnames(UC_REG_TA_logFC) = "Gene"
IDs_REG_TA_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "REG_TA"])
for (i in 1:20){
  per = sample(IDs_REG_TA_UC, size =224, replace=F) #size=75% of smallest group
  UCINI_REG_TA_try <- FindMarkers(data_UC[,per], subset.ident = "REG_TA", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_REG_TA_try) <- paste(i, colnames(UCINI_REG_TA_try), sep = "_") 
  UCINI_REG_TA_try$Gene <- rownames(UCINI_REG_TA_try)
  UC_REG_TA_Pvalues <- merge(UC_REG_TA_Pvalues, UCINI_REG_TA_try[,c(5,6)], by = "Gene")
  UC_REG_TA_logFC <- merge(UC_REG_TA_logFC, UCINI_REG_TA_try[,c(2,6)], by = "Gene")
}
write.csv(UC_REG_TA_Pvalues, file="Results/DE_2/UC_REG_TA_Pvalues.csv")
write.csv(UC_REG_TA_logFC, file="Results/DE_2/UC_REG_TA_logFC.csv")

# cycling TA
PSC_cyclingTA_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_cyclingTA_Pvalues) = "Gene"
PSC_cyclingTA_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_cyclingTA_logFC) = "Gene"
IDs_cyclingTA_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Cycling_TA"])
for (i in 1:20){
  per = sample(IDs_cyclingTA_PSC, size =636, replace=F) #size=75% of smallest group
  PSCINI_cyclingTA_try <- FindMarkers(data_PSC[,per], subset.ident = "Cycling_TA", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_cyclingTA_try) <- paste(i, colnames(PSCINI_cyclingTA_try), sep = "_") 
  PSCINI_cyclingTA_try$Gene <- rownames(PSCINI_cyclingTA_try)
  PSC_cyclingTA_Pvalues <- merge(PSC_cyclingTA_Pvalues, PSCINI_cyclingTA_try[,c(5,6)], by = "Gene")
  PSC_cyclingTA_logFC <- merge(PSC_cyclingTA_logFC, PSCINI_cyclingTA_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_cyclingTA_Pvalues, file="Results/DE_2/PSC_cyclingTA_Pvalues.csv")
write.csv(PSC_cyclingTA_logFC, file="Results/DE_2/PSC_cyclingTA_logFC.csv")

UC_cyclingTA_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_cyclingTA_Pvalues) = "Gene"
UC_cyclingTA_logFC <-data.frame(rownames(data_UC))
colnames(UC_cyclingTA_logFC) = "Gene"
IDs_cyclingTA_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Cycling_TA"])
for (i in 1:20){
  per = sample(IDs_cyclingTA_UC, size =636, replace=F) #size=75% of smallest group
  UCINI_cyclingTA_try <- FindMarkers(data_UC[,per], subset.ident = "Cycling_TA", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_cyclingTA_try) <- paste(i, colnames(UCINI_cyclingTA_try), sep = "_") 
  UCINI_cyclingTA_try$Gene <- rownames(UCINI_cyclingTA_try)
  UC_cyclingTA_Pvalues <- merge(UC_cyclingTA_Pvalues, UCINI_cyclingTA_try[,c(5,6)], by = "Gene")
  UC_cyclingTA_logFC <- merge(UC_cyclingTA_logFC, UCINI_cyclingTA_try[,c(2,6)], by = "Gene")
}
write.csv(UC_cyclingTA_Pvalues, file="Results/DE_2/UC_cyclingTA_Pvalues.csv")
write.csv(UC_cyclingTA_logFC, file="Results/DE_2/UC_cyclingTA_logFC.csv")

# stem cell
PSC_stem_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_stem_Pvalues) = "Gene"
PSC_stem_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_stem_logFC) = "Gene"
IDs_stem_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "stem"])
for (i in 1:20){
  per = sample(IDs_stem_PSC, size =725, replace=F) #size=75% of smallest group
  PSCINI_stem_try <- FindMarkers(data_PSC[,per], subset.ident = "stem", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_stem_try) <- paste(i, colnames(PSCINI_stem_try), sep = "_") 
  PSCINI_stem_try$Gene <- rownames(PSCINI_stem_try)
  PSC_stem_Pvalues <- merge(PSC_stem_Pvalues, PSCINI_stem_try[,c(5,6)], by = "Gene")
  PSC_stem_logFC <- merge(PSC_stem_logFC, PSCINI_stem_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_stem_Pvalues, file="Results/DE_2/PSC_stem_Pvalues.csv")
write.csv(PSC_stem_logFC, file="Results/DE_2/PSC_stem_logFC.csv")

UC_stem_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_stem_Pvalues) = "Gene"
UC_stem_logFC <-data.frame(rownames(data_UC))
colnames(UC_stem_logFC) = "Gene"
IDs_stem_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "stem"])
for (i in 1:20){
  per = sample(IDs_stem_UC, size =725, replace=F) #size=75% of smallest group
  UCINI_stem_try <- FindMarkers(data_UC[,per], subset.ident = "stem", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_stem_try) <- paste(i, colnames(UCINI_stem_try), sep = "_") 
  UCINI_stem_try$Gene <- rownames(UCINI_stem_try)
  UC_stem_Pvalues <- merge(UC_stem_Pvalues, UCINI_stem_try[,c(5,6)], by = "Gene")
  UC_stem_logFC <- merge(UC_stem_logFC, UCINI_stem_try[,c(2,6)], by = "Gene")
}
write.csv(UC_stem_Pvalues, file="Results/DE_2/UC_stem_Pvalues.csv")
write.csv(UC_stem_logFC, file="Results/DE_2/UC_stem_logFC.csv")

# enteroendocrine
PSC_enteroendocrine_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_enteroendocrine_Pvalues) = "Gene"
PSC_enteroendocrine_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_enteroendocrine_logFC) = "Gene"
IDs_enteroendocrine_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "enteroendocrine"])
for (i in 1:20){
  per = sample(IDs_enteroendocrine_PSC, size =103, replace=F) #size=75% of smallest group
  PSCINI_enteroendocrine_try <- FindMarkers(data_PSC[,per], subset.ident = "enteroendocrine", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_enteroendocrine_try) <- paste(i, colnames(PSCINI_enteroendocrine_try), sep = "_") 
  PSCINI_enteroendocrine_try$Gene <- rownames(PSCINI_enteroendocrine_try)
  PSC_enteroendocrine_Pvalues <- merge(PSC_enteroendocrine_Pvalues, PSCINI_enteroendocrine_try[,c(5,6)], by = "Gene")
  PSC_enteroendocrine_logFC <- merge(PSC_enteroendocrine_logFC, PSCINI_enteroendocrine_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_enteroendocrine_Pvalues, file="Results/DE_2/PSC_enteroendocrine_Pvalues.csv")
write.csv(PSC_enteroendocrine_logFC, file="Results/DE_2/PSC_enteroendocrine_logFC.csv")

UC_enteroendocrine_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_enteroendocrine_Pvalues) = "Gene"
UC_enteroendocrine_logFC <-data.frame(rownames(data_UC))
colnames(UC_enteroendocrine_logFC) = "Gene"
IDs_enteroendocrine_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "enteroendocrine"])
for (i in 1:20){
  per = sample(IDs_enteroendocrine_UC, size =103, replace=F) #size=75% of smallest group
  UCINI_enteroendocrine_try <- FindMarkers(data_UC[,per], subset.ident = "enteroendocrine", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_enteroendocrine_try) <- paste(i, colnames(UCINI_enteroendocrine_try), sep = "_") 
  UCINI_enteroendocrine_try$Gene <- rownames(UCINI_enteroendocrine_try)
  UC_enteroendocrine_Pvalues <- merge(UC_enteroendocrine_Pvalues, UCINI_enteroendocrine_try[,c(5,6)], by = "Gene")
  UC_enteroendocrine_logFC <- merge(UC_enteroendocrine_logFC, UCINI_enteroendocrine_try[,c(2,6)], by = "Gene")
}
write.csv(UC_enteroendocrine_Pvalues, file="Results/DE_2/UC_enteroendocrine_Pvalues.csv")
write.csv(UC_enteroendocrine_logFC, file="Results/DE_2/UC_enteroendocrine_logFC.csv")

# Immature goblet
PSC_Immature_goblet_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_Immature_goblet_Pvalues) = "Gene"
PSC_Immature_goblet_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_Immature_goblet_logFC) = "Gene"
IDs_Immature_goblet_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Immature_goblet"])
for (i in 1:20){
  per = sample(IDs_Immature_goblet_PSC, size =219, replace=F) #size=75% of smallest group
  PSCINI_Immature_goblet_try <- FindMarkers(data_PSC[,per], subset.ident = "Immature_goblet", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_Immature_goblet_try) <- paste(i, colnames(PSCINI_Immature_goblet_try), sep = "_") 
  PSCINI_Immature_goblet_try$Gene <- rownames(PSCINI_Immature_goblet_try)
  PSC_Immature_goblet_Pvalues <- merge(PSC_Immature_goblet_Pvalues, PSCINI_Immature_goblet_try[,c(5,6)], by = "Gene")
  PSC_Immature_goblet_logFC <- merge(PSC_Immature_goblet_logFC, PSCINI_Immature_goblet_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_Immature_goblet_Pvalues, file="Results/DE_2/PSC_Immature_goblet_Pvalues.csv")
write.csv(PSC_Immature_goblet_logFC, file="Results/DE_2/PSC_Immature_goblet_logFC.csv")

UC_Immature_goblet_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_Immature_goblet_Pvalues) = "Gene"
UC_Immature_goblet_logFC <-data.frame(rownames(data_UC))
colnames(UC_Immature_goblet_logFC) = "Gene"
IDs_Immature_goblet_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Immature_goblet"])
for (i in 1:20){
  per = sample(IDs_Immature_goblet_UC, size =219, replace=F) #size=75% of smallest group
  UCINI_Immature_goblet_try <- FindMarkers(data_UC[,per], subset.ident = "Immature_goblet", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_Immature_goblet_try) <- paste(i, colnames(UCINI_Immature_goblet_try), sep = "_") 
  UCINI_Immature_goblet_try$Gene <- rownames(UCINI_Immature_goblet_try)
  UC_Immature_goblet_Pvalues <- merge(UC_Immature_goblet_Pvalues, UCINI_Immature_goblet_try[,c(5,6)], by = "Gene")
  UC_Immature_goblet_logFC <- merge(UC_Immature_goblet_logFC, UCINI_Immature_goblet_try[,c(2,6)], by = "Gene")
}
write.csv(UC_Immature_goblet_Pvalues, file="Results/DE_2/UC_Immature_goblet_Pvalues.csv")
write.csv(UC_Immature_goblet_logFC, file="Results/DE_2/UC_Immature_goblet_logFC.csv")

# goblet
PSC_Goblet_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_Goblet_Pvalues) = "Gene"
PSC_Goblet_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_Goblet_logFC) = "Gene"
IDs_Goblet_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Goblet"])
for (i in 1:20){
  per = sample(IDs_Goblet_PSC, size =41, replace=F) #size=75% of smallest group
  PSCINI_Goblet_try <- FindMarkers(data_PSC[,per], subset.ident = "Goblet", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_Goblet_try) <- paste(i, colnames(PSCINI_Goblet_try), sep = "_") 
  PSCINI_Goblet_try$Gene <- rownames(PSCINI_Goblet_try)
  PSC_Goblet_Pvalues <- merge(PSC_Goblet_Pvalues, PSCINI_Goblet_try[,c(5,6)], by = "Gene")
  PSC_Goblet_logFC <- merge(PSC_Goblet_logFC, PSCINI_Goblet_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_Goblet_Pvalues, file="Results/DE_2/PSC_Goblet_Pvalues.csv")
write.csv(PSC_Goblet_logFC, file="Results/DE_2/PSC_Goblet_logFC.csv")

UC_Goblet_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_Goblet_Pvalues) = "Gene"
UC_Goblet_logFC <-data.frame(rownames(data_UC))
colnames(UC_Goblet_logFC) = "Gene"
IDs_Goblet_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Goblet"])
for (i in 1:20){
  per = sample(IDs_Goblet_UC, size =41, replace=F) #size=75% of smallest group
  UCINI_Goblet_try <- FindMarkers(data_UC[,per], subset.ident = "Goblet", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_Goblet_try) <- paste(i, colnames(UCINI_Goblet_try), sep = "_") 
  UCINI_Goblet_try$Gene <- rownames(UCINI_Goblet_try)
  UC_Goblet_Pvalues <- merge(UC_Goblet_Pvalues, UCINI_Goblet_try[,c(5,6)], by = "Gene")
  UC_Goblet_logFC <- merge(UC_Goblet_logFC, UCINI_Goblet_try[,c(2,6)], by = "Gene")
}
write.csv(UC_Goblet_Pvalues, file="Results/DE_2/UC_Goblet_Pvalues.csv")
write.csv(UC_Goblet_logFC, file="Results/DE_2/UC_Goblet_logFC.csv")

# Tuft
PSC_Tuft_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_Tuft_Pvalues) = "Gene"
PSC_Tuft_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_Tuft_logFC) = "Gene"
IDs_Tuft_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Tuft"])
for (i in 1:20){
  per = sample(IDs_Tuft_PSC, size =179, replace=F) #size=75% of smallest group
  PSCINI_Tuft_try <- FindMarkers(data_PSC[,per], subset.ident = "Tuft", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_Tuft_try) <- paste(i, colnames(PSCINI_Tuft_try), sep = "_") 
  PSCINI_Tuft_try$Gene <- rownames(PSCINI_Tuft_try)
  PSC_Tuft_Pvalues <- merge(PSC_Tuft_Pvalues, PSCINI_Tuft_try[,c(5,6)], by = "Gene")
  PSC_Tuft_logFC <- merge(PSC_Tuft_logFC, PSCINI_Tuft_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_Tuft_Pvalues, file="Results/DE_2/PSC_Tuft_Pvalues.csv")
write.csv(PSC_Tuft_logFC, file="Results/DE_2/PSC_Tuft_logFC.csv")

UC_Tuft_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_Tuft_Pvalues) = "Gene"
UC_Tuft_logFC <-data.frame(rownames(data_UC))
colnames(UC_Tuft_logFC) = "Gene"
IDs_Tuft_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Tuft"])
for (i in 1:20){
  per = sample(IDs_Tuft_UC, size =179, replace=F) #size=75% of smallest group
  UCINI_Tuft_try <- FindMarkers(data_UC[,per], subset.ident = "Tuft", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_Tuft_try) <- paste(i, colnames(UCINI_Tuft_try), sep = "_") 
  UCINI_Tuft_try$Gene <- rownames(UCINI_Tuft_try)
  UC_Tuft_Pvalues <- merge(UC_Tuft_Pvalues, UCINI_Tuft_try[,c(5,6)], by = "Gene")
  UC_Tuft_logFC <- merge(UC_Tuft_logFC, UCINI_Tuft_try[,c(2,6)], by = "Gene")
}
write.csv(UC_Tuft_Pvalues, file="Results/DE_2/UC_Tuft_Pvalues.csv")
write.csv(UC_Tuft_logFC, file="Results/DE_2/UC_Tuft_logFC.csv")

# Mast cell
PSC_MAST_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_MAST_Pvalues) = "Gene"
PSC_MAST_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_MAST_logFC) = "Gene"
IDs_MAST_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "MAST"])
for (i in 1:20){
  per = sample(IDs_MAST_PSC, size =53, replace=F) #size=75% of smallest group
  PSCINI_MAST_try <- FindMarkers(data_PSC[,per], subset.ident = "MAST", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_MAST_try) <- paste(i, colnames(PSCINI_MAST_try), sep = "_") 
  PSCINI_MAST_try$Gene <- rownames(PSCINI_MAST_try)
  PSC_MAST_Pvalues <- merge(PSC_MAST_Pvalues, PSCINI_MAST_try[,c(5,6)], by = "Gene")
  PSC_MAST_logFC <- merge(PSC_MAST_logFC, PSCINI_MAST_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_MAST_Pvalues, file="Results/DE_2/PSC_MAST_Pvalues.csv")
write.csv(PSC_MAST_logFC, file="Results/DE_2/PSC_MAST_logFC.csv")

UC_MAST_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_MAST_Pvalues) = "Gene"
UC_MAST_logFC <-data.frame(rownames(data_UC))
colnames(UC_MAST_logFC) = "Gene"
IDs_MAST_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "MAST"])
for (i in 1:20){
  per = sample(IDs_MAST_UC, size =53, replace=F) #size=75% of smallest group
  UCINI_MAST_try <- FindMarkers(data_UC[,per], subset.ident = "MAST", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_MAST_try) <- paste(i, colnames(UCINI_MAST_try), sep = "_") 
  UCINI_MAST_try$Gene <- rownames(UCINI_MAST_try)
  UC_MAST_Pvalues <- merge(UC_MAST_Pvalues, UCINI_MAST_try[,c(5,6)], by = "Gene")
  UC_MAST_logFC <- merge(UC_MAST_logFC, UCINI_MAST_try[,c(2,6)], by = "Gene")
}
write.csv(UC_MAST_Pvalues, file="Results/DE_2/UC_MAST_Pvalues.csv")
write.csv(UC_MAST_logFC, file="Results/DE_2/UC_MAST_logFC.csv")

# glia
PSC_glia_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_glia_Pvalues) = "Gene"
PSC_glia_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_glia_logFC) = "Gene"
IDs_glia_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "glia"])
for (i in 1:20){
  per = sample(IDs_glia_PSC, size =49, replace=F) #size=75% of smallest group
  PSCINI_glia_try <- FindMarkers(data_PSC[,per], subset.ident = "glia", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_glia_try) <- paste(i, colnames(PSCINI_glia_try), sep = "_") 
  PSCINI_glia_try$Gene <- rownames(PSCINI_glia_try)
  PSC_glia_Pvalues <- merge(PSC_glia_Pvalues, PSCINI_glia_try[,c(5,6)], by = "Gene")
  PSC_glia_logFC <- merge(PSC_glia_logFC, PSCINI_glia_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_glia_Pvalues, file="Results/DE_2/PSC_glia_Pvalues.csv")
write.csv(PSC_glia_logFC, file="Results/DE_2/PSC_glia_logFC.csv")

UC_glia_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_glia_Pvalues) = "Gene"
UC_glia_logFC <-data.frame(rownames(data_UC))
colnames(UC_glia_logFC) = "Gene"
IDs_glia_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "glia"])
for (i in 1:20){
  per = sample(IDs_glia_UC, size =49, replace=F) #size=75% of smallest group
  UCINI_glia_try <- FindMarkers(data_UC[,per], subset.ident = "glia", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_glia_try) <- paste(i, colnames(UCINI_glia_try), sep = "_") 
  UCINI_glia_try$Gene <- rownames(UCINI_glia_try)
  UC_glia_Pvalues <- merge(UC_glia_Pvalues, UCINI_glia_try[,c(5,6)], by = "Gene")
  UC_glia_logFC <- merge(UC_glia_logFC, UCINI_glia_try[,c(2,6)], by = "Gene")
}
write.csv(UC_glia_Pvalues, file="Results/DE_2/UC_glia_Pvalues.csv")
write.csv(UC_glia_logFC, file="Results/DE_2/UC_glia_logFC.csv")

# WNT2B+
PSC_WNT2B_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_WNT2B_Pvalues) = "Gene"
PSC_WNT2B_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_WNT2B_logFC) = "Gene"
IDs_WNT2B_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "WNT2B+"])
for (i in 1:20){
  per = sample(IDs_WNT2B_PSC, size =191, replace=F) #size=75% of smallest group
  PSCINI_WNT2B_try <- FindMarkers(data_PSC[,per], subset.ident = "WNT2B+", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_WNT2B_try) <- paste(i, colnames(PSCINI_WNT2B_try), sep = "_") 
  PSCINI_WNT2B_try$Gene <- rownames(PSCINI_WNT2B_try)
  PSC_WNT2B_Pvalues <- merge(PSC_WNT2B_Pvalues, PSCINI_WNT2B_try[,c(5,6)], by = "Gene")
  PSC_WNT2B_logFC <- merge(PSC_WNT2B_logFC, PSCINI_WNT2B_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_WNT2B_Pvalues, file="Results/DE_2/PSC_WNT2B_Pvalues.csv")
write.csv(PSC_WNT2B_logFC, file="Results/DE_2/PSC_WNT2B_logFC.csv")

UC_WNT2B_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_WNT2B_Pvalues) = "Gene"
UC_WNT2B_logFC <-data.frame(rownames(data_UC))
colnames(UC_WNT2B_logFC) = "Gene"
IDs_WNT2B_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "WNT2B+"])
for (i in 1:20){
  per = sample(IDs_WNT2B_UC, size =191, replace=F) #size=75% of smallest group
  UCINI_WNT2B_try <- FindMarkers(data_UC[,per], subset.ident = "WNT2B+", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_WNT2B_try) <- paste(i, colnames(UCINI_WNT2B_try), sep = "_") 
  UCINI_WNT2B_try$Gene <- rownames(UCINI_WNT2B_try)
  UC_WNT2B_Pvalues <- merge(UC_WNT2B_Pvalues, UCINI_WNT2B_try[,c(5,6)], by = "Gene")
  UC_WNT2B_logFC <- merge(UC_WNT2B_logFC, UCINI_WNT2B_try[,c(2,6)], by = "Gene")
}
write.csv(UC_WNT2B_Pvalues, file="Results/DE_2/UC_WNT2B_Pvalues.csv")
write.csv(UC_WNT2B_logFC, file="Results/DE_2/UC_WNT2B_logFC.csv")

# WNT5B+
PSC_WNT5B_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_WNT5B_Pvalues) = "Gene"
PSC_WNT5B_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_WNT5B_logFC) = "Gene"
IDs_WNT5B_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "WNT5B+"])
for (i in 1:20){
  per = sample(IDs_WNT5B_PSC, size =241, replace=F) #size=75% of smallest group
  PSCINI_WNT5B_try <- FindMarkers(data_PSC[,per], subset.ident = "WNT5B+", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_WNT5B_try) <- paste(i, colnames(PSCINI_WNT5B_try), sep = "_") 
  PSCINI_WNT5B_try$Gene <- rownames(PSCINI_WNT5B_try)
  PSC_WNT5B_Pvalues <- merge(PSC_WNT5B_Pvalues, PSCINI_WNT5B_try[,c(5,6)], by = "Gene")
  PSC_WNT5B_logFC <- merge(PSC_WNT5B_logFC, PSCINI_WNT5B_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_WNT5B_Pvalues, file="Results/DE_2/PSC_WNT5B_Pvalues.csv")
write.csv(PSC_WNT5B_logFC, file="Results/DE_2/PSC_WNT5B_logFC.csv")

UC_WNT5B_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_WNT5B_Pvalues) = "Gene"
UC_WNT5B_logFC <-data.frame(rownames(data_UC))
colnames(UC_WNT5B_logFC) = "Gene"
IDs_WNT5B_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "WNT5B+"])
for (i in 1:20){
  per = sample(IDs_WNT5B_UC, size =241, replace=F) #size=75% of smallest group
  UCINI_WNT5B_try <- FindMarkers(data_UC[,per], subset.ident = "WNT5B+", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_WNT5B_try) <- paste(i, colnames(UCINI_WNT5B_try), sep = "_") 
  UCINI_WNT5B_try$Gene <- rownames(UCINI_WNT5B_try)
  UC_WNT5B_Pvalues <- merge(UC_WNT5B_Pvalues, UCINI_WNT5B_try[,c(5,6)], by = "Gene")
  UC_WNT5B_logFC <- merge(UC_WNT5B_logFC, UCINI_WNT5B_try[,c(2,6)], by = "Gene")
}
write.csv(UC_WNT5B_Pvalues, file="Results/DE_2/UC_WNT5B_Pvalues.csv")
write.csv(UC_WNT5B_logFC, file="Results/DE_2/UC_WNT5B_logFC.csv")

# RSPO3+
PSC_RSPO3_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_RSPO3_Pvalues) = "Gene"
PSC_RSPO3_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_RSPO3_logFC) = "Gene"
IDs_RSPO3_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "RSPO3+"])
for (i in 1:20){
  per = sample(IDs_RSPO3_PSC, size =153, replace=F) #size=75% of smallest group
  PSCINI_RSPO3_try <- FindMarkers(data_PSC[,per], subset.ident = "RSPO3+", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_RSPO3_try) <- paste(i, colnames(PSCINI_RSPO3_try), sep = "_") 
  PSCINI_RSPO3_try$Gene <- rownames(PSCINI_RSPO3_try)
  PSC_RSPO3_Pvalues <- merge(PSC_RSPO3_Pvalues, PSCINI_RSPO3_try[,c(5,6)], by = "Gene")
  PSC_RSPO3_logFC <- merge(PSC_RSPO3_logFC, PSCINI_RSPO3_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_RSPO3_Pvalues, file="Results/DE_2/PSC_RSPO3_Pvalues.csv")
write.csv(PSC_RSPO3_logFC, file="Results/DE_2/PSC_RSPO3_logFC.csv")

UC_RSPO3_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_RSPO3_Pvalues) = "Gene"
UC_RSPO3_logFC <-data.frame(rownames(data_UC))
colnames(UC_RSPO3_logFC) = "Gene"
IDs_RSPO3_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "RSPO3+"])
for (i in 1:20){
  per = sample(IDs_RSPO3_UC, size =153, replace=F) #size=75% of smallest group
  UCINI_RSPO3_try <- FindMarkers(data_UC[,per], subset.ident = "RSPO3+", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_RSPO3_try) <- paste(i, colnames(UCINI_RSPO3_try), sep = "_") 
  UCINI_RSPO3_try$Gene <- rownames(UCINI_RSPO3_try)
  UC_RSPO3_Pvalues <- merge(UC_RSPO3_Pvalues, UCINI_RSPO3_try[,c(5,6)], by = "Gene")
  UC_RSPO3_logFC <- merge(UC_RSPO3_logFC, UCINI_RSPO3_try[,c(2,6)], by = "Gene")
}
write.csv(UC_RSPO3_Pvalues, file="Results/DE_2/UC_RSPO3_Pvalues.csv")
write.csv(UC_RSPO3_logFC, file="Results/DE_2/UC_RSPO3_logFC.csv")

# Inflammatory_fibroblast
PSC_inflamfibro_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_inflamfibro_Pvalues) = "Gene"
PSC_inflamfibro_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_inflamfibro_logFC) = "Gene"
IDs_inflamfibro_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Inflammatory_fibroblast"])
for (i in 1:20){
  per = sample(IDs_inflamfibro_PSC, size =176, replace=F) #size=75% of smallest group
  PSCINI_inflamfibro_try <- FindMarkers(data_PSC[,per], subset.ident = "Inflammatory_fibroblast", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_inflamfibro_try) <- paste(i, colnames(PSCINI_inflamfibro_try), sep = "_") 
  PSCINI_inflamfibro_try$Gene <- rownames(PSCINI_inflamfibro_try)
  PSC_inflamfibro_Pvalues <- merge(PSC_inflamfibro_Pvalues, PSCINI_inflamfibro_try[,c(5,6)], by = "Gene")
  PSC_inflamfibro_logFC <- merge(PSC_inflamfibro_logFC, PSCINI_inflamfibro_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_inflamfibro_Pvalues, file="Results/DE_2/PSC_inflamfibro_Pvalues.csv")
write.csv(PSC_inflamfibro_logFC, file="Results/DE_2/PSC_inflamfibro_logFC.csv")

UC_inflamfibro_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_inflamfibro_Pvalues) = "Gene"
UC_inflamfibro_logFC <-data.frame(rownames(data_UC))
colnames(UC_inflamfibro_logFC) = "Gene"
IDs_inflamfibro_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Inflammatory_fibroblast"])
for (i in 1:20){
  per = sample(IDs_inflamfibro_UC, size =176, replace=F) #size=75% of smallest group
  UCINI_inflamfibro_try <- FindMarkers(data_UC[,per], subset.ident = "Inflammatory_fibroblast", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_inflamfibro_try) <- paste(i, colnames(UCINI_inflamfibro_try), sep = "_") 
  UCINI_inflamfibro_try$Gene <- rownames(UCINI_inflamfibro_try)
  UC_inflamfibro_Pvalues <- merge(UC_inflamfibro_Pvalues, UCINI_inflamfibro_try[,c(5,6)], by = "Gene")
  UC_inflamfibro_logFC <- merge(UC_inflamfibro_logFC, UCINI_inflamfibro_try[,c(2,6)], by = "Gene")
}
write.csv(UC_inflamfibro_Pvalues, file="Results/DE_2/UC_inflamfibro_Pvalues.csv")
write.csv(UC_inflamfibro_logFC, file="Results/DE_2/UC_inflamfibro_logFC.csv")

# myofibroblasts
PSC_myofibro_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_myofibro_Pvalues) = "Gene"
PSC_myofibro_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_myofibro_logFC) = "Gene"
IDs_myofibro_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "myofibroblasts"])
for (i in 1:20){
  per = sample(IDs_myofibro_PSC, size =76, replace=F) #size=75% of smallest group
  PSCINI_myofibro_try <- FindMarkers(data_PSC[,per], subset.ident = "myofibroblasts", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_myofibro_try) <- paste(i, colnames(PSCINI_myofibro_try), sep = "_") 
  PSCINI_myofibro_try$Gene <- rownames(PSCINI_myofibro_try)
  PSC_myofibro_Pvalues <- merge(PSC_myofibro_Pvalues, PSCINI_myofibro_try[,c(5,6)], by = "Gene")
  PSC_myofibro_logFC <- merge(PSC_myofibro_logFC, PSCINI_myofibro_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_myofibro_Pvalues, file="Results/DE_2/PSC_myofibro_Pvalues.csv")
write.csv(PSC_myofibro_logFC, file="Results/DE_2/PSC_myofibro_logFC.csv")

UC_myofibro_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_myofibro_Pvalues) = "Gene"
UC_myofibro_logFC <-data.frame(rownames(data_UC))
colnames(UC_myofibro_logFC) = "Gene"
IDs_myofibro_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "myofibroblasts"])
for (i in 1:20){
  per = sample(IDs_myofibro_UC, size =76, replace=F) #size=75% of smallest group
  UCINI_myofibro_try <- FindMarkers(data_UC[,per], subset.ident = "myofibroblasts", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_myofibro_try) <- paste(i, colnames(UCINI_myofibro_try), sep = "_") 
  UCINI_myofibro_try$Gene <- rownames(UCINI_myofibro_try)
  UC_myofibro_Pvalues <- merge(UC_myofibro_Pvalues, UCINI_myofibro_try[,c(5,6)], by = "Gene")
  UC_myofibro_logFC <- merge(UC_myofibro_logFC, UCINI_myofibro_try[,c(2,6)], by = "Gene")
}
write.csv(UC_myofibro_Pvalues, file="Results/DE_2/UC_myofibro_Pvalues.csv")
write.csv(UC_myofibro_logFC, file="Results/DE_2/UC_myofibro_logFC.csv")

# pericytes
PSC_pericytes_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_pericytes_Pvalues) = "Gene"
PSC_pericytes_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_pericytes_logFC) = "Gene"
IDs_pericytes_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "pericytes"])
for (i in 1:20){
  per = sample(IDs_pericytes_PSC, size =131, replace=F) #size=75% of smallest group
  PSCINI_pericytes_try <- FindMarkers(data_PSC[,per], subset.ident = "pericytes", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_pericytes_try) <- paste(i, colnames(PSCINI_pericytes_try), sep = "_") 
  PSCINI_pericytes_try$Gene <- rownames(PSCINI_pericytes_try)
  PSC_pericytes_Pvalues <- merge(PSC_pericytes_Pvalues, PSCINI_pericytes_try[,c(5,6)], by = "Gene")
  PSC_pericytes_logFC <- merge(PSC_pericytes_logFC, PSCINI_pericytes_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_pericytes_Pvalues, file="Results/DE_2/PSC_pericytes_Pvalues.csv")
write.csv(PSC_pericytes_logFC, file="Results/DE_2/PSC_pericytes_logFC.csv")

UC_pericytes_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_pericytes_Pvalues) = "Gene"
UC_pericytes_logFC <-data.frame(rownames(data_UC))
colnames(UC_pericytes_logFC) = "Gene"
IDs_pericytes_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "pericytes"])
for (i in 1:20){
  per = sample(IDs_pericytes_UC, size =131, replace=F) #size=75% of smallest group
  UCINI_pericytes_try <- FindMarkers(data_UC[,per], subset.ident = "pericytes", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_pericytes_try) <- paste(i, colnames(UCINI_pericytes_try), sep = "_") 
  UCINI_pericytes_try$Gene <- rownames(UCINI_pericytes_try)
  UC_pericytes_Pvalues <- merge(UC_pericytes_Pvalues, UCINI_pericytes_try[,c(5,6)], by = "Gene")
  UC_pericytes_logFC <- merge(UC_pericytes_logFC, UCINI_pericytes_try[,c(2,6)], by = "Gene")
}
write.csv(UC_pericytes_Pvalues, file="Results/DE_2/UC_pericytes_Pvalues.csv")
write.csv(UC_pericytes_logFC, file="Results/DE_2/UC_pericytes_logFC.csv")

# endothelial
PSC_endothelial_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_endothelial_Pvalues) = "Gene"
PSC_endothelial_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_endothelial_logFC) = "Gene"
IDs_endothelial_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "endothelial"])
for (i in 1:20){
  per = sample(IDs_endothelial_PSC, size =241, replace=F) #size=75% of smallest group
  PSCINI_endothelial_try <- FindMarkers(data_PSC[,per], subset.ident = "endothelial", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_endothelial_try) <- paste(i, colnames(PSCINI_endothelial_try), sep = "_") 
  PSCINI_endothelial_try$Gene <- rownames(PSCINI_endothelial_try)
  PSC_endothelial_Pvalues <- merge(PSC_endothelial_Pvalues, PSCINI_endothelial_try[,c(5,6)], by = "Gene")
  PSC_endothelial_logFC <- merge(PSC_endothelial_logFC, PSCINI_endothelial_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_endothelial_Pvalues, file="Results/DE_2/PSC_endothelial_Pvalues.csv")
write.csv(PSC_endothelial_logFC, file="Results/DE_2/PSC_endothelial_logFC.csv")

UC_endothelial_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_endothelial_Pvalues) = "Gene"
UC_endothelial_logFC <-data.frame(rownames(data_UC))
colnames(UC_endothelial_logFC) = "Gene"
IDs_endothelial_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "endothelial"])
for (i in 1:20){
  per = sample(IDs_endothelial_UC, size =241, replace=F) #size=75% of smallest group
  UCINI_endothelial_try <- FindMarkers(data_UC[,per], subset.ident = "endothelial", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_endothelial_try) <- paste(i, colnames(UCINI_endothelial_try), sep = "_") 
  UCINI_endothelial_try$Gene <- rownames(UCINI_endothelial_try)
  UC_endothelial_Pvalues <- merge(UC_endothelial_Pvalues, UCINI_endothelial_try[,c(5,6)], by = "Gene")
  UC_endothelial_logFC <- merge(UC_endothelial_logFC, UCINI_endothelial_try[,c(2,6)], by = "Gene")
}
write.csv(UC_endothelial_Pvalues, file="Results/DE_2/UC_endothelial_Pvalues.csv")
write.csv(UC_endothelial_logFC, file="Results/DE_2/UC_endothelial_logFC.csv")

# IgG plasma
PSC_IgG_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_IgG_Pvalues) = "Gene"
PSC_IgG_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_IgG_logFC) = "Gene"
IDs_IgG_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "IgG_plasma"])
for (i in 1:20){
  per = sample(IDs_IgG_PSC, size =1621, replace=F) #size=75% of smallest group
  PSCINI_IgG_try <- FindMarkers(data_PSC[,per], subset.ident = "IgG_plasma", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_IgG_try) <- paste(i, colnames(PSCINI_IgG_try), sep = "_") 
  PSCINI_IgG_try$Gene <- rownames(PSCINI_IgG_try)
  PSC_IgG_Pvalues <- merge(PSC_IgG_Pvalues, PSCINI_IgG_try[,c(5,6)], by = "Gene")
  PSC_IgG_logFC <- merge(PSC_IgG_logFC, PSCINI_IgG_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_IgG_Pvalues, file="Results/DE_2/PSC_IgG_Pvalues.csv")
write.csv(PSC_IgG_logFC, file="Results/DE_2/PSC_IgG_logFC.csv")

UC_IgG_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_IgG_Pvalues) = "Gene"
UC_IgG_logFC <-data.frame(rownames(data_UC))
colnames(UC_IgG_logFC) = "Gene"
IDs_IgG_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "IgG_plasma"])
for (i in 1:20){
  per = sample(IDs_IgG_UC, size =1621, replace=F) #size=75% of smallest group
  UCINI_IgG_try <- FindMarkers(data_UC[,per], subset.ident = "IgG_plasma", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_IgG_try) <- paste(i, colnames(UCINI_IgG_try), sep = "_") 
  UCINI_IgG_try$Gene <- rownames(UCINI_IgG_try)
  UC_IgG_Pvalues <- merge(UC_IgG_Pvalues, UCINI_IgG_try[,c(5,6)], by = "Gene")
  UC_IgG_logFC <- merge(UC_IgG_logFC, UCINI_IgG_try[,c(2,6)], by = "Gene")
}
write.csv(UC_IgG_Pvalues, file="Results/DE_2/UC_IgG_Pvalues.csv")
write.csv(UC_IgG_logFC, file="Results/DE_2/UC_IgG_logFC.csv")

# IgA plasma
PSC_IgA_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_IgA_Pvalues) = "Gene"
PSC_IgA_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_IgA_logFC) = "Gene"
IDs_IgA_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "IgA_plasma"])
for (i in 1:20){
  per = sample(IDs_IgA_PSC, size =2372, replace=F) #size=75% of smallest group
  PSCINI_IgA_try <- FindMarkers(data_PSC[,per], subset.ident = "IgA_plasma", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_IgA_try) <- paste(i, colnames(PSCINI_IgA_try), sep = "_") 
  PSCINI_IgA_try$Gene <- rownames(PSCINI_IgA_try)
  PSC_IgA_Pvalues <- merge(PSC_IgA_Pvalues, PSCINI_IgA_try[,c(5,6)], by = "Gene")
  PSC_IgA_logFC <- merge(PSC_IgA_logFC, PSCINI_IgA_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_IgA_Pvalues, file="Results/DE_2/PSC_IgA_Pvalues.csv")
write.csv(PSC_IgA_logFC, file="Results/DE_2/PSC_IgA_logFC.csv")

UC_IgA_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_IgA_Pvalues) = "Gene"
UC_IgA_logFC <-data.frame(rownames(data_UC))
colnames(UC_IgA_logFC) = "Gene"
IDs_IgA_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "IgA_plasma"])
for (i in 1:20){
  per = sample(IDs_IgA_UC, size =2372, replace=F) #size=75% of smallest group
  UCINI_IgA_try <- FindMarkers(data_UC[,per], subset.ident = "IgA_plasma", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_IgA_try) <- paste(i, colnames(UCINI_IgA_try), sep = "_") 
  UCINI_IgA_try$Gene <- rownames(UCINI_IgA_try)
  UC_IgA_Pvalues <- merge(UC_IgA_Pvalues, UCINI_IgA_try[,c(5,6)], by = "Gene")
  UC_IgA_logFC <- merge(UC_IgA_logFC, UCINI_IgA_try[,c(2,6)], by = "Gene")
}
write.csv(UC_IgA_Pvalues, file="Results/DE_2/UC_IgA_Pvalues.csv")
write.csv(UC_IgA_logFC, file="Results/DE_2/UC_IgA_logFC.csv")

# IgM plasma
PSC_IgM_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_IgM_Pvalues) = "Gene"
PSC_IgM_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_IgM_logFC) = "Gene"
IDs_IgM_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "IgM_plasma"])
for (i in 1:20){
  per = sample(IDs_IgM_PSC, size =257, replace=F) #size=75% of smallest group
  PSCINI_IgM_try <- FindMarkers(data_PSC[,per], subset.ident = "IgM_plasma", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_IgM_try) <- paste(i, colnames(PSCINI_IgM_try), sep = "_") 
  PSCINI_IgM_try$Gene <- rownames(PSCINI_IgM_try)
  PSC_IgM_Pvalues <- merge(PSC_IgM_Pvalues, PSCINI_IgM_try[,c(5,6)], by = "Gene")
  PSC_IgM_logFC <- merge(PSC_IgM_logFC, PSCINI_IgM_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_IgM_Pvalues, file="Results/DE_2/PSC_IgM_Pvalues.csv")
write.csv(PSC_IgM_logFC, file="Results/DE_2/PSC_IgM_logFC.csv")

UC_IgM_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_IgM_Pvalues) = "Gene"
UC_IgM_logFC <-data.frame(rownames(data_UC))
colnames(UC_IgM_logFC) = "Gene"
IDs_IgM_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "IgM_plasma"])
for (i in 1:20){
  per = sample(IDs_IgM_UC, size =257, replace=F) #size=75% of smallest group
  UCINI_IgM_try <- FindMarkers(data_UC[,per], subset.ident = "IgM_plasma", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_IgM_try) <- paste(i, colnames(UCINI_IgM_try), sep = "_") 
  UCINI_IgM_try$Gene <- rownames(UCINI_IgM_try)
  UC_IgM_Pvalues <- merge(UC_IgM_Pvalues, UCINI_IgM_try[,c(5,6)], by = "Gene")
  UC_IgM_logFC <- merge(UC_IgM_logFC, UCINI_IgM_try[,c(2,6)], by = "Gene")
}
write.csv(UC_IgM_Pvalues, file="Results/DE_2/UC_IgM_Pvalues.csv")
write.csv(UC_IgM_logFC, file="Results/DE_2/UC_IgM_logFC.csv")

# cyclingB plasma
PSC_cyclingB_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_cyclingB_Pvalues) = "Gene"
PSC_cyclingB_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_cyclingB_logFC) = "Gene"
IDs_cyclingB_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Cycling_B"])
for (i in 1:20){
  per = sample(IDs_cyclingB_PSC, size =43, replace=F) #size=75% of smallest group
  PSCINI_cyclingB_try <- FindMarkers(data_PSC[,per], subset.ident = "Cycling_B", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_cyclingB_try) <- paste(i, colnames(PSCINI_cyclingB_try), sep = "_") 
  PSCINI_cyclingB_try$Gene <- rownames(PSCINI_cyclingB_try)
  PSC_cyclingB_Pvalues <- merge(PSC_cyclingB_Pvalues, PSCINI_cyclingB_try[,c(5,6)], by = "Gene")
  PSC_cyclingB_logFC <- merge(PSC_cyclingB_logFC, PSCINI_cyclingB_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_cyclingB_Pvalues, file="Results/DE_2/PSC_cyclingB_Pvalues.csv")
write.csv(PSC_cyclingB_logFC, file="Results/DE_2/PSC_cyclingB_logFC.csv")

UC_cyclingB_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_cyclingB_Pvalues) = "Gene"
UC_cyclingB_logFC <-data.frame(rownames(data_UC))
colnames(UC_cyclingB_logFC) = "Gene"
IDs_cyclingB_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Cycling_B"])
for (i in 1:20){
  per = sample(IDs_cyclingB_UC, size =43, replace=F) #size=75% of smallest group
  UCINI_cyclingB_try <- FindMarkers(data_UC[,per], subset.ident = "Cycling_B", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_cyclingB_try) <- paste(i, colnames(UCINI_cyclingB_try), sep = "_") 
  UCINI_cyclingB_try$Gene <- rownames(UCINI_cyclingB_try)
  UC_cyclingB_Pvalues <- merge(UC_cyclingB_Pvalues, UCINI_cyclingB_try[,c(5,6)], by = "Gene")
  UC_cyclingB_logFC <- merge(UC_cyclingB_logFC, UCINI_cyclingB_try[,c(2,6)], by = "Gene")
}
write.csv(UC_cyclingB_Pvalues, file="Results/DE_2/UC_cyclingB_Pvalues.csv")
write.csv(UC_cyclingB_logFC, file="Results/DE_2/UC_cyclingB_logFC.csv")

# GC 
PSC_GC_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_GC_Pvalues) = "Gene"
PSC_GC_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_GC_logFC) = "Gene"
IDs_GC_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "GC"])
for (i in 1:20){
  per = sample(IDs_GC_PSC, size =71, replace=F) #size=75% of smallest group
  PSCINI_GC_try <- FindMarkers(data_PSC[,per], subset.ident = "GC", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_GC_try) <- paste(i, colnames(PSCINI_GC_try), sep = "_") 
  PSCINI_GC_try$Gene <- rownames(PSCINI_GC_try)
  PSC_GC_Pvalues <- merge(PSC_GC_Pvalues, PSCINI_GC_try[,c(5,6)], by = "Gene")
  PSC_GC_logFC <- merge(PSC_GC_logFC, PSCINI_GC_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_GC_Pvalues, file="Results/DE_2/PSC_GC_Pvalues.csv")
write.csv(PSC_GC_logFC, file="Results/DE_2/PSC_GC_logFC.csv")

UC_GC_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_GC_Pvalues) = "Gene"
UC_GC_logFC <-data.frame(rownames(data_UC))
colnames(UC_GC_logFC) = "Gene"
IDs_GC_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "GC"])
for (i in 1:20){
  per = sample(IDs_GC_UC, size =71, replace=F) #size=75% of smallest group
  UCINI_GC_try <- FindMarkers(data_UC[,per], subset.ident = "GC", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_GC_try) <- paste(i, colnames(UCINI_GC_try), sep = "_") 
  UCINI_GC_try$Gene <- rownames(UCINI_GC_try)
  UC_GC_Pvalues <- merge(UC_GC_Pvalues, UCINI_GC_try[,c(5,6)], by = "Gene")
  UC_GC_logFC <- merge(UC_GC_logFC, UCINI_GC_try[,c(2,6)], by = "Gene")
}
write.csv(UC_GC_Pvalues, file="Results/DE_2/UC_GC_Pvalues.csv")
write.csv(UC_GC_logFC, file="Results/DE_2/UC_GC_logFC.csv")


# Follicular 
PSC_Follicular_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_Follicular_Pvalues) = "Gene"
PSC_Follicular_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_Follicular_logFC) = "Gene"
IDs_Follicular_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Follicular_B"])
for (i in 1:20){
  per = sample(IDs_Follicular_PSC, size =1058, replace=F) #size=75% of smallest group
  PSCINI_Follicular_try <- FindMarkers(data_PSC[,per], subset.ident = "Follicular_B", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_Follicular_try) <- paste(i, colnames(PSCINI_Follicular_try), sep = "_") 
  PSCINI_Follicular_try$Gene <- rownames(PSCINI_Follicular_try)
  PSC_Follicular_Pvalues <- merge(PSC_Follicular_Pvalues, PSCINI_Follicular_try[,c(5,6)], by = "Gene")
  PSC_Follicular_logFC <- merge(PSC_Follicular_logFC, PSCINI_Follicular_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_Follicular_Pvalues, file="Results/DE_2/PSC_Follicular_Pvalues.csv")
write.csv(PSC_Follicular_logFC, file="Results/DE_2/PSC_Follicular_logFC.csv")

UC_Follicular_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_Follicular_Pvalues) = "Gene"
UC_Follicular_logFC <-data.frame(rownames(data_UC))
colnames(UC_Follicular_logFC) = "Gene"
IDs_Follicular_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Follicular_B"])
for (i in 1:20){
  per = sample(IDs_Follicular_UC, size =1058, replace=F) #size=75% of smallest group
  UCINI_Follicular_try <- FindMarkers(data_UC[,per], subset.ident = "Follicular_B", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_Follicular_try) <- paste(i, colnames(UCINI_Follicular_try), sep = "_") 
  UCINI_Follicular_try$Gene <- rownames(UCINI_Follicular_try)
  UC_Follicular_Pvalues <- merge(UC_Follicular_Pvalues, UCINI_Follicular_try[,c(5,6)], by = "Gene")
  UC_Follicular_logFC <- merge(UC_Follicular_logFC, UCINI_Follicular_try[,c(2,6)], by = "Gene")
}
write.csv(UC_Follicular_Pvalues, file="Results/DE_2/UC_Follicular_Pvalues.csv")
write.csv(UC_Follicular_logFC, file="Results/DE_2/UC_Follicular_logFC.csv")

# CD8T 
PSC_CD8T_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_CD8T_Pvalues) = "Gene"
PSC_CD8T_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_CD8T_logFC) = "Gene"
IDs_CD8T_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "CD8T"])
for (i in 1:20){
  per = sample(IDs_CD8T_PSC, size =579, replace=F) #size=75% of smallest group
  PSCINI_CD8T_try <- FindMarkers(data_PSC[,per], subset.ident = "CD8T", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_CD8T_try) <- paste(i, colnames(PSCINI_CD8T_try), sep = "_") 
  PSCINI_CD8T_try$Gene <- rownames(PSCINI_CD8T_try)
  PSC_CD8T_Pvalues <- merge(PSC_CD8T_Pvalues, PSCINI_CD8T_try[,c(5,6)], by = "Gene")
  PSC_CD8T_logFC <- merge(PSC_CD8T_logFC, PSCINI_CD8T_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_CD8T_Pvalues, file="Results/DE_2/PSC_CD8T_Pvalues.csv")
write.csv(PSC_CD8T_logFC, file="Results/DE_2/PSC_CD8T_logFC.csv")

UC_CD8T_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_CD8T_Pvalues) = "Gene"
UC_CD8T_logFC <-data.frame(rownames(data_UC))
colnames(UC_CD8T_logFC) = "Gene"
IDs_CD8T_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "CD8T"])
for (i in 1:20){
  per = sample(IDs_CD8T_UC, size =579, replace=F) #size=75% of smallest group
  UCINI_CD8T_try <- FindMarkers(data_UC[,per], subset.ident = "CD8T", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_CD8T_try) <- paste(i, colnames(UCINI_CD8T_try), sep = "_") 
  UCINI_CD8T_try$Gene <- rownames(UCINI_CD8T_try)
  UC_CD8T_Pvalues <- merge(UC_CD8T_Pvalues, UCINI_CD8T_try[,c(5,6)], by = "Gene")
  UC_CD8T_logFC <- merge(UC_CD8T_logFC, UCINI_CD8T_try[,c(2,6)], by = "Gene")
}
write.csv(UC_CD8T_Pvalues, file="Results/DE_2/UC_CD8T_Pvalues.csv")
write.csv(UC_CD8T_logFC, file="Results/DE_2/UC_CD8T_logFC.csv")

# actT_FOSlow 
PSC_actT_FOSlow_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_actT_FOSlow_Pvalues) = "Gene"
PSC_actT_FOSlow_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_actT_FOSlow_logFC) = "Gene"
IDs_actT_FOSlow_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "T_activated_FOS_low"])
for (i in 1:20){
  per = sample(IDs_actT_FOSlow_PSC, size =413, replace=F) #size=75% of smallest group
  PSCINI_actT_FOSlow_try <- FindMarkers(data_PSC[,per], subset.ident = "T_activated_FOS_low", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_actT_FOSlow_try) <- paste(i, colnames(PSCINI_actT_FOSlow_try), sep = "_") 
  PSCINI_actT_FOSlow_try$Gene <- rownames(PSCINI_actT_FOSlow_try)
  PSC_actT_FOSlow_Pvalues <- merge(PSC_actT_FOSlow_Pvalues, PSCINI_actT_FOSlow_try[,c(5,6)], by = "Gene")
  PSC_actT_FOSlow_logFC <- merge(PSC_actT_FOSlow_logFC, PSCINI_actT_FOSlow_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_actT_FOSlow_Pvalues, file="Results/DE_2/PSC_actT_FOSlow_Pvalues.csv")
write.csv(PSC_actT_FOSlow_logFC, file="Results/DE_2/PSC_actT_FOSlow_logFC.csv")

UC_actT_FOSlow_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_actT_FOSlow_Pvalues) = "Gene"
UC_actT_FOSlow_logFC <-data.frame(rownames(data_UC))
colnames(UC_actT_FOSlow_logFC) = "Gene"
IDs_actT_FOSlow_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "T_activated_FOS_low"])
for (i in 1:20){
  per = sample(IDs_actT_FOSlow_UC, size =413, replace=F) #size=75% of smallest group
  UCINI_actT_FOSlow_try <- FindMarkers(data_UC[,per], subset.ident = "T_activated_FOS_low", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_actT_FOSlow_try) <- paste(i, colnames(UCINI_actT_FOSlow_try), sep = "_") 
  UCINI_actT_FOSlow_try$Gene <- rownames(UCINI_actT_FOSlow_try)
  UC_actT_FOSlow_Pvalues <- merge(UC_actT_FOSlow_Pvalues, UCINI_actT_FOSlow_try[,c(5,6)], by = "Gene")
  UC_actT_FOSlow_logFC <- merge(UC_actT_FOSlow_logFC, UCINI_actT_FOSlow_try[,c(2,6)], by = "Gene")
}
write.csv(UC_actT_FOSlow_Pvalues, file="Results/DE_2/UC_actT_FOSlow_Pvalues.csv")
write.csv(UC_actT_FOSlow_logFC, file="Results/DE_2/UC_actT_FOSlow_logFC.csv")

# actT_FOShi 
PSC_actT_FOShi_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_actT_FOShi_Pvalues) = "Gene"
PSC_actT_FOShi_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_actT_FOShi_logFC) = "Gene"
IDs_actT_FOShi_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "T_activated_FOS_high"])
for (i in 1:20){
  per = sample(IDs_actT_FOShi_PSC, size =347, replace=F) #size=75% of smallest group
  PSCINI_actT_FOShi_try <- FindMarkers(data_PSC[,per], subset.ident = "T_activated_FOS_high", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_actT_FOShi_try) <- paste(i, colnames(PSCINI_actT_FOShi_try), sep = "_") 
  PSCINI_actT_FOShi_try$Gene <- rownames(PSCINI_actT_FOShi_try)
  PSC_actT_FOShi_Pvalues <- merge(PSC_actT_FOShi_Pvalues, PSCINI_actT_FOShi_try[,c(5,6)], by = "Gene")
  PSC_actT_FOShi_logFC <- merge(PSC_actT_FOShi_logFC, PSCINI_actT_FOShi_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_actT_FOShi_Pvalues, file="Results/DE_2/PSC_actT_FOShi_Pvalues.csv")
write.csv(PSC_actT_FOShi_logFC, file="Results/DE_2/PSC_actT_FOShi_logFC.csv")

UC_actT_FOShi_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_actT_FOShi_Pvalues) = "Gene"
UC_actT_FOShi_logFC <-data.frame(rownames(data_UC))
colnames(UC_actT_FOShi_logFC) = "Gene"
IDs_actT_FOShi_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "T_activated_FOS_high"])
for (i in 1:20){
  per = sample(IDs_actT_FOShi_UC, size =347, replace=F) #size=75% of smallest group
  UCINI_actT_FOShi_try <- FindMarkers(data_UC[,per], subset.ident = "T_activated_FOS_high", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_actT_FOShi_try) <- paste(i, colnames(UCINI_actT_FOShi_try), sep = "_") 
  UCINI_actT_FOShi_try$Gene <- rownames(UCINI_actT_FOShi_try)
  UC_actT_FOShi_Pvalues <- merge(UC_actT_FOShi_Pvalues, UCINI_actT_FOShi_try[,c(5,6)], by = "Gene")
  UC_actT_FOShi_logFC <- merge(UC_actT_FOShi_logFC, UCINI_actT_FOShi_try[,c(2,6)], by = "Gene")
}
write.csv(UC_actT_FOShi_Pvalues, file="Results/DE_2/UC_actT_FOShi_Pvalues.csv")
write.csv(UC_actT_FOShi_logFC, file="Results/DE_2/UC_actT_FOShi_logFC.csv")

# memT 
PSC_memT_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_memT_Pvalues) = "Gene"
PSC_memT_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_memT_logFC) = "Gene"
IDs_memT_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "CD4_T_memory"])
for (i in 1:20){
  per = sample(IDs_memT_PSC, size =480, replace=F) #size=75% of smallest group
  PSCINI_memT_try <- FindMarkers(data_PSC[,per], subset.ident = "CD4_T_memory", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_memT_try) <- paste(i, colnames(PSCINI_memT_try), sep = "_") 
  PSCINI_memT_try$Gene <- rownames(PSCINI_memT_try)
  PSC_memT_Pvalues <- merge(PSC_memT_Pvalues, PSCINI_memT_try[,c(5,6)], by = "Gene")
  PSC_memT_logFC <- merge(PSC_memT_logFC, PSCINI_memT_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_memT_Pvalues, file="Results/DE_2/PSC_memT_Pvalues.csv")
write.csv(PSC_memT_logFC, file="Results/DE_2/PSC_memT_logFC.csv")

UC_memT_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_memT_Pvalues) = "Gene"
UC_memT_logFC <-data.frame(rownames(data_UC))
colnames(UC_memT_logFC) = "Gene"
IDs_memT_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "CD4_T_memory"])
for (i in 1:20){
  per = sample(IDs_memT_UC, size =480, replace=F) #size=75% of smallest group
  UCINI_memT_try <- FindMarkers(data_UC[,per], subset.ident = "CD4_T_memory", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_memT_try) <- paste(i, colnames(UCINI_memT_try), sep = "_") 
  UCINI_memT_try$Gene <- rownames(UCINI_memT_try)
  UC_memT_Pvalues <- merge(UC_memT_Pvalues, UCINI_memT_try[,c(5,6)], by = "Gene")
  UC_memT_logFC <- merge(UC_memT_logFC, UCINI_memT_try[,c(2,6)], by = "Gene")
}
write.csv(UC_memT_Pvalues, file="Results/DE_2/UC_memT_Pvalues.csv")
write.csv(UC_memT_logFC, file="Results/DE_2/UC_memT_logFC.csv")

# Treg 
PSC_Treg_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_Treg_Pvalues) = "Gene"
PSC_Treg_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_Treg_logFC) = "Gene"
IDs_Treg_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Treg"])
for (i in 1:20){
  per = sample(IDs_Treg_PSC, size =349, replace=F) #size=75% of smallest group
  PSCINI_Treg_try <- FindMarkers(data_PSC[,per], subset.ident = "Treg", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_Treg_try) <- paste(i, colnames(PSCINI_Treg_try), sep = "_") 
  PSCINI_Treg_try$Gene <- rownames(PSCINI_Treg_try)
  PSC_Treg_Pvalues <- merge(PSC_Treg_Pvalues, PSCINI_Treg_try[,c(5,6)], by = "Gene")
  PSC_Treg_logFC <- merge(PSC_Treg_logFC, PSCINI_Treg_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_Treg_Pvalues, file="Results/DE_2/PSC_Treg_Pvalues.csv")
write.csv(PSC_Treg_logFC, file="Results/DE_2/PSC_Treg_logFC.csv")

UC_Treg_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_Treg_Pvalues) = "Gene"
UC_Treg_logFC <-data.frame(rownames(data_UC))
colnames(UC_Treg_logFC) = "Gene"
IDs_Treg_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Treg"])
for (i in 1:20){
  per = sample(IDs_Treg_UC, size =349, replace=F) #size=75% of smallest group
  UCINI_Treg_try <- FindMarkers(data_UC[,per], subset.ident = "Treg", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_Treg_try) <- paste(i, colnames(UCINI_Treg_try), sep = "_") 
  UCINI_Treg_try$Gene <- rownames(UCINI_Treg_try)
  UC_Treg_Pvalues <- merge(UC_Treg_Pvalues, UCINI_Treg_try[,c(5,6)], by = "Gene")
  UC_Treg_logFC <- merge(UC_Treg_logFC, UCINI_Treg_try[,c(2,6)], by = "Gene")
}
write.csv(UC_Treg_Pvalues, file="Results/DE_2/UC_Treg_Pvalues.csv")
write.csv(UC_Treg_logFC, file="Results/DE_2/UC_Treg_logFC.csv")

# cyclingT 
PSC_cyclingT_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_cyclingT_Pvalues) = "Gene"
PSC_cyclingT_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_cyclingT_logFC) = "Gene"
IDs_cyclingT_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Cycling_T"])
for (i in 1:20){
  per = sample(IDs_cyclingT_PSC, size =16, replace=F) #size=75% of smallest group
  PSCINI_cyclingT_try <- FindMarkers(data_PSC[,per], subset.ident = "Cycling_T", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_cyclingT_try) <- paste(i, colnames(PSCINI_cyclingT_try), sep = "_") 
  PSCINI_cyclingT_try$Gene <- rownames(PSCINI_cyclingT_try)
  PSC_cyclingT_Pvalues <- merge(PSC_cyclingT_Pvalues, PSCINI_cyclingT_try[,c(5,6)], by = "Gene")
  PSC_cyclingT_logFC <- merge(PSC_cyclingT_logFC, PSCINI_cyclingT_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_cyclingT_Pvalues, file="Results/DE_2/PSC_cyclingT_Pvalues.csv")
write.csv(PSC_cyclingT_logFC, file="Results/DE_2/PSC_cyclingT_logFC.csv")

UC_cyclingT_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_cyclingT_Pvalues) = "Gene"
UC_cyclingT_logFC <-data.frame(rownames(data_UC))
colnames(UC_cyclingT_logFC) = "Gene"
IDs_cyclingT_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Cycling_T"])
for (i in 1:20){
  per = sample(IDs_cyclingT_UC, size =16, replace=F) #size=75% of smallest group
  UCINI_cyclingT_try <- FindMarkers(data_UC[,per], subset.ident = "Cycling_T", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_cyclingT_try) <- paste(i, colnames(UCINI_cyclingT_try), sep = "_") 
  UCINI_cyclingT_try$Gene <- rownames(UCINI_cyclingT_try)
  UC_cyclingT_Pvalues <- merge(UC_cyclingT_Pvalues, UCINI_cyclingT_try[,c(5,6)], by = "Gene")
  UC_cyclingT_logFC <- merge(UC_cyclingT_logFC, UCINI_cyclingT_try[,c(2,6)], by = "Gene")
}
write.csv(UC_cyclingT_Pvalues, file="Results/DE_2/UC_cyclingT_Pvalues.csv")
write.csv(UC_cyclingT_logFC, file="Results/DE_2/UC_cyclingT_logFC.csv")

# APC 
PSC_APC_Pvalues <-data.frame(rownames(data_PSC))
colnames(PSC_APC_Pvalues) = "Gene"
PSC_APC_logFC <-data.frame(rownames(data_PSC))
colnames(PSC_APC_logFC) = "Gene"
IDs_APC_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "APC"])
for (i in 1:20){
  per = sample(IDs_APC_PSC, size =510, replace=F) #size=75% of smallest group
  PSCINI_APC_try <- FindMarkers(data_PSC[,per], subset.ident = "APC", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_APC_try) <- paste(i, colnames(PSCINI_APC_try), sep = "_") 
  PSCINI_APC_try$Gene <- rownames(PSCINI_APC_try)
  PSC_APC_Pvalues <- merge(PSC_APC_Pvalues, PSCINI_APC_try[,c(5,6)], by = "Gene")
  PSC_APC_logFC <- merge(PSC_APC_logFC, PSCINI_APC_try[,c(2,6)], by = "Gene")
}
write.csv(PSC_APC_Pvalues, file="Results/DE_2/PSC_APC_Pvalues.csv")
write.csv(PSC_APC_logFC, file="Results/DE_2/PSC_APC_logFC.csv")

UC_APC_Pvalues <-data.frame(rownames(data_UC))
colnames(UC_APC_Pvalues) = "Gene"
UC_APC_logFC <-data.frame(rownames(data_UC))
colnames(UC_APC_logFC) = "Gene"
IDs_APC_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "APC"])
for (i in 1:20){
  per = sample(IDs_APC_UC, size =510, replace=F) #size=75% of smallest group
  UCINI_APC_try <- FindMarkers(data_UC[,per], subset.ident = "APC", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_APC_try) <- paste(i, colnames(UCINI_APC_try), sep = "_") 
  UCINI_APC_try$Gene <- rownames(UCINI_APC_try)
  UC_APC_Pvalues <- merge(UC_APC_Pvalues, UCINI_APC_try[,c(5,6)], by = "Gene")
  UC_APC_logFC <- merge(UC_APC_logFC, UCINI_APC_try[,c(2,6)], by = "Gene")
}
write.csv(UC_APC_Pvalues, file="Results/DE_2/UC_APC_Pvalues.csv")
write.csv(UC_APC_logFC, file="Results/DE_2/UC_APC_logFC.csv")

# Vulcanoplot activated T FOS low
PSC_actT_FOSlow_Pvalues$PSC_mean_Pvalue <- rowMeans(PSC_actT_FOSlow_Pvalues[,-1])
PSC_actT_FOSlow_logFC$PSC_mean_logFC <- rowMeans(PSC_actT_FOSlow_logFC[,-1])
UC_actT_FOSlow_Pvalues$UC_mean_Pvalue <- rowMeans(UC_actT_FOSlow_Pvalues[,-1])
UC_actT_FOSlow_logFC$UC_mean_logFC <- rowMeans(UC_actT_FOSlow_logFC[,-1])

actT_FOSlow <- merge(PSC_actT_FOSlow_Pvalues[,c(1,22)], PSC_actT_FOSlow_logFC[,c(1,22)], by = "Gene")
actT_FOSlow <- merge(actT_FOSlow, UC_actT_FOSlow_Pvalues[,c(1,22)], by = "Gene")
actT_FOSlow <- merge(actT_FOSlow, UC_actT_FOSlow_logFC[,c(1,22)], by = "Gene")

actT_FOSlow$both_sign <- ifelse((actT_FOSlow$UC_mean_Pvalue < 0.05 & actT_FOSlow$PSC_mean_Pvalue < 0.05), "yes", "no")

ggplot(data=actT_FOSlow, aes(x=UC_mean_logFC, y=-log10(UC_mean_Pvalue), col=both_sign, alpha = both_sign, label=Gene)) + 
  geom_text(data=subset(actT_FOSlow, both_sign == "no" & (UC_mean_Pvalue < 0.05 | PSC_mean_Pvalue < 0.05), overlap = T)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("grey0","red3")) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  ggplot(data=actT_FOSlow, aes(x=PSC_mean_logFC, y=-log10(PSC_mean_Pvalue), col=both_sign,  alpha = both_sign, label=Gene)) + 
  geom_text(data=subset(actT_FOSlow, both_sign == "no" & (UC_mean_Pvalue < 0.05 | PSC_mean_Pvalue < 0.05), overlap = T)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("grey0","red3")) +
  scale_alpha_discrete(range = c(0.3, 1))

EnhancedVolcano(actT_FOSlow, lab = actT_FOSlow$Gene, x = 'UC_mean_logFC', y = 'UC_mean_Pvalue', pCutoff = 0.05, FCcutoff = 0.5) +
  EnhancedVolcano(actT_FOSlow, lab = actT_FOSlow$Gene, x = 'PSC_mean_logFC', y = 'PSC_mean_Pvalue', pCutoff = 0.05, FCcutoff = 0.5)

# Vulcanoplot IgG
PSC_IgG_Pvalues$PSC_mean_Pvalue <- rowMeans(PSC_IgG_Pvalues[,-1])
PSC_IgG_logFC$PSC_mean_logFC <- rowMeans(PSC_IgG_logFC[,-1])
UC_IgG_Pvalues$UC_mean_Pvalue <- rowMeans(UC_IgG_Pvalues[,-1])
UC_IgG_logFC$UC_mean_logFC <- rowMeans(UC_IgG_logFC[,-1])

IgG <- merge(PSC_IgG_Pvalues[,c(1,22)], PSC_IgG_logFC[,c(1,22)], by = "Gene")
IgG <- merge(IgG, UC_IgG_Pvalues[,c(1,22)], by = "Gene")
IgG <- merge(IgG, UC_IgG_logFC[,c(1,22)], by = "Gene")

IgG$both_sign <- ifelse((IgG$UC_mean_Pvalue < 0.05 & IgG$PSC_mean_Pvalue < 0.05), "yes", "no")

ggplot(data=IgG, aes(x=UC_mean_logFC, y=-log10(UC_mean_Pvalue), col=both_sign, alpha = both_sign, label=Gene)) + 
  geom_text(data=subset(IgG, UC_mean_Pvalue < 1e-10 | PSC_mean_Pvalue < 1e-10), check_overlap = F, size = 2.5) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("grey0","red3")) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  xlim(-2, 2) +
  ylim(-1,35) +
  ggplot(data=IgG, aes(x=PSC_mean_logFC, y=-log10(PSC_mean_Pvalue), col=both_sign,  alpha = both_sign, label=Gene)) + 
  geom_text(data=subset(IgG, UC_mean_Pvalue < 1e-10 | PSC_mean_Pvalue < 1e-10), check_overlap = T, size = 2.5) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("grey0","red3")) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  xlim(-2, 2) +
  ylim(-1,35)

ggsave("Results/Figures/Vulcano_IgG.pdf")

EnhancedVolcano(IgG, lab = IgG$Gene, x = 'UC_mean_logFC', y = 'UC_mean_Pvalue', pCutoff = 0.05, FCcutoff = 0.5) +
  EnhancedVolcano(IgG, lab = IgG$Gene, x = 'PSC_mean_logFC', y = 'PSC_mean_Pvalue', pCutoff = 0.05, FCcutoff = 0.5)

# Vulcanoplot IgA
PSC_IgA_Pvalues$PSC_mean_Pvalue <- rowMeans(PSC_IgA_Pvalues[,-1])
PSC_IgA_logFC$PSC_mean_logFC <- rowMeans(PSC_IgA_logFC[,-1])
UC_IgA_Pvalues$UC_mean_Pvalue <- rowMeans(UC_IgA_Pvalues[,-1])
UC_IgA_logFC$UC_mean_logFC <- rowMeans(UC_IgA_logFC[,-1])

IgA <- merge(PSC_IgA_Pvalues[,c(1,22)], PSC_IgA_logFC[,c(1,22)], by = "Gene")
IgA <- merge(IgA, UC_IgA_Pvalues[,c(1,22)], by = "Gene")
IgA <- merge(IgA, UC_IgA_logFC[,c(1,22)], by = "Gene")

IgA$both_sign <- ifelse((IgA$UC_mean_Pvalue < 0.05 & IgA$PSC_mean_Pvalue < 0.05), "yes", "no")

ggplot(data=IgA, aes(x=UC_mean_logFC, y=-log10(UC_mean_Pvalue), col=both_sign, alpha = both_sign, label=Gene)) + 
  geom_text(data=subset(IgA, UC_mean_Pvalue < 1e-10 | PSC_mean_Pvalue < 1e-10), check_overlap = F, size = 2.5) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("grey0","red3")) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  xlim(-2, 2) +
  ylim(-1,40) +
  ggplot(data=IgA, aes(x=PSC_mean_logFC, y=-log10(PSC_mean_Pvalue), col=both_sign,  alpha = both_sign, label=Gene)) + 
  geom_text(data=subset(IgA, UC_mean_Pvalue < 1e-10 | PSC_mean_Pvalue < 1e-10), check_overlap = T, size = 2.5) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("grey0","red3")) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  xlim(-2, 2) +
  ylim(-1,40)

ggsave("Results/Figures/Vulcano_IgA.pdf")

EnhancedVolcano(IgA, lab = IgA$Gene, x = 'UC_mean_logFC', y = 'UC_mean_Pvalue', pCutoff = 0.05, FCcutoff = 0.5) +
  EnhancedVolcano(IgA, lab = IgA$Gene, x = 'PSC_mean_logFC', y = 'PSC_mean_Pvalue', pCutoff = 0.05, FCcutoff = 0.5)

# Vulcanoplot Follicular B cells
PSC_Follicular_Pvalues$PSC_mean_Pvalue <- rowMeans(PSC_Follicular_Pvalues[,-1])
PSC_Follicular_logFC$PSC_mean_logFC <- rowMeans(PSC_Follicular_logFC[,-1])
UC_Follicular_Pvalues$UC_mean_Pvalue <- rowMeans(UC_Follicular_Pvalues[,-1])
UC_Follicular_logFC$UC_mean_logFC <- rowMeans(UC_Follicular_logFC[,-1])

Follicular <- merge(PSC_Follicular_Pvalues[,c(1,22)], PSC_Follicular_logFC[,c(1,22)], by = "Gene")
Follicular <- merge(Follicular, UC_Follicular_Pvalues[,c(1,22)], by = "Gene")
Follicular <- merge(Follicular, UC_Follicular_logFC[,c(1,22)], by = "Gene")

Follicular$both_sign <- ifelse((Follicular$UC_mean_Pvalue < 0.05 & Follicular$PSC_mean_Pvalue < 0.05), "yes", "no")

ggplot(data=Follicular, aes(x=UC_mean_logFC, y=-log10(UC_mean_Pvalue), col=both_sign, alpha = both_sign, label=Gene)) + 
  geom_text(data=subset(Follicular, UC_mean_Pvalue < 1e-3 | PSC_mean_Pvalue < 1e-3), check_overlap = F, size = 2.5) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("grey0","red3")) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  xlim(-1, 1) +
  ylim(-1,15) +
  ggplot(data=Follicular, aes(x=PSC_mean_logFC, y=-log10(PSC_mean_Pvalue), col=both_sign,  alpha = both_sign, label=Gene)) + 
  geom_text(data=subset(Follicular, UC_mean_Pvalue < 1e-3 | PSC_mean_Pvalue < 1e-3), check_overlap = T, size = 2.5) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("grey0","red3")) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  xlim(-1, 1) +
  ylim(-1,15)

ggsave("Results/Figures/Vulcano_Follicular.pdf")


