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
PSCINI_IgG <- FindMarkers(data_PSC, subset.ident = "IgG_plasma", group.by = "inflammation", test.use = "MAST", ident.1 = "I", ident.2 = "NI")
PSCINI_IgG <- filter(PSCINI_IgG, PSCINI_IgG$p_val_adj < 0.05)
PSCINI_IgG$gene = rownames(PSCINI_IgG)
UCINI_IgG <- FindMarkers(data_UC, subset.ident = "IgG_plasma", group.by = "inflammation", test.use = "MAST", ident.1 = "I", ident.2 = "NI")
UCINI_IgG <- filter(UCINI_IgG, UCINI_IgG$p_val_adj < 0.05)
UCINI_IgG$gene = rownames(UCINI_IgG)

Genes_IgG = unique(c(UCINI_IgG$gene, PSCINI_IgG$gene))
PSCINI_IgG_selected <- FindMarkers(data_PSC, subset.ident = "IgG_plasma", features = Genes_IgG, group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
UCINI_IgG_selected <- FindMarkers(data_UC, subset.ident = "IgG_plasma", features = Genes_IgG, group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
colnames(UCINI_IgG_selected) <- paste("UC", colnames(UCINI_IgG_selected), sep = "_")
colnames(PSCINI_IgG_selected) <- paste("PSC", colnames(PSCINI_IgG_selected), sep = "_")
IgG_selected <- merge(UCINI_IgG_selected, PSCINI_IgG_selected, by = "row.names")
colnames(IgG_selected)[1] = "Gene"

IgG_selected$both_sign = ifelse((IgG_selected$UC_p_val_adj < 0.05 & IgG_selected$PSC_p_val_adj < 0.05), "yes", "no")
IgG_selected$both_sign = as.factor(IgG_selected$both_sign)

# IgG all genes
PSCINI_IgG_all <- FindMarkers(data_PSC, subset.ident = "IgG_plasma", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
UCINI_IgG_all <- FindMarkers(data_UC, subset.ident = "IgG_plasma", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
colnames(UCINI_IgG_all) <- paste("UC", colnames(UCINI_IgG_all), sep = "_")
colnames(PSCINI_IgG_all) <- paste("PSC", colnames(PSCINI_IgG_all), sep = "_")
IgG_all <- merge(UCINI_IgG_all, PSCINI_IgG_all, by = "row.names")
colnames(IgG_all)[1] = "Gene"

IgG_all$both_sign = ifelse((IgG_all$UC_p_val_adj < 0.05 & IgG_all$PSC_p_val_adj < 0.05), "yes", "no")
IgG_all$both_sign = as.factor(IgG_all$both_sign)

# Treg all genes
PSCINI_Treg_all <- FindMarkers(data_PSC, subset.ident = "Treg", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
UCINI_Treg_all <- FindMarkers(data_UC, subset.ident = "Treg", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
colnames(UCINI_Treg_all) <- paste("UC", colnames(UCINI_Treg_all), sep = "_")
colnames(PSCINI_Treg_all) <- paste("PSC", colnames(PSCINI_Treg_all), sep = "_")
Treg_all <- merge(UCINI_Treg_all, PSCINI_Treg_all, by = "row.names")
colnames(Treg_all)[1] = "Gene"

Treg_all$both_sign = ifelse((Treg_all$UC_p_val_adj < 0.05 & Treg_all$PSC_p_val_adj < 0.05), "yes", "no")
Treg_all$both_sign = as.factor(Treg_all$both_sign)

# Absorptive enterocyte all genes
PSCINI_absentero_all <- FindMarkers(data_PSC, subset.ident = "Absorptive_enterocyte", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
UCINI_absentero_all <- FindMarkers(data_UC, subset.ident = "Absorptive_enterocyte", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
colnames(UCINI_absentero_all) <- paste("UC", colnames(UCINI_absentero_all), sep = "_")
colnames(PSCINI_absentero_all) <- paste("PSC", colnames(PSCINI_absentero_all), sep = "_")
absentero_all <- merge(UCINI_absentero_all, PSCINI_absentero_all, by = "row.names")
colnames(absentero_all)[1] = "Gene"

absentero_all$both_sign = ifelse((absentero_all$UC_p_val_adj < 0.05 & absentero_all$PSC_p_val_adj < 0.05), "yes", "no")
absentero_all$both_sign = as.factor(absentero_all$both_sign)

# Absorptive enterocyte all genes
PSCINI_APC_all <- FindMarkers(data_PSC, subset.ident = "APC", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
UCINI_APC_all <- FindMarkers(data_UC, subset.ident = "APC", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
colnames(UCINI_APC_all) <- paste("UC", colnames(UCINI_APC_all), sep = "_")
colnames(PSCINI_APC_all) <- paste("PSC", colnames(PSCINI_APC_all), sep = "_")
APC_all <- merge(UCINI_APC_all, PSCINI_APC_all, by = "row.names")
colnames(APC_all)[1] = "Gene"

APC_all$both_sign = ifelse((APC_all$UC_p_val_adj < 0.05 & APC_all$PSC_p_val_adj < 0.05), "yes", "no")
APC_all$both_sign = as.factor(APC_all$both_sign)


# correlation plot
ggplot(IgG_all, aes(x=-log10(UC_p_val_adj), y=-log10(PSC_p_val_adj))) + 
  geom_point() +
  geom_abline(slope=1, intercept=0) + 
  theme_classic()

# vulcano plot
EnhancedVolcano(absentero_all, lab = absentero_all$Gene, x = 'UC_avg_logFC', y = 'UC_p_val_adj') +
EnhancedVolcano(absentero_all, lab = absentero_all$Gene, x = 'PSC_avg_logFC', y = 'PSC_p_val_adj')

ggplot(data=absentero_all, aes(x=UC_avg_logFC, y=-log10(UC_p_val_adj), col=both_sign, alpha = both_sign, label=Gene)) + 
  geom_text(data=subset(absentero_all, both_sign == "no" & (UC_p_val_adj < 0.05 | PSC_p_val_adj < 0.05), overlap = T)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("grey0","red3")) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  ggplot(data=absentero_all, aes(x=PSC_avg_logFC, y=-log10(PSC_p_val_adj), col=both_sign,  alpha = both_sign, label=Gene)) + 
  geom_text(data=subset(absentero_all, both_sign == "no" & (UC_p_val_adj < 0.05 | PSC_p_val_adj < 0.05), overlap = T)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("grey0","red3")) +
  scale_alpha_discrete(range = c(0.3, 1))

# permutation
PSC_absentero <-data.frame(rownames(data_PSC))
colnames(PSC_absentero) = "Gene"
IDs_absentero_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Absorptive_enterocyte"])
for (i in 1:100){
  per = sample(IDs_absentero_PSC, size =400, replace=F)
  PSCINI_absentero_try <- FindMarkers(data_PSC[,per], subset.ident = "Absorptive_enterocyte", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_absentero_try) <- paste(i, colnames(PSCINI_absentero_try), sep = "_") 
  PSCINI_absentero_try$Gene <- rownames(PSCINI_absentero_try)
  PSC_absentero <- merge(PSC_absentero, PSCINI_absentero_try[,c(2,5,6)], by = "Gene")
}

UC_absentero <-data.frame(rownames(data_UC))
colnames(UC_absentero) = "Gene"
IDs_absentero_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Absorptive_enterocyte"])
for (i in 1:100){
  per = sample(IDs_absentero_UC, size =400, replace=F)
  UCINI_absentero_try <- FindMarkers(data_UC[, per], subset.ident = "Absorptive_enterocyte", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_absentero_try) <- paste(i, colnames(UCINI_absentero_try), sep = "_") 
  UCINI_absentero_try$Gene <- rownames(UCINI_absentero_try)
  UC_absentero <- merge(UC_absentero, UCINI_absentero_try[,c(2,5,6)], by = "Gene")
}

# permutation
PSC_follicular <-data.frame(rownames(data_PSC))
colnames(PSC_follicular) = "Gene"
IDs_follicular_PSC=colnames(data_PSC[,data_PSC@meta.data[["celltypes"]] == "Follicular_B"])
for (i in 1:100){
  per = sample(IDs_follicular_PSC, size =1200, replace=F)
  PSCINI_follicular_try <- FindMarkers(data_PSC[,per], subset.ident = "Follicular_B", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(PSCINI_follicular_try) <- paste(i, colnames(PSCINI_follicular_try), sep = "_") 
  PSCINI_follicular_try$Gene <- rownames(PSCINI_follicular_try)
  PSC_follicular <- merge(PSC_follicular, PSCINI_follicular_try[,c(2,5,6)], by = "Gene")
}

UC_follicular <-data.frame(rownames(data_UC))
colnames(UC_follicular) = "Gene"
IDs_follicular_UC=colnames(data_UC[,data_UC@meta.data[["celltypes"]] == "Follicular_B"])
for (i in 1:100){
  per = sample(IDs_follicular_UC, size =1200, replace=F)
  UCINI_follicular_try <- FindMarkers(data_UC[, per], subset.ident = "Follicular_B", group.by = "inflammation", test.use = "MAST", ident.1 = "NI", ident.2 = "I", logfc.threshold = 0, min.pct = 0)
  colnames(UCINI_follicular_try) <- paste(i, colnames(UCINI_follicular_try), sep = "_") 
  UCINI_follicular_try$Gene <- rownames(UCINI_follicular_try)
  UC_follicular <- merge(UC_follicular, UCINI_follicular_try[,c(2,5,6)], by = "Gene")
}

# permutation vulcanoplot
p_values_columns = grepl("p_val_adj",names(PSC_absentero))
p_values_columns[1] = TRUE
PSC_absentero_Pvalues=PSC_absentero[,p_values_columns]
logFC_columns = grepl("logFC",names(PSC_absentero))
logFC_columns[1] = TRUE
PSC_absentero_logFC=PSC_absentero[,logFC_columns]

PSC_absentero_logFC$PSC_mean_logFC <- rowMeans(PSC_absentero_logFC[,-1])
PSC_absentero_logFC = transform(PSC_absentero_logFC, PSC_SD_logFC=apply(PSC_absentero_logFC[,-1],1, sd, na.rm = TRUE))

PSC_absentero_Pvalues$PSC_mean_Pvalues <- rowMeans(PSC_absentero_Pvalues[,-1])
PSC_absentero_Pvalues = transform(PSC_absentero_Pvalues, PSC_SD_Pvalues=apply(PSC_absentero_Pvalues[,-1],1, sd, na.rm = TRUE))

PSC_absentero_overview = merge(PSC_absentero_Pvalues[,c(1,102,103)], PSC_absentero_logFC[,c(1,102,103)], by = "Gene") 

p_values_columns = grepl("p_val_adj",names(UC_absentero))
p_values_columns[1] = TRUE
UC_absentero_Pvalues=UC_absentero[,p_values_columns]
logFC_columns = grepl("logFC",names(UC_absentero))
logFC_columns[1] = TRUE
UC_absentero_logFC=UC_absentero[,logFC_columns]

UC_absentero_logFC$UC_mean_logFC <- rowMeans(UC_absentero_logFC[,-1])
UC_absentero_logFC = transform(UC_absentero_logFC, UC_SD_logFC=apply(UC_absentero_logFC[,-1],1, sd, na.rm = TRUE))

UC_absentero_Pvalues$UC_mean_Pvalues <- rowMeans(UC_absentero_Pvalues[,-1])
UC_absentero_Pvalues = transform(UC_absentero_Pvalues, UC_SD_Pvalues=apply(UC_absentero_Pvalues[,-1],1, sd, na.rm = TRUE))

UC_absentero_overview = merge(UC_absentero_Pvalues[,c(1,102,103)], UC_absentero_logFC[,c(1,102,103)], by = "Gene") 

Absentero_overview = merge(UC_absentero_overview, PSC_absentero_overview, by = "Gene", all = T)

Absentero_overview$both_sign = ifelse((Absentero_overview$UC_mean_Pvalues < 0.05 & Absentero_overview$PSC_mean_Pvalues < 0.05), "yes", "no")
Absentero_overview$both_sign = as.factor(Absentero_overview$both_sign)

