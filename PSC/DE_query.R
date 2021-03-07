#### implicated genes literature

genes<-read.csv("/genes_junctions_extended.csv", sep=";")

# complication: name needs to be overlapping with ~/Downloads/DE/*.csv name
setwd("/DE/")
temp = list.files(pattern="^PSC*")
i<-temp[1]
x<-read.csv(i)
x$filename<-i
x<-merge(genes, x, by="gene", all=F)

for(i in temp){
y<-read.csv(i)
y$filename<-i
y<-merge(genes, y, by="gene", all=F)
x<-rbind(x,y)}
x<-x[-c(1:4),]
write.csv(x, "/DE/Results/interesting_genes_PSC.csv")

temp = list.files(pattern="^UC*")
i<-temp[1]
x<-read.csv(i)
x$filename<-i
x<-merge(genes, x, by="gene", all=F)

for(i in temp){
  y<-read.csv(i)
  y$filename<-i
  y<-merge(genes, y, by="gene", all=F)
  x<-rbind(x,y)}
write.csv(x, "/DE/Results/interesting_genes_UC.csv")

temp = list.files(pattern="^HC*")
i<-temp[1]
x<-read.csv(i)
x$filename<-i
x<-merge(genes, x, by="gene", all=F)

for(i in temp){
  y<-read.csv(i)
  y$filename<-i
  y<-merge(genes, y, by="gene", all=F)
  x<-rbind(x,y)}
write.csv(x, "/DE/Results/interesting_genes_HC.csv")

temp = list.files(pattern="^IBD*")
i<-temp[1]
x<-read.csv(i)
x$filename<-i
x<-merge(genes, x, by="gene", all=F)

for(i in temp){
  y<-read.csv(i)
  y$filename<-i
  y<-merge(genes, y, by="gene", all=F)
  x<-rbind(x,y)}
x<-x[-c(1:2),]
write.csv(x, "/DE/Results/interesting_genes_IBD.csv")


#####
# Drug Targets

DefaultAssay(data)<-"RNA"
dim(data)
allgenes <- data.frame(matrix(ncol = 2, nrow = 33538))
colnames(allgenes)[1]<-"Gene"
allgenes$Gene<-rownames(data@assays$RNA)

targets<-read.csv("~/Desktop/PSC/targets_associated_with_inflammatory_bowel_disease_dd17012021.csv")
colnames(targets)[1]<-"Gene"
targets_incuded<-merge(targets, allgenes, by="Gene", all=F) #4332/4407 included in our dataset

targets_incuded<-targets_incuded[-12]
# filter those celltypes with >100 cells per group
Idents(data) <- "celltypes"
celltypes <- levels(data$celltypes)
table(data$celltypes, data$state)
celltypes_allcomparisons <- celltypes[-c(2:4,6:9,12:14,17:23, 25:27, 33,34, 36, 38:43)] #filter out celltypes with too little cells
celltypes_PSC_and_UCNI_vs_HC <- celltypes[-c(2:4,6:9,12:14,17:23, 25:27, 33,34, 36, 38, 40:43)] #filter out celltypes with too little cells
celltypes_UC_and_PSCNI_vs_HC <- celltypes[-c(2:4,6:9,12:14,17:23, 25:27, 33, 36, 39:43)] #filter out celltypes with too little cells

DefaultAssay(data)<-"RNA"

PSC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(PSC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")

for(i in 1:16){
  DEgenes <- FindMarkers(data, subset.ident = celltypes_UC_and_PSCNI_vs_HC[i], group.by = "state", test.use = "MAST", ident.1 = "PSC-NI", ident.2 = "HC-NI")
  DEgenes$celltype <- celltypes_UC_and_PSCNI_vs_HC[i]
  DEgenes$Gene <- rownames(DEgenes)
  DEgenes <- merge(DEgenes, targets_incuded, by = "Gene")
  PSC <- rbind(PSC, DEgenes)
}

UC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(UC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")
for(i in 1:16){
  DEgenes <- FindMarkers(data, subset.ident = celltypes_UC_and_PSCNI_vs_HC[i], group.by = "state", test.use = "MAST", ident.1 = "UC-NI", ident.2 = "HC-NI")
  DEgenes$celltype <- celltypes_UC_and_PSCNI_vs_HC[i]
  #DEgenes1 <- DEgenes[DEgenes$p_val_adj < 0.05,]
  DEgenes$Gene <- rownames(DEgenes)
  DEgenes <- merge(DEgenes, targets_incuded, by = "Gene")
  UC <- rbind(UC, DEgenes)
}

UC_NI<-UC
PSC_NI<-PSC
write.csv(UC, "~/Desktop/PSC/202103_analyses/DE/Results/UC_drugtarget_genes_celltypes.csv")
write.csv(PSC, "~/Desktop/PSC/202103_analyses/DE/Results/PSC_drugtarget_genes_celltypes.csv")

PSC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(PSC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")
for(i in 1:15){
  DEgenes <- FindMarkers(data, subset.ident = celltypes_PSC_and_UCNI_vs_HC[i], group.by = "state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "HC-NI")
  DEgenes$celltype <- celltypes_PSC_and_UCNI_vs_HC[i]
  #DEgenes1 <- DEgenes[DEgenes$p_val_adj < 0.05,]
  DEgenes$Gene <- rownames(DEgenes)
  DEgenes <- merge(DEgenes, targets_incuded, by = "Gene")
  PSC <- rbind(PSC, DEgenes)
}

UC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(UC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")
for(i in 1:16){
  DEgenes <- FindMarkers(data, subset.ident = celltypes_UC_and_PSCNI_vs_HC[i], group.by = "state", test.use = "MAST", ident.1 = "UC-I", ident.2 = "HC-NI")
  DEgenes$celltype <- celltypes_UC_and_PSCNI_vs_HC[i]
  #DEgenes1 <- DEgenes[DEgenes$p_val_adj < 0.05,]
  DEgenes$Gene <- rownames(DEgenes)
  DEgenes <- merge(DEgenes, targets_incuded, by = "Gene")
  UC <- rbind(UC, DEgenes)
}




UC_I<-UC
PSC_I<-PSC

PSC_I$p_val_adj<-PSC_I$p_val*4332
PSC_NI$p_val_adj<-PSC_NI$p_val*4332

UC_I$p_val_adj<-UC_I$p_val*4332
UC_NI$p_val_adj<-UC_NI$p_val*4332

write.csv(UC_I, "~/Desktop/PSC/202103_analyses/DE/Results/UC_I_bonferroni_drugtarget_genes_celltypes.csv")
write.csv(PSC_I, "~/Desktop/PSC/202103_analyses/DE/Results/PSC_I_bonferroni_drugtarget_genes_celltypes.csv")
write.csv(UC_NI, "~/Desktop/PSC/202103_analyses/DE/Results/UC_NI_bonferroni_drugtarget_genes_celltypes.csv")
write.csv(PSC_NI, "~/Desktop/PSC/202103_analyses/DE/Results/PSC_NI_bonferroni_drugtarget_genes_celltypes.csv")

#####
#PSC risk genes


data<- readRDS("~/Desktop/PSC_processed_march_2021.rds")
DefaultAssay(data) = "RNA"
UC_NI<-subset(data, subset = state == "UC-NI")
saveRDS(UC_NI, "~/Desktop/PSC/202101analyses/UC_NI.rds")
UC_risk<-read.csv("~/Desktop/PSC/202101analyses/ucputativeriskgenes_smillie.csv", sep=";")
genes<-read.csv("~/Desktop/PSC/202010_analyses/genes_junctions_extended.csv", sep=";")
PSC_risk<-genes[genes$function. == "PSC_and_suggestive_risk_genes",]
DefaultAssay(UC_NI) = "RNA"
DotPlot(UC_NI, features = UC_risk$Gene) # 3 genes not found. save image 8*44inch pdf
DotPlot(UC_NI, features = PSC_risk$gene) # 5 genes not found save image 8*30inch pdf

PSC_NI<-subset(data, subset = state == "PSC-NI")
saveRDS(PSC_NI, "~/Desktop/PSC/202101analyses/PSC_NI.rds")
DefaultAssay(PSC_NI) = "RNA"
DotPlot(PSC_NI, features = UC_risk$Gene) # 3 genes not found. save image 8*44inch pdf
DotPlot(PSC_NI, features = PSC_risk$gene) # 5 genes not found save image 8*30inch pdf

HC_NI<-subset(data, subset = state == "HC-NI")
saveRDS(HC_NI, "~/Desktop/PSC/202101analyses/HC_NI.rds")
DefaultAssay(HC_NI) = "RNA"
DotPlot(HC_NI, features = UC_risk$Gene) # 3 genes not found. save image 8*44inch pdf
DotPlot(HC_NI, features = PSC_risk$gene) # 5 genes not found save image 8*30inch pdf

PSC_I<-subset(data, subset = state == "PSC-I")
saveRDS(PSC_I, "~/Desktop/PSC/202101analyses/PSC_I.rds")
DefaultAssay(PSC_I) = "RNA"
DotPlot(PSC_I, features = UC_risk$Gene) # 3 genes not found. save image 8*44inch pdf
DotPlot(PSC_I, features = PSC_risk$gene) # 5 genes not found save image 8*30inch pdf

UC_I<-subset(data, subset = state == "UC-I")
saveRDS(UC_I, "~/Desktop/PSC/202101analyses/UC_I.rds")
DefaultAssay(UC_I) = "RNA"
DotPlot(UC_I, features = UC_risk$Gene) # 3 genes not found. save image 8*44inch pdf
DotPlot(UC_I, features = PSC_risk$gene) # 5 genes not found save image 8*30inch pdf

inflamm_fibro<-subset(data, subset = celltypes == "Inflammatory_fibroblast")
DefaultAssay(inflamm_fibro) = "RNA"
VlnPlot(inflamm_fibro, split.by="state", c("FAP","TWIST1", "WNT2"), split.plot = TRUE)

#script amber

Idents(data) <- "celltypes"
celltypes <- levels(data$celltypes)
table(data$celltypes, data$state)
celltypes_allcomparisons <- celltypes[-c(2:4,6:9,12:14,17:23, 25:27, 33,34, 36, 38:43)] #filter out celltypes with too little cells
celltypes_PSC_and_UCNI_vs_HC <- celltypes[-c(2:4,6:9,12:14,17:23, 25:27, 33,34, 36, 38, 40:43)] #filter out celltypes with too little cells
celltypes_UC_and_PSCNI_vs_HC <- celltypes[-c(2:4,6:9,12:14,17:23, 25:27, 33, 36, 39:43)] #filter out celltypes with too little cells

DefaultAssay(data)<-"RNA"

PSC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(PSC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")
colnames(PSC_risk)[1]<-"Gene"
PSC_risk$Gene<-as.character(PSC_risk$Gene)


for(i in 1:16){
  try <- FindMarkers(data, subset.ident = celltypes_UC_and_PSCNI_vs_HC[i], group.by = "state", test.use = "MAST", ident.1 = "PSC-NI", ident.2 = "HC-NI")
  try$celltype <- celltypes_UC_and_PSCNI_vs_HC[i]
  try1 <- try[try$p_val_adj < 0.05,]
  try1$Gene <- rownames(try1)
  try1 <- merge(try1, PSC_risk, by = "Gene")
  PSC <- rbind(PSC, try1)
  }

UC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(UC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")
for(i in 1:16){
  try <- FindMarkers(data, subset.ident = celltypes_PSC_and_UCNI_vs_HC[i], group.by = "state", test.use = "MAST", ident.1 = "UC-NI", ident.2 = "HC-NI")
  try$celltype <- celltypes_PSC_and_UCNI_vs_HC[i]
  try1 <- try[try$p_val_adj < 0.05,]
  try1$Gene <- rownames(try1)
  try1 <- merge(try1, UC_risk, by = "Gene")
  UC <- rbind(UC, try1)
}

UC_genes<-unique(UC$Gene)
PSC_genes<-unique(PSC$Gene)

DotPlot(data, group.by = "state", features = UC_genes)

DotPlot(data, group.by = "celltypes", split.by = "state", features = UC_genes, cols = c("red", "green", "blue", "orange", "yellow", split.plot=T))
write.csv(UC, "/DE/Results/UC_sign_risk_genes_celltypes.csv")
write.csv(PSC, "/DE/Results/PSC_sign_risk_genes_celltypes.csv")



#####
allgenes <- data.frame(matrix(ncol = 2, nrow = 33538))
colnames(allgenes)[1]<-"Gene"
allgenes$Gene<-rownames(data@assays$RNA)
PSC_risk_in_set<-merge(allgenes, PSC_risk, by="Gene")
UC_risk_in_set<-merge(allgenes, UC_risk, by="Gene")
write.csv(UC_risk_in_set, "/DE/Results/ucputativeriskgenes_smillie_inset.csv")
write.csv(PSC_risk_in_set, "/DE/Results/PSCgenes_inset.csv")



PSC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(PSC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")

for(i in 1:16){
  try <- FindMarkers(data, subset.ident = celltypes_UC_and_PSCNI_vs_HC[i], group.by = "state", test.use = "MAST", ident.1 = "PSC-NI", ident.2 = "HC-NI")
  try$celltype <- celltypes_UC_and_PSCNI_vs_HC[i]
  #try1 <- try[try$p_val_adj < 0.05,]
  try$Gene <- rownames(try)
  try <- merge(try, PSC_risk, by = "Gene")
  PSC <- rbind(PSC, try)
}

UC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(UC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")
for(i in 1:16){
  try <- FindMarkers(data, subset.ident = celltypes_UC_and_PSCNI_vs_HC[i], group.by = "state", test.use = "MAST", ident.1 = "UC-NI", ident.2 = "HC-NI")
  try$celltype <- celltypes_UC_and_PSCNI_vs_HC[i]
  #try1 <- try[try$p_val_adj < 0.05,]
  try$Gene <- rownames(try)
  try <- merge(try, UC_risk, by = "Gene")
  UC <- rbind(UC, try)
}

UC_NI<-UC
PSC_NI<-PSC
write.csv(UC, "/DE/Results/UC_sign_risk_genes_16celltypes.csv")
write.csv(PSC, "/DE/Results/PSC_sign_risk_genes_16celltypes.csv")

PSC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(PSC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")
for(i in 1:15){
  try <- FindMarkers(data, subset.ident = celltypes_PSC_and_UCNI_vs_HC[i], group.by = "state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "HC-NI")
  try$celltype <- celltypes_PSC_and_UCNI_vs_HC[i]
  #try1 <- try[try$p_val_adj < 0.05,]
  try$Gene <- rownames(try)
  try <- merge(try, PSC_risk, by = "Gene")
  PSC <- rbind(PSC, try)
}

UC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(UC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")
for(i in 1:16){
  try <- FindMarkers(data, subset.ident = celltypes_UC_and_PSCNI_vs_HC[i], group.by = "state", test.use = "MAST", ident.1 = "UC-I", ident.2 = "HC-NI")
  try$celltype <- celltypes_UC_and_PSCNI_vs_HC[i]
  #try1 <- try[try$p_val_adj < 0.05,]
  try$Gene <- rownames(try)
  try <- merge(try, UC_risk, by = "Gene")
  UC <- rbind(UC, try)
}


write.csv(UC, "/DE/Results/UC_I_sign_risk_genes_16celltypes.csv")
write.csv(PSC, "/DE/Results/PSC_I_sign_risk_genes_15celltypes.csv")

UC_I<-UC
PSC_I<-PSC

PSC_significance_level_bonferroni<-0.05/40 #divided by nr of genes in set
UC_significance_level_bonferroni<-0.05/54 #divided by nr of genes in set

UC_NI_sign<-UC_NI[UC_NI$p_val < UC_significance_level_bonferroni,]
PSC_NI_sign<-PSC_NI[PSC_NI$p_val < PSC_significance_level_bonferroni,]

UC_I_sign<-UC_I[UC_I$p_val < UC_significance_level_bonferroni,]
PSC_I_sign<-PSC_I[PSC_I$p_val < PSC_significance_level_bonferroni,]

PSC_I$p_val_adj<-PSC_I$p_val*40
PSC_NI$p_val_adj<-PSC_NI$p_val*40

UC_I$p_val_adj<-UC_I$p_val*54
UC_NI$p_val_adj<-UC_NI$p_val*54

write.csv(UC_I_sign, "/DE/Results/UC_I_sign_bonferroni_risk_genes_celltypes.csv")
write.csv(PSC_I_sign, "/DE/Results/PSC_I_sign_bonferroni_risk_genes_celltypes.csv")

write.csv(UC_NI_sign, "/DE/Results/UC_NI_sign_bonferroni_risk_genes_celltypes.csv")
write.csv(PSC_NI_sign, "/DE/Results/SC_NI_sign_bonferroni_risk_genes_celltypes.csv")

write.csv(UC_I, "/DE/Results/UC_I_sign_bonferroni_risk_genes_celltypes.csv")
write.csv(PSC_I, "/DE/Results/PSC_I_sign_bonferroni_risk_genes_celltypes.csv")

write.csv(UC_NI, "/DE/Results/UC_NI_bonferroni_risk_genes_celltypes.csv")
write.csv(PSC_NI, "/DE/Results/PSC_NI_bonferroni_risk_genes_celltypes.csv")
library(viridis)

ggplot(PSC_NI, aes(x=celltype, y = Gene, color = avg_logFC , size = -log(p_val))) + 
  geom_point() + scale_color_viridis(option = "D") 


UC_I<-read.csv("/DE/Results/UC_I_sign_bonferroni_risk_genes_celltypes.csv")
PSC_I<-read.csv("/DE/Results/PSC_I_sign_bonferroni_risk_genes_celltypes.csv")

UC_NI<-read.csv("/DE/Results/UC_NI_bonferroni_risk_genes_celltypes.csv", )
PSC_NI<-read.csv("/DE/Results/PSC_NI_bonferroni_risk_genes_celltypes.csv")
library(viridis)
library(ggplot2)
ggplot(PSC_NI, aes(x=celltype, y = Gene, color = avg_logFC , size = -log(p_val))) + 
  geom_point()  + scale_color_viridis(option = "D") +theme(title = element_text()) +
  labs(title =  "PSC risk genes in PSC_NI vs HC") # saved PSC_riskgenes_PSC_NI_vs_HC

ggplot(PSC_I, aes(x=celltype, y = Gene, color = avg_logFC , size = -log(p_val))) + 
  geom_point()  + scale_color_viridis(option = "D") +theme(title = element_text()) +
  labs(title =  "PSC risk genes in PSC_I vs HC") # saved PSC_riskgenes_PSC_I_vs_HC


ggplot(UC_NI, aes(x=celltype, y = Gene, color = avg_logFC , size = -log(p_val))) + 
  geom_point()  + scale_color_viridis(option = "D") +theme(title = element_text()) +
  labs(title =  "UC risk genes in UC_NI vs HC") # saved UC_riskgenes_UC_NI_vs_HC

ggplot(UC_I, aes(x=celltype, y = Gene, color = avg_logFC , size = -log(p_val))) + 
  geom_point()  + scale_color_viridis(option = "D") +theme(title = element_text()) +
  labs(title =  "UC risk genes in UC_I vs HC") # saved UC_riskgenes_UC_I_vs_HC

##
PSC_risk<-read.csv("/DE/Results/PSC_genes.csv", sep=";")
colnames(PSC_risk)[1]<-"Gene"

try <- FindMarkers(data, subset.ident = celltypes_PSC_and_UCNI_vs_HC[1], group.by = "state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "HC-NI")
try$celltype <- celltypes_PSC_and_UCNI_vs_HC[1]
try$Gene <- rownames(try)
PSC <- merge(try, PSC_risk, by = "Gene")

for(i in 2:15){
  try <- FindMarkers(data, subset.ident = celltypes_PSC_and_UCNI_vs_HC[i], group.by = "state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "HC-NI")
  try$celltype <- celltypes_PSC_and_UCNI_vs_HC[i]
  try$Gene <- rownames(try)
  try1 <- merge(try, PSC_risk, by = "Gene")
  PSC <- rbind(PSC, try1)
}

PSC_I<-PSC

try <- FindMarkers(data, subset.ident = celltypes_UC_and_PSCNI_vs_HC[1], group.by = "state", test.use = "MAST", ident.1 = "PSC-NI", ident.2 = "HC-NI")
try$celltype <- celltypes_UC_and_PSCNI_vs_HC[1]
try$Gene <- rownames(try)
PSC <- merge(try, PSC_risk, by = "Gene")

for(i in 2:16){
  try <- FindMarkers(data, subset.ident = celltypes_UC_and_PSCNI_vs_HC[i], group.by = "state", test.use = "MAST", ident.1 = "PSC-NI", ident.2 = "HC-NI")
  try$celltype <- celltypes_UC_and_PSCNI_vs_HC[i]
  try$Gene <- rownames(try)
  try1 <- merge(try, PSC_risk, by = "Gene")
  PSC <- rbind(PSC, try1)
}

PSC_NI<-PSC
x<-merge(PSC_risk, allgenes, by="Gene", all=F) # 71 risk genes in set 
write.csv(x, "/DE/Results/allPSCgenes_in_dataset.csv")
PSC_I$p_val_adj<-PSC_I$p_val*71
PSC_NI$p_val_adj<-PSC_NI$p_val*71

write.csv(PSC_I, "/DE/Results/PSC_I_sign_bonferroni_allrisk_genes_celltypes.csv")
write.csv(PSC_NI, "/DE/Results/PSC_NI_sign_bonferroni_allrisk_genes_celltypes.csv")

PSC_I<-PSC_I[PSC_I$p_val_adj <0.05,]
PSC_NI<-PSC_NI[PSC_NI$p_val_adj <0.05,]

x<-unique(PSC_I$Gene) #13 unique genes
x<-unique(PSC_NI$Gene) #6 unique genes

x<-unique(PSC_I$celltype) #13 unique celltypes
x<-unique(PSC_NI$celltype) #6 unique celltypes
x

data<-ScaleData(data)
DoHeatmap(data, features = x)










