#analysis UC and PSC risk/implicated genes and their differential expression in different health statuses

data<- readRDS("~/PSC_processed.rds")
DefaultAssay(data) = "RNA"
UC_risk<-read.csv("~/ucputativeriskgenes.csv", sep=";")
genes<-read.csv("~/genes_junctions.csv", sep=";")
PSC_risk<-genes[genes$function. == "PSC_and_suggestive_risk_genes",]

#identifyng risk gene expression
Idents(data) <- "celltypes"
celltypes <- levels(data$celltypes)
table(data$celltypes, data$state)
celltypes <- celltypes[-c(2,3,4,6,7,8,9,10,13,14,15,18,19,20,21,22,24,25,26,32,35,39,40,41,42)] #filter out celltypes with too little (<100) cells per group
DefaultAssay(data)<-"RNA"

#####
allgenes <- data.frame(matrix(ncol = 2, nrow = 33538))
colnames(allgenes)[1]<-"Gene"
allgenes$Gene<-rownames(data@assays$RNA)
PSC_risk_in_set<-merge(allgenes, PSC_risk, by="Gene")
UC_risk_in_set<-merge(allgenes, UC_risk, by="Gene")
write.csv(UC_risk_in_set, "~/ucputativeriskgenes_smillie_inset.csv")
write.csv(PSC_risk_in_set, "~/PSCgenes_inset.csv")

PSC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(PSC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")

for(i in 1:17){
  DEgenes <- FindMarkers(data, subset.ident = celltypes[i], group.by = "state", test.use = "MAST", ident.1 = "PSC-NI", ident.2 = "HC-NI")
  DEgenes$celltype <- celltypes[i]
  DEgenes$Gene <- rownames(DEgenes)
  DEgenes <- merge(DEgenes, PSC_risk, by = "Gene")
  PSC <- rbind(PSC, DEgenes)
}

UC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(UC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")
for(i in 1:17){
  DEgenes <- FindMarkers(data, subset.ident = celltypes[i], group.by = "state", test.use = "MAST", ident.1 = "UC-NI", ident.2 = "HC-NI")
  DEgenes$celltype <- celltypes[i]
  DEgenes$Gene <- rownames(DEgenes)
  DEgenes <- merge(DEgenes, UC_risk, by = "Gene")
  UC <- rbind(UC, DEgenes)
}

UC_NI<-UC
PSC_NI<-PSC

PSC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(PSC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")
for(i in 1:17){
  DEgenes <- FindMarkers(data, subset.ident = celltypes[i], group.by = "state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "HC-NI")
  DEgenes$celltype <- celltypes[i]
  DEgenes$Gene <- rownames(DEgenes)
  DEgenes <- merge(DEgenes, PSC_risk, by = "Gene")
  PSC <- rbind(PSC, DEgenes)
}

UC <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(UC) <- c("p_val", "avg_logFC","pct.1","pct2","p_val_adj","celltype")
for(i in 1:17){
  DEgenes <- FindMarkers(data, subset.ident = celltypes[i], group.by = "state", test.use = "MAST", ident.1 = "UC-I", ident.2 = "HC-NI")
  DEgenes$celltype <- celltypes[i]
  DEgenes$Gene <- rownames(DEgenes)
  DEgenes <- merge(DEgenes, UC_risk, by = "Gene")
  UC <- rbind(UC, DEgenes)
}

UC_I<-UC
PSC_I<-PSC

PSC_significance_level_bonferroni<-0.05/40 #40 genes in total
UC_significance_level_bonferroni<-0.05/54 #54 genes in total

PSC_I$p_val_adj<-PSC_I$p_val*40
PSC_NI$p_val_adj<-PSC_NI$p_val*40

UC_I$p_val_adj<-UC_I$p_val*54
UC_NI$p_val_adj<-UC_NI$p_val*54

write.csv(UC_I, "~/UC_I_bonferroni_risk_genes_celltypes.csv")
write.csv(PSC_I, "~/PSC_I_bonferroni_risk_genes_celltypes.csv")

write.csv(UC_NI, "~/UC_NI_bonferroni_risk_genes_celltypes.csv")
write.csv(PSC_NI, "~/PSC_NI_bonferroni_risk_genes_celltypes.csv")

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

