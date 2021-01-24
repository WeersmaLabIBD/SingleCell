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
celltypes <- celltypes[-c(2,3,4,6,7,8,9,10,13,14,15,18,19,20,21,22,24,25,26,32,35,39,40,41,42)] #filter out celltypes with too few (<100) cells per group
DefaultAssay(data)<-"RNA"

#####
allgenes <- data.frame(matrix(ncol = 2, nrow = 33538))
colnames(allgenes)[1]<-"Gene"
allgenes$Gene<-rownames(data@assays$RNA)
PSC_risk_in_set<-merge(allgenes, PSC_risk, by="Gene")
UC_risk_in_set<-merge(allgenes, UC_risk, by="Gene")
write.csv(UC_risk_in_set, "~/ucputativeriskgenes_smillie_inset.csv")
write.csv(PSC_risk_in_set, "~/PSCgenes_inset.csv")

PSC_risk<-read.csv("~/PSC_genes.csv", sep=";")
colnames(PSC_risk)[1]<-"Gene"

DEgenes <- FindMarkers(data, subset.ident = celltypes[1], group.by = "state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "HC-NI")
DEgenes$celltype <- celltypes[1]
DEgenes$Gene <- rownames(DEgenes)
PSC <- merge(DEgenes, PSC_risk, by = "Gene")

for(i in 2:17){
  DEgenes <- FindMarkers(data, subset.ident = celltypes[i], group.by = "state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "HC-NI")
  DEgenes$celltype <- celltypes[i]
  DEgenes$Gene <- rownames(DEgenes)
  DEgenes1 <- merge(DEgenes, PSC_risk, by = "Gene")
  PSC <- rbind(PSC, DEgenes1)
}

PSC_I<-PSC

DEgenes <- FindMarkers(data, subset.ident = celltypes[1], group.by = "state", test.use = "MAST", ident.1 = "PSC-NI", ident.2 = "HC-NI")
DEgenes$celltype <- celltypes[1]
DEgenes$Gene <- rownames(DEgenes)
PSC <- merge(DEgenes, PSC_risk, by = "Gene")

for(i in 2:17){
  DEgenes <- FindMarkers(data, subset.ident = celltypes[i], group.by = "state", test.use = "MAST", ident.1 = "PSC-NI", ident.2 = "HC-NI")
  DEgenes$celltype <- celltypes[i]
  DEgenes$Gene <- rownames(DEgenes)
  DEgenes1 <- merge(DEgenes, PSC_risk, by = "Gene")
  PSC <- rbind(PSC, DEgenes1)
}

PSC_NI<-PSC

x<-merge(PSC_risk, allgenes, by="Gene", all=F) # 71 PSC risk genes in set 
PSC_I$p_val_adj<-PSC_I$p_val*71
PSC_NI$p_val_adj<-PSC_NI$p_val*71

write.csv(PSC_I, "~/PSC_I_bonferroni_allrisk_genes_celltypes.csv")
write.csv(PSC_NI, "~/PSC_NI_onferroni_allrisk_genes_celltypes.csv")

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
UC_significance_level_bonferroni<-0.05/54 #54 UCrisk genes in total in set

UC_I$p_val_adj<-UC_I$p_val*54
UC_NI$p_val_adj<-UC_NI$p_val*54

write.csv(UC_I, "~/UC_I_bonferroni_risk_genes_celltypes.csv")
write.csv(UC_NI, "~/UC_NI_bonferroni_risk_genes_celltypes.csv")

# plot
library(ggplot2)
PSC_I <- filter(PSC_I, PSC_I$p_val_adj <= 0.05)
UC_I <- filter(UC_I, UC_I$p_val_adj <= 0.05)
UC_NI <- filter(UC_NI, UC_NI$p_val_adj <=0.05)
PSC_NI <- filter(PSC_NI, PSC_NI$p_val_adj <=0.05)

UC_I$Adj_p_value <- ifelse(UC_I$p_val_adj <= 0.001, "< 0.001", ifelse(UC_I$p_val_adj <= 0.01, "< 0.01", "< 0.05"))
UC_I$Adj_p_value <- factor(UC_I$Adj_p_value, levels = c("< 0.05", "< 0.01","< 0.001"))
UC_NI$Adj_p_value <- ifelse(UC_NI$p_val_adj <= 0.001, "< 0.001", ifelse(UC_NI$p_val_adj <= 0.01, "< 0.01", "< 0.05"))
UC_NI$Adj_p_value <- factor(UC_NI$Adj_p_value, levels = c("< 0.05", "< 0.01","< 0.001"))
PSC_I$Adj_p_value <- ifelse(PSC_I$p_val_adj <= 0.001, "< 0.001", ifelse(PSC_I$p_val_adj <= 0.01, "< 0.01", "< 0.05"))
PSC_I$Adj_p_value <- factor(PSC_I$Adj_p_value, levels = c("< 0.05", "< 0.01","< 0.001"))
PSC_NI$Adj_p_value <- ifelse(PSC_NI$p_val_adj <= 0.001, "< 0.001", ifelse(PSC_NI$p_val_adj <= 0.01, "< 0.01", "< 0.05"))
PSC_NI$Adj_p_value <- factor(PSC_NI$Adj_p_value, levels = c("< 0.05", "< 0.01","< 0.001"))
UC_I$direction <- ifelse(UC_I$avg_logFC > 0, "up", "down")
UC_NI$direction <- ifelse(UC_NI$avg_logFC > 0, "up", "down")
PSC_I$direction <- ifelse(PSC_I$avg_logFC > 0, "up", "down")
PSC_NI$direction <- ifelse(PSC_NI$avg_logFC > 0, "up", "down")

ggplot(PSC_I, aes(x = Gene, y = celltype, color = direction)) + 
  geom_point(aes(alpha = Adj_p_value), size = 5) +
  scale_alpha_manual(values =c(0.2, 0.5, 1)) +
  scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
  theme_bw() +
  labs(title =  "PSC risk genes in PSC_I vs HC")

ggsave("Results/Figures/Riskgenes_PSCI.pdf", height = 7)

ggplot(PSC_NI, aes(x = Gene, y = celltype, color = direction)) + 
  geom_point(aes(alpha = Adj_p_value), size = 5) +
  scale_alpha_manual(values =c(1)) +
  scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
  theme_bw() +
  labs(title =  "PSC risk genes in PSC_NI vs HC")

ggsave("Results/Figures/Riskgenes_PSCNI.pdf", height = 7)

ggplot(UC_I, aes(x = Gene, y = celltype, color = direction)) + 
  geom_point(aes(alpha = Adj_p_value), size = 5) +
  scale_alpha_manual(values =c(0.2, 0.5, 1)) +
  scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
  theme_bw()+
  labs(title =  "UC risk genes in UC_I vs HC")

ggsave("Results/Figures/Riskgenes_UCI.pdf", height = 7)

ggplot(UC_NI, aes(x = Gene, y = celltype, color = direction)) + 
  geom_point(aes(alpha = Adj_p_value), size = 5) +
  scale_alpha_manual(values =c(0.5, 1)) +
  scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
  theme_bw()+
  labs(title =  "UC risk genes in UC_NI vs HC")

ggsave("Results/Figures/Riskgenes_UCNI.pdf", height = 7)

# old
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

