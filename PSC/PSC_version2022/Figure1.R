library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyr)
library(readxl)
library(openxlsx)
library(reshape2)

# load in dataset
epi <- readRDS("Nieuw/epi_azimuth_duox2.rds")
imm <- readRDS("Nieuw/imm_azimuth_with_plasma.rds")
str <- readRDS("Nieuw/stro_azimuth.rds")

DefaultAssay(epi) = "RNA"
DefaultAssay(imm) = "RNA"
DefaultAssay(str) = "RNA"
Idents(epi) <- "celltype.final"
Idents(imm) <- "predicted.cell_type.pred"
Idents(str) <- "predicted.cell_type.pred"

# make full umap
  # --> ask Paola


# make umaps per subtype 

DimPlot(str, pt.size = 0,   cols = c("#FF5733", "#D95F0E", "#B38CC6",
                                    "#8c564b", "#FADADD", "#7f7f7f", "#d62728", "#FF6F61",
                                    "#aec7e8", "#ffbb78","#C7A995", "#98df8a", "#ff7f0e")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) + custom_theme +
  ggtitle("Stromal cells") +   xlab("UMAP 1") +
  ylab("UMAP 2")

ggsave("Results/Figures_new/Dimplot_str.pdf", width = 8, height = 6)

DimPlot(str, pt.size = 0,   cols = c("#FF5733", "#D95F0E", "#B38CC6",
                                     "#8c564b", "#FADADD", "#7f7f7f", "#d62728", "#FF6F61",
                                     "#aec7e8", "#ffbb78","#C7A995", "#98df8a", "#ff7f0e")) + 
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) + NoLegend() + theme(axis.title.x = element_blank(),
                                                                                         axis.title.y = element_blank(),
                                                                                         axis.text.x = element_blank(),
                                                                                         axis.text.y = element_blank(),
                                                                                         axis.ticks = element_blank())
ggsave("Results/Figures_new/Dimplot_str_notext.pdf", width = 6, height = 6)

DimPlot(imm, pt.size = 0,   cols = c("#1f77b4", "#aec7e8", "#800080","#6baed6", "#08519c", "#555599",
                                    "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
                                    "#8c564b", "#c49c94", "#e377c2","#5E00EB", "#f7b6d2", "#7f7f7f",
                                    "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5",
                                    "#FDE1E2", "#ffbb78", "#ff7f0e", "#2ca02c")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) + custom_theme +
  ggtitle("Immune cells") +   xlab("UMAP 1") +
  ylab("UMAP 2")

ggsave("Results/Figures_new/Dimplot_imm.pdf", width = 8, height = 6)

DimPlot(imm, pt.size = 0,   cols = c("#1f77b4", "#aec7e8", "#800080","#6baed6", "#08519c", "#555599",
                                     "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
                                     "#8c564b", "#c49c94", "#e377c2","#5E00EB", "#f7b6d2", "#7f7f7f",
                                     "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5",
                                     "#FDE1E2", "#ffbb78", "#ff7f0e", "#2ca02c")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) + NoLegend() + theme(axis.title.x = element_blank(),
                                                                                         axis.title.y = element_blank(),
                                                                                         axis.text.x = element_blank(),
                                                                                         axis.text.y = element_blank(),
                                                                                         axis.ticks = element_blank())

ggsave("Results/Figures_new/Dimplot_imm_notext.pdf", width = 6, height = 6)

DimPlot(epi, pt.size = 0,   cols = c("#2D6100", "#d9ef8b", "#ffffbf",
                                     "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#E8E23A", "#66bd63", "#1f9915",
                                     "#1f77b4", "#bcbd22", "#b2df8a", "#8c564b", "#5E00EB", "#9e0142")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) + custom_theme +
  ggtitle("Immune cells") +   xlab("UMAP 1") +
  ylab("UMAP 2")

ggsave("Results/Figures_new/Dimplot_epi.pdf", width = 8, height = 6)



DimPlot(epi, pt.size = 0,   cols = c("#2D6100", "#d9ef8b", "#ffffbf",
                                     "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#E8E23A", "#66bd63", "#1f9915",
                                     "#1f77b4", "#bcbd22", "#b2df8a", "#8c564b", "#5E00EB", "#9e0142")) + NoLegend() + theme(axis.title.x = element_blank(),
                                                                                                                             axis.title.y = element_blank(),
                                                                                                                             axis.text.x = element_blank(),
                                                                                                                             axis.text.y = element_blank(),
                                                                                                                             axis.ticks = element_blank())

ggsave("Results/Figures_new/Dimplot_epi_notext.pdf", width = 6, height = 6)

#make frequency plots
counts=as.data.frame.matrix(table(imm$Final_HTO, imm$predicted.cell_type.pred))
total_row = apply(counts, 1, sum)
pcts = lapply(counts, function(x) {
  x / total_row
})
frequencies = as.data.frame(pcts)
rownames(counts) = rownames(frequencies)
frequencies$disease =  sapply(strsplit(rownames(frequencies),"-"), `[`, 1)
frequencies$inflammation =  sapply(strsplit(rownames(frequencies),"-"), `[`, 2)
frequencies$state =  apply( frequencies[ ,c(26,27)] , 1 , paste , collapse = "-" )
frequencies = frequencies[,-c(26,27)]
df = melt(frequencies)
df$state = factor(df$state, levels = c("HC-NI", "UC-NI", "UC-I", "PSC-NI", "PSC-I"))

df1 <- df[df$variable == "Inflammatory.Monocytes",]
custom_colors <- c("orange", "pink", "purple", "lightgreen", "darkgreen")
ggplot(df1, aes(x=state, y=value, fill=state)) + 
  geom_violin(alpha = 0.2) +
  geom_boxplot(width=0.2, outlier.size = 0.5) +
  #facet_wrap( ~ variable, scales="free") +
  theme_bw(base_size = 14) + 
  scale_fill_manual(values = custom_colors) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), aspect.ratio = 1) + 
  ggtitle("Inflammatory Monocytes") +
  xlab("") +
  ylab("proportion")
ggsave("Results/Figures_new/inflammono.pdf", width = 6)

# figure 2
DimPlot(epi, label = T, repel = T) +NoLegend()
FeaturePlot(epi, features = "DUOX2")
ggsave("Results/Figures_new/Featureplot_DUOX2.pdf", width = 8, height = 6)

# make umaps colored by inflammation
Idents(imm) <- "inflammation"
DimPlot(imm, pt.size = 0,   cols = c("#9467bd", "#17becf")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) + custom_theme
ggsave("Results/Figures_new/Dimplot_infl_imm.pdf", width = 8, height = 6)

DimPlot(imm, pt.size = 0,   cols = c("#9467bd", "#17becf")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1))  +NoLegend()+ theme(axis.title.x = element_blank(),
                                                                             axis.title.y = element_blank(),
                                                                             axis.text.x = element_blank(),
                                                                             axis.text.y = element_blank(),
                                                                             axis.ticks = element_blank()) 
ggsave("Results/Figures_new/Dimplot_infl_imm_notext.pdf", width = 6, height = 6)

Idents(epi) <- "inflammation"
DimPlot(epi, pt.size = 0,   cols = c("#1f9915", "#E8E23A")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1))  + NoLegend() + theme(axis.title.x = element_blank(),
                                                                                          axis.title.y = element_blank(),
                                                                                          axis.text.x = element_blank(),
                                                                                          axis.text.y = element_blank(),
                                                                                          axis.ticks = element_blank()) 
                                                                                           
ggsave("Results/Figures_new/Dimplot_infl_epi_notext.pdf", width = 6, height = 6)

Idents(str) <- "inflammation"
DimPlot(str, pt.size = 0,   cols = c("#ffbb78", "#d62728")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) + custom_theme
ggsave("Results/Figures_new/Dimplot_infl_str.pdf", width = 8, height = 6)

DimPlot(str, pt.size = 0,   cols = c("#ffbb78", "#d62728")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1))+ NoLegend() + theme(axis.title.x = element_blank(),
                                                                           axis.title.y = element_blank(),
                                                                           axis.text.x = element_blank(),
                                                                           axis.text.y = element_blank(),
                                                                           axis.ticks = element_blank()) 
ggsave("Results/Figures_new/Dimplot_infl_str_notext.pdf", width = 6, height = 6)

# stacked barplots
imm@meta.data$groups[imm@meta.data$predicted.cell_type.pred == c("IgA", "IgG", "IgM", "Ig_negative")] <- "Plasma cells"
imm@meta.data$groups[imm@meta.data$predicted.cell_type.pred == c("Follicular", "GC", "Cycling B")] <- "B cells"
imm@meta.data$groups[imm@meta.data$predicted.cell_type.pred == c("CD4+ Activated Fos-hi", "CD4+ Activated Fos-lo", "CD4+ Memory", "CD4+PD1+", "Tregs" )] <- "CD4+ T cells"
imm@meta.data$groups[imm@meta.data$predicted.cell_type.pred == c("CD8+ IELs", "CD8+ IL17+", "CD8+ LP" )] <- "CD8+ T cells"
imm@meta.data$groups[imm@meta.data$predicted.cell_type.pred == c("Cycling T")] <- "Cycling T"
imm@meta.data$groups[imm@meta.data$predicted.cell_type.pred == c("NKs")] <- "NK cells"
imm@meta.data$groups[imm@meta.data$predicted.cell_type.pred == c("ILCs")] <- "ILCs"
imm@meta.data$groups[imm@meta.data$predicted.cell_type.pred == c("CD69- Mast", "CD69+ Mast","Cycling Monocytes","DC1","DC2","Inflammatory Monocytes","Macrophages" )] <- "Myeloid cells"
table(imm$groups)
Idents(imm) <- "groups"

str@meta.data$groups[str@meta.data$predicted.cell_type.pred == c("Endothelial", "Microvascular","Pericytes","Post-capillary Venules")] <- "Endothelial cells"
str@meta.data$groups[str@meta.data$predicted.cell_type.pred == c("Myofibroblasts","RSPO3+","WNT2B+ Fos-hi","WNT2B Fos-lo 1", "WNT2B Fos-lo 2","WNT5B+ 1","WNT5B+ 2")] <- "Fibroblasts"
str@meta.data$groups[str@meta.data$predicted.cell_type.pred == c("Inflammatory Fibroblasts")] <- "Inflammatory Fibroblasts"
str@meta.data$groups[str@meta.data$predicted.cell_type.pred == c("Glia")] <- "Glial cells"
table(str$groups)
Idents(str) <- "groups"

epi@meta.data$groups[epi@meta.data$celltype.final == c("Best4+ Enterocytes", "Enterocytes")] <- "Enterocytes"
epi@meta.data$groups[epi@meta.data$celltype.final == c("DUOX2 enterocytes")] <- "DUOX2+ enterocytes"
epi@meta.data$groups[epi@meta.data$celltype.final == c("Cycling TA", "Secretory TA", "Stem", "TA 1", "TA 2")] <- "Undifferentiated"
epi@meta.data$groups[epi@meta.data$celltype.final == c("Enterocyte Progenitors", "Immature Enterocytes 1", "Immature Enterocytes 2")] <- "Immature absorptive"
epi@meta.data$groups[epi@meta.data$celltype.final == c("Enteroendocrine", "Goblet", "Immature Goblet", "M cells", "Tuft")] <- "Secretory"
table(epi$groups)
Idents(epi) <- "groups"

# create df with state - group - value, with proportion of the group from state total
epi_counts=as.data.frame.matrix(table(epi$state, epi$groups))

total_row = apply(epi_counts, 1, sum)
pcts = lapply(epi_counts, function(x) {
  x / total_row
})
epi_frequencies = as.data.frame(pcts)
rownames(epi_frequencies) = rownames(epi_counts)

epi_frequencies$state <- rownames(epi_frequencies)

epi_frequencies <- epi_frequencies %>%
  gather(group, value, -state) %>%
  arrange(state, group)

epi_frequencies$group <- as.factor(epi_frequencies$group)
epi_frequencies$group <- factor(epi_frequencies$group, levels = c("Enterocytes","DUOX2..enterocytes", "Immature.absorptive", "Secretory","Undifferentiated"))
epi_frequencies$state <- as.factor(epi_frequencies$state)
levels(epi_frequencies$state)
epi_frequencies$state <- factor(epi_frequencies$state, levels = c("HC-NI", "UC-NI", "UC-I", "PSC-NI", "PSC-I"))

data=epi_frequencies
ggplot(data, aes(fill=group, y=value, x=state)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#2D6100",  "#E8B094", "#d9ef8b","#fdae61", "#ffec80"))+
  theme_void()
ggsave("Results/Figures_new/Bar_epi.pdf", width = 6, height = 6)

ggplot(data, aes(fill=group, y=value, x=state)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#2D6100",  "#E8B094", "#d9ef8b","#fdae61", "#ffec80"))+
  theme_void() +NoLegend()
ggsave("Results/Figures_new/Bar_epi_notext.pdf", width = 4, height = 6)


# do imm
# create df with state - group - value, with proportion of the group from state total
imm_counts=as.data.frame.matrix(table(imm$state, imm$groups))

total_row = apply(imm_counts, 1, sum)
pcts = lapply(imm_counts, function(x) {
  x / total_row
})
imm_frequencies = as.data.frame(pcts)
rownames(imm_frequencies) = rownames(imm_counts)

imm_frequencies$state <- rownames(imm_frequencies)

imm_frequencies <- imm_frequencies %>%
  gather(group, value, -state) %>%
  arrange(state, group)

imm_frequencies$group <- as.factor(imm_frequencies$group)
levels(imm_frequencies$group)
imm_frequencies$group <- factor(imm_frequencies$group, levels = c("Myeloid.cells", "Plasma.cells", "B.cells", "CD8..T.cells", "CD4..T.cells", "Cycling.T", "NK.cells", "ILCs"))

imm_frequencies$state <- as.factor(imm_frequencies$state)
levels(imm_frequencies$state)
imm_frequencies$state <- factor(imm_frequencies$state, levels = c("HC-NI", "UC-NI", "UC-I", "PSC-NI", "PSC-I"))

data=imm_frequencies
ggplot(data, aes(fill=group, y=value, x=state)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#9edae5","#08519c","#17becf","#1f77b4" ,"#c7c7c7","#800080","#ff9896","#d62728"))+
  theme_void()
ggsave("Results/Figures_new/Bar_imm.pdf", width = 6, height = 6)

ggplot(data, aes(fill=group, y=value, x=state)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#9edae5","#08519c","#17becf","#1f77b4" ,"#c7c7c7","#800080","#ff9896","#d62728"))+
  theme_void() + NoLegend()
ggsave("Results/Figures_new/Bar_imm_notext.pdf", width = 4, height = 6)

# do str
# create df with state - group - value, with proportion of the group from state total
str_counts=as.data.frame.matrix(table(str$state, str$groups))

total_row = apply(str_counts, 1, sum)
pcts = lapply(str_counts, function(x) {
  x / total_row
})
str_frequencies = as.data.frame(pcts)
rownames(str_frequencies) = rownames(str_counts)

str_frequencies$state <- rownames(str_frequencies)

str_frequencies <- str_frequencies %>%
  gather(group, value, -state) %>%
  arrange(state, group)

str_frequencies$group <- as.factor(str_frequencies$group)
levels(str_frequencies$group)
str_frequencies$group <- factor(str_frequencies$group, levels = c("Inflammatory.Fibroblasts","Fibroblasts","Glial.cells","Endothelial.cells"))
str_frequencies$state <- as.factor(str_frequencies$state)
levels(str_frequencies$state)
str_frequencies$state <- factor(str_frequencies$state, levels = c("HC-NI", "UC-NI", "UC-I", "PSC-NI", "PSC-I"))

data=str_frequencies
ggplot(data, aes(fill=group, y=value, x=state)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#c5b0d5","#d62728", "#ff7f0e","#ffbb78"))+
  theme_void()
ggsave("Results/Figures_new/Bar_str.pdf", width = 6, height = 6)

ggplot(data, aes(fill=group, y=value, x=state)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#c5b0d5","#d62728", "#ff7f0e","#ffbb78"))+
  theme_void() +NoLegend()
ggsave("Results/Figures_new/Bar_str_notext.pdf", width = 4, height = 6)
