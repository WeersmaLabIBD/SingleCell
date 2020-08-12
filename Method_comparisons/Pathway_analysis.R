install.packages("ggpubr")
library(ggpubr)
library(ggplot2)

methodspaper@meta.data$alldatasets<-methodspaper@meta.data$Sample
methodspaper@meta.data$alldatasets[methodspaper@meta.data$dataset == "HC_2"]<-"HC_2"
methodspaper@meta.data$alldatasets[methodspaper@meta.data$dataset == "HC_6"]<-"HC_6"
methodspaper@meta.data$alldatasets[methodspaper@meta.data$dataset == "HC_9"]<-"HC_9"
methodspaper@meta.data$alldatasets[methodspaper@meta.data$dataset == "HC_3002"]<-"HC_3002"
methodspaper@meta.data$alldatasets[methodspaper@meta.data$dataset == "HC_3030"]<-"HC_3030"
methodspaper@meta.data$alldatasets[methodspaper@meta.data$dataset == "HC_3034"]<-"HC_3034"
methodspaper@meta.data$alldatasets[methodspaper@meta.data$dataset == "HC_3037"]<-"HC_3037"
methodspaper@meta.data$alldatasets[methodspaper@meta.data$dataset == "HC_3049"]<-"HC_3049"
methodspaper@meta.data$alldatasets[methodspaper@meta.data$dataset == "HC_3083"]<-"HC_3083"
methodspaper@meta.data$alldatasets[methodspaper@meta.data$dataset == "HC_3296"]<-"HC_3296"

methodspaper@meta.data$allsamples<-methodspaper@meta.data$alldatasets
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N10.EpiA"]<-"N10.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N10.LPA"]<-"N10.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N10.EpiB"]<-"N10.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N10.LPB"]<-"N10.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N11.EpiA"]<-"N11.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N11.LPA"]<-"N11.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N11.EpiB"]<-"N11.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N11.LPB"]<-"N11.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N13.EpiA"]<-"N13.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N13.LPA"]<-"N13.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N13.EpiB"]<-"N13.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N13.LPB"]<-"N13.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N15.EpiA"]<-"N15.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N15.LPA"]<-"N15.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N15.EpiB"]<-"N15.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N15.LPB"]<-"N15.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N16.EpiA"]<-"N16.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N16.LPA"]<-"N16.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N16.EpiB"]<-"N16.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N16.LPB"]<-"N16.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N17.EpiA"]<-"N17.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N17.LPA"]<-"N17.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N17.EpiB"]<-"N17.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N17.LPB"]<-"N17.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N18.EpiA"]<-"N18.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N18.LPA"]<-"N18.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N18.EpiB"]<-"N18.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N18.LPB"]<-"N18.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N20.EpiA"]<-"N20.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N20.LPA"]<-"N20.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N20.EpiB"]<-"N20.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N20.LPB"]<-"N20.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N21.EpiA"]<-"N21.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N21.LPA"]<-"N21.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N21.EpiB"]<-"N21.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N21.LPB"]<-"N21.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N46.EpiA"]<-"N46.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N46.LPA"]<-"N46.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N46.EpiB"]<-"N46.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N46.LPB"]<-"N46.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N51.EpiA"]<-"N51.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N51.LPA"]<-"N51.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N51.EpiB"]<-"N51.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N51.LPB"]<-"N51.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N8.EpiA"]<-"N8.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N8.LPA"]<-"N8.A"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N8.EpiB"]<-"N8.B"
methodspaper@meta.data$allsamples[methodspaper@meta.data$alldatasets == "N8.LPB"]<-"N8.B"

table(methodspaper@meta.data$allsamples)

dataframe_proportions_sample_celltype_method<-data.frame(prop.table(table(methodspaper@meta.data$Predicted_all, methodspaper@meta.data$allsamples, methodspaper@meta.data$method ), 2)*100)

dataframe_proportions_layers_celltype_method<-data.frame(prop.table(table(methodspaper@meta.data$Predicted_all, methodspaper@meta.data$alldatasets, methodspaper@meta.data$method ), 2)*100)


my_comparisions=list( c("splitcollagenase","splitprotease"), c("splitcollagenase","wholecollagenase"), c("splitprotease","wholecollagenase"))

ggplot (dataframe_proportions_sample_celltype_method, aes (Var3,Freq,fill=Var3)) + geom_violin() + geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Frequency") + xlab ("Protocol")

ggplot (dataframe_proportions_layers_celltype_method, aes (Var3,Freq,fill=Var3)) + geom_violin() + geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Frequency") + xlab ("Protocol")

dataframe_proportions_layers_celltype_method$relative_proportion[dataframe_proportions_layers_celltype_method$Var2 == "splitprotease"]<-dataframe_proportions_layers_celltype_method$Freq[,]/3

#sum(dataframe_proportions_layers_celltype_method$Freq[dataframe_proportions_layers_celltype_method$Var2 == "splitprotease"])

dataframe_proportions_layers_celltype_method[dataframe_proportions_layers_celltype_method$Var3=="splitprotease",]$relative_proportion=(dataframe_proportions_layers_celltype_method$Freq[dataframe_proportions_layers_celltype_method$Var3 == "splitprotease"])/sum(dataframe_proportions_layers_celltype_method$Freq[dataframe_proportions_layers_celltype_method$Var3 == "splitprotease"])
dataframe_proportions_layers_celltype_method[dataframe_proportions_layers_celltype_method$Var3=="splitcollagenase",]$relative_proportion=(dataframe_proportions_layers_celltype_method$Freq[dataframe_proportions_layers_celltype_method$Var3 == "splitcollagenase"])/sum(dataframe_proportions_layers_celltype_method$Freq[dataframe_proportions_layers_celltype_method$Var3 == "splitcollagenase"])
dataframe_proportions_layers_celltype_method[dataframe_proportions_layers_celltype_method$Var3=="wholecollagenase",]$relative_proportion=(dataframe_proportions_layers_celltype_method$Freq[dataframe_proportions_layers_celltype_method$Var3 == "wholecollagenase"])/sum(dataframe_proportions_layers_celltype_method$Freq[dataframe_proportions_layers_celltype_method$Var3 == "wholecollagenase"])

ggplot(dataframe_proportions_layers_celltype_method, aes(Var3,relative_proportion, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))

ggplot(dataframe_proportions_layers_celltype_method, aes(Var3,relative_proportion, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))


whole_coll<-dataframe_proportions_layers_celltype_method[dataframe_proportions_layers_celltype_method$Var3 == "wholecollagenase",]
split_coll<-dataframe_proportions_layers_celltype_method[dataframe_proportions_layers_celltype_method$Var3 == "splitcollagenase",]
split_prot<-dataframe_proportions_layers_celltype_method[dataframe_proportions_layers_celltype_method$Var3 == "splitprotease",]
whole_coll<-whole_coll[whole_coll$Freq != 0,]
split_coll<-split_coll[split_coll$Freq != 0,]

split_coll_epi<-split_coll[grepl(split_coll$Var2, pattern = "Epi"),]
split_coll_imm<-split_coll[grepl(split_coll$Var2, pattern = "LP"),]

split_prot<-split_prot[split_prot$Freq != 0,]
ggplot(whole_coll, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))

ggplot(split_coll, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))

ggplot(split_coll_epi, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() + scale_fill_manual(values=c("firebrick2","gray88","gold","deepskyblue2", "hotpink2"))

ggplot(split_coll_imm, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))


ggplot(split_prot, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","purple4","deepskyblue2", "hotpink2"))

datasets_finecelltypes<-data.frame(prop.table(table(methodspaper@meta.data$celltype,methodspaper@meta.data$alldatasets, methodspaper@meta.data$method), 2))
samples_finecelltypes<-data.frame(prop.table(table(methodspaper@meta.data$celltype,methodspaper@meta.data$allsamples, methodspaper@meta.data$method), 2))


finecelltype_whole_coll<-datasets_finecelltypes[datasets_finecelltypes$Var3 == "wholecollagenase",]
finecelltype_split_coll<-datasets_finecelltypes[datasets_finecelltypes$Var3 == "splitcollagenase",]
finecelltype_split_prot<-datasets_finecelltypes[datasets_finecelltypes$Var3 == "splitprotease",]
finecelltype_whole_coll<-finecelltype_whole_coll[finecelltype_whole_coll$Freq != 0,]
finecelltype_split_coll<-finecelltype_split_coll[finecelltype_split_coll$Freq != 0,]
finecelltype_split_prot<-finecelltype_split_prot[finecelltype_split_prot$Freq != 0,]

finecelltype_split_coll_epi<-finecelltype_split_coll[grepl(finecelltype_split_coll$Var2, pattern = "Epi"),]
finecelltype_split_coll_imm<-finecelltype_split_coll[grepl(finecelltype_split_coll$Var2, pattern = "LP"),]



ggplot(finecelltype_whole_coll, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() 

ggplot(finecelltype_split_coll_epi, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() 

ggplot(finecelltype_split_coll_imm, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() 

ggplot(finecelltype_split_prot, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() 

table(methodspaper@meta.data$celltype, methodspaper@meta.data$Location)

methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "0"]<-"IgAplasma1"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "1"]<-"T"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "2"]<-"immature_progenitor_enterocyte" 
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "3"]<-"IgGplasma"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "4"]<-"follicularB"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "5"]<-"mitoHi"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "6"]<-"cycling_TA"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "7"]<-"goblet"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "8"]<-"absorptive_enterocyte"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "9"]<-"stem"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "10"]<-"riboHi"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "11"]<-"MAST/fibroblasts"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "12"]<-"IgAplasma2"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "13"]<-"APC"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "14"]<-"Tuft"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "15"]<-"endothelium"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "16"]<-"BEST4_enterocyte"
methodspaper2@meta.data$celltype[methodspaper2@meta.data$seurat_clusters == "17"]<-"GC"
proportion_celltypes2<-data.frame(prop.table(table(methodspaper2@meta.data$celltype, methodspaper2@meta.data$method), 2))
round(prop.table(table(methodspaper@meta.data$celltype, methodspaper@meta.data$method), 2)*100)
markers_methodspaper2<-FindAllMarkers(methodspaper2, test.use = "MAST", only.pos = T, min.pct = 0.25)
cluster0<-subset(methodspaper2, subset = seurat_clusters == "0")

DefaultAssay(cluster0)<-"integrated"
cluster0<-FindNeighbors(cluster0, dims=1:30)
cluster0<-FindClusters(cluster0, resolution = 0.2)
markers_cluster0_methodspaper2<-FindAllMarkers(cluster0, test.use = "MAST", only.pos = T, min.pct = 0.25)

cluster0@meta.data$celltype[cluster0@meta.data$seurat_clusters == "0"]<-"epithelium"
cluster0@meta.data$celltype[cluster0@meta.data$seurat_clusters == "1"]<-"B"
cluster0@meta.data$celltype[cluster0@meta.data$seurat_clusters == "2"]<-"mitoHi"


cluster1<-subset(methodspaper2, subset = seurat_clusters == "1")
DefaultAssay(cluster1)<-"integrated"
cluster1<-FindNeighbors(cluster1, dims=1:30)
cluster1<-FindClusters(cluster1, resolution = 0.2)
DimPlot(cluster1)
DefaultAssay(cluster1)<-"RNA"
FeaturePlot(cluster1, c("CD8A", "CD4", "SELL"))
markers_cluster1_methodspaper2<-FindAllMarkers(cluster1, test.use = "MAST", only.pos = T, min.pct = 0.25)

cluster1@meta.data$celltype[cluster1@meta.data$seurat_clusters == "0"]<-"CD4_mem"
cluster1@meta.data$celltype[cluster1@meta.data$seurat_clusters == "1"]<-"CD4_activated"
cluster1@meta.data$celltype[cluster1@meta.data$seurat_clusters == "2"]<-"CD8_T"


CellsMeta<-methodspaper2@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-cluster0@meta.data
y<-cluster1@meta.data
x$NAME<-rownames(x)
y$NAME<-rownames(y)
row.names(CellsMeta)=NULL
row.names(x)<-NULL
row.names(y)<-NULL
x<-x[,c(20,42)]
y<-y[,c(20,42)]
z<-rbind(x,y)
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=z, by = "NAME", all=T)

CellsMeta$celltype.x[CellsMeta$celltype.y == "CD4_mem"]<-"CD4_memory"
colnames(CellsMeta)[42]<-"celltype_subdivided"
CellsMeta<-CellsMeta[c(1,42)]
rownames(CellsMeta)<-CellsMeta$NAME
methodspaper2<-AddMetaData(methodspaper2, CellsMeta)
DimPlot(methodspaper2, group.by="celltype_subdivided")


methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Best4+ Enterocytes"]<-"epithelium"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Enterocytes"]<-"epithelium"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Enterocyte Progenitors"]<-"epithelium"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Immature Enterocytes 1"]<-"epithelium"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Immature Enterocytes 2"]<-"epithelium"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "TA 1"]<-"epithelium"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "TA 2"]<-"epithelium"
                      
                      # epi general
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "M cells"]<-"epithelium"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Cycling TA"]<-"epithelium"
                      
                      # epi secretory
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Goblet"]<-"epithelium"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Enteroendocrine"]<-"epithelium"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Immature Goblet"]<-"epithelium"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Secretory TA"]<-"epithelium"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Stem"]<-"epithelium"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Tuft"]<-"epithelium"
                      
                      
                      #Tcell
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "CD4+ Activated Fos-lo"]<-"Tcell"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "CD4+ Activated Fos-hi"]<-"Tcell"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "CD4+ Memory"]<-"Tcell"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "CD4+ PD1+"]<-"Tcell"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "CD8+ IELs"]<-"Tcell"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "CD8+ IL17+"]<-"Tcell"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "CD8+ LP"]<-"Tcell"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "MT-hi"]<-"Tcell"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Tregs"]<-"Tcell"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Cycling T"]<-"Tcell"
                      
                      
                      # ILS?NK
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "ILCs"]<-"Lymphoid"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "NKs"]<-"Lymphoid"
                      
                      
                      
                      
                      #B
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Cycling B"]<-"Bcell"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "GC"]<-"Bcell"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Plasma"]<-"Bcell"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Follicular"]<-"Bcell"
                      
                      #Myeloid
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "CD69+ Mast"]<-"myeloid"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "CD69- Mast"]<-"myeloid"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Inflammatory Monocytes"]<-"myeloid"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Macrophages"]<-"myeloid"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Cycling Monocytes"]<-"myeloid"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "DC1"]<-"myeloid"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "DC2"]<-"myeloid"
                      
                      #mesenchymal
                      
                      
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Glia"]<-"Glia"
                      
                      
                      
                      # fibroblasts
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Inflammatory Fibroblasts"]<-"fibroblast"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "RSPO3+"]<-"fibroblast"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "WNT2B+ Fos-hi"]<-"fibroblast"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "WNT2B+ Fos-lo 1"]<-"fibroblast"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "WNT2B+ Fos-lo 2"]<-"fibroblast"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "WNT5B+ 1"]<-"fibroblast"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "WNT5B+ 2"]<-"fibroblast"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "WNT2B+ Fos-lo 1"]<-"fibroblast"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Myofibroblasts"]<-"fibroblast"
                      
                      #endothelial
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Post-capillary Venules"]<-"endothelial"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Pericytes"]<-"endothelial"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Endothelial"]<-"endothelial"
                      methodspaper2@meta.data$Predicted_all[methodspaper2@meta.data$Cluster == "Microvascular"]<-"endothelial"
                      



Idents(methodspaper2)<-"Predicted_all"
DefaultAssay(methodspaper2)<-"RNA"
Tcell<-subset(methodspaper2, subset = Predicted_all == "Tcell" )
Idents(Tcell)<-"method"
T_method_markers<-FindAllMarkers(object = Tcell, test.use = "MAST", min.pct = 0.25)
write.csv(T_method_markers, "

/T_method_markers_MAST.csv")
Bcell<-subset(methodspaper2, subset = Predicted_all == "Bcell" )
Idents(Bcell)<-"method"
B_method_markers<-FindAllMarkers(object = Bcell, test.use = "MAST", min.pct = 0.25)
write.csv(B_method_markers, "~/Desktop/methods_paper/methodspaper2/B_method_markers_MAST.csv")

endothelial<-subset(methodspaper2, subset = Predicted_all == "endothelial" )
Idents(endothelial)<-"method"
endothelial_method_markers<-FindAllMarkers(object = endothelial, test.use = "MAST", min.pct = 0.25)
write.csv(endothelial_method_markers, "~/Desktop/methods_paper/methodspaper2/endothelial_method_markers_MAST.csv")

epithelium<-subset(methodspaper2, subset = Predicted_all == "epithelium" )
Idents(epithelium)<-"method"
epithelium_method_markers<-FindAllMarkers(object = epithelium, test.use = "MAST", min.pct = 0.25)
write.csv(epithelium_method_markers, "~/Desktop/methods_paper/methodspaper2/epithelium_method_markers_MAST.csv")

fibroblast<-subset(methodspaper2, subset = Predicted_all == "fibroblast" )
Idents(fibroblast)<-"method"
fibroblast_method_markers<-FindAllMarkers(object = fibroblast, test.use = "MAST", min.pct = 0.25)
write.csv(fibroblast_method_markers, "~/Desktop/methods_paper/methodspaper2/fibroblast_method_markers_MAST.csv")

Glia<-subset(methodspaper2, subset = Predicted_all == "Glia" )
Idents(Glia)<-"method"
Glia_method_markers<-FindAllMarkers(object = Glia, test.use = "MAST", min.pct = 0.25)
write.csv(Glia_method_markers, "~/Desktop/methods_paper/methodspaper2/Glia_method_markers_MAST.csv")

myeloid<-subset(methodspaper2, subset = Predicted_all == "myeloid" )
Idents(myeloid)<-"method"
myeloid_method_markers<-FindAllMarkers(object = myeloid, test.use = "MAST", min.pct = 0.25)
write.csv(myeloid_method_markers, "~/Desktop/methods_paper/methodspaper2/myeloid_method_markers_MAST.csv")

Lymphoid<-subset(methodspaper2, subset = Predicted_all == "Lymphoid" )
Idents(Lymphoid)<-"method"
Lymphoid_method_markers<-FindAllMarkers(object = Lymphoid, test.use = "MAST", min.pct = 0.25)
write.csv(Lymphoid_method_markers, "~/Desktop/methods_paper/methodspaper2/Lymphoid_method_markers_MAST.csv")

saveRDS(methodspaper2, "~/Desktop/methods_paper/2nd_SCT_integration_7perdc_smillie.rds")

collagenase_genes<-read.csv("~/Desktop/methods_paper/collagenase_genes_from_article_wtc.csv", sep=";")
collagenase_genes<-collagenase_genes[c(2,3)]
colnames(collagenase_genes)[1]<-"gene"
T_cells_collagenase<-merge(collagenase_genes, T_method_markers, all=F)
write.csv(T_cells_collagenase, "~/Desktop/methods_paper/methodspaper2/T_method_collagenase.csv")

fibroblast_collagenase<-merge(collagenase_genes, fibroblast_method_markers, all=F)
write.csv(fibroblast_collagenase, "~/Desktop/methods_paper/methodspaper2/fibroblast_method_collagenase.csv")

epithelium_collagenase<-merge(collagenase_genes, epithelium_method_markers, all=F)
write.csv(epithelium_collagenase, "~/Desktop/methods_paper/methodspaper2/epithelium_method_collagenase.csv")

endothelial_collagenase<-merge(collagenase_genes, endothelial_method_markers, all=F)
write.csv(endothelial_collagenase, "~/Desktop/methods_paper/methodspaper2/endothelial_method_collagenase.csv")

myeloid_collagenase<-merge(collagenase_genes, myeloid_method_markers, all=F)
write.csv(myeloid_collagenase, "~/Desktop/methods_paper/methodspaper2/myeloid_method_collagenase.csv")

Lymphoid_collagenase<-merge(collagenase_genes, Lymphoid_method_markers, all=F)
write.csv(Lymphoid_collagenase, "~/Desktop/methods_paper/methodspaper2/Lymphoid_method_collagenase.csv")

Glia_collagenase<-merge(collagenase_genes, Glia_method_markers, all=F)
write.csv(Glia_collagenase, "~/Desktop/methods_paper/methodspaper2/Glia_method_collagenase.csv")
 
## methodspaper2



methodspaper2@meta.data$alldatasets<-methodspaper2@meta.data$Sample
methodspaper2@meta.data$alldatasets[methodspaper2@meta.data$dataset == "HC_2"]<-"HC_2"
methodspaper2@meta.data$alldatasets[methodspaper2@meta.data$dataset == "HC_6"]<-"HC_6"
methodspaper2@meta.data$alldatasets[methodspaper2@meta.data$dataset == "HC_9"]<-"HC_9"
methodspaper2@meta.data$alldatasets[methodspaper2@meta.data$dataset == "HC_3002"]<-"HC_3002"
methodspaper2@meta.data$alldatasets[methodspaper2@meta.data$dataset == "HC_3030"]<-"HC_3030"
methodspaper2@meta.data$alldatasets[methodspaper2@meta.data$dataset == "HC_3034"]<-"HC_3034"
methodspaper2@meta.data$alldatasets[methodspaper2@meta.data$dataset == "HC_3037"]<-"HC_3037"
methodspaper2@meta.data$alldatasets[methodspaper2@meta.data$dataset == "HC_3049"]<-"HC_3049"
methodspaper2@meta.data$alldatasets[methodspaper2@meta.data$dataset == "HC_3083"]<-"HC_3083"
methodspaper2@meta.data$alldatasets[methodspaper2@meta.data$dataset == "HC_3296"]<-"HC_3296"

methodspaper2@meta.data$allsamples<-methodspaper2@meta.data$alldatasets
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N10.EpiA"]<-"N10.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N10.LPA"]<-"N10.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N10.EpiB"]<-"N10.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N10.LPB"]<-"N10.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N11.EpiA"]<-"N11.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N11.LPA"]<-"N11.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N11.EpiB"]<-"N11.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N11.LPB"]<-"N11.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N13.EpiA"]<-"N13.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N13.LPA"]<-"N13.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N13.EpiB"]<-"N13.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N13.LPB"]<-"N13.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N15.EpiA"]<-"N15.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N15.LPA"]<-"N15.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N15.EpiB"]<-"N15.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N15.LPB"]<-"N15.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N16.EpiA"]<-"N16.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N16.LPA"]<-"N16.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N16.EpiB"]<-"N16.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N16.LPB"]<-"N16.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N17.EpiA"]<-"N17.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N17.LPA"]<-"N17.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N17.EpiB"]<-"N17.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N17.LPB"]<-"N17.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N18.EpiA"]<-"N18.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N18.LPA"]<-"N18.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N18.EpiB"]<-"N18.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N18.LPB"]<-"N18.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N20.EpiA"]<-"N20.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N20.LPA"]<-"N20.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N20.EpiB"]<-"N20.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N20.LPB"]<-"N20.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N21.EpiA"]<-"N21.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N21.LPA"]<-"N21.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N21.EpiB"]<-"N21.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N21.LPB"]<-"N21.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N46.EpiA"]<-"N46.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N46.LPA"]<-"N46.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N46.EpiB"]<-"N46.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N46.LPB"]<-"N46.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N51.EpiA"]<-"N51.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N51.LPA"]<-"N51.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N51.EpiB"]<-"N51.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N51.LPB"]<-"N51.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N8.EpiA"]<-"N8.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N8.LPA"]<-"N8.A"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N8.EpiB"]<-"N8.B"
methodspaper2@meta.data$allsamples[methodspaper2@meta.data$alldatasets == "N8.LPB"]<-"N8.B"

datasets_finecelltypes<-data.frame(prop.tablepro(table(methodspaper2@meta.data$celltype,methodspaper2@meta.data$alldatasets, methodspaper2@meta.data$method), 2))
samples_finecelltypes<-data.frame(prop.table(table(methodspaper2@meta.data$celltype,methodspaper2@meta.data$allsamples, methodspaper2@meta.data$method), 2))


finecelltype_whole_coll<-datasets_finecelltypes[datasets_finecelltypes$Var3 == "wholecollagenase",]
finecelltype_split_coll<-datasets_finecelltypes[datasets_finecelltypes$Var3 == "splitcollagenase",]
finecelltype_split_prot<-datasets_finecelltypes[datasets_finecelltypes$Var3 == "splitprotease",]
finecelltype_whole_coll<-finecelltype_whole_coll[finecelltype_whole_coll$Freq != 0,]
finecelltype_split_coll<-finecelltype_split_coll[finecelltype_split_coll$Freq != 0,]
finecelltype_split_prot<-finecelltype_split_prot[finecelltype_split_prot$Freq != 0,]

finecelltype_split_coll_epi<-finecelltype_split_coll[grepl(finecelltype_split_coll$Var2, pattern = "Epi"),]
finecelltype_split_coll_imm<-finecelltype_split_coll[grepl(finecelltype_split_coll$Var2, pattern = "LP"),]


ggplot(finecelltype_whole_coll, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() 

ggplot(finecelltype_split_coll_epi, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() 

ggplot(finecelltype_split_coll_imm, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() 

ggplot(finecelltype_split_prot, aes(Var2,Freq, fill=Var1 ))+
  geom_bar(position="stack", width = 0.8, stat = "identity", color="black") + theme_classic() 

dataframe_proportions_sample_celltype_method<-data.frame(prop.table(table(methodspaper2@meta.data$Predicted_all, methodspaper2@meta.data$allsamples, methodspaper2@meta.data$method ), 2)*100)

prop_samples_splitcoll<-dataframe_proportions_sample_celltype_method[dataframe_proportions_sample_celltype_method$Var3 == "splitcollagenase",]
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Bcell",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Tcell",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "endothelial",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "epithelium",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "fibroblast",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Glia",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Lymphoid",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "myeloid",]$Freq)


prop_samples_wholecoll<-dataframe_proportions_sample_celltype_method[dataframe_proportions_sample_celltype_method$Var3 == "wholecollagenase",]
prop_samples_wholecoll<-prop_samples_wholecoll[grepl(prop_samples_wholecoll$Var2, pattern = "HC_3"),]

median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Bcell",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Tcell",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "endothelial",]$Freq)
median(prop_samples_wholecoll[prop_samples_wholecoll$Var1 == "epithelium",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "fibroblast",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Glia",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Lymphoid",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "myeloid",]$Freq)

prop_samples_splitprot<-dataframe_proportions_sample_celltype_method[dataframe_proportions_sample_celltype_method$Var3 == "splitprotease",]
prop_samples_splitprot<-prop_samples_splitprot[grepl(prop_samples_splitprot$Var2, pattern = "HC_2|HC_6|HC_9"),]

median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Bcell",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Tcell",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "endothelial",]$Freq)
median(prop_samples_splitprot[prop_samples_splitprot$Var1 == "epithelium",]$Freq)
median(prop_samples_splitprot[prop_samples_splitprot$Var1 == "fibroblast",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Glia",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Lymphoid",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "myeloid",]$Freq)

prop_layer_splitcoll<-dataframe_proportions_layers_celltype_method[dataframe_proportions_layers_celltype_method$Var3 == "splitcollagenase",]


median(prop_layer_splitcoll[prop_layer_splitcoll$Var1 == "Bcell",]$Freq)
median(prop_layer_splitcoll[prop_layer_splitcoll$Var1 == "Tcell",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "endothelial",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "epithelium",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "fibroblast",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Glia",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "Lymphoid",]$Freq)
median(prop_samples_splitcoll[prop_samples_splitcoll$Var1 == "myeloid",]$Freq)

prop_layer_splitcoll_epi<-prop_layer_splitcoll[grepl(prop_layer_splitcoll$Var2, pattern = "Epi"),]
prop_layer_splitcoll_imm<-prop_layer_splitcoll[grepl(prop_layer_splitcoll$Var2, pattern = "LP"),]

x<-rbind(prop_layer_splitcoll_imm, prop_layer_splitcoll_epi)


median(x[x$Var1 == "epithelium",]$Freq)
median(x[x$Var1 == "fibroblast",]$Freq)

median(prop_layer_splitcoll_epi[prop_layer_splitcoll_epi$Var1 == "epithelium",]$Freq)
median(prop_layer_splitcoll_epi[prop_layer_splitcoll_epi$Var1 == "fibroblast",]$Freq)

median(prop_layer_splitcoll_imm[prop_layer_splitcoll_imm$Var1 == "epithelium",]$Freq)
median(prop_layer_splitcoll_imm[prop_layer_splitcoll_imm$Var1 == "fibroblast",]$Freq)

whole_coll_epi<-subset(methodspaper2, Location == "Epi")
whole_coll_LP<-subset(methodspaper2, Location == "LP")

mean(whole_coll_epi@meta.data$percent.mt)
mean(whole_coll_LP@meta.data$percent.mt)
mean(wholecoll@meta.data$percent.mt)
#####
whole_coll_epi<-subset(splitcoll, Location == "Epi")
whole_coll_LP<-subset(splitcoll, Location == "LP")

proportions_split_prot<-data.frame(prop.table(table(splitprot@meta.data$Predicted_all, splitprot@meta.data$dataset),2))
proportions_whole_coll<-data.frame(prop.table(table(wholecoll@meta.data$Predicted_all, wholecoll@meta.data$dataset),2))
proportions_split_coll_epi<-data.frame(prop.table(table(whole_coll_epi@meta.data$Predicted_all, whole_coll_epi@meta.data$Sample),2))
proportions_split_coll_LP<-data.frame(prop.table(table(whole_coll_LP@meta.data$Predicted_all, whole_coll_LP@meta.data$Sample),2))

proportions_split_coll_epi$method<-"split_coll_epi"
proportions_split_coll_LP$method<-"split_coll_LP"
proportions_split_prot$method<-"split_prot"
proportions_whole_coll$method<-"whole_coll"

proportions_all<-rbind(proportions_split_coll_epi, proportions_split_coll_LP)
proportions_all<-rbind(proportions_all, proportions_whole_coll)
proportions_all<-rbind(proportions_all, proportions_split_prot)

ggplot(proportions_whole_coll, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2")) 
# and save Barplot_freqs_allsamples_wholecollagenase 8*10inch as pdf 
ggplot(proportions_split_coll_epi, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","gray88","gold", "deepskyblue2", "hotpink2")) 
# and save Barplot_freqs_allsamples_splicollagenase_epi 8*40inch as pdf 
ggplot(proportions_split_coll_LP, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2")) 
# and save Barplot_freqs_allsamples_splitcollagenase_LP 8*40inch as pdf 
ggplot(proportions_split_prot, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","deepskyblue2", "hotpink2")) 
# and save Barplot_freqs_allsamples_splitprot 8*8inch as pdf 


  
proportions_B<-proportions_all[proportions_all$Var1 == "Bcell",]
proportions_T<-proportions_all[proportions_all$Var1 == "Tcell",]
proportions_epi<-proportions_all[proportions_all$Var1 == "epithelium",]
proportions_endo<-proportions_all[proportions_all$Var1 == "endothelial",]
proportions_myeloid<-proportions_all[proportions_all$Var1 == "myeloid",]
proportions_fib<-proportions_all[proportions_all$Var1 == "fibroblast",]



ggplot(proportions_B, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","deepskyblue2", "hotpink2")) + stat_compare_means(aes(group = Var2), method="kruskal.test", label = "p.format")

p <- ggboxplot(proportions_B, x = "method", y = "Freq",
               color = "method", palette = "jco",
               add = "jitter") + stat_compare_means(aes(group = method), method="wilcox.test", paired = F, label = "p.format")
p + theme(axis.text.x = element_text(angle = 45, hjust = 1))

my_comparisions=list( c("split_coll_epi","split_coll_LP"), c("split_coll_epi","whole_coll"), c("split_coll_epi","split_prot"),  c("split_coll_LP","whole_coll"), c("split_coll_LP","split_prot"), c("whole_coll","split_prot"))

ggplot (proportions_B, aes (method,Freq,fill=method)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol")

ggplot (proportions_T, aes (method,Freq,fill=method)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol")

ggplot (proportions_epi, aes (method,Freq,fill=method)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol")

alt_comparisions=list(   c("split_coll_LP","whole_coll"), c("split_coll_LP","split_prot"), c("whole_coll","split_prot"))

ggplot (proportions_endo, aes (method,Freq,fill=method)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = alt_comparisions, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol") 

ggplot (proportions_myeloid, aes (method,Freq,fill=method)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol") 


ggplot (proportions_fib, aes (method,Freq,fill=method)) + geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = alt_comparisions, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol") 


smillie_meta<-read.table("~/Desktop/smillie_2019/all.meta2.txt", sep="\t", header=T)
smillie_meta<-smillie_meta[-1,]

a<-read.csv("~/Desktop/methods_paper/test_integrations/epi_healthy_sub_table_7pct4.csv")
b<-read.csv("~/Desktop/methods_paper/test_integrations/epi_healthy_sub_table_7pct5.csv")
x<-read.csv("~/Desktop/methods_paper/test_integrations/epi_healthy_sub_table_7pct1.csv")
y<-read.csv("~/Desktop/methods_paper/test_integrations/epi_healthy_sub_table_7pct2.csv")
z<-read.csv("~/Desktop/methods_paper/test_integrations/epi_healthy_sub_table_7pct3.csv")

x$X<-"first"
y$X<-"second"
z$X<-"third"
a$X<-"fourth"
b$X<-"fifth"

c<-rbind(a,b)
c<-rbind(c,x)
c<-rbind(c,y)
c<-rbind(c,z)

d <- dcast(c, X~Var1)

row.names(d)<-d$X
d<-d[-1]


results<-matrix(nrow = 52, ncol = 2)
colnames(results)<-c("mean", "sd")
for( i in 1:ncol(d)){
results[i,1]=mean(d[,i])
results[i,2]=sd(d[,i])
}


celltype_sample_smillie_full<-data.frame(prop.table(table(smillie_meta_full$Predicted_all,smillie_meta_full$Sample, smillie_meta_full$Location),2))




rownames(results)<-colnames(d)

results[1,]=mean(d[,1])

celltype_sample_smillie_full<-data.frame(prop.table(table(smillie_meta_full$Predicted_all,smillie_meta_full$Sample),2))
celltype_sample_smillie_full<-celltype_sample_smillie_full[celltype_sample_smillie_full$Var2 != "group",]

celltype_sample_smillie_full_epi<-celltype_sample_smillie_full[celltype_sample_smillie_full$Var3 == "Epi",]
ggplot(celltype_sample_smillie_full_epi, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2")) 

celltype_sample_smillie_full_LP<-celltype_sample_smillie_full[celltype_sample_smillie_full$Var3 == "LP",]
ggplot(celltype_sample_smillie_full_LP, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2")) 



celltype_sample_smillie_full_epi<-celltype_sample_smillie_full[celltype_sample_smillie_full$Var3 == "Epi",]
celltype_sample_smillie_full_epi<-celltype_sample_smillie_full[celltype_sample_smillie_full$Var3 == "Epi",]

celltype_sample_smillie_full_epi$Var2 = with(celltype_sample_smillie_full_epi, reorder(Var2, Freq))

celltype_sample_smillie_full_epi<-droplevels(celltype_sample_smillie_full_epi)
celltype_sample_smillie_full_epi$Genotype <- factor(celltype_sample_smillie_full_epi$Genotype, levels = c("Genotype 2", "Genotype 3", "Genotype 1")).
ggplot(celltype_sample_smillie_full_epi, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2")) 

celltype_sample_smillie_full_LP<-celltype_sample_smillie_full[celltype_sample_smillie_full$Var3 == "LP",]
ggplot(celltype_sample_smillie_full_LP, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2")) 


celltype_sample_smillie_full_epi<-celltype_sample_smillie_full[grepl(celltype_sample_smillie_full$Var2, pattern = "Epi"),]
celltype_sample_smillie_full_LP<-celltype_sample_smillie_full[grepl(celltype_sample_smillie_full$Var2, pattern = "LP"),]

ggplot(celltype_sample_smillie_full_epi, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2")) 
ggplot(celltype_sample_smillie_full_LP, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2")) 


celltype_sample_smillie_full_epi <- celltype_sample_smillie_full_epi[order(celltype_sample_smillie_full_epi$Freq, decreasing = TRUE),]
ggplot(celltype_sample_smillie_full_epi, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2")) 


ggplot(proportions_whole_coll, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2")) 

ggplot(proportions_split_prot, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","deepskyblue2", "hotpink2")) 

## look into epi only protease samples

hc1<-readRDS("~/Desktop/methods_paper/Sangerdata_202004/OTARscRNA8355919/HC_Sanger1_dataset_sct.rds")
# etc for HC1,3,4,5,7,8

hc1<-data.frame(prop.table(table(hc1@meta.data$Predicted_all)))
hc1$Var2<-"HC1"
hc3<-data.frame(prop.table(table(hc3@meta.data$Predicted_all)))
hc3$Var2<-"HC3"
hc4<-data.frame(prop.table(table(hc4@meta.data$Predicted_all)))
hc4$Var2<-"HC4"
hc5<-data.frame(prop.table(table(hc5@meta.data$Predicted_all)))
hc5$Var2<-"HC5"
hc7<-data.frame(prop.table(table(hc7@meta.data$Predicted_all)))
hc7$Var2<-"HC7"
hc8<-data.frame(prop.table(table(hc8@meta.data$Predicted_all)))
hc8$Var2<-"HC8"

split_prot_epi<-rbind(hc1, hc3)
split_prot_epi<-rbind(split_prot_epi, hc4)
split_prot_epi<-rbind(split_prot_epi, hc5)
split_prot_epi<-rbind(split_prot_epi, hc7)
split_prot_epi<-rbind(split_prot_epi, hc8)

split_prot_epi$Var1 = with(split_prot_epi, reorder(Var1))
ggplot(split_prot_epi, aes(Var2,Freq, fill=Var1))+
  geom_bar(position="stack", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4", "deepskyblue2", "hotpink2")) 

split_prot_epi$protocol<-"protease_epi"
proportions_split_prot$protocol<-"protease_mix"
proportions_whole_coll$protocol<-"whole_collagenase"
celltype_sample_smillie_full_epi$protocol<-"collagenase_epi"
celltype_sample_smillie_full_LP$protocol<-"collagenase_LP"

all_samples_proportions<-rbind(split_prot_epi, proportions_split_prot)
all_samples_proportions<-rbind(all_samples_proportions, proportions_whole_coll)
all_samples_proportions<-rbind(all_samples_proportions, celltype_sample_smillie_full_epi)
all_samples_proportions<-rbind(all_samples_proportions, celltype_sample_smillie_full_LP)




proportions_B<-all_samples_proportions[all_samples_proportions$Var1 == "Bcell",]
proportions_T<-all_samples_proportions[all_samples_proportions$Var1 == "Tcell",]
proportions_epi<-all_samples_proportions[all_samples_proportions$Var1 == "epithelium",]
proportions_endo<-all_samples_proportions[all_samples_proportions$Var1 == "endothelial",]
proportions_myeloid<-all_samples_proportions[all_samples_proportions$Var1 == "myeloid",]
proportions_fib<-all_samples_proportions[all_samples_proportions$Var1 == "fibroblast",]
proportions_glia<-all_samples_proportions[all_samples_proportions$Var1 == "Glia",]
proportions_lymphoid<-all_samples_proportions[all_samples_proportions$Var1 == "Lymphoid",]

my_comparisions=list( c("collagenase_epi","collagenase_LP"), c("collagenase_epi","whole_collagenase"), c("collagenase_epi","protease_mix"),  c("collagenase_LP","whole_collagenase"), c("collagenase_LP","protease_mix"), c("whole_collagenase","protease_mix"), c("protease_epi", "collagenase_epi"),c("protease_epi", "protease_mix"))
ggplot (proportions_B, aes (protocol,Freq,fill=Var1)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol")

ggplot (proportions_T, aes (protocol,Freq,fill=Var1)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol")
ggplot (proportions_epi, aes (protocol,Freq,fill=Var1)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol")


ggplot (proportions_endo, aes (protocol,Freq,fill=Var1)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol")

ggplot (proportions_myeloid, aes (protocol,Freq,fill=Var1)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol")

ggplot (proportions_fib, aes (protocol,Freq,fill=Var1)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol")

my_comparisions_glia=list( c("collagenase_epi","collagenase_LP"), c("collagenase_epi","whole_collagenase"),   c("collagenase_LP","whole_collagenase"),  c("protease_epi", "collagenase_epi"))

ggplot (proportions_glia, aes (protocol,Freq,fill=Var1)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = my_comparisions_glia, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol")

my_comparisions_lymph=list( c("collagenase_epi","collagenase_LP"), c("collagenase_epi","whole_collagenase"),  c("collagenase_LP","whole_collagenase"))

ggplot (proportions_lymphoid, aes (protocol,Freq,fill=Var1)) +  geom_boxplot(width=0.1,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var1, nrow=2) + theme_bw() +scale_fill_manual(values=c("snow","firebrick2","blue2", "grey", "yellow", "pink"))  + stat_compare_means(comparisons = my_comparisions_lymph, method = "wilcox.test" ) + ylab ("Fraction") + xlab ("Protocol")


##  smillie several proportions comparison

prop.table(table(methodspaper2@meta.data$seurat_clusters,methodspaper2@meta.data$method ), 1)
smillie_25<-read.csv("~/Desktop/methods_paper/test_integrations/SCT_integration_25perc_smillie_meta.csv")
prop.table(table(smillie_25$seurat_clusters, smillie_25$method), 1)


# 
library(reactome.db)
library(clusterProfiler)
library(ReactomePA)

# findallmarkers on 100% smillie integrated dataset
Tcell_markers<-read.csv("~/Desktop/methods_paper/test_integrations/Tcell_markers.csv")
Tcell_pos<-Tcell_markers[Tcell_markers$avg_logFC > 0,]
Tcell_neg<-Tcell_markers[Tcell_markers$avg_logFC < 0,]



# select wholecollagenase up genes
wholecollagenase_Tcell_pos_genes<-subset(Tcell_pos, (Tcell_pos$cluster == "wholecollagenase")) 
wholecollagenase_Tcell_pos_genes<-wholecollagenase_Tcell_pos_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_Tcell_pos_genes = bitr(wholecollagenase_Tcell_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_Tcell_pos_entrez <- wholecollagenase_Tcell_pos_genes$ENTREZID
# do pathway analysis
wholecollagenase_Tcell_pos_pathways <- enrichPathway(gene=wholecollagenase_Tcell_pos_genes_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_Tcell_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Tcell_wholecollagenase_upregulated_pathways.csv")


# select splitcollagenase up genes
splitcollagenase_epi_Tcell_pos_genes<-subset(Tcell_pos, (Tcell_pos$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_Tcell_pos_genes<-splitcollagenase_epi_Tcell_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_Tcell_pos_genes = bitr(splitcollagenase_epi_Tcell_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_Tcell_pos_entrez <- splitcollagenase_epi_Tcell_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_Tcell_pos_pathways <- enrichPathway(gene=splitcollagenase_epi_Tcell_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_Tcell_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Tcell_splitcollagenase_epi_upregulated_pathways.csv")


##
# select splitcollagenase up genes
splitcollagenase_lp_Tcell_pos_genes<-subset(Tcell_pos, (Tcell_pos$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_Tcell_pos_genes<-splitcollagenase_lp_Tcell_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_Tcell_pos_genes = bitr(splitcollagenase_lp_Tcell_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_Tcell_pos_entrez <- splitcollagenase_lp_Tcell_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_Tcell_pos_pathways <- enrichPathway(gene=splitcollagenase_lp_Tcell_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_Tcell_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Tcell_splitcollagenase_lp_upregulated_pathways.csv")


###
# select splitcollagenase up genes
splitprotease_epi_Tcell_pos_genes<-subset(Tcell_pos, (Tcell_pos$cluster == "splitprotease_epi")) 
splitprotease_epi_Tcell_pos_genes<-splitprotease_epi_Tcell_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_Tcell_pos_genes = bitr(splitprotease_epi_Tcell_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_Tcell_pos_entrez <- splitprotease_epi_Tcell_pos_genes$ENTREZID
# do pathway analysis
splitprotease_epi_Tcell_pos_pathways <- enrichPathway(gene=splitprotease_epi_Tcell_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_Tcell_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Tcell_splitprotease_epi_upregulated_pathways.csv")

###
# select splitcollagenase up genes
splitprotease_Tcell_pos_genes<-subset(Tcell_pos, (Tcell_pos$cluster == "splitprotease")) 
splitprotease_Tcell_pos_genes<-splitprotease_Tcell_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_Tcell_pos_genes = bitr(splitprotease_Tcell_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_Tcell_pos_entrez <- splitprotease_Tcell_pos_genes$ENTREZID
# do pathway analysis
splitprotease_Tcell_pos_pathways <- enrichPathway(gene=splitprotease_Tcell_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_Tcell_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Tcell_splitprotease_upregulated_pathways.csv")


######
# select wholecollagenase down genes
wholecollagenase_Tcell_neg_genes<-subset(Tcell_neg, (Tcell_neg$cluster == "wholecollagenase")) 
wholecollagenase_Tcell_neg_genes<-wholecollagenase_Tcell_neg_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_Tcell_neg_genes = bitr(wholecollagenase_Tcell_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_Tcell_neg_entrez <- wholecollagenase_Tcell_neg_genes$ENTREZID
# do pathway analysis
wholecollagenase_Tcell_neg_pathways <- enrichPathway(gene=wholecollagenase_Tcell_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_Tcell_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Tcell_wholecollagenase_downregulated_pathways.csv")


# select splitcollagenase down genes
splitcollagenase_epi_Tcell_neg_genes<-subset(Tcell_neg, (Tcell_neg$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_Tcell_neg_genes<-splitcollagenase_epi_Tcell_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_Tcell_neg_genes = bitr(splitcollagenase_epi_Tcell_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_Tcell_neg_entrez <- splitcollagenase_epi_Tcell_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_Tcell_neg_pathways <- enrichPathway(gene=splitcollagenase_epi_Tcell_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_Tcell_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Tcell_splitcollagenase_epi_downregulated_pathways.csv")


##
# select splitcollagenase down genes
splitcollagenase_lp_Tcell_neg_genes<-subset(Tcell_neg, (Tcell_neg$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_Tcell_neg_genes<-splitcollagenase_lp_Tcell_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_Tcell_neg_genes = bitr(splitcollagenase_lp_Tcell_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_Tcell_neg_entrez <- splitcollagenase_lp_Tcell_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_Tcell_neg_pathways <- enrichPathway(gene=splitcollagenase_lp_Tcell_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_Tcell_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Tcell_splitcollagenase_lp_downregulated_pathways.csv")


###
# select splitcollagenase down genes
splitprotease_epi_Tcell_neg_genes<-subset(Tcell_neg, (Tcell_neg$cluster == "splitprotease_epi")) 
splitprotease_epi_Tcell_neg_genes<-splitprotease_epi_Tcell_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_Tcell_neg_genes = bitr(splitprotease_epi_Tcell_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_Tcell_neg_entrez <- splitprotease_epi_Tcell_neg_genes$ENTREZID
# do pathway analysis
splitprotease_epi_Tcell_neg_pathways <- enrichPathway(gene=splitprotease_epi_Tcell_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_Tcell_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Tcell_splitprotease_epi_downregulated_pathways.csv")

###
# select splitcollagenase down genes
splitprotease_Tcell_neg_genes<-subset(Tcell_neg, (Tcell_neg$cluster == "splitprotease")) 
splitprotease_Tcell_neg_genes<-splitprotease_Tcell_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_Tcell_neg_genes = bitr(splitprotease_Tcell_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_Tcell_neg_entrez <- splitprotease_Tcell_neg_genes$ENTREZID
# do pathway analysis
splitprotease_Tcell_neg_pathways <- enrichPathway(gene=splitprotease_Tcell_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_Tcell_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Tcell_splitprotease_downregulated_pathways.csv")





Bcell_markers<-read.csv("~/Desktop/methods_paper/test_integrations/Bcell_markers.csv")
Bcell_pos<-Bcell_markers[Bcell_markers$avg_logFC > 0,]
Bcell_neg<-Bcell_markers[Bcell_markers$avg_logFC < 0,]



# select wholecollagenase up genes
wholecollagenase_Bcell_pos_genes<-subset(Bcell_pos, (Bcell_pos$cluster == "wholecollagenase")) 
wholecollagenase_Bcell_pos_genes<-wholecollagenase_Bcell_pos_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_Bcell_pos_genes = bitr(wholecollagenase_Bcell_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_Bcell_pos_entrez <- wholecollagenase_Bcell_pos_genes$ENTREZID
# do pathway analysis
wholecollagenase_Bcell_pos_pathways <- enrichPathway(gene=wholecollagenase_Bcell_pos_genes_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_Bcell_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Bcell_wholecollagenase_upregulated_pathways.csv")


# select splitcollagenase up genes
splitcollagenase_epi_Bcell_pos_genes<-subset(Bcell_pos, (Bcell_pos$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_Bcell_pos_genes<-splitcollagenase_epi_Bcell_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_Bcell_pos_genes = bitr(splitcollagenase_epi_Bcell_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_Bcell_pos_entrez <- splitcollagenase_epi_Bcell_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_Bcell_pos_pathways <- enrichPathway(gene=splitcollagenase_epi_Bcell_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_Bcell_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Bcell_splitcollagenase_epi_upregulated_pathways.csv")


##
# select splitcollagenase up genes
splitcollagenase_lp_Bcell_pos_genes<-subset(Bcell_pos, (Bcell_pos$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_Bcell_pos_genes<-splitcollagenase_lp_Bcell_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_Bcell_pos_genes = bitr(splitcollagenase_lp_Bcell_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_Bcell_pos_entrez <- splitcollagenase_lp_Bcell_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_Bcell_pos_pathways <- enrichPathway(gene=splitcollagenase_lp_Bcell_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_Bcell_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Bcell_splitcollagenase_lp_upregulated_pathways.csv")


###
# select splitcollagenase up genes
splitprotease_epi_Bcell_pos_genes<-subset(Bcell_pos, (Bcell_pos$cluster == "splitprotease_epi")) 
splitprotease_epi_Bcell_pos_genes<-splitprotease_epi_Bcell_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_Bcell_pos_genes = bitr(splitprotease_epi_Bcell_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_Bcell_pos_entrez <- splitprotease_epi_Bcell_pos_genes$ENTREZID
# do pathway analysis
splitprotease_epi_Bcell_pos_pathways <- enrichPathway(gene=splitprotease_epi_Bcell_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_Bcell_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Bcell_splitprotease_epi_upregulated_pathways.csv")

###
# select splitcollagenase up genes
splitprotease_Bcell_pos_genes<-subset(Bcell_pos, (Bcell_pos$cluster == "splitprotease")) 
splitprotease_Bcell_pos_genes<-splitprotease_Bcell_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_Bcell_pos_genes = bitr(splitprotease_Bcell_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_Bcell_pos_entrez <- splitprotease_Bcell_pos_genes$ENTREZID
# do pathway analysis
splitprotease_Bcell_pos_pathways <- enrichPathway(gene=splitprotease_Bcell_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_Bcell_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Bcell_splitprotease_upregulated_pathways.csv")


######
# select wholecollagenase down genes
wholecollagenase_Bcell_neg_genes<-subset(Bcell_neg, (Bcell_neg$cluster == "wholecollagenase")) 
wholecollagenase_Bcell_neg_genes<-wholecollagenase_Bcell_neg_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_Bcell_neg_genes = bitr(wholecollagenase_Bcell_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_Bcell_neg_entrez <- wholecollagenase_Bcell_neg_genes$ENTREZID
# do pathway analysis
wholecollagenase_Bcell_neg_pathways <- enrichPathway(gene=wholecollagenase_Bcell_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_Bcell_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Bcell_wholecollagenase_downregulated_pathways.csv")


# select splitcollagenase down genes
splitcollagenase_epi_Bcell_neg_genes<-subset(Bcell_neg, (Bcell_neg$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_Bcell_neg_genes<-splitcollagenase_epi_Bcell_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_Bcell_neg_genes = bitr(splitcollagenase_epi_Bcell_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_Bcell_neg_entrez <- splitcollagenase_epi_Bcell_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_Bcell_neg_pathways <- enrichPathway(gene=splitcollagenase_epi_Bcell_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_Bcell_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Bcell_splitcollagenase_epi_downregulated_pathways.csv")


##
# select splitcollagenase down genes
splitcollagenase_lp_Bcell_neg_genes<-subset(Bcell_neg, (Bcell_neg$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_Bcell_neg_genes<-splitcollagenase_lp_Bcell_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_Bcell_neg_genes = bitr(splitcollagenase_lp_Bcell_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_Bcell_neg_entrez <- splitcollagenase_lp_Bcell_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_Bcell_neg_pathways <- enrichPathway(gene=splitcollagenase_lp_Bcell_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_Bcell_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Bcell_splitcollagenase_lp_downregulated_pathways.csv")


###
# select splitcollagenase down genes
splitprotease_epi_Bcell_neg_genes<-subset(Bcell_neg, (Bcell_neg$cluster == "splitprotease_epi")) 
splitprotease_epi_Bcell_neg_genes<-splitprotease_epi_Bcell_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_Bcell_neg_genes = bitr(splitprotease_epi_Bcell_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_Bcell_neg_entrez <- splitprotease_epi_Bcell_neg_genes$ENTREZID
# do pathway analysis
splitprotease_epi_Bcell_neg_pathways <- enrichPathway(gene=splitprotease_epi_Bcell_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_Bcell_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Bcell_splitprotease_epi_downregulated_pathways.csv")

###
# select splitcollagenase down genes
splitprotease_Bcell_neg_genes<-subset(Bcell_neg, (Bcell_neg$cluster == "splitprotease")) 
splitprotease_Bcell_neg_genes<-splitprotease_Bcell_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_Bcell_neg_genes = bitr(splitprotease_Bcell_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_Bcell_neg_entrez <- splitprotease_Bcell_neg_genes$ENTREZID
# do pathway analysis
splitprotease_Bcell_neg_pathways <- enrichPathway(gene=splitprotease_Bcell_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_Bcell_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/Bcell_splitprotease_downregulated_pathways.csv")


epitheliums_markers<-read.csv("~/Desktop/methods_paper/test_integrations/epitheliums_markers.csv")
epitheliums_pos<-epitheliums_markers[epitheliums_markers$avg_logFC > 0,]
epitheliums_neg<-epitheliums_markers[epitheliums_markers$avg_logFC < 0,]



# select wholecollagenase up genes
wholecollagenase_epitheliums_pos_genes<-subset(epitheliums_pos, (epitheliums_pos$cluster == "wholecollagenase")) 
wholecollagenase_epitheliums_pos_genes<-wholecollagenase_epitheliums_pos_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_epitheliums_pos_genes = bitr(wholecollagenase_epitheliums_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_epitheliums_pos_entrez <- wholecollagenase_epitheliums_pos_genes$ENTREZID
# do pathway analysis
wholecollagenase_epitheliums_pos_pathways <- enrichPathway(gene=wholecollagenase_epitheliums_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_epitheliums_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/epitheliums_wholecollagenase_upregulated_pathways.csv")


# select splitcollagenase up genes
splitcollagenase_epi_epitheliums_pos_genes<-subset(epitheliums_pos, (epitheliums_pos$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_epitheliums_pos_genes<-splitcollagenase_epi_epitheliums_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_epitheliums_pos_genes = bitr(splitcollagenase_epi_epitheliums_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_epitheliums_pos_entrez <- splitcollagenase_epi_epitheliums_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_epitheliums_pos_pathways <- enrichPathway(gene=splitcollagenase_epi_epitheliums_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_epitheliums_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/epitheliums_splitcollagenase_epi_upregulated_pathways.csv")


##
# select splitcollagenase up genes
splitcollagenase_lp_epitheliums_pos_genes<-subset(epitheliums_pos, (epitheliums_pos$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_epitheliums_pos_genes<-splitcollagenase_lp_epitheliums_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_epitheliums_pos_genes = bitr(splitcollagenase_lp_epitheliums_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_epitheliums_pos_entrez <- splitcollagenase_lp_epitheliums_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_epitheliums_pos_pathways <- enrichPathway(gene=splitcollagenase_lp_epitheliums_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_epitheliums_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/epitheliums_splitcollagenase_lp_upregulated_pathways.csv")


###
# select splitcollagenase up genes
splitprotease_epi_epitheliums_pos_genes<-subset(epitheliums_pos, (epitheliums_pos$cluster == "splitprotease_epi")) 
splitprotease_epi_epitheliums_pos_genes<-splitprotease_epi_epitheliums_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_epitheliums_pos_genes = bitr(splitprotease_epi_epitheliums_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_epitheliums_pos_entrez <- splitprotease_epi_epitheliums_pos_genes$ENTREZID
# do pathway analysis
splitprotease_epi_epitheliums_pos_pathways <- enrichPathway(gene=splitprotease_epi_epitheliums_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_epitheliums_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/epitheliums_splitprotease_epi_upregulated_pathways.csv")

###
# select splitcollagenase up genes
splitprotease_epitheliums_pos_genes<-subset(epitheliums_pos, (epitheliums_pos$cluster == "splitprotease")) 
splitprotease_epitheliums_pos_genes<-splitprotease_epitheliums_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epitheliums_pos_genes = bitr(splitprotease_epitheliums_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epitheliums_pos_entrez <- splitprotease_epitheliums_pos_genes$ENTREZID
# do pathway analysis
splitprotease_epitheliums_pos_pathways <- enrichPathway(gene=splitprotease_epitheliums_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epitheliums_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/epitheliums_splitprotease_upregulated_pathways.csv")


######
# select wholecollagenase down genes
wholecollagenase_epitheliums_neg_genes<-subset(epitheliums_neg, (epitheliums_neg$cluster == "wholecollagenase")) 
wholecollagenase_epitheliums_neg_genes<-wholecollagenase_epitheliums_neg_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_epitheliums_neg_genes = bitr(wholecollagenase_epitheliums_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_epitheliums_neg_entrez <- wholecollagenase_epitheliums_neg_genes$ENTREZID
# do pathway analysis
wholecollagenase_epitheliums_neg_pathways <- enrichPathway(gene=wholecollagenase_epitheliums_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_epitheliums_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/epitheliums_wholecollagenase_downregulated_pathways.csv")


# select splitcollagenase down genes
splitcollagenase_epi_epitheliums_neg_genes<-subset(epitheliums_neg, (epitheliums_neg$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_epitheliums_neg_genes<-splitcollagenase_epi_epitheliums_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_epitheliums_neg_genes = bitr(splitcollagenase_epi_epitheliums_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_epitheliums_neg_entrez <- splitcollagenase_epi_epitheliums_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_epitheliums_neg_pathways <- enrichPathway(gene=splitcollagenase_epi_epitheliums_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_epitheliums_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/epitheliums_splitcollagenase_epi_downregulated_pathways.csv")


##
# select splitcollagenase down genes
splitcollagenase_lp_epitheliums_neg_genes<-subset(epitheliums_neg, (epitheliums_neg$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_epitheliums_neg_genes<-splitcollagenase_lp_epitheliums_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_epitheliums_neg_genes = bitr(splitcollagenase_lp_epitheliums_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_epitheliums_neg_entrez <- splitcollagenase_lp_epitheliums_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_epitheliums_neg_pathways <- enrichPathway(gene=splitcollagenase_lp_epitheliums_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_epitheliums_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/epitheliums_splitcollagenase_lp_downregulated_pathways.csv")


###
# select splitcollagenase down genes
splitprotease_epi_epitheliums_neg_genes<-subset(epitheliums_neg, (epitheliums_neg$cluster == "splitprotease_epi")) 
splitprotease_epi_epitheliums_neg_genes<-splitprotease_epi_epitheliums_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_epitheliums_neg_genes = bitr(splitprotease_epi_epitheliums_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_epitheliums_neg_entrez <- splitprotease_epi_epitheliums_neg_genes$ENTREZID
# do pathway analysis
splitprotease_epi_epitheliums_neg_pathways <- enrichPathway(gene=splitprotease_epi_epitheliums_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_epitheliums_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/epitheliums_splitprotease_epi_downregulated_pathways.csv")

###
# select splitcollagenase down genes
splitprotease_epitheliums_neg_genes<-subset(epitheliums_neg, (epitheliums_neg$cluster == "splitprotease")) 
splitprotease_epitheliums_neg_genes<-splitprotease_epitheliums_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epitheliums_neg_genes = bitr(splitprotease_epitheliums_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epitheliums_neg_entrez <- splitprotease_epitheliums_neg_genes$ENTREZID
# do pathway analysis
splitprotease_epitheliums_neg_pathways <- enrichPathway(gene=splitprotease_epitheliums_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epitheliums_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/epitheliums_splitprotease_downregulated_pathways.csv")

endothelials_markers<-read.csv("~/Desktop/methods_paper/test_integrations/endothelials_markers.csv")
endothelials_pos<-endothelials_markers[endothelials_markers$avg_logFC > 0,]
endothelials_neg<-endothelials_markers[endothelials_markers$avg_logFC < 0,]



# select wholecollagenase up genes
wholecollagenase_endothelials_pos_genes<-subset(endothelials_pos, (endothelials_pos$cluster == "wholecollagenase")) 
wholecollagenase_endothelials_pos_genes<-wholecollagenase_endothelials_pos_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_endothelials_pos_genes = bitr(wholecollagenase_endothelials_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_endothelials_pos_entrez <- wholecollagenase_endothelials_pos_genes$ENTREZID
# do pathway analysis
wholecollagenase_endothelials_pos_pathways <- enrichPathway(gene=wholecollagenase_endothelials_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_endothelials_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/endothelials_wholecollagenase_upregulated_pathways.csv")


# select splitcollagenase up genes
splitcollagenase_epi_endothelials_pos_genes<-subset(endothelials_pos, (endothelials_pos$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_endothelials_pos_genes<-splitcollagenase_epi_endothelials_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_endothelials_pos_genes = bitr(splitcollagenase_epi_endothelials_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_endothelials_pos_entrez <- splitcollagenase_epi_endothelials_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_endothelials_pos_pathways <- enrichPathway(gene=splitcollagenase_epi_endothelials_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_endothelials_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/endothelials_splitcollagenase_epi_upregulated_pathways.csv")


##
# select splitcollagenase up genes
splitcollagenase_lp_endothelials_pos_genes<-subset(endothelials_pos, (endothelials_pos$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_endothelials_pos_genes<-splitcollagenase_lp_endothelials_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_endothelials_pos_genes = bitr(splitcollagenase_lp_endothelials_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_endothelials_pos_entrez <- splitcollagenase_lp_endothelials_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_endothelials_pos_pathways <- enrichPathway(gene=splitcollagenase_lp_endothelials_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_endothelials_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/endothelials_splitcollagenase_lp_upregulated_pathways.csv")


###
# select splitcollagenase up genes
splitprotease_epi_endothelials_pos_genes<-subset(endothelials_pos, (endothelials_pos$cluster == "splitprotease_epi")) 
splitprotease_epi_endothelials_pos_genes<-splitprotease_epi_endothelials_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_endothelials_pos_genes = bitr(splitprotease_epi_endothelials_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_endothelials_pos_entrez <- splitprotease_epi_endothelials_pos_genes$ENTREZID
# do pathway analysis
splitprotease_epi_endothelials_pos_pathways <- enrichPathway(gene=splitprotease_epi_endothelials_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_endothelials_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/endothelials_splitprotease_epi_upregulated_pathways.csv")

###
# select splitcollagenase up genes
splitprotease_endothelials_pos_genes<-subset(endothelials_pos, (endothelials_pos$cluster == "splitprotease")) 
splitprotease_endothelials_pos_genes<-splitprotease_endothelials_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_endothelials_pos_genes = bitr(splitprotease_endothelials_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_endothelials_pos_entrez <- splitprotease_endothelials_pos_genes$ENTREZID
# do pathway analysis
splitprotease_endothelials_pos_pathways <- enrichPathway(gene=splitprotease_endothelials_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_endothelials_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/endothelials_splitprotease_upregulated_pathways.csv")


######
# select wholecollagenase down genes
wholecollagenase_endothelials_neg_genes<-subset(endothelials_neg, (endothelials_neg$cluster == "wholecollagenase")) 
wholecollagenase_endothelials_neg_genes<-wholecollagenase_endothelials_neg_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_endothelials_neg_genes = bitr(wholecollagenase_endothelials_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_endothelials_neg_entrez <- wholecollagenase_endothelials_neg_genes$ENTREZID
# do pathway analysis
wholecollagenase_endothelials_neg_pathways <- enrichPathway(gene=wholecollagenase_endothelials_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_endothelials_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/endothelials_wholecollagenase_downregulated_pathways.csv")


# select splitcollagenase down genes
splitcollagenase_epi_endothelials_neg_genes<-subset(endothelials_neg, (endothelials_neg$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_endothelials_neg_genes<-splitcollagenase_epi_endothelials_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_endothelials_neg_genes = bitr(splitcollagenase_epi_endothelials_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_endothelials_neg_entrez <- splitcollagenase_epi_endothelials_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_endothelials_neg_pathways <- enrichPathway(gene=splitcollagenase_epi_endothelials_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_endothelials_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/endothelials_splitcollagenase_epi_downregulated_pathways.csv")


##
# select splitcollagenase down genes
splitcollagenase_lp_endothelials_neg_genes<-subset(endothelials_neg, (endothelials_neg$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_endothelials_neg_genes<-splitcollagenase_lp_endothelials_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_endothelials_neg_genes = bitr(splitcollagenase_lp_endothelials_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_endothelials_neg_entrez <- splitcollagenase_lp_endothelials_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_endothelials_neg_pathways <- enrichPathway(gene=splitcollagenase_lp_endothelials_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_endothelials_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/endothelials_splitcollagenase_lp_downregulated_pathways.csv")


###
# select splitcollagenase down genes
splitprotease_epi_endothelials_neg_genes<-subset(endothelials_neg, (endothelials_neg$cluster == "splitprotease_epi")) 
splitprotease_epi_endothelials_neg_genes<-splitprotease_epi_endothelials_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_endothelials_neg_genes = bitr(splitprotease_epi_endothelials_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_endothelials_neg_entrez <- splitprotease_epi_endothelials_neg_genes$ENTREZID
# do pathway analysis
splitprotease_epi_endothelials_neg_pathways <- enrichPathway(gene=splitprotease_epi_endothelials_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_endothelials_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/endothelials_splitprotease_epi_downregulated_pathways.csv")

###
# select splitcollagenase down genes
splitprotease_endothelials_neg_genes<-subset(endothelials_neg, (endothelials_neg$cluster == "splitprotease")) 
splitprotease_endothelials_neg_genes<-splitprotease_endothelials_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_endothelials_neg_genes = bitr(splitprotease_endothelials_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_endothelials_neg_entrez <- splitprotease_endothelials_neg_genes$ENTREZID
# do pathway analysis
splitprotease_endothelials_neg_pathways <- enrichPathway(gene=splitprotease_endothelials_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_endothelials_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/endothelials_splitprotease_downregulated_pathways.csv")







fibroblasts_markers<-read.csv("~/Desktop/methods_paper/test_integrations/fibroblasts_markers.csv")
fibroblasts_pos<-fibroblasts_markers[fibroblasts_markers$avg_logFC > 0,]
fibroblasts_neg<-fibroblasts_markers[fibroblasts_markers$avg_logFC < 0,]



# select wholecollagenase up genes
wholecollagenase_fibroblasts_pos_genes<-subset(fibroblasts_pos, (fibroblasts_pos$cluster == "wholecollagenase")) 
wholecollagenase_fibroblasts_pos_genes<-wholecollagenase_fibroblasts_pos_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_fibroblasts_pos_genes = bitr(wholecollagenase_fibroblasts_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_fibroblasts_pos_entrez <- wholecollagenase_fibroblasts_pos_genes$ENTREZID
# do pathway analysis
wholecollagenase_fibroblasts_pos_pathways <- enrichPathway(gene=wholecollagenase_fibroblasts_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_fibroblasts_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/fibroblasts_wholecollagenase_upregulated_pathways.csv")


# select splitcollagenase up genes
splitcollagenase_epi_fibroblasts_pos_genes<-subset(fibroblasts_pos, (fibroblasts_pos$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_fibroblasts_pos_genes<-splitcollagenase_epi_fibroblasts_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_fibroblasts_pos_genes = bitr(splitcollagenase_epi_fibroblasts_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_fibroblasts_pos_entrez <- splitcollagenase_epi_fibroblasts_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_fibroblasts_pos_pathways <- enrichPathway(gene=splitcollagenase_epi_fibroblasts_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_fibroblasts_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/fibroblasts_splitcollagenase_epi_upregulated_pathways.csv")


##
# select splitcollagenase up genes
splitcollagenase_lp_fibroblasts_pos_genes<-subset(fibroblasts_pos, (fibroblasts_pos$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_fibroblasts_pos_genes<-splitcollagenase_lp_fibroblasts_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_fibroblasts_pos_genes = bitr(splitcollagenase_lp_fibroblasts_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_fibroblasts_pos_entrez <- splitcollagenase_lp_fibroblasts_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_fibroblasts_pos_pathways <- enrichPathway(gene=splitcollagenase_lp_fibroblasts_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_fibroblasts_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/fibroblasts_splitcollagenase_lp_upregulated_pathways.csv")


###
# select splitcollagenase up genes
splitprotease_epi_fibroblasts_pos_genes<-subset(fibroblasts_pos, (fibroblasts_pos$cluster == "splitprotease_epi")) 
splitprotease_epi_fibroblasts_pos_genes<-splitprotease_epi_fibroblasts_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_fibroblasts_pos_genes = bitr(splitprotease_epi_fibroblasts_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_fibroblasts_pos_entrez <- splitprotease_epi_fibroblasts_pos_genes$ENTREZID
# do pathway analysis
splitprotease_epi_fibroblasts_pos_pathways <- enrichPathway(gene=splitprotease_epi_fibroblasts_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_fibroblasts_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/fibroblasts_splitprotease_epi_upregulated_pathways.csv")

###
# select splitcollagenase up genes
splitprotease_fibroblasts_pos_genes<-subset(fibroblasts_pos, (fibroblasts_pos$cluster == "splitprotease")) 
splitprotease_fibroblasts_pos_genes<-splitprotease_fibroblasts_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_fibroblasts_pos_genes = bitr(splitprotease_fibroblasts_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_fibroblasts_pos_entrez <- splitprotease_fibroblasts_pos_genes$ENTREZID
# do pathway analysis
splitprotease_fibroblasts_pos_pathways <- enrichPathway(gene=splitprotease_fibroblasts_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_fibroblasts_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/fibroblasts_splitprotease_upregulated_pathways.csv")


######
# select wholecollagenase down genes
wholecollagenase_fibroblasts_neg_genes<-subset(fibroblasts_neg, (fibroblasts_neg$cluster == "wholecollagenase")) 
wholecollagenase_fibroblasts_neg_genes<-wholecollagenase_fibroblasts_neg_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_fibroblasts_neg_genes = bitr(wholecollagenase_fibroblasts_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_fibroblasts_neg_entrez <- wholecollagenase_fibroblasts_neg_genes$ENTREZID
# do pathway analysis
wholecollagenase_fibroblasts_neg_pathways <- enrichPathway(gene=wholecollagenase_fibroblasts_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_fibroblasts_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/fibroblasts_wholecollagenase_downregulated_pathways.csv")


# select splitcollagenase down genes
splitcollagenase_epi_fibroblasts_neg_genes<-subset(fibroblasts_neg, (fibroblasts_neg$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_fibroblasts_neg_genes<-splitcollagenase_epi_fibroblasts_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_fibroblasts_neg_genes = bitr(splitcollagenase_epi_fibroblasts_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_fibroblasts_neg_entrez <- splitcollagenase_epi_fibroblasts_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_fibroblasts_neg_pathways <- enrichPathway(gene=splitcollagenase_epi_fibroblasts_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_fibroblasts_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/fibroblasts_splitcollagenase_epi_downregulated_pathways.csv")


##
# select splitcollagenase down genes
splitcollagenase_lp_fibroblasts_neg_genes<-subset(fibroblasts_neg, (fibroblasts_neg$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_fibroblasts_neg_genes<-splitcollagenase_lp_fibroblasts_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_fibroblasts_neg_genes = bitr(splitcollagenase_lp_fibroblasts_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_fibroblasts_neg_entrez <- splitcollagenase_lp_fibroblasts_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_fibroblasts_neg_pathways <- enrichPathway(gene=splitcollagenase_lp_fibroblasts_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_fibroblasts_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/fibroblasts_splitcollagenase_lp_downregulated_pathways.csv")


###
# select splitcollagenase down genes
splitprotease_epi_fibroblasts_neg_genes<-subset(fibroblasts_neg, (fibroblasts_neg$cluster == "splitprotease_epi")) 
splitprotease_epi_fibroblasts_neg_genes<-splitprotease_epi_fibroblasts_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_fibroblasts_neg_genes = bitr(splitprotease_epi_fibroblasts_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_fibroblasts_neg_entrez <- splitprotease_epi_fibroblasts_neg_genes$ENTREZID
# do pathway analysis
splitprotease_epi_fibroblasts_neg_pathways <- enrichPathway(gene=splitprotease_epi_fibroblasts_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_fibroblasts_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/fibroblasts_splitprotease_epi_downregulated_pathways.csv")

###
# select splitcollagenase down genes
splitprotease_fibroblasts_neg_genes<-subset(fibroblasts_neg, (fibroblasts_neg$cluster == "splitprotease")) 
splitprotease_fibroblasts_neg_genes<-splitprotease_fibroblasts_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_fibroblasts_neg_genes = bitr(splitprotease_fibroblasts_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_fibroblasts_neg_entrez <- splitprotease_fibroblasts_neg_genes$ENTREZID
# do pathway analysis
splitprotease_fibroblasts_neg_pathways <- enrichPathway(gene=splitprotease_fibroblasts_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_fibroblasts_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/fibroblasts_splitprotease_downregulated_pathways.csv")


glias_markers<-read.csv("~/Desktop/methods_paper/test_integrations/glias_markers.csv")
glias_pos<-glias_markers[glias_markers$avg_logFC > 0,]
glias_neg<-glias_markers[glias_markers$avg_logFC < 0,]



# select wholecollagenase up genes
wholecollagenase_glias_pos_genes<-subset(glias_pos, (glias_pos$cluster == "wholecollagenase")) 
wholecollagenase_glias_pos_genes<-wholecollagenase_glias_pos_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_glias_pos_genes = bitr(wholecollagenase_glias_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_glias_pos_entrez <- wholecollagenase_glias_pos_genes$ENTREZID
# do pathway analysis
wholecollagenase_glias_pos_pathways <- enrichPathway(gene=wholecollagenase_glias_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_glias_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/glias_wholecollagenase_upregulated_pathways.csv")


# select splitcollagenase up genes
splitcollagenase_epi_glias_pos_genes<-subset(glias_pos, (glias_pos$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_glias_pos_genes<-splitcollagenase_epi_glias_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_glias_pos_genes = bitr(splitcollagenase_epi_glias_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_glias_pos_entrez <- splitcollagenase_epi_glias_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_glias_pos_pathways <- enrichPathway(gene=splitcollagenase_epi_glias_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_glias_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/glias_splitcollagenase_epi_upregulated_pathways.csv")


##
# select splitcollagenase up genes
splitcollagenase_lp_glias_pos_genes<-subset(glias_pos, (glias_pos$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_glias_pos_genes<-splitcollagenase_lp_glias_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_glias_pos_genes = bitr(splitcollagenase_lp_glias_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_glias_pos_entrez <- splitcollagenase_lp_glias_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_glias_pos_pathways <- enrichPathway(gene=splitcollagenase_lp_glias_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_glias_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/glias_splitcollagenase_lp_upregulated_pathways.csv")


###
# select splitcollagenase up genes
splitprotease_epi_glias_pos_genes<-subset(glias_pos, (glias_pos$cluster == "splitprotease_epi")) 
splitprotease_epi_glias_pos_genes<-splitprotease_epi_glias_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_glias_pos_genes = bitr(splitprotease_epi_glias_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_glias_pos_entrez <- splitprotease_epi_glias_pos_genes$ENTREZID
# do pathway analysis
splitprotease_epi_glias_pos_pathways <- enrichPathway(gene=splitprotease_epi_glias_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_glias_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/glias_splitprotease_epi_upregulated_pathways.csv")

###
# select splitcollagenase up genes
splitprotease_glias_pos_genes<-subset(glias_pos, (glias_pos$cluster == "splitprotease")) 
splitprotease_glias_pos_genes<-splitprotease_glias_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_glias_pos_genes = bitr(splitprotease_glias_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_glias_pos_entrez <- splitprotease_glias_pos_genes$ENTREZID
# do pathway analysis
splitprotease_glias_pos_pathways <- enrichPathway(gene=splitprotease_glias_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_glias_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/glias_splitprotease_upregulated_pathways.csv")


######
# select wholecollagenase down genes
wholecollagenase_glias_neg_genes<-subset(glias_neg, (glias_neg$cluster == "wholecollagenase")) 
wholecollagenase_glias_neg_genes<-wholecollagenase_glias_neg_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_glias_neg_genes = bitr(wholecollagenase_glias_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_glias_neg_entrez <- wholecollagenase_glias_neg_genes$ENTREZID
# do pathway analysis
wholecollagenase_glias_neg_pathways <- enrichPathway(gene=wholecollagenase_glias_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_glias_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/glias_wholecollagenase_downregulated_pathways.csv")


# select splitcollagenase down genes
splitcollagenase_epi_glias_neg_genes<-subset(glias_neg, (glias_neg$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_glias_neg_genes<-splitcollagenase_epi_glias_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_glias_neg_genes = bitr(splitcollagenase_epi_glias_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_glias_neg_entrez <- splitcollagenase_epi_glias_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_glias_neg_pathways <- enrichPathway(gene=splitcollagenase_epi_glias_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_glias_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/glias_splitcollagenase_epi_downregulated_pathways.csv")


##
# select splitcollagenase down genes
splitcollagenase_lp_glias_neg_genes<-subset(glias_neg, (glias_neg$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_glias_neg_genes<-splitcollagenase_lp_glias_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_glias_neg_genes = bitr(splitcollagenase_lp_glias_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_glias_neg_entrez <- splitcollagenase_lp_glias_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_glias_neg_pathways <- enrichPathway(gene=splitcollagenase_lp_glias_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_glias_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/glias_splitcollagenase_lp_downregulated_pathways.csv")


###
# select splitcollagenase down genes
splitprotease_epi_glias_neg_genes<-subset(glias_neg, (glias_neg$cluster == "splitprotease_epi")) 
splitprotease_epi_glias_neg_genes<-splitprotease_epi_glias_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_glias_neg_genes = bitr(splitprotease_epi_glias_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_glias_neg_entrez <- splitprotease_epi_glias_neg_genes$ENTREZID
# do pathway analysis
splitprotease_epi_glias_neg_pathways <- enrichPathway(gene=splitprotease_epi_glias_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_glias_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/glias_splitprotease_epi_downregulated_pathways.csv")

###
# select splitcollagenase down genes
splitprotease_glias_neg_genes<-subset(glias_neg, (glias_neg$cluster == "splitprotease")) 
splitprotease_glias_neg_genes<-splitprotease_glias_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_glias_neg_genes = bitr(splitprotease_glias_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_glias_neg_entrez <- splitprotease_glias_neg_genes$ENTREZID
# do pathway analysis
splitprotease_glias_neg_pathways <- enrichPathway(gene=splitprotease_glias_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_glias_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/glias_splitprotease_downregulated_pathways.csv")


lymphoids_markers<-read.csv("~/Desktop/methods_paper/test_integrations/lymphoids_markers.csv")
lymphoids_pos<-lymphoids_markers[lymphoids_markers$avg_logFC > 0,]
lymphoids_neg<-lymphoids_markers[lymphoids_markers$avg_logFC < 0,]



# select wholecollagenase up genes
wholecollagenase_lymphoids_pos_genes<-subset(lymphoids_pos, (lymphoids_pos$cluster == "wholecollagenase")) 
wholecollagenase_lymphoids_pos_genes<-wholecollagenase_lymphoids_pos_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_lymphoids_pos_genes = bitr(wholecollagenase_lymphoids_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_lymphoids_pos_entrez <- wholecollagenase_lymphoids_pos_genes$ENTREZID
# do pathway analysis
wholecollagenase_lymphoids_pos_pathways <- enrichPathway(gene=wholecollagenase_lymphoids_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_lymphoids_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/lymphoids_wholecollagenase_upregulated_pathways.csv")


# select splitcollagenase up genes
splitcollagenase_epi_lymphoids_pos_genes<-subset(lymphoids_pos, (lymphoids_pos$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_lymphoids_pos_genes<-splitcollagenase_epi_lymphoids_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_lymphoids_pos_genes = bitr(splitcollagenase_epi_lymphoids_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_lymphoids_pos_entrez <- splitcollagenase_epi_lymphoids_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_lymphoids_pos_pathways <- enrichPathway(gene=splitcollagenase_epi_lymphoids_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_lymphoids_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/lymphoids_splitcollagenase_epi_upregulated_pathways.csv")


##
# select splitcollagenase up genes
splitcollagenase_lp_lymphoids_pos_genes<-subset(lymphoids_pos, (lymphoids_pos$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_lymphoids_pos_genes<-splitcollagenase_lp_lymphoids_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_lymphoids_pos_genes = bitr(splitcollagenase_lp_lymphoids_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_lymphoids_pos_entrez <- splitcollagenase_lp_lymphoids_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_lymphoids_pos_pathways <- enrichPathway(gene=splitcollagenase_lp_lymphoids_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_lymphoids_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/lymphoids_splitcollagenase_lp_upregulated_pathways.csv")


###
# select splitcollagenase up genes
splitprotease_epi_lymphoids_pos_genes<-subset(lymphoids_pos, (lymphoids_pos$cluster == "splitprotease_epi")) 
splitprotease_epi_lymphoids_pos_genes<-splitprotease_epi_lymphoids_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_lymphoids_pos_genes = bitr(splitprotease_epi_lymphoids_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_lymphoids_pos_entrez <- splitprotease_epi_lymphoids_pos_genes$ENTREZID
# do pathway analysis
splitprotease_epi_lymphoids_pos_pathways <- enrichPathway(gene=splitprotease_epi_lymphoids_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_lymphoids_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/lymphoids_splitprotease_epi_upregulated_pathways.csv")

###
# select splitcollagenase up genes
splitprotease_lymphoids_pos_genes<-subset(lymphoids_pos, (lymphoids_pos$cluster == "splitprotease")) 
splitprotease_lymphoids_pos_genes<-splitprotease_lymphoids_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_lymphoids_pos_genes = bitr(splitprotease_lymphoids_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_lymphoids_pos_entrez <- splitprotease_lymphoids_pos_genes$ENTREZID
# do pathway analysis
splitprotease_lymphoids_pos_pathways <- enrichPathway(gene=splitprotease_lymphoids_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_lymphoids_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/lymphoids_splitprotease_upregulated_pathways.csv")


######
# select wholecollagenase down genes
wholecollagenase_lymphoids_neg_genes<-subset(lymphoids_neg, (lymphoids_neg$cluster == "wholecollagenase")) 
wholecollagenase_lymphoids_neg_genes<-wholecollagenase_lymphoids_neg_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_lymphoids_neg_genes = bitr(wholecollagenase_lymphoids_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_lymphoids_neg_entrez <- wholecollagenase_lymphoids_neg_genes$ENTREZID
# do pathway analysis
wholecollagenase_lymphoids_neg_pathways <- enrichPathway(gene=wholecollagenase_lymphoids_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_lymphoids_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/lymphoids_wholecollagenase_downregulated_pathways.csv")


# select splitcollagenase down genes
splitcollagenase_epi_lymphoids_neg_genes<-subset(lymphoids_neg, (lymphoids_neg$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_lymphoids_neg_genes<-splitcollagenase_epi_lymphoids_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_lymphoids_neg_genes = bitr(splitcollagenase_epi_lymphoids_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_lymphoids_neg_entrez <- splitcollagenase_epi_lymphoids_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_lymphoids_neg_pathways <- enrichPathway(gene=splitcollagenase_epi_lymphoids_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_lymphoids_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/lymphoids_splitcollagenase_epi_downregulated_pathways.csv")


##
# select splitcollagenase down genes
splitcollagenase_lp_lymphoids_neg_genes<-subset(lymphoids_neg, (lymphoids_neg$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_lymphoids_neg_genes<-splitcollagenase_lp_lymphoids_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_lymphoids_neg_genes = bitr(splitcollagenase_lp_lymphoids_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_lymphoids_neg_entrez <- splitcollagenase_lp_lymphoids_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_lymphoids_neg_pathways <- enrichPathway(gene=splitcollagenase_lp_lymphoids_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_lymphoids_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/lymphoids_splitcollagenase_lp_downregulated_pathways.csv")


###
# select splitcollagenase down genes
splitprotease_epi_lymphoids_neg_genes<-subset(lymphoids_neg, (lymphoids_neg$cluster == "splitprotease_epi")) 
splitprotease_epi_lymphoids_neg_genes<-splitprotease_epi_lymphoids_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_lymphoids_neg_genes = bitr(splitprotease_epi_lymphoids_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_lymphoids_neg_entrez <- splitprotease_epi_lymphoids_neg_genes$ENTREZID
# do pathway analysis
splitprotease_epi_lymphoids_neg_pathways <- enrichPathway(gene=splitprotease_epi_lymphoids_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_lymphoids_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/lymphoids_splitprotease_epi_downregulated_pathways.csv")

###
# select splitcollagenase down genes
splitprotease_lymphoids_neg_genes<-subset(lymphoids_neg, (lymphoids_neg$cluster == "splitprotease")) 
splitprotease_lymphoids_neg_genes<-splitprotease_lymphoids_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_lymphoids_neg_genes = bitr(splitprotease_lymphoids_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_lymphoids_neg_entrez <- splitprotease_lymphoids_neg_genes$ENTREZID
# do pathway analysis
splitprotease_lymphoids_neg_pathways <- enrichPathway(gene=splitprotease_lymphoids_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_lymphoids_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/lymphoids_splitprotease_downregulated_pathways.csv")



myeloids_markers<-read.csv("~/Desktop/methods_paper/test_integrations/myeloids_markers.csv")
myeloids_pos<-myeloids_markers[myeloids_markers$avg_logFC > 0,]
myeloids_neg<-myeloids_markers[myeloids_markers$avg_logFC < 0,]



# select wholecollagenase up genes
wholecollagenase_myeloids_pos_genes<-subset(myeloids_pos, (myeloids_pos$cluster == "wholecollagenase")) 
wholecollagenase_myeloids_pos_genes<-wholecollagenase_myeloids_pos_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_myeloids_pos_genes = bitr(wholecollagenase_myeloids_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_myeloids_pos_entrez <- wholecollagenase_myeloids_pos_genes$ENTREZID
# do pathway analysis
wholecollagenase_myeloids_pos_pathways <- enrichPathway(gene=wholecollagenase_myeloids_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_myeloids_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/myeloids_wholecollagenase_upregulated_pathways.csv")


# select splitcollagenase up genes
splitcollagenase_epi_myeloids_pos_genes<-subset(myeloids_pos, (myeloids_pos$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_myeloids_pos_genes<-splitcollagenase_epi_myeloids_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_myeloids_pos_genes = bitr(splitcollagenase_epi_myeloids_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_myeloids_pos_entrez <- splitcollagenase_epi_myeloids_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_myeloids_pos_pathways <- enrichPathway(gene=splitcollagenase_epi_myeloids_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_myeloids_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/myeloids_splitcollagenase_epi_upregulated_pathways.csv")


##
# select splitcollagenase up genes
splitcollagenase_lp_myeloids_pos_genes<-subset(myeloids_pos, (myeloids_pos$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_myeloids_pos_genes<-splitcollagenase_lp_myeloids_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_myeloids_pos_genes = bitr(splitcollagenase_lp_myeloids_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_myeloids_pos_entrez <- splitcollagenase_lp_myeloids_pos_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_myeloids_pos_pathways <- enrichPathway(gene=splitcollagenase_lp_myeloids_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_myeloids_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/myeloids_splitcollagenase_lp_upregulated_pathways.csv")


###
# select splitcollagenase up genes
splitprotease_epi_myeloids_pos_genes<-subset(myeloids_pos, (myeloids_pos$cluster == "splitprotease_epi")) 
splitprotease_epi_myeloids_pos_genes<-splitprotease_epi_myeloids_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_myeloids_pos_genes = bitr(splitprotease_epi_myeloids_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_myeloids_pos_entrez <- splitprotease_epi_myeloids_pos_genes$ENTREZID
# do pathway analysis
splitprotease_epi_myeloids_pos_pathways <- enrichPathway(gene=splitprotease_epi_myeloids_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_myeloids_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/myeloids_splitprotease_epi_upregulated_pathways.csv")

###
# select splitcollagenase up genes
splitprotease_myeloids_pos_genes<-subset(myeloids_pos, (myeloids_pos$cluster == "splitprotease")) 
splitprotease_myeloids_pos_genes<-splitprotease_myeloids_pos_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_myeloids_pos_genes = bitr(splitprotease_myeloids_pos_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_myeloids_pos_entrez <- splitprotease_myeloids_pos_genes$ENTREZID
# do pathway analysis
splitprotease_myeloids_pos_pathways <- enrichPathway(gene=splitprotease_myeloids_pos_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_myeloids_pos_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/myeloids_splitprotease_upregulated_pathways.csv")


######
# select wholecollagenase down genes
wholecollagenase_myeloids_neg_genes<-subset(myeloids_neg, (myeloids_neg$cluster == "wholecollagenase")) 
wholecollagenase_myeloids_neg_genes<-wholecollagenase_myeloids_neg_genes$gene

# convert gene symbols to entrezid for reactome
wholecollagenase_myeloids_neg_genes = bitr(wholecollagenase_myeloids_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
wholecollagenase_myeloids_neg_entrez <- wholecollagenase_myeloids_neg_genes$ENTREZID
# do pathway analysis
wholecollagenase_myeloids_neg_pathways <- enrichPathway(gene=wholecollagenase_myeloids_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(wholecollagenase_myeloids_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/myeloids_wholecollagenase_downregulated_pathways.csv")


# select splitcollagenase down genes
splitcollagenase_epi_myeloids_neg_genes<-subset(myeloids_neg, (myeloids_neg$cluster == "split_collagenase_epi")) 
splitcollagenase_epi_myeloids_neg_genes<-splitcollagenase_epi_myeloids_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_epi_myeloids_neg_genes = bitr(splitcollagenase_epi_myeloids_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_epi_myeloids_neg_entrez <- splitcollagenase_epi_myeloids_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_epi_myeloids_neg_pathways <- enrichPathway(gene=splitcollagenase_epi_myeloids_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_epi_myeloids_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/myeloids_splitcollagenase_epi_downregulated_pathways.csv")


##
# select splitcollagenase down genes
splitcollagenase_lp_myeloids_neg_genes<-subset(myeloids_neg, (myeloids_neg$cluster == "split_collagenase_lp")) 
splitcollagenase_lp_myeloids_neg_genes<-splitcollagenase_lp_myeloids_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitcollagenase_lp_myeloids_neg_genes = bitr(splitcollagenase_lp_myeloids_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitcollagenase_lp_myeloids_neg_entrez <- splitcollagenase_lp_myeloids_neg_genes$ENTREZID
# do pathway analysis
splitcollagenase_lp_myeloids_neg_pathways <- enrichPathway(gene=splitcollagenase_lp_myeloids_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitcollagenase_lp_myeloids_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/myeloids_splitcollagenase_lp_downregulated_pathways.csv")


###
# select splitcollagenase down genes
splitprotease_epi_myeloids_neg_genes<-subset(myeloids_neg, (myeloids_neg$cluster == "splitprotease_epi")) 
splitprotease_epi_myeloids_neg_genes<-splitprotease_epi_myeloids_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_epi_myeloids_neg_genes = bitr(splitprotease_epi_myeloids_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_epi_myeloids_neg_entrez <- splitprotease_epi_myeloids_neg_genes$ENTREZID
# do pathway analysis
splitprotease_epi_myeloids_neg_pathways <- enrichPathway(gene=splitprotease_epi_myeloids_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_epi_myeloids_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/myeloids_splitprotease_epi_downregulated_pathways.csv")

###
# select splitcollagenase down genes
splitprotease_myeloids_neg_genes<-subset(myeloids_neg, (myeloids_neg$cluster == "splitprotease")) 
splitprotease_myeloids_neg_genes<-splitprotease_myeloids_neg_genes$gene

# convert gene symbols to entrezid for reactome
splitprotease_myeloids_neg_genes = bitr(splitprotease_myeloids_neg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
splitprotease_myeloids_neg_entrez <- splitprotease_myeloids_neg_genes$ENTREZID
# do pathway analysis
splitprotease_myeloids_neg_pathways <- enrichPathway(gene=splitprotease_myeloids_neg_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(splitprotease_myeloids_neg_pathways)
write.csv(y, "~/Desktop/methods_paper/test_integrations/myeloids_splitprotease_downregulated_pathways.csv")












# barplot pathways
#barplot(TregQuiescent_mucosa_down_pathways, showCategory=15)

# dotplot enrichment
#dotplot(TregQuiescent_mucosa_down_pathways, showCategory=15)

# enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
#cnetplot(TregQuiescent_mucosa_down_pathways, categorySize="pvalue", foldChange=geneList)

##### type the fresh vs 

metadata_methodspaper<-read.csv("~/Desktop/methods_paper/metadata_100perc_smillie_20200628.csv")


proportion_celltypes1<-data.frame(prop.table(table(metadata_methodspaper$dataset, metadata_methodspaper$method, metadata_methodspaper$Predicted_all), 1))

proportion_celltypes2<-data.frame(prop.table(table(metadata_methodspaper$Sample, metadata_methodspaper$Location, metadata_methodspaper$Predicted_all), 1))

proportion_celltypes<-rbind(proportion_celltypes2, proportion_celltypes1)

proportion_celltypes<-proportion_celltypes[proportion_celltypes$Freq > 0,]

proportion_celltypes<-proportion_celltypes[proportion_celltypes$Var2 != "splitcollagenase",]


x<-proportion_celltypes

library(stringr)
x$Var2 <- gsub('Epi', 'splitcollagenase_epi', x$Var2)
x$Var2 <- gsub('LP', 'splitcollagenase_lp', x$Var2)



colnames(x)[3]<-"celltypes"
x$celltypes <- gsub('Lymphoid', 'NK_ILC', x$celltypes)
x$celltypes <- gsub('endothelial', 'Endothelial', x$celltypes)
x$celltypes <- gsub('epithelium', 'Epithelial', x$celltypes)
x$celltypes <- gsub('fibroblast', 'Fibroblast', x$celltypes)
x$celltypes <- gsub('myeloid', 'Myeloid', x$celltypes)



ggplot (proportion_celltypes, aes (Var2,Freq,fill=Var3))  +  geom_boxplot( alpha=0.9, outlier.shape = NA) + geom_jitter(position=position_jitter(0.2), cex=0.5) +theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))  + ylab ("Frequency") + xlab ("Protocol") +  facet_wrap(~ Var3, nrow=1, scales = "free_y")
ggplot (x, aes (Var2,Freq,fill=celltypes))  +  geom_boxplot( alpha=0.9, outlier.shape = NA) + geom_jitter(position=position_jitter(0.2), cex=0.5) +theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))  + ylab ("Frequency") + xlab ("Protocol") +  facet_wrap(~ celltypes, nrow=1, scales = "free_y")


ggplot(proportion_celltypes, aes(Var2,Freq, fill=Var3))+
  geom_bar(position="fill", width = 0.8, stat = "identity") + theme_classic() + scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))




library(ggpubr)
my_comparisions=list( c("Epi","splitprotease_epi"), c("LP","wholecollagenase"), c("splitprotease","wholecollagenase"))

ggplot (proportion_celltypes, aes (Var3,Freq,fill=Var3)) + geom_violin() + geom_boxplot(width=0.06,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var2, nrow=2) + theme_bw() +scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Frequency") + xlab ("Protocol")

# best option
ggplot (proportion_celltypes, aes (Var2,Freq,fill=Var3)) + geom_jitter(aes(colour = Var3), position=position_jitter(0.2), cex=1.2)  + theme_bw() +scale_color_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))  + ylab ("Frequency") + xlab ("Protocol") +  facet_wrap(~ Var3, nrow=1, scales = "free_y")

ggplot (proportion_celltypes, aes (Var2,Freq,fill=Var3))  +  geom_boxplot( alpha=0.9, outlier.shape = NA) + geom_jitter(position=position_jitter(0.2), cex=0.5) +theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))  + ylab ("Frequency") + xlab ("Protocol") +  facet_wrap(~ Var3, nrow=1, scales = "free_y")

ggplot (proportion_celltypes, aes (Var3,Freq,fill=Var3)) + geom_violin() + geom_boxplot(width=0.06,position=position_dodge(1), alpha=0.2) + facet_wrap(~ Var2, nrow=2) + theme_bw() +scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))  + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test" ) + ylab ("Frequency") + xlab ("Protocol")

ggplot (proportion_celltypes, aes (Var3,Freq,fill=Var2)) + geom_boxplot( alpha=1)  + theme_bw() +scale_fill_manual(values=c("firebrick2","blue2","gray88","tan3","purple4","gold", "deepskyblue2", "hotpink2"))  + ylab ("Frequency") + xlab ("Protocol")


# als je een bepaalde string (hier "unassigned") in kolom 45 wil verwisselen voor een string in kolom 47 (let op, kolom 43 moet wel character zijn, geen factor) 
metadata<-metadata_methodspaper
metadata$method_layers<-metadata$method


metadata$method_layers<-as.character(metadata$method_layers)
metadata$Location<-as.character(metadata$Location)

for (i in 1:nrow(metadata)){ 
  if (metadata[i,49] == "splitcollagenase"){
    metadata[i,49]<-metadata[i,40] 
  } else {
    metadata[i,49]<-metadata[i,49]
  }
}

metadata$method_layers <- gsub('Epi', 'splitcollagenase_epi', metadata$method_layers)
metadata$method_layers <- gsub('LP', 'splitcollagenase_lp', metadata$method_layers)

x$var4=paste(x$celltypes, x$Var2, sep = "_")

x$Var2<-as.character(x$Var2)
epithelial<-x[x$celltypes == "Epithelial",]
tcell<-x[x$celltypes == "Tcell",]
bcell<-x[x$celltypes == "Bcell",]
endothelial<-x[x$celltypes == "Endothelial",]
glia<-x[x$celltypes == "Glia",]
nk_ilc<-x[x$celltypes == "NK_ILC",]
fibroblast<-x[x$celltypes == "Fibroblast",]
myeloid<-x[x$celltypes == "Myeloid",]


my_results=melt(as.matrix(pairwise.wilcox.test(myeloid$Freq,myeloid$Var2, p.adjust.method = "bonferroni")$p.value))


my_results<-my_results[!is.na(my_results$value),]
write.csv(my_results, "~/Desktop/methods_paper/new_celltype_imputation/results_myeloid_wilcox_bonferronicorr.csv")

write.csv(x, "~/Desktop/methods_paper/new_celltype_imputation/proportions_celltypes.csv" )
