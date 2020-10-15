library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyr)
library(MAST)
library(readxl)
library(openxlsx)

data<- readRDS("Data/PSC_202002_integrated_v2_noribo.rds")
data=UpdateSeuratObject(data)
DimPlot(data, label = T) + NoLegend()

data@meta.data$disease = sapply(strsplit(data$Final_HTO,"-"), `[`, 1)

# recluster subset
# epi = average expression EPCAM > 5
epi<-subset(data, idents = c("6","7","17","19","0","14","10","3","21"))
DefaultAssay(epi)<-"integrated"
epi<- RunPCA(epi)
epi <- RunUMAP(epi, dims = 1:30)
epi<-FindNeighbors(epi, dims = 1:30)
epi<-FindClusters(epi, resolution = 0.2)
DimPlot(epi, label=T)
DefaultAssay(epi)<-"RNA"
markers_epi <- FindAllMarkers(epi, only.pos = TRUE)
  selected_markers_epi <- markers_epi %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

# epi2 = epi cluster 2 and 5  
epi2<-subset(epi, idents = c("2","5"))
DefaultAssay(epi2)<-"integrated"
epi2<- RunPCA(epi2)
epi2 <- RunUMAP(epi2, dims = 1:30)
epi2<-FindNeighbors(epi2, dims = 1:30)
epi2<-FindClusters(epi2, resolution = 0.2)
DimPlot(epi2, label=T)
DefaultAssay(epi2)<-"RNA"
markers_epi2 <- FindAllMarkers(epi2, only.pos = TRUE)
  
# leuko = average expression PTPRC > 5
leuko<-subset(data, idents = c("2","11","18","24","27","15","4","26"))
DefaultAssay(leuko)<-"integrated"
leuko<- RunPCA(leuko)
leuko <- RunUMAP(leuko, dims = 1:30)
leuko<-FindNeighbors(leuko, dims = 1:30)
leuko<-FindClusters(leuko, resolution = 0.4)
DimPlot(leuko, label=T)
DefaultAssay(leuko)<-"RNA"
markers_leuko <- FindAllMarkers(leuko, only.pos = TRUE)
#selected_markers_entero <- markers_entero %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

# stromal =  average expression THY1 > 0.4 (fibroblast), SOX10 > 2 (glia), MADCAM1 > 1 (endotheel)
stromal<-subset(PSC_data, idents = c("9", "13", "20", "22", "25"))
DefaultAssay(stromal)<-"integrated"
stromal<- RunPCA(stromal)
stromal <- RunUMAP(stromal, dims = 1:30)
stromal <-FindNeighbors(stromal, dims = 1:30)
stromal <-FindClusters(stromal, resolution = 0.4)
DimPlot(stromal, label=T)
DefaultAssay(stromal)<-"RNA"
markers_stromal <- FindAllMarkers(stromal, only.pos = TRUE)

new.idents<-c("WNT5B+", "Inflammatory_fibroblast", "RSPO3+", "WNT2B+", "endothelial", "endothelial", "Mt-Hi_stromal", "pericytes", "glia", "myofibroblasts", "endothelial", "endo-pericytes")
names(new.idents) <- levels(stromal)
stromal=RenameIdents(stromal,new.idents)
DimPlot(stromal, label=T)
stromal[["celltypes"]] <- Idents(object = stromal)
DimPlot(stromal)
meta_stromal<-stromal@meta.data

# plasma = remaining cells
plasma<-subset(data, idents = c("1", "16", "12", "23", "5", "8"))
DefaultAssay(plasma)<-"integrated"
plasma<- RunPCA(plasma)
plasma <- RunUMAP(plasma, dims = 1:30)
plasma <-FindNeighbors(plasma, dims = 1:30)
plasma <-FindClusters(plasma, resolution = 0.25)
DimPlot(plasma)

new.idents<-c("IgA_plasma", "IgA_plasma", "IgG_plasma", "IgA_plasma", "IgG_plasma", "MT-Hi_plasma", "IgM_plasma", "IgA_plasma", "IgM_plasma")
names(new.idents) <- levels(plasma)
plasma=RenameIdents(plasma,new.idents)
DimPlot(plasma, label=T)
plasma[["celltypes"]] <- Idents(object = plasma)
DimPlot(plasma)

meta_plasma<-plasma@meta.data

colnames(meta_stromal)
colnames(meta_leuko)
colnames(meta_plasma)
meta_plasma<-meta_plasma[-c(25,26,27,28,29)]
meta_leuko<-meta_leuko[-c(25,26,27)]
meta_stromal<-meta_stromal[-c(24,26,27)]
x<-rbind(meta_stromal, meta_leuko)
x<-rbind(x, meta_plasma)
meta_epi<-read.csv("~/Desktop/PSC/202010_analyses/metadata.epi.csv")
colnames(meta_epi)
meta_epi<-meta_epi[-c(26,27,28,29)]
row.names(meta_epi)<-meta_epi$X
meta_epi<-meta_epi[-1]
x<-rbind(x, meta_epi)
write.csv(x, "~/Desktop/PSC/202010_analyses/meta_all.csv")

#add celltypes to PSC main file
x$NAME<-row.names(x)
CellsMeta<-PSC_data@meta.data
colnames(CellsMeta)
CellsMeta$NAME<-rownames(CellsMeta)
x<-x[,c(25,26)]
row.names(CellsMeta)=NULL
row.names(x)<-NULL
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 
CellsMeta<-keeping.order(CellsMeta, merge, y=x, by = "NAME", all=T)
CellsMeta<-CellsMeta[,c(1,25)]
row.names(CellsMeta)<-CellsMeta$NAME
PSC_data<-AddMetaData(PSC_data, CellsMeta)



#Split HTO naam in disease
data@meta.data$disease = sapply(strsplit(data$Final_HTO,"-"), `[`, 1)

vignette pseudotime: 
https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html

write.csv(markers_leuko, "~/Desktop/PSC/202010_analyses/markers_leuko.csv")
new.idents<-c("Follicular_B", "CD4_T_activated", "Treg", "CD4_T_memory", "MT_Hi_T", "CD8T", "APC", "CD8T", "MAST_leuko", "Follicular_B", "APC", "MT_Hi_B", "Doublet_T_plasma", "GC", "MT_Hi_T", "Cycling_B", "APC", "Cycling_T", "Doublet_T_B")
names(new.idents) <- levels(leuko)
leuko=RenameIdents(leuko,new.idents)
DimPlot(leuko, label=T)
leuko[["celltypes"]] <- Idents(object = leuko)
DimPlot(leuko)

#### calculate difference in cell type composition with clean dataset, without sample 3269 and doublets
PSC_data<-readRDS("~/Desktop/PSC/202010_analyses/PSC_processed_oct.rds")

install.packages("provenance")
library(provenance)
sample_table<-data.frame(table(PSC_data@meta.data$celltypes, PSC_data@meta.data$Final_HTO))
sample_table$Var2<-as.character(sample_table$Var2)
sample_table$disease = sapply(strsplit(sample_table$Var2,"-"), `[`, 1)
sample_table$status = sapply(strsplit(sample_table$Var2,"-"), `[`, 2)
sample_table$disease2="HC"
sample_table$disease2[sample_table$disease!="HC"]="IBD"
sample_table$disease3="non-PSC"
sample_table$disease3[sample_table$disease=="PSC"]="PSC"
library(reshape2)
prop2CLR=dcast(Var2+disease+status+disease2+disease3~Var1, value.var = "Freq" , data = sample_table)
rownames(prop2CLR)<-prop2CLR$Var2
prop2CLR<-prop2CLR[-6,]
prop2CLR_numeric<-prop2CLR
rownames(prop2CLR_numeric)<-prop2CLR_numeric$Var2
prop2CLR_numeric<-prop2CLR_numeric[-c(1:5)]

###
transform_and_filter_taxa=function(x, samples_row=T, method="asin", missing_filter=0){
  x[x=="NA"]=0
  x[x=="NaN"]=0
  #if samples are in columns transpose
  if (!samples_row){
    
    x=as.data.frame(t(x))
    
  } 
  #Exclude/keep columns that pass the missigness threshold
  if (missing_filter>100){
    
    stop("\n Hey! \n Values should be a proportion of missing values allowed per column: a value from 0 to 100") 
    
  }
  
  x_filt=x[,((colSums(x !=0) / nrow(x)) *100 )>missing_filter]
  my_num_removed=ncol(x)-ncol(x_filt)
  print (paste(my_num_removed, "species removed due to many missing values"))
  
  if (method=="asin"){
    print ("ASIN")
    x_filt=x_filt/100
    x_filt=asin(sqrt(x_filt))
  } else if (method=="log"){
    print ("LOG10")
    #replace 0 by the half of the smallest value observed
    my_min=min(x_filt[x_filt>0])/2
    x_filt=x_filt+my_min
    x_filt=log10(x_filt)
  }else if (method=="clr"){
    print ("CLR")
    #Adapted from Alexander Kurilshikov 
    #x_filt=x_filt/100
    #replace 0 by the half of the smallest value observed
    my_min=min(x[x>0])/2
    x=x+my_min
    #Calculate geometric mean
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x), na.rm=na.rm) / length(x))
    }
    Gmean_core = apply(x, 1, gm_mean)
    data_prepared = cbind(Gmean_core,x)
    d <- t(apply(data_prepared, 1, function(b) {
      log(b/b[1])[-1]
    }))
    x=d
    x_filt=x[,colnames(x) %in%colnames(x_filt)]
  }
  return(as.data.frame(x_filt))
}
###


x<-transform_and_filter_taxa(prop2CLR_numeric, method="clr")
y<-prop2CLR[c(2:5)]
x<-cbind(y, x)
prop2<-x


flag=1
for(i in 5:ncol(prop2)){
  my_cell_type=colnames(prop2)[i]
  result_IBD=wilcox.test(prop2[,i]~prop2$disease2, p.adjust.method = "none")$p.value
  result_PSC=wilcox.test(prop2[,i]~prop2$disease3, p.adjust.method = "none")$p.value
  result_INF=wilcox.test(prop2[,i]~prop2$status, p.adjust.method = "none")$p.value
  my_result=data.frame("Cell_type"=my_cell_type, "HC_vs_IBD"=result_IBD, "non-PSC_vs_PSC"= result_PSC, "I_vs_NI"=result_INF)
  if(flag==1){
    results_prop_cells=my_result
    flag=5
  }else{
    results_prop_cells=rbind(results_prop_cells,my_result)
  }
}

results_prop_cells_fdr=melt(results_prop_cells)
results_prop_cells_fdr$Bonf=p.adjust(results_prop_cells_fdr$value, method = "bonferroni")
results_prop_cells_fdr=dcast(results_prop_cells_fdr, formula=Cell_type~variable, value.var="Bonf")

library(ggplot2)
prop2$sample<-rownames(prop2)
x<-melt(prop2)
ggplot(x, aes(status,value, fill=status)) + geom_boxplot() + theme_bw() + facet_wrap(variable~.)
ggplot(x, aes(disease2,value, fill=disease2)) + geom_boxplot() + theme_bw() + facet_wrap(variable~.)
ggplot(x, aes(disease3,value, fill=disease3)) + geom_boxplot() + theme_bw() + facet_wrap(variable~.)


