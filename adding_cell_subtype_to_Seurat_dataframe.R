# Adding new celltype from a subset of cells to a Seurat dataframe

colonocytes<-subset(PSC_dataset, idents = 16)
DefaultAssay(colonocytes)<-"integrated"
colonocytes<-FindNeighbors(colonocytes, dims = 1:30)
colonocytes<-FindClusters(colonocytes, resolution = 0.1)
DimPlot(colonocytes)
table(colonocytes@meta.data$integrated_snn_res.0.1)


colonocytes@meta.data$best4_sub[colonocytes@meta.data$seurat_clusters == "0"]<-"Best4_enterocytes" 
colonocytes@meta.data$best4_sub[colonocytes@meta.data$seurat_clusters == "1"]<-"Best4_enterocytes"
colonocytes@meta.data$best4_sub[colonocytes@meta.data$seurat_clusters == "2"]<-"M_cells" 


CellsMeta<-data.frame(PSC_dataset@meta.data)
x<-data.frame(colonocytes@meta.data)
x$NAME<-row.names(x)
CellsMeta$NAME<-row.names(CellsMeta)
x<-x[,c(17,35)]
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
CellsMeta<-CellsMeta[,c(1,34)]
rownames(CellsMeta)<-CellsMeta$NAME
PSC_dataset<-AddMetaData(PSC_dataset, CellsMeta)

PSC_dataset@meta.data$celltype_subsets[PSC_dataset@meta.data$best4_sub == "Best4_enterocytes"]<-"Best4_enterocytes" 
PSC_dataset@meta.data$celltype_subsets[PSC_dataset@meta.data$best4_sub == "M_cells"]<-"M_cells"

# check whether it worked
table(PSC_dataset@meta.data$best4_sub)
DimPlot(PSC_dataset, group.by="celltype_subsets", label=T, repel=T)


