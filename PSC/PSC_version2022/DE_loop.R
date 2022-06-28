#######################################
# Generate DE list per cell type      #
#######################################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyr)
library(MAST)
library(readxl)
library(openxlsx)
library(enrichR)

#load in datasets
data <- readRDS("XXXX")
data = UpdateSeuratObject(object = data)

#set idents to final celltyping and assay to RNA
Idents(data) <- "XXXX"             
DefaultAssay(data) = "RNA"

#for loop for DE analysis ##paths and celltype specification should be changed to your likings
celltypes <- unique(data$celltypes_azi_final)
for (x in celltypes) {
  # safe pathnames
  pathPSCI <- paste("/Users/amberbangma/Documents/R/PSC/Results/",x,"_PSC_I_up.csv",sep="")
  pathUCI <- paste("/Users/amberbangma/Documents/R/PSC/Results/",x,"_UC_I_up.csv",sep="")
  pathPSCI_GO <- paste("/Users/amberbangma/Documents/R/PSC/Results/",x,"_PSC_I_up_GO.csv",sep="")
  pathUCI_GO <- paste("/Users/amberbangma/Documents/R/PSC/Results/",x,"_UC_I_up_GO.csv",sep="")
  # find markers
  x <- FindMarkers(data, subset.ident = x, group.by = "state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "UC-I")
  x <- filter(x, x$p_val_adj < 0.05)
  if(nrow(x) == 0)
  {
    #skip iteration because of zero significant genes
    next
  }
  #rest of iteration for case of no error
  # create a 'state' column
  x$state = NA
  x = as.data.frame(x)
  x$state = ifelse(x$avg_log2FC > 0,"PSC_I","UC_I")
  # create a 'gene' column
  x$gene = NA
  x$gene = rownames(x)
  # separate PSC_I up gene list and UC_I up gene list
  PSC_I_up_x = filter(x, x$state == "PSC_I")  
  UC_I_up_x = filter(x, x$state == "UC_I")
  # GO terms for each gene list 
  GO_PSCIup_x <- enrichr(PSC_I_up_x$gene, databases = "GO_Biological_Process_2018")$GO_Biological_Process
  GO_UCIup_x <- enrichr(UC_I_up_x$gene, databases = "GO_Biological_Process_2018")$GO_Biological_Process
  # save dataframes
  write.csv(PSC_I_up_x, pathPSCI)
  write.csv(UC_I_up_x, pathUCI)
  write.csv(GO_PSCIup_x, pathPSCI_GO)
  write.csv(GO_UCIup_x, pathUCI_GO)
}

## Findmarkers parameters and script can be tweeked for different DE analyses (for example UC-I vs UC-NI)
## check DE_analysis_202103.R for examples















  
