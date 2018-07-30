# Author: WTC
# Date: July 2018
# Script that merges DE genes for CD vs controls from Holgersen paper (JCC 2015) with cell type markers of CD gut mucosacells

## match holgersen with risk genes and DE genes
holgersen<-read.csv("~/Desktop/Single_cell/final_data_paper/holgersen/holgersen.csv", sep=";") # load only the 62 genes that are significant in CD 
mucosamarkers<-read.table("~/Desktop/Single_cell/final_data_paper/mucosacells_eightcelltypes_DEmarkers_1percMAST.txt") # load cell type markers within mucosa
mucosamarkers<-mucosamarkers[mucosamarkers$p_val_adj<0.05,]
holgersen_mucosamarkers<-merge(mucosamarkers, holgersen, by="gene", all=F)
write.csv(holgersen_mucosamarkers,"~/Desktop/Single_cell/final_data_paper/holgersen/overlap_mucosamarkers_cd_holgersen.csv")
