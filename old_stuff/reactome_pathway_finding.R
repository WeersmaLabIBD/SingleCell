# Author: WTC
# Year: 2018
# Extract enriched pathways from DE genes, using reactomePA package

# open markergenes file
celltype_markergenes<-read.csv("~/Desktop/Single_cell/final_data_paper/DE/allcells_pos_neg_DE_markers.csv")
# select only significant genes
celltype_markergenes<-subset(celltype_markergenes, celltype_markergenes$p_val_adj < 0.05)
# select all upregultated genes
markergenes_up<-celltype_markergenes[celltype_markergenes$avg_logFC >0,]
markergenes_down<-celltype_markergenes[celltype_markergenes$avg_logFC <0,]
write.csv(markergenes_down, "~/Desktop/Single_cell/final_data_paper/DE/downregulated_eight_cell_types_genes_1percMAST.csv")
# select all downreagulated genes

# libraries
library(reactome.db)
library(clusterProfiler)
library(ReactomePA)

# select TregQuiescent_mucosa_down genes
TregQuiescent_mucosa_down_genes<-subset(markergenes_down, (markergenes_down$cluster == "Treg/Quiescent_mucosa")) 
TregQuiescent_mucosa_down_genes<-TregQuiescent_mucosa_down_genes$gene

# convert gene symbols to entrezid for reactome
TregQuiescent_mucosa_down_genes = bitr(TregQuiescent_mucosa_down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# take only Entrez-id column in list
TregQuiescent_mucosa_down_entrez <- TregQuiescent_mucosa_down_genes$ENTREZID
# do pathway analysis
TregQuiescent_mucosa_down_pathways <- enrichPathway(gene=TregQuiescent_mucosa_down_entrez,pvalueCutoff=0.05, readable=T)
# write dataframe with results
y<-as.data.frame(TregQuiescent_mucosa_down_pathways)
write.csv(y, "~/Desktop/Single_cell/final_data_paper/DE/TregQuiescent_mucosa_downregulated_pathways.csv")

# barplot pathways
barplot(TregQuiescent_mucosa_down_pathways, showCategory=15)

# dotplot enrichment
dotplot(TregQuiescent_mucosa_down_pathways, showCategory=15)

# enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
cnetplot(TregQuiescent_mucosa_down_pathways, categorySize="pvalue", foldChange=geneList)


