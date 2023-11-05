library(CellChat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyr)
library(MAST)
library(readxl)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(ggalluvial)
options(stringsAsFactors = FALSE)

data<- readRDS("Data/PSC_processed_oct.rds")
DefaultAssay(data) = "RNA"

# Paired test PSC
Idents(data) = "sample"
data_PSC= subset(data, ident = c("X019","X069","X086","X147","X191","X317","X325"))
data_UC= subset(data, ident = c("X006","X059","X076","X085","X107","X125","X249"))
Idents(data_PSC) = "inflammation"
data_PSC_NI= subset(data_PSC, ident = "NI")
data_PSC_I= subset(data_PSC, ident = "I")
Idents(data_UC) = "inflammation"
data_UC_NI= subset(data_UC, ident = "NI")
data_UC_I= subset(data_UC, ident = "I")

# create random subset of data
cells.to.sample = 2000
set.seed(111)
sampled.cells <- sample(x = colnames(data_PSC_NI), size = cells.to.sample, replace = F)
data_PSC_NI.sub <- subset(data_PSC_NI, cells = sampled.cells)
data_PSC_NI.sub = ScaleData(data_PSC_NI.sub)
Idents(data_PSC_NI.sub) <- "celltypes"

cells.to.sample = 2000
set.seed(111)
sampled.cells <- sample(x = colnames(data_PSC_I), size = cells.to.sample, replace = F)
data_PSC_I.sub <- subset(data_PSC_I, cells = sampled.cells)
ata_PSC_I.sub = ScaleData(data_PSC_I.sub)
Idents(data_PSC_I.sub) <- "celltypes"

# create cellchat objects
data.input <- GetAssayData(data_PSC_I.sub, assay = "RNA", slot = "data") # normalized data matrix
labels_PSCI <- Idents(data_PSC_I.sub)
identity_PSCI <- data.frame(group = labels_PSCI, row.names = names(labels_PSCI)) # create a dataframe of the cell labels
cellchat_PSCI <- createCellChat(data = data.input)

data.input <- GetAssayData(data_PSC_NI.sub, assay = "RNA", slot = "data") # normalized data matrix
labels_PSCNI <- Idents(data_PSC_NI.sub)
identity_PSCNI <- data.frame(group = labels_PSCNI, row.names = names(labels_PSCNI)) # create a dataframe of the cell labels
cellchat_PSCNI <- createCellChat(data = data.input)

# add metadata slot
cellchat_PSCI <- addMeta(cellchat_PSCI, meta = identity_PSCI, meta.name = "labels")
cellchat_PSCI <- setIdent(cellchat_PSCI, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_PSCI@idents) # show factor levels of the cell labels
groupSize_PSCI <- as.numeric(table(cellchat_PSCI@idents))

cellchat_PSCNI <- addMeta(cellchat_PSCNI, meta = identity_PSCNI, meta.name = "labels")
cellchat_PSCNI <- setIdent(cellchat_PSCNI, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_PSCNI@idents) # show factor levels of the cell labels
groupSize_PSCNI <- as.numeric(table(cellchat_PSCNI@idents))

# set up receptor-ligand interaction databse
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat_PSCI@DB <- CellChatDB.use # set the used database in the object
cellchat_PSCNI@DB <- CellChatDB.use # set the used database in the object

cellchat_PSCI <- subsetData(cellchat_PSCI) # subset the expression data of signaling genes for saving computation cost
cellchat_PSCNI <- subsetData(cellchat_PSCNI) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel

# compute interactions
cellchat_PSCNI <- identifyOverExpressedGenes(cellchat_PSCNI)
cellchat_PSCNI <- identifyOverExpressedInteractions(cellchat_PSCNI)
cellchat_PSCNI <- projectData(cellchat_PSCNI, PPI.human)
cellchat_PSCNI <- computeCommunProb(cellchat_PSCNI)
cellchat_PSCNI <- computeCommunProbPathway(cellchat_PSCNI)
cellchat_PSCNI <- aggregateNet(cellchat_PSCNI)


cellchat_PSCI <- identifyOverExpressedGenes(cellchat_PSCI)
cellchat_PSCI <- identifyOverExpressedInteractions(cellchat_PSCI)
cellchat_PSCI <- projectData(cellchat_PSCI, PPI.human)
cellchat_PSCI <- computeCommunProb(cellchat_PSCI)
cellchat_PSCI <- computeCommunProbPathway(cellchat_PSCI)
cellchat_PSCI <- aggregateNet(cellchat_PSCI)

# merge and compute stuff
cellchat_PSCmerged <- mergeCellChat(list(cellchat_PSCI, cellchat_PSCNI), add.names = c("PSCI","PSCNI"))
cellchat_PSCmerged <- computeNetSimilarityPairwise(cellchat_PSCmerged,  type = "functional")
cellchat_PSCmerged <- netEmbedding(cellchat_PSCmerged,  type = "functional")
cellchat_PSCmerged <- netClustering(cellchat_PSCmerged,  type = "functional")

netVisual_embeddingPairwise(cellchat_PSCmerged, type = "functional")
netVisual_embeddingPairwiseZoomIn(cellchat_PSCmerged, type = "functional")

rankSimilarity(cellchat_PSCmerged, type = "functional")
ggsave("Results/Figures/Pathwaydistance_PSC.pdf")
rankNet(cellchat_PSCmerged, mode = "comparison")

# FOR UC
# create random subset of data
cells.to.sample = 2000
set.seed(111)
sampled.cells <- sample(x = colnames(data_UC_NI), size = cells.to.sample, replace = F)
data_UC_NI.sub <- subset(data_UC_NI, cells = sampled.cells)
data_UC_NI.sub = ScaleData(data_UC_NI.sub)
Idents(data_UC_NI.sub) <- "celltypes"

data_UC_I <-data_UC_I[,data_UC_I$celltypes != "Cycling_B"]
cells.to.sample = 2000
set.seed(111)
sampled.cells <- sample(x = colnames(data_UC_I), size = cells.to.sample, replace = F)
data_UC_I.sub <- subset(data_UC_I, cells = sampled.cells)
data_UC_I.sub = ScaleData(data_UC_I.sub)
Idents(data_UC_I.sub) <- "celltypes"

# create cellchat objects
data.input <- GetAssayData(data_UC_I.sub, assay = "RNA", slot = "data") # normalized data matrix
labels_UCI <- Idents(data_UC_I.sub)
identity_UCI <- data.frame(group = labels_UCI, row.names = names(labels_UCI)) # create a dataframe of the cell labels
cellchat_UCI <- createCellChat(data = data.input)

data.input <- GetAssayData(data_UC_NI.sub, assay = "RNA", slot = "data") # normalized data matrix
labels_UCNI <- Idents(data_UC_NI.sub)
identity_UCNI <- data.frame(group = labels_UCNI, row.names = names(labels_UCNI)) # create a dataframe of the cell labels
cellchat_UCNI <- createCellChat(data = data.input)

# add metadata slot
cellchat_UCI <- addMeta(cellchat_UCI, meta = identity_UCI, meta.name = "labels")
cellchat_UCI <- setIdent(cellchat_UCI, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_UCI@idents) # show factor levels of the cell labels
groupSize_UCI <- as.numeric(table(cellchat_UCI@idents))

cellchat_UCNI <- addMeta(cellchat_UCNI, meta = identity_UCNI, meta.name = "labels")
cellchat_UCNI <- setIdent(cellchat_UCNI, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_UCNI@idents) # show factor levels of the cell labels
groupSize_UCNI <- as.numeric(table(cellchat_UCNI@idents))

# set up receptor-ligand interaction databse
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat_UCI@DB <- CellChatDB.use # set the used database in the object
cellchat_UCNI@DB <- CellChatDB.use # set the used database in the object

cellchat_UCI <- subsetData(cellchat_UCI) # subset the expression data of signaling genes for saving computation cost
cellchat_UCNI <- subsetData(cellchat_UCNI) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel

# compute interactions
cellchat_UCNI <- identifyOverExpressedGenes(cellchat_UCNI)
cellchat_UCNI <- identifyOverExpressedInteractions(cellchat_UCNI)
cellchat_UCNI <- projectData(cellchat_UCNI, PPI.human)
cellchat_UCNI <- computeCommunProb(cellchat_UCNI)
cellchat_UCNI <- computeCommunProbPathway(cellchat_UCNI)
cellchat_UCNI <- aggregateNet(cellchat_UCNI)

cellchat_UCI <- identifyOverExpressedGenes(cellchat_UCI)
cellchat_UCI <- identifyOverExpressedInteractions(cellchat_UCI)
cellchat_UCI <- projectData(cellchat_UCI, PPI.human)
cellchat_UCI <- computeCommunProb(cellchat_UCI)
cellchat_UCI <- computeCommunProbPathway(cellchat_UCI)
cellchat_UCI <- aggregateNet(cellchat_UCI)

# merge and compute stuff
cellchat_UCmerged <- mergeCellChat(list(cellchat_UCI, cellchat_UCNI), add.names = c("UCI","UCNI"))
cellchat_UCmerged <- computeNetSimilarityPairwise(cellchat_UCmerged,  type = "functional")
cellchat_UCmerged <- netEmbedding(cellchat_UCmerged,  type = "functional")
cellchat_UCmerged <- netClustering(cellchat_UCmerged,  type = "functional")

netVisual_embeddingPairwise(cellchat_UCmerged, type = "functional")
netVisual_embeddingPairwiseZoomIn(cellchat_UCmerged, type = "functional")

rankSimilarity(cellchat_UCmerged, type = "functional")
ggsave("Results/Figures/Pathwaydistance_UC.pdf")
rankNet(cellchat_UCmerged, mode = "comparison")

saveRDS(cellchat_UCI, file = "Data/cellchat_UCI.rds")
saveRDS(cellchat_UCNI, file = "Data/cellchat_UCNI.rds")
saveRDS(cellchat_UCmerged, file = "Data/cellchat_UCmerged.rds")

# visualize interactions
cellchat_PSCmerged<- readRDS("Data/cellchat_PSCmerged.rds")
cellchat_PSCI<- readRDS("Data/cellchat_PSCI.rds")
cellchat_PSCNI<- readRDS("Data/cellchat_PSCNI.rds")

pathways.show <- c("IFN-II") 
vertex.receiver = c(1,24) # a numeric vector
groupSize_PSCI <- as.numeric(table(cellchat_PSCI@idents))
groupSize_PSCNI <- as.numeric(table(cellchat_PSCNI@idents))
groupSize_UCI <- as.numeric(table(cellchat_UCI@idents))
groupSize_UCNI <- as.numeric(table(cellchat_UCNI@idents))

# Hierarchy plot
netVisual_aggregate(cellchat_PSCI, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize_PSCI)

netVisual_aggregate(cellchat_PSCNI, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "circle", vertex.size = groupSize_PSCNI)
ggsave("Results/Figures/TNF_PSCNI.pdf")
netVisual_aggregate(cellchat_PSCI, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "circle", vertex.size = groupSize_PSCI)
ggsave("Results/Figures/TNF_PSCI.pdf")

netVisual_aggregate(cellchat_UCNI, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "circle", vertex.size = groupSize_UCNI)
ggsave("Results/Figures/TNF_UCNI.pdf")
netVisual_aggregate(cellchat_UCI, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "circle", vertex.size = groupSize_UCI)
ggsave("Results/Figures/TNF_UCI.pdf")

netAnalysis_contribution(cellchat, signaling = pathways.show)

cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP") 

