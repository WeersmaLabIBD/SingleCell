
#
# CellChat for PSC data - functions courtesy of Roy Oelen
#
#.libPaths("/home/xxx/R/x86_64-pc-linux-gnu-library/3.6/") # on cluster

library(reticulate)
use_python("/usr/local/opt/python@3.9/bin/python3.9")
#reticulate::py_install(packages = 'umap-learn'))
library(Seurat)
PSC<-readRDS("~/xxx/PSC_processed_march16_2021.rds") # on PC
#PSC<-readRDS("/groups/xxx/PSC_processed_march16_2021.rds") # on cluster
# subset Seurat object in order to reduce calculation time
PSC<-subset(PSC, downsample=250) # only for use on PC
PSC<-UpdateSeuratObject(PSC)
PSC@meta.data$celltypes<-as.factor(PSC@meta.data$celltypes)
table(PSC@meta.data$celltypes)
PSC_sub<-PSC

# load the libraries

library(Matrix)
library(CellChat)
library(patchwork)
library(ggpubr)
library(dplyr)

# set some options
options(stringsAsFactors = FALSE)

# set functions including default metadata for comparison (adjust if wishful)

########
init_cellchat_object <- function(seurat_object, assay = 'SCT', slot = 'data', ident = 'celltypes'){
  # set the default assay
  DefaultAssay(seurat_object) <- assay
  # extract the data
  data.input <- GetAssayData(seurat_object, assay = assay, slot = slot)
  meta <- seurat_object@meta.data
  # create the object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = ident)
  return(cellchat)
}

preprocess_cellchat_object <- function(chat_object, nthreads=8){
  # set the database
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  chat_object@DB <- CellChatDB.use
  chat_object <- subsetData(chat_object) # This step is necessary even if using the whole database
  # set multithreading options
  future::plan("multiprocess", workers = nthreads) # do parallel
  # get genes
  chat_object <- identifyOverExpressedGenes(chat_object)
  chat_object <- identifyOverExpressedInteractions(chat_object)
  # project gene expression data onto PPI network (optional)
  chat_object <- projectData(chat_object, PPI.human)
  return(chat_object)
}

inference_communication_network <- function(chat_object, min.cells=10, thresh=1){
  # Compute the communication probability and infer cellular communication network
  chat_object <- computeCommunProb(chat_object)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  chat_object <- filterCommunication(chat_object, min.cells = min.cells)
  # Infer the cell-cell communication at a signaling pathway level
  chat_object <- computeCommunProbPathway(chat_object, thresh = thresh)
  # Calculate the aggregated cell-cell communication network
  chat_object <- aggregateNet(chat_object, thresh = thresh)
  return(chat_object)
}

do_default_cellchat_workflow <- function(seurat_object, assay = 'SCT', slot = 'data', ident = 'celltypes', min.cells=10, nthreads=8, thresh=0.05){
  # go through the steps
  chat_object <- init_cellchat_object(seurat_object, assay = assay, slot = slot, ident = ident)
  chat_object <- preprocess_cellchat_object(chat_object, nthreads = nthreads)
  chat_object <- inference_communication_network(chat_object, min.cells = min.cells, thresh = thresh)
  return(chat_object)
}

do_default_cellchat_workflow_per_timepoint <- function(seurat_object, timepoint.column='timepoint.final', assay = 'SCT', slot = 'data', ident = 'celltypes', min.cells=10, nthreads=8, thresh=0.05){
  # save the objects in a list
  chat_per_timepoint <- list()
  # check each cell type
  for(timepoint in unique(seurat_object@meta.data[[timepoint.column]])){
    # subset to that timepoint
    seurat_timepoint <- seurat_object[, seurat_object@meta.data[[timepoint.column]] == timepoint]
    # go through the work flow
    chat_timepoint <- do_default_cellchat_workflow(seurat_timepoint, assay = assay, slot = slot, ident = ident, min.cells=min.cells, nthreads=nthreads, thresh = thresh)
    # add to list
    chat_per_timepoint[[timepoint]] <- chat_timepoint
  }
  return(chat_per_timepoint)
}


plot_communication_network <- function(chat_object){
  groupSize <- as.numeric(table(chat_object@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(chat_object@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(chat_object@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
}

plot_communication_network_separate <- function(chat_object, slot = 'weight'){
  # grab the number of groups
  groupSize <- as.numeric(table(chat_object@idents))
  # grab the right slot
  mat <- chat_object@net[[slot]]
  # make a perfect square of the plot
  par(mfrow = c(ceiling(sqrt(nrow(mat))),ceiling(sqrt(nrow(mat)))), xpd=TRUE)
  # make the plots
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
}

plot_all_communications_networks_separate <- function(chat_object_per_timepoint_and_chem, output_loc='', file_type='pdf', chems=c('V2', 'V3'), width=10, height=10, slot = 'weight'){
  # plot each timepoint
  for(timepoint in names(chat_object_per_timepoint_and_chem[[chems[[1]]]])){
    # plot each chem
    for(chem in chems){
      # paste the full output loc
      output_loc_full <- paste(output_loc, 'communications_separate_', timepoint, '_', chem, '_', slot, '.', file_type, sep = '')
      print(output_loc_full)
      # init where we will save
      if(file_type == 'pdf'){
        pdf(output_loc_full, width = width, height = height)
      }
      else if(file_type == 'png'){
        png(output_loc_full, width = width, height = height)
      }
      else{
        print('unknown file type, doing pdf instead')
        pdf(output_loc_full, width = width, height = height)
      }
      try({
      # grab the number of groups
      groupSize <- as.numeric(table(chat_object_per_timepoint_and_chem[[chem]][[timepoint]]@idents))
      # grab the right slot
      mat <- chat_object_per_timepoint_and_chem[[chem]][[timepoint]]@net[[slot]]
      # make a perfect square of the plot
      par(mfrow = c(ceiling(sqrt(nrow(mat))),ceiling(sqrt(nrow(mat)))), xpd=TRUE)
      # make the plots
      for (i in 1:nrow(mat)) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = paste(rownames(mat)[i], 'in', chem, timepoint))
      }
      })
      dev.off()
    }
  }
  par(mfrow=c(1,1))
}

plot_all_communications_networks <- function(chat_object_per_timepoint_and_chem, output_loc='', file_type='pdf', chems=c('V2', 'V3'), width=10, height=10){
  # plot each timepoint
  for(timepoint in names(chat_object_per_timepoint_and_chem[[chems[[1]]]])){
    # paste the full output loc
    output_loc_full <- paste(output_loc, 'communications_', timepoint, '.', file_type, sep = '')
    print(output_loc_full)
    # init where we will save
    if(file_type == 'pdf'){
      pdf(output_loc_full, width = width, height = height)
    }
    else if(file_type == 'png'){
      png(output_loc_full, width = width, height = height)
    }
    else{
      print('unknown file type, doing pdf instead')
      pdf(output_loc_full, width = width, height = height)
    }
    try({
      # plot each chem
      for(chem in chems){
        groupSize <- as.numeric(table(chat_object_per_timepoint_and_chem[[chem]][[timepoint]]@idents))
        netVisual_circle(chat_object_per_timepoint_and_chem[[chem]][[timepoint]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste("Number of interactions", 'in', chem, timepoint))
        netVisual_circle(chat_object_per_timepoint_and_chem[[chem]][[timepoint]]@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste("Interaction weights/strength", 'in', chem, timepoint))
      }
    })
    dev.off()
  }
}

plot_all_aggregate_pathways <- function(chat_object_per_timepoint_and_chem, output_loc='', file_type='pdf', chems=c('V2', 'V3'), width=10, height=10){
  # plot each timepoint
  for(timepoint in names(chat_object_per_timepoint_and_chem[[chems[[1]]]])){
    # plot each chem
    for(chem in chems){
      # paste the full output loc
      output_loc_full <- paste(output_loc, 'pathways_', timepoint, '_', chem, '.', file_type, sep = '')
      print(output_loc_full)
      # init where we will save
      if(file_type == 'pdf'){
        pdf(output_loc_full, width = width, height = height)
      }
      else if(file_type == 'png'){
        png(output_loc_full, width = width, height = height)
      }
      else{
        print('unknown file type, doing pdf instead')
        pdf(output_loc_full, width = width, height = height)
      }
      try({
        par(mfrow=c(3,2))
        for(pathway in chat_object_per_timepoint_and_chem[[chem]][[timepoint]]@netP$pathways){
          netVisual_aggregate(chat_object_per_timepoint_and_chem[[chem]][[timepoint]], signaling = c(pathway), layout = "circle")
        }
      })
      dev.off()
    }
  }
}

compare_conditions <- function(all_conditions_list, condition.1, condition.2){
  # subset to only the conditions we want
  conditions.list <- list(condition.1 = all_conditions_list[[condition.1]], condition.2 = all_conditions_list[[condition.2]])
  conditions <- mergeCellChat(conditions.list, add.names = names(conditions.list))
  # compare interactions and strength numbers
  gg1 <- compareInteractions(conditions, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(conditions, show.legend = F, group = c(1,2), measure = "weight")
  # check strengths and weights as circle plots
  #gg3 <- plot(netVisual_diffInteraction(conditions, weight.scale = T))
  #gg4 <- plot(netVisual_diffInteraction(conditions, weight.scale = T, measure = "weight"))
  # check strengths and weights as heatmaps
  gg5 <- netVisual_heatmap(conditions)
  gg6 <- netVisual_heatmap(conditions, measure = "weight")
  # compute net neutrality score to show differences
  conditions.list[[1]] <- netAnalysis_computeCentrality(conditions.list[[1]])
  conditions.list[[2]] <- netAnalysis_computeCentrality(conditions.list[[2]])
  # visualize in 2d space, #first get the number of links in both sets
  conditions.list.num.link <- sapply(conditions.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  # set a max and min, so the plots will have the same scale
  conditions.list.weight.MinMax <- c(min(conditions.list.num.link), max(conditions.list.num.link))
  # make the plots
  gg7 <- netAnalysis_signalingRole_scatter(conditions.list[[1]], title = names(conditions.list)[1], weight.MinMax = conditions.list.weight.MinMax)
  gg8 <- netAnalysis_signalingRole_scatter(conditions.list[[2]], title = names(conditions.list)[2], weight.MinMax = conditions.list.weight.MinMax)
  # identify signalling groups based on structure similarity
  conditions <- computeNetSimilarityPairwise(conditions, type = "functional")
  # Compute signaling network similarity for datasets 1 2
  conditions <- netEmbedding(conditions, type = "functional")
  # Manifold learning of the signaling networks for datasets 1 2
  conditions <- netClustering(conditions, type = "functional")
  # Visualization in 2D-space
  gg9 <- netVisual_embeddingPairwise(conditions, type = "functional", label.size = 3.5)
  # same, but on structure similarity this time
  conditions <- computeNetSimilarityPairwise(conditions, type = "structural")
  conditions <- netEmbedding(conditions, type = "structural")
  conditions <- netClustering(conditions, type = "structural")
  gg10 <- netVisual_embeddingPairwise(conditions, type = "structural", label.size = 3.5)
  # put everything in a nice list
  all_plots <- list('interactions_number' = gg1,
                    'interactions_weight' = gg2,
                    #'circle_number' = gg3,
                    #'circle_weight' = gg4,
                    'heatmap_number' = gg5,
                    'heatmap_weight' = gg6,
                    'scatter_cond1' = gg7,
                    'scatter_cond2' = gg8,
                    'clust_functional' = gg9,
                    'clust_structural' = gg10)
  return(all_plots)
}


do_compare_conditions <- function(all_conditions_list, condition.1, condition.2){
  # subset to only the conditions we want
  conditions.list <- list(condition.1 = all_conditions_list[[condition.1]], condition.2 = all_conditions_list[[condition.2]])
  conditions <- mergeCellChat(conditions.list, add.names = names(conditions.list))
  # compute net neutrality score to show differences
  conditions.list[[1]] <- netAnalysis_computeCentrality(conditions.list[[1]])
  conditions.list[[2]] <- netAnalysis_computeCentrality(conditions.list[[2]])
  # identify signalling groups based on structure similarity
  conditions <- computeNetSimilarityPairwise(conditions, type = "functional")
  # Compute signaling network similarity for datasets 1 2
  conditions <- netEmbedding(conditions, type = "functional")
  # Manifold learning of the signaling networks for datasets 1 2
  conditions <- netClustering(conditions, type = "functional")
  # Visualization in 2D-space
  conditions <- computeNetSimilarityPairwise(conditions, type = "structural")
  conditions <- netEmbedding(conditions, type = "structural")
  conditions <- netClustering(conditions, type = "structural")
  return(conditions)
}

chat_result_to_plots <-function(chat_plots, output_loc){
  # setup the output
  pdf(output_loc)
  #par(mfrow = c(1,2))
  #chat_plots[['circle_number']]
  #chat_plots[['circle_weight']]
  print(
    ggarrange(
      plotlist = list(
        chat_plots[['interactions_number']], 
        chat_plots[['interactions_weight']] 
      ), ncol = 1, nrow = 2
    )
  )
  print(chat_plots[['heatmap_number']] + chat_plots[['heatmap_weight']])
  print(
    ggarrange(
      plotlist = list(
        chat_plots[['scatter_cond1']], 
        chat_plots[['scatter_cond2']], 
        chat_plots[['clust_functional']], 
        chat_plots[['clust_structural']]
      ), ncol = 2, nrow = 2
    )
  )
  dev.off()
  # back to default
}

#######

# create cell categories if wishful

#########
PSC@meta.data$cell_cat<-NULL

PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "BEST4_enterocyte"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Immature_enterocyte"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "PLCG2_TA"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "REG_TA"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Ribo_TA"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Absorptive_enterocyte"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "DUOX2_enterocyte"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Enteroendocrine"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Absorptive_TA"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "MT-Hi-enterocyte"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Cycling_TA "]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Stem"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Goblet"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Tuft"]<-"epithelium"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Immature_goblet"]<-"epithelium"


#Tcell
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Memory_CD4T"]<-"Tcell"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Cycling_T"]<-"Tcell"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Treg"]<-"Tcell"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "CD8T"]<-"Tcell"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "MAIT"]<-"Tcell"

#mast
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "MAST"]<-"MASTcell"


#B
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Activated_B"]<-"Bcell"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Activated_cycling_B"]<-"Bcell"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Mature_B"]<-"Bcell"

#plasma
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "IgA_plasma"]<-"Plasma"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "IgM_plasma"]<-"Plasma"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "IgG_plasma"]<-"Plasma"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "MT-Hi_plasma"]<-"Plasma"


#Myeloid
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Macrophage"]<-"Myeloid"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "DC"]<-"Myeloid"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Inflammatory_monocyte"]<-"Myeloid"

#mesenchymal
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Glia"]<-"Glia"

# fibroblasts
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "WNT5B_fibroblast"]<-"fibroblast"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Inflammatory_fibroblast"]<-"fibroblast"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "RSPO3_fibroblast"]<-"fibroblast"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "WNT2B_fibroblast"]<-"fibroblast"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Mt-Hi_stromal"]<-"fibroblast"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Myofibroblasts"]<-"fibroblast"

#endothelial
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Pericytes"]<-"endothelial"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Endothelial"]<-"endothelial"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Mt"]<-"endothelial"
PSC@meta.data$cell_cat[PSC@meta.data$celltypes == "Doublets"]<-"endothelial"

#############

# select only cell types of interest (otherwise plot not readable)
# DUOX2 vs B cell subtypes
PSC <- PSC_sub[ , PSC_sub@meta.data$celltypes %in% c('DUOX2_enterocyte',"Activated_B", "Mature_B","Activated_cycling_B")]
PSC@meta.data$celltypes <- droplevels(PSC@meta.data$celltypes)

#try alternative cell type interactions (see below)
############ 
#DUOX2 vs T cell subtypes
#PSC <- PSC_sub[ , PSC_sub@meta.data$celltypes %in% c('DUOX2_enterocyte', 'CD8T', 'MAIT', 'Treg',"Cycling_T", "Memory_CD4T")]
#PSC@meta.data$celltypes <- droplevels(PSC@meta.data$celltypes)

# BEST4 vs T cell subtypes
#PSC <- PSC_sub[ , PSC_sub@meta.data$celltypes %in% c('BEST4_enterocyte', 'CD8T', 'MAIT', 'Treg',"Cycling_T", "Memory_CD4T")]
#PSC@meta.data$celltypes <- droplevels(PSC@meta.data$celltypes)

# BEST4 vs B cell subtypes
#PSC <- PSC_sub[ , PSC_sub@meta.data$celltypes %in% c('BEST4_enterocyte',"Activated_B", "Mature_B","Activated_cycling_B")]
#PSC@meta.data$celltypes <- droplevels(PSC@meta.data$celltypes)

# BEST4 vs APC subtypes
#PSC <- PSC_sub[ , PSC_sub@meta.data$celltypes %in% c('BEST4_enterocyte',"Inflammatory_monocyte", "Macrophage","DC")]
#PSC@meta.data$celltypes <- droplevels(PSC@meta.data$celltypes)
##########

# test interactions per required disease/HC and inflammation status (below all possible disease/inflammation status combinations)

# subset to specific health status and inflammation status
PSC_Healthy_NonInflamed <- PSC[, PSC@meta.data$disease == 'HC' & PSC@meta.data$inflammation == 'NI']
# create cell type
PSC_Healthy_NonInflamed.chat <- init_cellchat_object(PSC_Healthy_NonInflamed)
PSC_Healthy_NonInflamed.chat <- preprocess_cellchat_object(PSC_Healthy_NonInflamed.chat)
PSC_Healthy_NonInflamed.chat <- inference_communication_network(PSC_Healthy_NonInflamed.chat)
# plot the interactions
plot_communication_network_separate(PSC_Healthy_NonInflamed.chat)


# subset to specific health status and inflammation status
PSC_PSC_NonInflamed <- PSC[, PSC@meta.data$disease == 'PSC' & PSC@meta.data$inflammation == 'NI']
# create cell type
PSC_PSC_NonInflamed.chat <- init_cellchat_object(PSC_PSC_NonInflamed)
PSC_PSC_NonInflamed.chat <- preprocess_cellchat_object(PSC_PSC_NonInflamed.chat)
PSC_PSC_NonInflamed.chat <- inference_communication_network(PSC_PSC_NonInflamed.chat)
# plot the interactions
plot_communication_network_separate(PSC_PSC_NonInflamed.chat)


# subset to specific health status and inflammation status
PSC_PSC_Inflamed <- PSC[, PSC@meta.data$disease == 'PSC' & PSC@meta.data$inflammation == 'I']
# create cell type
PSC_PSC_Inflamed.chat <- init_cellchat_object(PSC_PSC_Inflamed)
PSC_PSC_Inflamed.chat <- preprocess_cellchat_object(PSC_PSC_Inflamed.chat)
PSC_PSC_Inflamed.chat <- inference_communication_network(PSC_PSC_Inflamed.chat)
# plot the interactions
plot_communication_network_separate(PSC_PSC_Inflamed.chat)

# subset to specific health status and inflammation status
PSC_UC_NonInflamed <- PSC[, PSC@meta.data$disease == 'UC' & PSC@meta.data$inflammation == 'NI']
# create cell type
PSC_UC_NonInflamed.chat <- init_cellchat_object(PSC_UC_NonInflamed)
PSC_UC_NonInflamed.chat <- preprocess_cellchat_object(PSC_UC_NonInflamed.chat)
PSC_UC_NonInflamed.chat <- inference_communication_network(PSC_UC_NonInflamed.chat)
# plot the interactions
plot_communication_network_separate(PSC_UC_NonInflamed.chat)


# subset to specific health status and inflammation status
PSC_UC_Inflamed <- PSC[, PSC@meta.data$disease == 'UC' & PSC@meta.data$inflammation == 'I']
# create cell type
PSC_UC_Inflamed.chat <- init_cellchat_object(PSC_UC_Inflamed)
PSC_UC_Inflamed.chat <- preprocess_cellchat_object(PSC_UC_Inflamed.chat)
PSC_UC_Inflamed.chat <- inference_communication_network(PSC_UC_Inflamed.chat)
# plot the interactions
plot_communication_network_separate(PSC_UC_Inflamed.chat)


# comparing the conditions
# put everything in the list
PSC.all.list <- list('PSC_Inflamed' = PSC_PSC_Inflamed.chat, 'PSC_NonInflamed' = PSC_PSC_NonInflamed.chat)
# do analyses we care about
PSC.chatresult.PSC_Inflamed_NonInflamed_DUOX2_Tcells <- compare_conditions(PSC.all.list,  'PSC_Inflamed', 'PSC_NonInflamed')
# save as plots
chat_result_to_plots(PSC.chatresult.PSC_Inflamed_NonInflamed_DUOX2_Tcells, "~/xxx/PSC_Inflamed_NonInflamed_DUOX2_Tcells.pdf")

# plot 1 celltype comparisons
par(mfrow=c(1,2))
netVisual_chord_gene(PSC_PSC_Inflamed.chat, sources.use = 4, targets.use = c(3), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(PSC_PSC_NonInflamed.chat, sources.use = 4, targets.use = c(3), lab.cex = 0.5,legend.pos.y = 30)
#chat_result_to_plots(PSC.chatresult.PSC_NonInflamed_Inflamed, '/groups/xxx/PSC_NonInflamed_Inflamed.pdf')
# we might also want the objects
PSC.chat.PSC_Inflamed_NonInflamed_DUOX2_Tcells <- do_compare_conditions(PSC.all.list, 'PSC_Inflamed', 'PSC_NonInflamed')
saveRDS(PSC.chat.PSC_Inflamed_NonInflamed_DUOX2_Tcells, '/groups/xxx/PSC_Inflamed_NonInflamed_DUOX2_Tcells.rds')

