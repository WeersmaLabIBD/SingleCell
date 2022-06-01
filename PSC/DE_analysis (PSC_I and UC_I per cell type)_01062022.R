#######################################
# Generate DE list per cell type      #
# between PSC-I and UC-I              #
#######################################



#######################################
# library                             #
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


#######################################
# function                            #
#######################################

DE_list <- function(celltype, cell_type, path1, path2){
  celltype <- FindMarkers(PSC, subset.ident = cell_type, group.by = "state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "UC-I")
  celltype <- filter(celltype, celltype$p_val_adj < 0.05)
  # create a 'state' column
  celltype$state = NA
  celltype = as.data.frame(celltype)
  celltype$state = ifelse(celltype$avg_log2FC > 0,"PSC_I","UC_I")
  # create a 'gene' column
  celltype$gene = NA
  celltype$gene = rownames(celltype)
  # separate PSC_I up gene list and UC_I up gene list
  PSC_I_up_celltype = filter(celltype, celltype$state == "PSC_I")  
  UC_I_up_celltype = filter(celltype, celltype$state == "PSC_I")
  write.csv(PSC_I_up_celltype, path1)
  write.csv(UC_I_up_celltype, path2)
}


#########################
# main code             #
#########################

PSC <- readRDS("/Applications/Chapters/4 PSC_Werna/PSC_processed_march16_2021.rds")

# how many cell types in PSC dataset?
unique(PSC$celltypes)

#[1] CD8T                    Cycling_TA              IgA_plasma              Absorptive_enterocyte  
#[5] Immature_enterocyte     Glia                    IgG_plasma              MAIT                   
#[9] Myofibroblasts          DUOX2_enterocyte        Mature_B                PLCG2_TA               
#[13] REG_TA                  RSPO3_fibroblast        Stem                    IgM_plasma             
#[17] DC                      Treg                    Activated_B             Macrophage             
#[21] Tuft                    Enteroendocrine         Memory_CD4T             BEST4_enterocyte       
#[25] Endothelial             WNT5B_fibroblast        Pericytes               Inflammatory_fibroblast
#[29] WNT2B_fibroblast        Activated_cycling_B     Inflammatory_monocyte   Absorptive_TA          
#[33] Ribo_TA                 Cycling_T               Goblet                  Immature_goblet        
#[37] MAST 

# PSC-I vs UC-I
Idents(PSC) <- "celltypes"

# Apply each cell type to function

CD8T <- DE_list(CD8T, "CD8T", "/Users/s.qs/Desktop/DE list/PSC_I_up_CD8T.csv", "/Users/s.qs/Desktop/DE list/UC_I_up_CD8T.csv")
Cycling_TA <- DE_list(Cycling_TA, "Cycling_TA", "/Users/s.qs/Desktop/DE list/PSC_I_up_Cycling_TA.csv", "/Users/s.qs/Desktop/DE list/UC_I_up_Cycling_TA.csv")
IgA_plasma <- DE_list(IgA_plasma, "IgA_plasma", "/Users/s.qs/Desktop/DE list/PSC_I_up_IgA_plasma.csv", "/Users/s.qs/Desktop/DE list/UC_I_up_IgA_plasma.csv")
Absorptive_enterocyte <- DE_list(Absorptive_enterocyte, "Absorptive_enterocyte", "/Users/s.qs/Desktop/DE list/PSC_I_up_Absorptive_enterocyte.csv", "/Users/s.qs/Desktop/DE list/UC_I_up_Absorptive_enterocyte.csv")
Immature_enterocyte <- DE_list(Immature_enterocyte, "Immature_enterocyte", "/Users/s.qs/Desktop/DE list/PSC_I_up_Immature_enterocyte.csv", "/Users/s.qs/Desktop/DE list/UC_I_up_Immature_enterocyte.csv")
Glia <- DE_list(Glia, "Glia", "/Users/s.qs/Desktop/DE list/PSC_I_up_Glia.csv", "/Users/s.qs/Desktop/DE list/UC_I_up_Glia.csv")
IgG_plasma <- DE_list(IgG_plasma, "IgG_plasma", "/Users/s.qs/Desktop/DE list/PSC_I_up_IgG_plasma.csv", "/Users/s.qs/Desktop/DE list/UC_I_up_IgG_plasma.csv")
MAIT <- DE_list(MAIT, "MAIT", "/Users/s.qs/Desktop/DE list/MAIT/PSC_I_up_MAIT.csv", "/Users/s.qs/Desktop/DE list/MAIT/UC_I_up_MAIT.csv")
Myofibroblasts <- DE_list(Myofibroblasts, "Myofibroblasts", "/Users/s.qs/Desktop/DE list/Myofibroblasts/PSC_I_up_Myofibroblasts.csv", "/Users/s.qs/Desktop/DE list/Myofibroblasts/UC_I_up_Myofibroblasts.csv")
DUOX2_enterocyte <- DE_list(DUOX2_enterocyte, "DUOX2_enterocyte", "/Users/s.qs/Desktop/DE list/DUOX2_enterocyte/PSC_I_up_DUOX2_enterocyte.csv", "/Users/s.qs/Desktop/DE list/DUOX2_enterocyte/UC_I_up_DUOX2_enterocyte.csv")
Mature_B <- DE_list(Mature_B, "Mature_B", "/Users/s.qs/Desktop/DE list/Mature_B/PSC_I_up_Mature_B.csv", "/Users/s.qs/Desktop/DE list/Mature_B/UC_I_up_Mature_B.csv")
PLCG2_TA <- DE_list(PLCG2_TA, "PLCG2_TA", "/Users/s.qs/Desktop/DE list/PLCG2_TA/PSC_I_up_PLCG2_TA.csv", "/Users/s.qs/Desktop/DE list/PLCG2_TA/UC_I_up_PLCG2_TA.csv")
REG_TA <- DE_list(REG_TA, "REG_TA", "/Users/s.qs/Desktop/DE list/REG_TA/PSC_I_up_REG_TA.csv", "/Users/s.qs/Desktop/DE list/REG_TA/UC_I_up_REG_TA.csv")
RSPO3_fibroblast <- DE_list(RSPO3_fibroblast, "RSPO3_fibroblast", "/Users/s.qs/Desktop/DE list/RSPO3_fibroblast/PSC_I_up_RSPO3_fibroblast.csv", "/Users/s.qs/Desktop/DE list/RSPO3_fibroblast/UC_I_up_RSPO3_fibroblast.csv")
Stem <- DE_list(Stem, "Stem", "/Users/s.qs/Desktop/DE list/Stem/PSC_I_up_Stem.csv", "/Users/s.qs/Desktop/DE list/Stem/UC_I_up_Stem.csv")
IgM_plasma <- DE_list(IgM_plasma, "IgM_plasma", "/Users/s.qs/Desktop/DE list/IgM_plasma/PSC_I_up_IgM_plasma.csv", "/Users/s.qs/Desktop/DE list/IgM_plasma/UC_I_up_IgM_plasma.csv")
DC <- DE_list(DC, "DC", "/Users/s.qs/Desktop/DE list/DC/PSC_I_up_DC.csv", "/Users/s.qs/Desktop/DE list/DC/UC_I_up_DC.csv")
Treg <- DE_list(Treg, "Treg", "/Users/s.qs/Desktop/DE list/Treg/PSC_I_up_TregC.csv", "/Users/s.qs/Desktop/DE list/Treg/UC_I_up_Treg.csv")
Activated_B <- DE_list(Activated_B, "Activated_B", "/Users/s.qs/Desktop/DE list/Activated_B/PSC_I_up_Activated_B.csv", "/Users/s.qs/Desktop/DE list/Activated_B/UC_I_up_Activated_B.csv")
Macrophage <- DE_list(Macrophage, "Macrophage", "/Users/s.qs/Desktop/DE list/Macrophage/PSC_I_up_Macrophage.csv", "/Users/s.qs/Desktop/DE list/Macrophage/UC_I_up_Macrophage.csv")
Tuft <- DE_list(Tuft, "Tuft", "/Users/s.qs/Desktop/DE list/Tuft/PSC_I_up_Tuft.csv", "/Users/s.qs/Desktop/DE list/Tuft/UC_I_up_Tuft.csv")
Enteroendocrine <- DE_list(Enteroendocrine, "Enteroendocrine", "/Users/s.qs/Desktop/DE list/Enteroendocrine/PSC_I_up_Enteroendocrine.csv", "/Users/s.qs/Desktop/DE list/Enteroendocrine/UC_I_up_Enteroendocrine.csv")
Memory_CD4T <- DE_list(Memory_CD4T, "Memory_CD4T", "/Users/s.qs/Desktop/DE list/Memory_CD4T/PSC_I_up_Memory_CD4T.csv", "/Users/s.qs/Desktop/DE list/Memory_CD4T/UC_I_up_Memory_CD4T.csv")
BEST4_enterocyte <- DE_list(BEST4_enterocyte, "BEST4_enterocyte", "/Users/s.qs/Desktop/DE list/BEST4_enterocyte/PSC_I_up_BEST4_enterocyte.csv", "/Users/s.qs/Desktop/DE list/BEST4_enterocyte/UC_I_up_BEST4_enterocyte.csv")
Endothelial <- DE_list(Endothelial, "Endothelial", "/Users/s.qs/Desktop/DE list/Endothelial/PSC_I_up_Endothelial.csv", "/Users/s.qs/Desktop/DE list/Endothelial/UC_I_up_Endothelial.csv")
WNT5B_fibroblast <- DE_list(WNT5B_fibroblast, "WNT5B_fibroblast", "/Users/s.qs/Desktop/DE list/WNT5B_fibroblast/PSC_I_up_WNT5B_fibroblast.csv", "/Users/s.qs/Desktop/DE list/WNT5B_fibroblast/UC_I_up_WNT5B_fibroblast.csv")
Pericytes <- DE_list(Pericytes, "Pericytes", "/Users/s.qs/Desktop/DE list/Pericytes/PSC_I_up_Pericytes.csv", "/Users/s.qs/Desktop/DE list/Pericytes/UC_I_up_Pericytes.csv")
Inflammatory_fibroblast <- DE_list(Inflammatory_fibroblast, "Inflammatory_fibroblast", "/Users/s.qs/Desktop/DE list/Inflammatory_fibroblast/PSC_I_up_Inflammatory_fibroblast.csv", "/Users/s.qs/Desktop/DE list/Inflammatory_fibroblast/UC_I_up_Inflammatory_fibroblast.csv")
WNT2B_fibroblast <- DE_list(WNT2B_fibroblast, "WNT2B_fibroblast", "/Users/s.qs/Desktop/DE list/WNT2B_fibroblast/PSC_I_up_WNT2B_fibroblast.csv", "/Users/s.qs/Desktop/DE list/WNT2B_fibroblast/UC_I_up_WNT2B_fibroblast.csv")
Activated_cycling_B <- DE_list(Activated_cycling_B, "Activated_cycling_B", "/Users/s.qs/Desktop/DE list/Activated_cycling_B/PSC_I_up_Activated_cycling_B.csv", "/Users/s.qs/Desktop/DE list/Activated_cycling_B/UC_I_up_Activated_cycling_B.csv")
Inflammatory_monocyte <- DE_list(Inflammatory_monocyte, "Inflammatory_monocyte", "/Users/s.qs/Desktop/DE list/Inflammatory_monocyte/PSC_I_up_Inflammatory_monocyte.csv", "/Users/s.qs/Desktop/DE list/Inflammatory_monocyte/UC_I_up_Inflammatory_monocyte.csv")
Absorptive_TA <- DE_list(Absorptive_TA, "Absorptive_TA", "/Users/s.qs/Desktop/DE list/Absorptive_TA/PSC_I_up_Absorptive_TA.csv", "/Users/s.qs/Desktop/DE list/Absorptive_TA/UC_I_up_Absorptive_TA.csv")
Ribo_TA <- DE_list(Ribo_TA, "Ribo_TA", "/Users/s.qs/Desktop/DE list/Ribo_TA/PSC_I_up_Ribo_TA.csv", "/Users/s.qs/Desktop/DE list/Ribo_TA/UC_I_up_Ribo_TA.csv")
Cycling_T <- DE_list(Cycling_T, "Cycling_T", "/Users/s.qs/Desktop/DE list/Cycling_T/PSC_I_up_Cycling_T.csv", "/Users/s.qs/Desktop/DE list/Cycling_T/UC_I_up_Cycling_T.csv")
Goblet <- DE_list(Goblet, "Goblet", "/Users/s.qs/Desktop/DE list/Goblet/PSC_I_up_Goblet.csv", "/Users/s.qs/Desktop/DE list/Goblet/UC_I_up_Goblet.csv")
Immature_goblet <- DE_list(Immature_goblet, "Immature_goblet", "/Users/s.qs/Desktop/DE list/Immature_goblet/PSC_I_up_Immature_goblet.csv", "/Users/s.qs/Desktop/DE list/Immature_goblet/UC_I_up_Immature_goblet.csv")
MAST <- DE_list(MAST, "MAST", "/Users/s.qs/Desktop/DE list/MAST/PSC_I_up_MAST.csv", "/Users/s.qs/Desktop/DE list/MAST/UC_I_up_MAST.csv")





















