# from: Gene expression profiles during human CD4+ T cell differentiation https://academic.oup.com/intimm/article/16/8/1109/865613/Gene-expression-profiles-during-human-CD4-T-cell

# cell states
FeaturePlot( b, c("RAG1", "DNTT", "TOP2A" , "TYMS", "RAD51AP1", "MCM2", "MCM10")) # DNA replication, recombination, repair (7:0) 
FeaturePlot( b, c("KIF11", "CCNB2", "BUB1B", "NCAPG", "CCNB1" , "KIF15", "CDK1", "CDKN3", "PRC1", "AURKA", "NUF2", "CENPA", "PTTG1", "NEK2", "NDC80", "CDC25A", "UBE2C", "CKS2", "TTK", "KIF2C", "DLGAP5" , "CYLD")) # Cell cycle regulation, progression, mitosis (21:1) 
FeaturePlot( b, c("CD1E", "CD1A", "HMMR" , "TSHR", "CD99", "TNFRSF21", "GPR171", "CCR7", "IL6R", "GPR183", "SEMA4C",  "IL27RA", "LRRN3", "EMP3",  "IFITM2",  "S1PR1",  "IL4R", "CD44", "CD27")) # Intercellular communication, receptor (7:18) 
FeaturePlot( b, c("GNA15", "SYK",  "RCAN1", "MPP1", "CYTIP", "TBC1D4"))  # Intracellular signaling (5:2) 
FeaturePlot( b, c("AEBP1", "UHRF1", "TFDP2", "MYB", "GFI1", "FOXO1" , "KLF2", "ID2", "SCML1"))  # Transcriptional regulation (5:4)
FeaturePlot( b, c("TUBB2A" ,"ACTN1", "TUBA4A")) # Cytoskeletal component (1:2) 
FeaturePlot( b, c("CPVL"))  # Protein degradation (1:0) 
FeaturePlot( b, c("SLC1A4", "KIF20A", "LDLRAP1" , "SLC2A3", "ANK3")) # Molecular transport (2:3) 
FeaturePlot( b, c("PON1", "GALNT7", "ADA" , "GGH", "RRM2", "MAN1C1", "CHST6", "PTS", "RNASET2")) # Metabolic enzyme (5:4) 
FeaturePlot( b, c("ADCYAP1", "FAIM", "BIRC5", "STK17A", "TNFRSF25", "TNFSF10")) # Apoptosis (3:3) 
# FeaturePlot( b, c("NUSAP1", "RAMP1" , "MELK", "PBK", "AKAP2",  "WHSC1", "LRR1", "CDCA3", "UBE2T", "HMGB3", "NELL2", "P51", "GBP2", "HRMT1L1", "SAMHD1" )) # Miscellaneous (14:6), Unknown (28:8) and EST (19:11)
