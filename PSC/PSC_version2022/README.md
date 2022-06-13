This folder contains code and information on the processing and analysis of data for the PSC project. Information on wet lab processes can be found upon publication online. Find below the workflow of the computational analyses.

Count matrices were made directly from raw data using cellranger v3.1.0 with cellranger-GRCh38-3.0.0 as alignment reference
Demultiplexing was done using souporcell (https://github.com/wheaton5/souporcell)
Each 10x lane was preprocessed (HTO demultiplexing, souporcell demultiplexing, QC filtering, SCtransformation) seperately following preprocess.R
Lanes were integrated into one set using SCT_integration_noribo_v2.R following standardized Seurat SCT integration workflow (https://satijalab.org/seurat/v3.0/integration.html)
Post-qc filtering and processing of samples was done following postqc_filtering.R
Doublet removal was done using Scrublet following XXXX
Celltype annotation was done using Azimuth following XXX

Analyses:

Proportion_analyses.R was used for differential abundance analysis
DE_analysis.R was used to identify differentially expressed genes for several comparisons, and GO terms associated
DE_analysis_final_hopelijk.R was used to identify differentially expressed genes, but corrected for number of cells using permutation analysis
Cellchat.R was used to perform cell-cell interaction analysis (https://github.com/sqjin/CellChat). This was done on a subset of cells for efficiency.
riskgene_expression_analysis.R was used to assess PSC and UC risk genes differentially expressed in different health and disease states, for each of the different cell types
