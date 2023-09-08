This folder contains code and information on the processing and analysis of data for the PSC project. Information on wet lab processes can be found upon publication online. Find below the workflow of the computational analyses.

1. Count matrices were made directly from raw data using cellranger v3.1.0 with cellranger-GRCh38-3.0.0 as alignment reference
2. Demultiplexing was done using souporcell (https://github.com/wheaton5/souporcell)
3. Each 10x lane was preprocessed (HTO demultiplexing, souporcell demultiplexing, QC filtering, SCtransformation) seperately following PSC_preprep.R
4. Lanes were integrated into one set using SCT_integration_noribo_v2.R following standardized Seurat SCT integration workflow (https://satijalab.org/seurat/v3.0/integration.html)
5. Post-qc filtering of samples was done following postqc_filtering.R
6. Celltype annotation was done using Azimuth following XXX


Further analyses:
1. DE_loop.R can be used to generate DE genes and enriched GO processes

Addition 2023:
1. DUOX2.R is gebruikt voor de DUOX2 analyse
2. Figure1.R is gebruikt voor figuur 1. De stacked barplot.R ook.
3. basicDE_2023.R is gebruikt voor figuur 3
4. PSC_riskgenes.R is gebruikt voor figuur 4
5. DE_analysis_2023 is gebaseerd op DE_loop.R en is gebruikt om de DE en GO tabellen te genereren
