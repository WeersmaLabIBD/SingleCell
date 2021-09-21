#!/bin/bash
#SBATCH --job-name=200611_lane1
#SBATCH --output=/groups/umcg-weersma/tmp01/Emilia/output/200611_lane1.out
#SBATCH --error=/groups/umcg-weersma/tmp01/Emilia/output/200611_lane1.err
#SBATCH --time=6-23:59:00
#SBATCH --cpus-per-task=22
#SBATCH --mem=100gb
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --tmp=1000gb
mkdir -p ${TMPDIR}/cellranger/200611_lane1
cd ${TMPDIR}/cellranger/200611_lane1
/groups/umcg-weersma/tmp01/singelcell/cellranger-3.1.0/cellranger count \
--id=cellranger_200611_lane1 \
--transcriptome=/groups/umcg-weersma/tmp01/Amber/refdata-cellranger-GRCh38-3.0.0 \
--libraries=/groups/umcg-weersma/tmp01/Emilia/libraries/library_200611_lane1.csv  \
--feature-ref=/groups/umcg-weersma/tmp01/Emilia/featurerefs/feature_ref_200611_lane1.csv \
--localcores=22

mv ${TMPDIR}/cellranger/200611_lane1 /groups/umcg-weersma/tmp01/Emilia/final_output/
