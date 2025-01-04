#!/bin/bash
#SBATCH --job-name=adj2mod
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/adj2mod_%A_%a.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/adj2mod_%A_%a.err  
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=24:00:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1
#SBATCH --array=0-2 

SUBTYPES=("ER" "TNBC" "HER2")
SUBTYPE=${SUBTYPES[$SLURM_ARRAY_TASK_ID]}

singularity exec \
  -B /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/grn:/mnt \
  /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/grn/aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
  pyscenic ctx \
  /mnt/Networks/${SUBTYPE}_network.tsv \
  /mnt/Common_files/databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  /mnt/Common_files/databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
  /mnt/Common_files/databases/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  /mnt/Common_files/databases/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
  --annotations_fname /mnt/Common_files/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname /mnt/exp_matrices/${SUBTYPE}_expression_matrix.loom \
  --output /mnt/regulones/${SUBTYPE}_regulons.csv \
  --num_workers 12
