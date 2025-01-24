#!/bin/bash
#SBATCH --job-name=adj2mod
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/adj2mod_%A_%a.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/adj2mod_%A_%a.err  
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=24:00:00              
#SBATCH --partition=short
#SBATCH --cpus-per-task=36        
#SBATCH --ntasks=1
#SBATCH --array=0-2 

SUBTYPES=("TNBC", "HER2", "TNBC_BRCA")
SUBTYPE=${SUBTYPES[$SLURM_ARRAY_TASK_ID]}

singularity exec \
  -B /home/igarzonalva/Proyecto_SC_TNBC:/home/igarzonalva/Proyecto_SC_TNBC \
  /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/grn/aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
  pyscenic ctx \
  /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/Networks/network_1_${SUBTYPE}.tsv \
  /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/Common_files/databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/Common_files/databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
  /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/Common_files/databases/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/Common_files/databases/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
  --annotations_fname /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/Common_files/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/updated_workflow/loom_files/${SUBTYPE}_expression_matrix.loom \
  --output /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/updated_workflow/Motifs/${SUBTYPE}_motifs.csv \
  --num_workers 34

