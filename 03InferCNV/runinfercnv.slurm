#!/bin/bash
#SBATCH -J icnv         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/infercnv_%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/infercnv_%A.err       
#SBATCH --time=72:00:00              
#SBATCH --partition=medium           
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo


# If script1.py was successful, run script2.py

cd /home/igarzonalva/Proyecto_SC_TNBC

module load R/4.4.0-gfbf-2023a
module load R-bundle-CRAN/2024.06-foss-2023a-R-4.4.0
module load JAGS


# Paths for split dataset inputs and outputs
METADATA_FILE="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV/metadata_infercnv_$SLURM_ARRAY_TASK_ID.tsv"
COUNTS_MATRIX="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV/counts_matrix_infercnv_$SLURM_ARRAY_TASK_ID.tsv"
GENE_ORDER_FILE="//home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV/your_gen_pos.txt"
OUTPUT_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV/split_$SLURM_ARRAY_TASK_ID"

vmstat 5 > memory_usage_$SLURM_ARRAY_TASK_ID.log &

Rscript --verbose /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/03InferCNV/RunInferCNV.R \
  $METADATA_FILE \
  $COUNTS_MATRIX \
  $GENE_ORDER_FILE \
  $OUTPUT_DIR
