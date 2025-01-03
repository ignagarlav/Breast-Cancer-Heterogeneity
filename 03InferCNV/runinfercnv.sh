#!/bin/bash
#SBATCH -J icnv-dock         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/infercnv_%A_%a.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/infercnv_%A_%a.err       
#SBATCH --time=23:00:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo

# Paths for split dataset inputs and outputs
METADATA_FILE="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV/metadata_infercnv_$SLURM_ARRAY_TASK_ID.tsv"
COUNTS_MATRIX="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV/counts_matrix_infercnv_$SLURM_ARRAY_TASK_ID.tsv"
GENE_ORDER_FILE="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV/your_gen_pos.txt"
OUTPUT_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV/split_$SLURM_ARRAY_TASK_ID"

SCRIPT_PATH="/home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/03InferCNV/RunInferCNV.R"

cd /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV 

singularity exec \
  --writable-tmpfs \
  -e \
  -B /home/igarzonalva/Proyecto_SC_TNBC:/home/igarzonalva/Proyecto_SC_TNBC \
  infercnv-1.20.0.simg \
  Rscript --verbose $SCRIPT_PATH \
    $METADATA_FILE \
    $COUNTS_MATRIX \
    $GENE_ORDER_FILE \
    $OUTPUT_DIR


