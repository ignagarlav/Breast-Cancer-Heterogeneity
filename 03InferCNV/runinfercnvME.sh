#!/bin/bash
#SBATCH -J MEcnv         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/infercnv_%A_%a.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/infercnv_%A_%a.err       
#SBATCH --time=35:00:00              
#SBATCH --partition=medium           
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo

# Define your tumor types in an array
TUMOR_TYPES=("ER" "TNBC" "HER2")

# Select the tumor based on the array index
TUMOR=${TUMOR_TYPES[$SLURM_ARRAY_TASK_ID]}
echo "Tumor type: $TUMOR"

# Paths for split dataset inputs and outputs
METADATA_FILE="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV/meta_mcroenv_cnv_$TUMOR.tsv"
COUNTS_MATRIX="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV/cnt_mt_cnv_mcoenv$TUMOR.tsv"


GENE_ORDER_FILE="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV/your_gen_pos.txt"
OUTPUT_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV/Microenv_$TUMOR"

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


