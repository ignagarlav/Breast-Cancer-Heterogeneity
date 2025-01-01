#!/bin/bash
#SBATCH -J icnvs         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529_176078/infercnv_%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529_176078/infercnv_%A.err      
#SBATCH --cpus-per-task=64    
#SBATCH --time=24:00:00              
#SBATCH --partition=medium           
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo


cd /home/igarzonalva/Proyecto_SC_TNBC

module load R/4.4.0-gfbf-2023a
module load R-bundle-CRAN/2024.06-foss-2023a-R-4.4.0
module load JAGS


# Paths for split dataset inputs and outputs
METADATA_FILE="/home/igarzonalva/Proyecto_SC_TNBC/data/GSE161529/alltypes/InferCNV/metadata_infercnv_$SLURM_ARRAY_TASK_ID.tsv"
COUNTS_MATRIX="/home/igarzonalva/Proyecto_SC_TNBC/data/GSE161529/alltypes/InferCNV/counts_matrix_infercnv_$SLURM_ARRAY_TASK_ID.tsv"
GENE_ORDER_FILE="/home/igarzonalva/Proyecto_SC_TNBC/data/GSE162929_176078/InferCNV/your_gen_pos.txt"
OUTPUT_DIR="/home/igarzonalva/Proyecto_SC_TNBC/data/GSE161529/alltypes/InferCNV/split_$SLURM_ARRAY_TASK_ID"


Rscript --verbose /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/InferCNV/RunInferCNV.R \
  $METADATA_FILE \
  $COUNTS_MATRIX \
  $GENE_ORDER_FILE \
  $OUTPUT_DIR