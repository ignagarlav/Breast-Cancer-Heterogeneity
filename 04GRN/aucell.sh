#!/bin/bash
#SBATCH --job-name=aucell
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/aucell_%A_%a.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/aucell_%A_%a.err  
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
  pyscenic aucell \
  /mnt/exp_matrices/${SUBTYPE}_expression_matrix.loom \
    /mnt/regulones/${SUBTYPE}_regulons.csv \
  --output /mnt/aucell/${SUBTYPE}_aucell.loom \
  --num_workers 20
