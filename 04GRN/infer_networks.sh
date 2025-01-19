#!/bin/bash
#SBATCH --job-name=adj2mod
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/adj2mod_%A_%a.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/adj2mod_%A_%a.err  
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=24:00:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36

#SUBTYPES=("ER" "TNBC" "HER2", "TNBC_BRCA")
#SUBTYPE=${SUBTYPES[$SLURM_ARRAY_TASK_ID]}

singularity exec \
  -B /home/igarzonalva/Proyecto_SC_TNBC/GSE161529:/mnt \
  /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/grn/aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
  pyscenic grn \
  /mnt/04_grn/Jan19Try/TNBC_data.loom\
  /mnt/04_grn/Common_files/TF_names_v_1.01.txt \
  --output /mnt/04_grn/Jan19Try/TNBC_network.tsv \
  --num_workers 36

