#!/bin/bash
#SBATCH --job-name=loom
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/loom%A_%a.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/lomm%A_%a.err  
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --time=00:30:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1


singularity exec \
  -B /home/igarzonalva/Proyecto_SC_TNBC:/mnt \
  /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/grn/aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
  python /mnt/repos/Breast-Cancer-Heterogeneity/04GRN/CreateLoomFiles.py \
    --data_dir /mnt/GSE161529/Integration/adata/adata_GenAnno.h5ad \
    --output_dir /mnt/GSE161529/grn/exp_matrices/ \
