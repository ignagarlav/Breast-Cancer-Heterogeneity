#!/bin/bash
#SBATCH --job-name=ad2loom
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/ad2loom%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/ad2loom%A.err  
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --time=05:00:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1

singularity exec \
  -B /home/igarzonalva/Proyecto_SC_TNBC:/home/igarzonalva/Proyecto_SC_TNBC \
  /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/grn/aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif \
  python /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/04GRN/updated_workflow/03CreateLoomFiles1.py \
    --data_dir /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/03_LabelTransfer/adatas/general_adata_postLT.h5ad \
    --output_dir /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/updated_workflow/loom_files \


