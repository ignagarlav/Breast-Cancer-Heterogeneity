#!/bin/bash
#SBATCH -J Arboreto        
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/arboreto%A_%a.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/arboreto%A._%a.err      
#SBATCH --time=24:00:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   
#SBATCH --mail-type=END,FAIL         
#SBATCH --array=0-2%1

# Change to the working directory
cd /home/igarzonalva/Proyecto_SC_TNBC

# Initialize Conda
source /beegfs/easybuild/common/software/Anaconda3/2022.10/etc/profile.d/conda.sh

# Activate the Conda environment
conda activate arboreto-env
conda info
conda list

# Check if the Conda environment exists
ENV_NAME="arboreto-env"
if ! conda env list | grep -q "^${ENV_NAME}\s"; then
  echo "Error: Conda environment '${ENV_NAME}' does not exist."
  exit 1
fi

echo "Python version: $(python --version)"
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python path: $(which python)"

python --version

TUMORS=("HER2" "TNBC_BRCA" "ER")

TUMOR=${TUMORS[$SLURM_ARRAY_TASK_ID]}
echo "Running array job index: $SLURM_ARRAY_TASK_ID"
echo "Processing tumor: $TUMOR"

export TUMOR_TYPE="$TUMOR"

export NETWORK_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/Networks"
export ADATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/02_Integration/adata/adata_scanvi_cuda_refinement.h5ad"

python /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/04GRN/arboreto_net_local.py