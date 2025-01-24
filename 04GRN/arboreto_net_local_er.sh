#!/bin/bash
#SBATCH -J Arboreto        
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/arboreto%A_%a.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/arboreto%A._%a.err      
#SBATCH --time=24:00:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   
#SBATCH --mail-type=END,FAIL         


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

export TUMOR_TYPE="ER"

export NETWORK_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/Networks"
export ADATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/02_Integration/adata/adata_scanvi_cuda_refinement.h5ad"

python /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/04GRN/arboreto_net_localer.py