#!/bin/bash
#SBATCH -J myelarches         
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/myelarches%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/myelarches%A.err      
#SBATCH --time=20:00:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   
#SBATCH --mail-type=END,FAIL         

# Change to the working directory
cd /home/igarzonalva/Proyecto_SC_TNBC

# Initialize Conda
source /beegfs/easybuild/common/software/Anaconda3/2022.10/etc/profile.d/conda.sh

# Activate the Conda environment
conda activate scarches

conda info
conda list

# Check if the Conda environment exists
ENV_NAME="scarches"
if ! conda env list | grep -q "^${ENV_NAME}\s"; then
  echo "Error: Conda environment '${ENV_NAME}' does not exist."
  exit 1
fi

echo "Python version: $(python --version)"
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python path: $(which python)"

python --version


export ADATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/LabelTransfer/adatas/Myeloid"
export MODEL_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/LabelTransfer/models/Myeloid"


python /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/03LabelTransfer/MyeloidTransfer.py

