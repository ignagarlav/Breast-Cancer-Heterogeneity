#!/bin/bash
#SBATCH -J Liana         
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/Liana%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/Liana%A.err      
#SBATCH --time=23:00:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1 
#SBATCH --nodelist=nodo10
#SBATCH --mem=50G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   
#SBATCH --mail-type=END,FAIL         

# Change to the working directory
cd /home/igarzonalva/Proyecto_SC_TNBC

# Initialize Conda
source /beegfs/easybuild/common/software/Anaconda3/2022.10/etc/profile.d/conda.sh

# Activate the Conda environment
conda activate cell2cell_gpu_env

conda info
conda list

# Check if the Conda environment exists
ENV_NAME="cell2cell_gpu_env"
if ! conda env list | grep -q "^${ENV_NAME}\s"; then
  echo "Error: Conda environment '${ENV_NAME}' does not exist."
  exit 1
fi

echo "Python version: $(python --version)"
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python path: $(which python)"

python --version


export LIANA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/liana_data"
export ADATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/03_LabelTransfer/adatas/general_adata_postLT_allvars"

python /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/cell2cell/liana_cell2cell.py

