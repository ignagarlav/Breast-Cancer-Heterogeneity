#!/bin/bash
#SBATCH --job-name=locgrn
#SBATCH --output=/home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/loc_dask_grnboost2_%A_%a.out
#SBATCH --error=/home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/loc_dask_grnboost2_%A_%a.err
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --partition=short

# Change to the working directory
cd /home/igarzonalva/Proyecto_SC_TNBC

# Initialize Conda
source /beegfs/easybuild/common/software/Anaconda3/2022.10/etc/profile.d/conda.sh

# Activate the Conda environment
conda activate pyscenic-git-env


conda info
conda list

# Check if the Conda environment exists
ENV_NAME="pyscenic-git-env"
if ! conda env list | grep -q "^${ENV_NAME}\s"; then
  echo "Error: Conda environment '${ENV_NAME}' does not exist."
  exit 1
fi

echo "Python version: $(python --version)"
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python path: $(which python)"

python --version


export DATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/Integration/adata/adata_GenAnno.h5ad"
export TF_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/grn/Common_files/TF_names_v_1.01.txt"
export NETWORK_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/grn"

echo $DATA_DIR
echo $TF_DIR
echo $NETWORK_DIR


python /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/04GRN/grn01.py 

