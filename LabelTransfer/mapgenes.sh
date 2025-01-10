#!/bin/bash
#SBATCH -J mapgen         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/mapgen_myeloid%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/mapgen_myeloid%A.err      
#SBATCH --time=1:00:00              # Tiempo máximo
#SBATCH --partition=short           # Partición
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo

# Change to the working directory
cd /home/igarzonalva/Proyecto_SC_TNBC

# Initialize Conda
source /beegfs/easybuild/common/software/Anaconda3/2022.10/etc/profile.d/conda.sh

# Activate the Conda environment
conda activate scvi-env

conda info
conda list

# Check if the Conda environment exists
ENV_NAME="scvi-env"
if ! conda env list | grep -q "^${ENV_NAME}\s"; then
  echo "Error: Conda environment '${ENV_NAME}' does not exist."
  exit 1
fi

echo "Python version: $(python --version)"
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python path: $(which python)"

python --version


export ADATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/LabelTransfer/adatas"
export QUERY_ADATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/Integration/adata/adata_scanvi_predictions.h5ad"

python /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/LabelTransfer/MapGenes.py

