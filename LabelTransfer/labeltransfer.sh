#!/bin/bash
#SBATCH -J LTscarches         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/LTscarches_%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/LTscarches_%A.err      
#SBATCH --time=24:00:00              # Tiempo máximo
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


export ADATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/LabelTransfer/adatas"
export MODEL_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/LabelTransfer/models"
export QUERY_ADATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/Integration/adata/adata_scanvi_predictions.h5ad"

python /home/igarzonalva/Proyecto_SC_TNBC/repos/LabelTransfer/labeltransfer.py

