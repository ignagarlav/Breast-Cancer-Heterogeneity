#!/bin/bash
#SBATCH -J cnv-dat         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/infercnv_dat_%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/infercnv_dat_%A.err      
#SBATCH --time=6:00:00              
#SBATCH --partition=medium           
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo

# Change to the working directory
cd /home/igarzonalva/Proyecto_SC_TNBC

# Load Anaconda module
module load Anaconda3/2022.10

# Initialize Conda
source $(conda info --base)/etc/profile.d/conda.sh

conda activate pyscenic-git-env

# Verify that the environment is activated
if [[ "$CONDA_DEFAULT_ENV" != "pyscenic-git-env" ]]; then
  echo "Error: Failed to activate Conda environment 'pyscenic-git-env'."
  exit 1
fi

# Check if required tools or packages are available
if ! command -v python &> /dev/null; then
  echo "Error: Python is not available in the environment."
  exit 1
fi

# Print Python version and environment information for debugging
echo "Python version: $(python --version)"
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python path: $(which python)"

export ADATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/Integration/adata"
export INFER_CNV_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV"

export PYTHON_BIN="/home/igarzonalva/.conda/envs/pyscenic-git-env/bin/python"
$PYTHON_BIN --version

$PYTHON_BIN /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/03InferCNV/PrepareData.py

conda deactivate 