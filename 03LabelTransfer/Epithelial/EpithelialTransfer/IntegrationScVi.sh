#!/bin/bash
#SBATCH -J CudScVi         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/ScViCuda%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/ScViCuda%A.err      
#SBATCH --time=10:00:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1 
#SBATCH --nodelist=nodo10
#SBATCH --mem=32G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo

# Change to the working directory
cd /home/igarzonalva/Proyecto_SC_TNBC

# Initialize Conda
source /beegfs/easybuild/common/software/Anaconda3/2022.10/etc/profile.d/conda.sh

# Activate the Conda environment
conda activate scvi_cuda_env3

conda info
conda list

# Check if the Conda environment exists
ENV_NAME="scvi_cuda_env3"
if ! conda env list | grep -q "^${ENV_NAME}\s"; then
  echo "Error: Conda environment '${ENV_NAME}' does not exist."
  exit 1
fi

echo "Python version: $(python --version)"
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python path: $(which python)"

python --version


export DATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/02_Integration/adata"
export MODEL_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/02_Integration/models"

python /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/02Integration/IntegrationScvi.py

