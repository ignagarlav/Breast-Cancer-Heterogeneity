#!/bin/bash
#SBATCH -J DatcnvEpi         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/DatcnvEpi%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/DatcnvEpi%A.err      
#SBATCH --time=6:00:00              
#SBATCH --partition=medium           
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo

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


export ADATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/Integration/adata"
export INFER_CNV_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV"


python /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/03InferCNV/PrepareDataEpi.py
