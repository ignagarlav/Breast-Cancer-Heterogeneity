#!/bin/bash
#SBATCH -J PropTest         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/PropTest_%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/PropTest_%A.err      
#SBATCH --time=6:00:00              # Tiempo máximo
#SBATCH --partition=short           # Partición
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo



# Definir directorios como variables de entorno
export DATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/data/"
export OUTPUT_DIR="/home/igarzonalva/Proyecto_SC_TNBC/data/GSE161529/alltypes/proportiontest"


# Crear directorios si no existen
mkdir -p $OUTPUT_DIR

# Change to the working directory
cd /home/igarzonalva/Proyecto_SC_TNBC

# Load Anaconda module
module load Anaconda3/2022.10

# Initialize Conda
source $(conda info --base)/etc/profile.d/conda.sh

# Check if the environment exists
if ! conda env list | grep -q "grn-env"; then
  echo "Error: Conda environment 'grn-env' does not exist."
  exit 1
fi

# Activate the Conda environment
conda activate grn-env

# Verify that the environment is activated
if [[ "$CONDA_DEFAULT_ENV" != "grn-env" ]]; then
  echo "Error: Failed to activate Conda environment 'grn-env'."
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

python /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/Proportions/TestProporcion.py


echo "Script completed successfully."