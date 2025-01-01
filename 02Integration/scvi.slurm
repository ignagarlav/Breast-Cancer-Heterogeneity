#!/bin/bash
#SBATCH -J scvi-in         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/scvi_in_%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/scvi_in_%A.err      
#SBATCH --time=24:00:00              # Tiempo máximo
#SBATCH --partition=medium           # Partición
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo

export DATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/Preprocessing"
export MODEL_DIR="/home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/Integration/models"
export OUTPUT_DIR="/home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/Integration"
export PYTHONPATH=/home/igarzonalva/Proyecto_SC_TNBC/SCPipeline:$PYTHONPATH



#export PYTHONPATH=/home/igarzonalva/Proyecto_SC_TNBC/SCPipeline:$PYTHONPATH

export PYTHON_BIN="/home/igarzonalva/.conda/envs/scvi-env/bin/python"


# Activar el entorno de conda
module load Anaconda3/2022.10
source activate scvi-env
cd /home/igarzonalva/Proyecto_SC_TNBC
# Correr el script ----------------------------------------
$PYTHON_BIN --version
$PYTHON_BIN /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/Integration/IntegrationScvi.py

