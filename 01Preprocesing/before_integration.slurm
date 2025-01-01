#!/bin/bash
#SBATCH -J scvi-prep         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/scvi_prep_%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/scvi_prep_%A.err      
#SBATCH --time=2:00:00              # Tiempo máximo
#SBATCH --partition=short           # Partición
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo


export DATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/Preprocessing"
#export PYTHONPATH=/home/igarzonalva/Proyecto_SC_TNBC/SCPipeline:$PYTHONPATH

export PYTHON_BIN="/home/igarzonalva/.conda/envs/scvi-env/bin/python"


# Activar el entorno de conda
module load Anaconda3/2022.10
source activate scvi-env
cd /home/igarzonalva/Proyecto_SC_TNBC
# Correr el script ----------------------------------------
$PYTHON_BIN --version
$PYTHON_BIN /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/Preprocessing/BeforeIntegration.py
