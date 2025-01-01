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
export PYTHONPATH=/home/igarzonalva/Proyecto_SC_TNBC/SCPipeline:$PYTHONPATH

# Crear directorios si no existen
mkdir -p $OUTPUT_DIR

# Activar el entorno de conda
module load Anaconda3/2022.10
source activate grn-env

cd /home/igarzonalva/Proyecto_SC_TNBC
# Correr el script ----------------------------------------
#/home/igarzonalva/.conda/envs/SC310_env/bin/python /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529_176078/SCPreprocessing.py

python /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/Proportions/TestProporcion.py
