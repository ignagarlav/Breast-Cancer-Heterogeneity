#!/bin/bash
#SBATCH -J SCPrepro         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529_176078/scPP_%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529_176078/scPP_%A.err      
#SBATCH --time=4:00:00              # Tiempo máximo
#SBATCH --partition=short           # Partición
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo



# Definir directorios como variables de entorno
export LOG_DIR="/home/igarzonalva/Proyecto_SC_TNBC/output/GSE162929_176078/tnbc"
export DATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/data/GSE162929_176078/raw"
export OUTPUT_DIR="/home/igarzonalva/Proyecto_SC_TNBC/data/GSE161529/alltypes/without_normal"
export PYTHONPATH=/home/igarzonalva/Proyecto_SC_TNBC/SCPipeline:$PYTHONPATH


# Crear directorios si no existen
mkdir -p $LOG_DIR
mkdir -p $DATA_DIR
mkdir -p $OUTPUT_DIR

# Activar el entorno de conda
module load Anaconda3/2022.10
source activate SC310_env

cd /home/igarzonalva/Proyecto_SC_TNBC
# Correr el script ----------------------------------------
#/home/igarzonalva/.conda/envs/SC310_env/bin/python /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529_176078/SCPreprocessing.py

python /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/SCPreprocessingNoNormal.py
