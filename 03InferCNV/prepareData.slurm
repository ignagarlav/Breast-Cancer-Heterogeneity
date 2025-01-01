#!/bin/bash
#SBATCH -J cnvdat         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/infercnv_dat_%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/infercnv_dat_%A.err      
#SBATCH --cpus-per-task=24    
#SBATCH --time=2:00:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo

export PYTHONPATH=/home/igarzonalva/Proyecto_SC_TNBC/SCPipeline:$PYTHONPATH

# Activar el entorno de conda
module load Anaconda3/2022.10
source activate SC310_env

cd /home/igarzonalva/Proyecto_SC_TNBC
# Correr el script ----------------------------------------
#/home/igarzonalva/.conda/envs/SC310_env/bin/python /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529_176078/SCPreprocessing.py

python /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/InferCNV/PrepareData.py
