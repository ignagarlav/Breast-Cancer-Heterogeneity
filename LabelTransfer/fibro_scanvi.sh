#!/bin/bash
#SBATCH -J fibanvi         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/fibanvi%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/fibanvi%A.err      
#SBATCH --time=6:00:00              # Tiempo máximo
#SBATCH --partition=short           # Partición
#SBATCH --mem=16G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo

# Change to the working directory
cd /home/igarzonalva/Proyecto_SC_TNBC
# Singularity image path
SIF_PATH="/home/igarzonalva/Proyecto_SC_TNBC/singularity_images/iga_scvi_env_v1.0.1.sif"

# Set environment variables that the script needs
export DATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/LabelTransfer/adatas/Fibroblast"
export MODEL_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/LabelTransfer/models/fibroblast"

# Run Python script with Singularity
singularity exec --bind /home/igarzonalva/Proyecto_SC_TNBC:/home/igarzonalva/Proyecto_SC_TNBC "${SIF_PATH}" python /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/LabelTransfer/Fibro_scanvi.py
