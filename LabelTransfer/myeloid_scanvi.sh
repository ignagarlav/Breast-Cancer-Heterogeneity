#!/bin/bash
#SBATCH -J myedScanvi         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/myedScanvi%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/myedScanvi%A.err      
#SBATCH --time=40:00:00              # Tiempo máximo
#SBATCH --partition=medium           # Partición
#SBATCH --mem=64G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo

# Change to the working directory
cd /home/igarzonalva/Proyecto_SC_TNBC
# Singularity image path
SIF_PATH="/home/igarzonalva/Proyecto_SC_TNBC/singularity_images/iga_scvi_env_v1.0.1.sif"

# Set environment variables that the script needs
export DATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/LabelTransfer/adatas"
export MODEL_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/LabelTransfer/models/myeloid"

# Run Python script with Singularity
singularity exec --bind /home/igarzonalva/Proyecto_SC_TNBC:/home/igarzonalva/Proyecto_SC_TNBC "${SIF_PATH}" python /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/LabelTransfer/myeloid_scanvi.py
