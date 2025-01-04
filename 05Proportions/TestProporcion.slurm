#!/bin/bash
#SBATCH -J PropTest         # Nombre del trabajo
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/PropTest_%A_%a.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/PropTest_%A_%a.err      
#SBATCH --time=24:00:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL         # Notificaciones por correo
#SBATCH --array=0-2

# Change to the working directory
cd /home/igarzonalva/Proyecto_SC_TNBC


############################### CONDA SETUP #################################

source /beegfs/easybuild/common/software/Anaconda3/2022.10/etc/profile.d/conda.sh
conda activate grn-env

conda info
conda list
echo "Python version: $(python --version)"
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python path: $(which python)"

python --version

############################ DIRECTORY SETUP ###############################

export DATA_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/Integration/adata"
export OUTPUT_DIR="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/proportiontest"

pairs=(
    "TNBC ER"
    "TNBC HER2"
    "HER2 ER"
)

# Get the current task ID (array index)
task_id=$SLURM_ARRAY_TASK_ID

# Retrieve the pair corresponding to the current task
pair=${pairs[$task_id]}

# Split the pair into input and output variables
group1=$(echo $pair | awk '{print $1}')
group2=$(echo $pair | awk '{print $2}')

echo "Running proportion test for groups $group1 and $group2"

python /home/igarzonalva/Proyecto_SC_TNBC/repos/Breast-Cancer-Heterogeneity/Proportions/TestProporcion1.py "$group1" "$group2"


