#!/bin/bash
#SBATCH -J scanvi         
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/pruebaPython_%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/pruebaPython_%A.err      
#SBATCH --time=0:30:00              
#SBATCH --partition=short           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-user=igarzonalva@alumni.unav.es   
#SBATCH --mail-type=END,FAIL        


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

export EXPORT_1="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/"
export EXPORT_2="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/"

echo $EXPORT_1
echo $EXPORT_2

# [Add your Python script execution here]
# Example:
# python /path/to/script.py

# Deactivate the Conda environment
conda deactivate

# End of script
echo "End of script"

