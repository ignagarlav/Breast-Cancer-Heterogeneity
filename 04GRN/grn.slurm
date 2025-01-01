#!/bin/bash
#SBATCH --job-name=pyscenic_grn          # Job name
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/grn_pyscenic_infer_%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/grn_pyscenic_infer_%A.err    
#SBATCH --ntasks=1                       # Number of tasks (processes)

#SBATCH --mem=128G                        # Memory per node
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --partition=medium       # Replace with your partition name
#SBATCH --mail-user=igarzonalva@alumni.unav.es   # Dirección de correo
#SBATCH --mail-type=END,FAIL 

cd /home/igarzonalva/Proyecto_SC_TNBC

module load Anaconda3/2022.10
# source /home/igarzonalva/.conda/envs/pyscenic-git-env/bin/activate
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pyscenic-git-env

# Print Python version and environment information for debugging
echo "Python version: $(python --version)"
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python path: $(which python)"

# Run the pyscenic grn command
pyscenic grn --sparse /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/grn/expression_adata.h5ad /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/grn/grn_common_files/TF_names_v_1.01.txt -o /home/igarzonalva/Proyecto_SC_TNBC/scripts/GSE161529/alltypes/grn/adyacencias.tsv

# Deactivate the virtual environment (if used)
deactivate