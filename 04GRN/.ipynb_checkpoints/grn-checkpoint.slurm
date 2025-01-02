#!/bin/bash
#SBATCH --job-name=pyscenic_grn          # Job name
#SBATCH -o /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/grn_pyscenic_infer_%A.out       
#SBATCH -e /home/igarzonalva/Proyecto_SC_TNBC/SlurmOutput/GSE161529/grn_pyscenic_infer_%A.err    
#SBATCH --ntasks=1                       # Number of tasks (processes)
#SBATCH --mem=128G                        # Memory per node
#SBATCH --time=72:00:00                  # Time limit hrs:min:sec
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

INPUT_FILE="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/grn/expression_adata_$SLURM_ARRAY_TASK_ID.h5ad"

OUTPUT_FILE="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/grn/adyacencias_$SLURM_ARRAY_TASK_ID.tsv"

TF_FILE="/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/grn/Common_files/TF_names_v_1.01.txt"

echo "Processing file $INPUT_FILE"

# Run the pyscenic grn command
pyscenic grn --sparse $INPUT_FILE $TF_FILE -o $OUTPUT_FILE

echo "Adyacencies file saved in $OUTPUT_FILE"
# Deactivate the virtual environment (if used)
source deactivate

echo "Done."