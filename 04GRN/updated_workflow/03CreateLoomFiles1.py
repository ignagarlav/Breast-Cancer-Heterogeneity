import scanpy as sc
import argparse 
import loompy
import os 
import numpy as np

# Parse arguments
parser = argparse.ArgumentParser(description="Creates loom files.")
parser.add_argument("--data_dir", required=True, help="Path to the .h5ad file")
parser.add_argument("--output_dir", required=True, help="To save loom file")
args = parser.parse_args()

adata = sc.read_h5ad(args.data_dir)

for tumor in adata.obs['subtype'].cat.categories:
    
    adata_tumor = adata[adata.obs['subtype'] == tumor,:].copy()

    floom_file = os.path.join(args.output_dir, f"{tumor}_expression_matrix.loom")
    
    row_attrs = {
        'Gene': np.array(adata_tumor.var_names)
    }
    col_attrs = {
        'CellID': np.array(adata_tumor.obs_names) ,
    }

    # Create Loom file
    loompy.create(floom_file, adata_tumor.X.transpose(), row_attrs, col_attrs)

    print(f"Loom file saved at {floom_file}")

