import scanpy as sc
from scipy.sparse import csc_matrix
import pandas as pd
import argparse 
import loompy
import os 

# Parse arguments
parser = argparse.ArgumentParser(description="Creates loom files.")
parser.add_argument("--data_dir", required=True, help="Path to the .h5ad file")
parser.add_argument("--output_dir", required=True, help="To save loom file")
args = parser.parse_args()

def fetch_exp_matrix(adata, tumor):
    adata_tumor = adata[adata.obs['subtype'] == tumor,:].copy()
    
    return(adata_tumor.X.toarray(), adata_tumor.obs_names.values,adata_tumor.var_names.values)

adata = sc.read_h5ad(args.data_dir)

for tumor in adata.obs['subtype'].unique():
    
    loom_file = os.path.join(args.output_dir, f"{tumor}_expression_matrix.loom")
    exp_mat, cells, genes = fetch_exp_matrix(adata, tumor)
    
    expression = exp_mat.T
    # Prepare row and column attributes
    row_attrs = {
        'Gene': genes
    }
    col_attrs = {
        'CellID': cells
    }

    # Create Loom file
    loompy.create(loom_file, expression, row_attrs=row_attrs, col_attrs=col_attrs)

    print(f"Loom file saved at {loom_file}")