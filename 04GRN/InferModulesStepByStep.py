import scanpy as sc
from scipy.sparse import csc_matrix
import pandas as pd
import pickle
import argparse
from pyscenic.utils import modules_from_adjacencies

# Parse arguments
parser = argparse.ArgumentParser(description="Process adjacencies and expression matrix and generates modules.")
parser.add_argument("--data_dir", required=True, help="Path to the .h5ad file")
parser.add_argument("--adjacencies_fname", required=True, help="Path to the adjacencies file")
parser.add_argument("--modules_fname", required=True, help="Output path for modules")
parser.add_argument("--tumor_type", required=True, help="Output directory")
args = parser.parse_args()


adata = sc.read_h5ad(args.data_dir)


def fetch_exp_matrix(adata, tumor):
    adata_tumor = adata[adata.obs['subtype'] == tumor,:].copy()
    
    return(pd.DataFrame(
        data = adata_tumor.X.toarray(),
        index = adata_tumor.obs_names.values,
        columns = adata_tumor.var_names.values))


exp_matrix = fetch_exp_matrix(adata, 'ER')
adjacencies = pd.read_csv(args.adjacencies_fname, sep='\t', names = ['TF','target','importance'])


modules = list(modules_from_adjacencies(adjacencies, exp_matrix))
with open(args.modules_fname, 'wb') as f:
    pickle.dump(modules, f)

