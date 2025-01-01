import pandas as pd
import numpy as np
import os
import scanpy as sc
import loompy

data_dir = "/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/data/BC_All_WithoutNormal"
adata = sc.read_h5ad(os.path.join(data_dir,"adata_integrated.h5ad"))

grn_dir = "/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/CyclingEpithelial"

expression_matrix = adata.raw.X
cell_names = adata.obs_names.to_series()
gene_names = adata.raw.var_names.to_series()

new_adata = sc.AnnData(
    X=expression_matrix,  # Raw expression matrix
    obs=pd.DataFrame(index=cell_names),  # Cell names as observations
    var=pd.DataFrame(index=gene_names)   # Gene names as variables
)

new_adata.write(os.path.join(grn_dir, "expression_data.h5ad"))

