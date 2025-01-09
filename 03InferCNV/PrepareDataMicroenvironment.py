'''The counts matrix can be generated using any conventional single cell transcriptome quantification pipeline, yielding a matrix of genes (rows) vs. cells (columns) containing assigned read counts'''

import scanpy as sc 
import pandas as pd
import os
import numpy as np



data_dir = os.getenv("ADATA_DIR")
adata_orig = sc.read_h5ad(os.path.join(data_dir,"adata_scanvi_predictions.h5ad"))

infer_CNV_dir = os.getenv("INFER_CNV_DIR")
os.makedirs(infer_CNV_dir, exist_ok=True)

for tumor in ["ER", "HER2", "TNBC"]:

    # Filter data for the current tumor type
    adata_tumor = adata_orig[adata_orig.obs.subtype == tumor]
    
    # Select only microenvironment cell populations
    mask = adata_tumor.obs['GennAnno_ScAnvi'].isin(['Cycling', 'Epithelial'])
    subset = adata_tumor[~mask].copy()
   
        
    # Fetch metadata 
    metadata = subset.obs.GennAnno_ScAnvi
    metadata_file = os.path.join(infer_CNV_dir, f"meta_mcroenv_cnv_{tumor}.tsv")
    metadata.to_csv(metadata_file, sep="\t", index=True, header=False)
    print(f"Saved metadata to {metadata_file}")
    
    
    countmatrix = subset.layers["counts"].todense().T
    print(f"Dimensions of count matrix: {countmatrix.shape}")

    counts_file = os.path.join(infer_CNV_dir,f"cnt_mt_cnv_mcoenv{tumor}.tsv")
    pd.DataFrame(countmatrix,
            index=subset.var_names, 
            columns=subset.obs_names).to_csv(
                counts_file, 
                sep="\t", 
                header=True, 
                index=True
                )
    print(f"Saved counts matrix to {counts_file}")
