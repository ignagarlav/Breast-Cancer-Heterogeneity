'''The counts matrix can be generated using any conventional single cell transcriptome quantification pipeline, yielding a matrix of genes (rows) vs. cells (columns) containing assigned read counts'''

import scanpy as sc 
import pandas as pd
import os
import numpy as np



data_dir = os.getenv("ADATA_DIR")
adata = sc.read_h5ad(os.path.join(data_dir,"adata_GenAnno.h5ad"))


infer_CNV_dir = os.getenv("INFER_CNV_DIR")
file_idx = 1
for tumor in ["ER", "HER2", "TNBC"]:

    # Filter data for the current tumor type
    adata_tumor = adata[adata.obs['subtype'] == tumor]
    
    # Fetch metadata before splitting
    metadata = adata_tumor.obs['GenAnno']
    metadata_file = os.path.join(infer_CNV_dir, f"metadata_infercnv_{file_idx}.tsv")
    metadata.to_csv(metadata_file, sep="\t", index=True, header=False)
    print(f"Saved metadata to {metadata_file}")

    # Split the data into 3 subsets
    adata_tumor = adata_tumor[np.random.permutation(adata_tumor.obs_names)]
    split_size = adata_tumor.n_obs // 3
    subset1 = adata_tumor[:split_size].copy()
    subset2 = adata_tumor[split_size:2*split_size].copy()
    subset3 = adata_tumor[2*split_size:].copy()

    # Fetch data
    for subset in [subset1, subset2, subset3]:
        
        countmatrix = subset.layers["counts"].todense().T
        print(f"Dimensions of count matrix: {countmatrix.shape}")

        counts_file = os.path.join(infer_CNV_dir,f"counts_matrix_infercnv_{file_idx}.tsv")
        pd.DataFrame(countmatrix,
                index=subset.var_names, 
                columns=subset.obs_names).to_csv(
                    counts_file, 
                    sep="\t", 
                    header=True, 
                    index=True
                    )
        print(f"Saved counts matrix to {counts_file}")
        file_idx += 1
