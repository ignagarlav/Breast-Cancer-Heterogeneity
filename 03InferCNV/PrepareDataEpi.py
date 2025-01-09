'''The counts matrix can be generated using any conventional single cell transcriptome quantification pipeline, yielding a matrix of genes (rows) vs. cells (columns) containing assigned read counts'''

import scanpy as sc 
import pandas as pd
import os
import numpy as np
import anndata as ad



data_dir = os.getenv("ADATA_DIR")
adata_orig = sc.read_h5ad(os.path.join(data_dir,"adata_scanvi_predictions.h5ad"))


infer_CNV_dir = os.getenv("INFER_CNV_DIR")
os.makedirs(infer_CNV_dir, exist_ok=True)

for tumor in ["ER", "HER2", "TNBC"]:

    # Filter data for the current tumor type
    adata_tumor = adata_orig[adata_orig.obs.subtype == tumor,:].copy()

    ec_mask = adata_tumor.obs.GennAnno_ScAnvi.isin(['Epithelial', 'Cycling'])
    t_mask = adata_tumor.obs.GennAnno_ScAnvi == 'T Cells'

    epi_subset = adata_tumor[ec_mask].copy()
    t_subset = adata_tumor[t_mask].copy()

    
    split_size = epi_subset.n_obs // 2
    split_1 = epi_subset[:split_size].copy()
    split_2 = epi_subset[split_size:].copy()

    print(split_size)
        
    split_1_ref = ad.concat([split_1,t_subset])
    split_2_ref = ad.concat([split_2, t_subset])
    
    print(len(split_1_ref.obs_names), len(split_2_ref.obs_names))
        
    

    for i, adata in enumerate([split_1_ref, split_2_ref], start=1):
    
        metadata_file = os.path.join(infer_CNV_dir, f"metadata_epi_infercnv_{tumor}_{i}.tsv")
        adata.obs.GennAnno_ScAnvi.to_csv(metadata_file, sep="\t", index=True, header=False)
               
    
        counts_file = os.path.join(infer_CNV_dir,f"cnt_mat_epi_cnv_{tumor}_{i}.tsv")    
        countmatrix = adata.layers["counts"].todense().T
        print(f"Dimensions of count matrix: {countmatrix.shape}")   

        pd.DataFrame(countmatrix,
                index=adata.var_names, 
                columns=adata.obs_names).to_csv(
                    counts_file, 
                    sep="\t", 
                    header=True, 
                    index=True
                    )
        print(f"Saved counts matrix to {counts_file}")
        