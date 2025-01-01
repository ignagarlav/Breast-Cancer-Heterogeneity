'''The counts matrix can be generated using any conventional single cell transcriptome quantification pipeline, yielding a matrix of genes (rows) vs. cells (columns) containing assigned read counts'''

import scanpy as sc 
import pandas as pd
import os

infer_CNV_dir = "/home/igarzonalva/Proyecto_SC_TNBC/data/GSE161529/alltypes/InferCNV"
data_dir = "/home/igarzonalva/Proyecto_SC_TNBC/data/GSE161529/alltypes"

adata = sc.read_h5ad(os.path.join(data_dir,"adata_GenAnno.h5ad"))


for i, tumor in enumerate(["ER", "HER2", "TNBC"], start=1):

    countmatrix = adata[adata.obs['subtype'] == tumor].layers["counts"].todense().T
    print(f"Dimensions of count matrix: {countmatrix.shape}")

    counts_file = os.path.join(infer_CNV_dir,f"counts_matrix_infercnv_{i}.tsv")
    pd.DataFrame(countmatrix,
                index=adata[adata.obs['subtype'] == tumor].var_names, 
                columns=adata[adata.obs['subtype'] == tumor].obs_names).to_csv(
                    counts_file, 
                    sep="\t", 
                    header=True, 
                    index=True
                    )
    print(f"Saved counts matrix to {counts_file}")

    metadata = adata[adata.obs['subtype'] == tumor].obs['GenAnno']
    metadata_file = os.path.join(infer_CNV_dir, f"metadata_infercnv_{i}.tsv")

    
    metadata.to_csv(metadata_file, sep="\t", index=True, header=False)
    print(f"Saved metadata to {metadata_file}")