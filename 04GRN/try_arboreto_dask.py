import scanpy as sc
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
from scipy.sparse import csc_matrix
import os 
from dask.distributed import Client
import pandas as pd

if __name__ == '__main__':

    DATA_DIR = "/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/02_Integration/adata/adata_scanvi_cuda_refinement.h5ad"
    TF_FNAME = "/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/Common_files/TF_names_v_1.01.txt"
    TUMOR = 'TNBC'

    NETWORK_FAME = f'/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/{TUMOR}_network.csv'


    adata = sc.read_h5ad(DATA_DIR)
    adata_raw = adata.raw.to_adata()
    adata_tumor = adata_raw[adata_raw.obs.subtype == TUMOR,:].copy()

    expression_mat = pd.DataFrame(adata_tumor.X.toarray(), index=adata_tumor.obs_names, columns=adata_tumor.var_names)

    
    tf_names = load_tf_names(TF_FNAME)
    client = Client('tcp://159.237.145.53:46800')
    network = grnboost2(
                expression_data=expression_mat,
                tf_names=tf_names,
                client_or_address=client,
                verbose = True,)


    network.to_csv(NETWORK_FAME, sep='\t', header=False, index=False)

