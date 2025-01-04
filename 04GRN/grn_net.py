import sys
import scanpy as sc
import os
import numpy as np
from scipy.sparse import csc_matrix
import pandas as pd

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from dask.distributed import Client, LocalCluster
import subprocess as sp


def main():
    
    ###################### Load adata ############################
    data_dir = os.getenv("DATA_DIR")
    adata = sc.read_h5ad(data_dir)

    def fetch_adata(adata):
        return csc_matrix(adata.X).toarray(), adata.var_names.values, adata.obs_names.values
    
    ###################### Load TF names ############################
    tf_dir = os.getenv("TF_DIR")
    tf_names = load_tf_names(tf_dir)
    
    ######################### Filter out adata #########################
    for tumor in adata.obs['subtype'].unique():

        print(f"Processing tumor subtype: {tumor}")
        if "subtype" not in adata.obs:
            raise ValueError("Column 'subtype' not found in adata.obs.")

        if tumor not in adata.obs['subtype'].unique():
            raise ValueError(f"Tumor subtype '{tumor}' not found in the data.")

        adata_tumor = adata[adata.obs['subtype'] == tumor].copy()


        ######################## Fetch data from filtered adata #######################
        mat, genes, cells = fetch_adata(adata_tumor)


        ####################### Set up Dask ############################
        #portdash = 40748
        #cluster = SLURMCluster(queue = "short", cores=8, processes=1, 
                            #memory="16GB", walltime="00:30:00",
                            #scheduler_options={"dashboard_address": f":{portdash}"})

        # Scale up to 4 workers total (cores=8 * processes=1 each = 8 cores)
        #cluster.scale(4) 
        # Ask for 4 workers (4 slurm jobs with 8 cores each = 32 cores total)
        # Ask for 16GB of memory per worker (4 workers * 16GB = 64GB total)

        cluster = LocalCluster(n_workers=8, threads_per_worker=1)
        client = Client(cluster)


        # # ------------------ Step 3: Launch Dask Client ------------------ #

        #client = Client(cluster)
        print(f"Dask client created: {client}")
        print(f"Dashboard link: {client.dashboard_link}")


        #################### Run grnboost2 with the Dask client ####################
        try:
            network = grnboost2(
                            expression_data=mat,
                            gene_names=genes,
                            tf_names=tf_names,
                            client_or_address=client)
            print(len(network))

            network_dir = os.getenv("NETWORK_DIR")
            network_file = os.path.join(network_dir, f"{tumor}_slurm_network.tsv")
            network.to_csv(network_file, sep='\t', header=False, index=False)

        finally:
            if client:
                client.close()
            if cluster:
                cluster.close()


if __name__ == '__main__':
    main()