from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from distributed import LocalCluster, Client
import scanpy as sc
import numpy as np
from dask_jobqueue import SLURMCluster
import os
import logging

# Configure the logging level and format
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

current_node = os.environ.get("SLURMD_NODENAME")
logging.info(f"Current node:{current_node}")

NETWORK_DIR = os.getenv("NETWORK_DIR")
TUMOR_TYPE = os.getenv("TUMOR_TYPE")
ADATA_DIR = os.getenv("ADATA_DIR")
if __name__ == '__main__':
      
    logging.info('loading data')  
    adata = sc.read_h5ad(ADATA_DIR)
    adata_sub = adata[adata.obs['subtype'] == TUMOR_TYPE,:].copy()
    
    
    logging.info('Creating expression matrix')
    mat = adata_sub.X.toarray()
    logging.info(f'expression matrix shape: {mat.shape}')
    genes = adata_sub.var_names.values
    logging.info(f'Number of genes: {(len(genes))}')

    assert mat.shape[1] == len(genes), "Number of genes does not match the number of columns in the expression matrix"
    

    logging.info('loading TF names')
    tf_names = load_tf_names("/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/04_grn/Common_files/TF_names_v_1.01.txt")

    
    

    logging.info('setting up cluster')
    cluster = SLURMCluster(queue = "short", cores=6,     processes=1, memory="16GB", walltime="24:00:00", scheduler_options={"host": current_node})
    logging.info(f'cluster: {cluster}')    
    
    num_workers = 8
    logging.info(f'Setting up {num_workers} workers')
    cluster.scale(num_workers)
    
    logging.info(f'Waiting for the workers to be ready...')
    # Block until the cluster reports exactly num_workers workers
    custom_client = Client(cluster)
    custom_client.wait_for_workers(num_workers)

    logging.info('Workers are ready. Resuming the workflow...')
    logging.info(f'client: {custom_client}')
    logging.info(f"Dask Scheduler is listening on: {cluster.scheduler_address}")


    network = grnboost2(expression_data=mat,
                            gene_names=genes,
                            tf_names=tf_names,
                            client_or_address=custom_client,  
                            seed=666, verbose=True)



    # close the Client and LocalCluster after use
    custom_client.close()
    cluster.close()
    
    network_name = os.path.join(NETWORK_DIR,f'network_{TUMOR_TYPE}.tsv')
    network.to_csv(network_name, sep='\t', index=False, header=False)