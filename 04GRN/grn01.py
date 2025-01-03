import sys
import scanpy as sc
import os
import numpy as np


from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2


from dask_jobqueue import SLURMCluster
from dask.distributed import Client


######################### Filter out adata #########################
tumor = sys.argv[1] 

data_dir = os.getenv("DATA_DIR")
adata = sc.read_h5ad(data_dir)

if "subtype" not in adata.obs:
    raise ValueError("Column 'subtype' not found in adata.obs.")

if tumor not in adata.obs['subtype'].unique():
    raise ValueError(f"Tumor subtype '{tumor}' not found in the data.")

adata_tumor = adata[adata.obs['subtype'] == tumor]


######################## Fetch data from filtered adata #######################
def fetch_adata(adata):
    return adata.X, adata.var_names.values, adata.obs_names.values
mat, genes, cells = fetch_adata(adata_tumor)

tf_dir = os.getenv("TF_DIR")
tf_names = load_tf_names(tf_dir)


####################### Set up Dask ############################
dashboard_port = 50000 + int(os.getenv('SLURM_ARRAY_TASK_ID', '0'))

cluster = SLURMCluster(
    # Each worker gets one process (worker) with one core (thread) in this setup
    cores=1,              
    processes=1,          
    memory='16GB',         # Memory per worker. With 4 workers, total ~64GB
    walltime='08:00:00',   
    job_extra=[
        '--partition=short',
        # Additional SLURM settings can be appended here, e.g.:
        # '--ntasks=1',
        # '--cpus-per-task=4',
        # '--mem=64G',
    ],
    scheduler_options={"dashboard_address": f":{dashboard_port}"},
)

# Scale up to 4 workers total (cores=1 * processes=1 each = 4 cores)
cluster.scale(4)

# ------------------ Step 3: Launch Dask Client ------------------ #
client = Client(cluster)
print(f"Dask client created: {client}")
print(f"Dashboard link: {client.dashboard_link}")


#################### Run grnboost2 with the Dask client ####################
try:
    network = grnboost2(
        expression_data=mat,
        gene_names=genes,
        tf_names=tf_names,
        client_or_address=client,
    )
    print(len(network))

    network_dir = os.getenv("NETWORK_DIR")
    network_file = os.path.join(network_dir, f"{tumor}_network.tsv")
    network.to_csv(network_file, sep='\t', header=False, index=False)
finally:
    if client:
        client.close()
    if cluster:
        cluster.close()
