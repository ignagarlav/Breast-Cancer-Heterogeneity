'''https://ccc-protocols.readthedocs.io/en/latest/notebooks/ccc_python/03-Generate-Tensor.html'''

################################### Imports ###################################
import time
start_time = time.time()
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",)

import pandas as pd
import scanpy as sc
import plotnine as p9
import liana as li
import cell2cell as c2c
import decoupler as dc # needed for pathway enrichment
import warnings
warnings.filterwarnings('ignore')
from collections import defaultdict
import numpy as np
import pickle
import os 
import torch
import tensorly as tl

device = "cuda" if torch.cuda.is_available() else "cpu"
if device == "cuda":
    tl.set_backend('pytorch')
else:
    device = "cpu"
logging.info(f"Device: {device}")
logging.info(f"TensorLy backend: {tl.get_backend()}")

logging.info(f'Pandas version used: {pd.__version__}') # to load 

##################### Paths and variables ###################################
LIANA_DIR = os.getenv("LIANA_DIR")
ADATA_DIR = os.getenv('ADATA_DIR')

sample_key = 'context'
condition_key = 'subtype'
groupby = 'IGA_LT_Anno'

logging.info(f'Saving files to: {LIANA_DIR}')

############################# Data ###################################
adata = sc.read_h5ad(ADATA_DIR) 
adata.X = adata.layers['scanvi_batch_corrected_counts'].copy()

################## Ligand-Receptor Inference by Sample ####################
li.mt.rank_aggregate.by_sample(
    adata,
    groupby=groupby,
    resource_name='consensus',
    sample_key=sample_key, 
    use_raw=False,
    verbose='full', 
    n_perms=1000, 
    return_all_lrs=True, # return all LR values
    )

############# Generate 3D tensor from interactions ####################
'''We will transform the structure of the communication scores from a set of 2D-matrices for each sample into a 3D Tensor where the third dimension is sample/context.'''
tensor = li.multi.to_tensor_c2c(adata,
                                sample_key=sample_key,
                                score_key='magnitude_rank', 
                                how='outer', 
                                outer_fraction=1/3,
                                inverse_fun=lambda x: 1 - x,
                                lr_fill=np.nan, 
                                cell_fill = np.nan, 
                                )

c2c.io.export_variable_with_pickle(tensor, os.path.join(LIANA_DIR,'initial_tensorV1.pkl'))

tam = tensor.shape
miss_frac = tensor.missing_fraction()
ten_spar = tensor.sparsity_fraction()

logging.info(f'Tensor shape (number of elements in each tensor dimension: (Contexts, LR pairs, Sender cells, Receiver cells)): {tam}')
logging.info(f'Tensor missing fraction (i.e. values that are missing. In this case, missing values are combinations of contexts x LR pairs x Sender cells x Receiver cells that did not have a communication score or were missing in the dataframes.): {miss_frac}')
logging.info(f'Tensor sparsity fraction (raction of values that are a real zero (excluding the missing values): {ten_spar}')

################## Generate tensor metadata ####################
context_dict = adata.obs[[sample_key, condition_key]].drop_duplicates()
context_dict = dict(zip(context_dict[sample_key], context_dict[condition_key]))
context_dict = defaultdict(lambda: 'Unknown', context_dict) 
tensor_meta = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor,
                                                  metadata_dicts=[context_dict, None, None, None],
                                                  fill_with_order_elements=True
                                                  )
c2c.io.export_variable_with_pickle(tensor_meta, os.path.join(LIANA_DIR,'tensor_metaV1.pkl'))

################## Run tensor cell2cell analysis ####################
tensor_cell2cell = c2c.analysis.run_tensor_cell2cell_pipeline(tensor,
                                                    tensor_meta,
                                                    copy_tensor=True, 
                                                    rank=None, # elbow 
                                                    tf_optimization='robust', 
                                                    random_state=0, 
                                                    device='cuda', 
                                                    elbow_metric='error', 
                                                    smooth_elbow=False, 
                                                    upper_rank=20, 
                                                    tf_init='random',
                                                    tf_svd='numpy_svd', 
                                                    cmaps=None,
                                                    sample_col='Element',
                                                    group_col='Category', 
                                                    output_fig=False,
                                                    )

c2c.io.export_variable_with_pickle(tensor_cell2cell, os.path.join(LIANA_DIR,'cell2cell_tensorV1.pkl'))

adata.write_h5ad(os.path.join(LIANA_DIR, 'liana_adata_post_allvars.h5ad'))

end_time = time.time()
elapsed_time = end_time - start_time
logging.info(f"Elapsed time: {elapsed_time} seconds")