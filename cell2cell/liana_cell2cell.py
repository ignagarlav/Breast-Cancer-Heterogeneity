################################### Imports ###################################
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
# Set the backend to PyTorch
use_gpu = True

if use_gpu:
    import torch
    import tensorly as tl

    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        tl.set_backend('pytorch')
else:
    device = "cpu"


LIANA_DIR = os.getenv("LIANA_DIR")
################################### Data ###################################
adata_fname = os.path.join(LIANA_DIR, 'liana_adata.h5ad')
adata_analysis = sc.read_h5ad(adata_fname) 

sample_key = 'context'
condition_key = 'subtype'
groupby = 'post_lt_anno'


###################### Ligand-Receptor Inference by Sample ####################
li.mt.rank_aggregate.by_sample(
    adata_analysis,
    groupby=groupby,
    resource_name='consensus',
    sample_key=sample_key, # sample key by which we which to loop
    use_raw=False,
    verbose='full', # use 'full' to show all verbose information
    n_perms=1000, 
    return_all_lrs=True, # return all LR values
    )


########################## Building a Tensor ##################################

tensor = li.multi.to_tensor_c2c(adata_analysis,
                                sample_key=sample_key,
                                score_key='magnitude_rank', 
                                how='outer_cells', # union
                                non_negative = True,
                                inverse_fun=lambda x: 1 - x,
                                non_expressed_fill=0,
                                outer_fraction=1/6,
                                lr_fill=np.nan, 
                                cell_fill = np.nan, # 
                                )
""" 
Fraction of samples as threshold to include cells and LR pairs.
Como tengo 4 grupos pero un poco descompensados y quiero mantener interacciones que son propias de un solo grupo (por ejemplo una interacción que sea de sólo TNBC o sólo de HER2), pondré 1/6, teniendo en cuenta que algunas muestras dentro de un grupo pueden ser outliers. """



c2c.io.export_variable_with_pickle(tensor, os.path.join(LIANA_DIR,'initial_tensorV1.pkl'))


context_dict = adata_analysis.obs[[sample_key, condition_key]].drop_duplicates()
context_dict = dict(zip(context_dict[sample_key], context_dict[condition_key]))
context_dict = defaultdict(lambda: 'Unknown', context_dict) # add unkown to missing keys

c2c.io.export_variable_with_pickle(context_dict, os.path.join(LIANA_DIR,'context_dictV1.pkl'))


tensor_meta = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor,
                                                  metadata_dicts=[context_dict, None, None, None],
                                                  fill_with_order_elements=True
                                                  )

c2c.io.export_variable_with_pickle(tensor_meta, os.path.join(LIANA_DIR,'tensor_metaV1.pkl'))

tensor_cell2cell = c2c.analysis.run_tensor_cell2cell_pipeline(tensor,
                                                    tensor_meta,
                                                    copy_tensor=True, # Whether to output a new tensor or modifying the original
                                                    rank=None, # Number of factors to perform the factorization. If None, it is automatically determined by an elbow analysis. Here, it was precomuputed.
                                                    tf_optimization='robust', # To define how robust we want the analysis to be.
                                                    random_state=0, # Random seed for reproducibility
                                                    device='cuda', # Device to use. If using GPU and PyTorch, use 'cuda'. For CPU use 'cpu'
                                                    elbow_metric='error', # Metric to use in the elbow analysis.
                                                    smooth_elbow=False, # Whether smoothing the metric of the elbow analysis.
                                                    upper_rank=20, # Max number of factors to try in the elbow analysis
                                                    tf_init='random', # Initialization method of the tensor factorization
                                                    tf_svd='numpy_svd', # Type of SVD to use if the initialization is 'svd'
                                                    cmaps=None, # Color palettes to use in color each of the dimensions. Must be a list of palettes.
                                                    sample_col='Element', # Columns containing the elements in the tensor metadata
                                                    group_col='Category', # Columns containing the major groups in the tensor metadata
                                                    output_fig=False, # Whether to output the figures. If False, figures won't be saved a files if a folder was passed in output_folder.
                                                    )

c2c.io.export_variable_with_pickle(tensor_cell2cell, os.path.join(LIANA_DIR,'cell2cell_tensorV1.pkl'))

adata_analysis.write_h5ad(os.path.join(LIANA_DIR, 'liana_adata_postV1.h5ad'))