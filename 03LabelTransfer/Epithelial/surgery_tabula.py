import scanpy as sc
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import numpy as np
import gdown
import os 


################################## Environment ################################
# Get environment variables
ADATA_DIR = os.getenv('ADATA_DIR')
MODEL_DIR = os.getenv('MODEL_DIR')



target_adata_fname = os.path.join(ADATA_DIR,'target_adata.h5ad')
target_adata = sc.read_h5ad(target_adata_fname)

ref_path = os.path.join(MODEL_DIR,'ref_model')


################# Perform surgery on reference model and train #################
model = sca.models.SCANVI.load_query_data(
    target_adata,
    ref_path,
)
model._unlabeled_indices = np.arange(target_adata.n_obs)
model._labeled_indices = []
print("Labelled Indices: ", len(model._labeled_indices))
print("Unlabelled Indices: ", len(model._unlabeled_indices))

model.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0), check_val_every_n_epoch=10,)


# Save the surgery model 
surgery_path = os.path.join(MODEL_DIR,'surgery_model')
model.save(surgery_path, overwrite=True)
