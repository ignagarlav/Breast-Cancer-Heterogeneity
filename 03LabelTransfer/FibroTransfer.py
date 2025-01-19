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


CONDITION_KEY = "condition"
CELL_TYPE_KEY = "CAF"
UNLABELED_CATEGORY = "unknown"

############################### Data Loading ##############################    
source_adata_fname = os.path.join(ADATA_DIR,'ref_adata.h5ad')
source_adata = sc.read_h5ad(source_adata_fname)

target_adata_fname = os.path.join(ADATA_DIR,'target_adata.h5ad')
target_adata = sc.read_h5ad(target_adata_fname)

source_adata = remove_sparsity(source_adata)
target_adata = remove_sparsity(target_adata)

################ Create and train ref model on fully labelled data ###########
sca.models.SCVI.setup_anndata(source_adata, 
                              batch_key=CONDITION_KEY, 
                              labels_key=CELL_TYPE_KEY)

vae = sca.models.SCVI(
    source_adata,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)

vae.train()

scanvae = sca.models.SCANVI.from_scvi_model(vae, 
unlabeled_category = UNLABELED_CATEGORY,labels_key=CELL_TYPE_KEY,)

print("Labelled Indices: ", len(scanvae._labeled_indices))
print("Unlabelled Indices: ", len(scanvae._unlabeled_indices))

scanvae.train(max_epochs=20,n_samples_per_label=100)

# Save the reference model
ref_path = os.path.join(MODEL_DIR,'ref_model')
scanvae.save(ref_path, overwrite=True)


################# Perform surgery on reference model and train #################
model = sca.models.SCANVI.load_query_data(
    target_adata,
    ref_path,
    freeze_dropout=True,
)
model._unlabeled_indices = np.arange(target_adata.n_obs)
model._labeled_indices = []
print("Labelled Indices: ", len(model._labeled_indices))
print("Unlabelled Indices: ", len(model._unlabeled_indices))

model.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0), check_val_every_n_epoch=10,)


# Save the surgery model 
surgery_path = os.path.join(MODEL_DIR,'surgery_model')
model.save(surgery_path, overwrite=True)
