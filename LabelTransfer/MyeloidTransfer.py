import scanpy as sc
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import numpy as np
import gdown
import os 
import mygene

# Get environment variables
adata_dir = os.getenv('ADATA_DIR')
model_dir = os.getenv('MODEL_DIR')


condition_key = "harm_study"
cell_type_key = "author_cell_type"


## Load reference dataset
reference_adata_fname = os.path.join(adata_dir,'MyeloidAdata_mapped.h5ad')
raw_adata = sc.read_h5ad(reference_adata_fname)


### Load the query dataset
query_adata_fname = os.path.join(adata_dir,'ScanviPredictionsAdata_mapped_to_myeloid.h5ad')
target_adata = sc.read_h5ad(query_adata_fname)

## Subset query to only myeloid cells
target_adata = target_adata[target_adata.obs.GennAnno_ScAnvi.isin(['Myeloid','Mast cells']),:].copy()


# Compute HVGs on raw data of reference dataset
sc.pp.highly_variable_genes(raw_adata, flavor='seurat_v3', n_top_genes=3000)

# Slice raw and target adatas to keep HVGs
source_adata = raw_adata[:, raw_adata.var['highly_variable']].copy()
target_adata = target_adata[:, raw_adata.var['highly_variable']].copy()

source_adata = remove_sparsity(source_adata)
target_adata = remove_sparsity(target_adata)


###### Create SCANVI model and train it on fully labelled reference dataset ####

sca.models.SCVI.setup_anndata(source_adata, batch_key=condition_key, labels_key=cell_type_key)

vae = sca.models.SCVI(
    source_adata,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)

vae.train()

scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown",labels_key=cell_type_key,)
print("Labelled Indices: ", len(scanvae._labeled_indices))
print("Unlabelled Indices: ", len(scanvae._unlabeled_indices))
scanvae.train(max_epochs=20,n_samples_per_label=100)


# Save the reference model
ref_path = os.path.join(model_dir,'ref_model/')
scanvae.save(ref_path, overwrite=True)


################# Perform surgery on reference model and train #################
###################### on query dataset without cell type labels ###########

############################################################################
# if there is no ‘.obs’ in the anndata for cell type labels (e.g. the data is unlabeled), one can only use scANVI in an unsupervised manner during surgery due to the nature of the classifier.
############################################################################


# Prepare target adata
target_adata.obs[condition_key] = "study1" 
target_adata.obs[cell_type_key] = "Unknown"


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
surgery_path = os.path.join(model_dir,'surgery_model')
model.save(surgery_path, overwrite=True)

# Get the latent representation of the query dataset
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = target_adata.obs[cell_type_key].tolist()
query_latent.obs['batch'] = target_adata.obs[condition_key].tolist()
query_latent.obs['predictions'] = model.predict()
query_latent.write_h5ad(os.path.join(adata_dir,'query_latent_myeloid.h5ad'))

####### EMBBED BOTH THE REFERENCE AND SURGERY MODEL IN THE SAME SPACE #########

adata_full = source_adata.concatenate(target_adata)
full_latent = sc.AnnData(model.get_latent_representation(adata=adata_full))
full_latent.obs['cell_type'] = adata_full.obs[cell_type_key].tolist()
full_latent.obs['batch'] = adata_full.obs[condition_key].tolist()
full_latent.write_h5ad(os.path.join(adata_dir,'full_latent_myeloid.h5ad'))