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
query_adata_dir = os.getenv('QUERY_ADATA_DIR')


condition_key = "harm_study"
cell_type_key = "author_cell_type"


vae_epochs = 500
scanvi_epochs = 200
surgery_epochs = 500


## Load reference dataset
reference_adata_dir = os.path.join(adata_dir,'MyeloidAdata.h5ad')
orig_adata = sc.read_h5ad(reference_adata_dir)
raw_adata = orig_adata.raw.to_adata()

## Convert ensembl ids to gene symbols
mg = mygene.MyGeneInfo()
ensembl_ids = list(raw_adata.var_names)
results = mg.querymany(
    ensembl_ids, 
    scopes="ensembl.gene", 
    fields="symbol", 
    species="human"
)
id_map = {}
for entry in results:
    if not entry.get("notfound", False):
        ens_id = entry["query"]
        symbol = entry.get("symbol")
        if symbol:
            id_map[ens_id] = symbol

raw_adata.var_names = [id_map.get(ens_id, ens_id) for ens_id in raw_adata.var_names]
raw_adata.var_names_make_unique()


# Load the query dataset
target_adata = sc.read_h5ad(query_adata_dir)
target_adata.X = target_adata.layers["counts"]

common_genes = [gene for gene in raw_adata.var_names if gene in target_adata.var_names]
print(f"Number of common genes: {len(common_genes)}")

target_adata = target_adata[:, common_genes].copy()
raw_adata = raw_adata[:, common_genes].copy()


# Compute HVGs on raw data 
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

scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")
print("Labelled Indices: ", len(scanvae._labeled_indices))
print("Unlabelled Indices: ", len(scanvae._unlabeled_indices))
scanvae.train(max_epochs=20)


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
    freeze_dropout = True,
)
model._unlabeled_indices = np.arange(target_adata.n_obs)
model._labeled_indices = []
print("Labelled Indices: ", len(model._labeled_indices))
print("Unlabelled Indices: ", len(model._unlabeled_indices))

model.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))


# Save the surgery model 
surgery_path = os.path.join(model_dir,'surgery_model')
model.save(surgery_path, overwrite=True)

# Get the latent representation of the query dataset
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = target_adata.obs[cell_type_key].tolist()
query_latent.obs['batch'] = target_adata.obs[condition_key].tolist()
query_latent.obs['predictions'] = model.predict()
query_latent.write_h5ad(os.path.join(adata_dir,'query_latent.h5ad'))

####### EMBBED BOTH THE REFERENCE AND SURGERY MODEL IN THE SAME SPACE #########

adata_full = source_adata.concatenate(target_adata)
full_latent = sc.AnnData(model.get_latent_representation(adata=adata_full))
full_latent.obs['cell_type'] = adata_full.obs[cell_type_key].tolist()
full_latent.obs['batch'] = adata_full.obs[condition_key].tolist()
full_latent.write_h5ad(os.path.join(adata_dir,'full_latent.h5ad'))