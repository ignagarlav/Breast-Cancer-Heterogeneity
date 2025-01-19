################################## Imports ################################
import scanpy as sc
import scvi
import os
import torch
import numpy as np

torch.manual_seed(311224)
np.random.seed(311224)

############################### Data Loading ###############################
data_dir = os.getenv("DATA_DIR")
model_dir = os.getenv("MODEL_DIR")  

SCANVI_LATENT_KEY = "X_scANVI"

adata = sc.read_h5ad(os.path.join(data_dir,"Fibro_inital_scarches_pred.h5ad"))

############################### SCVI #######################################

if not os.path.exists(model_dir):
    os.makedirs(model_dir, exist_ok=True)

scvi.model.SCVI.setup_anndata(adata, 
                            batch_key="batch", 
                            layer="counts")

model = scvi.model.SCVI(adata, 
                        n_layers=4, 
                        n_latent=80, 
                        gene_likelihood="nb")


model.train(max_epochs=400, 
             early_stopping=True,           
            plan_kwargs={"lr": 1e-3})


# guardar modelo de scvi
scvi_ref_path = os.path.join(model_dir, "scvi_model")
model.save(scvi_ref_path, overwrite=True)

#################### SCANVI ############################

# Load the model architecture (e.g., SCVI)
adata.obs['cell_type'] = adata.obs['scArches_Initial'].copy()
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
print(f"Unique cell types: {adata.obs['cell_type'].unique()}")

scanvi_model = scvi.model.SCANVI.from_scvi_model(model, 
                                      unlabeled_category="unknown",
                                      labels_key='cell_type') 


scanvi_model.train(max_epochs=40, n_samples_per_label=100)

# guardar modelo de scanvi
scanvi_ref_path = os.path.join(model_dir, "scanvi_model")
scanvi_model.save(scanvi_ref_path, overwrite=True)
print(f"Model saved in {scanvi_ref_path}")
