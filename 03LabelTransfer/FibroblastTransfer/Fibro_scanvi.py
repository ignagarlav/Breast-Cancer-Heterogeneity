################################## Imports ################################
import scanpy as sc
import scvi
import os
import torch
import numpy as np

torch.manual_seed(311224)
np.random.seed(311224)


print(f"Using torch version {torch.__version__}")
print(f"Using scvi version {scvi.__version__}")
print(f"CUDA: {torch.cuda.is_available()}")
print(f"CUDA version: {torch.version.cuda}")

if torch.cuda.is_available():
    accelerator = "gpu"
else:
    accelerator = "cpu"

############################### Data Loading ###############################
DATA_DIR = os.getenv("DATA_DIR")
MODEL_DIR = os.getenv("MODEL_DIR")  
CELLTYPE = os.getenv("CELL_TYPE") # column name for cell type

if not os.path.exists(MODEL_DIR):
    os.makedirs(MODEL_DIR, exist_ok=True)


adata = sc.read_h5ad(os.path.join(DATA_DIR,"post_scarches_post_scanvi_adata.h5ad"))

############################### SCVI #######################################
scvi.model.SCVI.setup_anndata(adata, 
                            batch_key="batch", 
                            layer="counts")

model = scvi.model.SCVI(adata, 
                        n_layers=3, 
                        n_latent=80, 
                        gene_likelihood="nb")


model.train(max_epochs=400, 
             early_stopping=True,           
            plan_kwargs={"lr": 1e-3})


# guardar modelo de scvi
scvi_ref_path = os.path.join(MODEL_DIR, "scvi_model")
model.save(scvi_ref_path, overwrite=True)

#################### SCANVI ############################

adata.obs[CELLTYPE] = adata.obs[CELLTYPE].astype('category')
print(f"Unique cell types: {adata.obs[CELLTYPE].unique()}")

mask = adata.obs[CELLTYPE] != "unknown"
n_label = adata.obs.loc[mask,CELLTYPE].count()
n_unlabel = adata.obs.loc[~mask,CELLTYPE].count()
print(f'Labelled cells: {n_label}')
print(f'Unlabelled cells: {n_unlabel}')


# Load the model architecture (e.g., SCVI)
scanvi_model = scvi.model.SCANVI.from_scvi_model(model, 
                                      unlabeled_category="unknown",
                                      labels_key=CELLTYPE) 


scanvi_model.train(max_epochs=40, n_samples_per_label=100)

# guardar modelo de scanvi
scanvi_ref_path = os.path.join(MODEL_DIR, "scanvi_model_2")
scanvi_model.save(scanvi_ref_path, overwrite=True)
print(f"Model saved in {scanvi_ref_path}")