'''El modelo de scanvi está pensado para refinar el espacio latente generado por scvi. 
El objeto adata se espera que esté integrado por scvi y que haya sido anotado con el pipeline de General Anno

El pipeline completo scvi-scanvi se puede ver en 

https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/harmonization.html


Se espera que previamente se hayan ejecutado: 

1. SCPreprocessing
2. Before Integration 

3. IntegrationScvi
4. GeneralAnnoBeforeScAnvi

This script IntegrationScAnvi is the 5th step of the pipeline. It refines the latent space generated by scvi using scanvi.


'''
import scanpy as sc
import scvi
import os
import torch

data_dir = os.getenv("DATA_DIR")
adata = sc.read_h5ad(os.path.join(data_dir,"adata_GenAnno.h5ad"))

#################### Crear modelo scanvi ############################

model_dir = os.getenv("MODEL_DIR")
scvi_ref_path = os.path.join(model_dir, "scvi_model")


# Load the model architecture (e.g., SCVI)
adata_hvg = adata[:, adata.var.highly_variable].copy()
model = scvi.model.SCVI.load(scvi_ref_path,adata_hvg)


adata.obs['cell_type'] = adata.obs['GenAnno']
print(f"Unique cell types: {adata.obs['cell_type'].unique()}")

# No hace falta aquí setup andata porque lo hace SCANVI 
scanvi_model = scvi.model.SCANVI.from_scvi_model(model, 
                                      unlabeled_category="Unknown",
                                      labels_key='cell_type') 


scanvi_model.train(max_epochs=40, n_samples_per_label=100)

# guardar modelo de scanvi
scanvi_ref_path = os.path.join(model_dir, "scanvi_model")
scanvi_model.save(scanvi_ref_path, overwrite=True)
print(f"Model saved in {scanvi_ref_path}")

###################### Transferir espacio latente al objeto adata ###########
SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()
sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY)
sc.tl.leiden(adata)

if torch.cuda.is_available():
    accelerator = "gpu"
else:
    accelerator = "cpu"

adata.obsm["X_scanvi_MDE"] = scvi.model.utils.mde(adata.obsm[SCANVI_LATENT_KEY], accelerator=accelerator)

# Guardar el objeto
adata.write_h5ad(os.path.join(data_dir,"adata_scanvi.h5ad"))