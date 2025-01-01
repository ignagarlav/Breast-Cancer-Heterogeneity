'''
Se espera un objeto directamente proveniente del paso "Before Integration", 
con adata.X conteniendo los conteos normalizados.

En ese mismo paso se espera que se haya realizado la selección de genes altamente variables.

Para el entrenamiento del modelo se usan sólo los hvgs, y luego se transfieren los resultados al objeto adata original.
'''
import scanpy as sc
import scvi
import os
import torch
import numpy as np

torch.manual_seed(311224)
np.random.seed(311224)


############### Cargar datos recén preprocesados ###############################
data_dir = os.getenv("DATA_DIR")

adata = sc.read_h5ad(os.path.join(data_dir,"adata_normalized.h5ad"))
adata_hvgs = adata[:,adata.var['highly_variable']].copy()

############################### Integar con scVI ##############################
model_dir = os.getenv("MODEL_DIR")  

if not os.path.exists(model_dir):
    os.makedirs(model_dir, exist_ok=True)


scvi.model.SCVI.setup_anndata(adata_hvgs, 
                            batch_key="batch", #  Eliminará los efectos técnicos entre las pacientes. 
                            labels_key="subtype", # Preservará las diferencias biológicas relevantes entre los subtipos tumorales (TNBC, ER+, HER2+).
                            layer="counts")


model = scvi.model.SCVI(adata_hvgs, 
                        n_layers=4, 
                        n_latent=80, 
                        gene_likelihood="nb")


model.train(max_epochs=400, 
             early_stopping=True,           
            plan_kwargs={"lr": 1e-3})


# guardar modelo de scvi
scvi_ref_path = os.path.join(model_dir, "scvi_model")
model.save(scvi_ref_path, overwrite=True)

##################### Transferir datos al objeto adata ###########
SCVI_LATENT_KEY = "X_scVI"

adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata,resolution=0.5, flavor="igraph", n_iterations=2)

adata.obsm["X_scvi_MDE"] = scvi.model.utils.mde(adata.obsm[SCVI_LATENT_KEY], accelerator="cpu")

# Guardar adata
output_dir = os.getenv("OUTPUT_DIR")
adata.write_h5ad(os.path.join(output_dir,"adata_scvi.h5ad"))
