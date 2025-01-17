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


data_dir = os.getenv("DATA_DIR") # where should be the data be saved
model_dir = os.getenv("MODEL_DIR") # where should the model be saved

if not os.path.exists(data_dir):
    os.makedirs(data_dir, exist_ok=True)
if not os.path.exists(model_dir):
    os.makedirs(model_dir, exist_ok=True)


print(f"Using torch version {torch.__version__}")
print(f"Using scvi version {scvi.__version__}")
print(f"CUDA: {torch.cuda.is_available()}")
print(f"CUDA version: {torch.version.cuda}")

if torch.cuda.is_available():
    accelerator = "gpu"
else:
    accelerator = "cpu"

############### Cargar datos recién preprocesados ##############################
adata = sc.read_h5ad('/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/01_Preprocessing/adata_post_qc.h5ad')
adata_hvgs = adata[:,adata.var['highly_variable']].copy()

############################### Integar con scVI ##############################
scvi.model.SCVI.setup_anndata(adata_hvgs, 
                            batch_key="batch", # each patient
                            layer="counts")


model = scvi.model.SCVI(adata_hvgs, 
                        n_layers=4, 
                        n_latent=80, 
                        gene_likelihood="nb")


model.train(max_epochs=400, 
            early_stopping=True,           
            plan_kwargs={"lr": 1e-3})


# Save ScVi model
scvi_ref_path = os.path.join(model_dir, "scvi_model_cuda")
model.save(scvi_ref_path, overwrite=True)

##################### Transferir datos al objeto adata ###########
SCVI_LATENT_KEY = "X_scVI"

adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata,resolution=0.5, flavor="igraph", n_iterations=2)
sc.tl.umap(adata)
# Save Adata
adata.write_h5ad(os.path.join(data_dir,"adata_scvi_cuda.h5ad"))
