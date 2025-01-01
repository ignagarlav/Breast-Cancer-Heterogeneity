from scvi.model import SCANVI
import scvi
import scanpy as sc

import scanpy as sc
import os 
base_dir = "/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/LatentSpace"
adata = sc.read_h5ad(os.path.join(base_dir,"adatascanvi1.h5ad")) 
adata_hvg = adata[:, adata.var.highly_variable].copy()

model_dir = "/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/LatentSpace/models/scvi_model/scanvi_model"

# Load the saved model
model = SCANVI.load(model_dir, adata=adata_hvg)


adata.obs['predicted_labels'] = model.predict()

sc.pl.embedding(
    adata,
    basis="X_scanvi_MDE",
    color=["predicted_labels"],
    frameon=False,
    ncols=1)