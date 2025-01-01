import os
import tempfile

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
import torch

print("Last run with scvi-tools version:", scvi.__version__)

base_dir = "/Users/joseignaciogarzonalvarez/Proyectos-R/SCRNASEQ/GSE162929_176078/data/"


####################### Generar la referencia ##########

## Cargar el objecto de referencia 
adata_path = os.path.join(base_dir, "PaperBC.h5ad")
adata_ref = sc.read(adata_path)


sc.pp.normalize_total(adata_ref, target_sum=1e4)
sc.pp.log1p(adata_ref)
adata_ref.raw = adata_ref  # keep full dimension safe
sc.pp.highly_variable_genes(
    adata_ref,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
    batch_key="sample",
    subset=True,
)

adata_ref.obs["batch"] = adata_ref.obs["batch"].astype("category")

scvi.model.SCVI.setup_anndata(adata_ref, 
                              batch_key="batch", 
                              layer="counts")

## entrenar aquí scvi de referencia 

'''We train the reference using the standard SCVI workflow, except we add a few non-default parameters that were identified to work well with scArches. It is essential to encode covariates here as this allows scArches to map new batches in the encoder to the existing data and thereby provides batch integration'''

scvi_ref = scvi.model.SCVI(
    adata_ref,
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    n_layers=3
)
scvi_ref.train()
# ahí tengo el modelo de referencia 

scvi_ref_path = os.path.join(base_dir, "scvi_ref")
scvi_ref.save(scvi_ref_path, overwrite=True)

############################## Comprobar modelo de referencia ###########
# En realidad opcional 

## ver con scanpy 
SCVI_LATENT_KEY = "X_scVI"

adata_ref.obsm[SCVI_LATENT_KEY] = scvi_ref.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata_ref)
sc.tl.umap(adata_ref)
sc.pl.umap(
    adata_ref,
    color=["sample"],
    frameon=False,
    ncols=1,
)


################ Meter query al modelo de SCVI ############################
## Cargar el query
adata_query = sc.read("/Users/joseignaciogarzonalvarez/Proyectos-R/SCRNASEQ/GSE162929_176078/data/FibroblastData/Adata.h5ad")
scvi.model.SCVI.prepare_query_anndata(adata_query, scvi_ref_path)

scvi_query = scvi.model.SCVI.load_query_data(
    adata_query,
    scvi_ref_path,
)
scvi_query.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})

adata_query.obsm[SCVI_LATENT_KEY] = scvi_query.get_latent_representation()

sc.pp.neighbors(adata_query, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata_query)
sc.tl.umap(adata_query)

sc.pl.umap(
    adata_query,
    color=["batch"],
    frameon=False,
    ncols=1,
)

#################### Crear modelo scanvi referencia ############################

## crear modelo scanvi desde el modelo scvi
SCANVI_LABELS_KEY = "labels_scanvi"

adata_ref.obs[SCANVI_LABELS_KEY] = adata_ref.obs["CAFtype"]
np.unique(adata_ref.obs[SCANVI_LABELS_KEY], return_counts=True)

scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_ref, 
                                      unlabeled_category="Unknown", # lo que usaremos luego para indicar cuales son las células no anotadas
                                      labels_key=SCANVI_LABELS_KEY) # la clave de las etiquetas en adata_ref.obs que se usará para entrenar el modelo
scanvi_model.train(max_epochs=100, n_samples_per_label=100)

# guardar modelo de scanvi
scanvi_ref_path = os.path.join(base_dir, "scanvi_ref")
scanvi_model.save(scanvi_ref_path, overwrite=True)

###################### ver el modelo de referencia 
SCANVI_LATENT_KEY = "X_scANVI"
adata_ref.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep=SCANVI_LATENT_KEY)
sc.tl.leiden(adata_ref)
sc.tl.umap(adata_ref)
sc.pl.umap(
    adata_ref,
    color=["sample"],
    frameon=False,
    ncols=1,
)
############################ usar el modelo ############################
#### Aquí es donde hacemos la transferencia de etiquetas

## Cargar el query
adata_query = sc.read("/Users/joseignaciogarzonalvarez/Proyectos-R/SCRNASEQ/GSE162929_176078/data/FibroblastData/Adata.h5ad")

scvi.model.SCANVI.prepare_query_anndata(adata_query, scanvi_ref_path)

## Cargar el query al modelo 
scanvi_query = scvi.model.SCANVI.load_query_data(adata_query, scanvi_ref_path)


# Entrenar el modelo
surgery_epochs = 500
train_kwargs_surgery = {
    "early_stopping": True,
    "early_stopping_monitor": "elbo_train",
    "early_stopping_patience": 10,
    "early_stopping_min_delta": 0.001,
    "plan_kwargs": {"weight_decay": 0.0},
}
scanvi_query.train(max_epochs=surgery_epochs, **train_kwargs_surgery)

# Guardar el modelo
query_model_path = os.path.join(base_dir, "query_model")
scanvi_query.save(query_model_path, overwrite=True)

################## Ver las predicciones sobre la query #########################

SCANVI_PREDICTIONS_KEY = "predictions_scanvi"

adata_query.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation()
adata_query.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query.predict()


sc.pp.neighbors(adata_query, use_rep=SCANVI_LATENT_KEY, n_neighbors=10)
sc.tl.leiden(adata_query)
sc.tl.umap(adata_query,min_dist=0.1,  spread=1)
sc.pl.umap(
    adata_query,
    color=[SCANVI_PREDICTIONS_KEY],
    frameon=False,
    ncols=1,
)


adata_query.write_h5ad("/Users/joseignaciogarzonalvarez/Proyectos-R/SCRNASEQ/GSE162929_176078/data/FibroblastData/FibroblastAnnotated.h5ad")

pd.DataFrame(adata_query.obsm[SCANVI_LATENT_KEY]).to_csv("/Users/joseignaciogarzonalvarez/Proyectos-R/SCRNASEQ/GSE162929_176078/data/FibroblastData/scanvi_latent.csv")


pd.DataFrame(adata_query.obs[SCANVI_PREDICTIONS_KEY]).to_csv("/Users/joseignaciogarzonalvarez/Proyectos-R/SCRNASEQ/GSE162929_176078/data/FibroblastData/scanvi_predictions.csv")

pd.DataFrame(adata_query.obsm[SCVI_LATENT_KEY]).to_csv("/Users/joseignaciogarzonalvarez/Proyectos-R/SCRNASEQ/GSE162929_176078/data/FibroblastData/scvi_latent.csv")