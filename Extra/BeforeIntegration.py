import scanpy as sc
import scvi
import torch
import os 
torch.manual_seed(311224)

data_dir = os.getenv("DATA_DIR")
# Cargar datos
adata = sc.read_h5ad(os.path.join(data_dir,"adata.h5ad"))

########################## Ajustar los metadatos #############################
her_samples = ['AH0308-', 'MH0031-', 'MH0069-', 'MH0161-', 'MH0176-', 'PM0337-']
tncb_samples_1 = ['B1-MH0131-', 'B1-MH0177-', 'B1-MH4031-', 'B1-Tum0554-',  'MH0114-T2-', 'MH0126-', 'MH0135-', 'SH0106-']
er_samples = ['AH0319-', 'MH0001-', 'MH0025-', 'MH0032-', 'MH0040-', 'MH0042-', 'MH0043-T-', 'MH0114-T3-', 'MH0125-', 'MH0151-', 'MH0163-', 'MH0167-T-', 'PM0360-']

def classify_sample(batch):
    if batch in tncb_samples_1:
        return 'TNBC'
    elif batch in er_samples:
        return 'ER'
    elif batch in her_samples:
        return 'HER2'
    else:
        return 'TNBC'

adata.obs['subtype'] = adata.obs['batch'].apply(classify_sample).astype('category')
adata.obs['batch'] = adata.obs['batch'].astype('category')

########################## Normalizar ########################################
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata, n_comps=30)
adata.obsm["X_pca_MDE"] = scvi.model.utils.mde(adata.obsm["X_pca"], accelerator="cpu")

adata.raw = adata # guardar datos normalizados (data, no counts)

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=4000,
    layer="counts",
    batch_key="sample",
    subset=False) 

########################### Guardar ##############################
adata.write_h5ad(os.path.join(data_dir,"adata_normalized.h5ad"))