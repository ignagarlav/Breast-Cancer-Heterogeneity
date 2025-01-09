'''
Se carga el objecto adata_scvi.
Después de hacer la integración con scvi, se realiza una anotación de las células usando las predicciones de celltypist y marcadores canónicos. 
Se guarda la anotación en "GenAnno". 
Esta anotación servirá a scanvi para mejora el espacio latente- '''


import scanpy as sc
import os


data_dir = "/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/LatentSpace"

adata = sc.read_h5ad(os.path.join(data_dir,"adata_scvi_complete.h5ad"))





# "EPCAM","CDH1" epithelial
# "MKI67" proliferative
# "CD3D", T cells
# "MS4A1" B cells
# JCHAIN plasma cells
# "CD68","CD14", "ITGAX" macrophages
# "PECAM1","VWF" endothelial
# "PDGFRB","FAP" fibroblasts
# MS4A2 Mast/Basophiphs cells

general_markers  = ["EPCAM","CDH1", "MKI67", "CD3D", "MS4A1", "CD68","CD14", "ITGAX","JCHAIN", "PECAM1","VWF", "PDGFRB","FAP","MS4A2","ACTA2"]


sc.pl.dotplot(adata=adata, var_names=general_markers, groupby="leiden")


rename_dict = {"0": "Myeloid",
                "1": "Unknown",#ActivePlasmaCells  
               "2": "Fibroblasts",
               "3": "Fibroblasts",
               "4": "B Cells",
               "5": "Epithelial",
               "6": "PlasmaCells",
               "7": "Epithelial",
               "8": "CyclingEpithelial",
               "9": "Epithelial",
               "10": "Epithelial",
               "11": "Epithelial",
               "12": "Epithelial",
               "13": "Epithelial",
               "14": "TCells",
               "15": "Unknown", # ACTA2 high
               "16": "MastCells",
               "17": "Endothelial",
               "18": "Unknown",
               "19": "Myeloid"}


# Renombrar los valores en la columna deseada
adata.obs['GenAnno'] = adata.obs['leiden'].replace(rename_dict)
sc.pl.embedding(
    adata,
    basis="X_scvi_MDE",
    color=["leiden"],
    frameon=False, 
    #legend_loc="on data",
    ncols=1,
)
sc.pl.dotplot(adata=adata, var_names=general_markers, groupby="GenAnno")

sc.pl.embedding(
    adata,
    basis="X_scvi_MDE",
    color=["GenAnno"],
    frameon=False, 
    #legend_loc="on data",
    ncols=1,
)


adata.write_h5ad(os.path.join(data_dir,"adata_GenAnno.h5ad"))
