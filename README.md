# Breast-Cancer-Heterogeneity
Source code for the analysis of the heterogenicity of ER/TNBC/HER2 tumors 

# Single Cell Data 
Raw data was retrieved from GEO accesion viewer [GEO](https://www.ncbi.nlm.nih.gov/geo/) under accession number [GSE161529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529)

# Single cell Data Preprocessing

Single cell data was preprocessed using the preprocessingV3 module of the SCPipeline package. 

Normalization: adata_normalized

# Integration 

Integration with scvi adata_scvi
General annotation: adata_GenAnno
Integration with scanvi: adata_scanvi 
With scanvi predictions: adata_scanvi_predictions

# Proportion test

Using adata_scanvi_predictions Proportion test was run
