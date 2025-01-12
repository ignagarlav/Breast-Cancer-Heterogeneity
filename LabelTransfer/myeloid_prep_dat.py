################################## Imports ################################
import scanpy as sc 
import pandas as pd
import pybiomart as pbm
import os 

############################## Helper Functions ################################
ensembl_version_dict = {'105': 'http://www.ensembl.org',

                        '104': 'http://may2021.archive.ensembl.org/',

                        '103': 'http://feb2021.archive.ensembl.org/',

                        '102': 'http://nov2020.archive.ensembl.org/',

                        '101': 'http://aug2020.archive.ensembl.org/',

                        '100': 'http://apr2020.archive.ensembl.org/',

                        '99': 'http://jan2020.archive.ensembl.org/',

                        '98': 'http://sep2019.archive.ensembl.org/',

                        '97': 'http://jul2019.archive.ensembl.org/',

                        '96': 'http://apr2019.archive.ensembl.org/',

                        '95': 'http://jan2019.archive.ensembl.org/',

                        '94': 'http://oct2018.archive.ensembl.org/',

                        '93': 'http://jul2018.archive.ensembl.org/',

                        '92': 'http://apr2018.archive.ensembl.org/',

                        '91': 'http://dec2017.archive.ensembl.org/',

                        '90': 'http://aug2017.archive.ensembl.org/',

                        '89': 'http://may2017.archive.ensembl.org/',

                        '88': 'http://mar2017.archive.ensembl.org/',

                        '87': 'http://dec2016.archive.ensembl.org/',

                        '86': 'http://oct2016.archive.ensembl.org/',

                        '80': 'http://may2015.archive.ensembl.org/',

                        '77': 'http://oct2014.archive.ensembl.org/',

                        '75': 'http://feb2014.archive.ensembl.org/',

                        '54': 'http://may2009.archive.ensembl.org/'}

def test_ensembl_host(adata, host):

    dataset = pbm.Dataset(name='hsapiens_gene_ensembl', host=host)
    
    qp = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
    print(f'Nº of genes in db: {len(qp)}')
    
    df_var = pd.DataFrame(adata.var_names.values, columns=["Gene stable ID"])
    print(f'Nº of genes in dataset: {len(df_var)}')
    
    df_mapped = df_var.merge(qp, how="inner", on="Gene stable ID")
    print(f'Nº of mapped genes to db: {len(df_mapped)}')
    
    df_valid = df_mapped[~df_mapped["Gene name"].isna()]
    ov = len(df_valid)
    print(f'Nº of valid mapped genes to db: {ov}')
    
    return(ov)

def map_genes(adata, host):

    dataset = pbm.Dataset(name='hsapiens_gene_ensembl', host=host)
    
    qp = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
    print(f'Nº of genes in db: {len(qp)}')
    
    df_var = pd.DataFrame(adata.var_names.values, columns=["Gene stable ID"])
    print(f'Nº of genes in dataset: {len(df_var)}')
    
    df_mapped = df_var.merge(qp, how="inner", on="Gene stable ID")
    print(f'Nº of mapped genes to db: {len(df_mapped)}')
    
    df_valid = df_mapped[~df_mapped["Gene name"].isna()]

    
    gene_order = df_valid["Gene stable ID"].tolist()

    adata = adata[:, gene_order].copy() 

    adata.var_names = df_valid['Gene name'].values
    
    return(adata)   


################################## Data Loading ################################
# Load reference adata
adata_dir = os.getenv('ADATA_DIR')
ref_adata_fname = os.path.join(adata_dir, 'MyeloidAdata.h5ad')
ref_adata = sc.read_h5ad(ref_adata_fname)
ref_adata = ref_adata.raw.to_adata()

# Load query adata
query_adata_fname = os.getenv('QUERY_ADATA_DIR')
query_adata = sc.read_h5ad(query_adata_fname)
# Make sure to work with raw data (counts)
query_adata.X = query_adata.layers["counts"]
# Subset to desired cell types
query_annotation = "GennAnno_ScAnvi"
query_adata = query_adata[query_adata.obs[query_annotation].isin(['Myeloid','Mast cells']),:].copy()

################################## Gene Mapping ################################
n_overlap = {}
for version in ensembl_version_dict.keys():

    print(f'host: {version}')

    try:

        n_overlap[version] =  test_ensembl_host(ref_adata, ensembl_version_dict[version])

    except:

        print('Host not reachable')
v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
host_to_use = ensembl_version_dict[v]
print(f"version: {v} has the largest overlap, use {host_to_use} as biomart host")

# Use ensemble biomart to map genes in the reference dataset  
ref_adata_renamed = map_genes(ref_adata, host_to_use)


# Subset query and reference with common genes
gcom = [g for g in ref_adata_renamed.var_names if g in query_adata.var_names]

query_matched = query_adata[:,gcom].copy()
ref_adata_renamed.var_names_make_unique()
ref_adata_matched = ref_adata_renamed[:,gcom].copy()    

# Make sure that the genes are unique
query_matched.var_names_make_unique()
ref_adata_matched.var_names_make_unique()

nmatch = len(set(query_matched.var_names).intersection(set(ref_adata_matched.var_names)))
print(f'Number of genes in common after filtering finishes: {nmatch}')

############################## Data Harmonization ##############################

# Compute HVGs on raw data of reference dataset
sc.pp.highly_variable_genes(ref_adata_matched, flavor='seurat_v3', n_top_genes=3000)

# Slice raw and target adatas to keep HVGs
ref_adata_matched = ref_adata_matched[:, ref_adata_matched.var['highly_variable']].copy()
query_matched = query_matched[:, ref_adata_matched.var['highly_variable']].copy()

############################## Data Saving ##############################

query_matched_fname = os.path.join(adata_dir, 'ScanviPredictionsAdata_mapped_to_myeloid.h5ad')
query_matched.write_h5ad(query_matched_fname)


ref_adata_matched_fname = os.path.join(adata_dir, 'MyeloidAdata_mapped.h5ad')
ref_adata_matched.write_h5ad(ref_adata_matched_fname)