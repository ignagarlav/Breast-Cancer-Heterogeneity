import scanpy as sc
import pandas as pd
import pybiomart as pbm
import os 

adata_dir = os.getenv('ADATA_DIR')

orig_adata_fname = os.path.join(adata_dir, 'MyeloidAdata.h5ad')
orig_adata = sc.read_h5ad(orig_adata_fname)
orig_adata = orig_adata.raw.to_adata()

query_adata_fname = os.getenv('QUERY_ADATA_DIR')
query_adata = sc.read_h5ad(query_adata_fname)


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

n_overlap = {}

for version in ensembl_version_dict.keys():

    print(f'host: {version}')

    try:

        n_overlap[version] =  test_ensembl_host(orig_adata, ensembl_version_dict[version])

    except:

        print('Host not reachable')

v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]

host_to_use = ensembl_version_dict[v]
print(f"version: {v} has the largest overlap, use {host_to_use} as biomart host")



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
    
mod_adata = map_genes(orig_adata, host_to_use)



ncommon = len(set(query_adata.var_names).intersection(set(mod_adata.var_names)))
print(f'Number of common genes: {ncommon}')


common_genes = [gene for gene in mod_adata.var_names if gene in query_adata.var_names]

query_matched = query_adata[:,common_genes].copy()

mod_adata.var_names_make_unique()
mod_adata_matched = mod_adata[:,common_genes].copy()    

query_matched.X = query_matched.layers["counts"]

query_matched.var_names_make_unique()
mod_adata_matched.var_names_make_unique()


nmatch = len(set(query_matched.var_names).intersection(set(mod_adata_matched.var_names)))
nquerybefore = len(query_adata.var_names)
nqueryafter = len(query_matched.var_names)
print(f'Number of genes in query before: {nquerybefore}')
print(f'Number of genes in query after: {nqueryafter}')
print(f'Number of genes in common after filtering finishes: {nmatch}')



query_matched_fname = os.path.join(adata_dir, 'ScanviPredictionsAdata_mapped_to_myeloid.h5ad')
query_matched.write_h5ad(query_matched_fname)


mod_adata_fname = os.path.join(adata_dir, 'MyeloidAdata_mapped.h5ad')
mod_adata_matched.write_h5ad(mod_adata_fname)



