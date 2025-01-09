'''PyScenic pipeline https://www.nature.com/articles/s41596-020-0336-2

Pre-processing
The SCENIC workflow starts with an expression matrix capturing the abundance of each gene’s transcript in every cell interrogated in an scRNA-seq experiment. Typically, each cell is represented as a separate row in this matrix, whereas genes are depicted as columns.

Network inference
In a first step, given a predefined list of TFs, regulatory interactions between these factors and putative target genes are inferred via regression-based network inference10 from the expression or count matrix

'''

# Import libraries 
import os
from scipy.io import mmread
import pandas as pd
from os import listdir
import pyarrow
from collections import defaultdict
import operator as op
from IPython.display import HTML, display
import numpy as np
from cytoolz import compose
import scanpy as sc

# --------------- Auxiliar functions --------------- #
# STEP 3. Create a function that displays the logos of the top enriched motifs.

BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"
COLUMN_NAME_LOGO = "MotifLogo"
COLUMN_NAME_MOTIF_ID = "MotifID"
COLUMN_NAME_TARGETS = "TargetGenes"

def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
    """
    :param df:
    :param base_url:
    """
    # Make sure the original dataframe is not altered.
    df = df.copy()
    
    # Add column with URLs to sequence logo.
    def create_url(motif_id):
        return '<img src="{}{}.png" style="max-height:124px;"></img>'.format(base_url, motif_id)
    df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    
    # Truncate TargetGenes.
    def truncate(col_val):
        return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
    df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))
    
    MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
    pd.set_option('display.max_colwidth', None)
    display(HTML(df.head().to_html(escape=False)))



# Step 4. Filter regulons
databases_directory = "/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/databases"
cis_target_data = {}
for file in listdir(databases_directory):
    if file.endswith(".feather"):
        cis_target_data[file] = pd.read_feather(os.path.join(databases_directory,file))  

databases_names = list(cis_target_data.keys())
databases_names = [os.path.splitext(filename)[0] for filename in databases_names] 

def derive_regulons(motifs, db_names=databases_names):
    if isinstance(motifs.columns, pd.MultiIndex):
        motifs.columns = motifs.columns.droplevel(0)

    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f


    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motifs for the modules that were created using the method 'weight>50.0%' (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.

    # This filters the motifs DataFrame, keeping only the rows whose Context value matches all 3 conditions:
    # It does not include 'weight>50.0%'
    # It must contain at least one of the strings in db_names
    # It must contain the string 'activating'


    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains('weight>50.0%')), motifs.Context), dtype=np.bool) & \
        np.fromiter(map(contains(*db_names), motifs.Context), dtype=np.bool) & \
        np.fromiter(map(contains('activating'), motifs.Context), dtype=np.bool)]

    # 1. We build regulons only using enriched motifs with a NES of 3.0 or higher; 
    # 
    # 2. We take only either 
    # directly annotated TFs 
    # TF annotated for an orthologous gene into account; 
    # 
    # 3. We only keep regulons with at least 10 genes.
    
    regulons = list(filter(lambda r: len(r) >= 10, 
    
    df2regulons(motifs[(motifs['NES'] >= 1.0)  & 
    ((motifs['Annotation'] == 'gene is directly annotated') | 
    
    (motifs['Annotation'].str.startswith('gene is orthologous to') & motifs['Annotation'].str.endswith('which is directly annotated for motif'))) 
        ])))
    
    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))


# Step 6. 
def savesvg(fname: str, fig, folder: str=FIGURES_FOLDERNAME) -> None:
    """
    Save figure as vector-based SVG image format.
    """
    fig.tight_layout()
    fig.savefig(os.path.join(folder, fname), format='svg')

# --------------- Pipeline --------------- # 

# ------------------- STEP 0: Load data ---------------------------- #
adata = sc.read_h5ad("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/data/BC_All_WithoutNormal/adata_integrated.h5ad")

adata_subset = adata[(adata.obs['GenAnno'].isin(['Cycling Epithelial']))&(adata.obs['subtype'].isin(['ER','HER2']))].copy()


expression_matrix = adata_subset.X
gene_names = adata_subset.var_names.to_list()

# ----------- STEP 1: Network inference based on GRNBoost2 ---------- #
####################################################################
#Candidate regulatory modules are inferred from coexpression patterns between genes 
####################################################################
'''pyscenic grn expression_data.csv TF_names_v_1.01.txt -o adjacencies2.tsv
'''


from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

tf_names = load_tf_names("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/grn_files/TF_names_v_1.01.txt")

# network = grnboost2(expression_data=expression_matrix, gene_names=gene_names, tf_names=tf_names)

# ---------- STEP 2-3: Regulon prediction aka cisTarget  ---------- #
####################################################################
# Derive enriched motifs associated with TFs
# Coexpression modules are refined by the elimination of indirect targets using TF motif information
# Incorporating motif information enhances the biological plausibility of inferred TF-target interactions by ensuring that TFs have binding motifs near target genes.
####################################################################

'''!pyscenic ctx {ADJACENCIES_FNAME} {DBS_PARAM} \ --annotations_fname {MOTIF_ANNOTATIONS_FNAME} \ --expression_mtx_fname {EXP_MTX_QC_FNAME} \ --output {MOTIFS_FNAME} \ --num_workers 18'''

# Motifs_fname = regulons 
from pyscenic.utils import load_motifs

# Load the motifs.csv file created by the ctx command.
motifs_df = load_motifs('/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/grn_files/regulons.csv')
display_logos(motifs_df.head())

# ---------- STEP 4: Cellular enrichment aka AUCell  ---------- #
####################################################################
# the activity of these discovered regulons is measured in each individual cell and used for clustering 
# Scores cells for the activity of regulons.
####################################################################

## Step 4.1: REGULON CREATION
from pyscenic.transform import df2regulons
regulons = derive_regulons(motifs_df)

# Pickle these regulons.
import pickle
REGULONS_DAT_FNAME = "/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/grn_files/regulons_new.p"
with open(REGULONS_DAT_FNAME, 'wb') as f:
    pickle.dump(regulons, f)

## Step 4.2: AUCell
from pyscenic.aucell import aucell

expression_data = pd.read_csv("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/grn_files/expression_data.csv", index_col=0)
auc_mtx = aucell(expression_data, regulons)

# ----------  STEP 5 - Regulon activity binarization  ---------- #
from pyscenic.binarization import binarize
import seaborn as sns

bin_mtx, thresholds = binarize(auc_mtx) 
bin_mtx.to_csv(BIN_MTX_FNAME) 
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv(THR_FNAME)

# bin_mtx = pd.read_csv(BIN_MTX_FNAME, index_col=0)
# thresholds = pd.read_csv(THR_FNAME, index_col=0).threshold


### Histogram of binarized results 
import matplotlib as mpl

def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f

N_COLORS = len(adata.obs.cell_type.dtype.categories)
COLORS = [color['color'] for color in mpl.rcParams["axes.prop_cycle"]]
cell_type_color_lut = dict(zip(adata.obs.cell_type.dtype.categories, COLORS))
#cell_type_color_lut = dict(zip(adata.obs.cell_type.dtype.categories, adata.uns['cell_type_colors']))
cell_id2cell_type_lut = df_metadata.set_index('cell_id').cell_type.to_dict()
bw_palette = sns.xkcd_palette(["white", "black"])
sns.set_theme()
sns.set_style("whitegrid")
fig = palplot(bw_palette, ['OFF', 'ON'], ['k', 'w'])
savesvg('legend - GSE115978 - on_off.svg', fig)

sns.set_theme()
sns.set_theme(font_scale=0.8)
fig = palplot(sns.color_palette(COLORS), adata.obs.cell_type.dtype.categories, size=1.0)
savesvg('legend - GSE115978 - cell_type_colors.svg', fig)

sns.set_theme()
sns.set_theme(font_scale=1.0)
sns.set_style("ticks", {"xtick.minor.size": 1, "ytick.minor.size": 0.1})
g = sns.clustermap(bin_mtx.T, 
               col_colors=auc_mtx.index.map(cell_id2cell_type_lut).map(cell_type_color_lut),
               cmap=bw_palette, figsize=(20,20))
g.ax_heatmap.set_xticklabels([])
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_xlabel('Cells')
g.ax_heatmap.set_ylabel('Regulons')
g.ax_col_colors.set_yticks([0.5])
g.ax_col_colors.set_yticklabels(['Cell Type'])
g.cax.set_visible(False)
g.fig.savefig(os.path.join(FIGURES_FOLDERNAME, 'clustermap - GSE115978.png'), format='png')

# Save results to excel

bin_mtx_clustered = bin_mtx.T.copy()
bin_mtx_clustered.rename(columns=df_annotations.set_index('cell_id')['cell_type'].to_dict(), inplace=True) 
bin_mtx_clustered.iloc[g.dendrogram_row.reordered_ind, g.dendrogram_col.reordered_ind].to_excel(os.path.join(RESULTS_FOLDERNAME, 'GSE115978 - Binarized regulon activity.xlsx'))




# ----------- STEP 6: Creation of adata --------------- #
from pyscenic.export import add_scenic_metadata

add_scenic_metadata(adata, auc_mtx, regulons)



import pandas as pd

from os import listdir
import pyarrow
from collections import defaultdict

directory = "/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/databases"



cis_target_data = {}
for file in listdir(directory):
    if file.endswith(".feather"):
        cis_target_data[file] = pd.read_feather(os.path.join(directory,file))  








db = cis_target_data['hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather']

motif_df = db.iloc[:5,:]

def create_motif_interaction_dict(dictionary):

    motif_interactions = defaultdict(list)
    

    
    for index,row in motif_df.iterrows():
        motif_name = row["motifs"]
        gene_columns = motif_df.columns[:-1]

        for gene_name in gene_columns:
                    #Check if the value is non-empty. 
                    if not pd.isna(row[gene_name]):
                        motif_interactions[motif_name].append(gene_name)        
        return motif_interactions

motif_interactions = create_motif_interaction_dict(directory)



from scipy.io import mmwrite
from scipy.sparse import csr_matrix
common_genes_matrix = adata[adata.obs['subtype'].isin(['ER','HER2']), adata.var_names.isin(common_genes)].X

mmwrite('common_genes.mtx',common_genes_matrix)

adata[adata.obs['subtype'].isin(['ER','HER2']),:].obs_names.to_series().to_csv('cells.txt', index=False, header=False)
adata[:,adata.var_names.isin(common_genes)].var_names.to_series().to_csv('genes.txt', index=False, header=False)



# Crear la matriz de expresión 
matriz = mmread('common_genes.mtx')

ADJACENCIES_FNAME = os.path.join(base_dir, "adjacencies.tsv")
MODULES_FNAME = os.path.join(base_dir, "modules.p")


gene_names = pd.read_csv("genes.txt", header=None, sep="\t",names=["Gene"])["Gene"].tolist()

cell_names = pd.read_csv("cells.txt", header=None, names=["Cell"])["Cell"].tolist()


expression_data = pd.DataFrame.sparse.from_spmatrix(
    matriz,
    index=cell_names,      # Nombres de las células
    columns=gene_names     # Nombres de los genes
)


expression_data.to_csv("./grn_files/expression_data.csv")

tf_names = load_tf_names("TF_names_v_1.01.txt")
tf_list_filtered = [tf for tf in tf_names if tf in expression_data.columns]

network = grnboost2(expression_data=expression_data,
                    tf_names=tf_names)

network.to_csv(ADJACENCIES_FNAME, index=False, sep='\t')
adjacencies = pd.read_csv(ADJACENCIES_FNAME, sep='\t')

# Derive potential regulomes from these co-expression modules
import pandas as pd
adjacencies = pd.read_table('/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/grn_files/adjacencies1.tsv', sep='\t')

expression_data = pd.read_csv("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929_176078/alltypes/grn_files/expression_data.csv", index_col=0)

modules = list(modules_from_adjacencies(adjacencies, expression_data))
import pickle
with open("modules.p", "wb") as f:
    pickle.dump(modules, f)



genes_in_adjacencies = set(adjacencies["target"].unique())
print(f"Number of unique genes in adjacencies: {len(genes_in_adjacencies)}")

# Obtener los nombres de las columnas de la matriz de expresión
genes_in_expression = set(expression_data.columns)
print(f"Number of genes in expression matrix: {len(genes_in_expression)}")

# Identificar genes que están en adjacencies pero no en expression_data
missing_genes = genes_in_adjacencies - genes_in_expression



with open(MODULES_FNAME, 'wb') as f:
    pickle.dump(modules, f)

with open(MODULES_FNAME, 'rb') as f:
    modules = pickle.load(f)

from pyscenic.utils import modules_to_df
modules_df = modules_to_df(modules)

modules_df.to_csv("modules.csv", index=False, sep="\t")



# Phase II: Prune modules for targets with cis regulatory footprints (aka RcisTarget)



import os
import glob
from pyscenic.ranking import FeatherRankingDatabase

# Define la ruta donde están las bases de datos


DATABASES_GLOB = "/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929_176078/alltypes/grn_files/*.feather"

# Función para extraer el nombre base del archivo
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]

# Encuentra todos los archivos .feather en la ruta
db_fnames = glob.glob(DATABASES_GLOB)

# Carga las bases de datos en objetos FeatherRankingDatabase
dbs = [FeatherRankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

# Imprime los objetos cargados para verificar
print(dbs)




import pandas as pd
import scanpy as sc

adata = sc.read_h5ad("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/data/BC_All_WithoutNormal/adata_integrated.h5ad")
results = pd.read_csv("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/grn_files/aucell_results.csv", index_col=0)




adyacencias = pd.read_csv("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/grn_files/adjacencies1.tsv", sep="\t")
enpp1_adjacencies = adyacencias[adyacencias['target'] == 'ENPP1']
enpp1_tfs = enpp1_adjacencies['TF'].unique().tolist()

relevant_tfs = [tf + "(+)" for tf in enpp1_tfs]
available_tfs = [tf for tf in relevant_tfs if tf in results.columns]


from sklearn.cluster import KMeans

cell_cycling = adata.obs_names[(adata.obs['GenAnno']== "Cycling Epithelial")&(adata.obs['subtype'].isin(['ER','HER2']))]

results_cycling = results[results.index.isin(cell_cycling)]

auc_data = results_cycling[available_tfs].values

# Realizar clusterización con K-means (e.g., 2 clusters para empezar)
n_clusters = 2
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
clusters = kmeans.fit_predict(auc_data)
results_cycling['Cluster'] = clusters
results_cycling['Cluster'].value_counts()


adata_cycling = adata[(adata.obs['GenAnno']== "Cycling Epithelial")&(adata.obs['subtype'].isin(['ER','HER2'])),:]


adata_cycling.obs['Cluster'] = results_cycling['Cluster'].astype('category')


sc.pl.embedding(adata_cycling, basis='X_scvi_MDE', color='Cluster', ncols=1, frameon=False)

sc.pl.embedding(adata_cycling, basis='X_scvi_MDE', color='ENPP1', ncols=1, frameon=False)

sc.pl.violin(adata_cycling, keys='ENPP1', groupby='Cluster')

aucell_results = pd.read_csv("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/grn_files/aucell_results.csv", index_col=0)


tfs = aucell_results.columns.values
threshold = 0.05
fraction_active = (aucell_results[tfs] > threshold).mean().sort_values(ascending=False)

top_regulons = fraction_active[fraction_active > 0.20].index

import matplotlib.pyplot as plt
# Crear gráficos de violín para estos regulones
fig, axes = plt.subplots(len(top_regulons), 1, figsize=(8, len(top_regulons) * 4), sharex=True)
for i, regulon in enumerate(top_regulons):
    axes[i].violinplot(aucell_results[regulon], showmeans=True)
    axes[i].set_title(f"Activity Distribution: {regulon}")
    axes[i].set_ylabel("AUC Score")
    axes[i].grid(axis='y')

plt.tight_layout()
plt.show()

for regulon in ['SPDEF(+)', 'AR(+)', 'TFAP2B(+)']:
    adata_fil.obs[regulon] = aucell_results[regulon]

sc.pl.embedding(
    adata_fil[adata_fil.obs['GenAnno'] == "Cycling Epithelial",:], basis='X_scvi_MDE', 
    color=['SPDEF(+)', 'AR(+)', 'TFAP2B(+)','ENPP1'], 
    ncols=1, frameon=False
)


# Seleccionar los regulones de interés
regulons_of_interest = ['SPDEF', 'AR', 'TFAP2B']

# Filtrar genes objetivo de estos regulones
target_genes = adyacencias[adyacencias['TF'].isin(regulons_of_interest)]['target'].unique()

# Mostrar los genes asociados a los regulones de interés
target_genes_list = target_genes.tolist()

from gprofiler import GProfiler
gp = GProfiler(return_dataframe=True)
enrichment_results = gp.profile(organism='hsapiens', query=target_genes_list,no_evidences=False)

# Filtrar los resultados más relevantes (e.g., FDR < 0.05)
significant_enrichments = enrichment_results[enrichment_results['p_value'] < 0.01].sort_values(['source','intersection_size'], ascending=False)


significant_enrichments.to_csv("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/grn_files/significant_enrichments.csv", index=False)



import numpy as np
np.mean(results, axis=0).sort_values(ascending=False).head(10)  






cell_cycle_genes = [x.strip() for x in open('/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/CellCycle/regev_lab_cell_cycle_genes.txt')]
print(len(cell_cycle_genes))

# Split into 2 lists
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
print(len(cell_cycle_genes))

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pl.violin(adata, ['S_score', 'G2M_score'],
             jitter=0.4, groupby = 'GenAnno', rotation=90)




import loompy as lp



row_attrs = { 
    "Gene": np.array(adata.var.index) ,
}
col_attrs = { 
    "CellID":  np.array(adata.obs.index) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}

#f_loom_path_unfilt = "/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/grn_files/loomUnfiltered.loom"
#lp.create( f_loom_path_unfilt, adata.X.transpose(), row_attrs, col_attrs )


lf = lp.connect("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/grn_files/loomUnfiltered.loom", mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)


regulons = pd.read_csv("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/GRN-/grn_files/regulons.csv",       skiprows=2,  # Omite las dos primeras filas
    names=["TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", 
           "OrthologousIdentity", "Annotation", "Context", "TargetGenes", "RankAtMax"],
    usecols=["AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity", 
             "Annotation", "Context", "TargetGenes", "RankAtMax"]
)
regulons = regulons.dropna(how="all")

from pyscenic.export import add_scenic_metadata

adata_fil = adata[adata.obs['subtype'].isin(['ER','HER2']),:]
add_scenic_metadata(adata_fil,aucell_results,regulons)