'''Use scvi-env to run this script'''


import scProportionTest as pt
import pandas as pd

results = pd.read_csv("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/PropTest/proportiontest/scPropTest_TNBCER_results.csv")

results1 = pd.read_csv("/Users/joseignaciogarzonalvarez/proyectosPython/singlecelltnbc/GSE162929/alltypes/PropTest/proportiontest/scPropTest_TNBCHER2_results.csv")

# En ambos resultados, TNBC es la referencia 

plot = pt.point_range_plot(results, figsize=(6,4), fold_difference = 1,alpha = 0.001,ascending = False,plot_title = "ER vs TNBC ")
plot = pt.point_range_plot(results1, figsize=(6,4),fold_difference = 1,alpha = 0.001,ascending = False, plot_title = "HER2 vs TNBC")

'''
AntigenMyeloid_TNBC          487
B Cells_ER                   475
Endothelial_TNBC             474
Plasmablasts_ER              240
Mast/Basophiph_HER2          142
Mast/Basophiph_ER            128
Mast/Basophiph_TNBC          113
AntigenMyeloid_HER2          104
AntigenMyeloid_ER             53


Cycling Epithelial_HER2     2249


(adata.obs['GenAnno'].value_counts() / adata.obs['GenAnno'].count())*100

GenAnno
Epithelial            51.090803 %
T Cells               13.455971 %
Cycling Epithelial    11.470664 %
Myeloid               10.072283 %
Fibroblasts            6.985996 %
Plasmablasts           3.200520 %
B Cells                1.744365 %
Endothelial            1.305156 %
AntigenMyeloid         0.422797 %
Mast/Basophiph         0.251446 %
Unknown                0.000000 %
Name: count, dtype: float64


'''

