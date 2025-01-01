## Comprobar cómo varía la proporción de células en los distintos tumores
import scanpy as sc
import scProportionTest as pt
import os 
import pandas as pd

def main():
    data_dir = os.getenv("DATA_DIR", "./data")  
    output_dir = os.getenv("OUTPUT_DIR", "./data")

    adata = sc.read_h5ad(os.path.join(data_dir,'adata.h5ad'))
    results = pt.permutation_test(adata,
                                   'TNBC',
                                   'ER',
                                   group_col='subtype',
                                   cell_type_col='GenAnno',
                                   nperm=10000,
                                   alpha=0.05,
                                   n_bootstrap=10000,
                                   verbose=True)
    results.to_csv(os.path.join(output_dir,"scPropTest_TNBCER_results.csv"), index=False)
    results1 = pt.permutation_test(adata,
                               'TNBC',
                               'HER2',
                               group_col='subtype',
                               cell_type_col='GenAnno',
                               nperm=10000,
                               alpha=0.05,
                               n_bootstrap=10000,
                               verbose=True)
    results1.to_csv(os.path.join(output_dir,"scPropTest_TNBCHER2_results.csv"), index=False)
    results2 = pt.permutation_test(adata,
                                'ER',
                                'HER2',
                                group_col='subtype',
                                cell_type_col='GenAnno',
                                nperm=10000,
                                alpha=0.05,
                                n_bootstrap=10000,
                                verbose=True)

    # Suponiendo que los resultados están en un DataFrame llamado `results`
    results2.to_csv(os.path.join(output_dir,"scPropTest_ERHER2_results.csv"), index=False)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error inesperado: {e}")