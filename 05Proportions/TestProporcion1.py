## Comprobar cómo varía la proporción de células en los distintos tumores
import scanpy as sc
import scProportionTest as pt
import os 
import pandas as pd
import sys

def main():
    
    group1 = sys.argv[1]      
    group2 = sys.argv[2]      
    data_dir = os.getenv("DATA_DIR")  
    output_dir = os.getenv("OUTPUT_DIR")

    adata_dir = os.path.join(data_dir,'adata_scanvi_predictions.h5ad')
    print(f"Reading adata from {adata_dir}")
    adata = sc.read_h5ad(adata_dir)

    print(f"Calculating proportions {group2} vs the reference {group1}")
    results = pt.permutation_test(adata,
                                   group1,
                                   group2,
                                   group_col='subtype',
                                   cell_type_col='predicted_labels',
                                   nperm=10000,
                                   alpha=0.05,
                                   n_bootstrap=10000,
                                   verbose=True)
    print(f"Saving results to {output_dir}")
    results.to_csv(os.path.join(output_dir,f"scPropTest_{group2}_vs_ref_{group1}_results.csv"), index=False)
    

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error inesperado: {e}")