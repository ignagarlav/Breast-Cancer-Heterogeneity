###################### ver el modelo de referencia 


sc.pl.embedding(
    adata,
    basis="X_scanvi_MDE",
    color=["batch","leiden","celltypist_labels_Immune_All_High","GenAnno"],
    frameon=False,
    ncols=1)

sc.pl.embedding(
    adata,
    basis="X_scanvi_MDE",
    color=["ENPP1","GenAnno", "subtype"],
    frameon=False,
    ncols=1)
