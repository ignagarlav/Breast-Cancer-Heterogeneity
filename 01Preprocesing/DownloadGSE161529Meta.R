BiocManager::install("GEOquery")
library(tidyverse)
library(GEOquery)
gse_object <- getGEO("GSE161529", GSEMatrix = TRUE)
sample_metadata <- pData(gse_object[[1]])


metadata <- sample_metadata %>% 
  dplyr::select(c("geo_accession", "characteristics_ch1.1","cancer type:ch1","cell population:ch1","supplementary_file_1","supplementary_file_2"))

row.names(metadata) <- NULL


metadata <- metadata %>% 
  dplyr::rename(
    "gender" = characteristics_ch1.1,
    "cancer_type" = `cancer type:ch1`,
    "cell_population"= `cell population:ch1`,
    "barcodes_file" = `supplementary_file_1`,
    "matrix_file" = supplementary_file_2
  )



metadata$gender <-  gsub('gender: ', "",metadata$gender)
metadata$cell_population <- gsub("Involved lymph node", "LN", metadata$cell_population)
metadata$cancer_type <- gsub(pattern = " |-","_",metadata$cancer_type)
metadata$cancer_type <- gsub(pattern = "\\+","",metadata$cancer_type)
metadata$cancer_type <- gsub("Triple_negative_tumour", "TNBC",metadata$cancer_type)
metadata$cancer_type <- gsub("Triple_negative_BRCA1_tumour", "TNBC_BRCA",metadata$cancer_type)
metadata$cancer_type <- gsub("_tumour","",metadata$cancer_type)
metadata$features_file <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161529/suppl/GSE161529_features.tsv.gz"


################################## Metadata of interest ##################################

meta_of_interest <- metadata %>% 
  filter(cell_population != "LN",
         !(cancer_type %in% c("Normal","BRCA1_pre_neoplastic","PR")),
           gender != "Male")
         

tmp_file <- file.path(tempdir(), "my_data.csv")
write.table(meta_of_interest, file = tmp_file, sep = "\t", row.names = F, col.names = T)
system(paste("scp", tmp_file, "hpclogin:/home/igarzonalva/Proyecto_SC_TNBC/GSE161529/01_Preprocessing/meta_of_interest.txt"))
file.remove(tmp_file)
