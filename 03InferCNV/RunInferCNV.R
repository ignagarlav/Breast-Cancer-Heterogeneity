library(infercnv)
library(rjags)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
metadata_file_dir <- args[1]      # First argument: Metadata file (annotations)
mat <- args[2]                    # Second argument: Counts matrix
gene_ordering <- args[3]          # Third argument: Gene ordering file
output_dir <- args[4]             # Fourth argument: Output directory

############## Run InferCNV ###########################################

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=mat,
                                    annotations_file=metadata_file_dir,
                                    delim="\t",
                                    gene_order_file=gene_ordering,
                                    ref_group_names=c('TCells','B Cells'))


options(scipen = 100)

infercnv_obj = infercnv::run(infercnv_obj,
                            cutoff=0.1,  
                            out_dir=output_dir,  
                            analysis_mode = "samples",
                            cluster_by_groups=T, 
                            plot_steps = F,
                            no_plot = T,
                            denoise=T,
                            HMM=T,
                            save_rds = TRUE,
                            save_final_rds = TRUE
)

