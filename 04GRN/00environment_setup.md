From: [pyscenic github page](https://github.com/aertslab/pySCENIC/blob/master/docs/installation.rst)



It is important  to load Squashfs if it is not by default: module load Squashfs/4.3
Singularity also needed (loaded by default in my case)


[Singularity/Apptainer images can be build from the Docker Hub image as source:](https://github.com/aertslab/pySCENIC/blob/master/docs/installation.rst#singularityapptainer)

# pySCENIC CLI version.
singularity build aertslab-pyscenic-0.12.1.sif docker://aertslab/pyscenic:0.12.1
apptainer build aertslab-pyscenic-0.12.1.sif docker://aertslab/pyscenic:0.12.1

# pySCENIC CLI version + ipython kernel + scanpy.
singularity build aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif docker://aertslab/pyscenic_scanpy:0.12.1_1.9.1
apptainer build aertslab-pyscenic-0.12.1-1.9.1.sif docker://aertslab/pyscenic_scanpy:0.12.1_1.9.1


In my case I run: singularity build aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif docker://aertslab/pyscenic_scanpy:0.12.1_1.9.1


To run pySCENIC with Singularity/Apptainer, the usage is very similar to that of Docker/Podman.



singularity run aertslab-pyscenic-0.12.1.sif \
    pyscenic grn \
        -B /data:/data
        --num_workers 6 \
        -o /data/expr_mat.adjacencies.tsv \
        /data/expr_mat.tsv \
        /data/allTFs_hg38.txt

    



----

Vignettes: 
https://github.com/aertslab/pySCENIC/blob/master/notebooks/pySCENIC%20-%20Full%20pipeline.ipynb
https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb



# WorkFlow 

## 1. Infer Networks using arboreto grnboost2: NetworkInference.ipynb
It is done in a jupyter notebook to leverage Dask paralelization. 

The extra file grn_net may potentially be used, but I have problems with Dask inside the cluster. 

In this section, we extract pairs of TF-gene that are co-expressed across the whole sample. 
The degree of coexpression is measured as "importance"

Results are saved into 01Networks folder

## 2. Create loom files with the expression data 

The next steps would leverage the cli version of pyscenic. 
A loom file would be the count matrix input. 

Results are saved into 02ExpMatrix

## 3. Infer Modules from the inferred networks
Cis target motifs (coming from .feather databases) are used 
The infer_modules.sh script is used, that leverages pyscenic ctx command. 
It creates motifs that are saved into 03Motifs 

From this motifs regulons are easily derived with the df2regulons command of pyscenic in a jupyter notebook

## 4. Infer the regulon activity for each cell

Results are saved into 04AUCell



--------------

# Analysis 

Relevant cell populations are obtained from cell phone DB analysis. 
The objective is to see what is different (unique, characteristic) and what is shared between twose populations 
across the three tumor types, integrating gene regulatory network data with cell-cell comunication data. 

07/01/25 00:56

Viendo los factores de transcripción activos por tipo celular me doy cuenta de que sólo aparecen poblaciones que son mas o menos pequeñas y con escasa variablilidad, como endoteliales o células plasmáticas. 
Por lo que antes de correr aucell anotar extensivamente las poblaciones de interés. 

servirá luego en cellphonedb



singularity exec --cleanenv --bind /home/igarzonalva/Proyecto_SC_TNBC/GSE161529/InferCNV:/data  my_infercnv_latest.sif jupyter notebook --KernelSpecManager.whitelist="['ir']"  --ip=0.0.0.0 --port=8888 --no-browser

