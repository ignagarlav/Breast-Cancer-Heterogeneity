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

# WorkFlow 

1. Infer Networks using arboreto grnboost2: NetworkInference.ipynb
It is done in a jupyter notebook to leverage Dask paralelization. 

The extra file grn_net may potentially be used, but I have problems with Dask inside the cluster. 

2. Create loom files with the expression data (would serve )

3. Infer Modules from the inferred networks
Cis target motifs (coming from .feather databases) are used 
I have used the infer_modules.sh script, that leverages pyscenic ctx command

