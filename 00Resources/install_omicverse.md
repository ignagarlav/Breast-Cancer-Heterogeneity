# Pull the docker image in "sandbox"

Install the non-cpu version as it is 1.5G (compared with 4.5G) and our HPC doesn't has GPU. 

```
singularity build --sandbox omicverse_sandbox docker://starlitnightly/omicverse:py310_cpu_latest
```


# Open the sandbox in interactive mode 

```
cd directorio_madre_del_sandbox
singularity shell --writable omicverse_sandbox
cd omicverse_sandbox # to move inside the folder in within the shell session
pip install --no-index --find-links=/packages torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 && pip install --no-index --find-links=/packages torch_geometric && pip install --no-index --find-links=/packages pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv && pip install --no-index --find-links=/packages omicverse && pip install --no-index --find-links=/packages dgl && pip install --no-index --find-links=/packages tangram pot tangram cvxpy libpysal gudhi openai patsy combat pydeseq2==0.4.1 pymde opencv-python scikit-image && pip install --no-index --find-links=/packages harmonypy intervaltree fbpca mofax metatime s-gd2 mellon scvelo cellrank && pip install --no-index --find-links=/packages numpy==1.23.5 llvmlite==0.41.1 numba==0.58.1 && jupyter lab --ip=0.0.0.0 --port=8899 --no-browser --allow-root 
exit
```

# Build final .sif image from the sandbox 

You should be in the parent directory of the sandbox
```
(base) [igarzonalva@hpc-2 singularity_images]$ ls
iga_scvi_env.sif         iga_scvi_env_v1.0.sif  omicverse_sandbox        scarches.sif
iga_scvi_env_v1.0.1.sif  images                 scarches_iga_latest.sif
```

Load **Squashfs** (*REALLY IMPORTANT*)

```
module load Squashfs/4.3
```

**Build the final .sif file**
```
singularity build omicverse_gpu_custom.sif omicverse_gpu_sandbox

# Run jupyter 

```
singularity exec omicverse_custom.sif jupyter lab --no-browser --port=8899
```
```

# Why is all this needed

In the [docker webpage](https://hub.docker.com/layers/starlitnightly/omicverse/py310_cpu_latest/images/sha256-ffa2a155024e73982a0f02dabd14e191fa1164632045f090a465e39adf4553c9)

in step 19 of the build: 

```
CMD ["sh" "-c" "pip install --no-index --find-links=/packages torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 && pip install --no-index --find-links=/packages torch_geometric && pip install --no-index --find-links=/packages pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv && pip install --no-index --find-links=/packages omicverse && pip install --no-index --find-links=/packages dgl && pip install --no-index --find-links=/packages tangram pot tangram cvxpy libpysal gudhi openai patsy combat pydeseq2==0.4.1 pymde opencv-python scikit-image && pip install --no-index --find-links=/packages harmonypy intervaltree fbpca mofax metatime s-gd2 mellon scvelo cellrank && pip install --no-index --find-links=/packages numpy==1.23.5 llvmlite==0.41.1 numba==0.58.1 && jupyter lab --ip=0.0.0.0 --port=8899 --no-browser --allow-root"]
```


the command is run through CMD, through which cached versions of the packages are installed (they are located in .whl)
Nonetheless this command is NOT run when using singularity, and we should do it by hand. 


-----------

[omicverse installation guide](https://starlitnightly.github.io/omicverse/Installation_guild/#docker)
[omicverse docker images](https://hub.docker.com/r/starlitnightly/omicverse/tags)
[omicverse github](https://github.com/Starlitnightly/omicverse)