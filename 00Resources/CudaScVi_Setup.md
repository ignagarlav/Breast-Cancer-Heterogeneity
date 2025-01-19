[Source](https://github.com/scverse/scvi-tools/issues/2333)

frac2738 on Jan 8, 2024 · edited by frac2738 said

```
conda create -n scvi-env python=3.9
conda activate scvi-env

# check cuda version
nvidia-smi
+---------------------------------------------------------------------------------------+
| NVIDIA-SMI 545.23.08              Driver Version: 545.23.08    CUDA Version: 12.3     |
|-----------------------------------------+----------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |         Memory-Usage | GPU-Util  Compute M. |
|                                         |                      |               MIG M. |
|=========================================+======================+======================|
|   0  Quadro P2000                   On  | 00000000:09:00.0  On |                  N/A |
| 65%   69C    P0              52W /  75W |   1165MiB /  5120MiB |     87%      Default |
|                                         |                      |                  N/A |
+-----------------------------------------+----------------------+----------------------+

pip install pandas PyYAML scipy
conda install jaxlib=*=*cuda* jax cuda-nvcc -c conda-forge -c nvidia

# some things I tried but all failed
# pip install --upgrade "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
# pip install -U --pre jaxlib -f https://storage.googleapis.com/jax-releases/jaxlib_nightly_cuda12_releases.html
# pip3 install torch torchvision torchaudio

pip install scvi-tools
pip install 'scanpy[leiden]'

```


What has specifically worked for me: 

````
conda create -n scvi-cuda_env3 python=3.12
conda activate scvi-cuda_env3
nvidia-smi
pip install pandas PyYAML scipy
conda install jaxlib=*=*cuda* jax cuda-nvcc -c conda-forge -c nvidia
pip install scvi-tools
pip install 'scanpy[leiden]'

