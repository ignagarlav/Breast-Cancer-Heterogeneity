
Primero module load Singularity/3.8.5 ––> luego Singularity/3.8.5-GCC-10.2.0
Luego  module load Squashfs/4.3 (importante para hacer el build)
Por último singularity pull nombre_sif.sif  docker://ignagaralv/scarches-iga:latest


Local SSH to HPC Login Node

On your local machine, run:
bash
Copiar código
ssh -L 9999:localhost:9999 igarzonalva@hpc-2.unav.es
This forwards local port 9999 to port 9999 on the HPC login node.
Allocate a Node with salloc

From the HPC login node (after the above SSH), request an interactive session:
bash
salloc --nodes=1 --ntasks=1 --time=02:00:00
Wait until a compute node is allocated (e.g., node123).
SSH to the Allocated Node (Second Tunnel)

Still inside the HPC login node, run:
bash
Copiar código
ssh -L 9999:localhost:8888 node123
Here, we’re linking port 9999 on the login node to port 8888 on node123.
Start Jupyter in Singularity

To start jupyter: 

(base) [igarzonalva@nodo06 Proyecto_SC_TNBC]$ singularity exec --bind /home/igarzonalva/Proyecto_SC_TNBC/GSE161529:/data   /home/igarzonalva/Proyecto_SC_TNBC/singularity_images/scarches_iga_latest.sif  conda run -n scarches jupyter notebook     --KernelSpecManager.whitelist="['scarches']"  --ip=0.0.0.0 --port=8888 --no-browser

Open in local computer: localhost:9999/token

