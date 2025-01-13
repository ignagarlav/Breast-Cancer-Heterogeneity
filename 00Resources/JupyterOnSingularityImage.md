
Primero module load Singularity/3.8.5 ––> luego Singularity/3.8.5-GCC-10.2.0
Luego  module load Squashfs/4.3 (importante para hacer el build)
Por último singularity pull nombre_sif.sif  docker://ignagaralv/scarches-iga:latest


Local SSH to HPC Login Node

1) Step 1
On your local machine, run:
ssh -L 9999:localhost:9999 igarzonalva@hpc-2.unav.es
This forwards local port 9999 to port 9999 on the HPC login node.

2) Step 2: Allocate a Node with salloc

From the HPC login node request an interactive session:
salloc --nodes=1 --ntasks=1 --time=02:00:00
Wait until a compute node is allocated (e.g., nodo07).

3) Step3: SSH to the Allocated Node (Second Tunnel)

Still inside the HPC login node, (in the terminal of step 1) run:

ssh -L 9999:localhost:8888 node123
Here, we’re linking port 9999 on the login node to port 8888 on node123.

4) Start Jupyter in Singularity (or through conda)

To start jupyter: 

singularity exec --bind    /home/igarzonalva/Proyecto_SC_TNBC/singularity_images/scarches_iga_latest.sif jupyter notebook --KernelSpecManager.whitelist=
"['scarches']"  --ip=0.0.0.0 --port=8888 --no-browser

Through conda would be: 
conda activate env
env jupyter notebook --no-browser --port=8888

Open in local computer: localhost:9999/token
(9999 or the port fordwarded in step 1)