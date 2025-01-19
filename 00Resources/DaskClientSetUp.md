
Set up dask 

portdash = 40748
cluster = SLURMCluster(queue = "short", cores=16, processes=1, 
                   memory="16GB", walltime="05:00:00",
                   scheduler_options={"dashboard_address": f":{portdash}", 'host':'nodo05'})
cluster.scale(2)
client = Client(cluster)

this is the output: 

Worker: SLURMCluster-1
Comm: tcp://172.16.2.7:46182	Total threads: 16
Dashboard: http://172.16.2.7:39084/status	Memory: 14.90 GiB
Nanny: tcp://172.16.2.7:37627	
Local directory: /tmp/dask-scratch-space/worker-ssuapi85
Tasks executing:	Tasks in memory:
Tasks ready:	Tasks in flight:
CPU usage: 13.2%	Last seen: Just now
Memory usage: 128.15 MiB	Spilled bytes: 0 B



3. Create a Second Set of Tunnels for Port 40748
Now you want to “chain” local → HPC login → nodo05:40748.
Open a new local terminal (or you can combine them, but we’ll do it step by step):

From your local machine to HPC login node:

bash
Copiar
ssh -N -L 40748:localhost:40748 hpclogin
This sets local port 40748 to forward to hpclogin:40748.
Leave this running.
From HPC login node to nodo05:

bash
Copiar
ssh -N -L 40748:localhost:40748 igarzonalva@nodo05
This sets HPC login node’s :40748 to forward to nodo05:40748.
So any request that hits “login_node:40748” is relayed to “nodo05:40748.”
Put another way:

Local port 40748 → HPC login node:40748
HPC login node:40748 → nodo05:40748
