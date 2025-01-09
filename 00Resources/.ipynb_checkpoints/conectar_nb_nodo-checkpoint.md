Conectarme desde local al root node del cluster: 
ssh -L 8888:localhost: 9999 hpclogin
bind [127.0.0.11:8888: Address already in use
Last login: Tue Jan
7 22:54:21 2025 from 10.71.192.34

Una vez dentro, conectarme al nodo correspondiente. 
(base) ssh -L 9999: localhost: 8888 igarzonalva@nodo08
Last login: Tue Jan 7 22:51:47 2025 from hpc-2

Para ver en qué nodo tengo la alocación, correr echo $HOSTNAME después de salloc --nodes=1 --ntasks=8 --mem=32G --time=5:00:00 --partition=short

