INSTALLATION & LAUNCHING:
The software is provided as an executable jar file. No installation is required, just place these files in a folder where you have write access to. You can launch the software either by double-clicking it or by using the command line:
java -jar SynTReN.jar

GENERATED FILES:
all generated files are saved in ./data/results
and apply to the following template:
nn<v>_nbgr<w>_hop<x>_bionoise<y>_expnoise<z>_(neigh|clust)Add_(dataset|external|network).(txt|sif|xml)

v: the number of nodes of the 'foreground' network
w: the number of nodes of the 'background' network
x: the probability that an interaction is 'complex'
y: the amount of biological noise
z: the amount of experimental noise
(neigh|clust): neighbour addition or cluster addition method
(dataset|external|network): the dataset (tab-delimited), the set of external nodes and complete description of the network (both in xml and sif)

MOST RECENT SOFTWARE VERSION:
http://homes.esat.kuleuven.be/~kmarchal/SyNTReN

REQUIREMENTS:
Java 5