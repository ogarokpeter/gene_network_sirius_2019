# RUN WITH /usr/bin/python3 minet.py (python 3.6)

import rpy2
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as rn
import numpy as np
import xml.dom.minidom as md
from sklearn.metrics import roc_curve, auc

datadirname = "/home/user/Sirius/gene_network_sirius_2019/Data"
datafilename = datadirname + "/{0}/{0}_data.txt"
graphfilename = datadirname + "/{0}/{0}_graph.xml"
matricesdirname = "/home/user/Sirius/gene_network_sirius_2019/Matrices"
predictedfilename = matricesdirname + "/{1}_{0}_predicted.txt"
truefilename = matricesdirname + "/{1}_{0}_true.txt"
datalist = ['exps_10', 'exps_10_2', 'exps_10_bgr', 'exps_50', 'exps_50_2', 'exps_50_bgr', 'exps_100', 'exps_100_2', 'exps_100_bgr', 'genes_200_exps_10_bgr', 'genes_400_exps_10_bgr', 'genes_600_exps_10_bgr', 'genes_700_exps_10_bgr', 'genes_1000_exps_10_bgr']
algolist = ['clr', 'aracne', 'mrnet', 'mrnetb']

def run_minet(filename, algo):

    rn.activate()

    code = """library(minet)
    filename <- '""" + filename + """'
    first <- readLines(filename, n=1)
    names <- strsplit(first, '\t')
    names <- unlist(names, use.names=FALSE)
    d <- read.table(filename, skip=1, col.names = names)

    mim <- build.mim(d, estimator = "mi.empirical", disc = "equalfreq")

    weight_adjacency_matrix <- minet(mim, method='""" + algo + """', estimator="mi.empirical", disc="equalfreq");

    weight_adjacency_matrix;
    """

    f = ro.r(code)

    weight_adjacency_matrix = np.array(f)

    long_array = weight_adjacency_matrix[np.triu_indices(weight_adjacency_matrix.shape[0])]
    return long_array, weight_adjacency_matrix



def xml_graph_to_adjacency_matrix(filename):
    dom = md.parse(filename)

    # print(dom.toprettyxml())

    nodes = dom.getElementsByTagName("Node")
    ids = [int(a.getAttribute('id')) for a in nodes]
    # parserint(ids)
    adjacency_matrix = np.zeros(shape=(len(ids), len(ids)))
    edges = dom.getElementsByTagName("Edge")
    for e in edges:
        source = int(e.getElementsByTagName('from')[0].firstChild.nodeValue)
        target = int(e.getElementsByTagName('to')[0].firstChild.nodeValue)
        adjacency_matrix[source][target] = 1
        adjacency_matrix[target][source] = 1
    return adjacency_matrix

results = np.zeros(shape=(len(datalist), len(algolist)))

for i, dataname in enumerate(datalist):

    for j, algo in enumerate(algolist):

        long_array, predicted_matrix = run_minet(datafilename.format(dataname), algo)
        true_matrix = xml_graph_to_adjacency_matrix(graphfilename.format(dataname))
        true_array = true_matrix[np.triu_indices(true_matrix.shape[0])]

        np.savetxt(predictedfilename.format(dataname, algo), predicted_matrix, delimiter='\t')
        np.savetxt(truefilename.format(dataname, algo), true_matrix, delimiter='\t')

        roc_auc = 0
        try:
            fpr, tpr, thresholds = roc_curve(true_array, long_array)
            roc_auc = auc(fpr, tpr)
        except:
            print("error", dataname, algo)
        results[i][j] = roc_auc
        print("done", dataname, algo, results[i][j])

print(results)

