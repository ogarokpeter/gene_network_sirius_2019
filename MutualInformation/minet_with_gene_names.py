# RUN WITH /usr/bin/python3 minet.py (python 3.6)

import rpy2
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as rn
import numpy as np
import pandas as pd
import xml.dom.minidom as md
from sklearn.metrics import roc_curve, auc

datadirname = "/home/user/Sirius/gene_network_sirius_2019/Data"
datafilename = datadirname + "/{0}/{0}_data.txt"
graphfilename = datadirname + "/{0}/{0}_graph.xml"
matricesdirname = "/home/user/Sirius/gene_network_sirius_2019/Matrices_1"
predictedfilename = matricesdirname + "/{1}_{0}_predicted.txt"
truefilename = matricesdirname + "/{1}_{0}_true.txt"
# datalist = ['exps_10', 'exps_10_2', 'exps_10_bgr', 'exps_50', 'exps_50_2', 'exps_50_bgr', 'exps_100', 'exps_100_2', 'exps_100_bgr', 'genes_200_exps_10_bgr', 'genes_400_exps_10_bgr', 'genes_600_exps_10_bgr', 'genes_700_exps_10_bgr', 'genes_1000_exps_10_bgr']
# datalist = ['genes_200_exps_20_bgr', 'genes_200_exps_40_bgr', 'genes_400_exps_40_bgr', 'genes_400_exps_80_bgr', 'genes_500_exps_50_bgr', 'genes_500_exps_100_bgr']
datalist = ['genes_500_exps_10_bgr']
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

    return weight_adjacency_matrix

    # long_array = weight_adjacency_matrix[np.triu_indices(weight_adjacency_matrix.shape[0])]
    # return long_array, weight_adjacency_matrix



def xml_graph_to_adjacency_matrix(filename):
    dom = md.parse(filename)

    # print(dom.toprettyxml())

    nodes = dom.getElementsByTagName("Node")
    ids = [int(a.getAttribute('id')) for a in nodes]
    genes = [None for id in ids]
    for node in nodes:
        id = int(node.getAttribute('id'))
        name = str(node.getAttribute('name'))
        genes[id] = name
    # print(genes)
    adjacency_matrix = np.zeros(shape=(len(ids), len(ids)))
    edges = dom.getElementsByTagName("Edge")
    for e in edges:
        source = int(e.getElementsByTagName('from')[0].firstChild.nodeValue)
        target = int(e.getElementsByTagName('to')[0].firstChild.nodeValue)
        adjacency_matrix[source][target] = 1
        adjacency_matrix[target][source] = 1
    adjacency_df = pd.DataFrame(adjacency_matrix, index=genes, columns=genes)
    return adjacency_df, genes


results = np.zeros(shape=(len(datalist), len(algolist)))

for i, dataname in enumerate(datalist):

    for j, algo in enumerate(algolist):

        predicted_matrix = run_minet(datafilename.format(dataname), algo)
        true_df, genes = xml_graph_to_adjacency_matrix(graphfilename.format(dataname))
        predicted_df = pd.DataFrame(predicted_matrix, index=genes, columns=genes)

        predicted_df.to_csv(predictedfilename.format(dataname, algo), index=True, header=True, sep='\t')
        true_df.to_csv(truefilename.format(dataname, algo), index=True, header=True, sep='\t')

        true_array = true_df.values[np.triu_indices(true_df.values.shape[0])]
        predicted_array = predicted_df.values[np.triu_indices(predicted_df.values.shape[0])]

        roc_auc = 0
        try:
            fpr, tpr, thresholds = roc_curve(true_array, predicted_array)
            roc_auc = auc(fpr, tpr)
        except:
            print("error", dataname, algo)
        results[i][j] = roc_auc
        print("done", dataname, algo, results[i][j])
        
        # print("done", dataname, algo)

print(results)

