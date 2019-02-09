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
datalist = ['exps_10_bgr', 'exps_50_bgr', 'exps_100_bgr', 'genes_200_exps_10_bgr', 'genes_400_exps_10_bgr', 'genes_500_exps_10_bgr']
# datalist = ['genes_700_exps_10_bgr', 'genes_1000_exps_10_bgr']
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
    # return long_array



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


def adjacency_matrix_to_array_of_degrees(adjacency_matrix):
    degrees = [[] for i in range(len(adjacency_matrix))]
    for i in range(len(adjacency_matrix)):
        deg_i = 0
        for j in range(len(adjacency_matrix[i])):
            if adjacency_matrix[i][j] == 1:
                deg_i += 1
        degrees[deg_i].append(i)
    return degrees


results = np.zeros(shape=(len(datalist), len(algolist)))

for ii, dataname in enumerate(datalist):

    for jj, algo in enumerate(algolist):

        matrix = run_minet(datafilename.format(dataname), algo)
        true_matrix = xml_graph_to_adjacency_matrix(graphfilename.format(dataname))

        # print(matrix, true_matrix)

        degrees = adjacency_matrix_to_array_of_degrees(true_matrix)

        # print(degrees)

        aucs_mean = np.zeros(shape=(len(degrees)))
        aucs_var = np.zeros(shape=(len(degrees)))
        for i, deg in enumerate(degrees):
            if len(deg) != 0:
                aucs_d = np.zeros(shape=(len(deg)))
                for j, ind in enumerate(deg):
                    true_array = true_matrix[ind]
                    array = matrix[ind]
                    # print(array, true_array)
                    roc_auc = 0
                    try:
                        fpr, tpr, thresholds = roc_curve(true_array, array)
                        roc_auc = auc(fpr, tpr)
                    except:
                        print("error", dataname, algo)
                    aucs_d[j] = roc_auc
                aucs_mean[i] = aucs_d.mean()
                aucs_var[i] = aucs_d.var()
        with open("result_analysis.txt", "a") as f:
            print(dataname, algo)
            f.write(dataname + " " + algo + ":" + '\n')
            for i in range(len(aucs_mean)):
                if aucs_mean[i] != 0:
                    f.write(str(i) + " " + str(aucs_mean[i]) + " " + str(aucs_var[i]) + '\n')
