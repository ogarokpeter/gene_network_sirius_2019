# RUN WITH /usr/bin/python3 minet.py (python 3.6)

import sys
import rpy2
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as rn
import numpy as np
import xml.dom.minidom as md
from sklearn.metrics import roc_curve, auc
from scipy.stats import rankdata
import pandas as pd

def run_ra(filename, n_genes):

    rn.activate()

    code = """library(RankAggreg)
d <- read.table('""" + filename + """')
d <- as.matrix(d)
CES <- RankAggreg(d, """ + str(n_genes) + """, method="GA", maxIter=100, verbose=FALSE)
res <- CES$top.list;
res;
    """
    f = ro.r(code)

    res = list(f)
    
    return res


def make_ranked_lists(matrix):
    mat = matrix.copy()
    lst = [list(mat[i].sort_values(ascending=False).index) for i in mat]
    # print(lst)
    return lst


def reconsrtuct_df(lst, genes):
    gns = genes.tolist()
    n_genes = len(gns)
    # print(gns, n_genes)
    # print(lst)
    zero_data = np.zeros(shape=(n_genes,n_genes))
    d = pd.DataFrame(zero_data, columns=genes, index=genes)
    for i, gene in enumerate(gns):
        d_gene = np.zeros(shape=(n_genes))
        # print(lst[i])
        d = d.reindex(lst[i])
        # d.index = lst[i]
        for i in range(n_genes):
            d_gene[n_genes - i - 1] = 2 * (i + 1) / (n_genes * n_genes + n_genes)
        d[gene] = d_gene
        # print(d)
        # break
    # d = pd.DataFrame(np.array([0]))
    # print(d.index)
    # print(lst[0])
    d = d.reindex(genes)
    # print(d.index)
    # print(d)
    return d


def compute_aggregated_matrix(matrixfiles_num, matrixfiles, savematrixfile, saveresultfile):
    # matrixfiles_num = int(sys.argv[1])
    # matrixfiles = [sys.argv[i] for i in range(2, matrixfiles_num + 2)]
    # savematrixfile = sys.argv[matrixfiles_num + 2]
    # saveresultfile = sys.argv[matrixfiles_num + 3]
    matrices = [pd.read_csv(f, index_col=0, sep='\t') for f in matrixfiles]
    genes = matrices[0].index
    # print(genes)

    # print(matrices)
    sz = len(matrices[0])
    for matrix in matrices:
        assert len(matrix) == sz
    # # make_ranked_lists(matrices[0])
    # print(matrices[0])
    # for i in matrices[0]:
    #     # print(matrix.loc[i])
    #     print(matrix[i])
    rankedmatrices = [make_ranked_lists(matrix) for matrix in matrices]
    # print(rankedmatrices[0][0])
    # print(sz)
    res = [None for i in range(sz)]
    for s in range(sz):
        print(s)
        with open(savematrixfile, 'w') as f:
            f.write('\n'.join(['\t'.join(matrix[s]) for matrix in rankedmatrices]) + '\n')
        # print('ready ra')
        res[s] = run_ra(savematrixfile, sz)
    # print(res)

    result_df = reconsrtuct_df(res, genes)
    result_df.to_csv(saveresultfile, index=True, header=True, sep='\t')
    # print(result_df)
    return result_df


matricesdirname = "/home/user/Sirius/gene_network_sirius_2019/Matrices_1"
savematricesdirname = "/home/user/Sirius/gene_network_sirius_2019/Matrices_2"
predictedfilename = matricesdirname + "/{1}_{0}_predicted.txt"
truefilename = matricesdirname + "/{1}_{0}_true.txt"
savematricesfilename = savematricesdirname + "/{0}_predicted.txt"
datalist = ['exps_10', 'exps_10_2', 'exps_10_bgr', 'exps_50', 'exps_50_2', 'exps_50_bgr', 'exps_100', 'exps_100_2', 'exps_100_bgr', 'genes_200_exps_10_bgr', 'genes_400_exps_10_bgr', 'genes_600_exps_10_bgr', 'genes_700_exps_10_bgr', 'genes_1000_exps_10_bgr']
algolist = ['clr', 'aracne', 'mrnet', 'mrnetb']

tmpfile = "/home/user/Sirius/gene_network_sirius_2019/RankAggregation/data/tmp1.txt"


if __name__ == "__main__":
    results = np.zeros(shape=(len(datalist)))

    for i, dataname in enumerate(datalist):

        true_df = pd.read_csv(truefilename.format(dataname, algolist[1]), index_col=0, sep='\t')
        predicted_df = compute_aggregated_matrix(len(algolist), [predictedfilename.format(dataname, algo) for algo in algolist], tmpfile, savematricesfilename.format(dataname))
        true_df.to_csv(savematricesdirname + "/{0}_true.txt".format(dataname), index=True, header=True, sep='\t')
        # print(true_df)

        true_array = true_df.values[np.triu_indices(true_df.values.shape[0])]
        predicted_array = predicted_df.values[np.triu_indices(predicted_df.values.shape[0])]
        
        roc_auc = 0
        # try:
        #     fpr, tpr, thresholds = roc_curve(true_array, predicted_array)
        #     roc_auc = auc(fpr, tpr)
        # except:
        #     print("error", dataname, algo)
        fpr, tpr, thresholds = roc_curve(true_array, predicted_array)
        roc_auc = auc(fpr, tpr)
        results[i] = roc_auc
        print("done", dataname, results[i])
        
        # print("done", dataname, algo)

    print(results)


