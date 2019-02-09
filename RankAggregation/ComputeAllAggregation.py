# RUN WITH /usr/bin/python3 minet.py (python 3.6)

import sys
import numpy as np
from sklearn.metrics import roc_curve, auc
import pandas as pd


def compute_aggregated_matrix(matrixfiles_num, matrixfiles_dima, matrixfiles_petr, savematrixfile, saveresultfile):
    coeffs = [1 for i in range(matrixfiles_num)]
    # matrixfiles_num = int(sys.argv[1])
    # matrixfiles = [sys.argv[i] for i in range(2, matrixfiles_num + 2)]
    # savematrixfile = sys.argv[matrixfiles_num + 2]
    # saveresultfile = sys.argv[matrixfiles_num + 3]
    matrices = [pd.read_csv(f, index_col=0, sep=',').T.abs() for f in matrixfiles_dima] + [pd.read_csv(f, index_col=0, sep='\t').abs() for f in matrixfiles_petr]
    genes = matrices[0].index
    # print(genes)

    # print(matrices)
    sz = len(matrices[0])
    
    for matrix in matrices:
        try:
            assert len(matrix) == sz
        except:
            print(sz, len(matrix))
            # print(matrix)
            break


    for matrix in matrices:
        for column in matrix:
            temp = matrix[column].argsort()
            ranks = np.empty_like(temp)
            ranks[temp] = np.arange(len(matrix[column]))
            matrix[column] = ranks
    
    res = np.zeros(shape=(sz, sz))
    for s in range(sz):
        for i, matrix in enumerate(matrices):
            res[s] += matrix.iloc[:, s].values * coeffs[i]
        res[s] /= len(matrices)

    for row in res:
        row /= row.sum()

    result_df = pd.DataFrame(res, columns=genes, index=genes)
    
    result_df.to_csv(saveresultfile, index=True, header=True, sep='\t')
    # print(result_df)
    return result_df


matricesdirname = "/home/user/Sirius/dimasdata"
truematricesdirname = "/home/user/Sirius/gene_network_sirius_2019/Matrices_1"
savematricesdirname = "/home/user/Sirius/gene_network_sirius_2019/Matrices_5"
predictedfilename_dima = matricesdirname + "/{0}/{0}_{1}.csv"
predictedfilename_petr = truematricesdirname + "/{1}_{0}_predicted.txt"
truefilename = truematricesdirname + "/{1}_{0}_true.txt"
datalist = ['exps_10', 'exps_10_2', 'exps_10_bgr', 'exps_50', 'exps_50_2', 'exps_50_bgr', 'exps_100', 'exps_100_2', 'exps_100_bgr', 'genes_200_exps_10_bgr', 'genes_400_exps_10_bgr', 'genes_600_exps_10_bgr', 'genes_700_exps_10_bgr', 'genes_1000_exps_10_bgr']
algolist_dima = ['BR', 'ET', 'GB', 'Lasso_0.001', 'Lasso_0.0001', 'Lasso_0.00001', 'Lasso_0.000001', 'Lasso_0.0000001', 'RF', 'XGB']
algolist_petr = ['clr', 'aracne', 'mrnet', 'mrnetb']
saveresultsfile = "/home/user/Sirius/gene_network_sirius_2019/RankAggregation/res5.txt"
tmpfile = "/home/user/Sirius/dimasdata/tmp5.txt"
savematricesfilename = savematricesdirname + "/{0}_predicted.txt"


def matrix_stringify(m):
    """Transform top-triangle matrix to one-dimensional array
    Parameters
    ----------
    m: np.ndarray
        two-dimensional matrix to transform
    Returns
    -------
        one-dimensional np.array
    """
    arr = np.array([])
    for i in range(m.shape[0] - 1):
        arr = np.concatenate([arr, m[i, i+1:]])
    return arr


if __name__ == "__main__":
    results = np.zeros(shape=(len(datalist)))

    for i, dataname in enumerate(datalist):

        true_df = pd.read_csv(truefilename.format(dataname, 'aracne'), index_col=0, sep='\t')
        predicted_df = compute_aggregated_matrix(len(algolist_petr) + len(algolist_dima), [predictedfilename_dima.format(dataname, algo) for algo in algolist_dima], [predictedfilename_petr.format(dataname, algo) for algo in algolist_petr], tmpfile, savematricesfilename.format(dataname))
        true_df.to_csv(savematricesdirname + "/{0}_true.txt".format(dataname), index=True, header=True, sep='\t')
        # print(true_df)

        true_array = true_df.values[np.triu_indices(true_df.values.shape[0], k=1)]
        predicted_array = predicted_df.values[np.triu_indices(predicted_df.values.shape[0], k=1)]
        
        roc_auc = 0
        # try:
        #     fpr, tpr, thresholds = roc_curve(true_array, predicted_array)
        #     roc_auc = auc(fpr, tpr)
        # except:
        #     print("error", dataname, algo)
        fpr, tpr, thresholds = roc_curve(true_array, predicted_array)
        roc_auc = auc(fpr, tpr)
        results[i] = roc_auc

        with open(savematricesdirname + "/{0}_auc.txt".format(dataname), 'w') as f:
            f.write(str(roc_auc) + '\n')
        print("done", dataname, results[i])
        with open(saveresultsfile, "a") as f:
            f.write("done " + dataname + str(results[i]) + '\n')
        
        # print("done", dataname, algo)

    print(results)
