# RUN WITH /usr/bin/python3 minet.py (python 3.6)

import sys
import numpy as np
from sklearn.metrics import roc_curve, auc
import pandas as pd


matricesdirname = "/home/user/Sirius/dimasdata"
truematricesdirname = "/home/user/Sirius/gene_network_sirius_2019/Matrices_1"
predictedfilename = matricesdirname + "/{0}/{0}_{1}.csv"
truefilename = truematricesdirname + "/{1}_{0}_true.txt"
datalist = ['exps_10', 'exps_10_2', 'exps_10_bgr', 'exps_50', 'exps_50_2', 'exps_50_bgr', 'exps_100', 'exps_100_2', 'exps_100_bgr', 'genes_200_exps_10_bgr', 'genes_400_exps_10_bgr', 'genes_600_exps_10_bgr', 'genes_700_exps_10_bgr', 'genes_1000_exps_10_bgr']
# algolist = ['BR', 'ET', 'GB', 'Lasso_0.001', 'Lasso_0.0001', 'Lasso_0.00001', 'Lasso_0.000001', 'Lasso_0.0000001', 'RF', 'XGB']
algolist = ['BR', 'ET', 'GB', 'Lasso_0.001', 'Lasso_0.0001', 'Lasso_0.00001', 'Lasso_0.000001', 'Lasso_0.0000001', 'RF', 'XGB']
saveresultsfile = "/home/user/Sirius/dimasdata/res2.txt"



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
    results = np.zeros(shape=(len(datalist), len(algolist)))

    for i, dataname in enumerate(datalist):
        for j, algo in enumerate(algolist):

            true_df = pd.read_csv(truefilename.format(dataname, 'aracne'), index_col=0, sep='\t')
            predicted_df = pd.read_csv(predictedfilename.format(dataname, algo), index_col=0, sep=',').abs().fillna(0)
            # print(predicted_df)
            # print(1/0)
            # true_df.to_csv(savematricesdirname + "/{0}_true.txt".format(dataname), index=True, header=True, sep='\t')
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
            fpr1, tpr1, thresholds1 = roc_curve(matrix_stringify(true_df.values), matrix_stringify(predicted_df.values))
            roc_auc1 = auc(fpr1, tpr1)

            print("ME AND DIMA", roc_auc, roc_auc1)
            results[i][j] = roc_auc

            print("done", dataname, algo, results[i][j])
            with open(saveresultsfile, "a") as f:
                f.write("done " + dataname + " " + algo + " " + str(results[i][j]) + '\n')
            
            # print("done", dataname, algo)

    print(results)


