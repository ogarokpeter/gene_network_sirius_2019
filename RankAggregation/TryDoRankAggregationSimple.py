# RUN WITH /usr/bin/python3 minet.py (python 3.6)

import sys
import numpy as np
from sklearn.metrics import roc_curve, auc
import pandas as pd


def compute_aggregated_matrix(matrixfiles_num, matrixfiles, savematrixfile, saveresultfile, coeffs=[1, 1, 1, 1]):
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


matricesdirname = "/home/user/Sirius/gene_network_sirius_2019/Matrices_1"
savematricesdirname = "/home/user/Sirius/gene_network_sirius_2019/Matrices_4"
predictedfilename = matricesdirname + "/{1}_{0}_predicted.txt"
truefilename = matricesdirname + "/{1}_{0}_true.txt"
savematricesfilename = savematricesdirname + "/{0}_predicted.txt"
datalist = ['exps_10', 'exps_10_2', 'exps_10_bgr', 'exps_50', 'exps_50_2', 'exps_50_bgr', 'exps_100', 'exps_100_2', 'exps_100_bgr', 'genes_200_exps_10_bgr', 'genes_400_exps_10_bgr', 'genes_600_exps_10_bgr', 'genes_700_exps_10_bgr', 'genes_1000_exps_10_bgr']
algolist = ['clr', 'aracne', 'mrnet', 'mrnetb']
saveresultsfile = "/home/user/Sirius/gene_network_sirius_2019/RankAggregation/res2.txt"
tmpfile = "/home/user/Sirius/gene_network_sirius_2019/RankAggregation/data/tmp4.txt"



if __name__ == "__main__":
    results = np.zeros(shape=(len(datalist)))

    for i, dataname in enumerate(datalist):

        true_df = pd.read_csv(truefilename.format(dataname, algolist[1]), index_col=0, sep='\t')
        predicted_df = compute_aggregated_matrix(len(algolist), [predictedfilename.format(dataname, algo) for algo in algolist], tmpfile, savematricesfilename.format(dataname))
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
            f.write("done " + dataname + str(results[i]))
        
        # print("done", dataname, algo)

    print(results)


