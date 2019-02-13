import sys
import numpy as np
from sklearn.metrics import roc_curve, auc
import pandas as pd
from itertools import combinations
import matplotlib.pyplot as plt


def compute_matrix(matrixfiles_num, matrixfiles_sergey):
    if matrixfiles_num == 0:
        return pd.DataFrame([0])
    coeffs = [1 for i in range(matrixfiles_num)]
    # matrixfiles_num = int(sys.argv[1])
    # matrixfiles = [sys.argv[i] for i in range(2, matrixfiles_num + 2)]
    # savematrixfile = sys.argv[matrixfiles_num + 2]
    # saveresultfile = sys.argv[matrixfiles_num + 3]
    matrices = [pd.read_csv(f, sep=',').abs() for f in matrixfiles_sergey]
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


    # for matrix in matrices:
    #     for column in matrix:
    #         temp = matrix[column].argsort()
    #         ranks = np.empty_like(temp)
    #         ranks[temp] = np.arange(len(matrix[column]))
    #         matrix[column] = ranks
    
    res = np.zeros(shape=(sz, sz))
    for s in range(sz):
        for i, matrix in enumerate(matrices):
            res[s] += matrix.iloc[:, s].values * coeffs[i]
        res[s] /= len(matrices)

    for row in res:
        if row.sum() != 0:
            row /= row.sum() 

    result_df = pd.DataFrame(res, columns=genes, index=genes)
    
    # result_df.to_csv(saveresultfile, index=True, header=True, sep='\t')
    # print(result_df)
    return result_df


def compute_aggregated_matrix(matrixfiles_num, matrixfiles_sergey):
    if matrixfiles_num == 0:
        return pd.DataFrame([0])
    coeffs = [1 for i in range(matrixfiles_num)]
    # matrixfiles_num = int(sys.argv[1])
    # matrixfiles = [sys.argv[i] for i in range(2, matrixfiles_num + 2)]
    # savematrixfile = sys.argv[matrixfiles_num + 2]
    # saveresultfile = sys.argv[matrixfiles_num + 3]
    matrices = [pd.read_csv(f, sep=',').abs() for f in matrixfiles_sergey]
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
        if row.sum() != 0:
            row /= row.sum() 

    result_df = pd.DataFrame(res, columns=genes, index=genes)
    
    # result_df.to_csv(saveresultfile, index=True, header=True, sep='\t')
    # print(result_df)
    return result_df


algolist_sergey = ['']

datafolder_sergey = "/home/user/Sirius/sergeysdata"
truedatafolder = "/home/user/Sirius/gene_network_sirius_2019/Matrices_1"

# datalist = ['genes_200_exps_10_bgr', 'genes_200_exps_20_bgr', 'genes_200_exps_40_bgr']
datalist = ['exps_100_2']

truefilename = truedatafolder + "/{1}_{0}_true.txt"
dataname_sergey = datafolder_sergey + "/{0}_{1}.csv"

# savematricesdirname = "/home/user/Sirius/Matrices_Subsets"
# savematricesfilename = savematricesdirname + "/{0}_predicted.txt"
saveresultsfile = "/home/user/Sirius/gene_network_sirius_2019/Graphics/result.txt"


if __name__ == "__main__":
    results = []

    
    for algo in algolist_sergey:
        for i, dataname in enumerate(datalist):
            true_df = pd.read_csv(truefilename.format(dataname, 'aracne'), index_col=0, sep='\t')
            predicted_df = compute_matrix(1, [dataname_sergey.format(dataname, algo)])
            predicted_agg_df = compute_aggregated_matrix(1, [dataname_sergey.format(dataname, algo)])

            true_array = true_df.values[np.triu_indices(true_df.values.shape[0], k=1)]
            predicted_array = predicted_df.values[np.triu_indices(predicted_df.values.shape[0], k=1)]
            predicted_agg_array = predicted_agg_df.values[np.triu_indices(predicted_agg_df.values.shape[0], k=1)]
            
            fpr, tpr, thresholds = roc_curve(true_array, predicted_array)
            roc_auc = auc(fpr, tpr)
            print("plotting...", algo, dataname)
            plt.plot(fpr, tpr)
            plt.show()
            fpr1, tpr1, thresholds1 = roc_curve(true_array, predicted_agg_array)
            roc_auc1 = auc(fpr1, tpr1)
            print("plotting...", algo, dataname)
            plt.plot(fpr1, tpr1)
            plt.show()
            print("done", algo, dataname, roc_auc, roc_auc1)
            