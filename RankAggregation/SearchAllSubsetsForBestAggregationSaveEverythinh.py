import sys
import numpy as np
from sklearn.metrics import roc_curve, auc
import pandas as pd
from itertools import combinations


def compute_aggregated_matrix(matrixfiles_num, matrixfiles_dima, matrixfiles_petr, matrixfiles_sergey, saveresultfile):
    if matrixfiles_num == 0:
        return pd.DataFrame([0])
    coeffs = [1 for i in range(matrixfiles_num)]
    # matrixfiles_num = int(sys.argv[1])
    # matrixfiles = [sys.argv[i] for i in range(2, matrixfiles_num + 2)]
    # savematrixfile = sys.argv[matrixfiles_num + 2]
    # saveresultfile = sys.argv[matrixfiles_num + 3]
    matrices = [pd.read_csv(f, index_col=0, sep=',').T.abs() for f in matrixfiles_dima] + \
               [pd.read_csv(f, index_col=0, sep='\t').abs() for f in matrixfiles_petr] + \
               [pd.read_csv(f, sep=',').abs() for f in matrixfiles_sergey]
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


# algolist_dima = ['BR', 'ET', 'GB', 'Lasso_0.001', 'Lasso_0.0001', 'Lasso_0.00001', 'Lasso_0.000001', 'Lasso_0.0000001', 'RF', 'XGB']
algolist_dima = ['ET', 'GB', 'Lasso_0.000001', 'RF', 'XGB']
algolist_petr = ['aracne', 'mrnet', 'mrnetb']
algolist_sergey = ['']

datafolder_dima = "/home/user/Sirius/dimasdata"
datafolder_petr = "/home/user/Sirius/gene_network_sirius_2019/Matrices_1"
datafolder_sergey = "/home/user/Sirius/sergeysdata"
truedatafolder = "/home/user/Sirius/gene_network_sirius_2019/Matrices_1"

datalist = ['exps_10', 'exps_10_2', 'exps_10_bgr', 'exps_50', 'exps_50_2', 'exps_50_bgr', 'exps_100', 'exps_100_2', 'exps_100_bgr', 'genes_200_exps_10_bgr', 'genes_200_exps_20_bgr', 'genes_200_exps_40_bgr', 'genes_400_exps_10_bgr', 'genes_500_exps_10_bgr', 'genes_500_exps_50_bgr', 'genes_600_exps_10_bgr', 'genes_700_exps_10_bgr', 'genes_1000_exps_10_bgr']

truefilename = truedatafolder + "/{1}_{0}_true.txt"
dataname_dima = datafolder_dima + "/{0}/{0}_{1}.csv"
dataname_petr = datafolder_petr + "/{1}_{0}_predicted.txt"
dataname_sergey = datafolder_sergey + "/{0}_{1}.csv"

tmpfile = "/home/user/Sirius/gene_network_sirius_2019/RankAggregation/tmp_Subsets_1.txt"
savematricesdirname = "/home/user/Sirius/Matrices_Subsets_2"
savematricesfilename = savematricesdirname + "/{0}_predicted.txt"
saveresultsfile = "/home/user/Sirius/gene_network_sirius_2019/RankAggregation/res_aggregated_full1.txt"


if __name__ == "__main__":
    results = []

    for d in range(len(algolist_dima) + 1):
        for d_t in combinations(algolist_dima, d):
            for p in range(len(algolist_petr) + 1):
                for p_t in combinations(algolist_petr, p):
                    for s in range(len(algolist_sergey) + 1):
                        for s_t in combinations(algolist_sergey, s):
                            if d + p + s == 0:
                                continue
                            auc_s = 0.0
                            auc_b = 0.0
                            snum = len([a for a in datalist if a.startswith("exps")])
                            bnum = len(datalist) - snum
                            aucs = []
                            for i, dataname in enumerate(datalist):
                                true_df = pd.read_csv(truefilename.format(dataname, 'aracne'), index_col=0, sep='\t')
                                predicted_df = compute_aggregated_matrix(len(d_t) + len(p_t) + len(s_t), 
                                    [dataname_dima.format(dataname, algo) for algo in d_t], 
                                    [dataname_petr.format(dataname, algo) for algo in p_t], 
                                    [dataname_sergey.format(dataname, algo) for algo in s_t], 
                                    savematricesfilename.format(dataname + "_" + str(d_t) + "_" + str(p_t) + "_" + str(s_t)))

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
                                aucs.append(roc_auc)
                                if dataname.startswith("exps"):
                                    auc_s += roc_auc
                                else:
                                    auc_b += roc_auc
                            auc_s /= snum
                            auc_b /= bnum
                            results.append((auc_b, auc_s, d_t, p_t, s_t, aucs))
                            print("done", results[len(results) - 1])
    results.sort()
    with open(saveresultsfile, 'a') as f:
        for t in results:
            f.write(str(t[0]) + " " + str(t[1]) + " " + str(t[2]) + " " + str(t[3]) + " " + str(t[4]) + "\n")
