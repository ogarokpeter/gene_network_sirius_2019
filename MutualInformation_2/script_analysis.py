# Run: python3 scripy_analysis.py file_with_true_data file_with_predicted_data file_to_save_results


import sys
import numpy as np
from sklearn.metrics import roc_curve, auc


def read_adjacency_matrix(filename):
    matrix = np.loadtxt(filename)
    return matrix


def adjacency_matrix_to_array_of_degrees(adjacency_matrix):
    degrees = [[] for i in range(len(adjacency_matrix))]
    for i in range(len(adjacency_matrix)):
        deg_i = 0
        for j in range(len(adjacency_matrix[i])):
            if adjacency_matrix[i][j] == 1:
                deg_i += 1
        degrees[deg_i].append(i)
    return degrees


if __name__ == "__main__":
    truefilename = sys.argv[1]
    filename = sys.argv[2]
    resultfilename = sys.argv[3]


    # matrix = read_adjacency_matrix(filename)
    matrix = np.array([[0.08737, 0.7681], [0.234142, 0.123234]])
    # true_matrix = read_adjacency_matrix(truefilename)
    true_matrix = np.array([[0., 1.], [0., 0.]])

    degrees = adjacency_matrix_to_array_of_degrees(true_matrix)

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
    with open(resultfilename, "a") as f:
        print(truefilename, filename)
        f.write(truefilename + " " + filename + ":" + '\n')
        for i in range(len(aucs_mean)):
            if aucs_mean[i] != 0:
                f.write(str(i) + " " + str(aucs_mean[i]) + " " + str(aucs_var[i]) + '\n')
